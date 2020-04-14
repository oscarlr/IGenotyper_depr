#!/bin/env python
import os
import sys # temporary
import pysam
import subprocess
import pybedtools
import networkx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import namedtuple
from networkx.algorithms.components.connected import connected_components

from ..common import *
from ..command_line import *

def not_overlapping(alignment):
    if not (int(alignment.qstart) < 500 and (int(alignment.slen) - int(alignment.send)) < 500):
        if not (int(alignment.sstart) < 500 and (int(alignment.qlen) - int(alignment.qend)) < 500):
            return True
    return False

def get_regions_with_hap(bed_file,hap):
    regions = []
    for interval in bed_file:
        chrom = show_value(interval.chrom)
        start = show_value(interval.start)
        end = show_value(interval.end)
        regions.append([chrom,start,end,hap])
    return regions

def get_assemble_regions(phased_regions,regions):
    phased_regions_bed = pybedtools.BedTool(phased_regions)
    regions_bed = pybedtools.BedTool(regions)
    phased_regions_bed = regions_bed.intersect(phased_regions_bed)
    non_phased_regions_bed = regions_bed.subtract(phased_regions_bed)
    regions = get_regions_with_hap(phased_regions_bed,"1")
    regions += get_regions_with_hap(phased_regions_bed,"2")
    regions += get_regions_with_hap(non_phased_regions_bed,"0")
    return regions

def merge_small_regions(regions,length=1000):
    new_regions = []
    regions = sorted(regions,key=lambda x: x[1])    
    for i,(chrom,start,end,hap) in enumerate(regions):
        new_region = [chrom,start,end,hap]
        if (end - start) < length:
            for chrom2,start2,end2,hap2 in regions[i:]:
                if hap != hap2:
                    continue
                if end == start2:
                    new_region = [chrom,start,end2,hap]
                    break
        new_regions.append(new_region)
    return new_regions            

def show_value(s):
    if sys.version_info.major == 2:
        if isinstance(s, unicode):
            return str(s)
    return s

def distance_between_regions(q_assembled_region,s_assembled_region):
    q_start = int(q_assembled_region[1])
    q_end = int(q_assembled_region[2])
    s_start = int(s_assembled_region[1])
    s_end = int(s_assembled_region[2])
    distance = max(0,max(q_start, s_start) - min(q_end, s_end))
    return distance

def pass_filters(alignment,hap,pairs,errors=5,length=5000):
    pident_min = 95.0
    mismatch_max = errors
    length_min = length 
    q_assembled_region = assembly_location(alignment.qseqid)
    q_haplotype = get_haplotype(alignment.qseqid)
    s_assembled_region = assembly_location(alignment.sseqid)
    s_haplotype = get_haplotype(alignment.sseqid)
    if distance_between_regions(q_assembled_region,s_assembled_region) > 100:
        return False
    if int(alignment.length) < length_min:
        return False
    if alignment.qseqid == alignment.sseqid:
        return False
    if float(alignment.pident) < pident_min:
        return False
    if float(alignment.mismatch) > mismatch_max:
        return False
    if alignment.sstrand != "plus":
        return False
    if not_overlapping(alignment):
        return False
    if int(alignment.sstart) > int(alignment.qstart):
        return False
    if sorted((alignment.qseqid,alignment.sseqid)) in pairs:
        return False
    return True

def to_graph(l):
    G = networkx.Graph()
    for part in l:
        G.add_nodes_from(part)
        G.add_edges_from(to_edges(part))
    return G

def to_edges(l):
    it = iter(l)
    last = next(it)
    for current in it:
        yield last, current
        last = current

def group_alignments(alignments,contig_names):
    contig_groupings = []
    for alignment in alignments:
        if alignment.qseqid in contig_names:
            contig_names.remove(alignment.qseqid)
        if alignment.sseqid in contig_names:
            contig_names.remove(alignment.sseqid)
        contig_groupings.append([alignment.qseqid,alignment.sseqid])
    G = to_graph(contig_groupings)
    groupings = {}
    i = 0
    for group in connected_components(G):
        groupings[i] = group
        i += 1
    if len(groupings.keys()) != 0:
        max_group = max(groupings.keys())
    else:
        max_group = 0
    for i,contig in enumerate(contig_names):
        groupings[i + max_group + 1] = [contig]
    return groupings

def create_directed_graph(alignments,contigs):
    G = networkx.DiGraph()
    for alignment in alignments:
        if alignment.qseqid not in contigs:
            continue
        if alignment.sseqid not in contigs:
            continue
        weight = int(alignment.length) + (int(alignment.slen) - int(alignment.length)) + (int(alignment.qlen) - int(alignment.length))
        G.add_edge(alignment.qseqid, alignment.sseqid, weight=weight)
    return G

def completely_overlapping(alignment):
    q_start = int(alignment.qstart) == 1
    q_end = (int(alignment.qlen) - int(alignment.qend)) == 0
    s_start = int(alignment.sstart) == 1
    s_end = (int(alignment.slen) - int(alignment.send)) == 0
    if (q_start and q_end) or (s_start and s_end):
        return True
    return False

def group_size_2_and_encaps(group,alignments):
    if len(group) == 2:
        for alignment in alignments:
            if (alignment.qseqid == group[0] and alignment.sseqid == group[1]) or \
               (alignment.qseqid == group[1] and alignment.sseqid == group[0]):
                if completely_overlapping(alignment):
                    return True
    return False

def get_contig_coords(alignments,max_path):
    if len(max_path) == 1:
        return [[max_path[0],1,-1]]
    coords = []
    i = 0
    for start_contig, end_contig in zip(max_path,max_path[1:]):
        for alignment in alignments:
            if alignment.qseqid != start_contig:
                continue
            if alignment.sseqid != end_contig:
                continue
            merge_coords = [alignment.qseqid,1,int(alignment.qstart) - 1]
            coords.append(merge_coords)
            i += 1
            if i == len(max_path) - 1:
                merge_coords = [alignment.sseqid,1,alignment.slen]
                coords.append(merge_coords)
    return coords

def valid_path(contig_path):
    previous_hap = None
    current_hap = None
    switch = False
    valid = True
    for contig in contig_path:
        hap = get_haplotype(contig)
        if current_hap == None:
            current_hap = hap
            continue
        previous_hap = current_hap
        current_hap = hap
        if not switch:
            if current_hap != previous_hap:
                switch = True
        else:            
            valid = False
            return valid
    return valid

def get_larger_contig(group,alignments):
    larger_contig = None
    for alignment in alignments:
        if (alignment.qseqid == group[0] and alignment.sseqid == group[1]) or \
           (alignment.qseqid == group[1] and alignment.sseqid == group[0]):
            if int(alignment.qlen) > int(alignment.slen):
                larger_contig = alignment.qseqid
            else:
                larger_contig = alignment.sseqid
            break
    assert larger_contig != None
    return larger_contig

def group_merging_contigs(alignments,contig_names,hap):
    groupings = {}
    groups = group_alignments(alignments,contig_names)
    for group in groups:
        G = create_directed_graph(alignments,groups[group])
        all_paths_weighted = dict(networkx.all_pairs_dijkstra_path_length(G))
        all_paths =  networkx.shortest_path(G)
        max_path = None
        num_contigs_in_max_path = 0
        for start in all_paths_weighted:
            for end in all_paths_weighted[start]:                
                # if not valid_path(all_paths[start][end]):
                #     continue
                length = all_paths_weighted[start][end]
                if length > num_contigs_in_max_path:
                    num_contigs_in_max_path = length
                    max_path = all_paths[start][end]
        if group_size_2_and_encaps(list(groups[group]),alignments):
            max_path = None
            groups[group] = [get_larger_contig(list(groups[group]),alignments)]
        if max_path == None:
            max_path = groups[group]
        contigs_coordinates = get_contig_coords(alignments,max_path)
        groupings[group] = contigs_coordinates
    return groupings

def get_merged_sequence_name(names,group,in_hap):
    starts = []
    ends = []
    haps = []
    for name in names:
        name_loc = assembly_location(name)
        starts.append(int(name_loc[1]))
        ends.append(int(name_loc[2]))
        haps.append(get_haplotype(name))
    if "1" in haps:
        hap = "1"
        assert "2" not in haps
    elif "2" in haps:
        hap = "2"
        assert "1" not in haps
    else:
        hap = "0"
    merged_name="c=igh:%s-%s_h=%s_g=%s_ih=%s_/0/0_0" % (min(starts),max(ends),hap,group,in_hap)
    return merged_name

class AssemblyRun():
    def __init__(self,Sample):
        self.sample = Sample
        self.command_line_tools = CommandLine(Sample)

    def load_whatshap_blocks(self,min_length=500,min_variants=2):
        if not non_emptyfile(self.sample.haplotype_blocks):
            self.command_line_tools.get_phased_blocks()
        blocks = []
        Block = namedtuple('Block',['sample','chrom','start_1',
                                    'start','end','num_variants'])
        with open(self.sample.haplotype_blocks,'r') as fh:
            header = fh.readline()
            for line in fh:
                line = line.rstrip().split('\t')
                block = Block._make(line)
                if int(block.num_variants) < min_variants:
                    continue
                if (int(block.end) - int(block.start)) < min_length:
                    continue
                blocks.append([block.chrom,int(block.start),int(block.end)])
        return sorted(blocks,key=lambda x: x[1])

    def get_phased_regions_to_assemble(self):
        whatshap_blocks = pybedtools.BedTool(self.load_whatshap_blocks())
        if not self.sample.dont_split:
            sv_regions = load_bed_regions(self.sample.sv_regions)
            non_sv_regions = load_bed_regions(self.sample.non_sv_regions)
            sv_regions_to_assemble = get_assemble_regions(whatshap_blocks,sv_regions)
            non_sv_regions_to_assemble = get_assemble_regions(whatshap_blocks,non_sv_regions)
            regions = sv_regions_to_assemble + non_sv_regions_to_assemble
        else:            
            non_dup_regions = load_bed_regions(self.sample.non_dup_regions)
            dup_regions = load_bed_regions(self.sample.dup_regions)
            non_dup_regions_assemble = get_assemble_regions(whatshap_blocks,non_dup_regions)
            dup_regions_assemble = get_assemble_regions(whatshap_blocks,dup_regions)
            regions = non_dup_regions_assemble + dup_regions_assemble
            #igh_region = load_bed_regions(self.sample.igh_coords)
            #regions = get_assemble_regions(whatshap_blocks,igh_region)
        regions = merge_small_regions(regions)        
        pybedtools.BedTool(regions).saveas(self.sample.regions_to_assemble)

    def create_assembly_scripts(self,flank=1000):     
        min_coverage = 10
        regions_to_assemble = pybedtools.BedTool(self.sample.regions_to_assemble)
        bashfiles = []
        for interval in regions_to_assemble:
            chrom = show_value(interval.chrom)
            start = show_value(interval.start)
            end = show_value(interval.end)
            hap = show_value(interval.name)
            if get_phased_coverage(self.sample.phased_subreads_mapped_reads,chrom,start,end,hap) < min_coverage:
                continue
            # if self.sample.add_unphased_reads and hap != "0":
            #     samtools_hap = "-r %s -r 0" % hap            
            # else:
            samtools_hap = "-r %s" % hap            
            directory = "%s/assembly/%s/%s_%s/%s" % (self.sample.outdir,chrom,start,end,hap)
            create_directory(directory)
            if os.path.isfile("%s/done" % directory):
                continue
            bashfile = "%s/assemble.sh" % directory
            params = {
                "hap": samtools_hap,
                "ccs_to_ref": self.sample.phased_ccs_mapped_reads,
                "chrom": chrom,
                "start": max(0,int(start) - flank),
                "end": int(end) + flank,
                "output": directory,
                "threads": self.sample.threads,
                "size": (int(end) - int(start)) + flank*2,
                "subreads": self.sample.input_bam,
                "subreads_to_ref": self.sample.phased_subreads_mapped_reads,
                "python_scripts": self.sample.python_scripts,
                "ref": self.sample.blasr_ref
                }
            write_to_bashfile(self.sample.assembly_script,bashfile,params)
            bashfiles.append(bashfile)
        return bashfiles

    def combine_sequence(self,type_,outfile,file_name="merged_contigs_quivered"):
        seqs = []
        bedfh = open(self.sample.regions_to_assemble, 'r')
        for line in bedfh:
            line = line.rstrip().split('\t')
            chrom,start,end,hap = line
            directory = "%s/assembly/%s/%s_%s/%s" % (self.sample.outdir,chrom,start,end,hap)
            contig = "%s/%s.%s" % (directory,file_name,type_)
            if os.path.isfile(contig):
                contigs = list(SeqIO.parse(contig,type_))
                total_contigs = len(contigs)
                for i,record in enumerate(contigs):
                    record.id = "c=%s:%s-%s_h=%s_i=%s_t=%s_/0/0_0" % (chrom,start,end,hap,i,total_contigs)
                    record.description = ""
                    seqs.append(record)
        SeqIO.write(seqs,outfile,type_)
        bedfh.close()
        
    def combine_sequences(self):
        self.combine_sequence("fastq",self.sample.locus_fastq)
        self.combine_sequence("fasta",self.sample.locus_fasta)
        self.combine_sequence("fasta",self.sample.locus_fasta_unquivered,file_name="merged_contigs")
        

    def filter_alignments_hap_specific(self,hap):
        pairs = []
        alignments = []
        columns = ["length","pident","nident","mismatch","gapopen","gaps","qseqid",
                   "qstart","qend","qlen","sseqid","sstart","send","slen","sstrand"]
        contig_names = set()
        Alignment = namedtuple('Alignment',columns)
        with open(self.sample.contigs_to_contigs_blast,'r') as fh:
            for line in fh:
                line = line.rstrip().split('\t')
                alignment = Alignment._make(line)
                q_haplotype = get_haplotype(alignment.qseqid)
                s_haplotype = get_haplotype(alignment.sseqid)
                hap_to_ignore = None
                if hap == "1":
                    hap_to_ignore = "2"
                else:
                    hap_to_ignore = "1"
                if q_haplotype == hap_to_ignore:
                    continue
                if s_haplotype == hap_to_ignore:
                    continue
                contig_names.add(alignment.qseqid)
                contig_names.add(alignment.sseqid)
                if pass_filters(alignment,hap,pairs):
                    pairs.append(sorted((alignment.qseqid,alignment.sseqid)))
                    alignments.append(alignment)
        return alignments,contig_names

    def get_hap_merges(self):
        hap_merges = {}
        for hap in ["1","2"]:
            hap_alignments,hap_contig_names = self.filter_alignments_hap_specific(hap)
            hap_groupings = group_alignments(hap_alignments,hap_contig_names)
            hap_groupings_with_contigs_to_merge = group_merging_contigs(hap_alignments,hap_contig_names,hap)
            hap_merges[hap] = [hap_alignments,hap_contig_names,hap_groupings,hap_groupings_with_contigs_to_merge]
        return hap_merges

    def save_alignments(self,alignments):
        with open(self.sample.contigs_to_contigs_blast_edited,'w') as fh:
            for alignment in alignments:
                fh.write("%s\n" % "\t".join(map(str,alignment)))

    def read_merging_instructions(self):
        infasta = SeqIO.to_dict(SeqIO.parse(self.sample.locus_fasta,"fasta"))
        outsequences = {}
        with open(self.sample.merge_alignments_instructions,'r') as fh:
            columns = ["group","hap","contig","start","end"]
            Line = namedtuple('Line',columns)
            for line in fh:
                line = line.rstrip().split('\t')
                line = Line._make(line)
                if line.hap not in outsequences:
                    outsequences[line.hap] = {}
                if line.group not in outsequences[line.hap]:
                    outsequences[line.hap][line.group] = {}
                    outsequences[line.hap][line.group]["seq"] = []
                    outsequences[line.hap][line.group]["names"] = []
                outsequences[line.hap][line.group]["seq"].append(str(infasta[line.contig].seq[int(line.start):int(line.end)]))
                outsequences[line.hap][line.group]["names"].append(line.contig)
        return outsequences

    def write_sequences(self):
        records = []
        hap0_sequences_seen = []
        outsequences = self.read_merging_instructions()
        for hap in outsequences:
            for group in outsequences[hap]:
                sequence = "".join(outsequences[hap][group]["seq"])
                if sequence in hap0_sequences_seen:
                    continue
                sequence_name = get_merged_sequence_name(outsequences[hap][group]["names"],group,hap)
                record = SeqRecord(Seq(sequence,"fasta"),id=sequence_name,name="",description="")
                if "h=0" in sequence_name:
                    hap0_sequences_seen.append(sequence)
                records.append(record)
        SeqIO.write(records,self.sample.merged_contigs,"fasta")
        
    def merge_sequences(self):
        self.command_line_tools.run_blast(self.sample.locus_fasta,self.sample.contigs_to_contigs_blast)
        hap_merges = self.get_hap_merges()
        alignments = set(hap_merges["1"][0]).union(set(hap_merges["2"][0]))
        self.save_alignments(alignments)
        with open(self.sample.contigs_grouped,'w') as fh:
            for hap in hap_merges:
                hap_group = hap_merges[hap][2]
                for group in hap_group:
                    for contig_name in hap_group[group]:
                        fh.write("%s\t%s\t%s\n" % (group,contig_name,hap))
        with open(self.sample.merge_alignments_instructions,'w') as fh:
            for hap in hap_merges:
                hap_group = hap_merges[hap][3]
                for group in hap_group:
                    for contig_coord in hap_group[group]:
                        output = [group,hap] + contig_coord
                        fh.write("%s\n" % "\t".join(map(str,output)))
        self.write_sequences()
        
    def __call__(self):
        self.get_phased_regions_to_assemble()
        assembly_scripts = self.create_assembly_scripts()
        run_assembly_scripts(assembly_scripts,self.sample.cluster,self.sample.cluster_walltime,
                             self.sample.threads,self.sample.cluster_mem,self.sample.cluster_queue)
        self.combine_sequences()
        self.command_line_tools.map_locus()
        self.command_line_tools.map_unquivered_locus()
        self.merge_sequences()
        self.command_line_tools.map_merged_locus()

def assemble_reads(self):
    assembly_runner = AssemblyRun(self)
    assembly_runner()

    
