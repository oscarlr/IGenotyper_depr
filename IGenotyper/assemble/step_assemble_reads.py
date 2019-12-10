#!/bin/env python
import os
import sys # temporary
import pysam
import subprocess
import pybedtools
from Bio import SeqIO
from collections import namedtuple

from ..common import *
from ..command_line import *

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

def pass_filters(hap,errors=5,length=5000):
    pident_min = 95.0
    mismatch_max = errors
    length_min = length 
    q_assembled_region = assembly_location(alignment.qseqid)
    q_haplotype = get_haplotype(alignment.qseqid)
    s_assembled_region = assembly_location(alignment.sseqid)
    s_haplotype = get_haplotype(alignment.sseqid)
    hap_to_ignore = None
    if hap == "1":
        hap_to_ignore = "2"
    else:
        hap_to_ignore = "1"
    if q_haplotype == hap_to_ignore:
        return False
    if s_haplotype == hap_to_ignore:
        return False
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

def group_alignments(alignments,fosmid_names):
    fosmid_groupings = []
    for alignment in alignments:
        if alignment.qseqid in fosmid_names:
            fosmid_names.remove(alignment.qseqid)
        if alignment.sseqid in fosmid_names:
            fosmid_names.remove(alignment.sseqid)
        fosmid_groupings.append([alignment.qseqid,alignment.sseqid])
    G = to_graph(fosmid_groupings)
    groupings = {}
    i = 0
    for group in connected_components(G):
        groupings[i] = group
        i += 1
    if len(groupings.keys()) != 0:
        max_group = max(groupings.keys())
    else:
        max_group = 0
    for i,fosmid in enumerate(fosmid_names):
        groupings[i + max_group + 1] = [fosmid]
    return groupings

def create_directed_graph(alignments,fosmids):
    G = networkx.DiGraph()
    for alignment in alignments:
        if alignment.qseqid not in fosmids:
            continue
        if alignment.sseqid not in fosmids:
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

def get_fosmid_coords(alignments,max_path):
    if len(max_path) == 1:
        return [[max_path[0],1,-1]]
    coords = []
    i = 0
    for start_fosmid, end_fosmid in zip(max_path,max_path[1:]):
        for alignment in alignments:
            if alignment.qseqid != start_fosmid:
                continue
            if alignment.sseqid != end_fosmid:
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
    for contig in contig_paths:
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

def group_merging_fosmids(alignments,fosmid_names,hap):
    groupings = {}
    groups = group_alignments(alignments,fosmid_names)
    for group in groups:
        G = create_directed_graph(alignments,groups[group])
        all_paths_weighted = dict(networkx.all_pairs_dijkstra_path_length(G))
        all_paths =  networkx.shortest_path(G)
        max_path = None
        num_fosmids_in_max_path = 0
        for start in all_paths_weighted:
            for end in all_paths_weighted[start]:                
                if not valid_path(all_paths[start][end]):
                    continue
                length = all_paths_weighted[start][end]
                if length > num_fosmids_in_max_path:
                    num_fosmids_in_max_path = length
                    max_path = all_paths[start][end]
        if group_size_2_and_encaps(list(groups[group]),alignments):
            max_path = None
            groups[group] = [get_larger_contig(list(groups[group]),alignments)]
        if max_path == None:
            max_path = groups[group]
        fosmids_coordinates = get_fosmid_coords(alignments,max_path)
        groupings[group] = fosmids_coordinates
    return groupings

class AssemblyRun():
    def __init__(self,Sample):
        self.sample = Sample
        self.command_line_tools = CommandLine(self)

    def load_whatshap_blocks(self,min_length=500,min_variants=2):
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
        sv_regions = load_bed_regions(self.sample.sv_regions)
        non_sv_regions = load_bed_regions(self.sample.non_sv_regions)
        sv_regions_to_assemble = get_assemble_regions(whatshap_blocks,sv_regions)
        non_sv_regions_to_assemble = get_assemble_regions(whatshap_blocks,non_sv_regions)
        regions = sv_regions_to_assemble + non_sv_regions_to_assemble
        regions = merge_small_regions(regions)        
        pybedtools.BedTool(regions).saveas(self.sample.regions_to_assemble)

    def create_assembly_scripts(self,flank=2000):     
        regions_to_assemble = pybedtools.BedTool(self.sample.regions_to_assemble)
        bashfiles = []
        for interval in regions_to_assemble:
            chrom = show_value(interval.chrom)
            start = show_value(interval.start)
            end = show_value(interval.end)
            hap = show_value(interval.name)
            if self.sample.add_unphased_reads and hap != "0":
                samtools_hap = "-r %s -r 0" % hap            
            else:
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

    def combine_sequence(self,type_,outfile):
        seqs = []
        bedfh = open(self.sample.regions_to_assemble, 'r')
        for line in bedfh:
            line = line.rstrip().split('\t')
            chrom,start,end,hap = line
            directory = "%s/assembly/%s/%s_%s/%s" % (self.sample.outdir,chrom,start,end,hap)
            contig = "%s/merged_contigs_quivered.%s" % (directory,type_)
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

    def filter_alignments_hap_specific(self,hap):
        pairs = []
        alignments = []
        columns = ["length","pident","nident","mismatch","gapopen","gaps","qseqid",
                   "qstart","qend","qlen","sseqid","sstart","send","slen","sstrand"]
        fosmid_names = set()
        Alignment = namedtuple('Alignment',columns)
        with open(self.sample.contig_to_contig_blast,'r') as fh:
            for line in fh:
                line = line.rstrip().split('\t')
                alignment = Alignment._make(line)
                fosmid_names.add(alignment.qseqid)
                fosmid_names.add(alignment.sseqid)
                if pass_filters(alignment,hap):
                    pairs.append(sorted((alignment.qseqid,alignment.sseqid)))
                    alignments.append(alignment)
        return alignments,fosmid_names

    def get_hap_merges(self):
        hap1_alignments,hap1_fosmid_names = self.filter_alignments_hap_specific("1")
        hap1_groupings = group_alignments(hap1_alignments,hap1_fosmid_names)
        hap1_groupings_with_fosmids_to_merge = group_merging_fosmids(hap1_alignments,hap1_fosmid_names,"1")
        hap2_alignments,hap2_fosmid_names = self.filter_alignments_hap_specific("2")


    def merge_sequences(self):
        self.run_blast(self.sample.locus_fasta,self.sample.contigs_to_contigs_blast)
        hap_merges = self.get_hap_merges()

        self.contigs_to_contigs_blast = None
        self.contigs_to_contigs_blast_edited = None
        self.contigs_grouped = None
        self.merge_alignments_instructions = None
        self.merged_contigs = None
        self.merged_contigs_to_ref = None
        get_possible_merges(fosmdis_to_fosmids_blast,blast_filtered,groupingsfn,fosmids_coords_to_mergefn,infosmids_fasta,outfosmids_fasta,length,errors)

        
    def __call__(self):
        self.get_phased_regions_to_assemble()
        assembly_scripts = self.create_assembly_scripts()
        run_assembly_scripts(assembly_scripts,self.sample.cluster,self.sample.cluster_walltime,
                             self.sample.cluster_threads,self.sample.cluster_mem,self.sample.cluster_queue)
        self.combine_sequences()
        self.command_line_tools.map_locus()
        self.merge_sequences()

def assemble_reads(self):
    assembly_runner = AssemblyRun(self)
    assembly_runner()

    
