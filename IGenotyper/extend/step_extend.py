#!/bin/env python
from collections import namedtuple
from ..command_line import *
from ..common import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_blast_on_contigs(locus_fasta,locus_to_locus_blast):    
    run_blast(locus_fasta,locus_fasta,blast_output)

def not_overlapping(alignment):
    if not (int(alignment.qstart) < 100 and (int(alignment.slen) - int(alignment.send)) < 100):
        if not (int(alignment.sstart) < 100 and (int(alignment.qlen) - int(alignment.qend)) < 100):
            return True
    return False

def get_entry_location(entry):
    hap = entry.split("_")[1].split('=')[1]
    loc = entry.split("_")[0].split('=')[1]
    start = int(loc.split(":")[1].split("-")[0])
    end = int(loc.split(":")[1].split("-")[1])
    return (start,end,hap)

    
def diff_hap(alignment):
    q_start,q_end,q_hap = get_entry_location(alignment.qseqid)
    s_start,s_end,s_hap = get_entry_location(alignment.sseqid)
    if (q_end == s_start) or (s_end == q_start):
        if q_hap == "1" and s_hap == "2":
            return True
        if q_hap == "2" and s_hap == "1":
            return True
    return False

def same_coord_diff_hap(alignment):
    q_start,q_end,q_hap = get_entry_location(alignment.qseqid)
    s_start,s_end,s_hap = get_entry_location(alignment.sseqid)
    if q_start == s_start and q_end == s_end:
        if q_hap != s_hap:
            return True
    return False
   
def filter_alignments(locus_to_locus_blast):
    pairs = []
    alignments = []
    pident_min = 95.0            
    mismatch_max = 3
    columns = ["length","pident","nident","mismatch","gapopen","gaps","qseqid",
               "qstart","qend","qlen","sseqid","sstart","send","slen","sstrand"]
    Alignment = namedtuple('Alignment',columns)
    with open(locus_to_locus_blast,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            alignment = Alignment._make(line)            
            if alignment.qseqid == alignment.sseqid:
                continue
            if float(alignment.pident) < pident_min:
                continue
            if float(alignment.mismatch) > mismatch_max:
                continue
            if alignment.sstrand != "plus":
                continue            
            if not_overlapping(alignment):
                continue
            if diff_hap(alignment):
                continue
            if same_coord_diff_hap(alignment):
                continue            
            if int(alignment.sstart) > int(alignment.qstart):
                continue
            if sorted((alignment.qseqid,alignment.sseqid)) in pairs:
                continue
            pairs.append(sorted((alignment.qseqid,alignment.sseqid)))
            q_start,q_end,q_hap = get_entry_location(alignment.qseqid)
            s_start,s_end,s_hap = get_entry_location(alignment.sseqid)
            output = ["igh",q_start,q_end,q_hap,"igh",s_start,s_end,s_hap] + line
            alignments.append(output)
    return sorted(alignments,key = lambda x: x[1])

def contig_by_contig_identity_matrix():
    pass

def get_adjacent_haplotypes(haplotype_blocks):
    regions = {}
    with open(haplotype_blocks,'r') as fh:
        for line in fh:
            line = line.strip().split('\t')
            region = [line[0],line[1],line[2]]
            hap = line[3]
            if region in regions:
                regions[region] = []
            regions[region].append(hap)
    regions_output = []
    for region in regions:
        region_with_hap = region
        haps = []
        for hap in regions[region]:
            haps.append(hap)
        region_with_hap.append(haps)
        regions_output.append(region_with_hap)
    return sorted(regions_output,key = lambda x: x[1])

def save_contigs_with_no_alignments(alignments,single_contigs_to_add,locus_fasta):
    contigs = SeqIO.to_dict(SeqIO.parse(locus_fasta,"fasta"))
    alignment_contigs = [i[14] for i in alignments] + [i[18] for i in alignments]
    with open(single_contigs_to_add,'w') as fh:
        for contig in contigs:
            if contig in alignment_contigs:
                continue
            fh.write("%s\n" % contig)

def save_alignments(alignments,alignments_file):
    with open(alignments_file,'w') as fh:
        for alignment in alignments:
            fh.write("%s\n" % "\t".join(map(str,alignment)))

def read_merge_instructions(merge_alignments_instructions):
    columns = ["contig_group","igh_1","qseqid_start_1","qseqid_end_1","qseqid_hap_1",
               "igh_2","sseqid_start_1","sseqid_end_1","sseqid_hap_1",
               "length","pident","nident","mismatch","gapopen","gaps",
               "qseqid","qstart","qend","qlen",
               "sseqid","sstart","send","slen",
               "sstrand"]
    Alignment = namedtuple('Alignment',columns)
    groupings = {}
    with open(merge_alignments_instructions,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            alignment = Alignment._make(line)            
            if alignment.contig_group not in groupings:
                groupings[alignment.contig_group] = []
                groupings[alignment.contig_group].append([alignment.qseqid,1,alignment.qlen,alignment.qseqid_hap_1])
                groupings[alignment.contig_group].append([alignment.sseqid,alignment.send,alignment.slen,alignment.sseqid_hap_1])
                continue
            assert alignment.qseqid == groupings[alignment.contig_group][-1][0]
            groupings[alignment.contig_group].append([alignment.sseqid,alignment.send,alignment.slen,alignment.sseqid_hap_1])
    return groupings

def merge_contigs(groupings,locus_fasta,merged_contigs,single_contigs_to_add):
    contigs = SeqIO.to_dict(SeqIO.parse(locus_fasta,"fasta"))
    contigs_to_output = []
    contigs_to_be_merged = [contigs[i][0] for i in contigs]    
    for group in groupings:
        hap=None
        starts=[]
        ends=[]
        sequences = []
        for contig_name,contig_start,contig_end,contig_hap in groupings[group]:            
            if hap == None or hap == "0":
                hap = contig_hap            
            locus_start,locus_end,c_hap = get_entry_location(contig_name)
            starts.append(int(locus_start))
            ends.append(int(locus_end))
            sequence = str(contigs[contig_name].seq[int(contig_start):int(contig_end)])
            sequences.append(sequence)
        sequence = "".join(sequences)
        sequence_name = "coord=igh:%s-%s_hap=%s_index=%s_total=merged_/0/0_0" % (min(starts),max(ends),hap,group)
        record = SeqRecord(Seq(sequence,"fasta"),id=sequence_name,name="",description="")
        contigs_to_output.append(record)
    with open(single_contigs_to_add,'r') as fh:
        for line in fh:
            contig = line.rstrip()
            contigs_to_output.append(contigs[contig])
    SeqIO.write(contigs_to_output,merged_contigs,"fasta")

def get_possible_merges(self):
    run_blast(self.locus_fasta,self.locus_to_locus_blast)
    alignments = filter_alignments(self.locus_to_locus_blast)
    save_alignments(alignments,self.contig_alignments)
    save_contigs_with_no_alignments(alignments,self.single_contigs_to_add,self.locus_fasta)
    if non_emptyfile(self.merge_alignments_instructions):
        groupings = read_merge_instructions(self.merge_alignments_instructions)
        merge_contigs(groupings,self.locus_fasta,self.merged_contigs,self.single_contigs_to_add)
        map_reads("1",self.pbmm2_ref,self.merged_contigs,self.merged_contigs_to_ref,"ccs")
