#!/bin/env python
import os
import sys
import pysam
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..command_line import *
from ..common import *

def load_het_snvs(snvs):
    positions = []
    with open(snvs,'r') as fh:
        for line in fh:
            if "#" in line:
                continue
            if "0/1" in line or "0|1" in line or "1|0" in line:
                line = line.rstrip().split('\t')
                positions.append(int(line[1]))
    return sorted(positions)

def read_is_unphased(read):
    haplotype = read.get_tag("RG",True)[0]
    if haplotype == "0":
        return True
    return False

def read_overlap_hets(read,het_snvs):
    het_snvs = het_snvs
    for snv in het_snvs:
        if snv < read.reference_start:
            continue
        if snv > read.reference_end:
            return False
        return True
    return False

# def get_secondary_alignments(phased_ccs_reads):
#     secondary_alignments = {}
#     samfile = pysam.AlignmentFile(phased_ccs_reads)
#     for read in samfile:
#         if read.is_secondary:
#             if read.query_name not in secondary_alignments:
#                 secondary_alignments[read.query_name] = []
#             secondary_alignments[read.query_name].append(read)
#     return secondary_alignments

def change_read(primary_read,secondary_reads):
    secondary_to_primary_read = None
    secondary_to_primary_flag = None
    secondary_to_primary_mapq = None
    index_to_ignore = None
    secondary_to_primary_score = 0
    for i,secondary_read in enumerate(secondary_reads):
        secondary_read_score = secondary_read.get_tag("AS",True)[0]
        if secondary_read_score < secondary_to_primary_score:
            secondary_to_primary_read = secondary_read
            secondary_to_primary_score = secondary_read_score
            secondary_to_primary_flag = secondary_read.flag
            secondary_to_primary_mapq = secondary_read.mapping_quality
            index_to_ignore = i
    secondary_to_primary_read.flag = primary_read.flag
    secondary_to_primary_read.mapping_quality = primary_read.mapping_quality
    primary_read.flag = secondary_to_primary_flag
    primary_read.mapping_quality = secondary_to_primary_mapq
    reads_to_return = [primary_read,secondary_to_primary_read]
    for i,read in enumerate(secondary_reads):
        if i == index_to_ignore:
            continue
        reads_to_return.append(read)
    return reads_to_return

# def change_primary_alignments(phased_ccs_reads,output_bamfile,het_snvs):
#     secondary_alignments = get_secondary_alignments(phased_ccs_reads)
#     samfile = pysam.AlignmentFile(phased_ccs_reads)
#     changed_reads = set()
#     for read in samfile.fetch("igh",476712,571415):
#         if read.is_secondary:
#             continue
#         if read.is_supplementary:
#             continue
#         if read.mapping_quality > 40:
#             continue
#         if read_is_unphased(read): #and read_overlap_hets(read,het_snvs):
#             if read.query_name in secondary_alignments:
#                 secondary_reads = secondary_alignments[read.query_name]
#                 output_reads = change_read(read,secondary_reads)
#                 for output_read in output_reads:
#                     changed_reads.add(output_read.query_name)
#                     output_bamfile.write(output_read)
#     samfile.close()
#     return changed_reads

# def fix_alignments(tmp_dir,phased_ccs_reads,snvs):
#     # 1. Find the unphased reads that overlap heterozygous SNPs
#     # 2. Check if the secondary alignment of the unphased read is phased
#     # 3. Change the secondary alignment to the primary alignment
#     # 4. If the scondary alignment is still not phased, move it to the primary alignment
#     het_snvs = load_het_snvs(snvs)
#     changed_bamfile = "%s/changed_alignments.sam" % tmp_dir
#     samfile = pysam.AlignmentFile(phased_ccs_reads)
#     output_bamfile = pysam.AlignmentFile(changed_bamfile,"w",template=samfile)
#     samfile.close()
#     samfile = pysam.AlignmentFile(phased_ccs_reads)
#     changed_reads = change_primary_alignments(phased_ccs_reads,output_bamfile,het_snvs)
#     for read in samfile:
#         if read.query_name in changed_reads:
#             continue
#         output_bamfile.write(read)
#     samfile.close()        
#     return changed_bamfile

# def save_previous_files(index,outdir):
#     src = "%s/variants/from_reads" % outdir
#     dst = "%s/variants/from_reads_%s" % (outdir,index)
#     os.mkdir(dst)
#     os.rename(src,dst)
#     file_names = ["ccs_to_ref_phased.sorted.bam","ccs_to_ref_phased.sorted.bam.bai",
#                   "ccs_to_ref.sorted.bam","ccs_to_ref.sorted.bam.bai"]
#     os.mkdir("%s/alignments/run_%s" % (outdir,index))
#     for file_name in file_names:
#         src = "%s/alignments/%s" % (outdir,file_name)
#         dst = "%s/alignments/run_%s/%s" % (outdir,index,file_name)
#         os.rename(src,dst)
#     src = "%s/tmp/changed_alignments.sorted.bam" % outdir
#     dst = "%s/alignments/ccs_to_ref.sorted.bam" % (outdir,index)
#     os.rename(src,dst)
#     src = "%s/tmp/changed_alignments.sorted.bam.bai" % outdir
#     dst = "%s/alignments/ccs_to_ref.sorted.bam.bai" % outdir
#     os.rename(src,dst)

class PhaseRun():
    def __init__(self,Sample):
        self.sample = Sample
        self.command_line_tools = CommandLine(self.sample)

    def get_initial_phasing(self):
        self.command_line_tools.get_ccs_reads()
        self.command_line_tools.turn_ccs_reads_to_fastq()
        self.command_line_tools.map_ccs_reads()
        self.command_line_tools.map_subreads()
        self.command_line_tools.phase_snps()
        self.command_line_tools.phase_ccs_reads()                

    def get_secondary_alignments(self,bam):
        secondary_alignments = {}
        samfile = pysam.AlignmentFile(bam)
        for read in samfile:
            if read.is_secondary:
                if read.query_name not in secondary_alignments:
                    secondary_alignments[read.query_name] = []
                secondary_alignments[read.query_name].append(read)
        return secondary_alignments

    def change_primary_alignments(self,output_bamfile,bam):
        start = 476712
        end = 571415
        secondary_alignments = self.get_secondary_alignments(bam)
        samfile = pysam.AlignmentFile(bam)
        changed_reads = set()
        for read in samfile.fetch("igh",start,end):
            if read.is_secondary:
                continue
            if read.is_supplementary:
                continue
            if read.mapping_quality > 40:
                continue
            if read_is_unphased(read):
                if read.query_name in secondary_alignments:
                    secondary_reads = secondary_alignments[read.query_name]
                    output_reads = change_read(read,secondary_reads)
                    for output_read in output_reads:
                        changed_reads.add(output_read.query_name)
                        output_bamfile.write(output_read)
        samfile.close()
        return changed_reads

    def get_changed_bamfile(self,bam):
        changed_bamfile = "%s/changed_alignments.sam" % self.sample.tmp_dir
        samfile = pysam.AlignmentFile(bam)
        output_bamfile = pysam.AlignmentFile(changed_bamfile,"w",template=samfile)
        samfile.close()        
        return output_bamfile

    def fix_alignments(self,bam):
        # 1. Find the unphased reads
        # 2. Check if the secondary alignment of the unphased read is phased
        # 3. Change the secondary alignment to the primary alignment
        # 4. If the scondary alignment is still not phased, move it to the primary alignment
        changed_bamfile = "%s/changed_alignments.sam" % self.sample.tmp_dir
        output_bamfile = self.get_changed_bamfile(bam)
        samfile = pysam.AlignmentFile(bam)
        changed_reads = self.change_primary_alignments(output_bamfile,bam)
        for read in samfile:
            if read.query_name in changed_reads:
                continue
            output_bamfile.write(read)
        samfile.close()        
        prefix = changed_bamfile[:-4]
        self.command_line_tools.sam_to_sorted_bam(prefix,"%s.sorted.bam" % prefix)

    def save_previous_variants(self,index):
        src = "%s/variants/from_reads" % self.sample.outdir
        dst = "%s/variants/from_reads_%s" % (self.sample.outdir,index)
        os.mkdir(dst)
        os.rename(src,dst)
        os.mkdir("%s/variants/from_reads" % self.sample.outdir)
        
    def save_previous_ccs_mappings(self,index):
        file_names = ["ccs_to_ref_phased.sorted.bam","ccs_to_ref_phased.sorted.bam.bai",
                      "ccs_to_ref.sorted.bam","ccs_to_ref.sorted.bam.bai"]
        os.mkdir("%s/alignments/run_%s" % (self.sample.outdir,index))
        for file_name in file_names:
            src = "%s/alignments/%s" % (self.sample.outdir,file_name)
            dst = "%s/alignments/run_%s/%s" % (self.sample.outdir,index,file_name)
            os.rename(src,dst)        
        src = "%s/tmp/changed_alignments.sorted.bam" % self.sample.outdir
        dst = "%s/alignments/ccs_to_ref.sorted.bam" % self.sample.outdir
        os.rename(src,dst)
        src = "%s/tmp/changed_alignments.sorted.bam.bai" % self.sample.outdir
        dst = "%s/alignments/ccs_to_ref.sorted.bam.bai" % self.sample.outdir
        os.rename(src,dst)        

    def remove_subreads_alignments(self,index):
        src = "%s/tmp/changed_alignments.sorted.bam" % self.sample.outdir
        dst = "%s/alignments/subreads_to_ref.sorted.bam" % self.sample.outdir
        os.rename(src,dst)
        src = "%s/tmp/changed_alignments.sorted.bam.bai" % self.sample.outdir
        dst = "%s/alignments/subreads_to_ref.sorted.bam.bai" % self.sample.outdir
        os.rename(src,dst)        

    def save_previous_files(self,index,bam):
        self.save_previous_variants(index)
        if bam == self.sample.phased_ccs_mapped_reads:
            self.save_previous_ccs_mappings(index)
        else:
            self.remove_subreads_alignments(index)

    def rephase(self):
        for iter in ["0","1"]:
            iter_dir = "%s/variants/from_reads_%s" % (self.sample.outdir,iter)
            if os.path.isdir(iter_dir):
                continue
            for bam in [self.sample.phased_ccs_mapped_reads,
                        self.sample.phased_subreads_mapped_reads]:
                self.fix_alignments(bam)
                self.save_previous_files(iter,bam)
            self.command_line_tools.phase_snps()
            self.command_line_tools.phase_ccs_reads()    
            self.command_line_tools.phase_subreads()

    def __call__(self):
        self.get_initial_phasing()
        self.rephase()
        self.command_line_tools.get_phased_blocks()

def phase_mapped_reads(self):
    phase_runner = PhaseRun(self)
    phase_runner()

