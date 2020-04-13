#!/bin/env python
import os
import sys
import pysam
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..command_line import *
from ..common import *

def read_is_unphased(read):
    haplotype = read.get_tag("RG",True)[0]
    if haplotype == "0":
        return True
    return False

def supplementary_score_diff(read,secondary_reads,thres):
    diffs = [] 
    primary_read_score = float(read.get_tag("AS",True)[0])
    for secondary_read in secondary_reads:
        secondary_read_score = float(secondary_read.get_tag("AS",True)[0])
        diffs.append(abs(primary_read_score-secondary_read_score))
    if min(diffs) > thres:
        return True
    return False

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
    # Change flag to primary alignment
    #secondary_to_primary_read.flag = primary_read.flag
    secondary_to_primary_read.flag = secondary_to_primary_read.flag - 256
    secondary_to_primary_read.mapping_quality = primary_read.mapping_quality
    # Change flag to secondary alignment
    #primary_read.flag = secondary_to_primary_flag
    primary_read.flag = primary_read.flag + 256
    primary_read.mapping_quality = secondary_to_primary_mapq
    #reads_to_return = [primary_read,secondary_to_primary_read]
    reads_to_return = [secondary_to_primary_read]
    for i,read in enumerate(secondary_reads):
        if i == index_to_ignore:
            continue
        reads_to_return.append(read)
    return reads_to_return

class PhaseRun():
    def __init__(self,Sample):
        self.sample = Sample
        self.command_line_tools = CommandLine(self.sample)

    def reads_phased(self):
        if non_emptyfile("%s/phasing.txt" % self.sample.tmp_dir):
            return True
        return False

    def get_initial_phasing(self):
        self.command_line_tools.get_ccs_reads()
        self.command_line_tools.turn_ccs_reads_to_fastq()
        self.command_line_tools.map_ccs_reads()
        self.command_line_tools.map_subreads()
        self.command_line_tools.phase_snps()
        self.command_line_tools.phase_ccs_reads()   
        self.command_line_tools.phase_subreads()             

    def get_secondary_alignments(self,bam):
        secondary_alignments = {}
        samfile = pysam.AlignmentFile(bam)
        for read in samfile.fetch("igh"):
            if read.is_secondary:
                if read.query_name not in secondary_alignments:
                    secondary_alignments[read.query_name] = []
                secondary_alignments[read.query_name].append(read)
        samfile.close()
        return secondary_alignments

    def change_primary_alignments(self,output_bamfile,bam):
        secondary_alignments = self.get_secondary_alignments(bam)
        samfile = pysam.AlignmentFile(bam)
        changed_reads = set()
        for read in samfile.fetch("igh"):
            if read.is_secondary:
                continue
            if read.is_supplementary:
                continue
            if read_is_unphased(read):
                if read.query_name in secondary_alignments:
                    secondary_reads = secondary_alignments[read.query_name]
                    # Always accept mapping quality less than 40
                    if read.mapping_quality > 40:
                        # If the scores is similar,accept                                                
                        if supplementary_score_diff(read,secondary_reads,self.sample.secondary_read_score):
                            continue
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
        changed_reads = self.change_primary_alignments(output_bamfile,bam)
        samfile = pysam.AlignmentFile(bam)
        for read in samfile:
            if read.query_name in changed_reads:
                continue
            output_bamfile.write(read)
        samfile.close()        
        output_bamfile.close()
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
        fn = "%s/alignments/subreads_to_ref_phased.sorted.bam" % self.sample.outdir
        os.remove(fn)
        os.remove("%s.bai" % fn)

    def save_previous_files(self,index,bam):
        if bam == self.sample.phased_ccs_mapped_reads:
            self.save_previous_variants(index)
            self.save_previous_ccs_mappings(index)
        else:
            self.remove_subreads_alignments(index)

    def rephase(self):
        for iter_ in ["0","1"]:
            iter_dir = "%s/variants/from_reads_%s" % (self.sample.outdir,iter_)
            if os.path.isdir(iter_dir):
                continue
            for bam in [self.sample.phased_ccs_mapped_reads,
                        self.sample.phased_subreads_mapped_reads]:
                self.fix_alignments(bam)
                self.save_previous_files(iter_,bam)
            self.command_line_tools.phase_snps()
            self.command_line_tools.phase_ccs_reads()    
            self.command_line_tools.phase_subreads()

    def check_phasing(self):
        finished = True
        fns = [self.sample.phased_ccs_mapped_reads,
               self.sample.phased_subreads_mapped_reads,
               self.sample.phased_variants_vcf,self.sample.variants_vcf]
        for fn in fns:
            if not non_emptyfile(fn):
                finished = False
        if finished:
            with open("%s/phasing.txt" % self.sample.tmp_dir,"w") as fh:
                fh.write("done")

    def clean_up(self):
        files_in_tmp_dir = ["ccs.fastq","ccs.fastq.gz","ccs_to_ref_all.bam","ccs_to_ref_all.sam",
                            "ccs_to_ref_all.sorted.bam","ccs_to_ref_all.sorted.bam.bai","changed_alignments.bam",
                            "changed_alignments.sam","subreads_to_ref_all.bam","subreads_to_ref_all.sam",
                            "subreads_to_ref_all.sorted.bam","subreads_to_ref_all.sorted.bam.bai"]
        files_in_alignment_dir = ["ccs_to_ref.sorted.bam","ccs_to_ref.sorted.bam.bai",
                                  "subreads_to_ref.sorted.bam","subreads_to_ref.sorted.bam.bai"]
        files_in_read_variants_dir = ["snp_candidates_filtered.vcf","snp_candidates.vcf"]
        dirs_in_alignment_dir = ["run_0","run_1"]
        dirs_in_variants_dir = ["from_reads_0","from_reads_1"]
        remove_files(self.sample.tmp_dir,files_in_tmp_dir)
        remove_files("%s/alignments/" % self.sample.outdir,files_in_alignment_dir)
        remove_files("%s/variants/from_reads" % self.sample.outdir,files_in_read_variants_dir)
        remove_dirs("%s/alignments" % self.sample.outdir,dirs_in_alignment_dir)
        remove_dirs("%s/variants" % self.sample.outdir,dirs_in_variants_dir)

    def __call__(self):
        if not self.reads_phased():
            self.get_initial_phasing()
            self.rephase()
            self.command_line_tools.get_phased_blocks()
            self.check_phasing()
            if not self.sample.keep:
                self.clean_up()
        
def phase_mapped_reads(self):
    phase_runner = PhaseRun(self)
    phase_runner()

