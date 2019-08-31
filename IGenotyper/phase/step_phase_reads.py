#!/bin/env python
import os
import sys
import pysam
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..command_line import *
from ..common import *
        
def phase_mapped_reads(self):
    ccs_reads = "%s/ccs.bam" % self.tmp_dir
    get_ccs_reads(self.input_bam,self.tmp_dir,self.threads)     
    
    # Map subreads ands ccs reads
    map_reads(self.threads,self.pbmm2_ref,ccs_reads,
              self.ccs_mapped_reads,"ccs")
    map_reads(self.threads,self.pbmm2_ref,self.input_bam,
              self.subreads_mapped_reads,"subreads")  
    
    # Phase snps, and mapped subreads and ccs reads
    phase_snps(self.pbmm2_ref,self.snp_candidates,self.variants_vcf,
               self.phased_variants_vcf,self.ccs_mapped_reads,
               self.haplotype_blocks,self.phasing_stats)    
    phase_reads(self.phased_variants_vcf,self.ccs_mapped_reads,
                self.phased_ccs_mapped_reads,self.phased_vcf_file_sample_name)
    phase_reads(self.phased_variants_vcf,self.subreads_mapped_reads,
                self.phased_subreads_mapped_reads,self.phased_vcf_file_sample_name)
    
    # Get haplotype blocks
    get_phased_blocks(self.phased_variants_vcf,self.phasing_stats)