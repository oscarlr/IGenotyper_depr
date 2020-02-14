#!/bin/env python
from read_stats import stats_on_reads
from assembly_stats import stats_on_assembly

def generate_stats(self):
    ccs_bam = "%s/ccs.bam" % self.tmp_dir
    stats_on_reads(self.input_bam,ccs_bam,self.phased_ccs_mapped_reads,
                   self.phased_subreads_mapped_reads,self.plots_dir,
                   self.bedgraph_dir,self.tables_dir,self.gene_coordinates,
                   self.sv_regions,self.region_types)
    stats_on_assembly(self.locus_fasta,self.tables_dir,self.mapped_locus)

    
