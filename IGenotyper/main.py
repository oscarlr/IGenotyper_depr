#!/bin/env python
import os
import sys
import argparse

from common import *
from file_names import FileNames

class Sample(object):
    from phase.step_phase_reads import phase_mapped_reads
    from assemble.step_assemble_reads import assemble_reads
    from extend.step_extend import get_possible_merges
    from detect.step_detect import detect_variants
    from stats.step_stats import generate_stats
    # from step_phase_stats import plot_phase_stats
    # from step_assemble_stats import plot_assemble_stats

    def __init__(self):
        self.input_bam = None
        self.outdir = None

        # Steps
        self.phase = None
        self.assemble = None
        self.extend_assembly = None
        self.detect = None
        self.stats = None

        # All
        self.pbmm2_ref = None
        self.tmp_dir = None
        
        # CPU Params
        self.threads = None        
        self.cluster = None
        self.cluster_queue = None
        self.cluster_threads = None
        self.cluster_walltime = None
        self.cluster_mem = None

        # Other
        self.haploid = None
        self.sv_regions = None
        self.non_sv_regions = None
        self.introns = None
        self.lpart1 = None
        self.rss = None
        self.region_types = None
        self.regions_to_ignore = None
        self.python_scripts = None

        # Phasing
        self.ccs_mapped_reads = None
        self.subreads_mapped_reads = None
        self.phased_ccs_mapped_reads = None
        self.phased_subreads_mapped_reads = None
        self.variants_vcf = None
        self.phased_variants_vcf = None
        self.haplotype_blocks = None
        self.snp_candidates = None
        self.snp_candidates_filtered = None
        #self.phasing_stats = None
        self.phased_vcf_file_sample_name = None

        # Assembly
        self.regions_to_assemble = None
        self.locus_fasta = None
        self.locus_fastq = None
        self.add_unphased_reads = None
        self.assembly_script = None
        self.phased_regions_with_coverage = None

        # Extend assembly
        self.locus_to_locus_blast = None
        self.contig_alignments = None
        self.merge_alignments_instructions = None
        self.merged_contigs = None
        self.merged_contigs_to_ref = None
        self.single_contigs_to_add = None
        
        # Detect
        self.alleles = None
        self.mapped_locus = None
        self.snps_in_sv_regions = None
        self.snps_not_in_sv_regions = None
        self.svs_genotyped = None
        self.gene_coordinates = None
        self.genes_from_assembly = None
        self.genes_from_reads = None
        self.allele_database = None
        self.novel_alleles = None
        self.genes_with_allele_assignment = None
        self.indels = None
        self.sv_signature = None
        self.sv_vcf = None
        
        # Stats
        self.stats_dir = None
        self.plots_dir = None
        self.bedgraph_dir = None
        self.tables_dir = None

    command_line_args_to_attrs = [
        ("input_bam","input_bam"),
        ("outdir","outdir"),
        ("phase","phase"),        
        ("assemble","assemble"),
        ("extend_assembly","extend_assembly"),
        ("detect","detect"),        
        ("stats","stats"),
        ("tmp_dir","tmp_dir"),
        ("threads","threads"),        
        ("cluster","cluster"),
        ("cluster_queue","cluster_queue"),
        ("cluster_threads","cluster_threads"),
        ("cluster_walltime","cluster_walltime"),
        ("cluster_mem","cluster_mem"),
        ("haploid","haploid"),
        ("phased_vcf_file","phased_variants_vcf"),
        ("phased_vcf_file_sample_name","phased_vcf_file_sample_name"),
        ("add_unphased_reads","add_unphased_reads")
        ]
    
    def load_attrs(self):
        file_names = FileNames(self.outdir)
        for attr in self.__dict__.keys():
            if getattr(self,attr) != None:
                continue
            setattr(self,attr,getattr(file_names,attr))

    @classmethod
    def load_args(cls,args):
        sample = cls()
        # Load command line options
        for arg, attr in sample.command_line_args_to_attrs:
            setattr(sample,attr,getattr(args,arg))         
        sample.load_attrs() # Load package options
        #sample.overwrite_params()
        return sample
        
    def __call__(self):
        create_directory(self.outdir)
        if self.phase:
            self.phase_mapped_reads()
        elif self.assemble:
            self.assemble_reads()
        elif self.extend_assembly:            
            self.get_possible_merges()
        elif self.detect:
            self.detect_variants()
        elif self.stats:
            self.generate_stats()
        elif self.phase_stats:
            self.plot_phase_stats()
        elif self.assemble_stats:
            self.plot_assemble_stats()
        else:
            sys.exit("Incorrect command")

def main():
    parser = argparse.ArgumentParser(description='Process IGH capture data')
    parser.add_argument('input_bam', metavar='input_bam', type=check_file_exist,
                        help='bam file containing raw reads')
    parser.add_argument('outdir', metavar='outdir', 
                        help='output directory')
    parser.add_argument('--phase', action='store_true', default=False,
                        help='Map and phase reads')
    parser.add_argument('--assemble', action='store_true', default=False,
                        help='Only assemble reads')
    parser.add_argument('--extend_assembly', action='store_true', default=False,
                        help='Extend assemblies')
    parser.add_argument('--detect', action='store_true', default=False,
                        help='Detect variants')
    parser.add_argument('--stats', action='store_true', default=False,
                        help='Generate stats')
    parser.add_argument('--tmp_dir', metavar='tmp_dir',
                        help='temporary directory')
    parser.add_argument('--threads', metavar='threads', default=1,
                        help='Number of threads for everything')
    parser.add_argument('--cluster', action='store_true', default=False,
                        help='Use cluster')
    parser.add_argument('--cluster_queue', metavar="", default="express",
                        help='Queue for cluster')
    parser.add_argument('--cluster_threads', metavar='', default=4,
                        help='Number of threads for cluster jobs')
    parser.add_argument('--cluster_walltime', metavar="", default=2,
                        help='Walltime for cluster')
    parser.add_argument('--cluster_mem', metavar="", default=5,
                        help='memory for cluster')
    parser.add_argument('--haploid',action='store_true', default=False,
                        help='Run Quiver in haploid mode')
    parser.add_argument('--phased_vcf_file', type=check_file_exist,
                        help='Run IG with this phased VCF file')
    parser.add_argument('--phased_vcf_file_sample_name',default="sample",
                        help='Sample name in phased VCF file')
    parser.add_argument('--add_unphased_reads',action='store_true', default=False,
                        help='Add unphased reads to phased region')
    args = parser.parse_args()
    sample = Sample.load_args(args)
    return sample()

