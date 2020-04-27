#!/bin/env python
import os
import sys
import argparse
import shutil

from common import non_emptyfile,create_directory
from file_names import FileNames
from command_line import CommandLine

class LoadTool(object):
    def __init__(self):
        self.command_line_args_to_attrs = []

    def load_args(self,args):
        # Load command line options
        for arg, attr in self.command_line_args_to_attrs:
            setattr(self,attr,getattr(args,arg))         

class FileManager(LoadTool):
    # The job of this class is to manage the files
    def __init__(self):
        # Inputs
        self.input_bam = None
        self.outdir = None

        # References
        self.pbmm2_ref = None
        self.blasr_ref = None

        # Directories
        self.tmp_dir = None

        # Package files
        self.sv_regions = None
        self.non_sv_regions = None
        self.introns = None
        self.lpart1 = None
        self.rss = None
        self.region_types = None
        self.regions_to_ignore = None
        self.python_scripts = None
        self.r_scripts = None
        self.igh_coords = None
        self.igh_fasta = None
        self.igh_fasta_fai = None
        self.non_dup_regions = None
        self.dup_regions = None
        self.aims = None
        self.target_regions = None

        # Phasing files
        self.ccs_reads = None
        self.ccs_mapped_reads = None
        self.subreads_mapped_reads = None
        self.phased_ccs_mapped_reads = None
        self.phased_subreads_mapped_reads = None
        self.variants_vcf = None
        self.phased_variants_vcf = None
        self.haplotype_blocks = None
        self.snp_candidates = None
        self.snp_candidates_filtered = None

        # Assembly files
        self.locus_fasta = None
        self.locus_fasta_unquivered = None
        self.locus_fastq = None
        #self.add_unphased_reads = None
        self.assembly_script = None
        self.regions_to_assemble = None
        self.phased_regions_with_coverage = None
        ## Merging
        self.contigs_to_contigs_blast = None
        self.contigs_to_contigs_blast_edited = None
        self.contigs_grouped = None
        self.merge_alignments_instructions = None
        self.merged_contigs = None
        self.merged_contigs_to_ref = None
        self.locus_fasta_unquivered_to_ref = None

        # Detect files
        self.alleles = None
        self.mapped_locus = None
        self.assembly_snps = None
        self.assembly_indels = None
        self.assembly_svs = None
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
        self.reassembly_gene_script = None
        self.sv_calling_script = None

        # Report files
        self.html_report = None
        self.report_template = None

        # Stats files
        self.stats_dir = None
        self.plots_dir = None
        self.bedgraph_dir = None
        self.tables_dir = None

        self.command_line_args_to_attrs = [
            ("input_bam","input_bam"),
            ("outdir","outdir"),
            ("tmp_dir","tmp_dir"),
            ("phased_vcf_file","phased_variants_vcf")
        ]
        
    def load_args(self,args):
        # Load command line options
        for arg, attr in self.command_line_args_to_attrs:
            setattr(self,attr,getattr(args,arg))         
        create_directory(self.outdir)
        file_names = FileNames(self.outdir)
        for attr in self.__dict__.keys():
            if getattr(self,attr) != None:
                continue
            setattr(self,attr,getattr(file_names,attr))

class CpuManager(LoadTool):
    # The job of this class is cpu management
    # How does it handle different cpu requirements for 
    # different steps?
    def __init__(self):
        self.threads = None        
        self.cluster = None
        self.cluster_queue = None
        #self.cluster_threads = None
        self.cluster_walltime = None
        self.cluster_mem = None

        self.command_line_args_to_attrs = [
            ("threads","threads"),        
            ("cluster","cluster"),
            ("cluster_queue","cluster_queue"),
            ("cluster_walltime","cluster_walltime"),
            ("cluster_mem","cluster_mem")
        ]

class Step(LoadTool):
    def __init__(self, file_manager, cpu_manager, command_line_tools):
        self.file_manager = file_manager
        self.cpu_manager = cpu_manager
        self.command_line_tools = command_line_tools

        self.done_file = None
        self.files_to_check = None

        # Parameters
        self.pacbio_data_type = None
        self.keep = None
        self.sample_name = None
        self.secondary_read_score = None
        self.split = None
        self.add_hom_ref_genotype = None
        
        self.command_line_args_to_attrs = [
            ("pacbio_data_type","pacbio_data_type"),
            ("keep","keep"),
            ("sample_name","sample_name"),        
            ("secondary_read_score","secondary_read_score"),
            ("split","split"),
            ("add_hom_ref_genotype","add_hom_ref_genotype")
        ]

    def step_complete(self):
        complete = False
        assert self.done_file != None
        if non_emptyfile(self.done_file):
            complete = True
        return complete            

    def step_completed(self):
        finished = True
        for fn in self.files_to_check:
            if not non_emptyfile(fn):
                finished = False
        return finished

    def clean_up(self):
        if not self.keep:
            if os.path.exists(self.file_manager.tmp_dir):
                shutil.rmtree(self.file_manager.tmp_dir)
    

    def __call__(self):
        if not self.step_complete():
            self.run()
        print "Step completed"
