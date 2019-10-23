#!/bin/env python
import os

from common import *

class FileNames():
    def __init__(self,directory):

        # Package data
        self.pbmm2_ref = None
        self.gene_coordinates = None
        self.sv_regions = None
        self.non_sv_regions = None
        self.allele_database = None
        self.introns = None
        self.lpart1 = None
        self.rss = None
        self.region_types = None
        self.regions_to_ignore = None

        # Alignments
        self.ccs_mapped_reads = None
        self.subread_mapped_reads = None
        self.phased_ccs_mapped_reads = None
        self.phased_subreads_mapped_reads = None
        self.mapped_locus = None

        # Assembly 
        self.locus_fasta = None
        self.locus_fastq = None
        self.regions_to_assemble = None
        self.haplotype_blocks = None
        self.phased_regions_with_coverage = None

        # Extend assembly
        self.locus_to_locus_blast = None
        self.contig_alignments = None        
        self.merge_alignments_instructions = None
        self.merged_contigs = None
        self.merged_contigs_to_ref = None
        self.single_contigs_to_add = None

        # Variants
        self.snp_candidates = None
        self.snp_candidates_filtered = None
        self.variants_vcf = None
        self.phased_variants_vcf = None
        self.snps_in_sv_regions = None
        self.snps_not_in_sv_regions = None
        self.indels = None
        self.svs_genotyped = None
        self.sv_signature = None
        self.sv_vcf = None

        # Alleles
        self.alleles = None
        self.genes_with_allele_assignment = None
        self.genes_from_assembly = None
        self.genes_from_reads = None
        self.novel_alleles = None

        # Stats
        #self.phasing_stats = None
        self.stats_dir = None
        self.plots_dir = None
        self.bedgraph_dir = None
        self.tables_dir = None

        # Tmp
        self.tmp_dir = None

        self.package_data()
        self.set_alignments(directory)
        self.set_assembly(directory)
        self.set_extend_assembly(directory)
        self.set_variants(directory)
        self.set_alleles(directory)
        self.set_stats(directory)
        self.set_tmp(directory)

    def package_data(self):
        directory = os.path.dirname(os.path.abspath(__file__))
        self.pbmm2_ref = "%s/data/pbmm2_index/reference.fasta" % directory
        self.gene_coordinates = "%s/data/gene_coords.bed" % directory
        self.sv_regions = "%s/data/sv_coords.bed" % directory
        self.non_sv_regions = "%s/data/non_sv_coords.bed" % directory
        self.allele_database = "%s/data/vdj_alleles.fasta" % directory
        self.introns = "%s/data/introns.bed" % directory
        self.lpart1 = "%s/data/lpart1.bed" % directory
        self.rss = "%s/data/rss.bed" % directory
        self.region_types = "%s/data/regions.bed" % directory
        self.regions_to_ignore = "%s/data/regions_to_ignore.bed" % directory
        self.assembly_script = "%s/bash_scripts/assemble.sh" % directory
        self.python_scripts = "%s/python_scripts" % directory

    def set_alignments(self,directory):
        folder_name = "alignments"
        create_directory("%s/%s" % (directory,folder_name))
        self.ccs_mapped_reads = "%s/%s/ccs_to_ref.sorted.bam" % (directory,folder_name)
        self.subreads_mapped_reads = "%s/%s/subreads_to_ref.sorted.bam" % (directory,folder_name)
        self.phased_ccs_mapped_reads = "%s/%s/ccs_to_ref_phased.sorted.bam" % (directory,folder_name)
        self.phased_subreads_mapped_reads = "%s/%s/subreads_to_ref_phased.sorted.bam" % (directory,folder_name)
        self.mapped_locus = "%s/%s/locus_to_ref.sorted.bam" % (directory,folder_name)

    def set_assembly(self,directory):
        folder_name = "assembly"
        create_directory("%s/%s" % (directory,folder_name))
        self.locus_fasta = "%s/%s/locus.fasta" % (directory,folder_name)
        self.locus_fastq = "%s/%s/locus.fastq" % (directory,folder_name)
        self.regions_to_assemble = "%s/%s/ccs_variants_phased.bed" % (directory,folder_name)
        self.haplotype_blocks = "%s/%s/ccs_variants_phased.tab" % (directory,folder_name)
        self.phased_regions_with_coverage = "%s/%s/phased_regions_with_coverage.bed" % (directory,folder_name)

    def set_extend_assembly(self,directory):
        folder_name = "extend"
        create_directory("%s/%s" % (directory,folder_name))
        self.locus_to_locus_blast = "%s/extend/locus_to_locus.blast" % directory
        self.contig_alignments = "%s/extend/self_locus_blast.txt" % directory
        self.merge_alignments_instructions = "%s/extend/self_locus_blast_edited.txt" % directory
        self.merged_contigs = "%s/extend/merged_locus.fasta" % directory
        self.merged_contigs_to_ref = "%s/extend/merged_locus_to_ref.sorted.bam" % directory
        self.single_contigs_to_add = "%s/extend/singleton_contigs.txt" % directory

    def set_variants(self,directory):
        folder_name = "variants"
        create_directory("%s/%s" % (directory,folder_name))
        create_directory("%s/%s/from_reads" % (directory,folder_name))
        create_directory("%s/%s/from_assembly" % (directory,folder_name))
        self.snp_candidates = "%s/%s/from_reads/snp_candidates.vcf" % (directory,folder_name)
        self.snp_candidates_filtered = "%s/%s/from_reads/snp_candidates_filtered.vcf" % (directory,folder_name)
        self.variants_vcf = "%s/%s/from_reads/ccs_variants.vcf" % (directory,folder_name)
        self.phased_variants_vcf = "%s/%s/from_reads/ccs_phased_variants.vcf" % (directory,folder_name)
        self.snps_in_sv_regions = "%s/%s/from_assembly/snps_in_svs.vcf" % (directory,folder_name)
        self.snps_not_in_sv_regions = "%s/%s/from_assembly/snps_not_in_svs.vcf" % (directory,folder_name)
        self.indels = "%s/%s/from_assembly/indels" % (directory,folder_name)
        self.svs_genotyped = "%s/%s/from_assembly/svs.txt" % (directory,folder_name)
        self.sv_signature = "%s/%s/from_assembly/sv.svsig.gz" % (directory,folder_name)
        self.sv_vcf = "%s/%s/from_assembly/sv.vcf" % (directory,folder_name)
        
    def set_alleles(self,directory):
        folder_name = "alleles"
        create_directory("%s/%s" % (directory,folder_name))
        self.alleles = "%s/%s/alleles.tab" % (directory,folder_name)
        self.genes_with_allele_assignment = "%s/%s/genes_assigned_to_alleles.txt" % (directory,folder_name)
        self.genes_from_assembly = "%s/%s/genes_from_assembly.fasta" % (directory,folder_name)
        self.genes_from_reads = "%s/%s/genes_from_ccs_reads.fasta" % (directory,folder_name)
        self.novel_alleles = "%s/%s/novel_alleles.txt" % (directory,folder_name)

    def set_stats(self,directory):
        folder_name = "stats"
        self.stats_dir = "%s/%s" % (directory,folder_name)
        self.plots_dir = "%s/%s/plots" % (directory,folder_name)
        self.bedgraph_dir = "%s/%s/bedgraphs" % (directory,folder_name)
        self.tables_dir = "%s/%s/tables" % (directory,folder_name)
        create_directory("%s/%s" % (directory,folder_name))
        create_directory(self.plots_dir)
        create_directory(self.bedgraph_dir)
        create_directory(self.tables_dir)
        #self.phasing_stats = "%s/%s/phasing_stats.out" % (directory,folder_name)
        #self.assemble_stats_dir = "%s/%s/assembly" % (directory,folder_name)

    def set_tmp(self,directory):
        folder_name = "tmp"
        create_directory("%s/%s" % (directory,folder_name))
        self.tmp_dir = "%s/%s" % (directory,folder_name)
