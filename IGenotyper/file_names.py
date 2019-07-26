#!/bin/env python

class FileNames():
    def __init__(self,directory):

        # Package data
        self.ref = None
        self.genes = None
        self.sv_regions = None
        self.non_sv_regions = None
        
        # Alignments
        self.ccs_aligned_bam = None
        self.subread_aligned_bam = None
        self.phased_ccs_aligned_bam = None
        self.phased_subreads_aligned_bam = None
        self.mapped_locus = None

        # Assembly 
        self.locus_fasta = None
        self.locus_fastq = None
        self.regions_to_assemble = None
        self.haplotype_blocks = None
        
        # Variants
        self.snp_candidates = None
        self.ccs_variants_vcf = None
        self.phased_ccs_variants_vcf = None
        self.snp_in_sv_vcf = None
        self.snp_not_in_sv_vcf = None
        self.svs_txt = None
        self.sv_sig = None
        self.sv_vcf = None

        # Alleles
        self.alleles = None
        self.allele_assignment = None
        self.genes_fasta = None
        self.novel_alleles = None

        # Stats
        self.phasing_stats = None

    def package_data(self):
        directory = os.path.dirname(os.path.abspath(__file__))
        self.pbmm2_ref = "%s/data/pbmm2_index/reference.fasta" % directory
        self.genes = "%s/data/gene_coords.bed" % directory
        self.sv_regions = "%s/data/sv_coords.bed" % directory
        self.non_sv_regions = "%s/data/non_sv_coords.bed" % directory

    def set_alignments(self,directory):
        folder_name = "alignments"
        self.ccs_aligned_bam = "%s/%s/ccs_to_ref.sorted.bam" % (directory,folder_name)
        self.subreads_aligned_bam = "%s/%s/subreads_to_ref.sorted.bam" % (directory,folder_name)
        self.phased_ccs_aligned_bam = "%s/%s/ccs_to_ref_phased.sorted.bam" % (directory,folder_name)
        self.phased_subreads_aligned_bam = "%s/%s/subreads_to_ref_phased.sorted.bam" % (directory,folder_name)
        self.mapped_locus = "%s/%s/locus_to_ref.sorted.bam" % (directory,folder_name)

    def set_assembly(self,directory):
        folder_name = "assembly"
        self.locus_fasta = "%s/%s/locus.fasta" % (directory,folder_name)
        self.locus_fastq = "%s/%s/locus.fastq" % (directory,folder_name)
        self.regions_to_assemble = "%s/%s/ccs_variants_phased.bed" % (directory,folder_name)
        self.haplotype_blocks = "%s/%s/ccs_variants_phased.tab" % (directory,folder_name)

    def set_variants(self,directory):
        folder_name = "variants"
        self.snp_candidates = "%s/%s/from_reads/snp_candidates.vcf" % (directory,folder_name)
        self.ccs_variants_vcf = "%s/%s/from_reads/ccs_variants.vcf" % (directory,folder_name)
        self.phased_ccs_variants_vcf = "%s/%s/from_reads/ccs_phased_variants.vcf" % (directory,folder_name)
        self.snp_in_sv_vcf = "%s/%s/from_assembly/snps_in_svs.vcf" % (directory,folder_name)
        self.snp_not_in_sv_vcf = "%s/%s/from_assembly/snps_not_in_svs.vcf" % (directory,folder_name)
        self.svs_txt = "%s/%s/from_assembly/svs.txt" % (directory,folder_name)
        self.sv_sig = "%s/%s/from_assembly/sv.svsig.gz" % (directory,folder_name)
        self.sv_vcf = "%s/%s/from_assembly/sv.vcf" % (directory,folder_name)
        
    def alleles(self,directory):
        folder_name = "alleles"
        self.alleles = "%s/%s/alleles.tab" % (directory,folder_name)
        self.allele_assignment = "%s/%s/genes_assigned_to_alleles.txt" % (directory,folder_name)
        self.genes_fasta = "%s/%s/genes_from_assembly.fasta" % (directory,folder_name)
        self.novel_alleles = "%s/%s/novel_alleles.txt" % (directory,folder_name)

    def stats(self,directory):
        folder_name = "stats"
        self.phasing_stats = "%s/%s/phasing_stats.out" % (directory,folder_name)
        
