#!/bin/env python
from step_detect_snps import detect_variants_type_snps
from step_detect_indels import detect_variants_type_indels
from step_detect_svs import detect_variants_type_svs
from step_detect_alleles import assign_alleles_to_genes

def detect_variants(self):
    detect_variants_type_snps(self.sv_regions,self.non_sv_regions,self.mapped_locus,
                              self.pbmm2_ref,self.snp_candidates,self.snps_in_sv_regions,
                              self.snps_not_in_sv_regions,self.introns,self.lpart1,self.rss,
                              self.gene_coordinates,self.phased_variants_vcf,self.haplotype_blocks)

    assign_alleles_to_genes(self.mapped_locus,self.gene_coordinates,self.genes_from_assembly,
                            self.allele_database,self.novel_alleles,self.genes_with_allele_assignment,
                            self.ccs_mapped_reads,self.genes_from_reads)

    indel_dir = "%s/indel" % self.tmp_dir
    detect_variants_type_indels(self.mapped_locus,indel_dir,self.pbmm2_ref,self.indels,
                                self.sv_regions,self.non_sv_regions,self.python_scripts)

    ccs = "%s/tmp/ccs.bam" % self.outdir
    detect_variants_type_svs(self.mapped_locus,self.svs_genotyped,self.phased_ccs_mapped_reads,
                             self.pbmm2_ref,ccs,self.ccs_mapped_reads,self.sv_signature,self.sv_vcf)

    