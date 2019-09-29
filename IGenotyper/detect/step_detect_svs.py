#!/bin/env python
import pysam
import sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from ..common import *
from ..command_line import *


def merge_intervals(intervals):
    intervals.sort()
    merged_intervals = []
    while len(intervals) > 0:
        if len(intervals) == 1:
            merged_intervals.append(intervals[0])
            intervals.pop(0)
            continue
        if intervals[0][1] >= intervals[1][0]:
            tmp = [intervals[0][0],max(intervals[0][1],intervals[1][1])]
            intervals[0] = tmp
            intervals.pop(1)
            continue
        merged_intervals.append(intervals[0])
        intervals.pop(0)
    return merged_intervals

def proportion_of_region_covered(alignment,ref,start,end):
    samfile = pysam.AlignmentFile(alignment)
    intervals = []
    for read in samfile.fetch(ref,start,end):
        if read.is_unmapped:
            continue
        if read.is_supplementary:
            continue
        if read.is_secondary:
            continue
        aligned_blocks = read.get_blocks()
        for aligned_block in aligned_blocks:
            start_block = aligned_block[0]
            end_block = aligned_block[1]
            start_interval = max(start,start_block)
            end_interval = min(end,end_block)
            if end_interval < start_interval:
                continue
            intervals.append([start_interval,end_interval])
    merged_intervals = merge_intervals(intervals)
    sum_of_coverage = 0.0
    for  start_interval, end_interval in merged_intervals:
        sum_of_coverage += (end_interval - start_interval)
    return sum_of_coverage/(end-start)

def number_haplotype_reads(phased_ccs_reads,start,end):
    haplotype_count = {"0": 1.0, "1": 1.0, "2": 1.0}
    samfile = pysam.AlignmentFile(phased_ccs_reads)
    for read in samfile.fetch("igh",start,end):
        if read.is_unmapped:
            continue
        if read.is_supplementary:
            continue
        if read.is_secondary:
            continue
        haplotype = read.get_tag("RG",True)[0]
        haplotype_count[haplotype] += 1.0
    return haplotype_count

def haplotype_coverage(phased_ccs_reads,start,end):    
    haplotype_count = number_haplotype_reads(phased_ccs_reads,start,end)
    total_reads = sum([haplotype_count[h] for h in haplotype_count])
    if total_reads > 10:
        hap_0 = haplotype_count["0"]/total_reads
        hap_1 = haplotype_count["1"]/total_reads
        hap_2 = haplotype_count["2"]/total_reads
    else:
        hap_0 = 0
        hap_1 = 0
        hap_2 = 0
    return (hap_0,hap_1,hap_2)

def get_haplotype_combination(h1_genes,h2_genes):
    gene_order = ["IGHV4-28","IGHV4-30-2","IGHV3-30-3","IGHV4-30-4","IGHV3-30-5","IGHV4-31"]
    haps = []
    for gene in gene_order:
        if gene in h1_genes and gene in h2_genes:
            combinations = ["HOM","HET"]
        elif gene in h1_genes:
            combinations = ["HOM"]
        elif gene in h2_genes:
            combinations = ["HOM"]
        else:
            combinations = ["NP"]
        haps.append(combinations)
    return haps

def hg38_possible_genes():
    genes = ["IGHV4-28","IGHV4-30-2"]
    return genes

def hg19_possible_genes():
    genes = ["IGHV4-28","IGHV4-31"]
    return genes

def na18555_hap_a_genes():
    genes = ["IGHV4-28","IGHV4-30-2","IGHV3-30-3","IGHV4-30-4","IGHV3-30-5","IGHV4-31"]
    return genes

def na18555_hap_b_genes():
    genes = ["IGHV4-28"]
    return genes

def na18507_genes():
    genes = ["IGHV3-30-3","IGHV4-30-4"]
    return genes

def hap_combinations():
    genes_in_haps = {
        "GRCh38":  hg38_possible_genes(),
        "GRCh37": hg19_possible_genes(),
        "NA18555 Hap A": na18555_hap_a_genes(),
        "NA18555 Hap B": na18555_hap_b_genes(),
        "NA18507": na18507_genes()
        }
    hap_combos = {}
    for hap1 in genes_in_haps:
        genes1 = genes_in_haps[hap1]
        for hap2 in genes_in_haps:
            genes2 = genes_in_haps[hap2]
            if (hap2,hap1) in hap_combos:
                continue
            hap_combos[(hap1,hap2)] = get_haplotype_combination(genes1,genes2)
    return hap_combos

def get_possible_haps(genotypes):
    hap_combos = hap_combinations()
    possible_haps = []
    for hap in hap_combos:
        hap_genotypes = hap_combos[hap]
        add = True
        for genotype,hap_genotype in zip(genotypes,hap_genotypes):
            if genotype not in hap_genotype:
                add = False
                break
        if add:
            possible_haps.append(hap)
    return possible_haps
            
def genotype_7_4_1(mapped_locus,phased_reads):    
    coverage_threshold = .90    
    hap_threshold = .30
    sv_start = 157688
    sv_end = 167669
    region_start = 163368
    region_end = 165900
    genotype = "None"
    prop_sv_covered =  proportion_of_region_covered(mapped_locus,"igh",sv_start,sv_end)
    hap_0,hap_1,hap_2 = haplotype_coverage(phased_reads,region_start,region_end)
    if prop_sv_covered > coverage_threshold:
        if min([hap_1,hap_2]) < hap_threshold:
            genotype = "IGHV7-4-1 insertion/No insertion"
        else:
            genotype = "IGHV7-4-1 insertion/IGHV7-4-1 insertion"
    else:
        genotype = "No insertion/No insertion"
    return genotype
        
def genotype_5_10_1(mapped_locus,phased_reads):
    coverage_threshold = .90    
    hap_0_threshold = .90
    sv_start = 210158
    sv_end = 257471
    region_start = 223784
    region_end = 249767
    genotype = "None"
    prop_sv_covered =  proportion_of_region_covered(mapped_locus,"igh",sv_start,sv_end)
    hap_0,hap_1,hap_2 = haplotype_coverage(phased_reads,region_start,region_end)
    if prop_sv_covered > coverage_threshold:
        if hap_0 > hap_0_threshold:
            genotype = "IGHV5-10-1;IGHV3-64D/No haplotype"
        else:
            genotype = "IGHV5-10-1;IGHV3-64D/IGHV5-10-1;IGHV3-64D"
    else:
        genotype = "No haplotype/No haplotype"
    return genotype

def get_read_coverage(alignment,chrom,start,end):
    samfile = pysam.AlignmentFile(alignment)
    read_bases = 0.0
    for read in samfile.fetch(chrom,start,end):
        if read.is_unmapped:
            continue
        if read.is_supplementary:
            continue
        if read.is_secondary:
            continue
        read_bases += read.query_length
    coverage = read_bases/(end-start)
    return coverage

def genotype_3_23(mapped_locus,phased_reads):
    coverage_threshold = .80    
    hap_threshold = .20
    sv_start = 412003
    sv_end = 414775
    region_start = 412003
    region_end = 414775
    genotype = "None"
    prop_sv_covered =  proportion_of_region_covered(mapped_locus,"igh",sv_start,sv_end)
    hap_0,hap_1,hap_2 = haplotype_coverage(phased_reads,region_start,region_end)
    print prop_sv_covered
    print hap_0,hap_1,hap_2     
    if prop_sv_covered > coverage_threshold:
        if min([hap_1,hap_2]) < hap_threshold:
            genotype = "IGHV3-23D duplication/IGHV3-23"
        else:
            genotype = "IGHV3-23D duplication/IGHV3-23D duplication"
    else:
        genotype = "IGHV3-23/IGHV3-23"
    return genotype

def genotype_region(phased_ccs_reads,start,end):
    genotype = None
    haplotype_count = number_haplotype_reads(phased_ccs_reads,start,end)
    total_reads = sum([haplotype_count[h] for h in haplotype_count])
    if total_reads < 20:
        genotype = "NP"
    else:
        hap_0 = haplotype_count["0"]/total_reads
        if hap_0 > .70:
            genotype = "HOM"
        else:
            genotype = "HET"
    return genotype

def genotype_3_30(mapped_locus,phased_ccs_reads):
    V4_28_gene = genotype_region(phased_ccs_reads,474597,474891)
    V4_30_2_gene = genotype_region(phased_ccs_reads,499626,499923)
    V3_30_3_gene = genotype_region(phased_ccs_reads,509356,509650)
    V4_30_4_gene = genotype_region(phased_ccs_reads,523554,523851)
    V3_30_5_gene = genotype_region(phased_ccs_reads,529654,530768)
    V4_31_gene = genotype_region(phased_ccs_reads,548556,548853)
    # V4_28_gene_coverage = get_read_coverage(phased_ccs_reads,"igh",474597,474891)
    # V4_30_2_gene_coverage = get_read_coverage(phased_ccs_reads,"igh",499626,499923)
    # V3_30_3_gene_coverage = get_read_coverage(phased_ccs_reads,"igh",509356,509650)
    # V4_30_4_gene_coverage = get_read_coverage(phased_ccs_reads,"igh",523554,523851)
    # V3_30_5_gene_coverage = get_read_coverage(phased_ccs_reads,"igh",529654,530768)
    # V4_31_gene_coverage = get_read_coverage(phased_ccs_reads,"igh",548556,548853)
    # print V4_28_gene_coverage
    # print V4_30_2_gene_coverage
    # print V3_30_3_gene_coverage
    # print V4_30_4_gene_coverage
    # print V3_30_5_gene_coverage
    # print V4_31_gene_coverage
    genotypes = [V4_28_gene,V4_30_2_gene,V3_30_3_gene,V4_30_4_gene,V3_30_5_gene,V4_31_gene]
    genes = ",".join(["IGHV4-28","IGHV4-30-2","IGHV3-30-3","IGHV4-30-4","IGHV3-30-5","IGHV4-31"])
    gene_genotypes = ",".join([V4_28_gene,V4_30_2_gene,V3_30_3_gene,V4_30_4_gene,V3_30_5_gene,V4_31_gene])
    possible_haps = ";".join(["/".join(i) for i in get_possible_haps(genotypes)])
    return (genes,gene_genotypes,possible_haps)

def genotype_4_38_2(mapped_locus,phased_reads):
    coverage_threshold = .70    
    hap_0_threshold = .90
    sv_start = 616403
    sv_end = 679227
    region_start = 634584
    region_end = 660749
    genotype = "None"
    prop_sv_covered =  proportion_of_region_covered(mapped_locus,"igh",sv_start,sv_end)
    hap_0,hap_1,hap_2 = haplotype_coverage(phased_reads,region_start,region_end)
    if prop_sv_covered > coverage_threshold:
        if hap_0 > hap_0_threshold:
            genotype = "IGHV4-38-2 to IGHV1-38-4 insertion/No insertion"
        else:
            genotype = "IGHV4-38-2 to IGHV1-38-4 insertion/IGHV4-38-2 to IGHV1-38-4 insertion"
    else:
        genotype = "No insertion/No insertion"
    return genotype
    
def genotype_1_69(mapped_locus,phased_reads):
    ## Genotype 1-69 and 2-70 to check for homozygosity
    coverage_threshold = .70
    hap_0_threshold = .70
    sv_start = 966204
    sv_end = 1028309
    region_start = 994215
    region_end = 1000406
    genotype = "None"
    prop_sv_covered =  proportion_of_region_covered(mapped_locus,"igh",sv_start,sv_end)
    hap_0,hap_1,hap_2 = haplotype_coverage(phased_reads,region_start,region_end)
    if prop_sv_covered > coverage_threshold:
        if hap_0 > hap_0_threshold:
            genotype = "IGHV1-69D to IGHV2-70D insertion/IGHV1-69 to IGHV2-70"
        else:
            genotype = "IGHV1-69D to IGHV2-70D insertion/IGHV1-69D to IGHV2-70D insertion"
    else:
        genotype = "IGHV1-69 to IGHV2-70/IGHV1-69 to IGHV2-70"
    return genotype
    
def genotype_1_8(mapped_locus,phased_reads):
    coverage_threshold = .90    
    hap_0_threshold = .90
    sv_start = 1149070
    sv_end = 1193130
    region_start = 1158302
    region_end = 1162596
    genotype = "None"
    prop_sv_covered =  proportion_of_region_covered(mapped_locus,"igh",sv_start,sv_end)
    hap_0,hap_1,hap_2 = haplotype_coverage(phased_reads,region_start,region_end)
    if prop_sv_covered > coverage_threshold:
        if hap_0 > hap_0_threshold:
            genotype = "IGHV1-8;IGHV3-9/No haplotype"
        else:
            genotype = "IGHV1-8;IGHV3-9/IGHV1-8;IGHV3-9"
    else:
        genotype = "No haplotype/No haplotype"
    return genotype

def svs_per_hap(mapped_locus,svs_genotyped,ccs_to_ref):
    output = [["SV","Genotype"]]
    sv_7_4_1_genotype = genotype_7_4_1(mapped_locus,ccs_to_ref)
    sv_5_10_1_genotype = genotype_5_10_1(mapped_locus,ccs_to_ref)
    sv_3_23_genotype = genotype_3_23(mapped_locus,ccs_to_ref)
    sv_3_30_genes,sv_3_30_gene_genotypes,sv_3_30_genotype =  genotype_3_30(mapped_locus,ccs_to_ref)
    sv_3_30_output = "haplotypes=%s;genes=%s;gene_haplotypes=%s" % (sv_3_30_genotype,sv_3_30_genes,sv_3_30_gene_genotypes)
    sv_4_38_2_genotype = genotype_4_38_2(mapped_locus,ccs_to_ref)
    sv_1_69_genotype = genotype_1_69(mapped_locus,ccs_to_ref)    
    sv_1_8_genotype = genotype_1_8(mapped_locus,ccs_to_ref)
    output.append(["IGHV7-4-1 insertion: ",sv_7_4_1_genotype])
    output.append(["IGHV5-10-1/IGHV3-64D haplotype: ",sv_5_10_1_genotype])
    output.append(["IGHV3-23D duplication: ",sv_3_23_genotype])
    output.append(["IGHV4-28 to IGHV3-33 complex event: ",sv_3_30_output])
    output.append(["IGHV4-38-2 to IGHV1-38-4 insertion: ",sv_4_38_2_genotype])
    output.append(["IGHV1-69D to IGHV2-70D insertion: ",sv_1_69_genotype])
    output.append(["IGHV1-8/IGHV3-9 haplotype: ",sv_1_8_genotype])
    with open(svs_genotyped,'w') as fh:
        for out in output:
            fh.write("%s\n" % "\t".join(out))

def detect_variants_type_svs(mapped_locus,svs_genotyped,phase_aligned_bam,ref,ccs,aligned_bam,sv_signature,sv_vcf):
    svs_per_hap(mapped_locus,svs_genotyped,phase_aligned_bam)
