#!/bin/env python
import subprocess
import datetime
import numpy as np
import pysam
import vcf
import sys

def read_phased_snps(vcffile,sample_name):
    phased_snps = {}
    vcf_reader = vcf.Reader(open(vcffile, 'r'))
    for record in vcf_reader:
        if record.var_subtype == "deletion":
            continue
        if record.var_subtype == "insertion":
            continue
        phased_snps.setdefault(record.CHROM, {})
        alleles  = record.genotype(sample_name)['GT'].split("|")
        if len(alleles) == 2 and alleles[0] != alleles[1]:
            alleles = map(int, alleles)
            if max(alleles) != 1:
                continue
            allele_bases = [record.REF,record.ALT[0]]
            sample_bases = map(lambda x: allele_bases[x], alleles)
            phased_snps[record.CHROM][record.POS-1] = sample_bases
    return phased_snps

def create_header(unphased_bam):
    unphased_tag = 0
    haplotype1_tag = 1
    haplotype2_tag = 2
    readgroup_unphased = { "ID": unphased_tag }
    readgroup_hap1 = { "ID": haplotype1_tag }
    readgroup_hap2 = { "ID": haplotype2_tag }
    if "RG" in unphased_bam.header:
        for group in unphased_bam.header["RG"]:
            for item in group:
                if item != "ID":
                    readgroup_unphased[item] = group[item]
                    readgroup_hap1[item] = group[item]
                    readgroup_hap2[item] = group[item]
    phased_bam_header = unphased_bam.header.to_dict()
    phased_bam_header["RG"] = [readgroup_unphased,readgroup_hap1,readgroup_hap2]
    return phased_bam_header

def create_tag(hap):
    haptag = ("RG", str(hap), "Z")
    return haptag

def calculate_no_prob(basetups,hap_dict): #,threshold=0.9,lr_threshold=10):
    '''
    Returns a 0, 1 or -1 as the value of a read
    '''
    unphased_tag = 0
    haplotype1_tag = 1
    haplotype2_tag = 2
    b1_bases = 0.0
    b0_bases = 0.0
    for basetup in basetups:
        rpos, b, qual = basetup
        b0, b1 = hap_dict[rpos]
        if b == b1 or b == b0: # check to make sure the base is called
            if b == b0:
                b0_bases += 1.0
            else:
                b1_bases += 1.0              
    if (b1_bases + b0_bases) == 0:
        return create_tag(unphased_tag)
    b0_frac = b0_bases/(b1_bases + b0_bases)
    b1_frac = b1_bases/(b1_bases + b0_bases)
    if b0_frac >= .9: # .8 sovles 1-69
        return create_tag(haplotype1_tag)
    elif b1_frac >= .9:
        return create_tag(haplotype2_tag)
    else:
        return create_tag(unphased_tag)

def phase_read(read,phased_snps,chrom,snps_to_ignore):
    unphased_tag = 0
    haplotype1_tag = 1
    haplotype2_tag = 2
    base_tuples = []
    if chrom in phased_snps:
        for read_pos, ref_pos in read.get_aligned_pairs():
            if ref_pos is None:
                continue
            if ref_pos in phased_snps[chrom]:
                if ref_pos in snps_to_ignore:
                    continue
                if read_pos is not None:
                    qual = 8
                    if read.query_qualities != None:
                        qual = read.query_qualities[read_pos]
                    base = read.query_sequence[read_pos]
                    base_tuples.append((ref_pos,base,qual))
    if len(base_tuples) > 0:
        #print "1"
        read_group_tag = calculate_no_prob(base_tuples,phased_snps[chrom])
    else:
        #print "2"
        read_group_tag = create_tag(unphased_tag)
    #print read_group_tag
    read_tags = read.get_tags()
    tags_to_add = []
    for tag in read_tags:
        if tag[0] != "RG":
            tags_to_add.append(tag)
    tags_to_add.append(read_group_tag)
    read.set_tags(tags_to_add)
    return read

def read_bases_overlapping_snps(read,phased_snps,chrom):
    base_tuples = []
    if chrom in phased_snps:
        for read_pos, ref_pos in read.get_aligned_pairs():
            if ref_pos is None:
                continue
            if ref_pos in phased_snps[chrom]:
                if read_pos is not None:
                    qual = read.query_qualities[read_pos]
                    base = read.query_sequence[read_pos]
                    base_tuples.append((ref_pos,base,qual))
    return base_tuples

def main():
    phased_variants_vcf = sys.argv[1]
    mapped_sequencing_data = sys.argv[2]
    phased_mapped_sequencing_data = sys.argv[3]
    sample_name = sys.argv[4]
    phased_snps = read_phased_snps(phased_variants_vcf,sample_name)
    unphased_bam = pysam.AlignmentFile(mapped_sequencing_data,'rb')    
    filtered_phased_snps = []
    phased_bam_header = create_header(unphased_bam)
    phased_bam = pysam.AlignmentFile(phased_mapped_sequencing_data,'wb',header=phased_bam_header)
    for read in unphased_bam.fetch(): 
        if read.is_unmapped:
            continue
        chrom = unphased_bam.get_reference_name(read.reference_id)
        tagged_read = phase_read(read,phased_snps,chrom,filtered_phased_snps)
        phased_bam.write(tagged_read)
    unphased_bam.close()
    phased_bam.close()

if __name__ == "__main__":
    main()
