#!/bin/env python
import sys
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..common import *
from ..command_line import *
from ..assemble.get_assembly_regions import get_phased_regions

def assembly_location(read_name):
    '''
    Return a list of the chrom, start and end 
    of regions that was assembled.
    '''
    read_origin = read_name.split("_")[0].split('=')[1]
    chrom = read_origin.split(":")[0]
    start = int(read_origin.split(":")[1].split("-")[0])
    end = int(read_origin.split(":")[1].split("-")[1])
    return [chrom,start,end]

def is_overlapping(a, b):
    if a[0] != b[0]: # chrom not matching
        return False
    overlapping = False
    num_overlapping = max(0, min(a[2], b[2]) - max(a[1], b[1]))
    if num_overlapping > 0:
        overlapping = True
    return overlapping

def get_haplotype(read_name):
    return read_name.split("_")[1].split('=')[1]

def snps_per_hap(bamfile,reffn,filter_on_region=True):
    samfile = pysam.AlignmentFile(bamfile)
    ref = pysam.FastaFile(reffn)
    regions = {}
    for read in samfile:
        if read.is_unmapped:
            continue
        if read.is_supplementary:
            continue
        if read.is_secondary:
            continue        
        if read.mapping_quality < 10:
            continue
        assembled_region = assembly_location(read.query_name)
        mapped_chrom = samfile.get_reference_name(read.reference_id)
        mapped_start = int(read.reference_start)
        mapped_end = int(read.reference_end)
        mapped_region = [mapped_chrom,mapped_start,mapped_end]
        if filter_on_region:
            if not is_overlapping(assembled_region,mapped_region):
                continue
        haplotype = get_haplotype(read.query_name)
        for read_pos, ref_pos in read.get_aligned_pairs():
            if read_pos == None:
                continue
            if ref_pos == None:
                continue
            if filter_on_region:
                if ref_pos < int(assembled_region[1]):
                    continue
                if ref_pos > int(assembled_region[2]):
                    continue
            ref_base = ref.fetch(mapped_chrom,ref_pos,ref_pos + 1).upper()
            read_base = read.query_sequence[read_pos].upper()
            if ref_base != read_base:
                read_qual = read.query_qualities[read_pos]
                if mapped_chrom not in regions:
                    regions[mapped_chrom] = {}
                if ref_pos not in regions[mapped_chrom]:
                    regions[mapped_chrom][ref_pos] = {}                
                regions[mapped_chrom][ref_pos][haplotype] = (ref_base,read_base,read_qual,read.query_name)
    return regions

def add_genotype_to_snps(haplotype_snps):
    output = []
    for chrom in haplotype_snps:
        for pos in haplotype_snps[chrom]:
            if "haploid" in haplotype_snps[chrom][pos]:
                ref,read_base,read_qual,read_name = haplotype_snps[chrom][pos]["haploid"]
                line = [chrom,int(pos),".",ref,read_base,"./.",read_qual,read_name]
                output.append(line)
            else:
                if "0" in haplotype_snps[chrom][pos]:
                    ref,read_base,read_qual,read_name = haplotype_snps[chrom][pos]["0"]
                    line = [chrom,int(pos),".",ref,read_base,"1/1",read_qual,read_name]
                    output.append(line)
                else:
                    if len(haplotype_snps[chrom][pos]) == 1:
                        hap = None
                        if "1" in haplotype_snps[chrom][pos]:
                            hap = "1"
                        if "2" in haplotype_snps[chrom][pos]:
                            hap = "2"
                        ref,read_base,read_qual,read_name = haplotype_snps[chrom][pos][hap]
                        line = [chrom,int(pos),".",ref,read_base,"0/1",read_qual,read_name]
                        output.append(line)
                    else:
                        hap_1_ref,hap_1_read_base,hap_1_read_qual,hap_1_read_name = haplotype_snps[chrom][pos]["1"]
                        hap_2_ref,hap_2_read_base,hap_2_read_qual,hap_2_read_name = haplotype_snps[chrom][pos]["2"]
                        if hap_1_read_base == hap_2_read_base:
                            genotype = "1/1"
                            alt_base = hap_1_read_base
                        else:
                            genotype = "1/2"
                            alt_base = "%s,%s" % (hap_1_read_base,hap_2_read_base)
                        read_names = "%s,%s" % (hap_1_read_name,hap_2_read_name)
                        alt_qual = "%s" % ((hap_1_read_qual + hap_2_read_qual)/2)
                        line = [chrom,int(pos),".",hap_1_ref,alt_base,genotype,alt_qual,read_names]
                        output.append(line)
    return output

def get_region_type(contigs,sv_regions,non_sv_regions):
    region_types = set()
    for contig in contigs.split(","):
        chrom,start,end = assembly_location(contig)
        for chrom_region,start_region,end_region,region_name in sv_regions:
            if int(start) >= int(start_region):
                if int(end) <= int(end_region):
                    region_types.add(region_name)
        for chrom_region,start_region,end_region in non_sv_regions:
            if int(start) >= int(start_region):
                if int(end) <= int(end_region):
                    region_types.add("nonsv_region")
    return list(region_types)

def snp_in_igh_region(snp_position):
    j_region = ["igh",1,5062,"IGHJ"]
    d_region = ["igh",5062,79255,"IGHD"]
    v_region = ["igh",79255,1193129,"IGHV"]
    regions = [j_region,d_region,v_region]
    for chrom,start,end,region in regions:
        if snp_position > int(start) and snp_position <= int(end):
            return region
    assert snp_position == 1
    return "IGHJ"

def snp_in_gene_feature(snp_position,features):
    in_feature = "No"
    for chrom_region,start_region,end_region,gene in features:
        if snp_position > int(start_region) and snp_position <= int(end_region):
            return gene
    return in_feature

def write_snps_output_to_vcf(genotyped_snps,regions,snps_in_regions,sv_regions,non_sv_regions,snps_from_reads,introns,lpart1_genes,rss_genes,gene_coords,phased_genotypes,hap_blocks):
    print "Calling SNPs %s.." % snps_in_regions
    output_lines = vcf_header()
    genotyped_snps.sort(key=lambda x: x[1])
    with open(snps_in_regions,"w") as vcf_fh:
        vcf_fh.write("%s\n" % "\n".join(output_lines))
        for genotyped_snp in genotyped_snps:
            in_region = False
            for region in regions:
                if int(genotyped_snp[1]) >= int(region[1]):
                    if int(genotyped_snp[1]) <= int(region[2]):
                        in_region = True
                        break
            if not in_region:
                continue
            genotype_with_qual = [":".join(map(str,genotyped_snp[5:-1]))]            
            region_type = get_region_type(genotyped_snp[-1],sv_regions,non_sv_regions)
            if type(genotyped_snp[-2]) != int:
                average_quality = sum(map(int,genotyped_snp[-2].split(",")))/2
            else:
                average_quality = genotyped_snp[-2]
            supported_by_reads = "No"
            if int(genotyped_snp[1]) in snps_from_reads:
                supported_by_reads = "Yes"
            intron = snp_in_gene_feature(genotyped_snp[1],introns)
            lpart1 = snp_in_gene_feature(genotyped_snp[1],lpart1_genes)
            rss = snp_in_gene_feature(genotyped_snp[1],rss_genes)
            in_gene = snp_in_gene_feature(genotyped_snp[1],gene_coords)
            igh_region = snp_in_igh_region(genotyped_snp[1])
            phased_genotype = "."
            if genotyped_snp[1] in phased_genotypes:
                phased_genotype = phased_genotypes[genotyped_snp[1]]
            haplotype_block = snp_in_gene_feature(genotyped_snp[1],hap_blocks)
            if haplotype_block == "No":
                haplotype_block = "."
            info_field = "contig=%s;region=%s;read_support=%s;intronic=%s;LP1=%s;RSS=%s;gene=%s;igh_region=%s;phased_genotype=%s,haplotype_block=%s" % \
                (genotyped_snp[-1],",".join(region_type),supported_by_reads,intron,lpart1,rss,in_gene,igh_region,phased_genotype,haplotype_block)
            append_to_line = [average_quality,"PASS",info_field,"GT:GQ"]
            outline = genotyped_snp[:5] + append_to_line + genotype_with_qual
            vcf_fh.write("%s\n" % "\t".join(map(str,outline)))

def read_vcf_genotype(vcf_file):
    positions = {}
    with open(vcf_file,'r') as vcf_fh:
        for line in vcf_fh:
            if line.startswith("#"):
                continue
            line = line.rstrip().split("\t")
            if line[0] not in ["igh","chr14","14"]:
                continue
            position = line[1]
            if "|" not in line[9]:
                continue
            genotype = line[9].split(":")[0]
            positions[int(position) - 1] = genotype
    return positions

def read_vcf_snps(vcf_file):
    positions = []
    with open(vcf_file,'r') as vcf_fh:
        for line in vcf_fh:
            if line.startswith("#"):
                continue
            line = line.rstrip().split("\t")
            if line[0] not in ["igh","chr14","14"]:
                continue
            positions.append(int(line[1]) - 1)
    return positions

def labeled_hap_blocks(haplotype_blocks):
    phased_regions = get_phased_regions(haplotype_blocks)
    labeled_phased_region = []
    label = 0
    for chrom,start,end in phased_regions:
        labeled_phased_region.append([chrom,start,end,label])
        label += 1
    return labeled_phased_region

def detect_variants_type_snps(sv_regions,non_sv_regions,mapped_locus,pbmm2_ref,snp_candidates,snps_in_sv_regions,snps_in_non_sv_regions,introns,lpart1,rss,gene_coords,phased_variants_vcf,haplotype_blocks):
    sv_regions = load_bed_regions(sv_regions,True)
    non_sv_regions = load_bed_regions(non_sv_regions)
    introns = load_bed_regions(introns,True)
    lpart1 = load_bed_regions(lpart1,True)
    rss = load_bed_regions(rss,True)
    gene_coords = load_bed_regions(gene_coords,True)
    haplotype_snps = snps_per_hap(mapped_locus,pbmm2_ref)
    genotyped_snps = add_genotype_to_snps(haplotype_snps)
    snps_from_reads = read_vcf_snps(snp_candidates)
    phased_snps_from_reads = read_vcf_genotype(phased_variants_vcf)
    haplotype_blocks = labeled_hap_blocks(haplotype_blocks)
    write_snps_output_to_vcf(genotyped_snps,sv_regions,snps_in_sv_regions,
                             sv_regions,non_sv_regions,snps_from_reads,introns,lpart1,rss,
                             gene_coords,phased_snps_from_reads,haplotype_blocks)
    write_snps_output_to_vcf(genotyped_snps,non_sv_regions,snps_in_non_sv_regions,
                             sv_regions,non_sv_regions,snps_from_reads,introns,lpart1,rss,
                             gene_coords,phased_snps_from_reads,haplotype_blocks)
