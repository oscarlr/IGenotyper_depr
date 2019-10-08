#!/bin/env python
import os
import datetime
import argparse
from string import Template


def non_emptyfile(checkfile):
    return os.path.isfile(checkfile) and os.path.getsize(checkfile) > 0

def check_file_exist(file_):
    if not non_emptyfile(file_):
        msg = "The file %s does not exist" % file_
        raise argparse.ArgumentTypeError(msg)
    return file_

def create_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def write_to_bashfile(template_bash,bashfile,params):
    filein = open(template_bash)
    src = Template(filein.read())
    output_lines = src.safe_substitute(params)
    bashfh = open(bashfile,'w')
    bashfh.write(output_lines)
    filein.close()
    bashfh.close()

def load_bed_regions(bedfile,add_fourth=False):
    bed_regions = []
    with open(bedfile,'r') as bedfh:
        for line in bedfh:
            line = line.rstrip().split('\t')
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            if add_fourth:
                annotation = line[3]
                bed_regions.append([chrom,start,end,annotation])
            else:
                bed_regions.append([chrom,start,end])
    return bed_regions

def read_bedfile(bedfile):
    regions = {}
    with open(bedfile,'r') as bedfh:
        for line in bedfh:
            line = line.rstrip().split('\t')
            chrom = line[0]
            coord = (int(line[1]),int(line[2]))
            if chrom not in regions:
                regions[chrom] = []
            regions[chrom].append(coord)
    return regions

def vcf_header(sample_name="sample"):
    i = datetime.datetime.now()
    line = [ "##fileformat=VCFv4.2",
             "##fileDate=%s%s%s" % (i.year,i.month,i.day),
             "##source=IGenotyper",
             "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
             "##INFO=<ID=contig,Number=2,Type=String,Description=\"Contig containing SNP\">",
             "##INFO=<ID=region,Number=1,Type=String,Description=\"Type of region\">",
             "##INFO=<ID=read_support,Number=1,Type=String,Description=\"Support from CCS reads\">",
             "##INFO=<ID=intronic,Number=1,Type=String,Description=\"SNP in intron of gene\">",
             "##INFO=<ID=LP1,Number=1,Type=String,Description=\"SNP in leader part 1 sequence of gene\">",
             "##INFO=<ID=RSS,Number=1,Type=String,Description=\"SNP in recombination signal sequence of gene\">",
             "##INFO=<ID=gene,Number=1,Type=String,Description=\"SNP in gene\">",
             "##INFO=<ID=igh_region,Number=1,Type=String,Description=\"SNP in IGHV, IGHD or IGHJ\">",
             "##INFO=<ID=phased_genotype,Number=1,Type=String,Description=\"Phased genotype only in phase in specified haplotype block\">",
             "##INFO=<ID=haplotype_block,Number=1,Type=String,Description=\"Haplotype block containing phased SNPs\">",
             "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",             
             "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % sample_name]
    return line

