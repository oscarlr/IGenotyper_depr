#!/bin/env python
import os
import shutil
import argparse
import datetime
import pysam
from string import Template


def remove_files(dir_,fns):
    for fn in fns:
        fn_path = "%s/%s" % (dir_,fn)
        if os.path.exists(fn_path):
            os.remove(fn_path)

def remove_dirs(dir_,dirs):
    for dirs_ in dirs:
        dir_path = "%s/%s" % (dir_,dirs_)        
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)

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

def get_haplotype(read_name):
    return read_name.split("_")[1].split('=')[1]

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

def get_gene_names(gene_coords):
    gene_coords = load_bed_regions(gene_coords,True)
    gene_names = [genes[3] for genes in gene_coords]
    return gene_names

def get_window_size(window_size,step,inbed,outbed):
    args = [inbed,window_size,step,outbed]
    command = ("bedtools makewindows -b %s -w %s -s %s > %s \n" % tuple(args))
    if not non_emptyfile("%s" % outbed):
        os.system(command)

def haplotype_coverage(bam,windows_bed,output_coverage):
    args = [bam,windows_bed,output_coverage,
            bam,windows_bed,output_coverage,
            bam,windows_bed,output_coverage]
    command = ("samtools view -Sbh %s -r 0 | bedtools coverage -counts -a %s -b stdin | awk '{ print $0\"\t0\"}' > %s \n"
               "samtools view -Sbh %s -r 1 | bedtools coverage -counts -a %s -b stdin | awk '{ print $0\"\t1\"}' >> %s \n"
               "samtools view -Sbh %s -r 2 | bedtools coverage -counts -a %s -b stdin | awk '{ print $0\"\t2\"}' >> %s \n" % tuple(args))
    if not non_emptyfile("%s" % output_coverage):
        os.system(command)

def read_genotype_svs(svs_genotyped_fn):
    sv_genotypes = {}
    with open(svs_genotyped_fn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            if line[0] == "chrom":
                continue
            sv_name = line[3]
            genotype = line[5]
            sv_genotypes[sv_name] = genotype
    return sv_genotypes

def is_overlapping(a, b):
    if a[0] != b[0]:
        return False
    overlapping = False
    num_overlapping = max(0, min(a[2], b[2]) - max(a[1], b[1]))
    if num_overlapping > 0:
        overlapping = True
    return overlapping

def get_phased_coverage(bamfile,chrom,start,end,read_group):
    coverage = get_coverage(bamfile,chrom,start,end,read_group)    
    return coverage

def get_coverage(bamfile,chrom,start,end,hap=None):
    samfile = pysam.AlignmentFile(bamfile)
    read_bases = 0.0
    for read in samfile.fetch(chrom,start,end):
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        if read.is_unmapped:
            continue
        if hap != None:
            read_group = read.get_tag("RG",True)[0]
            if hap != read_group:
                continue
        read_bases += float(read.query_length)
    coverage = read_bases/(end - start)
    return coverage

def check_if_step_completed(fns,outfn):
    finished = True
    for fn in fns:
        if not non_emptyfile(fn):
            finished = False
    if finished:
        with open(outfn,"w") as fh:
            fh.write("done\n")
