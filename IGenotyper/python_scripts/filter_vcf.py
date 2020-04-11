#!/bin/env python
import sys

region_1_69_to_ignore = ["igh",977000,1022947]

def read_bedfile(bedfile,add_1_69_region):
    regions = []
    with open(bedfile,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            region = [line[0],int(line[1]),int(line[2])]
            regions.append(region)
    if add_1_69_region:
        regions.append(region_1_69_to_ignore)
    return regions

def filter_line(line,regions_to_ignore):
    line = line.rstrip().split('\t')
    line_pos = int(line[1])
    filter_ = False
    for chrom,start,end in regions_to_ignore:
        if line_pos > start and line_pos < end:
            filter_ = True
    return filter_

def filter_vcf(input_variants_vcf,output_variants_vcf,regions_to_ignore_bed,add_1_69_region):
    regions_to_ignore = read_bedfile(regions_to_ignore_bed,add_1_69_region)
    with open(output_variants_vcf,'w') as output_variants_vcf_fh:
        with open(input_variants_vcf,'r') as input_variants_vcf_fh:
            for line in input_variants_vcf_fh:
                if "#" in line:
                    output_variants_vcf_fh.write(line)
                else:
                    if not filter_line(line,regions_to_ignore):
                        output_variants_vcf_fh.write(line)

def get_1_69_hets(input_variants_vcf):
    flank = 2000
    ighv_1_69_2_start = 997521 - flank
    ighv_1_69_2_end = 997815 + flank
    count = 0
    with open(input_variants_vcf,'r') as input_variants_vcf_fh:
        for line in input_variants_vcf_fh:
            if "#" in line:
                continue
            line_columns = line.rstrip().split('\t')
            line_pos = int(line_columns[1])
            if line_pos >  ighv_1_69_2_start and  line_pos < ighv_1_69_2_end:
                if "0/1" in line:
                    count += 1
                elif "1/0" in line:
                    count += 1
    return count

def main():
    input_variants_vcf = sys.argv[1]
    output_variants_vcf = sys.argv[2]
    regions_to_ignore_bed = sys.argv[3]
    num_hets = get_1_69_hets(input_variants_vcf)
    if num_hets < 10:
        add_1_69_region = True
    else:
        add_1_69_region = False
    filter_vcf(input_variants_vcf,output_variants_vcf,regions_to_ignore_bed,add_1_69_region)

if __name__ == "__main__":
    main()
