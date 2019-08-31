#!/bin/env python
import os
import sys
import pysam
import pybedtools

from ..common import *
from ..command_line import *


def get_phased_regions(whatshap_blocks):
    min_length = 0
    min_variants = 2
    regions = []
    with open(whatshap_blocks,'r') as fh:
        header = fh.readline()
        for line in fh:
            line = line.rstrip().split('\t')
            chrom = line[1]
            start = int(line[3])
            end = int(line[4])
            num_variants = int(line[5])
            if num_variants < min_variants:
                continue
            if (end - start) < min_length:
                continue
            regions.append([chrom,start,end])
    return sorted(regions,key=lambda x: x[1])

def load_phased_regions(phased_blocks):
    regions = []
    if os.path.isfile(phased_blocks):
        regions = get_phased_regions(phased_blocks)
    return regions

def show_value(s):
    if sys.version_info.major == 2:
        if isinstance(s, unicode):
            return str(s)
    return s

def get_intervals_with_coverage(chrom,start,end,chrom_hap_with_coverage):
    feature = pybedtools.BedTool([(chrom, start, end)])
    coverage_regions = []
    for s,e in chrom_hap_with_coverage:
        coverage_regions.append((chrom,s,e))
    coverage_regions = pybedtools.BedTool(coverage_regions)
    feature_with_coverage = feature.intersect(coverage_regions)
    return feature_with_coverage

def get_regions_with_hap(bed_file,hap,regions_with_coverage):
    regions = []
    for interval in bed_file:
        chrom = show_value(interval.chrom)
        start = show_value(interval.start)
        end = show_value(interval.end)
        chrom_hap_with_coverage = regions_with_coverage[interval.chrom][hap]
        intervals_with_coverage = get_intervals_with_coverage(chrom,start,end,chrom_hap_with_coverage)
        for interval_with_coverage in intervals_with_coverage:
            regions.append([show_value(interval_with_coverage.chrom),
                            show_value(interval_with_coverage.start),
                            show_value(interval_with_coverage.end),
                            hap])              
    return regions

def add_unphased_regions(regions_with_coverage):
    regions = []
    for s,e in regions_with_coverage["igh"]["0"]:
        regions.append(["igh",s,e,"0"])
    return regions

def get_assemble_regions(phased_regions_bed,regions,regions_with_coverage):
    regions_bed = pybedtools.BedTool(regions)
    phased_regions_bed = regions_bed.intersect(phased_regions_bed)
    non_phased_regions_bed = regions_bed.subtract(phased_regions_bed)
    regions = get_regions_with_hap(phased_regions_bed,"1",regions_with_coverage)
    regions += get_regions_with_hap(phased_regions_bed,"2",regions_with_coverage)
    #regions += get_regions_with_hap(non_phased_regions_bed,"0",regions_with_coverage)    
    regions += add_unphased_regions(regions_with_coverage)
    return regions

def get_read_groups_regions(phased_ccs_mapped_reads):
    regions = {"0":{},"1":{},"2":{}}
    samfile = pysam.AlignmentFile(phased_ccs_mapped_reads)
    for read in samfile:
        read_group = read.get_tag("RG")    
        chrom = samfile.get_reference_name(read.reference_id)
        ref_start = read.reference_start
        ref_end = read.reference_end
        if chrom not in regions[read_group]:
            regions[read_group][chrom] = []
        regions[read_group][chrom].append((ref_start,ref_end))
    return regions

def merge_regions(regions):
    sorted_by_lower_bound = sorted(regions, key=lambda tup: tup[0])
    merged_regions = []
    for higher in sorted_by_lower_bound:
        if len(merged_regions) == 0:
            merged_regions.append(higher)
            continue
        lower = merged_regions[-1]
        if higher[0] <= lower[1]:
            upper_bound = max(lower[1], higher[1])
            merged_regions[-1] = (lower[0], upper_bound) 
        else:
            merged_regions.append(higher)
    return merged_regions

def is_overlapping(a, b):
    overlapping = False
    num_overlapping = max(0, min(a[1], b[1]) - max(a[0], b[0]))
    if num_overlapping > 0:
        overlapping = True
    return overlapping

def get_region_coverage(regions,start,end):
    sorted_by_lower_bound = sorted(regions, key=lambda tup: tup[0])
    sum_ = 0.0
    for ref_start,ref_end in sorted_by_lower_bound:
        if not is_overlapping([start,end],[ref_start,ref_end]):
            continue
        sum_ += (ref_end - ref_start)
    return sum_/(end - start)

def get_phased_coverage(phased_ccs_mapped_reads,phased_regions_with_coverage):
    phase_merged_regions = {}
    read_group_regions = get_read_groups_regions(phased_ccs_mapped_reads)
    with open(phased_regions_with_coverage,'w') as fh:
        for read_group in read_group_regions:
            for chrom in read_group_regions[read_group]:
                merged_regions = merge_regions(read_group_regions[read_group][chrom])
                for start,end in merged_regions:
                    if get_region_coverage(read_group_regions[read_group][chrom],start,end) < 10:
                        continue
                    if chrom not in phase_merged_regions:
                        phase_merged_regions[chrom] = {}
                    if read_group not in phase_merged_regions[chrom]:
                        phase_merged_regions[chrom][read_group] = []
                    phase_merged_regions[chrom][read_group].append((start,end))
                    out = [chrom,start,end,read_group]
                    fh.write("%s\n" % "\t".join(map(str,out)))
    return phase_merged_regions

def get_regions_to_assemble(phased_ccs_mapped_reads,phased_regions_with_coverage,phased_regions,sv_regions,non_sv_regions):
    phase_merged_regions = get_phased_coverage(phased_ccs_mapped_reads,phased_regions_with_coverage)
    phased_regions = load_phased_regions(phased_regions)
    sv_regions = load_bed_regions(sv_regions)
    non_sv_regions = load_bed_regions(non_sv_regions)
    if len(phased_regions) == 0:
        return sv_regions + non_sv_regions
    phased_regions_bed = pybedtools.BedTool(phased_regions)
    sv_regions_with_haps = get_assemble_regions(phased_regions_bed,sv_regions,phase_merged_regions)
    non_sv_regions_with_haps = get_assemble_regions(phased_regions_bed,non_sv_regions,phase_merged_regions)
    return sv_regions_with_haps + non_sv_regions_with_haps

