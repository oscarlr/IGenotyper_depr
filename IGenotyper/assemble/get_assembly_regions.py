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

def get_regions_with_hap(bed_file,hap):
    regions = []
    for interval in bed_file:
        chrom = show_value(interval.chrom)
        start = show_value(interval.start)
        end = show_value(interval.end)
        regions.append([chrom,start,end,hap])
    return regions

def get_assemble_regions(phased_regions,regions):
    phased_regions_bed = pybedtools.BedTool(phased_regions)
    regions_bed = pybedtools.BedTool(regions)
    phased_regions_bed = regions_bed.intersect(phased_regions_bed)
    non_phased_regions_bed = regions_bed.subtract(phased_regions_bed)
    regions = get_regions_with_hap(phased_regions_bed,"1")
    regions += get_regions_with_hap(phased_regions_bed,"2")
    regions += get_regions_with_hap(non_phased_regions_bed,"0")
    return regions

def get_regions_to_assemble(phased_regions,sv_regions,non_sv_regions):
    phased_regions = load_phased_regions(phased_regions)
    sv_regions = load_bed_regions(sv_regions)
    non_sv_regions = load_bed_regions(non_sv_regions)
    if len(phased_regions) == 0:
        return sv_regions + non_sv_regions
    phased_regions_bed = pybedtools.BedTool(phased_regions)
    sv_regions_to_assemble = get_assemble_regions(phased_regions_bed,sv_regions)
    non_sv_regions_to_assemble = get_assemble_regions(phased_regions_bed,non_sv_regions)
    return sv_regions_to_assemble + non_sv_regions_to_assemble

