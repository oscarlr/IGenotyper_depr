#!/bin/env python
import os
import sys # temporary
import pysam
import subprocess
import pybedtools
from Bio import SeqIO
from collections import namedtuple

from ..common import *
from ..command_line import *

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

def merge_small_regions(regions,length=1000):
    new_regions = []
    regions = sorted(regions,key=lambda x: x[1])    
    for i,(chrom,start,end,hap) in enumerate(regions):
        new_region = [chrom,start,end,hap]
        if (end - start) < length:
            for chrom2,start2,end2,hap2 in regions[i:]:
                if hap != hap2:
                    continue
                if end == start2:
                    new_region = [chrom,start,end2,hap]
                    break
        new_regions.append(new_region)
    return new_regions            

def show_value(s):
    if sys.version_info.major == 2:
        if isinstance(s, unicode):
            return str(s)
    return s

class AssemblyRun():
    def __init__(self,Sample):
        self.sample = Sample

    def load_whatshap_blocks(self,min_length=500,min_variants=2):
        blocks = []
        Block = namedtuple('Block',['sample','chrom','start_1',
                                    'start','end','num_variants'])
        with open(self.sample.haplotype_blocks,'r') as fh:
            header = fh.readline()
            for line in fh:
                line = line.rstrip().split('\t')
                block = Block._make(line)
                if int(block.num_variants) < min_variants:
                    continue
                if (int(block.end) - int(block.start)) < min_length:
                    continue
                blocks.append([block.chrom,int(block.start),int(block.end)])
        return sorted(blocks,key=lambda x: x[1])

    def get_phased_regions_to_assemble(self):
        whatshap_blocks = pybedtools.BedTool(self.load_whatshap_blocks())
        sv_regions = load_bed_regions(self.sample.sv_regions)
        non_sv_regions = load_bed_regions(self.sample.non_sv_regions)
        sv_regions_to_assemble = get_assemble_regions(whatshap_blocks,sv_regions)
        non_sv_regions_to_assemble = get_assemble_regions(whatshap_blocks,non_sv_regions)
        regions = sv_regions_to_assemble + non_sv_regions_to_assemble
        regions = merge_small_regions(regions)        
        pybedtools.BedTool(regions).saveas(self.sample.regions_to_assemble)

    def create_assembly_scripts(self,flank=2000):     
        regions_to_assemble = pybedtools.BedTool(self.sample.regions_to_assemble)
        bashfiles = []
        for interval in regions_to_assemble:
            chrom = show_value(interval.chrom)
            start = show_value(interval.start)
            end = show_value(interval.end)
            hap = show_value(interval.name)
            if self.sample.add_unphased_reads and hap != "0":
                samtools_hap = "-r %s -r 0" % hap            
            else:
                samtools_hap = "-r %s" % hap            
            directory = "%s/assembly/%s/%s_%s/%s" % (self.sample.outdir,chrom,start,end,hap)
            create_directory(directory)
            if os.path.isfile("%s/done" % directory):
                continue
            bashfile = "%s/assemble.sh" % directory
            params = {
                "hap": samtools_hap,
                "ccs_to_ref": self.sample.phased_ccs_mapped_reads,
                "chrom": chrom,
                "start": max(0,int(start) - flank),
                "end": int(end) + flank,
                "output": directory,
                "threads": self.sample.threads,
                "size": (int(end) - int(start)) + flank*2,
                "subreads": self.sample.input_bam,
                "subreads_to_ref": self.sample.phased_subreads_mapped_reads,
                "python_scripts": self.sample.python_scripts,
                "ref": self.sample.blasr_ref
                }
            write_to_bashfile(self.sample.assembly_script,bashfile,params)
            bashfiles.append(bashfile)
        return bashfiles

    def combine_sequence(self,type_,outfile):
        seqs = []
        bedfh = open(self.sample.regions_to_assemble, 'r')
        for line in bedfh:
            line = line.rstrip().split('\t')
            chrom,start,end,hap = line
            directory = "%s/assembly/%s/%s_%s/%s" % (self.sample.outdir,chrom,start,end,hap)
            contig = "%s/merged_contigs_quivered.%s" % (directory,type_)
            if os.path.isfile(contig):
                contigs = list(SeqIO.parse(contig,type_))
                total_contigs = len(contigs)
                for i,record in enumerate(contigs):
                    record.id = "c=%s:%s-%s_h=%s_i=%s_t=%s_/0/0_0" % (chrom,start,end,hap,i,total_contigs)
                    record.description = ""
                    seqs.append(record)
        SeqIO.write(seqs,outfile,type_)
        bedfh.close()
        
    def combine_sequences(self):
        self.combine_sequence("fastq",self.sample.locus_fastq)
        self.combine_sequence("fasta",self.sample.locus_fasta)

    def __call__(self):
        self.get_phased_regions_to_assemble()
        assembly_scripts = self.create_assembly_scripts()
        run_assembly_scripts(assembly_scripts,self.sample.cluster,self.sample.cluster_walltime,
                             self.sample.cluster_threads,self.sample.cluster_mem,self.sample.cluster_queue)
        self.combine_sequences()

def assemble_reads(self):
    assembly_runner = AssemblyRun(self)
    assembly_runner()
    command_line_tools = CommandLine(self)
    command_line_tools.map_locus()
    
