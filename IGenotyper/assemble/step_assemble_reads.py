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
# from get_assembly_regions import get_regions_to_assemble

# def create_assembly_scripts(regions_to_assemble,assembly_dir,mapped_bam,reads_bam,
#                             threads,raw_mapped_reads,add_unphased_reads,assembly_script,python_scripts,ref):
#     flank = 2000
#     bashfiles = []
#     for region in regions_to_assemble:
#         chrom = region[0]
#         start = region[1]
#         end = region[2]
#         if len(region) == 4:
#             hap = region[3]
#             if add_unphased_reads and hap != "0":
#                 samtools_hap = "-r %s -r 0" % hap            
#             else:
#                 samtools_hap = "-r %s" % hap            
#         if len(region) == 3:
#             hap = "haploid"
#             samtools_hap = ""        
#         directory = "%s/%s/%s_%s/%s" % (assembly_dir,chrom,start,end,hap)
#         create_directory(directory)
#         if os.path.isfile("%s/done" % directory):
#             continue
#         template_bash = assembly_script
#         bashfile = "%s/assemble.sh" % directory
#         params = {
#             "hap": samtools_hap,
#             "subreads_to_ref": mapped_bam,
#             "chrom": chrom,
#             "start": max(0,int(start) - flank),
#             "end": int(end) + flank,
#             "output": directory,
#             "threads": threads,
#             "size": (int(end) - int(start)) + flank*2,
#             "subreads": reads_bam,
#             "raw_subreads_to_ref": raw_mapped_reads,
#             "python_scripts": python_scripts,
#             "ref": ref
#             }
#         write_to_bashfile(template_bash,bashfile,params)
#         bashfiles.append(bashfile)
#     return bashfiles

# def combine_sequence(outdir,bedfile,outfasta,outfastq):
#     fasta_seqs = []
#     fastq_seqs = []
#     bedfh = open(bedfile, 'r')
#     for line in bedfh:
#         line = line.rstrip().split('\t')
#         if len(line) == 4:
#             chrom,start,end,hap = line
#         else:
#             chrom,start,end = line
#             hap = "haploid"
#         directory = "%s/%s/%s_%s/%s" % (outdir,chrom,start,end,hap)
#         contig_fasta = "%s/merged_contigs_quivered.fasta" % directory
#         if os.path.isfile(contig_fasta):
#             contigs = list(SeqIO.parse(contig_fasta,"fasta"))
#             total_contigs = len(contigs)
#             for i,record in enumerate(contigs):
#                 record.id = "coord=%s:%s-%s_hap=%s_index=%s_total=%s_/0/0_0" % (chrom,start,end,hap,i,total_contigs)
#                 record.description = ""
#                 fasta_seqs.append(record)
#         contig_fastq = "%s/merged_contigs_quivered.fastq" % directory
#         if os.path.isfile(contig_fastq):
#             contigs = list(SeqIO.parse(contig_fastq,"fastq"))
#             total_contigs = len(contigs)
#             for i,record in enumerate(contigs):
#                 record.id = "coord=%s:%s-%s_hap=%s_index=%s_total=%s_/0/0_0" % (chrom,start,end,hap,i,total_contigs)
#                 record.description = ""
#                 fastq_seqs.append(record)
#     SeqIO.write(fasta_seqs, outfasta, "fasta")
#     SeqIO.write(fastq_seqs, outfastq, "fastq")
#     bedfh.close()

# def combine_headers(headers,ref):
#     inref = pysam.FastaFile(ref)
#     inref_lengths = inref.lengths
#     inref_references = inref.references
#     inref_references_lengths = {}
#     for length, reference in zip(inref_lengths,inref_references):
#         inref_references_lengths[reference] = length
#     SQ_line = headers[0]["SQ"]
#     SQ_line[0]['SN'] = SQ_line[0]['SN'].split(":")[0]
#     SQ_line[0]['LN'] = inref_references_lengths[SQ_line[0]['SN']]
#     if 'M5' in SQ_line[0]:
#         del SQ_line[0]['M5']
#     del headers[0]["RG"][0]["PU"]
#     outheader = {}
#     outheader["SQ"] = SQ_line
#     outheader["HD"] = headers[0]["HD"]
#     outheader["RG"] = headers[0]["RG"]    
#     outheader["PG"] = []
#     for i,header in enumerate(headers):
#         add_header = header["PG"][0]
#         add_header["ID"] = i
#         outheader["PG"].append(add_header)
#     return outheader

# def combine_alignment(outdir,bedfile,outbam,ref):
#     outbam_tmp = "%s/locus_to_ref.bam" % outdir
#     headers = []
#     bedfh = open(bedfile, 'r')
#     for line in bedfh:
#         line = line.rstrip().split('\t')
#         if len(line) == 4:
#             chrom,start,end,hap = line
#         else:
#             chrom,start,end = line
#             hap = "haploid"
#         directory = "%s/%s/%s_%s/%s" % (outdir,chrom,start,end,hap)
#         insam = "%s/contig_after_filter_to_ref.bam" % directory
#         if not os.path.isfile(insam):
#             continue
#         sam = pysam.AlignmentFile(insam)
#         headers.append(dict(sam.header))    
#     bedfh.close()
#     outheader = combine_headers(headers,ref)
#     outbamfh = pysam.AlignmentFile(outbam_tmp,'wb',header=outheader)
#     outbam_header = outbamfh.header    
#     alignments = []
#     bedfh = open(bedfile, 'r')
#     for line in bedfh:
#         line = line.rstrip().split('\t')
#         if len(line) == 4:
#             chrom,start,end,hap = line
#         else:
#             chrom,start,end = line
#             hap = "haploid"
#         directory = "%s/%s/%s_%s/%s" % (outdir,chrom,start,end,hap)
#         insam = "%s/contig_after_filter_to_ref.bam" % directory
#         if not os.path.isfile(insam):
#             continue
#         contig_fasta = "%s/merged_contigs_quivered.fasta" % directory
#         name_mapping = {}
#         if os.path.isfile(contig_fasta):
#             contigs = list(SeqIO.parse(contig_fasta,"fasta"))
#             total_contigs = len(contigs)
#             for i,record in enumerate(contigs):
#                 name_mapping[record.id] = "coord=%s:%s-%s_hap=%s_index=%s_total=%s_/0/0_0" % (chrom,start,end,hap,i,total_contigs)
#         print contig_fasta
#         print insam
#         sam = pysam.AlignmentFile(insam,header=outheader)
#         for alignment in sam.fetch():
#             out_alignment = alignment.to_dict()
#             offset = int(out_alignment["ref_name"].split(":")[1].split("-")[0])
#             next_ref_pos = max(1,int(out_alignment["next_ref_pos"]) + offset)
#             ref_pos = max(1,int(out_alignment["ref_pos"]) + offset)
#             if int(start) != 1:
#                 next_ref_pos = next_ref_pos - 1
#                 ref_pos = ref_pos - 1
#             out_alignment["next_ref_pos"] = str(next_ref_pos)
#             out_alignment["ref_pos"] = str(ref_pos)
#             out_alignment["ref_name"] = out_alignment["ref_name"].split(":")[0]
#             print out_alignment["name"].split("/")
#             print name_mapping
#             out_alignment["name"] = name_mapping[out_alignment["name"].split("/")[0]]
#             out_alignment = pysam.AlignedSegment.from_dict(out_alignment,outbam_header)
#             outbamfh.write(out_alignment)
#     outbamfh.close()
#     bedfh.close()
#     pysam.sort("-o",outbam,outbam_tmp)
#     pysam.index(outbam)
#     outbamfh.close()

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
            directory = "%s/assembly/%s/%s_%s/%s" % (self.outdir,chrom,start,end,hap)
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
            write_to_bashfile(self.assembly_script,bashfile,params)
            bashfiles.append(bashfile)
        return bashfiles

    def __call__(self):
        self.get_phased_regions_to_assemble()
        assembly_scripts = self.create_assembly_scripts()
        run_assembly_scripts(assembly_scripts,self.cluster,self.cluster_walltime,
                             self.cluster_threads,self.cluster_mem,self.cluster_queue)
        

def assemble_reads(self):
    assembly_runner = AssemblyRun(self)
    assembly_runner()

 #    assembly_dir = "%s/assembly" % self.outdir    
#     regions_to_assemble = get_regions_to_assemble(self.haplotype_blocks,self.sv_regions,self.non_sv_regions)

#     pybedtools.BedTool(regions_to_assemble).saveas(self.regions_to_assemble)    
#     assembly_scripts = create_assembly_scripts(regions_to_assemble,assembly_dir,self.phased_ccs_mapped_reads,
#                                                self.input_bam,self.threads,self.phased_subreads_mapped_reads,
#                                                self.add_unphased_reads,self.assembly_script,self.python_scripts,self.pbmm2_ref)    
#     combine_sequence(assembly_dir,self.regions_to_assemble,self.locus_fasta,self.locus_fastq)
#     map_reads(1,self.igh_fasta,self.locus_fastq,self.mapped_locus,"ccs")
    #align_assembly()
    #combine_alignment(assembly_dir,self.regions_to_assemble,self.mapped_locus,self.pbmm2_ref)
    
