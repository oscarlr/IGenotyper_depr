#!/bin/env python
import pysam
import sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from ..common import *
from ..command_line import *


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

def get_region_type(contigs,sv_regions,non_sv_regions):
    region_types = set()
    for contig in contigs.split(","):
        chrom,start,end = assembly_location(contig)
        for chrom_region,start_region,end_region,region_type in sv_regions:
            if int(start) >= int(start_region):
                if int(end) <= int(end_region):
                    region_types.add(region_type)
        for chrom_region,start_region,end_region in non_sv_regions:
            if int(start) >= int(start_region):
                if int(end) <= int(end_region):
                    region_types.add("nonsv_region")
    return list(region_types)

def indel_in_igh_region(position):
    j_region = ["igh",0,5062,"IGHJ"]
    d_region = ["igh",5062,79255,"IGHD"]
    v_region = ["igh",79255,1193129,"IGHV"]
    regions = [j_region,d_region,v_region]
    for chrom,start,end,region in regions:
        if position > int(start) and position <= int(end):
            return region
    assert position == 0
    return "IGHJ"

def split_mapped_hap0_locus_into_haps(mapped_locus,indel_dir):
    samfile = pysam.AlignmentFile(mapped_locus)
    hap1_bam = "%s/hap1_to_ref.sorted.bam" % indel_dir
    hap2_bam = "%s/hap2_to_ref.sorted.bam" % indel_dir
    hap1_sam = pysam.AlignmentFile(hap1_bam,"wb", template=samfile)
    hap2_sam = pysam.AlignmentFile(hap2_bam,"wb", template=samfile)
    for read in samfile:
        if "hap=0" in read.query_name:
            hap1_sam.write(read)
            hap2_sam.write(read)
    hap1_sam.close()
    hap2_sam.close()
    pysam.index(hap1_bam)
    pysam.index(hap2_bam)
 
def split_mapped_locus_into_haps(mapped_locus,indel_dir):
    samfile = pysam.AlignmentFile(mapped_locus)
    hap1_bam = "%s/hap1_to_ref.sorted.bam" % indel_dir
    hap2_bam = "%s/hap2_to_ref.sorted.bam" % indel_dir
    hap1_sam = pysam.AlignmentFile(hap1_bam,"wb", template=samfile)
    hap2_sam = pysam.AlignmentFile(hap2_bam,"wb", template=samfile)
    for read in samfile:
        if "haploid" in read.query_name:
            hap1_sam.write(read)
            hap2_sam.write(read)
        if "hap=1" in read.query_name:
            hap1_sam.write(read)
        if "hap=2" in read.query_name:
            hap2_sam.write(read)
    hap1_sam.close()
    hap2_sam.close()
    pysam.index(hap1_bam)
    pysam.index(hap2_bam)

def get_hap_sequence(hap_bamfile,regions):
    samfile = pysam.AlignmentFile(hap_bamfile,'rb')
    hap_sequences = {}
    for chrom in regions:
        for start,end in regions[chrom]:
            for contig in samfile.fetch(chrom,start,end):
                if contig.is_unmapped:
                    continue
                if contig.is_secondary:
                    continue
                if contig.is_supplementary:
                    continue
                if contig.reference_start > start:
                    continue
                if contig.reference_end < end:
                    continue
                aligned_pairs = contig.get_aligned_pairs()
                query_start = None
                query_end = None
                matched_ref_start = None
                matched_ref_end = None
                for query_pos, ref_pos in aligned_pairs:
                    if query_pos == None:
                        continue
                    if ref_pos == None:
                        continue
                    if int(ref_pos) <= int(start):
                        query_start = query_pos
                        matched_ref_start = ref_pos
                    query_end = query_pos
                    matched_ref_end = ref_pos
                    if int(ref_pos) > int(end):
                        break
                assert query_start != None
                assert query_end != None
                hap_sequence = contig.query_sequence[query_start:query_end]
                sequence_qual = contig.query_qualities[query_start:query_end]
                assert len(sequence_qual) == len(hap_sequence)
                hap_sequences[(chrom,start,end)] = (hap_sequence,sequence_qual,contig.query_name)
    return hap_sequences

def extract_msa_sequence(indel_dir,ref):
    msa_coordsfn = "%s/msa_coords.bed" % indel_dir
    msa_coords = read_bedfile(msa_coordsfn)
    hap1_bamfn = "%s/hap1_to_ref.sorted.bam" % indel_dir
    hap1_sequence = get_hap_sequence(hap1_bamfn,msa_coords)
    hap2_bamfn = "%s/hap2_to_ref.sorted.bam" % indel_dir
    hap2_sequence = get_hap_sequence(hap2_bamfn,msa_coords)
    fasta = pysam.FastaFile(ref)
    for chrom in msa_coords:
        for i,(start,end) in enumerate(msa_coords[chrom]):
            ref_seq = fasta.fetch(reference=chrom,start=max(1,start),end=end)
            h1_seq,h1_qual,h1_contig = hap1_sequence[(chrom,start,end)]
            h2_seq,h2_qual,h2_contig = hap2_sequence[(chrom,start,end)]
            directory = "%s/%s/%s_%s" % (indel_dir,chrom,start,end)
            create_directory(directory)
            outseqfn = "%s/seq.fa" % directory
            outqualfn = "%s/seq.qual" % directory
            with open(outseqfn,'w') as outfh:
                outfh.write(">ref\n%s\n" % ref_seq)
                outfh.write(">hap1\n%s\n" % h1_seq)
                outfh.write(">hap2\n%s\n" % h2_seq)
            with open(outqualfn,'w') as outqualfh:
                outqualfh.write(">hap1_%s\n%s\n" % (h1_contig,",".join(map(str,h1_qual))))
                outqualfh.write(">hap2_%s\n%s\n" % (h2_contig,",".join(map(str,h2_qual))))

def filter_indel(indel_chrom,indel_start,indel_end,indel_genotype,hap_1_region_type,hap_2_region_type,sv_genotypes):
    filter_ = False
    hap1_svs = hap_1_region_type.split(',')
    hap2_svs = hap_2_region_type.split(',')
    if len(set(hap1_svs)) == 1 and len(set(hap2_svs)) == 1:
        if "nonsv_region" in hap1_svs and "nonsv_region" in hap2_svs:
            return filter_
    hap_svs = hap1_svs + hap2_svs
    for sv in hap_svs:
        if sv_genotypes[sv] == "0/0":
            if sv == "IGHV3-23.region.ABC9":
                if is_overlapping([indel_chrom,indel_start,indel_end],["igh",410643,426627]):
                    filter_ = True
            if sv == "IGHV1-69.region.CH17":
                if is_overlapping([indel_chrom,indel_start,indel_end],["igh",983547,1023544]):
                    filter_ = True
        if sv_genotypes[sv] == "0/1":
            if sv == "IGHV1-69.region.CH17":
                if is_overlapping([indel_chrom,indel_start,indel_end],["igh",983547,1023544]):
                    if indel_genotype == "0/1":
                        filter_ = True
    return filter_

def merge_indels(indel_dir,indels,sv_regions,non_sv_regions,sv_genotypes):
    msa_coordsfn = "%s/msa_coords.bed" % indel_dir
    #msa_coords = read_bedfile(msa_coordsfn)    
    with open(indels,'w') as fh:
        with open(msa_coordsfn) as msa_coordsfh:
            for line in msa_coordsfh:
                line = line.rstrip().split('\t')
                chrom = line[0]
                start = line[1]
                end = line[2]
                hap_1_contig = line[3]
                hap_2_chrom = line[4]
                hap_2_start = line[5]
                hap_2_end = line[6]
                hap_2_contig = line[7]
                hap_1_region_type = ",".join(get_region_type(hap_1_contig,sv_regions,non_sv_regions))
                hap_2_region_type = ",".join(get_region_type(hap_2_contig,sv_regions,non_sv_regions))
                igh_region = indel_in_igh_region(int(start))
                directory = "%s/%s/%s_%s" % (indel_dir,chrom,start,end)    
                infile = "%s/svs.bed" % directory
                with open(infile) as infilefh:
                    for inline in infilefh:
                        inline = inline.rstrip().split('\t')
                        inline += [hap_1_contig,hap_2_contig,hap_1_region_type,hap_2_region_type,igh_region]
                        indel_chrom = inline[0]
                        indel_start = int(inline[1])
                        indel_end = int(inline[2])
                        indel_genotype = inline[4]
                        if filter_indel(indel_chrom,indel_start,indel_end,indel_genotype,hap_1_region_type,hap_2_region_type,sv_genotypes):
                            continue                
                        fh.write("%s\n" % "\t".join(inline))

def detect_indel_in_haps(mapped_locus,indel_dir,ref,package_python_directory,indels,sv_regions,non_sv_regions,sv_genotypes,hap0=False):
    create_directory(indel_dir)
    if hap0:
        split_mapped_hap0_locus_into_haps(mapped_locus,indel_dir)
    else:
        split_mapped_locus_into_haps(mapped_locus,indel_dir)
    get_msa_coordinates(indel_dir,ref,package_python_directory)
    extract_msa_sequence(indel_dir,ref)
    calls_svs_from_msa(indel_dir,package_python_directory)
    merge_indels(indel_dir,indels,sv_regions,non_sv_regions,sv_genotypes)    

def add_to_indels_fh(fh,infile):
    with open(infile,'r') as infile_fh:
        for inline in infile_fh:
            inline = inline.strip().split('\t')
            sv_chrom = inline[0]
            sv_start = int(inline[1])
            sv_end = int(inline[2])
            contig_chrom,contig_start,contig_end = assembly_location(inline[6])                
            if sv_start < contig_start:
                continue
            if sv_end > contig_end:
                continue
            contig_chrom,contig_start,contig_end = assembly_location(inline[7])
            if sv_start < contig_start:
                continue
            if sv_end > contig_end:
                continue
            fh.write("%s\n" % "\t".join(inline))    

def combine_indels(haps_indels,hap0_indels,indels): #,haps_coords):
    #hap_coords = load_bed_regions(haps_coords)
    with open(indels,'w') as indels_fh:
        add_to_indels_fh(indels_fh,haps_indels)
        add_to_indels_fh(indels_fh,hap0_indels)

def detect_variants_type_indels(mapped_locus,indel_dir,ref,indels_prefix,sv_regions,non_sv_regions,python_scripts,svs_genotyped_fn):
    sv_regions = load_bed_regions(sv_regions,True)
    non_sv_regions = load_bed_regions(non_sv_regions)
    sv_genotypes = read_genotype_svs(svs_genotyped_fn)
    ###
    haps_indel_dir = "%s/haps" % indel_dir
    haps_indels = "%s_in_haps.bed" % indels_prefix    
    detect_indel_in_haps(mapped_locus,haps_indel_dir,ref,python_scripts,haps_indels,sv_regions,non_sv_regions,sv_genotypes)
    ###
    hap0_indel_dir = "%s/hap0" % indel_dir
    hap0_indels = "%s_in_hap0.bed" % indels_prefix
    detect_indel_in_haps(mapped_locus,hap0_indel_dir,ref,python_scripts,hap0_indels,sv_regions,non_sv_regions,sv_genotypes,True)
    ###
    indels = "%s.bed" % indels_prefix
    haps_coords = "%s/msa_coords.bed" % haps_indel_dir
    combine_indels(haps_indels,hap0_indels,indels)#,haps_coords)
