#!/bin/env python
import pysam
from Bio import SeqIO

def get_n50(fasta_file):
    igh_size = 1193129.0
    assembly_lengths = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        assembly_lengths.append(len(record.seq))
    assembly_lengths = sorted(assembly_lengths,reverse=True)
    sum_of_lengths = 0.0
    n50_length = None
    for length in assembly_lengths:
        sum_of_lengths += length
        if sum_of_lengths >= (igh_size/2):
            n50_length = length
            break
    return n50_length

def get_assembly_size(fasta_file):
    assembly_size = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        assembly_size += len(record.seq)
    return assembly_size

def get_number_of_contigs(fasta_file):
    num_contigs = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        num_contigs += 1
    return num_contigs

def merge(intervals):
    if not intervals:
        return []
    data = []
    for interval in intervals:
        data.append((interval[0], 0))
        data.append((interval[1], 1))
    data.sort()
    merged = []
    stack = [data[0]]
    for i in xrange(1, len(data)):
        d = data[i]
        if d[1] == 0:
            stack.append(d)
        elif d[1] == 1:
            if stack:
                start = stack.pop()
            if len(stack) == 0:
                merged.append( (start[0], d[0]))
    return merged

def get_assembled_locus_coverage(mapped_locus):
    igh_size = 1193129.0
    contig_coords = []
    samfile = pysam.AlignmentFile(mapped_locus)
    for contig in samfile:
        if contig.is_supplementary:
            continue
        if contig.is_secondary:
            continue
        if contig.is_unmapped:
            continue
        ref = samfile.get_reference_name(contig.reference_id)
        if ref == "igh":
            contig_coords.append([contig.reference_start,contig.reference_end])
    contig_coords = merge(contig_coords)
    covered = 0.0
    for c1, c2  in contig_coords:
        covered += (c2 - c1)
    return round((covered/igh_size)*100,2)

def stats_on_assembly(fasta_file,tables_dir,mapped_locus):
    assembly_size = get_assembly_size(fasta_file)
    number_contigs = get_number_of_contigs(fasta_file)
    n50 = get_n50(fasta_file)
    locus_coverage = get_assembled_locus_coverage(mapped_locus)
    data = [["Assembly size",assembly_size],
            ["Number of contigs",number_contigs],
            ["N50",n50],
            ["% of IGH covered",locus_coverage]]
    assemby_stats = "%s/assembly_stats.txt" % tables_dir
    with open(assemby_stats,'w') as fh:
        for column1, column2 in data:
            fh.write("%s\t%s\n" % (column1,column2))                 
