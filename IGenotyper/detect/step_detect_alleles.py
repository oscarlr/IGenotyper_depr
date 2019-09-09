#!/bin/env python
import pysam
import sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from ..common import *
from ..command_line import *

def extract_sequence_from(read,chrom,start,end):
    read_start = None
    read_end = None
    aligned_pairs = read.get_aligned_pairs()
    for query_pos, ref_pos in aligned_pairs:
        if query_pos == None:
            continue
        if ref_pos == None:
            continue
        if ref_pos <= start:
            read_start = query_pos
        read_end = query_pos
        if ref_pos > end:
            break
    assert read_start != None
    assert read_end != None
    return read.query_sequence[read_start:read_end].upper()

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

def write_gene_sequence_to_file(sequences,gene_sequence):
    out_sequences = []
    for gene in sequences:
        for haplotype in sequences[gene]:
            total = len(sequences[gene][haplotype])
            for i,seq in enumerate(sequences[gene][haplotype]):
                name = "gene=%s.hap=%s.index=%s.total=%s" % (gene,haplotype,i,total)
                record = SeqRecord(Seq(seq),id=name,name=name,description="")
                out_sequences.append(record)
    SeqIO.write(out_sequences,gene_sequence,"fasta")            
    
def is_overlapping(a, b):
    if a[0] != b[0]:
        return False
    overlapping = False
    num_overlapping = max(0, min(a[2], b[2]) - max(a[1], b[1]))
    if num_overlapping > 0:
        overlapping = True
    return overlapping

def read_overlap_region(read,gene_coord):
    assembled_region = assembly_location(read.query_name)
    if is_overlapping(assembled_region,gene_coord):
        return True
    return False

def extract_genes_from_assembly(mapped_locus,gene_coords,gene_sequencefn):
    sequences = {}
    gene_coords = load_bed_regions(gene_coords,True)
    samfile = pysam.AlignmentFile(mapped_locus)
    for chrom,start,end,gene in gene_coords:
        start = int(start)
        end = int(end)
        for read in samfile.fetch(chrom,start,end):
            mapping_cord = [chrom,read.reference_start,read.reference_end]
            if not read_overlap_region(read,mapping_cord):
                continue
            gene_coord = [chrom,start,end]
            if not read_overlap_region(read,gene_coord):
                continue
            gene_sequence = extract_sequence_from(read,chrom,start,end)
            haplotype = get_haplotype(read.query_name)
            if gene not in sequences:
                sequences[gene] = {}
            if haplotype not in sequences[gene]:
                sequences[gene][haplotype] = []
            sequences[gene][haplotype].append(gene_sequence)
    write_gene_sequence_to_file(sequences,gene_sequencefn)

def get_matches(query_entries,database_entries):
    matches = {}
    for query_entry in query_entries:
        query_seq = query_entries[query_entry].seq.upper()
        matches[(query_entry,query_seq)] = set()
        for database_entry in database_entries:
            database_seq = database_entries[database_entry].seq.upper()
            if query_seq.count(database_seq) > 0:
                matches[(query_entry,query_seq)].add(database_entry)
            if query_seq.count(database_seq.reverse_complement()) > 0:
                matches[(query_entry,query_seq)].add(database_entry)
    return matches

def match_gene_to_allele_db(gene_sequence,database):
    query_entries = SeqIO.to_dict(SeqIO.parse(gene_sequence,"fasta"))
    database_entries = SeqIO.to_dict(SeqIO.parse(database,"fasta"))
    matches = get_matches(query_entries,database_entries)
    novel_genes = []
    genes_to_alleles = {}
    for gene_extraction_name, gene_seq in matches:
        gene_name = gene_extraction_name.split(".")[0].split("=")[1]
        haplotype = gene_extraction_name.split(".")[1].split("=")[1]
        if gene_name not in genes_to_alleles:
            genes_to_alleles[gene_name] = {}
        if haplotype not in genes_to_alleles[gene_name]:
            genes_to_alleles[gene_name][haplotype] = set()
        if len(matches[(gene_extraction_name, gene_seq)]) == 0:
            novel_genes.append([gene_name,gene_seq])
        else:
            for allele_hit in matches[(gene_extraction_name, gene_seq)]:
                genes_to_alleles[gene_name][haplotype].add(allele_hit)    
    return genes_to_alleles,novel_genes

def get_gene_names(gene_coords):
    gene_coords = load_bed_regions(gene_coords,True)
    gene_names = [genes[3] for genes in gene_coords]
    return gene_names

def get_output_assignment(genes_to_alleles,gene_names):
    output_lines = []
    header = ["gene_name","haplotype_0","haplotype_1","haplotype_2"]
    output_lines.append(header)
    for gene_name in genes_to_alleles:
        if gene_name in gene_names:
            gene_names.remove(gene_name)
        if "haploid" in genes_to_alleles[gene_name]:            
            alleles = ",".join(genes_to_alleles[gene_name]["haploid"])
            output_lines.append([gene_name,alleles])
        else:
            out = [gene_name]
            for i in ["0","1","2"]:
                if i in genes_to_alleles[gene_name]:
                    alleles = "Novel"
                    if len(list(genes_to_alleles[gene_name][i])) != 0:
                        alleles = ",".join(list(genes_to_alleles[gene_name][i]))
                    out.append(alleles)
                else:
                    out.append(".")
            output_lines.append(out)
    for gene_name in gene_names:
        out = [gene_name,"NA","NA","NA"]
        output_lines.append(out)
    return output_lines

def assign_alleles_to_genes(mapped_locus,gene_coords,gene_sequence,database,novel,assignment):
    # Save gene sequence from assembly into a fasta file
    extract_genes_from_assembly(mapped_locus,gene_coords,gene_sequence)
    # Assign gene to allele
    genes_to_alleles,novel_genes = match_gene_to_allele_db(gene_sequence,database)
    # Get gene names
    gene_names = get_gene_names(gene_coords)
    # Write allele assignment to file
    output_lines = get_output_assignment(genes_to_alleles,gene_names)
    with open(assignment,'w') as fh:
        for line in output_lines:
            fh.write("%s\n" % "\t".join(line))
    with open(novel,'w') as fh:        
        for gene_name,gene_seq in novel_genes:        
            fh.write("%s\t%s\n" % (gene_name,gene_seq))
