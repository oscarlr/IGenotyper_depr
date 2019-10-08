
#!/bin/env python
import sys
import copy
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
from Bio.SeqRecord import SeqRecord

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
    # assert read_start != None
    # assert read_end != None
    if read_start == None:
        return ""
    if read_end == None:
        return ""
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

def extract_genes_from_sequence(mapped_locus,gene_coords,gene_sequencefn,assembly=False):
    sequences = {}
    gene_coords = load_bed_regions(gene_coords,True)
    samfile = pysam.AlignmentFile(mapped_locus)
    for chrom,start,end,gene in gene_coords:
        start = int(start)
        end = int(end)
        for read in samfile.fetch(chrom,start,end):
            mapping_cord = [chrom,read.reference_start,read.reference_end]
            if assembly:
                if not read_overlap_region(read,mapping_cord):
                    continue
            gene_coord = [chrom,start,end]
            if assembly:
                if not read_overlap_region(read,gene_coord):
                    continue
            gene_sequence = extract_sequence_from(read,chrom,start,end)
            if gene_sequence == "":
                continue
            if assembly:
                haplotype = get_haplotype(read.query_name)
            else:
                haplotype = "reads"
            if gene not in sequences:
                sequences[gene] = {}
            if haplotype not in sequences[gene]:
                sequences[gene][haplotype] = []
            sequences[gene][haplotype].append(gene_sequence)
    write_gene_sequence_to_file(sequences,gene_sequencefn)

def extract_genes_from_assembly(mapped_locus,gene_coords,gene_sequencefn):
    extract_genes_from_sequence(mapped_locus,gene_coords,gene_sequencefn,True)

def extract_genes_from_ccs_reads(mapped_ccs,gene_coords,genes_seq_from_ccs):
    extract_genes_from_sequence(mapped_ccs,gene_coords,genes_seq_from_ccs)

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
            genes_to_alleles[gene_name][haplotype] = []
        if len(matches[(gene_extraction_name, gene_seq)]) == 0:
            novel_genes.append([gene_name,gene_seq])
        else:
            for allele_hit in matches[(gene_extraction_name, gene_seq)]:
                if "IGHV1-69D" != gene_name:
                    if allele_hit.split("_")[0].split("=")[1] != gene_name:
                        continue
                genes_to_alleles[gene_name][haplotype].append(allele_hit)    
    return genes_to_alleles,novel_genes

def get_gene_names(gene_coords):
    gene_coords = load_bed_regions(gene_coords,True)
    gene_names = [genes[3] for genes in gene_coords]
    return gene_names

def get_output_assignment(assembly,ccs_reads,gene_names):
    output_lines = []
    ccs_read_thres = 5
    output_order = copy.deepcopy(gene_names)
    header = ["gene_name","haplotype_0","haplotype_1","haplotype_2","ccs_reads"]
    output_lines.append(header)
    for gene_name in assembly:
        if gene_name in gene_names:
            gene_names.remove(gene_name)
            out = [gene_name]
            for i in ["0","1","2"]:
                if i in assembly[gene_name]:
                    alleles = "."
                    if len(list(assembly[gene_name][i])) != 0:
                        alleles = ",".join(list(set(assembly[gene_name][i])))                    
                    out.append(alleles)
                else:
                    out.append(".")
            alleles_from_reads = "."
            if gene_name in ccs_reads:
                allele_hits_from_reads = Counter(ccs_reads[gene_name]["reads"])
                alleles_from_reads = []
                for allele in list(allele_hits_from_reads):
                    if allele_hits_from_reads[allele] > ccs_read_thres:
                        alleles_from_reads.append(allele)
                if len(alleles_from_reads) != 0:
                    alleles_from_reads = ",".join(alleles_from_reads)
                else:
                    alleles_from_reads = "."
            out.append(alleles_from_reads)
            output_lines.append(out)
    for gene_name in gene_names:
        alleles_from_reads = "."
        if gene_name in ccs_reads:
            alleles_from_reads = []
            allele_hits_from_reads = Counter(ccs_reads[gene_name]["reads"])
            for allele in list(allele_hits_from_reads):
                if allele_hits_from_reads[allele] > ccs_read_thres:
                    alleles_from_reads.append(allele)
            if len(alleles_from_reads) != 0:
                alleles_from_reads = ",".join(alleles_from_reads)        
            else:
                alleles_from_reads = "."
        out = [gene_name,"NA","NA","NA",alleles_from_reads]
        output_lines.append(out)
    reordered_output_lines = []
    for gene in output_order:
        for line in output_lines:
            if line[0] == gene:
                reordered_output_lines.append(line)
    return reordered_output_lines

def assign_alleles_to_genes(mapped_locus,gene_coords,gene_sequence,database,novel,assignment,mapped_ccs,genes_seq_from_ccs):
    # Save gene sequence from assembly and ccs reads into a fasta file
    extract_genes_from_assembly(mapped_locus,gene_coords,gene_sequence)
    extract_genes_from_ccs_reads(mapped_ccs,gene_coords,genes_seq_from_ccs)
    # Assign gene to allele
    assembly_genes_to_alleles,assembly_novel_genes = match_gene_to_allele_db(gene_sequence,database)
    ccs_genes_to_alleles,ccs_novel_genes = match_gene_to_allele_db(genes_seq_from_ccs,database)
    # Get gene names
    gene_names = get_gene_names(gene_coords)
    # Write allele assignment to file
    count_novel_genes_in_ccs_reads = Counter(["%s_%s" % (i,j) for i,j in ccs_novel_genes])
    with open(novel,'w') as fh:        
        header = ["gene_name","gene_seq","count_in_ccs_reads","origin"]
        fh.write("%s\n" % "\t".join(header))
        for gene_name,gene_seq in assembly_novel_genes:     
            count_in_ccs_reads = count_novel_genes_in_ccs_reads["%s_%s" % (gene_name,gene_seq)]
            fh.write("%s\t%s\t%s\tassembly\n" % (gene_name,gene_seq,count_in_ccs_reads))
        for gene_name_seq,count in count_novel_genes_in_ccs_reads.most_common():
            if count > 5:
                gene_name,gene_seq = gene_name_seq.split("_")
                fh.write("%s\t%s\t%s\tccs_reads\n" % (gene_name,gene_seq,count))
            else:
                break

    output_lines = get_output_assignment(assembly_genes_to_alleles,ccs_genes_to_alleles,gene_names)
    with open(assignment,'w') as fh:
        for line in output_lines:
            fh.write("%s\n" % "\t".join(line))
