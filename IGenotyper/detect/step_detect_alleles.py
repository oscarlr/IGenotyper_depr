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
from step_detect_snps import labeled_hap_blocks,snp_in_gene_feature

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
            if read.is_secondary:
                continue
            if read.is_supplementary:
                continue            
            if read.is_unmapped:
                continue                
            mapping_cord = [chrom,read.reference_start,read.reference_end]
            if assembly:
                if not read_overlap_region(read,mapping_cord):
                    continue
            gene_coord = [chrom,start,end]
            # if assembly:
            #     if not read_overlap_region(read,gene_coord):
            #         continue
            gene_sequence = extract_sequence_from(read,chrom,start,end)
            if gene_sequence == "":
                continue
            if len(gene_sequence) < (end - start - 10):
                continue
            if len(gene_sequence) > ((end - start) + 10):
                continue
            if assembly:
                haplotype = get_haplotype(read.query_name)
            else:
                haplotype = read.get_tag("RG",True)[0]
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
    novel_genes = {}
    genes_to_alleles = {}
    potentially_novel_genes = []
    for gene_extraction_name, gene_seq in matches:
        gene_name = gene_extraction_name.split(".")[0].split("=")[1]
        haplotype = gene_extraction_name.split(".")[1].split("=")[1]
        if len(matches[(gene_extraction_name, gene_seq)]) == 0:
            if gene_name not in novel_genes:
                novel_genes[gene_name] = {}
            if haplotype not in novel_genes[gene_name]:
                novel_genes[gene_name][haplotype] = []
            novel_genes[gene_name][haplotype].append(gene_seq)
        else:
            for allele_hit in matches[(gene_extraction_name, gene_seq)]:
                if allele_hit.split("_")[0].split("=")[1] != gene_name:
                    potentially_novel_genes.append(gene_extraction_name)
                    continue
                if gene_name not in genes_to_alleles:
                    genes_to_alleles[gene_name] = {}
                if haplotype not in genes_to_alleles[gene_name]:
                    genes_to_alleles[gene_name][haplotype] = []
                genes_to_alleles[gene_name][haplotype].append(allele_hit)    
    # Could have hits but to the wrong gene
    for potentially_novel_gene in potentially_novel_genes:
        for gene_extraction_name, gene_seq in matches:
            if potentially_novel_gene != gene_extraction_name:
                continue
            gene_name = gene_extraction_name.split(".")[0].split("=")[1]
            haplotype = gene_extraction_name.split(".")[1].split("=")[1]            
            if gene_name in genes_to_alleles:
                if haplotype in genes_to_alleles[gene_name]:
                    continue
            if gene_name not in novel_genes:
                novel_genes[gene_name] = {}
            if haplotype not in novel_genes[gene_name]:
                novel_genes[gene_name][haplotype] = []
            novel_genes[gene_name][haplotype].append(gene_seq)
    return genes_to_alleles,novel_genes

def get_gene_names(gene_coords):
    gene_coords = load_bed_regions(gene_coords,True)
    gene_names = [genes[3] for genes in gene_coords]
    return gene_names

def get_genes_per_sv():
    # "IGHV4-28 to IGHV3-33 complex event": ["IGHV4-28","IGHV3-30","IGHV4-30-2","IGHV3-30-3",
    #                                        "IGHV4-30-4","IGHV3-30-5","IGHV4-31","IGHV3-33"],
    # "IGHV4-38-2 to IGHV1-38-4 insertion": ["IGHV4-38-2","IGHV3-43D","IGHV3-38-3","IGHV1-38-4"],
    genes_per_sv = {"IGHV7-4-1 insertion": ["IGHV7-4-1"],
                    "IGHV5-10-1/IGHV3-64D haplotype": ["IGHV5-10-1","IGHV3-64D"],
                    "IGHV3-23D duplication": ["IGHV3-23D","IGHV3-23"],
                    "IGHV4-38-2 to IGHV1-38-4 insertion": ["IGHV4-38-2","IGHV3-43D","IGHV3-38-3","IGHV1-38-4"],
                    "IGHV1-69D to IGHV2-70D insertion": ["IGHV1-69","IGHV2-70D","IGHV1-69-2","IGHV1-69D","IGHV2-70"],
                    "IGHV1-8/IGHV3-9 haplotype": ["IGHV1-8","IGHV3-9"]}    
    return genes_per_sv

def get_svs_per_gene():
    sv_per_gene = {
        "IGHV7-4-1": "IGHV7-4-1.CH17",
        "IGHV5-10-1": "IGHV5-10-1.CH17",
        "IGHV3-64D": "IGHV5-10-1.CH17",
        "IGHV3-23D": "IGHV3-23.region.ABC9",
        "IGHV4-38-2": "IGHV4-38-2.region.mixFosmids",
        "IGHV3-43D": "IGHV4-38-2.region.mixFosmids",
        "IGHV3-38-3": "IGHV4-38-2.region.mixFosmids",
        "IGHV1-38-4": "IGHV4-38-2.region.mixFosmids",
        "IGHV2-70D": "IGHV1-69.region.CH17",
        "IGHV1-69-2": "IGHV1-69.region.CH17",
        "IGHV1-69D": "IGHV1-69.region.CH17",
        "IGHV1-8": "IGHV1-8.GRCh37",
        "IGHV3-9": "IGHV1-8.GRCh37"
    }
    return sv_per_gene

def get_genes_in_sv():
    genes = []
    genes_per_sv = get_genes_per_sv()
    for sv in genes_per_sv:
        for gene in genes_per_sv[sv]:
            genes.append(gene)
    return genes

def ccs_reads_per_allele(ccs_reads,gene,allele):
    count = 0
    if gene in ccs_reads:
        alleles = []
        for hap in ccs_reads[gene]:
            for read_allele in ccs_reads[gene][hap]:
                alleles.append(read_allele)
        allele_counts = Counter(alleles)
        for allele_count in allele_counts:
            allele_num = allele_count.split("=")[-1]
            if allele_num == allele:
                count = allele_counts[allele_count]
    return count

def get_haplotype_allele_from_assembly(assembly_allele,haplotype,sv_gene):
    allele = None
    if "0" in assembly_allele:
        assert 1 not in assembly_allele
        assert 2 not in assembly_allele
        allele = assembly_allele["0"][0].split('=')[-1]
        if haplotype == "2":
            if sv_gene:
                allele = "Deleted"
    else:
        if haplotype in assembly_allele:
            allele = assembly_allele[haplotype][0].split('=')[-1]
    return allele

def get_haplotype_allele_from_ccs_reads(ccs_reads_allele,gene,haplotype):
    allele = None
    thres = 5
    if gene in ccs_reads_allele:
        if haplotype in ccs_reads_allele:
            allele_counts = Counter(ccs_reads_allele[haplotype])
            for allele_count in allele_counts:
                if allele_counts[allele_count] > thres:
                    allele = allele_count.split('=')[-1]
        if allele == None:
            if haplotype != "0":
                haplotype = "0"
            if haplotype in ccs_reads_allele:
                allele_counts = Counter(ccs_reads_allele[haplotype])
                for allele_count in allele_counts:
                    if allele_counts[allele_count] > thres:
                        allele = allele_count.split('=')[-1]            
    return allele

def get_novel_seq(assembly_novel_genes,gene_name,hap):
    seq = None
    if gene_name in assembly_novel_genes:
        if hap in assembly_novel_genes[gene_name]:
            seq = str(assembly_novel_genes[gene_name][hap][0])
        if seq == None:
            hap = "0"
            if hap in assembly_novel_genes[gene_name]:
                seq = str(assembly_novel_genes[gene_name][hap][0])
    return seq

def count_novel_ccs_reads(ccs_novel_genes,gene,assembly_seq):
    count = 0
    if gene not in ccs_novel_genes:
        return count
    for hap in ccs_novel_genes[gene]:
        for ccs_seq in ccs_novel_genes[gene][hap]:
            if ccs_seq == assembly_seq:
                count += 1
    return count


def check_novel_or_deleted(novel_alleles,gene,hap,ccs_reads,sv_gene):
    # Is this allele novel?
    haplotype_novel_seq = get_novel_seq(novel_alleles,gene,hap)
    if haplotype_novel_seq == None:
        # Is this is allele present in CCS reads?
        haplotype_allele = get_haplotype_allele_from_ccs_reads(ccs_reads,gene,hap)
        if haplotype_allele == None:
            haplotype_allele = "ND"
            if sv_gene:
                haplotype_allele = "Deleted"
    else:
        haplotype_allele = "Novel"
    return haplotype_allele

def coverage_in_hap_gene(hap,gene,gene_coverage):
    thres = 10
    coverage = 0
    with open(gene_coverage,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            ingene = line[3]
            if gene != ingene:
                continue
            if hap == "1":
                coverage = float(line[5])
            if hap == "2":
                coverage = float(line[6])
    if coverage > thres:
        return True
    return False

def get_gene_coords(gene_coords,gene):
    for chrom,start,end,ingene in gene_coords:
        if gene  == ingene:
            return (chrom,start,end)
    sys.exit("SHOULD NOT HAPPEN")

def reassemble_gene(hap,gene,tmp,ccs_to_ref,raw_subreads_to_ref,python_scripts,subreads,assembly_script,cluster,cluster_walltime,cluster_threads,cluster_mem,cluster_queue,gene_coords):
    chrom, start, end = get_gene_coords(gene_coords,gene)
    threads = cluster_threads
    size = 2000
    output = "%s/%s/%s" % (tmp,gene,hap)
    ccs_to_ref = ccs_to_ref 
    raw_subreads_to_ref = raw_subreads_to_ref
    python_scripts = python_scripts
    subreads = subreads
    create_directory(output)
    template_bash = assembly_script
    bashfile = "%s/assemble.sh" % output
    params = {
        "chrom": chrom,
        "start": start,
        "end": end,
        "hap": "-r %s" % hap,
        "threads": threads,
        "size": size,
        "output": output,
        "ccs_to_ref": ccs_to_ref,
        "raw_subreads_to_ref": raw_subreads_to_ref,
        "python_scripts": python_scripts,
        "subreads": subreads,
    }
    contig_fasta = "%s/contig.fasta" % output
    if not os.path.isfile(contig_fasta):
        write_to_bashfile(template_bash,bashfile,params)
        run_assembly_scripts([bashfile],cluster,cluster_walltime,cluster_threads,cluster_mem,cluster_queue)
    if os.path.isfile(contig_fasta):
        return contig_fasta
    return None

def regenotype_gene(hap,gene,gene_coverage,tmp,ccs_to_ref,raw_subreads_to_ref,python_scripts,subreads,assembly_script,cluster,cluster_walltime,cluster_threads,cluster_mem,cluster_queue,database,gene_coords):
    if coverage_in_hap_gene(hap,gene,gene_coverage):
        sequencefn = reassemble_gene(hap,gene,tmp,ccs_to_ref,raw_subreads_to_ref,python_scripts,subreads,assembly_script,cluster,cluster_walltime,cluster_threads,cluster_mem,cluster_queue,gene_coords)
        if sequencefn != None:
            query_entries = SeqIO.to_dict(SeqIO.parse(sequencefn,"fasta"))
            database_entries = SeqIO.to_dict(SeqIO.parse(database,"fasta"))
            matches = get_matches(query_entries,database_entries)
            for gene_extraction_name, gene_seq in matches:
                for allele_hit in matches[(gene_extraction_name, gene_seq)]:
                    if allele_hit.split("_")[0].split("=")[1] != gene:
                        continue
                    allele = allele_hit.split("_")[1].split("=")[1]
                    return allele
    return "ND"

def get_alleles_per_hap(assembly,novel_alleles,ccs_reads,genes,gene_coverage,tmp,ccs_to_ref,raw_subreads_to_ref,python_scripts,subreads,assembly_script,cluster,cluster_walltime,cluster_threads,cluster_mem,cluster_queue,database,gene_coords):
    sv_genes = get_genes_in_sv()
    genes_with_alleles = {}
    for gene in genes:
        sv_gene = False
        if gene in sv_genes:
            sv_gene = True
        if gene in assembly:
            haplotype_1_allele = get_haplotype_allele_from_assembly(assembly[gene],"1",sv_gene)
            haplotype_2_allele = get_haplotype_allele_from_assembly(assembly[gene],"2",sv_gene)
            if haplotype_1_allele == None:
                haplotype_1_allele = check_novel_or_deleted(novel_alleles,gene,"1",ccs_reads,sv_gene)
            elif haplotype_2_allele == None:
                haplotype_2_allele = check_novel_or_deleted(novel_alleles,gene,"2",ccs_reads,sv_gene)
        else:
            haplotype_1_allele = get_haplotype_allele_from_ccs_reads(ccs_reads,gene,"1")
            haplotype_2_allele = get_haplotype_allele_from_ccs_reads(ccs_reads,gene,"2")
            if haplotype_1_allele == None:
                haplotype_1_allele = check_novel_or_deleted(novel_alleles,gene,"1",ccs_reads,sv_gene)
            if haplotype_2_allele == None:
                haplotype_2_allele = check_novel_or_deleted(novel_alleles,gene,"2",ccs_reads,sv_gene)
        assert haplotype_1_allele != None
        assert haplotype_2_allele != None
        if haplotype_1_allele == "ND" and gene not in ["IGHV3-23","IGHV3-23D"]:
            haplotype_1_allele = regenotype_gene("1",gene,gene_coverage,tmp,ccs_to_ref,raw_subreads_to_ref,python_scripts,subreads,assembly_script,cluster,cluster_walltime,cluster_threads,cluster_mem,cluster_queue,database,gene_coords)
        if haplotype_2_allele == "ND" and gene not in ["IGHV3-23","IGHV3-23D"]:
            haplotype_2_allele = regenotype_gene("2",gene,gene_coverage,tmp,ccs_to_ref,raw_subreads_to_ref,python_scripts,subreads,assembly_script,cluster,cluster_walltime,cluster_threads,cluster_mem,cluster_queue,database,gene_coords)
        genes_with_alleles[gene] = [haplotype_1_allele,haplotype_2_allele]
    return genes_with_alleles

def get_deleted_allele(alleles):
    hap1_allele = alleles[0]
    hap2_allele = alleles[1]
    deleted_allele = None
    if hap1_allele == "Deleted":    
        deleted_allele = "1"
    if hap2_allele == "Deleted":    
        deleted_allele = "2"
    return deleted_allele

def fix_1_69_alleles(genes_with_alleles,hap1_allele,hap2_allele):
    v2_70d = get_deleted_allele(genes_with_alleles["IGHV2-70D"])
    v1_69_2 = get_deleted_allele(genes_with_alleles["IGHV1-69-2"])
    v1_69d = get_deleted_allele(genes_with_alleles["IGHV1-69D"])
    if [v2_70d,v1_69_2,v1_69d].count("1") > 1:
        hap1_allele = "Deleted"
    if [v2_70d,v1_69_2,v1_69d].count("2") > 1:
        hap2_allele = "Deleted"                        
    return hap1_allele,hap2_allele
    

def print_gene_alleles(genes_with_alleles,assembly_novel_genes,gene_coords,ccs_genes_to_alleles,ccs_novel_genes,haplotype_blocks,svs_genotyped,outputfh):
    header = ["chrom","start","end","gene_name",
              "haplotype_1_allele","haplotype_2_allele",
              "num_of_ccs_reads_support_for_haplotype_1","num_of_ccs_reads_support_for_haplotype_2",
              "haplotype_block",
              "haplotype_1_novel_sequence","haplotype_2_novel_sequence"]
    phased_blocks = labeled_hap_blocks(haplotype_blocks)
    output = [header]
    svs_per_gene = get_svs_per_gene()
    sv_genotypes = read_genotype_svs(svs_genotyped)
    for chrom,start,end,gene in gene_coords:
        hap1_gene_ccs_support = ccs_reads_per_allele(ccs_genes_to_alleles,gene,genes_with_alleles[gene][0])
        hap2_gene_ccs_support = ccs_reads_per_allele(ccs_genes_to_alleles,gene,genes_with_alleles[gene][1])
        haplotype_1_novel_sequence = get_novel_seq(assembly_novel_genes,gene,"1")
        haplotype_2_novel_sequence = get_novel_seq(assembly_novel_genes,gene,"2")
        if haplotype_1_novel_sequence != None:
            hap1_gene_ccs_support = count_novel_ccs_reads(ccs_novel_genes,gene,haplotype_1_novel_sequence)
        if haplotype_2_novel_sequence != None:
            hap2_gene_ccs_support = count_novel_ccs_reads(ccs_novel_genes,gene,haplotype_2_novel_sequence)
        haplotype_block = snp_in_gene_feature(start,phased_blocks)
        hap1_allele = genes_with_alleles[gene][0]
        hap2_allele = genes_with_alleles[gene][1]        
        if gene in svs_per_gene:
            sv = svs_per_gene[gene]
            if sv_genotypes[sv] == "0/0":
                hap1_allele = "Deleted"
                hap2_allele = "Deleted"
            if sv == "IGHV1-69.region.CH17" and sv_genotypes[sv] == "0/1":
                if "Deleted" not in [hap1_allele,hap2_allele]:
                    hap1_allele, hap2_allele = fix_1_69_alleles(genes_with_alleles,hap1_allele,hap2_allele)
        # 3-23 will always come first
        if gene == "IGHV3-23":
            ighv_3_23_hap1_allele = genes_with_alleles["IGHV3-23"][0]
            ighv_3_23_hap2_allele = genes_with_alleles["IGHV3-23"][1]
            ighv_3_23d_hap1_allele = genes_with_alleles["IGHV3-23D"][0]
            ighv_3_23d_hap2_allele = genes_with_alleles["IGHV3-23D"][1]
            if ighv_3_23_hap1_allele == "Deleted" and ighv_3_23d_hap1_allele != "Deleted":
                hap1_allele = ighv_3_23d_hap1_allele
                genes_with_alleles["IGHV3-23D"][0] = "Deleted"
            if ighv_3_23_hap2_allele == "Deleted" and ighv_3_23d_hap2_allele != "Deleted":
                hap2_allele = ighv_3_23d_hap2_allele
                genes_with_alleles["IGHV3-23D"][1] = "Deleted"
        out = [chrom,start,end,gene,
               hap1_allele,hap2_allele,
               hap1_gene_ccs_support,hap2_gene_ccs_support,
               haplotype_block,
               haplotype_1_novel_sequence,haplotype_2_novel_sequence]
        output.append(out)
    with open(outputfh,'w') as fh:
        for line in output:
            fh.write("%s\n" % "\t".join(map(str,line)))               

def assign_alleles_to_genes(mapped_locus,gene_coords,gene_sequence,database,novel,assignment,mapped_ccs,genes_seq_from_ccs,haplotype_blocks,pacbio_data_type,sv_genotypes,
                            stats_dir,tmp,raw_subreads_to_ref,python_scripts,subreads,assembly_script,
                            cluster,cluster_walltime,cluster_threads,cluster_mem,cluster_queue):
    # Save gene sequence from assembly and ccs reads into a fasta file
    extract_genes_from_assembly(mapped_locus,gene_coords,gene_sequence)
    extract_genes_from_ccs_reads(mapped_ccs,gene_coords,genes_seq_from_ccs)
    # Assign gene to allele
    assembly_genes_to_alleles,assembly_novel_genes = match_gene_to_allele_db(gene_sequence,database)
    ccs_genes_to_alleles,ccs_novel_genes = match_gene_to_allele_db(genes_seq_from_ccs,database)
    # Get gene names
    gene_names = get_gene_names(gene_coords)
    gene_coords = load_bed_regions(gene_coords,True)
    gene_coverage = "%s/tables/gene_coverage.txt" % stats_dir
    #genes_with_alleles = get_alleles_per_hap(assembly_genes_to_alleles,assembly_novel_genes,ccs_genes_to_alleles,gene_names)
    genes_with_alleles = get_alleles_per_hap(assembly_genes_to_alleles,assembly_novel_genes,ccs_genes_to_alleles,gene_names,
                                             gene_coverage,tmp,mapped_ccs,raw_subreads_to_ref,python_scripts,subreads,assembly_script,
                                             cluster,cluster_walltime,cluster_threads,cluster_mem,cluster_queue,database,gene_coords)
    print_gene_alleles(genes_with_alleles,assembly_novel_genes,gene_coords,ccs_genes_to_alleles,ccs_novel_genes,haplotype_blocks,sv_genotypes,assignment)
