#!/bin/env python
import textwrap
from Bio import SeqIO
from tabulate import tabulate


def get_coverage(read_type,output_file):
    coverage = None
    with open(output_file,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            if line[3] == read_type:
                if line[1] == "IGH coverage":
                    coverage = line[2]
    return round(float(coverage),2)
      
def get_region_coverage(region,output_file):
    coverage = None
    with open(output_file,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            if line[3] == region:
                coverage = sum(map(float,line[4:]))
    return round(coverage,2)

def get_SNV_count(vcf_file,info_field,value,read_support=True):
    count = 0
    string_to_search = "%s=%s" % (info_field,value)
    with open(vcf_file,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            if "#" in line[0]:
                continue
            if string_to_search in line[7]:
                if read_support:
                    if "read_support=Yes" in line[7]:
                        count += 1
                else:
                    count += 1
    return count

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

def get_indel_count(infile,variant_type):
    count = 0
    size_thres = 2
    with open(infile,'r') as infile_fh:
        for inline in infile_fh:
            inline = inline.strip().split('\t')
            variant = inline[3]
            size = int(inline[5])
            if size < size_thres:
                continue
            if variant == variant_type:
                count += 1
    return count

def get_alleles(genes_with_allele_assignment,gene_type):
    alleles = set()
    with open(genes_with_allele_assignment,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            if gene_type in line[0]:
                for column in range(1,5):
                    if line[column] == ".":
                        continue
                    for gene_alleles in line[column].split(","):
                        if "gene" not in gene_alleles:
                            continue
                        gene = gene_alleles.split("=")[1].split("_")[0]
                        allele = gene_alleles.split("=")[2]
                        gene_allele = "%s*%s" % (gene,allele)
                        alleles.add(gene_allele)
    return "\n".join(textwrap.wrap(",".join(list(alleles))))

def get_novel_alleles(novel_alleles_fn):
    header = ["gene_name","gene_seq","count_in_ccs_reads","origin"]
    novel_allele_seq = {}
    with open(novel_alleles_fn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            gene_seq = line[1]
            gene_name = line[0]
            ccs_read_support = line[2]
            origin = line[3]
            if gene_seq not in novel_allele_seq:
                novel_allele_seq[gene_seq] = {}
            novel_allele_seq[gene_seq]["gene_name"] = gene_name
            novel_allele_seq[gene_seq]["origin"] = []
            novel_allele_seq[gene_seq]["ccs_read_support"] = ccs_read_support
            novel_allele_seq[gene_seq]["origin"].append(origin)
            output = []
            for gene_seq in novel_allele_seq:
                gene_output = [
                    "Gene: %s" % novel_allele_seq[gene_seq]["gene_name"],
                    "Sequence: %s" % '\n'.join(textwrap.wrap(gene_seq)),
                    "# of CCS support: %s" % novel_allele_seq[gene_seq]["ccs_read_support"],
                    "Found in: %s" % ",".join(novel_allele_seq[gene_seq]["origin"]),
                    "-----"
                ]
            output.append("\n".join(gene_output))
        return "\n".join(output)
        
def get_SV_genotypes():
    pass

def write_report(self):
    phasing_stats_output = "%s/phasing_stats.txt" % self.tables_dir
    region_coverage_output = "%s/region_coverage.txt" % self.tables_dir
    indels = "%s.bed" % self.indels
    snv_total_count = get_SNV_count(self.assembly_snps,"read_support","Yes")    
    report = [
        ["IGH Coverage (CCS)", get_coverage("ccs",phasing_stats_output)],
        ["IGH Coverage (subreads)", get_coverage("subreads",phasing_stats_output)],
        ["IGHJ Coverage", get_region_coverage("IGHJ",region_coverage_output)],
        ["IGHD Coverage", get_region_coverage("IGHD",region_coverage_output)],
        ["IGHV Coverage", get_region_coverage("IGHV",region_coverage_output)],
        ["IGH assembly size (bp)",get_assembly_size(self.locus_fasta)],
        ["IGH assembly number of contigs",get_number_of_contigs(self.locus_fasta)],
        ["# of SNVs in IGHJ region",get_SNV_count(self.assembly_snps,"igh_region","IGHJ")],
        ["# of SNVs in IGHD region",get_SNV_count(self.assembly_snps,"igh_region","IGHD")],
        ["# of SNVs in IGHV region",get_SNV_count(self.assembly_snps,"igh_region","IGHV")],
        ["# of SNVs in RSS",snv_total_count - get_SNV_count(self.assembly_snps,"RSS","No")],
        ["# of SNVs in LP1",snv_total_count - get_SNV_count(self.assembly_snps,"LP1","No")],
        ["# of SNVs in introns",snv_total_count - get_SNV_count(self.assembly_snps,"intronic","No")],
        ["# of deletions (>3bps)",get_indel_count(indels,"DEL")],
        ["# of insertions (>3bps)",get_indel_count(indels,"INS")],
        ["IGHJ alleles",get_alleles(self.genes_with_allele_assignment,"IGHJ")],
        ["IGHD alleles",get_alleles(self.genes_with_allele_assignment,"IGHD")],
        ["IGHV alleles",get_alleles(self.genes_with_allele_assignment,"IGHV")],
        ["Novel alleles",get_novel_alleles(self.novel_alleles)]
    ]
    with open(self.report_file,'w') as fh:
        fh.write("%s\n" % tabulate(report,tablefmt="plain"))
