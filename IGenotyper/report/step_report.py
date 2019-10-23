#!/bin/env python

def get_coverage():
	pass
  
def get_region_coverage(region):
 	pass

def get_SNV_count(vcf_file,info_field,value,read_support=True):
	count = 0
	string_to_search = "%s=%s" % (info_field,value)
	with open(vcf_file,'r') as fh:
		for line in fh:
			line = line.rstrip().split('\t')
			if line[0] == "#":
				continue
			if string_to_search in line[XXX]:
				if read_support:
					if "read_support=Yes" in line[XXX]:
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
			variant = inline[XXXX]
			size = int(inline[XXXX])
			if size < size_thres:
				continue
			if variant == variant_type:
				count += 1
	return count
	
def get_alleles(genes_with_allele_assignment,gene_type):
	header = ["gene_name","haplotype_0","haplotype_1","haplotype_2","ccs_reads"]
	alleles = set()
	with open(genes_with_allele_assignment,'r') as fh:
		for line in fh:
			line = line.rstrip().split('\t')				
			if gene_type in line[0]:
				for column in range(1,5):
					if line[column] == ".":
						continue
					gene_allele = "*".join(line[1].split("="))
					alleles.add(gene_allele)
	return ",".join(list(alleles))
	
def get_novel_alleles():
	pass
	
def get_SV_genotypes():
	pass

def write_report(self):
	snv_total_count = get_SNV_count(self.phased_variants_vcf,"read_support","Yes")
	report = [
		["IGH Coverage (CCS)",get_coverage("ccs")],
		["IGH Coverage (subreads)",get_coverage("subreads")],
		["IGHJ Coverage",get_region_coverage("IGHJ")],
		["IGHD Coverage",get_region_coverage("IGHD")],
		["IGHV Coverage",get_region_coverage("IGHV")],
		["IGH assembly size (bp)",get_assembly_size(self.locus_fasta)],
		["IGH assembly number of contigs",get_number_of_contigs(self.locus_fasta)],
		["# of SNVs in IGHJ region",get_SNV_count(self.phased_variants_vcf,"region","IGHJ")],
		["# of SNVs in IGHD region",get_SNV_count(self.phased_variants_vcf,"region","IGHD")],
		["# of SNVs in IGHV region",get_SNV_count(self.phased_variants_vcf,"region","IGHV")],
		["# of SNVs in RSS",snv_total_count - get_SNV_count(self.phased_variants_vcf,"RSS","None")],
		["# of SNVs in LP1",snv_total_count - get_SNV_count(self.phased_variants_vcf,"LP1","None")],
		["# of SNVs in introns",snv_total_count - get_SNV_count(self.phased_variants_vcf,"intronic","None")],
		["# of deletions (>3bps)",get_indel_count(self.indels,"DEL")],
		["# of insertions (>3bps)",get_indel_count(self.indels,"INS")],
		["IGHJ alleles",get_alleles(self.genes_with_allele_assignment,"IGHJ")],
		["IGHD alleles",get_alleles(self.genes_with_allele_assignment,"IGHD")],
		["IGHV alleles",get_alleles(self.genes_with_allele_assignment,"IGHV")],
		["Novel alleles",get_novel_alleles()],
		get_SV_genotypes()
		]
	]
