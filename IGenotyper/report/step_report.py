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
	
def get_assembly_size():
	pass
	
def get_number_of_contigs():
	pass
	
def get_RSS_SNV_count():
	pass
	
def get_leader_part_SNV_count():
	pass
 
def get_intron_SNV_count():
	pass

def get_deletion_count():
	pass
	
def get_insertion_count():
	pass
	
def get_alleles(gene_type):
	pass
	
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
		["IGH assembly size (bp)",get_assembly_size()],
		["IGH assembly number of contigs",get_number_of_contigs()],
		["# of SNVs in IGHJ region",get_SNV_count(self.phased_variants_vcf,"region","IGHJ")],
		["# of SNVs in IGHD region",get_SNV_count(self.phased_variants_vcf,"region","IGHD")],
		["# of SNVs in IGHV region",get_SNV_count(self.phased_variants_vcf,"region","IGHV")],
		["# of SNVs in RSS",snv_total_count - get_SNV_count(self.phased_variants_vcf,"RSS","None")],
		["# of SNVs in LP1",snv_total_count - get_SNV_count(self.phased_variants_vcf,"LP1","None")],
		["# of SNVs in introns",snv_total_count - get_SNV_count(self.phased_variants_vcf,"intronic","None")],
		["# of deletions (>3bps)",get_deletion_count()],
		["# of insertions (>3bps)",get_insertion_count()],
		["IGHJ alleles",get_alleles("IGHJ")],
		["IGHD alleles",get_alleles("IGHD")],
		["IGHV alleles",get_alleles("IGHV")],
		["Novel alleles",get_novel_alleles()],
		get_SV_genotypes()
		]
	]
