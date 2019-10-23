#!/bin/env python

def get_coverage():
	pass
  
def get_region_coverage(region):
 	pass

def get_region_SNV_count(region):
	pass
	
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
	report = [
		["IGH Coverage (CCS)",get_coverage("ccs")],
		["IGH Coverage (subreads)",get_coverage("subreads")],
		["IGHJ Coverage",get_region_coverage("IGHJ")],
		["IGHD Coverage",get_region_coverage("IGHD")],
		["IGHV Coverage",get_region_coverage("IGHV")],
		["IGH assembly size (bp)",get_assembly_size()],
		["IGH assembly number of contigs",get_number_of_contigs()],
		["# of SNVs in IGHJ region",get_SNV_count("IGHJ")],
		["# of SNVs in IGHD region",get_SNV_count("IGHD")],
		["# of SNVs in IGHV region",get_SNV_count("IGHV")],
		["# of SNVs in RSS",get_SNV_count("IGHV")],
		["# of SNVs in LP1",get_SNV_count("IGHV")],
		["# of SNVs in introns",get_SNV_count("IGHV")]
	]
