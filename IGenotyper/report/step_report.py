#!/bin/env python
from Bio import SeqIO
from string import Template

from ..common import *

def get_locus_coverage(read_type,output_file):
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

def get_SNV_count(vcf_file,info_fields,values,read_support=True):
    count = 0    
    strings_to_search = []
    for info_field,value in zip(info_fields,values):       
        strings_to_search.append("%s=%s" % (info_field,value))
    with open(vcf_file,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            if "#" in line[0]:
                continue
            strings_found = 0
            for string_to_search in strings_to_search:
                if string_to_search in line[7]:
                    strings_found += 1
            if len(strings_to_search) == strings_found:
                if read_support:
                    if "read_support=Yes" in line[7]:
                        count += 1
                else:
                    count += 1
    return count

def get_assembly_stats(output_file,stat):
    value = None
    with open(output_file,'r') as fh:
        for line in fh:
            line =  line.rstrip().split('\t')
            if line[0] == stat:
                value = line[1]
    return value

def get_assembly_size(output_file):
    size = get_assembly_stats(output_file,"Assembly size")
    return size

def get_number_of_contigs(output_file):
    num_contigs = get_assembly_stats(output_file,"Number of contigs")
    return num_contigs

def get_n50(output_file):
    n50 = get_assembly_stats(output_file,"N50")
    return n50

def get_locus_covered(output_file):
    percent_covered = get_assembly_stats(output_file,"% of IGH covered")
    return percent_covered

def get_indel_sv_count(infile,variant_type,region_type,min_size,max_size):
    count = 0
    with open(infile,'r') as infile_fh:
        for inline in infile_fh:
            inline = inline.strip().split('\t')
            variant = inline[3]
            size = int(inline[5])
            region = inline[20]
            if size < min_size:
                continue
            if size > max_size:
                continue
            if region_type != region:
                continue
            if variant != variant_type:
                continue
            count += 1
    return count

def get_indel_count(infile,variant_type,region_type):
    return get_indel_sv_count(infile,variant_type,region_type,2,50)

def get_sv_count(infile,variant_type,region_type):
    return get_indel_sv_count(infile,variant_type,region_type,49,1000000)
        
def get_SV_genotypes():
    pass
    
def get_coverage(phasing_stats_output,region_coverage_output):
    params = {
        'igh_coverage_ccs': get_locus_coverage("ccs",phasing_stats_output),
        'igh_coverage_subreads': get_locus_coverage("subreads",phasing_stats_output),
        'ighv_coverage_ccs': get_region_coverage("IGHV",region_coverage_output),
        'ighd_coverage_ccs': get_region_coverage("IGHD",region_coverage_output),
        'ighj_coverage_ccs': get_region_coverage("IGHJ",region_coverage_output),
    }
    return params

def get_assembly(assembly_stats_output):
    params = {
        'assembly_size': get_assembly_size(assembly_stats_output),
        'num_contigs': get_number_of_contigs(assembly_stats_output),
        'n50': get_n50(assembly_stats_output),
        'locus_coverage': get_locus_covered(assembly_stats_output)
    }
    return params

def get_variants_per_region(indels,assembly_snps,region):
    total_count = get_SNV_count(assembly_snps,["read_support","igh_region"],["Yes",region])
    segments = total_count - get_SNV_count(assembly_snps,["gene","igh_region"],["No",region])
    introns = total_count - get_SNV_count(assembly_snps,["intronic","igh_region"],["No",region])
    lp1 = total_count - get_SNV_count(assembly_snps,["LP1","igh_region"],["No",region])
    intergenic = total_count - segments - introns - lp1
    rss = total_count - get_SNV_count(assembly_snps,["RSS","igh_region"],["No",region])
    indel_dels = get_indel_count(indels,"DEL",region)
    indel_ins = get_indel_count(indels,"INS",region)
    sv_dels = get_sv_count(indels,"DEL",region)
    sv_ins = get_sv_count(indels,"INS",region)
    params = {
        '%s_snvs_total' % region: total_count,
        '%s_snvs_segment' % region: segments,
        '%s_snvs_introns' % region: introns,
        '%s_snvs_intergenic' % region: intergenic,
        '%s_snvs_lp1' % region: lp1,
        '%s_snvs_rss' % region: rss,
        '%s_indels_ins' % region: indel_ins,
        '%s_indels_del' % region: indel_dels,
        '%s_sv_ins' % region: sv_ins,
        '%s_sv_del' % region: sv_dels
    }
    return params

def get_variants(indels,assembly_snps):
    ighv = get_variants_per_region(indels,assembly_snps,"IGHV")
    ighd = get_variants_per_region(indels,assembly_snps,"IGHD")
    ighj = get_variants_per_region(indels,assembly_snps,"IGHJ")
    variants = {}
    variants.update(ighv)
    variants.update(ighd)
    variants.update(ighj)
    return variants

def get_alleles(allele_assignments):
    params = {}
    with open(allele_assignments,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            gene = line[3]
            allele = "%s,%s" % (line[4],line[5])
            params[gene.replace("-","_")] = allele
    return params

def get_novel_alleles_region(allele_assignments,region):
    params = {}
    novel_seqs = set()
    html_line = ""
    with open(allele_assignments,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            gene = line[3]
            if region not in gene:
                continue            
            if line[4] == "Novel":
                if line[9] in novel_seqs:
                    continue
                novel_seqs.add(line[9])
                html_line = "%s\n<tr><td>%s</td><td>%s</td></tr>" % (html_line,gene,line[9])
            if line[5] == "Novel":
                if line[10] in novel_seqs:
                    continue
                novel_seqs.add(line[10])
                html_line = "%s\n<tr><td>%s</td><td>%s</td></tr>" % (html_line,gene,line[10])
    params["%s_novel_alleles" % region] = html_line
    return params

def get_novel_alleles(allele_assignments):
    ighv = get_novel_alleles_region(allele_assignments,"IGHV")
    ighd = get_novel_alleles_region(allele_assignments,"IGHD")
    ighj = get_novel_alleles_region(allele_assignments,"IGHJ")
    novel_alleles = {}
    novel_alleles.update(ighv)
    novel_alleles.update(ighd)
    novel_alleles.update(ighj)
    return novel_alleles

def get_haplotype_coverage(igh_bed,igh_fai,ccs_phased_bam,tmp):
    step_size = 10
    window_size = 100 
    windows_bed = "%s/window.bed" % tmp
    output_coverage = "%s/coverage_per_window.bed" % tmp
    get_window_size(window_size,step_size,igh_bed,windows_bed)
    haplotype_coverage(ccs_phased_bam,windows_bed,output_coverage)

def get_phased_snvs(phased_ccs_snvs,tmp):
    snvs_per_pos = "%s/phased_snvs.txt" % tmp
    with open(snvs_per_pos,'w') as outfh:
        with open(phased_ccs_snvs,'r') as infh:
            for line in infh:
                line = line.rstrip().split('\t')
                if line[0][0] == "#":
                    continue
                if line[0] != "igh":
                    continue
                pos = line[1]
                gt = line[9].split(":")[0]
                outfh.write("%s\t%s\n" % (pos,gt))

def get_haplotype_blocks(regions_to_assemble,tmp):
    hap_blocks = "%s/hap_blocks.txt" % tmp
    with open(hap_blocks,'w') as outfh:
        with open(regions_to_assemble,'r') as infh:
            for line in infh:
                line = line.rstrip().split('\t')
                if line[3] != "1":
                    continue
                start = int(line[1])
                end = int(line[2])
                y = 1
                outfh.write("%s\t%s\t%s\tHaplotype\n" % (start,y,end))

def get_gene_coverage(ccs_phased_bam,gene_coords,tmp):
    output_coverage = "%s/gene_coverage.bed" % tmp
    haplotype_coverage(ccs_phased_bam,gene_coords,output_coverage)

def get_data_files(igh_bed,igh_fai,ccs_phased_bam,phased_ccs_snvs,regions_to_assemble,gene_coords,tmp):
    get_haplotype_coverage(igh_bed,igh_fai,ccs_phased_bam,tmp)
    get_phased_snvs(phased_ccs_snvs,tmp)
    get_haplotype_blocks(regions_to_assemble,tmp)
    #get_gene_coverage(ccs_phased_bam,gene_coords,tmp)


def create_plots(igh_bed,igh_fai,ccs_phased_bam,phased_ccs_snvs,regions_to_assemble,gene_coords,sv_coverage,allele_assignments,gene_coverage_output,plots,r_scripts,tmp):
    get_data_files(igh_bed,igh_fai,ccs_phased_bam,phased_ccs_snvs,regions_to_assemble,gene_coords,tmp)
    genes_coords_file = gene_coords
    coverage_per_window_file = "%s/coverage_per_window.bed" % tmp
    phased_snvs_file = "%s/phased_snvs.txt" % tmp
    hap_blocks_file = "%s/hap_blocks.txt" % tmp
    sv_coverage_file = sv_coverage
    gene_coverage_file = gene_coverage_output
    haplotype_alleles_file = allele_assignments
    locus_phasing_plot = "%s/locus_phasing.png" % plots
    sv_coverage_plot = "%s/sv_coverage.png" % plots
    gene_coverage_plot = "%s/gene_coverage.png" % plots
    haplotype_allele_plot = "%s/haplotype_alleles.png" % plots
    args = [r_scripts,genes_coords_file,coverage_per_window_file,phased_snvs_file,hap_blocks_file,
            sv_coverage_file,gene_coverage_file,haplotype_alleles_file,locus_phasing_plot,
            sv_coverage_plot,gene_coverage_plot,haplotype_allele_plot]
    command = ("Rscript %s/plot.R %s %s %s %s %s %s %s %s %s %s %s" % tuple(args))
    os.system(command)

def write_report(self):
    phasing_stats_output = "%s/phasing_stats.txt" % self.tables_dir
    region_coverage_output = "%s/region_coverage.txt" % self.tables_dir
    gene_coverage_output = "%s/gene_coverage.txt" % self.tables_dir
    assembly_stats_output = "%s/assembly_stats.txt" % self.tables_dir
    sv_coveragefn = "%s/sv_coverage.txt" % self.tables_dir
    indels = "%s.bed" % self.indels
    snvs = self.assembly_snps
    allele_assignments = self.genes_with_allele_assignment    
    create_plots(self.igh_coords,self.igh_fasta_fai,self.phased_ccs_mapped_reads,
                 self.phased_variants_vcf,self.regions_to_assemble,self.gene_coordinates,
                 sv_coveragefn,allele_assignments,gene_coverage_output,self.plots_dir,self.r_scripts,self.tmp_dir)
    all_params = [ get_coverage(phasing_stats_output,region_coverage_output),
                   get_assembly(assembly_stats_output),
                   get_variants(indels,snvs),
                   get_alleles(allele_assignments),
                   get_novel_alleles(allele_assignments)]
    params = {}
    for param in all_params:
        params.update(param)
    html_template = self.report_template
    filein = open(html_template)
    src = Template(filein.read())
    output_lines = src.safe_substitute(params)
    html = open(self.html_report,'w')
    html.write(output_lines)
    filein.close()
    html.close()
