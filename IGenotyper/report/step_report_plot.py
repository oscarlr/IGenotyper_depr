#!/bin/env python
import os


###

def non_emptyfile(checkfile):
    return os.path.isfile(checkfile) and os.path.getsize(checkfile) > 0

def get_window_size(window_size,step,inbed,outbed):
    args = [inbed,window_size,step,outbed]
    command = ("bedtools makewindows -b %s -w %s -s %s > %s \n" % tuple(args))
    if not non_emptyfile("%s" % outbed):
        os.system(command)

def haplotype_coverage(bam,windows_bed,output_coverage):
    args = [bam,windows_bed,output_coverage,
            bam,windows_bed,output_coverage,
            bam,windows_bed,output_coverage]
    command = ("samtools view -Sbh %s -r 0 | bedtools coverage -counts -a %s -b stdin | awk '{ print $0\"\t0\"}' > %s \n"
               "samtools view -Sbh %s -r 1 | bedtools coverage -counts -a %s -b stdin | awk '{ print $0\"\t1\"}' >> %s \n"
               "samtools view -Sbh %s -r 2 | bedtools coverage -counts -a %s -b stdin | awk '{ print $0\"\t2\"}' >> %s \n" % tuple(args))
    if not non_emptyfile("%s" % output_coverage):
        os.system(command)
###

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
                print line
                start = int(line[1])
                end = int(line[2])
                y = 1
                outfh.write("%s\t%s\t%s\tHaplotype\n" % (start,y,end))

def get_gene_coverage(ccs_phased_bam,gene_coords,tmp):
    output_coverage = "%s/gene_coverage.bed" % tmp
    haplotype_coverage(ccs_phased_bam,gene_coords,output_coverage)

def get_data_files():
    igh_bed = "/sc/orga/work/rodrio10/software/in_github/IGenotyper_clean_version/IGenotyper/data/igh.bed"
    igh_fai = "/sc/orga/work/rodrio10/software/in_github/IGenotyper_clean_version/IGenotyper/data/igh.fasta.fai"
    ccs_phased_bam = "/sc/hydra/work/rodrio10/projects/Human_IGH_1KG_fosmids/data/project_scratch/2019-10-03_redo_flu_samples/flu_4_4_10_run_merge/alignments/ccs_to_ref_phased.sorted.bam"
    phased_ccs_snvs = "/sc/hydra/work/rodrio10/projects/Human_IGH_1KG_fosmids/data/project_scratch/2019-10-03_redo_flu_samples/flu_4_4_10_run_merge/variants/from_reads/ccs_phased_variants.vcf"
    regions_to_assemble = "/sc/hydra/work/rodrio10/projects/Human_IGH_1KG_fosmids/data/project_scratch/2019-10-03_redo_flu_samples/flu_4_4_10_run_merge/assembly/ccs_variants_phased.bed"
    gene_coords = "/sc/orga/work/rodrio10/software/in_github/IGenotyper_clean_version/IGenotyper/data/gene_coords.bed"
    tmp = "tmp"
    get_haplotype_coverage(igh_bed,igh_fai,ccs_phased_bam,tmp)
    get_phased_snvs(phased_ccs_snvs,tmp)
    get_haplotype_blocks(regions_to_assemble,tmp)
    get_gene_coverage(ccs_phased_bam,gene_coords,tmp)

get_data_files()
