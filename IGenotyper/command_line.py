#!/bin/env python
import os
from common import *
from lsf.lsf import Lsf

def get_bedgraph(mapped_reads,bedgraph,bam_filter):
    args = [mapped_reads,bam_filter,bedgraph]
    command = ("samtools view -F 3844 -Sbh %s %s | "
               "bedtools genomecov -bg -ibam stdin > %s" % tuple(args))
    if not non_emptyfile(bedgraph):
        os.system(command)

def get_ccs_reads(input_bam,dir_,threads):
    print "Running CCS..."
    command = ("ccs "
               "--numThreads %s "
               "--minPasses 2 "               
               "%s "
               "%s/ccs.bam " % (threads,input_bam,dir_))
    if not non_emptyfile("%s/ccs.bam.pbi" % dir_):
        os.system(command)

def turn_reads_to_fastq(dir_,reads_in_bam):
    print "Turnning reads to fastq..."
    args = [dir_,dir_,dir_,dir_,reads_in_bam,dir_]
    command = ("bam2fastq "
               "-o %s/ccs.reads.tmp %s/ccs.bam\n"
               "gunzip %s/ccs.reads.tmp.fastq.gz\n" 
               "bam2fastq "
               "-o %s/subreads.reads %s\n"
               "gunzip %s/subreads.reads.fastq.gz" % tuple(args))
    subreads_file = "%s/subreads.reads.fastq" % dir_ # Change back to subreads
    if not non_emptyfile(subreads_file):
        os.system(command)

def map_reads(threads,ref,reads,sorted_bamfile,type_of_read):
    print "Mapping %s..." % reads
    opts=""
    if type_of_read == "ccs":
        opts="--preset CCS"        
    args = [ref,reads,sorted_bamfile,threads,threads,opts]
    command = ("pbmm2 align "
               "%s "
               "%s "
               "%s "
               "--sort "
               "-j %s -J %s "
               "%s " % tuple(args))
    if not non_emptyfile("%s.bai" % sorted_bamfile):
        print command
        os.system(command)

def find_snp_candidates(ref,mapped_reads,snp_candidates):
    args = [ref,mapped_reads,snp_candidates]
    command = ("source activate whatshap-tool \n"
               "whatshap find_snv_candidates "
               "%s "
               "%s "
               "--pacbio "
               "-o %s \n" % tuple(args))
    if not non_emptyfile(snp_candidates):
        os.system(command)

def phase_snps(ref,snp_candidates,vcf,phased_vcf,ccs_mapped_reads,snp_candidates_filtered,regions_to_ignore):
    print "Calling and phasing SNPs..."
    args = [ref,ccs_mapped_reads,snp_candidates,
            snp_candidates,snp_candidates_filtered,regions_to_ignore,
            ref,vcf,snp_candidates_filtered,ccs_mapped_reads,
            ref,phased_vcf,vcf,ccs_mapped_reads]
    command = ("source activate whatshap-tool \n"
               "whatshap find_snv_candidates "
               "%s "
               "%s "
               "--pacbio "
               "-o %s \n"
               "conda deactivate \n"
               "IG-filter-vcf %s %s %s\n "
               "source activate whatshap-tool \n"
               "whatshap genotype "
               "--chromosome igh "
               "--sample sample "
               "--ignore-read-groups "
               "--reference %s "
               "-o %s "
               "%s "
               "%s \n"
               "whatshap phase "
               "--sample sample "
               "--reference %s "
               "--ignore-read-groups "
               "--distrust-genotypes "
               "-o %s "
               "%s "
               "%s \n"
               "conda deactivate" % tuple(args))
    if not non_emptyfile(phased_vcf):
        print command
        os.system(command)    

def get_phased_blocks(phased_vcf,hap_blocks):
    args = [hap_blocks,phased_vcf]
    command = ("source activate whatshap-tool \n"
               "whatshap stats "
               "--block-list %s "
               "%s \n"
               "conda deactivate whatshap-tool" % tuple(args))
    if not non_emptyfile(hap_blocks):
        os.system(command)    

def phase_reads(phased_vcf,mapped_reads,phased_mapped_reads,sample_name):
    print "Phasing reads in %s..." % mapped_reads
    args = [phased_vcf,mapped_reads,phased_mapped_reads,sample_name,phased_mapped_reads]
    command = ("IG-phase-reads %s %s %s %s \n"
               "samtools index %s" % tuple(args))
    if not non_emptyfile("%s.bai" % phased_mapped_reads):
        os.system(command)

def run_assembly_scripts(assembly_scripts,cluster,walltime,core,mem,queue):
    if not cluster:
        for script in assembly_scripts:
            command = "sh %s" % script
            os.system(command)
    else:
        hpc = Lsf()        
        for job in assembly_scripts:
            hpc.config(cpu=core,walltime=walltime,memory=mem,queue=queue)
            hpc.submit("%s" % job)
        hpc.wait()

def map_locus(threads,ref,tmp_dir,locus,sorted_bamfile):
    print "Mapping %s..." % locus
    threads = threads
    args = [ref,locus,sorted_bamfile,threads,threads]
    command = ("pbmm2 align "
               "%s "
               "%s "
               "%s "
               "--sort "
               "-j %s -J %s "
               "--preset CCS " % tuple(args))
    if not non_emptyfile("%s.bai" % sorted_bamfile):
        os.system(command)

def get_msa_coordinates(indel_dir,ref,package_python_directory):
    template_bash = "/sc/orga/work/rodrio10/software/in_github/IGenotyper/IGenotyper/get_msa_coords.sh"
    bashfile = "%s/indel_dir.sh" % indel_dir
    params = {
        'python_packages': package_python_directory,
        'ref': ref,
        'sv_calling_dir': indel_dir
        }
    write_to_bashfile(template_bash,bashfile,params)
    coords = "%s/msa_coords.bed" % indel_dir
    if not non_emptyfile(coords):
        os.system("sh %s" % bashfile)

def calls_svs_from_msa(indel_dir,package_python_directory):
    msa_coordsfn = "%s/msa_coords.bed" % indel_dir
    msa_coords = read_bedfile(msa_coordsfn)
    template_bash = "/sc/orga/work/rodrio10/software/in_github/IGenotyper/IGenotyper/sv_calling.sh"
    for chrom in msa_coords:
        for i,(start,end) in enumerate(msa_coords[chrom]):
            directory = "%s/%s/%s_%s" % (indel_dir,chrom,start,end)
            bashfile = "%s/sv_calling.sh" % directory
            params = {
                'dir': directory,
                'python_scripts': "/sc/orga/work/rodrio10/software/in_github/IGenotyper/IGenotyper",
                'chrom': chrom,
                'start': start,
                'end': end
                }
            write_to_bashfile(template_bash,bashfile,params)
            if not non_emptyfile("%s/svs.bed"):
                os.system("sh %s" % bashfile)

def call_svs_using_pbsv(ref,ref_index,ccs,aligned_bam,sample_name,sv_signature,sv_vcf,threads):
    args = [threads,threads,ref_index,ccs,aligned_bam,sample_name,
            aligned_bam,sv_signature,
            ref,sv_signature,sv_vcf,threads]
    command = ("#pbmm2 align -j %s -J %s %s %s %s --sort --preset CCS --sample %s \n "
               "pbsv discover %s %s \n "
               "pbsv call %s %s %s --ccs -j %s -m 10 " % tuple(args))
    print command
    if not non_emptyfile("%s" % sv_vcf):
        os.system(command)

def plot_bedgraph(bedgraphfn,tmp,bedgraph_plot):
    ini_file = "%s/out.ini" % tmp
    args = [bedgraphfn,ini_file,
            ini_file,bedgraph_plot]
    command = ("make_tracks_file --trackFiles %s -o %s \n "
               "pyGenomeTracks --tracks %s --region igh --outFileName %s " % tuple(args))
    if not non_emptyfile("%s" % bedgraph_plot):
        os.system(command)
    
