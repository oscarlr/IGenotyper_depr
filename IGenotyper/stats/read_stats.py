#!/bin/env python
import os
import pysam
import pandas
import matplotlib
import pybedtools
from array import *
matplotlib.use('Agg')
from matplotlib import pyplot

from ..common import *
from ..command_line import *

def plot_length_distribution(lengths,type_,plotfn):
    values = pandas.Series(lengths)
    fig   = pyplot.figure()
    axes  = values.hist(color='gray', bins=100)
    fig   = pyplot.gcf()
    title = 'Histogram of %s lengths' % type_
    axes.set_title(title)
    axes.set_xlabel('Length')
    axes.set_ylabel('Count')
    axes.xaxis.grid(False)
    axes.axvline(values.mean(), color='k', linestyle='dashed', linewidth=1)
    _, max_ = pyplot.ylim()
    axes.text(values.mean() + values.mean()/10,
              max_ - max_/10,
              'Mean: {:.2f}'.format(values.mean()))
    axes.text(values.mean() + values.mean()/10,
              max_ - (max_/10)*2,
              'Std: {:.2f}'.format(values.std()))  
    axes.text(values.mean() + values.mean()/10,
              max_ - (max_/10)*3,
              '.5 q: {:.2f}'.format(values.quantile(.5)))  
    axes.text(values.mean() + values.mean()/10,
              max_ - (max_/10)*4,
              '.75 q: {:.2f}'.format(values.quantile(.75)))
    fig.savefig(plotfn, format='png')

def get_lengths(bam_file):
    lengths = []
    samfile = pysam.AlignmentFile(bam_file, "rb",check_sq=False)
    for read in samfile:
        lengths.append(read.query_length)
    return lengths

def plot_subread_length_distribution(subreads_bam,out_dir):
    plotfn = "%s/subreads_lengths.png" % out_dir
    if not non_emptyfile(plotfn):
        lengths = get_lengths(subreads_bam)
        plot_length_distribution(lengths,"subreads",plotfn)
    
def plot_ccs_length_distribution(ccs_bam,out_dir):
    plotfn = "%s/ccs_lengths.png" % out_dir
    #if not non_emptyfile(plotfn):
    lengths = get_lengths(ccs_bam)
    plot_length_distribution(lengths,"ccs",plotfn)

def generate_bedgraphs(mapped_phased_ccs,bedgraph_dir):
    for i in range(0,3):
        bedgraph = "%s/hap%s.bedGraph" % (bedgraph_dir,i)
        bam_filter = "-r %s" % i
        get_bedgraph(mapped_phased_ccs,bedgraph,bam_filter)

def phasing_counts(bam_file):
    hap_read_count = {}
    samfile = pysam.AlignmentFile(bam_file, "rb")
    references = samfile.references
    for ref in references:
        if ref not in hap_read_count:
            hap_read_count[ref] = {}
        for read in samfile.fetch(ref):
            if read.is_secondary:
                continue
            if read.is_unmapped:
                continue
            if read.is_supplementary:
                continue
            hap = read.get_tag("RG",True)[0]
            if hap not in hap_read_count[ref]:
                hap_read_count[ref][hap] = [0,0,0]
            hap_read_count[ref][hap][0] += 1
            hap_read_count[ref][hap][1] += read.query_length
            if read.query_qualities != None:
                hap_read_count[ref][hap][2] += sum(map(float,list(read.query_qualities)))/read.query_length
    return hap_read_count

def get_phasing_counts(mapped_phased_reads):
    number_reads = 0.0
    igh_hap_0_reads = 0.0
    igh_hap_1_reads = 0.0
    igh_hap_2_reads = 0.0
    igh_read_bases = 0.0
    igh_reads = 0.0
    target_rate = 0.0
    igh_coverage = 0.0
    igh_read_quality = 0.0
    igh_locus_length = 1193129.0
    read_counts = phasing_counts(mapped_phased_reads)
    for ref in read_counts:
        for hap in read_counts[ref]:
            number_reads += read_counts[ref][hap][0]
            if ref == "igh":
                igh_read_bases += read_counts[ref][hap][1]
                igh_read_quality += read_counts[ref][hap][2]
            if ref == "igh" and hap == "0":
                igh_hap_0_reads += read_counts[ref][hap][0]
            if ref == "igh" and hap == "1":
                igh_hap_1_reads += read_counts[ref][hap][0]
            if ref == "igh" and hap == "2":
                igh_hap_2_reads += read_counts[ref][hap][0]
    igh_reads = igh_hap_0_reads + igh_hap_1_reads + igh_hap_2_reads
    target_rate = igh_reads/number_reads
    igh_coverage = igh_read_bases/igh_locus_length
    average_read_length = igh_read_bases/igh_reads
    average_read_quality = igh_read_quality/igh_reads
    phasing_rate = (igh_hap_1_reads + igh_hap_2_reads)/igh_reads
    output = [["stat-0","Number of reads",number_reads],
              ["stat-0","IGH reads",igh_reads],
              ["stat-0","Percent of phased reads",phasing_rate],
              ["stat-0","IGH hap 0 reads",igh_hap_0_reads],
              ["stat-0","IGH hap 1 reads",igh_hap_1_reads],
              ["stat-0","IGH hap 2 reads",igh_hap_2_reads],
              ["stat-0","IGH on-target rate",target_rate],
              ["stat-0","IGH average read length",average_read_length],
              ["stat-0","IGH average read quality",average_read_quality],
              ["stat-0","IGH coverage",igh_coverage]]
    for ref in read_counts:
        for hap in read_counts[ref]:
            num_reads = read_counts[ref][hap][0]
            out = ["stat-1","chrom",ref,hap,num_reads]
            output.append(out)
    return output

def phasing_read_counts(mapped_phased_ccs,mapped_phased_subreads,tables_dir):
    phasing_stats_output = "%s/phasing_stats.txt" % tables_dir        
    if not non_emptyfile(phasing_stats_output):        
        ccs_stats = get_phasing_counts(mapped_phased_ccs)    
        subreads_stats = get_phasing_counts(mapped_phased_subreads)       
        with open(phasing_stats_output,'w') as fh:
            for line in ccs_stats:
                line.append("ccs")
                fh.write("%s\n" % "\t".join(map(str,line)))
            for line in subreads_stats:
                line.append("subreads")
                fh.write("%s\n" % "\t".join(map(str,line)))

def get_coverage(coords,bedgraph_dir,coveragefn,type_):
    if not non_emptyfile(coveragefn):
        total_coverage = {}
        lengths = {}
        for i in range(0,3):
            total_coverage[i] = {}
            bedgraph = "%s/hap%s.bedGraph" % (bedgraph_dir,i)
            hap_coverage = pybedtools.BedTool(bedgraph)
            entries = pybedtools.BedTool(coords)
            base_coverage = entries.intersect(hap_coverage,wao=True)
            for interval in base_coverage:
                name = interval[3]
                length = int(interval[2]) - int(interval[1])
                coverage = 0
                if int(interval[8]) != 0:
                    coverage = int(interval[7]) * int(interval[8])
                if name not in total_coverage[i]:
                    total_coverage[i][name] = 0
                if name not in lengths:
                    lengths[name] = length
                total_coverage[i][name] += coverage        
        entries = pybedtools.BedTool(coords)
        with open(coveragefn,'w') as fh:        
            header = ["chrom","start","end",type_,"0","1","2"]
            fh.write("%s\n" % "\t".join(header))
            for interval in entries:
                chrom = interval[0]
                start = interval[1]
                end = interval[2]
                entry = interval[3]
                out = [chrom,start,end,entry]
                for hap in range(0,3):
                    hap_coverage = total_coverage[hap][entry]/float(lengths[entry])
                    out.append(hap_coverage)
                fh.write("%s\n" % "\t".join(map(str,out)))

def get_genes_coverage(gene_coords,bedgraph_dir,tables_dir):
    coveragefn = "%s/gene_coverage.txt" % tables_dir
    get_coverage(gene_coords,bedgraph_dir,coveragefn,"gene")

def get_sv_coverage(sv_coords,bedgraph_dir,tables_dir):
    coveragefn = "%s/sv_coverage.txt" % tables_dir
    get_coverage(sv_coords,bedgraph_dir,coveragefn,"sv")

def get_region_coverage(region_coords,bedgraph_dir,tables_dir):
    coveragefn = "%s/region_coverage.txt" % tables_dir
    get_coverage(region_coords,bedgraph_dir,coveragefn,"sv")

def stats_on_reads(subreads,ccs,mapped_phased_ccs,mapped_phased_subreads,plots_dir,bedgraph_dir,tables_dir,gene_coords,sv_coords,region_coords):
    plot_subread_length_distribution(subreads,plots_dir)
    plot_ccs_length_distribution(ccs,plots_dir)
    generate_bedgraphs(mapped_phased_ccs,bedgraph_dir)
    phasing_read_counts(mapped_phased_ccs,mapped_phased_subreads,tables_dir)
    get_genes_coverage(gene_coords,bedgraph_dir,tables_dir)
    get_sv_coverage(sv_coords,bedgraph_dir,tables_dir)
    get_region_coverage(region_coords,bedgraph_dir,tables_dir)
    
