#!/bin/env python
import os
import sys
import argparse

from load import FileManager, CpuManager
from common import check_file_exist
from command_line import CommandLine
from phase.step_phase_reads import Phase
from assemble.step_assemble_reads import Assemble

def load_managers(args):
    file_manager = FileManager()
    cpu_manager = CpuManager()
    for manager in [file_manager,cpu_manager]:
        manager.load_args(args)
    return (file_manager,cpu_manager)
    
def run(args):
    file_manager,cpu_manager = load_managers(args)
    command_line_tools = CommandLine(file_manager,cpu_manager,args.sample_name)
    if args.phase:
        step = Phase(file_manager,cpu_manager,command_line_tools)
    elif args.assemble:
        step = Assemble(file_manager,cpu_manager,command_line_tools)
    else:
        sys.exit("Select phase, assemble, detect, stats or report")
    step.load_args(args)
    return step()
        
def main():
    parser = argparse.ArgumentParser(description='Process IGH capture data')
    parser.add_argument('input_bam', metavar='input_bam', type=check_file_exist,
                        help='pacbio subreads bam file')
    parser.add_argument('outdir', metavar='outdir', 
                        help='output directory')
    parser.add_argument('--phase', action='store_true', default=False,
                        help='Map and phase reads')
    parser.add_argument('--assemble', action='store_true', default=False,
                        help='Assemble reads')
    parser.add_argument('--extend_assembly', action='store_true', default=False,
                        help='Extend assemblies. Still testing. Dont use.')
    parser.add_argument('--detect', action='store_true', default=False,
                        help='Detect variants')
    parser.add_argument('--stats', action='store_true', default=False,
                        help='Generate stats')
    parser.add_argument('--report', action='store_true', default=False,
                        help='Generate report')
    parser.add_argument('--tmp_dir', metavar='',
                        help='temporary directory')
    parser.add_argument('--threads', metavar='', default=1,
                        help='Number of threads for everything')
    parser.add_argument('--cluster', action='store_true', default=False,
                        help='Use cluster')
    parser.add_argument('--cluster_queue', metavar="", default="express",
                        help='Queue for cluster')
    # parser.add_argument('--cluster_threads', metavar='', default=4,
    #                     help='Number of threads for cluster jobs')
    parser.add_argument('--cluster_walltime', metavar="", default=2,
                        help='Walltime for cluster')
    parser.add_argument('--cluster_mem', metavar="", default=5,
                        help='memory for cluster')
    # parser.add_argument('--haploid',action='store_true', default=False,
    #                     help='Run Quiver in haploid mode')
    parser.add_argument('--phased_vcf_file', metavar="", type=check_file_exist,
                        help='Run IG with this phased VCF file')
    parser.add_argument('--pacbio_data_type', metavar="", default="RS",
                        help='Pacbio data type (RS or Sequel), either "RS" or "SQ"')
    # parser.add_argument('--phased_vcf_file_sample_name',metavar="",default="sample",
    #                     help='Sample name in phased VCF file')
    parser.add_argument('--sample_name',metavar="",default="sample",
                        help='Sample name')
    # parser.add_argument('--add_unphased_reads',action='store_true', default=False,
    #                     help='Add unphased reads to phased region')
    parser.add_argument('--split',action='store_true', default=False,
                        help='Split assembly regions into SV/non-SV regions')
    parser.add_argument('--secondary_read_score',metavar="",default=500,type=int,
                        help='Min secondary read score to move')
    parser.add_argument('--keep',action='store_true',default=False,
                        help='Keep intermediate files')
    args = parser.parse_args()
    return run(args)
