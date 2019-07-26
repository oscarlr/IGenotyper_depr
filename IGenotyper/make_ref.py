#!/bin/env python
import os
import pysam
import argparse

from common import *

def create_reference(hg19_reffn,igh_sequence,igh_specific_reference):
    print "Creating IGH specific reference ..."
    hg19 = pysam.FastaFile(hg19_reffn)
    igh = pysam.FastaFile(igh_sequence)
    
    igh_sequence = igh.fetch("igh")

    with open(igh_specific_reference,'w') as fh:
        for chromosome in hg19.references:
            if chromosome == "chr14":
                fh.write(">igh\n")
                fh.write("%s\n" % igh_sequence)
                fh.write(">chr14\n")
                fh.write("%s\n" % hg19.fetch("chr14",1,106326710))
            else:
                fh.write(">%s\n" % chromosome)
                fh.write("%s\n" % hg19.fetch(chromosome))

def create_blasr_index(blasr_index_directory,igh_specific_reference):
    blasr_index_ref = "%s/reference.fasta" % blasr_index_directory
    os.symlink(igh_specific_reference, blasr_index_ref)

def create_pbmm2_index(pbmm2_index_directory,igh_specific_reference):
    print "Indexing IGH specific reference ..."
    pbmm2_index_ref = "%s/reference.fasta" % pbmm2_index_directory
    os.symlink(igh_specific_reference, pbmm2_index_ref)
    pbmm2_index = "%s/reference.mmi" % pbmm2_index_directory
    command = ("pbmm2 index %s %s" % (pbmm2_index_ref,pbmm2_index))
    if not non_emptyfile(pbmm2_index):
        os.system(command)

def main():
    parser = argparse.ArgumentParser(description='Create IGH specific reference')
    parser.add_argument('hg19', help='Path to hg19 reference')
    parser.add_argument('--directory', help='Write to this directory instead of package directory')
    args = parser.parse_args()

    output_directory = args.directory
    data_directory = os.path.dirname(os.path.abspath(__file__))
    if output_directory is None:
        output_directory = "%s/data" % data_directory

    igh_sequence = "%s/data/igh.fasta" % data_directory
    igh_specific_reference = "%s/reference.fasta" % output_directory
    blasr_index_directory = "%s/blasr_index" % output_directory
    pbmm2_index_directory = "%s/pbmm2_index" % output_directory

    #create_directory(blasr_index_directory)
    create_directory(pbmm2_index_directory)

    create_reference(args.hg19,igh_sequence,igh_specific_reference)
    #create_blasr_index(blasr_index_directory,igh_specific_reference)
    create_pbmm2_index(pbmm2_index_directory,igh_specific_reference)
