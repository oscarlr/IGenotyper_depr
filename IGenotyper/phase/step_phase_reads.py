#!/bin/env python
import os
import sys
import pysam
import numpy
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


from ..common import read_is_unphased,file_paths,load_bed_regions
from ..load import Step

def supplementary_score_diff(read,secondary_reads,thres):
    diffs = [] 
    primary_read_score = float(read.get_tag("AS",True)[0])
    for secondary_read in secondary_reads:
        secondary_read_score = float(secondary_read.get_tag("AS",True)[0])
        diffs.append(abs(primary_read_score-secondary_read_score))
    if min(diffs) > thres:
        return True
    return False

def change_read(primary_read,secondary_reads):
    secondary_to_primary_read = None
    secondary_to_primary_flag = None
    secondary_to_primary_mapq = None
    index_to_ignore = None
    secondary_to_primary_score = 0
    for i,secondary_read in enumerate(secondary_reads):
        secondary_read_score = secondary_read.get_tag("AS",True)[0]
        if secondary_read_score < secondary_to_primary_score:
            secondary_to_primary_read = secondary_read
            secondary_to_primary_score = secondary_read_score
            secondary_to_primary_flag = secondary_read.flag
            secondary_to_primary_mapq = secondary_read.mapping_quality
            index_to_ignore = i
    # Change flag to primary alignment
    #secondary_to_primary_read.flag = primary_read.flag
    secondary_to_primary_read.flag = secondary_to_primary_read.flag - 256
    secondary_to_primary_read.mapping_quality = primary_read.mapping_quality
    # Change flag to secondary alignment
    #primary_read.flag = secondary_to_primary_flag
    primary_read.flag = primary_read.flag + 256
    primary_read.mapping_quality = secondary_to_primary_mapq
    #reads_to_return = [primary_read,secondary_to_primary_read]
    reads_to_return = [secondary_to_primary_read]
    for i,read in enumerate(secondary_reads):
        if i == index_to_ignore:
            continue
        reads_to_return.append(read)
    return reads_to_return

class Phase(Step):
    def __init__(self,file_manager,cpu_manager,command_line_tools):
        super(Phase,self).__init__(file_manager,cpu_manager,command_line_tools)
        self.done_file = "%s/phasing.txt" % self.file_manager.outdir
        self.files_to_check = [self.file_manager.phased_ccs_mapped_reads,
                               self.file_manager.phased_subreads_mapped_reads,
                               self.file_manager.phased_variants_vcf,
                               self.file_manager.variants_vcf]
        self.number_of_on_target_subreads = 0
        self.number_of_on_target_ccs_reads = 0
        self.igh_ccs_coverage = 0
        self.igh_subread_coverage = 0
        self.phased_gene_coverage = {}
        self.phased_sv_coverage = {}
        self.phased_region_coverage = {}
        self.aim_coverage = {}
        self.output = []

    def get_initial_phasing(self):
        self.command_line_tools.get_ccs_reads()
        self.command_line_tools.turn_ccs_reads_to_fastq()
        self.command_line_tools.map_ccs_reads()
        self.command_line_tools.map_subreads()
        self.command_line_tools.phase_snps()
        self.command_line_tools.phase_ccs_reads()   
        self.command_line_tools.phase_subreads()             

    def get_secondary_alignments(self,bam):
        secondary_alignments = {}
        samfile = pysam.AlignmentFile(bam)
        for read in samfile.fetch("igh"):
            if read.is_secondary:
                if read.query_name not in secondary_alignments:
                    secondary_alignments[read.query_name] = []
                secondary_alignments[read.query_name].append(read)
        samfile.close()
        return secondary_alignments

    def change_primary_alignments(self,output_bamfile,bam):
        secondary_alignments = self.get_secondary_alignments(bam)
        samfile = pysam.AlignmentFile(bam)
        changed_reads = set()
        for read in samfile.fetch("igh"):
            if read.is_secondary:
                continue
            if read.is_supplementary:
                continue
            if not read_is_unphased(read):
                continue
            if read.query_name in secondary_alignments:
                secondary_reads = secondary_alignments[read.query_name]
                # Always accept mapping quality less than 40
                if read.mapping_quality > 40:
                    # If the scores is similar,accept                                                
                    if supplementary_score_diff(read,secondary_reads,self.secondary_read_score):
                        continue
                output_reads = change_read(read,secondary_reads)
                for output_read in output_reads:
                    changed_reads.add(output_read.query_name)
                    output_bamfile.write(output_read)
        samfile.close()
        return changed_reads

    def get_changed_bamfile(self,bam):
        changed_bamfile = "%s/changed_alignments.sam" % self.file_manager.tmp_dir
        samfile = pysam.AlignmentFile(bam)
        output_bamfile = pysam.AlignmentFile(changed_bamfile,"w",template=samfile)
        samfile.close()        
        return output_bamfile

    def fix_alignments(self,bam):
        # 1. Find the unphased reads
        # 2. Check if the secondary alignment of the unphased read is phased
        # 3. Change the secondary alignment to the primary alignment
        # 4. If the scondary alignment is still not phased, move it to the primary alignment
        changed_bamfile = "%s/changed_alignments.sam" % self.file_manager.tmp_dir
        output_bamfile = self.get_changed_bamfile(bam)
        changed_reads = self.change_primary_alignments(output_bamfile,bam)
        samfile = pysam.AlignmentFile(bam)
        for read in samfile:
            if read.query_name in changed_reads:
                continue
            output_bamfile.write(read)
        samfile.close()        
        output_bamfile.close()
        prefix = changed_bamfile[:-4]
        self.command_line_tools.sam_to_sorted_bam(prefix,"%s.sorted.bam" % prefix)

    def save_previous_variants(self,index):
        src = "%s/variants/from_reads" % self.file_manager.outdir
        dst = "%s/tmp/variants/from_reads_%s" % (self.file_manager.outdir,index)
        os.mkdir(dst)
        os.rename(src,dst)
        os.mkdir("%s/variants/from_reads" % self.file_manager.outdir)
        
    def save_previous_ccs_mappings(self,index):
        file_names = ["ccs_to_ref_phased.sorted.bam","ccs_to_ref_phased.sorted.bam.bai",
                      "ccs_to_ref.sorted.bam","ccs_to_ref.sorted.bam.bai"]
        os.mkdir("%s/tmp/alignments/run_%s" % (self.file_manager.outdir,index))
        for file_name in file_names:
            src = "%s/alignments/%s" % (self.file_manager.outdir,file_name)
            dst = "%s/tmp/alignments/run_%s/%s" % (self.file_manager.outdir,index,file_name)
            os.rename(src,dst)        
        src = "%s/tmp/changed_alignments.sorted.bam" % self.file_manager.outdir
        dst = "%s/alignments/ccs_to_ref.sorted.bam" % self.file_manager.outdir
        os.rename(src,dst)
        src = "%s/tmp/changed_alignments.sorted.bam.bai" % self.file_manager.outdir
        dst = "%s/alignments/ccs_to_ref.sorted.bam.bai" % self.file_manager.outdir
        os.rename(src,dst)        

    def remove_subreads_alignments(self,index):
        src = "%s/tmp/changed_alignments.sorted.bam" % self.file_manager.outdir
        dst = "%s/alignments/subreads_to_ref.sorted.bam" % self.file_manager.outdir
        os.rename(src,dst)
        src = "%s/tmp/changed_alignments.sorted.bam.bai" % self.file_manager.outdir
        dst = "%s/alignments/subreads_to_ref.sorted.bam.bai" % self.file_manager.outdir
        os.rename(src,dst)        
        fn = "%s/alignments/subreads_to_ref_phased.sorted.bam" % self.file_manager.outdir
        os.remove(fn)
        os.remove("%s.bai" % fn)

    def save_previous_files(self,index,bam):
        if bam == self.file_manager.phased_ccs_mapped_reads:
            self.save_previous_variants(index)
            self.save_previous_ccs_mappings(index)
        else:
            self.remove_subreads_alignments(index)

    def rephase(self):
        for iter_ in ["0","1"]:
            print "Rephasing iteration %s" % iter_
            iter_dir = "%s/tmp/variants/from_reads_%s" % (self.file_manager.outdir,iter_)
            if os.path.isdir(iter_dir):
                continue
            for bam in [self.file_manager.phased_ccs_mapped_reads,
                        self.file_manager.phased_subreads_mapped_reads]:
                self.fix_alignments(bam)
                self.save_previous_files(iter_,bam)
            self.command_line_tools.phase_snps()
            self.command_line_tools.phase_ccs_reads()    
            self.command_line_tools.phase_subreads()

    def get_read_lengths(self,bam_file):
        lengths = []
        samfile = pysam.AlignmentFile(bam_file, "rb",check_sq=False)
        for read in samfile:
            lengths.append(read.query_length)
        return lengths

    def get_reads_stats(self):
        subread_lengths = self.get_read_lengths(self.file_manager.input_bam)
        ccs_lengths = self.get_read_lengths(self.file_manager.ccs_reads)        
        self.output.append(["reads","Number of subreads",len(subread_lengths)])
        self.output.append(["reads","Number of ccs reads",len(ccs_lengths)])
        self.output.append(["reads","Average subread length",numpy.mean(subread_lengths)])
        self.output.append(["reads","Average ccs length",numpy.mean(ccs_lengths)])
        self.output.append(["reads","Median subread length",numpy.median(subread_lengths)])
        self.output.append(["reads","Median ccs length",numpy.median(ccs_lengths)])

    def get_variants_stats(self):
        number_of_het_snvs = 0
        number_of_hom_snvs = 0
        number_of_phased_snvs = 0
        number_of_ref_bp_phased = 0
        with open(self.file_manager.variants_vcf,'r') as fh:
            for line in fh:
                line = line.rstrip()
                if line[0] == "#":
                    continue
                if "1/0" in line:
                    number_of_het_snvs += 1
                if "0/1" in line:
                    number_of_het_snvs += 1
                if "1/1" in line:
                    number_of_hom_snvs += 1
        with open(self.file_manager.phased_variants_vcf,'r') as fh:
            for line in fh:
                line = line.rstrip()
                if line[0] == "#":
                    continue
                if "0|1" in line:
                    number_of_phased_snvs += 1
                if "1|0" in line:
                    number_of_phased_snvs += 1
        with open(self.file_manager.haplotype_blocks,'r') as fh:
            for line in fh:
                line = line.rstrip().split('\t')
                if "#" in line[0]:
                    continue
                number_of_ref_bp_phased += (int(line[4]) - int(line[3])) 
        self.output.append(["snv","Number of heterozygous snv",number_of_het_snvs])
        self.output.append(["snv","Number of homozygous snv",number_of_hom_snvs])
        self.output.append(["snv","Number of phased heterozygous snv",number_of_phased_snvs])
        self.output.append(["snv","Number of references bases phased",number_of_ref_bp_phased])

        
    def get_phased_reads_from(self,bam_file,chrom,start,end):
        read_lengths = {}
        samfile = pysam.AlignmentFile(bam_file, "rb")
        for read in samfile.fetch(chrom,start,end):
            if read.is_secondary:
                continue
            if read.is_unmapped:
                continue
            if read.is_supplementary:
                continue
            hap = read.get_tag("RG",True)[0]
            if hap not in read_lengths:
                read_lengths[hap] = []
            length = max(0, min(end, int(read.reference_end)) - max(start, int(read.reference_start)))
            read_lengths[hap].append(length)
        return read_lengths

    def get_igh_phased_reads(self,bam_file):
        chrom = "igh"
        start = 1
        end = 1193129
        read_lengths = self.get_phased_reads_from(bam_file,chrom,start,end)
        return read_lengths
        
    def on_target_reads(self,read_lengths):
        reads = 0
        for hap in read_lengths:
            reads += len(read_lengths[hap])
        return reads

    def get_igh_coverage(self,read_lengths):
        sum_read_lengths = 0
        for hap in read_lengths:
            sum_read_lengths += numpy.sum(read_lengths[hap])
        return sum_read_lengths/1193129.0

    def get_phased_region_coverage(self,bedfile,name):
        phased_coverage = {}
        features = load_bed_regions(bedfile,True)
        for chrom,start,end,feature in features:
            reads = self.get_phased_reads_from(self.file_manager.phased_ccs_mapped_reads,chrom,start,end)
            output_line = [name,chrom,start,end,feature]
            for hap in ["0","1","2"]:
                coverage = 0
                if hap in reads:
                    coverage = numpy.sum(reads[hap])/(end - start)
                output_line.append(coverage)
            self.output.append(output_line)

    def get_aligment_stats(self):
        igh_subreads = self.get_igh_phased_reads(self.file_manager.phased_subreads_mapped_reads)
        igh_ccs = self.get_igh_phased_reads(self.file_manager.phased_ccs_mapped_reads)
        self.output.append(["alignment","Number of on-target subreads",self.on_target_reads(igh_subreads)])
        self.output.append(["alignment","Number of on-target ccs reads",self.on_target_reads(igh_ccs)])
        self.output.append(["alignment","IGH ccs coverage",self.get_igh_coverage(igh_ccs)])
        self.output.append(["alignment","IGH subread coverage",self.get_igh_coverage(igh_subreads)])
        self.get_phased_region_coverage(self.file_manager.gene_coordinates,"gene")
        self.get_phased_region_coverage(self.file_manager.sv_regions,"sv")
        self.get_phased_region_coverage(self.file_manager.region_types,"region")
        self.get_phased_region_coverage(self.file_manager.aims,"aims")
        
    def get_phase_stats(self):
        if self.step_completed():
            self.get_reads_stats()
            self.get_variants_stats()
            self.get_aligment_stats()
            with open(self.done_file,"w") as fh:                
                for output in self.output:
                    fh.write("%s\n" % "\t".join(map(str,output)))

    def run(self):
        self.get_initial_phasing()
        self.rephase()
        self.command_line_tools.get_phased_blocks()
        self.get_phase_stats()
        self.clean_up()
        
