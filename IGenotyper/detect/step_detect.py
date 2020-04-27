#!/bin/env python
import pysam
import pybedtools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import namedtuple

from ..load import Step
from ..common import assembly_location,is_overlapping,get_haplotype,vcf_header,load_bed_regions
from ..common import show_value,create_directory,write_to_bashfile,extract_sequence_from,read_overlap_region
from ..command_line import run_assembly_scripts

class Variant(object):
    def __init__(self,file_manager):
        self.file_manager = file_manager

    def detect(self):
        pass

    def write(self):
        pass
        
    def labeled_hap_blocks(self,min_length=500,min_variants=2):
        labeled_phased_region = []
        label = 0
        Block = namedtuple('Block',['sample','chrom','start_1','start','end','num_variants'])
        with open(self.file_manager.haplotype_blocks,'r') as fh:
            header = fh.readline()
            for line in fh:
                line = line.rstrip().split('\t')
                block = Block._make(line)
                if int(block.num_variants) < min_variants:
                    continue
                if (int(block.end) - int(block.start)) < min_length:
                    continue
                labeled_phased_region.append([block.chrom,int(block.start),int(block.end),label])
                label += 1
        return labeled_phased_region

    def pos_in_bed_region(self,position,features):
        out_feature = "None"
        position = int(position)
        for chrom_region,start_region,end_region,feature in features:
            if position > int(start_region) and position <= int(end_region):
                out_feature = feature
        return out_feature

    def in_feature(self,position,bed_file):
        features = load_bed_regions(bed_file,True)
        return self.pos_in_bed_region(position,features)


class SV_type(object):
    def __init__(self,ct,ht,sv_start,sv_end,gt_start,gt_end,file_manager):
        self.coverage_threshold = ct
        self.hap_0_threshold = ht
        self.sv_start = sv_start
        self.sv_end = sv_end
        self.genotype_start = gt_start
        self.genotype_end = gt_end
        self.genotype = "1/1" # Not present
        self.file_manager = file_manager
        self.ref = "igh" # bad
        self.name = None

    def merge_intervals(self,intervals):
        intervals.sort()
        merged_intervals = []
        while len(intervals) > 0:
            if len(intervals) == 1:
                merged_intervals.append(intervals[0])
                intervals.pop(0)
                continue
            if intervals[0][1] >= intervals[1][0]:
                tmp = [intervals[0][0],max(intervals[0][1],intervals[1][1])]
                intervals[0] = tmp
                intervals.pop(1)
                continue
            merged_intervals.append(intervals[0])
            intervals.pop(0)
        return merged_intervals

    def genomic_regions_with_assembly(self):
        intervals = []
        samfile = pysam.AlignmentFile(self.file_manager.mapped_locus)
        for read in samfile.fetch(self.ref,self.sv_start,self.sv_end):
            if read.is_unmapped:
                continue
            if read.is_supplementary:
                continue
            if read.is_secondary:
                continue
            aligned_blocks = read.get_blocks()
            for aligned_block in aligned_blocks:
                start_block = aligned_block[0]
                end_block = aligned_block[1]
                start_interval = max(self.sv_start,start_block)
                end_interval = min(self.sv_end,end_block)
                if end_interval < start_interval:
                    continue
                intervals.append([start_interval,end_interval])
        merged_intervals = self.merge_intervals(intervals)
        return merged_intervals

    def get_prop_covered(self):
        intervals = self.genomic_regions_with_assembly()
        sum_of_coverage = 0.0
        for  start_interval, end_interval in intervals:
            sum_of_coverage += (end_interval - start_interval)
        prop_covered = sum_of_coverage/(self.sv_end - self.sv_start)
        return prop_covered

    def number_haplotype_reads(self):
        haplotype_count = {"0": 1.0, "1": 1.0, "2": 1.0}
        samfile = pysam.AlignmentFile(self.file_manager.phased_ccs_mapped_reads)
        for read in samfile.fetch(self.ref,self.genotype_start,self.genotype_end):
            if read.is_unmapped:
                continue
            if read.is_supplementary:
                continue
            if read.is_secondary:
                continue
            haplotype = read.get_tag("RG",True)[0]
            haplotype_count[haplotype] += 1.0
        return haplotype_count

    def haplotype_coverage(self):
        haplotype_count = self.number_haplotype_reads()
        total_reads = sum([haplotype_count[h] for h in haplotype_count])
        if total_reads > 10:
            hap_0 = haplotype_count["0"]/total_reads
            hap_1 = haplotype_count["1"]/total_reads
            hap_2 = haplotype_count["2"]/total_reads
        else:
            hap_0 = 0
            hap_1 = 0
            hap_2 = 0
        return (hap_0,hap_1,hap_2)

    def genotype_sv(self):
        prop_of_bp_covered = self.get_prop_covered()
        hap_0,hap_1,hap_2 = self.haplotype_coverage()
        if prop_of_bp_covered > self.coverage_threshold:
            if hap_0 > self.hap_0_threshold:
                self.genotype = "0/1" # het deletion --> homozygous reads
            else:
                self.genotype = "0/0" # hom presence --> heterzygous reads

class SV_1(SV_type):
    # 7-4-1
    def __init__(self,ct,ht,sv_start,sv_end,gt_start,gt_end,file_manager):
        super(SV_1,self).__init__(ct,ht,sv_start,sv_end,gt_start,gt_end,file_manager)
        self.name = "SV_1"

class SV_2A(SV_type):
    # V3-64D, V5-10-1
    def __init__(self,ct,ht,sv_start,sv_end,gt_start,gt_end,file_manager):
        super(SV_2A,self).__init__(ct,ht,sv_start,sv_end,gt_start,gt_end,file_manager)
        self.name = "SV_2A"

class SV_3(SV_type):
    # V3-23, V3-23D
    def __init__(self,ct,ht,sv_start,sv_end,gt_start,gt_end,file_manager):
        super(SV_3,self).__init__(ct,ht,sv_start,sv_end,gt_start,gt_end,file_manager)
        self.name = "SV_3"

class SV_4(SV_type):
    # V4-28,V3-30,V4-30-2,V3-30-3,V4-30-4,V3-30-5, V4-31,V3-33,V4-34
    def __init__(self):
        super(SV_4,self).__init__()
        self.name = "SV_4"

class SV_5(SV_type):
    # V3-38,V4-38-2,V3-43D,V3-38-3,V4-39,V1-38-4
    def __init__(self,ct,ht,sv_start,sv_end,gt_start,gt_end,file_manager):
        super(SV_5,self).__init__(ct,ht,sv_start,sv_end,gt_start,gt_end,file_manager)
        self.name = "SV_5"

class SV_6(SV_type):
    # V1-69,V2-70D,V1-69-2,V2-69D,V2-70
    def __init__(self,ct,ht,sv_start,sv_end,gt_start,gt_end,file_manager):
        super(SV_6,self).__init__(ct,ht,sv_start,sv_end,gt_start,gt_end,file_manager)
        self.name = "SV_6"

class SV_2b(SV_type):
    # V1-8,V3-9
    def __init__(self,ct,ht,sv_start,sv_end,gt_start,gt_end,file_manager):
        super(SV_2b,self).__init__(ct,ht,sv_start,sv_end,gt_start,gt_end,file_manager)
        self.name = "SV_2b"

class SV(Variant):
    def __init__(self,file_manager):
        super(SV,self).__init__(file_manager)
        self.sv_1 = SV_1(.9,.7,157688,167669,163368,165900,self.file_manager)
        self.sv_2a = SV_2A(.9,.9,210158,257471,223784,249767,self.file_manager)
        self.sv_3 = SV_3(.8,.8,412003,414775,412003,414775,self.file_manager)
        self.sv_4 = None
        self.sv_5 = SV_5(.7,.9,616403,679227,634584,660749,self.file_manager)
        self.sv_6 = SV_6(.7,.7,966204,1028309,989752,1001191,self.file_manager)
        self.sv_2b = SV_2b(.9,.9,1149070,1193130,1158302,1162596,self.file_manager)
        self.svs = [self.sv_1,self.sv_2a,self.sv_3,self.sv_5,self.sv_6,self.sv_2b]

    def detect(self):
        for sv in self.svs:
            sv.genotype_sv()

    def write(self):
        with open(self.file_manager.svs_genotyped,'w') as fh:
            for sv in self.svs:
                fh.write("%s\t%s\n" % (sv.name,sv.genotype))

class Snv(Variant):
    def __init__(self,file_manager,add_hom_ref_genotype):
        super(Snv,self).__init__(file_manager)
        self.hap_snps = {}
        self.genotypes = {}
        self.quality_scores = {}
        self.alt_base = {}
        self.alt_contigs = {}
        self.ref_base = {}
        self.snps_from_reads = {}
        self.add_hom_ref_genotype = add_hom_ref_genotype

    def snp_position_per_hap(self):
        samfile = pysam.AlignmentFile(self.file_manager.mapped_locus)
        ref = pysam.FastaFile(self.file_manager.blasr_ref)
        for read in samfile:
            if read.is_unmapped:
                continue
            if read.is_supplementary:
                continue
            if read.is_secondary:
                continue
            assembled_region = assembly_location(read.query_name)
            mapped_chrom = samfile.get_reference_name(read.reference_id)
            mapped_start = int(read.reference_start)
            mapped_end = int(read.reference_end)
            mapped_region = [mapped_chrom,mapped_start,mapped_end]
            if not is_overlapping(assembled_region,mapped_region):
                continue
            haplotype = get_haplotype(read.query_name)
            for read_pos, ref_pos in read.get_aligned_pairs():
                if read_pos == None:
                    continue
                if ref_pos == None:
                    continue
                ref_base = ref.fetch(mapped_chrom,ref_pos,ref_pos + 1).upper()
                read_base = read.query_sequence[read_pos].upper()
                read_qual = read.query_qualities[read_pos]
                if mapped_chrom not in self.hap_snps:
                    self.hap_snps[mapped_chrom] = {}
                if ref_pos not in self.hap_snps[mapped_chrom]:
                    self.hap_snps[mapped_chrom][ref_pos] = {}
                self.hap_snps[mapped_chrom][ref_pos][haplotype] = (ref_base,read_base,read_qual,read.query_name)

    def get_snp_genotype(self,ref,alt_bases,haps):
        if len(haps) == 1:
            if haps[0] == "0":
                if ref != alt_bases[0]:
                    genotype = "1/1"
                else:
                    genotype = "0/0"
            else:
                if ref != alt_bases[0]:
                    genotype = "./1"
                else:
                    genotype = "./0"
        else:
            if alt_bases[0] == alt_bases[1]:
                if alt_bases[0] == ref:
                    genotype = "0/0"
                else:
                    genotype = "1/1"
            else:
                if alt_bases[0] == ref:
                    genotype = "0/1"
                elif alt_bases[1] == ref:
                    genotype = "0/1"
                else:
                    genotype = "1/2"
        return genotype

    def get_alt_base(self,ref,alt_bases):
        alt_base = None
        if len(alt_bases) == 1:
            alt_base = alt_bases[0]
        else:
            if alt_bases[0] == alt_bases[1]:
                if alt_bases[0] == ref:
                    alt_base = ref
                else:
                    alt_base = alt_bases[0]
            else:
                if alt_bases[0] == ref:
                    alt_base = alt_bases[1]
                elif alt_bases[1] == ref:
                    alt_base = alt_bases[0]
                else:
                    alt_base = "%s,%s" % (alt_bases[0],alt_bases[1])
        return alt_base

    def get_snp_data(self):
        for chrom in self.hap_snps:
            if chrom not in self.genotypes:
                self.genotypes[chrom] = {}
            if chrom not in self.quality_scores:
                self.quality_scores[chrom] = {}
            if chrom not in self.alt_base:
                self.alt_base[chrom] = {}
            if chrom not in self.alt_contigs:
                self.alt_contigs[chrom] = {}
            if chrom not in self.ref_base:
                self.ref_base[chrom] = {}
            for pos in self.hap_snps[chrom]:
                haps = None
                alt_bases = None
                if "0" in self.hap_snps[chrom][pos]:
                    ref,alt_base,alt_qual,alt_name = self.hap_snps[chrom][pos]["0"]
                    haps = ["0"]
                    alt_bases = [alt_base]
                else:
                    if len(self.hap_snps[chrom][pos]) == 1:
                        hap = None
                        if "1" in self.hap_snps[chrom][pos]:
                            hap = "1"
                        if "2" in self.hap_snps[chrom][pos]:
                            hap = "2"
                        ref,alt_base,alt_qual,alt_name = self.hap_snps[chrom][pos][hap]
                        haps = [hap]
                        alt_bases = [alt_base]
                    else:
                        ref,hap_1_read_base,hap_1_read_qual,hap_1_read_name = self.hap_snps[chrom][pos]["1"]
                        ref,hap_2_read_base,hap_2_read_qual,hap_2_read_name = self.hap_snps[chrom][pos]["2"]
                        haps = ["1","2"]      
                        alt_bases = [hap_1_read_base,hap_2_read_base]
                        alt_qual = (hap_1_read_qual + hap_2_read_qual)/2.0                        
                        alt_name = "%s,%s" % (hap_1_read_name,hap_2_read_name)
                self.genotypes[chrom][pos] = self.get_snp_genotype(ref,alt_bases,haps)
                self.quality_scores[chrom][pos] = alt_qual
                self.alt_base[chrom][pos] = self.get_alt_base(ref,alt_bases)
                self.alt_contigs[chrom][pos] = alt_name
                self.ref_base[chrom][pos] = ref                        
    
    def get_snps_from_reads(self):
        with open(self.file_manager.phased_variants_vcf,'r') as vcf_fh:
            for line in vcf_fh:
                if line.startswith("#"):
                    continue
                line = line.rstrip().split("\t")
                if line[0] not in ["igh"]:
                    continue
                position = line[1]
                genotype = line[9].split(":")[0]
                self.snps_from_reads[int(position) - 1] = genotype
        return self.snps_from_reads

    def detect(self):
        self.snp_position_per_hap()
        self.get_snp_data()
        self.get_snps_from_reads()

    def get_info_field(self,chrom,pos,haplotype_blocks,genotype):
        if int(pos) in self.snps_from_reads:
            read_support = "Yes"
            phased_genotype = self.snps_from_reads[int(pos)]
            haplotype_block = self.pos_in_bed_region(pos,haplotype_blocks)
        else:
            read_support = "No"
            phased_genotype = "."
            haplotype_block = "."
        if genotype == "0/0" or genotype == "./0":
            read_support = "Yes" # fix
        contig = "contig=%s" % self.alt_contigs[chrom][pos]
        region = "region=%s" % self.in_feature(pos,self.file_manager.sv_regions)
        read_support = "read_support=%s" % read_support
        intronic = "intronic=%s" % self.in_feature(pos,self.file_manager.introns)
        lp1 = "LP1=%s" % self.in_feature(pos,self.file_manager.lpart1)
        rss = "RSS=%s" % self.in_feature(pos,self.file_manager.rss)
        gene = "gene=%s" % self.in_feature(pos,self.file_manager.gene_coordinates)
        igh_region = "igh_region=%s" % self.in_feature(pos,self.file_manager.region_types)
        phased_genotype = "phased_genotype=%s" % phased_genotype
        haplotype_block = "haplotype_block=%s" % haplotype_block
        sv_filter = "sv_filter=."
        info = [contig,region,read_support,intronic,lp1,rss,gene,igh_region,phased_genotype,haplotype_block,sv_filter]
        info = ";".join(info)
        return info

    def write(self):
        print "Calling SNPs %s.." % self.file_manager.assembly_snps
        haplotype_blocks = self.labeled_hap_blocks()
        output_lines = vcf_header()        
        with open(self.file_manager.assembly_snps,"w") as vcf_fh:
            vcf_fh.write("%s\n" % "\n".join(output_lines))
            for chrom in self.hap_snps:
                sorted_positions = self.hap_snps[chrom].keys()
                sorted_positions.sort()
                for pos in sorted_positions:
                    genotype = self.genotypes[chrom][pos]        
                    if genotype == "0/0" and (not self.add_hom_ref_genotype):
                        continue
                    if genotype == "./0" and (not self.add_hom_ref_genotype):
                        continue
                    rs_id= "."
                    ref = self.ref_base[chrom][pos]
                    alt_base = self.alt_base[chrom][pos]
                    alt_qual = self.quality_scores[chrom][pos]
                    info = self.get_info_field(chrom,pos,haplotype_blocks,genotype)
                    out = [chrom,pos,rs_id,ref,alt_base,alt_qual,"PASS",info,"GT",genotype]
                    vcf_fh.write("%s\n" % "\t".join(map(str,out)))

class IndelScript(object):
    def __init__(self,interval,file_manager):
        self.chrom = interval[0]
        self.start = interval[1]
        self.end = interval[2]
        self.file_manager = file_manager

    def create_script(self):
        directory = "%s/%s/%s_%s" % (self.file_manager.indels,self.chrom,self.start,self.end)
        bashfile = "%s/assemble.sh" % directory
        params = {
            'dir': directory,
            'python_scripts': self.file_manager.python_scripts,
            'chrom': self.chrom,
            'start': self.start,
            'end': self.end
        }
        write_to_bashfile(self.file_manager.sv_calling_script,bashfile,params)
        return bashfile

class Indels(Variant):
    def __init__(self,file_manager,cpu_manager):
        super(Indels,self).__init__(file_manager)
        self.cpu_manager = cpu_manager
        self.length_thres = 100
        self.igh_region_coord = "igh 1 1193129"
        self.regions_to_detect_indels = {}
        self.contig_coord = {}
        self.scripts = []

    def get_contig_coord(self):
        samfile = pysam.AlignmentFile(self.file_manager.mapped_locus)
        for contig in samfile:
            if contig.is_secondary:
                continue
            if contig.is_supplementary:
                continue
            if contig.is_unmapped:
                continue
            ref = samfile.get_reference_name(contig.reference_id)
            start = contig.reference_start
            end = contig.reference_end
            hap = get_haplotype(contig.query_name)
            if hap not in self.contig_coord:
                self.contig_coord[hap] = {}
            self.contig_coord[hap][(ref,start,end)] = contig        
                
    def add_region_to_detect(self,hap,region):
        self.regions_to_detect_indels[hap] = []
        coords = pybedtools.BedTool(self.contig_coord[hap].keys()).intersect(region)
        coords_coverage = coords.genomecov(bg=True,g="%s" % self.file_manager.igh_fasta_fai)
        for chrom,start,end,coverage in coords_coverage:
            if int(show_value(end)) - int(show_value(start)) < self.length_thres:
                continue
            if int(show_value(coverage)) > 1:
                continue
            self.regions_to_detect_indels[hap].append((chrom,start,end))

    def get_overlapping_regions(self):
        hap1 = pybedtools.BedTool(self.regions_to_detect_indels["1"])
        hap2 = pybedtools.BedTool(self.regions_to_detect_indels["2"])
        overlapping_regions = hap1.intersect(hap2)
        for hap in ["1","2"]:
            self.regions_to_detect_indels[hap] = []
            for chrom,start,end in overlapping_regions:
                overlapping_region = (show_value(chrom),show_value(start),show_value(end))
                self.regions_to_detect_indels[hap].append(overlapping_region)

    def get_coordinates_for_msa(self):
        haplotype_blocks = self.labeled_hap_blocks()        
        locus = pybedtools.BedTool(self.igh_region_coord,from_string=True)
        phased_regions = pybedtools.BedTool(haplotype_blocks)
        unphased_regions = locus.subtract(phased_regions)
        self.add_region_to_detect("0",unphased_regions)
        self.add_region_to_detect("1",phased_regions)
        self.add_region_to_detect("2",phased_regions)
        self.get_overlapping_regions()
        # write coordinates to file  

    def get_assembly_sequence(self,chrom,start,end,hap):
        hap_sequences = []
        hap_quals = []
        contigs_names = []
        samfile = pysam.AlignmentFile(self.file_manager.mapped_locus,'rb')
        for contig in samfile.fetch(chrom,start,end):
            if contig.is_unmapped:
                continue
            if contig.is_secondary:
                continue
            if contig.is_supplementary:
                continue
            if contig.reference_start > start:
                continue
            if contig.reference_end < end:
                continue
            contig_hap = get_haplotype(contig.query_name)
            if hap == "0":
                if contig_hap != "0":
                    continue
            else:
                if contig_hap == "0":
                    continue
            aligned_pairs = contig.get_aligned_pairs()
            query_start = None
            query_end = None
            matched_ref_start = None
            matched_ref_end = None
            for query_pos, ref_pos in aligned_pairs:
                if query_pos == None:
                    continue
                if ref_pos == None:
                    continue
                if int(ref_pos) <= int(start):
                    query_start = query_pos
                    matched_ref_start = ref_pos
                query_end = query_pos
                matched_ref_end = ref_pos
                if int(ref_pos) > int(end):
                    break
            assert query_start != None
            assert query_end != None
            hap_sequence = contig.query_sequence[query_start:query_end]
            hap_sequences.append(hap_sequence)
            sequence_qual = contig.query_qualities[query_start:query_end]
            hap_quals.append(sequence_qual)
            contigs_names.append(contig.query_name)
            assert len(sequence_qual) == len(hap_sequence)
        return (hap_sequences,hap_quals,contigs_names)

    def get_sequences_for_msa(self):
        fasta = pysam.FastaFile(self.file_manager.blasr_ref)
        for hap in self.regions_to_detect_indels:
            if hap == "2":
                continue # only need to go through hap1 or hap2, not both    
            for chrom,start,end in self.regions_to_detect_indels[hap]:
                start = int(start)
                end = int(end)
                assembly_sequences,hap_quals,contigs_names = self.get_assembly_sequence(chrom,start,end,hap)
                ref_seq = fasta.fetch(reference=chrom,start=max(1,start),end=end)
                hap1 = assembly_sequences[0]
                qual1 = hap_quals[0]
                h1_contig = contigs_names[0]
                if hap == "0":
                    assert len(assembly_sequences) == 1
                    hap2 = assembly_sequences[0]
                    qual2 = hap_quals[0]
                    h2_contig = contigs_names[0]
                else:
                    assert len(assembly_sequences) == 2
                    hap2 = assembly_sequences[1]
                    qual2 = hap_quals[1]
                    h2_contig = contigs_names[1]
                directory = "%s/%s/%s_%s" % (self.file_manager.indels,chrom,start,end)
                create_directory(directory)
                outseqfn = "%s/seq.fa" % directory
                with open(outseqfn,'w') as outfh:
                    outfh.write(">ref\n%s\n" % ref_seq)
                    outfh.write(">hap1\n%s\n" % hap1)
                    outfh.write(">hap2\n%s\n" % hap2)
                outqualfn = "%s/seq.qual" % directory
                with open(outqualfn,'w') as outqualfh:
                    outqualfh.write(">hap1_%s\n%s\n" % (h1_contig,",".join(map(str,qual1))))
                    outqualfh.write(">hap2_%s\n%s\n" % (h2_contig,",".join(map(str,qual2))))


    def run_msa(self):
        for hap in self.regions_to_detect_indels:
            if hap == "2":
                continue # only need to go through hap1 or hap2, not both    
            for interval in self.regions_to_detect_indels[hap]:
                indel_script = IndelScript(interval,self.file_manager)
                self.scripts.append(indel_script.create_script())
        run_assembly_scripts(self.scripts,self.cpu_manager.cluster,
                             self.cpu_manager.cluster_walltime,
                             self.cpu_manager.threads,self.cpu_manager.cluster_mem,
                             self.cpu_manager.cluster_queue)        

    def write(self):
        variants = []
        for hap in self.regions_to_detect_indels:
            if hap == "2":
                continue # only need to go through hap1 or hap2, not both    
            for chrom,start,end in self.regions_to_detect_indels[hap]:
                variant_bed = "%s/%s/%s_%s/svs.bed" % (self.file_manager.indels,chrom,start,end)
                with open(variant_bed,'r') as variant_bedfh:
                    for line in variant_bedfh:
                        line = line.strip().split('\t')
                        line[1] = int(line[1])
                        variants.append(line)
        variants.sort(key=lambda x: x[1])
        assembly_indels_fh = open(self.file_manager.assembly_indels,'w')
        assembly_svs_fh = open(self.file_manager.assembly_svs,'w')
        for variant in variants:
            size = int(variant[5])
            if size < 51:
                assembly_indels_fh.write("%s\n" % "\t".join(map(str,variant)))
            else:
                assembly_svs_fh.write("%s\n" % "\t".join(map(str,variant)))
        assembly_indels_fh.close()
        assembly_svs_fh.close()

    def detect(self):
        self.get_contig_coord()
        self.get_coordinates_for_msa()
        self.get_sequences_for_msa()        
        self.run_msa()

class AllelesFrom(object):
    def __init__(self,gene_coordinates,database,input_alignment,output_gene_sequence,is_assembly=False):
        self.gene_coordinates = gene_coordinates
        self.allele_database = database
        self.input_alignment = input_alignment
        self.output_gene_sequence = output_gene_sequence
        self.is_assembly = is_assembly
        self.sequence = {}
        self.matches = {}
        self.novel_genes = {}
        self.genes_to_alleles = {}

    def write_gene_sequence_to_file(self): 
        out_sequences = []
        for gene in self.sequence:
            for haplotype in self.sequence[gene]:
                total = len(self.sequence[gene][haplotype])
                for i,seq in enumerate(self.sequence[gene][haplotype]):
                    name = "gene=%s.hap=%s.index=%s.total=%s" % (gene,haplotype,i,total)
                    record = SeqRecord(Seq(seq),id=name,name=name,description="")
                    out_sequences.append(record)
        SeqIO.write(out_sequences,self.output_gene_sequence,"fasta")

    def extract_genes_from_sequence(self):
        gene_coordinates = load_bed_regions(self.gene_coordinates,True)
        samfile = pysam.AlignmentFile(self.input_alignment)
        for chrom,start,end,gene in gene_coordinates:
            start = int(start)
            end = int(end)
            for read in samfile.fetch(chrom,start,end):
                if read.is_secondary:
                    continue
                if read.is_supplementary:
                    continue
                if read.is_unmapped:
                    continue
                mapping_cord = [chrom,read.reference_start,read.reference_end]
                if self.is_assembly:
                    if not read_overlap_region(read,mapping_cord):
                        continue
                gene_sequence = extract_sequence_from(read,chrom,start,end)
                if gene_sequence == "":
                    continue
                if len(gene_sequence) < (end - start - 10):
                    continue
                if len(gene_sequence) > ((end - start) + 10):
                    continue
                if self.is_assembly:
                    haplotype = get_haplotype(read.query_name)
                else:
                    haplotype = read.get_tag("RG",True)[0]
                if gene not in self.sequence:
                    self.sequence[gene] = {}
                if haplotype not in self.sequence[gene]:
                    self.sequence[gene][haplotype] = []
                self.sequence[gene][haplotype].append(gene_sequence)
        self.write_gene_sequence_to_file() #self.sequence,self.output_gene_sequence)

    def get_matches(self):
        query_entries = SeqIO.to_dict(SeqIO.parse(self.output_gene_sequence,"fasta"))
        database_entries = SeqIO.to_dict(SeqIO.parse(self.allele_database,"fasta"))
        for query_entry in query_entries:
            query_seq = query_entries[query_entry].seq.upper()
            self.matches[(query_entry,query_seq)] = set()
            for database_entry in database_entries:
                database_seq = database_entries[database_entry].seq.upper()
                if query_seq.count(database_seq) > 0:
                    self.matches[(query_entry,query_seq)].add(database_entry)
                if query_seq.count(database_seq.reverse_complement()) > 0:
                    self.matches[(query_entry,query_seq)].add(database_entry)

    def is_novel_allele(self,gene_extraction_name,gene_seq):
        if len(self.matches[(gene_extraction_name, gene_seq)]) == 0:
            return True
        return False

    def add_novel_allele(self,gene_name,haplotype,gene_seq):
        if gene_name not in self.novel_genes:
            self.novel_genes[gene_name] = {}
        if haplotype not in self.novel_genes[gene_name]:
            self.novel_genes[gene_name][haplotype] = []
        self.novel_genes[gene_name][haplotype].append(gene_seq)        

    def get_assigned_allele(self):
        for gene_extraction_name, gene_seq in self.matches:
            gene_name = gene_extraction_name.split(".")[0].split("=")[1]
            haplotype = gene_extraction_name.split(".")[1].split("=")[1]
            for allele_hit in self.matches[(gene_extraction_name, gene_seq)]:
                allele_hit_gene = allele_hit.split("_")[0].split("=")[1]
                if allele_hit_gene == gene_name:
                    allele_hit = allele_hit.split("_")[1].split("=")[1]
                if gene_name not in self.genes_to_alleles:
                    self.genes_to_alleles[gene_name] = {}
                if haplotype not in self.genes_to_alleles[gene_name]:
                    self.genes_to_alleles[gene_name][haplotype] = []
                self.genes_to_alleles[gene_name][haplotype].append(allele_hit)        

    def get_novel_alleles(self):
        for gene_extraction_name, gene_seq in self.matches:
            gene_name = gene_extraction_name.split(".")[0].split("=")[1]
            haplotype = gene_extraction_name.split(".")[1].split("=")[1]
            if self.is_novel_allele(gene_extraction_name,gene_seq):
                self.add_novel_allele(gene_name,haplotype,gene_seq)        
        
    def match_gene_to_allele_db(self):
        self.get_matches()
        self.get_novel_alleles()
        self.get_assigned_allele()

    def get_allele_count(self,gene,hap,allele):
        count = 0
        for location in [self.genes_to_alleles,self.novel_genes]:
            if gene not in location:
                continue
            if hap not in location[gene]:                    
                continue
            for allele_hit in location[gene][hap]:
                if allele == allele_hit:
                    count += 1
        return count

    def get_allele(self,gene,hap):
        alleles = None
        if gene in self.genes_to_alleles:
            if hap in self.genes_to_alleles[gene]:
                alleles = self.genes_to_alleles[gene][hap]
        if alleles == None:
            if gene in self.novel_genes:
                if hap in self.novel_genes[gene]:
                    alleles = [str(self.novel_genes[gene][hap][0])]
        if alleles == None:
            alleles = ["NotDetermined"]
        return list(set(alleles))

class Alleles(Variant):
    def __init__(self,file_manager):
        super(Alleles,self).__init__(file_manager)
        self.alleles_from_assembly = AllelesFrom(self.file_manager.gene_coordinates,
                                                 self.file_manager.allele_database,
                                                 self.file_manager.mapped_locus,
                                                 self.file_manager.genes_from_assembly,
                                                 True)
        self.alleles_from_ccs_reads = AllelesFrom(self.file_manager.gene_coordinates,
                                                  self.file_manager.allele_database,
                                                  self.file_manager.phased_ccs_mapped_reads,
                                                  self.file_manager.genes_from_reads)

    def detect(self):
        for sequence_type in [self.alleles_from_assembly,self.alleles_from_ccs_reads]:
            sequence_type.extract_genes_from_sequence()
            sequence_type.match_gene_to_allele_db()                                                
            
    def get_haplotype_alleles(self,chrom,start,end,gene):
        output = [chrom,int(start),int(end),gene]
        gene_in_sv = self.in_feature(int(start),self.file_manager.sv_regions)            
        haplotype_blocks = self.labeled_hap_blocks()
        haplotype_block = self.pos_in_bed_region(int(start),haplotype_blocks)
        if haplotype_block != "None":
            haps = ["1","2"]
        else:
            haps = ["0"]
        for hap in haps:
            alleles = self.alleles_from_assembly.get_allele(gene,hap)
            if alleles[0] == "NotDetermined":
                alleles = self.alleles_from_ccs_reads.get_allele(gene,hap)
            allele_read_count = []
            for allele in alleles:
                allele_read_count.append(self.alleles_from_ccs_reads.get_allele_count(gene,hap,allele))
            alleles = ",".join(alleles)
            allele_read_count = ",".join(map(str,allele_read_count))
            output += [alleles,allele_read_count]
            if hap == "0" and gene_in_sv == "None":
                output += [alleles,allele_read_count]
            if hap == "0" and gene_in_sv != "None":
                output += ["Deleted","."]
        return output
        
    def write(self):
        header = ["chrom",
                  "start",
                  "end",
                  "gene_name",
                  "haplotype_1_allele",
                  "num_of_ccs_reads_support_for_haplotype_1",
                  "haplotype_2_allele",
                  "num_of_ccs_reads_support_for_haplotype_2",
                  "haplotype_block"]
        genes = load_bed_regions(self.file_manager.gene_coordinates,True)
        outputs = []
        for chrom,start,end,gene in genes:
            output = self.get_haplotype_alleles(chrom,start,end,gene)
            outputs.append(output)
        outputs.sort(key = lambda x: x[1])
        with open(self.file_manager.genes_with_allele_assignment,'w') as fh:
            fh.write("%s\n" % "\t".join(header))
            for output in outputs:
                fh.write("%s\n" % "\t".join(map(str,output)))
        
        
class DetectVariants(Step):
    def __init__(self,file_manager,cpu_manager,command_line_tools):
        super(DetectVariants,self).__init__(file_manager,cpu_manager,command_line_tools)
        self.done_file = "%s/detect.txt" % self.file_manager.outdir        

    def run(self):
        svs = SV(self.file_manager)
        snvs = Snv(self.file_manager,self.add_hom_ref_genotype)
        indels = Indels(self.file_manager,self.cpu_manager)
        alleles = Alleles(self.file_manager)
        for variant in [svs,snvs,indels,alleles]:
            variant.detect()
            variant.write()

