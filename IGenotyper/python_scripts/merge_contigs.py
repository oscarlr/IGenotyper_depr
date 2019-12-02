#!/bin/env python
import sys
import pysam

def score_contigs(bamfile,c1,c2,ref,start,end):
    samfile = pysam.AlignmentFile(bamfile)
    mismatches = 0.0
    comparisons = 0.0
    for pileupcolumn in samfile.pileup(ref, start, end):
        c1_base = None
        c2_base = None
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                if pileupread.alignment.query_name == c1:
                    c1_base = pileupread.alignment.query_sequence[pileupread.query_position]
                if pileupread.alignment.query_name == c2:
                    c2_base = pileupread.alignment.query_sequence[pileupread.query_position]
                if None in [c1_base,c2_base]:
                    continue
                if c1_base != c2_base:
                    mismatches += 1.0
                comparisons += 1.0
    #print (comparisons-mismatches)/comparisons,mismatches,comparisons
    if comparisons > 100:
        return mismatches,(comparisons-mismatches)/comparisons
    return None,None

def get_contained_contigs(contig_coords,bamfile):
    contained_contigs = []
    for contig1 in contig_coords:
        c1_ref,c1_start,c1_end = contig_coords[contig1]
        for contig2 in contig_coords:
            if contig1 == contig2:
                continue
            c2_ref,c2_start,c2_end = contig_coords[contig2]
            if c2_start > c1_start and c2_end < c1_end:
                mismatches,alignment_score = score_contigs(bamfile,contig1,contig2,c1_ref,c1_start,c1_end)
                if alignment_score == None:
                    continue
                if alignment_score == 1 or mismatches < 3:
                    contained_contigs.append(contig2)
    return contained_contigs

def is_overlapping(a, b):
    if a[0] != b[0]: # chrom not matching
        return False
    overlapping = False
    num_overlapping = max(0, min(a[2], b[2]) - max(a[1], b[1]))
    if num_overlapping > 0:
        overlapping = True
    return overlapping

def get_overlapping_contigs(contig_coords,contained_contigs):
    overlapping_contigs = []
    for contig1 in contig_coords:
        if contig1 in contained_contigs:
            continue
        for contig2 in contig_coords:
            if contig1 == contig2:
                continue
            if sorted([contig1,contig2]) in overlapping_contigs:
                continue
            if contig2 in contained_contigs:
                continue
            if not is_overlapping(contig_coords[contig1],contig_coords[contig2]):
                continue                
            start = max(contig_coords[contig1][1],contig_coords[contig2][1])
            end = min(contig_coords[contig1][2],contig_coords[contig2][2])
            ref = contig_coords[contig1][0]
            mismatches,alignment_score = score_contigs(bamfile,contig1,contig2,ref,start,end)
            if alignment_score == None:
                continue
            if alignment_score == 1 or mismatches < 3:
                overlapping_contigs.append(sorted([contig1,contig2]))
    return overlapping_contigs

def group_overlapping_contigs(overlapping_contigs):
    groups = {}
    current_group = 0
    for c1,c2 in overlapping_contigs:
        c1_or_c2_in_groups = []
        for group in groups:
            if c1 in groups[group] or c2 in groups[group]:
                c1_or_c2_in_groups.append(group)                
                merge_group = group
        if len(c1_or_c2_in_groups) == 1:
            groups[merge_group].add(c1)
            groups[merge_group].add(c2)
        elif len(c1_or_c2_in_groups) > 1:
            groups[current_group] = set()
            groups[current_group].add(c1)
            groups[current_group].add(c2)
            for group_to_be_merged in c1_or_c2_in_groups:
                for c in groups[group_to_be_merged]:
                    groups[current_group].add(c)
                del groups[group_to_be_merged]
            current_group += 1
        else:
            groups[current_group] = set()
            groups[current_group].add(c1)
            groups[current_group].add(c2)
            current_group += 1
    return groups

def get_group_ref_coords(bamfile,contigs):
    samfile = pysam.AlignmentFile(bamfile)
    contig_coords = []
    for contig in samfile:
        if contig.is_secondary:
            continue
        if contig.is_unmapped:
            continue
        if contig.is_supplementary:
            continue
        if contig.query_name not in contigs:
            continue
        ref_start = contig.reference_start
        ref_end = contig.reference_end                
        contig_coords.append([contig.query_name,ref_start,ref_end])
    contig_coords.sort(key=lambda x: x[1])
    return contig_coords

def get_group_contig_coords_extract(contig_coords):
    contig_coords_with_ref_pos = []
    for i in range(0,len(contig_coords)):
        contig_name = contig_coords[i][0]
        if i == 0:
            contig_ref_start = contig_coords[0][1]
        else:
            contig_ref_start = contig_coords[i-1][2]
        contig_ref_end = contig_coords[i][2]
        contig_coords_with_ref_pos.append([contig_name,contig_ref_start,contig_ref_end])
    return contig_coords_with_ref_pos
   
def get_query_pos(bamfile,contig_name,ref_pos):
    samfile = pysam.AlignmentFile(bamfile)
    query_pos = None
    for contig in samfile:
        if contig.is_secondary:
            continue
        if contig.is_unmapped:
            continue
        if contig.is_supplementary:
            continue
        if contig.query_name != contig_name:
            continue
        for q,r in contig.get_aligned_pairs():
            if r == ref_pos:
                query_pos = q
    assert query_pos != None
    return query_pos

def get_query_seq(bamfile,contig_name,query_start,query_end):
    samfile = pysam.AlignmentFile(bamfile)
    sequence = None
    for contig in samfile:
        if contig.is_secondary:
            continue
        if contig.is_unmapped:
            continue
        if contig.is_supplementary:
            continue
        if contig.query_name != contig_name:
            continue
        sequence = contig.query_sequence[query_start:query_end]
    assert sequence != None
    return sequence

def get_group_sequence(bamfile,group_contigs):
    grouped_sequences = {}
    for group in group_contigs:
        grouped_sequences[group] = []
        contig_coords = get_group_ref_coords(bamfile,group_contigs[group])
        contig_coords_with_ref_pos = get_group_contig_coords_extract(contig_coords)
        for i in range(0,len(contig_coords_with_ref_pos)):
            contig_name = contig_coords_with_ref_pos[i][0]
            ref_start = contig_coords_with_ref_pos[i][1]
            ref_end = contig_coords_with_ref_pos[i][2]
            if i == 0:
                query_start = 0
            else:                
                query_start = get_query_pos(bamfile,contig_name,ref_start)
            if i == len(contig_coords_with_ref_pos):
                query_end = -1
            else:
                query_end = get_query_pos(bamfile,contig_name,ref_end-1)
            sequence = get_query_seq(bamfile,contig_name,query_start,query_end)
            grouped_sequences[group].append(sequence)            
        grouped_sequences[group] = "".join(grouped_sequences[group])
    return grouped_sequences

def get_contig_ref_coords(bamfile):
    coords = {}
    samfile = pysam.AlignmentFile(bamfile)
    for contig in samfile:
        if contig.is_secondary:
            continue
        if contig.is_unmapped:
            continue
        if contig.is_supplementary:
            continue
        coords[contig.query_name] = [samfile.get_reference_name(contig.reference_id),contig.reference_start,contig.reference_end]
    return coords

def return_sequence(bamfile,contained_contigs,group_contigs,group_sequence):
    grouped_contigs = []
    sequence = []
    for group in group_contigs:
        for contig in group_contigs[group]:
            grouped_contigs.append(contig)
    samfile = pysam.AlignmentFile(bamfile)
    for contig in samfile:
        if contig.is_secondary:
            continue
        if contig.is_unmapped:
            continue
        if contig.is_supplementary:
            continue
        if contig.query_name in contained_contigs:
            continue
        if contig.query_name in grouped_contigs:
            continue
        sequence.append([contig.query_name,contig.query_sequence])
    for group in group_sequence:
        sequence.append([group,group_sequence[group]])
    return sequence

bamfile = sys.argv[1]
contig_coords = get_contig_ref_coords(bamfile)
contained_contigs = get_contained_contigs(contig_coords,bamfile)
overlapping_contigs = get_overlapping_contigs(contig_coords,contained_contigs)
group_contigs = group_overlapping_contigs(overlapping_contigs)
group_sequence = get_group_sequence(bamfile,group_contigs)
sequences = return_sequence(bamfile,contained_contigs,group_contigs,group_sequence)
for name,sequence in sequences:
    print ">%s\n%s" % (name,sequence) 
