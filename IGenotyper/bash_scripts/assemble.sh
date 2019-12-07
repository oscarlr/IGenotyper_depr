#!/bin/bash
set -x #-e

### IGNORE
# hap=$1
# ccs_reads_to_ref=$2
# chrom=$3
# start=$4
# end=$5
# output=$6
# threads=$7
# size=$8
# subreads=${9}
# raw_subreads_to_ref=${10}
### IGNORE

if [ -s ${output}/contig_after_filter_to_ref.bam ]
then
    if [ -s ${output}/contig_after_filter.fastq ]
    then
	if [ -s ${output}/contig_after_filter.fasta ]
	then
	    if [ -s ${output}/merged_contigs.fasta ]
	    then
		echo "" > ${output}/done
		exit 0
	    fi
	fi
    fi
fi

### 1. Assemble with higher coverage
if [ ! -s ${output}/canu/canu.contigs.fasta ]
then    
    samtools view ${ccs_to_ref} ${hap} ${chrom}:${start}-${end} | awk '{ print ">"$1"\n"$10}' > ${output}/reads.fasta
    canu \
	-p canu \
	-d ${output}/canu \
	corOutCoverage=200 \
	minThreads=${threads} \
	genomeSize=${size} \
	useGrid=0 \
	contigFilter="2 0 1.0 0.5 0" \
	-pacbio-raw ${output}/reads.fasta
fi


### 2. If no contigs are assembled then assemble with lower coverage
if [ ! -s ${output}/canu/canu.contigs.fasta ]
then    
    rm -fr ${output}/canu
    samtools view ${subreads_to_ref} ${hap} ${chrom}:${start}-${end} | awk '{ print ">"$1"\n"$10}' > ${output}/reads.fasta
    canu \
	-p canu \
	-d ${output}/canu \
	minThreads=${threads} \
	genomeSize=${size} \
	useGrid=0 \
	contigFilter="2 0 1.0 0.5 0" \
	-pacbio-raw ${output}/reads.fasta
fi

### 3. Clean contigs if contigs are assembled
if [ -s ${output}/canu/canu.contigs.fasta ]
then
    if [ ! -s ${output}/subreads.bam.pbi ]
    then
	samtools view ${subreads_to_ref} ${hap} ${chrom}:${start}-${end} | awk '{ print $1 }' | sort | uniq > ${output}/reads.names
	python ${python_scripts}/extract_reads.py \
	    -b ${subreads} -n ${output}/reads.names -o ${output}/subreads.bam
	pbindex ${output}/subreads.bam
    fi
    if [ ! -s ${output}/canu/reads_to_canu_contigs.sorted.bam.pbi ]
    then
	pbmm2 align \
	    --sort \
	    -J ${threads} \
	    -j ${threads} \
	    ${output}/canu/canu.contigs.fasta \
	    ${output}/subreads.bam \
	    ${output}/canu/reads_to_canu_contigs.sorted.bam
	pbindex ${output}/canu/reads_to_canu_contigs.sorted.bam
    fi
    if [ ! -s ${output}/contig.fastq ]
    then
	samtools faidx ${output}/canu/canu.contigs.fasta
	variantCaller \
	    --referenceFilename ${output}/canu/canu.contigs.fasta \
	    -j ${threads} \
	    -o ${output}/contig.fastq \
	    -o ${output}/contig.fasta \
	    ${output}/canu/reads_to_canu_contigs.sorted.bam
    fi
fi

### 4. Filter contigs based on quality scores and overlapping contigs
if [ ! -s ${output}/contig_after_filter.fastq ]
then
    python ${python_scripts}/average_quality_per_contig.py \
	${output}/contig.fastq > \
	${output}/contig.qual_values
    
    samtools faidx ${ref} ${chrom}:${start}-${end} > \
	${output}/ref.fa
    
    pbmm2 align \
	--sort \
	-J ${threads} \
	-j ${threads} \
	${output}/ref.fa \
	${output}/contig.fasta \
	${output}/contig_to_ref.bam \
	--preset "CCS"
    
    python ${python_scripts}/mapping_coords.py ${output}/contig_to_ref.bam > \
	${output}/contig_to_ref.bed
    
    bedtools intersect -wb -f .9 \
	-a ${output}/contig_to_ref.bed \
	-b ${output}/contig_to_ref.bed | awk '$4 != $10' > \
	${output}/contig_overlap.bed
    
    cat ${output}/contig.qual_values | awk '$2 < 20' | cut -f1 > ${output}/contig_to_filter.txt
    
    cat ${output}/contig_overlap.bed | cut -f4 | cut -f1 -d/ >> ${output}/contig_to_filter.txt
    
    samtools faidx ${output}/contig.fasta
    cat ${output}/contig.fasta.fai | cut -f1 | grep -v -f ${output}/contig_to_filter.txt | while read contig_name
    do
	samtools faidx ${output}/contig.fasta ${contig_name}
    done > ${output}/contig_after_filter.fasta
fi

### 5. Filter fastq file as well 
if [ ! -s ${output}/contig_after_filter.fastq ]
then
    cat ${output}/contig.fasta.fai | cut -f1 | grep -v -f ${output}/contig_to_filter.txt | while read contig_name
    do
	samtools faidx ${output}/contig.fastq ${contig_name} --fastq --length 1000000000
    done > ${output}/contig_after_filter.fastq
fi

### 6. If all the contigs were filtered out then re-assemble with lower coverage
if [ ! -s ${output}/contig_after_filter.fastq ]
then
    rm -fr ${output}/c*
    
    if [ ! -s ${output}/canu/canu.contigs.fasta ]
    then    
	samtools view ${subreads_to_ref} ${hap} ${chrom}:${start}-${end} | awk '{ print ">"$1"\n"$10}' > ${output}/reads.fasta
	canu \
	    -p canu \
	    -d ${output}/canu \
	    minThreads=${threads} \
	    genomeSize=${size} \
	    useGrid=0 \
	    contigFilter="2 0 1.0 0.5 0" \
	    -pacbio-raw ${output}/reads.fasta
    fi
    
    if [ -s ${output}/canu/canu.contigs.fasta ]
    then
	if [ ! -s ${output}/subreads.bam.pbi ]
	then
	    samtools view ${subreads_to_ref} ${hap} ${chrom}:${start}-${end} | awk '{ print $1 }' | sort | uniq > ${output}/reads.names
	    python ${python_scripts}/extract_reads.py \
		-b ${subreads} -n ${output}/reads.names -o ${output}/subreads.bam
	    pbindex ${output}/subreads.bam
	fi
	if [ ! -s ${output}/canu/reads_to_canu_contigs.sorted.bam.pbi ]
	then
	    pbmm2 align \
		--sort \
		-J ${threads} \
		-j ${threads} \
		${output}/canu/canu.contigs.fasta \
		${output}/subreads.bam \
		${output}/canu/reads_to_canu_contigs.sorted.bam
	    pbindex ${output}/canu/reads_to_canu_contigs.sorted.bam
	fi
	if [ ! -s ${output}/contig.fastq ]
	then
	    samtools faidx ${output}/canu/canu.contigs.fasta
	    variantCaller \
		--referenceFilename ${output}/canu/canu.contigs.fasta \
		-j ${threads} \
		-o ${output}/contig.fastq \
		-o ${output}/contig.fasta \
		${output}/canu/reads_to_canu_contigs.sorted.bam
	fi
    fi
    
    if [ ! -s ${output}/contig_after_filter.fastq ]
    then
	python ${python_scripts}/extract_reads.py \
	    ${output}/contig.fastq > \
	    ${output}/contig.qual_values
	
	samtools faidx ${ref} ${chrom}:${start}-${end} > \
	    ${output}/ref.fa
	
	pbmm2 align \
	    --sort \
	    -J ${threads} \
	    -j ${threads} \
	    ${output}/ref.fa \
	    ${output}/contig.fasta \
	    ${output}/contig_to_ref.bam \
	    --preset "CCS"
	
	python ${python_scripts}/mapping_coords.py ${output}/contig_to_ref.bam > \
	    ${output}/contig_to_ref.bed
	
	bedtools intersect -wb -f .9 \
	    -a ${output}/contig_to_ref.bed \
	    -b ${output}/contig_to_ref.bed | awk '$4 != $10' > \
	    ${output}/contig_overlap.bed
	
	cat ${output}/contig.qual_values | awk '$2 < 20' | cut -f1 > ${output}/contig_to_filter.txt
	
	cat ${output}/contig_overlap.bed | cut -f4 | cut -f1 -d/ >> ${output}/contig_to_filter.txt	
	
	samtools faidx ${output}/contig.fasta
	cat ${output}/contig.fasta.fai | cut -f1 | grep -v -f ${output}/contig_to_filter.txt | while read contig_name
	do
	    samtools faidx ${output}/contig.fasta ${contig_name}
	done > ${output}/contig_after_filter.fasta
    fi
    
    cat ${output}/contig.fasta.fai | cut -f1 | grep -v -f ${output}/contig_to_filter.txt | while read contig_name
    do
	samtools faidx ${output}/contig.fastq ${contig_name} --fastq --length 1000000000
    done > ${output}/contig_after_filter.fastq    
fi

samtools faidx ${ref} ${chrom}:${start}-${end} > \
    ${output}/ref.fa

### 7. Align to the reference
pbmm2 align \
    --sort \
    -J ${threads} \
    -j ${threads} \
    ${output}/ref.fa \
    ${output}/contig_after_filter.fastq \
    ${output}/contig_after_filter_to_ref_before_merge.bam \
    --preset "CCS"

python ${python_scripts}/merge_contigs.py \
    ${output}/contig_after_filter_to_ref_before_merge.bam > ${output}/merged_contigs.fasta

samtools faidx ${output}/merged_contigs.fasta

pbmm2 align \
    --sort \
    -J ${threads} \
    -j ${threads} \
    ${output}/merged_contigs.fasta \
    ${output}/subreads.bam \
    ${output}/reads_to_merged_contigs.sorted.bam

pbindex ${output}/reads_to_merged_contigs.sorted.bam

samtools faidx ${output}/merged_contigs.fasta

variantCaller \
    --referenceFilename ${output}/merged_contigs.fasta  \
    -j ${threads} \
    -o ${output}/merged_contigs_quivered.fastq \
    -o ${output}/merged_contigs_quivered.fasta \
    ${output}/reads_to_merged_contigs.sorted.bam

pbmm2 align \
    --sort \
    -J ${threads} \
    -j ${threads} \
    ${output}/ref.fa \
    ${output}/merged_contigs_quivered.fastq \
    ${output}/contig_after_filter_to_ref.bam \
    --preset "CCS"

echo "" > ${output}/done
