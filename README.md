# IG_clean

[Introduction](#introduction)  
[Tool requirements](#tool-requirements)  
[Installation](#installation)  
[Creating IGH specific reference](#creating-igh-specific-reference)<br>
[Quick start](#quick-start)<br>
[Explanation of steps](#explanation-of-steps)<br>
[Output directories](#output-directories)<br>
[Notes](#notes)

## Introduction
IGenotyper (or IG) was developed for PacBio capture data to assemble the Immunoglobulin Heavy Chain locus (IGH), genotype the IGH genes, and identify SNPs and SVs within the IGH locus.

## Tool requirements
1. Linux operating system
2. [Conda package](https://conda.io/en/latest/)
3. [cluster python package](https://github.com/oscarlr/cluster)

## Installation
```
### Installing IGenotyper and it's dependencies
git clone https://github.com/oscarlr/IG_clean.git
cd IG_clean
conda env create -f environment.yml 
conda activate IG_clean
python setup.py install

### Installing cluster package that's needed
cd ..
git clone https://github.com/oscarlr/cluster.git
cd cluster
python setup.py install
```

## Creating IGH specific reference
IG uses a specific reference. To create this reference, run the command `IG-make-ref`. The input to `IG-make-ref` is the path to the hg19 reference fasta file. `IG-make-ref` will create the reference and index the reference.
```
# Example of running IG-make-ref
IG-make-ref reference/hg19.fasta
```
## Quick start
```
IG --phase <pacbio bam file> <output> 
IG --assemble <pacbio bam file> <output> 
IG --detect <pacbio bam file> <output> 
IG --stats <pacbio bam file> <output> 
```
## Explanation of steps
### Phase
In the first step `--phase`, the subreads and CCS reads are phased and aligned to the IGH specific reference. Each read has a read group annotation. A read group annotation of 1 and 2 corresponds to haplotype 1 and 2. The read group annotation of 0 corresponds to unassignable reads. In IGV, you can seperate these reads by left clicking and selecting group by read group.

### Assemble
In the second step `assembly`, the haplotypes are assembled. During this process folders will be created for each region/haplotype block. Within each folder there is a bash script that runs the assembly process. These can be submitted as a single job into the cluster (this speeds up the process).

### Detect
### Stats
## Output directories
| Directories            | Description                                          |
|------------------------|------------------------------------------------------|
| `<output>/alignments`  | Alignments of CCS and subreads (phased and unphased) |
| `<output>/assembly`    | Assembly of IGH locus                                |
| `<output>/variants`    | SNVs, indels and SVs                                 |
| `<output>/alleles`     | Alleles in sample                                    |
| `<output>/stats`       | Statistics from different IG steps                   |
| `<output>/tmp`         | Temporary files. Could be deleted.                   |

## Notes
1. The most computationaly expensive step is turning the PacBio subreads to CCS reads. This usually takes 1 day with 10+ cores. Of course, this depends on the coverage of the sample.
