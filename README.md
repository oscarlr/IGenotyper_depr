# IG_clean

[Introduction](#introduction)  
[Tool requirements](#tool-requirements)  
[Installation](#installation)  
[Creating IGH specific reference](#creating-igh-specific-reference)<br>
[Quick start](#quick-start)<br>
[Explanation of steps](#explanation-of-steps)<br>
[Output directories](#output-directories)

## Introduction
## Installation
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
```
## Explanation of steps

## Output directories
| Directories            | Description                                          |
|------------------------|------------------------------------------------------|
| `<output>/alignments`  | Alignments of CCS and subreads (phased and unphased) |
| `<output>/assembly`    | Assembly of IGH locus                                |
| `<output>/variants`    | SNVs, indels and SVs                                 |
| `<output>/alleles`     | Alleles in sample                                    |
| `<output>/stats`       | Statistics from different IG steps                   |
| `<output>/tmp`         | Temporary files. Could be deleted.                   |
