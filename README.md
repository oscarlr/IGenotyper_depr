# IG_clean

[Introduction](#introduction)  
[Tool requirements](#tool-requirements)  
[Installation](#installation)  
[Creating IGH specific reference](#creating-igh-specific-reference)
[Output directories](#output-directories)

## Introduction
## Installation
## Creating IGH specific reference
IG uses a specific reference. To create this reference, run the command `IG-make-ref`. The input to `IG-make-ref` is the path to the hg19 reference fasta file. `IG-make-ref` will create the reference and index the reference.
```
# Example of running IG-make-ref
IG-make-ref reference/hg19.fasta
```
## Output directories
| File names   | Description   |
|--------------|---------------|
| alignments/  |               |
| assembly/    |               |
| variants/    |               |
| alleles/     |               |
| stats/       |               |
| tmp/         |               |
