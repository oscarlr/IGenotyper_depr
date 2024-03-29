# IGenotyper

[Introduction](#introduction)  
[Requirements](#requirements)  
[Installation](#installation)  
[Creating IGH specific reference](#creating-igh-specific-reference)<br>
[Testing IGenotyper installation](#testing-igenotyper-installation)<br>
[Running IGenotyper](#running-igenotyper)<br>
[Explanation of steps](#explanation-of-steps)<br>
[Output directories](#output-directories)<br>
[Output files](#output-files)<br>
[Notes](#notes)

## Introduction
IGenotyper (or IG) was developed for PacBio capture data to assemble the Immunoglobulin Heavy Chain locus (IGH), genotype the IGH genes, and identify SNPs, indels and SVs within the IGH locus.

## Requirements
### Tool requirements
1. Linux operating system
    1. Built on CentOS.
2. [Conda package](https://conda.io/en/latest/)
    1. Built using python2.7 verion. [Download link](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh) 
3. [cluster python package](https://github.com/oscarlr/cluster)
4. gcc v4.8+

## Installation
### Install whatshap conda environment
```
conda create -n whatshap-latest python=3.6
conda activate whatshap-latest
pip install git+https://bitbucket.org/whatshap/whatshap
conda deactivate
```

### Installing IGenotyper
```
git clone https://github.com/oscarlr/IGenotyper.git
cd IGenotyper
conda env create -f environment.yml 
conda activate IG
python setup.py install
```

### Installing cluster package
```
cd ..
git clone https://github.com/oscarlr/cluster.git
cd cluster
python setup.py install
# To run IGenotyper locally run this command
export SJOB_DEFALLOC=NONE
# Or to run IGenotyper on a cluster set SJOB_DEFALLOC to the allocation account
```

## Creating IGH specific reference or download IGH specific reference
IG uses a specific reference fasta file. This reference contains the hg38 (chr1-22, X and Y) chromosomes with hg38/chr14 IGH locus removed and an expanded modified IGH lcous, and two ~20KB insertions TRB sequences as alternate sequences. The reference can be obtained as passed to IG using the following commands:

```
wget https://watsonlab.blob.core.windows.net/public/reference/hg38_igh_trab/reference.fasta
wget https://watsonlab.blob.core.windows.net/public/reference/hg38_igh_trab/reference.fasta.fai

IG-make-ref reference.fasta
```
### Using AzCopy for super fast downloads
```
# https://docs.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10

azcopy copy https://watsonlab.blob.core.windows.net/public/reference/hg38_igh_trab/reference.fasta .
azcopy copy https://watsonlab.blob.core.windows.net/public/reference/hg38_igh_trab/reference.fasta.fai .

IG-make-ref reference.fasta
```

## Testing IGenotyper installation
```
## HAVE TO REDO SINCE SERVERS WERE CHANGED
```

## Running IGenotyper
```
IG-make-ref hg19_no_alts.fa
IG --phase <pacbio bam file> <output> 
IG --assemble <pacbio bam file> <output> 
IG --detect <pacbio bam file> <output> 
```

## Usage
```
(IG)[oscarlr]$ IG --help
usage: IG [-h] [--phase] [--assemble] [--extend_assembly] [--detect] [--stats]
          [--report] [--tmp_dir] [--threads] [--cluster] [--cluster_queue]
          [--cluster_walltime] [--cluster_mem] [--phased_vcf_file]
          [--pacbio_data_type] [--phased_vcf_file_sample_name] [--dont_split]
          [--secondary_read_score] [--keep]
          input_bam outdir

Process IGH capture data

positional arguments:
  input_bam             pacbio subreads bam file
  outdir                output directory

optional arguments:
  -h, --help            show this help message and exit
  --phase               Map and phase reads
  --assemble            Assemble reads
  --extend_assembly     Extend assemblies. Still testing. Dont use.
  --detect              Detect variants
  --stats               Generate stats
  --report              Generate report
  --tmp_dir             temporary directory
  --threads             Number of threads for everything
  --cluster             Use cluster
  --cluster_queue       Queue for cluster
  --cluster_walltime    Walltime for cluster
  --cluster_mem         memory for cluster
  --phased_vcf_file     Run IG with this phased VCF file
  --pacbio_data_type    Pacbio data type (RS or Sequel), either "RS" or "SQ"
  --phased_vcf_file_sample_name 
                        Sample name in phased VCF file
  --dont_split          Do not split assembly regions into SV/non-SV regions
  --secondary_read_score 
                        Min secondary read score to move
  --keep                Keep intermediate files
```

## Explanation of steps
### Phase
In the first step `--phase`, the subreads and CCS reads are phased and aligned to the IGH specific reference. Each read has a read group annotation. A read group annotation of 1 and 2 corresponds to haplotype 1 and 2. The read group annotation of 0 corresponds to unassignable reads. In IGV, you can seperate these reads by left clicking and selecting group by read group.

### Assemble
In the second step `--assembly`, the haplotypes are assembled. During this process folders will be created for each region/haplotype block. Within each folder there is a bash script that runs the assembly process. These can be submitted as a single job into the cluster (this speeds up the process).

### Detect
In the third step `--detect`, SNVs, indels, SVs and gene/alleles are genotyped. A VCF file is created for the SNVs, a BED file for the indels and SVs, a TAB-delimited file for the gene/alleles calls.  

## Output directories
alignments  alleles  assembly  logs  plots  preprocessed  report.html  tmp  variants
| Directories            | Description                                          |
|------------------------|------------------------------------------------------|
| `<output>/alignments`  | Alignments of CCS, subreads and contigs (phased and unphased) |
| `<output>/assembly`    | Assembly of IGH locus                                |
| `<output>/variants`    | SNVs, indels and SVs                                 |
| `<output>/alleles`     | Alleles in sample                                    |
| `<output>/logs`        | Log files with input parameters                   |
| `<output>/tmp`         | Temporary files. Could be deleted.                   |

## Output files
1. `alignments/`
    1. `ccs_to_ref*`: CCS reads aligned to reference
    2. `contigs_to_ref*`: All assembled contigs aligned to refernece
    3. `igh_contigs_to_ref*`: IGH assembled contigs aligned to igh reference
2. `assembly/`
    1. `contigs.fasta`: All assembled contigs
    2. `igh_contigs.fasta`: IGH assembled contigs
3. `alleles/`
    1. `assembly_alleles.bed`: Alleles extracted from the assembly for each gene
    2. `assembly_genes.fasta`: Fasta sequence from the assembly for each gene/allele
    3. `ccs_alleles.bed`: Alleles extracted from the CCS reads for each gene
    4. `ccs_genes.fasta`: Fasta sequence from the CCS reads for each gene/allele
4. `logs/`
    1. `gene_cov.txt`: Haplotype coverage for each gene
5. `variants/`
    1. `snvs_phased_from_ccs.vcf`: Phased SNVs detected from the CCS reads
    2. `snvs_assembly.vcf`: SNVS detected from the assembly
    3. `indel_assembly.bed`: Indels detected from the assembly
    4. `sv_assembly.bed`: SVs detected from the assembly
    5. `phased_blocks.txt`: Phased haplotype blocks
