# IGenotyper

[Introduction](#introduction)  
[Requirements](#requirements)  
[Installation](#installation)  
[Creating IGH specific reference](#creating-igh-specific-reference)<br>
[Testing IGenotyper installation](#testing-igenotyper-installation)<br>
[Running IGenotyper](#running-igenotyper)<br>
[Explanation of steps](#explanation-of-steps)<br>
[Output directories](#output-directories)<br>
[Notes](#notes)

## Introduction
IGenotyper (or IG) was developed for PacBio capture data to assemble the Immunoglobulin Heavy Chain locus (IGH), genotype the IGH genes, and identify SNPs and SVs within the IGH locus.

## Requirements
### Tool requirements
1. Linux operating system
    1. Built on CentOS.
2. [Conda package](https://conda.io/en/latest/)
    1. Built using python2.7 verion. [Download link](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh) 
3. [cluster python package](https://github.com/oscarlr/cluster)

### CPU requirements
1. At least 14 GBs for the `phase` step

## Installation
### Install whatshap conda environment
```
conda create -n whatshap-latest python=3.6
conda activate whatshap-latest
pip install git+https://bitbucket.org/whatshap/whatshap
conda deactivate
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

### Installing IGenotyper
```
git clone https://github.com/oscarlr/IGenotyper.git
cd IGenotyper
conda env create -f environment.yml 
conda activate IG
python setup.py install
```

## Creating IGH specific reference or download IGH specific reference
IG uses a specific reference. To create this reference, run the command `IG-make-ref`. The input to `IG-make-ref` is the path to the hg19 reference fasta file.`IG-make-ref` will create the reference and index the reference. The reference fasta file should only contain chr1-22, X and Y. No alternate sequences should be in the reference fasta file.

### Create IGH specific reference
You can create the IGH specific reference using the commands below or download the reference as shown in the following section.
```
### Downloading hg19 reference and creating hg19 reference without alternate sequences
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
samtools faidx hg19.fa
for i in {1..22} X Y
do
    samtools faidx hg19.fa chr${i}
done > hg19_no_alts.fa
samtools faidx hg19_no_alts.fa

# Example of running IG-make-ref
IG-make-ref hg19_no_alts.fa
```
### Download IGH reference
```
wget https://users.hpc.mssm.edu/~rodrio10/public/IGenotyper/ref/reference.fasta
wget https://users.hpc.mssm.edu/~rodrio10/public/IGenotyper/ref/reference.fasta.sa 
wget https://users.hpc.mssm.edu/~rodrio10/public/IGenotyper/ref/reference.mmi

IG-make-ref reference.fasta --sa reference.fasta.sa --mmi reference.mmi
```

## Testing IGenotyper installation
```
wget https://users.hpc.mssm.edu/~rodrio10/public/IGenotyper/test_data/NA19240_subreads.bam
wget hhttps://users.hpc.mssm.edu/~rodrio10/public/IGenotyper/test_data/NA19240_subreads.bam.bpi

IG --phase NA19240_subreads.bam test_dir
```

## Running IGenotyper
```
IG --phase <pacbio bam file> <output> 
IG --assemble <pacbio bam file> <output> 
IG --detect <pacbio bam file> <output> 
IG --stats <pacbio bam file> <output> 
```

## Usage
```
(IG)[oscarlr]$ IG --help
usage: IG [-h] [--phase] [--assemble] [--extend_assembly] [--detect] [--stats]
          [--report] [--tmp_dir tmp_dir] [--threads threads] [--cluster]
          [--cluster_queue] [--cluster_threads] [--cluster_walltime]
          [--cluster_mem] [--haploid] [--phased_vcf_file PHASED_VCF_FILE]
          [--pacbio_data_type PACBIO_DATA_TYPE]
          [--phased_vcf_file_sample_name PHASED_VCF_FILE_SAMPLE_NAME]
          [--add_unphased_reads] [--dont_split]
          [--secondary_read_score SECONDARY_READ_SCORE]
          input_bam outdir

Process IGH capture data

positional arguments:
  input_bam             bam file containing raw reads
  outdir                output directory

optional arguments:
  -h, --help            show this help message and exit
  --phase               Map and phase reads
  --assemble            Only assemble reads
  --extend_assembly     Extend assemblies
  --detect              Detect variants
  --stats               Generate stats
  --report              Generate report
  --tmp_dir tmp_dir     temporary directory
  --threads threads     Number of threads for everything
  --cluster             Use cluster
  --cluster_queue       Queue for cluster
  --cluster_threads     Number of threads for cluster jobs
  --cluster_walltime    Walltime for cluster
  --cluster_mem         memory for cluster
  --haploid             Run Quiver in haploid mode
  --phased_vcf_file PHASED_VCF_FILE
                        Run IG with this phased VCF file
  --pacbio_data_type PACBIO_DATA_TYPE
                        Pacbio data type (RS or Sequel), either "RS" or "SQ"
  --phased_vcf_file_sample_name PHASED_VCF_FILE_SAMPLE_NAME
                        Sample name in phased VCF file
  --add_unphased_reads  Add unphased reads to phased region
  --dont_split          Do not split assembly regions into SV/non-SV regions
  --secondary_read_score SECONDARY_READ_SCORE
                        Min secondary read score to move
```

## Explanation of steps
### Phase
In the first step `--phase`, the subreads and CCS reads are phased and aligned to the IGH specific reference. Each read has a read group annotation. A read group annotation of 1 and 2 corresponds to haplotype 1 and 2. The read group annotation of 0 corresponds to unassignable reads. In IGV, you can seperate these reads by left clicking and selecting group by read group.

### Assemble
In the second step `--assembly`, the haplotypes are assembled. During this process folders will be created for each region/haplotype block. Within each folder there is a bash script that runs the assembly process. These can be submitted as a single job into the cluster (this speeds up the process).

### Detect
In the third step `--detect`, SNVs, indels, SVs and gene/alleles are genotyped. A VCF file is created for the SNVs, a BED file for the indels and SVs, a TAB-delimited file for the gene/alleles calls.  

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
