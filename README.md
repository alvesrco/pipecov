# PiPeCOV
This project proposes to contribute to fill the knowledge gap about Covid-19SARS-CoV-2 in Brazil and with the global knowledge about the pathogen. 

The scientific community and governments of other countries are looking for actions along these lines. The United Kingdom recently established a research network for genomic studies of SARS-CoV-2Covid-19 with a contribution of GBP 20 million (https://www.gov.uk/government/news/uk-launches-whole- genome-sequence-alliance-to-map-spread-of-coronavirus).

The PiPeCOV pipeline can handle quality assesment, assembly and annotation of SARS-CoV-2 genomes sequenced by Illumina (Amplicon & mNGS). Once running PiPeCOV you obtain the assembled and annotated SARS-CoV-2 genomes at the end of the procedure.

PipeCoV is free to use for non-commercial users, under a GPLv3 License.

The PiPeCOV workflow:
![Screenshot](pipecov.png)

All the steps and commands used in the pipeline for Quality Assesment and Mapping, Assembly, Annotation, and Phylogenetic Assignment of lineages of the SARS-CoV-2 Genomes are encapsulated in images and Dockers containers. The user just needs to have Docker installed on his machine, without worrying about installing all the tools used in the pipelines.

The docker images used in the pipeline can be found at (https://hub.docker.com/u/itvds)

PiPeCOV must be downloaded from this repo (https://github.com/alvesrco/covid19_itvds)

All Dockerfiles, and pipes repo are developped by the Covid19 Project Network @ ITVDS.

## Sample file

Samples in .fastq.gz format. They should be in the pattern: EC114_S15_R1_001.fastq.gz

## How To

**: Quality Assesment :**
```
$ ./qc_docker.sh -i illumina -1 SAMPLE_R1.fastq -2 SAMPLE_R2.fastq -a adapters.txt -q 20 -l 50 -o output_qc -t 24
```
> Parameters
- i illumina		[Sequencing platform]
- q 20		[Minimum PHRED quality for trimming and filtering. Default: 20]
- l 50		[Minimum size of post-trimming sequences. Default: 50]
- t 24		[Number of threads to be used. Default: 1]
> Input
- 1 SAMPLE_R1.fastq	[Forward strings in the original raw format]
- 2 SAMPLE_R2.fastq	[Reverse strings in the original raw format]
- a adapters.txt		[File with sequence adapters that must be removed]
> Output
- o output_qc		[Folder where the results will be saved. Default: “output”]

**: Genome Assembly, Annotation and Phylogenetic Assignment :**
```
$ ./assembly_docker.sh -i illumina -1 output_qc/SRR11587600_good.pair1.truncated -2 output_qc/SRR11587600_good.pair2.truncated -r sars-cov-2_MN908947.fasta -k 31 -m 2 -l 100 -c 10 -o output_assembly -t 24 -s illumina_rtpcr
```
> Parameters
- i illumina		[Sequencing platform]
- k 31		[Size of the kmer in the decontamination step. Default: 31]
- m 2		[Maximum mismatch to be accepted in kmers. Default: 2]
- l 100		[Minimum contig size. Default: 100]
- c 10		[Minimum contig coverage. Default: 10]
- m 2		[Maximum mismatch to be accepted in kmers. Default: 2]
- s SAMPLE_NAME		[Sample name. Default: “sample”]
- t 24		[Number of threads to be used. Default: 1]
- g 80		[Maximum of memory in Gigabytes to use in decontamination step. Defaul: 80]
> Input
- 1 SAMPLE_good_R1.fastq	[Forward sequences after quality treatment]
- 2 SAMPLE_good_R2.fastq	[Reverse sequences after quality treatment]
- r reference.fasta		[Fasta file with the reference (s) to be used]
> Output
- o output_assembly		[Folder where the results will be saved. Default: “output”]
