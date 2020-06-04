#Author: Renato Oliveira
#version: 1.0
#Date: 04-06-2020

###    Copyright (C) 2020  Renato Oliveira
###
###    This program is free software: you can redistribute it and/or modify
###    it under the terms of the GNU General Public License as published by
###    the Free Software Foundation, either version 3 of the License, or
###    any later version.
###
###    This program is distributed in the hope that it will be useful,
###    but WITHOUT ANY WARRANTY; without even the implied warranty of
###    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###    GNU General Public License for more details.
###
###    You should have received a copy of the GNU General Public License
###    along with this program.  If not, see <http://www.gnu.org/licenses/>.

###Contacts:
###    ronnie.alves@itv.org
###    renato.renison@gmail.com or renato.oliveira@pq.itv.org

#!/bin/bash
#usage for Illumina Data: ./qc_docker.sh illumina <forward_reads> <reverse_reads> <Adapters.txt> <min_quality> <min_len> <output_folder> <threads>
#<forward_reads> = path to the forward reads.
#<reverse_reads> = path to the reverse reads.
#<Adapters.txt> = txt file with the pairs of adapters separated by tab, to remove from rawdata.
#<min_quality> = minimum quality PHRED value for the trimming and filtering steps.
#<min_len> = minimum length for the trimming steps. Reads smaller than min_len will be discarded.
#<output_folder> = folder where all the results will be saved.
#<threads> = number of threads to use in the analisys.
#Example = ./qc_docker.sh illumina rawdata/SRR11587600_1.fastq rawdata/SRR11587600_2.fastq adapters.txt 20 50 output_qc 24