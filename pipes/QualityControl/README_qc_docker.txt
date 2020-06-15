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
#usage for Illumina Data: ./qc_docker.sh -i illumina -1 <forward_reads> -2 <reverse_reads> -a <Adapters.txt> -q <min_quality> -l <min_len> -o <output_folder> -t <threads>
#-i = Sequencing platform. Either "illumina" or "pacbio"
#-1 = path to the forward reads.
#-2 = path to the reverse reads.
#-a = txt file with the pairs of adapters separated by tab, to remove from rawdata.
#-q = minimum quality PHRED value for the trimming and filtering steps. Default is 20.
#-l = minimum length for the trimming steps. Reads smaller than min_len will be discarded. Default is 50.
#-o = folder where all the results will be saved. Default is "output"
#-t = number of threads to use in the analisys.
#Example = ./qc_docker.sh -i illumina -1 rawdata/SRR11587600_1.fastq -2 rawdata/SRR11587600_2.fastq -a rawdata/adapters.txt -q 20 -l 50 -o output_qc -t 24
