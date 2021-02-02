#Author: Renato Oliveira
#version: 1.0
#Date: 27-07-2020

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
#usage for Illumina Data: ./assembly_docker.sh -i illumina -1 <forward_reads> -2 <reverse_reads> -r <ref_fasta> -k <kmer_decontamination> -m <max_mismatch> -o <output_folder> -s <sample_name> -t <threads> -g <max_memory>
#-i = Sequencinf platform. Either "illumina" or "pacbio"
#-1 = path to the forward reads.
#-2 = path to the reverse reads.
#-k = Kmer to use in decontamination when comparing to the reference fasta. Default is 31.
#-m = Maximum of mismatch to accept in the kmers. Default is 2.
#-l = Minimum length for sequence filtering. Default is 100.
#-c = Minimum coverage for sequence filtering. Default is 10.
#-r = Reference file to use in decontamination
#-o = folder where all the results will be saved. Default is "output"
#-s = Sample name. Default is "sample"
#-t = number of threads to use in the analisys. Default is 1.
#-g = Maximum of memory in Gigabytes to use in decontamination step. Default is 80.
#Example = ./assembly_docker.sh -i illumina -1 output_qc/SRR11587600_good.pair1.truncated -2 output_qc/SRR11587600_good.pair2.truncated -r sars-cov-2_MN908947.fasta -k 31 -m 2 -o output_assembly -t 24 -s illumina_rtpcr -g 80
