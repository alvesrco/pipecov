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
SEQUENCER=$1
ADAPTERS=$4
MIN_QUAL=$5
MIN_LEN=$6
OUTPUT=$7
THREADS=$8
CURRENT_PATH=$(pwd)

if [ ! -d $OUTPUT ]; then
	mkdir $OUTPUT
fi

if [ $SEQUENCER = "illumina" ];
then
	FORWARD_READS=$2
	REVERSE_READS=$3

	DIR_NAME_FOR=$(dirname $FORWARD_READS)
	cd $DIR_NAME_FOR
	FULL_PATH_FOR=$(pwd)
	cd $CURRENT_PATH

	DIR_NAME_REV=$(dirname $REVERSE_READS)
	cd $DIR_NAME_REV
	FULL_PATH_REV=$(pwd)
	cd $CURRENT_PATH

	DIR_NAME_ADAP=$(dirname $ADAPTERS)
	cd $DIR_NAME_ADAP
	FULL_PATH_ADAP=$(pwd)
	cd $CURRENT_PATH

	COMMON_PATH=$({ echo $FULL_PATH_FOR; echo $FULL_PATH_REV; echo $FULL_PATH_ADAP;} | sed -e 'N;s/^\(.*\).*\n\1.*$/\1\n\1/;D')
	newfile="$(basename $FORWARD_READS)"
	newfile="${newfile%%_*}"
	newfile="${newfile%%.*}"

	#Checking rawdata quality with FastQC
	echo "Creating a FastQC Container: "
	docker run -id -v $COMMON_PATH:/fastqc/ -v $CURRENT_PATH:/output_qc/ --name fastqc itvds/covid19_fastqc:v0.11.9
	
	echo "Running FastQC Container on rawdata: "
	docker exec -i fastqc /bin/bash -c "fastqc /fastqc/$(echo ${FULL_PATH_FOR#"$COMMON_PATH"})/$(basename $FORWARD_READS) /fastqc/$(echo ${FULL_PATH_REV#"$COMMON_PATH"})/$(basename $REVERSE_READS) -t $THREADS --nogroup -o /output_qc/${OUTPUT}; chmod -R 777 /output_qc/${OUTPUT}"

	#Removing adapters and filtering sequences by quality
	echo "Creating an AdapterRemoval Container: "
	docker run -id -v $COMMON_PATH:/adapter_removal/ -v $CURRENT_PATH:/output/ --name adapter_removal itvds/covid19_adapterremoval:v.2.2.2

	echo "Running the AdapterRemoval Container: "
	docker exec -i adapter_removal  /bin/bash -c "cd /output/${OUTPUT}; AdapterRemoval --file1 /adapter_removal/$(echo ${FULL_PATH_FOR#"$COMMON_PATH"})/$(basename $FORWARD_READS) --file2 /adapter_removal/$(echo ${FULL_PATH_REV#"$COMMON_PATH"})/$(basename $REVERSE_READS) --threads $THREADS --mate-separator ' ' --adapter-list /adapter_removal/$(echo ${FULL_PATH_ADAP#"$COMMON_PATH"})/$(basename $ADAPTERS) --trimwindows 10 --minquality $MIN_QUAL --minlength $MIN_LEN --qualitymax 64 --basename ${newfile}_good --mm 5; chmod 777 *"

	#Checking quality for good quality reads
	echo "Running FastQC Container on good data: "
	docker exec -i fastqc /bin/bash -c "fastqc /output_qc/${OUTPUT}/${newfile}_good.pair1.truncated /output_qc/${OUTPUT}/${newfile}_good.pair2.truncated -t $THREADS --nogroup -o /output_qc/${OUTPUT}; chmod -R 777 /output_qc/${OUTPUT}"



fi

echo "Stopping Containeres: "
docker stop fastqc
docker stop adapter_removal

echo "Removing Containeres: "
docker rm fastqc
docker rm adapter_removal

echo "Done!"