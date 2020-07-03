#Author: Renato Oliveira
#version: 1.0
#Date: 13-06-2020

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
#usage for Illumina Data: ./assembly_docker.sh -i illumina -1 <forward_reads> -2 <reverse_reads> -r <ref_fasta> -k <kmer_decontamination> -m <max_mismatch> -o <output_folder> -s <sample_name> -t <threads>
#-i = Sequencinf platform. Either "illumina" or "pacbio"
#-1 = path to the forward reads.
#-2 = path to the reverse reads.
#-k = Kmer to use in decontamination when comparing to the reference fasta. Default is 31.
#-m = Maximum of mismatch to accept in the kmers. Default is 2.
#-r = Reference file to use in decontamination
#-o = folder where all the results will be saved. Default is "output"
#-s = Sample name. Default is "sample"
#-t = number of threads to use in the analisys. Default is 1.
#Example = 

KMER="31"
MAX_MISMATCH="2"
OUTPUT="output"
SAMPLE_NAME="sample"
THREADS="1"

while getopts "i:1:2:k:r:s:o:t:m:" opt; do
	case $opt in
		i) SEQUENCER="$OPTARG"
		;;
		1) FORWARD_READS="$OPTARG"
		;;
		2) REVERSE_READS="$OPTARG"
		;;
		k) KMER="$OPTARG"
		;;
		m) MAX_MISMATCH="$OPTARG"
		;;
		r) REF="$OPTARG"
		;;
		s) SAMPLE_NAME="$OPTARG"
		;;
		o) OUTPUT="$OPTARG"
		;;
		t) THREADS="$OPTARG"
		;;
		\?) echo "Invalid option -$OPTARG" >&2
    	;;
	esac
done

CURRENT_PATH=$(pwd)

if [ ! -d $OUTPUT ]; then
	mkdir $OUTPUT
fi

if [ $SEQUENCER = "illumina" ];
then

	DIR_NAME_FOR=$(dirname $FORWARD_READS)
	cd $DIR_NAME_FOR
	FULL_PATH_FOR=$(pwd)
	cd $CURRENT_PATH

	DIR_NAME_REV=$(dirname $REVERSE_READS)
	cd $DIR_NAME_REV
	FULL_PATH_REV=$(pwd)
	cd $CURRENT_PATH

	DIR_NAME_REF=$(dirname $REF)
	cd $DIR_NAME_REF
	FULL_PATH_REF=$(pwd)
	cd $CURRENT_PATH

	COMMON_PATH=$({ echo $FULL_PATH_FOR; echo $FULL_PATH_REV; echo $FULL_PATH_REF;} | sed -e 'N;s/^\(.*\).*\n\1.*$/\1\n\1/;D')
	FORWARD_READS=$(echo ${FULL_PATH_FOR#"$COMMON_PATH"})/$(basename $FORWARD_READS)
	REVERSE_READS=$(echo ${FULL_PATH_REV#"$COMMON_PATH"})/$(basename $REVERSE_READS)
	REF=$(echo ${FULL_PATH_REF#"$COMMON_PATH"})/$(basename $REF)
	newfile="$(basename $FORWARD_READS)"
	newfile="${newfile%%_*}"
	newfile="${newfile%%.*}"

	#Filtering contaminants 
#	echo "Creating a BBDuk Container: "
#	docker run -id -v $COMMON_PATH:/bbduk/ -v $CURRENT_PATH:/output/ --name bbduk itvds/covid19_bbtools:BBMap.38.76

	matched_fastq1='/output/'$OUTPUT'/decontamination_out/'$SAMPLE_NAME'_matched_r1.fastq'
	matched_fastq2='/output/'$OUTPUT'/decontamination_out/'$SAMPLE_NAME'_matched_r2.fastq'
	unmatched_fastq1='/output/'$OUTPUT'/decontamination_out/'$SAMPLE_NAME'_unmatched_r1.fastq' 
	unmatched_fastq2='/output/'$OUTPUT'/decontamination_out/'$SAMPLE_NAME'_unmatched_r2.fastq' 
	filter_stats='/output/'$OUTPUT'/decontamination_out/'$SAMPLE_NAME'_stats_filtering.txt'
		
#	echo "Running BBDuk Container: "
#	docker exec -i bbduk /bin/bash -c "mkdir /output/'$OUTPUT'; mkdir /output/'$OUTPUT'/decontamination_out; bbduk.sh -Xmx80g in1=/bbduk/'$FORWARD_READS' in2=/bbduk/'$REVERSE_READS' out1=$unmatched_fastq1 out2=$unmatched_fastq2 outm1=$matched_fastq1 outm2=$matched_fastq2 ref=/bbduk/'$REF' k=$KMER hdist=$MAX_MISMATCH stats=$filter_stats overwrite=TRUE t=$THREADS ;chmod -R 777 /output/'$OUTPUT'/decontamination_out"

	#Alignin with Bowtie
	echo "Creating a Bowtie Container: "
	docker run -id -v $COMMON_PATH:/bowtie/ -v $CURRENT_PATH:/output/ --name bowtie itvds/covid19_bowtie2:v2.3.5.1

	echo "Running the Bowtie Container: "
	docker exec -i bowtie  /bin/bash -c "mkdir /output/'$OUTPUT'/bowtie_out; cd /output/'$OUTPUT'/bowtie_out ; \
			bowtie2-build /bowtie/'$REF' /bowtie/'$REF'; \
			bowtie2 -x /bowtie/'$REF'-1 /bowtie/'$FORWARD_READS' -2 /bowtie/'$REVERSE_READS' -p $THREADS > bowtie_output; \
			chmod -R 777 /output/'$OUTPUT'/bowtie_out"

	#Starting with SAMTOOLS
	echo "Creating a SAMTOOLS Container: "
	docker run -id -v $COMMON_PATH:/samtools/ -v $CURRENT_PATH:/output/ --name samtools covid19_samtools

	mappedtoref_bam='/output/'$OUTPUT'/samtools_out/'$SAMPLE_NAME'.bam'
	samtools view -bS - > $mappedtoref_bam
	samtools sort -@ $THREADS -o '/output/'$OUTPUT'/samtools_out/'$SAMPLE_NAME'.sorted.bam' $mappedtoref_bam 
	rm $mappedtoref_bam 
	mv './mapped_reads/'$sampname'.sorted.bam' $mappedtoref_bam 

	echo "Running the SAMTOOLS Container: "
	docker exec -i samtools  /bin/bash -c "mkdir /output/'$OUTPUT'/samtools_out; cd /output/'$OUTPUT'/samtools_out ; \
			samtools view -bS -@ $THREADS /output/'$OUTPUT'/bowtie_out > $mappedtoref_bam; \
			samtools sort -@ $THREADS -o /output/'$OUTPUT'/samtools_out/'$SAMPLE_NAME'.sorted.bam $mappedtoref_bam; \
			rm $mappedtoref_bam; mv /output/'$OUTPUT'/samtools_out/'$SAMPLE_NAME'.sorted.bam $mappedtoref_bam; \
			chmod -R 777 /output/'$OUTPUT'/samtools_out"

	#Checking quality for good quality reads
	#echo "Running FastQC Container on good data: "
	#docker exec -i fastqc /bin/bash -c "fastqc /output_qc/${OUTPUT}/${newfile}_good.pair1.truncated /output_qc/${OUTPUT}/${newfile}_good.pair2.truncated -t $THREADS --nogroup -o /output_qc/${OUTPUT}; chmod -R 777 /output_qc/${OUTPUT}"



fi

echo "Stopping Containeres: "
docker stop bbduk
#docker stop adapter_removal

echo "Removing Containeres: "
docker rm bbduk
#docker rm adapter_removal

echo "Done!"