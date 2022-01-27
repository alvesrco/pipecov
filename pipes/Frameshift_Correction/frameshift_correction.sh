#Author: Renato Oliveira
#version: 1.0
#Date: 24-01-2022

###    Copyright (C) 2022  Renato Oliveira
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
#usage : ./frameshift_correction.sh -i <input_fasta> -r <ref_fasta> -t <threads> -o <output>
#-i = Fasta file with genome sequence to be corrected.
#-r = Reference file to aid in the frameshift correction.
#-t = number of threads to use in the analisys. Default is 1.
#-o = folder where all the results will be saved. Default is "output"
#Example = ./frameshift_correction.sh -i sample_genome.fasta -r sars-cov-2_MN908947.fasta -t 24 -o output

OUTPUT="output"
THREADS="1"
TIMESTAMP=$(date +"%Y%m%d%H%M%S")

while getopts "i:r:o:t:" opt; do
	case $opt in
		i) INPUT="$OPTARG"
		;;
		r) REF="$OPTARG"
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


DIR_NAME_INPUT=$(dirname $INPUT)
cd $DIR_NAME_INPUT
FULL_PATH_INPUT=$(pwd)
cd $CURRENT_PATH

DIR_NAME_REF=$(dirname $REF)
cd $DIR_NAME_REF
FULL_PATH_REF=$(pwd)
cd $CURRENT_PATH

COMMON_PATH=$({ echo $FULL_PATH_INPUT; echo $FULL_PATH_REF;} | sed -e 'N;s/^\(.*\).*\n\1.*$/\1\n\1/;D')
INPUT=$(echo ${FULL_PATH_INPUT#"$COMMON_PATH"})/$(basename $INPUT)
REF=$(echo ${FULL_PATH_REF#"$COMMON_PATH"})/$(basename $REF)
newfile="$(basename $INPUT)"
newfile="${newfile%%_*}"
newfile="${newfile%%.*}"

# echo "Current path: '$CURRENT_PATH'"
# echo "Common path: '$COMMON_PATH'"


#Aligning genomes
echo "Creating a MAFFT Container: "
docker run -id -v $COMMON_PATH:/common/ -v $CURRENT_PATH:/output/ --name mafft_$TIMESTAMP staphb/mafft:latest
	
echo "Running MAFFT Container: "
docker exec -i mafft_$TIMESTAMP /bin/bash -c "mkdir -p /output/'$OUTPUT'/mafft_results; cd /output/'$OUTPUT'/mafft_results; cat /common/'$INPUT' /common/'$REF' > '$newfile'_REF.fasta;\
 		mafft --thread '$THREADS' --auto '$newfile'_REF.fasta > '$newfile'_aligned.fasta ;\
 		chmod -R 777 /output/'$OUTPUT'"

#Correcting Frameshifts
echo "Creating a Frameshift_Correction Container: "
docker run -id -v $COMMON_PATH:/common/ -v $CURRENT_PATH:/output/ --name framecorrect_$TIMESTAMP itvdsbioinfo/frameshift_correction:v1

echo "Running the Frameshift_Correction Container v1: "
docker exec -i framecorrect_$TIMESTAMP  /bin/bash -c "cd /output/'$OUTPUT' ; \
		python3.8 /framecorrect/correct_frameshift.py mafft_results/'$newfile'_aligned.fasta '$newfile'; \
		chmod -R 777 /output/'$OUTPUT'"


echo "Stopping Containeres: "
docker stop mafft_$TIMESTAMP
docker stop framecorrect_$TIMESTAMP

echo "Removing Containeres: "
docker rm mafft_$TIMESTAMP
docker rm framecorrect_$TIMESTAMP

echo "Done!"