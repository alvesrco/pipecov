# HSV script but works with any viral sequence: This script makes a new reference sequence from de novo assembled scaffolds
# Pavitra Roychoudhury
# Sep 2017
 
# Built to be called from hsv_wgs_pipeline.sh with input arguments specifying input filename
# Requires wgs_functions.R which contains several utility scripts plus multiple R packages listed below

rm(list=ls()); 
sessionInfo();
library(Rsamtools);
library(GenomicAlignments);
library(Biostrings);
library(RCurl);
library(parallel);

#Get latest stable version of wgs_functions.R from github
source('/Scripts/wgs_functions.R'); #or locally
#script<-getURL('https://raw.githubusercontent.com/proychou/ViralWGS/master/wgs_functions.R', sl.verifypeer=FALSE)
# eval(parse(text=script));

#Get args from command line 
args<-(commandArgs(TRUE));
if(length(args)==0){
	print("No arguments supplied.")
}else{
	for(i in 1:length(args)){
		eval(parse(text=args[[i]]))
		print(args[[i]])
	}
}

#For testing
# sampname='GH120016_CGATCCAC-TCGCGCAT_L002'
# scaffname=paste0('/fh/fast/jerome_k/RSV_WGS/contigs/',sampname,'/scaffolds.fasta')
# reffname='./refs/NC_045512.2.fasta'
# ncores

#First import scaffolds and filter by length (>200) and coverage (>10x)
newref<-make_ref_from_assembly(bamfname,reffname)

if(newref==FALSE) print('Failed to generate consensus from scaffolds')
