# RSV : This script imports bam files and makes a consensus sequence
# Pavitra Roychoudhury
# Adapted from hsv_generate_consensus.R on 6-Mar-19

# Built to be called from hhv6_wgs_pipeline.sh with input arguments specifying input filename
# Requires wgs_functions.R which contains several utility scripts plus multiple R packages listed below

rm(list=ls()); 
sessionInfo();
library(Rsamtools);
library(GenomicAlignments);
library(ShortRead);
library(Biostrings);
library(RCurl);

#Get latest stable version of wgs_functions.R from github
#source('/Scripts/wgs_functions.R'); #or locally
script<-getURL('https://raw.githubusercontent.com/proychou/ViralWGS/master/wgs_functions.R',
							 ssl.verifypeer=FALSE)
eval(parse(text=script));

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

#For testing (these args should come from command line)
# sampname='2016-01040_S451_L001'
# ref='NC_016842' 
# remapped_bamfname

#Files, directories, target site
mapped_reads_folder<-'./mapped_reads/';

#Make consensus sequence--returns TRUE if this worked
conseq<-clean_consensus_hcov(sampname,remapped_bamfname,mappedtoref_bamfname,ref);

#Prepare seqs for annotation -- will make separate folders for A and B
if(conseq==TRUE){
	if(!dir.exists('./annotations_prokka')) dir.create('./annotations_prokka');

  #Write consensus seq to folder for prokka
  fname<-paste('./consensus_seqs/',sampname,'.fasta',sep='')
  con_seq<-readDNAStringSet(fname);
  names(con_seq)<-substring(names(con_seq),1,20); #prokka needs contig name to be <=20 chars long
  sampdir<-paste('./annotations_prokka/',sampname,sep='');
  if(!dir.exists(sampdir)) dir.create(sampdir); #create folder for the sample
  writeXStringSet(con_seq,file=paste(sampdir,'/',sampname,'.fa',sep=''),format='fasta');
	
}else{
	print('Failed to generate consensus sequences.')
}

