#Collection of functions for working with WGS data
#Pavitra Roychoudhury
#Aug 2017

#Return the number of mapped reads in a bam file
n_mapped_reads<-function(bamfname){
  require(Rsamtools)
  indexBam(bamfname)
  if(file.exists(bamfname)&class(try(scanBamHeader(bamfname),silent=T))!='try-error'){
    return(idxstatsBam(bamfname)$mapped)
  }else{
    return(NA)
  }
}



#Make a new reference from scaffolds
make_ref_from_assembly<-function(bamfname,reffname){
	require(Rsamtools);
	require(GenomicAlignments);
	require(parallel)
	ncores<-detectCores();
	
	#Read reference sequence
	ref_seq<-readDNAStringSet(reffname);
	
	if(!is.na(bamfname)&class(try(scanBamHeader(bamfname),silent=T))!='try-error'){
		
		#Index bam if required
		if(!file.exists(paste(bamfname,'.bai',sep=''))){
			baifname<-indexBam(bamfname); 
		}else{
			baifname<-paste(bamfname,'.bai',sep='');
		}
		
		#Import bam file
		params<-ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
												 what=c('qname','rname','strand','pos','qwidth','mapq','cigar','seq'));
		gal<-readGAlignments(bamfname,index=baifname,param=params);

		#Remove any contigs with width <200 bp
		#gal<-gal[width(gal)>200];
		
		#First lay contigs on reference space--this removes insertions and produces a seq of the same length as ref
		qseq_on_ref<-sequenceLayer(mcols(gal)$seq,cigar(gal),from="query",to="reference");
		qseq_on_ref_aligned<-stackStrings(qseq_on_ref,1,max(mcols(gal)$pos+qwidth(gal)-1,width(ref_seq)),
																			shift=mcols(gal)$pos-1,Lpadding.letter='N',Rpadding.letter='N');
		
		#Make a consensus matrix and get a consensus sequence from the aligned scaffolds
		cm<-consensusMatrix(qseq_on_ref_aligned,as.prob=T,shift=0)[c('A','C','G','T','N','-'),];
		# cm[c('N','-'),]<-0;
		cm['N',]<-0;
		cm<-apply(cm,2,function(x)if(all(x==0))return(x) else return(x/sum(x)));
		cm['N',colSums(cm)==0]<-1;
		con_seq<-DNAStringSet(gsub('\\?','N',consensusString(cm,threshold=0.25)));
		con_seq<-DNAStringSet(gsub('\\+','N',con_seq));
		
		
		#Now fill in the Ns with the reference
		temp<-as.matrix(con_seq);
		temp[temp=='N']<-as.matrix(ref_seq)[temp=='N'];
		con_seq<-DNAStringSet(paste0(temp,collapse=''));
		names(con_seq)<-sub('.bam','_consensus',basename(bamfname));
		
		#Look for insertions in bam cigar string
		cigs_ref<-cigarRangesAlongReferenceSpace(cigar(gal),with.ops=F,ops='I',
																						 reduce.ranges=T,drop.empty.ranges=F,
																						 pos=mcols(gal)$pos);
		cigs_query<-cigarRangesAlongQuerySpace(cigar(gal),ops='I',with.ops=F,
																					 reduce.ranges=T,drop.empty.ranges=F);
		all_ins<-mclapply(c(1:length(cigs_query)),function(i)
		  extractAt(mcols(gal)$seq[i],cigs_query[[i]])[[1]]);
		
		#Merge all insertions
		all_ins_merged<-do.call('rbind',mclapply(c(1:length(cigs_ref)),function(i)
		  return(data.frame(
		    start_ref=start(cigs_ref[[i]]),end_ref=end(cigs_ref[[i]]),
		    start_query=start(cigs_query[[i]]),end_query=end(cigs_query[[i]]),
		    ins_seq=all_ins[[i]],width_ins=width(all_ins[[i]]))),
		  mc.cores=ncores));
		all_ins_merged<-all_ins_merged[order(all_ins_merged$end_ref),];
		
		# write.csv(all_ins_merged,'./testing/all_ins.csv',row.names=F);
		
		#TO DO: Check for overlaps--should be minimal since scaffolds don't usually overlap that much
		if(any(table(all_ins_merged$start_ref)>1)){
		  print('Overlapping insertions')
		  #not the best way, but just pick the first for now
		  all_ins_merged<-all_ins_merged[!duplicated(all_ins_merged[,c('start_ref','end_ref')]),];
		}
		
		#Now the beauty part of inserting the strings back in
		#Split ref seq by the insert positions
		if(nrow(all_ins_merged)!=0){
		  new_strs<-DNAStringSet(rep('',nrow(all_ins_merged)+1))
		  for(i in 1:nrow(all_ins_merged)){
		    if(i==1){
		      new_strs[i]<-paste0(extractAt(con_seq,IRanges(start=1,end=all_ins_merged$end_ref[i]))[[1]],
		                          all_ins_merged$ins_seq[i]);
		    }else{
		      new_strs[i]<-paste0(extractAt(con_seq,IRanges(start=all_ins_merged$start_ref[i-1],
		                                                    end=all_ins_merged$end_ref[i]))[[1]],
		                          all_ins_merged$ins_seq[i]);
		    }
		  }
		  
		  #Last bit
		  new_strs[i+1]<-paste0(extractAt(con_seq,IRanges(start=all_ins_merged$start_ref[i],
		                                                  end=width(con_seq)))[[1]])
		  temp_str<-paste0(as.character(new_strs),collapse='');
		  
		  #Remove gaps to get final sequence
		  con_seq_final<-DNAStringSet(gsub('-','',temp_str));
		
		#No insertions
		}else{
		  con_seq_final<-con_seq;
		}
		names(con_seq_final)<-sub('.bam','_consensus',basename(bamfname));
		
		if(!dir.exists('./ref_for_remapping')) dir.create('./ref_for_remapping');
		writeXStringSet(con_seq_final,
		                paste0('./ref_for_remapping/',names(con_seq_final),'.fasta'));
		
		#Delete bai file
		file.remove(baifname);
		
	}else{
		print('Bam file could not be opened.')
		return(NA)
	}
}



#Takes in a bam file, produces consensus sequence
generate_consensus<-function(bamfname){
  require(Rsamtools)
  require(GenomicAlignments)
  require(Biostrings)
	require(parallel)
	ncores<-detectCores()
	
	#for testing this function--comment out or remove
	# bamfname<-'./testing/ABI-HHV6A_S385_L001_A.sorted.bam'
	
  if(!is.na(bamfname)&class(try(scanBamHeader(bamfname),silent=T))!='try-error'){
    
  	#Index bam if required
    if(!file.exists(paste(bamfname,'.bai',sep=''))){
      baifname<-indexBam(bamfname); 
    }else{
      baifname<-paste(bamfname,'.bai',sep='');
    }
    
    #Import bam file
    params<-ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
                         what=c('qname','rname','strand','pos','qwidth','mapq','cigar','seq'));
    gal<-readGAlignments(bamfname,index=baifname,param=params);
    # summary(gal);
    
    #Remove any contigs with mapq <2 -- this leads to a loss of a lot of the DR seqs even though there are reads there
    # gal<-gal[mcols(gal)$mapq>2];
    
    #First lay reads on reference space--this doesn't include insertions
    qseq_on_ref<-sequenceLayer(mcols(gal)$seq,cigar(gal),from="query",to="reference");

    #Make a consensus matrix and get a consensus sequence from the aligned scaffolds
    # cm<-consensusMatrix(qseq_on_ref,as.prob=T,shift=start(gal)-1,width=seqlengths(gal))[c('A','C','G','T','N','-'),];
    # cm['N',colSums(cm)==0]<-1;
    
    #Edit to include a coverage threshold
    cm<-consensusMatrix(qseq_on_ref,as.prob=F,shift=start(gal)-1,width=seqlengths(gal))[c('A','C','G','T','N','-'),];
    poor_cov<-which(colSums(cm)<0);
    cm<-apply(cm,2,function(x)x/sum(x));
    cm[,poor_cov]<-0;
    cm['N',poor_cov]<-1;
    
    tmp_str<-strsplit(consensusString(cm,ambiguityMap='?',threshold=0.5),'')[[1]];
    ambig_sites<-which(tmp_str=='?');
    ambig_bases<-unlist(lapply(ambig_sites,function(i){mixedbase<-paste(names(cm[,i])[cm[,i]>0],collapse=''); 
    			 if(mixedbase%in%IUPAC_CODE_MAP) return(names(IUPAC_CODE_MAP)[IUPAC_CODE_MAP==mixedbase]) 
    else return('N')}));
    tmp_str[ambig_sites]<-ambig_bases
    con_seq<-DNAStringSet(paste0(tmp_str,collapse='')); 
    names(con_seq)<-sub('.bam','_consensus',basename(bamfname));
    rm(tmp_str);
    
    #Remove gaps and leading and trailing Ns to get final sequence
    con_seq_trimmed<-DNAStringSet(gsub("N*N$",'',gsub("^N*",'',as.character(con_seq))));
    con_seq_final<-DNAStringSet(gsub('-','',as.character(con_seq_trimmed)));
    names(con_seq_final)<-sub('.bam','_consensus',basename(bamfname));
    
    #Delete bai file
    file.remove(baifname);
    
    return(con_seq_final);
  }else{
    return(NA)
  }
}

clean_consensus_hsv<-function(sampname,merged_bam_folder,mapped_reads_folder){
  require(Rsamtools); 
  require(GenomicAlignments);
  require(Biostrings);
	sampname<-paste0(sampname,'_');
  mapping_stats<-data.frame(ref=c('hsv1_ref','hsv2_sd90e','hsv2_ref_hg52'),
                            bamfname_merged=c(grep(sampname,list.files(merged_bam_folder,'_hsv1_ref.*bam$',full.names=T),value=T),
                                              grep(sampname,list.files(merged_bam_folder,'_hsv2_sd90e.*bam$',full.names=T),value=T),
                                              grep(sampname,list.files(merged_bam_folder,'_hsv2_ref_hg52.*bam$',full.names=T),value=T)),
  													bamfname_mapped=c(grep(sampname,list.files(mapped_reads_folder,'_hsv1_ref.*bam$',full.names=T),value=T),
  																						grep(sampname,list.files(mapped_reads_folder,'_hsv2_sd90e.*bam$',full.names=T),value=T),
  																						grep(sampname,list.files(mapped_reads_folder,'_hsv2_ref_hg52.*bam$',full.names=T),value=T)),
  													mapped_reads_ref=0,mapped_reads_assemblyref=0,perc_Ns=0,num_Ns=0,width=0,
                            stringsAsFactors=F);
  
  #Import mapped reads + assembly and generate consensus
  con_seqs<-lapply(mapping_stats$bamfname_merged,generate_consensus);
  if(!dir.exists('./consensus_seqs_all')) dir.create('./consensus_seqs_all');
  dummyvar<-lapply(con_seqs,function(x)
  	writeXStringSet(x,file=paste('./consensus_seqs_all/',names(x),'.fasta',sep=''),format='fasta'));
  rm(dummyvar)
  
  #Compute #mapped reads and %Ns
  mapping_stats$mapped_reads_ref<-unlist(lapply(mapping_stats$bamfname_mapped,n_mapped_reads));
  mapping_stats$mapped_reads_assemblyref<-unlist(lapply(mapping_stats$bamfname_merged,n_mapped_reads));
  mapping_stats$num_Ns<-unlist(lapply(con_seqs,function(x)sum(letterFrequency(x,c('N','+')))));
  mapping_stats$width<-unlist(lapply(con_seqs,width));
  mapping_stats$perc_Ns<-100*mapping_stats$num_Ns/mapping_stats$width;
  if(!dir.exists('./stats/')) dir.create('./stats/');
  write.csv(mapping_stats,file=paste('./stats/',sampname,'_mappingstats.csv',sep=''),row.names=F);
  
  return(TRUE)
}

clean_consensus_hhv6<-function(sampname,merged_bam_folder,mapped_reads_folder){
	require(Rsamtools); 
	require(GenomicAlignments);
	require(Biostrings);
	mapping_stats<-data.frame(
		ref=c('hhv6A_ref_U1102','hhv6B_ref_z29'),
		bamfname_merged=c(grep(paste0('\\/',sampname,'_'),list.files(merged_bam_folder,"_hhv6A_ref_U1102.*bam$",full.names=T),value=T),
											grep(paste0('\\/',sampname,'_'),list.files(merged_bam_folder,'_hhv6B_ref_z29.*bam$',full.names=T),value=T)),
		bamfname_mapped=c(grep(paste0('\\/',sampname,'_'),list.files(mapped_reads_folder,'_hhv6A_ref_U1102.*bam$',full.names=T),value=T),
											grep(paste0('\\/',sampname,'_'),list.files(mapped_reads_folder,'_hhv6B_ref_z29.*bam$',full.names=T),value=T)),
		mapped_reads_ref=0,mapped_reads_assemblyref=0,perc_Ns=0,num_Ns=0,width=0,
		stringsAsFactors=F);
	
	#Import mapped reads + assembly and generate consensus
	con_seqs<-lapply(mapping_stats$bamfname_merged,generate_consensus);
	if(!dir.exists('./consensus_seqs_all')) dir.create('./consensus_seqs_all');
	dummyvar<-lapply(con_seqs,function(x)
		writeXStringSet(x,file=paste('./consensus_seqs_all/',names(x),'.fasta',sep=''),format='fasta'));
	rm(dummyvar)
	
	#Compute #mapped reads and %Ns
	mapping_stats$mapped_reads_ref<-unlist(lapply(mapping_stats$bamfname_mapped,n_mapped_reads));
	mapping_stats$mapped_reads_assemblyref<-unlist(lapply(mapping_stats$bamfname_merged,n_mapped_reads));
	mapping_stats$num_Ns<-unlist(lapply(con_seqs,function(x)sum(letterFrequency(x,c('N','+')))));
	mapping_stats$width<-unlist(lapply(con_seqs,width));
	mapping_stats$perc_Ns<-100*mapping_stats$num_Ns/mapping_stats$width;
	if(!dir.exists('./stats/')) dir.create('./stats/');
	write.csv(mapping_stats,file=paste('./stats/',sampname,'_mappingstats.csv',sep=''),row.names=F);
	
	return(TRUE)
}



clean_consensus_hiv<-function(sampname,merged_bam_folder,mapped_reads_folder){
	require(Rsamtools); 
	require(GenomicAlignments);
	require(Biostrings);
	mapping_stats<-data.frame(
		ref=c('hiv_hxb2_ref'),
		bamfname_merged=c(grep(sampname,list.files(merged_bam_folder,"_hiv_hxb2_ref.*bam$",full.names=T),value=T)),
		bamfname_mapped=c(grep(sampname,list.files(mapped_reads_folder,'_hiv_hxb2_ref.*bam$',full.names=T),value=T)),
		mapped_reads_ref=0,mapped_reads_assemblyref=0,perc_Ns=0,num_Ns=0,width=0,
		stringsAsFactors=F);
	
	#Import mapped reads + assembly and generate consensus
	con_seqs<-lapply(mapping_stats$bamfname_merged,generate_consensus);
	if(!dir.exists('./consensus_seqs_all')) dir.create('./consensus_seqs_all');
	dummyvar<-lapply(con_seqs,function(x)
		writeXStringSet(x,file=paste('./consensus_seqs_all/',names(x),'.fasta',sep=''),format='fasta'));
	rm(dummyvar)
	
	#Compute #mapped reads and %Ns
	mapping_stats$mapped_reads_ref<-unlist(lapply(mapping_stats$bamfname_mapped,n_mapped_reads));
	mapping_stats$mapped_reads_assemblyref<-unlist(lapply(mapping_stats$bamfname_merged,n_mapped_reads));
	mapping_stats$num_Ns<-unlist(lapply(con_seqs,function(x)sum(letterFrequency(x,c('N','+')))));
	mapping_stats$width<-unlist(lapply(con_seqs,width));
	mapping_stats$perc_Ns<-100*mapping_stats$num_Ns/mapping_stats$width;
	if(!dir.exists('./stats/')) dir.create('./stats/');
	write.csv(mapping_stats,file=paste('./stats/',sampname,'_mappingstats.csv',sep=''),row.names=F);
	
	return(TRUE)
}




clean_consensus_hhv8<-function(sampname,merged_bam_folder,mapped_reads_folder){
  require(Rsamtools); 
  require(GenomicAlignments);
  require(Biostrings);
  sampname<-paste0(sampname,'_');
  mapping_stats<-data.frame(ref='hhv8_ref',
                            bamfname_merged=grep(sampname,list.files(merged_bam_folder,'_hhv8.*bam$',full.names=T),value=T),
                            bamfname_mapped=grep(sampname,list.files(mapped_reads_folder,'_hhv8.*bam$',full.names=T),value=T),
                            mapped_reads_ref=0,mapped_reads_assemblyref=0,perc_Ns=0,num_Ns=0,width=0,
                            stringsAsFactors=F);
  
  #Import mapped reads + assembly and generate consensus
  con_seqs<-lapply(mapping_stats$bamfname_merged,generate_consensus);
  if(!dir.exists('./consensus_seqs_all')) dir.create('./consensus_seqs_all');
  dummyvar<-lapply(con_seqs,function(x)
    writeXStringSet(x,file=paste('./consensus_seqs_all/',names(x),'.fasta',sep=''),format='fasta'));
  rm(dummyvar)
  
  #Compute #mapped reads and %Ns
  mapping_stats$mapped_reads_ref<-unlist(lapply(mapping_stats$bamfname_mapped,n_mapped_reads));
  mapping_stats$mapped_reads_assemblyref<-unlist(lapply(mapping_stats$bamfname_merged,n_mapped_reads));
  mapping_stats$num_Ns<-unlist(lapply(con_seqs,function(x)sum(letterFrequency(x,c('N','+')))));
  mapping_stats$width<-unlist(lapply(con_seqs,width));
  mapping_stats$perc_Ns<-100*mapping_stats$num_Ns/mapping_stats$width;
  if(!dir.exists('./stats/')) dir.create('./stats/');
  write.csv(mapping_stats,file=paste('./stats/',sampname,'_mappingstats.csv',sep=''),row.names=F);
  
  return(TRUE)
}

clean_consensus_rsv<-function(sampname,merged_bam_folder,mapped_reads_folder){
	require(Rsamtools); 
	require(GenomicAlignments);
	require(Biostrings);
	sampname<-paste0(sampname,'_');
	mapping_stats<-data.frame(ref=c('rsvA_ref','rsvB_ref'),
														bamfname_merged=c(grep(sampname,list.files(merged_bam_folder,'_rsvA_ref.*bam$',full.names=T),value=T),
																							grep(sampname,list.files(merged_bam_folder,'_rsvB_ref.*bam$',full.names=T),value=T)),
														bamfname_mapped=c(grep(sampname,list.files(mapped_reads_folder,'_rsvA_ref.*bam$',full.names=T),value=T),
																							grep(sampname,list.files(mapped_reads_folder,'_rsvB_ref.*bam$',full.names=T),value=T)),
														mapped_reads_ref=0,mapped_reads_assemblyref=0,perc_Ns=0,num_Ns=0,width=0,
														stringsAsFactors=F);
	
	#Import mapped reads + assembly and generate consensus
	con_seqs<-lapply(mapping_stats$bamfname_merged,generate_consensus);
	if(!dir.exists('./consensus_seqs_all')) dir.create('./consensus_seqs_all');
	dummyvar<-lapply(con_seqs,function(x)
		writeXStringSet(x,file=paste('./consensus_seqs_all/',names(x),'.fasta',sep=''),format='fasta'));
	rm(dummyvar)
	
	#Compute #mapped reads and %Ns
	mapping_stats$mapped_reads_ref<-unlist(lapply(mapping_stats$bamfname_mapped,n_mapped_reads));
	mapping_stats$mapped_reads_assemblyref<-unlist(lapply(mapping_stats$bamfname_merged,n_mapped_reads));
	mapping_stats$num_Ns<-unlist(lapply(con_seqs,function(x)sum(letterFrequency(x,c('N','+')))));
	mapping_stats$width<-unlist(lapply(con_seqs,width));
	mapping_stats$perc_Ns<-100*mapping_stats$num_Ns/mapping_stats$width;
	if(!dir.exists('./stats/')) dir.create('./stats/');
	write.csv(mapping_stats,file=paste('./stats/',sampname,'_mappingstats.csv',sep=''),row.names=F);
	
	return(TRUE)
}

#Measles (added Aug 2019)
clean_consensus_measles<-function(sampname,merged_bam_folder,mapped_reads_folder){
	require(Rsamtools); 
	require(GenomicAlignments);
	require(Biostrings);
	sampname<-paste0(sampname,'_');
	mapping_stats<-data.frame(ref='measles_ref',
														bamfname_merged=grep(sampname,list.files(merged_bam_folder,'_measles_ref.*bam$',full.names=T),value=T),
														bamfname_mapped=grep(sampname,list.files(mapped_reads_folder,'_measles_ref.*bam$',full.names=T),value=T),
														mapped_reads_ref=0,mapped_reads_assemblyref=0,perc_Ns=0,num_Ns=0,width=0,
														stringsAsFactors=F);
	
	#Import mapped reads + assembly and generate consensus
	con_seqs<-lapply(mapping_stats$bamfname_merged,generate_consensus);
	if(!dir.exists('./consensus_seqs_all')) dir.create('./consensus_seqs_all');
	dummyvar<-lapply(con_seqs,function(x)
		writeXStringSet(x,file=paste('./consensus_seqs_all/',names(x),'.fasta',sep=''),format='fasta'));
	rm(dummyvar)
	
	#Compute #mapped reads and %Ns
	mapping_stats$mapped_reads_ref<-unlist(lapply(mapping_stats$bamfname_mapped,n_mapped_reads));
	mapping_stats$mapped_reads_assemblyref<-unlist(lapply(mapping_stats$bamfname_merged,n_mapped_reads));
	mapping_stats$num_Ns<-unlist(lapply(con_seqs,function(x)sum(letterFrequency(x,c('N','+')))));
	mapping_stats$width<-unlist(lapply(con_seqs,width));
	mapping_stats$perc_Ns<-100*mapping_stats$num_Ns/mapping_stats$width;
	if(!dir.exists('./stats/')) dir.create('./stats/');
	write.csv(mapping_stats,file=paste('./stats/',sampname,'_mappingstats.csv',sep=''),row.names=F);
	
	return(TRUE)
}

#Treponema (added Dec 2019)
clean_consensus_tp<-function(sampname,merged_bam_folder,mapped_reads_folder,ref){
  require(Rsamtools); 
  require(GenomicAlignments);
  require(Biostrings);
  mapping_stats<-data.frame(ref=ref,
                            bamfname_merged=grep(sampname,list.files(merged_bam_folder,'*.bam$',full.names=T),value=T),
                            bamfname_mapped=grep(sampname,list.files(mapped_reads_folder,'*.bam$',full.names=T),value=T),
                            mapped_reads_ref=0,mapped_reads_assemblyref=0,perc_Ns=0,num_Ns=0,width=0,
                            stringsAsFactors=F);
  
  #Import mapped reads + assembly and generate consensus
  con_seqs<-lapply(mapping_stats$bamfname_merged,generate_consensus);
  if(!dir.exists('./consensus_seqs_all')) dir.create('./consensus_seqs_all');
  dummyvar<-lapply(con_seqs,function(x)
    writeXStringSet(x,file=paste('./consensus_seqs_all/',names(x),'.fasta',sep=''),format='fasta'));
  rm(dummyvar)
  
  #Compute #mapped reads and %Ns
  mapping_stats$mapped_reads_ref<-unlist(lapply(mapping_stats$bamfname_mapped,n_mapped_reads));
  mapping_stats$mapped_reads_assemblyref<-unlist(lapply(mapping_stats$bamfname_merged,n_mapped_reads));
  mapping_stats$num_Ns<-unlist(lapply(con_seqs,function(x)sum(letterFrequency(x,c('N','+')))));
  mapping_stats$width<-unlist(lapply(con_seqs,width));
  mapping_stats$perc_Ns<-100*mapping_stats$num_Ns/mapping_stats$width;
  if(!dir.exists('./stats/')) dir.create('./stats/');
  write.csv(mapping_stats,file=paste('./stats/',sampname,'_mappingstats.csv',sep=''),row.names=F);
  
  return(TRUE)
}

#hCoV (added Mar 2020)
clean_consensus_hcov<-function(sampname,remapped_bamfname,mappedtoref_bamfname,ref){
  require(Rsamtools); 
  require(GenomicAlignments);
  require(Biostrings);
  mapping_stats<-data.frame(ref=ref,
                            remapped_bam=remapped_bamfname,
                            mappedtoref_bam=mappedtoref_bamfname,
                            mapped_reads_ref=0,mapped_reads_assemblyref=0,perc_Ns=0,num_Ns=0,width=0,
                            stringsAsFactors=F);
  
  #Import mapped reads + assembly and generate consensus
  con_seq<-generate_consensus(mapping_stats$remapped_bam);
  if(!dir.exists('./consensus_seqs')) dir.create('./consensus_seqs');
  writeXStringSet(con_seq,file=paste('./consensus_seqs/',sampname,'.fasta',sep=''),format='fasta');
  
  #Compute #mapped reads and %Ns
  mapping_stats$mapped_reads_ref<-unlist(lapply(mapping_stats$mappedtoref_bam,n_mapped_reads));
  mapping_stats$mapped_reads_assemblyref<-unlist(lapply(mapping_stats$remapped_bam,n_mapped_reads));
  mapping_stats$num_Ns<-sum(letterFrequency(con_seq,c('N','+')));
  mapping_stats$width<-width(con_seq);
  mapping_stats$perc_Ns<-100*mapping_stats$num_Ns/mapping_stats$width;
  if(!dir.exists('./stats/')) dir.create('./stats/');
  write.csv(mapping_stats,file=paste('./stats/',sampname,'_mappingstats.csv',sep=''),row.names=F);
  
  return(TRUE)
}



#Find coverage at each position in the alignment
cov_by_pos<-function(bamfname){
	require(Rsamtools);
	require(GenomicAlignments);
	
	if(file.exists(bamfname)&class(try(scanBamHeader(bamfname),silent=T))!='try-error'){
		#Import alignment
		if(file.exists(paste(bamfname,'.bai',sep='')))
			file.remove(paste(bamfname,'.bai',sep='')); #remove any old index files
		baifname<-indexBam(bamfname); #Make an index file
		params<-ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
												 what=c('qname','rname','strand','pos','qwidth','mapq','cigar','seq'));
		gal<-readGAlignments(bamfname,index=baifname,param=params);
		cov<-coverage(gal);
		file.remove(baifname);
		return(cov)
		
	}else{
		return(NA)
	}
}



#Compute coverage stats
get_coverage<-function(bamfname){
 	if(length(bamfname)==0){
 		mapped<-NA; avg_cov<-NA;
 		sd_cov<-NA; min_cov<-NA; max_cov<-NA;
 	}else if(file.exists(bamfname)&class(try(scanBamHeader(bamfname),silent=T))!='try-error'){
  	require(Rsamtools);
  	require(GenomicAlignments);
    #Import alignment
    if(file.exists(paste(bamfname,'.bai',sep='')))
      file.remove(paste(bamfname,'.bai',sep='')); #remove any old index files
    baifname<-indexBam(bamfname); #Make an index file
    params<-ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
                         what=c('qname','rname','strand','pos','qwidth','mapq','cigar','seq'));
    gal<-readGAlignments(bamfname,index=baifname,param=params);
    # summary(gal);
    cov<-coverage(gal);
    mapped<-length(gal);
    avg_cov<-mean(cov);
    sd_cov<-sd(cov);
    min_cov<-min(cov);
    max_cov<-max(cov);
    file.remove(baifname);
  }else{
    mapped<-NA; avg_cov<-NA;
    sd_cov<-NA; min_cov<-NA; max_cov<-NA;
  }
  return(data.frame(mapped,avg_cov,sd_cov,min_cov,max_cov))
}

#Extracts number of reads and read widths from html report generated by fastqc
fastqc_readstats<-function(fname){
	require(rvest)
  if(file.exists(fname)){
    tmp_fastqc<-read_html(fname);
    tmp_table<-html_table(tmp_fastqc)[[1]];
    fastq_reads<-as.numeric(tmp_table[tmp_table$Measure=='Total Sequences','Value']);
    fastq_width<-tmp_table[tmp_table$Measure=='Sequence length','Value']; #returns single number for raw reads and range for trimmed
    gc<-as.numeric(tmp_table[tmp_table$Measure=='%GC','Value']);
  }else{
    fastq_reads<-NA;
    fastq_width<-NA;
    gc<-NA;
  }
  return(data.frame(fastq_reads,fastq_width,gc,stringsAsFactors=F));
}


#Compute stats on a consensus seq (or really any fasta file)
conseq_stats<-function(fname){
  require(Biostrings)
  if(length(fname)==0){
  	width<-NA; Ns<-NA; percNs<-NA;
  }else if(file.exists(fname)){
  	conseq<-readDNAStringSet(fname,format='fasta');
  	width<-width(conseq);
  	Ns<-sum(letterFrequency(conseq,c('N','+')));
  	percNs<-100*Ns/width;
  }else{
  	width<-NA; Ns<-NA; percNs<-NA;
  }
  return(data.frame(width,Ns,percNs));
}



#VCF to data frame for a vcf generated by Lofreq
vcf_to_df<-function(vcf_fname,sampid){
	require(VariantAnnotation);
	vcf<-readVcf(vcf_fname);
	results<-data.frame(samp_id=sampid,pos=start(rowRanges(vcf)),af=info(vcf)$AF,dp=info(vcf)$DP,ref=ref(vcf),
											alt=unlist(alt(vcf)),stringsAsFactors=F);
	results$snpid<-paste(results$ref,'_',results$pos,'_',results$alt,sep='');
	results$major_af<-unlist(lapply(results$af,function(x)max(x,1-x)));
	results$minor_af<-unlist(lapply(results$af,function(x)min(x,1-x)));
	return(results)
}


#Extract VRC samp year and ID from the fastq file name
get_year<-function(in_string){
	yr<-strsplit(in_string,"-")[[1]][1];
	if(grepl("^[0-9]{1,2}_(19[0-9][0-9]|20[0,1][0-9])",yr)){
		return(strsplit(yr,'_')[[1]][2]);
	}else if(!grepl("19[0-9][0-9]|20[0,1][0-9]",yr)){
		return(NA)
	}else{
		return(yr)
	}
}
get_sampid<-function(in_string){
	if(!is.na(get_year(in_string))){
		return(strsplit(strsplit(in_string,'-')[[1]][2],'_')[[1]][1]);
	}else{
		return(NA);
	}
}



