#Script written by Kanishk Asthana kasthana@eng.ucsd.edu
  
require('intervals')
macsData=read.table('macs2Out_peaks.xls', comment.char='#', blank.lines.skip=TRUE);
macs_column_labels<-unlist(macsData[1,])
macsData<-macsData[2:nrow(macsData),]
names(macsData)<-macs_column_labels

homerData=read.table('homerPeaks.txt',comment.char="#", blank.lines.skip=TRUE);
names(homerData)<-c("PeakID","chr","start","end","strand","Normalized_Tag_Count","focus_ratio","findPeaks_Score","Fold_ChangevsLocal","p-valuevsLocal","ClonalFoldChange")

chromosomes_macs=unique(macsData$chr)
chromosomes_homer=unique(homerData$chr)

chromosomes=chromosomes_homer

#Choosing the larger of the two chromosome lists
if(length(chromosomes_homer)<length(chromosomes_macs)){
  chromosomes=chromosomes_macs
}

  mergedPeaksDataFrame<-sapply(chromosomes, function(chr_name){
  chr_name=as.character(chr_name)
  
  macs_chromosome_peaks=macsData[macsData$chr==chr_name,]
  homer_chromosome_peaks=homerData[homerData$chr==chr_name,]
  
  #Matrix of Start and End positions for peaks in a given chromosome for Homer peaks
  homer_chromosome_matrix=as.matrix(cbind(homer_chromosome_peaks$start,homer_chromosome_peaks$end));
  
  #Matrix of Start and End positions for peaks
  macs_chromosome_matrix=as.matrix(cbind(macs_chromosome_peaks$start,macs_chromosome_peaks$end));
  
  #Creatging Interval objects for start and end values computed using both homer and macs
  homer_interval=Intervals(homer_chromosome_matrix);
  macs_interval=Intervals(macs_chromosome_matrix);

  homer_overlapping_row_numbers=interval_overlap(homer_interval,macs_interval);
  
});
