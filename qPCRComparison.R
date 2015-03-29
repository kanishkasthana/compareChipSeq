#Script Written By Kanishk Asthana
require('intervals');
qPCRData=read.table("gabp_qPCR_hg19.bed");
names(qPCRData)<-c("chr","start","end");

qPCR_matrix=as.matrix(cbind(qPCRData$start,qPCRData$end))
qPCR_interval=Intervals(qPCR_matrix);
macsData=read.table('GABP_peaks.xls',comment.char="#",header=TRUE);
homerData=read.table('homerGABPPeaks.txt',comment.char="#", blank.lines.skip=TRUE);
names(homerData)<-c("PeakID","chr","start","end","strand","Normalized_Tag_Count","focus_ratio","findPeaks_Score","Fold_ChangevsLocal","p-valuevsLocal","ClonalFoldChange")


#Sorting for better comparison of column numbers
homerData=homerData[order(homerData$chr),]
macsData=macsData[order(macsData$chr),]


#Getting List of all Chromosomes to Iterate through for Peak Comparison

chromosomes_macs=unique(macsData$chr)
chromosomes_homer=unique(homerData$chr)

chromosomes=chromosomes_homer

#Choosing the larger of the two chromosome lists
if(length(chromosomes_homer)<length(chromosomes_macs)){
  chromosomes=chromosomes_macs
}
#ordering peaks in decreasing order
homerPeaksInOrder=order(-homerData$findPeaks_Score);
#print(length(homerPeaksInOrder))
#getting array for number of qPCR detected sites 
range=20000

qPCRPositives=sapply(seq(1,range,by=250),function(peakNumber)
{
  
  homerPeaks=homerData[homerPeaksInOrder[1:peakNumber],]
  
  
  chr_values_mat=sapply(chromosomes, function(chr_name){
  
    chr_name=as.character(chr_name)
    qPCRPeaks=qPCRData[qPCRData$chr==chr_name,]
    homer_chromosome_peaks=homerPeaks[homerPeaks$chr==chr_name,]
  
    #Matrix of Start and End positions for peaks in a given chromosome for Homer peaks
    homer_chromosome_matrix=as.matrix(cbind(homer_chromosome_peaks$start,homer_chromosome_peaks$end));
  
    #Matrix of Start and End positions for peaks
    qPCRmatrix=as.matrix(cbind(qPCRPeaks$start,qPCRPeaks$end));
  
    #Creatging Interval objects for start and end values computed using both homer and macs
    homer_interval=Intervals(homer_chromosome_matrix);
    qPCR_interval=Intervals(qPCRmatrix);
  
    PeaksCovered=length(unlist(interval_overlap(homer_interval,qPCR_interval)))
  
});

return(sum(chr_values_mat));

});

plot(seq(1,range,by=250),qPCRPositives,type='o',col="green",xlab="Peak Rank",ylab="Number of qPCR Positives Detected")
title(main="Comparison of Homer vs Macs in detecting qPCR verified positives for GABP")

macsPeaksInOrder=order(-macsData$X.log10.pvalue.)

qPCRPositives=sapply(seq(1,range,by=250),function(peakNumber)
{
  
  macsPeaks=macsData[macsPeaksInOrder[1:peakNumber],]
  
  chr_values_mat=sapply(chromosomes, function(chr_name){
  
    chr_name=as.character(chr_name)
    qPCRPeaks=qPCRData[qPCRData$chr==chr_name,]
    macs_chromosome_peaks=macsPeaks[macsPeaks$chr==chr_name,]
  
    #Matrix of Start and End positions for peaks in a given chromosome for Homer peaks
    macs_chromosome_matrix=as.matrix(cbind(macs_chromosome_peaks$start,macs_chromosome_peaks$end));
  
    #Matrix of Start and End positions for peaks
    qPCRmatrix=as.matrix(cbind(qPCRPeaks$start,qPCRPeaks$end));
  
    #Creatging Interval objects for start and end values computed using both homer and macs
    macs_interval=Intervals(macs_chromosome_matrix);
    qPCR_interval=Intervals(qPCRmatrix);
  
    PeaksCovered=length(unlist(interval_overlap(macs_interval,qPCR_interval)))
  
  });

return(sum(chr_values_mat));

});

lines(seq(1,range,by=250),qPCRPositives,type='o',col="red");