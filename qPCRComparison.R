#Script Written By Kanishk Asthana
require('intervals');
qPCRData=read.table("gabp_qPCR_hg19.bed");
names(qPCRData)<-c("chr","start","end");

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

homerPeaksInOrder=order(homerData$findPeaks_Score);

