#Script Written by Kanishk Asthana kasthana@eng.ucsd.edu
require('intervals')

macsData=read.table('macs2Out_peaks.xls',comment.char="#",header=TRUE);
homerData=read.table("homerPeaks.txt",comment.char="#", blank.lines.skip=TRUE);
names(homerData)<-c("PeakID","chr","start","end","strand","Normalized_Tag_Count","focus_ratio","findPeaks_Score","Fold_ChangevsLocal","p-valuevsLocal","ClonalFoldChange")


#Sorting for better comparison of column numbers
homerData=homerData[order(homerData$chr),]
macsData=macsData[order(macsData$chr),]


#Formatting macsData for Homer Peak Finder
macsData=macsData[,c(1,2,3,10)]
unused=rep("unused",nrow(macsData));
#setting all strands as + : although it might seem like i am making up information homer also predicts peaks with the 
#majority of or all of peaks on the + strand, so this should not affect the analysis
Strand=rep("+",nrow(macsData));
macsData=cbind(macsData,unused,Strand);

#Getting List of all Chromosomes to Iterate through for Peak Comparison

chromosomes_macs=unique(macsData$chr)
chromosomes_homer=unique(homerData$chr)

chromosomes=chromosomes_homer

#Choosing the larger of the two chromosome lists
if(length(chromosomes_homer)<length(chromosomes_macs)){
  chromosomes=chromosomes_macs
}


mat=lapply(chromosomes, function(chr_name){
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
  
  #Getting Intersections
  presentInHomerButNotInMacs=homer_chromosome_peaks[-1*unlist(interval_overlap(macs_interval,homer_interval)),]
  
  presentInMacsButNotInHomer=macs_chromosome_peaks[-1*unlist(interval_overlap(homer_interval,macs_interval)),]
  

  if(nrow(presentInHomerButNotInMacs)>0){
    write.table(presentInHomerButNotInMacs,file="hmacsinter.txt",append=TRUE, quote=FALSE, col.names=FALSE,row.names=FALSE, sep="\t")
  }
  
  if(nrow(presentInMacsButNotInHomer)>0){
    write.table(presentInMacsButNotInHomer,file="mhomerinter.bed",append=TRUE, quote=FALSE, col.name=FALSE,row.names=FALSE, sep="\t")
  }
  
});



