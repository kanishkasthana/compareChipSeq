#Script Written by Kanishk Asthana kasthana@eng.ucsd.edu
require('intervals')
macsData=read.table('macsEncodePeaks_Filtered.txt', header=TRUE)
homerData=read.table('homerEncodePeaks_Filtered.txt', header=TRUE )

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

chr_values_mat=sapply(chromosomes, function(chr_name){
  chr_name=as.character(chr_name)
  macs_chromosome_peaks=macsData[macsData$chr==chr_name,]
  column1=nrow(macs_chromosome_peaks)
  homer_chromosome_peaks=homerData[homerData$chr==chr_name,]
  column2=nrow(homer_chromosome_peaks)
  
  #Matrix of Start and End positions for peaks in a given chromosome for Homer peaks
  homer_chromosome_matrix=as.matrix(cbind(homer_chromosome_peaks$start,homer_chromosome_peaks$end));
  
  #Matrix of Start and End positions for peaks
  macs_chromosome_matrix=as.matrix(cbind(macs_chromosome_peaks$start,macs_chromosome_peaks$end));
  
  #Creatging Interval objects for start and end values computed using both homer and macs
  homer_interval=Intervals(homer_chromosome_matrix);
  macs_interval=Intervals(macs_chromosome_matrix);
  
  #Storing the number of values that are common for each chromosome
  #NOTE: Hope this is an efficient implementation or else it can take a long time and will have to implement interval trees in 
  #Java
  column3=length(unlist(interval_overlap(homer_interval,macs_interval)))
  
  return(cbind(column1,column2,column3))
});

#Labelling Comparison Data frame
homerVsMacs=data.frame(t(chr_values_mat));
homerVsMacs=cbind(homerVsMacs,chromosomes);
colnames(homerVsMacs)<-c("Number of MACs peaks", "Number of Homer Peaks", "Shared Peaks","Chromosome Number")

require('xtable')
sink("HomerVSMacs_Encode_Filtered.html")

print(xtable(homerVsMacs),type="html")

#Total Values for data
total_sums=apply(chr_values_mat,1,sum)
print("<br>Percentage of Peaks in MACS shared with HOMER:<br>")
print(100*total_sums[3]/total_sums[1])
print("<br>Percentage of Peaks in HOMER shared with MACS:<br>")
print(100*total_sums[3]/total_sums[2])


