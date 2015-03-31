#Script Written by Kanishk Asthana kasthana@eng.ucsd.edu

#Script formats MACS data in the format required for HOMER peak calling.
macsData=read.table('macs2Out_peaks.xls',comment.char="#",header=TRUE);

#Formatting macsData for Homer Peak Finder
macsData=macsData[,c(1,2,3,10)]
unused=rep("unused",nrow(macsData));
#setting all strands as + : although it might seem like i am making up information homer also predicts peaks with the 
#majority of or all of peaks on the + strand, so this should not affect the analysis
Strand=rep("+",nrow(macsData));
macsData=cbind(macsData,unused,Strand);

if(nrow(macsData)>0){
  write.table(macsData,file="macsData.bed",append=TRUE, quote=FALSE, col.name=FALSE,row.names=FALSE, sep="\t")
}