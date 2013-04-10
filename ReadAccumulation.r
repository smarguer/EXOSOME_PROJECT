
ReadAccumulation=function(c1=chr1,c2=chr2,c3=chr3,annot=gff,my.col=1,what=c("cs","gff"),cs.type="proximal_cs",length=200)
{
 c1=toLINEAR(c1[,my.col])
 c2=toLINEAR(c2[,my.col])
 c3=toLINEAR(c3[,my.col])
print(cs.type)
 if(what=="gff")
 {
  annotg=annot[which(annot$feature=="gene"),]
  annotg=annotg[which(annotg$chr < 4),]
  out=rep(0,((length*2)+1))
  for (i in 1:nrow(annotg))
  {
   chr=get(paste("c",annotg$chr[i],sep=''))

   if(annotg$strand[i]=="+")
   {
    range=c((annotg$end[i]-length),(annotg$end[i]+length))
    out=cbind(out,chr[range[1]:range[2]])
    out=rowSums(out,na.rm=T)
   }
   else if (annotg$strand[i]=="-")
   {
    range=c((annotg$start[i]-length),(annotg$start[i]+length))
    range=range+length(chr)
    out=cbind(out,rev(chr[range[1]:range[2]]))
    out=rowSums(out,na.rm=T)
   }
  } 
 }

 else if(what=="cs")
 {
  annotg=annot[which(annot$main_utr > 0),]
  #annotg=annot
  out=rep(0,((length*2)+1))
  #print("2")
  
  for (i in 1:nrow(annotg))
  {
  #print(i)
   chr=get(paste("c",annotg$chr[i],sep=''))

   if(annotg$strand[i]=="+")
   {
    range=c((annotg[i,cs.type]-length),(annotg[i,cs.type]+length))
    out=cbind(out,chr[range[1]:range[2]])
    out=rowSums(out,na.rm=T)
   }
   else if (annotg$strand[i]=="-")
   {
    range=c((annotg[i,cs.type]-length),(annotg[i,cs.type]+length))
    range=range+length(chr)
    out=cbind(out,rev(chr[range[1]:range[2]]))
    out=rowSums(out,na.rm=T)
   }
  } 
 }
 return(round(out,2))
}





