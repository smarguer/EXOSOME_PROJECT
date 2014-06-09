
####
plotTTS=function(data,my.ylim=c(0.8,1.5),my.col=NULL,norm=NULL,plot=TRUE,my.log="",y_norm=NULL,my.colour=c("black","black","black","dark blue","dark blue","dark green","dark green","red","red","orange","orange"))
{
 if(is.null(my.col)==TRUE)
 {
  #my.col=c("black","red","blue","green","orange","purple")
##col for RNA exosome figures##
  #my.col=c("dark grey","dark grey","dark blue","dark blue","dark green","dark green","olivedrab","olivedrab","orange","orange")
##col for chip-seq 1##
  #my.col=c("dark grey","dark blue","dark green","orange")
##col for chip-seq Tfs1##
  #my.col=c("dark grey","dark blue","red","purple")
  #my.col=c("grey","grey","grey","grey","light blue","blue","dark blue","purple")
##FIGURES FINAL
  #my.col=c("black","black","dark blue","dark blue","dark green","dark green","red","red","orange","orange")
my.col=my.colour
 }
 if(is.null(norm)==F)
 {
  for(i in 1:length(data))
  {
   if(norm=="median")
   {
    data[[i]]=data[[i]]/median(data[[i]],na.rm=T)
   }
   else if (norm=="fac")
   {
    data[[i]]=data[[i]]/10000
   }
   else
   {
    data[[i]]=data[[i]]/data[[i]][norm]
   }
  }
 }
 X_axis=seq(-floor(length(data[[1]])/2),floor(length(data[[1]])/2),1)

 out=matrix(0,length(data),2)
 row.names(out)=names(data)
 colnames(out)=c("RT","Max")
 
 for(i in 1:length(data))
 {
  out[i,1]=tail(X_axis[which(data[[i]]==max(data[[i]]))],n=1)
  out[i,2]=max(data[[i]])
 }
 
 if(plot==TRUE)
 {
#  plot(X_axis,data[[1]],type="l",lwd=6,col=my.col[1],ylab=expression(paste("Normalised hits per nucleotide [x",10^4,"]"),sep=''),xlab="position realtive to wt main poly(A) cleavage site",ylim=my.ylim,log=my.log,cex.axis=1.8,cex.lab=2)
  plot(X_axis,data[[1]],type="l",lwd=6,col=my.col[1],ylab="Normalised hits per nucleotide",xlab="position realtive to wt main poly(A) cleavage site",ylim=my.ylim,log=my.log,cex.axis=1.8,cex.lab=2)
  #abline(v=X_axis[which(data[[1]]==max(data[[1]]))],col=my.col[1])
  #abline(v=X_axis[which(data[[1]]>quantile(data[[1]],0.98))],col=my.col[1])
  if(length(data) > 1)
  {
   for(i in 2:length(data))
   {
    lines(X_axis,data[[i]],type="l",lwd=6,col=my.col[i])
    #abline(v=X_axis[which(data[[i]]==max(data[[i]]))],col=my.col[i])
    #abline(v=X_axis[which(data[[i]]>quantile(data[[i]],0.98))],col=my.col[i])
   }
  }
  abline(v=0,lwd=3,lty=2)
 }
 return(out)
}

####
ReadAccumulation=function(c1=chr1,c2=chr2,c3=chr3,annot=gff,my.col=NULL,what=c("cs","gff"),cs.type="proximal_cs",length=200,log=T,li=NULL,smooth=0.5,heatmap=FALSE)
{
 
 if(log==T)
 {
  if(is.null(my.col)==TRUE)
  {
   c1=toLINEAR(c1)
   c2=toLINEAR(c2)
   c3=toLINEAR(c3)
  }
  else
  {
   c1=toLINEAR(c1[,my.col])
   c2=toLINEAR(c2[,my.col])
   c3=toLINEAR(c3[,my.col])
  }
 }
 else
 {
  if(is.null(my.col)==FALSE)
  {
   c1=c1[,my.col]
   c2=c2[,my.col]
   c3=c3[,my.col]
  }
 }
 x=seq(1,((length*2)+1),1)
####
 if(what=="gff")
 {
  print("gff")
  annotg=annot[which(annot$feature=="gene"),]
  annotg=annotg[which(annotg$chr < 4),]
  out=rep(0,((length*2)+1))
  out1=rep(0,nrow(annotg))
  out2=matrix(0,nrow(annotg),((length*2)+1))
  for (i in 1:nrow(annotg))
  #for (i in 1:10)
  {
   print(annotg$Name[i])
   chr=get(paste("c",annotg$chr[i],sep=''))
   if(annotg$strand[i]=="+")
   {
    print("p")
    range=c((annotg$end[i]-length),(annotg$end[i]+length))
    out=cbind(out,chr[range[1]:range[2]])
    hold=out[,2]
    print(hold[1])
    out=rowSums(out,na.rm=T)
   }
   else if (annotg$strand[i]=="-")
   {
    print("m")
    range=c((annotg$start[i]-length),(annotg$start[i]+length))
    range=range+(length(chr)/2)
    out=cbind(out,rev(chr[range[1]:range[2]]))
    hold=out[,2]
    print(hold[1])
    out=rowSums(out,na.rm=T)
   }
   if(heatmap==TRUE)
   {
    if(is.null(smooth)==F)
    {
     hold=smooth.spline(x,hold, spar=smooth)$y
    }
    out2[i,]=hold
   }
  } 
 }
####
 else if(what=="cs")
 {
  print(cs.type)
  annotg=annot[which(annot$main_utr > 0),]
  if(is.null(li)==FALSE)
  {
   annotg=annotg[which(annotg$name %in% li),]
  }
  out=rep(0,((length*2)+1))
  out1=rep(0,nrow(annotg))
  out2=matrix(0,nrow(annotg),((length*2)+1)) 
  for (i in 1:nrow(annotg))
  {
   chr=get(paste("c",annotg$chr[i],sep=''))

   if(annotg$strand[i]=="+")
   {
    range=c((annotg[i,cs.type]-length),(annotg[i,cs.type]+length))
    out=cbind(out,chr[range[1]:range[2]])
    hold=out[,2]
    out=rowSums(out,na.rm=T)
   }
   else if (annotg$strand[i]=="-")
   {
    range=c((annotg[i,cs.type]-length),(annotg[i,cs.type]+length))
    range=range+(length(chr)/2)
    out=cbind(out,rev(chr[range[1]:range[2]]))
    hold=out[,2]
    out=rowSums(out,na.rm=T)
   }
   if(heatmap==TRUE)
   {
    if(is.null(smooth)==F)
    {
     hold=smooth.spline(x,hold, spar=smooth)$y
    }
    out2[i,]=hold
   }
  }
 }
 if(is.null(smooth)==F)
 {
  out=smooth.spline(x,out, spar=smooth)$y
 }
 if(heatmap==TRUE)
 {
  out=out2 
 }
 return(round(out,2))
}

#######
#######
####
Screen_ReadAccumulation=function(c1=chr1,c2=chr2,c3=chr3,annot=gff,my.col=NULL,what=c("cs","gff"),cs.type="proximal_cs",length=200,log=T,li=NULL,smooth=0.5,heatmap=FALSE)
{
 
 if(log==T)
 {
  if(is.null(my.col)==TRUE)
  {
   c1=toLINEAR(c1)
   c2=toLINEAR(c2)
   c3=toLINEAR(c3)
  }
  else
  {
   c1=toLINEAR(c1[,my.col])
   c2=toLINEAR(c2[,my.col])
   c3=toLINEAR(c3[,my.col])
  }
 }
 else
 {
  if(is.null(my.col)==FALSE)
  {
   c1=c1[,my.col]
   c2=c2[,my.col]
   c3=c3[,my.col]
  }
 }
 x=seq(1,((length*2)+1),1)
####
 if(what=="gff")
 {
  #print("gff")
  annotg=annot[which(annot$feature=="gene"),]
  annotg=annotg[which(annotg$chr < 4),]
  out3=data.frame()
  for (i in 1:nrow(annotg))
  #for (i in 1:10)
  {
   print(annotg$Name[i])
   chr=get(paste("c",annotg$chr[i],sep=''))
   if(annotg$strand[i]=="+")
   {
    #print("p")
    range=c((annotg$end[i]-length),(annotg$end[i]+length))
    out3[i,1]=annotg$Name[i]
    out3[i,2]=sum(chr[range[1]:range[2]],na.rm=T)
    print(sum(chr[range[1]:range[2]],na.rm=T))
   }
   else if (annotg$strand[i]=="-")
   {
    #print("m")
    range=c((annotg$start[i]-length),(annotg$start[i]+length))
    range=range+(length(chr)/2)
    out3[i,1]=annotg$Name[i]
    out3[i,2]=sum(chr[range[1]:range[2]],na.rm=T)
    print(sum(chr[range[1]:range[2]],na.rm=T))
   }
  } 
 }
####
 return((out3))
}





