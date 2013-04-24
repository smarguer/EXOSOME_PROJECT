
####
plotTTS=function(data,my.ylim=c(0.8,1.5),my.col=NULL,norm=NULL,plot=TRUE)
{
 if(is.null(my.col)==TRUE)
 {
  my.col=c("black","red","blue","green","orange","purple")
 }
 if(is.null(norm)==F)
 {
  for(i in 1:length(data))
  {
   data[[i]]=data[[i]]/data[[i]][norm]
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
  plot(X_axis,data[[1]],type="l",lwd=3,col=my.col[1],ylab="A.U.",xlab="position realtive to TTS",ylim=my.ylim)
  #abline(v=X_axis[which(data[[1]]==max(data[[1]]))],col=my.col[1])
  abline(v=X_axis[which(data[[1]]>quantile(data[[1]],0.98))],col=my.col[1])
  if(length(data) > 1)
  {
   for(i in 2:length(data))
   {
    lines(X_axis,data[[i]],type="l",lwd=3,col=my.col[i])
    #abline(v=X_axis[which(data[[i]]==max(data[[i]]))],col=my.col[i])
    abline(v=X_axis[which(data[[i]]>quantile(data[[i]],0.98))],col=my.col[i])
   }
  }
  abline(v=0,lwd=3,lty=2)
 }
 return(out)
}

####
ReadAccumulation=function(c1=chr1,c2=chr2,c3=chr3,annot=gff,my.col=1,what=c("cs","gff"),cs.type="proximal_cs",length=200,log=T,li=NULL,smooth=0.5,heatmap=FALSE)
{
 if(log==T)
 {
  c1=toLINEAR(c1[,my.col])
  c2=toLINEAR(c2[,my.col])
  c3=toLINEAR(c3[,my.col])
 }
 else
 {
  c1=c1[,my.col]
  c2=c2[,my.col]
  c3=c3[,my.col]
 }
 x=seq(1,((length*2)+1),1)
 if(what=="gff")
 {
  print("gff")
  annotg=annot[which(annot$feature=="gene"),]
  annotg=annotg[which(annotg$chr < 4),]
  out=rep(0,((length*2)+1))
  out1=rep(0,nrow(annotg))
  out2=matrix(0,nrow(annotg),((length*2)+1))
str(out2)
  for (i in 1:nrow(annotg))
  {
   chr=get(paste("c",annotg$chr[i],sep=''))
   if(annotg$strand[i]=="+")
   {
    range=c((annotg$end[i]-length),(annotg$end[i]+length))
    out=cbind(out,chr[range[1]:range[2]])
    hold=out[,2]
    out=rowSums(out,na.rm=T)
   }
   else if (annotg$strand[i]=="-")
   {
    range=c((annotg$start[i]-length),(annotg$start[i]+length))
    range=range+length(chr)
    out=cbind(out,rev(chr[range[1]:range[2]]))
    hold=out[,2]
    out=rowSums(out,na.rm=T)
   }
   if(heatmap==TRUE)
   {
    out2[i,]=hold
   }
  } 
 }

 else if(what=="cs")
 {
  print(cs.type)
  annotg=annot[which(annot$main_utr > 0),]
  if(is.null(li)==FALSE)
  {
   annotg=annotg[which(annotg$name %in% li),]
  }
  out=rep(0,((length*2)+1))
  
  for (i in 1:nrow(annotg))
  {
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





