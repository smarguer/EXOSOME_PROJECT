 cs_gff=read.delim("pA_analysis/0304_1_RITA.pA.map",stringsAsFactors=F)
 #cs_gff=read.delim("EXOSOME_local/0304_1_RITA.pA.map",stringsAsFactors=F)
 type="cs"
 cs="distal_cs"
 an=cs_gff

 test=list()
 for(i in 1:11)
 {
   test[[i]]=ReadAccumulation(my.col=i,log=F,what=type,length=400,annot=an,cs.type=cs,smooth=0.8)
 }
 names(test)=colnames(chr1)
 #stop("so far so good...")
 pdf("figure_TTS_180413.pdf",width=30,height=10)
 par(mfrow=c(2,4))
 ct=c(1,3,5,5,5,5,10)
 ex=c(2,4,6,7,8,9,11)
 for(i in 1:7)
 {
  print(i)
  if(i==6)
  {
   plotTTS(test[c(ct[i],ex[i])],my.ylim=c(500,2000))
  }
  else
  {
   plotTTS(test[c(ct[i],ex[i])],my.ylim=c(500,1500))
  }
  legend(x="topleft",legend=names(test)[ex[i]],bty="n",cex=3)
 }
 
 #plotTTS(test[c(3,4,5,6)],my.ylim=c(500,1500),my.col=c("black","red","black","red"))
 out=plotTTS(data=test,norm=1,plot=F)
 plot(out,pch=20,col="red",cex=3,xlim=c(0,300),ylim=c(1,2))
 points(out[c(1,3,5,10),],pch=20,col="black",cex=3)
 textxy(out[,1],out[,2],labs=rownames(out),cx=1.2)
 dev.off()
 rm(an,cs,ct,ex,i,type)
