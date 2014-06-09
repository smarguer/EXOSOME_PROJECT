
#load("NORM_0304_ALL_280512.rda")
#load("../EXOSOME_PROJECT/RT_ann_070612.rda")
#source("toLINEAR.r")




RT_signal=function(annot,c1=chr1,c2=chr2,c3=chr3,cut=1,inter=50,win=20000)
{
 out=data.frame(stringsAsFactors=F)
 j=0

#for(i in 1:30)
 for(i in 1:nrow(annot))
 {
#print(i)
  out[i,1]=row.names(annot)[i]
  out[i,2]=NA
  out[i,3]=NA
  out[i,4]=NA
  out[i,5]=NA
  out[i,6]=NA
  out[i,7]=NA
 
  if (annot[i,"chr"]<4)
  { 
   chr=get(paste("c",annot[i,"chr"],sep='')) 
   end=annot[i,"end"]
   start=annot[i,"start"] 
  #annot[i,which(annot[i,] < 0)]=0
 
   if(annot[i,2] == "-")
   {
   #print("in")
    chr=chr[((length(chr)/2)+1):length(chr)]
    end=annot[i,"start"]
    start=annot[i,"end"]  
   }
#print(end)
#print(row.names(annot)[i]) 
#print(i)
#if (i==j+1000)
#{
# print(i)
# j=i
#}
   out[i,1]=row.names(annot)[i]
   out[i,2]=sum(chr[annot[i,"start"]:annot[i,"end"]],na.rm=T)/(annot[i,"end"]-annot[i,"start"])
   out[i,3]=sum(chr[annot[i,"prom_start"]:annot[i,"prom_end"]],na.rm=T)/(annot[i,"prom_end"]-annot[i,"prom_start"])

   if(is.na(annot[i,"utr5_start"])==FALSE)
   {
    out[i,4]=sum(chr[annot[i,"utr5_start"]:annot[i,"utr5_end"]],na.rm=T)/(annot[i,"utr5_end"]-annot[i,"utr5_start"])
   }
  
   else
   {
    out[i,4]=NA
   }
   
   out[i,5]=sum(chr[annot[i,"RT_start"]:annot[i,"RT_end"]],na.rm=T)/(annot[i,"RT_end"]-annot[i,"RT_start"])
 
   if(is.na(annot[i,"utr3_start"])==FALSE)
   { 
    out[i,6]=sum(chr[annot[i,"utr3_start"]:annot[i,"utr3_end"]],na.rm=T)/(annot[i,"utr3_end"]-annot[i,"utr3_start"])
   }
   
   else
   {
    out[i,6]
   }


####
#extend from orf end
####
   if(annot[i,2] == "+")
   {
   #print("plus")
    if((end+win)>length(chr))
    {
     enddown=length(chr)
    }
    else
    {
     enddown=end+win
    }
    if((start-win)<1)
    {
     startup=1
    }
    else
    {
     startup=start-win
    }
    chr_3ext=chr[end:enddown]
    out[i,7]=which.min(chr_3ext)
    out[i,8]=chr_3ext[which.min(chr_3ext)]
#interval cut off
    chr_3str=ceiling(chr_3ext)
    chr_3str[which(chr_3str > 9)]=9
    out[i,9]=regexpr(paste(rep(0,inter),collapse=''),paste(chr_3str,collapse=''))[1]

    chr_5ext=chr[startup:start]
    chr_5ext=rev(chr_5ext)
    out[i,10]=which.min(chr_5ext)
    out[i,11]=chr_5ext[which.min(chr_5ext)]
#interval cut off
    chr_5str=ceiling(chr_5ext)
    chr_5str[which(chr_5str > 9)]=9
    out[i,12]=regexpr(paste(rep(0,inter),collapse=''),paste(chr_5str,collapse=''))[1]
   }
  
   else if(annot[i,2] == "-")
   {
   #print("minus")
    if((end-win) < 0)
    {
     enddown=1
    }
    else
    {
     enddown=end-win
    } 
    if((start+win)>length(chr))
    {
     startup=length(chr)
    }
    else
    {
     startup=start+win
    }
   
    chr_3ext=chr[enddown:end]
    chr_3ext=rev(chr_3ext)
    out[i,7]=which.min(chr_3ext)
    out[i,8]=chr_3ext[which.min(chr_3ext)]
#interval cut off
    chr_3str=ceiling(chr_3ext)
    chr_3str[which(chr_3str > 9)]=9
    out[i,9]=regexpr(paste(rep(0,inter),collapse=''),paste(chr_3str,collapse=''))[1]

    chr_5ext=chr[start:startup]
    out[i,10]=which.min(chr_5ext)
    out[i,11]=chr_5ext[which.min(chr_5ext)]
#interval cut off
    chr_5str=ceiling(chr_5ext)
    chr_5str[which(chr_5str > 9)]=9
    out[i,12]=regexpr(paste(rep(0,inter),collapse=''),paste(chr_5str,collapse=''))[1]
   }
####
  }
 }
 out[,2]=round(out[,2],3)
 out[,3]=round(out[,3],3)
 out[,4]=round(out[,4],3)
 out[,5]=round(out[,5],3)
 out[,6]=round(out[,6],3)
 out[which(out[,8]>=cut),7]=win
 out[which(out[,11]>=cut),10]=win
 colnames(out)=c("name","gene","prom","utr5","down","utr3","ext3","ext3.min","match3","ext5","ext5.min","match5")
 return(out)
}
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################

#for (i in c(1,7))
#for (i in 1:20)
#for (i in 1:41)
#for (i in 1)
#{
# print(i)
# assign(paste("RT",i,sep=''),RT_signal(c1=toLINEAR(chr1[,i]),c2=toLINEAR(chr2[,i]),c3=toLINEAR(chr3[,i]),annot=ann,inter=25))
# if (i==1)
# {
#  save(RT1,file="RT1_AVERAGE_260612.rda")
# }
#}

#par(mfrow=c(3,7))
#names=colnames(chr1)
#out1=data.frame(stringsAsFactors=F)
PLOT=FALSE
names=colnames(chr1)
if(PLOT==TRUE)
{
 par(mfrow=c(4,5))
 for(i in 1:20)
 {
 a=get(paste("RT",i,sep=''))
 b=get(paste("RT",i,sep=''))
 what="match3"
 # out1[i,1]=names[i]
 # out1[i,2]=cor(log2(a[which((a[,"gene"]*a[,"down"]) !=0),"gene"]),log2(a[which((a[,"gene"]*a[,"down"]) !=0),"down"]),use="complete")
 ##
 # #hist(log2(a[,"down"])-log2(a[,"gene"]),breaks=100)
 plot(log2(RT2[,what]),log2(a[,what]),xlim=c(0,17),ylim=c(0,17),main=names[i],pch=20,cex=0.5)
 #lines(log2(RT1[grep('SPNCRNA',a[,1]),"gene"]),log2(a[grep('SPNCRNA',a[,1]),"gene"]),pch=20,cex=0.5,type="p",col="green")
 #plot(log2(a[,"ext3"]),log2(a[,"match"]),xlim=c(0,17),ylim=c(0,17),main=names[i],pch=20,cex=0.5)
 #lines(log2(a[grep('SPNCRNA',a[,1]),"gene"]),log2(a[grep('SPNCRNA',a[,1]),"prom"]),pch=20,cex=0.5,type="p",col="green")
 #legend(x="bottomright",legend=names[i],cex=2)
 abline(0,1,col="red")
 abline(1,1,lty=2,col="red")
 abline(-1,1,lty=2,col="red")
 }
}

#save(list=c("RT1","RT2","RT3","RT4","RT5","RT6","RT7","RT8","RT9","RT10","RT11","RT12","RT13","RT14","RT15","RT16","RT17","RT18","RT19","RT20","RT21","RT22","RT23","RT24","RT25","RT26","RT27","RT28","RT29","RT30","RT31","RT32","RT33","RT34","RT35","RT36","RT37","RT38","RT39","RT40","RT41"),file="RT_140612.rda")
#save(list=c("RT1","RT2","RT3","RT4","RT5","RT6","RT7","RT8","RT9","RT10","RT11","RT12","RT13","RT14","RT15","RT16","RT17","RT18","RT19","RT20"),file="RT_AVERAGE_260612.rda")
format=FALSE
if(format==TRUE)
{
RT_match=data.frame(stringsAsFactors=F,
RT1[,"match3"],
RT2[,"match3"],
RT3[,"match3"],
RT4[,"match3"],
RT5[,"match3"],
RT6[,"match3"],
RT7[,"match3"],
RT8[,"match3"],
RT9[,"match3"],
RT10[,"match3"],
RT11[,"match3"],
RT12[,"match3"],
RT13[,"match3"],
RT14[,"match3"],
RT15[,"match3"],
RT16[,"match3"],
RT17[,"match3"],
RT18[,"match3"],
RT19[,"match3"],
RT20[,"match3"],
RT1[,"match5"],
RT2[,"match5"],
RT3[,"match5"],
RT4[,"match5"],
RT5[,"match5"],
RT6[,"match5"],
RT7[,"match5"],
RT8[,"match5"],
RT9[,"match5"],
RT10[,"match5"],
RT11[,"match5"],
RT12[,"match5"],
RT13[,"match5"],
RT14[,"match5"],
RT15[,"match5"],
RT16[,"match5"],
RT17[,"match5"],
RT18[,"match5"],
RT19[,"match5"],
RT20[,"match5"])
colnames(RT_match)=c(paste(colnames(chr1),"_3",sep=''),paste(colnames(chr1),"_5",sep=''))

for(i in 1:40)
{
 RT_match[which(RT_match[,i]<1),i]=1
 RT_match[which(is.na(RT_match[,i])==T),i]=1
}

}










