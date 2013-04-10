

#polishCHIP=function(IP,IN,win=10,log=TRUE)
#{
# if(log==TRUE)
# {
#  IP=toLINEAR(IP)
#  IN=toLINEAR(IN)
# }
##print(i)
# out=rep(0,length(IP))
# for (i in 1:((length(IP)/2)))
# {
#  fac=mean(IN[i:(i+win)])
# 
##print(i)
#  if(fac!=0)
#  {
#   out[i]=mean(IP[i:(i+win)])/fac
#  }
# }
# for(i in ((length(IP)/2)+1):(length(IP)-win))
# {
#  fac=mean(IN[i:(i+win)])
##print(i)
#  if(fac!=0)
#  {
#   out[i+win]=mean(IP[i:(i+win)])/fac
#  }
# }
#
# return(out)
#}


polishCHIP.2=function(IP,IN,win=11,log=TRUE)
{
 if(log==TRUE)
 {
  IP=toLINEAR(IP)
  IN=toLINEAR(IN)
 }
 out=rep(0.00,length(IP))

 for (i in ((1+floor(win/2)):((length(IP)/2)-floor(win/2))))
 {
   fac=mean(IN[(i-floor(win/2)):(i+floor(win/2))],na.rm=T)
#print(i)
  if(fac!=0)
  {
   out[i]=mean(IP[(i-floor(win/2)):(i+floor(win/2))],na.rm=T)/fac
  }
 }
 out[((length(out)/2)+1):length(out)]=out[1:(length(out)/2)]
 out=out/median(out)
 out=round(out,digits=2)
 return(out)
}
