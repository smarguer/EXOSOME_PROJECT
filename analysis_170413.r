
f1=FALSE

if(f1=TRUE)
{
par(mfrow=c(2,6))
boxplot(test[which(row.names(test) %in% names[1:500]),])
boxplot(test[which(row.names(test) %in% names[500:1000]),])
boxplot(test[which(row.names(test) %in% names[1000:1500]),])
boxplot(test[which(row.names(test) %in% names[1500:2000]),])
boxplot(test[which(row.names(test) %in% names[2000:2500]),])
boxplot(test[which(row.names(test) %in% names[2500:3000]),])
boxplot(test[which(row.names(test) %in% names[3000:3500]),])
boxplot(test[which(row.names(test) %in% names[3500:4000]),])
boxplot(test[which(row.names(test) %in% names[4000:4500]),])
boxplot(test[which(row.names(test) %in% names[4500:5000]),])
boxplot(test[which(row.names(test) %in% names[5000:length(names)]),])
} 

load("CHIP_EXO_100413.rda")
chr1C=chr1[,c("wt-1","dis3-1","wt-2","dis3-2")]
chr2C=chr2[,c("wt-1","dis3-1","wt-2","dis3-2")]
chr3C=chr3[,c("wt-1","dis3-1","wt-2","dis3-2")]
load("ALL_080612.rda")
chr1R=chr1[,c("wtTHY-1","dis3-1","dis3-2")]
chr2R=chr2[,c("wtTHY-1","dis3-1","dis3-2")]
chr3R=chr3[,c("wtTHY-1","dis3-1","dis3-2")]
rm(chr1,chr2,chr3)
gc()
colnames(chr1C)=c("CH_wt-1","CH_dis3-1","CH_wt-2","CH_dis3-2")
colnames(chr2C)=c("CH_wt-1","CH_dis3-1","CH_wt-2","CH_dis3-2")
colnames(chr3C)=c("CH_wt-1","CH_dis3-1","CH_wt-2","CH_dis3-2")
colnames(chr1R)=c("R_wtTHY-1","R_dis3-1","R_dis3-2")
colnames(chr2R)=c("R_wtTHY-1","R_dis3-1","R_dis3-2")
colnames(chr3R)=c("R_wtTHY-1","R_dis3-1","R_dis3-2")

chr1=cbind(chr1C,chr1R)
chr2=cbind(chr2C,chr2R)
chr3=cbind(chr3C,chr3R)
chr1[,1:4]=toLOG(chr1[,1:4])
chr2[,1:4]=toLOG(chr2[,1:4])
chr3[,1:4]=toLOG(chr3[,1:4])
chr1=round(chr1,2)
chr2=round(chr2,2)
chr3=round(chr3,2)

