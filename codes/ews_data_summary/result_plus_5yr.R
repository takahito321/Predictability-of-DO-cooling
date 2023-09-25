library(palinsol)
library(Hmisc)
library(misty)
library(Kendall)
library(nonlinearTseries)
library(fields)
library(RColorBrewer)
library(viridis)   


tab1<-read.xlsx("../Rasmussen_et_al_2014_QSR_Table_2.xlsx",sheet=1,skip=28)
tab1<-tab1[,1:7]
i1<-c(9,13,15,17,19,21,23,25,29,33,35,37,39,43,53,55,57,61,63,67,69,71,73,75,79,85,87,97,99,103,105,109,111)
i2<-c(2,12,14,16,18,20,22,24,26,30,34,36,38,40,44,54,56,58,62,64,68,70,72,74,76,80,86,88,98,100,104,106,110)
t1<-tab1[i1,3] # GI-start
t2<-tab1[i2,3] # GI-end
u1<-tab1[i1,4] # GI-start uncertainty
u2<-tab1[i2,4] # GI-end uncertainty
t1<-t1-2*4*(u1=="±4")-2*20*(u1=="a")-2*60*(u1=="b")-2*200*(u1=="c")-2*100*(u1=="d")-2*40*(u1=="f")
t2<-t2-2*4*(u2=="±4")+2*20*(u2=="a")+2*60*(u2=="b")+2*200*(u2=="c")+2*100*(u2=="d")+2*40*(u2=="f")
GI<-tab1[i1,6]
GI2<-tab1[i1,7] 
    
t1<--rev(t1)   # reverse
t2<--rev(t2)
GI<-rev(GI) 
GI2<-rev(GI2)
dur<-t2-t1
    
p<-matrix(NA,length(t1),16)

for(k in 1:6){
    filename<-paste("result",k,".dat",sep="")
    tab<-read.table(filename)
    if(k==1) p[,3]<-tab[,2]
    if(k==1) p[,4]<-tab[,3]
    if(k==2) p[,5]<-tab[,2]
    if(k==2) p[,6]<-tab[,3]
    if(k==3) p[,7]<-tab[,2]
    if(k==3) p[,8]<-tab[,3]
    if(k==4) p[,11]<-tab[,2]
    if(k==4) p[,12]<-tab[,3]
    if(k==5) p[,13]<-tab[,2]
    if(k==5) p[,14]<-tab[,3]
    if(k==6) p[,15]<-tab[,2]
    if(k==6) p[,16]<-tab[,3]
}

for(k in 1:2){
    filename<-paste("../ews_data_summary_5yr/result",k,".dat",sep="")
    tab<-read.table(filename)
    if(k==1) p[,1]<-tab[,2]
    if(k==1) p[,2]<-tab[,3]
    if(k==2) p[,9]<-tab[,2]
    if(k==2) p[,10]<-tab[,3]
}

p[,1:16]<-p[,16:1]

ii<-c(1, 2, 3, 4, 6, 8, 9, 10, 16, 19, 20, 24, 33) # NGRIP d18O
postscript(file="result_matrix.eps",encoding="WinAnsi.enc",horizontal=TRUE)
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
par(mai = c(0.4, 0.4, 0.1, 0.2))
laby<-c("AC","Var")
laby<-rep(laby,8)
image.plot(1:length(ii),1:16,p[ii,],zlim=c(0.6,100),xlab="",ylab="",col=rev(magma(100)),xaxt="n",yaxt="n",legend.lab="", legend.width=1, legend.cex=2, axis.args=list(cex.axis=1.2))
axis(1, at=1:length(ii), labels=GI2[ii], cex.axis=1) 
axis(2, at=1:16, labels=laby, cex.axis=1)
rect(0.5, 14.5, 8.5, 16.5, density = 50, col = "gray", angle = -30, border = "transparent")
rect(0.5, 8.5, 4.5, 12.5, density = 50, col = "gray", angle = -30, border = "transparent")
rect(0.5, 6.5, 8.5, 8.5, density = 50, col = "gray", angle = -30, border = "transparent")
rect(0.5, 4.5, 3.5, 6.5, density = 50, col = "gray", angle = -30, border = "transparent")
rect(0.5, 0.5, 4.5, 4.5, density = 50, col = "gray", angle = -30, border = "transparent")
dev.off()

postscript(file="result_matrix_full.eps",encoding="WinAnsi.enc",horizontal=TRUE)
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
par(mai = c(0.5, 0.4, 0.1, 0.2))
laby<-c("AC","Var")
laby<-rep(laby,8)
image.plot(1:33,1:16,p,zlim=c(0.6,100),xlab="",ylab="",col=rev(magma(100)),xaxt="n",yaxt="n",legend.lab="", legend.width=1, legend.cex=2, axis.args=list(cex.axis=1.2))
axis(1, at=1:33, labels=GI2, cex.axis=1) 
axis(2, at=1:16, labels=laby, cex.axis=1)
rect(0.5, 14.5, 13.5, 16.5, density = 50, col = "gray", angle = -30, border = "transparent")
rect(0.5, 8.5, 4.5, 12.5, density = 50, col = "gray", angle = -30, border = "transparent")
rect(0.5, 6.5, 13.5, 8.5, density = 50, col = "gray", angle = -30, border = "transparent")
rect(0.5, 4.5, 3.5, 6.5, density = 50, col = "gray", angle = -30, border = "transparent")
rect(0.5, 0.5, 4.5, 4.5, density = 50, col = "gray", angle = -30, border = "transparent")
dev.off()


postscript(file="result_barplot.eps",encoding="WinAnsi.enc",height=2,width=12)
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
par(mai = c(0.4, 1, 0.1, 0.2))
barplot(dur[ii]/1000,ylab="Duration (kyr)", cex.axis=1.5, cex.lab=1.5)
dev.off()
