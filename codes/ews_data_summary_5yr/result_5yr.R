library(palinsol)
library(Hmisc)
library(misty)
library(Kendall)
library(nonlinearTseries)
library(fields)
library(RColorBrewer)
library(viridis)   


tab1<-read.xlsx("../Rasmussen_et_al_2014_QSR_Table_2.xlsx",sheet=1,skip=28)
tab1<-tab1[,1:8]
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
GI3<-tab1[i1,8] 
    
t1<--rev(t1)   # reverse
t2<--rev(t2)
GI<-rev(GI) 
GI2<-rev(GI2)
GI3<-rev(GI3)
dur<-t2-t1
    
p<-matrix(NA,length(t1),4)

for(k in 1:2){
    filename<-paste("result",k,".dat",sep="")
    tab<-read.table(filename)
    if(k==1) p[,1]<-tab[,2]
    if(k==1) p[,2]<-tab[,3]
    if(k==2) p[,3]<-tab[,2]
    if(k==2) p[,4]<-tab[,3]
}

p[,1:4]<-p[,4:1]

ii<-c(14, 16, 17, 19, 20, 21, 22, 24, 25, 27, 33)  
cl<-c(brewer.pal(n = 9, name = "Blues")[2:6],brewer.pal(n = 9, name = "RdPu")[5:9])
postscript(file="result_matrix_5yr.eps",encoding="WinAnsi.enc",horizontal=TRUE)
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
par(mai = c(0.5, 1.9, 0.2, 0))
laby<-c(expression(paste("NGRIP Ca"^{"2+"}," AC")),
        expression(paste("NGRIP Ca"^{"2+"},"  Var")),
        expression(paste("NGRIP ",delta^{18},"O AC")),
        expression(paste("NGRIP ",delta^{18},"O Var")))
image.plot(1:length(ii),1:4,p[ii,],zlim=c(1.51,28.51),xlab="",ylab="",col=cl,xaxt="n",yaxt="n",legend.lab="", legend.width=1, legend.cex=2, axis.args=list(cex.axis=1.3))
axis(1, at=1:length(ii), labels=GI2[ii], cex.axis=1.3) 
axis(2, at=1:4, labels=laby, cex.axis=1.3, las = 1)
rect(2.5, 0.5, 3.5, 2.5, density = 50, col = "gray", angle = -30, border = "transparent")
rect(7.5, 0.5, 9.5, 2.5, density = 50, col = "gray", angle = -30, border = "transparent")
text(3,1.5,"Gaps in data",cex=1.5,srt=90)
text(8.5,1.5,"Gaps in data",cex=1.5,srt=90)
dev.off()

# the number of robust EWS
# sum(p[ii,]>15)  = 34
# sum(p[ii,]>15)/(12*12-28) = 0.29 (probability)
# sum(p[ii,]>=15) = 36
# sum(p[ii,]>=15)/(12*12-28) = 0.3103448 (probability)

postscript(file="result_barplot_5yr.eps",encoding="WinAnsi.enc",height=2,width=12)
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
par(mai = c(0.4, 1.9, 0.1, 0.1))
barplot(dur[ii]/1000,ylab="Duration (kyr)", cex.axis=1.5, cex.lab=1.5, col="orange")
abline(h=1.5,lty=2,lwd=2,col=2)
dev.off()
