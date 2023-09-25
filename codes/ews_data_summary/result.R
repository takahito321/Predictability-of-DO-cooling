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
    
p<-matrix(NA,length(t1),12)

for(k in 1:6){
    filename<-paste("result",k,".dat",sep="")
    tab<-read.table(filename)
    if(k==1) p[,1]<-tab[,2]
    if(k==1) p[,2]<-tab[,3]
    if(k==2) p[,3]<-tab[,2]
    if(k==2) p[,4]<-tab[,3]
    if(k==3) p[,5]<-tab[,2]
    if(k==3) p[,6]<-tab[,3]
    if(k==4) p[,7]<-tab[,2]
    if(k==4) p[,8]<-tab[,3]
    if(k==5) p[,9]<-tab[,2]
    if(k==5) p[,10]<-tab[,3]
    if(k==6) p[,11]<-tab[,2]
    if(k==6) p[,12]<-tab[,3]
}

p[,1:12]<-p[,12:1]

ii<-c(2, 3, 4, 6, 8, 9, 10, 16, 19, 20, 24, 33) # NGRIP d18O
postscript(file="result_matrix.eps",encoding="WinAnsi.enc",horizontal=TRUE)
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
par(mai = c(0.4, 1.9, 0.1, 0))
laby<-c(expression(paste("GISP2 Ca"^{"2+"}," AC")),
        expression(paste("GISP2 Ca"^{"2+"}," Var")),
        expression(paste("GRIP Ca"^{"2+"}," AC")),
        expression(paste("GRIP Ca"^{"2+"}," Var")),
        expression(paste("NGRIP Ca"^{"2+"}," AC")),
        expression(paste("NGRIP Ca"^{"2+"},"  Var")),
        expression(paste("GISP2 ",delta^{18},"O AC")),
        expression(paste("GISP2 ",delta^{18},"O Var")),
        expression(paste("GRIP ",delta^{18},"O AC")),
        expression(paste("GRIP ",delta^{18},"O Var")),
        expression(paste("NGRIP ",delta^{18},"O AC")),
        expression(paste("NGRIP ",delta^{18},"O Var")))
image.plot(1:length(ii),1:12,p[ii,],zlim=c(0.6,30),xlab="",ylab="",col=rev(magma(100)),xaxt="n",yaxt="n",legend.lab="", legend.width=1, legend.cex=2, axis.args=list(cex.axis=1.3))
axis(1, at=1:length(ii), labels=GI2[ii], cex.axis=1.3) 
axis(2, at=1:12, labels=laby, cex.axis=1.3, las = 1)
rect(0.5, 6.5, 3.5, 10.5, density = 50, col = "gray", angle = -30, border = "transparent")
rect(0.5, 4.5, 2.5, 6.5, density = 50, col = "gray", angle = -30, border = "transparent")
rect(0.5, 0.5, 3.5, 4.5, density = 50, col = "gray", angle = -30, border = "transparent")
text(2,2.5,"Unpublished",cex=1.5)
dev.off()

# the number of robust EWS
# sum(p[ii,]>15)  = 34
# sum(p[ii,]>15)/(12*12-28) = 0.29 (probability)
# sum(p[ii,]>=15) = 36
# sum(p[ii,]>=15)/(12*12-28) = 0.3103448 (probability)

postscript(file="result_barplot.eps",encoding="WinAnsi.enc",height=2,width=12)
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
par(mai = c(0.4, 1.9, 0.1, 0.1))
barplot(dur[ii]/1000,ylab="Duration (kyr)", cex.axis=1.5, cex.lab=1.5, col="rosybrown1")
dev.off()
