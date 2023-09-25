library(Hmisc)
library(misty)
library(Kendall)
library(nonlinearTseries) # FFTsurrogate
library(smoother)
library(fields)
library(RColorBrewer)

tab1<-read.xlsx("Rasmussen_et_al_2014_QSR_Table_2.xlsx",sheet=1,skip=28)
tab1<-tab1[,1:8] # 1:8!
i1<-c(9,13,15,17,19,21,23,25,29,33,35,37,39,43,53,55,57,61,63,67,69,71,73,75,79,85,87,97,99,103,105,109,111) # GI-start : row index in tab1
i2<-c(2,12,14,16,18,20,22,24,26,30,34,36,38,40,44,54,56,58,62,64,68,70,72,74,76,80,86,88,98,100,104,106,110) # GI-end : row index in tab1
t1<-tab1[i1,3] # GI-start yr (positive)
t2<-tab1[i2,3] # GI-end yr (positive)
u1<-tab1[i1,4] # GI-start uncertainty: a=20yr (1sigma), b=40-60, c=200, d=100, f=20--40
u2<-tab1[i2,4] # GI-end uncertainty
t1<-t1-2*4*(u1=="±4")-2*20*(u1=="a")-2*60*(u1=="b")-2*200*(u1=="c")-2*100*(u1=="d")-2*40*(u1=="f")
t2<-t2-2*4*(u2=="±4")+2*20*(u2=="a")+2*60*(u2=="b")+2*200*(u2=="c")+2*100*(u2=="d")+2*40*(u2=="f")
GI<-tab1[i1,6]   # label 
GI2<-tab1[i1,7]  # label2
GI3<-tab1[i1,8]  # new! GI-13-14, GI-23.1-22
    
t1<--rev(t1)   # reverse (negative)
t2<--rev(t2)   # reverse (negative)
GI<-rev(GI)    # reverse
GI2<-rev(GI2)  # reverse
GI3<-rev(GI3)  # simplified label

# Proxy time series data : seierstad.xlsx
tab<-read.xlsx("seierstad.xlsx",sheet=3,skip=51)
t<-tab[seq(1,nrow(tab),by=2),1]+10        # start t=30yr
d18O.NGRIP<-tab[seq(1,nrow(tab),by=2),5] 
d18O.GRIP<-tab[seq(1,nrow(tab),by=2),8]   # t=20yr d18O.GRIP=-35.13, t=40yr d18O.GRIP=-35.13
d18O.GISP2<-tab[seq(1,nrow(tab),by=2),11]
d18O.GISP2<-approx(t[!is.na(d18O.GISP2)],d18O.GISP2[!is.na(d18O.GISP2)],t)$y
ca.NGRIP<-tab[seq(1,nrow(tab),by=2),6]
ca.GRIP<-tab[seq(1,nrow(tab),by=2),9]
ca.GISP2<-tab[seq(1,nrow(tab),by=2),12]
ca.NGRIP<-approx(t[!is.na(ca.NGRIP)],ca.NGRIP[!is.na(ca.NGRIP)],t)$y
ca.GRIP<-approx(t[!is.na(ca.GRIP)],ca.GRIP[!is.na(ca.GRIP)],t)$y
ca.GISP2<-approx(t[!is.na(ca.GISP2)],ca.GISP2[!is.na(ca.GISP2)],t)$y
d18O.NGRIP<-as.numeric(d18O.NGRIP[(t>=10000) & (t<=123000)])
d18O.GRIP<-as.numeric(d18O.GRIP[(t>=10000) & (t<=123000)])
d18O.GISP2<-as.numeric(d18O.GISP2[(t>=10000) & (t<=123000)])
ca.NGRIP<-log10(as.numeric(ca.NGRIP[(t>=10000) & (t<=123000)]))
ca.GRIP<-log10(as.numeric(ca.GRIP[(t>=10000) & (t<=123000)]))
ca.GISP2<-log10(as.numeric(ca.GISP2[(t>=10000) & (t<=123000)]))
t<-t[(t>=10000) & (t<=123000)] # positive

t<--rev(t) # reverse (negative)
d18O.NGRIP<-rev(d18O.NGRIP)
d18O.GRIP<-rev(d18O.GRIP)
d18O.GISP2<-rev(d18O.GISP2)
ca.NGRIP<-rev(ca.NGRIP)
ca.GRIP<-rev(ca.GRIP)
ca.GISP2<-rev(ca.GISP2)

## main ###############################################################################################
k<-6        # which proxy?
set.seed(k) # set a seed for random values

# GI[ii] = interstadials
if(k==1) ii<-c(2, 3, 4, 6, 8, 9, 10, 16, 19, 20, 24, 33) # NGRIP d18O
if(k==2) ii<-c(6, 8, 9, 10, 16, 19, 20, 24, 33)    # GRIP d18O
if(k==3) ii<-c(6, 8, 9, 10, 16, 19, 20, 24, 33)    # GISP2 d18O
if(k==4) ii<-c(4, 6, 8, 9, 10, 16, 19, 20, 24, 33) # NGRIP ca2+
if(k==5) ii<-c(6, 8, 9, 10, 16, 19, 20, 24, 33)    # GRIP ca2+
if(k==6) ii<-c(6, 8, 9, 10, 16, 19, 20, 24, 33)    # GISP2 ca2+

# parameter setting
ns<-1000      # sample size for p values
#smoothing<-"gaussian"
smoothing<-"loess"
bandwidth<-c(20,30,40,50,60,70) # %
winsize<-c(20,30,40,50,60)      # %

trend1<-array(NA,dim=c(length(t1),length(winsize),length(bandwidth)))
trend2<-array(NA,dim=c(length(t1),length(winsize),length(bandwidth)))
#trend3<-array(NA,dim=c(length(t1),length(winsize),length(bandwidth)))
p1<-array(NA,dim=c(length(t1),length(winsize),length(bandwidth)))
p2<-array(NA,dim=c(length(t1),length(winsize),length(bandwidth)))
#p3<-array(NA,dim=c(length(t1),length(winsize),length(bandwidth)))

for(i in ii){  # for-loop of each GI
    
# DO number i
tx<-t[t>t1[i] & t<t2[i]]
if(k==1)  x<-d18O.NGRIP[t>t1[i] & t<t2[i]]
if(k==2)  x<-d18O.GRIP[t>t1[i] & t<t2[i]]
if(k==3)  x<-d18O.GISP2[t>t1[i] & t<t2[i]]
if(k==4)  x<-ca.NGRIP[t>t1[i] & t<t2[i]]
if(k==5)  x<-ca.GRIP[t>t1[i] & t<t2[i]]
if(k==6)  x<-ca.GISP2[t>t1[i] & t<t2[i]]
#x<-x[round(0.3*length(x)):length(x)]     # remove first % (not used)
#tx<-tx[round(0.3*length(tx)):length(tx)] # remove first % (not used)
n<-length(x)                   # data length
    
for(l in 1:length(winsize)){
for(m in 1:length(bandwidth)){

bw <- round(n * bandwidth[m]/100) # bandwidth     (actual data size)
wd <- round(n * winsize[l]/100)   # window length (actual data size)
span <- bandwidth[m]/100          # bandwidth for loess (ratio)
if(smoothing=="gaussian") xs <- ksmooth(1:n, x, kernel = "normal", bandwidth = bw, n.points=n, x.points=1:n)$y
if(smoothing=="loess")    xs<-predict(loess(x~tx, span=span, degree=1)) # linear                       
y<-x-xs                        # residuals
v<-rep(NA,length(y))           # variance
ac<-rep(NA,length(y))          # lag-1 autocorrelation
#lam<-rep(NA,length(y))         # lag-1 lambda
#dy<-y[1:n]-c(0,y[-n])          # dy(j)=y(j)-y(j-1). dy[1]=y[1]-0 (unnatural but unused), dy[2]=y[2]-y[1]
ix<-wd:n                   
for(j in ix){
    v[j]<-var(y[(j-wd+1):j])   # window length wd
    Z<-y[(j-wd+1):j]-mean(y[(j-wd+1):j])
    ac[j]<-sum(Z[-1]*Z[-wd])/sum(Z^2)
    #Y<-dy[(j-wd+2):j]-mean(dy[(j-wd+2):j])        # dy[2:wd]     (j=wd), dy[3:(wd+1)]
    #X<-y[(j-wd+1):(j-1)]-mean(y[(j-wd+1):(j-1)])  #  y[1:(wd-1)] (j=wd),  y[2:wd]
    #lam[j]<-sum(X*Y)/sum(X^2)                   
}
fit1<-lm(v[ix]~tx[ix])
fit2<-lm(ac[ix]~tx[ix])
#fit3<-lm(lam[ix]~tx[ix])
trend1[i,l,m]<-as.numeric(fit1$coef[2])
trend2[i,l,m]<-as.numeric(fit2$coef[2])
#trend3[i,l,m]<-as.numeric(fit3$coef[2])

sg<-FFTsurrogate(y, n.samples = ns) # surrogate data [instance x time]
sg<-sg-rowMeans(sg)                 # remove rowMeans
vs<-matrix(NA,ns,length(y))         # variance
acs<-matrix(NA,ns,length(y))        # lag-1 autocorrelation
#lams<-matrix(NA,ns,length(y))       # lag-1 lambda
#dys<-sg[,1:n]-cbind(matrix(0,ns,1),sg[,-n]) # dy(j)=y(j)-y(j-1)
for(j in ix){ 
    vs[,j]<-apply(sg[,(j-wd+1):j],1,var)
    Z<-sg[,(j-wd+1):j]-rowMeans(sg[,(j-wd+1):j])
    acs[,j]<-rowSums(Z[,-1]*Z[,-wd])/rowSums(Z^2)
    #Y<-dys[,(j-wd+2):j]-rowMeans(dys[,(j-wd+2):j])        # dy(j)=y(j)-y(j-1)~y(j-1)
    #X<-sg[,(j-wd+1):(j-1)]-rowMeans(sg[,(j-wd+1):(j-1)])  
    #lams[,j]<-rowSums(X*Y)/rowSums(X^2)
}
vs<-vs[,ix]
acs<-acs[,ix]
#lams<-lams[,ix]
vs<-vs-rowMeans(vs)
acs<-acs-rowMeans(acs)
#lams<-lams-rowMeans(lams)
tm<-matrix(tx[ix],dim(vs)[1],dim(vs)[2],byrow=TRUE)
tm<-tm-rowMeans(tm)
trend1.surrogate<-rowSums(vs*tm)/rowSums(tm^2)
trend2.surrogate<-rowSums(acs*tm)/rowSums(tm^2)
#trend3.surrogate<-rowSums(lams*tm)/rowSums(tm^2)
p1[i,l,m]<-length(trend1.surrogate[trend1.surrogate>trend1[i,l,m]])/ns
p2[i,l,m]<-length(trend2.surrogate[trend2.surrogate>trend2[i,l,m]])/ns
#p3[i,l,m]<-length(trend3.surrogate[trend3.surrogate>trend3[i,l,m]])/ns
}
}
}


#cols <- brewer.pal(3, "BuGn")
cols <- brewer.pal(3, "PRGn")
#cols <- rev(brewer.pal(3, "RdBu"))
pal <- colorRampPalette(cols)
## plot ############################################################################
for(iii in 1:2){
if(k==1 && iii==1) postscript(file="ews_data_summary/ews_data_summary_NGRIP_d18O_1.eps", onefile=FALSE, horizontal=FALSE,encoding="WinAnsi.enc")
if(k==1 && iii==2) postscript(file="ews_data_summary/ews_data_summary_NGRIP_d18O_2.eps", onefile=FALSE, horizontal=FALSE,encoding="WinAnsi.enc")
if(k==2 && iii==1) postscript(file="ews_data_summary/ews_data_summary_GRIP_d18O_1.eps", onefile=FALSE, horizontal=FALSE,encoding="WinAnsi.enc")
if(k==2 && iii==2) postscript(file="ews_data_summary/ews_data_summary_GRIP_d18O_2.eps", onefile=FALSE, horizontal=FALSE,encoding="WinAnsi.enc")
if(k==3 && iii==1) postscript(file="ews_data_summary/ews_data_summary_GISP2_d18O_1.eps", onefile=FALSE, horizontal=FALSE,encoding="WinAnsi.enc")
if(k==3 && iii==2) postscript(file="ews_data_summary/ews_data_summary_GISP2_d18O_2.eps", onefile=FALSE, horizontal=FALSE,encoding="WinAnsi.enc")
if(k==4 && iii==1) postscript(file="ews_data_summary/ews_data_summary_NGRIP_ca_1.eps", onefile=FALSE, horizontal=FALSE,encoding="WinAnsi.enc")
if(k==4 && iii==2) postscript(file="ews_data_summary/ews_data_summary_NGRIP_ca_2.eps", onefile=FALSE, horizontal=FALSE,encoding="WinAnsi.enc")
if(k==5 && iii==1) postscript(file="ews_data_summary/ews_data_summary_GRIP_ca_1.eps", onefile=FALSE, horizontal=FALSE,encoding="WinAnsi.enc")
if(k==5 && iii==2) postscript(file="ews_data_summary/ews_data_summary_GRIP_ca_2.eps", onefile=FALSE, horizontal=FALSE,encoding="WinAnsi.enc")
if(k==6 && iii==1) postscript(file="ews_data_summary/ews_data_summary_GISP2_ca_1.eps", onefile=FALSE, horizontal=FALSE,encoding="WinAnsi.enc")
if(k==6 && iii==2) postscript(file="ews_data_summary/ews_data_summary_GISP2_ca_2.eps", onefile=FALSE, horizontal=FALSE,encoding="WinAnsi.enc")

if(iii==1) ii2<-ii[1:round(length(ii)/2)]
if(iii==2) ii2<-ii[(round(length(ii)/2)+1):length(ii)]

layout( matrix(1:(2*length(ii2)), ncol=2) ) # ncol=6
par(mar=c(0,0,0,0))
par(oma=c(3,2,2,3))
par(mai = c(0.25, 0.65, 0.15, 0.4))

for(i in ii2){
zmax<-max(abs(1000*trend1[i,,]))    
image.plot(winsize,bandwidth,1000*trend1[i,,],zlim=c(-zmax,zmax),cex.lab=2,col=pal(20),xlab="",ylab="")
for(l in 1:length(winsize)){
    for(m in 1:length(bandwidth)){
        if(p1[i,l,m]<0.05){
            points(winsize[l],bandwidth[m],pch=4,cex=2,lwd=2)
        } else if(p1[i,l,m]<0.1){
            points(winsize[l],bandwidth[m],pch=1,cex=1,lwd=2)
        }   
    }
}
mtext("Smoothing span (%)", side = 2, line = 2.5, col=1)
mtext(paste("GI-",GI3[i],sep=""), side = 2, line = 4.5, col=1, lwd=2) #  change!
if(i==min(ii2)) mtext("Variance", side = 3, line = 1, col=1)
if(i==max(ii2)) mtext("Rolling window size (%)", side = 1, line = 2.8, col=1)
 
zmax<-max(abs(1000*trend2[i,,]))    
image.plot(winsize,bandwidth,1000*trend2[i,,],zlim=c(-zmax,zmax),cex.lab=2,col=pal(20),xlab="",ylab="")
for(l in 1:length(winsize)){
    for(m in 1:length(bandwidth)){
        if(p2[i,l,m]<0.05){
            points(winsize[l],bandwidth[m],pch=4,cex=2,lwd=2)
        } else if(p2[i,l,m]<0.1){
            points(winsize[l],bandwidth[m],pch=1,cex=1,lwd=2)
        }   
    }
}
if(i==min(ii2)) mtext("Lag-1 autocorrelation", side = 3, line = 1, col=1)
if(i==max(ii2)) mtext("Rolling window size (%)", side = 1, line = 2.8, col=1)
    
#image.plot(winsize,bandwidth,trend3[i,,],cex.lab=2,col=pal(20),xlab="",ylab="")
#for(l in 1:length(winsize)){
#    for(m in 1:length(bandwidth)){
#        if(p3[i,l,m]<0.05){
#            points(winsize[l],bandwidth[m],pch=4,cex=2,lwd=2)
#        } else if(p3[i,l,m]<0.1){
#            points(winsize[l],bandwidth[m],pch=1,cex=1,lwd=2)
#        }   
#    }
#}
#if(i==min(ii)) mtext(expression(lambda), side = 3, line = 1, col=1)
#if(i==max(ii)) mtext("Window size", side = 1, line = 2.8, col=1)
}
dev.off()
} 



## write out results ########################################################################
#if(FALSE){
dat<-matrix(0,length(t1),3) # 10 DO cases, 3 = i, var, ac  
for(i in ii){
    dat[i,]<-c(i,
      sum(p1[i,,]<0.05), #/length(winsize)/length(bandwidth)*100,
      sum(p2[i,,]<0.05)  #/length(winsize)/length(bandwidth)*100
      )
}
filename<-paste("ews_data_summary/result",k,".dat",sep="")
write.table(dat, file = filename, col.names=F, row.names=F) 
#}



## plot ############################################################################
if(k==1){
postscript(file="ews_data_summary/ews_data_summary_example.eps", onefile=FALSE, horizontal=TRUE,encoding="WinAnsi.enc")

ii2<-2

layout( matrix(1:(2*length(ii2)), ncol=2) ) # ncol=6
par(mar=c(0,0,0,0))
par(oma=c(4,0,4,3))
par(mai = c(1.2, 1.2, 0.9, 0.7))

for(i in ii2){
zmax<-max(abs(1000*trend1[i,,]))   
image.plot(winsize,bandwidth,1000*trend1[i,,],zlim=c(-zmax,zmax),cex.lab=2,col=pal(20),xlab="",ylab="",cex.axis=1.4)
for(l in 1:length(winsize)){
    for(m in 1:length(bandwidth)){
        if(p1[i,l,m]<0.05){
            points(winsize[l],bandwidth[m],pch=4,cex=2.5,lwd=2)
        } else if(p1[i,l,m]<0.1){
            points(winsize[l],bandwidth[m],pch=1,cex=1,lwd=2)
        }   
    }
}
mtext("Smoothing span (%)", side = 2, line = 3, col=1,cex=1.4)
#mtext(GI[i], side = 2, line = 5, col=1, lwd=2,cex=1.4)
if(i==min(ii2)) mtext("Variance", side = 3, line = 1, col=1,cex=1.4)
if(i==max(ii2)) mtext("Rolling window size (%)", side = 1, line = 2.8, col=1,cex=1.4)
    
zmax<-max(abs(1000*trend2[i,,]))   
image.plot(winsize,bandwidth,1000*trend2[i,,],zlim=c(-zmax,zmax),cex.lab=2,col=pal(20),xlab="",ylab="",cex.axis=1.4)
for(l in 1:length(winsize)){
    for(m in 1:length(bandwidth)){
        if(p2[i,l,m]<0.05){
            points(winsize[l],bandwidth[m],pch=4,cex=2.5,lwd=2)
        } else if(p2[i,l,m]<0.1){
            points(winsize[l],bandwidth[m],pch=1,cex=1,lwd=2)
        }   
    }
}
if(i==min(ii2)) mtext("Lag-1 autocorrelation", side = 3, line = 1, col=1,cex=1.4)
if(i==max(ii2)) mtext("Rolling window size (%)", side = 1, line = 2.8, col=1,cex=1.4)
    
#image.plot(winsize,bandwidth,trend3[i,,],cex.lab=2,col=pal(20),xlab="",ylab="")
#for(l in 1:length(winsize)){
#    for(m in 1:length(bandwidth)){
#        if(p3[i,l,m]<0.05){
#            points(winsize[l],bandwidth[m],pch=4,cex=2,lwd=2)
#        } else if(p3[i,l,m]<0.1){
#            points(winsize[l],bandwidth[m],pch=1,cex=1,lwd=2)
#        }   
#    }
#}
#if(i==min(ii)) mtext(expression(lambda), side = 3, line = 1, col=1)
#if(i==max(ii)) mtext("Window size", side = 1, line = 2.8, col=1)
}
dev.off()
}

