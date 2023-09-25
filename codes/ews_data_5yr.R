library(Hmisc)
library(misty)
library(Kendall)
library(nonlinearTseries) # FFTsurrogate
library(smoother)
library(RColorBrewer)
library(gdata)
cl<-brewer.pal(11,"RdBu")

tab1<-read.xlsx("Rasmussen_et_al_2014_QSR_Table_2.xlsx",sheet=1,skip=28)
tab1<-tab1[,1:8] # change also this!
trb<-0.5*c(110940+110640,89800+88920,79240+78740,78080+77760,74320+74100,2*56800,49280+49120,44560+44280,42240+42450,36860+36580) # times of rebound events: 25a, 22, 21.1c-b-a (x2), 20a, 16, 13a, 12a, 11, 8a
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
GI3<-tab1[i1,8] # new! GI-13-14, GI-23.1-22
    
t1<--rev(t1)   # reverse
t2<--rev(t2)
GI<-rev(GI) 
GI2<-rev(GI2)
GI3<-rev(GI3)  # simplified label

tab<-read.xls("NGRIP_d18O_and_dust_5cm.xls",sheet=3,header=T) # 60 ka
tr<--rev(tab[,4])
xr<-rev(tab[,2])
#spl<-spline(tr, xr, n = 1+(59940-9530)/5 , method = "fmm", xmax = -9530, xmin = -59940, ties = mean)
tt<-seq(-59940,-9530,by=5)
spl<-approx(tr, xr, tt)
t<-spl$x
d18O.NGRIP<-spl$y

xr2<-log10(rev(tab[,3]))
#spl<-spline(tr, xr2, n = 1+(59940-9530)/5 , method = "fmm", xmax = -9530, xmin = -59940, ties = mean)
spl<-approx(tr, xr2, tt)
tc<-spl$x
ca.NGRIP<-spl$y

# which((t2-t1)>=400)
# c(14, 16, 17, 19, 20, 21, 22, 24, 25, 27, 33) >300 yr
# c(16, 19, 20, 21, 22, 24, 25, 33) > 400 yr
# c(16, 19, 20, 21, 22, 24, 25, 33) > 500 yr


# d18O
for(i in c(14, 16, 17, 19, 20, 21, 22, 24, 25, 27, 33)){
    ts<-t[t>t1[i] & t<t2[i]]
    to<-tr[tr>t1[i] & tr<t2[i]] # cut from original
    xo<-xr[tr>t1[i] & tr<t2[i]] # cut from original
    to<-to[!is.na(xo)]          # omit NA
    xo<-xo[!is.na(xo)]          # omit NA
    gap<-ts[-1]-ts[-length(ts)]
    gapo<-to[-1]-to[-length(to)]
    print( c(as.numeric(GI3[i]), t2[i]-t1[i], max(gapo), length(ts), length(to), length(to)/length(ts)) )
}


# Ca2+
for(i in c(14, 16, 17, 19, 20, 21, 22, 24, 25, 27, 33)){
    ts<-t[t>t1[i] & t<t2[i]]
    to<-tr[tr>t1[i] & tr<t2[i]]  # cut from original
    xo<-xr2[tr>t1[i] & tr<t2[i]] # cut from original
    to<-to[!is.na(xo)]           # omit NA
    xo<-xo[!is.na(xo)]           # omit NA
    gap<-ts[-1]-ts[-length(ts)]
    gapo<-to[-1]-to[-length(to)]
    print( c(as.numeric(GI3[i]), t2[i]-t1[i], max(gapo), length(ts), length(to), length(to)/length(ts)) )
}

# Exclude GI-15.2 (i=17) of Ca2+ (total 320 yr, largest gap 42.9 yr) HIGH GAP RATIO 42.9/320
# Exclude GI-8 (i=24) of Ca2+ (total 1480 yr, largest gap 90.6 yr)   HIGH GAP RATIO
# Exclude GI-7 (i=25) of Ca2+ (total 580 yr, largest gap 217.9 yr)   HIGH GAP RATIO

ii<-c(14, 16, 17, 19, 20, 21, 22, 24, 25, 27, 33) # which(dur>=300)
#ii<-c(16, 19, 20, 21, 22, 24, 25, 33)
ii2<-c(14, 16, 19, 20, 21, 22, 27, 33)
## for look #################################################################################################
postscript(file="ews_data_5yr/NGRIP_time_series_5yr.eps", horizontal=TRUE, encoding="WinAnsi.enc")
par(mfrow=c(2,1))
par(mar=c(0,0,0,0))
par(oma=c(3,0,1,0.5))
par(mai = c(0.3, 1.0, 0.1, 0.1))
plot(-t/1000,d18O.NGRIP,type="l",col=4,xlim=c(60,10),bty = "n",cex.lab=1.3,cex.axis=1.3,ylab=expression(paste(delta^{18},"O (\u2030)")),xaxs="i")
for(i in ii) rect(-t1[i]/1000,-1000,-t2[i]/1000,1000,col="gray80",border=NA)
rect(90140/1000,-1000,90040/1000,1000,col="gray90",border=NA)  # GS-23.1
rect(49600/1000,-1000,49280/1000,1000,col="gray90",border=NA)  # GS-14
lines(-t/1000,d18O.NGRIP,type="l",col=4)
for(i in ii){
    if(i%in%c(16)){
        text((-t1[i]-t2[i])/2000,-33.8,GI2[i])
    } else{
        text((-t1[i]-t2[i])/2000,-32.7,GI2[i])
    }
}
#text((90040+87600)/2000,-32.7,22) # GI-22
#text((49280+48340)/2000,-32.7,13) # GI-22
#text(59.5,-32.7,"GI") # GI-22
#for(i in 1:length(trb)){
# arrows(trb[i]/1000, -35.25, trb[i]/1000, -36, length = 0.075, angle = 30, code = 2, col="magenta", lwd = 2)
#}
mtext("(a)",side=3,adj=-0.1,cex=1.3,line=-0.2,font=2,xpd=TRUE)

plot(-tc/1000,ca.NGRIP,type="l",col=4,xlim=c(60,10),ylim=c(6.3,3.5),bty = "n",cex.lab=1.3,cex.axis=1.3,ylab=expression(paste("log"[10],"(count ml"^{"-1"},")")),xaxs="i")
for(i in ii2) rect(-t1[i]/1000,-1000,-t2[i]/1000,1000,col="gray80",border=NA)
rect(90140/1000,-1000,90040/1000,1000,col="gray90",border=NA)  # GS-23.1
rect(49600/1000,-1000,49280/1000,1000,col="gray90",border=NA)  # GS-14
lines(-t/1000,ca.NGRIP,type="l",col=4)
minor.tick(nx=10, ny=1, tick.ratio=0.6)
for(i in ii2){
    if(i%in%c(16)){
        text((-t1[i]-t2[i])/2000,3.8,GI2[i])
    } else{
        text((-t1[i]-t2[i])/2000,3.6,GI2[i])
    }
}
#text((90040+87600)/2000,0.85,22) # GI-22
#text((49280+48340)/2000,0.85,13) # GI-22
#text(59.5,3.6,"GI") # GI-22
#for(i in 2:length(trb)){
# arrows(trb[i]/1000, 1.08, trb[i]/1000, 1.2, length = 0.075, angle = 30, code = 2, col="magenta", lwd = 2)
#}
mtext("(b)",side=3,adj=-0.1,cex=1.3,line=-0.2,font=2,xpd=TRUE)
mtext("Age (kyr b2k)",side=1,cex=1.3,line=2.5)
dev.off()


## each event ###############################################################################################
k<-4  # which proxy?
set.seed(k)
if(k==1) postscript(file="ews_data_5yr/ews_data_5yr_NGRIP_d18O_1.eps", onefile=FALSE, horizontal=TRUE,encoding="WinAnsi.enc")
if(k==2) postscript(file="ews_data_5yr/ews_data_5yr_NGRIP_d18O_2.eps", onefile=FALSE, horizontal=TRUE,encoding="WinAnsi.enc")
if(k==3) postscript(file="ews_data_5yr/ews_data_5yr_NGRIP_ca_1.eps", onefile=FALSE, horizontal=TRUE,encoding="WinAnsi.enc")
if(k==4) postscript(file="ews_data_5yr/ews_data_5yr_NGRIP_ca_2.eps", onefile=FALSE, horizontal=TRUE,encoding="WinAnsi.enc")

if(k==1) ii<-c(14, 16, 17, 19, 20)      # which(dur>=300)
if(k==2) ii<-c(21, 22, 24, 25, 27, 33)  # which(dur>=300)
if(k==3) ii<-c(14, 16, 19, 20)          # which(dur>=300)
if(k==4) ii<-c(21, 22, 27, 33)          # which(dur>=300)

layout( matrix(1:(4*length(ii)), nrow=4) ) # ncol=6
par(mar=c(0,0,0,0))
par(oma=c(2,2.5,1,1.5))
par(mai = c(0.3, 0.3, 0.1, 0.05))

for(i in ii){

# DO number i
tx<-t[t>t1[i] & t<t2[i]]
if(k==1)  x<-d18O.NGRIP[t>t1[i] & t<t2[i]]
if(k==2)  x<-d18O.NGRIP[t>t1[i] & t<t2[i]]
if(k==3)  x<-ca.NGRIP[t>t1[i] & t<t2[i]]
if(k==4)  x<-ca.NGRIP[t>t1[i] & t<t2[i]]

if(i%in%c(6)){
   t3<-t2[i]+0.15*(t2[i]-t1[i])
} else if(i%in%c(1,10,16,24)){
   t3<-t2[i]+0.25*(t2[i]-t1[i])
} else{
   t3<-t2[i]+0.2*(t2[i]-t1[i])
}
tex<-t[t>t1[i] & t<t3]
if(k==1) xex<-d18O.NGRIP[t>t1[i] & t<t3]
if(k==2) xex<-d18O.NGRIP[t>t1[i] & t<t3]
if(k==3) xex<-ca.NGRIP[t>t1[i] & t<t3]
if(k==4) xex<-ca.NGRIP[t>t1[i] & t<t3]
xl<-c(max(-tex)/1000,min(-tex)/1000) # time interval to plot (positive) 
    
# parameter setting
ns<-1000       # sample size for p values
#smoothing<-"gaussian"
smoothing<-"loess"
bandwidth<-50 # 0-100% for smoothing 
winsize<-50   # 0-100% (default 50%) for rolling statistics
#x<-x[round(0.2*length(x)):length(x)] # remove first %
#tx<-tx[round(0.2*length(tx)):length(tx)] # remove first %
th<-0.1
n<-length(x)                   # data length

bw <- round(n * bandwidth/100) # bandwidth
wd <- round(n * winsize/100)   # window length
span <- bandwidth/100            # bandwidth for loess
if(smoothing=="gaussian") xs <- ksmooth(1:n, x, kernel = "normal", bandwidth = bw, n.points=n, x.points=1:n)$y
if(smoothing=="loess")    xs <- predict(loess(x~tx, span=span, degree=1))             
y<-x-xs                        # residuals
v<-rep(NA,length(y))           # variance
ac<-rep(NA,length(y))          # lag-1 autocorrelation
ac2<-rep(NA,length(y))         # lag-2 autocorrelation
ac3<-rep(NA,length(y))         # lag-3 autocorrelation
#ac4<-rep(NA,length(y))         # lag-4 autocorrelation
lam<-rep(NA,length(y))         # lag-1 lambda
ACF<-matrix(NA,length(y),wd)       # lag-1 lambda
dy<-y[1:n]-c(0,y[-n])          # dy(j)=y(j)-y(j-1). dy[1]=y[1]-0 (unnatural but unused), dy[2]=y[2]-y[1]
ix<-wd:n                   
for(j in ix){
    v[j]<-var(y[(j-wd+1):j])   # window length wd
    Z<-y[(j-wd+1):j]-mean(y[(j-wd+1):j])
    ac[j]<-sum(Z[-1]*Z[-wd])/sum(Z^2)
    ac2[j]<-sum(Z[-c(1,2)]*Z[-c(wd-1,wd)])/sum(Z^2)
    ac3[j]<-sum(Z[-c(1,2,3)]*Z[-c(wd-2,wd-1,wd)])/sum(Z^2)
    #ac4[j]<-sum(Z[-c(1,2,3,4)]*Z[-c(wd-3,wd-2,wd-1,wd)])/sum(Z^2)
    ACF[j,]<-acf(y[(j-wd+1):j],lag.max=wd,plot=FALSE)$acf
    Y<-dy[(j-wd+2):j]-mean(dy[(j-wd+2):j])        # dy[2:wd]     (j=wd), dy[3:(wd+1)]
    X<-y[(j-wd+1):(j-1)]-mean(y[(j-wd+1):(j-1)])  #  y[1:(wd-1)] (j=wd),  y[2:wd]
    lam[j]<-sum(X*Y)/sum(X^2)                   
}
fit1<-lm(v[ix]~tx[ix])
fit2<-lm(ac[ix]~tx[ix])
fit3<-lm(lam[ix]~tx[ix])
trend1<-as.numeric(fit1$coef[2])
trend2<-as.numeric(fit2$coef[2])
trend3<-as.numeric(fit3$coef[2])
pred1<-predict(fit1)
pred2<-predict(fit2)
pred3<-predict(fit3)

sg<-FFTsurrogate(y, n.samples = ns) # surrogate data [instance x time]
sg<-sg-rowMeans(sg)                 # remove rowMeans
vs<-matrix(NA,ns,length(y))         # variance
acs<-matrix(NA,ns,length(y))        # lag-1 autocorrelation
lams<-matrix(NA,ns,length(y))       # lag-1 lambda
dys<-sg[,1:n]-cbind(matrix(0,ns,1),sg[,-n]) # dy(j)=y(j)-y(j-1)
for(j in ix){ 
    vs[,j]<-apply(sg[,(j-wd+1):j],1,var)
    Z<-sg[,(j-wd+1):j]-rowMeans(sg[,(j-wd+1):j])
    acs[,j]<-rowSums(Z[,-1]*Z[,-wd])/rowSums(Z^2)
    Y<-dys[,(j-wd+2):j]-rowMeans(dys[,(j-wd+2):j])        # dy(j)=y(j)-y(j-1)~y(j-1)
    X<-sg[,(j-wd+1):(j-1)]-rowMeans(sg[,(j-wd+1):(j-1)])  
    lams[,j]<-rowSums(X*Y)/rowSums(X^2)
}
vs<-vs[,ix]
acs<-acs[,ix]
lams<-lams[,ix]
vs<-vs-rowMeans(vs)
acs<-acs-rowMeans(acs)
lams<-lams-rowMeans(lams)
tm<-matrix(tx[ix],dim(vs)[1],dim(vs)[2],byrow=TRUE)
tm<-tm-rowMeans(tm)
trend1.surrogate<-rowSums(vs*tm)/rowSums(tm^2)
trend2.surrogate<-rowSums(acs*tm)/rowSums(tm^2)
trend3.surrogate<-rowSums(lams*tm)/rowSums(tm^2)
p1<-length(trend1.surrogate[trend1.surrogate>trend1])/ns
p2<-length(trend2.surrogate[trend2.surrogate>trend2])/ns
p3<-length(trend3.surrogate[trend3.surrogate>trend3])/ns
c(p1,p2,p3)
    

## plot ############################################################
plot(-tex/1000,xex,xlim=xl,type="l",col="gray",bty = "n",cex.lab=1.4,cex.axis=1.4,ylab="")
if(k<=2) lines(-tr[!is.na(xr)]/1000,xr[!is.na(xr)],col="darkorange",lwd=1)
if(k>=3) lines(-tr[!is.na(xr2)]/1000,xr2[!is.na(xr2)],col="darkorange",lwd=1)
lines(-tx/1000,x,col=4)
lines(-tx/1000,xs,col=2,lwd=2)
if(i==min(ii) && k<=2) mtext(expression(paste(delta^{18},"O (\u2030)")), side = 2, line = 2.8, col=1)
if(i==min(ii) && k>=3) mtext(expression(paste("log"[10],"(count ml"^{"-1"},")")), side = 2, line = 2.8, col=1)
title(paste("GI-",GI3[i],sep=""),cex.main=1.4,col=1,xpd = NA,line = 0.5)  #  change!
    
plot(-tx/1000,x-xs,xlim=xl,type="l",col=3,bty = "n",cex.lab=1.4,cex.axis=1.4,ylab="")
if(i==min(ii)) mtext(expression(paste("Residual (\u2030)")), side = 2, line = 2.8, col=1)

yl<-c(min(c(v,pred1),na.rm=TRUE),max(c(v,pred1),na.rm=TRUE))
plot(-tx/1000,v,xlim=xl,ylim=yl,type="l",col=1,bty = "n",cex.lab=1.4,cex.axis=1.4,ylab="")
p<-max(floor(p1*10^3)/10^3,0.001)
if(p<=0.001){
    lines(-tx[ix]/1000,pred1,col=2,lwd=2,lty=1)
    mtext("p<0.001", col=2, side=3, adj=0.07, line=-2, cex=0.8) 
} else if(p<0.05){
    lines(-tx[ix]/1000,pred1,col=2,lwd=2,lty=1)
    mtext(sprintf("p=%.3f",p), col=2, side=3, adj=0.07, line=-2, cex=0.8) 
} else if(p<0.1){
    lines(-tx[ix]/1000,pred1,col=2,lwd=2,lty=2)
    mtext(sprintf("p=%.2f",p), col=2, side=3, adj=0.07, line=-2, cex=0.8) 
} else{
    lines(-tx[ix]/1000,pred1,col=2,lwd=2,lty=3)
    mtext(sprintf("p=%.2f",p), col=2, side=3, adj=0.07, line=-2, cex=0.8) 
}   
if(i==min(ii) && k<=2) mtext(expression(paste("Variance (\u2030"^2,")")), side = 2, line = 2.8, col=1)
if(i==min(ii) && k>=3) mtext("Variance", side = 2, line = 2.8, col=1)
    
yl<-c(min(c(ac,pred2),na.rm=TRUE),max(c(ac,pred2),na.rm=TRUE))
plot(-tx/1000,ac,xlim=xl,ylim=yl,type="l",col=1,bty = "n",cex.lab=1.4,cex.axis=1.4,ylab="")
#lines(tx,ac2,col="dodgerblue")
#lines(tx,ac3,col="darkorange")
#lines(tx,ac4,col="gray")
#abline(h=0,lty=2)
p<-max(floor(p2*10^3)/10^3,0.001)
if(p<=0.001){
    lines(-tx[ix]/1000,pred2,col=2,lwd=2,lty=1)
    mtext("p<0.001", col=2, side=3, adj=0.07, line=-2, cex=0.8) 
} else if(p<0.05){
    lines(-tx[ix]/1000,pred2,col=2,lwd=2,lty=1)
    mtext(sprintf("p=%.3f",p), col=2, side=3, adj=0.07, line=-2, cex=0.8) 
} else if(p<0.1){
    lines(-tx[ix]/1000,pred2,col=2,lwd=2,lty=2)
    mtext(sprintf("p=%.2f",p), col=2, side=3, adj=0.07, line=-2, cex=0.8) 
} else{
    lines(-tx[ix]/1000,pred2,col=2,lwd=2,lty=3)
    mtext(sprintf("p=%.2f",p), col=2, side=3, adj=0.07, line=-2, cex=0.8) 
}
if(i==min(ii)) mtext("Lag-1 autocorrlation", side = 2, line = 2.8, col=1)
    
#yl<-c(min(c(lam,pred3),na.rm=TRUE),max(c(lam,pred3),na.rm=TRUE))
#plot(tx,lam,xlim=xl,ylim=yl,type="l",col=1,bty = "n",cex.lab=1.4,cex.axis=1.4,ylab="")
#p<-max(floor(p3*10^3)/10^3,0.001)
#if(p<=0.001){
#    lines(tx[ix],pred3,col=2,lwd=2,lty=1)
#    text(0.8*min(xl)+0.2*max(xl),max(yl)-0.1,"p<0.001",cex=1.4,col=4)
#} else if(p<th){
#    lines(tx[ix],pred3,col=2,lwd=2,lty=1)
#    text(0.8*min(xl)+0.2*max(xl),max(yl)-0.1,sprintf("p=%.3f",p),cex=1.4,col=4)
#} else{
#    lines(tx[ix],pred3,col=2,lwd=2,lty=2)
#    text(0.8*min(xl)+0.2*max(xl),max(yl)-0.1,sprintf("p=%.3f",p),cex=1.4,col=4)
#}
#if(i==min(ii)) mtext(expression(lambda), side = 2, line = 2.8, col=1)
    
#image(tx,1:wd,ACF,xlim=xl,ylim=c(1,round(wd/2)),zlim=c(-max(ACF[,-1],na.rm=TRUE),max(ACF[,-1],na.rm=TRUE)),col=rev(cl),bty = "n",cex.lab=1.4,cex.axis=1.4,ylab="")
#if(i==min(ii)) mtext("lag", side = 2, line = 2.8, col=1)
mtext("Age (kyr b2k)", side = 1, line = 2.8, col=1)
}

dev.off()




