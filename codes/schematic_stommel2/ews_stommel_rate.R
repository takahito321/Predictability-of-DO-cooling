#library(Hmisc)
#library(misty)
library(nonlinearTseries) # FFTsurrogate

dt<-0.001
T<-300 # 30 300
q<-1
N<-T/dt
t<-(1:N)*dt
a<-seq(0.1,-0.4,length=N)
skip<-100 # 100

set.seed(1)
s<-rep(1.1,N)
dw<-rnorm(N,0,0.03*sqrt(dt)) 
for(i in 1:(N-1)) s[i+1]<-s[i]+q*(abs(s[i])*(1-s[i])+a[i])*dt+dw[i]

i1<-min(which(s>0.95)) # 0.95
i2<-min(which(s<(-0.1) & t>t[i1])) # 0.2
i3<-min(which(s<0.5 & t>t[i1])) # 0.6

plot(t,s,type="l")
abline(h=0.3)
lines(t[seq(i1,i3,skip)],s[seq(i1,i3,skip)],col=4)
length(seq(i1,i3,skip))

## each event #########################################################################
postscript(file="ews_stommel_rate.eps", onefile=FALSE, horizontal=TRUE,encoding="WinAnsi.enc")

ii<-1:2

layout( matrix(1:(4*length(ii)), nrow=4) )
par(mar=c(0,0,0,0))
par(oma=c(2,2.5,1,1.5))
par(mai = c(0.3, 0.3, 0.1, 0.05))

for(i in ii){  # for-loop of each GI

if(i==1) T<-300 # 300 
if(i==2) T<-30  # 30
q<-1
N<-T/dt
t<-(1:N)*dt
a<-seq(0.1,-0.4,length=N)
    
set.seed(i+34) #15
s<-rep(1.1,N)
dw<-rnorm(N,0,0.03*sqrt(dt)) 
for(j in 1:(N-1)) s[j+1]<-s[j]+q*(abs(s[j])*(1-s[j])+a[j])*dt+dw[j]

j1<-min(which(s>0.95)) # 0.95
j2<-N #min(which(s<(-0.35) & t>t[j1])) # 0.2
j3<-min(which(s<0.5 & t>t[j1])) # 0.6

tx<-t[seq(j1,j3,skip)]
x<-s[seq(j1,j3,skip)]

# time series of extended interval xex (gray)
tex<-t[seq(j1,j2,skip)]
xex<-s[seq(j1,j2,skip)]
xeq1<-0.5*(1+sqrt(1+4*a[seq(j1,j2,skip)]))
xeq2<-0.5*(1-sqrt(1+4*a[seq(j1,j2,skip)]))
xeq2[xeq2<0]<-NA
xeq3<-0.5*(1-sqrt(1-4*a[seq(j1,j2,skip)]))
xeq3[xeq3>0]<-NA    
xl<-c(min(tex),max(tex)) # time interval to plot (positive) 
    
# parameter setting
ns<-1000      # sample size for p values
#smoothing<-"gaussian"
smoothing<-"loess"
bandwidth<-50 # 0-100% for smoothing 
winsize<-50   # 0-100% (default 50%) for rolling statistics
#x<-x[round(0.2*length(x)):length(x)]     # remove first % (not used)
#tx<-tx[round(0.2*length(tx)):length(tx)] # remove first % (not used)
th<-0.1
n<-length(x)                   # data length

bw <- round(n * bandwidth/100) # bandwidth     (actual data size)
wd <- round(n * winsize/100)   # window length (actual data size)
span <- bandwidth/100          # bandwidth for loess (ratio)
if(smoothing=="gaussian") xs <- ksmooth(1:n, x, kernel = "normal", bandwidth = bw, n.points=n, x.points=1:n)$y
if(smoothing=="loess")    xs <- predict(loess(x~tx, span=span, degree=1)) # linear
    
y<-x-xs                        # residuals
v<-rep(NA,length(y))           # variance
ac<-rep(NA,length(y))          # lag-1 autocorrelation
ix<-wd:n                   
for(j in ix){
    v[j]<-var(y[(j-wd+1):j])   # window length wd
    Z<-y[(j-wd+1):j]-mean(y[(j-wd+1):j])
    ac[j]<-sum(Z[-1]*Z[-wd])/sum(Z^2)
}
fit1<-lm(v[ix]~tx[ix])
fit2<-lm(ac[ix]~tx[ix])
trend1<-as.numeric(fit1$coef[2])
trend2<-as.numeric(fit2$coef[2])
pred1<-predict(fit1)
pred2<-predict(fit2)

sg<-FFTsurrogate(y, n.samples = ns) # surrogate data [instance x time]
sg<-sg-rowMeans(sg)                 # remove rowMeans
vs<-matrix(NA,ns,length(y))         # variance
acs<-matrix(NA,ns,length(y))        # lag-1 autocorrelation)
for(j in ix){ 
    vs[,j]<-apply(sg[,(j-wd+1):j],1,var)
    Z<-sg[,(j-wd+1):j]-rowMeans(sg[,(j-wd+1):j])
    acs[,j]<-rowSums(Z[,-1]*Z[,-wd])/rowSums(Z^2)
}
vs<-vs[,ix]
acs<-acs[,ix]
vs<-vs-rowMeans(vs)
acs<-acs-rowMeans(acs)
tm<-matrix(tx[ix],dim(vs)[1],dim(vs)[2],byrow=TRUE)
tm<-tm-rowMeans(tm)
trend1.surrogate<-rowSums(vs*tm)/rowSums(tm^2)
trend2.surrogate<-rowSums(acs*tm)/rowSums(tm^2)
p1<-length(trend1.surrogate[trend1.surrogate>trend1])/ns
p2<-length(trend2.surrogate[trend2.surrogate>trend2])/ns
c(p1,p2)
    

## plot ############################################################
plot(tex,xex,xlim=xl,ylim=c(-0.3,1.2),type="l",col="gray",bty = "n",cex.lab=1.7,cex.axis=1.7,ylab="")
lines(tx,x,col=4)
lines(tx,xs,col=2,lwd=2)
lines(tex,xeq1,col=1)
lines(tex,xeq2,col=1,lty=2)
lines(tex,xeq3,col=1)
if(i==1) mtext("x", side = 2, line = 2.8, col=1)
if(i==1) mtext("(a)",side=3,adj=0.04,cex=1.2,line=0,font=2,xpd=TRUE)
if(i==2) mtext("(b)",side=3,adj=0.04,cex=1.2,line=0,font=2,xpd=TRUE)
if(i==1) mtext("Slow parameter change",side=3,adj=0.5,cex=1.2,line=0,font=2,xpd=TRUE)
if(i==2) mtext("Fast parameter change",side=3,adj=0.5,cex=1.2,line=0,font=2,xpd=TRUE)
    
    
plot(tx,x-xs,xlim=xl,type="l",col=3,bty = "n",cex.lab=1.7,cex.axis=1.7,ylab="")
if(i==1) mtext(expression(paste("Residual")), side = 2, line = 2.8, col=1)
if(i==1) mtext("(c)",side=3,adj=0.04,cex=1.2,line=-2,font=2,xpd=TRUE)
if(i==2) mtext("(d)",side=3,adj=0.04,cex=1.2,line=-2,font=2,xpd=TRUE)

yl<-c(min(c(v,pred1),na.rm=TRUE),max(c(v,pred1),na.rm=TRUE))
plot(tx,v,xlim=xl,ylim=yl,type="l",col=1,bty = "n",cex.lab=1.7,cex.axis=1.7,ylab="")
p<-max(floor(p1*10^3)/10^3,0.001)
if(p<=0.001){
    lines(tx[ix],pred1,col=2,lwd=2,lty=1)
    mtext("p<0.001", col=2, side=3, adj=0.3, line=-4, cex=1.2) 
} else if(p<0.05){
    lines(tx[ix],pred1,col=2,lwd=2,lty=1)
    mtext(sprintf("p=%.3f",p), col=2, side=3, adj=0.3, line=-4, cex=1.2) 
} else if(p<0.1){
    lines(tx[ix],pred1,col=2,lwd=2,lty=2)
    mtext(sprintf("p=%.2f",p), col=2, side=3, adj=0.3, line=-4, cex=1.2) 
} else{
    lines(tx[ix],pred1,col=2,lwd=2,lty=3)
    mtext(sprintf("p=%.2f",p), col=2, side=3, adj=0.3, line=-4, cex=1.2) 
}
if(i==1) mtext("Variance", side = 2, line = 2.8, col=1)
if(i==1) mtext("(e)",side=3,adj=0.04,cex=1.2,line=-2,font=2,xpd=TRUE)
if(i==2) mtext("(f)",side=3,adj=0.04,cex=1.2,line=-2,font=2,xpd=TRUE)
    
yl<-c(min(c(ac,pred2),na.rm=TRUE),max(c(ac,pred2),na.rm=TRUE))
plot(tx,ac,xlim=xl,ylim=yl,type="l",col=1,bty = "n",cex.lab=1.7,cex.axis=1.7,ylab="")
p<-max(floor(p2*10^3)/10^3,0.001)
if(p<=0.001){
    lines(tx[ix],pred2,col=2,lwd=2,lty=1)
    mtext("p<0.001", col=2, side=3, adj=0.3, line=-4, cex=1.2) 
} else if(p<0.05){
    lines(tx[ix],pred2,col=2,lwd=2,lty=1)
    mtext(sprintf("p=%.3f",p), col=2, side=3, adj=0.3, line=-4, cex=1.2) 
} else if(p<0.1){
    lines(tx[ix],pred2,col=2,lwd=2,lty=2)
    mtext(sprintf("p=%.2f",p), col=2, side=3, adj=0.3, line=-4, cex=1.2) 
} else{
    lines(tx[ix],pred2,col=2,lwd=2,lty=3)
    mtext(sprintf("p=%.2f",p), col=2, side=3, adj=0.3, line=-4, cex=1.2) 
}
if(i==1) mtext("Lag-1 autocorrlation", side = 2, line = 2.8, col=1)
if(i==1) mtext("(g)",side=3,adj=0.04,cex=1.2,line=-2,font=2,xpd=TRUE)
if(i==2) mtext("(h)",side=3,adj=0.04,cex=1.2,line=-2,font=2,xpd=TRUE)
    
mtext("Time t", side = 1, line = 2.8, col=1)
}

dev.off()

