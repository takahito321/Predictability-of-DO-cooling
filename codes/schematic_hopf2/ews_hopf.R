library(Hmisc)
library(misty)
library(Kendall)
library(nonlinearTseries) # FFTsurrogate
library(smoother)
library(RColorBrewer)
cl<-brewer.pal(11,"RdBu")

dt<-0.001
T<-10
q<-100
th<-0.2
N<-T/dt
t<-(1:N)*dt
a<-seq(0.3,0.2,length=N)
ah<-((1-1/q)/2)^2

set.seed(2) # 2 13
s<-rep(0,N) # -0.3
y<-rep(0,N) # -0.15
dw1<-rnorm(N,0,0.05*sqrt(dt)) #0.05
dw2<-rnorm(N,0,0.0*sqrt(dt))
for(i in 1:(N-1)){
    s[i+1]<-s[i]+q*(abs(s[i])*(1-s[i])+y[i]+a[i])*dt+dw1[i]
    y[i+1]<-y[i]+(-s[i]-y[i])*dt+dw2[i]
}

i1<-min(which(s>1.1))+100
i2<-min(which(s<0.3 & t>t[i1]))
i3<-max(which(s>0.5 & t<t[i2]))

plot(t,s,type="l")
abline(h=0.3)
lines(t[i1:i3],s[i1:i3],col=4)



## each event ###############################################################################################
postscript(file="ews_hopf.eps", onefile=FALSE, horizontal=TRUE,encoding="WinAnsi.enc")

ii<-1:6

layout( matrix(1:(4*length(ii)), nrow=4) )
par(mar=c(0,0,0,0))
par(oma=c(2,2.5,1,1.5))
par(mai = c(0.3, 0.3, 0.1, 0.05))

for(i in ii){  # for-loop of each GI

dt<-0.001
T<-10
q<-100
N<-T/dt
t<-(1:N)*dt
a<-seq(0.3,0.2,length=N)

set.seed(i) #13
s<-rep(0,N)
y<-rep(0,N)
dw1<-rnorm(N,0,0.05*sqrt(dt)) # 0.2
dw2<-rnorm(N,0,0.0*sqrt(dt))
for(j in 1:(N-1)){
    s[j+1]<-s[j]+q*(abs(s[j])*(1-s[j])+y[j]+a[j])*dt+dw1[j]
    y[j+1]<-y[j]+(-s[j]-y[j])*dt+dw2[j]
}

j1<-min(which(s>1.1))+100
j2<-min(which(s<0.3 & t>t[j1]))
j3<-max(which(s>0.5 & t<t[j2]))

#plot(t[j1:j3],s[j1:j3],type="l")

# DO number i
tx<-t[seq(j1,j3,by=10)]
x<-s[seq(j1,j3,by=10)]

# time series of extended interval xex (gray)
tex<-t[seq(j1,(j2+200),by=10)]
xex<-s[seq(j1,(j2+200),by=10)]
xl<-c(min(tex),max(tex)) # time interval to plot (positive) 
    
# parameter setting
ns<-1000      # sample size for p values
#smoothing<-"gaussian"
smoothing<-"loess"
bandwidth<-20 # 0-100% for smoothing 
winsize<-30   # 0-100% (default 50%) for rolling statistics
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
ac2<-rep(NA,length(y))         # lag-2 autocorrelation
ac3<-rep(NA,length(y))         # lag-3 autocorrelation
#ac4<-rep(NA,length(y))        # lag-4 autocorrelation
lam<-rep(NA,length(y))         # lag-1 lambda
ACF<-matrix(NA,length(y),wd)   # lag-1 lambda
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
plot(tex,xex,xlim=xl,type="l",col="gray",bty = "n",cex.lab=1.4,cex.axis=1.4,ylab="")
lines(tx,x,col=4)
lines(tx,xs,col=2,lwd=2)
if(i==1) mtext("x", side = 2, line = 2.8, col=1)
    
plot(tx,x-xs,xlim=xl,type="l",col=3,bty = "n",cex.lab=1.4,cex.axis=1.4,ylab="")
if(i==1) mtext(expression(paste("Residual")), side = 2, line = 2.8, col=1)

yl<-c(min(c(v,pred1),na.rm=TRUE),max(c(v,pred1),na.rm=TRUE))
plot(tx,v,xlim=xl,ylim=yl,type="l",col=1,bty = "n",cex.lab=1.4,cex.axis=1.4,ylab="")
p<-max(floor(p1*10^3)/10^3,0.001)
if(p<=0.001){
    lines(tx[ix],pred1,col=2,lwd=2,lty=1)
    mtext("p<0.001", col=2, side=3, adj=0.07, line=-2, cex=0.8) 
} else if(p<0.05){
    lines(tx[ix],pred1,col=2,lwd=2,lty=1)
    mtext(sprintf("p=%.3f",p), col=2, side=3, adj=0.07, line=-2, cex=0.8) 
} else if(p<0.1){
    lines(tx[ix],pred1,col=2,lwd=2,lty=2)
    mtext(sprintf("p=%.2f",p), col=2, side=3, adj=0.07, line=-2, cex=0.8) 
} else{
    lines(tx[ix],pred1,col=2,lwd=2,lty=3)
    mtext(sprintf("p=%.2f",p), col=2, side=3, adj=0.07, line=-2, cex=0.8) 
}
if(i==1) mtext("Variance", side = 2, line = 2.8, col=1)
    
yl<-c(min(c(ac,pred2),na.rm=TRUE),max(c(ac,pred2),na.rm=TRUE))
plot(tx,ac,xlim=xl,ylim=yl,type="l",col=1,bty = "n",cex.lab=1.4,cex.axis=1.4,ylab="")
p<-max(floor(p2*10^3)/10^3,0.001)
if(p<=0.001){
    lines(tx[ix],pred2,col=2,lwd=2,lty=1)
    mtext("p<0.001", col=2, side=3, adj=0.07, line=-2, cex=0.8) 
} else if(p<0.05){
    lines(tx[ix],pred2,col=2,lwd=2,lty=1)
    mtext(sprintf("p=%.3f",p), col=2, side=3, adj=0.07, line=-2, cex=0.8) 
} else if(p<0.1){
    lines(tx[ix],pred2,col=2,lwd=2,lty=2)
    mtext(sprintf("p=%.2f",p), col=2, side=3, adj=0.07, line=-2, cex=0.8) 
} else{
    lines(tx[ix],pred2,col=2,lwd=2,lty=3)
    mtext(sprintf("p=%.2f",p), col=2, side=3, adj=0.07, line=-2, cex=0.8) 
}
if(i==1) mtext("Lag-1 autocorrlation", side = 2, line = 2.8, col=1)
    
mtext("Time t", side = 1, line = 2.8, col=1)
}

dev.off()

