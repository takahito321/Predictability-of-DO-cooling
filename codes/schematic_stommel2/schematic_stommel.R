dt<-0.001
T<-500
q<-1
N<-T/dt
t<-(1:N)*dt
a<-seq(0.1,-0.4,length=N)

set.seed(1)
x<-rep(1.1,N)
dw<-rnorm(N,0,0.03*sqrt(dt)) 
for(i in 1:(N-1)) x[i+1]<-x[i]+q*(abs(x[i])*(1-x[i])+a[i])*dt+dw[i]

a1<-seq(-0.25,0.1,by=0.001)
x1<-0.5+0.5*sqrt(1+4*a1)
a2<-seq(-0.25,0,by=0.001)
x2<-0.5-0.5*sqrt(1+4*a2)
a3<-seq(-0.5,0,by=0.001)
x3<-0.5-0.5*sqrt(1-4*a3)

postscript(file="schematic_stommel.eps", horizontal=TRUE, encoding="WinAnsi.enc")
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
par(oma=c(0,0,1,0))
par(mai = c(0.9, 1, 0.2, 0.2))
## plot ############################################################
plot(a[(1:N)%%100==0],x[(1:N)%%100==0],type="l",xlim=c(min(a),max(a)),col=4,ylab="",xlab="",bty = "n",axes=F,lwd=1.5)
axis(side=2, cex.axis=2)
axis(side=1, cex.axis=2)
mtext("x", side = 2, line = 3.5, cex=2.5)
mtext("p(t)", side = 1, line = 2.8, cex=2.5)
lines(a1,x1,col=3,lwd=2.5)
lines(a2,x2,col=3,lty=2,lwd=2.5)
lines(a3,x3,col=3,lwd=2.5)
#mtext("(a)",side=3,adj=0.02,cex=2.5,line=-0.5,font=2,xpd=TRUE)
arrows(-0.07, 1.08, -0.1, 1.04, length = 0.2, angle = 30, code = 2, col = 1, lty = 1, lwd = 1)
dev.off()
