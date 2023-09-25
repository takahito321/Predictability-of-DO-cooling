dt<-0.001
T<-20
a<-0.26
q<-100 # 100
beta<-1.1
n<-T/dt
t<-(1:n)*dt

set.seed(4) # 4
x<-rep(0.9,n)
y<-rep(-0.35,n)
dw1<-rnorm(n,0,0.4*sqrt(dt)) # 0.3
dw2<-rnorm(n,0,0.0*sqrt(dt))
for(i in 1:(n-1)){
    x[i+1]<-x[i]+q*(abs(x[i])*(1-x[i])+y[i]+a)*dt+dw1[i]
    y[i+1]<-y[i]+(-x[i]-y[i])*dt+dw2[i]
}

plot(t,x,type="l")

xs<-0.03*((-100):100)
ys<--abs(xs)*(1-xs)-a
yf<--xs

#postscript(file="schematic_fhn.eps", horizontal=TRUE, encoding="WinAnsi.enc")
par(mar=c(0,0,0,0))
par(oma=c(0,0,1,0))
par(mai = c(0.9, 1, 0.2, 0.2))
plot(t,x,xlim=c(0,max(T)),type="l",col=4,ylab="",xlab="",bty = "n",axes=F)
axis(side=2, cex.axis=2)
axis(side=1, cex.axis=2)
mtext("x", side = 2, line = 3.5, cex=2.5)
mtext("p", side = 1, line = 2.8, cex=2.5)
mtext("(b)",side=3,adj=0.02,cex=2.5,line=-0.5,font=2,xpd=TRUE)
#dev.off()

#postscript(file="fhn_2d.eps", horizontal=TRUE, encoding="WinAnsi.enc")
par(mar=c(0,0,0,0))
par(oma=c(0,0,1,1))
par(mai = c(0.9, 1, 0.2, 0.2))
plot(y,x,type="l",xlim=c(-0.55,-0.2),ylim=c(-0.3,1.1),col=4,ylab="",xlab="",bty = "n",xaxs="i",axes=F)
lines(ys,xs,lwd=2,col=3)
lines(yf,xs,lwd=2,col=6,lty=2)
axis(side=2, cex.axis=2)
axis(side=1, cex.axis=2)
mtext("x", side = 2, line = 3.5, cex=2.5)
mtext("y", side = 1, line = 2.8, cex=2.5)
mtext("(c)",side=3,adj=0.02,cex=2.5,line=-0.5,font=2,xpd=TRUE)
arrows(-0.36, 0.99, -0.39, 0.95, length = 0.2, angle = 30, code = 2, col = 1, lty = 1, lwd = 1)
arrows(-0.38,-0.23 , -0.35, -0.21, length = 0.2, angle = 30, code = 2, col = 1, lty = 1, lwd = 1)
#dev.off()
