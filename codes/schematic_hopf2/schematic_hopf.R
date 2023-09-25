dt<-0.001
T<-10
q<-100
th<-0.2
n<-T/dt
t<-(1:n)*dt
a<-seq(0.3,0.2,length=n)
ah<-((1-1/q)/2)^2

set.seed(2) # 2 13
x<-rep(0,n) # -0.3
y<-rep(0,n) # -0.15
dw1<-rnorm(n,0,0.05*sqrt(dt)) #0.05
dw2<-rnorm(n,0,0.0*sqrt(dt))
for(i in 1:(n-1)){
    x[i+1]<-x[i]+q*(abs(x[i])*(1-x[i])+y[i]+a[i])*dt+dw1[i]
    y[i+1]<-y[i]+(-x[i]-y[i])*dt+dw2[i]
}

plot(x,type="l")



x1<-sqrt(a)

dat<- data.frame(a, x, y)
write.table(dat, file = "a.dat", col.names=F, row.names=F) # -50 used   

x0<-sqrt(a[a>ah])
y0<--x0
dat<-data.frame(a[a>ah],x0,y0)
write.table(dat, file = "b.dat", col.names=F, row.names=F) # -50 used   

x1<-sqrt(a[a<ah])
y1<--x1
dat<-data.frame(a[a<ah],x1,y1)
write.table(dat, file = "c.dat", col.names=F, row.names=F) # -50 used   

dat<-numeric(0)
for(i in 1:100){
aa<-0.3-0.001*i
xs<-0.01*((-200):200)
ys<--abs(xs)*(1-xs)-aa
yf<--xs
dat<-rbind(dat,data.frame(xs, ys, yf, aa))
}
write.table(dat, file = "d.dat", col.names=F, row.names=F) # -50 used

postscript(file="schematic_hopf.eps", horizontal=TRUE, encoding="WinAnsi.enc")
par(mar=c(0,0,0,0))
par(oma=c(0,0,1,0))
par(mai = c(0.9, 1, 0.2, 0.2))
plot(a,x,xlim=c(0.3,0.2),type="l",col=4,ylab="",xlab="",bty = "n",axes=F,lwd=1.5)
lines(a[a>ah],x0,col=1,lwd=3)
lines(a[a<ah],x1,col=6,lty=2,lwd=3)
axis(side=2, cex.axis=2)
axis(side=1, cex.axis=2)
mtext("x", side = 2, line = 3.5, cex=2.5)
mtext("p(t)", side = 1, line = 2.8, cex=2.5)
#mtext("(d)",side=3,adj=0.02,cex=2.5,line=-0.5,font=2,xpd=TRUE)
dev.off()
