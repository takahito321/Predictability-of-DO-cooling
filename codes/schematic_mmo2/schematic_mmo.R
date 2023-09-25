dt<-0.001
T<-40
n<-T/dt
t<-(1:n)*dt
k<-0.8
p<-0.225

x<-rep(0.5,n)
y<-rep(-0.5,n)
z<-rep(-0.5,n)
#set.seed(20) #11
dw1<-rnorm(n,0,0.0*sqrt(dt)) 
dw2<-rnorm(n,0,0.0*sqrt(dt))
dw3<-rnorm(n,0,0.0*sqrt(dt))
for(i in 1:(n-1)){
    x[i+1]<-x[i]+50*(abs(x[i])*(1-x[i])+y[i]+p)*dt+dw1[i]
    y[i+1]<-y[i] +0.5*(-x[i]-y[i]+k*(z[i]-y[i]))*dt+dw2[i]
    z[i+1]<-z[i]+0.25*(-x[i]-z[i]+k*(y[i]-z[i]))*dt+dw3[i]
}

# saddle?

J<-matrix(1:9,3,3,byrow=TRUE)
J<-matrix(c((1-2*sqrt(p))*50, 50, 0, -0.5, -(1+k)*0.5, k*0.5, -1*0.25, k*0.25, -(1+k)*0.25),3,3,byrow=TRUE)
eigen(J)

plot(t,x,type="l",xlab="t",ylab="x,y",cex.lab=1.8,cex.axis=1.8,lwd=0.5,ylim=c(-0.2,max(c(x,y))), col=1)

dat<- data.frame(x[t>20], y[t>20], z[t>20])
write.table(dat, file = "a.dat", col.names=F, row.names=F) # -50 used   

dat<-numeric(0)
for(i in 1:100){
zz<--0.5+0.003*i
xs<-0.005*((-400):400)
ys<--abs(xs)*(1-xs)-p
dat<-rbind(dat,data.frame(xs, ys, zz))
}
write.table(dat, file = "b.dat", col.names=F, row.names=F) # -50 used

dat<-numeric(0)
v<-0.07*Re(eigen(J)$vectors[3,])
vx<-c(sqrt(p)-v[1],sqrt(p),sqrt(p)+v[1])
vy<-c(-sqrt(p)-v[2],-sqrt(p),-sqrt(p)+v[2])
vz<-c(-sqrt(p)-v[3],-sqrt(p),-sqrt(p)+v[3])
dat<-rbind(dat,data.frame(vx, vy, vz))
write.table(dat, file = "c.dat", col.names=F, row.names=F) # -50 used

dat<-numeric(0)
dat<-rbind(dat,data.frame(sqrt(p), sqrt(p), sqrt(p)))
write.table(dat, file = "p.dat", col.names=F, row.names=F) # -50 used

postscript(file="schematic_mmo.eps", horizontal=TRUE, encoding="WinAnsi.enc")
par(mar=c(0,0,0,0))
par(oma=c(0,0,1,0))
par(mai = c(0.9, 1, 0.2, 0.2))
plot(t-19,x,xlim=c(0,20),type="l",col=4,ylab="",xlab="",bty = "n",axes=F,lwd=1.5)
axis(side=2, cex.axis=2)
axis(side=1, cex.axis=2)
mtext("x", side = 2, line = 3.5, cex=2.5)
mtext("t", side = 1, line = 3, cex=2.5)
#mtext("(f)",side=3,adj=0.02,cex=2.5,line=-0.5,font=2,xpd=TRUE)
dev.off()

