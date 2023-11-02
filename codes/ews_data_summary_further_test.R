### Instruction ######################################################
# PROGRAM to calculate the probability of observing robust precursor signals
# for 5000 phase-randomized surrogtes
# in Appendix A (Mitsui and Boers, Climate of the Past).
# 
# 0. The outputs of this R-code (like proxy1_DO**_case*.dat) are generated in the directry ews_data_summary_further_test 
# 1. Set index 'i' below to specify which DO is examined: i=2 for DO-25 or i=20 for DO-12. Those two are considered in the paper.  
# 2. Set index 'case' below from 1 to 10. Then run this R-code for each case.
# Finally go to directry ews_data_summary_further_test. Then run result.R to obtain the probability of observing robust precursor signals for variance (>=15 times over 30 parameter sets), lag-1AC (>=15), variance (>15), and lag-1 AC (>15)
case<-1

library(misty)
library(nonlinearTseries) # FFTsurrogate
library(smoother)

tab1<-read.xlsx("Rasmussen_et_al_2014_QSR_Table_2.xlsx",sheet=1,skip=28)
tab1<-tab1[,1:8] # 1:8!
i1<-c(9,13,15,17,19,21,23,25,29,33,35,37,39,43,53,55,57,61,63,67,69,71,73,75,79,85,87,97,99,103,105,109,111) # GI-start
i2<-c(2,12,14,16,18,20,22,24,26,30,34,36,38,40,44,54,56,58,62,64,68,70,72,74,76,80,86,88,98,100,104,106,110) # GI-end
t1<-tab1[i1,3] # GI-start (positive)
t2<-tab1[i2,3] # GI-end   (positive)
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

# Data for seierstad.xlsx
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
k<-1         # which proxy?
i<-20         # which DO? DO-25 (i=2), DO-12 (i=20)
#case<-1     # sub-division from 1 to 20: input this at the top
ii<-1:500     # how many experiments?
set.seed(case)
event<-i

# parameter setting
ns<-1000     # sample size for p values
#smoothing<-"gaussian"
smoothing<-"loess"
bandwidth<-c(20,30,40,50,60,70) 
winsize<-c(20,30,40,50,60)

trend1<-array(NA,dim=c(max(ii),length(winsize),length(bandwidth)))
trend2<-array(NA,dim=c(max(ii),length(winsize),length(bandwidth)))
p1<-array(NA,dim=c(max(ii),length(winsize),length(bandwidth)))
p2<-array(NA,dim=c(max(ii),length(winsize),length(bandwidth)))
    
# DO number
tx<-t[t>t1[i] & t<t2[i]]
if(k==1)  xd<-d18O.NGRIP[t>t1[i] & t<t2[i]]
if(k==2)  xd<-d18O.NGRIP[t>t1[i] & t<t2[i]]
if(k==3)  xd<-d18O.GRIP[t>t1[i] & t<t2[i]]
if(k==4)  xd<-d18O.GRIP[t>t1[i] & t<t2[i]]
if(k==5)  xd<-d18O.GISP2[t>t1[i] & t<t2[i]]
if(k==6)  xd<-d18O.GISP2[t>t1[i] & t<t2[i]]
if(k==7)  xd<-ca.NGRIP[t>t1[i] & t<t2[i]]
if(k==8)  xd<-ca.NGRIP[t>t1[i] & t<t2[i]]
if(k==9)  xd<-ca.GRIP[t>t1[i] & t<t2[i]]
if(k==10) xd<-ca.GRIP[t>t1[i] & t<t2[i]]
if(k==11) xd<-ca.GISP2[t>t1[i] & t<t2[i]]
if(k==12) xd<-ca.GISP2[t>t1[i] & t<t2[i]]
n<-length(xd)                      # data length


## main ####################################################################################
for(i in ii){  # experiment loop, generating different surrogates of the original ice reocrd

x<-as.numeric( FFTsurrogate(xd, n.samples = 1) ) # random surrogate data [instance x time]

for(l in 1:length(winsize)){
for(m in 1:length(bandwidth)){

bw <- round(n * bandwidth[m]/100)  # bandwidth
wd <- round(n * winsize[l]/100)    # window length
span <- bandwidth[m]/100           # bandwidth for loess
if(smoothing=="gaussian") xs <- ksmooth(1:n, x, kernel = "normal", bandwidth = bw, n.points=n, x.points=1:n)$y
if(smoothing=="loess")    xs<-predict(loess(x~tx, span=span, degree=1))             
y<-x-xs                            # residuals

v<-rep(NA,length(y))               # variance
ac<-rep(NA,length(y))              # lag-1 autocorrelation
ix<-wd:n                   
for(j in ix){
    v[j]<-var(y[(j-wd+1):j])       # window length wd
    Z<-y[(j-wd+1):j]-mean(y[(j-wd+1):j])
    ac[j]<-sum(Z[-1]*Z[-wd])/sum(Z^2)       
}
fit1<-lm(v[ix]~tx[ix])
fit2<-lm(ac[ix]~tx[ix])
trend1[i,l,m]<-as.numeric(fit1$coef[2])
trend2[i,l,m]<-as.numeric(fit2$coef[2])
    
sg<-FFTsurrogate(y, n.samples = ns) # surrogate data [instance x time]
sg<-sg-rowMeans(sg)                 # remove rowMeans
vs<-matrix(NA,ns,length(y))         # variance
acs<-matrix(NA,ns,length(y))        # lag-1 autocorrelation
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
p1[i,l,m]<-length(trend1.surrogate[trend1.surrogate>trend1[i,l,m]])/ns
p2[i,l,m]<-length(trend2.surrogate[trend2.surrogate>trend2[i,l,m]])/ns
}
}
}

# p1[i,l,m] = p-values for smoothing and window (l,m) of surrogate record i
# robust EW signal for surrogate i if sum(p1[i,,]<0.05)>=15


count1<-numeric(dim(p1)[1])
for(i in 1:length(count1)){
   count1[i]<-sum(p1[i,,]<0.05) # variance 
}
hist(count1,breaks=100,xlab="Probability of p<0.05",main="Variance")
count2<-numeric(dim(p2)[1])
for(i in 1:length(count1)){
   count2[i]<-sum(p2[i,,]<0.05) # lag-1 AC
}


filename<-paste("ews_data_summary_further_test/proxy",k,"_DO",event,"_case",case,".dat",sep="")
dat<-c(sum(count1>=15),sum(count2>=15),sum(count1>15),sum(count2>15),max(ii)) 
write.table(dat, file = filename, col.names=F, row.names=F)






#postscript(file=filename, horizontal=TRUE, encoding="WinAnsi.enc")
#layout( matrix(1:2, ncol=2) ) # ncol=6
#par(mar=c(0,0,0,0))
#par(oma=c(1,1,1,1))
#par(mai = c(0.9, 0.8, 0.15, 0.1))
#count1<-numeric(dim(p1)[1])
#for(i in 1:length(count1)){
#   count1[i]<-sum(p1[i,,]<0.05) #/length(p1[i,,])
#}
#hist(count1,breaks=100,xlab="Probability of p<0.05",main="Variance")
#count2<-numeric(dim(p2)[1])
#for(i in 1:length(count1)){
#   count2[i]<-sum(p2[i,,]<0.05) #/length(p2[i,,])
#}
#hist(count2,breaks=100,xlab="Probability of p<0.05",main="Lag-1 autocorrelation")
#dev.off()

#filename<-paste("ews_data_summary_further_test/proxy",k,"_DO",event,".dat",sep="")
#dat<-c(sum(p1<0.05)/sum(length(p1)),
#sum(p2<0.05)/sum(length(p2)),
#(sum(p1<0.05)+sum(p2<0.05))/sum(length(p1)+length(p2)))
#write.table(dat, file = filename, col.names=F, row.names=F) 
