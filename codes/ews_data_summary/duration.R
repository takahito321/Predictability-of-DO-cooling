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

t1.o<-t1
t2.o<-t2

t1<-t1-2*4*(u1=="±4")-2*20*(u1=="a")-2*60*(u1=="b")-2*200*(u1=="c")-2*100*(u1=="d")-2*40*(u1=="f")
t2<-t2-2*4*(u2=="±4")+2*20*(u2=="a")+2*60*(u2=="b")+2*200*(u2=="c")+2*100*(u2=="d")+2*40*(u2=="f")

t1<--rev(t1)   # reverse
t2<--rev(t2)
dur<-t2-t1

t1.o<--rev(t1.o)   # reverse
t2.o<--rev(t2.o)
dur.o<-t2.o-t1.o

ii<-c(2, 3, 4, 6, 8, 9, 10, 16, 19, 20, 24, 33) # NGRIP d18O

dur[ii]   # duration after 2*sigma removing uncertainty ranges
dur.o[ii] # duration original
