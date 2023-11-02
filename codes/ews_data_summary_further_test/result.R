result<-numeric(0)
for(case in 1:10){
    filename<-paste("proxy1_DO2_case",case,".dat",sep="")
    tab<-read.table(filename)[,1]
    result<-rbind(result,tab)
}

# i=2 for DO-25: p for (>=15 times over 30 parameter sets), lag-1AC (>=15), variance (>15), and lag-1 AC (>15)
colMeans(result)/500


result<-numeric(0)
for(case in 1:10){
    filename<-paste("proxy1_DO20_case",case,".dat",sep="")
    tab<-read.table(filename)[,1]
    result<-rbind(result,tab)
}

# Case i=20 for DO-12:  p for (>=15 times over 30 parameter sets), lag-1AC (>=15), variance (>15), and lag-1 AC (>15)
colMeans(result)/500
