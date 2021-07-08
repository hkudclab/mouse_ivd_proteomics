pvec<-rep(NA,1e4)
for(i in 1:1e4){
	x<-rnorm(17)
	y<-rnorm(14)+2
	pvec[i]<-t.test(x,y)$p.value
}
table(p.adjust(pvec)<0.05)


pvec<-rep(NA,1e4)
for(i in 1:1e4){
	x<-rnorm(17)
	y<-rnorm(14)
	pvec[i]<-t.test(x,y)$p.value
}
table(p.adjust(pvec)<0.05)
