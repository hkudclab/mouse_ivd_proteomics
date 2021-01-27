PLOT<-function(indj){
	Compart<-factor(Compartments)
	YLIM<-range(exprBR[indj,],na.rm=T)
	agg1<-aggregate(exprBR[indj,],by=list(Compart,Levels),mean,na.rm=T)
	agg2<-aggregate(exprBR[indj,],by=list(Compart,Levels),sd,na.rm=T)
	agg2[["x"]][is.na(agg2[["x"]])]<-0
	x1<-agg1[["x"]][c(1,3)]
	x2<-agg2[["x"]][c(1,3)]

	x3<-agg1[["x"]][c(2,4)]
	x4<-agg2[["x"]][c(2,4)]

	MATP<-signif(ANOVA_P[indj,5],2)
	plot(x1,ylim=YLIM,xlim=c(0.75,2.5),main=paste0(indj,":  ",
			GeneSym[indj]," (F-test P=",MATP,"; FDR q=",signif(FDRall[indj],3),")"),
			ylab="heavy/light ratioss",type='n',axes=F)
	axis(2)
	axis(1,at=seq(2),labels=c("Lumbar","Tail"))

	polygon(c(seq(2),rev(seq(2))),c(x1-x2,rev(x1)+rev(x2)),col="grey")
	lines(c(1,2),x1,lwd=2,lty=2)
	polygon(c(seq(2),rev(seq(2))),c(x3-x4,rev(x3)+rev(x4)),col="#F8766D80")
	lines(c(1,2),x3,lty=2,lwd=2,col="#FE766D")
	LN<-split(exprBR[indj,Compart=="NP"],Levels[Compart=="NP"])
	LO<-split(exprBR[indj,Compart=="AF"],Levels[Compart=="AF"])
	#boxplot(LN,add=T)
	#boxplot(LO,add=T,col=2)
	res<-sapply(seq(2),function(i){
		points(jitter(rep(i,length(LN[[i]])),amount=0.1),LN[[i]],pch=17,cex=2,col="#F8766D")
		points(jitter(rep(i,length(LO[[i]])),amount=0.1),LO[[i]],pch=16,cex=2)
	})
	legend("top",legend=c("AF","NP"),pch=c(16,17),fill=c("grey","#F8766D"))

	MATP2<-signif(ANOVA_P[indj,2:4],2)
	legend("bottomleft",legend=paste0(names(MATP2),": ",MATP2))
}
PLOT_2<-function(indj){
	Compart<-factor(Compartments)
	YLIM<-range(exprBR[indj,],na.rm=T)
	agg1<-aggregate(exprBR[indj,],by=list(Compart,Levels),mean,na.rm=T)
	agg2<-aggregate(exprBR[indj,],by=list(Compart,Levels),sd,na.rm=T)
	agg2[["x"]][is.na(agg2[["x"]])]<-0
	x1<-agg1[["x"]][c(3,4)]
	x2<-agg2[["x"]][c(3,4)]

	x3<-agg1[["x"]][c(1,2)]
	x4<-agg2[["x"]][c(1,2)]

	MATP<-signif(ANOVA_P[indj,5],2)
	plot(x1,ylim=YLIM,xlim=c(0.75,2.5),main=paste0(indj,":  ",
			GeneSym[indj]," (F-test P=",MATP,"; FDR q=",signif(FDRall[indj],3),")"),
			ylab="heavy/light ratioss",type='n',axes=F)
	axis(2)
	axis(1,at=seq(2),labels=c("AF","NP"))

	polygon(c(seq(2),rev(seq(2))),c(x1-x2,rev(x1)+rev(x2)),col="grey")
	lines(c(1,2),x1,lwd=2,lty=2)

	polygon(c(seq(2),rev(seq(2))),c(x3-x4,rev(x3)+rev(x4)),col="#F8766D80")
	lines(c(1,2),x3,lty=2,lwd=2,col="#FE766D")

	LTail<-split(exprBR[indj,Levels=="T"],Compart[Levels=="T"])
	LLumbar<-split(exprBR[indj,Levels=="L"],Compart[Levels=="L"])


	res<-sapply(seq(2),function(i){
		points(jitter(rep(i,length(LTail[[i]])),amount=0.1),LTail[[i]],pch=17,cex=2)
		points(jitter(rep(i,length(LLumbar[[i]])),amount=0.1),LLumbar[[i]],pch=16,cex=2,col="#F8766D")
	})
	legend("top",legend=c("Tail","Lumbar"),pch=c(17,16),fill=c("grey","#F8766D"))

	MATP2<-signif(ANOVA_P[indj,2:4],2)
	legend("bottomleft",legend=paste0(names(MATP2),": ",MATP2))
}
