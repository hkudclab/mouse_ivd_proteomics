#step11.help.funcs.R


DEP<-function(pheno,phenoName){
	LEV<-levels(factor(pheno))
	print(table(as.character(pheno)))
	print(table(LEV))
	#print(table(as.character(pheno),LEV))
	LMEMB<-split(as.character(pheno),LEV)
	#sapply(LMEMB,function(feno){
		
	#})#seq(nrow(exprBR))
	print(str(LMEMB))
	ANOVP<-sapply(seq(nrow(exprBR)),function(i){
		x<-exprBR[i,]
		#print(x)
		#if(sum(!is.na(x))>10){
			y<-  tryCatch({
				fit1<-aov(x~pheno)
				aP<-summary(fit1)[[1]][[5]][1]
				tres<-TukeyHSD(fit1)
				indMAX<-which.max(abs(tres[[1]][,1]))
				COMPName<-rownames(tres[[1]])[indMAX]
				#print(COMPName)
				tmpi<-c(aP,COMPName,as.vector(tres[[1]][indMAX,]))
				#print(summary(fit1))
				tmpi
			}, warning = function(w) {
				rep(NA,6)
			}, error = function(e) {
				rep(NA,6)
			}, finally = {
				rep(NA,6)
			})
			return(y)
		#}
		#return(NA)
		#rep(NA,6)

	})
	t(ANOVP)
	
}

############################

DEP.t<-function(pheno,phenoName,indS=seq(31)){
	#pheno<-droplevels(pheno[indS])
	LEV<-levels(pheno)
	print(table(pheno))

	TTEST<-sapply(seq(nrow(exprBR)),function(i){
		x<-exprBR[i,indS]
		#print(x)
		#if(sum(!is.na(x))>10){
			y<-  tryCatch({
				fit1<-t.test(x~pheno)
				tmpi<-c(fit1$p.value,
					paste0(rev(levels(pheno)),collapse=" - "),
					diff(fit1$estimate))
				tmpi
			}, warning = function(w) {
				rep(NA,3)
			}, error = function(e) {
				rep(NA,3)
			}, finally = {
				rep(NA,3)
			})
			return(y)
		#}
		#return(NA)
		#rep(NA,6)

	})
	t(TTEST)
}