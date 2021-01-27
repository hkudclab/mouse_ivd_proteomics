library(scales)
library(stringr)
library(mice)
###################################################
load("../resources/Matrisome.mouse.RData")
load("Matt.n91.RData")
load("../resources/gencode.vM11.annotation.gene.gtf.RData")

GNames<-MATT.dat$Genenames
levels(GNames)[which(levels(GNames)=="1-Mar")]<-"Mtarc1"
levels(GNames)[grep("-Sep",levels(GNames))]<-paste0("Sept",
	gsub("-.*","",levels(GNames)[grep("-Sep",levels(GNames))]))

table(GNames%in%annot.g[,4])
setdiff(GNames,annot.g[,4])


GeneSym<-gsub(";.*","",GNames)
table(GeneSym%in%annot.g[,4])
setdiff(GeneSym,annot.g[,4])


LgMatrisome2<-LgMatrisome[c(1,4,2,5,7,8)]

NumProt<-colSums(!is.na(expr))
BioRep<- gsub("_.*","",HEADER[,5])
HEADERbioRep<-HEADER[match(unique(BioRep),BioRep),]
mouseIDs<-factor(paste0(HEADERbioRep[,4],"_",gsub("[TL].?.?","",HEADERbioRep[,3])))
levels(mouseIDs)<-paste0("mouse",seq(10))

BRID<-paste0("BR_",gsub("_.*","",HEADERbioRep[,5]))
HEADER_BR<-data.frame(HEADERbioRep,
		Compartment=gsub("[0-9]","",HEADERbioRep[,3]),
		BRID,mouseIDs)[order(HEADERbioRep[,4]),]

exprBR<-sapply(as.character(HEADER_BR$BRID),function(x){
	IDi<-gsub("BR_","",x)
	rowMeans(expr[,which(gsub("X|_.*","",colnames(expr))==IDi),drop=F],na.rm=T)
})

names(LgMatrisome2)[1:6]<-c("Collagens","Proteoglycans","Glycoproteins",
		"ECM affiliated","ECM regulators","Secreted factors")
LgMatrisome3<-LgMatrisome2

LgMatrisome3[["Core"]]<-LgMatrisome$corematrisome_mm
LgMatrisome3[["non-Core"]]<-setdiff(LgMatrisome$matrisome_mm_masterlist,LgMatrisome$corematrisome_mm)
LgMatrisome3[["Matrisome"]]<-LgMatrisome$matrisome_mm_masterlist
LgMatrisome3[["Others"]]<-setdiff(GeneSym,unlist(LgMatrisome2))


##############################################
Compartments<-HEADER_BR[,2]
Gender<-HEADER_BR[,4]
Levels<-factor(substr(HEADER_BR[,3],1,1))
Litter<-factor(substr(HEADER_BR[,3],4,4))
##############################################
COL_AFNP<-hue_pal()(2)[as.integer(Compartments)]
PCH_gender<-c(16,17)[as.integer(Gender)]
COL_levels<-seq(3)[as.integer(Levels)]
COL_Litter<-seq(5)[as.integer(Litter)]
##############################################
source("step11.help.funcs.R")
library(scales)

ANOVA_P<-matrix(NA,nrow(exprBR),5)
colnames(ANOVA_P)<-c("(Intercept)", "CompartmentsNP", "LevelsT",               
	"CompartmentsNP:LevelsT","FTEST")
ANOVA_COEF<-ANOVA_P

ANOVA_full<-t(sapply(seq(nrow(exprBR)),function(i){
	x<-exprBR[i,]
	y<-tryCatch({
			av1<-lm(x~Compartments*Levels)
			1
		}, warning = function(w) {
			rep(NA,1)
		}, error = function(e) {
			rep(NA,1)
		}, finally = {
			rep(NA,1)
	})
	if(!is.na(y)){
		sum1<-summary(av1)
		PVAL<-1-pf(sum1$fstatistic[1],sum1$fstatistic[2],sum1$fstatistic[3])
		coeffi<-coefficients(sum1)[,"Pr(>|t|)"]
		est<-coefficients(sum1)[,"Estimate"]
		indi<-match(names(coeffi),colnames(ANOVA_P))
		z<-rep(NA,5);names(z)<-colnames(ANOVA_P)
		z[indi]<-as.numeric(coeffi)
		z[5]<-PVAL
		ANOVA_COEF[i,names(est)]<<-est
		ANOVA_P[i,]<<-as.numeric(z)
	}
}))

colnames(ANOVA_COEF)<-gsub("Compartments","Comp",colnames(ANOVA_COEF))
colnames(ANOVA_P)<-colnames(ANOVA_COEF)
#############################################

FTESTP<-ANOVA_P[,5]
indnonNA<-which(!is.na(FTESTP))
FTEST_FDR<-p.adjust(FTESTP[indnonNA])

FDRall<-rep(1.5,nrow(exprBR))
FDRall[indnonNA]<-FTEST_FDR

source("step11.help.funcs-PLOT.R")
PLOT(1523)

#############################################
pdf("ANOVA/full-model-interactions.pdf",width=28,height=16)
	table(p.adjust(ANOVA_P[,4])<0.05,FDRall<0.05)
	table((ANOVA_P[,4])<0.05,FDRall<0.05)

	FLAG1<-(ANOVA_P[,4])<0.01
	indINT1<-which(FLAG1&FDRall<0.05)

	par(mfrow=c(1,1))
	plot(0,type='n',axes=F,xlab="",ylab="")
	text(1,0,paste0("Below are DE-proteins with interactions in a full two-way ANOVA model
		(total ",length(indINT1)," proteins)"))

	par(mfrow=c(4,6))
	for(k in indINT1[order(FDRall[indINT1])]){
		PLOT(k)	
		PLOT_2(k)
		boxplot(exprBR[k,]~HEADER_BR$Compartment)
	}
#############
	indINT1<-which(FLAG1&FDRall>0.05)

	par(mfrow=c(1,1))
	plot(0,type='n',axes=F,xlab="",ylab="")
	text(1,0,paste0("Below are non-DE-proteins with interactions in a full two-way ANOVA model
		(total ",length(indINT1)," proteins)"))

	par(mfrow=c(4,6))
	for(k in indINT1[order(FDRall[indINT1])]){
		PLOT(k)	
		PLOT_2(k)
		boxplot(exprBR[k,]~HEADER_BR$Compartment)
	}
#############
	indINT1<-which(p.adjust(ANOVA_P[,4])<0.05)

	par(mfrow=c(1,1))
	plot(0,type='n',axes=F,xlab="",ylab="")
	text(1,0,paste0("Below are DEP/non-DEPs with STRICT interactions in a full two-way ANOVA model
		(total ",length(indINT1)," proteins)"))

	par(mfrow=c(4,6))
	for(k in indINT1[order(ANOVA_P[,4][indINT1])]){
		PLOT(k)	
		PLOT_2(k)
		boxplot(exprBR[k,]~HEADER_BR$Compartment)
	}

	table((ANOVA_P[,4])<0.05,FDRall<0.05)
#############
	indINT1<-which(ANOVA_P[,4]<0.05&FDRall<0.05)

dev.off()
#############################################
UNIQ<-c("Lamc1", "Comp", "Thbs1", "Cilp", "Fgg", "Atp5b", "Acta1", "Myl1", "Serpina1b", "Ak1", "Acat1", "Map4", "Got2", "Col1a1", "Col1a2", "Col5a2", "Col15a1", "Vcan", "Fn1", "Pcolce", "Lgals3", "Dstn", "Hnrnpm", "Akr1a1", "Sod3", "Myoc")
indCandidate<-match(UNIQ,GeneSym)

#############################################
pdf("ANOVA/full-model-interactions-candidate-list.pdf",width=28,height=16)
	
	par(mfrow=c(1,1))
	plot(0,type='n',axes=F,xlab="",ylab="")
	text(1,0,paste0("Below are a list of candidate proteins with their
		 interactions in a full two-way ANOVA model
		(total ",length(indCandidate)," proteins)"))

	par(mfrow=c(4,6))
	for(k in indCandidate[order(ANOVA_P[indCandidate,4])]){
		PLOT(k)	
		PLOT_2(k)
		boxplot(exprBR[k,]~HEADER_BR$Compartment)
	}
dev.off()

#############


