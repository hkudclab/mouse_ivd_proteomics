library(scales)
library(stringr)
library(mice)
###################################################
load("../resources/matrisome-mit-mouse/Matrisome.mouse.RData")
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

LPheno<-list(Compartments=factor(Compartments),
		Gender=factor(Gender),
		Levels=Levels,
		Litter=Litter)
LDEPi<-rep(list())
LDEPTi<-LDEPi
for(i in 1:4)
	LDEPi[[i]]<-DEP(LPheno[[i]],names(LPheno)[i])
for(i in 1:3)
	LDEPTi[[i]]<-DEP.t(LPheno[[i]],names(LPheno)[i])

indNP<-which(Compartments=="NP")
indAF<-which(Compartments=="AF")
indLNP<-which(HEADER_BR$Compartment=="LNP")
indLAF<-which(HEADER_BR$Compartment=="LAF")
indTNP<-which(HEADER_BR$Compartment=="TNP")
indTAF<-which(HEADER_BR$Compartment=="TAF")


LIND<-rep(list(seq(31)),3)
LIND[[4]]<-indNP
LIND[[5]]<-indAF
LIND[[6]]<-indLNP
LIND[[7]]<-indLAF
LIND[[8]]<-indTNP
LIND[[9]]<-indTAF

LPHENO<-LPheno[1:3]
LPHENO[[4]]<-Levels[indNP]
levels(LPHENO[[4]])<-paste0("NP_",levels(LPHENO[[4]]))
LPHENO[[5]]<-Levels[indAF]
levels(LPHENO[[5]])<-paste0("AF_",levels(LPHENO[[5]]))

LPHENO[[6]]<-factor(paste0("LNP_",Gender[indLNP]))
LPHENO[[7]]<-factor(paste0("LAF_",Gender[indLAF]))
LPHENO[[8]]<-factor(paste0("TNP_",Gender[indTNP]))
LPHENO[[9]]<-factor(paste0("TAF_",Gender[indTAF]))

LDEPTi[[4]]<-DEP.t(LPHENO[[4]],"Levels at NP",indNP)
LDEPTi[[5]]<-DEP.t(LPHENO[[5]],"Levels at AF",indAF)

LDEPTi[[6]]<-DEP.t(LPHENO[[6]],"Genders at LNP",indLNP)
LDEPTi[[7]]<-DEP.t(LPHENO[[7]],"Genders at LAF",indLAF)
LDEPTi[[8]]<-DEP.t(LPHENO[[8]],"Genders at TNP",indTNP)
LDEPTi[[9]]<-DEP.t(LPHENO[[9]],"Genders at TAF",indTAF)

#############################################
library(gplots)

pdf("DEP/DEP.BioRep31.Jan26.ttest.pdf",width=9)

	for(i in 1:9){
		DEPi<-LDEPTi[[i]]
		indNNA<-which(!is.na(as.numeric(DEPi[,1])))

		INDSi<-LIND[[i]]
		pheno<-droplevels(LPHENO[[i]])
		tabFeno<-table(pheno)
		LEV<-levels(pheno)
		Num1<-rowSums(!is.na(exprBR[,INDSi][,pheno==LEV[1]]))
		Num2<-rowSums(!is.na(exprBR[,INDSi][,pheno==LEV[2]]))

		par(mfrow=c(1,1))
		L1<-list(GeneSym[Num1>=3|Num1>(tabFeno[LEV[1]]/2)],
			GeneSym[Num2>=3|Num2>(tabFeno[LEV[2]]/2)])
		L1<-list(GeneSym[Num1>0],
			GeneSym[Num2>0])

		names(L1)<-LEV
		venn(L1)
		L2<-list(GeneSym[Num1>0],
			GeneSym[Num2>0],
			GeneSym[indNNA])
		names(L2)<-c(LEV,"tested")
		venn(L2)

		
		indG1<-which(Num1>= max(Num1)/2 & Num2==0)
		indG2<-which(Num2>= max(Num2)/2 & Num1==0)
		G1<-GeneSym[indG1]
		G2<-GeneSym[indG2]
######################################
######################################
		if(0){
		if(length(c(indG1,indG2))>0&0){
			fpkm0<-exprBR[c(indG1,indG2),INDSi]
			dat.G12<-data.frame(GeneSym[c(indG1,indG2)],
				rep(LEV,c(length(indG1),length(indG2))),fpkm0)
			###################
			feno<-droplevels(data.frame(BRID=colnames(fpkm0),
					site=substr(HEADER_BR[INDSi,6],1,1),
					HEADER_BR[INDSi,c(2,6,4)]))
			colnames(feno)[5]<-"Gender"
			MAThead<-matrix("",ncol(feno),2)
			MAThead[1,]<-c("Proteins","geneGroup")
			header<-data.frame(MAThead,t(feno))
			header[,1]<- colnames(feno)

			FNOUT<-paste0("DEP/exclusive.",i,".DEP.",paste0(paste0(LEV,".n",c(length(indG1),length(indG2))),collapse=".vs."),".txt")
			cat("FNOUT:",FNOUT,"\n")
			write.table(header,file=FNOUT,quote=F,row.names=F,col.names=F,sep="\t")
			write.table(dat.G12,file=FNOUT,quote=F,row.names=F,col.names=F,sep="\t",append=T)
		}
		}
######################################
######################################
		indNNA12<-which(!is.na(as.numeric(DEPi[,1]))&Num1>2&Num2>2)

		DEPi2<-DEPi[indNNA12,]
		PVALi<-as.numeric(DEPi2[,1])
		DIFFi<-as.numeric(DEPi2[,3])
		FDRi<-p.adjust(PVALi)
		indSig<-which(FDRi<0.05)

		SIGNi<-sign(DIFFi[indSig])
		KOLi<-rev(hue_pal()(2))[(SIGNi+3)/2]

		DATOUTi<-data.frame(GeneSym[indNNA12],DEPi2,PVALi,DIFFi,FDRi,DEP=FDRi<0.05)[order(PVALi),]

		layout(t(matrix(seq(2))),widths=c(6,3))
		plot(DIFFi,-log10(PVALi),main=paste0("No.",i,":",names(LPHENO)[i]),
			pch=16,col="grey",xlab=DEPi2[1,2],ylab="-log10(P)")
		if(length(indSig)>0){
			points(DIFFi[indSig],-log10(PVALi)[indSig],
				pch=16,col=KOLi)
			text(DIFFi[indSig],-log10(PVALi)[indSig],
				GeneSym[indNNA12][indSig],xpd=T)
#################################
			gUp<-GeneSym[indNNA12][indSig][SIGNi>0]
			gDown<-GeneSym[indNNA12][indSig][SIGNi<0]


			source("../my-func-lib/CAT_DEGs3.R")
			CAT_DEG2(gUp,gDown,BSIZE=6,WIDTH=30,Xpos=0.8)
		}
		FNOUT2<-paste0("DEP/statistical.",i,".DEP.",paste0(paste0(LEV,".n",c(sum(SIGNi<0),sum(SIGNi>0))),collapse=".vs."),".txt")
		write.table(DATOUTi,file=FNOUT2,quote=F,row.names=F,col.names=T,sep="\t")

	}

dev.off()


########################
########################
	LINVAR<-rep(list(),3)
	for(i in 1:3){
		DEPi<-LDEPTi[[i]]
		indNNA<-which(!is.na(as.numeric(DEPi[,1])))
		DEPi2<-DEPi[indNNA,]
		PVALi<-as.numeric(DEPi2[,1])
		FDRi<-p.adjust(PVALi)
		DIFFi<-as.numeric(DEPi2[,3])
		LINVAR[[i]]<-indNNA[which(FDRi==1&abs(DIFFi)<1)]
	}
	V1<-venn(LINVAR)

	indABC<-as.integer(attr(V1,"intersections")[["A:B:C"]])
	indHalf<-which(rowSums(!is.na(exprBR))>15)
	ind27<-which(rowSums(!is.na(exprBR))>27)
	Invariome<-setdiff(unique(sort(GeneSym[intersect(indABC,indHalf)])),"")
	Invariome27<-setdiff(unique(sort(GeneSym[intersect(indABC,ind27)])),"")
	MatIN<-sapply(Invariome27,function(gene){
		sapply(LgMatrisome3,function(Set){
			gene%in%Set
		})
	})
apply(MatIN,1,function(x){
	paste0(sort(colnames(MatIN)[x]),collapse=", ")})

pdf("DEP/invariome.pdf")
	pie(rowSums(MatIN)[c(1:6,10)],col=c(1:6,"grey"))
dev.off()

pdf("DEP/invariome-2.pdf")
	pie(c(0,6,2,1,4,2),col=c(1:6))
dev.off()
########################
Human.invariome<-read.table("G:/kathy-cheah-ddd/vivian-tam/old-and-young/invariome/invariome.in.young.245.genes.txt",header=T)[,1]
V1<-venn(list(toupper(Invariome),Human.invariome))

attr(V1,"intersect")[["A:B"]]


indCD<-which(GeneSym%in%Invariome)
fpkm0<-exprBR[indCD,]
range(rowSums(!is.na(fpkm0)))
dat.CD<-data.frame(GeneSym[indCD],fpkm0)
###################
feno<-droplevels(data.frame(BRID=colnames(fpkm0),
		site=substr(HEADER_BR[,6],1,1),
		HEADER_BR[,c(2,6,4)]))
colnames(feno)[5]<-"Gender"
header<-data.frame(matrix("",ncol(feno),1),t(feno))
header[,1]<- colnames(feno)
header[1,1]<- "Genes"

write.table(header,file="invariome/exprtable.Matt.invariome.BioRep.n31.txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(dat.CD,file="invariome/exprtable.Matt.invariome.BioRep.n31.txt",quote=F,row.names=F,col.names=F,sep="\t",append=T)



