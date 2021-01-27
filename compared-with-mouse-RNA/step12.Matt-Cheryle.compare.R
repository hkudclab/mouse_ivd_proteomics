library(scales)
library(stringr)
library(mice)
###################################################
load("../resources/matrisome-mit-mouse/Matrisome.mouse.RData")
load("../resources/Matt.n91.RData")
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
load("../resources/Cheryle.Processed.RData")
library(scales)
library(gplots)

##############################################
SYBM<-intersect(GeneSym,geneSym)
indGPROT<-match(SYBM,GeneSym)
indGARRAY<-match(SYBM,geneSym)
#############################################
table(rowSums(!is.na(exprBR[indGPROT,]))>30)
ind10<-rowSums(!is.na(exprBR[indGPROT,]))>30

EXPRBR<-exprBR[indGPROT,][ind10,]
ARRAY<-expr[indGARRAY,][ind10,]
colnames(EXPRBR)<-paste0(HEADER_BR[,3],"_",HEADER_BR[,4]," (",HEADER_BR[,7],")")
colnames(ARRAY)<-pheno2[,2]
CORR<-cor(EXPRBR,ARRAY,method="spearman",use="pairwise.complete.obs")

DEP_AF<-c("Col2a1", "Col6a2", "Col6a1", "Col1a1", "Col1a2", "Col11a2", "Col11a1", "Col12a1", "Bgn", "Prelp", "Fmod", "Dcn", "Lum", "Vcan", "Kera", "Fn1", "Mfge8", "Sparc", "Comp", "Pcolce", "Thbs1", "Matn3", "Nid1", "Serpinf1", "Clec3a", "S100a13", "Clu", "Ahsg", "Apoe", "Atp5b", "Rpsa", "Myl1", "Mylpf", "Lect1", "Rps12", "Myl3", "Des", "Acat1", "Lect2", "Atp5d")
DEP_NP<-c("Lamp2", "Cd109", "Anxa2", "Lgals3", "Anxa1", "Anxa11", "Anxa7", "Anxa3", "S100a11", "Prdx2", "Hbb?b1", "Ldha", "Eef1a1", "Ca3", "Vim", "Lmna", "H2afj", "Krt19", "Krt8", "Tubb4b", "Sptan1", "Hspb1", "Tuba1b", "Actn4", "Flnb", "Kctd12", "Ptrf", "P4hb", "Atp6v1b2", "Vcp", "Rps8", "Cap1", "Ugdh", "Sptbn1", "Arpc4", "Cct3", "Stip1", "Tubb2a", "Tcp1", "Tubb5", "Hsp90ab1", "Rpl3", "Gstm2", "Ywhaq", "Vat1", "Atic", "Ddx3x")

pdf("cross-compare-with.cheryle.pdf",height=10)
	heatmap.2(CORR,col=bluered(50),margin=c(15,10))

	NPprot<-rowMeans(exprBR[indGPROT,Compartments=="NP"],na.rm=T)
	AFprot<-rowMeans(exprBR[indGPROT,Compartments=="AF"],na.rm=T)
	NParray<-rowMeans(expr[indGARRAY,Compartment=="NP"],na.rm=T)
	AFarray<-rowMeans(expr[indGARRAY,Compartment=="AF"],na.rm=T)

	PCC<-cor(NPprot-AFprot,NParray-AFarray,use="pairwise.complete.obs")
	PCC.p<-cor.test(NPprot-AFprot,NParray-AFarray,use="pairwise.complete.obs")$p.value
	SCC<-cor(NPprot-AFprot,NParray-AFarray,method="spearman",use="pairwise.complete.obs")
dev.off()

pdf("cross-compare-with.cheryle-scatter.pdf")

	Xi<-AFprot-NPprot
	Yi<-AFarray-NParray
	plot(Xi,Yi,pch=16,col="grey",main=paste0("PCC=",PCC,"; p=",PCC.p))
	abline(h=0)
	abline(v=0)
	abline(lm(Yi~Xi))
	inAFi<-which(SYBM%in%DEP_AF)
	inNPi<-which(SYBM%in%DEP_NP)
	points(Xi[inAFi],Yi[inAFi],pch=16,col="#F8766D")
	text(Xi[inAFi],Yi[inAFi],SYBM[inAFi],xpd=T)
	points(Xi[inNPi],Yi[inNPi],pch=16,col="#619CFF")
	text(Xi[inNPi],Yi[inNPi],SYBM[inNPi],xpd=T)


dev.off()









