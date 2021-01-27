library(scales)
library(stringr)
library(mice)
###################################################
load("../resources/Matrisome.mouse.RData")
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

colnames(exprBR)<-paste0(HEADER_BR[,3],"_",HEADER_BR[,4]," (",HEADER_BR[,7],")")

barplot(exprBR[which(GeneSym=="Acan"),],las=2)

##############################################
Compartments<-HEADER_BR[,2]
Gender<-HEADER_BR[,4]
Levels<-factor(substr(HEADER_BR[,3],1,1))
Litter<-factor(substr(HEADER_BR[,3],4,4))
##############################################
#######################################
oydata<-read.csv("../resources/old-and-young.aug21-2019.txt",sep="\t")

head(oydata)
geneSymb<-as.character(oydata[,3])
geneSymb2<-gsub(";.*$","",geneSymb)

nonMatrisome<-setdiff(geneSymb2,MATRISOME[,1])

data0<-as.matrix(oydata[,-seq(3)])

agegrp<-"young"
sampleIDS<-colnames(data0)
vLpheno<-strsplit(sampleIDS,"\\.")
vAges<-tolower(sapply(vLpheno,function(x)x[2]))

vDATA<-data0[,vAges%in%unlist(strsplit(agegrp,"_"))]
sampleIDS<-colnames(vDATA)
vLpheno<-strsplit(sampleIDS,"\\.")

vLevels<-sapply(vLpheno,function(x)x[1])
vAges<-tolower(sapply(vLpheno,function(x)x[2]))
vDirections<-sapply(vLpheno,function(x)x[3])
vCompart<-sapply(vLpheno,function(x)x[4])
vCompart[is.na(vCompart)]<-vDirections[is.na(vCompart)]

#######################################
#######################################
library(scales)
library(gplots)

##############################################
load("../resources/human.mouse.2020.RData")
GENESYMB<-as.character(human.mouse.2020[match(GeneSym,human.mouse.2020[,5]),6])

SYBM<-intersect(GENESYMB,geneSymb2)
indMATT<-match(SYBM,GENESYMB)
indVivian<-match(SYBM,geneSymb2)
#############################################
library(ggplot2)
SimilarDraw<-function(LCORR2,strMain){
	LCORR3<-c(LCORR2[[1]],LCORR2[[2]],LCORR2[[3]],LCORR2[[4]])

	fac<-rep(names(LCORR3),sapply(LCORR3,length))
	quant<-unlist(LCORR3)

	dat34<-data.frame(quant,fac)
	MED<-aggregate(quant,by=list(fac),median)
	print(MED$x)
	p34<- ggplot(dat34, aes(x=fac, y=quant)) + 
		geom_violin( scale = "width")+ 
		geom_boxplot(width=0.45,outlier.shape =NA)+
		geom_jitter(size=2,width=0.2,height=0,color="#808080")+
		geom_text(data=MED,aes(x = Group.1, y = max(quant), 
			label=signif(MED$x,3)),color="red",angle=90)

		tres<-summary(aov(quant~fac))
		p34<-p34+ggtitle(paste0(strMain,": p=",signif(as.matrix(tres)[[1]][1,5],3)))

	plot(p34)

}
#############################################
table(rowSums(!is.na(exprBR[indMATT,]))>20)
ind10<-rowSums(!is.na(exprBR[indMATT,]))>20

EXPRBR<-exprBR[indMATT,][ind10,]
VIVIAN<-vDATA[indVivian,][ind10,]
colnames(EXPRBR)<-paste0(HEADER_BR[,3],"_",HEADER_BR[,4]," (",HEADER_BR[,7],")")

CORR<-cor(EXPRBR,VIVIAN,method="spearman",use="pairwise.complete.obs")

DEP_AF<-c("Col2a1", "Col6a2", "Col6a1", "Col1a1", "Col1a2", "Col11a2", "Col11a1", "Col12a1", "Bgn", "Prelp", "Fmod", "Dcn", "Lum", "Vcan", "Kera", "Fn1", "Mfge8", "Sparc", "Comp", "Pcolce", "Thbs1", "Matn3", "Nid1", "Serpinf1", "Clec3a", "S100a13", "Clu", "Ahsg", "Apoe", "Atp5b", "Rpsa", "Myl1", "Mylpf", "Lect1", "Rps12", "Myl3", "Des", "Acat1", "Lect2", "Atp5d")
DEP_NP<-c("Lamp2", "Cd109", "Anxa2", "Lgals3", "Anxa1", "Anxa11", "Anxa7", "Anxa3", "S100a11", "Prdx2", "Hbb?b1", "Ldha", "Eef1a1", "Ca3", "Vim", "Lmna", "H2afj", "Krt19", "Krt8", "Tubb4b", "Sptan1", "Hspb1", "Tuba1b", "Actn4", "Flnb", "Kctd12", "Ptrf", "P4hb", "Atp6v1b2", "Vcp", "Rps8", "Cap1", "Ugdh", "Sptbn1", "Arpc4", "Cct3", "Stip1", "Tubb2a", "Tcp1", "Tubb5", "Hsp90ab1", "Rpl3", "Gstm2", "Ywhaq", "Vat1", "Atic", "Ddx3x")

pdf("cross-compare/with.Vivian.pdf",width=10,height=10)
	heatmap.2(CORR,col=bluered(50),margin=c(10,10))

	LCORR<-lapply(levels(HEADER_BR$Compartment),function(LEV){
		tmpi<-lapply(unique(vCompart),function(vCPMT){
			XY<-CORR[which(HEADER_BR$Compartment==LEV),vCompart==vCPMT]
			as.vector(XY)
		})
		names(tmpi)<-paste0("m",LEV,"_h",unique(vCompart))
		print(str(tmpi))
		tmpi
	})
	SimilarDraw(LCORR,"Compare with Vivian")

dev.off()

pdf("cross-compare.with.Vivian-scatter.pdf")

	NPprot<-rowMeans(exprBR[indMATT,Compartments=="NP"],na.rm=T)
	AFprot<-rowMeans(exprBR[indMATT,Compartments=="AF"],na.rm=T)
	NParray<-rowMeans(vDATA[indVivian,vCompart!="OAF"],na.rm=T)
	AFarray<-rowMeans(vDATA[indVivian,vCompart=="OAF"],na.rm=T)

	PCC<-cor(NPprot-AFprot,NParray-AFarray,use="pairwise.complete.obs")
	PCC.p<-cor.test(NPprot-AFprot,NParray-AFarray,use="pairwise.complete.obs")$p.value
	SCC<-cor(NPprot-AFprot,NParray-AFarray,method="spearman",use="pairwise.complete.obs")

	Xi<-AFprot-NPprot
	Yi<-AFarray-NParray
	plot(Xi,Yi,pch=16,col="grey",main=paste0("PCC=",PCC,"; p=",PCC.p))
	abline(h=0)
	abline(v=0)
	abline(lm(Yi~Xi))
	inAFi<-which(SYBM%in%toupper(DEP_AF))
	inNPi<-which(SYBM%in%toupper(DEP_NP))
	inAFi2<-which(SYBM%in%toupper(DEP_AF)&(Xi<0&Yi<0|Xi>0&Yi>0))
	inNPi2<-which(SYBM%in%toupper(DEP_NP)&(Xi<0&Yi<0|Xi>0&Yi>0))
	points(Xi[inAFi],Yi[inAFi],pch=16,col="#F8766D")
	text(Xi[inAFi2],Yi[inAFi2],SYBM[inAFi2],xpd=T)
	points(Xi[inNPi],Yi[inNPi],pch=16,col="#00BFC4")
	text(Xi[inNPi2],Yi[inNPi2],SYBM[inNPi2],xpd=T)


dev.off()









