library(lumi)
library(limma)
library(openxlsx)
library(ggsci)

filename <- "../data/GSE32707_non_normalized.txt"
x <- read.ilmn(filename,expr="AVG_Signal",probeid="PROBE_ID",other.columns="Detection Pval")
y <- neqc(x,detection.p="Detection Pval")
human.expr <- y$E
metadata <- read.xlsx("../data/GSE32707.samples.xlsx",sheet = 1)
metadata <- metadata[order(metadata$Group,decreasing = F),]
table(metadata$Group)
#  control se/ARDS  Sepsis 
#      34      31      58 

group.col = pal_nejm(c("default"))(3)
names(group.col) = unique(metadata$Group)
human.expr = human.expr[,match(metadata$Sample_description,colnames(human.expr))]
colnames(human.expr) = metadata$ID

pdf(file = "../figures/human.normalized.boxplot.pdf",width = 8,height = 3.5)
boxplot(data.frame(human.expr),outline = F,main = "Human",col = rep(group.col,as.numeric(table(metadata$Group))))
dev.off()


clinical <- openxlsx::read.xlsx("../data/GSE32707.samples.xlsx",sheet = 2)
clinical = clinical[,1:8]
clinical$Group = 0
clinical$Group[which(clinical$sepsis==0&clinical$ARDS==0)] = "Control"
clinical$Group[which(clinical$sepsis==1&clinical$ARDS==0)] = "Sepsis"
clinical$Group[which(clinical$sepsis==1&clinical$ARDS==1)] = "ARDS"
table(clinical$Group)
str(clinical)
mydata = clinical
mydata$Gender[which(mydata$Gender=="M")] = 1
mydata$Gender[which(mydata$Gender=="F")] = 2
mydata$Gender = as.numeric(mydata$Gender)
mydata$Gender = factor(mydata$Gender,levels = c(1,2),labels = c("Male","Female"))
mydata$Group[which(mydata$Group=="Control")] = 1
mydata$Group[which(mydata$Group=="Sepsis")] = 2
mydata$Group[which(mydata$Group=="ARDS")] = 3
mydata$Group = as.numeric(mydata$Group)
mydata$Group = factor(mydata$Group,levels = c(1,2,3),labels = c("Control","Sepsis","ARDS"))
mytable = table(mydata$Gender,mydata$Group)
addmargins(mytable,2)
print(prop.table(mytable, 1),digits = 3)
chisq.test(mytable)
mean(mydata$Age[which(mydata$Group=="Control")])
sd(mydata$Age[which(mydata$Group=="Control")])
mean(mydata$Age[which(mydata$Group=="Sepsis")])
sd(mydata$Age[which(mydata$Group=="Sepsis")])
mean(mydata$Age[which(mydata$Group=="ARDS")])
sd(mydata$Age[which(mydata$Group=="ARDS")])
shapiro.test(mydata$Age[which(mydata$Group=="Control")])
shapiro.test(mydata$Age[which(mydata$Group=="Sepsis")])
shapiro.test(mydata$Age[which(mydata$Group=="ARDS")])
summary(aov(Age~Group,data = mydata))



normadj <- (0.5+0.5*bicor(human.expr))^2
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers = (z.ku < -2.5)
table(outliers)
pdf(file = "../figures/outliers.detect.pdf",width = 8,height = 7)
par(mar = rep(6,4))
plot(1:length(z.ku),z.ku,col = metadata$color,pch=19, main = "Outlier Detection", xlab="", ylab="Standardized Network\nConnectivity (Z score)")
abline(h=-2.5, lty=2)
legend("bottomleft",pch = 16,col = group.col,legend = names(group.col),bty = "n")
dev.off()

keep.meta <- metadata[!outliers,]
keep.expr <- human.expr[,!outliers]

filters <- rownames(keep.expr)
str(filters)
annotation_illuminagene <- function(x){
  require("biomaRt")
  mart=useMart("ensembl",dataset="hsapiens_gene_ensembl")
  attributes=c("illumina_humanht_12_v4","entrezgene_id","hgnc_symbol","gene_biotype")
  genes = as.data.frame(getBM(attributes=attributes,filters=c("illumina_humanht_12_v4"),values=x,mart=mart))
  return(genes)
}

xv = annotation_illuminagene(filters)
xv = xv[!is.na(xv$entrezgene_id),]
xv = xv[which(xv$hgnc_symbol!=""),]
xv = subset(xv,xv$gene_biotype=="protein_coding")
table(duplicated(xv$illumina_humanht_12_v4)) # 3024
dup.probes = xv[duplicated(xv$illumina_humanht_12_v4),1]
xv = xv[!xv$illumina_humanht_12_v4%in%dup.probes,]
xv = xv[-grep("LINC",xv$hgnc_symbol),]
table(duplicated(xv$entrezgene_id))
redup.expr = keep.expr[match(xv$illumina_humanht_12_v4,rownames(keep.expr)),]
a = collapseRows(redup.expr,rowGroup = xv$entrezgene_id,rowID = rownames(redup.expr),method = "MaxMean")
redup.expr = redup.expr[a$selectedRow,] # 18066
rownames(redup.expr) = xv$hgnc_symbol[a$selectedRow]
human.anno <- xv
gsg = goodSamplesGenes(t(redup.expr), verbose = 3)
gsg$allOK
sampleTree = hclust(dist(t(redup.expr)), method = "average")
pdf(file = "../figures/sampleclustering outliers.pdf",width = 6,height = 4)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 0.8,cex.axis = 0.8, cex.main = 1,cex = 0.3)
abline(h = 180, col = "red");
dev.off()
clust = cutreeStatic(sampleTree, cutHeight = 180, minSize = 10)
keep.samples = clust==1
keep.meta = keep.meta[keep.samples,]
redup.expr = redup.expr[,keep.samples]
save(redup.expr,keep.meta,group.col,file = "human.filter.expr.RData")
