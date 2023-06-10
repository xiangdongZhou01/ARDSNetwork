library(WGCNA)
enableWGCNAThreads()
nGenes = nrow(redup.expr); nSamples = nrow(keep.meta)
corType = "cor"
networkType = "signed"
TOMType = "signed"
TOMDenom = "mean"

powers = c(seq(4,10,by=1), seq(12,30, by=2))
sft = pickSoftThreshold(t(redup.expr), powerVector = powers, verbose = 5,corFnc = corType,networkType = networkType)
sizeGrWindow(9, 5)
pdf(file = "../figures/scaleFreeAnalysis.pdf",width = 8,height = 6)
par(mfrow = c(1,2))
cex1 = 0.6
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
power = sft$powerEstimate
adj <- adjacency(t(redup.expr), power = power,type = networkType,corFnc = corType)
TOM <- TOMsimilarity(adj,TOMType = TOMType,TOMDenom = TOMDenom)
dissTOM <- 1-TOM
geneTree = hclust(as.dist(dissTOM),method="average")
minModSize <- 30 # Modules are at least 20 genes large
dthresh <- 0.25 # MEs are no more than 0.85 correlated, if they are then the modules are merged and the ME is re-calculated
ds <- 4 # deep split parameter to determine how finely to cut the tree
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = ds, pamRespectsDendro = FALSE,minClusterSize = minModSize,method = "hybrid")
merged <- mergeCloseModules(exprData = t(redup.expr),colors = dynamicMods,cutHeight = dthresh)
mColorh <- labels2colors(merged$colors)
table(mColorh)

MEs = moduleEigengenes(t(redup.expr), mColorh)$eigengenes
MEs = orderMEs(MEs)
moduleTraitCor = cor(MEs, keep.meta$Group.index+1, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), ", p = ",signif(moduleTraitPvalue, 1), sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf(file = "../figures/Module-Group correlation.pdf",width = 6,height = 10)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,xLabels = "Group",yLabels = names(MEs),ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,zlim = c(-1,1),main = paste("Module-Group relationships"))
dev.off()

which.max(abs(moduleTraitCor[-which(rownames(moduleTraitCor) == "MEgrey"),]))
# MEdarkgrey 
which.min(moduleTraitPvalue[-which(rownames(moduleTraitPvalue) == "MEgrey"),])
# MEdarkgrey 
GS1 <- as.numeric(cor(t(redup.expr),keep.meta$Group.index+1,use="p"))
GS.pvalue <- as.data.frame(corPvalueStudent(as.matrix(GS1), nSamples))
GeneSignificance <- abs(GS1)
ModuleSignificance <- tapply(GeneSignificance,mColorh,mean,na.rm=T)
which.max(ModuleSignificance[names(ModuleSignificance)!="grey"])
# darkgrey 
ind = which(mColorh=="grey")
pdf(file = "../figures/Gene significance across modules.pdf",width = 10,height = 6)
plotModuleSignificance(GeneSignificance[-ind],mColorh[-ind],main = "Gene significance across modules,",ylab = "Gene Significance",ylim = c(0,0.4),cex.lab = 1.5)
dev.off()
geneModuleMembership = as.data.frame(cor(t(redup.expr), MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
modNames = substring(names(MEs), 3)
pdf(file = "../figures/MMvsGS.pdf",width = 7,height = 6)
module = "darkgrey"
moduleGenes = mColorh==module
column = match(module, modNames)
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),abs(GS1[moduleGenes]),xlab = paste("Module Membership in", module, "module"),ylab = "Gene significance for group",main = paste("Module membership vs. gene significance\n"),cex.main = 0.8, cex.lab = 0.6, cex.axis = 0.7, col = module,pch = 16)
abline(lm(abs(GS1[moduleGenes])~abs(geneModuleMembership[moduleGenes, column])),col = "red")
dev.off()

KIM <- intramodularConnectivity(adjMat = adj,mColorh,scaleByMax = T)

library(gplots)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
plotTOM = dissTOM^power
diag(plotTOM) = NA
tiff(filename = "../figures/Network.heatmap.tiff",width = 5,height = 5,res = 300,units = "in")
TOMplot(plotTOM, geneTree, mColorh, main = "Network heatmap plot",col = myheatcol)
dev.off()

pdf(file = "../figures/MMvsGS2.pdf",width = 7,height = 6)
module = "darkmagenta"
moduleGenes = mColorh==module
column = match(module, modNames)
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),abs(GS1[moduleGenes]),xlab = paste("Module Membership in", module, "module"),ylab = "Gene significance for group",main = paste("Module membership vs. gene significance\n"),cex.main = 0.8, cex.lab = 0.6, cex.axis = 0.7, col = module,pch = 16)
abline(lm(abs(GS1[moduleGenes])~abs(geneModuleMembership[moduleGenes, column])),col = "red")
dev.off()


