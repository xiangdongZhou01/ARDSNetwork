library(limma)
rm(list = ls()) 
load("human.filter.expr.RData")
group <- as.factor(keep.meta$Group)
design <- model.matrix(~0+group)
colnames(design) = levels(group)
contrast.matrix <- makeContrasts(ARDS-Sepsis,ARDS-control,Sepsis-control, levels=design)
fit <- lmFit(redup.expr, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2,trend = T,robust = T)
human.res <- lapply(1:3,function(x){topTable(fit2,coef = x,sort.by = "none",adjust.method = "BH",n=Inf)})
names(human.res) = c("ARDS-Sepsis","ARDS-Control","Sepsis-Control")
human.DEGs <- lapply(human.res,function(x){subset(x,x$adj.P.Val<0.05)})
str(human.DEGs)
human.DEGs.gene <- lapply(human.DEGs,rownames)
str(human.DEGs.gene)
# List of 3
# $ ARDS-Sepsis   : chr [1:162] "AOC1" "ABCB6" "AADACL4" "ALKBH5" ...
# $ ARDS-control  : chr [1:150] "ACAP1" "ABCA7" "ADORA2A" "AFTPH" ...
# $ Sepsis-control: chr [1:180] "ACAP1" "AKAP3" "ADORA2A" "ADRA2C" ...

venn.plot <- venn.diagram(human.DEGs.gene,filename =NULL,height = 5,width = 4.5,resolution = 700,col = "black",fill = c("seagreen","purple","orange"),alpha = 0.50, lwd = 0.7,label.col = "black",cex = 1,fontfamily = "sans",cat.col = c("darkgreen","purple","darkorange"),cat.cex = 0.6,cat.pos = 0,cat.dist = 0.07,cat.fontfamily = "sans",cat.fontface="bold",margin = 0.2)
pdf(file = "../../figures/human.DEGs.venn.pdf",w = 5,h = 4.5)
grid.draw(venn.plot)
dev.off()
for(i in 1:3){
 write.table(data.frame(Gene = rownames(human.DEGs[[i]]),human.DEGs[[i]]),file = paste("DEGs/",names(human.DEGs)[i],".human.DEGs.txt",sep = ""),sep = "\t",col.names = T,row.names = F) 
}

library(ggplot2);library(ggpubr);
volcano.list <- p.list <- list()
ind = c(2,3,1)
for(i in 1:3){
  volcano.list[[i]] <- data.frame(human.res[[ind[i]]],logP = -log10(human.res[[ind[i]]]$adj.P.Val),sig = "N")
  volcano.list[[i]]$sig[which(volcano.list[[i]]$adj.P.Val<0.05)] = "Y"
  p.list[[i]] <- ggplot(volcano.list[[i]],aes(x = logFC,y = logP,color = sig)) + geom_point(size = 0.4,alpha = 0.6) + scale_color_manual(values = c("grey80","indianred")) + theme_bw() + labs(x = "log2 Fold-change",y = "-log10 adjusted p-value",title = names(human.res)[ind[i]],element_text(size = 2)) + theme(plot.title = element_text(hjust = 0.5,face = "bold"))
}
names(p.list) = names(human.DEGs.gene)[ind]
p <- ggarrange(plotlist = p.list,nrow = 1,ncol = 3)
ggsave(p,file = "../../figures/volcano.pdf",width = 180,height = 50,units = "mm")

#### GSEA analysis
rnk.list <- list()
for(i in 1:3){
  rnk.list[[i]] = data.frame(Gene = rownames(volcano.list[[i]]),volcano.list[[i]][,c(1,7)])
  rnk.list[[i]]$metrics = rnk.list[[i]]$logP*sign(rnk.list[[i]]$logFC)
  rnk.list[[i]] = rnk.list[[i]][,c(1,4)]
  write.table(rnk.list[[i]],file = paste(names(human.res)[i],".rnk",sep = ""),sep = "\t",col.names = T,row.names = F)
}

a <- dir("GSEA/")
rnk.files <- paste("GSEA/",a,sep = "")
outputDirectory = "GSEA/"
for(i in 1:3){
  enrichResult <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens", enrichDatabase="geneontology_Biological_Process_noRedundant", interestGeneFile=rnk.files[i], interestGeneType="genesymbol", referenceGeneType="genesymbol", sigMethod = "fdr",isOutput=TRUE, outputDirectory=outputDirectory,projectName=a[i],nThreads = 16,gseaPlotFormat = "svg",fdrThr = 0.25)
}

setwd("GSEA/")
gsea.go <- list(read.table("Project_ARDS_control_go/enrichment_results_ARDS_control_go.txt",header = T,sep = "\t",quote = ""),read.table("Project_ARDS_Sepsis_go/enrichment_results_ARDS_Sepsis_go.txt",header = T,sep = "\t",quote = ""),read.table("Project_Sepsis_control_go/enrichment_results_Sepsis_control_go.txt",header = T,sep = "\t",quote = ""))
names(gsea.go) = names(human.res)[c(2,1,3)]
gsea.go <- lapply(gsea.go, function(x){x[abs(x$normalizedEnrichmentScore)>1.5,]})
length(unique(c(gsea.go$`ARDS-control`$description,gsea.go$`Sepsis-control`$description,gsea.go$`ARDS-Sepsis`$description)))
# 54
terms <- unique(c(gsea.go$`ARDS-control`$description,gsea.go$`Sepsis-control`$description,gsea.go$`ARDS-Sepsis`$description))
a = list(gsea.go$`ARDS-control`$description,gsea.go$`Sepsis-control`$description)
names(a) = names(human.DEGs)[2:3]
venn.plot <- venn.diagram(a,filename =NULL,height = 850,width = 950,resolution =300,col = "black",fill = c("seagreen","orange"),alpha = 0.50, lwd=0.8,label.col = "black",cex = 0.8,fontfamily = "sans",cat.col = c("darkgreen","darkorange"),cat.cex = 0.6,cat.pos = 0,cat.dist = 0.07,cat.fontfamily = "sans",cat.fontface="bold",margin = 0.2)
pdf(file = "../../figures/gsea.go.pdf",w = 5,h = 4.5,onefile = F)
grid.draw(venn.plot)
dev.off()

terms <- unique(c(gsea.go$`ARDS-Control`$description,gsea.go$`Sepsis-Control`$description,gsea.go$`ARDS-Sepsis`$description))
ind <- lapply(gsea.go, function(x){match(terms,x$description)})
go.matrix <- data.frame(term = rep(terms,3),logfdr = -log10(c(gsea.go$`ARDS-Control`$FDR[ind[[1]]],gsea.go$`Sepsis-Control`$FDR[ind[[3]]],gsea.go$`ARDS-Sepsis`$FDR[ind[[2]]])),NES = c(gsea.go$`ARDS-Control`$normalizedEnrichmentScore[ind[[1]]],gsea.go$`Sepsis-Control`$normalizedEnrichmentScore[ind[[3]]],gsea.go$`ARDS-Sepsis`$normalizedEnrichmentScore[ind[[2]]]),group = rep(c("ARDS-Control","Sepsis-Control","ARDA-Sepsis"),each = length(terms)))
go.matrix$NES = abs(go.matrix$NES)
go.matrix$logfdr[go.matrix$logfdr==Inf] = 3
go.matrix$term = Hmisc::capitalize(as.character(go.matrix$term))
p1 <- ggplot(go.matrix,aes(x = group,y = term,color = logfdr,size = NES)) + geom_point() + scale_color_gradient(low = "skyblue",high = "midnightblue") + theme_bw() + theme(axis.text.x = element_text(angle = 35,vjust = 0.6,color = "black"),axis.text.y = element_text(color = "black"),plot.title = element_text(hjust = 0.5)) + labs(x = "",y = "",title = "GSEA")
ggsave(p1,filename = "../figures/gsea.go.point.pdf",width = 8,height = 12)