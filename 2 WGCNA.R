library(tinyarray)
gse = "GSE54129"
geo = geo_download(gse,destdir=".",colon_remove = T)
geo$gpl
# Group information
library(stringr)
k = str_detect(geo$pd$title,"normal");table(k) #不在title就在pd的其他列
Group = ifelse(k,"Normal","GC")
Group = factor(Group,levels = c("Normal","GC"))
Group
table(Group)
data.frame(geo$pd$title,Group)
# Obtain probe annotation
gpl_number=geo$gpl
pkg_all[pkg_all$gpl==gpl_number,2]
if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db",ask = F,update = F)
library(hgu133plus2.db)
ids <- toTable(hgu133plus2SYMBOL)
ids = ids[na.omit(ids$symbol!=""),]
exp = trans_array(geo$exp,ids)
pd = geo$pd
table(Group)
save(exp,Group,pd,file = "Dat.Rdata")

rm(list = ls())
library(WGCNA)
library(tinyarray)
load("Dat.Rdata")
range(exp)
boxplot(exp[,sample(1:132,20)],las = 2)
exp = limma::normalizeBetweenArrays(exp)
png("1.exp.png", width = 2000, height = 1200,res = 300)
boxplot(exp)
dev.off()

datExpr0 = t(exp[order(apply(exp, 1, mad), decreasing = T)[1:10000],])
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK 
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average")

png(file = "2.sampleClustering.png", width = 2000, height = 2000,res = 300)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

if(F){
  clust = cutreeStatic(sampleTree, cutHeight = 170, minSize = 10)
  table(clust) 
  keepSamples = (clust!=0)
  datExpr = datExpr0[keepSamples, ]
}else{
  datExpr = datExpr0
}

library(stringr)
traitData = data.frame(
  Normal = as.numeric(pd$tissue=="normal"),
  GC = as.numeric(pd$tissue=="gastric cancer tumor"))

datTraits = traitData
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
png(file = "3.Sample dendrogram and trait heatmap.png", width = 2000, height = 2000,res = 300)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Select soft-thresholding power
powers = c(1:10, seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft$powerEstimate

cex1 = 0.87
png(file = "4.Soft threshold.png", width = 2000, height = 1500,res = 300)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=cex1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

power = sft$powerEstimate
net = blockwiseModules(datExpr, power = power,
                       TOMType = "unsigned", 
                       minModuleSize =30, 
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25,
                       deepSplit = 2 ,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "testTOM",
                       verbose = 3)

mergedColors = labels2colors(net$colors)
png(file = "5.DendroAndColors.png", width = 2000, height = 1200,res = 300)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Save genes for each module
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
gm = data.frame(net$colors)
gm$color = moduleColors
head(gm)
genes = split(rownames(gm),gm$color)
save(genes,file = "genes.Rdata")

# Module–trait correlation
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
png(file = "6.labeledHeatmap.png", width = 2000, height = 2000,res = 300)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed (50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

i = 2 #⭐
module = "blue" # module = c("blue","turquoise"，"green"，"pink")
assign(colnames(traitData)[i],traitData[i])
instrait = eval(parse(text = colnames(traitData)[i]))
geneTraitSignificance = as.data.frame(cor(datExpr, instrait, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(instrait), sep="")
names(GSPvalue) = paste("p.GS.", names(instrait), sep="")
png(file = paste0("72.MM-GS-scatterplot.png"), width = 2000, height = 2000,res = 300)
column = match(module, modNames) 
moduleGenes = moduleColors==module 
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

# TOM plot
nSelect = 400
set.seed(101)
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 7)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
library(gplots)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
png(file = "8.Sub400-netheatmap.png", width = 2000, height = 2000,res = 300)
plotDiss = selectTOM^7
diag(plotDiss) = NA 
TOMplot(plotDiss, selectTree, selectColors, col=myheatcol,main = "Network heatmap plot, selected genes")
dev.off()