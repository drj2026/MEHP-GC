
options(timeout = 100000) 
options(scipen = 20)


library(GEOquery)
eSet = getGEO("GSE54129", destdir = '.', getGPL = F)
class(eSet)
length(eSet)
eSet = eSet[[1]] 
class(eSet)

# Extract expression matrix (exp)
exp <- exprs(eSet)
dim(exp)
range(exp)
boxplot(exp,las = 2) 
exp = limma::normalizeBetweenArrays(exp)
boxplot(exp[,sample(1:130,20)],las = 2)

# Extract clinical information
pd <- pData(eSet)

library(stringr)
k1 = str_detect(pd$title,"normal");table(k1)
k2 = str_detect(pd$title,"tumor");table(k2)
pd = pd[k1|k2,]
identical(rownames(pd),colnames(exp))

# Extract platform ID
gpl_number <- eSet@annotation;gpl_number
save(pd,exp,gpl_number,file = "step1output.Rdata")

# Experimental grouping
rm(list = ls())  
load(file = "step1output.Rdata")
library(stringr)
k = str_detect(pd$title,"normal");table(k) #不在title就在pd的其他列
Group = ifelse(k,"Normal","GC")
Group = factor(Group,levels = c("Normal","GC"))
data.frame(pd$title,Group)
# Obtain probe annotation
library(tinyarray)
gpl_number 
pkg_all[pkg_all$gpl==gpl_number,2]
if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db",ask = F,update = F)
library(hgu133plus2.db)
ids <- toTable(hgu133plus2SYMBOL) 
save(exp,Group,ids,file = "step2output.Rdata")

# PCA plot
rm(list = ls())  
load(file = "step2output.Rdata")
dat=as.data.frame(t(exp))
library(FactoMineR)
library(factoextra) 
dat.pca <- PCA(dat, graph = FALSE)#graph = FALSE默认参数不需要管
d=fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = Group, # color by groups
               palette = c("#016074", "#b02113"),
               addEllipses = T, # Concentration ellipses
               legend.title = "Groups")
ggsave(d,file="GC_PCA.png")

# Differential expression analysis
rm(list = ls()) 
load(file = "step2output.Rdata")
library(limma)
design = model.matrix(~Group)
fit = lmFit(exp,design)
fit = eBayes(fit)
deg = topTable(fit,coef = 2,number = Inf)#coef = 2统计design的第二列，两组不需要改动，Inf正无穷全部列出来的意思


library(dplyr)
deg = mutate(deg,probe_id = rownames(deg))
ids = distinct(ids,symbol,.keep_all = T)
deg = inner_join(deg,ids,by="probe_id")


# Add 'change' column to label up- and down-regulated genes
logFC_t = 0.5
p_t = 0.01
k1 = (deg$P.Value < p_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < p_t)&(deg$logFC > logFC_t)
deg = mutate(deg,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(deg$change)

# Volcano plot
library(ggplot2)
g=ggplot(data = deg, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, aes(color=change)) +
  scale_color_manual(values=c("#016074", "grey","#b02113"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",linewidth=0.8) +
  geom_hline(yintercept = -log10(p_t),lty=4,col="black",linewidth=0.8) +
  theme_bw()
ggsave(g,file="GC_vol.png")

# Save differentially expressed genes (DEGs)
GC_gene_diff = deg$symbol[deg$change != "stable"] 
write.csv(GC_gene_diff,
          file = "GC_gene_diff.csv",
          row.names = FALSE)

# Heatmap of differentially expressed genes (DEGs)
exp = exp[deg$probe_id,]
rownames(exp) = deg$symbol
n1 = deg[deg$change !="stable",]
n2 = n1[order(n1$P.Value),]
n3 = n2[1:50,]
diff_gene50 = n3$symbol
n = exp[diff_gene50,]
library(pheatmap)
annotation_col = data.frame(group = Group)
rownames(annotation_col) = colnames(n) 
ann_colors = list(
  group = c(Normal = "#016074", GC = "#b02113"))
g=pheatmap(n,
           show_colnames =F,
           show_rownames = T,
           fontsize_row = 7,
           cluster_cols=F,#不聚类
           annotation_col=annotation_col,
           annotation_colors = ann_colors,
           scale = "row", #按行标准化，只保留行内差别，不保留行间差别，会把数据范围缩放到大概-5~5之间
           breaks = seq(-3,3,length.out = 100), #设置色带分布范围为-3~3之间，超出此范围的数字显示极限颜色，#length.out = 100有100个渐变颜色
           color = colorRampPalette(c("#016074", "white", "#b02113"))(100))
ggsave(g,file="GC deg Heat Map50.png")