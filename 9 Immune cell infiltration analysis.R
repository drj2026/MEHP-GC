library(CIBERSORT)
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(reshape2)
library(tidyverse)
DEG_expr <- read.csv("DEG_expr.csv",row.names = 1)
group <- read.csv("group.csv",row.names = 1)

# Data Validation
boxplot(DEG_expr,outline=F, notch=F , las=2)
qx <- as.numeric(quantile(DEG_expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))  #数据的分布，样本分位数
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)#判断是否进行log的标准
LogC
if (LogC) { 
  DEG_expr[which(DEG_expr <= 0)] <- NaN
  DEG_expr <- log2(DEG_expr) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
}   

library(limma)
DEG_expr=normalizeBetweenArrays(DEG_expr)
boxplot(DEG_expr,outline=FALSE, notch=F , las=2)

LM22_local <- read.table("LM22.txt",header = T,row.names = 1,sep = "\t")
data(LM22)
all(LM22==LM22_local)
result <- cibersort(sig_matrix = LM22, mixture_file = DEG_expr, perm = 1, QN = TRUE)
result <- as.data.frame(result)
colnames(result)
write.csv(result,"cibersort_result.csv")
result1 <- result[,1:ncol(LM22)]
result1 <- result1[,apply(result1, 2, function(x){sum(x)>0})]#删除全是0的列

data3 <- cbind(result1,group)
colnames(data3)
data3 <- data3[,c(23,24,1:22)]
data3 <- pivot_longer(data = data3,
                      cols = 3:24,
                      names_to = "celltype",
                      values_to = "proportion")
# plot
pdf(file="immune infiltration.pdf",width = 10,height = 8)
ggboxplot(data = data3,
          x = "celltype",
          y = "proportion",
          combine = TRUE,
          merge = FALSE,
          color = "black",
          fill = "group",
          palette = c("#b02113","#016074"),
          title = "TME Cell composition",
          xlab = NULL,
          ylab = "Cell composition",
          bxp.errorbar = FALSE,
          bxp.errorbar.width = 0.2,
          facet.by = NULL,
          panel.labs = NULL,
          short.panel.labs = TRUE,
          linetype = "solid",
          size = NULL,
          width = 0.8,
          notch = FALSE,
          outlier.shape = NA,
          select = NULL,
          remove = NULL,
          order = NULL,
          error.plot = "pointrange",
          label = NULL,
          font.label = list(size = 12, color = "black"),
          label.select = NULL,
          repel = TRUE,
          label.rectangle = TRUE, ggtheme = theme_pubr())+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1)) +  
          stat_compare_means(label = "p.signif",method = "t.test",aes(group = group),hide.ns = T) 
dev.off()  

