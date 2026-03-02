
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(tidyverse)
library(TCGAbiolinks)
library(stringr)
library(SummarizedExperiment)
projs <- getGDCprojects()$project_id %>% 
  str_subset("TCGA")
projs
proj = "TCGA-GC"
load("stad_exp.Rdata")
exp = stad
library(tinyarray)
exp = trans_exp_new(exp)
exp[1:4,1:4]
nrow(exp)
exp = exp[apply(exp, 1, function(x) sum(x > 0) > 0.5*ncol(exp)), ]
nrow(exp)
"EDNRA" %in% rownames(exp)
Group = make_tcga_group(exp)
table(Group)
exp = exp[,Group=='tumor']
data = exp

gene="EDNRA"             
dataL=data[,data[gene,]<=median(data[gene,]),drop=F]
dataH=data[,data[gene,]>median(data[gene,]),drop=F]
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
logFC=sort(logFC,decreasing=T)
genes=names(logFC)

gmt=read.gmt("c2.cp.kegg.v7.4.symbols.gmt")
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
write.table(kkTab,file=paste0(gene,"_GSEA.result-KEGG.txt"),sep="\t",quote=F,row.names = F)

termNum=10  
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="")
  gseaplot[[1]]=gseaplot[[1]]+theme(legend.position="right",legend.direction="vertical")
  pdf(file=paste0(gene,"_GSEA-KEGG.pdf"), width=10, height=8)
  print(gseaplot)
  dev.off()
}
if(nrow(kkTab)<termNum){
  showTerm=row.names(kkTab)[1:nrow(kkTab)]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="")
  gseaplot[[1]]=gseaplot[[1]]+theme(legend.position="right",legend.direction="vertical")
  pdf(file=paste0(gene,"_GSEA-KEGG.pdf"), width=10, height=8)
  print(gseaplot)
  dev.off()
}


