library(clusterProfiler)
library(ggthemes)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
library(ggalluvial)
join=read.delim("join_target.txt" )
colnames(join) = "symbol"
join_gene = join
s2e = bitr(join_gene$symbol, 
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = org.Hs.eg.db)
nrow(join_gene) 
join_gene = inner_join(join_gene,s2e,by=c("symbol"="SYMBOL"))
nrow(join_gene)
gene_rich = join_gene$ENTREZID 

# GO enrichment analysis
ego_CC <- enrichGO(gene = gene_rich,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene_rich,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene_rich,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

display_number = c(10, 10, 10)
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
  pvalue = c(ego_result_BP$pvalue,ego_result_CC$pvalue,ego_result_MF$pvalue),
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("BP", display_number[1]), 
                rep("CC", display_number[2]),
                rep("MF", display_number[3])), 
              levels=c("BP", "CC","MF" )))
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
write.csv(go_enrich_df,file = "GO.csv")
COLS <- c("#ff642a", "#009c66", "#33337a")

P1 = ggplot(data=go_enrich_df, aes(x=type_order,y=-log10(pvalue), fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Enrichment Score") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()

P2 = P1 + theme(axis.text = element_text(size = 12,color="black"),
                panel.border = element_rect(color = "black", linewidth = 1, fill = NA)) 
P3 = P2 + scale_x_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(P3,file="GO.png",width = 12,height = 12)

# KEGG enrichment analysis
kk <- enrichKEGG(gene = gene_rich,keyType = "kegg",organism= "human", qvalueCutoff = 0.05, pvalueCutoff=0.05)
gg1 <- as.data.frame(kk)
rownames(gg1) <- 1:nrow(gg1)
gg1$order=factor(rev(as.integer(rownames(gg1))),labels = rev(gg1$Description))
G1 = ggplot(gg1,aes(y=order,x=-log10(pvalue)))+
  geom_point(aes(size=Count,color=-1*pvalue))+
  scale_color_gradient(low="green",high = "red")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Enrichment Score",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()
G2 = G1 +theme(axis.text = element_text(size = 12,color="black"),
               panel.border = element_rect(color = "black", linewidth = 1, fill = NA)) 
ggsave(G2,file="KEGG.png",width = 8,height = 6)
pathway = data.frame(pathnames = gg1$Description,
                     symbol = gg1$geneID)
write.csv(pathway,file="pathway.csv",row.names = F)
write.csv(join_gene,file="ENTREZID.csv",row.names = F)
