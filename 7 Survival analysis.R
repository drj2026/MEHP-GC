library(survival)
library(survminer)
load("TCGA-GC_sur_model.Rdata")
ls()
exprSet[1:4,1:4]
str(meta)
g = "EDNRA"
meta$gene2 = exprSet[g,]
cut = surv_cutpoint(meta, time = "time", event = "event", 
                    variables = "gene2")
m = cut[["cutpoint"]][1, 1]
meta$gene = ifelse(exprSet[g,]>m,'EDNRA-high','EDNRA-low')
table(meta$gene)
sfit=survfit(Surv(time, event)~gene, data=meta)
P = ggsurvplot(sfit,pval =TRUE, data = meta, risk.table = F,
               xlab = "Time(month)",legend.labs = c('EDNRA-high','EDNRA-low'),
               conf.int = F,
               palette = c("#b02113", "#016074"))
ggsave(filename = "KM_EDNRA.png", plot = P$plot, width = 6, height = 4, dpi = 300)
dev.off()
