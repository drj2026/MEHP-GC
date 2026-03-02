proj = "RF"
load(file = "step2output.Rdata")
ls()
exp[1:4,1:4]
library(dplyr)
exp = as.data.frame(exp)
exp = mutate(exp,probe_id = rownames(exp))
ids = distinct(ids,symbol,.keep_all = T)
exp = inner_join(exp,ids,by="probe_id")
rownames(exp)=exp$symbol
exp = select(exp,-c(probe_id,symbol))
exp = as.matrix(exp)
table(Group)
G=read.table("join.txt")
g = G$V1
range(exp)
exp = exp[g,]
exp = t(exp)

library(randomForest)
library(ROCR)
library(genefilter)
library(Hmisc)
library(limma)
library(ggpubr)
set.seed(1112)

RF = randomForest(Group~.,data = exp,ntree = 500)
png("RFpicture.png", width = 1200, height = 1200,res = 300)
plot(RF,main = "Random forest",lwd = 2)
dev.off()
optionTrees = which.min(RF$err.rate[,1])
optionTrees
RF2 = randomForest(Group~.,data = exp,ntree = optionTrees)
plot(RF2,main = "Random forest",lwd = 2)
importance = importance(x = RF2)
importance = as.data.frame(importance)
importance$size = rownames(importance)
importance = importance[,c(2,1)]
names(importance) = c("Gene","importance")
af = importance[order(importance$importance,decreasing = T),]
af = af[1:10,]

p1 = ggplot(af,aes(x = reorder(Gene,importance),
                   y = importance,fill = importance))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_fill_gradient(low = ggsci::pal_npg()(2)[1],high =ggsci::pal_npg()(2)[2])+
  labs(x = "Gene",y = "Importance",title = "Top 10 Genes by Importance")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p1
png("RFimportance.png", width = 1200, height = 1200,res = 300)
plot(p1)
dev.off()
RFgene = importance[order(importance$importance,decreasing = T),]
write.csv(RFgene,"RFgene.csv")