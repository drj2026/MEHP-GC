rm(list = ls()) 
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

x = t(exp)
y = ifelse(Group=="Normal",0,1)
library(glmnet)
set.seed(1001)
cv_fit <- cv.glmnet(x=x, y=y)
png("cv_fit.png", width = 2000, height = 1200,res = 300)
plot(cv_fit)
dev.off()
fit <- glmnet(x=x, y=y)
png("fit.png", width = 2000, height = 1200,res = 300)
plot(fit,xvar = "lambda")
dev.off()
model <- glmnet(x=x, y=y,lambda=cv_fit$lambda.min)

lassoGene=rownames(model$beta)[as.numeric(model$beta)!=0]
length(lassoGene)
write.csv(lassoGene,"lassoGene.csv")