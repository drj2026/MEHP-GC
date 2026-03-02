library(tidyverse)
library(glmnet)
source('msvmRFE.R') 
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)
library(ggplot2)


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
Group = ifelse(Group=="Normal",0,1)
input = cbind(Group,exp)

svmRFE(input, k = 5, halve.above = 100)
nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100)

top.features = WriteFeatures(results, input, save=F)
head(top.features)
write.csv(top.features,"feature_svm.csv")


featsweep = lapply(1:37, FeatSweep.wrap, results, input) 
save(featsweep,file = "result.RData")
load("result.RData")
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
PlotErrors(errors, no.info=no.info)
Plotaccuracy(1-errors,no.info=no.info) 
which.min(errors)
top<-top.features[1:which.min(errors), "FeatureName"]
write.csv(top,"SVM key genes.csv")



