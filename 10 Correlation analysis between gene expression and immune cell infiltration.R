library(ggcorrplot)
library(tidyr)

DEG_expr <- read.csv("DEG_expr.csv",row.names = 1)
cibersort_result <- read.csv("cibersort_result.csv",row.names = 1)
if (max(DEG_expr) < 50) {
  DEG_expr <- 2^DEG_expr - 1   
}
data <- t(DEG_expr)
cibersort_result <- cibersort_result[,1:22]
identical(rownames(data),rownames(cibersort_result))

data <- data[,c("BCL2L1","EDNRA")]
data <- as.data.frame(data)


calculate_correlation <- function(Gene_expr, cibersort_result) {
  cor_matrix <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:ncol(cibersort_result)) {
    result <- cor.test(Gene_expr, cibersort_result[, i], method = "pearson")
    new_row <- data.frame("Gene" = "Gene", "im_cell" = colnames(cibersort_result)[i], "Cor" = result$estimate, "p-value" = result$p.value)
    cor_matrix <- rbind(cor_matrix, new_row)
  }
  return(cor_matrix)
}


data_to_calculate <- data[,1:2]

results <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)

for (i in 1:ncol(data_to_calculate)) {
  print(i)
  gene_expr <- data_to_calculate[, i]
  corr_result <- calculate_correlation(gene_expr, cibersort_result)
  
  results <- rbind(results, corr_result)
}

colnames(data_to_calculate)

results$Gene <- c(rep("BCL2L1", 22),
                  rep("EDNRA", 22))
write.csv(results,file = "correlation_result.csv")

library(ggplot2)
data <- read.csv("correlation_result.csv", row.names = 1,header = T, sep = ",", stringsAsFactors = F)

target_gene <- "EDNRA"
mydata <- data[data$Gene == target_gene,]# 选择目标基因的数据
head(mydata)

mydata$im_cell <- reorder(mydata$im_cell, mydata$Cor)

plot <- ggplot(mydata, aes(x=Cor, y=im_cell)) +   # 绘制散点图
  geom_segment(aes(x = 0, xend = Cor, y = im_cell, yend = im_cell), color="black") + 
  geom_point(aes(size=abs(Cor), colour=p.value), alpha=0.5) +  
  
  scale_colour_gradient(low = "#016074", high = "#b02113") + 
  scale_size_continuous(range = c(2, 10)) + 
  theme_minimal() +  
  theme(axis.line = element_line(size = 1.0),
        axis.text = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(hjust = 1),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  
        axis.title = element_text(size = 14, face = "bold")) +
  labs(title = "EDNRA cor", xlab="Cor", ylab="im_cell")
plot + theme_grey() + ggtitle("theme_grey()")
ggsave("EDNRA.png", width = 6, height = 6, dpi = 300)