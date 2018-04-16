#清除变量空间
rm(list=ls())

##############################################################
#lesson 1: analyse .diff by boxplot and pheatmap
##############################################################

diff = read.table(file = "~/Desktop/share/RNA_Seq_learn/cuffdiff/gene_exp.diff",header = T,sep = "\t") #header表示是否有表头，sep表示分隔符是什么
filter_vector_p_value = diff$p_value < 0.05      #三种过滤条件筛选合适的数据
filter_vector_fc = abs(diff$log2.fold_change.) > 1
filter_vector_fpkm = diff$value_1>1 | diff$value_2>1

filter_vector_all = filter_vector_p_value & filter_vector_fc & filter_vector_fpkm  #条件汇总

diff.chosen = diff[filter_vector_all,]  #生成新表

boxplot(log2(diff.chosen$value_1+1), log2(diff.chosen$value_2+1),col = c("pink","orange"))#log将数据范围变小，图像更容易展现


install.packages("pheatmap")
library(pheatmap)
pheatmap(log2(diff.chosen[,c(8,9)]+1))