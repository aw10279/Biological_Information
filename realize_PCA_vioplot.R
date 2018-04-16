#清除变量空间
rm(list=ls())

##############################################################
#lesson 2: realize PCA & vioplot
##############################################################

#princomp()  rstudio自带的主成分分析函数，功能有限制
#install.packages("psych")
library(psych)

combine_fpkm_table = read.table(file="./nature_2014_single_seq/fpkm_result/E18_combine_fpkm.table",header = T,sep = "\t")
dim(combine_fpkm_table)  #dim计算行列的数量

input_matrix = combine_fpkm_table[,c(2:ncol(combine_fpkm_table))]

#进行PCA分析
pca_result = principal(input_matrix,nfactors = 3)  #nfactors表示要取几个主成分

# PCA分析的各种结果，存在不同的变量里面

pca_result$values  ##各主成分数值
pca_result$scores  #各主成分得分
pca_result$weights  #各主成分权重

# 直接使用PCA的结果，进行绘图
plot(pca_result$scores[,1],pca_result$scores[,2],xlim=c(0,50),ylim=c(0,50))  #xlim和ylim限制图表的两轴范围

#vioplot&boxplot   vioplot可以看出数据在范围内的具体分布，boxplot只是展示数据的平均范围
gene_fpkm_table = read.table(file="~/Desktop/share/RNA_Seq_learn/cufflinks/hela_ctrl_rep1/genes.fpkm_tracking",header = T,sep = "\t")

# 筛选表达的基因画boxplot
select_vector = gene_fpkm_table$FPKM > 0
gene_fpkm_table.select = gene_fpkm_table[select_vector,]
boxplot(log2(gene_fpkm_table.select$FPKM + 1),col="orange")

# 使用vioplot库，画violin plot
#install.packages("vioplot")
library(vioplot)
vioplot(log2(gene_fpkm_table.select$FPKM + 1),col="orange")