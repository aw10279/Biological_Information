#清除变量空间
rm(list=ls())

##############################################################
#lesson3：Basic R & statistical analysis & kinds of ploting
##############################################################

## 常用的赋值操作
a <- 1
2 -> a
a <- 1e5
a <- 1e-10
a <- log2(10)
a <- log10(10)
a <- log(10)
a <- pi
a <- exp(1) 

## 常用的方法生成向量
b <- c(1:10)
a <- 1 
c <- a + b

## 创建序列 
d <- seq(from=1,to=100,by=5)
e <- rep(c(1:10),5)
f <- rep(c(1:10),each=5)

## 用c函数可以合并2个向量
a = c(1:10)
b = rep(c(1:10),each=5)
c = rep(c(1:5),2)
d = c(a,b,c)

## 创建矩阵
g <- matrix(c(1:100), nrow = 10)
h <- matrix(c(1:15),5,5)

# 查看变量空间中都有哪些变量
ls()

g_up <- upper.tri(g, diag = T) #上三角全部置TRUE，得到布尔矩阵,如要留对角线：diag = T
g_low <- lower.tri(g) #下三角全部置TRUE，得到布尔矩阵
g[upper.tri(g)] #根据布尔矩阵得到上三角的具体数值
det(h) # 求方阵行列式
t(g) # 转置
eigen(g)  #特征值

## R中的for循环
d <- seq(from=1,to=100,by=5)
for(i in d){
  if(i >= 50){
    print(i)
  }
}

d[d>50 & d<80] #数组中可以直接用条件筛选

#对于数据框，也同样可以做筛选
gene_table
gene_table.select = gene_table[gene_table$gene_fpkm >= 10,]

## R中写自己的函数
my_ABS <- function(x){  #用function()直接给变量名赋值即可
  if(x <= 0){
    -x
  }else{
    x
  }
}

my_ABS(10)
my_ABS(-10)G852018031915561957204

?hist  #  ?+函数名可查使用参数，example(函数名)可查看实例

## 创建data.frame
gene_id = c(1:100)
gene_fpkm = abs(rnorm(100,10,5)) #rnorm 符合正态分布的随机数100个，均值为10，方差为5.可能出现负值，再取绝对值。
sample_id = round(runif(100,1,10))#round为取整，runif(100,1,10) 从1-10中平均地取100个随机数
gene_table = data.frame(gene_id,gene_fpkm,sample_id) #拼装成数表
colnames(gene_table) 
rownames(gene_table) #显示行名和列名
dim(gene_table) #显示表的大小
table(gene_table$sample_id) #对gene_id这一列进行数值统计
#gene_table$gene_id和gene_table[0,1]等价，都是指某一列，用列名的好处是不用专门数列数

barplot(table(gene_table$sample_id)) #barplot对sample_id的统计情况画直方图
hist(gene_table$sample_id) #实际所画图像和上面的barplot基本差不多，因为所实现功能基本相同
hist(gene_table$gene_fpkm,col = "#00A8E8", border = F) #hist 显示fpkm的分布情况统计并画直方图,border为是否要边线，颜色可自定义
abline(v = mean(gene_table$gene_fpkm),col = "#003459", lwd = 3,lty=3)
abline(v = median(gene_table$gene_fpkm))
#abline为画线工具，v表示竖线h表示横线，mean为平均数，median为中位数，lwd为线宽，lty为线的种类。

par(new = True)  #表示在原图上画新的图
plot(density(gene_table$gene_fpkm),xaxt="n",yaxt="n",bty="n") #density求出fpkm值(x)的概率分布(y),再用plot画出拟合曲线。xaxt和yaxt表示画不画横纵坐标，bty表示画框
plot(x=gene_table$gene_id,y=gene_table$gene_fpkm,type = "o")




############# R与统计分析 #############
case_1 <- rnorm(50,mean = 20,sd=5) #生成正态分布样本1，mean平均数，sd标准差
plot(case_1,main="case_1")
case_2 <- rnorm(40,mean = 25,sd=5) #生成正态分布样本2
plot(case_2,main="case_2")

hist(case_1,breaks = c(1:40)) #分布图
hist(case_2,breaks = c(1:40))

case_T.test <- t.test(case_1,case_2) #独立样本T检验
case_T.test$p.value


############# 使用R画boxplot #############
case_1 <- rnorm(10,5,1)
case_2 <- rnorm(10,6,1)
case_3 <- rnorm(10,7,1)
case_4 <- rnorm(10,8,1)
case_all <- cbind(case_1,case_2,case_3,case_4)  #cbind 将几组数组合成数表
boxplot(case_all)

col_1 <- "red"
col_2 <- rainbow(10)[5]
col_3 <- rgb(1,0,0,alpha = 0.5)
col_4 <- rgb(0,0,1,alpha = 0.5)
col_list <- c(col_1,col_2,col_3,col_4) #将上面4中颜色组合起来，一起使用在给case_all画图上，可以同时给几列数据颜色
boxplot(case_all,col = col_list)

col_list <- c(col_3,col_3,col_4,col_4)
boxplot(case_all,col = col_list)


############# 绘制RNA-Seq基因表达量的散点图 #############

cuffnorm_result = read.csv(file="~/Desktop/live_R_data/cuffnorm_genes_fpkm.csv")
head(cuffnorm_result)

x.vector = cuffnorm_result$Empty_KD_0
y.vector = cuffnorm_result$Empty_KD_1

plot(x=x.vector,y=y.vector,col="#0000FF11",pch=16,cex=0.5,xlim=c(0,10),ylim=c(0,10),xlab = "Repeat_1 Log2(FPKM)",ylab = "Repeat_2 Log2(FPKM)")
# plot中的pch即是画的符号的样式，可以搜索R pch直接查各种符号的代码，cex为符号大小
text(2,8,round(cor(x.vector.filter,y.vector.filter,method = "spearman"),4))
# cor是计算2组数据的相关性系数（原数据可取log），method是使用的统计学方法，text是将传入内容写在图表上，前面的2和8是书写的具体坐标位置
hist(x.vector.filter[x.vector.filter>=0],breaks = seq(0,15,0.1),border = F,col = "blue",xlab = "FPKM",ylab = "Frequency",main = "")
hist(y.vector.filter[y.vector.filter>=0],breaks = seq(0,15,0.1),border = F,col = "blue",xlab = "FPKM",ylab = "Frequency",main = "")