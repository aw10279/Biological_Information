#清除变量空间
rm(list=ls())

##############################################################
#lesson4：volcano & Hi-C heatmap by hand
##############################################################


#volcano plot

diff = read.table(file = "~/Desktop/share/R/20171203-Live-R_partII/20171203-Live-R_partII/data_file/gene_exp.diff",header = T)

fc = log2(diff$value_1 / diff$value_2) #计算横坐标初始数据
#log2_foldchange中的INF、-INF、NAN分别为正无穷（正数/0）、负无穷（负数/0 或 0/非0）、非数字（0/0）。
#所以直接用变量==0的方式筛选出这些数并直接赋值0
fc[diff$value_1 == 0] = 0
fc[diff$value_2 == 0] = 0
p_value = -1 * log10(diff$p_value) #计算纵坐标初始数据


#创建筛选器，将所有判断显著表达的基本要求用&连接，得到一个布尔矩阵
exp_filter = abs(fc)>=1 & diff$p_value<=0.05 & (diff$value_1>0)&(diff$value_2>0) & (diff$value_1>=1|diff$value_2>=1)

#用同一个过滤条件去掉0值数据，保证3种数据的数量完全一致，可一一对应，并创建颜色筛选器（也是布尔矩阵）
p_value_filter = p_value[p_value >= 0.001]
fc_filter = fc[p_value >= 0.001]
col_filter = exp_filter[p_value >= 0.001]

col_vector = rep(rgb(1,0,0,0.1), length(p_value_filter))#创建颜色向量，数量=p_value_filter
col_vector[col_filter] = rgb(0,0,1) #用颜色筛选器把筛出的点赋值另一种颜色
plot(x = fc_filter, y = p_value_filter, xlim = c(-10,10),ylim = c(0,4),col = col_vector, pch = 16)

abline(h=-1*log10(0.05),lwd=1,lty=6,col="#4C5B61")  #画出p_value值的基准线


#基础画图工具手绘heatmap plot
rm(list = ls())
input_matrix= matrix(c(1:30),5,6)
x_size = dim(input_matrix)[1]
y_size = dim(input_matrix)[2] #x和y的大小可以直接由输入矩阵的尺寸赋予

#用rect画方格，每个格子都由xleft、xright、ybottom、ytop四条直线决定。因为数值从上到下对应从小到大，故：
# 1st column
# xleft 0,0,0,0,0
# ybottom 4,3,2,1,0
# xright 1,1,1,1,1
# ytop 5,4,3,2,1

# 2nd column
# xleft 1,1,1,1,1
# ybottom 4,3,2,1,0
# xright 2,2,2,2,2
# ytop 5,4,3,2,1
#故：
my_xleft = rep(c(0:(y_size-1)),each = x_size)
my_xright = my_xleft + 1
my_ybottom = rep(c((x_size-1):0),y_size)
my_ytop = my_ybottom + 1

#确定颜色及变化规则
matrix_max = max(input_matrix)
col_rate = input_matrix / matrix_max
is.matrix(col_rate)
is.vector(col_rate)
col_matrix = rgb(1,0,0,as.vector(col_rate))

#开始绘图
plot(x=c(0:y_size), y=c(0:y_size), type="n", frame.plot=F, xaxt="n", yaxt="n",xlab="",ylab="")#画出尺寸范围，清空所有内容
rect(xleft=my_xleft, xright=my_xright, ybottom=my_ybottom, ytop=my_ytop, col=col_matrix, border=F)#画出热图


#Hi-C数据画图
rm(list = ls())

hic.raw = read.table(file = "~/Desktop/share/R/20171203-Live-R_partII/20171203-Live-R_partII/data_file/chr_16_100000_MAPQ20.txt", sep = ",", header = F)
input_matrix = as.matrix(hic.raw) #将读取数据矩阵化，下面操作才能进行

x_size = dim(input_matrix)[1]
y_size = dim(input_matrix)[2] #获取图像尺寸

my_xleft = rep(c(0:(y_size-1)),each = x_size)
my_xright = my_xleft + 1
my_ybottom = rep(c((x_size-1):0),y_size)
my_ytop = my_ybottom + 1  #确定每个方格的4条线

matrix_max = quantile(input_matrix,prob = 0.9) #用quantile计算矩阵的90%位数，作为最大值
col_rate = input_matrix / matrix_max #用最大值求得透明度比例矩阵
col_rate[col_rate>1] = 1 #将其中大于1的都赋值1
col_matrix = rgb(1,0,0,col_rate) #用该矩阵得到颜色数组（注意col_matrix并不是矩阵，但排列顺序和画图顺序完全一致，故可用）

plot(x=c(0:y_size), y=c(0:y_size), type="n", frame.plot=F, xaxt="n", yaxt="n",xlab="",ylab="")
rect(xleft=my_xleft, xright=my_xright, ybottom=my_ybottom, ytop=my_ytop, col=col_matrix, border=NA)

# tiff(file="~/test_hic.png",width = 2000,height = 2000)可直接生成图片
# dev.off() 关闭函数