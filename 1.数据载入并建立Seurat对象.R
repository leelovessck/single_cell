rm(list = ls())
setwd("E:/2023.10/PE scRNA（多数据库运行）")
getwd()

#打开必要的package
{
if(!require(multtest))BiocManager::install("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(tidyverse))install.packages("tidyverse")
}

#原始数据的答疑
{
#注1：
#如果下面涉及到原始数据不清楚从哪里来的，或者不了解分组情况
#参考“下载的原始数据”文件夹中的“使用的数据-总结.xlsx”

#注2：
#下列原始数据中有以下3种存储形式：
#（1）barcodes、features、matrix3个文件；
#（2）稀疏矩阵（txt或tsv）；
#（3）Rdata；
#分别用不同的方法载入后才能创建Seurat文件，标题中将体现原文件形式
#还有h5、h5ad等文件形式，这两种都是来自Python的对象，下面没有涉及

#注3
#rds文件是二进制文件
#在新的项目中，载入rds文件可以得到相同的对象
#例：保存Seurat的rds，再次载入可以在新的R里得直接到Seurat对象
#例：保存CellChat的rds，再次载入可以在新的R里直接得到CellChat对象

#注4
#在Seurat对象中，metadata可以理解成一个Excel表格
#一种修改方式：“Seurat_name@meta.data$orig.ident <- xxx”
#另外一种修改方式：
#先用“metadata <- Seurat_name@meta.data”提取出metadata
#然后按矩阵修改，或者保存导出后在excel中修改
#然后再用Seurat_name@meta.data <- metadata把修改好的metadata导入
  
#注5
#如果是表达矩阵形式的文件（txt或tsv），使用“read.table”函数读取
#如何同时读取一个目录下的多个表达矩阵，可能需要使用循环语句（技能待学习）
#下面给出的方法可能有点繁琐，但是能用
  
#注6
#GSE173193给出的是同时读取同一目录下多个cellranger输出文件的方法
#如果只有一个文件，可以简单通过以下两行代码完成
#in[1]：counts1 <- Read10X(data.dir = "D:/XX/文件目录/")  ##最后一个/必不可少
#in[2]：scRNA1 = CreateSeuratObject(counts1)
  
#注7
#另外还有文件格式为h5、h5ad格式的文件
#这种文件格式是从python导出的单细胞对象
#下面是h5文件的载入方法（使用hdf5r包）
#in[1]：devtools::install_github(repo = "hhoeflin/hdf5r")
#in[2]：library(hdf5r)
#in[3]：scRNA1 <-Read10X_h5("xxx.h5")
#下面是h5ad文件的载入方法（使用SeuratDisk包）
#in[1]：remotes::install_github("mojaveazure/seurat-disk", force = TRUE)
#in[2]：library(SeuratDisk)
#in[3]：setwd('D:/')  ##路径要是英文
#in[4]：Convert("data.h5ad", "h5seurat",overwrite = TRUE, assay = "RNA")
#in[5]：scRNA1 <- LoadH5Seurat("data.h5seurat")
}

#载入GSE173193（barcodes、features、matrix3个文件）
{
samples1 <- list.files("raw_data/GSE173193")  ##列出路径中的样本
samples1
dir1 <- file.path("raw_data/GSE173193",samples1)  ##得到样本路径
names(dir1) <- samples1
dir1
counts1 <- Read10X(data.dir = dir1)
scRNA1 = CreateSeuratObject(counts1)
unique(scRNA1$orig.ident)  ##变量就存在这里
  ##得到了control1-2、EOPE1-2的Seurat对象
saveRDS(scRNA1,"./Seurat_data/con1-2和EOPE1-2.rds")  ##保存rds文件
}

#载入GSE87720（稀疏矩阵，txt）
{
#首先是control3
{
counts2 <- read.table("./raw_data/GSE87720/control3/GSE87720_Placenta2.txt",
                      sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
counts2 <- counts2 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(counts2) <- counts2[,1]
  ##把gene设置为列名
counts2 <- counts2[,-1]
  ##删除数据集中gene这一行
  ##这里要保证，除了列名、行名外，数据集中全部是数字
counts2 <- as(as.matrix(counts2), "dgCMatrix")
scRNA2 <- CreateSeuratObject(counts = counts2)
unique(scRNA2$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNA2$orig.ident <- "control3"
unique(scRNA2$orig.ident)  ##变量就存在这里
  ##得到了control3的Seurat对象
saveRDS(scRNA2,"./Seurat_data/con3.rds")  ##保存rds文件
}
  
#下面是control4
{
counts3 <- read.table("./raw_data/GSE87720/control4/GSE87720_Placenta1.txt",sep = "\t",header = T)
  ##这里发现有重复的基因
counts3 <- counts3 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(counts3) <- counts3[,1]
  ##把gene设置为列名
counts3 <- counts3[,-1]
  ##删除数据集中gene这一行
counts3 <- as(as.matrix(counts3), "dgCMatrix")
scRNA3 <- CreateSeuratObject(counts = counts3)
unique(scRNA3$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNA3$orig.ident <- "control4"
unique(scRNA3$orig.ident)  ##变量就存在这里
  ##得到了control3的Seurat对象
saveRDS(scRNA3,"./Seurat_data/con4.rds")  ##保存rds文件
}
}

#载入GSE182381（Rda）
{
scRNA4 <- readRDS(file = "./raw_data/GSE182381/GSE182381.rda")
  ##这里可以直接得到Seurat对象
  ##原文里面的分组情况也在数据集中,下面把不需要的数据删除
metadata4 <- scRNA4@meta.data
metadata4 <- metadata4[,1:3]
scRNA4@meta.data <- metadata4
  ##删除了原作者的后续分析内容
unique(scRNA4@meta.data$orig.ident)
  ##这里发现运行结果和实验目的不符，下面手动修改
unique(metadata4[,1])
metadata4$orig.ident[metadata4$orig.ident == "kc.40.1"] <- "control5"
metadata4$orig.ident[metadata4$orig.ident == "kc.40.2"] <- "control6"
metadata4$orig.ident[metadata4$orig.ident == "kc.42.1"] <- "control7"
metadata4$orig.ident[metadata4$orig.ident == "kc.42.2"] <- "control8"
unique(metadata4[,1])
scRNA4@meta.data <- metadata4
unique(scRNA4@meta.data$orig.ident)  ##变量就存在这里
  ##得到了control5-8的Seurat对象
saveRDS(scRNA4,"./Seurat_data/con5-8.rds")  ##保存rds文件
}

#载入github数据（稀疏矩阵，txt）
{
counts5 <- read.table("./raw_data/github/control9/control9.txt",
                      sep = "\t",header = T,row.names = 1)
  ##这里没有重复的基因
counts5 <- as(as.matrix(counts5), "dgCMatrix")
scRNA5 <- CreateSeuratObject(counts = counts5)
unique(scRNA5$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNA5$orig.ident <- "control9"
unique(scRNA5$orig.ident)  ##变量就存在这里
  ##得到了control9的Seurat对象
saveRDS(scRNA5,"./Seurat_data/con9.rds")  ##保存rds文件

#control10
counts6 <- read.table("./raw_data/github/control10/control10.txt",
                      sep = "\t",header = T,row.names = 1)
##这里没有重复的基因
counts6 <- as(as.matrix(counts6), "dgCMatrix")
scRNA6 <- CreateSeuratObject(counts = counts6)
unique(scRNA6$orig.ident)  
##这里发现运行结果和实验目的不符，下面手动修改
scRNA6$orig.ident <- "control10"
unique(scRNA6$orig.ident)  ##变量就存在这里
##得到了control10的Seurat对象
saveRDS(scRNA6,"./Seurat_data/con10.rds")  ##保存rds文件

#control11
counts7 <- read.table("./raw_data/github/control11/control11.txt",
                      sep = "\t",header = T,row.names = 1)
##这里没有重复的基因
counts7 <- as(as.matrix(counts7), "dgCMatrix")
scRNA7 <- CreateSeuratObject(counts = counts7)
unique(scRNA7$orig.ident)  
##这里发现运行结果和实验目的不符，下面手动修改
scRNA7$orig.ident <- "control11"
unique(scRNA7$orig.ident)  ##变量就存在这里
##得到了control11的Seurat对象
saveRDS(scRNA7,"./Seurat_data/con11.rds")  ##保存rds文件

#LOPE1
counts8 <- read.table("./raw_data/github/LOPE1/LOPE1.txt",
                      sep = "\t",header = T,row.names = 1)
##这里没有重复的基因
counts8 <- as(as.matrix(counts8), "dgCMatrix")
scRNA8 <- CreateSeuratObject(counts = counts8)
unique(scRNA8$orig.ident)  
##这里发现运行结果和实验目的不符，下面手动修改
scRNA8$orig.ident <- "LOPE1"
unique(scRNA8$orig.ident)  ##变量就存在这里
##得到了LOPE1的Seurat对象
saveRDS(scRNA8,"./Seurat_data/LOPE1.rds")  ##保存rds文件

#LOPE2
counts9 <- read.table("./raw_data/github/LOPE2/LOPE2.txt",
                      sep = "\t",header = T,row.names = 1)
##这里没有重复的基因
counts9 <- as(as.matrix(counts9), "dgCMatrix")
scRNA9 <- CreateSeuratObject(counts = counts9)
unique(scRNA9$orig.ident)  
##这里发现运行结果和实验目的不符，下面手动修改
scRNA9$orig.ident <- "LOPE2"
unique(scRNA9$orig.ident)  ##变量就存在这里
##得到了LOPE2的Seurat对象
saveRDS(scRNA9,"./Seurat_data/LOPE2.rds")  ##保存rds文件

#LOPE3
counts10 <- read.table("./raw_data/github/LOPE3/LOPE3.txt",
                      sep = "\t",header = T,row.names = 1)
##这里没有重复的基因
counts10 <- as(as.matrix(counts10), "dgCMatrix")
scRNA10 <- CreateSeuratObject(counts = counts10)
unique(scRNA10$orig.ident)  
##这里发现运行结果和实验目的不符，下面手动修改
scRNA10$orig.ident <- "LOPE3"
unique(scRNA10$orig.ident)  ##变量就存在这里
##得到了LOPE3的Seurat对象
saveRDS(scRNA10,"./Seurat_data/LOPE3.rds")  ##保存rds文件

#LOPE4
counts11 <- read.table("./raw_data/github/LOPE4/LOPE4.txt",
                       sep = "\t",header = T,row.names = 1)
##这里没有重复的基因
counts11 <- as(as.matrix(counts11), "dgCMatrix")
scRNA11 <- CreateSeuratObject(counts = counts11)
unique(scRNA11$orig.ident)  
##这里发现运行结果和实验目的不符，下面手动修改
scRNA11$orig.ident <- "LOPE4"
unique(scRNA11$orig.ident)  ##变量就存在这里
##得到了LOPE3的Seurat对象
saveRDS(scRNA11,"./Seurat_data/LOPE4.rds")  ##保存rds文件
}

#载入GSE192693（barcodes、features、matrix3个文件）
{
samples2 <- list.files("raw_data/GSE192693")  ##列出路径中的样本
samples2
dir2 <- file.path("raw_data/GSE192693",samples2)  ##得到样本路径
names(dir2) <- samples2
dir2
counts12 <- Read10X(data.dir = dir2)
scRNA12 = CreateSeuratObject(counts12)
unique(scRNA12$orig.ident)  ##变量就存在这里
  ##得到了control12-15、EOPE3-8的Seurat对象
saveRDS(scRNA12,"./Seurat_data/con12-15和EOPE3-8.rds")  ##保存rds文件
}

#载入文章自带large data（稀疏矩阵，txt）
{
#下面使用的文件已经是经python拆分后的文件
#首先是control16
{
counts13 <- read.table("./raw_data/large data/control16/control16.txt",
                      sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(counts13)[1] <- "Gene"
counts13 <- counts13 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(counts13) <- counts13[,1]
  ##把gene设置为列名
counts13 <- counts13[,-1]
  ##删除数据集中gene这一行
counts13 <- counts13[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
class(counts13$AAACCCAAGCAGTAAT.1_48.1)
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(counts13,"./raw_data/large data（处理）/control16_num.csv")
counts13_num <- read.csv("./raw_data/large data（处理）/control16_num.csv")
class(counts13_num$AAACCCAAGCAGTAAT.1_48.1)
rownames(counts13_num) <- counts13_num[,1]
counts13_num <- counts13_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
counts13_num <- as(as.matrix(counts13_num), "dgCMatrix")
scRNA13 <- CreateSeuratObject(counts = counts13_num)
unique(scRNA13$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNA13$orig.ident <- "control16"
unique(scRNA13$orig.ident)  ##变量就存在这里
  ##得到了control16的Seurat对象
saveRDS(scRNA13,"./Seurat_data/con16.rds")  ##保存rds文件
}

#这里是control17
{
counts17 <- read.table("./raw_data/large data/control17/control17.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(counts17)[1] <- "Gene"
counts17 <- counts17 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(counts17) <- counts17[,1]
  ##把gene设置为列名
counts17 <- counts17[,-1]
  ##删除数据集中gene这一行
counts17 <- counts17[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
class(counts17$AAACCCAAGTGACACG.1_53.1)
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(counts17,"./raw_data/large data（处理）/control17_num.csv")
counts17_num <- read.csv("./raw_data/large data（处理）/control17_num.csv")
class(counts17_num$AAACCCAAGTGACACG.1_53.1)
rownames(counts17_num) <- counts17_num[,1]
counts17_num <- counts17_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
counts17_num <- as(as.matrix(counts17_num), "dgCMatrix")
scRNA17 <- CreateSeuratObject(counts = counts17_num)
unique(scRNA17$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNA17$orig.ident <- "control17"
unique(scRNA17$orig.ident)  ##变量就存在这里
  ##得到了control17的Seurat对象
saveRDS(scRNA17,"./Seurat_data/con17.rds")  ##保存rds文件
}

#这里是control18
{
counts18 <- read.table("./raw_data/large data/control18/control18.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(counts18)[1] <- "Gene"
counts18 <- counts18 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(counts18) <- counts18[,1]
  ##把gene设置为列名
counts18 <- counts18[,-1]
  ##删除数据集中gene这一行
counts18 <- counts18[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(counts18,"./raw_data/large data（处理）/control18_num.csv")
counts18_num <- read.csv("./raw_data/large data（处理）/control18_num.csv")
rownames(counts18_num) <- counts18_num[,1]
counts18_num <- counts18_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
counts18_num <- as(as.matrix(counts18_num), "dgCMatrix")
scRNA18 <- CreateSeuratObject(counts = counts18_num)
unique(scRNA18$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNA18$orig.ident <- "control18"
unique(scRNA18$orig.ident)  ##变量就存在这里
  ##得到了control17的Seurat对象
saveRDS(scRNA18,"./Seurat_data/con18.rds")  ##保存rds文件
}

#这里是control19
{
counts19 <- read.table("./raw_data/large data/control19/control19.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(counts19)[1] <- "Gene"
counts19 <- counts19 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(counts19) <- counts19[,1]
  ##把gene设置为列名
counts19 <- counts19[,-1]
  ##删除数据集中gene这一行
counts19 <- counts19[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(counts19,"./raw_data/large data（处理）/control19_num.csv")
counts19_num <- read.csv("./raw_data/large data（处理）/control19_num.csv")
rownames(counts19_num) <- counts19_num[,1]
counts19_num <- counts19_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
counts19_num <- as(as.matrix(counts19_num), "dgCMatrix")
scRNA19 <- CreateSeuratObject(counts = counts19_num)
unique(scRNA19$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNA19$orig.ident <- "control19"
unique(scRNA19$orig.ident)  ##变量就存在这里
  ##得到了control19的Seurat对象
saveRDS(scRNA19,"./Seurat_data/con19.rds")  ##保存rds文件
}

#这里是control20
{
counts20 <- read.table("./raw_data/large data/control20/control20.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(counts20)[1] <- "Gene"
counts20 <- counts20 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(counts20) <- counts20[,1]
  ##把gene设置为列名
counts20 <- counts20[,-1]
  ##删除数据集中gene这一行
counts20 <- counts20[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(counts20,"./raw_data/large data（处理）/control20_num.csv")
counts20_num <- read.csv("./raw_data/large data（处理）/control20_num.csv")
rownames(counts20_num) <- counts20_num[,1]
counts20_num <- counts20_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
counts20_num <- as(as.matrix(counts20_num), "dgCMatrix")
scRNA20 <- CreateSeuratObject(counts = counts20_num)
unique(scRNA20$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNA20$orig.ident <- "control20"
unique(scRNA20$orig.ident)  ##变量就存在这里
  ##得到了control20的Seurat对象
saveRDS(scRNA20,"./Seurat_data/con20.rds")  ##保存rds文件
}

#这里是control21
{
counts21 <- read.table("./raw_data/large data/control21/control21.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(counts21)[1] <- "Gene"
counts21 <- counts21 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(counts21) <- counts21[,1]
  ##把gene设置为列名
counts21 <- counts21[,-1]
  ##删除数据集中gene这一行
counts21 <- counts21[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(counts21,"./raw_data/large data（处理）/control21_num.csv")
counts21_num <- read.csv("./raw_data/large data（处理）/control21_num.csv")
rownames(counts21_num) <- counts21_num[,1]
counts21_num <- counts21_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
counts21_num <- as(as.matrix(counts21_num), "dgCMatrix")
scRNA21 <- CreateSeuratObject(counts = counts21_num)
unique(scRNA21$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNA21$orig.ident <- "control21"
unique(scRNA21$orig.ident)  ##变量就存在这里
  ##得到了control21的Seurat对象
saveRDS(scRNA21,"./Seurat_data/con21.rds")  ##保存rds文件
}

#这里是control22
{
counts22 <- read.table("./raw_data/large data/control22/control22.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(counts22)[1] <- "Gene"
counts22 <- counts22 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(counts22) <- counts22[,1]
  ##把gene设置为列名
counts22 <- counts22[,-1]
  ##删除数据集中gene这一行
counts22 <- counts22[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(counts22,"./raw_data/large data（处理）/control22_num.csv")
counts22_num <- read.csv("./raw_data/large data（处理）/control22_num.csv")
rownames(counts22_num) <- counts22_num[,1]
counts22_num <- counts22_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
counts22_num <- as(as.matrix(counts22_num), "dgCMatrix")
scRNA22 <- CreateSeuratObject(counts = counts22_num)
unique(scRNA22$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNA22$orig.ident <- "control22"
unique(scRNA22$orig.ident)  ##变量就存在这里
  ##得到了control22的Seurat对象
saveRDS(scRNA22,"./Seurat_data/con22.rds")  ##保存rds文件
}

#这里是control23
{
counts23 <- read.table("./raw_data/large data/control23/control23.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(counts23)[1] <- "Gene"
counts23 <- counts23 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(counts23) <- counts23[,1]
  ##把gene设置为列名
counts23 <- counts23[,-1]
  ##删除数据集中gene这一行
counts23 <- counts23[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(counts23,"./raw_data/large data（处理）/control23_num.csv")
counts23_num <- read.csv("./raw_data/large data（处理）/control23_num.csv")
rownames(counts23_num) <- counts23_num[,1]
counts23_num <- counts23_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
counts23_num <- as(as.matrix(counts23_num), "dgCMatrix")
scRNA23 <- CreateSeuratObject(counts = counts23_num)
unique(scRNA23$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNA23$orig.ident <- "control23"
unique(scRNA23$orig.ident)  ##变量就存在这里
  ##得到了control23的Seurat对象
saveRDS(scRNA23,"./Seurat_data/con23.rds")  ##保存rds文件
}

#这里是control24
{
counts24 <- read.table("./raw_data/large data/control24/control24.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(counts24)[1] <- "Gene"
counts24 <- counts24 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(counts24) <- counts24[,1]
  ##把gene设置为列名
counts24 <- counts24[,-1]
  ##删除数据集中gene这一行
counts24 <- counts24[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(counts24,"./raw_data/large data（处理）/control24_num.csv")
counts24_num <- read.csv("./raw_data/large data（处理）/control24_num.csv")
rownames(counts24_num) <- counts24_num[,1]
counts24_num <- counts24_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
counts24_num <- as(as.matrix(counts24_num), "dgCMatrix")
scRNA24 <- CreateSeuratObject(counts = counts24_num)
unique(scRNA24$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNA24$orig.ident <- "control24"
unique(scRNA24$orig.ident)  ##变量就存在这里
  ##得到了control24的Seurat对象
saveRDS(scRNA24,"./Seurat_data/con24.rds")  ##保存rds文件
}

#这里是EOPE9
{
countsE9 <- read.table("./raw_data/large data/EOPE9/EOPE9.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsE9)[1] <- "Gene"
countsE9 <- countsE9 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsE9) <- countsE9[,1]
  ##把gene设置为列名
countsE9 <- countsE9[,-1]
  ##删除数据集中gene这一行
countsE9 <- countsE9[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsE9,"./raw_data/large data（处理）/EOPE9_num.csv")
countsE9_num <- read.csv("./raw_data/large data（处理）/EOPE9_num.csv")
rownames(countsE9_num) <- countsE9_num[,1]
countsE9_num <- countsE9_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsE9_num <- as(as.matrix(countsE9_num), "dgCMatrix")
scRNAE9 <- CreateSeuratObject(counts = countsE9_num)
unique(scRNAE9$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAE9$orig.ident <- "EOPE9"
unique(scRNAE9$orig.ident)  ##变量就存在这里
  ##得到了EOPE9的Seurat对象
saveRDS(scRNAE9,"./Seurat_data/EOPE9.rds")  ##保存rds文件
}

#这里是EOPE10
{
countsE10 <- read.table("./raw_data/large data/EOPE10/EOPE10.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsE10)[1] <- "Gene"
countsE10 <- countsE10 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsE10) <- countsE10[,1]
  ##把gene设置为列名
countsE10 <- countsE10[,-1]
  ##删除数据集中gene这一行
countsE10 <- countsE10[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsE10,"./raw_data/large data（处理）/EOPE10_num.csv")
countsE10_num <- read.csv("./raw_data/large data（处理）/EOPE10_num.csv")
rownames(countsE10_num) <- countsE10_num[,1]
countsE10_num <- countsE10_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsE10_num <- as(as.matrix(countsE10_num), "dgCMatrix")
scRNAE10 <- CreateSeuratObject(counts = countsE10_num)
unique(scRNAE10$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAE10$orig.ident <- "EOPE10"
unique(scRNAE10$orig.ident)  ##变量就存在这里
  ##得到了EOPE10的Seurat对象
saveRDS(scRNAE10,"./Seurat_data/EOPE10.rds")  ##保存rds文件
}

#这里是EOPE11
{
countsE11 <- read.table("./raw_data/large data/EOPE11/EOPE11.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsE11)[1] <- "Gene"
countsE11 <- countsE11 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsE11) <- countsE11[,1]
  ##把gene设置为列名
countsE11 <- countsE11[,-1]
  ##删除数据集中gene这一行
countsE11 <- countsE11[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsE11,"./raw_data/large data（处理）/EOPE11_num.csv")
countsE11_num <- read.csv("./raw_data/large data（处理）/EOPE11_num.csv")
rownames(countsE11_num) <- countsE11_num[,1]
countsE11_num <- countsE11_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsE11_num <- as(as.matrix(countsE11_num), "dgCMatrix")
scRNAE11 <- CreateSeuratObject(counts = countsE11_num)
unique(scRNAE11$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAE11$orig.ident <- "EOPE11"
unique(scRNAE11$orig.ident)  ##变量就存在这里
  ##得到了EOPE11的Seurat对象
saveRDS(scRNAE11,"./Seurat_data/EOPE11.rds")  ##保存rds文件
}

#这里是EOPE12
{
countsE12 <- read.table("./raw_data/large data/EOPE12/EOPE12.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsE12)[1] <- "Gene"
countsE12 <- countsE12 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsE12) <- countsE12[,1]
  ##把gene设置为列名
countsE12 <- countsE12[,-1]
  ##删除数据集中gene这一行
countsE12 <- countsE12[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsE12,"./raw_data/large data（处理）/EOPE12_num.csv")
countsE12_num <- read.csv("./raw_data/large data（处理）/EOPE12_num.csv")
rownames(countsE12_num) <- countsE12_num[,1]
countsE12_num <- countsE12_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsE12_num <- as(as.matrix(countsE12_num), "dgCMatrix")
scRNAE12 <- CreateSeuratObject(counts = countsE12_num)
unique(scRNAE12$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAE12$orig.ident <- "EOPE12"
unique(scRNAE12$orig.ident)  ##变量就存在这里
  ##得到了EOPE12的Seurat对象
saveRDS(scRNAE12,"./Seurat_data/EOPE12.rds")  ##保存rds文件
}

#这里是EOPE13
{
countsE13 <- read.table("./raw_data/large data/EOPE13/EOPE13.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsE13)[1] <- "Gene"
countsE13 <- countsE13 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsE13) <- countsE13[,1]
  ##把gene设置为列名
countsE13 <- countsE13[,-1]
  ##删除数据集中gene这一行
countsE13 <- countsE13[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsE13,"./raw_data/large data（处理）/EOPE13_num.csv")
countsE13_num <- read.csv("./raw_data/large data（处理）/EOPE13_num.csv")
rownames(countsE13_num) <- countsE13_num[,1]
countsE13_num <- countsE13_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsE13_num <- as(as.matrix(countsE13_num), "dgCMatrix")
scRNAE13 <- CreateSeuratObject(counts = countsE13_num)
unique(scRNAE13$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAE13$orig.ident <- "EOPE13"
unique(scRNAE13$orig.ident)  ##变量就存在这里
  ##得到了EOPE13的Seurat对象
saveRDS(scRNAE13,"./Seurat_data/EOPE13.rds")  ##保存rds文件
}

#这里是EOPE14
{
countsE14 <- read.table("./raw_data/large data/EOPE14/EOPE14.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsE14)[1] <- "Gene"
countsE14 <- countsE14 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsE14) <- countsE14[,1]
  ##把gene设置为列名
countsE14 <- countsE14[,-1]
  ##删除数据集中gene这一行
countsE14 <- countsE14[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsE14,"./raw_data/large data（处理）/EOPE14_num.csv")
countsE14_num <- read.csv("./raw_data/large data（处理）/EOPE14_num.csv")
rownames(countsE14_num) <- countsE14_num[,1]
countsE14_num <- countsE14_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsE14_num <- as(as.matrix(countsE14_num), "dgCMatrix")
scRNAE14 <- CreateSeuratObject(counts = countsE14_num)
unique(scRNAE14$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAE14$orig.ident <- "EOPE14"
unique(scRNAE14$orig.ident)  ##变量就存在这里
  ##得到了EOPE14的Seurat对象
saveRDS(scRNAE14,"./Seurat_data/EOPE14.rds")  ##保存rds文件
}

#这里是EOPE15
{
countsE15 <- read.table("./raw_data/large data/EOPE15/EOPE15.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsE15)[1] <- "Gene"
countsE15 <- countsE15 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsE15) <- countsE15[,1]
  ##把gene设置为列名
countsE15 <- countsE15[,-1]
  ##删除数据集中gene这一行
countsE15 <- countsE15[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsE15,"./raw_data/large data（处理）/EOPE15_num.csv")
countsE15_num <- read.csv("./raw_data/large data（处理）/EOPE15_num.csv")
rownames(countsE15_num) <- countsE15_num[,1]
countsE15_num <- countsE15_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsE15_num <- as(as.matrix(countsE15_num), "dgCMatrix")
scRNAE15 <- CreateSeuratObject(counts = countsE15_num)
unique(scRNAE15$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAE15$orig.ident <- "EOPE15"
unique(scRNAE15$orig.ident)  ##变量就存在这里
  ##得到了EOPE15的Seurat对象
saveRDS(scRNAE15,"./Seurat_data/EOPE15.rds")  ##保存rds文件
}

#这里是EOPE16
{
countsE16 <- read.table("./raw_data/large data/EOPE16/EOPE16.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsE16)[1] <- "Gene"
countsE16 <- countsE16 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsE16) <- countsE16[,1]
  ##把gene设置为列名
countsE16 <- countsE16[,-1]
  ##删除数据集中gene这一行
countsE16 <- countsE16[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsE16,"./raw_data/large data（处理）/EOPE16_num.csv")
countsE16_num <- read.csv("./raw_data/large data（处理）/EOPE16_num.csv")
rownames(countsE16_num) <- countsE16_num[,1]
countsE16_num <- countsE16_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsE16_num <- as(as.matrix(countsE16_num), "dgCMatrix")
scRNAE16 <- CreateSeuratObject(counts = countsE16_num)
unique(scRNAE16$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAE16$orig.ident <- "EOPE16"
unique(scRNAE16$orig.ident)  ##变量就存在这里
  ##得到了EOPE16的Seurat对象
saveRDS(scRNAE16,"./Seurat_data/EOPE16.rds")  ##保存rds文件
}

#这里是EOPE17
{
countsE17 <- read.table("./raw_data/large data/EOPE17/EOPE17.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsE17)[1] <- "Gene"
countsE17 <- countsE17 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsE17) <- countsE17[,1]
  ##把gene设置为列名
countsE17 <- countsE17[,-1]
  ##删除数据集中gene这一行
countsE17 <- countsE17[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsE17,"./raw_data/large data（处理）/EOPE17_num.csv")
countsE17_num <- read.csv("./raw_data/large data（处理）/EOPE17_num.csv")
rownames(countsE17_num) <- countsE17_num[,1]
countsE17_num <- countsE17_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsE17_num <- as(as.matrix(countsE17_num), "dgCMatrix")
scRNAE17 <- CreateSeuratObject(counts = countsE17_num)
unique(scRNAE17$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAE17$orig.ident <- "EOPE17"
unique(scRNAE17$orig.ident)  ##变量就存在这里
  ##得到了EOPE17的Seurat对象
saveRDS(scRNAE17,"./Seurat_data/EOPE17.rds")  ##保存rds文件
}

#这里是EOPE18
{
countsE18 <- read.table("./raw_data/large data/EOPE18/EOPE18.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsE18)[1] <- "Gene"
countsE18 <- countsE18 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsE18) <- countsE18[,1]
  ##把gene设置为列名
countsE18 <- countsE18[,-1]
  ##删除数据集中gene这一行
countsE18 <- countsE18[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsE18,"./raw_data/large data（处理）/EOPE18_num.csv")
countsE18_num <- read.csv("./raw_data/large data（处理）/EOPE18_num.csv")
rownames(countsE18_num) <- countsE18_num[,1]
countsE18_num <- countsE18_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsE18_num <- as(as.matrix(countsE18_num), "dgCMatrix")
scRNAE18 <- CreateSeuratObject(counts = countsE18_num)
unique(scRNAE18$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAE18$orig.ident <- "EOPE18"
unique(scRNAE18$orig.ident)  ##变量就存在这里
  ##得到了EOPE18的Seurat对象
saveRDS(scRNAE18,"./Seurat_data/EOPE18.rds")  ##保存rds文件
}

#这里是EOPE19
{
countsE19 <- read.table("./raw_data/large data/EOPE19/EOPE19.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsE19)[1] <- "Gene"
countsE19 <- countsE19 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsE19) <- countsE19[,1]
  ##把gene设置为列名
countsE19 <- countsE19[,-1]
  ##删除数据集中gene这一行
countsE19 <- countsE19[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsE19,"./raw_data/large data（处理）/EOPE19_num.csv")
countsE19_num <- read.csv("./raw_data/large data（处理）/EOPE19_num.csv")
rownames(countsE19_num) <- countsE19_num[,1]
countsE19_num <- countsE19_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsE19_num <- as(as.matrix(countsE19_num), "dgCMatrix")
scRNAE19 <- CreateSeuratObject(counts = countsE19_num)
unique(scRNAE19$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAE19$orig.ident <- "EOPE19"
unique(scRNAE19$orig.ident)  ##变量就存在这里
  ##得到了EOPE19的Seurat对象
saveRDS(scRNAE19,"./Seurat_data/EOPE19.rds")  ##保存rds文件
}

#这里是EOPE20
{
countsE20 <- read.table("./raw_data/large data/EOPE20/EOPE20.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsE20)[1] <- "Gene"
countsE20 <- countsE20 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsE20) <- countsE20[,1]
  ##把gene设置为列名
countsE20 <- countsE20[,-1]
  ##删除数据集中gene这一行
countsE20 <- countsE20[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsE20,"./raw_data/large data（处理）/EOPE20_num.csv")
countsE20_num <- read.csv("./raw_data/large data（处理）/EOPE20_num.csv")
rownames(countsE20_num) <- countsE20_num[,1]
countsE20_num <- countsE20_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsE20_num <- as(as.matrix(countsE20_num), "dgCMatrix")
scRNAE20 <- CreateSeuratObject(counts = countsE20_num)
unique(scRNAE20$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAE20$orig.ident <- "EOPE20"
unique(scRNAE20$orig.ident)  ##变量就存在这里
  ##得到了EOPE20的Seurat对象
saveRDS(scRNAE20,"./Seurat_data/EOPE20.rds")  ##保存rds文件
}

#这里是LOPE5
{
countsL5 <- read.table("./raw_data/large data/LOPE5/LOPE5.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsL5)[1] <- "Gene"
countsL5 <- countsL5 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsL5) <- countsL5[,1]
  ##把gene设置为列名
countsL5 <- countsL5[,-1]
  ##删除数据集中gene这一行
countsL5 <- countsL5[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsL5,"./raw_data/large data（处理）/LOPE5_num.csv")
countsL5_num <- read.csv("./raw_data/large data（处理）/LOPE5_num.csv")
rownames(countsL5_num) <- countsL5_num[,1]
countsL5_num <- countsL5_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsL5_num <- as(as.matrix(countsL5_num), "dgCMatrix")
scRNAL5 <- CreateSeuratObject(counts = countsL5_num)
unique(scRNAL5$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAL5$orig.ident <- "LOPE5"
unique(scRNAL5$orig.ident)  ##变量就存在这里
  ##得到了LOPE5的Seurat对象
saveRDS(scRNAL5,"./Seurat_data/LOPE5.rds")  ##保存rds文件
}

#这里是LOPE6
{
countsL6 <- read.table("./raw_data/large data/LOPE6/LOPE6.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsL6)[1] <- "Gene"
countsL6 <- countsL6 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsL6) <- countsL6[,1]
  ##把gene设置为列名
countsL6 <- countsL6[,-1]
  ##删除数据集中gene这一行
countsL6 <- countsL6[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsL6,"./raw_data/large data（处理）/LOPE6_num.csv")
countsL6_num <- read.csv("./raw_data/large data（处理）/LOPE6_num.csv")
rownames(countsL6_num) <- countsL6_num[,1]
countsL6_num <- countsL6_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsL6_num <- as(as.matrix(countsL6_num), "dgCMatrix")
scRNAL6 <- CreateSeuratObject(counts = countsL6_num)
unique(scRNAL6$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAL6$orig.ident <- "LOPE6"
unique(scRNAL6$orig.ident)  ##变量就存在这里
  ##得到了LOPE6的Seurat对象
saveRDS(scRNAL6,"./Seurat_data/LOPE6.rds")  ##保存rds文件
}

#这里是LOPE7
{
countsL7 <- read.table("./raw_data/large data/LOPE7/LOPE7.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsL7)[1] <- "Gene"
countsL7 <- countsL7 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsL7) <- countsL7[,1]
  ##把gene设置为列名
countsL7 <- countsL7[,-1]
  ##删除数据集中gene这一行
countsL7 <- countsL7[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsL7,"./raw_data/large data（处理）/LOPE7_num.csv")
countsL7_num <- read.csv("./raw_data/large data（处理）/LOPE7_num.csv")
rownames(countsL7_num) <- countsL7_num[,1]
countsL7_num <- countsL7_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsL7_num <- as(as.matrix(countsL7_num), "dgCMatrix")
scRNAL7 <- CreateSeuratObject(counts = countsL7_num)
unique(scRNAL7$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAL7$orig.ident <- "LOPE7"
unique(scRNAL7$orig.ident)  ##变量就存在这里
  ##得到了LOPE7的Seurat对象
saveRDS(scRNAL7,"./Seurat_data/LOPE7.rds")  ##保存rds文件
}

#这里是LOPE8
{
countsL8 <- read.table("./raw_data/large data/LOPE8/LOPE8.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsL8)[1] <- "Gene"
countsL8 <- countsL8 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsL8) <- countsL8[,1]
  ##把gene设置为列名
countsL8 <- countsL8[,-1]
  ##删除数据集中gene这一行
countsL8 <- countsL8[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsL8,"./raw_data/large data（处理）/LOPE8_num.csv")
countsL8_num <- read.csv("./raw_data/large data（处理）/LOPE8_num.csv")
rownames(countsL8_num) <- countsL8_num[,1]
countsL8_num <- countsL8_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsL8_num <- as(as.matrix(countsL8_num), "dgCMatrix")
scRNAL8 <- CreateSeuratObject(counts = countsL8_num)
unique(scRNAL8$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAL8$orig.ident <- "LOPE8"
unique(scRNAL8$orig.ident)  ##变量就存在这里
  ##得到了LOPE8的Seurat对象
saveRDS(scRNAL8,"./Seurat_data/LOPE8.rds")  ##保存rds文件
}

#这里是LOPE9
{
countsL9 <- read.table("./raw_data/large data/LOPE9/LOPE9.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsL9)[1] <- "Gene"
countsL9 <- countsL9 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsL9) <- countsL9[,1]
  ##把gene设置为列名
countsL9 <- countsL9[,-1]
  ##删除数据集中gene这一行
countsL9 <- countsL9[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsL9,"./raw_data/large data（处理）/LOPE9_num.csv")
countsL9_num <- read.csv("./raw_data/large data（处理）/LOPE9_num.csv")
rownames(countsL9_num) <- countsL9_num[,1]
countsL9_num <- countsL9_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsL9_num <- as(as.matrix(countsL9_num), "dgCMatrix")
scRNAL9 <- CreateSeuratObject(counts = countsL9_num)
unique(scRNAL9$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAL9$orig.ident <- "LOPE9"
unique(scRNAL9$orig.ident)  ##变量就存在这里
  ##得到了LOPE9的Seurat对象
saveRDS(scRNAL9,"./Seurat_data/LOPE9.rds")  ##保存rds文件
}

#这里是LOPE10
{
countsL10 <- read.table("./raw_data/large data/LOPE10/LOPE10.txt",
                       sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsL10)[1] <- "Gene"
countsL10 <- countsL10 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsL10) <- countsL10[,1]
  ##把gene设置为列名
countsL10 <- countsL10[,-1]
  ##删除数据集中gene这一行
countsL10 <- countsL10[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsL10,"./raw_data/large data（处理）/LOPE10_num.csv")
countsL10_num <- read.csv("./raw_data/large data（处理）/LOPE10_num.csv")
rownames(countsL10_num) <- countsL10_num[,1]
countsL10_num <- countsL10_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsL10_num <- as(as.matrix(countsL10_num), "dgCMatrix")
scRNAL10 <- CreateSeuratObject(counts = countsL10_num)
unique(scRNAL10$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAL10$orig.ident <- "LOPE10"
unique(scRNAL10$orig.ident)  ##变量就存在这里
  ##得到了LOPE10的Seurat对象
saveRDS(scRNAL10,"./Seurat_data/LOPE10.rds")  ##保存rds文件
}

#这里是LOPE11
{
countsL11 <- read.table("./raw_data/large data/LOPE11/LOPE11.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsL11)[1] <- "Gene"
countsL11 <- countsL11 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsL11) <- countsL11[,1]
  ##把gene设置为列名
countsL11 <- countsL11[,-1]
  ##删除数据集中gene这一行
countsL11 <- countsL11[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsL11,"./raw_data/large data（处理）/LOPE11_num.csv")
countsL11_num <- read.csv("./raw_data/large data（处理）/LOPE11_num.csv")
rownames(countsL11_num) <- countsL11_num[,1]
countsL11_num <- countsL11_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsL11_num <- as(as.matrix(countsL11_num), "dgCMatrix")
scRNAL11 <- CreateSeuratObject(counts = countsL11_num)
unique(scRNAL11$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAL11$orig.ident <- "LOPE11"
unique(scRNAL11$orig.ident)  ##变量就存在这里
  ##得到了LOPE11的Seurat对象
saveRDS(scRNAL11,"./Seurat_data/LOPE11.rds")  ##保存rds文件
}

#这里是LOPE12
{
countsL12 <- read.table("./raw_data/large data/LOPE12/LOPE12.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsL12)[1] <- "Gene"
countsL12 <- countsL12 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsL12) <- countsL12[,1]
  ##把gene设置为列名
countsL12 <- countsL12[,-1]
  ##删除数据集中gene这一行
countsL12 <- countsL12[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsL12,"./raw_data/large data（处理）/LOPE12_num.csv")
countsL12_num <- read.csv("./raw_data/large data（处理）/LOPE12_num.csv")
rownames(countsL12_num) <- countsL12_num[,1]
countsL12_num <- countsL12_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsL12_num <- as(as.matrix(countsL12_num), "dgCMatrix")
scRNAL12 <- CreateSeuratObject(counts = countsL12_num)
unique(scRNAL12$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAL12$orig.ident <- "LOPE12"
unique(scRNAL12$orig.ident)  ##变量就存在这里
  ##得到了LOPE12的Seurat对象
saveRDS(scRNAL12,"./Seurat_data/LOPE12.rds")  ##保存rds文件
}

#这里是LOPE13
{
countsL13 <- read.table("./raw_data/large data/LOPE13/LOPE13.txt",
                        sep = "\t",header = T)
  ##这里发现有重复的基因
  ##从上一行代码里删除“,row.names = 1”，并进行以下操作
colnames(countsL13)[1] <- "Gene"
countsL13 <- countsL13 %>% distinct(Gene, .keep_all = T)
  ##去重操作
rownames(countsL13) <- countsL13[,1]
  ##把gene设置为列名
countsL13 <- countsL13[,-1]
  ##删除数据集中gene这一行
countsL13 <- countsL13[-c(1:22),]
  ##删除表中多余的数据
  ##这里要保证，除了列名、行名外，数据集中全部是数字
  ##发现数据类型是"character"
  ##这里不能直接创建表达矩阵，进行下面处理
write.csv(countsL13,"./raw_data/large data（处理）/LOPE13_num.csv")
countsL13_num <- read.csv("./raw_data/large data（处理）/LOPE13_num.csv")
rownames(countsL13_num) <- countsL13_num[,1]
countsL13_num <- countsL13_num[,-1]
  ##处理实质就是删除掉造成character的行后重新载入
countsL13_num <- as(as.matrix(countsL13_num), "dgCMatrix")
scRNAL13 <- CreateSeuratObject(counts = countsL13_num)
unique(scRNAL13$orig.ident)  
  ##这里发现运行结果和实验目的不符，下面手动修改
scRNAL13$orig.ident <- "LOPE13"
unique(scRNAL13$orig.ident)  ##变量就存在这里
  ##得到了LOPE13的Seurat对象
saveRDS(scRNAL13,"./Seurat_data/LOPE13.rds")  ##保存rds文件
}
}

#保存RData
save.image("./Rdata/1.数据载入并建立Seur对象.RData")
