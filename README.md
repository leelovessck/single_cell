# single_cell
#single-cell
这里面包括了单细胞原始数据、代码、结果及R环境
阅读时建议使用“code”界面，而不是“preview”
  "code"界面中排版更清晰


#raw_data文件夹：
这个文件夹里是初步处理过的单细胞数据，里面包括了有3种原始文件类型：
  （1）RDa
  （2）txt
  （3）cellranger（barcodes、features、matrix3个文件）
文件的名字就是该文件里面包含的样本

#原始数据的答疑
#注1：
  受github上传文件大小限制，部分原始数据无法上传，缺少的原始数据可以参考以下链接：
  

#注2：
  如果下面涉及到原始数据不清楚从哪里来的，或者不了解分组情况
  参考“下载的原始数据”文件夹中的“使用的数据-总结.xlsx”
    在这2个Excel文件中，可以得到以下内容 
      1.使用的数据集的下载链接或GEO ID
      2.与数据集相关的原文标题
      3.数据集构成的简介

#注3：
  单细胞的数据还有h5、h5ad等文件形式，这两种都是来自Python的对象
  #下面是h5文件的载入方法（使用hdf5r包）
    in[1]：devtools::install_github(repo = "hhoeflin/hdf5r")
    in[2]：library(hdf5r)
    in[3]：scRNA1 <-Read10X_h5("xxx.h5")
  #下面是h5ad文件的载入方法（使用SeuratDisk包）
    in[1]：remotes::install_github("mojaveazure/seurat-disk", force = TRUE)
    in[2]：library(SeuratDisk)
    in[3]：setwd('D:/')  ##路径要是英文
    in[4]：Convert("data.h5ad", "h5seurat",overwrite = TRUE, assay = "RNA")
    in[5]：scRNA1 <- LoadH5Seurat("data.h5seurat")

#注4：
  rds文件是二进制文件
    在新的项目中，载入rds文件可以得到相同的对象
      例：保存Seurat的rds，再次载入可以在新的R里得直接到Seurat对象
      例：保存CellChat的rds，再次载入可以在新的R里直接得到CellChat对象

#注5：
  在Seurat对象中，metadata可以理解成一个Excel表格
  #一种修改metadata的方式
    in[1]:Seurat_name@meta.data$orig.ident <- xxx
  #另外一种修改方式：
    in[1]:metadata <- Seurat_name@meta.data”  ##提取出metadata
    然后按矩阵修改，或者保存导出后在excel中修改
    in[2]:Seurat_name@meta.data <- metadata  ##把修改好的metadata导入
  
#注6：
  如果是表达矩阵形式的文件（txt或tsv），使用“read.table”函数读取
  如何同时读取一个目录下的多个表达矩阵，可能需要使用循环语句（技能待学习）
    下面给出的方法可能有点繁琐，但是能用
  
#注7：
  GSE173193给出的是同时读取同一目录下多个cellranger输出文件的方法
  #如果只有一个文件，可以简单通过以下两行代码完成
    in[1]：counts1 <- Read10X(data.dir = "D:/XX/文件目录/")  ##最后一个/必不可少
    in[2]：scRNA1 = CreateSeuratObject(counts1)
