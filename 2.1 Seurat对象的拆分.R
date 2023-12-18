rm(list = ls())
setwd("D:/2023.10/PE scRNA（多数据库运行）")
getwd()

#打开必要的package
{
if(!require(Seurat))install.packages("Seurat")
if(!require(multtest))BiocManager::install("multtest")
if(!require(dplyr))install.packages("dplyr")
}

#对象拆分的答疑
{
#注1：
#GSE173193和GSE182381分别使用了两种不同函数
#即SplitObject函数和subset函数
#两种方法都可以从大的Seurat对象中按条件筛选，从而得到需要的子对象
  #SplitObject函数使用"split.by = "限定筛选的对象
  #sbuset函数使用"=="左边是筛选的对象
#两种方法没有差异
}

#载入并拆分GSE173193
{
Seurat1 <- readRDS("./Seurat_data/con1-2和EOPE1-2.rds")
  ##这是一个多个对象构成的数据
  ##下面对它进行拆分和整合
unique(Seurat1@meta.data$orig.ident)
Seurat1.list <- SplitObject(Seurat1,split.by = "orig.ident") 
control1 <- Seurat1.list$control1
control1@meta.data$orig.ident <- "control1"
saveRDS(control1,"./Seurat_data（整理）/control/con1.rds")
  ##得到control1的Seurat对象
control2 <- Seurat1.list$control2
control2@meta.data$orig.ident <- "control2"
saveRDS(control2,"./Seurat_data（整理）/control/con2.rds")
  ##得到control2的Seurat对象
EOPE1 <- Seurat1.list$EOPE1
EOPE1@meta.data$orig.ident <- "EOPE1"
saveRDS(EOPE1,"./Seurat_data（整理）/EOPE/EOPE1.rds")
  ##得到EOPE1的Seurat对象
EOPE2 <- Seurat1.list$EOPE2
EOPE2@meta.data$orig.ident <- "EOPE2"
saveRDS(EOPE2,"./Seurat_data（整理）/EOPE/EOPE2.rds")
  ##得到EOPE2的Seurat对象
}

#载入并拆分GSE182381
{
Seurat2 <- readRDS("./Seurat_data/con5-8.rds")
  ##这是一个多个对象构成的数据
  ##下面对它进行拆分和整合
unique(Seurat2@meta.data$orig.ident)
control5 <- subset(Seurat2,orig.ident == "control5")
control5@meta.data$orig.ident <- "control5"
saveRDS(control5,"./Seurat_data（整理）/control/con5.rds")
  ##得到control5的Seurat对象
control6 <- subset(Seurat2,orig.ident == "control6")
control6@meta.data$orig.ident <- "control6"
saveRDS(control6,"./Seurat_data（整理）/control/con6.rds")
  ##得到control6的Seurat对象
control7 <- subset(Seurat2,orig.ident == "control7")
control7@meta.data$orig.ident <- "control7"
saveRDS(control7,"./Seurat_data（整理）/control/con7.rds")
  ##得到control7的Seurat对象
control8 <- subset(Seurat2,orig.ident == "control8")
control8@meta.data$orig.ident <- "control8"
saveRDS(control8,"./Seurat_data（整理）/control/con8.rds")
  ##得到control8的Seurat对象
}

#载入并拆分GSE192693
{
Seurat3 <- readRDS("./Seurat_data/con12-15和EOPE3-8.rds")
  ##这是一个多个对象构成的数据
  ##下面对它进行拆分和整合
unique(Seurat3@meta.data$orig.ident)
control12 <- subset(Seurat3,orig.ident == "control12")
control12@meta.data$orig.ident <- "control12"
saveRDS(control12,"./Seurat_data（整理）/control/con12.rds")
  ##得到control12的Seurat对象
control13 <- subset(Seurat3,orig.ident == "control13")
control13@meta.data$orig.ident <- "control13"
saveRDS(control13,"./Seurat_data（整理）/control/con13.rds")
  ##得到control13的Seurat对象
control14 <- subset(Seurat3,orig.ident == "control14")
control14@meta.data$orig.ident <- "control14"
saveRDS(control14,"./Seurat_data（整理）/control/con14.rds")
  ##得到control14的Seurat对象
control15 <- subset(Seurat3,orig.ident == "control15")
control15@meta.data$orig.ident <- "control15"
saveRDS(control15,"./Seurat_data（整理）/control/con15.rds")
  ##得到control15的Seurat对象
EOPE3 <- subset(Seurat3,orig.ident == "EOPE3")
EOPE3@meta.data$orig.ident <- "EOPE3"
saveRDS(EOPE3,"./Seurat_data（整理）/EOPE/EOPE3.rds")
  ##得到EOPE3的Seurat对象
EOPE4 <- subset(Seurat3,orig.ident == "EOPE4")
EOPE4@meta.data$orig.ident <- "EOPE4"
saveRDS(EOPE4,"./Seurat_data（整理）/EOPE/EOPE4.rds")
  ##得到EOPE4的Seurat对象
EOPE5 <- subset(Seurat3,orig.ident == "EOPE5")
EOPE5@meta.data$orig.ident <- "EOPE5"
saveRDS(EOPE5,"./Seurat_data（整理）/EOPE/EOPE5.rds")
  ##得到EOPE5的Seurat对象
EOPE6 <- subset(Seurat3,orig.ident == "EOPE6")
EOPE6@meta.data$orig.ident <- "EOPE6"
saveRDS(EOPE6,"./Seurat_data（整理）/EOPE/EOPE6.rds")
  ##得到EOPE6的Seurat对象
EOPE7 <- subset(Seurat3,orig.ident == "EOPE7")
EOPE7@meta.data$orig.ident <- "EOPE7"
saveRDS(EOPE7,"./Seurat_data（整理）/EOPE/EOPE7.rds")
  ##得到EOPE7的Seurat对象
EOPE8 <- subset(Seurat3,orig.ident == "EOPE8")
EOPE8@meta.data$orig.ident <- "EOPE8"
saveRDS(EOPE8,"./Seurat_data（整理）/EOPE/EOPE8.rds")
  ##得到EOPE8的Seurat对象
}

save.image("./Rdata/2.1 Seurat对象的拆分.RData")
