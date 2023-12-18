rm(list = ls())
setwd("E:/2023.10/PE scRNA（多数据库运行）")
getwd()

#打开必要的package
{
if(!require(Seurat))install.packages("Seurat")
if(!require(multtest))BiocManager::install("multtest")
if(!require(dplyr))install.packages("dplyr")
if(!require(harmony))install.packages("harmony")
}

#创建对象的答疑
{
#注1：
#这里的方法适用于有多个样本时
#只有1个样本时不需要进行ancher操作
}

#定义一个函数帮助NormalizeData（标准化表达矩阵）与FindVariableFeatures（高变基因计算）
Normol <- function(testA.seu){
  testA.seu <- NormalizeData(testA.seu, normalization.method = "LogNormalize", scale.factor = 10000)
  testA.seu <- FindVariableFeatures(testA.seu, selection.method = "vst", nfeatures = 2000)
  return(testA.seu)
}

#首先处理control
{
#载入control的所有Seurat对象
{
control1 <- readRDS(file = "./Seurat_data（整理）/control/con1.rds")
control2 <- readRDS(file = "./Seurat_data（整理）/control/con2.rds")
control3 <- readRDS(file = "./Seurat_data（整理）/control/con3.rds")
control4 <- readRDS(file = "./Seurat_data（整理）/control/con4.rds")
control5 <- readRDS(file = "./Seurat_data（整理）/control/con5.rds")
control6 <- readRDS(file = "./Seurat_data（整理）/control/con6.rds")
control7 <- readRDS(file = "./Seurat_data（整理）/control/con7.rds")
control8 <- readRDS(file = "./Seurat_data（整理）/control/con8.rds")
control9 <- readRDS(file = "./Seurat_data（整理）/control/con9.rds")
control10 <- readRDS(file = "./Seurat_data（整理）/control/con10.rds")
control11 <- readRDS(file = "./Seurat_data（整理）/control/con11.rds")
control12 <- readRDS(file = "./Seurat_data（整理）/control/con12.rds")
control13 <- readRDS(file = "./Seurat_data（整理）/control/con13.rds")
control14 <- readRDS(file = "./Seurat_data（整理）/control/con14.rds")
control15 <- readRDS(file = "./Seurat_data（整理）/control/con15.rds")
control16 <- readRDS(file = "./Seurat_data（整理）/control/con16.rds")
control17 <- readRDS(file = "./Seurat_data（整理）/control/con17.rds")
control18 <- readRDS(file = "./Seurat_data（整理）/control/con18.rds")
control19 <- readRDS(file = "./Seurat_data（整理）/control/con19.rds")
control20 <- readRDS(file = "./Seurat_data（整理）/control/con20.rds")
control21 <- readRDS(file = "./Seurat_data（整理）/control/con21.rds")
control22 <- readRDS(file = "./Seurat_data（整理）/control/con22.rds")
control23 <- readRDS(file = "./Seurat_data（整理）/control/con23.rds")
control24 <- readRDS(file = "./Seurat_data（整理）/control/con24.rds")
}

#将control的Seurat对象标准化并计算高变基因
{
control1 <- Normol(control1)
control2 <- Normol(control2)
control3 <- Normol(control3)
control4 <- Normol(control4)
control5 <- Normol(control5)
control6 <- Normol(control6)
control7 <- Normol(control7)
control8 <- Normol(control8)
control9 <- Normol(control9)
control10 <- Normol(control10)
control11 <- Normol(control11)
control12 <- Normol(control12)
control13 <- Normol(control13)
control14 <- Normol(control14)
control15 <- Normol(control15)
control16 <- Normol(control16)
control17 <- Normol(control17)
control18 <- Normol(control18)
control19 <- Normol(control19)
control20 <- Normol(control20)
control21 <- Normol(control21)
control22 <- Normol(control22)
control23 <- Normol(control23)
control24 <- Normol(control24)
}

#得到control的合并Seurat对象
{
#首先使用merge函数简单合并数据
  ##mergr只能简单的合并对象，不具有去除批次化效应的作用
{
control_raw1 <- merge(control1, 
                      y = c(control2,control3,
                            control4,control5,
                            control6,control7,
                            control8,control9), 
                      add.cell.ids = c("control1", 
                                       "control2","control3",
                                       "control4","control5",
                                       "control6","control7",
                                       "control8","control9"), 
                      project = "ALL")
control_raw2 <- merge(control10, 
                      y = c(control11,control12,
                            control13,control14,
                            control15,control16,
                            control17,control18), 
                      add.cell.ids = c("control10", 
                                       "control11","control12",
                                       "control13","control14",
                                       "control15","control16",
                                       "control17","control18"), 
                      project = "ALL")
control_raw3 <- merge(control19, 
                      y = c(control20,control21,
                            control22,control23,
                            control24), 
                      add.cell.ids = c("control19", 
                                       "control20","control21",
                                       "control22","control23",
                                       "control24"), 
                      project = "ALL")
control_raw <- merge(control_raw1, 
                      y = c(control_raw2,control_raw3), 
                      add.cell.ids = c("control_raw1",
                                       "control_raw2",
                                       "control_raw3"), 
                      project = "ALL")
}

#精简工作环境
{
rm(control1,control2,control3,control4,control5,
   control6,control7,control8,control9,control10,
   control11,control12,control13,control14,control15,
   control16,control17,control18,control19,control20,
   control21,control22,control23,control24)
gc()
}

#使用harmony去除批次效应
{
control_raw <- NormalizeData(control_raw,
                             normalization.method = "LogNormalize",scale.factor = 10000) 
  #NormalizeData进行标准化
control_raw <- FindVariableFeatures(control_raw,
                                    selection.method = "vst", nfeatures = 2000)
  #FindVariableFeatures计算高变基因
control_raw <- ScaleData(control_raw, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo",
                                                          "percent.mt", "nCount_RNA"),
                         verbose = T)
  #ScaleData选择矫正的因素
control_raw <- RunPCA(control_raw, features = VariableFeatures(object = control_raw))
control <- RunHarmony(control_raw, "orig.ident")
names(control@reductions)
control <- RunUMAP(control,  dims = 1:20, reduction = "harmony")
unique(control@meta.data$orig.ident)
control <- FindNeighbors(control, reduction = "harmony",dims = 1:20) 
} 

saveRDS(control,"./Seurat_data（整理）/control.rds")
}
}

#下面处理EOPE
{
#载入EOPE的所有Seurat对象
{
EOPE1 <- readRDS("./Seurat_data（整理）/EOPE/EOPE1.rds")
EOPE2 <- readRDS("./Seurat_data（整理）/EOPE/EOPE2.rds")
EOPE3 <- readRDS("./Seurat_data（整理）/EOPE/EOPE3.rds")
EOPE4 <- readRDS("./Seurat_data（整理）/EOPE/EOPE4.rds")
EOPE5 <- readRDS("./Seurat_data（整理）/EOPE/EOPE5.rds")
EOPE6 <- readRDS("./Seurat_data（整理）/EOPE/EOPE6.rds")
EOPE7 <- readRDS("./Seurat_data（整理）/EOPE/EOPE7.rds")
EOPE8 <- readRDS("./Seurat_data（整理）/EOPE/EOPE8.rds")
EOPE9 <- readRDS("./Seurat_data（整理）/EOPE/EOPE9.rds")
EOPE10 <- readRDS("./Seurat_data（整理）/EOPE/EOPE10.rds")
EOPE11 <- readRDS("./Seurat_data（整理）/EOPE/EOPE11.rds")
EOPE12 <- readRDS("./Seurat_data（整理）/EOPE/EOPE12.rds")
EOPE13 <- readRDS("./Seurat_data（整理）/EOPE/EOPE13.rds")
EOPE14 <- readRDS("./Seurat_data（整理）/EOPE/EOPE14.rds")
EOPE15 <- readRDS("./Seurat_data（整理）/EOPE/EOPE15.rds")
EOPE16 <- readRDS("./Seurat_data（整理）/EOPE/EOPE16.rds")
EOPE17 <- readRDS("./Seurat_data（整理）/EOPE/EOPE17.rds")
EOPE18 <- readRDS("./Seurat_data（整理）/EOPE/EOPE18.rds")
EOPE19 <- readRDS("./Seurat_data（整理）/EOPE/EOPE19.rds")
EOPE20 <- readRDS("./Seurat_data（整理）/EOPE/EOPE20.rds")
}

#将EOPE的Seurat对象标准化并计算高变基因
{
EOPE1 <- Normol(EOPE1)
EOPE2 <- Normol(EOPE2)
EOPE3 <- Normol(EOPE3)
EOPE4 <- Normol(EOPE4)
EOPE5 <- Normol(EOPE5)
EOPE6 <- Normol(EOPE6)
EOPE7 <- Normol(EOPE7)
EOPE8 <- Normol(EOPE8)
EOPE9 <- Normol(EOPE9)
EOPE10 <- Normol(EOPE10)
EOPE11 <- Normol(EOPE11)
EOPE12 <- Normol(EOPE12)
EOPE13 <- Normol(EOPE13)
EOPE14 <- Normol(EOPE14)
EOPE15 <- Normol(EOPE15)
EOPE16 <- Normol(EOPE16)
EOPE17 <- Normol(EOPE17)
EOPE18 <- Normol(EOPE18)
EOPE19 <- Normol(EOPE19)
EOPE20 <- Normol(EOPE20)
}

#得到EOPE的合并Seurat对象
{
#首先使用merge函数简单合并数据
  ##mergr只能简单的合并对象，不具有去除批次化效应的作用
{
EOPE_raw1 <- merge(EOPE1, 
                      y = c(EOPE2,EOPE3,
                            EOPE4,EOPE5,
                            EOPE6,EOPE7,
                            EOPE8,EOPE9,
                            EOPE10), 
                      add.cell.ids = c("EOPE2","EOPE3",
                                       "EOPE4","EOPE5",
                                       "EOPE6","EOPE7",
                                       "EOPE8","EOPE9",
                                       "EOPE10","EOPE1"), 
                      project = "ALL")
EOPE_raw2 <- merge(EOPE11, 
                   y = c(EOPE12,EOPE13,
                         EOPE14,EOPE15,
                         EOPE16,EOPE17,
                         EOPE18,EOPE19,
                         EOPE20), 
                   add.cell.ids = c("EOPE12","EOPE13",
                                    "EOPE14","EOPE15",
                                    "EOPE16","EOPE17",
                                    "EOPE18","EOPE19",
                                    "EOPE20","EOPE11"), 
                   project = "ALL")
EOPE_raw <- merge(EOPE_raw1, 
                  y = c(EOPE_raw2), 
                  add.cell.ids = c("EOPE_raw1","EOPE_raw2"), 
                  project = "ALL")
}

#精简工作环境
{
rm(EOPE1,EOPE2,EOPE3,EOPE4,EOPE5,
   EOPE6,EOPE7,EOPE8,EOPE9,EOPE10,
   EOPE11,EOPE12,EOPE13,EOPE14,EOPE15,
   EOPE16,EOPE17,EOPE18,EOPE19,EOPE20)
gc()
}
  
#使用harmony去除批次效应
{
EOPE_raw <- NormalizeData(EOPE_raw,
                          normalization.method = "LogNormalize",scale.factor = 10000) 
  #NormalizeData进行标准化
EOPE_raw <- FindVariableFeatures(EOPE_raw,
                                 selection.method = "vst", nfeatures = 2000)
  #FindVariableFeatures计算高变基因
EOPE_raw <- ScaleData(EOPE_raw, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo",
                                                    "percent.mt", "nCount_RNA"),
                      verbose = T)
  #ScaleData选择矫正的因素
EOPE_raw <- RunPCA(EOPE_raw, features = VariableFeatures(object = EOPE_raw))
EOPE <- RunHarmony(EOPE_raw, "orig.ident")
names(EOPE@reductions)
EOPE <- RunUMAP(EOPE,  dims = 1:20, reduction = "harmony")
unique(EOPE@meta.data$orig.ident)
EOPE <- FindNeighbors(EOPE, reduction = "harmony",dims = 1:20) 
} 
saveRDS(EOPE,"./Seurat_data（整理）/EOPE.rds")
}
}

#下面处理LOPE
{
#载入LOPE的所有Seurat对象
{
LOPE1 <- readRDS("./Seurat_data（整理）/LOPE/LOPE1.rds")
LOPE2 <- readRDS("./Seurat_data（整理）/LOPE/LOPE2.rds")
LOPE3 <- readRDS("./Seurat_data（整理）/LOPE/LOPE3.rds")
LOPE4 <- readRDS("./Seurat_data（整理）/LOPE/LOPE4.rds")
LOPE5 <- readRDS("./Seurat_data（整理）/LOPE/LOPE5.rds")
LOPE6 <- readRDS("./Seurat_data（整理）/LOPE/LOPE6.rds")
LOPE7 <- readRDS("./Seurat_data（整理）/LOPE/LOPE7.rds")
LOPE8 <- readRDS("./Seurat_data（整理）/LOPE/LOPE8.rds")
LOPE9 <- readRDS("./Seurat_data（整理）/LOPE/LOPE9.rds")
LOPE10 <- readRDS("./Seurat_data（整理）/LOPE/LOPE10.rds")
LOPE11 <- readRDS("./Seurat_data（整理）/LOPE/LOPE11.rds")
LOPE12 <- readRDS("./Seurat_data（整理）/LOPE/LOPE12.rds")
LOPE13 <- readRDS("./Seurat_data（整理）/LOPE/LOPE13.rds")
}

#将LOPE的Seurat对象标准化并计算高变基因
{
LOPE1 <- Normol(LOPE1)
LOPE2 <- Normol(LOPE2)
LOPE3 <- Normol(LOPE3)
LOPE4 <- Normol(LOPE4)
LOPE5 <- Normol(LOPE5)
LOPE6 <- Normol(LOPE6)
LOPE7 <- Normol(LOPE7)
LOPE8 <- Normol(LOPE8)
LOPE9 <- Normol(LOPE9)
LOPE10 <- Normol(LOPE10)
LOPE11 <- Normol(LOPE11)
LOPE12 <- Normol(LOPE12)
LOPE13 <- Normol(LOPE13)
}

#得到LOPE的合并Seurat对象
{
#首先使用merge函数简单合并数据
  ##mergr只能简单的合并对象，不具有去除批次化效应的作用
{  
LOPE_raw <- merge(LOPE1, 
                  y = c(LOPE2,LOPE3,
                        LOPE4,LOPE5,
                        LOPE6,LOPE7,
                        LOPE8,LOPE9,
                        LOPE10,LOPE11,
                        LOPE12,LOPE13), 
                  add.cell.ids = c("LOPE2","LOPE3",
                                   "LOPE4","LOPE5",
                                   "LOPE6","LOPE7",
                                   "LOPE8","LOPE9",
                                   "LOPE10","LOPE11",
                                   "LOPE12","LOPE13",
                                   "LOPE1"), 
                  project = "ALL")
}
  
#精简工作环境
{  
rm(LOPE2,LOPE3,LOPE4,LOPE5,LOPE6,LOPE7,LOPE8,LOPE9,
   LOPE10,LOPE11,LOPE12,LOPE13,LOPE1)
gc()
}  
  
#使用harmony去除批次效应
{  
LOPE_raw <- NormalizeData(LOPE_raw,
                          normalization.method = "LogNormalize",scale.factor = 10000) 
  #NormalizeData进行标准化
LOPE_raw <- FindVariableFeatures(LOPE_raw,
                                 selection.method = "vst", nfeatures = 2000)
  #FindVariableFeatures计算高变基因
LOPE_raw <- ScaleData(LOPE_raw, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo",
                                                    "percent.mt", "nCount_RNA"),
                      verbose = T)
  #ScaleData选择矫正的因素
LOPE_raw <- RunPCA(LOPE_raw, features = VariableFeatures(object = LOPE_raw))
LOPE <- RunHarmony(LOPE_raw, "orig.ident")
names(LOPE@reductions)
LOPE <- RunUMAP(LOPE,  dims = 1:20, reduction = "harmony")
unique(LOPE@meta.data$orig.ident)
LOPE <- FindNeighbors(LOPE, reduction = "harmony",dims = 1:20) 
}  
saveRDS(LOPE,"./Seurat_data（整理）/LOPE.rds")
}
}
