rm(list = ls())
setwd("D:/2023.10/PE scRNA（多数据库运行）")
getwd()

#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(multtest))BiocManager::install("multtest")
  if(!require(dplyr))install.packages("dplyr")
  if(!require(harmony))install.packages("harmony")
  if(!require(data.table))install.packages("data.table")
}

#定义一个函数帮助NormalizeData（标准化表达矩阵）与FindVariableFeatures（高变基因计算）
Normol <- function(testA.seu){
  testA.seu <- NormalizeData(testA.seu, normalization.method = "LogNormalize", scale.factor = 10000)
  testA.seu <- FindVariableFeatures(testA.seu, selection.method = "vst", nfeatures = 2000)
  return(testA.seu)
}

#首先对control进行质量控制
{
control <- readRDS("./Seurat_data（整理）/control.rds")
control@meta.data$group <- "control"
control <- Normol(control)

##查看质量控制前的各细胞RNA表达数、RNA表达的种类数、线粒体RNA所占的比例
control[["percent.mt"]] <- PercentageFeatureSet(control,pattern = "^MT-")
control$orig.ident <- factor(x = control$orig.ident,
                             levels = c("control1","control2","control3","control4","control5",
                                        "control6","control7","control8","control9","control10",
                                        "control11","control12","control13","control14","control15",
                                        "control16","control17","control18","control19","control20",
                                        "control21","control22","control23","control24"))
  ##做图发现小提琴图的x轴顺序不满意
  ##使用上面的步骤修改x轴的顺序
preQC_control <- VlnPlot(control, 
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                         raster=FALSE,
                         group.by = "orig.ident",
                         pt.size = 0)
preQC_control  ##得到1.control的原始质量图 -保存图片

##查看原始的细胞数目
control_cellnumber <- as.data.table(table(control$orig.ident))
write.csv(control_cellnumber,"./result/control质控前细胞数.csv")

##质量控制
control <- subset(control, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)

##查看质量控制后的各细胞RNA表达数、RNA表达的种类数、线粒体RNA所占的比例
afterQC_control <- VlnPlot(control, 
                           features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                           raster=FALSE,
                           group.by = "orig.ident",
                           pt.size = 0)
afterQC_control  ##得到2.control的质控后质量图 -保存图片

##查看质控后的细胞数目
control_cellnumber <- as.data.table(table(control$orig.ident))
write.csv(control_cellnumber,"./result/control质控后细胞数.csv")

saveRDS(control,"./Seurat_data（整理）/control（质控后）.rds")
}

#下面对EOPE进行质量控制
{
EOPE <- readRDS("./Seurat_data（整理）/EOPE.rds")
EOPE@meta.data$group <- "EOPE"
EOPE <- Normol(EOPE)

##查看质量控制前的各细胞RNA表达数、RNA表达的种类数、线粒体RNA所占的比例
EOPE[["percent.mt"]] <- PercentageFeatureSet(EOPE,pattern = "^MT-")
EOPE$orig.ident <- factor(x = EOPE$orig.ident,
                             levels = c("EOPE1","EOPE2","EOPE3","EOPE4","EOPE5",
                                        "EOPE6","EOPE7","EOPE8","EOPE9","EOPE10",
                                        "EOPE11","EOPE12","EOPE13","EOPE14","EOPE15",
                                        "EOPE16","EOPE17","EOPE18","EOPE19","EOPE20"))
  ##做图发现小提琴图的x轴顺序不满意
  ##使用上面的步骤修改x轴的顺序
preQC_EOPE <- VlnPlot(EOPE, 
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     raster=FALSE,
                     group.by = "orig.ident",
                     pt.size = 0)
preQC_EOPE  ##得到3.EOPE的原始质量图 -保存图片

##查看原始的细胞数目
EOPE_cellnumber <- as.data.table(table(EOPE$orig.ident))
write.csv(EOPE_cellnumber,"./result/EOPE质控前细胞数.csv")

##质量控制
EOPE <- subset(EOPE, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)

##查看质量控制后的各细胞RNA表达数、RNA表达的种类数、线粒体RNA所占的比例
afterQC_EOPE <- VlnPlot(EOPE, 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        raster=FALSE,
                        group.by = "orig.ident",
                        pt.size = 0)
afterQC_EOPE  ##得到4.EOPE的质控后质量图 -保存图片

##查看质控后的细胞数目
EOPE_cellnumber <- as.data.table(table(EOPE$orig.ident))
write.csv(EOPE_cellnumber,"./result/EOPE质控后细胞数.csv")

saveRDS(EOPE,"./Seurat_data（整理）/EOPE（质控后）.rds")
}

#下面对LOPE进行质量控制
{
LOPE <- readRDS("./Seurat_data（整理）/LOPE.rds")
LOPE@meta.data$group <- "LOPE"
LOPE <- Normol(LOPE)

##查看质量控制前的各细胞RNA表达数、RNA表达的种类数、线粒体RNA所占的比例
LOPE[["percent.mt"]] <- PercentageFeatureSet(LOPE,pattern = "^MT-")
LOPE$orig.ident <- factor(x = LOPE$orig.ident,
                          levels = c("LOPE1","LOPE2","LOPE3","LOPE4","LOPE5",
                                     "LOPE6","LOPE7","LOPE8","LOPE9","LOPE10",
                                     "LOPE11","LOPE12","LOPE13"))
  ##做图发现小提琴图的x轴顺序不满意
  ##使用上面的步骤修改x轴的顺序
preQC_LOPE <- VlnPlot(LOPE, 
                      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                      raster=FALSE,
                      group.by = "orig.ident",
                      pt.size = 0)
preQC_LOPE  ##得到5.LOPE的原始质量图 -保存图片

##查看原始的细胞数目
LOPE_cellnumber <- as.data.table(table(LOPE$orig.ident))
write.csv(LOPE_cellnumber,"./result/LOPE质控前细胞数.csv")

##质量控制
LOPE <- subset(LOPE, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)

##查看质量控制后的各细胞RNA表达数、RNA表达的种类数、线粒体RNA所占的比例
afterQC_LOPE <- VlnPlot(LOPE, 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        raster=FALSE,
                        group.by = "orig.ident",
                        pt.size = 0)
afterQC_LOPE  ##得到6.LOPE的质控后质量图 -保存图片

##查看质控后的细胞数目
LOPE_cellnumber <- as.data.table(table(LOPE$orig.ident))
write.csv(LOPE_cellnumber,"./result/LOPE质控后细胞数.csv")

saveRDS(LOPE,"./Seurat_data（整理）/LOPE（质控后）.rds")
}

#下面合并3个对象
control <- readRDS("./Seurat_data（整理）/control（质控后）.rds")
LOPE <- readRDS("./Seurat_data（整理）/LOPE（质控后）.rds")
EOPE <- readRDS("./Seurat_data（整理）/EOPE（质控后）.rds")

##首先使用merge函数简单合并数据
scRNA_raw <- merge(LOPE, 
                   y = c(EOPE,control), 
                   add.cell.ids = c("control","EOPE","LOPE"), 
                   project = "ALL")
saveRDS(scRNA_raw,"./Seurat_data（整理）/合并（raw）.rds")

##精简工作环境
rm(control,EOPE,LOPE)
gc()

##使用harmony去除批次效应
scRNA_raw <- NormalizeData(scRNA_raw,
                           normalization.method = "LogNormalize",scale.factor = 10000) 
  #NormalizeData进行标准化
scRNA_raw <- FindVariableFeatures(scRNA_raw,
                                  selection.method = "vst", nfeatures = 2000)
  #FindVariableFeatures计算高变基因
scRNA_raw <- ScaleData(scRNA_raw, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo",
                                                      "percent.mt", "nCount_RNA"),
                       verbose = T)
  #ScaleData选择矫正的因素
scRNA_raw <- RunPCA(scRNA_raw, features = VariableFeatures(object = scRNA_raw))
scRNA <- RunHarmony(scRNA_raw, "group",
                    verbose = TRUE)
names(scRNA@reductions)
scRNA <- RunUMAP(scRNA, dims = 1:20, reduction = "harmony",
                 verbose = TRUE)
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:20,
                       verbose = TRUE) 

saveRDS(scRNA,"./Seurat_data（整理）/合并.rds")
