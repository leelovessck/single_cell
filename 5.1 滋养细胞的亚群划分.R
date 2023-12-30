rm(list = ls())
setwd("E:/2023.10/PE scRNA（多数据库运行）")
getwd()

#打开必要的package
{
if(!require(Seurat))install.packages("Seurat")
if(!require(clustree))install.packages("clustree")
if(!require(scatterplot3d))install.packages("scatterplot3d")
if(!require(dplyr))install.packages("dplyr")
if(!require(ggplot2))install.packages("ggplot2")
if(!require(data.table))install.packages("data.table")
}

#载入seurat对象
object <- readRDS("./Seurat_data（整理）/合并（细胞注释）.rds")

#首先处理EVT
{
#筛选出EVT细胞
EVT <- subset(object,celltype == "EVT")
DimPlot(EVT,label = TRUE,raster=FALSE)|
  DimPlot(object,label = TRUE,raster=FALSE)

#EVT细胞的再次降维
EVT <- NormalizeData(EVT,
                     normalization.method = "LogNormalize",scale.factor = 10000) 
  ##NormalizeData进行标准化
EVT <- FindVariableFeatures(EVT,
                            selection.method = "vst", nfeatures = 2000)
  ##FindVariableFeatures计算高变基因
EVT <- ScaleData(EVT, features = rownames(EVT))
EVT <- RunPCA(EVT, features = VariableFeatures(object = EVT))
EVT <- RunUMAP(EVT, dims = 1:20, reduction = "harmony",
               verbose = TRUE)
EVT <- FindNeighbors(EVT, reduction = "harmony",dims = 1:20,
                     verbose = TRUE) 

#设置不同的分辨率，查看EVT聚类树
set.resolutions = seq(0, 1.2, by = 0.1)
EVT_for_clustree <- FindClusters(object = EVT,
                                 resolution = set.resolutions,
                                 verbose = TRUE) 
clustree(EVT_for_clustree)
  ##得到1.EVT细胞分群的聚类树
  ##根据以上结果，确定resolutions = 0.1比较合适

#按resolutions = 0.1对EVT进行分群
rm(EVT_for_clustree)
EVT <- FindClusters(EVT, resolution = 0.1)
EVT <- RunUMAP(EVT, dims = 1:10)
new.cluster.ids <- c("0" = "EVT-1",
                     "1" = "EVT-2",
                     "2" = "EVT-3",
                     "3" = "EVT-4",
                     "4" = "EVT-5")
names(new.cluster.ids) <- levels(EVT)
EVT <- RenameIdents(EVT, new.cluster.ids)
EVT@meta.data$EVTtype <- Idents(EVT)
DimPlot(EVT, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  ##得到2.resolutions = 0.1时的EVT细胞分群图
DimPlot(EVT, reduction = "umap", label = TRUE, pt.size = 0.5,
        split.by = "group",ncol = 3) + NoLegend()
  ##得到3.resolutions = 0.1时的EVT细胞分群图（各组）

#显示EVT细胞亚群3D分析结果
EVT <- RunTSNE(EVT, dims = 1:40,dim.embed = 3)
data <- as.data.frame(EVT@reductions$tsne@cell.embeddings)
celltype <- as.data.frame(EVT@meta.data$EVTtype)
df <- cbind(data, celltype)
colnames(df)[4] <- 'celltype'
  ##从Seurat对象中提取出用于3D绘图的数据
scatterplot3d(df[,1:3])  ##初步得到3D图像
unique(df$celltype) 
df <- df %>% mutate(colour = case_when(
  df$celltype == "EVT-1" ~ "#DC050C",
  df$celltype == "EVT-2" ~ "#FB8072",
  df$celltype == "EVT-3" ~ "#1965B0",
  df$celltype == "EVT-4" ~ "#7BAFDE",
  df$celltype == "EVT-5" ~ "#ffff00"
))
##设置3D图的颜色
scatterplot3d(df[,1:3],color=df$colour,
              pch = 16,angle=30,
              cex.symbols = 0.5,
              type="p",
              lty.hide=2,lty.grid = 2,
              grid=T, box=FALSE)
legend("right",c("EVT-1", "EVT-2", "EVT-3", "EVT-4", "EVT-5"),
       fill=c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#ffff00"),
       horiz = F)
  ##得到4.EVT细胞注释图（3D图像）
saveRDS(EVT,"./Seurat_data（亚群分析）/EVT.rds")

#把EVT的亚群分析写入object中
Idents(object, cells = colnames(EVT)) <- Idents(EVT)
DimPlot(object, reduction = "umap", label = TRUE, 
        pt.size = 0.5,raster=FALSE) +
  NoLegend()
}

#下面处理VCT
{
#筛选出VCT细胞
VCT <- subset(object,celltype == "VCT")
DimPlot(VCT,label = TRUE,raster=FALSE)|
  DimPlot(object,label = TRUE,raster=FALSE)

#VCT细胞的再次降维
VCT <- NormalizeData(VCT,
                     normalization.method = "LogNormalize",scale.factor = 10000) 
  ##NormalizeData进行标准化
VCT <- FindVariableFeatures(VCT,
                            selection.method = "vst", nfeatures = 2000)
  ##FindVariableFeatures计算高变基因
VCT <- ScaleData(VCT, features = rownames(VCT))
VCT <- RunPCA(VCT, features = VariableFeatures(object = VCT))
VCT <- RunUMAP(VCT, dims = 1:20, reduction = "harmony",
               verbose = TRUE)
VCT <- FindNeighbors(VCT, reduction = "harmony",dims = 1:20,
                     verbose = TRUE) 

#设置不同的分辨率，查看VCT聚类树
set.resolutions = seq(0, 1.2, by = 0.1)
VCT_for_clustree <- FindClusters(object = VCT,
                                 resolution = set.resolutions,
                                 verbose = TRUE) 
clustree(VCT_for_clustree)
  ##得到1.VCT细胞分群的聚类树
  ##根据以上结果，确定resolutions = 0.2比较合适

#按resolutions = 0.2进行分群
rm(VCT_for_clustree)
VCT <- FindClusters(VCT, resolution = 0.2)
VCT <- RunUMAP(VCT, dims = 1:10)
new.cluster.ids <- c("0" = "VCT-1",
                     "1" = "VCT-2",
                     "2" = "VCT-3",
                     "3" = "VCT-4",
                     "4" = "VCT-5",
                     "5" = "VCT-6")
names(new.cluster.ids) <- levels(VCT)
VCT <- RenameIdents(VCT, new.cluster.ids)
VCT@meta.data$VCTtype <- Idents(VCT)
DimPlot(VCT, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  ##得到2.resolutions = 0.2时的VCT细胞分群图
DimPlot(VCT, reduction = "umap", label = TRUE, pt.size = 0.5,
        split.by = "group",ncol = 3) + NoLegend()
  ##得到3.resolutions = 0.2时的VCT细胞分群图（各组）

#显示SCT细胞亚群分析3D结果
rm(new.cluster.ids,set.resolutions)
VCT <- RunTSNE(VCT, dims = 1:40,dim.embed = 3)
data <- as.data.frame(VCT@reductions$tsne@cell.embeddings)
celltype <- as.data.frame(VCT@meta.data$VCTtype)
df <- cbind(data, celltype)
colnames(df)[4] <- 'celltype'
  ##从Seurat对象中提取出用于3D绘图的数据
scatterplot3d(df[,1:3])  ##初步得到3D图像
unique(df$celltype) 
df <- df %>% mutate(colour = case_when(
  df$celltype == "VCT-1" ~ "#DC050C",
  df$celltype == "VCT-2" ~ "#FB8072",
  df$celltype == "VCT-3" ~ "#1965B0",
  df$celltype == "VCT-4" ~ "#7BAFDE",
  df$celltype == "VCT-5" ~ "#ffff00",
  df$celltype == "VCT-6" ~ "#56ff0d",
))
  ##设置3D图的颜色
scatterplot3d(df[,1:3],color=df$colour,
              pch = 16,angle=30,
              cex.symbols = 0.5,
              type="p",
              lty.hide=2,lty.grid = 2,
              grid=T, box=FALSE)
legend("right",c("VCT-1", "VCT-2", "VCT-3", "VCT-4", "VCT-5","VCT-6"),
       fill=c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#ffff00","#56ff0d"),
       horiz = F)
  ##得到4.VCT细胞注释图（3D图像）

saveRDS(VCT,"./Seurat_data（亚群分析）/VCT.rds")

#把VCT的亚群分析写入object中
Idents(object, cells = colnames(VCT)) <- Idents(VCT)
DimPlot(object, reduction = "umap", label = TRUE, 
        pt.size = 0.5,raster=FALSE) +
  NoLegend()
}

#下面处理SCT
{
#筛选出SCT细胞
SCT <- subset(object,celltype == "SCT")
DimPlot(SCT,label = TRUE,raster=FALSE)|
  DimPlot(object,label = TRUE,raster=FALSE)

#SCT细胞的再次降维
SCT <- NormalizeData(SCT,
                     normalization.method = "LogNormalize",scale.factor = 10000) 
##NormalizeData进行标准化
SCT <- FindVariableFeatures(SCT,
                            selection.method = "vst", nfeatures = 2000)
##FindVariableFeatures计算高变基因
SCT <- ScaleData(SCT, features = rownames(SCT))
SCT <- RunPCA(SCT, features = VariableFeatures(object = SCT))
SCT <- RunUMAP(SCT, dims = 1:20, reduction = "harmony",
               verbose = TRUE)
SCT <- FindNeighbors(SCT, reduction = "harmony",dims = 1:20,
                     verbose = TRUE) 

#设置不同的分辨率，查看VCT聚类树
set.resolutions = seq(0, 1.2, by = 0.1)
SCT_for_clustree <- FindClusters(object = SCT,
                                 resolution = set.resolutions,
                                 verbose = TRUE) 
clustree(SCT_for_clustree)
  ##得到1.SCT细胞分群的聚类树
  ##根据以上结果，确定resolutions = 0.2比较合适

#按resolutions = 0.2进行分群
rm(SCT_for_clustree)
SCT <- FindClusters(SCT, resolution = 0.2)
SCT <- RunUMAP(SCT, dims = 1:10)
new.cluster.ids <- c("0" = "SCT-1",
                     "1" = "SCT-2",
                     "2" = "SCT-3",
                     "3" = "SCT-4",
                     "4" = "SCT-5")
names(new.cluster.ids) <- levels(SCT)
SCT <- RenameIdents(SCT, new.cluster.ids)
SCT@meta.data$SCTtype <- Idents(SCT)
DimPlot(SCT, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  ##得到2.resolutions = 0.2时的SCT细胞分群图
DimPlot(SCT, reduction = "umap", label = TRUE, pt.size = 0.5,
        split.by = "group",ncol = 3) + NoLegend()
  ##得到3.resolutions = 0.2时的SCT细胞分群图（各组）

#显示SCT细胞亚群分析3D结果
rm(new.cluster.ids,set.resolutions)
SCT <- RunTSNE(SCT, dims = 1:40,dim.embed = 3)
data <- as.data.frame(SCT@reductions$tsne@cell.embeddings)
celltype <- as.data.frame(SCT@meta.data$SCTtype)
df <- cbind(data, celltype)
colnames(df)[4] <- 'celltype'
  ##从Seurat对象中提取出用于3D绘图的数据
scatterplot3d(df[,1:3])  ##初步得到3D图像
unique(df$celltype) 
df <- df %>% mutate(colour = case_when(
  df$celltype == "SCT-1" ~ "#DC050C",
  df$celltype == "SCT-2" ~ "#FB8072",
  df$celltype == "SCT-3" ~ "#1965B0",
  df$celltype == "SCT-4" ~ "#7BAFDE",
  df$celltype == "SCT-5" ~ "#ffff00"
))
  ##设置3D图的颜色
scatterplot3d(df[,1:3],color=df$colour,
              pch = 16,angle=30,
              cex.symbols = 0.5,
              type="p",
              lty.hide=2,lty.grid = 2,
              grid=T, box=FALSE)
legend("right",c("SCT-1", "SCT-2", "SCT-3", "SCT-4", "SCT-5"),
       fill=c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#ffff00"),
       horiz = F)
  ##得到4.SCT细胞注释图（3D图像）

saveRDS(SCT,"./Seurat_data（亚群分析）/SCT.rds")

#把VCT的亚群分析写入object中
Idents(object, cells = colnames(SCT)) <- Idents(SCT)
DimPlot(object, reduction = "umap", label = TRUE, 
        pt.size = 0.5,raster=FALSE) +
  NoLegend()
}

#创建滋养细胞的Seurat对象
Trophoblast <- subset(object,
                      celltype =="EVT"|celltype =="SCT"|celltype =="VCT")
Idents(Trophoblast)
Trophoblast@meta.data$celltype <- Idents(Trophoblast)
Idents(Trophoblast) <- Trophoblast@meta.data$celltype
DimPlot(Trophoblast, reduction = "umap", label = TRUE, 
        pt.size = 0.5,raster=FALSE) +
  NoLegend()
  ##得到1.滋养细胞的亚群分析图
DimPlot(Trophoblast, reduction = "umap", label = FALSE, 
        pt.size = 0.5,raster=FALSE)
  ##得到2.滋养细胞的亚群分析图（不带标签）
DimPlot(Trophoblast, reduction = "umap", label = TRUE, 
        pt.size = 0.5,split.by = "group",raster=FALSE) +
  NoLegend()
  ##得到3.滋养细胞的亚群分析图(各组)
DimPlot(Trophoblast, reduction = "umap", label = FALSE, 
        pt.size = 0.5,split.by = "group",raster=FALSE)
  ##得到4.滋养细胞的亚群分析图(各组+不带标签)
saveRDS(Trophoblast,"./Seurat_data（亚群分析）/滋养细胞.rds")
