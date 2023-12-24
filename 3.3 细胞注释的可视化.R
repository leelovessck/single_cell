rm(list = ls())
setwd("E:/2023.10/PE scRNA（多数据库运行）")
getwd()

#打开必要的package
{
if(!require(Seurat))install.packages("Seurat")
if(!require(scatterplot3d))install.packages("scatterplot3d")
if(!require(dplyr))install.packages("dplyr")
if(!require(ggplot2))install.packages("ggplot2")
}

#载入seurat对象
object <- readRDS("./Seurat_data（整理）/合并（细胞注释）.rds")
object <- readRDS("./Seurat_data（整理）/合并（注释+TSNE）.rds")

#显示细胞注释结果(平面图)
{
DimPlot(object, reduction = "umap",
        label = TRUE,
        pt.size = 0,
        cols = c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                 "#FF7F00", "#E78AC3", "#33A02C", "#B2DF8A", "#A6761D",
                 "#999999", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"),
        raster=FALSE)
  ##得到1.细胞注释图
DimPlot(object, reduction = "umap",
        label = FALSE,
        pt.size = 0,
        cols = c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                 "#FF7F00", "#E78AC3", "#33A02C", "#B2DF8A", "#A6761D",
                 "#999999", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"),
        raster=FALSE)
  ##得到2.细胞注释图（不带标签）
DimPlot(object, reduction = "umap",
        label = FALSE,
        pt.size = 0,
        raster=FALSE,
        cols = c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                 "#FF7F00", "#E78AC3", "#33A02C", "#B2DF8A", "#A6761D",
                 "#999999", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"),
        split.by = "orig.ident",
        ncol = 10)
  ##得到3.细胞注释图（按样本）
DimPlot(object, reduction = "umap",
        label = FALSE,
        pt.size = 0,
        raster=FALSE,
        cols = c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                 "#FF7F00", "#E78AC3", "#33A02C", "#B2DF8A", "#A6761D",
                 "#999999", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"),
        split.by = "group",
        ncol = 3)
  ##得到4.细胞注释图（按分组）
DimPlot(object, reduction = "umap",
        label = TRUE,
        pt.size = 0,
        raster=FALSE,
        cols = c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                 "#FF7F00", "#E78AC3", "#33A02C", "#B2DF8A", "#A6761D",
                 "#999999", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"),
        split.by = "group",
        ncol = 3)
##得到5.细胞注释图（按分组+标签）
}

#显示细胞注释结果(立体图)
{
object <- RunTSNE(object, dims = 1:40,dim.embed = 3)
#saveRDS(object,"./Seurat_data（整理）/合并（注释+TSNE）.rds")
data <- as.data.frame(object@reductions$tsne@cell.embeddings)
celltype <- as.data.frame(object@meta.data$celltype)
df <- cbind(data, celltype)
colnames(df)[4] <- 'celltype'
  ##从Seurat对象中提取出用于3D绘图的数据
scatterplot3d(df[,1:3])  ##初步得到3D图像
unique(df$celltype) 
df <- df %>% mutate(colour = case_when(
  df$celltype == "EVT" ~ "#DC050C",
  df$celltype == "SCT" ~ "#FB8072",
  df$celltype == "VCT" ~ "#1965B0",
  df$celltype == "Monocyte" ~ "#7BAFDE",
  df$celltype == "T cell" ~ "#882E72",
  df$celltype == "Endothelial cell" ~ "#ffff00",
  df$celltype == "Dendritic cell" ~ "#FF7F00",
  df$celltype == "Macrophage" ~ "#E78AC3",
  df$celltype == "Hofbauer cell" ~ "#33A02C",
  df$celltype == "Stromal cell" ~ "#B2DF8A",
  df$celltype == "RBC" ~ "#A6761D",
  df$celltype == "B cell" ~ "#999999",
  df$celltype == "Decidual cell" ~ "#1e90ff",
  df$celltype == "NK cell" ~ "#00bfff",
  df$celltype == "Smooth muscle cell" ~"#56ff0d" 
))
  ##设置3D图的颜色
scatterplot3d(df[,1:3],color=df$colour,
              pch = 16,angle=45,
              cex.symbols = 0.5,
              type="p",
              lty.hide=2,lty.grid = 2,
              grid=T, box=FALSE)
  ##得到6.细胞注释图（3D图像）
legend("right",c("EVT", "SCT", "VCT", "Monocyte", "T cell"),
       fill=c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72"),
       horiz = F)
  ##得到图例1
legend("right",c("Endothelial cell", "Dendritic cell", "Macrophage",
                 "Hofbauer cell","Stromal cell"),
       fill=c("#ffff00", "#FF7F00", "#E78AC3", "#33A02C", "#B2DF8A"),
       horiz = F)
  ##得到图例2
legend("right",c("RBC", "B cell",
                 "Decidual cell", "NK cell","Smooth muscle cell"),
       fill=c("#A6761D", "#999999", "#1e90ff", "#00bfff", "#56ff0d"),
       horiz = F)
  ##得到图例3
scatterplot3d(df[,1:3],highlight.3d = TRUE,
              pch = 16,angle=45,
              cex.symbols = 0.5,
              type="p",
              lty.hide=2,lty.grid = 2,
              grid=T, box=FALSE)
  ##得到7.细胞注释图（3D图像+渐变）
}

#显示细胞marker
{
markers <- c("HLA-G","EBI3","CYP19A1","LGALS16","ERVFRD-1",
             "HLA-A","HLA-B","HLA-DOB","CDH1","EGFR","CD14",
             "CD74","CD163","RNASE1","CD52","CD83","MYH11",
             "CDH5","IGFBP1","PRL","VIM","COL1A1","TAGLN",
             "LUM","CD79A","NKG7","CD2")
DotPlot(object,features = markers)+coord_flip()
  ##得到8.细胞marker在各种细胞中的表达（点图）
VlnPlot(object,
        features = markers,pt.size = 0,
        slot = "counts",log = TRUE,
        ncol = 5,raster=FALSE)
  ##得到9.细胞marker在各种细胞中的表达（小提琴图）
}
