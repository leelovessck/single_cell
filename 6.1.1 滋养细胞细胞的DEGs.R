rm(list = ls())
setwd("E:/2023.10/PE scRNA（多数据库运行）")
getwd()

#打开必要的package
{
if(!require(Seurat))install.packages("Seurat")
if(!require(EnhancedVolcano))BiocManager::install("EnhancedVolcano")
}

#载入seurat对象
object <- readRDS("./Seurat_data（亚群分析）/滋养细胞.rds")

#下面是control vs EOPE
{
DEG_cve_Trophoblast <- FindMarkers(object, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.10,
                                   group.by = "group",
                                   ident.1 = "EOPE",
                                   ident.2 = "control")
  ##这里ident.1是实验组，ident.2是对照组
write.csv(DEG_cve_Trophoblast,
          "./result/6.1.1 滋养细胞的DEG/1.滋养细胞的DEGs（control vs EOPE）.csv")
  ##得到滋养细胞的差异表达基因（control vs EOPE）
EnhancedVolcano(DEG_cve_Trophoblast,
                lab = rownames(DEG_cve_Trophoblast),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-5,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 6.0,
                title = "滋养细胞的DEGs（control vs EOPE）")
  ##得到2.滋养细胞的DEGs（control vs EOPE）（火山图）
}

#下面是control vs LOPE
{
DEG_cvl_Trophoblast <- FindMarkers(object, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.10,
                                   group.by = "group",
                                   ident.1 = "LOPE",
                                   ident.2 = "control")
  ##这里ident.1是实验组，ident.2是对照组
write.csv(DEG_cvl_Trophoblast,
          "./result/6.1.1 滋养细胞的DEG/3.滋养细胞的DEGs（control vs LOPE）.csv")
  ##得到滋养细胞的差异表达基因（control vs LOPE）
EnhancedVolcano(DEG_cvl_Trophoblast,
                lab = rownames(DEG_cvl_Trophoblast),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-5,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 6.0,
                title = "滋养细胞的DEGs（control vs LOPE）")
  ##得到4.滋养细胞的DEGs（control vs LOPE）（火山图）
}

#下面是EOPE vs LOPE
{
DEG_evl_Trophoblast <- FindMarkers(object, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.10,
                                   group.by = "group",
                                   ident.1 = "LOPE",
                                   ident.2 = "EOPE")
  ##这里ident.1是实验组，ident.2是对照组
write.csv(DEG_evl_Trophoblast,
          "./result/6.1.1 滋养细胞的DEG/5.滋养细胞的DEGs（EOPE vs LOPE）.csv")
  ##得到滋养细胞的差异表达基因（EOPE vs LOPE）
EnhancedVolcano(DEG_evl_Trophoblast,
                lab = rownames(DEG_evl_Trophoblast),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-5,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 6.0,
                title = "滋养细胞的DEGs（EOPE vs LOPE）")
  ##得到6.滋养细胞的DEGs（EOPE vs LOPE）（火山图）
}
