rm(list = ls())
setwd("E:/2023.10/PE scRNA（多数据库运行）")
getwd()

#打开必要的package
{
if(!require(org.Hs.eg.db))BiocManager::install("org.Hs.eg.db")
if(!require(dplyr))install.packages("dplyr")
if(!require(clusterProfiler))BiocManager::install("clusterProfiler")
if(!require(ggplot2))install.packages("ggplot2")
if(!require(R.utils))install.packages("R.utils")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(stringr))install.packages("stringr")
if(!require(enrichplot))install.packages("enrichplot")
if(!require(msigdbr))install.packages("msigdbr")
if(!require(GSVA))BiocManager::install("GSVA")
if(!require(pheatmap))install.packages("pheatmap")
if(!require(limma))install.packages("limma")
if(!require(BiocParallel))install.packages("BiocParallel")
}

#定义绘制气泡图的函数
enrichplot <- function(data4plot){
  library(ggplot2)
  data4plot <- data4plot[order(data4plot$qvalue,decreasing = F)[1:15],]
  
  data4plot <- separate(data4plot,
                        col = BgRatio,
                        into = c("BgRatio_1", "BgRatio_2"),
                        sep = "/")
  data4plot$BgRatio_1 <- as.numeric(data4plot$BgRatio_1)
  data4plot$BgRatio_2 <- as.numeric(data4plot$BgRatio_2)
  data4plot$BgRatio <- data4plot$BgRatio_1 / data4plot$BgRatio_2
  
  p <- ggplot(data4plot,aes(BgRatio,Description))
  p<-p + geom_point()
  
  pbubble <- p + geom_point(aes(size=Count,color=-1*log10(qvalue)))
  
  pr <- pbubble + scale_colour_gradient(low="#90EE90",high="red") + 
    labs(color=expression(-log[10](qvalue)),size="observed.gene.count", 
         x="Richfactor", y="term.description",title="Enrichment Process")
  
  pr <- pr + theme_bw()
  pr
}

#修改下载协议
R.utils::setOption("clusterProfiler.download.method","auto")

#载入数据
object <- readRDS("./Seurat_data（亚群分析）/滋养细胞.rds")
DEG <- read.csv("./result/6.1.1 滋养细胞的DEG/1.滋养细胞的DEGs（control vs EOPE）.csv")
rownames(DEG) <- DEG[,1]
colnames(DEG)[1] <- "SYMBOL"
DEG_rich <- DEG %>% top_n(n = -200, wt = p_val_adj) %>% rownames()

#利用org.Hs.eg.db整合基因名、 ENTREZID、ENSEMBL
DEG_df <- bitr(DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
               OrgDb = org.Hs.eg.db)
DEG_df <- DEG_df %>% distinct(SYMBOL, .keep_all = T)
  ##去重操作

#KEGG富集分析
kegg <- enrichKEGG(unique(DEG_df$ENTREZID), organism='hsa',
                          pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                          minGSSize=10,maxGSSize=500,use_internal_data=F)
kegg <- setReadable(kegg,"org.Hs.eg.db","ENTREZID")  
  ##把ENTREZID转换为SYMBOL
write.csv(kegg@result,
          "./result/6.1.2 滋养细胞的GO和KEGG（con vs EOPE）/1.滋养细胞的KEGG（con vs EOPE）.csv")
enrichplot(kegg@result)
  ##得到2.滋养细胞的KEGG气泡图（con vs EOPE）

#GO富集分析（molecular function）
goMF <- enrichGO(DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goMF@result,
          "./result/6.1.2 滋养细胞的GO和KEGG（con vs EOPE）/3.滋养细胞的GO_MF（con vs EOPE）.csv")
enrichplot(goMF@result)
  ##得到4.滋养细胞的GO_MF气泡图（con vs EOPE）

#GO富集分析（cell component）
goCC <- enrichGO(DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goCC@result,
          "./result/6.1.2 滋养细胞的GO和KEGG（con vs EOPE）/5.滋养细胞的GO_CC（con vs EOPE）.csv")
enrichplot(goCC@result)
  ##得到6.滋养细胞的GO_CC气泡图（con vs EOPE）

#GO富集分析（biological process）
goBP <- enrichGO(DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goBP@result,
          "./result/6.1.2 滋养细胞的GO和KEGG（con vs EOPE）/7.滋养细胞的GO_BP（con vs EOPE）.csv")
enrichplot(goBP@result)
  ##得到8.滋养细胞的GO_BP气泡图（con vs EOPE）
