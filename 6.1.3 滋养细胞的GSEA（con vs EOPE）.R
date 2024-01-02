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

#定义用于GSEA分析的函数
GSEA <- function(gene,
                 LogFC,
                 analysis=c('GO',"KEGG"),
                 package=c('clusterProfiler','fgsea'),
                 OrgDb=c("org.Hs.eg.db", "org.Mm.eg.db")){
  if(OrgDb=="org.Hs.eg.db"){
    organism = "hsa"
    species = "Homo sapiens"
  }
  
  if(OrgDb=="org.Mm.eg.db"){
    organism = "mmu"
    species = "Mus musculus"
  }
  
  if(package=="clusterProfiler"){
    entrezID <- bitr(gene, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = OrgDb)#genesymbol转化为ID
    
    genelist <- LogFC
    names(genelist) <- gene
    genelist <- genelist[names(genelist) %in% entrezID[,1]]
    names(genelist) <- entrezID[match(names(genelist),entrezID[,1]),2]
    genelist <- sort(genelist,decreasing = T)
    
    #install.packages('R.utils')
    R.utils::setOption( "clusterProfiler.download.method",'auto')
    
    if(analysis=="KEGG"){
      KEGG_gesa <- gseKEGG(geneList = genelist,
                           organism = organism,#不同物种查询：https://www.genome.jp/kegg/catalog/org_list.html
                           minGSSize = 10,
                           maxGSSize = 500,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           verbose = FALSE,
                           eps = 0)
      
      KEGG_gesa_result = setReadable(KEGG_gesa, OrgDb = OrgDb, keyType = "ENTREZID")
      
      KEGG_gesa_table <- KEGG_gesa_result@result
      write.csv(KEGG_gesa_table, file = './KEGG_gesa_table.csv')
      
      
      return(KEGG_gesa_result)
      
    }
    else{
      GO_gesa <- gseGO(geneList = genelist,
                       ont = "BP",
                       OrgDb=OrgDb,
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       verbose = FALSE)
      
      GO_gesa_result = setReadable(GO_gesa,OrgDb = OrgDb, keyType = "ENTREZID")
      
      GO_gesa_table <- GO_gesa_result@result 
      write.csv(GO_gesa_table, file = './GO_gesa_table.csv')
      
      return(GO_gesa_result)
      
    }
    
    
  }
  
  if(package=="fgsea"){
    if(analysis=="KEGG"){
      
      geneset_KEGG = msigdbr(species = species,#mouse:Mus musculus
                             category = "C2", 
                             subcategory = "CP:KEGG") %>% dplyr::select(gs_name,gene_symbol)
      geneset_KEGG$gs_name <- gsub('KEGG_','',geneset_KEGG$gs_name)#去除前缀KEGG_
      geneset_KEGG$gs_name <- tolower(geneset_KEGG$gs_name)#将大写换为小写
      geneset_KEGG$gs_name <- gsub('_',' ',geneset_KEGG$gs_name)#将_转化为空格
      library(Hmisc)
      geneset_KEGG$gs_name <- capitalize(geneset_KEGG$gs_name)#首字母大写
      GSEA_geneset <- geneset_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
      
    }else{
      
      geneset_GO = msigdbr(species = species,#mouse:Mus musculus
                           category = "C5", 
                           subcategory = "GO:BP") %>% dplyr::select(gs_name,gene_symbol)
      
      geneset_GO$gs_name <- gsub('GOBP_','',geneset_GO$gs_name)#去除前缀KEGG_
      geneset_GO$gs_name <- tolower(geneset_GO$gs_name)#将大写换为小写
      geneset_GO$gs_name <- gsub('_',' ',geneset_GO$gs_name)#将_转化为空格
      library(Hmisc)
      geneset_GO$gs_name <- capitalize(geneset_GO$gs_name)#首字母大写
      GSEA_geneset <- geneset_GO %>% split(x = .$gene_symbol, f = .$gs_name)
    }
    
    df <- data.frame(logFC = LogFC, gene = gene)
    df <- df[order(df$logFC,decreasing = T),]
    ranks <- df$logFC
    names(ranks) <- df$gene
    
    ## GSEA分析
    GSEA_df <- fgsea(pathways = GSEA_geneset, 
                     stats = ranks,
                     minSize=10,
                     maxSize=500,
                     eps=0.0)
    library(data.table)
    fwrite(GSEA_df, file="./GSEA_df.txt", sep="\t", sep2=c("", " ", ""))
    
    return(GSEA_df)
    
  }
  
}

#定义用于GSEA分析结果可视化的函数
GSEA_plot <- function(inputType=c('clusterProfiler','fgsea'),
                      analysis=c('GO',"KEGG"),
                      data,
                      term,
                      gene,
                      LogFC,
                      OrgDb){
  
  if(OrgDb=="org.Hs.eg.db"){
    species = "Homo sapiens"
  }
  
  if(OrgDb=="org.Mm.eg.db"){
    species = "Mus musculus"
  }
  
  
  
  if(inputType=='clusterProfiler'){
    #clusterprofile
    gseaScores <- getFromNamespace("gseaScores", "DOSE")
    
    # define function
    gsInfo <- function(object, terms) {
      geneList <- object@geneList
      
      gsea_result_table <- object@result
      site=which(gsea_result_table$Description==terms, arr.ind = TRUE)
      genesetid=gsea_result_table$ID[site]
      
      geneSetID <- object@result[genesetid, "ID"]
      
      geneSet <- object@geneSets[[geneSetID]]
      exponent <- object@params[["exponent"]]
      df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
      df$ymin <- 0
      df$ymax <- 0
      pos <- df$position == 1
      h <- diff(range(df$runningScore))/20
      df$ymin[pos] <- -h
      df$ymax[pos] <- h
      df$geneList <- geneList
      
      df$Description <- object@result[geneSetID, "Description"]
      return(df)
    }
    
    gsea_plotData = gsInfo(data, terms = term)
    
    gsea_plotData <- gsea_plotData %>%
      mutate("gene_name" = data@gene2Symbol) %>%
      filter(position == 1) 
    
    colnames(gsea_plotData)[2] <- 'y'
    
    test = gsea_plotData
    
    #NES和adjust Pvalue
    gesa_table <- data@result
    terms = gesa_table[gesa_table$Description==term,]
    labels_NES = round(terms$NES, 4)
    labels_FDR = format(terms$qvalue, scientific = T,digits = 2)
    
    
    
  }
  
  if(inputType=='fgsea'){
    if(analysis=="KEGG"){
      geneset_KEGG = msigdbr(species = species,#mouse:Mus musculus
                             category = "C2", 
                             subcategory = "CP:KEGG") %>% dplyr::select(gs_name,gene_symbol)
      geneset_KEGG$gs_name <- gsub('KEGG_','',geneset_KEGG$gs_name)#去除前缀KEGG_
      geneset_KEGG$gs_name <- tolower(geneset_KEGG$gs_name)#将大写换为小写
      geneset_KEGG$gs_name <- gsub('_',' ',geneset_KEGG$gs_name)#将_转化为空格
      library(Hmisc)
      geneset_KEGG$gs_name <- capitalize(geneset_KEGG$gs_name)#首字母大写
      GSEA_geneset <- geneset_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
      
    }else{
      geneset_GO = msigdbr(species = species,#mouse:Mus musculus
                           category = "C5", 
                           subcategory = "GO:BP") %>% dplyr::select(gs_name,gene_symbol)
      
      geneset_GO$gs_name <- gsub('GOBP_','',geneset_GO$gs_name)#去除前缀KEGG_
      geneset_GO$gs_name <- tolower(geneset_GO$gs_name)#将大写换为小写
      geneset_GO$gs_name <- gsub('_',' ',geneset_GO$gs_name)#将_转化为空格
      library(Hmisc)
      geneset_GO$gs_name <- capitalize(geneset_GO$gs_name)#首字母大写
      GSEA_geneset <- geneset_GO %>% split(x = .$gene_symbol, f = .$gs_name)
      
      
    }
    
    
    df <- data.frame(logFC = LogFC, gene = gene)
    df <- df[order(df$logFC,decreasing = T),]
    ranks <- df$logFC
    names(ranks) <- df$gene
    
    
    
    ToPlot  <- sapply(data$pathway, function(Pathway){
      pathway <- GSEA_geneset[[Pathway]]
      stats <- ranks
      rnk <- rank(-stats)
      ord <- order(rnk)
      statsAdj <- stats[ord]
      statsAdj <- sign(statsAdj)*(abs(statsAdj)^1)
      statsAdj <- statsAdj/max(abs(statsAdj))
      pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
      pathway <- sort(pathway)
      gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
      bottoms <- gseaRes$bottoms
      tops <- gseaRes$tops
      n <- length(statsAdj)
      xs <- as.vector(rbind(pathway - 1, pathway))
      ys <- as.vector(rbind(bottoms, tops))
      toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0), pathway = Pathway)
      return(list(toPlot))
    })
    
    #点图数据
    gsea_plotData <- ToPlot[[term]]
    
    test <-  gsea_plotData[1:floor(nrow(gsea_plotData)/length(1)),]
    #NES和adjust Pvalue
    terms = data[data$pathway==term,]
    labels_NES = round(terms$NES, 4)
    labels_FDR = format(terms$padj, scientific = T,digits = 2)
    
  }
  
  
  
  test$xend = test$x+1
  test$xend2 = ''
  test$xend3 = ''
  
  for (i in 1:nrow(test)) {
    test$xend2[i] = test$xend[i+1] - test$xend[i]
  }
  
  for (i in 1:ncol(test)){
    test[,i][is.na(test[,i])] <- 5
  } 
  
  test$xend2 <- as.numeric(test$xend2)
  
  for (i in 1:nrow(test)) {
    test$xend3[i] = (test$x[i] + test$xend2[i])
  }
  
  test$xend3 <- as.numeric(test$xend3)
  
  
  if(analysis=="KEGG"){
    title_terms = paste0("KEGG:",term)
  }
  
  if(analysis=="GO"){
    title_terms = paste0("GO:",term)
  }
  
  
  if(labels_NES>0){
    fill_color = "#90191B"
  }
  
  if(labels_NES<0){
    fill_color = "#3F90C9"
  }
  
  
  p=ggplot(test) + 
    geom_rect(aes(xmin = x-1,xmax = xend3, ymin = -0.04 , ymax = 0.04, fill=x), lwd=4)+
    scale_fill_gradientn(colours = colorRampPalette(c("#90191B","white","#3F90C9"))(100))+
    geom_rect(aes(xmin = x,xmax = xend, ymin = -0.04, ymax = 0.04, fill=x), color="black", size=0.5)+
    geom_point(data=gsea_plotData, aes(x = x, y = y),fill=fill_color, shape=21, size=4) + 
    geom_hline(yintercept = 0, linetype=3, lwd = 1) +
    scale_x_continuous(expand = c(0.01,0.01))+
    ylab("Enrichment Score") +
    xlab('')+
    labs(title=title_terms)+
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black",size = 1),
          axis.text.x=element_blank(),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 15,colour = 'black'),
          legend.position="none",
          plot.title = element_text(size=15,hjust =0.5, face = 'bold'),
          axis.ticks.x = element_blank())
  
  if(labels_NES>0){
    lable_y = max(gsea_plotData$y)-0.2
    
    a=nrow(gsea_plotData)
    
    lable_x = gsea_plotData$x[a]-gsea_plotData$x[ceiling(a/1.8)]
    
    ylim = expand_limits(y=c(-0.2, NA))
  }
  
  if(labels_NES<0){
    
    lable_y = min(gsea_plotData$y)+0.2
    a=nrow(gsea_plotData)
    
    lable_x = gsea_plotData$x[a]-gsea_plotData$x[ceiling(a/2)]
    ylim = expand_limits(y=c(NA, 0.2))
    
  }
  
  p+ylim+annotate(geom = 'text', label=paste("NES =",labels_NES, "\n", "FDR =", labels_FDR), 
                  x=lable_x, y=lable_y, size=4)
  
  
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

#KEGG的GSEA分析
KEGG_GSEA <- GSEA(gene = DEG$SYMBOL,
                  LogFC = DEG$avg_log2FC,
                  analysis = "KEGG",
                  package = 'clusterProfiler',
                  OrgDb = 'org.Hs.eg.db')
write.csv(KEGG_GSEA@result,
          "./result/6.1.3 滋养细胞的GSEA（con vs EOPE）/1.滋养细胞的GSEA_KEGG（con vs EOPE）.csv")
dotplot(KEGG_GSEA,split=".sign")+facet_wrap(~.sign,scales = "free")
  ##得到2.滋养细胞的GSEA_KEGG（con vs EOPE）（点图）
gseaplot2(KEGG_GSEA,c(1:10),color="green") 
  ##得到3.滋养细胞的GSEA_KEGG（总折线图）
for (i in 1:length(KEGG_GSEA@result$Description)) {
  plot1 <- gseaplot(KEGG_GSEA,geneSetID = i,
                    title = KEGG_GSEA$Description[i],
                    pvalue_table = T)
  ggsave(filename =  paste0(KEGG_GSEA$Description[i],".pdf"),
         plot = plot1[[1]]/plot1[[2]],height = 10,width = 10)
}
  ##得到4.滋养细胞的GSEA_KEGG（各折线图）

#GO的GSEA分析
GO_GSEA <- GSEA(gene = DEG$SYMBOL,
                LogFC = DEG$avg_log2FC,
                analysis = "GO",
                package = 'clusterProfiler',
                OrgDb = 'org.Hs.eg.db')
write.csv(GO_GSEA@result,
          "./result/6.1.3 滋养细胞的GSEA（con vs EOPE）/5.滋养细胞的GSEA_GO（con vs EOPE）.csv")
dotplot(GO_GSEA,split=".sign")+facet_wrap(~.sign,scales = "free")
  ##得到6.滋养细胞的GO_KEGG（con vs EOPE）（点图）
gseaplot2(GO_GSEA,c(1:10),color="green") 
  ##得到7.滋养细胞的GO_KEGG（总折线图）
for (i in 1:length(GO_GSEA@result$Description)) {
  plot1 <- gseaplot(GO_GSEA,geneSetID = i,
                    title = GO_GSEA$Description[i],
                    pvalue_table = T)
  ggsave(filename =  paste0(GO_GSEA$Description[i],".pdf"),
         plot = plot1[[1]]/plot1[[2]],height = 10,width = 10)
}
  ##得到8.滋养细胞的GO_KEGG（各折线图）
