rm(list = ls())
setwd("E:/2023.10/PE scRNA（多数据库运行）")
getwd()

#打开必要的package
{
if(!require(Seurat))install.packages("Seurat")
if(!require(scatterplot3d))install.packages("scatterplot3d")
if(!require(dplyr))install.packages("dplyr")
if(!require(ggplot2))install.packages("ggplot2")
if(!require(remotes))install.packages("remotes")  
if(!require(MySeuratWrappers))remotes::install_github("lyc-1995/MySeuratWrappers")
if(!require(ggsci))install.packages("ggsci")
if(!require(paletteer))install.packages("paletteer")
if(!require(ggtree))install.packages("ggtree")
if(!require(aplot)) install.packages("aplot")
if(!require(patchwork))install.packages("patchwork")
if(!require(FlexDotPlot))devtools::install_github("Simon-Leonard/FlexDotPlot")
if(!require(forcats))install.packages("forcats")
}

#载入seurat对象
object <- readRDS("./Seurat_data（整理）/合并（注释+TSNE）.rds")

#选择细胞marker
markers <- c("HLA-G","EBI3","CYP19A1","LGALS16","ERVFRD-1",
             "HLA-A","HLA-B","HLA-DOB","CDH1","EGFR","CD14",
             "CD74","CD163","RNASE1","CD52","CD83","MYH11",
             "CDH5","IGFBP1","PRL","VIM","COL1A1","TAGLN",
             "LUM","CD79A","NKG7","CD2")

DotPlot(object,features = markers)+coord_flip()
  ##得到1.细胞marker在各种细胞中的表达（点图1）

VlnPlot(object,
        features = markers,pt.size = 0,
        slot = "counts",log = TRUE,
        ncol = 5,raster=FALSE)
  ##得到2.细胞marker在各种细胞中的表达（小提琴图1）

my33colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#585658', '#9FA3A8', '#E0D4CA', 
               '#5F3D69', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35',
               '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',
               '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')
  ##设置颜色
VlnPlot(object, features = markers,
        stacked=T,pt.size=0,
        cols = my33colors,#颜色
        direction = "horizontal", 
        x.lab = '', y.lab = '')+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())
  ##得到3.细胞marker在各种细胞中的表达（小提琴图2）

DotPlot(object, features = markers)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#330066','#336699','#66CC66','#FFCC33'))
  ##得到4.细胞marker在各种细胞中的表达（点图2）

p <- DotPlot(object, features = markers)
exp <- p$data
p1 <- ggplot(exp,aes(x=features.plot,y=id))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  theme_bw()+coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(hjust = 1,vjust=0.5))+
  scale_color_gradient(low="lightgrey",high="blue")+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
df1 <- reshape2::dcast(exp,id~features.plot,value.var = "avg.exp.scaled")
rownames(df1) <- df1[,1] 
df1 <- df1[,-1]
df1 <- t(df1)
head(df1)
p2 <- ggtree(hclust(dist(df1)))
p2+geom_tiplab()
p3 <- ggtree(hclust(dist(t(df1))))+layout_dendrogram()
p3
p1%>%insert_left(p2)
p1%>%insert_top(p3,height = 0.15)%>%insert_left(p2,width = 0.3)
  ##得到5.细胞marker在各种细胞中的表达（聚类点图1）

dp = DotPlot(object, features = markers) + RotatedAxis()
dot_plot(dp$data[,c(3,4,1,2,5)],
         size_var = "pct.exp", 
         col_var = "avg.exp.scaled",
         size_legend = "Percent Expressed", 
         col_legend = "Average Expression",
         x.lab.pos = "bottom", 
         y.lab.pos = "right",
         display_max_sizes = F, 
         shape.scale=8,
         hclust_method = "ward.D2",
         dend_x_var = c("pct.exp", "avg.exp.scaled"),
         dend_y_var = c("pct.exp", "avg.exp.scaled"),
         text.size = 0.5,
         text.vjust = 0.5,
         size.breaks.number=6,
         color.breaks.number=4,
         x.lab.size.factor = 0.8,
         cols.use = c('#330066','#336699','#66CC66','#FFCC33'),
         y.lab.size.factor = 1)
  ##得到6.细胞marker在各种细胞中的表达（聚类点图2）

markers <- as.data.frame(markers)
markerdata <- ScaleData(object,
                        features = as.character(unique(markers$markers)),
                        assay = "RNA")
p <- DoHeatmap(markerdata,
               features = as.character(unique(markers$markers)),
               group.by = "celltype",
               assay = 'RNA')
DoHeatmap(markerdata,
          features = as.character(unique(markers$markers)),
          group.by = "celltype",
          assay = 'RNA',
          group.colors = c("#00BFC4","#AB82FF","#00CD00","#C77CFF"))+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
  ##得到7.细胞marker在各种细胞中的表达（热图）
