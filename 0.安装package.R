rm(list = ls())
setwd("D:/2023.10/PE scRNA（多数据库运行）/R（code）")
getwd()

#注：
#一些难装的package需要多行代码完成，新电脑第一次安装时可以在这里进行

#CellChat的安装
if(!require(NMF))install.packages('NMF')
if(!require(circlize))devtools::install_github("jokergoo/circlize")
if(!require(ComplexHeatmap))devtools::install_github("jokergoo/ComplexHeatmap") 
if(!require(CellChat))devtools::install_github("sqjin/CellChat")
