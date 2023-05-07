library(devtools)
library(reshape2)
library(ggplot2)
library(factoextra)
library(stringr)
library(ggsci)
library(ggbiplot)
setwd("F:/RDX/")
data=read.csv("RegLogEx.csv",header = T,row.names = 1)
data=data[rowSums(data==0)<round(ncol(data)*0.8),] 
exprData<-data
data_t <- t(data)
variableL <- ncol(data_t)
pca <- prcomp(data_t, scale=T)
print(str(pca))
fviz_eig(pca, addlabels = TRUE)
fviz_pca_ind(pca, repel=T) 
A=stringr::str_extract(rownames(data_t),"[A-Za-z0-9]*")
A[A=="Control"]="Ctrl"
#A[A=="R1"]="R1-05"
#A[A=="R25"]="R25-05"
A=factor(A,levels = c("Ctrl","SE","ME"))

p <- fviz_pca_ind(pca, col.ind=A, mean.point=F, addEllipses = T, legend.title="",
             geom = "point",palette= c("#008B00",'#104E8B',"#CD2626"),title = "")
p1 <- p + theme_bw() + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), legend.text = element_text(size=20))
ggsave(p1,filename = "PCA.pdf") 



