library(stringr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

NetP <- read.csv("resultTable.csv")
n <- nrow(NetP)

BoxplotForP <-  function(Pdata,outFiles){
  p2=ggplot(Pdata,aes(x=group,y=value,fill=group))+ 
    stat_boxplot(mapping=aes(x=group,y=value),
                 geom ="errorbar",                          
                 width=0.2,position=position_dodge(0.5))+     
    geom_boxplot(aes(fill=group),                            
                 position=position_dodge(0.8),                
                 width=0.55)+                    
    stat_summary(mapping=aes(group=group),                    
                 fun="mean",                                   
                 geom="point",shape=23,size=3,fill="white",    
                 position=position_dodge(0.8))+               
    scale_fill_manual(values=c('#104E8B',"#CD2626"))+theme_bw()+
    theme(panel.grid = element_blank(),axis.text=element_text(face = "bold",size = 14,colour = "black"), legend.position = 'none')+
    labs(x="",y="")+
    stat_compare_means(aes(group=group),method = "wilcox.test",label="p.signif",comparisons = list(c("0","1")),label.x.npc = 0.5,exact = FALSE)+
    scale_x_discrete(labels = c('SE','ME'))
  p2
  ggsave(p2,filename = outFiles,width = 5,height = 5)
  
}

#K
Pdata <- data.frame(value=c(NetP$sindegree_w,NetP$duodegree_w),group = c(rep("0",n),rep("1",n)))
BoxplotForP(Pdata,"K_compare.pdf")

#B
Pdata <- data.frame(value=c(NetP$sin_betweeness,NetP$duo_betweeness),group = c(rep("0",n),rep("1",n)))
BoxplotForP(Pdata,"B_compare.pdf")

#B
Pdata <- data.frame(value=c(NetP$sin_betweeness.0.1,NetP$duo_betweeness.0.1),group = c(rep("0",n),rep("1",n)))
BoxplotForP(Pdata,"B_compare_nor01.pdf")


#C
Pdata <- data.frame(value=c(NetP$sin_closeness.0.1,NetP$duo_closeness.0.1),group = c(rep("0",n),rep("1",n)))
BoxplotForP(Pdata,"C_compare_nor01.pdf")

#L
Pdata <- data.frame(value=c(NetP$sin_L,NetP$duo_L),group = c(rep("0",n),rep("1",n)))
BoxplotForP(Pdata,"L_compare_nor01.pdf")


#EC: Eigenvector centrality 
Pdata <- data.frame(value=c(NetP$sin_enig,NetP$duo_enig),group = c(rep("0",n),rep("1",n)))
BoxplotForP(Pdata,"EC_compare_nor01.pdf")




#enrichment analysis by clusterProfiler
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggnewscale)
library(ggupset)
library(R.utils)

KEGGFunction <- function(GeneList,outputFile,cutoff) {
  ConvertList <-bitr(GeneList, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db") #GeneList is t     he symbol0
  #KEGG
  R.utils::setOption("clusterProfiler.download.method",'auto')
  kk <- enrichKEGG(gene = ConvertList$ENTREZID, organism = "hsa", qvalueCutoff = cutoff)
  #set the gene readble
  edox <- setReadable(kk, OrgDb=org.Hs.eg.db, 'ENTREZID')
  write.csv(summary(edox),outputFile,row.names =FALSE)
  return(edox)
} 

#enrichment analysis for top 5% degree genes in SENet and MENet
cutoff <- 0.05

#1. pathways for SENet
GeneList_SE <- read.table(paste(getwd(),"/SE_K_TOP5.txt",sep=''),sep='\t',header = F)[,1]
outputFile_SE <- "KEGG_SENet_K_TOP5.csv"
edox_SE <- KEGGFunction(GeneList_SE,outputFile_SE,cutoff)

#2. pathways for MENet
GeneList_ME <- read.table(paste(getwd(),"/ME_K_TOP5.txt",sep=''),sep='\t',header = F)[,1]
outputFile_ME <- "KEGG_MENet_K_TOP5.csv"
edox_ME <- KEGGFunction(GeneList_ME,outputFile_ME,cutoff)

edoxDF_SE <- as.data.frame(edox_SE)
edoxDF_ME <- as.data.frame(edox_ME)

overlapPath <- intersect(edoxDF_SE$ID,edoxDF_ME$ID)#overlap pathways

edoxDF_SE_Unique <- edoxDF_SE[!(edoxDF_SE$ID %in% overlapPath),] # SE specific pathways
edoxDF_ME_Unique <- edoxDF_ME[!(edoxDF_ME$ID %in% overlapPath),] # ME specific pathways

#convert data.frame to enrichResult object
remotes::install_github("jmw86069/jamenrich")
library(multienrichjam)
edox_SE_Unique <- enrichDF2enrichResult(enrichDF = edoxDF_SE_Unique, keyColname = "ID", geneColname = "geneID", 
                                        pvalueColname = "pvalue", descriptionColname = "Description")
edox_ME_Unique <- enrichDF2enrichResult(enrichDF = edoxDF_ME_Unique, keyColname = "ID", geneColname = "geneID", 
                                        pvalueColname = "pvalue", descriptionColname = "Description")


write.csv(as.data.frame(edox_SE_Unique),file = "KEGG_SENet_K_TOP5_Unique.csv",row.names =FALSE)
write.csv(as.data.frame(edox_ME_Unique),file = "KEGG_MENet_K_TOP5_Unique.csv",row.names =FALSE)


