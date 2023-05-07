#2022-09-05
#ywy

setwd("/Users/ywy/weiyun/project2022/2_SongYiDan/Results/3_Modules/")
BiocManager::install("ggpubr")
install.packages("tidyverse")
library(stringr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

SEMP <- read.table("SEModule1",sep='\t',header = T)
MEMP <- read.table("MEModule1",sep='\t',header = T)

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
Pdata <- data.frame(value=c(SEMP$C,MEMP$C),group = c(rep("0",nrow(SEMP)),rep("1",nrow(MEMP))))
BoxplotForP(Pdata,"C_compare.pdf")


##Enrichment analysis

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


write.csv(as.data.frame(edox_SE_Unique),file = "KEGG_SENet.csv",row.names =FALSE)
write.csv(as.data.frame(edox_ME_Unique),file = "KEGG_MENet.csv",row.names =FALSE)



