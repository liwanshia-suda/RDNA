library(stringr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

SENetCoreP <- read.table("SENetCore",sep='\t',header = T)
MENetCoreP <- read.table("MENetCore",sep='\t',header = T)

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

Normal <- function(x){
  s <- (x-min(x))/(max(x)-min(x))
  return(s)
}

ScaleNor <- function(x){
  s  <- (x-min(x))/sd(x)
  return(s)
}

#EC
Pdata <- data.frame(value=c(SENetCoreP$EC,MENetCoreP$EC),group = c(rep("0",nrow(SENetCoreP)),rep("1",nrow(MENetCoreP))))
BoxplotForP(Pdata,"NetCore_EC_compare.pdf")

#SENetCoreP: EC mean:0.04017024; median:0.00000
#MENetCoreP: EC mean:0.044565  ; median:0.001426
#just C and EC have significantly difference. 


#RNs: not used; no significance

SENetCoreP$RNs2 <- SENetCoreP$K*SENetCoreP$Rs*Normal(SENetCoreP$C)
MENetCoreP$RNs2 <- MENetCoreP$K*MENetCoreP$Rs*Normal(MENetCoreP$C)
Pdata <- data.frame(value=c(SENetCoreP$RNs2,MENetCoreP$RNs2),group = c(rep("0",nrow(SENetCoreP)),rep("1",nrow(MENetCoreP))))
BoxplotForP(Pdata,"NetCore_RNs3_compare.pdf")


SENetCoreP$RNs2 <- ScaleNor(SENetCoreP$EC*SENetCoreP$Rs*SENetCoreP$C)
MENetCoreP$RNs2 <- ScaleNor(MENetCoreP$EC*MENetCoreP$Rs*MENetCoreP$C)

SENetCoreP$RNs2 <- log10(SENetCoreP$B + SENetCoreP$Rs*100 + SENetCoreP$EC*100)
MENetCoreP$RNs2 <- log10(MENetCoreP$B + MENetCoreP$Rs*100 + MENetCoreP$EC*100)

SENetCoreP$RNs2 <- log10(SENetCoreP$B + SENetCoreP$Rs*0.1) * (SENetCoreP$EC + 0.0001*SENetCoreP$Rs)
MENetCoreP$RNs2 <- log10(MENetCoreP$B + MENetCoreP$Rs*0.1) * (MENetCoreP$EC + 0.0001*MENetCoreP$Rs)

S_RNs2 <- log10(SENetCoreP$RNs2[which(SENetCoreP$RNs2 > 1)])
M_RNs2 <- log10(MENetCoreP$RNs2[which(MENetCoreP$RNs2 > 1)])

Pdata <- data.frame(value=c(S_RNs2,M_RNs2),group = c(rep("0",length(S_RNs2)),rep("1",length(M_RNs2))))
BoxplotForP(Pdata,"NetCore_RNs2_2_compare.pdf")

#######################################################
library(dplyr)#2022-09-13 add the DEGs and module information to network core
DESE <- read.csv("DEG_0.05_SE.csv")
SE_Mod1 <- read.table("/Users/ywy/weiyun/project2022/2_SongYiDan/Results/3_Modules/round2/SEModule1",sep='\t',header = T)
SENetCoreP2  <- left_join(SENetCoreP, DESE, by = "ID")
SENetCoreP2$Module1 <- SENetCoreP2$ID %in% SE_Mod1$ID
write.csv(SENetCoreP2,file = "SENetCoreNodeInfor.csv")

DEME <- read.csv("DEG_0.05_ME.csv")
ME_Mod1 <- read.table("/Users/ywy/weiyun/project2022/2_SongYiDan/Results/3_Modules/round2/MEModule1",sep='\t',header = T)
MENetCoreP2  <- left_join(MENetCoreP, DEME, by = "ID")
MENetCoreP2$Module1 <- MENetCoreP2$ID %in% ME_Mod1$ID
write.csv(MENetCoreP2,file = "MENetCoreNodeInfor.csv")

