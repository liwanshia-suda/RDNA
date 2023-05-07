library(stringr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

ourRMSD=function(data){
  RMSDmatrix<-matrix(0,nrow = ncol(data),ncol =ncol(data))
  rownames(RMSDmatrix)<-colnames(data)
  colnames(RMSDmatrix)<-colnames(data)
  for(i in 1:ncol(data)){
    for(j in 1:ncol(data)){
   RMSDmatrix[i,j]<-sqrt(sum((data[,i]-data[,j])^2)/nrow(data)) 
    } 
  }
  return(RMSDmatrix)
}


compareGroup <- function(dat,ourdata,
                         sampleNum=10,
                         #The number of random
                         randomNum=100,
                         #random seed
                         randomSeedFold=100,
                         #RMSD/cor
                         method="RMSD"){
  
  ourdataSampleNum <- sum(str_detect(colnames(ourdata),"sin|mul"))
  
  for(i in 1:randomNum){
    ##set random seed,make it Repeatable
    set.seed(i*randomSeedFold)
    
    if(sampleNum < ncol(dat)-1)
      data <- dat[,c(1,sample(2:ncol(dat),sampleNum,replace = F))]
    else
      data <- dat
    
    data <- inner_join(data,ourdata,by="gene")
    data <- data[,str_detect(colnames(data),"TCGA|_P")]
    data <- edgeR::cpm(data,log = T)
    
    if(method == "RMSD")dataCor <- ourRMSD(data)
    ##Pearson correlation
    if(method == "cor")dataCor <- cor(data)
    
    ggplot2table <- reshape2::melt(dataCor[(ncol(data)-ourdataSampleNum+1):ncol(data),1:(ncol(data)-ourdataSampleNum)])
    ggplot2table$group <- substr(ggplot2table$Var1,1,3)
    ggplot2table$group<-factor(ggplot2table$group,levels = c("sin","mul"))
    
    ggplot2table$ourp <- wilcox.test(ggplot2table[ggplot2table$group=="sin",3],ggplot2table[ggplot2table$group=="mul",3])$p.value
    ggplot2table$randomseed <- i*randomSeedFold
    cat(paste0("computing:",i,"...\r"))
    if(i==1) ResTable <- ggplot2table else ResTable <- rbind(ResTable,ggplot2table)
    if(sampleNum >= ncol(dat)-1) {
      cat("\nAll samples are used to compute!!!")
      return(ResTable)}
  }
  return(ResTable)
}

#=======================================
##Load TCGA data and matching normal data
data <- read.table(file = "TCGA-LUAD.htseq_counts.tsv.gz",header = T,row.names = 1,check.names = F)
data <- round(2^data-1)
pheotype <- read.delim(file = "TCGA-LUAD.GDC_phenotype.tsv.gz",header = T)
samples <- intersect(colnames(data),pheotype$submitter_id.samples)
pheotype <- pheotype[match(samples,pheotype$submitter_id.samples),]
data <-data[,samples]
table(colnames(data)==pheotype$submitter_id.samples)

##Load our data
ourdata <- read.csv(file = "count.csv")
#ourdata<-ourdata[,c("gene","sin_P10_1","sin_P10_2","sin_P10_3","mul_P10_1","mul_P10_2","mul_P10_3")]

dataannotation <- read.table(file = "gencode.v22.annotation.gene.probeMap",header = T)
data <- rownames_to_column(data,var="Ensembl_ID")
data <- inner_join(data,dataannotation,by=c("Ensembl_ID" = "id"))
data <- data[,c("gene",samples)]
pheotype$tumor_stage.diagnoses=gsub("[ab]$","",pheotype$tumor_stage.diagnoses)

##select the tumor stage
selectedColumn <- c("gene",pheotype$submitter_id.samples[pheotype$tumor_stage.diagnoses=="stage iv" &
                                                           substr(colnames(data)[-1],14,15)!="11"])
##select the Normal stage
#selectedColumn <- c("gene",pheotype$submitter_id.samples[substr(colnames(data)[-1],14,15)=="11"])

##select the compared samples
dataTemp <- data[,selectedColumn]

#set sampleNum as Inf and all samples will be used
ggplot2table <- compareGroup(dataTemp,ourdata,sampleNum=Inf,randomNum=100,randomSeedFold=100,method="cor")
#============================
##plot boxplot of sin and mul
ourp <- wilcox.test(ggplot2table[ggplot2table$group=="sin",3],ggplot2table[ggplot2table$group=="mul",3])$p.value
ourp


p2=ggplot(ggplot2table,aes(x=group,y=value,fill=group))+ 
  stat_boxplot(mapping=aes(x=group,y=value),
               geom ="errorbar",                          
               width=0.2,position=position_dodge(0.5))+     
  geom_boxplot(aes(fill=group),                            
               position=position_dodge(0.8),                
               width=0.55,                                   
               outlier.color = "white")+                    
  stat_summary(mapping=aes(group=group),                    
               fun="mean",                                   
               geom="point",shape=23,size=3,fill="white",    
               position=position_dodge(0.8))+               
  scale_fill_manual(values=c('#104E8B',"#CD2626"))+theme_bw()+
  theme(panel.grid = element_blank(),axis.text=element_text(face = "bold",size = 14,colour = "black"), legend.position = 'none')+
  labs(x="",y="")+
  stat_compare_means(aes(group=group),method = "wilcox.test",label="p.signif",comparisons = list(c("sin","mul")),label.x.npc = 0.5,exact = FALSE) +
  scale_x_discrete(labels = c('SE','ME'))
p2
ggsave(p2,filename = paste0("cor.pdf"),width = 5,height = 5)

write.csv(ggplot2table,file = "ggplot2table.csv",quote = F)

 