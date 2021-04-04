library(readxl)
#install.packages("BiocManager")
library("BiocManager")
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("RandomWalkRestartMH")
library(RandomWalkRestartMH)
library(igraph)
setwd("G:/colon_cancer")

PPI<-read_excel("41467_2019_9186_MOESM3_ESM.xlsx")
# transform id  
map1<- bitr(PPI$Protein_A_Entrez_ID, fromType = "ENTREZID",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
map2<- bitr(PPI$Protein_B_Entrez_ID, fromType = "ENTREZID",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
dt_merge <- merge(map1,PPI, by.y = "Protein_A_Entrez_ID", by.x = "ENTREZID")
PPI1<-merge(map2,dt_merge, by.y = "Protein_B_Entrez_ID", by.x = "ENTREZID")
PPI2<-PPI1[,c(2,4)]
#去除自己与自己相连的节点
NPPI<-PPI2[-which(PPI2$SYMBOL.x==PPI2$SYMBOL.y),]
####################################
PPI_table<-NPPI
seedgene<-read.table("FDEG&driver.txt",sep="\t",header=FALSE)
seedgene1<-as.matrix(seedgene)
specific<-PPI_table
#创建图
PPI_Network <- graph.data.frame(specific,directed=FALSE)
#删除多循环的边
PPI_Network <- igraph::simplify(PPI_Network, remove.multiple = TRUE, remove.loops = TRUE)
#创建复用网络
PPI_MultiplexObject <- create.multiplex(PPI_Network,Layers_Name=c("PPI"))
#计算复用网络的邻接矩阵
AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
#邻接矩阵规范化
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)

#########
#挖模块

j=1
Gene_Result<-matrix(data=NA,nrow = 46,ncol=195)
while(j<=46){
  n<-which(specific[,1:2]==seedgene1[j])
  if(length(n)>0)
  {
    SeedGene<-seedgene1[j]
    RWR_PPI_Results<-Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,PPI_MultiplexObject,SeedGene)
    result<-RWR_PPI_Results[[1]][which(RWR_PPI_Results[[1]][2]>0.001),] ##阈值可调
    l<-length(result$NodeNames)+1
    Gene_Result[j,1:l]<-c(seedgene1[j],result$NodeNames)
  }
  j<-j+1
}
###################
#计算10个簇中的基因
library(MCL)
F<-read.delim("G:/colon_cancer/Result/Sim_result_0.001.txt",header = FALSE)
G<-matrix(data=F$V3,nrow=46,ncol=46,byrow=T)
res<-mcl(G,inflation =1.3,addLoops=T,max.iter = 1000)
clu<-res$Cluster
clu<-as.factor(clu)
levels(clu)<-1:10
RR<-Gene_Result
RR<-as.data.frame(RR)
RR["group"]<-clu
Clu_Result<-matrix(data=NA,nrow = 10,ncol = 900)

i<-1
for (i in 1:46) {
  
  str<-min(which(is.na(Clu_Result[RR$group[i],])))
  end<-str+length(Gene_Result[i,which(!is.na(Gene_Result[i,]))])-1
  Clu_Result[RR$group[i],str:end]<-Gene_Result[i,which(!is.na(Gene_Result[i,]))]
  
}

##################
#10个簇的聚类个数
Clu_num<-rbind(8,2,13,5,2,2,2,6,3,3)
#富集分析Hallmark与10个簇
hallmark<-read.delim("h.all.v7.2.symbols.gmt",,sep="\t",header=FALSE)
hallmark<-as.matrix(hallmark)
Ver_Gene<-as.vector(Clu_Result)
Ver_Hallmark<-as.vector(hallmark[,-1])
N<-length(unique(c(Ver_Gene,Ver_Hallmark)))-2
Hall_Gene_Enrich<-NULL
i<-1
j<-1
for (i in 1:dim(Clu_Result)[1]) {
  gene<-unique(Clu_Result[i,which(!is.na(Clu_Result[i,]))])
  n<-length(gene)
  for (j in 1:dim(hallmark)[1]) {
    hall<-hallmark[j,2:max(which(hallmark[j,]!=""))]
    m<-length(hall)
    k<-length(intersect(gene,hall))
    p<-1-phyper(k-1,m,N-m,n)
    p.adjust<-p.adjust(p,method = "BH")
    if(p.adjust<=0.05){
      Hall_Gene_Enrich<-rbind(Hall_Gene_Enrich,c(i,hallmark[j,1],p,p.adjust,Clu_num[i]))
    }
    
  }
  
}

#富集分析miRNA与10个簇
miRNA<-read.delim("c3.mir.mirdb.v7.2.symbols.gmt",sep="\t",header=FALSE)
miRNA<-as.matrix(miRNA)
Ver_Gene<-as.vector(Clu_Result)
Ver_miRNA<-as.vector(miRNA[,-1])
N<-length(unique(c(Ver_Gene,Ver_miRNA)))-2
miRNA_Gene_Enrich<-NULL
i<-1
j<-1
for (i in 1:dim(Clu_Result)[1]) {
  gene<-unique(Clu_Result[i,which(!is.na(Clu_Result[i,]))])
  n<-length(gene)
  for (j in 1:dim(miRNA)[1]) {
    RNA<-miRNA[j,2:max(which(miRNA[j,]!=""))]
    m<-length(RNA)
    k<-length(intersect(gene,RNA))
    p<-1-phyper(k-1,m,N-m,n)
    p.adjust<-p.adjust(p,method = "BH")
    if(p.adjust<=0.001){
      miRNA_Gene_Enrich<-rbind(miRNA_Gene_Enrich,c(i,miRNA[j,1],p,p.adjust,Clu_num[i]))
    }
    
  }
  
}
setwd("G:/colon_cancer/Result")
write.table(Hall_Gene_Enrich,"Hall_Gene_Enrich_0.05.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(miRNA_Gene_Enrich,"miRNA_Gene_Enrich_0.001.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)

#Hall
library(grid)
library(stringr)
library(ggplot2)
Hall_plot<-as.data.frame(Hall_Gene_Enrich)
colnames(Hall_plot)<-c("Clu","Hall","P","P.adjust","Number")
Hall_plot$P<-as.numeric(Hall_plot$P)
Hall_plot$P.adjust<-as.numeric(Hall_plot$P.adjust)
Hall_plot$Clu<-as.factor(Hall_plot$Clu)
Hall_plot["Y"]<--log10(Hall_plot$P.adjust)
#0-1标准化x=x/sum(x)
Hall_plot$Y<-Hall_plot$Y/sum(Hall_plot$Y)
#将Hall拆分为多列
Hall_plot[,c("DAG","hall")] <- str_split_fixed(Hall_plot$Hall, "_", 2)
#将下划线用空格替换
Hall_plot$hall<-gsub(pattern="_", replacement=" ", Hall_plot$hall)
#设置类的名字
levels(Hall_plot$Clu)<-c("Angiogenesis","Apoptosis","Surviral","Growth",
                         "Immune","Metabolism","Metastasis",
                         "Proliferation","TGF-Beta","MYC")

#按类顺序画图
Hall_plot$Clu<- factor(Hall_plot$Clu, levels=c("Angiogenesis","Surviral","Growth",
                                               "Immune","Metabolism","Metastasis",
                                               "Proliferation","TGF-Beta","MYC","Apoptosis"), ordered=TRUE)
#画图
p<-ggplot(data = Hall_plot, 
          mapping = aes(x = hall, y =Y ,fill=Y))+
  labs(x=NULL, y=NULL) +
  scale_fill_gradient(low="#D7BDE2",high="#0E4E75")+
  geom_bar(stat='identity')+
  coord_flip()+
  facet_grid(~ Clu) +
  guides(fill=FALSE)+
  theme_bw()+
  scale_y_continuous(breaks=c(0,0.05))+
  theme(strip.text.x=element_text(size=rel(0.9)),
        axis.text.x = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black"),
        panel.background = element_rect(size = 1), 
        plot.background = element_rect())+
  labs(y="Enrichment Strength",subtitle = " Cluster1    Cluster2    Cluster3    Cluster4    Cluster5    Cluster6    Cluster7    Cluster8    Cluster9    Cluster10")
#设置分面小框的颜色
g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("#F5B7B1","#BB8FCE","#D7BDE2","#A9CCE3","#A3E4D7","#A9DFBF","#F9E79F","#F8C471","#F5CBA7","#B2BABB")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

tiff(filename = "G:/图片/6.tiff",width = 6050,
     height = 3000,res=600);  #保存图片
grid.draw(g)
dev.off();

#########################
#mirna
mir_plot<-as.data.frame(miRNA_Gene_Enrich)
colnames(mir_plot)<-c("Clu","Mir","P","P.adjust","Number")
#大写转为小写
mir_plot$Mir<-tolower(mir_plot$Mir)
#将mir拆分为多列
mir_plot[,c("mir", "DAG")] <- str_split_fixed(mir_plot$Mir, "_", 2)
mir_plot<-mir_plot[,-7]
mir_plot<-mir_plot[-12,]
#将r用R替换
mir_plot$mir<-gsub(pattern="r", replacement="R",mir_plot$mir)
mir_plot$P<-as.numeric(mir_plot$P)
mir_plot$P.adjust<-as.numeric(mir_plot$P.adjust)
mir_plot$Clu<-as.factor(mir_plot$Clu)
mir_plot["Y"]<--log10(mir_plot$P.adjust)


#设置类的名字
levels(mir_plot$Clu)<-c("miR-3158","miR-4802","miR-142","miR-200",
                        "miR-92","miR-5682","miR-6876",
                        "miR-12119")
#0-1标准化x=(x-min)/(max-min)
mir_plot$Y<-mir_plot$Y/sum(mir_plot$Y)
#按类顺序画图
mir_plot$Clu<- factor(mir_plot$Clu, levels=c("miR-4802","miR-142","miR-200",
                                             "miR-92","miR-5682","miR-6876",
                                             "miR-12119","miR-3158"), ordered=TRUE)



#画图
m<-ggplot(data = mir_plot, 
          mapping = aes(x = mir, y =Y ,fill=Y))+
  labs(x=NULL, y=NULL) +
  scale_fill_gradient(low="#D7BDE2",high="#0E4E75")+
  geom_bar(stat='identity')+
  coord_flip()+
  facet_grid(~ Clu) +
  guides(fill=FALSE)+
  theme_bw()+
  scale_y_continuous(breaks=c(0,0.02))+
  theme(strip.text.x=element_text(size=rel(1.5)),
        axis.text.x = element_text(colour = "black",size = 13), 
        axis.text.y = element_text(colour = "black",size = 13),
        panel.background = element_rect(size = 3), 
        plot.background = element_rect(),
        title=element_text(size=13))+
  labs(y="Enrichment Strength",subtitle = "     Cluster2         Cluster3          Cluster4          Cluster6          Cluster7           Cluster8          Cluster9         Cluster10")
#设置分面小框的颜色
mm <- ggplot_gtable(ggplot_build(m))
strip_both <- which(grepl('strip-', mm$layout$name))
fills <- c("#BB8FCE","#D7BDE2","#A9CCE3","#A9DFBF","#F9E79F","#F8C471","#F5CBA7","#B2BABB")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', mm$grobs[[i]]$grobs[[1]]$childrenOrder))
  mm$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(mm)
tiff(filename = "G:/图片/7.tiff",width = 6200,
     height = 5600,res=600);  #保存图片
grid.draw(mm)
dev.off();

