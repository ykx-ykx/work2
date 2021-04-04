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

NPPI<-PPI2[-which(PPI2$SYMBOL.x==PPI2$SYMBOL.y),]
####################################
PPI_table<-NPPI
seedgene<-read.table("DEG&driver.txt",sep="\t",header=F)
seedgene1<-as.matrix(seedgene)
specific<-PPI_table

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
setwd("G:/colon_cancer/Result-0.001")
j=1
while(j<=50){
  n<-which(specific[,1:2]==seedgene1[j])
  if(length(n)>0)
  {
    SeedGene<-seedgene1[j]
    RWR_PPI_Results<-Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,PPI_MultiplexObject,SeedGene)
    result<-RWR_PPI_Results[[1]][which(RWR_PPI_Results[[1]][2]>0.001),] ##阈值可调
    Result<-result
    #提取模块的基因集
    #resulti为包含种子节点的基因，用来计算相关性
    assign(paste("result",j, sep=""),c(result$NodeNames,seedgene1[j]))
    #Result为不包含种子节点的基因集，用来计算模块与模块之间的关系
    assign(paste("Result",j, sep=""),c(Result$NodeNames))
    write.table(result,paste(seedgene1[j], ".txt", sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
    sep=t(c("Gene",SeedGene,dim(result)[1]))
    write.table(sep,"RWR_result_0.001.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
  }
  j=j+1
}
###################
#计算模块之间的关系，相离，相交，包含，删除模块
Disjoint<-NULL
Inter<-NULL
Contain<-NULL
for(i in 1:49){
  for (j in (i+1):50) {
    k<-length(intersect(get(paste("Result",i,sep = "")),get(paste("Result",j,sep = ""))))
    gene1<-length(get(paste("Result",i,sep = "")))
    gene2<-length(get(paste("Result",j,sep = "")))
    MIN<-min(gene1,gene2)
    if(k==0){
      Disjoint<-rbind(Disjoint,c(seedgene1[i],gene1,seedgene1[j],gene2,k))
    }else if(k<MIN){
      Inter<-rbind(Inter,c(seedgene1[i],gene1,seedgene1[j],gene2,k))
    }else{
      Contain<-rbind(Contain,c(seedgene1[i],gene1,seedgene1[j],gene2,k))
    }
  }
}
setwd("G:/colon_cancer/Result")
write.table(Disjoint,"PDisjoint.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(Inter,"PInter.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(Contain,"PContain.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
#############
#计算模块相似性
#计算超几何
#PPI网络中不重复的基因作为背景
GENE<-c(t(NPPI$SYMBOL.x),t(NPPI$SYMBOL.y))
N<-length(unique(GENE))
i<-1
j<-2
for (i in 1:50) {
  a<-seedgene1[i]
  n<-length(get(paste("result",i,sep = "")))
  for (j in 1:50) {
    b<-seedgene1[j]
    m<-length(get(paste("result",j,sep = "")))
    k<-length(intersect(get(paste("result",i,sep = "")),get(paste("result",j,sep = ""))))
    p<-1-phyper(k-1,m,N-m,n)
    p.adjust<-p.adjust(p,method = "BH")
    sim<-length(intersect(get(paste("result",i,sep = "")),get(paste("result",j,sep = ""))))/
      length(union(get(paste("result",i,sep = "")),get(paste("result",j,sep = ""))))
    if(p.adjust<0.001){
      result=t(c(seedgene1[i],seedgene1[j],sim,p,p.adjust))
      write.table(result,"result_0.001.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
      sep=t(c(seedgene1[i],seedgene1[j],sim,p,p.adjust))
    }else{
      sep=t(c(seedgene1[i],seedgene1[j],0,p,p.adjust))
    }
    write.table(sep,"Sim_result_0.001.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
    }
}
