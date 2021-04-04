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
#实际基因集
seedgene<-read.table("DEG&driver.txt",sep="\t",header=F)
seedgene1<-as.matrix(seedgene)
##############
#生成不包含driver与差异基因的基因集合
driver<-read.csv("cancer_gene_census.csv")
driver<-driver$Gene.Symbol
DEG2<-read.table("DEG2.txt",sep="\t",header=TRUE)
DEG2<-DEG2$Gene.Symbol
Gene<-unique(c(NPPI$SYMBOL.x,NPPI$SYMBOL.y))
seedgene<-setdiff(Gene, driver)
seedgene<-setdiff(seedgene,DEG2)
######################################################

init.igraph<-function(dat,dir=FALSE,rem.multi=TRUE)
{
  labels<-union(unique(dat[,1]),unique(dat[,2]))
  ids<-1:length(labels);names(ids)<-labels
  from<-as.character(dat[,1]);to<-as.character(dat[,2])
  edges<-matrix(c(ids[from],ids[to]),nc=2)
  g<-graph.empty(directed = dir)
  g<-add.vertices(g,length(labels))
  V(g)$label=labels
  g<-add.edges(g,t(edges))
  if (rem.multi)
  {
    E(g)$weight<-count.multiple(g)
    g<-simplify(g,remove.multiple = TRUE,
                remove.loops = TRUE,edge.attr.comb = "mean")
  }
  g
}

#构建网络
dat<-as.data.frame(NPPI)
g <- init.igraph(dat)
#计算节点的度
num<-degree(g,mode="total") 
#提取节点的名字
node<-get.vertex.attribute(g)
node<-as.data.frame(node)
#将生成画图需要的文件，即节点，度
degree_res<-cbind(node,num)
###########################################
#计算最短路径
dis<-shortest.paths(g)
path<-dis
Distance_res<-NULL
for (i in 1:15899) {
  path[i,which(path[i,]==Inf)]<-NA
  r<-mean(path[i,],na.rm=TRUE)
  Distance_res<-rbind(Distance_res,c(node$label[i],r))
}
Distance_res<-as.data.frame(Distance_res)
Distance_res$V2<-as.numeric(Distance_res$V2)
###################################################
#计算随机1000次的平均最短路径与平均度
set.seed(123)
res_deg<-NULL
res_dis<-NULL
for (i in 1:1000) {
  #随机从集合中抽取50个基因做随机游走
  seedgene1<-as.data.frame(sample(seedgene,50))
  colnames(seedgene1)<-"gene"
  seed1<-merge(seedgene1,degree_res, by.x= "gene", by.y = "label")
  seed2<-merge(seedgene1,Distance_res, by.x= "gene", by.y = "V1")
  res_deg<-rbind(res_deg,c("Degree",mean(seed1$num)))
  res_dis<-rbind(res_dis,c("Distance",mean(seed2$V2)))
}
res_deg<-as.data.frame(res_deg)
res_dis<-as.data.frame(res_dis)
res_deg$V2<-as.numeric(res_deg$V2)
res_dis$V2<-as.numeric(res_dis$V2)
random<-data.frame(dregee=res_deg$V2,distance=res_dis$V2)

#读取实际种子基因
seed<-read.delim("G:/colon_cancer/DEG&driver.txt",sep="\t",header=FALSE)
seed<-merge(seed,degree_res, by.x= "V1", by.y = "label")
seed<-merge(seed,Distance_res, by.x= "V1", by.y = "V1")
colnames(seed)<-c("mir","degree","distance")
######################################
real<-data.frame(dregee=rep(mean(seed$degree),100),distance=rep(mean(seed$distance),100))

data<-data.frame(dregee=NA,distance =NA)
data[1,]<-colMeans(random)
data[2,]<-colMeans(real)
data["Group"]<-c("Random","Real")
data<-melt(data)
colnames(data)[3]<-"Mean"
data["SD"]<-c(sd(random[,1]),0,sd(random[,2]),0)
p<-data.frame(dregee=NA,distance =NA)
for (i in 1:2) {
  p[1,i]<-wilcox.test(random[,i],real[,i])$p.value
}

data<-data[order(data$Group,decreasing=T),]
data$Group<-factor(data$Group,levels=c("Real","Random"))
data1<-data[c(1,3),]
data2<-data[c(2,4),]
ggplot(data = data1, aes(x=Group, y = Mean, fill = Group)) + 
  geom_bar(position ="dodge",stat="identity", width = 0.4)+
  scale_fill_manual(values = c("#2E86C1","#A6ACAF"))+
  geom_errorbar(aes(x=Group, ymin=Mean-SD, ymax=Mean+SD),  # 添加误差线
                width=0.1, color='black', position = position_dodge(0.6),  
                size=0.6)+
  theme_bw() +
  scale_y_continuous(expand = c(0, 0),limits = c(0,75))+
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8),
        axis.title = element_text(colour = "black",size = 13),
        axis.text = element_text(colour = "black",size = 12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  labs(title ="P (dregee) < 5.6e-39" ,x=NULL,y=NULL)

ggsave("G:/图片/补度随机.tiff", dpi=600)

ggplot(data = data2, aes(x=Group, y = Mean, fill = Group)) + 
  geom_bar(position ="dodge",stat="identity", width = 0.4)+
  scale_fill_manual(values = c("#FF7F50","#A6ACAF"))+
  geom_errorbar(aes(x=Group, ymin=Mean-SD, ymax=Mean+SD),  # 添加误差线
                width=0.1, color='black', position = position_dodge(0.6),  
                size=0.6)+
  theme_bw() +
  scale_y_continuous(expand = c(0, 0),limits = c(0,3))+
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8),
        axis.title = element_text(colour = "black",size = 13),
        axis.text = element_text(colour = "black",size = 12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  labs(title ="P (distance) < 2.9e-33" ,x=NULL,y=NULL)

ggsave("G:/图片/补最短路径随机.tiff", dpi=600)  
  
  ################
  #########
#实际基因集
seedgene<-read.table("DEG&driver.txt",sep="\t",header=F)
seedgene1<-as.matrix(seedgene)
##############
#生成不包含driver与差异基因的基因集合
driver<-read.csv("cancer_gene_census.csv")
driver<-driver$Gene.Symbol
DEG2<-read.table("DEG2.txt",sep="\t",header=TRUE)
DEG2<-DEG2$Gene.Symbol
Gene<-unique(c(NPPI$SYMBOL.x,NPPI$SYMBOL.y))
seedgene<-setdiff(Gene, driver)
seedgene<-setdiff(seedgene,DEG2)

set.seed(123)
Res<-NULL
for (num in 1:1000) {
#随机从集合中抽取50个基因做随机游走
seedgene1<-sample(seedgene,50)
#########
#挖模块

j=1
Result<-list()
while(j<=50){
  n<-which(specific[,1:2]==seedgene1[j])
  if(length(n)>0)
  {
    SeedGene<-seedgene1[j]
    RWR_PPI_Results<-Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,PPI_MultiplexObject,SeedGene)
    result<-RWR_PPI_Results[[1]][which(RWR_PPI_Results[[1]][2]>0.001),] ##阈值可调
    #提取模块的基因集
    Result<-c(Result,list(c(SeedGene,t(result$NodeNames))))
  }
  j=j+1
}
###################
#计算模块之间的关系，相离，相交，包含，删除模块
Disjoint<-0
Inter<-0
Contain<-0
for(i in 1:49){
  for (j in (i+1):50) {
    k<-length(intersect(Result[[i]],Result[[j]]))
    gene1<-length(Result[[i]])
    gene2<-length(Result[[j]])
    MIN<-min(gene1,gene2)
    if(k==0){
      Disjoint<-Disjoint+1
    }else if(k<MIN){
      Inter<-Inter+1
    }else{
      Contain<-Contain+1
    }
  }
}
Res<-rbind(Res,c(Disjoint,Inter,Contain))
}
setwd("G:/")
colnames(Res)<-c("相离","相交","包含")
write.table(Res,"随机.txt",quote=FALSE,sep="\t",row.names=FALSE)
#############
library(reshape2)
library(ggplot2)
library(ggsignif)
random<-read.table("G:/随机.txt",header = T)
colnames(random)<-c("Disjoint","Intersect","Contain")
random<-random/1225
real<-data.frame(rep(238,1000),rep(978,1000),rep(9,1000))
colnames(real)<-c("Disjoint","Intersect","Contain")
real<-real/1225
data<-data.frame(Disjoint=NA,Intersect =NA,Contain=NA)
data[1,]<-colMeans(random)
data[2,]<-colMeans(real)
data["Group"]<-c("random","real")
data<-melt(data)
colnames(data)[3]<-"Mean"
data["SD"]<-c(sd(random[,1]),0,sd(random[,2]),0,sd(random[,3]),0)
p<-data.frame(Disjoint=NA,Intersect =NA,Contain=NA)
for (i in 1:3) {
p[1,i]<-wilcox.test(random[,i],real[,i])$p.value
}
data<-data[order(data$Group,decreasing=T),]
data$Group<-factor(data$Group,levels=c("real","random"))
ggplot(data = data, aes(x=variable, y = Mean, fill = Group,color=Group)) + 
  geom_bar(position =position_dodge(0.6),stat="identity", width = 0.4,alpha=0.3,size=0.8)+
  scale_fill_manual(values = c("#EC7063","#2874A6"))+
  scale_color_manual(values = c("#EC7063","#2874A6"))+
  geom_errorbar(aes(x=variable, ymin=Mean-SD, ymax=Mean+SD),  # 添加误差线
                width=0.1, color='black', position = position_dodge(0.6),  
                size=0.6)+
  geom_signif(y_position=c(0.85), xmin=c(1.85), xmax=c(2.14), color="black" ,
              annotation=c("P<e-10"), tip_length=c(0.05,0.02), size=0.6,
              textsize = 3.8, vjust = -0.5)  +
  geom_signif(y_position=c(0.38), xmin=c(0.85), xmax=c(1.14), color="black" ,
              annotation=c("P<e-14"), tip_length=c(0.2,0.03), size=0.6,
              textsize = 3.8, vjust = -0.5)  +
  geom_signif(y_position=c(0.05), xmin=c(2.85), xmax=c(3.14), color="black" ,
              annotation=c("P<e-41"), tip_length=c(0.02,0.03), size=0.6, 
              textsize = 3.8, vjust = -0.5)  +
theme_bw() +
  scale_y_continuous(expand = c(0, 0),limits = c(0,1))+
  theme(axis.title = element_text(colour = "black",size = 14),
        axis.text = element_text(colour = "black",size = 14),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x=NULL,y="Fraction")+
  theme(legend.position = c(0.15, 0.9),
  legend.direction = "vertical",
  legend.key.width=unit(.6,"inches"),
  legend.background=element_rect(colour="#566573",size=0.4),
  legend.key.height=unit(.2,"inches"),
  legend.text=element_text(colour="black",size=13),
  legend.title=element_blank()) 
ggsave("G:/图片/补1.tiff", dpi=600)  

