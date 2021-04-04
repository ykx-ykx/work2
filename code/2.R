library(igraph)
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

#读取种子基因
seed<-read.delim("G:/colon_cancer/DEG&driver.txt",sep="\t",header=FALSE)
#构建网络
dat<-as.data.frame(NPPI)
g <- init.igraph(dat)
#计算节点的度
num<-degree(g,mode="total") 
#提取节点的名字
node<-get.vertex.attribute(g)
node<-as.data.frame(node)
#将生成画图需要的文件，即节点，度，分组
degree_res<-cbind(node,num)
degree_res["Group"]<-"PPI_node"
colnames(degree_res)<-c("Gene","Degree","Group")
seed<-merge(seed,degree_res, by.y= "Gene", by.x = "V1")
seed$Group<-"Seed_node"
colnames(seed)<-c("Gene","Degree","Group")
degree_res<-rbind(degree_res,seed)
degree_res$Degree<-log10(degree_res$Degree)
degree_res$Group<-factor(degree_res$Group,levels=c("Seed_node","PPI_node"),ordered=TRUE)
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
Distance_res["group"]<-"ALL"
colnames(Distance_res)<-c("Gene","Distance","group")
seed<-read.delim("G:/colon_cancer/DEG&driver.txt",sep="\t",header=FALSE)
seed<-merge(seed,Distance_res, by.y= "Gene", by.x = "V1")
seed$group<-"SeedGene"
colnames(seed)<-c("Gene","Distance","group")
Distance_res<-rbind(Distance_res,seed)
Distance_res$Distance<-as.numeric(Distance_res$Distance)
#画图
#将数据做常用对数变换， 但是相应的坐标轴刻度的数值还标成变换之前的原始值
library(ggsignif)
library(ggplot2)

compaired <- list(c("Seed_node","PPI_node"))
ggplot(degree_res) +
  aes(x = Group, y = Degree, fill = Group) +
  scale_fill_manual(values = c("#EC7063","#839192"))+
  geom_boxplot(outlier.colour = NA,width=0.5)+
  geom_signif(comparisons = compaired,step_increase = 0.1,
              color="black",
              y_position = 3.6,x_position = c(1, 2),
              map_signif_level=function(p)sprintf("p = %.2g", p),
              test = wilcox.test)+###显著性
  theme_bw() +
  ylim(0,3.7)+
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8),
        axis.title = element_text(colour = "black",size = 13),
        axis.text = element_text(colour = "black",size = 12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = NULL,y="Log10(Degree)")
ggsave("G:/图片/2.tiff", dpi=600)
dev.off()
#画图 
colnames(Distance_res)<-c("Gene","Distance","Group")
Distance_res$Group<- factor(Distance_res$Group, levels=c("SeedGene","ALL"), ordered=TRUE)
levels(Distance_res$Group)<-c("Seed_node","PPI_node")
ggplot(Distance_res) +
  aes(x = Group, y = Distance, fill = Group) +
  scale_fill_manual(values = c("#5DADE2","#95A5A6"))+
  geom_boxplot(outlier.colour = NA,width=0.5)+
  geom_signif(comparisons = compaired,step_increase = 0.1,
              color="black",
              y_position = 5,x_position = c(1, 2),
              map_signif_level=function(p)sprintf("p = %.2g", p),
              test = wilcox.test)+###显著性
  theme_bw() +
  ylim(1,5.2)+
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8),
        axis.title = element_text(colour = "black",size = 13),
        axis.text = element_text(colour = "black",size = 12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = NULL,y="Distance")
ggsave("G:/图片/3.tiff", dpi=600)
dev.off()
