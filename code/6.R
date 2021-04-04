#BiocManager::install("impute")
library(impute)
library(WGCNA)
library(ggplot2)
library(WeightedCluster)
F<-read.delim("G:/colon_cancer/Result/Sim_result_0.001.txt",header = FALSE)
G<-matrix(data=F$V3,nrow=46,ncol=46,byrow=T)
#导入马尔可夫聚类的R包
library(MCL)

#尝试不同的inflation计算平均轮廓系数并画图
#导入WeightedCluster包计算轮廓系数
### 把邻接矩阵转换为拓扑重叠矩阵，以降低噪音和假相关，获得距离矩阵。
TOM = TOMsimilarity(G)
dissTOM = 1-TOM

i<-1
Result<-NULL
for(i in seq(1,1.9, by=0.05)){
  
  res<-mcl(G,inflation =i,addLoops=T,max.iter = 1000)
  
  qual <- wcClusterQuality(dissTOM, res$Cluster)
  
  Result<-rbind(Result,c(i,qual$stats[5],qual$stats[4]))
}
#画图
colnames(Result)<-c("inflation","ASWw","ASW")
Result<-as.data.frame(Result)
ggplot(Result, aes(x=inflation, y=ASWw)) + 
  geom_line() + 
  geom_point(size=1)+
  annotate(geom = "point",
           x = 1.3, y = 0.22861986,color="red",size=3)+
  annotate(geom = "line",
           x =1.3, y = c(0,0.22),linetype="dashed")+
  scale_x_continuous(breaks=c(1, 1.1, 1.2, 1.3, 1.4, 1.5,1.6,1.7,1.8))+
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x = "inflation factor",y="Average Silhouette width (weighted)")
ggsave("G:/图片/4.tiff", dpi=600)
dev.off()
##############
txt1<-1:10
txt1<-as.data.frame(txt1)
txt1["sum"]<-c(8,2,13,5,2,2,2,6,3,3)
txt1$txt1<-as.factor(txt1$txt1)
txt1$txt1<- factor(txt1$txt1, levels=c("3","1","8","4","9","10","2","5","6","7"), ordered=TRUE)
ggplot(txt1)+
  aes(x=txt1,fill=txt1,y=sum)+
  geom_bar(stat = "identity",width = 0.7)+
  scale_fill_manual(values = c("#D7BDE2","#F5B7B1","#F8C471","#A9CCE3","#F5CBA7","#B2BABB","#BB8FCE","#A3E4D7","#A9DFBF","#F9E79F"))+
  theme(panel.background = element_rect(fill = NA),
        legend.position = "none",
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8),
        axis.title = element_text(colour = "black",size = 14),
        axis.text = element_text(colour = "black",size = 14),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x = "Cluster",y="Number of class")+
  annotate(geom = "text",x = 1, y = 13.5, hjust = 0.5,size=4,label="13")+
  annotate(geom = "text",x = 2, y = 8.5, hjust = 0.5,size=4,label="8")+
  annotate(geom = "text",x = 3, y = 6.5, hjust = 0.5,size=4,label="6")+
  annotate(geom = "text",x = 4, y = 5.5, hjust = 0.5,size=4,label="5")+
  annotate(geom = "text",x = 5, y = 3.5, hjust = 0.5,size=4,label="3")+
  annotate(geom = "text",x = 6, y = 3.5, hjust = 0.5,size=4,label="3")+
  annotate(geom = "text",x = 7, y = 2.5, hjust = 0.5,size=4,label="2")+
  annotate(geom = "text",x = 8, y = 2.5, hjust = 0.5,size=4,label="2")+
  annotate(geom = "text",x = 9, y = 2.5, hjust = 0.5,size=4,label="2")+
  annotate(geom = "text",x = 10,y = 2.5, hjust = 0.5,size=4,label="2")+
  annotate(
    geom = "rect",
    ymin = -Inf, ymax = Inf, 
    xmin = c(0, 4.5),
    xmax = c(4.5, 11),
    fill = c("#F5B7B1", "#AED6F1"), alpha = 0.3)
ggsave("G:/图片/5.tiff", dpi=600)
dev.off()
