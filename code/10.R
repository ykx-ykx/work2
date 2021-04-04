library(corrplot)
library(clusterProfiler)
library(ggplot2)
library(stringr)
setwd("G:/colon")
miRNA_target<-read.delim("G:/colon_cancer/c3.mir.mirdb.v7.2.entrez.gmt",header = FALSE)
target<-miRNA_target[c(274,524,1104,760),]
mir_29a<-na.omit(unique(t(cbind(target[1,2:1071],target[2,2:1071]))))
mir_891a<-na.omit(t(target[3,2:1071]))
mir_548v<-na.omit(t(target[4,2:1071]))
mir_res<-unique(c(t(mir_29a),t(mir_891a),t(mir_548v)))

#富集
#GO富集分析
#这里只对BP进行富集
#对其他富集可以改ont = "BP"参数。
ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = mir_res, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析

#KEGG分析
ekk <- enrichKEGG(gene= mir_res,organism  = 'hsa',qvalueCutoff = 0.05)	 #KEGG富集分析
barplot(ekk,showCategory=10,drop=T)
barplot(ego,showCategory=10,drop=T)
#写出富集结果
write.csv(summary(ekk),"G:/KEGG.csv",row.names =FALSE)
write.csv(summary(ego),"G:/GO.csv",row.names =FALSE)

kegg <- read.csv("G:/KEGG.csv",stringsAsFactors=FALSE)
go<-read.csv("G:/GO.csv",stringsAsFactors=FALSE)
#对富集结果按照p.adjust进行从小到大排序，保证最显著的通路在前
kegg <- kegg[order(kegg$p.adjust),]

#这里画图只展示top10的通路
kegg <- kegg[1:11,]
kegg<-kegg[-5,]
#提取每条通路里面差异表达的基因数
top10 <- data.frame(kegg$Description,kegg$Count ,kegg$p.adjust)
top10<-top10[order(top10$kegg.p.adjust,decreasing=T),]
colnames(top10) <- c("Description","count","p.adj")
top10$Description<-factor(top10$Description,levels =top10$Description)

#fill=padj fill颜色填充，使用连续值padj
ggplot(data=top10,aes(x=Description,y=count,fill=p.adj))+
  geom_bar(stat="identity") + coord_flip()+
  theme(panel.background=element_rect(fill='transparent',color='black'),
        axis.text.y=element_text(color="black",size=13.5),
        axis.text.x =element_text(color="black",size=13.5),
        plot.title = element_text(hjust = 0.5))+ 
  scale_fill_gradient(low="red",high="blue")+
  #scale_x_discrete(limits=rev(top10[,1])) +
  labs(x="",y="",title="KEGG")+
  scale_x_discrete(labels=function(x) str_wrap(x, width=30))
ggsave("G:/图片/12.tiff", dpi=600)




go <- go[order(go$p.adjust),]

#这里画图只展示top10的通路
go <- go[1:10,]

#提取每条通路里面差异表达的基因数
top10 <- data.frame(go$Description,go$Count ,go$p.adjust)
top10<-top10[order(top10$go.p.adjust,decreasing=T),]
colnames(top10) <- c("Description","count","p.adj")
top10$Description<-factor(top10$Description,levels =top10$Description)

#fill=padj fill颜色填充，使用连续值padj
ggplot(data=top10,aes(x=Description,y=count,fill=p.adj))+
  geom_bar(stat="identity") + coord_flip()+
  theme(panel.background=element_rect(fill='transparent',color='black'),
        axis.text.y=element_text(color="black",size=13.5),
        axis.text.x =element_text(color="black",size=13.5),
        plot.title = element_text(hjust = 0.5))+ 
  scale_fill_gradient(low="red",high="blue")+
  #scale_x_discrete(limits=rev(top10[,1])) +
  labs(x="",y="",title="GO")+
  scale_x_discrete(labels=function(x) str_wrap(x, width=30))
ggsave("G:/图片/13.tiff", dpi=600)
