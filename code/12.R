library(NbClust)
library(ggplot2)
setwd("G:/colon")
#处理数据，对于dead样本，overall survival采用day_to_death
#对于alive样本，overall survival采用day_to_last_follow_up
#最后只保留样本和生存时间两列
phen<-read.delim("TCGA-COAD.GDC_phenotype.tsv")
phen<-phen[,c(1,93,39,75)]
#删除无用数据
phen<-phen[which(phen$pathologic_M!=""),]
phen<-phen[which(phen$pathologic_M!="MX"),]
#只保留3,4期数据
d1<-phen[which(phen$tumor_stage.diagnoses=="stage iv"),]
d2<-phen[which(phen$tumor_stage.diagnoses=="stage iva"),]
d3<-phen[which(phen$tumor_stage.diagnoses=="stage ivb"),]

d4<-phen[which(phen$tumor_stage.diagnoses=="stage iii"),]
d5<-phen[which(phen$tumor_stage.diagnoses=="stage iiia"),]
d6<-phen[which(phen$tumor_stage.diagnoses=="stage iiib"),]
d7<-phen[which(phen$tumor_stage.diagnoses=="stage iiic"),]


phen<-rbind(d1,d2,d3,d4,d5,d6,d7)
phen["Group"]<-NA
phen[which(phen$tumor_stage.diagnoses=="stage iv"),5]<-"Distance Metastases"
phen[which(phen$tumor_stage.diagnoses=="stage iva"),5]<-"Distance Metastases"
phen[which(phen$tumor_stage.diagnoses=="stage ivb"),5]<-"Distance Metastases"

phen[which(phen$tumor_stage.diagnoses=="stage iii"),5]<-"Lymph Node Metastasis"
phen[which(phen$tumor_stage.diagnoses=="stage iiia"),5]<-"Lymph Node Metastasis"
phen[which(phen$tumor_stage.diagnoses=="stage iiib"),5]<-"Lymph Node Metastasis"
phen[which(phen$tumor_stage.diagnoses=="stage iiic"),5]<-"Lymph Node Metastasis"


phen<-na.omit(phen)
phen$submitter_id.samples<-gsub(pattern="-", replacement=".", phen$submitter_id.samples)
#只保留3列显著的miRNA
miRNA<-read.delim("batch_mir.txt",header = T)
miRNA<-miRNA[c(567,209,756),]
miRNA<-as.data.frame(t(miRNA))
colnames(miRNA)<-miRNA[1,]
miRNA<-miRNA[-1,]
miRNA['Name']<-rownames(miRNA)
rownames(miRNA)<-1:428
miR_phen<-merge(miRNA,phen,by.x = "Name",by.y = "submitter_id.samples")
miR_phen<-miR_phen[,-c(1,5,6)]
miR_phen1<-as.data.frame(apply(miR_phen[,1:3], 2, as.numeric))
miR_phen[,1:3]<-miR_phen1

#层次聚类
df<-scale(miR_phen[,1:3])
dist<-dist(df,method = "euclidean")
hc<-hclust(dist,method = "ward.D2")
plot(hc,hang = -1,labels=F,
     xlab="3-miRNA expression(n=900)",ylab = "Similarity distance score")
nc<-NbClust(df,distance = "euclidean",min.nc = 2,max.nc = 15,method = "average")
par(mfrow=c(1,1))

tiff(filename = "G:/图片/15.tiff",width = 3400,height = 3000,res=600);  #保存图片
barplot(table(nc$Best.n[1,]), xlab="Number of Clusters", ylab="Number of Supporting", 
        main="Determine the Number of Clustering")
dev.off();
tiff(filename = "G:/图片/16.tiff",width = 4000,height = 3000,res=600) 

plot(hc,hang = -1,labels=F,mgp=c(1.2,0.5,-0.5),
     xlab="3-miRNA expression",ylab = "Similarity distance score")
rect.hclust(hc, k=2,border = "red")

dev.off()

clusters<-cutree(hc,k=2)
mir<-miR_phen
mir["cul"]<-clusters
colnames(mir)<-c("mir-548v","mir-29a","mir-891a","Status","Transfer","Group")
sub1<-mir[which(mir$Group==1),]
sub2<-mir[which(mir$Group==2),]
sub3<-mir[which(mir$Group==3),]
mir$Group<-as.factor(mir$Group)
levels(mir$Group)<-c("S1","S2","S3")
percent<-c("0%","25%","50%","75%","100%")
ggplot(data = mir, mapping = aes(x = Group, fill = Status)) +
  geom_bar(position = 'fill',width = 0.3)+
  scale_fill_manual(values = c("#85C1E9","#99A3A4"))+
  scale_y_continuous(expand = c(0,0),labels=percent)+
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8),
        axis.title = element_text(colour = "black",size = 14),
        axis.text = element_text(colour = "black",size = 14),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = rel(1.2)), 
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = "bottom",legend.direction = "horizontal")+
  labs(title="Survival state",y="Frequency")+
  guides(fill=guide_legend(title=NULL))+
  annotate(geom = "text",x = 2, y = 0.15,label="31%",size=6)+
  annotate(geom = "text",x = 2, y = 0.68,label="69%",size=6)+
  annotate(geom = "text",x = 1, y = 0.28,label="52%",size=6)+
  annotate(geom = "text",x = 1, y = 0.8,label="48%",size=6)

#ggsave("G:/图片/17.tiff", dpi=600)

mir<-mir[order(mir$Transfer,decreasing=T),]
mir$Transfer<-factor(mir$Transfer,levels =c("Lymph Node Metastasis","Distant Metastasis"))
p21<-ggplot(data = mir, mapping = aes(x = Group, fill = Transfer)) +
  geom_bar(position = 'fill',width = 0.3)+
  scale_fill_manual(values = c("#85C1E9","#99A3A4"))+
  scale_y_continuous(expand = c(0,0),labels=percent)+
  #scale_x_discrete(expand=c(0.6, 0))+
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8),
        axis.title = element_text(colour = "black",size = 14),
        axis.text = element_text(colour = "black",size = 14),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = rel(1)), 
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  labs(title="Metastatic state",y="Frequency")+
  guides(fill=guide_legend(title=NULL))+
  annotate(geom = "text",x = 2, y = 0.2,label="36%",size=6)+
  annotate(geom = "text",x = 2, y = 0.7,label="64%",size=6)+
  annotate(geom = "text",x = 1, y = 0.3,label="56%",size=6)+
  annotate(geom = "text",x = 1, y = 0.78,label="44%",size=6)

#ggsave("G:/图片/18.tiff", dpi=600)
tiff(filename = "G:/改/20.tiff",width = 6000,
     height = 3000,res=600);  #保存图片
ggarrange(p20,p21,nrow = 1, ncol = 2)
dev.off()