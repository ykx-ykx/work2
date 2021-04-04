library(ggplot2)
library(NbClust)
library(reshape2)
library(ggsignif)
library(dplyr)
#MIR靶基因
miRNA_target<-read.delim("G:/colon_cancer/c3.mir.mirdb.v7.2.symbols.gmt",header = FALSE)
target<-miRNA_target[c(274,524,1104,760),]
mir_29a<-na.omit(unique(t(cbind(target[1,2:1071],target[2,2:1071]))))
mir_29a<-as.data.frame(mir_29a[-435,])
mir_891a<-as.data.frame(t(target[3,2:128]))
mir_548v<-as.data.frame(t(target[4,2:138]))
colnames(mir_29a)<-"mir_29a"


#读取基因分类
Group<-read.delim("G:/Newcon/cancer_gene_census.txt")
#读取免疫基因
inna<-read.delim("G:/Newcon/innate.txt")
inna<-toupper(t(inna))
mir_29a1<-merge(mir_29a,Group,by.x = "mir_29a",by.y="Gene.Symbol")
mir_891a1<-merge(mir_891a,Group,by.x = "1104",by.y="Gene.Symbol")
mir_548v1<-merge(mir_548v,Group,by.x = "760",by.y="Gene.Symbol")

mir_29ai<-intersect(t(mir_29a),inna)
mir_891ai<-intersect(t(mir_891a),inna)
mir_548vi<-intersect(t(mir_548v),inna)

mir_29a1[56:98,1]<-mir_29ai
mir_29a1[56:98,2]<-"Immune Gene"
mir_548v1[7:11,1]<-mir_548vi
mir_548v1[7:11,2]<-"Immune Gene"
mir_891a1[11:17,1]<-mir_891ai
mir_891a1[11:17,2]<-"Immune Gene"

mir_29a1[,1]<-"mir-29a"
mir_548v1[,1]<-"mir-548v"
mir_891a1[,1]<-"mir-891a"
colnames(mir_29a1)<-c("Mir","Group")
colnames(mir_548v1)<-c("Mir","Group")
colnames(mir_891a1)<-c("Mir","Group")
res<-rbind(mir_29a1,mir_548v1,mir_891a1)
res[which(res$Group=="Immune Gene"),2]<-1
res[which(res$Group=="oncogene"),2]<-2
res[which(res$Group=="TSG"),2]<-3
res[which(res$Group=="fusion"),2]<-4
res$Mir<- factor(res$Mir, levels=c("mir-29a", "mir-891a", "mir-548v"), ordered=TRUE)
percent<-c("0%","25%","50%","75%","100%")
ggplot(data = res, mapping = aes(x = Mir, fill = Group)) +
  geom_bar(position = 'fill',width = 0.3)+
  scale_fill_manual(values = c("#85C1E9","#F1948A","#7DCEA0","#DCDCDC"),
                    breaks=c("1", "2", "3","4"),
                    labels=c("Immune Gene", "Oncogene", "TSG","Other Gene"))+
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
  labs(x="",y="Frequency of genes")+
  guides(fill=guide_legend(title=NULL))+
  scale_x_discrete(breaks=c("mir-29a", "mir-891a", "mir-548v"),
                   labels=c("miR-29a", "miR-891a", "miR-548v"))+
  annotate(geom = "text",x = 1, y = 0.8,label="44%",size=6)+
  annotate(geom = "text",x = 1, y = 0.51,label="10%",size=6)+
  annotate(geom = "text",x = 1, y = 0.39,label="15%",size=6)+
  annotate(geom = "text",x = 1, y = 0.16,label="31%",size=6)+
  annotate(geom = "text",x = 2, y = 0.8,label="41%",size=6)+
  annotate(geom = "text",x = 2, y = 0.53,label="12%",size=6)+
  annotate(geom = "text",x = 2, y = 0.25,label="47%",size=6)+
  annotate(geom = "text",x = 3, y = 0.8,label="45%",size=6)+
  annotate(geom = "text",x = 3, y = 0.5,label="10%",size=6)+
  annotate(geom = "text",x = 3, y = 0.37,label="18%",size=6)+
  annotate(geom = "text",x = 3, y = 0.15,label="27%",size=6)

ggsave("G:/图片/14.tiff", dpi=600)
###########3
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
miR_phen<-miR_phen[,-c(5,6)]
miR_phen1<-as.data.frame(apply(miR_phen[,2:4], 2, as.numeric))
miR_phen[,2:4]<-miR_phen1

#层次聚类
df<-scale(miR_phen[,2:4])
dist<-dist(df,method = "euclidean")
hc<-hclust(dist,method = "ward.D2")

nc<-NbClust(df,distance = "euclidean",min.nc = 2,max.nc = 15,method = "average")
par(mfcol=c(1,1))

clusters<-cutree(hc,k=2)
mir<-miR_phen
mir["cul"]<-clusters
colnames(mir)<-c("Name","mir-548v","mir-29a","mir-891a","Status","Transfer","Group")
sub1<-mir[which(mir$Group==1),]
sub2<-mir[which(mir$Group==2),]

mir$Group<-as.factor(mir$Group)
levels(mir$Group)<-c("S1","S2")

mir<-mir[order(mir$Transfer,decreasing=T),]
mir$Transfer<-factor(mir$Transfer,levels =c("Lymph Node Metastasis","Distance Metastases"))

#############################3
mir1<-mir[,c(1,7)]
score<-read.table("G:/colorectal_imm_sco.txt",header = T)
score$ID<-gsub(pattern="-", replacement=".", score$ID)
mir1$Name<-gsub('.{1}$', '', mir1$Name)
mir2<-merge(mir1,score,by.x = "Name",by.y = "ID")
mir2<-mir2[,-1]

sub1<-subset(mir2,Group%in%"S1")
sub2<-subset(mir2,Group%in%"S2")
p<-data.frame(Name="P",Stromal_score=NA,Immune_score =NA,ESTIMATE_score=NA)
for (i in 2:4) {
  p[1,i]<-wilcox.test(sub1[,i],sub2[,i])$p.value
}

c<-matrix(data=NA,ncol=2,nrow=3)
for(i in 2:4){
  m1<-median(sub1[,i])
  m2<-median(sub2[,i])
  
  c[i-1,1]<-colnames(sub1)[i]
  c[i-1,2]<-log2(m1/m2)
  
  
}
#log-fold-change
c<-as.data.frame(c)
colnames(c)<-c("imma","m1/m2")

mir2<-melt(mir2)
mir2$Group<-as.character(mir2$Group)

compaired<-list("Stromal_score")
ggplot(mir2) +
  aes(x = variable, y = value, fill = Group) +
  geom_boxplot(outlier.colour = NA,width=0.5)+
  scale_fill_manual(values = c("#2ECC71","#F1C40F"))+
  geom_jitter(aes(fill=Group),width =0.2,size=0.5)+
  geom_signif(y_position=c(1600), xmin=c(1.87), xmax=c(2.125), color="red" ,
              annotation=c("LFC= -2.74"), tip_length=c(0.04,0.15), size=0.5,
              textsize = 3.8, vjust = -0.3)  +
  geom_signif(y_position=c(1250), xmin=c(0.87), xmax=c(1.125), color="red" ,
              annotation=c("LFC= -0.310"), tip_length=c(0.05,0.02), size=0.5,
              textsize = 3.8, vjust = -0.3)  +
  geom_signif(y_position=c(2490), xmin=c(2.87), xmax=c(3.125), color="red" ,
              annotation=c("LFC= -0.72"), tip_length=c(0.1,0.03), size=0.5, 
              textsize = 3.8, vjust = -0.3)  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =14))+
  labs(x=NULL,y="value")+
  theme(axis.title = element_text(colour = "black",size = 14),
        axis.text = element_text(colour = "black",size = 14),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.direction = "vertical",
        #legend.key.width=unit(.6,"inches"),
        legend.background=element_rect(colour="#566573",size=0.4),
        legend.key.height=unit(.2,"inches"),
        legend.text=element_text(colour="black",size=13),
        legend.title=element_blank()) 
ggsave("G:/图片/19.tiff", dpi=600)  