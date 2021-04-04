library(reshape2)
library(ggplot2)
library(ggsignif)
library(gridExtra)
#读取数据
setwd("G:/colon")
#读取表型数据
phenotype<-read.delim("TCGA-COAD.GDC_phenotype.tsv")
phenotype<-phenotype[,c(1,93,39,40,71,75)]
colnames(phenotype)<-c("samples","tumor_stage","M","N","death","status")
#删除无用值
phenotype<-phenotype[which(phenotype$M!=""),]
phenotype<-phenotype[which(phenotype$M!="MX"),]
#分组
phenotype["Group"]<-NA
phenotype[which(phenotype$M=="M0"),7]="Primary"
phenotype[which(is.na(phenotype$Group)),7]="OM"
table(phenotype$Group)
####################################################
#将横杠用.替换
phenotype$samples<-gsub(pattern="-", replacement=".", phenotype$samples)
#读取miR数据
miR<-read.delim("batch_mir.txt",header = T)
mir<-read.delim("mir.txt",header = T)
mir$mir<-tolower(mir$mir)
miR<-miR[grep("mir",miR$miRNA_ID),]
#删除0超过一般的数据
miR<-miR[rowMeans(miR == 0)<0.8,]
#将宽矩阵变为长矩阵
miR<-melt(miR)
miR<-merge(miR,phenotype,by.x = "variable",by.y = "samples")
MIR<-merge(mir,miR,by.x = "mir",by.y = "miRNA_ID")
Name<-unique(MIR$mir)
#检验正态性
b<-matrix(data=NA,ncol=4,nrow=17)
for(i in 1:length(Name)){
  mir1<-MIR[which(MIR$mir==Name[i]),]
  b[i,1]<-Name[i]
  b[i,2]<-shapiro.test(mir1$value)$p.value
  b[i,3]<-shapiro.test(mir1[which(mir1$Group=="Primary"),3])$p.value
  b[i,4]<-shapiro.test(mir1[which(mir1$Group=="OM"),3])$p.value
  
  
}
#计算FC值
c<-matrix(data=NA,ncol=2,nrow=17)
for(i in 1:length(Name)){
  mir1<-MIR[which(MIR$mir==Name[i]),]
  m1<-mean(mir1[which(mir1$Group=="Primary"),3])
  m2<-mean(mir1[which(mir1$Group=="OM"),3])
  
  c[i,1]<-Name[i]
  c[i,2]<-log2(m2/m1)
  
  
}
#log-fold-change
c<-as.data.frame(c)
colnames(c)<-c("miR","m2/m1")

a<-list()
for (i in 1:length(Name)) {
  mir1<-MIR[which(MIR$mir==Name[i]),]
  mir1$Group<-as.factor(mir1$Group)
  mir1$Group <- factor(mir1$Group, levels=c("Primary", "OM"), ordered=TRUE)
  compaired <- list(c("Primary","OM"))
  name1<-unique(mir1$mir)
  name1<-gsub(pattern="hsa-mir", replacement="miR",name1)
  p<-ggplot(mir1) +
    aes(x = Group, y = value, fill = Group) +
    geom_boxplot(outlier.colour = NA,width=0.5) +
    scale_fill_manual(values = c("#52BE80","#F1C40F"))+
    geom_jitter(aes(fill=Group),width =0.2,size=0.5)+
    geom_signif(comparisons = compaired,step_increase = 0.1,
                color="#E74C3C",
                textsize =5 ,
                map_signif_level=function(p)sprintf("p = %.2g", p),
                test = t.test)+###显著性
    theme_bw() +
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(size = 0.8),
          axis.title = element_text(colour = "black",size = 14),
          axis.text = element_text(colour = "black",size = 14),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    labs(x=name1,y="Expression Level")
   # ggtitle(unique(mir1$mir)) +
   # theme(plot.title = element_text(hjust = 0.5,face="bold",size =14)) #设置标题居中
  a<-c(a,list(p))
}
grid.arrange(a[[1]],a[[2]],a[[3]],a[[4]],a[[5]],a[[6]],ncol=3,nrow=2)
grid.arrange(a[[7]],a[[8]],a[[9]],a[[10]],a[[11]],a[[12]],ncol=3,nrow=2)
grid.arrange(a[[13]],a[[14]],a[[15]],a[[16]],a[[17]],ncol=3,nrow=2)


ggarrange(a[[6]],a[[15]],a[[17]],nrow = 2, ncol = 2)%>%
  ggexport(filename = "G:/图片/9.tiff",width = 5700,
           height = 5700,res=600)
