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
#处理为分类器可用数据
MIR<-merge(miR,mir,by.y = "mir",by.x = "miRNA_ID")
MIR<-as.data.frame(t(MIR))
colnames(MIR)<-MIR[1,]
MIR<-MIR[-1,]
MIR["NAME"]<-rownames(MIR)
rownames(MIR)<-1:428
MIR<-merge(MIR,phenotype,by.x = "NAME",by.y = "samples")
MIR<-MIR[,-c(1,4,21,22,23,24)]
MIR1<-as.data.frame(lapply(MIR[,1:17],as.numeric))
MIR1['Group']<-as.factor(MIR$Group)
MIR1['Stage']<-MIR$tumor_stage

MIR1[which(MIR1$Stage=="not reported"),19]=1
MIR1[which(MIR1$Stage=="stage i"),19]=1
MIR1[which(MIR1$Stage=="stage ii"),19]=2
MIR1[which(MIR1$Stage=="stage iia"),19]=2
MIR1[which(MIR1$Stage=="stage iib"),19]=2
MIR1[which(MIR1$Stage=="stage iic"),19]=2
MIR1[which(MIR1$Stage=="stage iii"),19]=3
MIR1[which(MIR1$Stage=="stage iiia"),19]=3
MIR1[which(MIR1$Stage=="stage iiib"),19]=3
MIR1[which(MIR1$Stage=="stage iiic"),19]=3
MIR1[which(MIR1$Stage=="stage iv"),19]=4
MIR1[which(MIR1$Stage=="stage iva"),19]=4
MIR1[which(MIR1$Stage=="stage ivb"),19]=4

A<-MIR1[which(MIR1$Stage==1),]
B<-MIR1[which(MIR1$Stage==2),]
C<-MIR1[which(MIR1$Stage==3),]
D<-MIR1[which(MIR1$Stage==4),]
set.seed(12)
a<-A[sample(nrow(A), 32), ]
b<-B[sample(nrow(B), 32), ]
c<-C[sample(nrow(C), 32), ]
N_MIR<-rbind(a,b,c,D)

N_MIR<-N_MIR[,-19]

N_MIR3<-N_MIR[,c(15,6,17,18)]

###########
library(plyr)
library(randomForest)
library(pROC)
CVgroup <- function(k,datasize,seed){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]    #将数据分成K份，并生成的完成数据集n
  temp <- sample(n,datasize)   #把n打乱
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])  #dataseq中随机生成k个随机有序数据列
  return(cvlist)
}
Res<-NULL
for (a in 1:100) {
  
  
  k <- 5
  datasize <- nrow(N_MIR3)#换要需要的数据集
  cvlist <- CVgroup(k = k,datasize = datasize,seed = a)#1206
  cvlist
  data <- N_MIR3 #换要需要的数据集
  pred <- data.frame()   #存储预测结果
  
  m <- seq(60,500,by = 20)  #如果数据量大尽量间隔大点，间隔过小没有实际意义
  for(j in m){   #j指的是随机森林的数量
    progress.bar <- create_progress_bar("text")  #plyr包中的create_progress_bar函数创建一个进度条，
    progress.bar$init(k)   #设置上面的任务数，几折就是几个任务
    for (i in 1:k){
      train <- data[-cvlist[[i]],]  #刚才通过cvgroup生成的函数
      test <- data[cvlist[[i]],]
      model1 <-randomForest(Group~.,data = train,ntree = j)   #建模，ntree=j 指的树数
      prediction <- predict(model1,subset(test,select = -Group))   #预测
      roc<-multiclass.roc (as.ordered(test$Group) ,as.ordered(prediction))
      randomtree <- rep(j,length(prediction))   #随机森林树的数量
      kcross <- rep(i,length(prediction))   #i是第几次循环交叉，共K次
      A<-as.numeric(roc$auc)
      temp <- data.frame(cbind(subset(test,select = Group),prediction,randomtree,kcross,A))#真实值、预测值、随机森林树数、预测组编号捆绑在一起组成新的数据框tenp
      pred <- rbind(pred,temp)   #temp按行和pred合并
      print(paste("随机森林：",j))  #循环至树数j的随机森林模型
      progress.bar$step() #输出进度条。告知完成了这个任务的百分之几
    }
  }
  
  res<-cbind(a,pred[which.max(pred$A),])
  Res<-rbind(Res,res)
  
}

#############
##############
#MIR17,seed=89,AUC-0.845
library(pROC)
library(ROCR)
train0 <- data[-cvlist[[2]],]  #刚才通过cvgroup生成的函数
test0 <- data[cvlist[[2]],]

mod <-randomForest(Group~.,data = train0,ntree = 340)   #建模，ntree=j 指的树数

prediction <- predict(mod,subset(test0,select = -Group))   #预测
roc1<-roc(as.ordered(test0$Group) ,as.ordered(prediction))
plot(roc1,print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), grid.col=c("green","red"), 
     max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T)
table <- table(prediction,test0$Group)  
#预测准确率
sum(diag(table))/sum(table)  

#install.packages("ROCR")
#library(ROCR)
pred_out_0<-predict(mod,subset(test0,select = -Group),"prob")  
pred0<-prediction(pred_out_0[,2],test0$Group)
perf0<-performance(pred0,"tpr","fpr")
plot(perf0)
###################################################
#MIR3-seed=51 p=0.794
train1 <- data[-cvlist[[3]],]  #刚才通过cvgroup生成的函数
test1 <- data[cvlist[[3]],]

mod1 <-randomForest(Group~.,data = train1,ntree = 240)   #建模，ntree=j 指的树数

prediction1 <- predict(mod1,subset(test1,select = -Group))   #预测
roc2<-roc(as.ordered(test1$Group) ,as.ordered(prediction1))
plot(roc2,print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), grid.col=c("green","red"), 
     max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T)
table <- table(prediction1,test1$Group)  
#预测准确率
sum(diag(table))/sum(table)  

#install.packages("ROCR")
#library(ROCR)
pred_out_1<-predict(mod1,subset(test1,select = -Group),"prob")  
pred1<-prediction(pred_out_1[,2],test1$Group)
perf<-performance(pred1,"tpr","fpr")
plot(perf,colorize=TRUE)
##########################################################################
# Create a basic roc object

r1 <- roc(as.ordered(test0$Group) ,as.ordered(pred_out_0[,2]))
r1
r2 <- roc(as.ordered(test1$Group) ,as.ordered(pred_out_1[,2]))
r2
ggroc(r1,legacy.axes = TRUE,show_guide=FALSE,size=1,color="#5499C7")+
  theme_bw() + 
  #scale_color_manual(values ="#5499C7")+
  geom_abline(slope = 1, intercept = 0,color="grey",size=0.8) +
  #geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey",size=0.8)+
  #theme_classic()+
  theme(#panel.border=element_rect(color="black",size=0.8),
    #axis.line = element_line(size = 0.8),
    axis.title = element_text(colour = "black",size = 14),
    axis.text = element_text(colour = "black",size = 14),
    axis.title.x = element_text(vjust=-1),
    axis.title.y = element_text(vjust=2),
    plot.title = element_text(hjust = 0.5),
    
    #
    panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(title = "ROC Curve (17 miRNAs)")+
  labs(x = "False positive Rate",y="Ture positive Rate")+
  #annotate("rect", xmin=0.44, xmax=0.9, ymin=0.2 , ymax=0.35,alpha=0, color="grey")+
  #annotate(geom = "line", x = c(0.45,0.51),y = 0.31,size=1,color="#5499C7") +
  annotate(geom = "text",size=5,
           x = 0.7, y = 0.25,
           label = "AUC=0.8019")
#annotate(geom = "line",x = c(0.45,0.51), y = 0.24,size=1,color="#F39C12") +
#annotate(geom = "text",x = 0.7, y = 0.25,label = "miRNA: 3 (AUC:0.8333)")
ggsave("G:/ͼƬ/10.tiff", dpi=600)

ggroc(r2,legacy.axes = TRUE,show_guide=FALSE,size=1,color="#F39C12")+
  theme_bw() + 
  #scale_color_manual(values ="#5499C7")+
  geom_abline(slope = 1, intercept = 0,color="grey",size=0.8) +
  #geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey",size=0.8)+
  #theme_classic()+
  theme(#panel.border=element_rect(color="black",size=0.8),
    #axis.line = element_line(size = 0.8),
    axis.title = element_text(colour = "black",size = 14),
    axis.text = element_text(colour = "black",size = 14),
    axis.title.x = element_text(vjust=-1),
    axis.title.y = element_text(vjust=2),
    plot.title = element_text(hjust = 0.5),
    #
    panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(title = "ROC Curve (3 miRNAs)")+
  labs(x = "False positive Rate",y="Ture positive Rate")+
  #annotate("rect", xmin=0.44, xmax=0.9, ymin=0.2 , ymax=0.35,alpha=0, color="grey")+
  #annotate(geom = "line", x = c(0.45,0.51),y = 0.31,size=1,color="#5499C7") +
  annotate(geom = "text",size=5,
           x = 0.7, y = 0.25,
           label = "AUC=0.8377")
#annotate(geom = "line",x = c(0.45,0.51), y = 0.24,size=1,color="#F39C12") +
#annotate(geom = "text",x = 0.7, y = 0.25,label = "miRNA: 3 (AUC:0.8333)")
ggsave("G:/ͼƬ/11.tiff", dpi=600)

res1<-rbind(c("Accracy",0.875),c("Balanced Accuracy",0.846),c("Sensitivity",0.778),
            c("Specificity",0.913),c("Positive Prediction Value(PPV)",0.778),
            c("Negative Prediction Value(NPV)",0.913),c("Matthew`s Correlation Coefficient",0.691),
            c("AUC",0.802))
res1<-as.data.frame(res1)
colnames(res1)<-c("Preformance metrics","Values")

library(xtable)
library(flextable)
library(officer)
tt<-xtable_to_flextable(xtable(res1))
doc = read_docx()
doc = body_add_flextable(doc,tt)
print(doc,"G:/ͼƬ/ROC_17.docx.docx")

res2<-rbind(c("Accracy",0.844),c("Balanced Accuracy",0.794),c("Sensitivity",0.636),
            c("Specificity",0.952),c("Positive Prediction Value(PPV)",0.875),
            c("Negative Prediction Value(NPV)",0.833),c("Matthew`s Correlation Coefficient",0.646),
            c("AUC",0.838))
res2<-as.data.frame(res2)
colnames(res2)<-c("Preformance metrics","Values")

tt<-xtable_to_flextable(xtable(res2))
doc = read_docx()
doc = body_add_flextable(doc,tt)
print(doc,"G:/ͼƬ/ROC_3.docx")
