setwd("G:/colon")
library(sva)
phen <-read.delim("TCGA-COAD.GDC_phenotype.tsv")
phen<-phen[1:447,c(1,3,93)]
phen$submitter_id.samples<-gsub(pattern="-", replacement=".", phen$submitter_id.samples)
edata <- read.delim("TCGA-COAD.mirna.tsv",header = T)
edata<-as.data.frame(t(edata))
colnames(edata)<-edata[1,]
edata<-edata[-1,]
edata['Name']<-rownames(edata)
rownames(edata)<-1:461
res<-merge(edata,phen,by.x = "Name",by.y = "submitter_id.samples")
res1<-as.data.frame(apply(res[,2:1882], 2, as.numeric))
res2<-as.data.frame(t(res1))
colnames(res2)<-res$Name
res2<-res2[rowMeans(res2 == 0)<0.8,]
phe<-res[,c(1,1883,1884)]
batch<-as.factor(res$batch_number)
levels(batch)<-1:20
batch<-as.numeric(batch)
mod = model.matrix(~as.factor(tumor_stage.diagnoses),data=phe)
exp2 = ComBat(dat=res2, batch=batch, mod=mod, par.prior=TRUE, ref.batch="1")
write.table(exp2,"batch_mir.txt",quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
