
txt<-c(1132,127)
txt<-as.data.frame(txt)
txt["Name"]<-c("Up-regulated","Down-regulated")
txt$Name<-as.factor(txt$Name)

txt$Name<- factor(txt$Name, levels=c("Up-regulated","Down-regulated"), ordered=TRUE)
p<-ggplot(txt) +
  aes(x = Name, fill = Name, weight = txt) +
  geom_bar(width = 0.6)+
  scale_fill_manual(values = c("#C0392B","#2980B9"))+
  theme(panel.background = element_rect(fill = NA),
        legend.position = "none",
        panel.border = element_blank(),
        axis.line = element_line(size = 0.8),
        axis.title = element_text(colour = "black",size = 14),
        axis.text = element_text(colour = "black",size = 14),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x = "Oligometastatic DEGs",y="Number of gene")+
  annotate(geom = "text",x = 1, y = 1180, hjust = 0.5,size=5,label="1132")+
  annotate(geom = "text",x = 2, y = 180, hjust = 0.5,size=5,label="127")
ggsave("G:/图片/1.tiff", dpi=600)
dev.off()




