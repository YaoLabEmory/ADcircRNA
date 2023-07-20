library(pheatmap)
library(RColorBrewer)
library(ggplot2)

RPM<-read.table(file="../RPM.xls",header=T,row.names=1)
RPM$logFC_5m<-log((apply(RPM[,c(10:12)],1,mean)+0.1)/(apply(RPM[,c(7:9)],1,mean)+0.1),2)
RPM$logFC_7m<-log((apply(RPM[,c(16:18)],1,mean)+0.1)/(apply(RPM[,c(13:15)],1,mean)+0.1),2)
RPM$logFC_6w<-log((apply(RPM[,c(4:6)],1,mean)+0.1)/(apply(RPM[,c(1:3)],1,mean)+0.1),2)

logFC<-RPM

I<-logFC[logFC$logFC_5m>0 & logFC$logFC_7m>0,]
II<-logFC[logFC$logFC_5m>0 & logFC$logFC_7m<0,]
III<-logFC[logFC$logFC_5m<0 & logFC$logFC_7m>0,]
IV<-logFC[logFC$logFC_5m<0 & logFC$logFC_7m<0,]

I<-I[order(-I$logFC_5m),]
II<-II[order(-II$logFC_5m),]
III<-III[order(-III$logFC_5m),]
IV<-IV[order(-IV$logFC_5m),]

RPM<-rbind(I,II,III,IV)
data<-RPM[,c(7:18)]

pdf(file="7mFAD.pdf",3,8)
pheatmap(data,show_rownames=F,show_colnames=F,cluster_rows=F,cluster_cols=F,scale="row",legend = FALSE)
dev.off()

cols = colorRampPalette(rev(brewer.pal(n = 7, name ="PRGn")))(100)
logFC<-RPM[,c(20:21)]
range <- max(abs(c(logFC$logFC_5m,logFC$logFC_7m)))

pdf(file="logFC.pdf",3,8)
pheatmap(logFC,show_rownames=F,cluster_rows=F,show_colnames=F,cluster_cols=F,legend = FALSE,color = cols,breaks = seq(-range, range,
length.out =100))
dev.off()

logFC<-RPM[,c(19:21)]

Iup<-I[I$logFC_6w>0,]
Idown<-I[I$logFC_6w<0,]
value<-c(length(Iup[,1]),length(Idown[,1]))
value
group<-c("up","down")
data<-data.frame(group=group,value=value)
p<-ggplot(data, aes(x="", y=value, fill=group)) +
	scale_fill_manual(values=c("blue","red"))+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()+
  theme(legend.position="none")
pdf(file="I.pie.pdf")
p
dev.off()

IIup<-II[II$logFC_6w>0,]
IIdown<-II[II$logFC_6w<0,]
value<-c(length(IIup[,1]),length(IIdown[,1]))
value
group<-c("up","down")
data<-data.frame(group=group,value=value)
p<-ggplot(data, aes(x="", y=value, fill=group)) +
        scale_fill_manual(values=c("blue","red"))+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()+
  theme(legend.position="none")
pdf(file="II.pie.pdf")
p
dev.off()

IIIup<-III[III$logFC_6w>0,]
IIIdown<-III[III$logFC_6w<0,]
value<-c(length(IIIup[,1]),length(IIIdown[,1]))
value
group<-c("up","down")
data<-data.frame(group=group,value=value)
p<-ggplot(data, aes(x="", y=value, fill=group)) +
        scale_fill_manual(values=c("blue","red"))+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()+
  theme(legend.position="none")
pdf(file="III.pie.pdf")
p
dev.off()

IVup<-IV[IV$logFC_6w>0,]
IVdown<-IV[IV$logFC_6w<0,]
value<-c(length(IVup[,1]),length(IVdown[,1]))
value
group<-c("up","down")
data<-data.frame(group=group,value=value)
p<-ggplot(data, aes(x="", y=value, fill=group)) +
        scale_fill_manual(values=c("blue","red"))+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()+
  theme(legend.position="none")
pdf(file="IV.pie.pdf")
p
dev.off()
