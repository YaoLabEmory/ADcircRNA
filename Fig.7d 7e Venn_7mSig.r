firstup<-read.table(file="../firstup",header=F)[,1]
firstdown<-read.table(file="../firstdown",header=F)[,1]

length(firstup)
length(firstdown)

firstup<-toupper(firstup)
firstdown<-toupper(firstdown)

Upgene<-read.table(file="/projects/compbio/users/yli/RNAseq/ADmouseCircRNAAtailing/RNAseq/Cuffdiff/7m_5xFADvsWT_diffout/CaseUpgene",header=F)[,1]
Downgene<-read.table(file="/projects/compbio/users/yli/RNAseq/ADmouseCircRNAAtailing/RNAseq/Cuffdiff/7m_5xFADvsWT_diffout/CaseDowngene",header=F)[,1]
DEgene<-c(Upgene,Downgene)
DEgene<-toupper(DEgene)

length(intersect(firstup,DEgene))
length(intersect(firstdown,DEgene))

CPSF6KD_Up<-read.table(file="/projects/compbio/users/yli/RNAseq/CPSF6_KD/RNAseq/Cuffdif/CPSF6_KDvsWT_diffout/CaseUpgene",header=F)[,1]
CPSF6KD_Down<-read.table(file="/projects/compbio/users/yli/RNAseq/CPSF6_KD/RNAseq/Cuffdif/CPSF6_KDvsWT_diffout/CaseDowngene",header=F)[,1]

CPSF6OE_Up<-read.table(file="/projects/compbio/users/yli/RNAseq/CPSF6_KD/RNAseq/Cuffdif/CPSF6_OEvsWT_diffout/CaseUpgene",header=F)[,1]
CPSF6OE_Down<-read.table(file="/projects/compbio/users/yli/RNAseq/CPSF6_KD/RNAseq/Cuffdif/CPSF6_OEvsWT_diffout/CaseDowngene",header=F)[,1]

CPSF6Target<-c(CPSF6OE_Up,CPSF6OE_Down)

up<-intersect(firstup,CPSF6Target)
down<-intersect(firstdown,CPSF6Target)

length(unique(CPSF6Target))

length(up)
length(down)

up<-as.character(up)
down<-as.character(down)

DE7m<-read.table(file="/projects/compbio/users/yli/RNAseq/ADmouseCircRNAAtailing/RNAseq/Cuffdiff/7m_5xFADvsWT_diffout/gene_exp.diff",header=T)
DE5m<-read.table(file="/projects/compbio/users/yli/RNAseq/ADmouseCircRNAAtailing/RNAseq/Cuffdiff/5m_5xFADvsWT_diffout/gene_exp.diff",header=T)

DE7m<-DE7m[,c(1,10,13)]
DE5m<-DE5m[,c(1,10)]

DE<-cbind(DE5m,DE7m)
row.names(DE)<-toupper(DE[,1])
DE<-DE[,-c(1,3)]
names(DE)<-c("logFC5m","logFC7m","FDR")

up<-DE[up,]
down<-DE[down,]

up<-up[!is.na(up$logFC5m),]
down<-down[!is.na(down$logFC5m),]

up<-up[up$logFC7m*up$logFC5m>0 & abs(up$logFC7m)>abs(up$logFC5m) & up$FDR<0.05,]
down<-down[down$logFC7m*down$logFC5m>0 & abs(down$logFC7m)>abs(down$logFC5m) & down$FDR<0.05,]

up<-up[!is.na(up$logFC5m),]
down<-down[!is.na(down$logFC5m),]

up<-row.names(up)
down<-row.names(down)

up
down

library(ggplot2)
library(reshape2)

DE<-DE[c(up,down),]
write.table(DE,file="GeneListSelectedlogFCFDR.xls",row.names=T,col.names=T,sep="\t",quote=F)

DE$gene<-row.names(DE)
DE<-DE[,-3]
DE$gene<-factor(DE$gene,levels=DE[order(-DE$logFC7m),]$gene)
DE<-melt(DE)

names(DE)<-c("Gene","Age","LogFC")

p<-ggplot(DE, aes(x=Gene, y=LogFC, fill=Age)) +geom_bar(position="dodge", stat="identity",width = 0.5)+
        scale_fill_manual(values = c("logFC5m" = "blue", "logFC7m" = "red"))+
        scale_colour_manual(values = c("logFC5m" = "blue", "logFC7m" = "red"))
p<-p+theme()+
            xlab("Gene")+
            ylab("logFC")+
            theme_classic()+
            theme( axis.title.x = element_text(size = 0),
            axis.text.y = element_text(size = 20,colour="black"),
            axis.text.x = element_text(size = 20,colour="black"),
            axis.title.y = element_text(size = 20,colour="black"))+
            theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1))+
            theme(legend.text=element_text(size=20,colour="black"))+
            theme(legend.title=element_text(size=0))+
#            theme(legend.position='none')+
            theme(axis.line = element_line(colour = 'black', size = 1))

pdf(file="logFC.pdf",10,4)
p
dev.off()


DE<-read.table(file="GeneListSelectedlogFCFDR.xls",header=T,row.names=1)
DE$gene<-row.names(DE)
DE5m<-read.table(file="/projects/compbio/users/yli/RNAseq/ADmouseCircRNAAtailing/RNAseq/Cuffdiff/5m_5xFADvsWT_diffout/gene_exp.diff",header=T)
DE5m<-DE5m[,c(1,13)]
names(DE5m)<-c("gene","FDR5m")
DE5m$gene<-toupper(DE5m$gene)

DE7m5m<-merge(DE,DE5m,by="gene")

ADRisk<-read.csv(file="/projects/compbio/users/yli/Download/AD_RiskGene/ADriskGene.csv",header=T)
ADRisk<-ADRisk[,c(1,6)]
dim(ADRisk)
names(ADRisk)<-c("RBP","RiskScore")
ADRisk<-ADRisk[ADRisk$RiskScore>2,1]

intersect(DE7m5m$gene,ADRisk)

