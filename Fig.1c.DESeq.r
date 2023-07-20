#!/apps/R-4.0.2/bin/R

AllEqual <- structure(function(
        ##title<<
        ## Check if all values in a vector are the same
        ##description<<
        ## This function is used to check if all values in a vector are equal. It can be used for example to check if a time series contains only 0 or NA values.

        x
        ### numeric, character vector, or time series of type ts
) {
        res <- FALSE
        x <- na.omit(as.vector(x))
        if (length(unique(x)) == 1 | length(x) == 0) res <- TRUE
        return(res)
        ### The function returns TRUE if all values are equal and FALSE if it contains different values.
},ex=function(){
# check if all values are equal in the following vectors:
AllEqual(1:10)
AllEqual(rep(0, 10))
AllEqual(letters)
AllEqual(rep(NA, 10))
})

config<-read.table(file="../../bin/config")
length<-config[which(config[,1]=="length"),2]

ARS<-read.table(file="sample_names_AR.txt")[,1]
RMS<-read.table(file="sample_names_RM.txt")[,1]
RMJ<-read.table(file=paste0("../../MergeCirc/juncmap/",length,"bp/summary.txt"),row.names=1,header=T)
RMJ<-RMJ[RMS,2]
RMM<-read.table(file="../../RNAseq/summary.txt",header=T,row.name=1)
RMM<-RMM[RMS,2]
NF<-RMM/RMJ
NF<-NF/(min(NF))
NF

RC<-read.table(file="../../MergeCirc/juncmap/8bp/RC",header=T,row.names=1,check.names = FALSE)
RC<-RC[,ARS]
RC<-sweep(RC,MARGIN = 2,NF,FUN = "/")

norm<-apply(RC,2,sum)
norm<-norm/(min(norm))
RC<-sweep(RC,MARGIN = 2,norm,FUN = "/")

RC<-RC[(RC[,1]>2 & RC[,2]>2 & RC[,3]>2) | (RC[,4]>2 & RC[,5]>2 & RC[,6]>2),]

designs=c(rep(0,3),rep(1,3))
n0=sum(designs==0)
n1=sum(designs==1)
m0=apply(RC[,designs==0],1,mean)
m1=apply(RC[,designs==1],1,mean)
s0=apply(RC[,designs==0],1,sd)
s1=apply(RC[,designs==1],1,sd)

n=nrow(RC)
pval=rep(1,n)

for(i in 1:n){
   if(AllEqual(as.numeric(RC[i,c(1:6)]))){pval[i]=1}else{
#   else if(AllEqual(as.numeric(RC[i,c(1:3)])) & AllEqual(as.numeric(RC[i,c(4:6)]))){pval[i]=0} else {
   z=(m1[i]-m0[i])/sqrt(s1[i]**2/n1+s0[i]**2/n0)
   #z=(m1[i]-m0[i])/sqrt(m1[i]/n1+m0[i]/n0)
   pval[i]=2*pnorm(-abs(z))
   }
}
fdr=p.adjust(pval,method='fdr')

logFC=log((m1+0.1)/(m0+0.1),2)

res=cbind(RC=RC,logFC=logFC,pval=pval, fdr=fdr)
res<-as.data.frame(res)
res=res[order(abs(res$pval),decreasing = F),]

genename<-read.table(file="../../MergeCirc/circlist/circlist.bed.gene")
names(genename)<-c("chr","start","end","strand","isoform","score","gene")
row.names(genename)<-paste(genename$chr,":",genename$start,"-",genename$end,sep="")
res$gene<-genename[row.names(res),]$gene

write.table(res,file="DESeq.xls",col.names=T,row.names=T,sep="\t",quote=FALSE)

up<-res[res$fdr<0.05 & res$logFC>0,]
down<-res[res$fdr<0.05 & res$logFC<0,]
dim(up)
dim(down)
