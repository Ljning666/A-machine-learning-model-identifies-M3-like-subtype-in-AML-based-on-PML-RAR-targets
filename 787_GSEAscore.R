library(ggplot2)
library(ggpubr)
library(fgsea)
library(pROC)
##pml_rara_geneset
target_down<-read.table("E:\\Mirror\\Classifier\\data\\PML_RARA_target\\down.txt",stringsAsFactors = F)
target_up<-read.table("E:\\Mirror\\Classifier\\data\\PML_RARA_target\\up.txt",stringsAsFactors = F)
list_down<-list(target_down=target_down$V1)
list_up<-list(target_up=target_up$V1)
AML_Sample_info<-read.delim("E:\\Mirror\\Classifier\\data\\AML_Sample_info_sur.txt",stringsAsFactors = F)
case_data1_2<-AML_Sample_info[which(AML_Sample_info$Dataset!=3&AML_Sample_info$GSE!="own"),]
case_data1_2<-case_data1_2[which(case_data1_2$Tissue=="BM"&case_data1_2$FAB!="unknown"),]
control_data1_2<-AML_Sample_info[which(AML_Sample_info$Dataset!=3&AML_Sample_info$GSE!="own"),]
control_data1_2<-control_data1_2[which(control_data1_2$Tissue=="BM"&control_data1_2$Disease=="healthy"),]
GSE122511_Dataset_2<-read.delim("E:\\Mirror\\Classifier\\data\\GSE122511_Dataset_2.txt",check.names = F,stringsAsFactors = F)
###training cohort model########################################################
case_data2<-GSE122511_Dataset_2[,case_data1_2$GSM[which(case_data1_2$Dataset==2)]]
control_data2<-GSE122511_Dataset_2[,control_data1_2$GSM[which(control_data1_2$Dataset==2)]]
control<-control_data2
case<-case_data2
control<-2^control
case<-2^case
control<-apply(control,1,mean)
data2_logFC<-apply(case,2,function(x){return(log2(x/control))})
for(i in 1:1){
  ref<-list()
  ref[[1]]<-list_down[[i]]
  ref[[2]]<-list_up[[i]]
  names(ref)<-c("down","up")
  res2 <- apply(data2_logFC, 2, function(x){
    names(x) <- rownames(data2_logFC)
    x <-x[order(x,decreasing = T)]
    fgsea_res<-fgseaMultilevel(pathways= ref,
                               stats=x,
                               minSize = 1,maxSize = 2000,nproc  =  0,
                               gseaParam  =  1,
                               BPPARAM  =  NULL
    )
    fgsea_res2 <- c(fgsea_res$ES[1],fgsea_res$ES[2],fgsea_res$pval[1],fgsea_res$pval[2])
    return(fgsea_res2)
  }
  )
  res2<-t(res2)
  colnames(res2)<-c("down","up","down_pval","up_pval")  
  res2<-as.data.frame(res2)
  res2$score<-res2$down-res2$up
  res2$score_norm<-unlist(lapply(res2$score,function(x){
    return((x-min(res2$score))/(max(res2$score)-min(res2$score)))
  }))
  res2$FAB<-AML_Sample_info$FAB[match(rownames(res2),AML_Sample_info$GSM)] write.table(res2,"E:\\Mirror\\Classifier\\data\\PML_RARA_target\\GSEA_score\\data2_score.txt",quote = F,sep = "\t")
###AUC
data2_score<-read.delim("E:\\Mirror\\Classifier\\data\\PML_RARA_target\\GSEA_score\\data2_score.txt",stringsAsFactors = F)
  AUC_data<-data2_score[,6:7]
  AUC_data$FAB_type[which(AUC_data$FAB!="M3")]<-0
  AUC_data$FAB_type[which(AUC_data$FAB=="M3")]<-1
  colnames(AUC_data)<-c("score","FAB","FAB_type")
  p_AUC<-roc(AUC_data$FAB_type,AUC_data$score)
  pdf("E:\\Mirror\\Classifier\\pic\\gene787score\\data2_ROC.pdf",height = 6,width = 6)
  plot(p_AUC,print.auc=T,auc.polygon=T,
       grid=c(0.2,0.2),
       grid.col=c("green","red"),
       max.auc.polygon=T,
       legacy.axes = TRUE,
       auc.polygon.col="skyblue",
       print.thres=T,
       xlim=c(1,0),
       xlab = "1-Specificity", ylab = "Sensitivity",
       main="data2_787gene")
  dev.off() 
}
data2_score$type[which(data2_score[,6]<0.560)]<-"other"
data2_score$type[which(data2_score[,6]>=0.560)]<-"M3_like"
data<-as.data.frame(table(data2_score$FAB,data2_score$type))
p<-ggplot(data,aes(Var1,Freq,fill=Var2))+
  geom_bar(stat="identity",position="fill")+
  ylab("Ratio") +
  xlab("Type")+
  theme(panel.grid=element_blank())+
  scale_fill_manual(values = c("M3_like"="#E31A1C","other"="#1F78B4"),name="Var2")+
  coord_flip()+
  ggtitle("data2_787gene")
pdf("E:\\Mirror\\Classifier\\pic\\gene787score\\data2_787gene.pdf",width=8,height=6)
print(p)
dev.off()