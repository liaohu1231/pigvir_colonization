data=read.csv("all.MAGs_nucl_diversity.tsv",
              check.names = F,row.names = 1,sep='\t',header = 1)
metawrap_F3_2=data[grep(pattern = "metawrap_F3-2.bin.10",row.names(data)),]
group=read.table("../metadata.txt",header = 1,sep='\t')
metawrap_F3_2$contigs=rownames(metawrap_F3_2)
library(reshape2)
metawrap_F3_2=melt(metawrap_F3_2)
rownames(group)=group$Sample
metawrap_F3_2$group=group[as.character(metawrap_F3_2$variable),4]
for (i in 1:nrow(metawrap_F3_2)){
  if (metawrap_F3_2[i,3]==0){
    metawrap_F3_2[i,3]=NA
  }
}
metawrap_F3_2=na.omit(metawrap_F3_2)
metawrap_F3_2$group=factor(metawrap_F3_2$group,levels = c("PPF","T1","T3","T7","T14",
                                                          "F1","F3","F7","F14"))
library(ggpubr)
library(ggplot2)
ggplot(metawrap_F3_2,aes(group,value))+
  geom_point(stat = 'identity')+
  geom_boxplot(aes(fill=group))+theme_bw()+
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("PPF","T3"),
                                        c("T3","T7"),
                                        c("PPF","F1"),
                                        c("F1","F3"),
                                        c("F3","F14"),
                                        c("F3","T3")))

data$scaffolds=rownames(data)
data_average=matrix(data = NA,nrow = 10,ncol = ncol(data)-1)
for (i in 1:ncol(data_average)) {
  tmp=as.data.frame(data[data[,i]>0,i])
  tmp=na.omit(tmp)
  for(j in 1:10){
    samp=sample(tmp[,1],size = 500,replace = F)
    tmp2=mean(samp)
    data_average[j,i]=tmp2
  }
}
colnames(data_average)=colnames(data[,1:52])
total=melt(as.data.frame(data_average))
total=total[total$value>0,]
total$group=group[as.character(total$variable),4]
total$group=factor(total$group,levels = c("PPF","CK0","CK1","CK3","CK7","CK14","CK56",
                                          "T1","T3","T7","T14","T56",
                                          "F1","F3","F7","F14","F56"))
total$days=as.factor(group[as.character(total$variable),3])
total$days=factor(total$days,levels = c("0","1","3","7","14","56"))
total$treatment=group[as.character(total$variable),2]
total$treatment=factor(total$treatment,levels = c("PPF","CK","F","T"))
total=na.omit(total)
library(ggpubr)
library(ggplot2)
ggplot(total,aes(days,value,fill=group))+
  geom_boxplot(aes(fill=treatment),outlier.shape = NA,alpha=0.8,
               show.legend = T)+theme_bw()+
  stat_compare_means(method = "t.test",
                     comparisons = list(c("PPF","CK"),
                                        c("CK","F"),c("F","T")))+
  scale_fill_d3()


ggplot(total[grep("56",total$days),],aes(treatment,value,fill=treatment))+
  geom_boxplot(aes(fill=treatment),outlier.shape = NA,alpha=0.8,
               show.legend = T)+theme_bw()+
  stat_compare_means(method = "t.test",
                     comparisons = list(c("CK","T"),
                                        c("CK","F"),c("F","T")))+
  scale_fill_d3()
ggsave("micro56_mi.pdf",device = 'pdf',width = 3,height = 4)
###viral_diversity
viral_mi=read.csv("all_pigviruses_scaffold_nul.tsv",
                    header = T,row.names = 1,check.names = F,sep='\t')
viromes_metadata=data.frame(row.names = colnames(viral_mi),
                            days=rep(c("14","56","14","0","14","56"),
                                     c(3,1,3,5,3,2)))
tmp=strsplit(rownames(viromes_metadata),split = '-')
tmp=do.call(rbind,tmp)
viromes_metadata$group=tmp[,1]

library(reshape2)
viral_mi=melt(as.matrix(viral_mi))
viral_mi$group=viromes_metadata[as.character(viral_mi$Var2),2]
viral_mi=viral_mi[viral_mi$value>0,]
viral_mi=na.omit(viral_mi)
viral_mi=viral_mi[-grep("CK",viral_mi$group),]
viral_mi=viral_mi[-grep("PS0",viral_mi$group),]
viral_mi$group=factor(viral_mi$group,levels = c("PPF0","PF14","PT14","PT56"))
library(ggplot2)
library(ggpubr)
ggplot(viral_mi,aes(group,value,fill=group))+
  geom_boxplot(show.legend = F,alpha=0.8)+
  theme_bw()+
  labs(y="Microdiversity (pi)")+
  stat_compare_means(comparisons = list(c("PF14","PPF0"),
                                        c("PF14","PT14")))+
  scale_fill_d3()

Acine=read.table("pigvir_acine_phage.txt")
rownames(Acine)=Acine$V1
viral_mi$Acine=Acine[as.character(viral_mi$Var1),1]
viral_mi_acine=na.omit(viral_mi)
ggplot(viral_mi_acine,aes(group,value,fill=group))+
  geom_boxplot(alpha=0.8)+
  theme_bw()+
  labs(y="Microdiversity (pi)")+
  scale_fill_d3()+
  stat_compare_means(comparisons = list(c("PF14","PPF0"),
                                        c("PF14","PT14")))
for (i in 1:nrow(viral_mi)) {
  if(is.na(viral_mi$Acine[i])==TRUE){
    viral_mi$Acinetobacter[i]="Non_Acinetobacter"
  }else{
    viral_mi$Acinetobacter[i]="Acinetobacter"
  }
}
ggplot(viral_mi,aes(group,value,fill=Acinetobacter))+
  geom_boxplot()+
  theme_bw()+
  labs(y="Microdiversity (pi)")


###metagenome
meta_vir_m=read.table("all_pigviruses_metagenome_scaffold_nul.tsv",
                      row.names = 1,sep='\t',header = 1,check.names = F)
viral_mi=read.csv("all_pigviruses_scaffold_nul.tsv",
                  header = T,row.names = 1,check.names = F,sep='\t')
library(dplyr)
source('new.cbind.R')
meta_vir_m=new.cbind(meta_vir_m,viral_mi[,grep("PPF",colnames(viral_mi))])
meta_vir_m=melt(as.matrix(meta_vir_m))
meta_vir_m=na.omit(meta_vir_m)
meta_vir_m$group=group[as.character(meta_vir_m$Var2),2]
meta_vir_m$days=as.factor(group[as.character(meta_vir_m$Var2),3])
meta_vir_m=meta_vir_m[meta_vir_m$value>0,]
ggplot(meta_vir_m,aes(days,value,fill=group))+
  geom_boxplot()+
  theme_bw()+
  labs(y="Microdiversity (pi)")

Acine=read.table("pigvir_acine_phage.txt")
rownames(Acine)=Acine$V1
meta_vir_m$Acine=Acine[as.character(meta_vir_m$Var1),1]
for (i in 1:nrow(meta_vir_m)) {
  if(is.na(meta_vir_m$Acine[i])==TRUE){
    meta_vir_m$Acinetobacter[i]="Non_Acinetobacter"
  }else{
    meta_vir_m$Acinetobacter[i]="Acinetobacter"
  }
}
ggplot(meta_vir_m,aes(days,value,fill=Acinetobacter))+
  geom_boxplot()+
  theme_bw()+
  labs(y="Microdiversity (pi)")

Acine_abu=read.table("../pigvir/pigvir_acine_matrix.tsv",sep='\t',check.names = F)
Acine_abu=melt(as.matrix(Acine_abu))

meta_vir_m_merg=merge(meta_vir_m,Acine_abu[,c("Var1","Var2","value")],by=c("Var1","Var2"),
                      all.x = T)

ggplot(meta_vir_m[grep("3",meta_vir_m$days),],
       aes(Acinetobacter,value,fill=Acinetobacter))+
  geom_boxplot()+
  theme_bw()+
  labs(y="Microdiversity (pi)")

meta_vir_m_merg=na.omit(meta_vir_m_merg)
ggplot(meta_vir_m_merg,
       aes(value.y,value.x,fill=days))+
  geom_point(shape=21)+
  geom_line(aes(group=Acine),alpha=0.7)+
  theme_bw()+
  labs(y="Microdiversity (pi)")


soil_vir_m=read.table("all_soilviruses_metagenome_scaffold_nul.tsv",
                      row.names = 1,sep='\t',header = 1,check.names = F)
data_average=matrix(data = NA,nrow = 10,ncol = ncol(soil_vir_m))
for (i in 1:ncol(data_average)) {
  tmp=as.data.frame(soil_vir_m[soil_vir_m[,i]>0,i])
  tmp=na.omit(tmp)
  for(j in 1:10){
    samp=sample(tmp[,1],size = 500,replace = F)
    tmp2=mean(samp)
    data_average[j,i]=tmp2
  }
}
colnames(data_average)=colnames(soil_vir_m)

soil_vir_m=melt(as.matrix(data_average))
soil_vir_m$days=as.factor(group[as.character(soil_vir_m$Var2),3])
soil_vir_m$group=as.factor(group[as.character(soil_vir_m$Var2),2])
soil_vir_m=soil_vir_m[soil_vir_m$value>0,]
soil_vir_m=na.omit(soil_vir_m)

ggplot(soil_vir_m,aes(days,value,fill=group))+
  geom_boxplot()+
  theme_bw()+
  labs(y="Microdiversity (pi)")
ggplot(soil_vir_m[grep("56",soil_vir_m$days),],
       aes(group,value,fill=group))+
  geom_boxplot(show.legend = F)+
  theme_bw()+
  labs(y="Microdiversity (pi)")+
  stat_compare_means(comparisons = list(c("CK","F"),
                                        c("CK","T")))
ggsave("soilvir56_mi.pdf",device = 'pdf',width = 2,height = 4)

soil_vir_m_merge=merge(soil_vir_m,total[,c("variable","value")],
                       by.x = 'Var2',by.y = 'variable',all.x = T)

model=lm(formula = value.x~value.y,soil_vir_m_merge[grep("T",soil_vir_m_merge$group),])
r_squared <- summary(model)$r.squared
p_value <- summary(model)$coefficients[2, "Pr(>|t|)"]
ggplot(soil_vir_m_merge,aes(value.x,value.y))+
  geom_point(aes(fill=days),shape=21,alpha=0.5,size=0.8)+
  geom_smooth(aes(group=group,color=group),method='lm')+
  theme_bw()+
  labs(x="Average microdiversity of\n phage communities",
       y="Average microdiversity of\n microbial communities")+
  annotate("text", x = max(soil_vir_m_merge$value.x), 
           y = min(soil_vir_m_merge$value.y), 
           label = paste("R-squared:", round(r_squared, 3)), 
           vjust = -1, hjust = 1.1) +
  annotate("text", x = max(soil_vir_m_merge$value.x), 
           y = min(soil_vir_m_merge$value.y) - 0.05 * diff(range(soil_vir_m_merge$value.y)), 
           label = paste("P-value:", p_value),
           vjust = -1, hjust = 1.1)+
  scale_fill_aaas()


soil_vir_m$source="soil"
meta_vir_m$source="pig_manure"
mergered_mi=rbind(soil_vir_m,meta_vir_m)
library(dplyr)
mergered_mi_filtered=mergered_mi %>%
  filter(!source %in% outlet)
library(ggsci)
ggplot(mergered_mi,
       aes(source,value,fill=group))+
  geom_boxplot(notch = TRUE, notchwidth = 0.5,outlier.shape = NA,alpha=0.8)+
  theme_bw()+
  labs(y="Microdiversity (pi)")+ylim(0,0.03)+
  stat_compare_means(label.y = 0.027)+
  scale_fill_d3()


##soil_line
soil_vir_m=read.table("all_soilviruses_metagenome_scaffold_nul.tsv",
                      row.names = 1,sep='\t',header = 1,check.names = F)
soil_vir_m=as.data.frame(t(soil_vir_m))
soil_vir_m$group=group[rownames(soil_vir_m),4]
soil_vir_m=aggregate.data.frame(soil_vir_m[,1:ncol(soil_vir_m)-1],
                                list(soil_vir_m$group),median)
soil_vir_m=melt(soil_vir_m)
group2=group[!duplicated(group$Group),]
rownames(group2)=group2$Group
soil_vir_m$days=as.factor(group2[as.character(soil_vir_m$Group.1),3])
soil_vir_m$group=as.factor(group2[as.character(soil_vir_m$Group.1),2])
soil_vir_m=soil_vir_m[soil_vir_m$value>0,]
ssoil_vir_m=na.omit(soil_vir_m)
ggplot(soil_vir_m[grep("CK",soil_vir_m$Group.1),],
       aes(days,value,group=variable,fill=variable))+
  geom_line(stat = 'identity')+
  theme_bw()+
  labs(y="Microdiversity (pi)")


