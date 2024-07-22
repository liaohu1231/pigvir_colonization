type=read.table("all.normalize_16s.type.tab.txt",sep = '\t',header = 1,row.names = 1,check.names = F)
library(ggplot2)
library(reshape2)
metadata=read.table("../metadata.txt",sep='\t',header = 1,row.names = 1)
type=as.data.frame(t(type))
type=type[order(match(rownames(type),rownames(metadata))),]
average_type=rowsum(x = type,group = as.factor(metadata$Group))
average_type=average_type/3
average_type=average_type[order(match(rownames(average_type),metadata$Group)),]
library(pheatmap)
average_type=t(average_type)
pheatmap(average_type, color = colorRampPalette(c("white","black"))(100),
         border_color = "gray",cluster_rows = F,cluster_cols = F)


###human related args
human_related=read.table("all.human_related_ARGs.count",sep = '\t',header = T,row.names = 1,check.names = F)
human_related_sum=data.frame(sum=apply(human_related,2,sum))
human_related_sum[,2:3]=metadata[rownames(human_related_sum),1:2]
hmd=read.csv("arg_name_uniq_hmd.blastp",sep = '\t',header = F)
hmd2=hmd[!duplicated(hmd$V1),]
rownames(hmd2)=hmd2$V1
human_related$annotation=hmd2[as.character(rownames(human_related)),2]
tmp=strsplit(human_related$annotation,split = '\\|')
tmp=do.call(rbind,tmp)
human_related$class=tmp[,4]
human_related_class=matrix(data = NA,nrow = nrow(human_related[!duplicated(human_related$class),]),
                           ncol = 48)
for(i in 1:48){
  human_related_class[,i]=tapply(human_related[,i],human_related$class,sum)
}
rownames(human_related_class)=rownames(as.data.frame(tapply(human_related[,1],
                                                            human_related$class,sum)))
colnames(human_related_class)=colnames(human_related)[1:48]
library(reshape2)
human_related_class=as.data.frame(t(human_related_class))
human_related_class$group=metadata[rownames(human_related_class),3]
human_related_class_g=aggregate(human_related_class[,1:10],by = list(human_related_class$group),mean)
human_related_class_g=melt(human_related_class_g)
metadata_2=metadata[!duplicated(metadata$Group),]
rownames(metadata_2)=metadata_2$Group
human_related_class_g[,4:6]=metadata_2[as.character(human_related_class_g$Group.1),1:3]
human_related_class_g$Days=factor(human_related_class_g$Days,levels = c("0","1","3","7","14","56"))
human_related_class_g$Microcosm=factor(human_related_class_g$Microcosm,levels = c("CK","T","F"))
library(ggsci)
library(ggplot2)
ggplot(human_related_class_g,aes(x=Days,y = value,fill=variable))+
  geom_bar(stat = 'identity')+
  facet_wrap(~Microcosm,scales = 'free_x')+
  theme_bw()+scale_fill_aaas()

###human_related_sankey
human_sankey=data.frame(row.names = rownames(human_related),
                        args_class=human_related$class)
plasmid=read.csv("arg_name_uniq_human_related_id99.fasta.plasflowout",sep='\t',header = 1)
rownames(plasmid)=plasmid$contig_name
for (i in 1:nrow(human_sankey)){
  tmp=strsplit(rownames(human_sankey)[i],split = '_')
  tmp=as.data.frame(do.call(rbind,tmp))
  tmp=paste(tmp[,1:ncol(tmp)-1],collapse = '_')
  human_sankey$contigs[i]=tmp
}
human_sankey$plasflow=plasmid[as.character(human_sankey$contigs),6]
phage=read.csv("arg_name_uniq_carrier_phage.txt",header = F)
phage$phage="True"
rownames(phage)=phage$V1
human_sankey$phage=phage[as.character(human_sankey$contigs),2]
phage_host=read.csv("../all.soil.vira_meta_95-85.host.csv",header = F)
phage_host=phage_host[!duplicated(phage_host$V1),]
rownames(phage_host)=phage_host$V1
human_sankey$phage_host=phage_host[human_sankey$contigs,8]
plasmid_host=read.csv("arg_plasimd_Host_prediction_to_genus_m90.csv",check.names = F)
plasmid_host=plasmid_host[!duplicated(plasmid_host$Virus),]
rownames(plasmid_host)=plasmid_host$Virus
tmp=strsplit(plasmid_host$Host_genus,split = ";")
tmp=do.call(rbind,tmp)
plasmid_host$genus=tmp[,6]
human_sankey$plasmid_host=plasmid_host[human_sankey$contigs,6]
tmp=strsplit(human_sankey$plasflow,split = '\\.')
tmp=do.call(rbind,tmp)
human_sankey$plasflow=tmp[,1]
for(i in 1:nrow(human_sankey)){
  if(human_sankey$plasflow[i]=="chromosome" && is.na(human_sankey$phage[i])==TRUE){
    human_sankey$type[i]="chromosome"
  }
  if(human_sankey$plasflow[i]=="plasmid" && is.na(human_sankey$phage[i])==TRUE){
    human_sankey$type[i]="plasmid"
  }
  if(human_sankey$plasflow[i]=="plasmid" && is.na(human_sankey$phage[i])==FALSE){
    human_sankey$type[i]="phage-plasmids"
  }
  if(human_sankey$plasflow[i]=="unclassified" && is.na(human_sankey$phage[i])==FALSE){
    human_sankey$type[i]="phage"
  }
}
human_sankey_contig=human_sankey[!duplicated(human_sankey$contigs),]
table(human_sankey_contig$type)
pathogen=read.csv("arg_name_uniq_human_related_id99_pathogen_id99.blastx",sep='\t',header = F,row.names = 1)
human_sankey$pathogen=pathogen[human_sankey$contigs,1]
for(i in 1:nrow(human_sankey)){
  if(is.na(human_sankey$phage_host)[i]==FALSE){
    human_sankey$host[i]=human_sankey$phage_host[i]
  }
  if(is.na(human_sankey$phage_host)[i]==TRUE){
    human_sankey$host[i]=human_sankey$plasmid_host[i]
  }
}


write.table(human_sankey,file = "Supplementary_table_S1.csv",quote = F,sep='\t',row.names = T)

human_related$classification=human_sankey[rownames(human_related),7]
human_related_clssification=as.data.frame(matrix(data = NA,
                                                 nrow = nrow(human_related[!duplicated(human_related$classification),]),
                           ncol = 48))
for(i in 1:48){
  human_related_clssification[,i]=tapply(human_related[,i],human_related$classification,sum)
}
rownames(human_related_clssification)=rownames(as.data.frame(tapply(human_related[,1],
                                                            human_related$classification,sum)))
colnames(human_related_clssification)=colnames(human_related)[1:48]
library(reshape2)
human_related_class=as.data.frame(t(human_related_clssification))
human_related_class$group=metadata[rownames(human_related_class),3]
human_related_class_g=aggregate(human_related_class[,1:4],
                                by = list(human_related_class$group),mean)
human_related_class_g=melt(human_related_class_g)
metadata_2=metadata[!duplicated(metadata$Group),]
rownames(metadata_2)=metadata_2$Group
human_related_class_g[,4:6]=metadata_2[as.character(human_related_class_g$Group.1),1:3]
human_related_class_g$Days=factor(human_related_class_g$Days,levels = c("0","1","3","7","14","56"))
human_related_class_g$Microcosm=factor(human_related_class_g$Microcosm,levels = c("CK","T","F"))
library(ggsci)
library(ggplot2)
ggplot(human_related_class_g,aes(x=Days,y = value,fill=variable))+
  geom_bar(stat = 'identity')+
  facet_wrap(~Microcosm,scales = 'free_x')+
  theme_bw()+scale_fill_aaas()

##host
human_related$classification=human_sankey[rownames(human_related),9]
human_related_2=na.omit(human_related)
human_related_clssification=as.data.frame(matrix(data = NA,
                                                 nrow = nrow(human_related_2[!duplicated(human_related_2$classification),]),
                                                 ncol = 48))
for(i in 1:48){
  human_related_clssification[,i]=tapply(human_related_2[,i],human_related_2$classification,sum)
}
rownames(human_related_clssification)=rownames(as.data.frame(tapply(human_related_2[,1],
                                                                    human_related_2$classification,sum)))
colnames(human_related_clssification)=colnames(human_related_2)[1:48]
library(reshape2)
human_related_class=as.data.frame(t(human_related_clssification))
human_related_class$group=metadata[rownames(human_related_class),3]
human_related_class_g=aggregate(human_related_class[,1:4],
                                by = list(human_related_class$group),mean)
human_related_class_g=melt(human_related_class_g)
metadata_2=metadata[!duplicated(metadata$Group),]
rownames(metadata_2)=metadata_2$Group
human_related_class_g[,4:6]=metadata_2[as.character(human_related_class_g$Group.1),1:3]
human_related_class_g$Days=factor(human_related_class_g$Days,levels = c("0","1","3","7","14","56"))
human_related_class_g$Microcosm=factor(human_related_class_g$Microcosm,levels = c("CK","T","F"))
human_related_class_g$variable=as.character(human_related_class_g$variable)
for (i in 1:nrow(human_related_class_g)) {
  if(grepl("__",human_related_class_g$variable[i])==TRUE){
    tmp=strsplit(as.character(human_related_class_g$variable[i]),"__")
    tmp=do.call(rbind,tmp)
    human_related_class_g$variable[i]=tmp[,2]
  }else{
    human_related_class_g$variable[i]=human_related_class_g$variable[i]
  }
}
library(ggsci)
library(ggplot2)
ggplot(human_related_class_g,aes(x=Days,y = value,fill=variable))+
  geom_bar(stat = 'identity')+
  facet_wrap(~Microcosm,scales = 'free_x')+
  theme_bw()+scale_fill_aaas()

###viromes_human_related_ARGs
arg_viromes=read.table("all.human_related_viromes.count",sep = '\t',header = 1,
                       row.names = 1,check.names = F)

viromes_metadata=data.frame(row.names = colnames(arg_viromes),
                            days=rep(c("14","56","14","PPF","0","14","56"),
                                     c(3,2,3,3,3,3,2)))
tmp=strsplit(rownames(viromes_metadata),split = '-')
tmp=do.call(rbind,tmp)
viromes_metadata$group=tmp[,1]

arg_viromes_count=data.frame(count=apply(arg_viromes,2,sum))
arg_viromes_count$days=viromes_metadata$days
arg_viromes_count$group=viromes_metadata$group
arg_viromes_count_mean=data.frame(mean=tapply(arg_viromes_count$count,
                                              arg_viromes_count$group,mean))
viromes_metadata2=viromes_metadata[!duplicated(viromes_metadata$group),]
viromes_metadata2$microcosm=c("CK","CK","F","PV","CK","T","T")
rownames(viromes_metadata2)=viromes_metadata2$group
arg_viromes_count_mean$sd=tapply(arg_viromes_count$count,
                                 arg_viromes_count$group,sd)
arg_viromes_count_mean$days=viromes_metadata2[rownames(arg_viromes_count_mean),1]
arg_viromes_count_mean$microcosm=viromes_metadata2[rownames(arg_viromes_count_mean),3]
arg_viromes_count$days=factor(arg_viromes_count$days,
                              levels = as.factor(c("PPF0","0","14","56")))
arg_viromes_count_mean$group=rownames(arg_viromes_count_mean)

library(ggsci)
library(ggplot2)
ggplot(arg_viromes_count_mean,aes(x=days,y = mean,fill=microcosm))+
  geom_bar(stat = 'identity')+
  #facet_wrap(~microcosm,scales = 'free_x')+
  theme_bw()+scale_fill_aaas()+
  labs(y="Mapped reads count of \nhuman-related ARGs-carriers\nin viromes")
ggsave("human-related_ARGs_carriers_viromes.pdf",device = pdf,width = 4,height = 3.5)
