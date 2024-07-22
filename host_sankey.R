##host
host_phylum=read.table("host_class_poly_or_single.txt",sep='\t')
host_phylum=as.data.frame(table(host_phylum$V2))
host_phylum$var2=as.factor("total")
library(ggplot2)
ggplot(host_phylum, aes(x=var2, y=Freq,fill=Var1)) +
  geom_bar(stat="identity",width=1,color="black")+
  coord_polar(theta = "y",start=0)+
  theme_void()+scale_fill_d3()
ggsave("polyvalent_in_class.pdf",device = "pdf",width = 4,height = 3.5)

host=read.csv("all.pigviruses_host.csv",sep = ",",header = F)

library(ggalluvial)
library(tidyverse)
library(readr)
library(ggplot2)
library(ggalluvial)
library(latex2exp)
library(ggsci)
data01=read.csv("all.soil.vira_meta_95-85.host_sankey.csv",sep=',',header = F)
data01=data01[data01$V5>200,]
data01$V1=as.factor(data01$V1)
data01$V2=as.factor(data01$V2)
data01$V3=as.factor(data01$V3)
data01$V4=as.factor(data01$V4)
data01=data01[-grep("NAmissing",data01$V4),]
data01$V4=factor(data01$V4,levels = c(data01$V4))
cp2 <- c("#FC7945", "#4AB08C", "#482870","gray")
ggplot(data=data01,
       aes(axis1=V1,axis2=V2,axis3=V3,axis4=V4,y=V5))+
  geom_stratum(width=0.1)+
  geom_alluvium(aes(fill=V1,color=V1),size=0.8,width = 0.1,alpha=0.3,
                show.legend = FALSE,curve_type = "cubic")+
  scale_fill_d3(palette = "category20")+
  scale_color_d3()+
  geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_x_continuous(breaks = c(1,2),
                     labels = c("first","second"),
                     expand = expansion(mult = c(0,0)))+
  theme_void()+
  scale_x_continuous(breaks = 1:4, labels = c("phylum", "class","order","family"), expand = c(0.025,0.025))
  