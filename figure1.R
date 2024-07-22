data=data.frame(name=as.factor(c("Virues","Non_virues")),
                number=c(62,527))
data$var2=as.factor("total")
library(ggplot2)
ggplot(data, aes(x=var2, y=number,fill=name)) +
  geom_bar(stat="identity",width=1,color="black")+
  coord_polar(theta = "y",start=0)+
  theme_void()
data=data.frame(name=as.factor(c("Virues","Non_virues")),
                number=c(108663,190945-108663))
data$var2=as.factor("total")
library(ggplot2)
ggplot(data, aes(x=var2, y=number,fill=name)) +
  geom_bar(stat="identity",width=1,color="black")+
  coord_polar(theta = "y",start=0)+
  theme_void()
data=read.table("pigviromes_arg_name_uniq_viruses.txt",sep='\t',header = F,row.names = 1)
data=as.data.frame(table(data[,7]))
data$var2=as.factor("total")
library(ggplot2)
library(ggsci)
ggplot(data, aes(x=var2, y=Freq,fill=Var1)) +
  geom_bar(stat="identity",width=1,color="black")+
  coord_polar(theta = "y",start=0)+
  theme_void()+scale_fill_d3()

abu=read.table("pigvirues_viromes.count",header = 1,row.names = 1,sep = '\t',check.names = F)
abu=abu[,-c(9:11)]
number=data.frame(number=apply(abu,2,sum))
number$group=rep(c("CK14","CK56","F14","S0","T14","T56"),c(3,2,3,3,3,2))
number_mean=data.frame(mean=tapply(number$number, number$group, mean),
                       sd=tapply(number$number, number$group, sd))
number_mean$days=as.factor(c(14,56,14,0,14,56))
number_mean$group=as.factor(c("CK","CK","F","S0","T","T"))
dodge <- position_dodge(width=.9)
ggplot(number_mean,aes(days,mean))+
  geom_bar(aes(fill=group),
           stat="identity", position=dodge)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=group),
                stat="identity", position=dodge, width=.3)+
  scale_fill_lancet()+
  scale_color_discrete()+
  theme_bw()
