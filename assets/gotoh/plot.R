library(ggplot2)

t=read.delim("runtime.txt",header=F,sep=" ")                                                                         
names(t)=c("Init","Int","Variant","Runtime")
t$Init=factor(ifelse(t$Init==0,"Init once","Init always"))
t$Int=factor(ifelse(t$Int==0,"Float","Integer"))
t$Variant=as.factor(paste(t$Variant," "))

png('gotoh.png')
ggplot(t,aes(Variant,Runtime,colour=Variant))+geom_point(size=4,shape=4)+facet_grid(Init~Int)+scale_y_continuous(limits=c(0,250))+theme(text=element_text(size=22),legend.position="none")+ylab("Runtime [sec]")
dev.off()

