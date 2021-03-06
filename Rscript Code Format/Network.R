library(igraph)
library(RColorBrewer)
setwd("E:/学习/研究/protein corona/201908/factor.selection.10/class/function.prediction.network")

data1<-read.csv("proximity-class.a.csv")
data<-as.matrix(data1[, 2:568])
diag(data) <- 0
data[data < 4*mean(data)] <- 0
data[data >= 4*mean(data)] <- 1
net<-graph.adjacency(adjmatrix=data,mode="undirected")
set.seed(12345)
#color.list <- data1[, 180:185]
graph.density(net)

#------------------------------
#RPA.average=vertex.size
#------------------------------
#cores
plot(net,layout=  layout_with_kk,
     vertex.size=6,
     vertex.color=as.character(data1$ver.core),
     vertex.label="",
     vertex.label.size=0.5,
     vertex.label.cex=1,
     vertex.label.dist=0,
     vertex.label.color="black")

legend("topleft",c("ADM105","Ag","Au","Ca3(PO4)2", "CaCO3", "Liposome", "Fe3O4", 
                   "LPD", "CNT", "PS", "PSOSO3", "Si related", "TiO2", "zeolite"),
       pch=16, cex=1,
       col=c("#BF0B0B", "#FF3E3E", "#FF9F9F", "#D17777", "#B21C1C", "#FFFF28", "#183899", 
             "#4C6FD6", "#A6B7EB", "#27805E", "#4EFFBD", "#20CCCC","#FFAE0C", "#FFDF9E"))

#mo
plot(net,layout=  layout_with_kk,
     vertex.size=6,
     vertex.color=as.character(data1$mo),
     vertex.label="",
     vertex.label.size=0.5,
     vertex.label.cex=1,
     vertex.label.dist=0,
     vertex.label.color="black")

legend("topleft",c("none","NH2","COOH","PEG", "CIT", "Amino acid", "DDT", 
                   "EMT", "FAU", "PVP", "NT", "LA", "BPEI", "PVA", "Others"),
       pch=16, cex=1,
       col=c("#3098BF", "#FF8071", "#33B18F", "#7DC7FF", "#FFBFB8", "#99D8C7", "#A8D9FF", 
             "#26A1FF", "#FFBFB8", "#007756", "#66C5AB", "#FFDFDB","#D4ECFF", "#CCECE3", "#FF604D"))






ver.c <- data1$logp
ver.c[ver.c==1] <- "#43F4DB"
ver.c[ver.c==2] <- "#619EFF"
ver.c[ver.c==3] <- "#0CA792"
ver.c[ver.c==4] <- "#F99A85"
ver.c[ver.c==5] <- "#F46D43"

ver.c.type <- data1$core.type2
ver.c.type[ver.c.type==1] <- "#43F4DB"
ver.c.type[ver.c.type==2] <- "#619EFF"
ver.c.type[ver.c.type==3] <- "#0CA792"
ver.c.type[ver.c.type==4] <- "#F99A85"

#ver.m; the colors of modification types
#ver.core; the colors of cores
#pdf("total-proximity.pdf",family="GB1")

#ver.logp
plot(net,layout=  layout_with_kk,
     vertex.size=6,
     vertex.color=ver.c,
     vertex.label="",
     vertex.label.size=0.5,
     vertex.label.cex=1,
     vertex.label.dist=0,
     vertex.label.color="black")

legend("topleft",c("[-5, -1]","(-1, -0.2]","(-0.2, 0]","(0, 5]","(5, 10]"),
       pch=16, cex=1,
       col=c("#43F4DB","#619EFF","#0CA792","#F99A85","#F46D43"))

#ver.mofication.charge
plot(net,layout=  layout_with_kk,
     vertex.size=6,
     vertex.color=data1$ver.m,
     vertex.label="",
     vertex.label.size=0.5,
     vertex.label.cex=1,
     vertex.label.dist=0,
     vertex.label.color="black")

legend("topleft",c("Cationic","Neutral","Antionic"),
       pch=16, cex=1,
       col=c("#F46D43","#619EFF","#0CA792"))

#core.type
plot(net,layout=  layout_with_kk,
     vertex.size=6,
     vertex.color=ver.c.type,
     vertex.label="",
     vertex.label.size=0.5,
     vertex.label.cex=1,
     vertex.label.dist=0,
     vertex.label.color="black")

legend("topleft",c("carbon","liposome","metal","other"),
       pch=16, cex=1,
       col=c("#43F4DB","#619EFF","#0CA792","#F99A85"))




legend("topleft",c("Albumin","Apolipoprotein","Clusterin","Coagulation","Complement", "Immunity","Others"),
       pch=16, cex=1,
       col=c("#43F4DB","#0CA792","#619EFF","#F99A85","#F46D43","#48C192","#A73E1D"))
legend("bottomright",c("<0.2","[0.2, 0.5]","[0.5, 1.0]","[1.0, 1.5]","[1.5, 2.0]","2.0<"),pch=1, cex=c(1, 1, 1, 1, 1, 1))
#dev.off()


##subgroup
##pdf("similar-subgroup-1.pdf",family="GB1")
com<-walktrap.community(net,steps=6)
V(net)$sg=com$membership+1
V(net)$color=rainbow(max(V(net)$sg))[V(net)$sg]
par(mar=c(0,0,0,0))
set.seed(12345)
plot(net,layout=layout_with_kk	,vertex.size=5,
     vertex.color=V(net)$color,vertex.label=NA, edge.color=grey(0.5),
     edge.arrow.mode="-")
membership(com)
plot(com,net,layout=layout_with_kk,vertex.size=5,
     vertex.color=V(net)$color,vertex.label=NA, edge.color=grey(0.5),
     edge.arrow.mode="-")
##dev.off()

##dendPlot(com)  ##将subgroup节点融合过程以树形图展示出
##modularity(com,com[[1]])  ##计算各模块的模块度Q值

##subgroup
sub1<-induced.subgraph(net,V(net)[sg==2])
graph.density(sub1)
graph.density(net)

png(file="similarnetwork1.png",width=1080,height=1080,bg="white",res=160)
plot(net,layout=layout.fruchterman.reingold,
     vertex.size=0.2,
     vertex.label=data1$name,
     vertex.label.cex=0.5,
     vertex.label.dist=0,
     vertex.label.color="black")
dev.off()



#------------------------------
#RPA.average=vertex.size
#------------------------------
ver.c <- color.list$color.r2,
ver.c[ver.c==1] <- "#ABDDA4"
ver.c[ver.c==2] <- "#FEE08B"
ver.c[ver.c==3] <- "#FDAE61"
ver.c[ver.c==4] <- "#F46D43"
ver.c[ver.c==5] <- "#D53E4F"
#-------------------------------
#ver.function
#-------------------------------
ver.f <- color.list$color.function
ver.f[ver.f==1] <- "#43F4DB"
ver.f[ver.f==2] <- "#0CA792"
ver.f[ver.f==3] <- "#619EFF"
ver.f[ver.f==4] <- "#F99A85"
ver.f[ver.f==5] <- "#F46D43"
ver.f[ver.f==6] <- "#48C192"
ver.f[ver.f==7] <- "#A73E1D"

