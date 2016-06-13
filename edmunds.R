library(graph)
library(RBGL)

rm(list=ls())
test1<-as.matrix(read.table("matrix.csv",sep=","))
test1[test1==-Inf]<-0;test1[test1==0]<-NaN;test2=exp(test1);test2[test2=='NaN']<-0

# name row and columns with numbers
rownames(test2)=seq(1,dim(test2)[1],1)
colnames(test2)=seq(1,dim(test2)[1],1)

# Create a nel.graph object from a matrix
nel.graph <-as(test2, Class="graphNEL")

# Find an optimum using Edmonds algorithm
search <- edmondsOptimumBranching(nel.graph)

#Output transmission tree
mat=cbind(as.numeric(search$edgeList[1,]),as.numeric(search$edgeList[2,]),t(search$weights))
write.table(mat,'transtree.csv',quote=F,row.names = F,col.names = F,sep=',')

