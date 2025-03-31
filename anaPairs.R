#This code performs a TransPairs analysis

source('likelihoodSEIR.R')

#Plot the likelihood
mat=matrix(NA,30,30)
for (i in 1:30) for (j in 1:30) {
  mat[i,j]=likelihoodSEIR(0,i/365,j/365)
}
lattice::levelplot(mat,xlab='Sampling date of infector (days since MRCA)',ylab='Sampling date of infectee (days since MRCA)')

#Perform analysis for a given tree
library(ape)
t=rtree(10,br=function(x) runif(x,min=0,max=0.01))#replace this line with your own tree
n=Ntip(t)
m=mrca(t)
d=dist.nodes(t)
r=matrix(NA,n,n)
rownames(r)<-colnames(r)<-t$tip.label
for (i in 1:n) for (j in 1:n) {
  r[i,j]=likelihoodSEIR(0,d[m[i,j],i],d[m[i,j],j])
}
r=r/max(r)#normalise likelihoods
lattice::levelplot(r,xlab='Infector',ylab='Infectee')

#Compute transmission tree using Edmunds algorithm
library(graph)
library(RBGL)
nel.graph=as(r, Class="graphNEL")
search=edmondsOptimumBranching(nel.graph)
ttree=data.frame(infector=search$edgeList[1,],infectee=search$edgeList[2,],prob=t(search$weights))

