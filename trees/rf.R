#!/usr/bin/env Rscript
library('phangorn')
library('geiger')
library('xlsx')


args = commandArgs(trailingOnly=TRUE)


names.filename <- args[1]
trees.filename <- args[2]
out.txt <- args[3]
out.xlsx <- args[4]


file <- read.csv(names.filename,header=FALSE,sep=" ")
tree.names <- file$V2
trees <- read.tree(trees.filename,tree.names=tree.names)

RobFou <- RF.dist(trees,normalize=TRUE)

out <- as.matrix(RobFou)
write.table(data.frame(value=out[,1],name=tree.names),out.txt,quote=FALSE,sep=" ",col.names=FALSE,row.names=FALSE)
write.xlsx(out,out.xlsx)


