#!/usr/bin/env Rscript
library('phangorn')
library('geiger')
library('xlsx')


args = commandArgs(trailingOnly=TRUE)


names.filename <- args[1]
trees.filename <- args[2]
out.txt <- args[3]
out.xlsx <- args[4]

# Determine RF
file <- read.csv(names.filename,header=FALSE,sep=" ")
tree.names <- file$V2
trees <- read.tree(trees.filename,tree.names=tree.names)
RobFou <- RF.dist(trees,normalize=TRUE)

# Write matrix
out <- as.matrix(RobFou)
write.xlsx(out,out.xlsx)

# Write list
values=c()
names=c()
n=length(tree.names)
for (i in c(1:n)){
    for (j in c(1:n)){
        name1=tree.names[i];name2=tree.names[j]
        name=paste(name1,name2,sep="_")
        rf=round(out[i,j],digits=4)
        values=c(values,rf)
        names=c(names,name)
    }
}
df=data.frame(rf=values,name=names)
write.table(df,out.txt,quote=FALSE,sep=" ",col.names=FALSE,row.names=FALSE)

