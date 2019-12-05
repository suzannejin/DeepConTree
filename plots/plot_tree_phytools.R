library(phangorn)
library(geiger)
library(phytools)

args = commandArgs(trailingOnly=TRUE)
filename <- args[1]

tree <- read.tree(filename)

plotTree <- function(tree){
    plot(tree,use.edge.length=FALSE,cex=0.75,font=1)
    nodelabels(tree$node.label,cex=0.6)
}