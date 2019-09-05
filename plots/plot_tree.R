#!/usr/bin/env Rscript
library('phangorn')
library('geiger')


args = commandArgs(trailingOnly=TRUE)


filename <- args[1]
outname <- args[2]


tree <- read.tree(filename)


png(outname)
plot(tree,show.node.label=TRUE)
dev.off()

