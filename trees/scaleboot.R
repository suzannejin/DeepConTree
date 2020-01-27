library(scaleboot)

args=commandArgs(trailingOnly=TRUE)

# Perform test
mat=read.mt(args[1])
rell=relltest(mat)
summary(rell)

