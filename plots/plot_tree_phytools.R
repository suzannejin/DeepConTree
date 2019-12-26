library(phangorn)
library(geiger)
library(phytools)

args = commandArgs(trailingOnly=TRUE)
filename <- args[1]

tree <- read.tree(filename)

plotTree <- function(tree,color){
    plot(tree,use.edge.length=FALSE,cex=0.75,font=1)
    nodelabels(tree$node.label,cex=0.6,bg=color)
}


plotTreelabels <- function(tree,tree2,tree3,color1,color2,color3){
    boot_2 <- prop.clades(tree,tree2)
    boot_3 <- prop.clades(tree,tree3)
    
    plot(tree,use.edge.length=FALSE,cex=0.75,font=1,label.offset=1)
    nodelabels(boot_2,cex=0.5,font=2,adj=c(-1.8,-0.8),bg=color2)
    nodelabels(boot_3,cex=0.5,font=2,adj=c(-1.8,1.2),bg=color3)
    nodelabels(tree$node.label,cex=0.6,bg=color1)
}



plotTreelabels2 <- function(tree,tree2,tree3,color1,color2,color3){

    plot(tree,use.edge.length=FALSE,cex=0.75,font=1,label.offset=1)
    nodelabels(tree2$node.label,cex=0.5,font=2,adj=c(-1.2,-0.8),bg=color2)
    nodelabels(tree3$node.label,cex=0.5,font=2,adj=c(-1.2,1.2),bg=color3)
    nodelabels(tree$node.label,cex=0.6,bg=color1)
}


plotTreelabels3 <- function(tree,tree2,tree3,tree4,color2,color3,color4){

    plot(tree,use.edge.length=FALSE,cex=0.75,font=1,label.offset=0.25)
    nodelabels(prop.clades(tree,tree2,rooted=TRUE),cex=0.5,font=2,adj=c(-1.2,-1.5),bg=color2)
    nodelabels(prop.clades(tree,tree3,rooted=TRUE),cex=0.5,font=2,adj=c(-1.2,0.5),bg=color3)
    nodelabels(prop.clades(tree,tree4,rooted=TRUE),cex=0.5,font=2,adj=c(-1.2,2.5),bg=color4)
}







##
exp <- read.tree(paste(pfam,"/trees/experimental/out/tree3d/",pfam,"_unw_d4_tmalign4-2.mat_1.txt.nwk",sep=""))
exp_1d <- read.tree(paste(pfam,"/trees/experimental/out/tree1d/",pfam,".nwk",sep=""))
nat <- read.tree(paste(pfam,"/trees/native_models250_n1/out/tree3d/",pfam,"_unw_d4_tmalign4-2.mat_1.txt.nwk",sep=""))
nat_1d <- read.tree(paste(pfam,"/trees/native_models250_n1/out/tree1d/",pfam,".nwk",sep=""))
deep <- read.tree(paste(pfam,"/trees/deepconpred35_models250_n1/out/tree3d/",pfam,"_unw_d4_tmalign4-2.mat_1.txt.nwk",sep=""))
deep_1d <- read.tree(paste(pfam,"/trees/deepconpred35_models250_n1/out/tree1d/",pfam,".nwk",sep=""))


RF.dist(exp,deep,normalize=TRUE)
RF.dist(exp,exp_1d,normalize=TRUE)
RF.dist(deep,deep_1d,normalize=TRUE)
RF.dist(exp_1d,deep_1d,normalize=TRUE)




##OMA
tol <- read.tree(dir(path=paste("original/",pfam,sep=""),pattern="*_taxa_tree_decoded.nwk",full.names=TRUE))
exp <- read.tree(dir(path=paste("out/",pfam,"/trees/experimental/out/tree3d",sep=""),pattern="*_species_unw_d4_tmalign4-2.mat_1.txt.nwk",full.names=TRUE))
exp_1d <- read.tree(dir(path=paste("out/",pfam,"/trees/experimental/out/tree1d",sep=""),pattern="*_species.nwk",full.names=TRUE))
deep <- read.tree(dir(path=paste("out/",pfam,"/trees/deepconpred_score0-35_seqsep8_n1/out/tree3d",sep=""),pattern="*_species_unw_d4_tmalign4-2.mat_1.txt.nwk",full.names=TRUE))




# Booster
trees <- list(read.tree(paste(pfam,"/trees/booster_trees/exp3d_exp1d.nwk",sep="")),read.tree(paste(pfam,"/trees/booster_trees/exp3d_deep3d.nwk",sep="")))#,read.tree(paste(pfam,"/trees/booster_trees/exp3d_nat3d.nwk",sep="")))
#for (n in 1:length(trees)){trees[[n]]$node.label <- as.integer(as.numeric(trees[[n]]$node.label)*100)}

for (n in 1:length(trees)){print(mean(as.numeric(trees[[n]]$node.label,na.rm=TRUE),na.rm=TRUE))}