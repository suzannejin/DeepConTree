#!/usr/bin/env Rscript
library('phangorn')
args = commandArgs(trailingOnly=TRUE)

L=read.table(args[1])
list=L[,1]
dis=args[2]

for (fam in list){
  if (dis=='d1'){
    tree3d=read.tree(paste(fam,"/",fam,"_d1_ratio_tmalign1-2.mat_1.txt.nwk", sep=""))
    bootrees3d=read.tree(paste(fam,"/boot3d-d1/", fam, "_bootrees", sep=""))
    } else if (dis=='d4'){
    tree3d=read.tree(paste(fam,"/",fam,"_d4_tmalign4-2.mat_1.txt.nwk", sep=""))
    bootrees3d=read.tree(paste(fam,"/boot3d-d4/", fam, "_bootrees", sep=""))
    } else if (dis=='d5'){
    tree3d=read.tree(paste(fam,"/",fam,"_d5_tmalign5-2.mat_1.txt.nwk", sep=""))
    bootrees3d=read.tree(paste(fam,"/boot3d-d5/", fam, "_bootrees", sep=""))
    } else if (dis=='unw_d1'){
    tree3d=read.tree(paste(fam,"/",fam,"_unw_d1_ratio_tmalign1-2.mat_1.txt.nwk", sep=""))
    bootrees3d=read.tree(paste(fam,"/boot3d-d1-unw/", fam, "_bootrees", sep=""))
    } else if (dis=='unw_d4'){
    tree3d=read.tree(paste(fam,"/",fam,"_unw_d4_tmalign4-2.mat_1.txt.nwk", sep=""))
    bootrees3d=read.tree(paste(fam,"/boot3d-d4-unw/", fam, "_bootrees", sep=""))
    } else if (dis=='unw_d5'){
    tree3d=read.tree(paste(fam,"/",fam,"_unw_d5_tmalign5-2.mat_1.txt.nwk", sep=""))
    bootrees3d=read.tree(paste(fam,"/boot3d-d5-unw/", fam, "_bootrees", sep=""))
    }
  #tree3d=read.tree(paste(fam,"/",fam,"_d4-no_BS_SQRLEN_tmalign4-2.mat_1.txt.nwk", sep=""))
  #bootrees3d=read.tree(paste(fam,"/boot3d-d4/", fam, "_bootrees", sep="")) 
  b3d=prop.clades(tree3d, bootrees3d)
  b3d[is.na(b3d)]<-0
  tree3d$node.label=b3d
  if (dis=='d1'){
    write.tree(tree3d, file=paste(fam,"/",fam,"_d1_ratio_tmalign1-2.mat_1.txt.nwk", sep=""))
    } else if (dis=='d4'){
    write.tree(tree3d, file=paste(fam,"/",fam,"_d4_tmalign4-2.mat_1.txt.nwk", sep=""))
    } else if (dis=='d5'){
    write.tree(tree3d, file=paste(fam,"/",fam,"_d5_tmalign5-2.mat_1.txt.nwk", sep=""))
    } else if (dis=='unw_d1'){
    write.tree(tree3d, file=paste(fam,"/",fam,"_unw_d1_ratio_tmalign1-2.mat_1.txt.nwk", sep=""))
    } else if (dis=='unw_d4'){
    write.tree(tree3d, file=paste(fam,"/",fam,"_unw_d4_tmalign4-2.mat_1.txt.nwk", sep=""))
    } else if (dis=='unw_d5'){
    write.tree(tree3d, file=paste(fam,"/",fam,"_unw_d5_tmalign5-2.mat_1.txt.nwk", sep=""))
    }
  #write.tree(tree3d, file=paste(fam,"/",fam,"_d4-no_BS_SQRLEN_tmalign4-2.mat_1.txt.nwk", sep=""))
  }

#Rscript bootrees.R list
