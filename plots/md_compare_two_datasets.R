

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)


filename1 <- args[1]
filename2 <- args[2]
out <- args[3]
titl <- args[4]
name1 <- args[5]
name2 <- args[6]

file1 <- read.csv(filename1,header=TRUE,sep=" ")
file2 <- read.csv(filename2,header=TRUE,sep=" ")

len <- ifelse(file1$Length>=100,"Length>=100","Length<100")



# One method per graph
df <- data.frame( id=file1$PDB_ID, len=len, rmsd1=file1$RMSD, rmsd2=file2$RMSD, tmscore1=file1$TMscore, tmscore2=file2$TMscore )
g.rmsd <- ggplot(df,aes(x=rmsd1,y=rmsd2,color=len)) + geom_point() + geom_abline(intercept =0 , slope = 1,linetype="dotted") +
            xlim(3,16) + ylim(3,16) + labs(color="") + labs(title=titl) + xlab(paste("RMSD (",name1,")",sep="")) + ylab(paste("RMSD (",name2,")",sep="")) +
            theme(legend.position="bottom")
            
g.tmscore <- ggplot(df,aes(x=tmscore1,y=tmscore2,color=len)) + geom_point() + geom_abline(intercept =0 , slope = 1,linetype="dotted") +
            xlim(0.2,0.8) + ylim(0.2,0.8) + labs(color="") + labs(title=titl) + xlab(paste("TM-score (",name1,")",sep="")) + ylab(paste("TM-score (",name2,")",sep="")) +
            theme(legend.position="bottom")
            


# In facet wrap: three methods    
deepcov1 <- file1[which(file1$Method=="DeepCov"),]       
deepcov2 <- file2[which(file2$Method=="DeepCov"),]
deepconpred1 <- file1[which(file1$Method=="DeepConPred2"),]
deepconpred2 <- file2[which(file2$Method=="DeepConPred2"),]
deepcontact1 <- file1[which(file1$Method=="DeepContact"),]
deepcontact2 <- file2[which(file2$Method=="DeepContact"),]
len <- ifelse(deepcov1$Length>=100,"Length>=100","Length<100")
df <- data.frame(id=rep(unique(file1$PDB_ID),3),
                 len=rep(len,3),
                 rmsd1=c(deepcov1$RMSD,deepconpred1$RMSD,deepcontact1$RMSD),
                 rmsd2=c(deepcov2$RMSD,deepconpred2$RMSD,deepcontact2$RMSD),
                 tmscore1=c(deepcov1$TMscore,deepconpred1$TMscore,deepcontact1$TMscore),
                 tmscore2=c(deepcov2$TMscore,deepconpred2$TMscore,deepcontact2$TMscore),
                 method=c( rep("DeepCov",nrow(deepcov1)), rep("DeepConPred2",nrow(deepconpred1)), rep("DeepContact",nrow(deepcontact1)) ) )
df$method <- factor(df$method,levels=c("DeepCov","DeepConPred2","DeepContact"))
g.rmsd <- ggplot(df,aes(x=rmsd1,y=rmsd2,color=len)) + geom_point() + geom_abline(intercept =0 , slope = 1,linetype="dotted") + facet_wrap(~method) +
            xlim(3,16) + ylim(3,16) + labs(color="") + xlab(paste("RMSD (",name1,")",sep="")) + ylab(paste("RMSD (",name2,")",sep="")) +
            theme(legend.position="bottom")
g.tmscore <- ggplot(df,aes(x=tmscore1,y=tmscore2,color=len)) + geom_point() + geom_abline(intercept =0 , slope = 1,linetype="dotted") + facet_wrap(~method) +
            xlim(0.2,0.8) + ylim(0.2,0.8) + labs(color="") + xlab(paste("TM-score (",name1,")",sep="")) + ylab(paste("TM-score (",name2,")",sep="")) +
            theme(legend.position="bottom")
