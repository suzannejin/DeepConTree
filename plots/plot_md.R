
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

filename <- args[1]   # Input file with the following columns: PDB_ID RMSD TMscore Method
out <- args[2]  # Output plots filename prefix

out.rmsd <- paste(out,"_rmsd.png",sep="")
out.tmscore <- paste(out,"_tmscore.png",sep="")

file <- read.csv(filename,header=TRUE,sep=" ")

len <- ifelse(file$Length>=100,"Length>=100","Length<100")

df <- data.frame(id=file$PDB_ID, rmsd=file$RMSD, tmscore=file$TMscore, method=file$Method,len=len)
df$method <- factor(df$method,levels=c("DeepCov","DeepConPred2","DeepContact","Native"))

g.rmsd <- ggplot(df,aes(x=id,y=rmsd,color=method)) + geom_point(size=2.5) + facet_wrap(~len) +
              xlab("Protein") + ylab("RMSD") + labs(color="") + ylim(0,17) +
              theme(
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                legend.position="bottom",
                plot.title = element_text(size=18,face="bold"),
                legend.text=element_text(size=18),
                axis.title=element_text(size=18,face="bold"),
                axis.text=element_text(size=18)
              )
g.tmscore <- ggplot(df,aes(x=id,y=tmscore,color=method)) + geom_point(size=2.5) +
              xlab("Protein") + ylab("TM-score") + labs(color="") + ylim(0,1) +
              theme(
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                plot.title = element_text(size=18,face="bold"),
                legend.text=element_text(size=18),
                axis.title=element_text(size=18,face="bold"),
                axis.text=element_text(size=18)
              )

ggsave(filename=out.rmsd,g.rmsd,width=18,height=12,units="cm")
ggsave(filename=out.tmscore,g.tmscore,width=18,height=12,units="cm")



# ggplot(df,aes(x=id,y=rmsd,color=method,group=method,order=method)) + geom_line() + xlab("Protein") + ylab("RMSD") + labs(color="") + ylim(0,17) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom")
# ggplot(df,aes(x=id,y=tmscore,color=method,group=method)) + geom_line() + xlab("Protein") + ylab("TM-score") + labs(color="") + ylim(0,1) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom")
# ggplot(df,aes(x=id,y=precision,color=method,group=method,order=method)) + geom_line() + xlab("Protein") + ylab("Precision") + labs(color="") + ylim(0.5,1) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom")