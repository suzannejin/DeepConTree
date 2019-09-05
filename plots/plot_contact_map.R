#
# This is script is used to plot contact maps, given a contactbench file with 5 columns { res i, res j, distance in sequence, probability, label 0 or 1 } 
# The triangle at top is the predicted contact map, and the triangle at bottom is the native contact map.
#

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

filename <- args[1]  # Contactbench filename
out <- args[2]  # Output plot filename .png
id <- args[3]  # Protein id: plot title
plotname <- args[4]  # Plot subtitle


file <- read.csv(filename,header=TRUE,sep=" ")

df1 <- data.frame(i=file[,1],j=file[,2],probab=as.numeric(file[,4]))
df2 <- data.frame(i=file[,2],j=file[,1],probab=as.numeric(file[,5]))
df <- rbind(df1,df2)


cols <- rev(rainbow(7)[-7])

g <- ggplot(df,aes(i,j,fill = probab)) + xlab("Residue") + ylab("Residue") + labs(fill="") + labs(title = id, subtitle = plotname)+
  geom_raster() +
  scale_fill_gradientn(colours = cols,limits=c(0,1.0))+
  theme(
      plot.title = element_text(size=18,face="bold"),
      plot.subtitle = element_text(size=18),
      legend.text=element_text(size=18),
      axis.title=element_text(size=20,face="bold"),
      axis.text=element_text(size=18),
      legend.key.height = unit(2, "cm"))

  
ggsave(filename=out,g)




#for folder in $(echo "DEEPCOV DEEPCONPRED DEEPCONTACT/original"); do for id in $(echo "1AZBA-1 1BQKA-1 1CUOA-1 1DYZA-1 1IUZA-1 1JOIA-1 1KDIA-1 1MDAA-1 1PCSA-1 1PMYA-1 1QHQA-1 1RKRA-1 1VLXA-1 2AANA-1 2GIMA-1 2H3XC-1 2PLTA-1 3AY2A-1 3C75A-1 3EF4A-1 4BWUA-1 7PCYA-1"); do dir=plots/contact_map/PF00127/${id}; [[ ! -d $dir ]] && mkdir $dir; awk '$0~/[0-9]+ [0-9]+ [0-9]+/{if($5=="TRUE"){label=1}else{label=0};printf "%s %s %s %s %s\n",$1,$2,$3,$4,label}' BENCHFAM_28.0_2018/first15/${folder}/contactbench/cutoff8/PF00127/${id}.contactbench > a; [[ $folder == "DEEPCOV" ]] && method1=deepcov && method=DeepCov; [[ $folder == "DEEPCONPRED" ]] && method1=deepconpred && method2=DeepConPred2; [[ $folder == "DEEPCONTACT/original" ]] && method1=deepcontact && method2=DeepContact; Rscript BENCHFAM_28.0_2018/first15/scripts/plots/plot_contacts.R a ${dir}/${id}_${method1}_native.png "$id" "$method2 vs Native Contacts"; done; done