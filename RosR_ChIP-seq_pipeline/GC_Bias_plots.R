### GC Bias plots script  for use with 'RosR_ChIP-seq_pipeline'
# Cynthia L. Darnell, Rylee K. Hackley, and Amy K. Schmid
# https://github.com/amyschmid/rosr-chip-utils/tree/master/RosR_ChIP-seq_pipeline

library(ggplot2)

bins <- as.matrix(seq(0, 100, by=1/3))
log2ratio <- matrix()

files <- list.files("./aligned_data/GC_bias/")[endsWith(list.files("./aligned_data/GC_bias/"), "WCE_sorted_GC.txt")]


for (i in 1:length(files)){
  
  GCBias <- read.delim(paste("./aligned_data/GC_bias/", files[i], sep=""), header=FALSE, sep="")
  
  all_data <- cbind(bins, GCBias, log2ratio)
  colnames(all_data) <- c("GC%", "Observed_reads", "Expected_reads", "Observed/Expected", "Log2ratio")
  
  log2ratio <- (log2(all_data$"Observed/Expected"))
  all_data$Log2ratio <- log2ratio
  
  ggplot(all_data)+
    geom_line(aes(x=all_data$"GC%", y=all_data$"Observed_reads"), color="blue")+
    geom_line(aes(x=all_data$"GC%", y=all_data$"Expected_reads"), color="red", linetype=2)+
    theme_bw()+
    labs(x="%GC", y="Reads per bin")+
    scale_x_continuous(limits=c(35, 80))
  
  
  ggsave(paste("./aligned_data/GC_bias/", files[i], "_reads.png", sep = ""), width=5, height=5)
  
  ggplot(all_data)+
    geom_line(aes(x=all_data$"GC%", y=all_data$"Log2ratio"), color="blue")+
    theme_bw()+
    labs(x="%GC", y="Log2 ratio observed/expected reads per bin")+
    scale_x_continuous(limits=c(35, 80))+
    geom_hline(aes(yintercept=0))
  
  ggsave(paste("./aligned_data/GC_bias/", files[i], "_ratio.png", sep = ""), width=5, height=5)
  
}