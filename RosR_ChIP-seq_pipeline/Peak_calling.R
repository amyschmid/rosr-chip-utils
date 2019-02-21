### Peak calling script using MOSAiCS R package for use with 'RosR_ChIP-seq_pipeline'
# Cynthia L. Darnell, Rylee K. Hackley, and Amy K. Schmid
# https://github.com/amyschmid/rosr-chip-utils/tree/master/RosR_ChIP-seq_pipeline

# MOSAiCS 
# Dongjun Chung, Pei Fen Kuan, Rene Welch, Sunduz Keles
# https://bioconductor.org/packages/release/bioc/html/mosaics.html


library(mosaics)
library(hexbin)

# Construct bins for IPs

IP_files <- list.files("./aligned_data/")[endsWith(list.files("aligned_data/"), "IP_sorted.bam")]

for (i in 1:length(IP_files)){
  
  constructBins(infile=paste("./aligned_data/", IP_files[i], sep=""),
                fileFormat="bam",
                outfileLoc="./peak_analysis/MOSAiCS_results/bins/",
                byChr=FALSE,
                fragLen=300,
                binSize=200,
                capping=0,
                PET=FALSE)
}

# Construct bins for WCEs (GC corrected)

WCE_files <- list.files("./aligned_data/GC_bias/")[endsWith(list.files("aligned_data/GC_bias/"), "WCE_GCcorrected.bam")]

for (i in 1:length(WCE_files)){
  
  constructBins(infile=paste("./aligned_data/GC_bias/", WCE_files[i], sep=""),
                fileFormat="bam",
                outfileLoc="./peak_analysis/MOSAiCS_results/bins/",
                byChr=FALSE,
                fragLen=300,
                binSize=200,
                capping=0,
                PET=FALSE)
}

# Matching sample IP and WCE

pyrF1_minus <- c("./peak_analysis/MOSAiCS_results/bins/pyrF1_minus_IP_sorted.bam_fragL300_bin200.txt", "./peak_analysis/MOSAiCS_results/bins/pyrF1_minus_WCE_GCcorrected.bam_fragL300_bin200.txt")
pyrF2_minus <- c("./peak_analysis/MOSAiCS_results/bins/pyrF2_minus_IP_sorted.bam_fragL300_bin200.txt", "./peak_analysis/MOSAiCS_results/bins/pyrF2_minus_WCE_GCcorrected.bam_fragL300_bin200.txt")
trmBHA1_minus <- c("./peak_analysis/MOSAiCS_results/bins/trmBHA1_minus_IP_sorted.bam_fragL300_bin200.txt", "./peak_analysis/MOSAiCS_results/bins/trmBHA1_minus_WCE_GCcorrected.bam_fragL300_bin200.txt")
trmBHA1_plus <- c("./peak_analysis/MOSAiCS_results/bins/trmBHA1_plus_IP_sorted.bam_fragL300_bin200.txt", "./peak_analysis/MOSAiCS_results/bins/trmBHA1_plus_WCE_GCcorrected.bam_fragL300_bin200.txt")
trmBHA2_minus <- c("./peak_analysis/MOSAiCS_results/bins/trmBHA2_minus_IP_sorted.bam_fragL300_bin200.txt", "./peak_analysis/MOSAiCS_results/bins/trmBHA2_minus_WCE_GCcorrected.bam_fragL300_bin200.txt")
trmBHA2_plus <- c("./peak_analysis/MOSAiCS_results/bins/trmBHA2_plus_IP_sorted.bam_fragL300_bin200.txt", "./peak_analysis/MOSAiCS_results/bins/trmBHA2_plus_WCE_GCcorrected.bam_fragL300_bin200.txt")
trmBHA3_minus <- c("./peak_analysis/MOSAiCS_results/bins/trmBHA3_minus_IP_sorted.bam_fragL300_bin200.txt", "./peak_analysis/MOSAiCS_results/bins/trmBHA3_minus_WCE_GCcorrected.bam_fragL300_bin200.txt")
trmBHA4_minus <- c("./peak_analysis/MOSAiCS_results/bins/trmBHA4_minus_IP_sorted.bam_fragL300_bin200.txt", "./peak_analysis/MOSAiCS_results/bins/trmBHA4_minus_WCE_GCcorrected.bam_fragL300_bin200.txt")

# Running peak calling analysis

all_samples <- rbind(pyrF1_minus, pyrF2_minus, trmBHA1_minus, trmBHA1_plus, trmBHA2_minus, trmBHA2_plus, trmBHA3_minus, trmBHA4_minus)

for (i in 1:nrow(all_samples)){
  
  binTest <- readBins(type=c("chip", "input"), fileName=all_samples[i,])
  count_data <- hexbin (binTest@input, binTest@tagCount, xbins=100)
  control <- plot(count_data, trans=log, inv=exp, colramp=rainbow, xlab="WCE", ylab="ChIP", lcex=0.9)
  hexVP.abline(control$plot.vp, a=0, b=sum(binTest@tagCount)/sum(binTest@input), lwd=0.2)
  
  dev.copy(png, paste("./peak_analysis/MOSAiCS_results/plots/", rownames(all_samples)[i], "_counts.png", sep=""))
  dev.off()
  
  fitTest <- mosaicsFit(binTest, analysisType="IO", bgEst="rMOM")
  plot(fitTest)
  
  dev.copy(png, paste("./peak_analysis/MOSAiCS_results/plots/", rownames(all_samples)[i], "_fit.png", sep=""))
  dev.off()
  
  peakTest <- mosaicsPeak(fitTest, signalModel="2S", FDR=0.01)
  export(peakTest, type="bed", filename=paste("./peak_analysis/MOSAiCS_results/output/", rownames(all_samples)[i], "_GC.bed", sep=""))
  export(peakTest, type="txt", filename=paste("./peak_analysis/MOSAiCS_results/output/", rownames(all_samples)[i], "_GC.txt", sep=""))
  
}