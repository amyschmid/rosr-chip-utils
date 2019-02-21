### Generate wiggle files script using MOSAiCS R package for use with 'RosR_ChIP-seq_pipeline'
# Cynthia L. Darnell, Rylee K. Hackley, and Amy K. Schmid
# https://github.com/amyschmid/rosr-chip-utils/tree/master/RosR_ChIP-seq_pipeline

# MOSAiCS 
# Dongjun Chung, Pei Fen Kuan, Rene Welch, Sunduz Keles
# https://bioconductor.org/packages/release/bioc/html/mosaics.html


library(mosaics)

# for original bam files

bam_files <- list.files("./aligned_data/")[endsWith(list.files("aligned_data/"), "_sorted.bam")]


for (i in 1:length(bam_files)){
  
  generateWig(infile=paste("./aligned_data/", bam_files[i], sep=""),
              fileFormat="bam", 
              outfileLoc="./aligned_data/wiggle_files/",
              byChr=FALSE,
              PET=FALSE, 
              fragLen=200, 
              span=50, 
              capping=0, 
              normConst=1)
}

# for GC corrected files

GCbam_files <- list.files("./aligned_data/GC_bias/")[endsWith(list.files("./aligned_data/GC_bias/"), "_GCcorrected.bam")]

for (i in 1:length(GCbam_files)){
  
  generateWig(infile=paste("./aligned_data/GC_bias", bam_files[i], sep=""),
              fileFormat="bam", 
              outfileLoc="./aligned_data/wiggle_files/",
              byChr=FALSE,
              PET=FALSE, 
              fragLen=200, 
              span=50, 
              capping=0, 
              normConst=1)
}