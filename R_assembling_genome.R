#assembling a draft genome and analyzing structural variants with long-read sequencing technologies
#https://www.sciencedirect.com/science/article/pii/S2666166722003860
#I read through suggestions in paper and re-wrote this to work in R

#-----SETUP-----

#setwd
setwd("C:/Users/drewj/OneDrive/R_based_work/assembling_genome")


#-----get data-----

#setwd to different directory (so we don't store large data files in a github repo)
setwd("C:/Users/drewj/OneDrive/R_based_work/assembling_genome_nonrepo")

#download data to local
x <- 
'curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR130/025/SRR13070625/SRR13070625_1.fastq.gz -o SRR13070625_Nanopore_sequencing_of_Drosophila_melanogaster_whole_adult_flies_pooled_male_and_female_1.fastq.gz'
system(x)

x <- 
'curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/080/SRR12473480/SRR12473480_subreads.fastq.gz -o SRR12473480_Drosophila_PacBio_HiFi_UltraLow_subreads.fastq.gz'
system(x)

x <- 
  'curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/022/SRR12099722/SRR12099722_1.fastq.gz -o SRR12099722_WGS_Drosophila_melanogaster_adult_ISCs_1.fastq.gz'
system(x)

x <- 
  'curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/022/SRR12099722/SRR12099722_2.fastq.gz -o SRR12099722_WGS_Drosophila_melanogaster_adult_ISCs_2.fastq.gz'
system(x)

x <- 
  'curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/025/SRR11906525/SRR11906525_subreads.fastq.gz -o SRR11906525_WGS_of_drosophila_melanogaster_female_adult_subreads.fastq.gz'
system(x)

x <- 
  'curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/042/SRR15130842/SRR15130842_1.fastq.gz -o SRR15130842_GSM5452672_Control_CM2_Drosophila_melanogaster_RNA-Seq_1.fastq.gz'
system(x)

x <- 
  'curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/042/SRR15130842/SRR15130842_2.fastq.gz -o SRR15130842_GSM5452672_Control_CM2_Drosophila_melanogaster_RNA-Seq_2.fastq.gz'
system(x)

x <- 
  'curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/041/SRR15130841/SRR15130841_1.fastq.gz -o SRR15130841_GSM5452671_Control_CM1_Drosophila_melanogaster_RNA-Seq_1.fastq.gz'
system(x)

x <- 
  'curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/041/SRR15130841/SRR15130841_2.fastq.gz -o SRR15130841_GSM5452671_Control_CM1_Drosophila_melanogaster_RNA-Seq_2.fastq.gz'
system(x)

# Download Drosophila melanogaster genome version r6.44 (released Jan 2022)

x <- 
  'curl -L http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.44_FB2022_01/fasta/dmel-all-chromosome-r6.44.fasta.gz -o dmel-all-chromosome-r6.44.fasta.gz'
system(x)
rm(x)

#change wdir back to original location
setwd("C:/Users/drewj/OneDrive/R_based_work/assembling_genome")

#-----VISUALIZE READ LENGTH DISTRIBUTION-----

#setwd to different directory (so we don't store large data files in a github repo)
setwd("C:/Users/drewj/OneDrive/R_based_work/assembling_genome_nonrepo")

#Run your CLI command through R (using system())

#NOTE: R is installed on Windows, so it doesn't necessarily know to use WSL to run commands, which is where gunzip is installed; thus, any time you need to run a command through WSL, start your system() command with 'wsl'

#NOTE: to capture output, include 'intern = TRUE' in your system command

x <- "wsl gunzip -c SRR11906525_WGS_of_drosophila_melanogaster_female_adult_subreads.fastq.gz | awk 'NR%4==2{print length($0)}'"

lengths_1 <- system(x,intern=TRUE)

x <- "wsl gunzip -c SRR12473480_Drosophila_PacBio_HiFi_UltraLow_subreads.fastq.gz | awk 'NR%4==2{print length($0)}'"

lengths_2 <- system(x,intern=TRUE)

x <- "wsl gunzip -c SRR13070625_Nanopore_sequencing_of_Drosophila_melanogaster_whole_adult_flies_pooled_male_and_female_1.fastq.gz | awk 'NR%4==2{print length($0)}'"

lengths_3 <- system(x,intern=TRUE)

#change wdir back to original location
setwd("C:/Users/drewj/OneDrive/R_based_work/assembling_genome")

library(ggplot2)
library(cowplot)
library(plyr)
library(dplyr)

#specify Sequencing platform
test1 <- as.data.frame(lengths_1)
test1$platform <- "PacBio_CLR"
test1 <- dplyr::rename(test1,"length"="lengths_1")
test2 <- as.data.frame(lengths_2)
test2$platform <- "PacBio_HiFi"
test2 <- dplyr::rename(test2,"length"="lengths_2")
test3 <- as.data.frame(lengths_3)
test3$platform <- "ONT"
test3 <- dplyr::rename(test3,"length"="lengths_3")

#combine all dfs
read_length_df <- rbind(test1,test2,test3)
read_length_df$length <- as.numeric(read_length_df$length)
rm(test1,test2,test3)

#calculate average read-lengths for each platform
summary_df <- ddply(read_length_df,"platform",summarise,grp.mean=mean(length,na.rm=TRUE))

#graph: read-length distribution
total.length.plot <- ggplot(read_length_df, aes(x=length, fill=platform, color=platform)) +
  geom_histogram(binwidth=100, alpha=0.5, position="dodge") +
  geom_vline(data=summary_df, aes(xintercept=grp.mean, color=platform), linetype="dashed", linewidth=0.2) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Read length (bp)", y = "Count") +
  theme_bw()
#SOMETHING NOT RIGHT, GRAPH NOT DISPLAYING

#graph: read-length distribution plot for reads <=20kb in length
t20kb.length.plot <- ggplot(read_length_df, aes(x=length, fill=platform, color=platform)) +
  geom_histogram(binwidth=50, alpha=0.5, position="dodge") +
  geom_vline(data=summary_df, aes(xintercept=grp.mean, color=platform), linetype="dashed", size=0.2) +
  scale_x_continuous(labels = scales::comma, limit = c(0,20000)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Read length (bp)", y = "Count") +
  theme_bw()

#merge both plots
plot <- plot_grid(total.length.plot,t20kb.length.plot,ncol=1)
#save the figure using the file name "read.length.pdf
pdf("read.length.pdf",width=6,height=8,paper='special')
print(plot)
dev.off()

#FINISH#