#assembling a draft genome and analyzing strucutral variants with long-read sequencing technologies
#https://www.sciencedirect.com/science/article/pii/S2666166722003860

#-----SETUP-----

#setwd
setwd("C:/Users/drewj/OneDrive/R_based_work/assembling_genome")


#-----get data-----

#setwd to different directory (so we don't store large data files in a github repo)
setwd("C:/Users/drewj/OneDrive/R_based_work/assembling_genome_nonrepo")

#download data to local
system(
  "curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR130/025/SRR13070625/SRR13070625_1.fastq.gz -o SRR13070625_Nanopore_sequencing_of_Drosophila_melanogaster_whole_adult_flies_pooled_male_and_female_1.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/080/SRR12473480/SRR12473480_subreads.fastq.gz -o SRR12473480_Drosophila_PacBio_HiFi_UltraLow_subreads.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/022/SRR12099722/SRR12099722_1.fastq.gz -o SRR12099722_WGS_Drosophila_melanogaster_adult_ISCs_1.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/022/SRR12099722/SRR12099722_2.fastq.gz -o SRR12099722_WGS_Drosophila_melanogaster_adult_ISCs_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/025/SRR11906525/SRR11906525_subreads.fastq.gz -o SRR11906525_WGS_of_drosophila_melanogaster_female_adult_subreads.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/042/SRR15130842/SRR15130842_1.fastq.gz -o SRR15130842_GSM5452672_Control_CM2_Drosophila_melanogaster_RNA-Seq_1.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/042/SRR15130842/SRR15130842_2.fastq.gz -o SRR15130842_GSM5452672_Control_CM2_Drosophila_melanogaster_RNA-Seq_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/041/SRR15130841/SRR15130841_1.fastq.gz -o SRR15130841_GSM5452671_Control_CM1_Drosophila_melanogaster_RNA-Seq_1.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/041/SRR15130841/SRR15130841_2.fastq.gz -o SRR15130841_GSM5452671_Control_CM1_Drosophila_melanogaster_RNA-Seq_2.fastq.gz

# Download Drosophila melanogaster genome version r6.44 (released Jan 2022)

> wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.44_FB2022_01/fasta/dmel-all-chromosome-r6.44.fasta.gz"
)

#change wdir back to original location
setwd("C:/Users/drewj/OneDrive/R_based_work/assembling_genome")

#-----VISUALIZE READ LENGTH DISTRIBUTION-----

#setwd to different directory (so we don't store large data files in a github repo)
setwd("C:/Users/drewj/OneDrive/R_based_work/assembling_genome_nonrepo")

#run this through CLI
system(

#Create a new file and generate a header line
"echo "platform,length" -> length.csv"

#Add each read length into the length.csv file.
"bioawk -c fastx '{print "PacBio_CLR," length($seq)}' SRR11906525_WGS_of_drosophila_melanogaster_female_adult_subreads.fastq.gz >> length.csv"

bioawk -c fastx '{print "PacBio_HiFi," length($seq)}' SRR12473480_Drosophila_PacBio_HiFi_UltraLow_subreads.fastq.gz >> length.csv

bioawk -c fastx '{print "ONT," length($seq)}' SRR13070625_1.fastq.gz >> length.csv

)


#change wdir back to original location
setwd("C:/Users/drewj/OneDrive/R_based_work/assembling_genome")