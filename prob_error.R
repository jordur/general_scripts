setwd("/home/bec2-jcalvete/Feina_Jordi/MALPOLON/rawdata")
qual1 <- read.delim("150807_SND405_A_L003_GZX-17_R1.fastq.quality")
qual2 <-read.delim("150807_SND405_A_L003_GZX-17_R2.fastq.quality")
par(mfrow=c(2,1))
boxplot(t(qual1), col='light blue', ylim=c(0,.4), frame.plot=F, outline=F, xaxt = "n", ylab='Probability of nucleotide error', xlab='Nucleotide Position', main='Read1')
axis(1, at=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130), labels=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130))
boxplot(t(qual2), col='light blue', ylim=c(0,.4), frame.plot=F, outline=F, xaxt = "n", ylab='Probability of nucleotide error', xlab='Nucleotide Position', main='Read2')
axis(1, at=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130), labels=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130))
