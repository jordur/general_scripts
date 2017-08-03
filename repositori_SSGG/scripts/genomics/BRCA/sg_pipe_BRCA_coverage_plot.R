#!/share/apps/local/R/bin/Rscript

######################################################
# Copyright 2013 - Sistemas Genomicos S.L.           #
# @Desc: Script to generate depth of coverage        #
# plots for resequencing studies based on amplicons  #
# @Author: Sheila                                    #
# @Date: April 2013                                  #
# @Version: v0.1                                     #
######################################################




args<-commandArgs();
#Argument 6 is the first you can really pass to R, the other arguments are for internal use. To check the meaning of other arguments you can do print(args)
y<-args[6]; ## Amplicon counts ordered by GENE
x<-args[7]; ## Amplicon counts ordered by MULTIPLEX NUMBER
a<-args[8]; ## Amplicon counts ordered by GENE and normalized
b<-args[9]; ## Amplicon counts ordered by MULTIPLEX NUMBER and normalized
maxdepth <-as.numeric(args[10])
maxdepthnorm <- as.numeric(args[11])
#legendDepth <- maxdepth/2
legendDepthNorm <- maxdepthnorm/2
outpath <-args[12]
#global parameters
colours <- c("THISISNOTACOLOR","cadetblue4", "burlywood4", "chocolate", "darkblue", "darkgoldenrod3", "darkmagenta", "darkolivegreen", "deeppink4", "gold4", "firebrick", "dodgerblue4", "gray9", "mediumpurple4", "red4", "violetred4", "darkslategrey", "gray0", "orange", "aquamarine3","green3")
#colours <- c("THISISNOTACOLOR","aquamarine3","red","orange","green3","deeppink","blue","lightpink","black","blueviolet","darkorange","darkmagenta","goldenrod1","deeppink4","grey","brown","chartreuse2","wheat2","violet")
amplicountsGene<-read.table(y,sep="\t",header=TRUE);
amplicountsMultiplx<-read.table(x,sep="\t",header=TRUE);
amplicountsGeneNorm<-read.table(a,sep=" ",header=TRUE);
amplicountsMultiplxNorm<-read.table(b,sep=" ",header=TRUE);
numcols <- ncol(amplicountsGene)
samples <- colnames(amplicountsGene)
outfilegene<-paste(outpath,"/Amplicon_counts_by_gene_position.pdf",sep="");
pdf(outfilegene,width=14)
#o<-par(xpd=T, mar=c(6,4,4,8)+0.1)
o<-par(xpd=T,mar=c(5, 4, 4, 7)+0.1)
plot(amplicountsGene[[samples[2]]],type="o", col=colours[2],xaxt="n",ylab="Depth of coverage",main="Depth of coverage ordered by gene",cex=0.5,cex.axis=1,cex.lab=1,xlab="",ylim=c(0,maxdepth),lwd = 3)
# Places a legend at the appropriate place c("Health","Defense")
#legendxy<-c(0.000003,6000)
legend(98,maxdepth, legend = samples[2:length(samples)], col = colours[2:numcols],lwd = 3, bty = "n",cex=0.7)
par(o)
#It starts at 3 since the first column is the amplicon name and the second has already been plotted above.
for (i in 3:length(samples))
{
	lines(amplicountsGene[[samples[i]]], col=colours[i],type="o",cex.lab=1,cex=0.5, lwd = 2)
}
axis(1, at=1:93,labels=amplicountsGene$AMPLICON,las = 2,cex.axis=0.4)
dev.off()
outfilemultiplex<-paste(outpath,"/Amplicon_counts_by_multiplex.pdf",sep="");
pdf(outfilemultiplex,width=14)
o<-par(xpd=T,mar=c(5, 4, 4, 7)+0.1)
plot(amplicountsMultiplx[[samples[2]]],type="o", col=colours[2],xaxt="n",ylab="Depth of coverage",main="Depth of coverage ordered by multiplex",cex=0.5,cex.axis=1,cex.lab=1,xlab="",ylim=c(0,maxdepth),lwd = 3)
# Places a legend at the appropriate place c("Health","Defense")
legend(98,maxdepth, legend = samples[2:length(samples)], col = colours[2:numcols],lwd = 3, bty = "n",cex=0.7)
par(o)
#It starts at 3 since the first column is the amplicon name and the second has already been plotted above.
for (i in 3:length(samples))
{
    lines(amplicountsMultiplx[[samples[i]]], col=colours[i],type="o",cex.lab=1,cex=0.5,  lwd = 2)
}
axis(1, at=1:93,labels=amplicountsMultiplx$AMPLICON,las = 2,cex.axis=0.4)
dev.off()

#####
outfilegeneNorm<-paste(outpath,"/Amplicon_counts_by_gene_position_normalized.pdf",sep="");
pdf(outfilegeneNorm,width=14)
o<-par(xpd=T,mar=c(5, 4, 4, 7)+0.1)
plot(amplicountsGeneNorm[[samples[2]]],type="o", col=colours[2],xaxt="n",ylab="Depth of coverage",main="Normalized depth of coverage ordered by gene",cex=0.5,cex.axis=1,cex.lab=1,xlab="",ylim=c(0,maxdepthnorm),lwd = 3)
# Places a legend at the appropriate place
legend(98,maxdepthnorm, legend = samples[2:length(samples)], col = colours[2:numcols],lwd = 3, bty = "n",cex=0.7)
par(o)
#It starts at 3 since the first column is the amplicon name and the second has already been plotted above.
for (i in 3:length(samples))
{
    lines(amplicountsGeneNorm[[samples[i]]], col=colours[i],type="o",cex.lab=1,cex=0.5, lwd = 2)
}
axis(1, at=1:93,labels=amplicountsGeneNorm$AMPLICON,las = 2,cex.axis=0.4)
dev.off()
outfilemultiplexNorm<-paste(outpath,"/Amplicon_counts_by_multiplex_normalized.pdf",sep="");
pdf(outfilemultiplexNorm,width=14)
o<-par(xpd=T,mar=c(5, 4, 4, 7)+0.1)
plot(amplicountsMultiplxNorm[[samples[2]]],type="o", col=colours[2],xaxt="n",ylab="Depth of coverage",main="Normalized depth of coverage ordered by multiplex",cex=0.5,cex.axis=1,cex.lab=1,xlab="",ylim=c(0,maxdepthnorm),lwd = 3)
# Places a legend at the appropriate place c("Health","Defense")
legend(98,maxdepthnorm, legend = samples[2:length(samples)], col = colours[2:numcols],lwd = 3, bty = "n",cex=0.7)
par(o)
#It starts at 3 since the first column is the amplicon name and the second has already been plotted above.
for (i in 3:length(samples))
{
    lines(amplicountsMultiplxNorm[[samples[i]]], col=colours[i],type="o",cex.lab=1,cex=0.5, lwd = 2)
}
axis(1, at=1:93,labels=amplicountsMultiplxNorm$AMPLICON,las = 2,cex.axis=0.4)
dev.off()
