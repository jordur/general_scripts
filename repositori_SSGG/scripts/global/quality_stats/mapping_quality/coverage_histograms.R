#!/share/apps/local/R/bin/Rscript

args<-commandArgs()
if (is.na(as.character(args[7]))){
	cat("ERROR: Not enough arguments!!!!\n")
	cat("Script for obtaining the coverage_histogram of a sample\n")
	cat("usage:",as.character(args[4]),"sample.pileup sample_name\n")
	cat("\twhere sample.pileup is the pileup file coming from mapping stats\n")
	quit()
}

# Arguments
sample_file<-as.character(args[6])
sample_name<-as.character(args[7])

# Read data
cov_data<-read.table(sample_file,sep="\t")

# Obtain coverage values in vector format, so that histogram functions can be applied. Please note that the 4th column contains the coverage data
coverage<-as.vector(as.matrix(cov_data[4]))

# Plot histogram
png(paste(sample_name,"_histogram.png",sep=""),width=768,height=768)
h<-hist(coverage,breaks=200,col="red",xlab="coverage",main=paste("Coverage histogram of",sample_name))
dev.off()

# For zooming in and setting limits of the new plot
png(paste(sample_name,"_histogram_low_coverages.png",sep=""),width=768,height=768)
h<-hist(coverage,breaks=200,col="red",xlab="coverage",main=paste("Low coverages histogram of",sample_name),xlim=c(0,100),ylim=range(h$counts[1:100]))
dev.off()