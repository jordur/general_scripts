#!/share/apps/local/R/bin/Rscript

# Load libraries
library(car)
library(optparse)

option_list<-list(
  make_option(c("-v","--verbose"),action="store_true",help="Print extra output"),
  make_option(c("-q","--quietly"),action="store_false",default=FALSE,dest="verbose",help="Print little output [default]"),
  make_option(c("-f","--freqs_file"),type="character",action="store",help="Sample alternative allele frequencies file",metavar="file"),
  make_option(c("-o","--output"),type="character",default="freqs_histogram.png",help="Output name (default freqs_histogram.png)",metavar="path name")
)


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

cat("Script for creating histogram plots of alternative allelic frequencies\n")
args<-commandArgs()
if (is.null(opt$freqs_file)){
  cat("ERROR: Not enough arguments!!!!\n")
  cat("Call script with -h or --help parameter for getting usage instructions!!\n")
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