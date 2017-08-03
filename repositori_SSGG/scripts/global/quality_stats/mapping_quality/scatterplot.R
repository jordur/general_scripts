#!/share/apps/local/R/bin/Rscript

# Load libraries
library(car)
library(optparse)

option_list<-list(
	make_option(c("-v","--verbose"),action="store_true",help="Print extra output"),
	make_option(c("-q","--quietly"),action="store_false",default=FALSE,dest="verbose",help="Print little output [default]"),
	make_option(c("-f","--coverage_file1"),type="character",action="store",help="Region coverages file 1",metavar="file1"),
	make_option(c("-s","--coverage_file2"),type="character",help="Region coverages file 2",metavar="file2"),
	make_option(c("--name1"),type="character",default="coverage1",help="Sample name for coverage file 1",metavar="name1"),
	make_option(c("--name2"),type="character",default="coverage2",help="Sample name for coverage file 2",metavar="name2"),
	make_option(c("-o","--output"),type="character",default="scatterplot.png",help="Output name (default scatterplot.png)",metavar="path name")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

cat("Script for creating a scatterplot image of coverages_file2 vs coverages_file1\n")
args<-commandArgs()
if (is.null(opt$coverage_file1) | is.null(opt$coverage_file2)){
	cat("ERROR: Not enough arguments!!!!\n")
	cat("Call script with -h or --help parameter for getting usage instructions!!\n")
	quit()
}

# Arguments
coverage_file1<-opt$coverage_file1
coverage_file2<-opt$coverage_file2

# Read data
cov1_data<-read.table(coverage_file1,sep="\t")
cov2_data<-read.table(coverage_file2,sep="\t")

# Obtain coverage values in vector format, so that histogram functions can be applied. Please note that the 4th column contains the coverage data
cov1<-as.vector(as.matrix(cov1_data[6]))
cov2<-as.vector(as.matrix(cov2_data[6]))

# Plot histogram
png(opt$output,width=768,height=768)
sc<-scatterplot(cov2~cov1,col=c("red","green","blue"),main="Scatterplot of coverages",xlab=opt$name1,ylab=opt$name2)
abline(fit<-lm(cov2~cov1), col="black")
legend("topright", bty="n", legend=paste("R2=",format(summary(fit)$adj.r.squared, digits=4)))
dev.off()