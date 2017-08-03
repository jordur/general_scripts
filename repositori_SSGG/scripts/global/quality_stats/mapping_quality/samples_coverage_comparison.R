#!/share/apps/local/R/bin/Rscript

# Load libraries
library(car)
library(optparse)
library(aCGH)

option_list<-list(
	make_option(c("-v","--verbose"),action="store_true",help="Print extra output"),
	make_option(c("-q","--quietly"),action="store_false",default=FALSE,dest="verbose",help="Print little output [default]"),
	make_option(c("-s","--samples"),type="character",action="store",help="Region coverages file 1",metavar="comma separated files list"),
	make_option(c("-n","--names"),type="character",help="Sample names",metavar="comma separated names list"),
	make_option(c("-t","--thresholds"),type="character",default="1,10,20",help="Thresholds list",metavar="comma separated values"),
	make_option(c("-o","--output"),type="character",default="below_thresholds.png",help="Output name (default scatterplot.png)",metavar="path name")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

cat("Script for creating a comparison barplot image of percentage of regions below thresholds\n")
args<-commandArgs()
if (is.null(opt$samples)){
	cat("ERROR: Not enough arguments!!!!\n")
	cat("Call script with -h or --help parameter for getting usage instructions!!\n")
	quit()
}

# Arguments
samples<-unlist(strsplit(opt$samples,","))
thresholds<-as.numeric(unlist(strsplit(opt$thresholds,",")))
if (is.null(opt$names)){
  names<-c()
  for (sample in samples){
    names[length(names)+1]<-paste("s", length(names)+1,sep="")
  }
} else {
  names<-unlist(strsplit(opt$names,","))
}

# Read data
percentages<-c()
for (sample in samples){
  cov_data<-read.table(sample,sep="\t")
  
  # Create comparative vectors of acumulative percentage of mean coverage regions
  perc_under_threshold<-c()
  regions=length(as.matrix(cov_data[6]))
  for (thres in thresholds){
    perc_under_threshold[length(perc_under_threshold)+1] = length(which(as.matrix(cov_data[6]<thres)))/regions*100
  }
  
  # Join all vectors
  percentages<-rbind(percentages,perc_under_threshold)
}

# Barplot
png(opt$output,width=768,height=768)
bp<-barplot(percentages,col=c("darkblue","darkred"),beside=TRUE,main="Percentage of regions whose mean coverage is below threshold",
            xlab="Coverage thresholds",names.arg=thresholds,legend.text=names,args.legend=list(x="topleft"))
dev.off()