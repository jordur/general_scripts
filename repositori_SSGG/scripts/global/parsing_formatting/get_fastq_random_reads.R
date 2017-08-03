#!/share/apps/local/R/bin/Rscript
# author: Arbol <oscar.rodriguez@sistemasgenomicos.com>

# Load libraries
library(car)
library(optparse)
library(ShortRead)

option_list<-list(
	make_option(c("-v","--verbose"),action="store_true",help="Print extra output"),
	make_option(c("-q","--quietly"),action="store_false",default=FALSE,dest="verbose",help="Print little output [default]"),
	make_option(c("-f","--R1_fastq"),type="character",action="store",help="R1 fastq file (script also works with fastq.gz)",metavar="R1"),
	make_option(c("-r","--R2_fastq"),type="character",action="store",help="R2 fastq file (script also works with fastq.gz)",metavar="R2"),
	make_option(c("-n","--reads_number"),type="integer",default="1e6",help="Number of random sampling reads to get from pair-end reads files",metavar="n"),
	make_option(c("-o","--output"),type="character",default="./",action="store",help="Output name (default scatterplot.png)",metavar="path name to output"),
	make_option(c("--output_R1"),type="character",default="",action="store",help="R1 output fastq file (script also works with fastq.gz)",metavar="R1 output"),
	make_option(c("--output_R2"),type="character",default="",action="store",help="R2 output fastq file (script also works with fastq.gz)",metavar="R2 output")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

cat("Script for random sampling reads from fastq (or fastq.gz) pair-end reads files\n")
args<-commandArgs()
if (is.null(opt$R1_fastq) | is.null(opt$R2_fastq)){
	cat("ERROR: Not enough arguments!!!!\n")
	cat("Call script with -h or --help parameter for getting usage instructions!!\n")
	quit()
}

# read fastq or fastq.gz files
f1 <- FastqSampler(opt$R1_fastq, n=opt$reads_number)
f2 <- FastqSampler(opt$R2_fastq, n=opt$reads_number)
# get the random reads
set.seed(123L); p1 <- yield(f1)
set.seed(123L); p2 <- yield(f2)
# output results files
if ((opt$output_R1 != "") && (opt$output_R2 != "")){
  writeFastq(p1,file=file.path(opt$output_R1))
  writeFastq(p2,file=file.path(opt$output_R2))
} else {
  writeFastq(p1,file=file.path(opt$output,paste("R1_",opt$reads_number,".fastq",sep="")))
  writeFastq(p2,file=file.path(opt$output,paste("R2_",opt$reads_number,".fastq",sep="")))
}