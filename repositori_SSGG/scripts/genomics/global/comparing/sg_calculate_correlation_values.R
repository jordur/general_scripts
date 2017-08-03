#!/share/apps/local/R/bin/Rscript

######################################################
# Copyright 2013 - Sistemas Genomicos S.L.           #
# @Desc: Script to generate depth of coverage        #
# calculates correlation values from mapping         #
# results between samples                            #
# @Author: Sheila                                    #
# @Date: June 2013                                   #
# @Version: v0.1                                     #
######################################################




args<-commandArgs();
#Argument 6 is the first you can really pass to R, the other arguments are for internal use. To check the meaning of other arguments you can do print(args)
y<-args[6]; ## matrix file from the samtools depth command with depth per base
outpath <-args[7]
outfile<-paste(outpath,"/correlation_values_between_samples.tsv",sep="");
tab <- read.table(y,header=TRUE)
write.table(cor(tab),outfile)
