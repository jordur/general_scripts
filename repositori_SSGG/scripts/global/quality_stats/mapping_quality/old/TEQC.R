#!/share/apps/local/R/bin/Rscript

suppressWarnings(suppressPackageStartupMessages(library(TEQC)));

args<-commandArgs();
x<-as.character(args[6]);##BAM
y<-as.character(args[7]);#targets
z<-as.character(args[8]);##Chr
v<-as.character(args[9]);##THRESHOLDS VALUES
offset_value<-as.character(args[10]);##OFFSET VALUE
offset_value<-as.numeric(offset_value);

vv<-do.call(rbind, strsplit(v, ","));
result<-(dim(length(vv)));


reads<-get.reads(x,filetype="bam");
targets<-get.targets(y,sep="\t");
coverage<-suppressWarnings(coverage.target(reads,targets,perTarget=T,perBase=T));
data_frame<-as.data.frame(coverage$targetCoverages);
suppressWarnings(write.table(data_frame,file="coverage_targetCoverage.txt",sep="\t",append=TRUE,col.name=FALSE));
#c<-suppressWarnings(covered.k(coverage$coverageTarget,k=c(1,10,20)));
result[1]<-z;

for(i in 1:length(vv))
{
	j<-i+1;
	c<-suppressWarnings(covered.k(coverage$coverageTarget,k=i));
	result[j]<-c;
}

result<-as.vector(result);
tresult<-t(result);
###
#tresult2<-paste(z,tresult,sep="\t");
suppressWarnings(write.table(tresult,file="coverage_threshold.txt",sep="\t",append=TRUE));
###


frt<-fraction.reads.target(reads,targets);
frt<-paste(z,fraction.reads.target(reads,targets),sep="\t");
if(offset_value > 0)
{
	frt_offset<-paste(z,fraction.reads.target(reads,targets,Offset=offset_value),sep="\t");
#	frt_offset<-fraction.reads.target(reads,targets,Offset=offset_value);
	suppressWarnings(write.table(frt_offset,file="reads_target_offSet",sep="\t",append=TRUE));
}
suppressWarnings(write.table(frt,file="reads_target",sep="\t",append=TRUE));

