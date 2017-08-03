#!/share/apps/local/R/bin/Rscript

#suppressWarnings(suppressPackageStartupMessages(library(TEQC)));

library(TEQC)
library(biovizBase);

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
targets<-get.targets(y,sep="\t",header=FALSE,skip=0);

coverage<-suppressWarnings(coverage.target(reads,targets,perTarget=T,perBase=T));
data_frame<-as.data.frame(coverage$targetCoverages);

file_intervals<-paste(z,"coverage_targetCoverage.txt",sep="_");

gr<-GRanges(seqnames=targets$space,ranges=IRanges(targets$ranges),strand=rep("*"));
for(i in 1:length(gr))
{
	pileup<-pileupAsGRanges(x,gr[i]);
	pileup_OK<-as.data.frame(pileup);
	pileup_OK2<-paste(pileup_OK$seqnames,pileup_OK$start,pileup_OK$depth,sep="\t");
	
	string<-paste(z,i,"coverage_bases.txt",sep="_");
	suppressWarnings(write.table(pileup_OK2,file=string,sep="\t",append=TRUE,col.name=FALSE));
}

suppressWarnings(write.table(data_frame,file=file_intervals,sep="\t",append=TRUE,col.name=FALSE));
result[1]<-z;

for(i in 1:length(vv))
{
	j<-i+1;
	c<-suppressWarnings(covered.k(coverage$coverageTarget,k=i));
	result[j]<-c;
}

result<-as.vector(result);
tresult<-t(result);

file_threshold<-paste(z,"coverage_threshold.txt",sep="_");

suppressWarnings(write.table(tresult,file=file_threshold,sep="\t",append=TRUE));

frt<-fraction.reads.target(reads,targets);
frt<-paste(z,fraction.reads.target(reads,targets),sep="\t");
if(offset_value > 0)
{
	file_specificity_offset<-paste(z,"reads_target_offSet",sep="_");

	frt_offset<-paste(z,fraction.reads.target(reads,targets,Offset=offset_value),sep="\t");
	suppressWarnings(write.table(frt_offset,file=file_specificity_offset,sep="\t",append=TRUE));
}
file_specificity<-paste(z,"reads_target",sep="_");

suppressWarnings(write.table(frt,file=file_specificity,sep="\t",append=TRUE));