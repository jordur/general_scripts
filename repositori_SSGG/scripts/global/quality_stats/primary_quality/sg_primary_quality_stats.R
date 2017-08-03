#!/share/apps/local/R/bin/Rscript
library(lattice)

args<-commandArgs();
y<-as.character(args[6]);
x<-read.table(y,sep="\t");##QV PER POSITIONS
x_ylim=max(x$V3)+1;

yy<-as.character(args[7]);
z<-read.table(yy,sep="\t");

yyy<-as.character(args[8]);##MEANS SEQUENCE QUALITY
s<-read.table(yyy,sep="\t");
s1<-as.numeric(s$V1);
s2<-as.numeric(s$V2);
max_s2<-max(s2);

y4<-as.character(args[9]);
r<-read.table(y4,sep="\t");
r1<-as.numeric(r$V1);
r2<-as.numeric(r$V2);
max_r2<-max(r2);

y5<-as.character(args[10]);##MEAN INDETERMINATION
t<-read.table(y5,sep="\t");
t1<-as.numeric(t$V1)
t2<-as.numeric(t$V2);
max_t2<-max(t2);

name<-as.character(args[11]);


stringA<-paste(name,"_quality_value_per_position.pdf",sep="");

pdf(stringA);

plot(x$V3,cex=0.001,ylim=c(-1,x_ylim),xlab="positions in reads",ylab="QVs",col='#A54239',type="l",lwd=2,main="Quality Value Per Position");
abline(h=5,lwd=0.25,col='#949494');
abline(h=10,lwd=0.25,col='#949494');
abline(h=15,lwd=0.25,col='#949494');
abline(h=20,lwd=0.25,col='#949494');
abline(h=25,lwd=0.25,col='#949494');
abline(h=30,lwd=0.25,col='#949494');
lines(x$V2,col='#4A7BAD',lwd=2);
lines(x$V4,col='#8CAD4A',lwd=2);
lines(x$V5,col='#BDADCE',lwd=2);

dev.off()

stringB<-paste(name,"_mean_sequence_quality.pdf",seep="");
pdf(stringB);
plot(s2,cex=0.001,xlab="mean sequence quality (Phred score)",ylab="% sequences",col='#4A7BAD',lwd=2,type='l',main="Mean Quality Per Sequence");
dev.off()

stringC<-paste(name,"_mean_indetermination.pdf",sep="");
pdf(stringC);
plot(t2,cex=0.001,xlab="mean indetermination",ylab="% sequences",col='#4A7BAD',lwd=2,type='l',main="Mean Indetermination per Sequence");
dev.off();

stringD<-paste(name,"_quality_score_distribution_per_position.pdf",sep="");
pdf(stringD);
barplot(z$V2,col='#4A7BAD',width=0.7,xlab="QV",ylab="%",main="Quality Value Distribution Per Position");
axis(1,xlim=c(-1,x_ylim));
axis(2,ylim=c(0,max_s2));
#box();
par(new=TRUE);
F<-ecdf(z$V2);
plot(F,ylim=c(0,1),axes=FALSE,xlab="",ylab="",main="",do.points = TRUE,col='#A54239');
axis(4,ylim=c(0,1));
dev.off();

stringE<-paste(name,"_indetermination_distribution_per_position.pdf",sep="");
pdf(stringE);
barplot(r2,col='#A54239',width=0.7,xlab="Positions",ylab="%",main="Indetermination Distribution Per Position");
axis(1,xlim=c(-1,x_ylim));
axis(2,ylim=c(0,max_r2));
par(new=TRUE);
FF<-ecdf(r2);
plot(FF,ylim=c(0,1),axes=FALSE,xlab="",ylab="",main="",do.points = TRUE,col='#4A7BAD');
axis(4,ylim=c(0,1));

dev.off();
