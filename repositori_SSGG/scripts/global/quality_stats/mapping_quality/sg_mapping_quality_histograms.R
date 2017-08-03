#!/share/apps/local/R/bin/Rscript

args<-commandArgs();
x<-as.character(args[6]);
xx<-read.table(x,sep="\t");
y<-as.character(args[7]);
yy<-read.table(y,sep="\t");

mean_xx<-colMeans(xx);
mean_yy<-colMeans(yy);

max_xx<-max(xx);
max_yy<-max(yy);
max<-0;

density_xx<-density(xx$V1);
density_yy<-density(yy$V1);

max_xx<-as.numeric(max(density_xx$y));
max_yy<-as.numeric(max(density_yy$y));

if(max_xx > max_yy)
{
        max<-as.numeric(max_xx);
}else
{
        max<-as.numeric(max_yy);
}

pdf("mapping_quality_density.pdf");
plot(density_xx,main="Mapping Quality",xlab="QV Phred",col="red");
lines(density_yy,col="blue");
legend(0,max,c(x,y),lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red"),cex=0.6)

dev.off();

