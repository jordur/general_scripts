#!/share/apps/local/R/bin/Rscript

args<-commandArgs();
x<-as.character(args[6]);#VectorA
y<-as.character(args[7]);#VectorB

xx<-read.csv(x,sep=",",header=FALSE,row.names=1);
yy<-read.csv(y,sep=",",header=FALSE,row.names=1);

for(i in 1:24)
{
	a<-as.vector(xx[i,]);
	b<-as.vector(yy[i,]);

	a[is.na(a)] <- 0;
	b[is.na(b)] <- 0;

	c<-cor(as.numeric(a),as.numeric(b));
	print(paste(i,c,sep=" "));



}

