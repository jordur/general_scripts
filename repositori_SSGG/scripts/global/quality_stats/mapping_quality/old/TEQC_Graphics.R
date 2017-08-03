#!/share/apps/local/R/bin/Rscript

args<-commandArgs();
y<-as.character(args[6]);#COVERAGE
z<-as.character(args[7]);#THRESHOLD
d<-as.character(args[8]);#CROMOSOME
s<-as.character(args[9]);#SPECIFITY

offset_value<-as.character(args[10]);##OFFSET
offset_value<-as.numeric(offset_value);


th<-as.character(args[11]);##THRESHOlDS
threshold<-do.call(rbind, strsplit(th, ","));
cols_legend<-"";
###VER QUE EL NÚMERO DE THRESHOLDS NO SUPERE EL DE 10

cols<-(dim(10));
cols[2]<-'#4A7BAD';
cols[3]<-'#A54239';
cols[4]<-'#8CAD4A';
cols[5]<-'#FFD700';
cols[6]<-'#4B0082';
cols[7]<-'#000000';
cols[8]<-'#A9A9A9';
cols[9]<-'#FF69B4';
cols[10]<-'#008080';
##HACER MAS COLORES Y LINEAS

x<-read.table(y,sep="\t");
x<-as.matrix(x$V6);
H<-hist(x, breaks = 50, plot = FALSE);

x_limit<-max(H$intensities);
entry<-0;

tab<-table(x)
cs<-cumsum(rev(tab))
m<-max(H$counts)

cs.s<-cs*(m/sum(tab))

entry<-0;
for(i in 1:length(H$counts))
{
	if(H$counts[i] <=0 && entry==0)
	{
		entry<-1;
		mx<-H$breaks[i];
	}
}

pdf("coverage_distribution.pdf");
plot(H, freq = TRUE, xlab = "Coverage", main = "Coverage Distribution",col="#336699",xlim=c(0,mx))

lines(x=names(cs.s), y=cs.s, col="red", lwd=2)
axis(side = 4, at = seq(0, m, length.out = 11), labels = seq(0,1, by = 0.1))
dev.off()

k<-read.table(z,sep="\t");
min_k<-2;

for(i in 3:length(names(k)))
{
	if(min_k > min(k[i]))
	{
		min_k<-min(k[i]);
	}
}

chr<-read.table(d,sep="\n");
chr_coord<-row.names(chr);

pdf("coverage_threshold.pdf",width=14);
plot(k$V3,cex=0.001,col=cols[2],ylim=c(min_k,1.1),ylab="% Bases",xlab="Chromosome",cex.axis=0.00000001,lwd=2,type="n")

for(i in 3:length(names(k)))
{
	if(i <11) ##LIMIT IN 10 VALUES
	{
		j<-i-1;
		jj<-i-2;
		lines(k[i],col=cols[j],lwd=2,type="b")
		cols_legend[jj]<-cols[j];
	}
}
axis(1,at=chr_coord,labels=chr$V1,las=2)
axis(2)
legend(0.5,1.1,threshold,lty=c(1,1,1),lwd=c(2,2,2),col=cols_legend);
dev.off()

r<-read.table(s,sep="\t");
rr<-r[,-1];
rr<-as.matrix(t(rr*100));

pdf("specificity_plot.pdf");
if(offset_value > 0 )
{
	barplot(rr,xlab="Chromosome",ylab="% Bases",col=c('#4A7BAD','#A54239'))
}
if(offset_value ==0)
{
	barplot(rr,xlab="Chromosome",ylab="% Bases",col='#4A7BAD')
}
axis(1,at=chr_coord,labels=chr$V1,las=2)
axis(2)
#legend AÑADIR MAS TARDE
dev.off()

