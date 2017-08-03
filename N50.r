## Script per calcular N50 des del fitxer stats.txt de cada ensamblatge fet amb Oases ##
## NOTA: cal ajustar el paràmetre myKmerValue acada k-mer utilitzat en cada cas ###
myStatsTable<-read.table("stats.txt",header=TRUE)
contigs<-rev(sort(myStatsTable$lgth+myKmerValue-1))
n50<-contigs[cumsum(contigs) >= sum(contigs)/2][1]
