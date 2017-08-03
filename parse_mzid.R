## Script per parsejar resultats mzid de forma que apareguin en forma de Taul Excel ##
## Els noms dels espectres hande ser del tipus Css404_12c_snakes.mzid per a quie funcionen (separacio per "_" ###

source("https://bioconductor.org/biocLite.R")
biocLite("mzID")
#browseVignettes("mzID")
setwd("/home/bec2-jcalvete/Feina_Jordi/proteogenomics/Scutulatus")
library(mzID)
#mzResults <- mzID("MSGF_output.mzid")
myFiles <- dir(pattern="*mzid")

for (i in 1:length(myFiles)){
    tempo=strsplit(myFiles[i],"_")
    nomspec=paste(tempo[[1]][1],tempo[[1]][2], sep="_")
    mzResults <- mzID(myFiles[i])
    #mzResults <- mzID("Css404_13b_snakes.mzid")
    #mzResults
    flatResults <- flatten(mzResults)
    #names(flatResults)
###per obtindre els resultats d'una fila en concret
    #flatResults[1,]
    names<-c("Description,Rank,Charge,m/z,Seq,SpecEvalue,Accession")
    description<-paste(flatResults$description,flatResults$rank,flatResults$chargestate,flatResults$calculatedmasstocharge,flatResults$pepseq,flatResults$`ms-gf:specevalue`,flatResults$accession, sep=",")
    final<-c(names, description)
    write(final, file=paste(nomspec,".csv", sep=""))
}
