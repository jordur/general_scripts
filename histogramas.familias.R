
                                        #Una proteïna en diferents serp
pdf("Figure1")
lista <- dir (pattern = ".txt")
layout(matrix(c(1:12), 4, 3))
ntotal <-c(1131,26711,2614,9370,1269,8948,2778,5293)
table.fami  <- matrix(0, ncol = 8, nrow = length(lista))
fami <- rep(0,8)
colnames (table.fami) <- names(fami)
for (i in 1 : length (lista)) {
  fam <- read.delim (lista[i], header = F)
  names(fami) <- c("mid1", "mid2", "mid3", "mid4", "mid5", "mid6", "mid7", "mid8")
  fami[names(table(fam[,4]))] <- table(fam[,4])
  barplot(fami/ntotal, ylab = "abundance",xlab = "")
  table.fami[i,] <- fami/ntotal 
  title (sub(".class.txt", "", lista[i]))
  print(fami/ntotal)
}
dev.off()

#Distàncies entre serps
colnames (table.fami) <- names(fami)
rownames(table.fami) <- sapply(lista, function (x) {sub(".class.txt", "",x)})
#table.fami <- table.fami[,-1]
pca <- princomp(table.fami, cor = T)
pca$scores[,1] <- rescale (pca$scores[,1], newmin = min(pca$loadings[,1]), newmax=max(pca$loadings[,1]))
pca$scores[,2] <- rescale (pca$scores[,2], newmin = min(pca$loadings[,2]), newmax=max(pca$loadings[,2]))

explained <- pca$sdev/sum( pca$sdev) * 100
pdf("Figure4")
plot(pca$loadings[,1], pca$loadings[,2], type = "n", main = "Por su veneno las conocerás",
    xlab = paste("PC1:", round(explained[1],2), "% explained variance", sep = ""),
    ylab = paste("PC2:", round(explained[2],2), "% explained variance", sep = "") )
text(pca$loadings[,1], pca$loadings[,2], label = names (fami))
text(pca$scores[,1], pca$scores[,2], label = rownames (table.fami), col = "blue")
dev.off()

#Una proteïna en diferents serps
layout(matrix(c(1:12), 4, 3))
for ( i in 1 : nrow(table.fami)) {
   barplot(table.fami[i,], ylab = "abundance",xlab = "")
   title (rownames(table.fami)[i])
       }

#Les proteïnes de cada serp en pie
pdf("Figure2")
layout(matrix(c(1:8), 4, 2))
colnames (table.fami) <- names(fami)
rownames(table.fami) <- sapply(lista, function (x) {sub(".class.txt", "",x)})
for ( i in 1 : nrow(t(table.fami))) {
   pie(t(table.fami)[i,],labels = rownames(table.fami) , radius = 1.1)
   title (rownames(t(table.fami))[i])
       }
dev.off()


#Les proteïnes de cada serp en barres
pdf("Figure3")
lista <- dir (pattern = ".txt")
layout(matrix(c(1:8), 4, 2))
table.fami  <- matrix(0, ncol = 8, nrow = length(lista))
colnames (table.fami) <- names(fami)
ntotal <-c(1131,26711,2614,9370,1269,8948,2778,5293)
rownames(table.fami) <- sapply(lista, function (x) {sub(".class.txt", "",x)})
for (i in 1:length(lista)) {
  print(i)
   fami <- rep(0,8)
   names(fami) <- c("mid1", "mid2", "mid3", "mid4", "mid5", "mid6", "mid7", "mid8")
   fam <- read.delim (lista[i], header = F)
   fami[names(table(fam[,4]))] <- table(fam[,4])
   table.fami[i,] <- fami/ntotal
}

layout(matrix(c(1:8), 4, 2))
for ( i in 1 : 8) {
  barplot(t(table.fami)[i,], ylab = "abundance",xlab = "", las = 2)
  title (colnames(table.fami)[i])
}

dev.off()


source ("~/rescale.R")



