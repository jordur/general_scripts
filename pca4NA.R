#####################################################################
# This function fills in NAs by PCA				       	  #
#											  #
# Developed by Ana Conesa aconesa@cipf.es, 14-Aug-2007		  #
#											  #
#####################################################################


pca4NA <- function (X, fac = min(dim(X)), conver = 1e-07, max.iter = 1000) 
{
   if (any(is.na(X))) { 
	X <- as.matrix(X)
   	print("Missing values are taken care of by PCA")
	NA.position <- which(is.na(X)) # finds NAs of X
   	print(paste(length(NA.position), "missing values"))
   	NA.values <- rnorm(length(NA.position), mean(X, na.rm = T), 
                      sd(X, na.rm = T)) # creates random values
   	X[NA.position] <- NA.values # put random values in X
   	SST <- 0
   	for (it in 1:max.iter) { # start loop
     	 	SST.old <- SST
   		pca <- PCA.GENES(t(X)) # fits PCA
   		T <- as.matrix(pca$scores[,1:fac]) # scores 
   		P <- as.matrix(pca$loadings[,1:fac]) # loadings
   		Xe <- T %*% t(P) # model
		Xe <- t(Xe)
      	X[NA.position] <- Xe[NA.position] # put estimated values in X
      	SST <- sum(X^2)
   		# Convergence
   		if (abs(SST-SST.old)/SST < conver) break
   	}
   } else { print("No missing values, unchanged output") }
   return(X) # output
}