splitfasta <- function (file, by=6420) {
      name <-sub(".fasta","",substitute(file))
	fasta <- readLines (file)
	#fasta <- fasta[1:1000]
	num.seq <- grep(">",fasta)
	principio <- seq(1,length(num.seq), by=by)
	start <- num.seq[principio]
	end <- num.seq[principio]-1
	end <- end[-1]
	end <- c(end, length(fasta))
	for (i in 1:length(end)) {
		writeLines(fasta[start[i]:end[i]], con=paste(name,"_", i,".fasta",sep=""))
	}
}
