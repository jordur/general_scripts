awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' input.fastq 
# calculate mean read length in a fastq file
