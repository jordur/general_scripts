## FILE must be a one-line column contig length ##
perl -ne 'chomp(); push(@contigs,-s);+=-s;END{foreach(sort{<=>}@contigs){+=-s;=-s;if(>=*0.5){print TOTAL: nN50 : n;exit;} ;}}' FILE
