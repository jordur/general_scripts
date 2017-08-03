#! usr/bin/perl -w
use strict;

#############Input sub-brak.pl#########
# 0=variable
# 1=no variable
#####PLA2-mid1
#~ my $seq1="1001111111111111111111111111111111111111111111111111101100001010100001010110000100111111001101010000010001000000010000001";
while (my $seq1=<>){
 my @pos=split ("",$seq1);
 my $count0=0;
 my $count1=0;
 foreach my $a(@pos){
	 if ($a =~ /0/){
		 $count0++;
	 }else{
		 $count1++;
	 }
 }
	 print "Numero de residuos conservados:",$count1,"\n";
	 print "Numero de residuos variables:",$count0,"\n";
 }