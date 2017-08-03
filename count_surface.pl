#! usr/bin/perl -w
use strict;


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
	 print "Numero de superficie:",$count1,"\n";
	 print "Numero internos:",$count0,"\n";
 }