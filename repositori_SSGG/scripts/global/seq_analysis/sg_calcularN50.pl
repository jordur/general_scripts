#!/usr/bin/perl -w
use strict;
my ($len,$total)=(0,0);
my @x;
if (@ARGV == 0) {
   die "
       Please specify a fasta file
\n".localtime()."\n=================================================================\n\n";
}


while(<>){
	if(/^[\>\@]/){
		if($len>0){
			$total+=$len;
			push @x,$len;
		}
		$len=0;
	}
	else{
		s/\s//g;
		$len+=length($_);
	}
}
if ($len>0){
	$total+=$len;
	push @x,$len;
}
@x=sort{$b<=>$a} @x; 
my ($count,$half)=(0,0);
for (my $j=0;$j<@x;$j++){
	$count+=$x[$j];
	if (($count>=$total/2)&&($half==0)){
		print "$x[$j]\n";
		$half=$x[$j]
	}elsif ($count>=$total*0.9){
		print "N90: $x[$j]\n";
		exit;
	}
}
