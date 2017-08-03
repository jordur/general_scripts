#!/usr/bin/perl -w
use strict;

### Script per obtindre informació de seqüències úniques i compartides a partir del .clstr de cd-hit ### 

if (@ARGV != 1){
    print "Arguments: arxiu .clstr resultant de cd-hit\n";
}else{

my $countN = 0;
my $countA = 0;
chomp ($countA);
chomp($countN);

my $name;
my $name1;
my $name2;
my @newborn;
my $total_newborn = 0;
my @adult;
my $total_adult = 0;
my@shared;
my $total_shared = 0;
my @maxA;
my @maxN;
my $n=0;
my $a=0;
chomp;

my $in = $ARGV[0];

##### Subrutina uniq #########
sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}

###############
##Main script #
###############
# open OUT1, ">clusters_composition.txt";
# open OUT2, ">clusters_summary.txt";
# open OUT3, ">PPN_cluster-PPA_cluster.txt";
open OUT4, ">uniq_newborn_cluster.txt";
open OUT5, ">uniq_adult_cluster.txt";

#print OUT1 "ClusterID\t\tAdult counts\tNewborn counts\n";

while(<>){
    if (/^>Cluster (\d+)/){
	$name = $_;
	chomp();
	my $cluster = $1 - 1;
	if ($cluster >= 0){
#	    	print "Total counts Adult Cluster $cluster:".$a,"\n";
#	    	print "Total counts Newborn Cluster $cluster:".$n,"\n";
#	    	print "Cluster $cluster:\n";
	    push @maxN, $countN;
#	    	print @maxN,"\n";
	    push @maxA, $countA;
#	    	print @maxA,"\n";
	    # 	print "Cluster $cluster:\n";
	    # 	print "\tNewborn clusters: $countN\n";
	    # 	print "\tAdult clusters: $countA\n\n";  
	    
	    if($maxN[0] >= 1 && $maxA[0] == 0 ){  ### si el nombre maxim d'adults == 0, no n'hi ha. Per tant el cluster és newborn
#	     		push @newborn, $cluster;
	  		# print OUT1 "Cluster"." ".$cluster,"\t";
	  		# print OUT1 "newborn\n";
	 		print OUT4 $name1,"\n";
	 		$total_newborn = $total_newborn +1;
		 
	    }elsif($maxN[0] == 0 && $maxA[0] >= 1){  ### De la mateixa manera, si si el nombre maxim de newborn == 0 el cluster és adult
	 		# print OUT1 "Cluster"." ".$cluster,"\t";
	 		# print OUT1 "adult\n";
	 		$total_adult = $total_adult +1;
	 		print OUT5 $name2,"\n";

	    }elsif($maxN[0] >= 1 && $maxA[0] >=1){ ###
	 		# print OUT1 "Cluster"." ".$cluster,"\t";
	 		# print OUT1 "shared\t";
	 		# print OUT3 $name1,"\t",$name2,"\n";
	 		$total_shared = $total_shared +1;
	 		# print OUT1 $a,"\t";
	 		# print OUT1 $n,"\n";
	    }
	}
	$countN = 0;
	$countA = 0;
	@maxA =();
	@maxN = ();
	$n =0;
	$a=0;
    }else{
#		if(/.>lcl\|PPN(\d+)-(\d+)./){
        if (/.>PPN(\d+)-(\d+)./){
	    $name1="PPN".$1."-".$2;
	    $countN = $countN + 1;
#	    	$n = $n + $1;
#	    	print $1,"\n";
#	    	push @maxN, $countN;
#	    	print $countN,"\n";
#		}elsif(/.>lcl\|PPA(\d+)-(\d+)./){
	}elsif(/.>PPA(\d+)-(\d+)./){
	    $name2="PPA".$1."-".$2;
	    $countA = $countA + 1;
#	    	$a = $a + $1;
#	    	print $1,"\n";
#	    	push @maxA, $countA;
#	    	print $countA,"\n";
	}
    }
}

# print OUT2 "total newborn clusters: $total_newborn\n";
# print OUT2 "total adult clusters: $total_adult\n";
# print OUT2 "total shared clusters: $total_shared\n";
# close OUT1;
# close OUT2;
# close OUT3;
close OUT4;
close OUT5;
}
