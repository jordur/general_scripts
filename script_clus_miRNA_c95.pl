#!/usr/bin/perl -w
use strict;

### Script per parsejar un arxiu .clstr de cd-hit amb threshold del 95%. 		###
### OUTPUT: el fitxer cluster-hit.list -> CLusterID	Sequencia representativa 	###
###			el fitxer clustersum -> ClusterID	counts 							###

# my @files=<*>;
# my @files2 = grep (/clusters.clstr/, @files);



open OUT1, ">PPA_cluster-hit.list";
#open OUT2, ">miRNA11_collapsed_hit_100.list";
open OUT3, ">PPA_clustersum";
#foreach my $file(@files2){
#    open IN, $file;
my $countN = 0;
my $countA = 0;
my $clusters = 0;
chomp ($countA);
chomp($countN);

print OUT1 "\n";
print OUT3 "\n";

my $name;
my $clus_name;
my @clusters;
#chomp();


my $in = $ARGV[0];

##### Subrutina uniq #########
sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}

###############
##Main script #
###############

while(<>){
    if (/^>Clus/){
#        open OUT, ">>temp";
	print OUT3 $countN,"\n";
	chomp;
	$name = $_;
	print OUT3 $name,"\t";
	chomp($name);
	$countN = 0;
	$countA = 0;
	my @clusters;
    }else{
#	print OUT3 $name,"\n";
		if(/.PPA(\d+)-(\d+)./){
	    	$countN = $countN + $2;
#	    	print $name,"\n";
		}
		if(/.(>PPA\d+)-(\d+).+\*/){
	    	$clus_name = $1."-".$2;
	    	chomp;
	    	print OUT1 $name."\t".$clus_name,"\n";
		}
#	if ($countN >= 100){
#	    $clus_name=$name."\t".$countN."\n";
#	    print OUT2 $name,"\n";
#	    print $countN, "\n";
#	    push @clusters, $name,"\n";
#	    push @clusters, $clus_name;
#	    print $clus_name,"\n";

#	}
       
    }
}
print "DONE";
#print uniq(@clusters);
close OUT1;
#close OUT2;
close OUT3;
