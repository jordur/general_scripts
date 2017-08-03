#!/usr/bin/perl
sub usage
{
    print "\nPARA QUE SIRVE: Saca los indels comunes entre los ficheros HapMap y bioscope \n COMO SE USA:\n";
    print "Input: 1) Archivos de HapMap parseados de /data/results/Solid0065/Exoma-HapMap-FA104-FA287/HapMap/indels_HapMap/ \n";
    print "       2) fichero tabulado con los indels en el que la primera columna es el cromosoma\n";
    print " Output: fichero tabulado con los SNPs comunes \n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(HAPMAP,"<",@ARGV[0]);
open(INDELS,"<",@ARGV[1]);
my @bedfile;
my @indels;
my @split_new_indels;
my @split_hapmap;

while ($lines_hapmap=<HAPMAP>)
{
	chomp($lines_hapmap);
	@split_lines_indels=split ("\t",$lines_hapmap);
	$ref_length = length($split_lines_indels[2]);
	$cons_length = length ($split_lines_indels[3]);
	$inicio=$split_lines_indels[1]-2;
	$final = $split_lines_indels[1]+$split_lines_indels[6]+2;
	$final_inser= $split_lines_indels[1]+2;
	
	if ($ref_length > $cons_length)
	{
		push (@hapmap, $split_lines_indels[0]."\tdeletion\t".$inicio."\t".$final."\t".$split_lines_indels[2]."\t".$split_lines_indels[3]."\t".$split_lines_indels[4]."\t".$split_lines_indels[5]."\t".$split_lines_indels[6]."\t".$split_lines_indels[7]."\t".$split_lines_indels[8]."\t".$split_lines_indels[9]);
	}
	else
	{
		push (@hapmap, $split_lines_indels[0]."\tinsertion\t".$inicio."\t".$final_inser."\t".$split_lines_indels[2]."\t".$split_lines_indels[3]."\t".$split_lines_indels[4]."\t".$split_lines_indels[5]."\t".$split_lines_indels[6]."\t".$split_lines_indels[7]."\t".$split_lines_indels[8]."\t".$split_lines_indels[9]);
	}
}


close (HAPMAP);
#print join ("\n",@hapmap),"\n";

#chr11   100147006       C       CT      1       1       2
#chr11   insertion	100147006	100147008       C       CT      1       1       2
#chr1    deletion        864728  864728  1       T/      864726  TCTCC/TCCC      REF,5   4       HEMIZYGOUS      1.0000

my $l = 0;
while (my $indel=<INDELS>)
{
	chomp ($indel);
        @split_new_indels= split (/\t/,$indel);
#	print $split_new_indels[2],"---\n";
	for (my $j=0; $j<=$#hapmap;)
	{
		$j=$l;
		@split_hapmap = split (/\t/,$hapmap[$l]);
		#print $split_new_indels[2]," --- ", $split_new_indels[3]," --- ", $split_hapmap[2]," --- $split_hapmap[3]\n";
		if (($split_new_indels[2] <= $split_hapmap[2] ) && ($split_new_indels[3] >= $split_hapmap[2])|| (($split_new_indels[2] <= $split_hapmap[3]) && ($split_new_indels[3] >= $split_hapmap[3])) || (( $split_new_indels[2] >= $split_hapmap[2]) && ($split_new_indels[3] <= $split_hapmap[3])))
		{
			print $indel,"\n";
			$j++;
			$l=$l+1;
			last;
		}
		elsif($split_new_indels[3] < $split_hapmap[2])
		{	
			last;
		}
		else 
		{
			$j++;
			$l=$l+1;
		}
	}
}
close(INDELS);
