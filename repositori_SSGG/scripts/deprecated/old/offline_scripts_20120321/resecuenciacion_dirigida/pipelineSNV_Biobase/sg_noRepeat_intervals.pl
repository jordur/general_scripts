#!/usr/bin/perl

if($ARGV[0] eq "")
{
	die("\n\n\t\tPara que sirve: Parte de un fichero multifasta de sequencias enmascaradas donde N son zonas de baja complejidad. En la cabezera se encuentran los intervalos que correspondes.\n\t\toutput Obtiene los intervalos que no son de baja complejidad y su secuencia:\n\n\tEjemplo:\n>chr1_78398867_78398993_+\nNNNNNNNNNNNNNNNNNNNNNNTCCCAGTTTCTTCAAGAAATAGTCCCTTTTTAGTATGTGTAATTCTGGCCAGAGTGATAAAATAATTATTTTAAATAGGTAGTAGATGATGACTCCCCAGAGATG\n>chr1_78401408_78401470_+\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n>chr1_116268164_116268274_+\nGTAGAGTGGGTCTGGNNNNNNNNNNNNNNNNNNNNNNNNNNNGAGCTTCCTTGAGTAGGGATCACTGTGGCAGACACAGCAGCCACAGCATTATGAAATGTAAAATGAGAT\n\n");
}

my $input=$ARGV[0];

my @explode=();
my @sequence=();

my $entry=0;
my $start=0;
my $end=0;
my $startReal=0;
my $endReal=0;
my $count=0;

my $line="";
my $string="";

open(FILE,"<",$ARGV[0]);
while($line=<FILE>)
{
        chomp($line);

	if($line=~/^>/)
	{
		$start=0;
		$end=0;

		@explode=split("_",$line);
		$start=$explode[1];
		$end=$explode[2];
		if($start==0 && $end==0)
		{
			die("\n\n\t\t**Error en la cabezera, no leido los intervalos correctamente**\n\n$line\n\n");
		}
		$entry=1;
	}
	else
	{
		if($entry==1)
		{
			$string="";
			$count=0;
			$startReal=0;
			$endReal=0;
			@sequence=split("",$line);
			for(my $i=0;$i<scalar(@sequence);$i++)
			{
				if($sequence[$i] eq "N")
				{
					if($string ne "")
					{
						$startReal=$start;
						$endReal=$startReal+$count-1;
						print ">$startReal"."\t".$endReal."\n";
						print "$string\n";
						$string="";
						
						$start=$endReal+1;
						$count=0;
					}
					else
					{
						$start=$start+1;
					}
				}
				else
				{
					$string=$string."".$sequence[$i];
					$count++;
				}
			}
			if($string ne "")
			{
				$startReal=$start;
				$endReal=$end;
                                print ">$startReal"."\t".$endReal."\n";
                                print "$string\n";
			}
		}
		$entry=0;
	}
}

