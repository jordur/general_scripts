#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Extrae los SNPs de un fichero tabulado que se encuentran en las regiones indicadas por un archivo donde indicamos los intervalos de secuencias que estamos analizando \n COMO SE USA:\n";
    print "Input: 1) fichero tabulado con las regiones que estudiamos. La primera columna es el cromosoma, la segunda la coordenada de inicio y la tercera, la coordenada final. Debe estar ordenado de menor a mayor.\nejemplo: chr1    343432      343445\n";
    print "       2) fichero tabulado de SNPs or Indels en formato VCF\n";
    print " Output: fichero tabulado con los SNPs que se encuentran en las regiones indicadas en el fichero de input 1. Hay que hacerlo cromosoma por cromosoma y darle un fichero de salida. \n";
    print "\n";
    print "\n\n\t** Se debe indicar si se trata de un vcf para SNPs (\"snps\") o para Indels (\"indels\")**\t\n\n\n";

    exit(1);
}

if(scalar(@ARGV) < 2)
{
    usage();
}

if($ARGV[2] ne "snps" && $ARGV[2] ne "indels")
{
	usage();
}


my %hash=();

my @explode1=();
my @explode2=();
my @finalColumn=();
my @f=();

my $line="";
my $line2="";
my $chr="";

my $position=0;
my $start=0;
my $end=0;
my $scalar=0;
my $count=0;

my $identification="";

open(FILE1,"<",@ARGV[0]);
@f=<FILE1>;

if($ARGV[2] eq "snps")
{
	open(FILE2,"<",@ARGV[1]);
	while($line2=<FILE2>)
	{
		chomp($line2);
	
		if($line2!~/^#/)
		{
			@explode2=split("\t",$line2);
			$chr=$explode2[0];
			$chr=~s/chr//g;
			$position=$explode2[1];

                        $scalar=(scalar(@explode2)-1);
                        @finalColumn=split(":",$explode2[$scalar]);
                        $finalColumn[0]=~s/\W//g;

			if($finalColumn[0] ne "00")
			{
				$identification=$chr."".$position;

	
				for($k=$count;$k<scalar(@f);$k++)
				{
 					chomp($f[$k]);
	
					@explode=split("\t",$f[$k]);
					$explode[0]=~s/chr//g;
					if($explode[0] eq $chr)
					{
						if($position>=$explode[1] && $position <=$explode[2])
						{
							if($hash{$identification}!=1)
							{
								print "$line2\n";
								$count=$k;
								$hash{$identification}=1;
							}
						}	
					}
				}
			}
		}
	}
}
if($ARGV[2] eq "indels")
{
        open(FILE2,"<",@ARGV[1]);
        while($line2=<FILE2>)
        {
                chomp($line2);
                if($line2!~/^#/)
                {
			$start=0;
			$end=0;

                        @explode2=split("\t",$line2);
			$chr=$explode2[0];
			$chr=~s/chr//g;
			
			$scalar=(scalar(@explode2)-1);
			@finalColumn=split(":",$explode2[$scalar]);
			$finalColumn[0]=~s/\W//g;

			if($finalColumn[0] ne "00")
			{
				$chr=~s/chr//gi;

				if((length($explode2[3]) - length($explode2[4]))<0)
				{
					$start=$explode2[1];
					$end=($start+length($explode2[4]))-1;

				}
				if((length($explode2[3]) - length($explode2[4]))>0)
				{
					$start=$explode2[1];
					$end=($explode2[1]+length($explode2[3]))-1;
				}

				$identification=$chr."".$start."".$end;
   	        	        for($k=$count;$k<scalar(@f);$k++)
				{                        
               		              	chomp($f[$k]);

                                       	@explode=split("\t",$f[$k]);
					$explode[0]=~s/chr//gi;
					if($explode[0] eq $chr)
					{
						if($hash{$identification}!=1)
						{
							if($explode[1] <= $start && $explode[2] >=$end)
							{
								$count=$k;
								print "$line2\n";
								$hash{$identification}=1;
								last;
							}
							elsif($explode[1] >= $start && $explode[1] <= $end && $explode[2]>=$start && $explode[2]>=$end)
							{
								$count=$k;
								print "$line2\n";
								$hash{$identification}=1;
								last;
							}
							elsif($explode[1]<=$start && $explode[1]<=$end && $explode[2]>= $start &&  $explode[2] <=$end)
							{
								$count=$k;
								print "$line2\n";
								$hash{$identification}=1;
								last;
							}
							elsif($explode[1]>=$start && $explode[2]<=$end)
							{
								$count=$k;
								print "$line2\n";
								$hash{$identification}=1;
								last;
							}
						}
					}
				}
			}
		}
	}
}
