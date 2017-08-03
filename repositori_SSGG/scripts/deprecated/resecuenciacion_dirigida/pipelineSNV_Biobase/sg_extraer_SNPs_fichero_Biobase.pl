#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Extrae los SNPs de un fichero tabulado que se encuentran en las regiones indicadas por un archivo donde indicamos los intervalos de secuencias que estamos analizando \n COMO SE USA:\n";
    print "Input: 1) fichero tabulado con las regiones que estudiamos. La primera columna es el cromosoma, la segunda la coordenada de inicio y la tercera, la coordenada final. Debe estar ordenado de menor a mayor.\nejemplo: chr1    343432      343445\n";
    print "       2) Biobase por cromosoma\n";
    print " Output: fichero tabulado con los SNPs que se encuentran en las regiones indicadas en el fichero de input 1. Hay que hacerlo cromosoma por cromosoma y darle un fichero de salida. \n";
    print "\n";
    

    exit(1);
}

if(scalar(@ARGV) < 2)
{
    usage();
}

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

open(FILE2,"<",@ARGV[1]);
while($line2=<FILE2>)
{
	chomp($line2);
        
	$start=0;
	$end=0;

	@explode2=split("\t",$line2);

	if($explode[4] ne "-")
	{

       		$explode2[4]=~s/chr[A-Z]://;
		$explode2[4]=~s/chr[0-9]{1,2}://g;
		@coor=split("-",$explode2[4]);	

		$start=$coor[0];
		$end=$coor[1];

  		for($k=$count;$k<scalar(@f);$k++)
		{                        
			chomp($f[$k]);
	
       	        	@explode=split("\t",$f[$k]);




			if($explode[1] <= $start && $explode[2] >=$end)
			{
#				$count=$k;
				print "$line2\n";
				last;
			}
			elsif($explode[1] >= $start && $explode[1] <= $end && $explode[2]>=$start && $explode[2]>=$end)
			{
#				$count=$k;
				print "$line2\n";
				last;
			}
			elsif($explode[1]<=$start && $explode[1]<=$end && $explode[2]>= $start &&  $explode[2] <=$end)
			{
#				$count=$k;
				print "$line2\n";
				last;
			}
			elsif($explode[1]>=$start && $explode[2]<=$end)
			{
#				$count=$k;
				print "$line2\n";
				last;
			}
		}
	}
}
