#!/usr/bin/perl
use List::Util qw (sum);

sub usage
{
    print "\nPARA QUE SIRVE: Rescata los valores de calidad de un grupo de lecturas no mapeados\n";
    print " Input: 1) fichero con las lecturas no mapeadas en fasta 2) fichero con los valores de calidad sin cabecera con almohadillas, 3) F3 o R3 \n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}



open (FILE1, "<", @ARGV[0]);# no_mapean

while ($lineas=<FILE1>)
{
	chomp ($lineas);
    	$mayor = substr ($lineas, 0,1);
	if ($mayor eq ">")
	{
		$lineas=~s/>//;
		$lineas=~s/_F3//;
		push (@fichero1,$lineas);
	}
}
close (FILE1);
@fichero1=sort(@fichero1);

open (FILE2, "<", @ARGV[1]);# qual

while ($lines= <FILE2>)
{
	chomp ($lines);
	push (@fichero2,$lines);
}

close (FILE2);
$f3r3= @ARGV[2];
open (OUT1, ">", $f3r3.".qual") or die "No puedo abrir el fichero de salida\n";

                    
use constant VALUE_NOT_FOUND => -1;
use constant ARRAY_UNORDERED => -2;
use constant LIMITS_REVERSED => -3;
use constant VALUE_FOUND     =>  0;

my @list  = @fichero1;
#print @list,"\n";
my $lower = 0;
my $upper = scalar(@list) - 1;
my $found;

for ($i=0; $i<$#fichero2;) 
{ 
  $fichero2[$i]=~s/_F3//;	
  $fichero2[$i]=~s/>//; 
  $valor1 =$fichero2[$i]; 
  $found = binarySearch(\@list, $lower, $upper, $valor1);
  if ($found == VALUE_FOUND) 
  {
	print OUT1 ">",$fichero2[$i],"_",$f3r3,"\n";
	print OUT1 $fichero2[$i+1],"\n";
  } 
  
  $i = $i+2;
}

sub binarySearch () {
        my $array    = shift;
        my $lower    = shift;
        my $upper    = shift;
        my $value    = shift;
        my $range = $upper - $lower;
        if ($range < 0) {
                return LIMITS_REVERSED;
        } elsif ( ($range == 0) && ($${array[$lower]} != $value) ) {
                return VALUE_NOT_FOUND;
        }
        if ($array[$lower] > $${array[$upper]}) {
                return ARRAY_UNORDERED;
        }
        my $center = (int($range / 2)) + $lower;
        if ($value eq $${array[$center]}) {
                return VALUE_FOUND;
        } elsif ($value lt $${array[$center]}) {
                return binarySearch(
                                $array,
                                $lower,
                                $center - 1,
                                $value
                        );
        } else {
                return binarySearch(
                                $array,
                                $center + 1,
                                $upper,
                                $value
                        );
        }
}

close (OUT1);

print "\nOperacion finalizada.\n\n";
