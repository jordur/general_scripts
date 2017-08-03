#!/usr/bin/perl
use List::Util qw (sum);
#use Smart::Comments '###';

#Inputs:
#F3_R3.mates
#[Respuesta Si vamos a procesar un fichero F3 o un R3]
#Solid0065_20080722_1_MPControl_BACs_Wt_P2_Wt_MP_2to3Kb_005_F3.unique.csfasta.ma.25.2
#o
#Solid0065_20080722_1_MPControl_BACs_Wt_P2_Wt_MP_2to3Kb_005_R3.unique.csfasta.ma.25.2
#Outputs:
#F3_no_incluidas.unique.csfasta.ma.25.2, con las parejas nombre-lectura de F3 que no estén en F3_R3
#o
#R3_no_incluidas.unique.csfasta.ma.25.2, con las parejas nombre-lectura de R3 que no estén en F3_R3

print "\nCategoria F3 o R3 [F3 por defecto]: ";
$categoria = <STDIN>;
chop ($categoria);

if ($categoria eq "R3")
{
  $categoria = "R3";  
}else
{
  $categoria = "F3";
} 

open (FILE1, "<", @ARGV[0]);# F3_R3

print "\nCargando el fichero  <".$ARGV[0]."> en memoria ...\n";

for (my $i=1; my $lines=<FILE1>;$i++)
{	
  $invalidos = index($lines,"#");
  if (($invalidos == -1) and (length($lines) > 2))#Se descartan las lineas comentadas y de longitud menor a 2 (sup. vacias).
  {
    @fila = split (/\s+/,$lines);
    push (@fichero1,$fila[0]);
  }
}
close (FILE1);

@fichero1=sort(@fichero1);

open (FILE2, "<", @ARGV[1]);# F3 ó R3

print "\nCargando el fichero  <".$ARGV[1]."> en memoria ...\n";

for (my $i=1; my $lines=<FILE2>;$i++) { 
  $invalidos = index($lines,"#");
  if (($invalidos == -1) and (length($lines) > 2))#Se descartan las lineas comentadas y de longitud menor a 2 (sup. vacias). 
  {
    push (@fichero2,$lines);
  }
}

close (FILE2);

open (OUT1, ">", $categoria."_individual.unique.csfasta.ma.25.2") or die "No puedo abrir el fichero de salida\n";

print "\nCotejando ".$categoria." ...\n";
                    
use constant VALUE_NOT_FOUND => -1;
use constant ARRAY_UNORDERED => -2;
use constant LIMITS_REVERSED => -3;
use constant VALUE_FOUND     =>  0;

my @list  = @fichero1;
my $lower = 0;
my $upper = scalar(@list) - 1;
my $found;

my $i = 0;
for (@fichero2) { ### |...   | [%]
  @fila = split (/\s+/,$fichero2[$i]);
  $valor1 = $fila[0];

  $posicion = index($valor1,"_".$categoria);
  if ($posicion != -1)
  {
    $valor1 = substr($valor1,1,$posicion-1);
  }

  $found = binarySearch(\@list, $lower, $upper, $valor1);

  if ($found != VALUE_FOUND) 
  {
    print OUT1 $fichero2[$i];
    print OUT1 $fichero2[$i+1];
  } 
  
  $i = $i+2;
}#for

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
