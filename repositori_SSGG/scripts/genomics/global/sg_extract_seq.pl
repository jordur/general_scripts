#!/usr/bin/perl
sub usage {
    print "\nCOMO SE USA:\n";
    print "sg_extract_seq.pl 1 2 .\n";
    print "1) Archivo con las coordenadas separadas por cromosomas ej: nombre_coordenada<TAB>start<TAB>end<TAB>strand.\n";
    print "2) Secuencia de referencia.\n";
    print "\n";
    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

open(EXTRACT,"<",@ARGV[0]);
open(REFERENCIA,"<",@ARGV[1]);
while  (my $miRNA=<EXTRACT>)
{
        chomp($miRNA);
        push (@extraer,$miRNA);
}
while  (my $lineas=<REFERENCIA>)
{
        chomp($lineas);

	$inicio= substr ($lineas,0,1);
	if ($inicio ne ">")
	{
       		$secuencia.=$lineas; 
	}

}

for (my $i=0; $i<=$#extraer; $i++)
{
	@separar = split ("\t",$extraer[$i]);
	$substring = substr ($secuencia,$separar[1]-1,(($separar[2]+1)-$separar[1]));

	if ($separar[3] eq "+")
	{
		print ">".$separar[0],"\n";
		print $substring,"\n";
	}
	elsif ($separar[3] eq "-")
	{
		$substring_rev=reverse ($substring);
		$substring_rev=~ tr/AaTtGgCcNn/TtAaCcGcNn/;
		print ">".$separar[0],"\n";
		print $substring_rev,"\n";
	}
}

close (EXTRACT);
close (REFERENCIA);
