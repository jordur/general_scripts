#!/usr/bin/perl

#Primero hay que eliminar los espacioes con sed 's/ //gi'

sub usage
{
    print "\nPARA QUE SIRVE: hace el recuento de lecturas por gen a partir del fichero de conteo de Bioscope\n";
    print " Input: 1) fichero con los nombres de los genes únicos (utilizar awk para seleccionar columna 9 y sort -u para coger solo las columnas unicas), 2) fichero con el conteo sin espacios en la última columna , 3) lecturas únicas mapeables \n\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}


open (GEN, "<", @ARGV[0]);
open (CONTEO, "<", @ARGV[1]);
$lecturas_totales = @ARGV[2];
while ($lines=<GEN>)
{
  chomp ($lines);
  push (@gen,$lines);
  
}

close (GEN);

while ($lineas_conteo=<CONTEO>)
{
  chomp ($lineas_conteo);
  push (@conteo,$lineas_conteo);

}

close (CONTEO);


$contador=0;
$longitud_exones_gen=0;
print "GENE_ID\tTRANSCRIPT_ID\tTRANSCRIPT_LENGTH\tCOUNTS\tRPKM\n";
for ($i=0; $i<= $#gen ;$i++)
{
	
	for ($a=0; $a<=$#conteo; $a++)
	{
		@split_conteo = split ("\t", $conteo[$a]);
		if ($split_conteo[8] eq $gen[$i])
		{
			$contador = $contador + $split_conteo[5];
			if ($split_conteo[5]!=0)
			{
				$long_exon= $split_conteo[4]-$split_conteo[3];
				$longitud_exones_gen=$longitud_exones_gen+$long_exon;
			}
		}

	}	
	if ($contador>0)
        {
		#print "GENE_ID\tCOUNTS\n";
		$transcrito = $gen[$i];
		$transcrito =~s/gene_id".*;transcript_id"//;
		$transcrito =~s/";//;
		$gen[$i]=~s/gene_id"//;
		$gen[$i]=~s/";transcript_id".*$//;
		$rpkm = ((1000000000*$contador)/($longitud_exones_gen*$lecturas_totales));
		print $gen[$i],"\t", $transcrito,"\t",$longitud_exones_gen,"\t",$contador,"\t",$rpkm,"\n";
	}
	$longitud_exones_gen=0;
	$contador =0;
}
