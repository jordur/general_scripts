#!/usr/bin/perl 

##ESTE SCRIPT ES COMO PASO INTERMEDIO EN EL PROCESO DE ANOTACION PARA PANELES DE GENES(CARDIO ,CANCER etc) ,454 Y EXOMA. CAMBIA UN _ POR UN ESPACIO EN LA COLUMNA 7 DONDE ESTA DEFINIDO EL TIPO DE MUTTACION QUE ES (INTRONIC...). ESTA TRANSFORMACION ES FUNDAMENTAL PARA PROCESOS POSTERIORES EN EL PIPELINE

$input=$ARGV[0];

if($ARGV[0] eq "")
{
	die("\n\n\t\t ** I:variant_effect_output.txt **\n\n");
}


open(FILE,'<',$input);
while($line=<FILE>)
{
        chomp($line);
        @explode=split("\t",$line);
	for($i=0;$i<scalar(@explode);$i++)
	{
		if($i==0)
		{
			$explode[0]=~s/_/\t/gi;
			print "$explode[0]";
		}
		else
		{
			if($i==6 || $i==5 || $i==4)
			{
				$explode[$i]=~s/_/ /g;
			}
			if($i==12)
			{
				if($explode[$i]=~/Condel/)
				{
					$explode[$i]=~s/Condel=/\t/;
					$explode[$i]=~s/\(/ /;
					$explode[$i]=~s/\)//;
				}
				else
				{
					$explode[$i]=~s/_/-/gi;
					$explode[$i]=$explode[$i]."\t-";
				}
			}
			print "\t$explode[$i]";
		}
	}
	print "\n";
}
