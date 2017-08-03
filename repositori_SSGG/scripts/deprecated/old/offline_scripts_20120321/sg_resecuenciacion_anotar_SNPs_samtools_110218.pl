#!/usr/bin/perl

##PARA QUE SIRVE: EL PROGRAMA POSEE COMO INPUTS UN FICHERO DONDE SE ENCUENTRA LOS SNPs Y LAS COORDENADAS (EL PROGRAMA ESTA DISEÑADO EN EL FORMATO DONDE LAS CORDENADAS ESTA EN LA SEGUNDA COLUMNA. EL SIGUIENTE INPUT ES UN FICHERO GFF. EL PROGRAMA DETECTA QUE SNPs SE ENCUENTRAN ANOTADOS EN FICHERO GFF COMO CDS,tRNA,rRNA O PSEUDOGENE, CUALQUIER OTRA ANOTACION SERA IGNORADA.
##SALIDA.FICHERO SIMILAR AL INPUT 1 (SNPs) DONDE QUEDARA AÑADIDO EN LA LINEA LA ANOTACION DEL GFF (SI SE HA ENCONTRADO) Y QUEDARA IGUAL SI NO SE HAENCONTRADO UNA ANOTACION.


if($ARGV[0] eq "" ||  $ARGV[1] eq "")
{
	print "\n\n\t\t**INPUTS:\n\n";
	print "\tI:\tInput File 1 \n\tII:\tInput File 2 (gff format)\n\n";
}
else
{
	my $input=$ARGV[0];
	my $gff=$ARGV[1];

	my @SNP_coord=();
	my @SNP_line=();
	my @explode=();
	my @gff_file=();
	my @explode_gff=();

	my $count=0;
	my $exit=0;

	my $line="";
	my $description="";

	open(FILE2,"<",$gff);
 	@gff_file=<FILE2>;

	open(FILE,"<",$input);
	while($line=<FILE>)
	{
		chomp($line);
		@explode=split("\t",$line);

		if($explode[1] ne "")
		{
			$SNP_coord[$count]=$explode[1];
			$SNP_line[$count]=$line;

			$count=$count+1;
		}

	}

	for(my $i=0;$i<scalar(@SNP_coord);$i++)
	{
		$exit=0;
		$description="";

		for(my $k=0;$k<scalar(@gff_file);$k++)
		{
			if($exit==0)
			{
				@explode_gff=split("\t",$gff_file[$k]);
				if($explode_gff[2] eq "CDS" || $explode_gff[2] eq "tRNA" || $explode_gff[2] eq "rRNA" || $explode_gff[2] eq "pseudogene")
				{
					if($SNP_coord[$i]>= $explode_gff[3] && $SNP_coord[$i] <= $explode_gff[4])
					{
						$description=$explode_gff[8];
						$description=~s/;/\t/g;
						chomp($description);
						print "$SNP_line[$i]\t$description\n";
						$exit=1;
					}	
				}
			}
		}
		if($exit==0)
		{
			print "$SNP_line[$i]\n";
		}
	}
}
