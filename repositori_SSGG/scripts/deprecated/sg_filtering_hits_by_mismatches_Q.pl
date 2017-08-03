#!/usr/bin/perl -w

use strict;

##PARA QUE SIRVE: SIRVE PARA ELIMINAR LAS LECTURAS QUE NO CUMPLAN CON LOS REQUISITOS DE MISMATCHES Y CALIDAD IMPUESTOS POR EL USUARIO
##INPUTS: FICHERO .csfasta O .qual
##OUTPUT: FICHERO .csfasta O .qual DONDE SE HA ELIMINADO LOS ALINEAMINETOS DE LOS READS QUE NO HAN CUMPLIDO LOS REQUISITOS IMPUESTOS POR EL USUARIO


if ($ARGV[0] eq "" || $ARGV[1] eq "" || $ARGV[2] eq "")
{
        print"\n\n **INPUTS: **\n\t\tI:File in .csfasta or .qual format \n\t\tII:Mismatches Threshold (in this case always is >, not >=)\n\t\tIII:Quality Threshold (idem II)\n\n";
}
else
{
	my $input=$ARGV[0];	
	my $mismatches=$ARGV[1];
	my $quality=$ARGV[2];	

	my $line="";
	
	my @explode=();
	my @explode2=();
	my @explode3=();

	open(FILE,'<',$input);
	while($line=<FILE>)
	{
		chomp($line);

		if($line!~/^#/)
		{
			if($line=~/^>/)
			{
				@explode=split(",",$line);
				
				print "$explode[0]";

				for(my $i=1;$i<scalar(@explode);$i++)
				{
					if($explode[$i] ne "")
					{
						@explode2=split(":",$explode[$i]);

						$explode2[1]=~s/\./&/gi;

						@explode3=split("&",$explode2[1]);
						
						if($explode3[1]<=$mismatches)
						{
							$explode2[2]=~s/q//gi;
						
							if($explode2[2]>=$quality)
							{
								print ",$explode[$i]"
							}
						}
					}
				}
				print "\n";

			}
			else
			{
				print "$line\n";
			}
		}
		else
		{
			print "$line\n";
		}
	}
}

