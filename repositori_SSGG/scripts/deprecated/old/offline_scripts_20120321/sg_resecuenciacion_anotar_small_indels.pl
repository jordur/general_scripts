#!/usr/bin/perl
sub usage
{
        print "\nPARA QUE SIRVE: Anota el fichero gff de los small indels\n";
        print "Input: 1) Fichero gff con los small indels\n";
        print "Output: 1) Fichero gff modificado con la anotacion \n\n\n";
    exit(1);
}
#chr1    AB_SOLiD Small Indel Tool       insertion_site  12136   12136   1       .       .       ins_len=1;tigh_reads=30;coverage_ratio=0.4000;zygosity=HOMOZYGOUS;zygosity-score=0.0000;bead_ids=1201_1740_1753,1083

if(scalar(@ARGV) == 0)
{
    usage();
}

open(INDELS,"<",@ARGV[0]);
open(GFF,"<",@ARGV[1]);
while  (my $linea=<INDELS>)
{
        chomp($linea);
        push (@array,$linea);
}

close (INDELS);

while  (my $linea2=<GFF>)
{
        chomp($linea2);
        push (@gff,$linea2);
}

close (INDELS);
close (GFF);

#print join ("\t",@array2),"\n";
my $contador=0;
for (my $i=0; $i<= $#array ; $i++)
{
	@split_array = split ("\t", $array[$i]);
	@indels_columna9 = split (";",$split_array[8]);
	for (my $a=0; $a<=$#gff; $a++)
        {

                @split_gff = split (/\t/,$gff[$a]);
		@gff_columna9 = split (";",$split_gff[8]);
		#print  $split_gff[3],"\t",$split_gff[4],"\n";
                if ((($split_array[3] >= $split_gff[3]) && ($split_array[4] <= $split_gff[4]))|| (($split_array[3] >= $split_gff[3]) && ($split_array[3] <= $split_gff[4])) || (($split_array[4] >= $split_gff[3]) && ($split_array[4] <= $split_gff[4]))|| (($split_array[3] <= $split_gff[3]) && ($split_array[4] >= $split_gff[4])))
                {
			$indels_columna9[0]=~s/.*_len=//;
			
			$indels_columna9[3] =~s/no_nonred_reads=//;
                        $indels_columna9[5] =~s/zygosity=//;
			print $split_array[3],"\t",$split_array[4],"\t",$split_array[2],"\t",$indels_columna9[0],"\t",$indels_columna9[3],"\t",$split_gff[3],";",$split_gff[4],";";#$gff_columna9[0],";";
			for (my $o=0; $o <=$#gff_columna9 ; $o++)
			{
				#$indels_columna9[3] =~s/no_nonred_reads=//;
				#$indels_columna9[5] =~s/zygosity=//;
				#print "---------->",$o,"\n";

				if ($gff_columna9[$o] =~m/ID=/)
				{
					$gff_columna9[$o] =~s/ID=//;
					print $gff_columna9[$o],";";
					#$gff_columna9[5]=~s/product=//;
					#$gff_columna9[6]=~s/protein_id=//;
				}
				elsif ($gff_columna9[$o] =~m/product=/)
				{
					$gff_columna9[$o]=~s/product=//;
					print $gff_columna9[$o],";";
				}
				elsif ($gff_columna9[$o] =~m/protein_id=/)
				{
					$gff_columna9[$o]=~s/protein_id=//;
					print $gff_columna9[$o],"\n";
				}
#				print $split_array[3],"\t",$split_array[4],"\t",$split_array[2],"\t",$split_array[5],"\t",$indels_columna9[3],"\t",$indels_columna9[5],"\t",$split_gff[3],";",$split_gff[4],";",$gff_columna9[0],";",$gff_columna9[5],";",$gff_columna9[6],"\n";
				#$contador= $contador+1;
			}
			$contador= $contador+1;
		}
	}
	if ($contador ==0)
	{
		$indels_columna9[0]=~s/.*_len=//;
		$indels_columna9[3] =~s/no_nonred_reads=//;
                $indels_columna9[5] =~s/zygosity=//;
		print $split_array[3],"\t",$split_array[4],"\t",$split_array[2],"\t",$indels_columna9[0],"\t",$indels_columna9[3],"\n";
	}	
	$contador =0;
}
