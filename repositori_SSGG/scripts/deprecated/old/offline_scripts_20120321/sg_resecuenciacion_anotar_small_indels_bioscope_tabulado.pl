#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Anota los SNPs encontrados \n COMO SE USA:\n";
    print " Input: 1) fichero gff con la anotación de los CDS \n";
    print "        2) fichero InDels\n
	chr1    deletion        38650   38659   10      CGGGGTCAACGAA/CGG/TAA   38649   TCGGGGTCAACGAACG/TCGGCG/TTAACG  REF,1,1 2   HEMIZYGOUS       0.9387
chr1    deletion        63248   63248   1       G/      63246   GTGGTC/GTGTC/NO_CALL    REF,1,1 2       HEMIZYGOUS      0.9576
chr1    deletion        84607   84613   7       ACCGGTGGTCTA/CAAGG      84605   GCACCGGTGGTCTACGG/GCCAAGGCGG/NO_CALL    REF,13,2     15      HOMOZYGOUS      0.0000";
    print "\nOutput: fichero InDels  con una columna anadida que indica el gen en el que se encuentra  \n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(FILE1,"<",@ARGV[0]);
open(FILE2,"<",@ARGV[1]);
my @fichero1;
my @fichero2;
my @fila1;
my @fila2;
for (my $i=1; my $lines=<FILE1>; $i++)
{
  chomp($lines);
  push (@fichero1,$lines);
}

close (FILE1);


for (my $i=1; my $lines2=<FILE2>; $i++)
{
  chomp($lines2);
  push (@fichero2,$lines2);
}

close (FILE2);

my $contador=0;
for (my $i=0; $i<=$#fichero2; $i++)
{
        @fila2 = split (/\t/,$fichero2[$i]);
        for (my $i=0; $i<=$#fichero1; $i++)
        {

                @fila1 = split (/\t/,$fichero1[$i]);
                if (((@fila2[2] >= $fila1[3]) || (@fila2[3] >= $fila1[3])) && ((@fila2[2] <= $fila1[4]) || (@fila2[3] <= $fila1[3])))
                {
                        #@Gene_ID=split (";",$fila1[8]);
                        #$fila1[8]=~s/eC_number=.*//;
                        $fila1[8]=~s/Parent=.*product=/product=/;
                        $fila1[8]=~s/;transl_table=.*$//;
                        @Gene_ID=split (";",$fila1[8]);
                        print $fila2[0],"\t",$fila2[1],"\t",$fila2[2],"\t",$fila2[3],"\t",$fila2[4],"\t",$fila2[5],"\t",$fila2[6],"\t",$fila2[7],"\t",$fila2[8],"\t",$fila2[9],"\t",$fila2[10],"\t",$fila2[11],"\t",$fila1[3],"\t",$fila1[4],"\t",$Gene_ID[0],"\t",$Gene_ID[1],"\t",$Gene_ID[2],"\n";
			$contador=$contador+1;
                }
        }
	if ($contador == 0 )
	{
		print $fichero2[$i],"\n";
	}
	$contador=0;
}

