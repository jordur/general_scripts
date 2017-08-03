#!/usr/bin/perl

sub usage 
{
    print "\nCOMO SE USA:\n";
    print "EL ORDEN ES MUY IMPORTANTE \n\t ";
    print "Primero el fichero que contiene el conteo de cada secuencia individual  \n\t ";
    print "SEGUNDO FICHERO .gff con las secuencias en base-space> \n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(MA,"<",@ARGV[0]);
open(GFF,"<",@ARGV[1]);
#open(FASTAMIRNA,">", 'miRNA_igual_ref.fasta') or die "miRNA_igual_ref.fasta.\n";
open (CHR1, ">", 'chr1_extractseq') or die "No puedo abrir chr1_extractseq.\n";
open (CHR2, ">", 'chr2_extractseq') or die "No puedo abrir chr2_extractseq.\n";
open (CHR3, ">", 'chr3_extractseq') or die "No puedo abrir chr3_extractseq.\n";
open (CHR4, ">", 'chr4_extractseq') or die "No puedo abrir chr4_extractseq.\n";
open (CHR5, ">", 'chr5_extractseq') or die "No puedo abrir chr5_extractseq.\n";
open (CHR6, ">", 'chr6_extractseq') or die "No puedo abrir chr6_extractseq.\n";

my @ma;
my @gff;
my @ma_linea;
my @gff_linea;
while  (my $linea=<MA>)
{
  	chomp($linea);
  	push (@ma,$linea);
}

close (MA);

while ( my $linea2=<GFF>)
{
  chomp($linea2);
  push (@gff,$linea2);
}

close (GFF);


for (my $i=0; $i<=$#ma;$i=$i+2)
{
	$ma[$i]=~s />//;
	for (my $a=0; $a<=$#gff;) #recorro fichero gff
	{
		@splitgffportabas= split ("\t",$gff[$a]); #divido gff por tabuladores
		@splitgff= split ("b=",$gff[$a]); #divido gff por b= para coger la secuencia y poder compararla con el fichero fasta
		$secuenciagff = substr ($splitgff[1],0,length ($ma[$i+1])); #cojo substring del gff para compararla con la del fichero fasta
		$splitgffportabas[8]=~s/.*i=/i=/;
        	$splitgffportabas[8]=~s/;.*$/;/;
        	$poschr = index($splitgffportabas[8],"=");
        	$chr= substr ($splitgffportabas[8],$poschr+1,1);
         	
		if ($secuenciagff ne $ma[$i+1])
		{	
			$a++;
		}
		
		elsif ($secuenciagff eq $ma[$i+1])
		{
			$a++;		
			if ($splitgffportabas[6] eq "+")
			{
				if ($chr eq "1")
				{
					
					print CHR1 $ma[$i],"\t",$splitgffportabas[3],"\t",$splitgffportabas[3]+length($ma[$i+1]),"\t",$splitgffportabas[6],"\n";
					last;
				}
				elsif ($chr eq "2")
				{
					print CHR2 $ma[$i],"\t",$splitgffportabas[3],"\t",$splitgffportabas[3]+length($ma[$i+1]),"\t",$splitgffportabas[6],"\n";
					last;
				}
				elsif ($chr eq "3")
                                {
                                        print CHR3 $ma[$i],"\t",$splitgffportabas[3],"\t",$splitgffportabas[3]+length($ma[$i+1]),"\t",$splitgffportabas[6],"\n";
                                        last;
                                }
				elsif ($chr eq "4")
                                {
                                        print CHR4 $ma[$i],"\t",$splitgffportabas[3],"\t",$splitgffportabas[3]+length($ma[$i+1]),"\t",$splitgffportabas[6],"\n";
                                        last;
                                }
				elsif ($chr eq "5")
                                {
                                        print CHR5 $ma[$i],"\t",$splitgffportabas[3],"\t",$splitgffportabas[3]+length($ma[$i+1]),"\t",$splitgffportabas[6],"\n";
                                        last;
                                }
				elsif ($chr eq "6")
                                {
                                        print CHR6 $ma[$i],"\t",$splitgffportabas[3],"\t",$splitgffportabas[3]+length($ma[$i+1]),"\t",$splitgffportabas[6],"\n";
                                        last;
                                }
			}
			elsif ($splitgffportabas[6] eq "-")
			{
				if ($chr eq "1")
				{
					print CHR1 $ma[$i],"\t",$splitgffportabas[4]-length($ma[$i+1]),"\t",$splitgffportabas[4],"\t",$splitgffportabas[6],"\n";
					last;
				}
				elsif ($chr eq "2")
				{
					print CHR2 $ma[$i],"\t",$splitgffportabas[4]-length($ma[$i+1]),"\t",$splitgffportabas[4],"\t",$splitgffportabas[6],"\n";
					last;
				}
				elsif ($chr eq "3")
                                {
                                        print CHR3 $ma[$i],"\t",$splitgffportabas[4]-length($ma[$i+1]),"\t",$splitgffportabas[4],"\t",$splitgffportabas[6],"\n";
                                        last;
                                }
				elsif ($chr eq "4")
                                {
                                        print CHR4 $ma[$i],"\t",$splitgffportabas[4]-length($ma[$i+1]),"\t",$splitgffportabas[4],"\t",$splitgffportabas[6],"\n";
                                        last;
                                }
				elsif ($chr eq "5")
                                {
                                        print CHR5 $ma[$i],"\t",$splitgffportabas[4]-length($ma[$i+1]),"\t",$splitgffportabas[4],"\t",$splitgffportabas[6],"\n";
                                        last;
                                }
				elsif ($chr eq "6")
                                {
                                        print CHR6 $ma[$i],"\t",$splitgffportabas[4]-length($ma[$i+1]),"\t",$splitgffportabas[4],"\t",$splitgffportabas[6],"\n";
                                        last;
                                }
			}
		}
	}
}

close (CHR1);
close (CHR2);
close (CHR3);
close (CHR4);
close (CHR5);
close (CHR6);
