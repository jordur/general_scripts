#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Separa las lecturas mapeadas por cromosomas y las pone en el formado que necesita cal_miRNA_extract_seq.pl \n COMO SE USA:\n";
    print " Input: 1) fichero ma.csfasta  \n";
    print " Output: fichero tabulado que contiene el nombre de la secuencia, inicio, fin y orientacion  \n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open (FILE1, "<", @ARGV[0]);
open (CHR1, ">", 'CHR1.ma.tab') or die "No puedo abrir CHR1.ma.tab.\n";
open (CHR2, ">", 'CHR2.ma.tab') or die "No puedo abrir CHR2.ma.tab.\n";
open (CHR3, ">", 'CHR3.ma.tab') or die "No puedo abrir CHR3.ma.tab.\n";
open (CHR4, ">", 'CHR4.ma.tab') or die "No puedo abrir CHR4.ma.tab.\n";
open (CHR5, ">", 'CHR5.ma.tab') or die "No puedo abrir CHR5.ma.tab.\n";
open (CHR6, ">", 'CHR6.ma.tab') or die "No puedo abrir CHR6.ma.tab.\n";

my $sub;
my $pos;
my $sufijo;
my $split_linea;

while (my $lines=<FILE1>)
{		
	$sub = substr($lines,0,1);     
    	if ($sub eq ">")
    	{
      		$lines =~ s/>//;
		@splitpunto = split (/\./,$lines);
		@split_linea = split (",",$lines);
		$pos = index($lines,",");
		$pos_punto = index($lines,".");
		$posicion = substr($lines,$pos+3,($pos_punto-($pos+3)));
		$posicion_negativo = substr($lines,$pos+4,($pos_punto-($pos+4)));
		$sufijo = substr($lines,$pos+1,1);
		$orientacion = substr($lines,$pos+3,1);
            	if ($sufijo eq "1")          
            	{
              		#print CHR1 $split_linea[0],"_",$orientacion,$posicion,"\t";
			if ($orientacion eq "-")
			{
				print CHR1 $split_linea[0],"_",$orientacion,$posicion_negativo-$splitpunto[2],"_",$posicion_negativo,"\t";
				print CHR1 $posicion_negativo-$splitpunto[2],"\t";
				print CHR1 $posicion_negativo,"\t";
				print CHR1 "-","\n";
			}
			else
			{
				print CHR1 $split_linea[0],"_",$posicion,"_",$posicion+$splitpunto[2],"\t";
				print CHR1 $posicion,"\t";
				print CHR1 $posicion+$splitpunto[2],"\t";
				print CHR1 "+","\n";
			}
            	}
            	elsif ($sufijo eq "2")
            	{
			if ($orientacion eq "-")
                        {
                                print CHR2 $split_linea[0],"_",$orientacion,$posicion_negativo-$splitpunto[2],"_",$posicion_negativo,"\t";
				print CHR2 $posicion_negativo-$splitpunto[2],"\t";
                                print CHR2 $posicion_negativo,"\t";
                                print CHR2 "-","\n";
                        }
                        else
                        {
                                print CHR2 $split_linea[0],"_",$posicion."_",$posicion+$splitpunto[2],"\t";
				print CHR2 $posicion,"\t";
                                print CHR2 $posicion+$splitpunto[2],"\t";
                                print CHR2 "+","\n";
			}
            	}
            	elsif ($sufijo eq "3")
            	{
			if ($orientacion eq "-")
                        {
                                print CHR3 $split_linea[0],"_",$orientacion,$posicion_negativo-$splitpunto[2],"_",$posicion_negativo,"\t";
				print CHR3 $posicion_negativo-$splitpunto[2],"\t";
                                print CHR3 $posicion_negativo,"\t";
                                print CHR3 "-","\n";
                        }
                        else
                        {
                                print CHR3 $split_linea[0],"_",$posicion,"_",$posicion+$splitpunto[2],"\t";
				print CHR3 $posicion,"\t";
                                print CHR3 $posicion+$splitpunto[2],"\t";
                                print CHR3 "+","\n";
                        }

            	}
	        elsif ($sufijo eq "4")
            	{
			if ($orientacion eq "-")
                        {
                                print CHR4 $split_linea[0],"_",$orientacion,$posicion_negativo-$splitpunto[2],"_",$posicion_negativo,"\t";
				print CHR4 $posicion_negativo-$splitpunto[2],"\t";
                                print CHR4 $posicion_negativo,"\t";
                                print CHR4 "-","\n";
                        }
                        else
                        {
                                print CHR4 $split_linea[0],"_",$posicion,"_",$posicion+$splitpunto[2],"\t";
				print CHR4 $posicion,"\t";
                                print CHR4 $posicion+$splitpunto[2],"\t";
                                print CHR4 "+","\n";
                        }

            	}
		elsif ($sufijo eq "5")
                {
			if ($orientacion eq "-")
                        {
                                print CHR5 $split_linea[0],"_",$orientacion,$posicion_negativo-$splitpunto[2],"_",$posicion_negativo,"\t";
				print CHR5 $posicion_negativo-$splitpunto[2],"\t";
                                print CHR5 $posicion_negativo,"\t";
                                print CHR5 "-","\n";
                        }
                        else
                        {
                                print CHR5 $split_linea[0],"_",$posicion,"_",$posicion+$splitpunto[2],"\t";
				print CHR5 $posicion,"\t";
                                print CHR5 $posicion+$splitpunto[2],"\t";
                                print CHR5 "+","\n";
                        }

                }
		elsif ($sufijo eq "6")
                {
			if ($orientacion eq "-")
                        {
                                print CHR6 $split_linea[0],"_",$orientacion,$posicion_negativo-$splitpunto[2],"_",$posicion_negativo,"\t";
				print CHR6 $posicion_negativo-$splitpunto[2],"\t";
                                print CHR6 $posicion_negativo,"\t";
                                print CHR6 "-","\n";
                        }
                        else
                        {
                                print CHR6 $split_linea[0],"_",$posicion_negativo-$splitpunto[2],"_",$posicion_negativo,"\t";
				print CHR6 $posicion,"\t";
                                print CHR6 $posicion+$splitpunto[2],"\t";
                                print CHR6 "+","\n";
                        }

                }
    	}
}

close (FILE1);
close (CHR1);
close (CHR2);
close (CHR3);
close (CHR4);
close (CHR5);
close (CHR6);
