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
open (CHR7, ">", 'CHR7.ma.tab') or die "No puedo abrir CHR7.ma.tab.\n";
open (CHR8, ">", 'CHR8.ma.tab') or die "No puedo abrir CHR8.ma.tab.\n";
open (CHR9, ">", 'CHR9.ma.tab') or die "No puedo abrir CHR9.ma.tab.\n";
open (CHR10, ">", 'CHR10.ma.tab') or die "No puedo abrir CHR10.ma.tab.\n";
open (CHR11, ">", 'CHR11.ma.tab') or die "No puedo abrir CHR11.ma.tab.\n";
open (CHR12, ">", 'CHR12.ma.tab') or die "No puedo abrir CHR12.ma.tab.\n";
open (CHR13, ">", 'CHR13.ma.tab') or die "No puedo abrir CHR13.ma.tab.\n";
open (CHR14, ">", 'CHR14.ma.tab') or die "No puedo abrir CHR14.ma.tab.\n";
open (CHR15, ">", 'CHR15.ma.tab') or die "No puedo abrir CHR15.ma.tab.\n";


my $sub;
my $pos;
my $sufijo;
my $split_linea;
#>1_145_325_F3,12_19791085.2.20
while (my $lines=<FILE1>)
{		
	$sub = substr($lines,0,1);     
    	if ($sub eq ">")
    	{
      		$lines =~ s/>//;
		@splitpunto = split (/\./,$lines);
		@split_linea = split (",",$lines);
		#$pos = index($lines,",");
		$linea_baja= index ($split_linea[1],"_");
		#$pos_punto = index($lines,".");
		$pos_punto2 = index ($split_linea[1],".");
		$posicion = substr($split_linea[1],$linea_baja+1,$pos_punto2-$linea_baja-1);
		$posicion_negativo = substr($split_linea[1],$linea_baja+2,$pos_punto2-$linea_baja-1);
		#$linea_baja= index ($split_linea[1],"_");
		#$numero_posiciones= $linea_baja-$pos;
		$sufijo = substr($split_linea[1],0,$linea_baja);
		$orientacion = substr ($split_linea[1],$linea_baja+1,1);
		#$orientacion = substr($lines,$pos+3,1);
            	#print "secuencia ",$split_linea[0]," posicion \t",$posicion,"\n";
		if ($sufijo eq "1")          
            	{
              		#print CHR1 $split_linea[0],"_",$orientacion,$posicion,"\t";
			if ($orientacion eq "-")
			{
				print CHR1 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
				print CHR1 (($posicion_negativo+2)-$splitpunto[2]),"\t";
				print CHR1 $posicion_negativo+1,"\t";
				print CHR1 "-","\n";
			}
			else
			{
				print CHR1 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
				print CHR1 $posicion+1,"\t";
				print CHR1 $posicion+$splitpunto[2],"\t";
				print CHR1 "+","\n";
			}
            	}
            	elsif ($sufijo eq "2")
            	{
			if ($orientacion eq "-")
                        {
                                print CHR2 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
				print CHR2 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR2 $posicion_negativo+1,"\t";
                                print CHR2 "-","\n";
                        }
                        else
                        {
                                print CHR2 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
				print CHR2 $posicion+1,"\t";
                                print CHR2 $posicion+$splitpunto[2],"\t";
                                print CHR2 "+","\n";
			}
            	}
            	elsif ($sufijo eq "3")
            	{
			if ($orientacion eq "-")
                        {
                                print CHR3 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
				print CHR3 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR3 $posicion_negativo+1,"\t";
                                print CHR3 "-","\n";
                        }
                        else
                        {
                                print CHR3 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
				print CHR3 $posicion+1,"\t";
                                print CHR3 $posicion+$splitpunto[2],"\t";
                                print CHR3 "+","\n";
                        }

            	}
	        elsif ($sufijo eq "4")
            	{
			if ($orientacion eq "-")
                        {
                                print CHR4 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
				print CHR4 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR4 $posicion_negativo+1,"\t";
                                print CHR4 "-","\n";
                        }
                        else
                        {
                                print CHR4 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
				print CHR4 $posicion+1,"\t";
                                print CHR4 $posicion+$splitpunto[2],"\t";
                                print CHR4 "+","\n";
                        }

            	}
		elsif ($sufijo eq "5")
                {
			if ($orientacion eq "-")
                        {
                                print CHR5 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
				print CHR5 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR5 $posicion_negativo+1,"\t";
                                print CHR5 "-","\n";
                        }
                        else
                        {
                                print CHR5 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
				print CHR5 $posicion+1,"\t";
                                print CHR5 $posicion+$splitpunto[2],"\t";
                                print CHR5 "+","\n";
                        }

                }
		elsif ($sufijo eq "6")
                {
			if ($orientacion eq "-")
                        {
                                print CHR6 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
				print CHR6 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR6 $posicion_negativo+1,"\t";
                                print CHR6 "-","\n";
                        }
                        else
                        {
                                print CHR6 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
				print CHR6 $posicion+1,"\t";
                                print CHR6 $posicion+$splitpunto[2],"\t";
                                print CHR6 "+","\n";
                        }

                }
		elsif ($sufijo eq "7")
                {
                        if ($orientacion eq "-")
                        {
                                print CHR7 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
                                print CHR7 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR7 $posicion_negativo+1,"\t";
                                print CHR7 "-","\n";
                        }
                        else
                        {
                                print CHR7 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
                                print CHR7 $posicion+1,"\t";
                                print CHR7 $posicion+$splitpunto[2],"\t";
                                print CHR7 "+","\n";
                        }

                }
		elsif ($sufijo eq "8")
                {
                        if ($orientacion eq "-")
                        {
                                print CHR8 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
                                print CHR8 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR8 $posicion_negativo+1,"\t";
                                print CHR8 "-","\n";
                        }
                        else
                        {
                                print CHR8 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
                                print CHR8 $posicion+1,"\t";
                                print CHR8 $posicion+$splitpunto[2],"\t";
                                print CHR8 "+","\n";
                        }

                }
		elsif ($sufijo eq "9")
                {
                        if ($orientacion eq "-")
                        {
                                print CHR9 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
                                print CHR9 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR9 $posicion_negativo+1,"\t";
                                print CHR9 "-","\n";
                        }
                        else
                        {
                                print CHR9 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
                                print CHR9 $posicion+1,"\t";
                                print CHR9 $posicion+$splitpunto[2],"\t";
                                print CHR9 "+","\n";
                        }

                }
		elsif ($sufijo eq "10")
                {
                        if ($orientacion eq "-")
                        {
                                print CHR10 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
                                print CHR10 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR10 $posicion_negativo+1,"\t";
                                print CHR10 "-","\n";
                        }
                        else
                        {
                                print CHR10 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
                                print CHR10 $posicion+1,"\t";
                                print CHR10 $posicion+$splitpunto[2],"\t";
                                print CHR10 "+","\n";
                        }

                }
		elsif ($sufijo eq "11")
                {
                        if ($orientacion eq "-")
                        {
                                print CHR11 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
                                print CHR11 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR11 $posicion_negativo+1,"\t";
                                print CHR11 "-","\n";
                        }
                        else
                        {
                                print CHR11 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
                                print CHR11 $posicion+1,"\t";
                                print CHR11 $posicion+$splitpunto[2],"\t";
                                print CHR11 "+","\n";
                        }

                }
		elsif ($sufijo eq "12")
                {
                        if ($orientacion eq "-")
                        {
                                print CHR12 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
                                print CHR12 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR12 $posicion_negativo+1,"\t";
                                print CHR12 "-","\n";
                        }
                        else
                        {
                                print CHR12 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
                                print CHR12 $posicion+1,"\t";
                                print CHR12 $posicion+$splitpunto[2],"\t";
                                print CHR12 "+","\n";
                        }

                }
		elsif ($sufijo eq "13")
                {
                        if ($orientacion eq "-")
                        {
                                print CHR13 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
                                print CHR13 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR13 $posicion_negativo+1,"\t";
                                print CHR13 "-","\n";
                        }
                        else
                        {
                                print CHR13 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
                                print CHR13 $posicion+1,"\t";
                                print CHR13 $posicion+$splitpunto[2],"\t";
                                print CHR13 "+","\n";
                        }

                }
		elsif ($sufijo eq "14")
                {
                        if ($orientacion eq "-")
                        {
                                print CHR14 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
                                print CHR14 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR14 $posicion_negativo+1,"\t";
                                print CHR14 "-","\n";
                        }
                        else
                        {
                                print CHR14 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
                                print CHR14 $posicion+1,"\t";
                                print CHR14 $posicion+$splitpunto[2],"\t";
                                print CHR14 "+","\n";
                        }

                }
		elsif ($sufijo eq "15")
                {
                        if ($orientacion eq "-")
                        {
                                print CHR15 $split_linea[0],"_",$orientacion,(($posicion_negativo+2)-$splitpunto[2]),"_",$posicion_negativo+1,"\t";
                                print CHR15 (($posicion_negativo+2)-$splitpunto[2]),"\t";
                                print CHR15 $posicion_negativo+1,"\t";
                                print CHR15 "-","\n";
                        }
                        else
                        {
                                print CHR15 $split_linea[0],"_",$posicion+1,"_",$posicion+$splitpunto[2],"\t";
                                print CHR15 $posicion+1,"\t";
                                print CHR15 $posicion+$splitpunto[2],"\t";
                                print CHR15 "+","\n";
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
close (CHR7);
close (CHR8);
close (CHR9);
close (CHR10);
close (CHR11);
close (CHR12);
close (CHR13);
close (CHR14);
close (CHR15);
