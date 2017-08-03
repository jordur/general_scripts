#!/usr/bin/perl
use List::Util qw (sum);
use Smart::Comments;
use strict;
use warnings;
use diagnostics;

#Inputs:
#Fichero de exones tipo:
#chr16   refGene2gff     exon    1	40	.       +       .       ID=NLRC5;refSeq=NM_032206
#chr20   refGene2gff     exon    2000	3630	.       -       .       ID=KIAA1755;refSeq=NM_001029864
#Fichero de cromosomas tipo:
#>chr1
#TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC
#CCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCT
#ACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCTAACC
#>chr2
#TAACCCTAACCCTAACCCTAACCCTAACCCqrACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAT
#...
#Outputs:
#Un fichero por cromosoma con las parejas _b y _a que indican los fragmentos a concatenar y sus secuencias concatenadas tipo:
#>chr10_COL13A1_3b_5a
#GGTCATTTTCTTGCCATTGAGTTGTTTGATTTCCTTATATTTATTTCTTCAAAGGTAAATATATTTACTTAAAAATATGC        -
#>chr10_COL13A1_4b_7a
#GGTCATTTTCTTGCCATTGAGTTGTTTGATTTCCTTATATTCTCATTTTCATCTAAACCTTTTGATTCTTGAATCACTCT        -
#Para crear las combinaciones solo se tienen en cuenta las del mismo gen y solo genes con mas de 1 exon.
#Se considera tambien si es la cadena + (directa) o - (inversa) por separado. Para estas ultimas secuencias hay que realizar el reverso complementario.



sub fragmento
{
  my ($indice, $pos_rel, $posicion2, $j, $v_ref) = @_; 
  my $superfragmento;
  my $fragmento;
  my $posicion3;
  my $desplazamiento;
  
  $posicion3 = $posicion2;
  $desplazamiento = $j;
  $superfragmento = @{$v_ref}[$desplazamiento];

  while ( ($posicion3 - $indice < 40) )
  {
    $desplazamiento ++;
    if (defined(@{$v_ref}[$desplazamiento]))
    {
      $superfragmento .= @{$v_ref}[$desplazamiento];
      $posicion3 += length(@{$v_ref}[$desplazamiento]);
      #print " DEFINIDO DES $desplazamiento POS3 $posicion3 POS $indice \n";
    }else#CREO QUE AQUI YA NO SE ENTRA
    {
      #print "\nERROR: Posicion dentro de secuencia cuyo rango cae DESPUES, situacion imposible implica error conceptual.\n";
      #print " NO DEFINIDO DES $desplazamiento POS3 $posicion3 POS $indice \n";
      return "ERROR: RANGO CAE DESPUES.";
      last;
    }
  }

  #print " Superfragmento $superfragmento\n";
  #print " Indice ".($indice-1)."";
  $fragmento = substr($superfragmento, (($indice-1) - $pos_rel), 40);
  #print " Fragmento $fragmento\n";
  if (length($fragmento) < 40)
  {
    #print "\nERROR: Posicion dentro de secuencia cuyo rango cae ANTES, situacion imposible implica error conceptual.\n";
    return "ERROR: RANGO CAE ANTES.";
    last; 
  }
  return $fragmento;
}


                   
sub anade_fragmentos_iniciales
{
  my ($fichero1_ref, $fichero2_ref, $fichero1_fragmentos_ref) = @_;
  my $k;
  my $v;
  my $posicion2;
  my @fichero1_fragmentos;
  my $dimension;
  my $escapa;
  my $pos_rel;

  if (%{$fichero2_ref}) 
  {
    while (($k,$v) = each %{$fichero2_ref})
    {### |...   | [%]
      #print "\n---> k=$k v=$v ".@{$v}[0]." - ".@{$v}[1]."\n\n";
      #print "\n Cromosoma $k\n";
      if ( exists(${$fichero1_ref}{$k}) )
      {
        @{${$fichero1_ref}{$k}} = sort { $a->[0] <=> $b->[0] } @{${$fichero1_ref}{$k}};
        $posicion2 = 0;
      }else
      {
        next;#cuando ni siquiera hay array ignora el indice
      }
      #print "UNO - ".@{${$fichero1_ref}{$k}}[0]." DOS - ".@{$v}[0]."\n";
    
      $escapa = 0;
      for (my $j=0; $j<=$#{$v}; $j++)
      {
        $pos_rel = $posicion2;
        $posicion2 += length(@{$v}[$j]);           
        #print "\nPOSICION DOS ".$posicion2." - ".@{$v}[$j]."\n";
            
        $dimension = $#{${$fichero1_ref}{$k}};
        for (my $l=0; $l<=$dimension; $l++)          
        {
          #print " POSICION UNO ".${@{${$fichero1_ref}{$k}}}[0][0]." - ".${@{${$fichero1_ref}{$k}}}[0][1]."\n";
          #print " l=$l dim=$dimension -> ".(${@{${$fichero1_ref}{$k}}}[$l][0])." con ".$posicion2;
          
          if ( (${@{${$fichero1_ref}{$k}}}[$l][0]) <= $posicion2 )
          {
            #print " MENOR ";
            #print fragmento(${@{${$fichero1_ref}{$k}}}[$l][0],$pos_rel,$posicion2,$j,\@{$v});
            
            push @{${$fichero1_ref}{$k}[0]}, fragmento(${@{${$fichero1_ref}{$k}}}[$l][0],$pos_rel,$posicion2,$j,\@{$v});
            push @{${$fichero1_fragmentos_ref}{$k}}, ${$fichero1_ref}{$k}[0];
            #use Data::Dumper;
            #print Dumper \@{${$fichero1_ref}{$k}}; print "\n+++++\n";
            #print Dumper \@{${$fichero1_fragmentos_ref}{$k}}; print "\n+++++++\n";
            
            if ($dimension == 0)
            {
              #print " ENTRA1 \n";
              delete(${$fichero1_ref}{$k});
              $escapa = 1;
              last;
            }else
            {
              #print " ENTRA2 \n";
              shift(@{${$fichero1_ref}{$k}});
              $dimension--;
              redo;
            }
          }else
          {
            last;#cuando encuentra un valor del array mayor que la posicion pasa a la siguiente linea de DOS
          }

        }#end for 1

        if ($escapa == 1)
        {
          last;#cuando no hay mas valores en el array escapa
        }      
      }#end for 2
    }#end while 2
  }
}#end sub anade_fragmentos_iniciales


                                                                       
sub anade_fragmentos_finales
{
  my ($fichero1_ref, $fichero2_ref, $fichero1_fragmentos_ref) = @_;
  my $k;
  my $v;
  my $posicion2;
  my $dimension;
  my $escapa;
  my $pos_rel;

  if (%{$fichero2_ref}) 
  {
    while (($k,$v) = each %{$fichero2_ref})
    {### |...   | [%]
      #print "\n---> k=$k v=$v ".@{$v}[0]." - ".@{$v}[1]."\n\n";
      #print "\n Cromosoma $k\n";
      if ( exists(${$fichero1_ref}{$k}) )
      {
        @{${$fichero1_ref}{$k}} = sort { $a->[1] <=> $b->[1] } @{${$fichero1_ref}{$k}};
        $posicion2 = 0;
      }else
      {
        next;#cuando ni siquiera hay array ignora el indice
      }
      #print "UNO - ".@{${$fichero1_ref}{$k}}[0]." DOS - ".@{$v}[0]."\n";
    
      $escapa = 0;
      for (my $j=0; $j<=$#{$v}; $j++)
      {
        $pos_rel = $posicion2;
        $posicion2 += length(@{$v}[$j]);           
        #print "\nPOSICION DOS ".$posicion2." - ".@{$v}[$j]."\n";
            
        $dimension = $#{${$fichero1_ref}{$k}};
        for (my $l=0; $l<=$dimension; $l++)          
        {
          #print " POSICION UNO ".${@{${$fichero1_ref}{$k}}}[0][0]." - ".${@{${$fichero1_ref}{$k}}}[0][1]."\n";
          #print " l=$l dim=$dimension -> ".((${@{${$fichero1_ref}{$k}}}[$l][1])-39)." con ".$posicion2;
          
          if ( (${@{${$fichero1_ref}{$k}}}[$l][1]-39) <= $posicion2 )
          {
            #print " MENOR ";
            #print fragmento((${@{${$fichero1_ref}{$k}}}[$l][1]-39),$pos_rel,$posicion2,$j,\@{$v});

            #print ${${$fichero1_fragmentos_ref}{$k}}[1]." IGUAL ". ${$fichero1_ref}{$k}[0][1];
            foreach (@{${$fichero1_fragmentos_ref}{$k}})
            {
              #print ${$_}[0]." ";
              if ( (${$_}[0] == ${$fichero1_ref}{$k}[0][0]) && (${$_}[1] == ${$fichero1_ref}{$k}[0][1]) && (${$_}[2] eq ${$fichero1_ref}{$k}[0][2]) && (${$_}[3] eq ${$fichero1_ref}{$k}[0][3]) )
              {
                #
                #print ${$_}[1]." IGUAL ".${$fichero1_ref}{$k}[0][1]."\n";
                push @{$_}, fragmento((${@{${$fichero1_ref}{$k}}}[$l][1]-39),$pos_rel,$posicion2,$j,\@{$v});
                last;
              }
            }
            #use Data::Dumper;
            #print Dumper \@{${$fichero1_ref}{$k}}; print "\n+++++\n";
            #print Dumper \@{${$fichero1_fragmentos_ref}{$k}}; print "\n+++++++\n";
              
            if ($dimension == 0)
            {
              #print " ENTRA1 \n";
              delete(${$fichero1_ref}{$k});
              $escapa = 1;
              last;
            }else
            {
              #print " ENTRA2 \n";
              shift(@{${$fichero1_ref}{$k}});
              $dimension--;
              redo;
            }
          }else
          {
            last;#cuando encuentra un valor del array mayor que la posicion pasa a la siguiente linea de DOS
          }

        }#end for 1

        if ($escapa == 1)
        {
          last;#cuando no hay mas valores en el array escapa
        }      
      }#end for 2
    }#end while 2
  }
}#end sub anada_fragmentos_finales


                   
sub reverso_complementario
{
  my ($sentence_ref) = @_;
  my $sentence = ${$sentence_ref};  
  my $reverso = "";
  
  $sentence =~ tr/ATGC/TACG/;

  for (my $i=0; $i<40; $i++)
  {
    $reverso .= chop($sentence);
  }

  return $reverso;
}



sub proceso
{
  my ($chr, $exonA_pos_ref, $exonB_pos_ref, $exonA_neg_ref, $exonB_neg_ref) = @_;
  my $kB;
  my $vB;
  my @valorB;
  my @valorA;

  #print "\nCromosoma $chr ...\n";
  open (OUT_CHR, ">", 'out_'.$chr.'.txt') or die "No puedo abrir el fichero.\n";

  ordenar_array(\%{$exonA_pos_ref});
  #ordenar_array(\%{$exonB_pos_ref});
        
  if ( (%{$exonA_pos_ref}) && (%{$exonB_pos_ref}) )
  {                                                                                                                                               
    while (($kB,$vB) = each %{$exonB_pos_ref})
    { 
      @valorB = @{$vB};
      #print "\n---> k=$kB v=@valorB\n";
      #use Data::Dumper;
      #print Dumper \@valorB; print "\n";
      
      if ( exists(${$exonA_pos_ref}{$kB}) )
      {
        @valorA = @{${$exonA_pos_ref}{$kB}};       
        #print Dumper \@valorA; print "\n";
        #print "\nTamB=".scalar(@valorB)." TamA=".$#valorA."\n";
        
        for (my $i=0; $i<=$#valorB; $i++)
        {
          #print "\nvalorB=".$valorB[$i][0]." ".$valorB[$i][1]."\n";
                      
          for (my $j=0; $j<=$#valorA; $j++)
          {
            #print "   valorA=".$valorA[$j][0]." ".$valorA[$j][1]."\n";
            if ($valorB[$i][0] < $valorA[$j][0])
            {
              print OUT_CHR ">".$chr."_".$kB."_".($i+1)."b_".($j+1)."a\n";
              print OUT_CHR $valorB[$i][1]."".$valorA[$j][1]."\n";
            }
          }                 
        }
      }       
    } 
  }

  #ordenar_array(\%{$exonA_neg_ref});
  #ordenar_array(\%{$exonB_neg_ref});

  if ( (%{$exonA_neg_ref}) && (%{$exonB_neg_ref}) )
  {                                                                                                                                               
    while (($kB,$vB) = each %{$exonB_neg_ref})
    { 
      @valorB = @{$vB};
      #print "\n---> k=$kB v=@valorB\n";
      #use Data::Dumper;
      #print Dumper \@valorB; print "\n";
      
      if ( exists(${$exonA_neg_ref}{$kB}) )
      {
        @valorA = @{${$exonA_neg_ref}{$kB}};       
        #print Dumper \@valorA; print "\n";
        #print "\nTamB=".scalar(@valorB)." TamA=".$#valorA."\n";
        
        for (my $i=0; $i<=$#valorB; $i++)
        {
          #print "valorB=".$valorB[$i][0]." ".$valorB[$i][1]."\n";
          for (my $j=0; $j<=$#valorA; $j++)
          {
            #print "	valorA=".$valorA[$j][0]." ".$valorA[$j][1]."\n";            
            if ($valorB[$i][0] < $valorA[$j][0])
            {
              print OUT_CHR ">".$chr."_".$kB."_".($j+1)."a_".($i+1)."b\n";
              #print $valorA[$j][1]." ".$valorB[$i][1]."\n";
              print OUT_CHR reverso_complementario(\$valorA[$j][1])."".reverso_complementario(\$valorB[$i][1])."\n";
            }
          }          
        }
      }       
    } 
  }

  close (OUT_CHR);
}# end sub proceso  



sub ordenar_array
{
  my ($exonA_pos_ref) = @_;

  if (%{$exonA_pos_ref})
  {
    #print "El hash tiene valores \n";
    my $k;
    #print "tamano hash=".(keys %{$exonA_pos_ref})."\n";
    foreach $k (keys %{$exonA_pos_ref})
    {
      #print "key=".$k." valor hash=".${$exonA_pos_ref}{$k}." tamano array=".($#{${$exonA_pos_ref}{$k}}+1)."\n";

      @{${$exonA_pos_ref}{$k}} = sort { $a->[0] <=> $b->[0] } @{${$exonA_pos_ref}{$k}};

      #for (my $l=0; $l<=$#{${$exonA_pos_ref}{$k}}; $l++)
      #{
      #  print "key=$k l=($l+1) pos=@{${${$exonA_pos_ref}{$k}}[$l]}[0] sec=@{${${$exonA_pos_ref}{$k}}[$l]}[1] \n";
      #}
    }
  }                                                                              
  #print "\n";
}


                                                                       
#PROCESAMOS



print "\nExtrayendo la informacion seleccionada del fichero <".$ARGV[0]."> ...\n";

open (FILE1, "<", $ARGV[0]);#Fichero de entrada como parametro. I/O en ruta local usuario.

my @fila;
my %fichero1a;
my %fichero1b;
my $posicion;
my $valor;

for (my $i=1; my $lines=<FILE1>; $i++)
{
  @fila = split (/\s+/,$lines);
  $posicion = index(substr($fila[8],2),";");
  if ($posicion != -1)
  {
    $valor = substr(substr($fila[8],2),1,$posicion-1);
  }
  push @{$fichero1a{$fila[0]}}, [ int($fila[3]), int($fila[4]), $valor, $fila[6] ];    
  push @{$fichero1b{$fila[0]}}, [ int($fila[3]), int($fila[4]), $valor, $fila[6] ];
}

close (FILE1);#Cierro la entrada del fichero
                


print "\nExtrayendo la informacion seleccionada del fichero <".$ARGV[1]."> ...\n\n";

open (FILE2, "<", $ARGV[1]);#Fichero de entrada como parametro. I/O en ruta local usuario.

my %fichero2;
my $chr;

for (my $i=1; my $lines=<FILE2>; $i++)
{
  chomp($lines);
  if (substr($lines,0,1) eq ">")
  {
    $chr = substr($lines,1);
    print "  Cromosoma $chr\n";      
  }else
  {
    push @{$fichero2{$chr}}, $lines;
  }
}

close (FILE2);#Cierro la entrada del fichero



print "\nCruzando ficheros para componer fragmentos iniciales ...\n\n";

my %fichero1_fragmentos;

anade_fragmentos_iniciales(\%fichero1a,\%fichero2,\%fichero1_fragmentos);



print "\n\nCruzando ficheros para componer fragmentos finales ...\n\n";

anade_fragmentos_finales(\%fichero1b,\%fichero2,\%fichero1_fragmentos);
#use Data::Dumper;
#print Dumper \%fichero1_fragmentos;



print "\n\nCargando vectores y emparejando exones ...\n\n";

my %hashA_pos;
my %hashB_pos;
my %hashA_neg;
my %hashB_neg;
my $k;
my $v;

if (%fichero1_fragmentos)
{
  while (($k,$v) = each %fichero1_fragmentos)
  {### |...   | [%]
    #print "\n\n---> k=$k\n\n";
    if ( exists($fichero1_fragmentos{$k}) )
    {
      for (my $j=0; $j<=$#{$v}; $j++)
      {
        #print "\n $k $j - ".${@{$v}[$j]}[0]." ".${@{$v}[$j]}[1]." ".${@{$v}[$j]}[2]." ".${@{$v}[$j]}[3]." ".${@{$v}[$j]}[4]." ".${@{$v}[$j]}[5]."\n";
      
        if (${@{$v}[$j]}[3] eq "+")
        {
          push @{$hashA_pos{${@{$v}[$j]}[2]}}, [ ${@{$v}[$j]}[0], ${@{$v}[$j]}[4] ];
          push @{$hashB_pos{${@{$v}[$j]}[2]}}, [ ${@{$v}[$j]}[1], ${@{$v}[$j]}[5] ];
        }else
        {
          push @{$hashA_neg{${@{$v}[$j]}[2]}}, [ ${@{$v}[$j]}[0], ${@{$v}[$j]}[4] ];
          push @{$hashB_neg{${@{$v}[$j]}[2]}}, [ ${@{$v}[$j]}[1], ${@{$v}[$j]}[5] ];
        }
      }

      proceso($k,\%hashA_pos,\%hashB_pos,\%hashA_neg,\%hashB_neg);

      %hashA_pos = ();
      %hashB_pos = ();
      %hashA_neg = ();
      %hashB_neg = ();

    }
  }
}


                                                                                                                                                                                        
print "\n\nOperacion finalizada.\n\n";

