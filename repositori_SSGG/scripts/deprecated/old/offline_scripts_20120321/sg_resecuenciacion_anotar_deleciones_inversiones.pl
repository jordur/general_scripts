#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Anota las deñecopmes encontradas \n COMO SE USA:\n";
    print " Input: 1) fichero tabulado con 3 columnas, rango_inicio, rango_final y tamano aproximado de la delecion \n 71925-71927     73752-73754     1827 \n";
    print "        2) fichero gff con la anotacion\n";
    print " Output: fichero de las deleciones con informacion anadida sobre el gen en el que mapean  \n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(DELECIONES,"<",@ARGV[0]);
open(GFF,"<",@ARGV[1]);
my @deleciones;
my @gff;
for (my $i=1; my $lines=<DELECIONES>; $i++)
{
  chomp($lines);
  push (@deleciones,$lines);
}

close (DELECIONES);

#print join ("\n",@deleciones),"\n";
for (my $i=1; my $lines2=<GFF>; $i++)
{
  chomp($lines2);
  push (@gff,$lines2);
}

close (GFF);
#print join ("\n",@gff),"\n";
print "Start","\t","End","\t","Length","\t","CDS start","\t","CDS end","\t","Gene name","\t","Product","\t","Accession number","\n";
my $contador =0;
for (my $i=0; $i<=$#deleciones; $i++)
{
        @split_deleciones = split (/\t/,$deleciones[$i]);
	@split_inicio = split ("-",$split_deleciones[0]);
	@split_final = split ("-",$split_deleciones[1]);
	
        for (my $a=0; $a<=$#gff; $a++)
        {

                @split_gff = split (/\t/,$gff[$a]);
                if ((($split_deleciones[0] >= $split_gff[3]) && ($split_deleciones[0] <= $split_gff[4]))|| (($split_deleciones[1] >= $split_gff[3]) && ($split_deleciones[1] <= $split_gff[4])) || (($split_deleciones[0] >= $split_gff[3]) && ($split_deleciones[1] <= $split_gff[4]))|| (($split_deleciones[0] <= $split_gff[3]) && ($split_deleciones[1] >= $split_gff[4])))
                {
                        
			#NC_004431       GenBank CDS     190     255     .       +       .       "ID=thrL.p01;Parent=thrL.t01;Dbxref=GI:26245918,GeneID:1040154;Note=Thr operon attenuator%3B Escherichia coli K-12 ortholog: b0001%3B Escherichia coli O157:H7 ortholog: z0001;codon_start=1;function=leader%3B Amino acid biosynthesis: Threonine;gene=thrL;locus_tag=c5491;product=thr operon leader peptide;protein_id=NP_751957.1;transl_table=11"
			$split_gff[8]=~s/Parent=.*product=//;
                        #$split_gff[8]=~s/;transl_table=.*$//;
                        @Gene_ID=split (";",$split_gff[8]);
                        $Gene_ID[0]=~s/"ID=//; 
			$Gene_ID[2]=~s/protein_id=//;
			#print  $Gene_ID[2],"\n";
			print $deleciones[$i],";",$split_gff[3],";",$split_gff[4],";",$Gene_ID[0],";",$Gene_ID[1],";",$Gene_ID[2],"\n";
			$contador = $contador+1;
                }
        }
	
	if ($contador == 0)
	{
		print $deleciones[$i],"\n";
	}
	$contador=0;
}

