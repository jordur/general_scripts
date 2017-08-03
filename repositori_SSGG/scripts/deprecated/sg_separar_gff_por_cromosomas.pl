#!/usr/bin/perl

open (FILE1, "<", @ARGV[0]);#Fichero de entrada como parametro. I/O en ruta local usuario.
open (CHR1, ">", 'chr1.gff') or die "No puedo abrir chr1.gff.\n";
open (CHR2, ">", 'chr2.gff') or die "No puedo abrir chr2.gff.\n";
open (CHR3, ">", 'chr3.gff') or die "No puedo abrir chr3.gff.\n";
open (CHR4, ">", 'chr4.gff') or die "No puedo abrir chr4.gff.\n";
open (CHR5, ">", 'chr5.gff') or die "No puedo abrir chr5.gff.\n";
open (CHR6, ">", 'chr6.gff') or die "No puedo abrir chr6.gff.\n";

#my $almohadilla = substr($lines,0,1);
for (my $i=1; my $lines=<FILE1>;$i++)
{	
	
	$lines2 = $lines;
	$lines2 =~s/.*i=/i=/;
	$lines2 =~s/;.*$/;/;
	#print $lines2,"\n";	
#$lines =~s /;.*$//;

	$poschr = index($lines2,"=");
	$chr= substr ($lines2,$poschr+1,1);
	#print $chr,"\n";		
	if ($chr == "1")
	{
		print CHR1 $lines,"\n";

	}
	elsif ($chr == "2")
	{
		print CHR2 $lines,"\n";
	}
	elsif ($chr == "3")
	 {
                print CHR3 $lines,"\n";
        }
	 elsif ($chr == "4")
        {
                print CHR4 $lines,"\n";
        }
	 elsif ($chr == "5")
        {
                print CHR5 $lines,"\n";
        }
	 elsif ($chr == "6")
        {
                print CHR6 $lines,"\n";
        }


}

close (FILE1);
close (CHR1);
close (CHR2);
close (CHR3);
close (CHR4);
close (CHR5);
close (CHR6);
