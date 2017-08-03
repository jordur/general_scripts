#!usr/bin/perl -w
use strict;
######### script per navegar per diferents directoris i extreure informació de fitxers amb una determinada extensió #######
open OUT, ">num_seq";
my @dirs=qw /\/media\/disk\/Feina\/prova\/rm\/class_script\/final\/toxines_senseNA\/mids\/mid1\/reads\/blast 
			\/media\/disk\/Feina\/prova\/rm\/class_script\/final\/toxines_senseNA\/mids\/mid2\/reads\/blast
			\/media\/disk\/Feina\/prova\/rm\/class_script\/final\/toxines_senseNA\/mids\/mid3\/reads\/blast
			\/media\/disk\/Feina\/prova\/rm\/class_script\/final\/toxines_senseNA\/mids\/mid4\/reads\/blast
			\/media\/disk\/Feina\/prova\/rm\/class_script\/final\/toxines_senseNA\/mids\/mid5\/reads\/blast
			\/media\/disk\/Feina\/prova\/rm\/class_script\/final\/toxines_senseNA\/mids\/mid6\/reads\/blast
			\/media\/disk\/Feina\/prova\/rm\/class_script\/final\/toxines_senseNA\/mids\/mid7\/reads\/blast
			\/media\/disk\/Feina\/prova\/rm\/class_script\/final\/toxines_senseNA\/mids\/mid8\/reads\/blast/;
foreach my $direc (@dirs){
	chdir $direc;
#~ print "@dirs\n";
#~ foreach my $direc (@dirs){
#~ opendir(DIR, $dir1) or die "$!";;
my @files=glob "*_cobalt";
chomp(@files);
#~ print "@files\n";
foreach my $a (@files){
		open DAT, $a  or die "Can't open the file!";
		my @data=<DAT>;
		#~ print "@data\n";
		my $count=0;
		foreach my $line (@data){
			if ($line =~/^>[0-9]./){
				#~ print $line,"\n";
				$count++;
				
			}
			
		}
		print OUT $a,"\t",$count,"\n";
	}
print OUT "\n";
}
