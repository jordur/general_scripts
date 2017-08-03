use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Getopt::Std;
use Scalar::Util qw(looks_like_number);
use Cwd;
use FileHandle;

my @chr_list = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X);

#Create one biobase file for each chr
foreach my $chr (@chr_list){
	&convert_biobase_bed($chr);
}
#



sub get_bio_file_handle {
	my $chr = shift;
	my $bio_file_handle = new FileHandle;
	my $file = "/home/likewise-open/SGNET/gmarco/biobase/bioBaseSortUniques_chr$chr.txt";
	$bio_file_handle->open($file) or die("ERROR: Could not read from input file: $file\n");
	return $bio_file_handle;
}


sub convert_biobase_bed () {
	
    my ($chr) = $_[0];
    my $bio_file_handle = &get_bio_file_handle($chr);
    my $desc = "-\t-";
    open (OUTPUT, ">/home/likewise-open/SGNET/gmarco/sg_hgnc_no_cover_testing/bioBed/bioBed_chr$chr.bed") or die "aaaa";
    while (my $line = <$bio_file_handle>) {
        chomp $line;
        
        if ($line =~ /^Type/) {next;}
        
        my @data = split (/\t/, $line);
        
        my ($chr, $start, $end) = split (/[:-]/, $data[4]);
        
        if (!$end or $end eq '' or $end eq '+') { $end = $start + 1;}
        
        $desc = $chr."\t".$start."\t".$end."\t".$data[1]."|".$data[2]."|".$data[6]."|".$data[9]."|".$data[7]."|";
   
		print OUTPUT $desc,"\n";
    }
    close OUTPUT;
    
   	#Sort bioBed file
   	system("sortBed -i /home/likewise-open/SGNET/gmarco/sg_hgnc_no_cover_testing/bioBed/bioBed_chr$chr.bed > /home/likewise-open/SGNET/gmarco/sg_hgnc_no_cover_testing/bioBed/sorted/bioBed_sorted_chr$chr.bed");
   	
    return 1;
}