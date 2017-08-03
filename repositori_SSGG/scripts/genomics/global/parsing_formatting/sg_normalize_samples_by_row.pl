#!/usr/bin/perl -w
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script to normalize a matrix
# Date: May 2013
# @Author: Sheila
# @Version: 0.1
########################################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Cwd;
use Math::NumberCruncher;

# -------------------------
# -- Variable definition --
# -------------------------

my @split_lines;
my $sumRow=0;
my @normByrow;
my @split_mids;
my @midCols;
my @sampleNames;
my @normRow;
my @medianbyrow;
my $median;
# ---------------------------
# -- functions definitions --
# ---------------------------

sub normalize_by_row
{
    my ($file,$outfile,$samples) = @_;
	open (MATRIX,"<","$file");
	open (BYROW,">","$outfile"."_normalized_by_row.tab");
	while (my $lines=<MATRIX>)
    {
        chomp($lines);
		@split_lines=split("\t",$lines);
		@split_mids=split(",",$samples);
		if ($lines =~ /AMPLICON/)
        {
			for (my $a=1;$a<=$#split_lines;$a++)
			{
				for (my $e=0;$e<=$#split_mids;$e++)
            	{
					if ($split_lines[$a] eq "MID".$split_mids[$e])
					{
						push(@midCols,$a);
						push(@sampleNames,"MID".$split_mids[$e]);
					}
				}
			}
		}
#		print BYROW "AMPLICON\t",join("\t",@sampleNames),"\n";
		if ($lines !~ /AMPLICON/)
		{
			for (my $i=0;$i<=$#midCols;$i++)
			{
				#$sumRow=$sumRow+$split_lines[$midCols[$i]];
				push (@medianbyrow,$split_lines[$midCols[$i]]);
			}
			$median = Math::NumberCruncher::Median(\@medianbyrow,1);
			for (my $o=0;$o<=$#midCols;$o++)
    	    {
        	    #push (@normByrow,($split_lines[$midCols[$o]]/($sumRow/($#midCols+1))));
				if (($split_lines[$midCols[$o]] == 0) || ($median ==0))
				{
					push (@normByrow,"0");
				}
				else
				{
					push  (@normByrow,($split_lines[$midCols[$o]]/$median));
				}
    	    }
			push (@normRow,$split_lines[0]."\t".join ("\t",@normByrow));
			#print BYROW join ("\t",@normRow),"\n";
			@normByrow=();
			@medianbyrow=();
			#@normRow=();
			#$sumRow=0;
		}
	}
	print BYROW "AMPLICON\t",join("\t",@sampleNames),"\n";
	print BYROW join ("\n",@normRow),"\n";
	close(BYROW);
	close(MATRIX);
}


# ---------------------------------
# -- global & default parameters --
# ---------------------------------
my $input;
my $output;
my $mids;
# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"i=s"	=> \$input,
				"o=s"	=> \$output,
				"mid=s"   => \$mids,
);

print "\n=================================================================
$basename: is a script to normalize values giving a matrix. It's intendeed to be used for BRCA samples process with amplicon kits such as Multiplicom's\n";

if (not defined $input) {
	die "
			Options:
		-i Input file. A matrix with the following format:\n
			AMPLICON    MID1    MID2    MID3    MID4    MID5    MID6    MID7    MID8\n
			BRCA1_ex02_01   115 0   108 0   89  100 135 113\n
			NOTE: This is a simple example, the number of columns can vary as well as the name of the amplicon.
		-o Ouput path
		-mid A coma-separated list of MID (Multiplex Identifier) numbers to be analyzed. Ex: MID1,MID2,MID4
\n".localtime()."\n=================================================================\n\n";
}

print  "\nProcess started at... ".localtime()."\n";
&main ();
print  "\nProcess finished... ".localtime()."\n";

# ---------------
# -- functions --
# ---------------
sub main (){

	print "Running $basename with parameters f=$input, o=$output\n";
	&normalize_by_row($input,$output,$mids)

}

exit;
