#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genómicos S.L.
# @Desc: Script for generating the BioScope folders structures and launching the analyses
# @Author: Arbol
# @Contributors: Arbol
########################################################

# -----------------------
# -- Include libraries --
# -----------------------
use File::Basename;
use Cwd;
use Cwd qw( abs_path );
use strict;
use warnings;
use POSIX qw( strftime );
use XML::Simple;
use Data::Dumper;

# ---------------
# -- functions --
# ---------------

sub usage {
# Displays the script usage
	# Parameters:
	my ( $scriptname ) = @_;

	print "
$scriptname is a script for generating the BioScope analysis folders and configuration files. It even launches the analyses after creation

Usage:
  $scriptname <options> input

* Input (mandatory parameters): *
 -pipe name		Pipeline to run (mapping, smallindels, all, wt)
 -seq type		Type of sequencing probes (PE -pair ends-, MP -mate pairs-, F3)
 -csfasta file	File containing F3 reads (mandatory if performing mapping analysis)
 -qual file		File containing F3 quality values (mandatory if performing mapping analysis)
 -bam file		Bam alignment file, obtained from mapping analysis (mandatory if performing smallindels analysis)

* Options: *
 -outdir output_path    Path where the results will be output (default path is <cwd>/bioscope_YYYYMMDD, where <cwd> is the current working directory)
 -csfasta2 file         File containing F5 or R3 reads
 -qual2 file            File containing F5 or R3 quality values
 -reference file        Fasta reference file (for its default value check /share/apps/scripts/resecuenciacion_dirigida/bioscope_templates/config.xml. Please consider that at the moment there are two different references, one for whole-transcriptome analysis and other for the rest)
 -cmap file             Cmap reference file (for its default value check /share/apps/scripts/resecuenciacion_dirigida/bioscope_templates/config.xml)
 -gtf file              Gtf file for whole transcriptome analysis (for its default value check /share/apps/scripts/resecuenciacion_dirigida/bioscope_templates/config.xml)
 -filter file           Reference filter file for whole transcriptome analysis (for its default value check /share/apps/scripts/resecuenciacion_dirigida/bioscope_templates/config.xml)
 -F3RL ##               F3 read length (default 75, unless diferrent set in /share/apps/scripts/resecuenciacion_dirigida/bioscope_templates/config.xml. Allowed values 75,50)
 -F5RL ##               F5 (or R3 if mate-pairs) read length (default 35, unless different set in /share/apps/scripts/resecuenciacion_dirigida/bioscope_templates/config.xml. Allowed values 35,25)
 -config config.xml     Default configuration XML file, containing all the default parameters (default is /share/apps/scripts/resecuenciacion_dirigida/bioscope_templates/config.xml)
 -run                   Runs the analysis inmediately after creation of folders and conf files (default won't run, unless different set in /share/apps/scripts/resecuenciacion_dirigida/bioscope_templates/config.xml)

* Output: *
 In the output_path are created the folders (config, intermediate, log, output, tmp) and conf files (F3.ini, F5.ini, etc...)

* Examples of running: *
 \$$scriptname -pipe mapping -seq PE -csfasta F3.csfasta -qual F3_QV.qual -csfasta2 F5.csfasta -qual2 F5_QV.qual -outdir arbol_mola
 \n";
 
	exit(1);
}

sub check_params{
# Function to check the parameters set into the files to be parsed.
# Please note that all parameters must have "$" at the beginning and "/" at the end (or simply be at the end of the line)

	# Arguments (please notice that $csfasta may also be $bam when calling function):
	my ( $param, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $date_time) = @_;

	if ($param eq "outdir") {return abs_path($outdir);}
	elsif ($param eq "csfasta" || $param eq "bam") {return abs_path($csfasta);}
	elsif ($param eq "qual") {return abs_path($qual);}
	elsif ($param eq "csfasta2") {return abs_path($csfasta2);}
	elsif ($param eq "qual2") {return abs_path($qual2);}
	elsif ($param eq "reference") {return $reference;}
	elsif ($param eq "filter") {return $filter;}
	elsif ($param eq "gtf") {return $gtf;}
	elsif ($param eq "cmap") {return $cmap;}
	elsif ($param eq "F3RL") {return $F3RL;}
	elsif ($param eq "F5RL") {return $F5RL;}
	elsif ($param eq "date-time") {
		if ($pipe eq "mapping") {return $date_time . "_Paired_End_";}
		elsif ($pipe eq "smallindels") {return $date_time . "_Find_Small_InDels";}
		elsif ($pipe eq "wt") {return $date_time . "_WT_Paired_End";}
	}
	elsif ($param eq "csfasta_name.ma"){
		(my $basename, my $dir, my $ext) = fileparse($csfasta, qr/\.[^.]*/);
		return $basename . $ext . ".ma";
	}
	elsif ($param eq "csfasta2_name.ma"){
		(my $basename, my $dir, my $ext) = fileparse($csfasta2, qr/\.[^.]*/);
		return $basename . $ext . ".ma";
	}
	else{
		print "ERROR: Parameter " . $param . " not found in list of allowed parameters!!!\n";
		exit();
	}
}


sub parse_template{
# Function to parse the contents of a template file ($template) and create the destination file ($config)
	
	# Arguments (please notice that $csfasta may also be $bam when calling function):
	my ( $config, $template, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $date_time) = @_;
	
	# Open config file as output:
	open CONFIG, '>', $config or die ("ERROR: Impossible to create $config!! $!");
	
	# Parse the template files and substitute with the given values where needed
	open (TEMPLATE, "<", $template) or die ("ERROR: Impossible to open $template!! $!");
	while (<TEMPLATE>) {
		my $line = $_;
		chomp($line);
		my $dollars_in_odd = 1;
		if (substr($line,0,1) eq '$'){
			$dollars_in_odd = 0;
		}
		my @params = split(/\$/, $line);
		
		my $counter = 0;
		foreach (@params){
			my $param = $_;
			if ($counter >= $dollars_in_odd){
				my $last_with = 0;
				if (substr($param,0,1) eq "/"){print CONFIG "/";}
				if (substr($param,-1,1) eq "/"){$last_with = 1;}
				
				my @folders = split(/\//,$param);
	
				my $first = 0;
				foreach (@folders){
					my $folder = $_;
					if ($first == 0) {
						my $param_to_print = check_params($folder ,$outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $date_time);
						print CONFIG $param_to_print;
					}
					else {
						print CONFIG "/" . $folder;
					}
					$first ++;
				}
				if ($last_with == 1){
					print CONFIG "/";
				}
			}
			else{
				print CONFIG $param;
			}
			$counter ++;
		}
		print CONFIG "\n";
	}
	close CONFIG;
	close TEMPLATE;
}


# **************************************************************************************************
# **************************************************************************************************

# -----------------------
# -- global parameters --
# -----------------------
# our $fooGlobalVar = $ENV{"PATH"};

# ------------------------
# -- default parameters --
# ------------------------
# Load default configuration from XML file:
my $xml = new XML::Simple;
my $data = XMLin("/share/apps/scripts/utils/templates/bioscope_templates/config.xml");

# print Dumper($data); #uncomment for showing contents of configuration stored in $data

my $yyyymmdd = strftime("%Y%m%d", localtime);
my $set_reference = 1;
my $reference;
my $filter = $data->{filter};
my $gtf = $data->{gtf};
my $cmap = $data->{cmap};
my $outdir = cwd() . "/bioscope_" . $yyyymmdd;
my $pipe = "";
my $seq = "";
my $csfasta = "";
my $qual = "";
my $bam = "";
my $csfasta2 = "";
my $qual2 = "";
my $F3RL = $data->{F3RL};
my $F5RL = $data->{F5RL};
my $run = $data->{run};

my $template_F3;
my $template_F5;
my $template_paired;
my $template_smallindels;
my $template_wt_F3rescue;
my $template_wt_F5rescue;
my $template_run1;
my $template_run2;
my $template_plan;

# -----------------------------
# -- Name, version and usage --
# -----------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);
print "$basename$ext v0.2
Script for generating the BioScope analysis folders and configuration files. It even launches the analyses after creation · Copyright 2012 by Arbol 
***************************************************************************************************************************************************\n";

# If no arguments are passed to the script, the script usage gets displayed
if ($#ARGV+1 < 1) {
	usage( $basename.$ext );
}

# ------------------------------------------
# -- Setting input options and parameters --
# ------------------------------------------
while ($#ARGV >= 0) {
	my $opt   = shift @ARGV;
	if ($opt =~ /^\-/) {
		if ($opt eq "-outdir") {
			$outdir = shift @ARGV;
        }
		elsif ($opt eq "-pipe") {
			$pipe = shift @ARGV;
        }
        elsif ($opt eq "-seq") {
			$seq = shift @ARGV;
        }
        elsif ($opt eq "-csfasta") {
			$csfasta = shift @ARGV;
        }
        elsif ($opt eq "-qual") {
			$qual = shift @ARGV;
        }
        elsif ($opt eq "-bam") {
			$bam = shift @ARGV;
        }
        elsif ($opt eq "-csfasta2") {
			$csfasta2 = shift @ARGV;
        }
		elsif ($opt eq "-qual2") {
			$qual2 = shift @ARGV;
        }
        elsif ($opt eq "-reference") {
			$reference = shift @ARGV;
			$set_reference = 0;
        }
        elsif ($opt eq "-cmap") {
			$cmap = shift @ARGV;
        }
        elsif ($opt eq "-F3RL") {
			$F3RL = shift @ARGV;
        }
        elsif ($opt eq "-F5RL") {
			$F5RL = shift @ARGV;
        }
        elsif ($opt eq "-run") {
			$run = 1;
        }
		else {
			print "Bad option: $opt\n";
			print "Call $basename$ext with no parameters for help.\n";
			exit(1);
		}
	}
}
if ($pipe ne "mapping" && $pipe ne "smallindels" && $pipe ne "all" && $pipe ne "wt") {
	print "Not allowed value for parameter pipeline in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}
if ($seq ne "PE" && $seq ne "MP" && $seq ne "F3") {
	print "Not allowed value for parameter sequencing type in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}
if (($csfasta eq "" || $qual eq "") && $pipe eq "mapping") {
	print "Please be sure to define all the reads and quality values files in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}
if ($bam eq "" && $pipe eq "smallindels") {
	print "Please be sure to define the alignment file (.bam) in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}
if ((($seq eq "PE" || $seq eq "MP") && ($csfasta2 eq "" || $qual2 eq "")) && $pipe eq "mapping") {
	print "Please be sure to define all the reads and quality values files in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}
if (($F3RL != 75 && $F3RL != 50) || ($F5RL != 35 && $F5RL != 25)) {
	print "Unallowed values for read lengths in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}
if ($set_reference){
	if ($pipe eq "wt"){
		$reference = $data->{wt_reference};
	} else {
		$reference = $data->{reference};
	}
}
#if ($pipe eq "all"){
#	if ($run) {
#		my $command=`$basename.$ext -pipe mapping -seq $seq -csfasta $csfasta -qual $qual -csfasta2 $csfasta2 -qual2 $qual2 -F3RL $F3RL -F5RL $F5RL -cmap $cmap -reference $reference -outdir $outdir -run`;
#	}
#	else {
#		my $command=`$basename.$ext -pipe mapping -seq $seq -csfasta $csfasta -qual $qual -csfasta2 $csfasta2 -qual2 $qual2 -F3RL $F3RL -F5RL $F5RL -cmap $cmap -reference $reference -outdir $outdir`;
#	}
#	my $command=`$basename.$ext -pipe smallindels -seq $seq -bam $outdir/output/pairing/F3-F5-Paired.bam -F3RL $F3RL -F5RL $F5RL -cmap $cmap -reference $reference -outdir $outdir`;
#	exit(1);
#}

#################################################################################################3
# ---------------
# -- main body -- 
# ---------------
################################################################################################3

my $now = localtime;
print scalar($now) . "\n";

# Some further default configuration parameters are loaded hier:
if ($pipe eq "mapping"){
	$template_F3 = $data->{mapping}->{F3_ini};
	$template_F5 = $data->{mapping}->{F5_ini};
	$template_paired = $data->{mapping}->{Paired_End_ini};
	$template_plan = $data->{mapping}->{analysis_plan};
} elsif ($pipe eq "smallindels"){
	if ($seq eq "PE") {
		$template_smallindels = $data->{smallindels_PE}->{Find_Small_InDels_ini};
		$template_plan = $data->{smallindels_PE}->{analysis_plan};		
	}
	elsif ($seq eq "F3") {
		$template_smallindels = $data->{smallindels_F3}->{Find_Small_InDels_ini};
		$template_plan = $data->{smallindels_F3}->{analysis_plan};		
	}
} elsif ($pipe eq "wt"){
	if ($seq eq "PE") {
		$template_F3 = $data->{whole_transcriptome_PE}->{F3_ini};
		$template_F5 = $data->{whole_transcriptome_PE}->{F5_ini};
		$template_paired = $data->{whole_transcriptome_PE}->{Paired_End_ini};
		$template_plan = $data->{whole_transcriptome_PE}->{analysis_plan};
		$template_wt_F3rescue = $data->{whole_transcriptome_PE}->{F3_rescue};
		$template_wt_F5rescue = $data->{whole_transcriptome_PE}->{F5_rescue};
	}
	elsif ($seq eq "F3") {
		$template_F3 = $data->{whole_transcriptome_PE}->{Find_Small_InDels_ini};
		$template_plan = $data->{whole_transcriptome_PE}->{analysis_plan};		
	}
}
$template_run1 = $data->{general}->{run1};
$template_run2 = $data->{general}->{run2};

# Create output folders:
unless(-d $outdir){
	mkdir $outdir or die ("ERROR: Impossible to create folder $outdir!!");
}
my $folder = $outdir . "/config";
unless(-d $folder){
	mkdir $folder or die ("ERROR: Impossible to create folder $folder!!");
}
$folder = $outdir . "/intermediate";
unless(-d $folder){
	mkdir $folder or die ("ERROR: Impossible to create folder $folder!!");
}
$folder = $outdir . "/log";
unless(-d $folder){
	mkdir $folder or die ("ERROR: Impossible to create folder $folder!!");
}
$folder = $outdir . "/output";
unless(-d $folder){
	mkdir $folder or die ("ERROR: Impossible to create folder $folder!!");
}
$folder = $outdir . "/tmp";
unless(-d $folder){
	mkdir $folder or die ("ERROR: Impossible to create folder $folder!!");
}
my $hhmmss = strftime("%H%M%S", localtime);
if ($seq eq "PE" && ( $pipe eq "mapping" or $pipe eq "wt")){
	if ($pipe eq "mapping"){
		$folder = $outdir . "/config/" . $yyyymmdd . "-" . $hhmmss . "_Paired_End_";
	} else {
		$folder = $outdir . "/config/" . $yyyymmdd . "-" . $hhmmss . "_WT_Paired_End";
	}
	unless(-d $folder){
		mkdir $folder or die ("ERROR: Impossible to create folder $folder!!");
	}
	unless(-d $folder . "/F3"){
		mkdir $folder . "/F3" or die ("ERROR: Impossible to create folder $folder!!");
	}
	unless(-d $folder . "/F5"){
		mkdir $folder . "/F5" or die ("ERROR: Impossible to create folder $folder!!");
	}
}
elsif ($pipe eq "smallindels"){
	$folder = $outdir . "/config/" . $yyyymmdd . "-" . $hhmmss . "_Find_Small_InDels";
	unless(-d $folder){
		mkdir $folder or die ("ERROR: Impossible to create folder $folder!!");
	}
}
else { $folder = $outdir . "/tmp";}

if ($pipe eq "mapping"){
	parse_template($folder . "/F3/F3.ini", $template_F3, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
	parse_template($folder . "/F5/F5.ini", $template_F5, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
	parse_template($folder . "/Paired_End_.ini", $template_paired, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
	parse_template($folder . "/analysis.plan", $template_plan, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
} elsif ($pipe eq "smallindels"){
	parse_template($folder . "/Find_Small_InDels.ini", $template_smallindels, $outdir, $bam, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
	parse_template($folder . "/analysis.plan", $template_plan, $outdir, $bam, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
} elsif ($pipe eq "wt"){
	parse_template($folder . "/F3/F3.ini", $template_F3, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
	parse_template($folder . "/F5/F5.ini", $template_F5, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
	parse_template($folder . "/WT_Paired_End.ini", $template_paired, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
	parse_template($folder . "/analysis.plan", $template_plan, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
	parse_template($folder . "/F3rescue.ini", $template_wt_F3rescue, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
	parse_template($folder . "/F5rescue.ini", $template_wt_F5rescue, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
}

parse_template($folder . "/run1.sh", $template_run1, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
chmod 0744, $folder . "/run1.sh" or die "Couldn't chmod " . $folder . "/run1.sh: $!";
parse_template($folder . "/run2.sh", $template_run2, $outdir, $csfasta, $qual, $csfasta2, $qual2, $reference, $cmap, $filter, $gtf, $F3RL, $F5RL, $pipe, $yyyymmdd . "-" . $hhmmss);
chmod 0744, $folder . "/run2.sh" or die "Couldn't chmod " . $folder . "/run2.sh: $!";

# Run the analysis if parameter $run==1:
if ($run == 1){
	chdir (abs_path($folder));
	my $run_analysis = `./run1.sh`;
}

$now = localtime;
print "Script finished at $now\n";
#******************************************************************* 