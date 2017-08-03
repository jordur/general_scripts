#!/usr/bin/perl -w

=head1 LICENSE

  This script was developed by Sistemas Genomicos SL

=head1 CONTACT

  Please email comments or questions to:
  jm.rosa@sistemasgenomicos.com
  oscar.rodriguez@sistemasgenomicos.com

=cut

=head1 NAME

Customized Variant Effect Predictor - a script to predict the consequences of genomic variants

Version 1.0

by Juan Manuel Rosa (jm.rosa@sistemasgenomicos.com) & Arbol (oscar.rodriguez@sistemasgenomicos.com)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use Parallel::ForkManager;


sub debug ($);
sub get_out_file_handle ($);
sub getting_annotation ($$);
sub getting_completation ($$);


# configure from command line opts
my $config = &configure(scalar @ARGV);

# run the main sub routine
&main($config);

# this is the main sub-routine - it needs the configured $config hash

sub main {
	my $config = shift;
	my @chrs = split (/,/, &get_number_chr ($config));
	my $pm = new Parallel::ForkManager($config->{threads});
	$pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; } );
	$pm->run_on_start( sub { my ($pid,$ident)=@_; } );
	$pm->run_on_wait( sub {}, 0.5 );

	foreach my $chr (@chrs) {
		my $pid = $pm->start($chr) and next;
		print "Annotating variants for chr$chr... ".localtime()."\n";
		my $anno_job_id = &getting_annotation ($config, $chr);
		print "Completing information for chr$chr... ".localtime()."\n";
		my $comp_job_id = &getting_completation ($config, $chr);
		$pm->finish;
	}
	$pm->wait_all_children;
}

debug("Finished!") if defined $config->{verbose};
exit;

sub getting_completation ($$) {
	my ($config, $chr) = @_;

	my $species = $config->{species};
	my $input = $config->{tmp_file}.".$chr";
	my $vcf = $config->{input_file};
	my $output = $config->{out_file}.".$chr";
	my $host = $config->{host};
	my $port = $config->{port};
	my $user = $config->{user};
	my $password = $config->{password};
	my $pass = '';
	my $dbversion = $config->{db_version};
	
	my $output_file = &get_out_file_handle($output.".job");
	
	if ($password) { 
		$pass = "--password $password";
	}

	print $output_file "variant_information_completer_64.pl -i $input -o $output --vcf $vcf --species $species --host $host --port $port --user $user $pass --db_version $dbversion --samples $config->{samples} --freq $config->{freq} $config->{type} > $output.complete.log\n"; 
	print "Launching command:\n$output.job\n".localtime()."\n";

	my $job = `bash $output.job`;
	wait;
	return 1;
}

sub getting_annotation ($$) {
	my ($config, $chr) = @_;
	
	my $species = $config->{species};
	my $input = "<(awk \'\$1 == \"chr$chr\"\' $config->{input_file})";
	my $output = $config->{tmp_file}.".$chr";
	my $host = $config->{host};
	my $port = $config->{port};
	my $user = $config->{user};
	my $password = $config->{password};
	my $pass = '';
	my $dbversion = $config->{db_version};
	
	my $output_file = &get_out_file_handle("chr$chr.job");
	
	if ($password) {
		$pass = "--password $password";
	}

	if ($species eq 'human' || $species eq 'homo_sapiens') {
		print $output_file "variant_effect_predictor_64.pl -i $input --format vcf -o $output --sift b --polyphen b --condel b --regulatory --hgnc  --hgvs --no_intergenic --force_overwrite --species $species --host $host --port $port --user $user $pass $config->{check} --db_version $dbversion > $output.anno.log\n"; 
		print "Launching command:\nchr$chr.job\n".localtime()."\n";
	}
	else {
        print $output_file "variant_effect_predictor_64.pl -i $input --format vcf -o $output --regulatory --hgnc  --hgvs --no_intergenic --force_overwrite --species $species --host $host --port $port --user $user $pass $config->{check} --db_version $dbversion > $output.anno.log\n"; 
		print "Launching command:\nchr$chr.job\n".localtime()."\n";
	}
	my $job = `bash chr$chr.job`;
	wait;
	return 1;
}

sub get_number_chr (){
	my ($config) = @_;
	my $last_chr = 'N\A';
	my $file_handle = $config->{in_file_handle};
	my $chrs;
	while (my $line = <$file_handle>) {
		chomp $line;
		my $chr = (split (/\s+/, $line))[0];
		if ($chr eq "#CHROM") { next;}
		$chr =~ s/chr//;
		if ($chr eq $last_chr) { next; }
		else {
			$last_chr = $chr;
			if (defined $chrs) {
				$chrs .= ",$chr";
			}
			else {
				$chrs = $chr;
			}
		}
	}
	return $chrs;
}

sub configure {
	my $args = shift;
	
	my $config = {};
	
	GetOptions(
		$config,
		'help',
		
		# input options,
		'config=s',
		'input_file=s',
		'bed_file=s',
		'samples=s',
		
		# DB options
		'species=s',
		'host=s',
		'user=s',
		'port=s',
		'password=s',
		'pair',
		'threads=i',
		'db_version=i',
		
		# runtime options
		'check_existing',
		'freq=f',
		
		# output options
		'output_file=s',
		'verbose',
		'quiet',
		'hgnc',
		'hgvs',
		'sift=s',
		'polyphen=s',
		'condel=s',
	);
	
	# print usage message if requested or no args supplied
	if(defined($config->{help}) || !$args) {
		&usage;
		exit(0);
	}
	
	# summarise options if verbose
	if(defined $config->{verbose}) {
		my $header =<<INTRO;
#---------------------------------------------#
# CUSTOMIZED ENSEMBL VARIANT EFFECT PREDICTOR #
#---------------------------------------------#

version 1.0

By sistemas Genomicos

Configuration options:

INTRO
		print $header;
		
		my $max_length = (sort {$a <=> $b} map {length($_)} keys %$config)[-1];
		
		foreach my $key(sort keys %$config) {
			print $key.(' ' x (($max_length - length($key)) + 4)).$config->{$key}."\n";
		}
		
		print "\n".("-" x 20)."\n\n";
	}
	
	# set defaults
	
	$config->{species} ||= "homo_sapiens";
	$config->{host}    ||= 'ensembldb.ensembl.org';
	$config->{port}    ||= 5306;
	$config->{user}    ||= 'anonymous';
	$config->{db_version}	||= 64;

	if (defined $config->{pair}) {
		$config->{type} = '--pair';
	}
	else {
		$config->{type} = '';
	}
	if (not defined $config->{samples}) {
		&usage;
		exit(0);
	}
	
	if (defined $config->{check_existing}){
		$config->{check} = '--check_existing';
	}
	else {
		$config->{check} = '';
	}
	
	$config->{freq} ||= 1;
	
	
	$config->{threads}    ||= 2;

	my $input = $config->{input_file};
	$input =~ s/\.\.\/+//g;
	my @output_name = split (/\./, $input); 
	$config->{output_file} ||= $output_name[0]."_annotated";
	$config->{out_file} = $config->{output_file};
	$config->{tmp_file} = $config->{output_file}.".tmp";
	$config->{log_file} = $config->{output_file}.".log";

	# get input file handles
	$config->{in_file_handle} = &get_in_file_handle($config);

	return $config;
}

# gets file handle for input
sub get_in_file_handle {
    my $config = shift;

    # define the filehandle to read input from
    my $in_file_handle = new FileHandle;
    
    if(defined($config->{input_file})) {
        
        # check defined input file exists
        die("ERROR: Could not find input file ", $config->{input_file}, "\n") unless -e $config->{input_file};
        
        if($config->{input_file} =~ /\.gz$/){
            $in_file_handle->open($config->{compress}." ". $config->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
        else {
            $in_file_handle->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
    }
    
    # no file specified - try to read data off command line
    else {
        $in_file_handle = 'STDIN';
        debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
    }
    
    return $in_file_handle;
}



sub usage {
	my $usage =<<END;
#---------------------------------------------#
# CUSTOMIZED ENSEMBL VARIANT EFFECT PREDICTOR #
#---------------------------------------------#

version 1.0

By Sistemas Genomicos

Usage:
perl variant_effect_predictor_Ensembl64.pl [arguments]

Options
=======

--help                 Display this message and quit

-i | --input_file      Input file - vcf format.

-o | --output_file     Output file. Write to STDOUT by specifying -o STDOUT - this
                       will force --quiet [default: "variant_effect_output.txt"]
					   
--samples				Samples ID (Comma separated for multiple). Mandatory
                       
--species [species]    Species to use [default: "human"]

--pair					Tumour/normal pair (default independant sample)
--threads				Number of threads to use (default 2)
                                     
--check_existing		Reporting existing variantions
--freq					Frequency value to filter known alleles (default "1" -> no filter)
                               
--db_version			Ensembl version (default 64)
--host                 Manually define database host [default: "ensembldb.ensembl.org"]
--user 		          Database username [default: "anonymous"]
--port                 Database port [default: 5306]
--password             Database password [default: no password]

--bed_file		   BED file with transcripts ID

END

	print $usage;
}

# gets time
sub get_time() {
	my @time = localtime(time());

	# increment the month (Jan = 0)
	$time[4]++;

	# add leading zeroes as required
	for my $i(0..4) {
		$time[$i] = "0".$time[$i] if $time[$i] < 10;
	}

	# put the components together in a string
	my $time =
 		($time[5] + 1900)."-".
 		$time[4]."-".
 		$time[3]." ".
		$time[2].":".
		$time[1].":".
		$time[0];

	return $time;
}

# prints debug output with time
sub debug ($) {
	my $text = (@_ ? (join "", @_) : "No message");
	my $time = get_time;
	
	print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}

#get file for output
sub get_out_file_handle ($) {
	my $file = shift;
	
	# define filehandle to write to
	my $out_file_handle = new FileHandle;
	$out_file_handle->open(">".$file) or die("ERROR: Could not write to output file ", $file, "\n");

	return $out_file_handle;
}
