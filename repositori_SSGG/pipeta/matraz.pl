#!/usr/bin/env perl
##############################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Matraz v0.3 JSON generator for pipeta
# @Author: gmarco
# @Contributors: Arbol, Jordi
# @Last Modification: 23/01/2014
###############################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use SOAP::Lite;
use Term::ANSIColor;
use file_handle;
use json_handle;
use Data::Dumper;

# Global variables
my $gluster = $ENV{'PROJECTS'};
my $gluster2 = $ENV{'PROJECTS2'};
my $config_json = $ENV{'MATRAZ_CONFIG'};
my ($bf,$analysis_option, $lv_id, $user_aligner, $nocheck, $pipeta);

# WS variables
my $ws_user = 'SGNET.LOCAL\\auren';
my $ws_pass = 'C0intec2012';
my $ws_host = 'sg02.sgnet.local:7047';	# always include the port

main();

sub main{
	my ($json_config_filename, $trimming_string);
	
	#Get user options, BF_dir and and target analysis type
	($bf,$analysis_option, $lv_id, $user_aligner, $nocheck)  = menu();
	
	#Check if project is located on gluster or gluster2.
	my $dir = check_project_location($gluster, $gluster2, $bf);
	
	#Get all information from matraz target json config file
	my ($chr_split, $target, $target_chrs, $capture, $three_prime_primers, $five_prime_primers, $haloplex_primer1, $haloplex_primer2, $haloplex_primer1_rc, $haloplex_primer2_rc, $seq_platform, $sequencer, $size, $capture_system, $prefix, $pipeline, $aligner) = config_reader($dir, $bf, $analysis_option, $lv_id);
	
	#Get hash with samples
	my %all_samples = dir_explorer($dir, $bf, $prefix);
	
	#Get sample string separated by comma
	my $samples = sample_list(%all_samples);
	
	print "▀▄▀▄▀▄ мαтяαz ν0.3 ▄▀▄▀▄▀\n";
	print "\nMatraz is working...\n";
	print "BF folder: $dir/$bf\n";
	print "Analysis type: $analysis_option\n";
	print "Pipeline: $pipeline\n";
	print "Mapping module: $aligner\n";
	
	if ($lv_id ne ""){
		print "Filtering by $lv_id\n";
	}
	
	#Remove config directory and re-create it.
	#`rm -rf $dir/$bf/config`;
	`mkdir -p $dir/$bf/config`;
	
	#Open output json config file
	$json_config_filename = "$dir/$bf/config/$bf"."_$analysis_option";
	if ($lv_id ne ""){
		$json_config_filename .= "_$lv_id";
	}
	my %json;
	
	#Choose aligner
	$aligner = choose_aligner($user_aligner, $aligner);

	# Define json fields:
	header(\%json, $dir, $bf);
	essay(\%json, $dir, $bf, $analysis_option, $pipeline, $aligner, $samples, $chr_split, $target, $target_chrs, $capture, $three_prime_primers, $five_prime_primers, $haloplex_primer1, $haloplex_primer2, $haloplex_primer1_rc, $haloplex_primer2_rc, $trimming_string, $lv_id);
	samples(\%json, $dir, $bf, $seq_platform, $sequencer, $size, $capture_system, %all_samples);
	
	# Create json
	create_json(\%json, $json_config_filename);
	
	print "#########################################\n";
	print "DONE: JSON Config file created for project $bf in: $json_config_filename.json\n";
	
	if ($pipeta){
		print "\n";
		print "INFO: Launching PipeTA...\n";
		print "pipeta.pl $json_config_filename.json\n";
		`pipeta.pl $json_config_filename.json`;
	}
}

sub menu {
	my $available_analysis_list;
	$lv_id = "";
	$user_aligner = "";
	$nocheck = 0;
	$pipeta = 0;
	
	GetOptions ("bf=s" => \$bf, 'a=s' => \$analysis_option, 'm:s' => \$user_aligner, 'lv:s' => \$lv_id, 'nocheck!' => \$nocheck, 'pipeta!' => \$pipeta);
	
	if ($lv_id){
		#If LV user input starts with a number add LV to string.
		if  ($lv_id =~ /^\d/){
			$lv_id = "LV$lv_id";
		}
	}
	
	if (not defined $bf or not defined $analysis_option){
		$available_analysis_list = analisis_available_list();
		help($available_analysis_list);
	}
	
	#Menu options to lower case
	$analysis_option = lc $analysis_option;
	$user_aligner = lc $user_aligner;
	
	#Check if mapping module is correct.
	available_list_aligners($user_aligner);
	
	return $bf,$analysis_option, $lv_id, $user_aligner, $nocheck, $pipeta;
}

sub help {
	my $available_analysis_list = $_[0];
	print "
мαтяαz ν0.3
           .   Before launching Matraz check:
      o .  .      - BF directory exists.
       . O o .    - rawdata folder exists inside the specified BF directory.
     O  .  .      - Each sample has its own directory inside rawdata directory.
       * O.
       o * .      Expected directory structure:
       o O.       /share/gluster/OnGoing/BF123_example
       o o .                     └── rawdata
       O o .                         ├── Sample_EX5109_0075
       o o .                         └── Sample_EX5165-III-1           
     aaaaaaaa   
     \"8. o 8\"     You must specify the following options:
      8 O .8         -bf BF Bioinformatics code (directory name) ie: -bf BF123_example
      8 o .8         -a  Analisis type: ie: -a Exome51
      8. O 8 ";
	print colored ("         Available analysis: $available_analysis_list\n", 'bold green');
	print "      8. O 8         -m mapping module (default is BWA) Be sure that aligner is compatible with pipeline ! ie: -m bwa is not compatible with DNA-reseq_trimming !
      8. O 8         -nocheck    Use this option for not performing the sample name integrity checks (useful for data coming from outside SSGG)
      8.o  8         -pipeta     Use this option for directly launching the relating bioinformatics analysis by means of PipeTA
      8. 0 8      
      8. O 8        LV filtering for panels (optional):
      8. O 8         -lv lv_number ie: -lv 1601 or -lv LV1601
      8. O 8 
      8 o. 8        Usage example:
   ,adP O .Yba,     matraz.pl -bf BF123_my_project -a onco2
  dP\". O  o  \"Yb
 dP' O . o .O `Yb\n";
 print " 8^^^^^^^^^^^^^^8   ";
 print colored "MATRAZ WILL VERIFY SAMPLE TYPE BASED ON DIRECTORY PREFIX.\n", 'bold red';
 print " 8    EnJoY !   8   ";
 print colored "HOWEVER YOU MUST CHECK JSON CONFIG BEFORE LAUNCHING PIPETA !\n", 'bold red';
 print " Yb,          ,dP
  \"Ya,______,aP\"
    `\"\"\"\"\"\"\"\"'
";

exit;

}

sub check_project_location{
	my ($gluster, $gluster2, $bf) = @_;
	my $dir;
	if(-d($gluster."/".$bf)){
		return $dir = $gluster;
	}
	if (-d($gluster2."/".$bf)){
		return $dir = $gluster2;
	}
	die "No $bf project found in Gluster or Gluster2!\n";
}

sub available_list_aligners {
	my $mapping_module = shift;
	if ($mapping_module !~ /bwa|bowtie2|novoalign||/){
		die "ERROR: $mapping_module is not a valid aligner. Available aligners are: bwa, bowtie2 & novoalign.\n"; 
	}
}
 
sub analisis_available_list {
	my (@available_list, $available_list, @gene_filter_list, %gene_filter_list);
	my %analysis_type = get_from_json($config_json);
	
	#Obtain a list of available target regions in case user is not selecting any valid option.
	foreach my $key (keys %analysis_type ) {
		foreach my $target_available (sort keys %{$analysis_type{$key}}) {
			push (@available_list, $target_available);
		}
	}
	
	$available_list = join(", ",@available_list);
	return $available_list;
}


sub config_reader {
	my ($dir, $bf, $option, $lv_id) = @_;
	my (@available_list, $available_list, @gene_filter_list, %gene_filter_list);
	my ($chr_split, $target, $target_chrs, $capture, $prefix, $panel_capture, 
	    $target_genes_dir, $three_prime_primers, $five_prime_primers,
	    $haloplex_primer1, $haloplex_primer2, $haloplex_primer1_rc, $haloplex_primer2_rc,
	    $seq_platform, $sequencer, $size, $capture_system, $pipeline, $aligner);
	my %analysis_type = get_from_json($config_json);
	
	$three_prime_primers = "";
	$five_prime_primers = "";
	$haloplex_primer1 = "";
	$haloplex_primer2 = "";
	$haloplex_primer1_rc = "";
	$haloplex_primer2_rc = "";
	
	#Obtain a list of available target regions in case user is not selecting any valid option.
	foreach my $key (keys %analysis_type ) {
		foreach my $target_available (sort keys %{$analysis_type{$key}}) {
			push (@available_list, $target_available);
		}
	}
	
	$available_list = join(", ",@available_list);
	die "ERROR: $option doesn't exist in Matraz JSON config file !\nAvailable analisis list: $available_list.\n" if not exists ${analysis_type{'target_reference'}{$option}};

	#If no lv_id is specified
	if ($lv_id eq ""){
		$chr_split = ${analysis_type{'target_reference'}{$option}{chr_split}};
		$target = ${analysis_type{'target_reference'}{$option}{target}};
		$target_chrs = ${analysis_type{'target_reference'}{$option}{target_chrs}};
		$capture = ${analysis_type{'target_reference'}{$option}{capture}};
		$prefix = ${analysis_type{'target_reference'}{$option}{prefix}};
		$capture_system = ${analysis_type{'target_reference'}{$option}{capture_system}};
		$seq_platform = ${analysis_type{'target_reference'}{$option}{seq_platform}};
		$sequencer = ${analysis_type{'target_reference'}{$option}{sequencer}};
		$size = ${analysis_type{'target_reference'}{$option}{size}};
		$pipeline = ${analysis_type{'target_reference'}{$option}{pipeline}};
		$aligner = ${analysis_type{'target_reference'}{$option}{aligner}};
		
		#Check if multiplicom primers defined.
		if (defined ${analysis_type{'target_reference'}{$option}{three_prime_primers}}){
			$three_prime_primers = ${analysis_type{'target_reference'}{$option}{three_prime_primers}};
		}
		if (defined ${analysis_type{'target_reference'}{$option}{five_prime_primers}}){
			$five_prime_primers = ${analysis_type{'target_reference'}{$option}{five_prime_primers}};
		}
		
		#Check if haloplex adapters defined.
		if (defined ${analysis_type{'target_reference'}{$option}{haloplex_primer1}}){
			$haloplex_primer1 = ${analysis_type{'target_reference'}{$option}{haloplex_primer1}};
		}
		if (defined ${analysis_type{'target_reference'}{$option}{haloplex_primer2}}){
			$haloplex_primer2 = ${analysis_type{'target_reference'}{$option}{haloplex_primer2}};
		}
		if (defined ${analysis_type{'target_reference'}{$option}{haloplex_primer1_rc}}){
			$haloplex_primer1_rc = ${analysis_type{'target_reference'}{$option}{haloplex_primer1_rc}};
		}
		if (defined ${analysis_type{'target_reference'}{$option}{haloplex_primer2_rc}}){
			$haloplex_primer2_rc = ${analysis_type{'target_reference'}{$option}{haloplex_primer2_rc}};
		}
		
		
		
	}
	
	#If LV is specified
	else{
		die "ERROR: You can't filter by LV if you're not using a panel design!\n" if ${analysis_type{'target_reference'}{$option}{panel_id}} eq "no";
		
		%gene_filter_list = get_lv_genes($lv_id);
		
		#Check that panel code from LV correspond to panel cod in matraz_config.json
		foreach my $gene_id (sort keys %gene_filter_list) {
			#die "ERROR: Panel ID doesn't match with your LV Panel ID.\n" if "$gene_filter_list{$gene_id}" ne "${analysis_type{'target_reference'}{$option}{panel_id}}";
			push (@gene_filter_list, $gene_id);
		}
		
		$target_genes_dir = ${analysis_type{'target_reference'}{$option}{bed_gene_dir}};
		
		#Get BED target regions for LV genes.
		$target = lv_target_region($lv_id, $dir, $bf, $option, $target_genes_dir, @gene_filter_list);
		
		#Get BED capture regions intersecting panel capture regions with target regions for LV.
		$panel_capture = ${analysis_type{'target_reference'}{$option}{capture}};
		$capture = lv_capture_region($lv_id, $dir, $bf, $option, $target, $panel_capture);
		
		$chr_split = ${analysis_type{'target_reference'}{$option}{chr_split}};
		$target_chrs = ${analysis_type{'target_reference'}{$option}{target_chrs}};
		
		$prefix = ${analysis_type{'target_reference'}{$option}{prefix}};
		$capture_system = ${analysis_type{'target_reference'}{$option}{capture_system}};
		$seq_platform = ${analysis_type{'target_reference'}{$option}{seq_platform}};
		$sequencer = ${analysis_type{'target_reference'}{$option}{sequencer}};
		$size = ${analysis_type{'target_reference'}{$option}{size}};
		$pipeline = ${analysis_type{'target_reference'}{$option}{pipeline}};
		$aligner = ${analysis_type{'target_reference'}{$option}{aligner}};
		
		#Check if multiplicom primers defined.
		if (defined ${analysis_type{'target_reference'}{$option}{three_prime_primers}}){
			$three_prime_primers = ${analysis_type{'target_reference'}{$option}{three_prime_primers}};
		}
		if (defined ${analysis_type{'target_reference'}{$option}{five_prime_primers}}){
			$five_prime_primers = ${analysis_type{'target_reference'}{$option}{five_prime_primers}};
		}
		
		#Check if haloplex adapters defined.
		if (defined ${analysis_type{'target_reference'}{$option}{haloplex_primer1}}){
			$haloplex_primer1 = ${analysis_type{'target_reference'}{$option}{haloplex_primer1}};
		}
		if (defined ${analysis_type{'target_reference'}{$option}{haloplex_primer2}}){
			$haloplex_primer2 = ${analysis_type{'target_reference'}{$option}{haloplex_primer2}};
		}
		if (defined ${analysis_type{'target_reference'}{$option}{haloplex_primer1_rc}}){
			$haloplex_primer1_rc = ${analysis_type{'target_reference'}{$option}{haloplex_primer1_rc}};
		}
		if (defined ${analysis_type{'target_reference'}{$option}{haloplex_primer2_rc}}){
			$haloplex_primer2_rc = ${analysis_type{'target_reference'}{$option}{haloplex_primer2_rc}};
		}
		

	}
	
	return $chr_split, $target, $target_chrs, $capture, $three_prime_primers, $five_prime_primers, $haloplex_primer1, $haloplex_primer2, $haloplex_primer1_rc, $haloplex_primer2_rc, $seq_platform, $sequencer, $size, $capture_system, $prefix, $pipeline, $aligner;
}

sub is_folder_empty {
	my $dirname = shift;
	opendir(my $dh, $dirname) or die "ERROR: $dirname is not a directory or not exsists !!\n";
	return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}

sub get_lv_genes {

	my $lv = shift;
	my (@genes, @lv_genes);
	my %genes = ();
	my $result = SOAP::SOM->new;
	my $endpoint2 = "http://sg02.sgnet.local:7047/DynamicsNAV/WS/Sistemas%20Gen%C3%B3micos,%20S.A/Page/LVGEN";
	
	my $soap = SOAP::Lite
		-> uri ('urn:microsoft-dynamics-schemas/page/lvgen')
		-> on_action( sub { return 'urn:microsoft-dynamics-schemas/page/lvgen:ReadMultiple' } )
		-> proxy($endpoint2, keep_alive => 1, credentials => [$ws_host, '', $ws_user, $ws_pass]);
	
	
	$result = $soap->ReadMultiple(SOAP::Data->name('filter' => \SOAP::Data->value(
					SOAP::Data->new(name => 'Field', value => 'Cod_item', type => 'tns:LVGEN_Fields'),
					SOAP::Data->new(name => 'Criteria', value => $lv, type => 'xsd:string')))->type('tns:LVGEN_Filter'),
					
		SOAP::Data->new(name => 'setSize', value => '0', type => 'xsd:int')
		)->result;
	
	#If WS returns no genes, die.
	die "ERROR: $lv doesn't exist or has no genes !\n" if $result eq "";
	
	#Check if there's only one gene in return from WS	
	if (ref $result->{LVGEN} eq 'HASH'){
		die "ERROR: GENE ID MISSING !! NO HGNC for $lv\n" if not exists ($result->{LVGEN}{'Cod_gen'});
		die "ERROR: GENE ID MISSING !! NO HGNC for $lv\n" if ($result->{LVGEN}{'Cod_gen'} eq "");
		$genes{$result->{LVGEN}{'Cod_gen'}} = $result->{LVGEN}{'Cod_gen'};
		#$genes{$result->{LVGEN}{'Cod_gen'}} = $result->{LVGEN}{'Panel_Prod'};
	}
	
	else{
		@lv_genes = @{$result->{LVGEN}};
		foreach my $j ( @lv_genes ) {
			die "ERROR: GENE ID MISSING !! NO HGNC for $lv\n" if not exists ($j->{'Cod_gen'});
			die "ERROR: GENE ID MISSING !! NO HGNC for $lv\n" if ($j->{'Cod_gen'} eq "");
			#die "ERROR: One or more genes in $lv have no associated panel code !\n" if not defined $j->{'Panel_Prod'};
			#$genes{$j->{'Cod_gen'}} = $j->{'Panel_Prod'};
			$genes{$j->{'Cod_gen'}} = $j->{'Cod_gen'};
		}
	} 
	#Return a hash with gene name and panel ID associated to gene.
	return %genes;
}

sub get_lv_validation{

	my $lv = shift;
	my $result = SOAP::SOM->new;
	my $endpoint2 = "http://sg02.sgnet.local:7047/DynamicsNAV/WS/Sistemas%20Gen%C3%B3micos,%20S.A/Page/Productos";
	
	my $soap = SOAP::Lite
		-> uri ('urn:microsoft-dynamics-schemas/page/productos')
		-> on_action( sub { return 'urn:microsoft-dynamics-schemas/page/productos:Read' } )
		-> proxy($endpoint2, keep_alive => 1, credentials => [$ws_host, '', $ws_user, $ws_pass]);
	
	$result = $soap->Read(SOAP::Data->name('Codigo_Producto' => \SOAP::Data->value(
				SOAP::Data->new(name => 'Codigo_Producto', value => $lv, type => 'xsd:string')))
	)->result;
}

sub lv_target_region {
	my ($lv_id, $dir, $bf, $option, $target_genes_dir, @genes) = @_;
	my @gene_bed_files;
	my $target_sequence_gene_merge;
	my $bedops_merge = "bedops --merge ";
	
	foreach my $gene (@genes){
		my $gene_file = "$target_genes_dir/$gene/target.bed";
		print "gene_file: $gene_file\n";
		die "ERROR: Missing gene target file for gene: $gene\nProbably the LV you're trying to filter doesn't belong to this panel !\n" if (!-e $gene_file);
		push (@gene_bed_files, $gene_file);
		$bedops_merge .= "$gene_file ";
	}
	
	#Clean target_seq for lv and re-create it.
	`rm -rf $dir/$bf/lv_regions`;
	`mkdir -p $dir/$bf/lv_regions`;
	
	my $lv_target_file = "$dir/$bf/lv_regions/$option";
	if ($lv_id ne ""){
		$lv_target_file .= "_$lv_id.csv";
	}
	open LV_TARGET_INFO, ">$lv_target_file";
	print LV_TARGET_INFO join("\n",@genes);
	close LV_TARGET_INFO;
	
	$target_sequence_gene_merge = "$dir/$bf/lv_regions/$option"."_$lv_id"."_target.bed";
	
	$bedops_merge .= "> $target_sequence_gene_merge.unsorted";
	print "bedops_merge: $bedops_merge\n";
	`$bedops_merge`;
	
	#After merge we sort file and remove tmp.
	`sed -e 's,chrX,chr23,' -e 's,chrY,chr24,' -e 's,chrM,chr25,' $target_sequence_gene_merge.unsorted | sort -k 1.4,1n -k 2,2n -k 3,3n | sed -e 's,chr23,chrX,' -e 's,chr24,chrY,' -e 's,chr25,chrM,' > $target_sequence_gene_merge`;
	`rm $target_sequence_gene_merge.unsorted`;
	
	#`sort-bed $target_sequence_gene_merge > $target_sequence_gene_merge.final`;
	#`mv $target_sequence_gene_merge.final $target_sequence_gene_merge`;
	
	return $target_sequence_gene_merge;
}

sub lv_capture_region {
	my ($lv_id, $dir, $bf, $option, $target_lv, $panel_capture) = @_;
	my $capture_lv;
	
	die "ERROR: Missing target file for $lv_id\n" if !-e $target_lv;
	$capture_lv = "$dir/$bf/lv_regions/$option"."_$lv_id"."_capture.bed";
	`bedops --intersect $target_lv $panel_capture > $capture_lv`;
	return $capture_lv;
}

sub dir_explorer{
	 my ($dir, $bf, $prefix) = @_;
	 my (%HoSoL, @samples, @sample_files, @lane_samples);
	 my ($lane_actual, $lane_ant, $read_ant, $read_actual);
	 my $work_dir = "$dir/$bf/rawdata";
	 $lane_actual = $lane_ant = $read_ant = $read_actual = "";	

	 
	 if(is_folder_empty($work_dir)){
		die "ERROR:$work_dir is empty !\n";	
	 }
	 
	 opendir (DIR, $work_dir) or die $!;
	 
	 while (my $file = readdir(DIR)) {	 	
	 	#Skip . and .. directories
	 	if ($file eq "." or $file eq ".."){
	 		next;
	 	}
	 	
	 	if ($nocheck == 0 and $file =~ /^(?!Sample_$prefix)/i){
	 		next;
	 	}
	 	# A file test to check that it is a directory
	 	# Use -f to test for a file
	 	next unless (-d "$work_dir/$file");
	 	#print "$file\n";
#	 	if ($file =~ /^Sample_EX71/){
#	 		print "yes\n";
#	 	}
	 	push (@samples, $file);
	 }
	 
	closedir(DIR);
	
	die "ERROR: Analysis configuration file prefix and real file prefix doesn't match.\nCheck rawdata file names or specify a file prefix!\n" if (scalar @samples == 0);
	
	foreach my $sample (@samples){
		my $sample_name = $sample;
		$sample_name =~ s/Sample_//g;
		$HoSoL{$sample_name} = ();
		#Foreach sample name, its directory will be open and the files in it listed, skipping .csv
		opendir (DIR_SAMPLE, "$work_dir/$sample");
		
		if(is_folder_empty("$work_dir/$sample")){
			die "ERROR:$work_dir/$sample is empty !\n";
		}
		
		#Obtain the file list
		while (my $file = readdir(DIR_SAMPLE)) {
			#Skip directories
			next unless (!-d "$file");
			
			my $ext = extension($file);
			
			#Just keep .fastq.gz files
			next unless ($ext eq "fastq.gz");
			
			#@sample_files contains the file list for this sample
			push (@sample_files, "$work_dir/$sample/$file");
			
			#Sort file sample order
			@sample_files = sort @sample_files;
		}
		
		#Check if @sample_files is empty
		if (!@sample_files){
			die "ERROR: No .fastq.gz files found in $work_dir/$sample\n";
		}
		
		#For the file list of the current sample, split it between lanes.
		foreach my $file (@sample_files){
			
			if ($lane_actual ne "") {
				$lane_ant = $lane_actual;
		   	}
		   	
		   	if ($read_actual ne ""){
		   		$read_ant = $read_actual;
		   	}
		   	
		   	#Regex for lane number
			$file =~ /(L\d\d\d)/;
			$lane_actual = $1;
			
			#Regex for read number
			$file =~ /_(R\d)_/;
			$read_actual = $1;
			
			#Detect lane change
			if ($lane_ant ne "" && $lane_actual ne $lane_ant){
				push @{ $HoSoL{$sample_name}{$lane_ant}{$read_ant} }, @lane_samples;
				@lane_samples = ();
			}
			
			#Detect read change
			if ($read_ant ne "" && $read_actual ne $read_ant){
				push @{ $HoSoL{$sample_name}{$lane_ant}{$read_ant} }, @lane_samples;
				@lane_samples = ();
			}
		
		#Save lane files
		push (@lane_samples, $file);
		}
		
	#Save sample, lane, files
	push @{ $HoSoL{$sample_name}{$lane_actual}{$read_actual} }, @lane_samples;
	
	#Reset vars
	@lane_samples = @sample_files = ();
	$lane_actual = $lane_ant = "";
	$read_actual = $read_ant = "";
	
   }
   return %HoSoL;
}

sub extension {
	#Returns a string with file extension. ie: "fastq.gz"
	my $path = shift;
	my $ext = (fileparse($path,'\..*'))[2];
	$ext =~ s/^\.//;
	return $ext;
}

sub sample_list {
	#Returns a string with sample names separated by commas
	my %all_samples = @_;
	my (@samples, $samples);
	for my $sample (keys %all_samples) {
		push (@samples, $sample);
	}
	# TODO: do not output "samples" prefix
	$samples = join(",",@samples);
	
	return $samples;
}

sub choose_aligner {
	my ($user_aligner, $aligner) = @_;

	if ($user_aligner ne "" && $user_aligner ne $aligner){
		return $user_aligner;
	}
	return $aligner;
}

sub header {
	my ($json, $dir, $bf) = @_;
	$$json{"name"} = "$bf";
	$$json{"path"} = "$dir/$bf";
	$$json{"references"} = {};
	$$json{"references"}{"genomes"} = {};
	$$json{"references"}{"genomes"}{"human"} = {};
	$$json{"references"}{"genomes"}{"human"}{"fasta"} = "/share/references/genomes/human/hg19/reference/human_hg19.fa";
	$$json{"references"}{"genomes"}{"human"}{"bwa_fasta"} = "/share/references/genomes/human/hg19/reference/bwa/hg19.fa";
	$$json{"references"}{"genomes"}{"human"}{"bowtie2_reference"} = "/share/references/genomes/human/hg19/reference/bowtie2/hg19";
	$$json{"references"}{"genomes"}{"human"}{"novoalign_index"} = "/share/references/genomes/human/hg19/reference/novoindex/hg19.ndx";
	$$json{"references"}{"genomes"}{"human"}{"vcf"} = "/share/references/realign_recalibrate_hapmap/common_all.vcf";
}

sub essay {
	my ($json, $dir, $bf, $analysis_option, $pipeline, $aligner, $samples, $chr_split, $target, $target_chrs, $capture, $three_prime_primers, $five_prime_primers, $haloplex_primer1, $haloplex_primer2, $haloplex_primer1_rc, $haloplex_primer2_rc, $trimming_hash, $lv_id) = @_;
	my $essay_name = "$analysis_option";
	if ($lv_id ne ""){
		$essay_name .= "_$lv_id";
	}
	$$json{"essays"} = {};
	$$json{"essays"}{$essay_name} = {};
	$$json{"essays"}{$essay_name}{"pipeline"} = $pipeline;
	$$json{"essays"}{$essay_name}{"target_reference"} = {};
	$$json{"essays"}{$essay_name}{"target_reference"}{"name"} = $analysis_option;
	$$json{"essays"}{$essay_name}{"target_reference"}{"chr_split"} = $chr_split;
	$$json{"essays"}{$essay_name}{"target_reference"}{"target"} = $target;
	$$json{"essays"}{$essay_name}{"target_reference"}{"target_chrs"} = $target_chrs;
	$$json{"essays"}{$essay_name}{"target_reference"}{"capture"} = $capture;
	
	if ($three_prime_primers ne "" && $five_prime_primers ne ""){
		$$json{"essays"}{$essay_name}{"target_reference"}{"five_prime_primers"} = $five_prime_primers;
		$$json{"essays"}{$essay_name}{"target_reference"}{"three_prime_primers"} = $three_prime_primers;
	}
	
	if ($haloplex_primer1 ne "" && $haloplex_primer2 ne "" && $haloplex_primer1_rc ne "" && $haloplex_primer2_rc ne ""){
		$$json{"essays"}{$essay_name}{"target_reference"}{"haloplex_primer1"} = $haloplex_primer1;
		$$json{"essays"}{$essay_name}{"target_reference"}{"haloplex_primer2"} = $haloplex_primer2;
		$$json{"essays"}{$essay_name}{"target_reference"}{"haloplex_primer1_rc"} = $haloplex_primer1_rc;
		$$json{"essays"}{$essay_name}{"target_reference"}{"haloplex_primer2_rc"} = $haloplex_primer2_rc;
	}
	
	$$json{"essays"}{$essay_name}{"samples"} = "$samples";
	$$json{"essays"}{$essay_name}{"modules"} = {};
	if ($pipeline =~ /trimming/i){
		$$json{"essays"}{$essay_name}{"modules"}{"trimming"} = {};
	}
	$$json{"essays"}{$essay_name}{"modules"}{"mapping"} = {};
	$$json{"essays"}{$essay_name}{"modules"}{"mapping"}{"parameters"} = {};
	$$json{"essays"}{$essay_name}{"modules"}{"mapping"}{"parameters"}{"aligner"} = $aligner;
	$$json{"essays"}{$essay_name}{"modules"}{"variant_calling"} = {};
	$$json{"essays"}{$essay_name}{"modules"}{"variant_calling"}{"parameters"} = {};
	$$json{"essays"}{$essay_name}{"modules"}{"variant_calling"}{"parameters"}{"filter"} = "yes";
	$$json{"essays"}{$essay_name}{"modules"}{"annotation"} = {};
	$$json{"essays"}{$essay_name}{"modules"}{"annotation"}{"parameters"} = {};
	$$json{"essays"}{$essay_name}{"modules"}{"annotation"}{"parameters"}{"filter"} = "no";
	$$json{"essays"}{$essay_name}{"modules"}{"mapping_stats"} = {};
	$$json{"essays"}{$essay_name}{"modules"}{"mapping_stats"}{"parameters"} = {};
	$$json{"essays"}{$essay_name}{"modules"}{"mapping_stats"}{"parameters"}{"bamstats"} = {};
	$$json{"essays"}{$essay_name}{"modules"}{"mapping_stats"}{"parameters"}{"bamstats"}{"depth_thresholds"} = {};
	$$json{"essays"}{$essay_name}{"modules"}{"mapping_stats"}{"parameters"}{"bamstats"}{"depth_thresholds"}{"graph_coverages_values"} = "1,10,20";
	$$json{"essays"}{$essay_name}{"modules"}{"mapping_stats"}{"parameters"}{"bamstats"}{"depth_thresholds"}{"region_distribution"} = "0,20,30,40,50,60,70,100";
	$$json{"essays"}{$essay_name}{"modules"}{"mapping_stats"}{"parameters"}{"NC_stats"} = {};
	$$json{"essays"}{$essay_name}{"modules"}{"mapping_stats"}{"parameters"}{"NC_stats"}{"depth_thresholds"} = "10,20,30";
}

sub samples {
		my ($json, $dir, $bf, $seq_platform, $sequencer, $size, $capture_system, %all_samples) = @_;
		my ($R1, $R2);
		my $sample_counter = 0;
		my $num_samples = scalar(keys %all_samples);
		my $sample_lane_ant = "";
		$R1 = $R2 = "";
		
		$$json{"samples"} = {};
		
		for my $sample (keys %all_samples) {
			$$json{"samples"}{$sample} = {};
			$$json{"samples"}{$sample}{"specie"} = "human";
			$$json{"samples"}{$sample}{"replicates"} = {};
			$$json{"samples"}{$sample}{"replicates"}{"replicate1"} = {};
			$$json{"samples"}{$sample}{"replicates"}{"replicate1"}{"seq_platform"} = $seq_platform;
			$$json{"samples"}{$sample}{"replicates"}{"replicate1"}{"sequencer"} = $sequencer;
			$$json{"samples"}{$sample}{"replicates"}{"replicate1"}{"size"} = $size;
			$$json{"samples"}{$sample}{"replicates"}{"replicate1"}{"capture_system"} = $capture_system;
			$$json{"samples"}{$sample}{"replicates"}{"replicate1"}{"lanes"} = {};
 					
			for my $lane (keys %{ $all_samples{$sample} }){
				my $lane_num_counter = 1;
				my $num_files = scalar (@{$all_samples{$sample}{$lane}{'R1'}});

				for (my $i=0; $i<$num_files; $i++){
					$R1 = "$all_samples{$sample}{$lane}{'R1'}[$i]";
					$R2 = "$all_samples{$sample}{$lane}{'R2'}[$i]";
					if ($R1 eq "" or $R2 eq ""){
						die "ERROR: You're missing R1 or R2 fastq.gz files in $sample";
					}
					$$json{"samples"}{$sample}{"replicates"}{"replicate1"}{"lanes"}{"$lane-$lane_num_counter"} = {};
					$$json{"samples"}{$sample}{"replicates"}{"replicate1"}{"lanes"}{"$lane-$lane_num_counter"}{"type"} = "PE";
					$$json{"samples"}{$sample}{"replicates"}{"replicate1"}{"lanes"}{"$lane-$lane_num_counter"}{"data"} = {};
					$$json{"samples"}{$sample}{"replicates"}{"replicate1"}{"lanes"}{"$lane-$lane_num_counter"}{"data"}{"R1"} = $all_samples{$sample}{$lane}{'R1'}[$i];
					$$json{"samples"}{$sample}{"replicates"}{"replicate1"}{"lanes"}{"$lane-$lane_num_counter"}{"data"}{"R2"} = $all_samples{$sample}{$lane}{'R2'}[$i];
					$lane_num_counter++;
				}#loop num_files
				$sample_lane_ant = $sample;
			}#loop lane
			$sample_counter++;
		}#loop sample
}