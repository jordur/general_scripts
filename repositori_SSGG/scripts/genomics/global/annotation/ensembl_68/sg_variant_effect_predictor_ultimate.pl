#!/usr/bin/perl

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.
  

=cut

=head1 NAME

SG Variant Effect Predictor Ultimate - a script to predict the consequences of genomic variants
163352
Version 1.0

by Juan Manuel Rosa (jm.rosa@sistemasgenomicos.com)

Based on the version 2.6 from Ensembl

=cut

use strict;
use Getopt::Long;
use FileHandle;
#use FindBin qw($Bin);
#use lib $Bin;
use Data::Dumper;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use VEP_SG qw(
    parse_line
    vf_to_consequences
    validate_vf
    convert_to_vcf
    load_dumped_adaptor_cache
    dump_adaptor_cache
    get_all_consequences
    get_slice
    build_full_cache
    read_cache_info
    get_time
    debug
    @OUTPUT_COLS
    @REG_FEAT_TYPES
    %FILTER_SHORTCUTS
);

# global vars
my $VERSION = '1.0';

 
# define headers that would normally go in the extra field
# keyed on the config parameter used to turn it on
my %extra_headers = (
    protein         => ['ENSP'],
    hgvs            => ['HGVSc','HGVSp'],
    hgnc            => ['HGNC'],
    sift            => ['SIFT'],
    polyphen        => ['PolyPhen'],
    domains         => ['DOMAINS'],
    regulatory      => ['MOTIF_NAME','HIGH_INF_POS','MOTIF_SCORE_CHANGE'],
    gmaf            => ['GMAF'],
);

my %extra_descs = (
    'HGNC'         => 'HGNC gene identifier',
    'ENSP'         => 'Ensembl protein identifer',
    'HGVSc'        => 'HGVS coding sequence name',
    'HGVSp'        => 'HGVS protein sequence name',
    'SIFT'         => 'SIFT prediction',
    'PolyPhen'     => 'PolyPhen prediction',
    'DOMAINS'      => 'The source and identifer of any overlapping protein domains',
    'MOTIF_NAME'   => 'The source and identifier of a transcription factor binding profile (TFBP) aligned at this position',
    'HIGH_INF_POS' => 'A flag indicating if the variant falls in a high information position of the TFBP',
    'MOTIF_SCORE_CHANGE' => 'The difference in motif score of the reference and variant sequences for the TFBP',
    'GMAF'         => 'Minor allele and frequency of existing variation in 1000 Genomes Phase 1',
 );

# set output autoflush for progress bars
$| = 1;

# configure from command line opts
my $config = &configure(scalar @ARGV);

# run the main sub routine
&main($config);

# this is the main sub-routine - it needs the configured $config hash
sub main {
    my $config = shift;
    
    #debug("Starting...") unless defined $config->{quiet};
    
    $config->{start_time} = time();
    $config->{last_time} = time();
    
    my $tr_cache = {};
    my $rf_cache = {};
    
    # create a hash to hold slices so we don't get the same one twice
    my %slice_cache = ();
    
    my @vfs;    
    my ($vf_count, $total_vf_count);
    my $in_file_handle = $config->{in_file_handle};
    
    # initialize line number in config
    $config->{line_number} = 0;
    
    # read the file
    while(<$in_file_handle>) {
        chomp;
        
        $config->{line_number}++;
        
        # header line?
        if(/^\#/) { next;}
        
        # configure output file
        $config->{out_file_handle} ||= &get_out_file_handle($config);
        
        # get variant object
        my $vf = &parse_line($config, $_);
            
		# now get the slice
		while (!defined($vf->{slice})) {
			#print "Getting slice...\n";
			my $slice;
			my $start = $vf->{start} - 1;
			my $end = $vf->{end} + 1;	
			my $region = "\'$vf->{chr}\', ".($vf->{start} - 1).", ".($vf->{end} + 1);
			
			#print "$region\n";
			# check if we have fetched this slice already
			if(defined $slice_cache{$region}) {
				$slice = $slice_cache{$region};
			}
			
			# if not create a new one
			else {
				$slice = &get_slice($config, $vf->{chr}, $start, $end);
				
				# if failed, warn and skip this line
				if(!defined($slice)) {
					warn("WARNING: Could not fetch slice named ".$vf->{chr}." on line ".$config->{line_number}."\n") unless defined $config->{quiet};
					next;
				}    
				
				# store the hash
				$slice_cache{$region} = $slice;
			}
			$vf->{slice} = $slice;
		}
		#print Dumper ($vf);
		
		my @consequences = &vf_to_consequences($config, $vf);
		foreach my $consequence (@consequences) {
			#print Dumper ($consequence);
			#print "$consequence->{existing}\n";
		}
		#print_line($config, $_) foreach @{vf_to_consequences($config, $vf)};
		$vf_count++;
		$total_vf_count++;
		debug("Processed $vf_count variants") if $vf_count =~ /0$/ && defined($config->{verbose});
	}

    #debug("Executed ", defined($Bio::EnsEMBL::DBSQL::StatementHandle::count_queries) ? $Bio::EnsEMBL::DBSQL::StatementHandle::count_queries : 'unknown number of', " SQL statements") if defined($config->{count_queries}) && !defined($config->{quiet});
    
    #debug("Finished!") unless defined $config->{quiet};
}

# sets up configuration hash that is used throughout the script
sub configure {
    my $args = shift;
    
    my $config = {};
    
    GetOptions(
        $config,
        'help',                    # displays help message
        
        # input options,
        'config=s',                # config file name
        'input_file|i=s',          # input file name
        'format=s',                # input file format
        
        # DB options
        'species=s',               # species e.g. human, homo_sapiens
        'registry=s',              # registry file
        'host=s',                  # database host
        'port=s',                  # database port
        'user=s',                  # database user name
        'password=s',              # database password
        'db_version=i',            # Ensembl database version to use e.g. 62
        'genomes',                 # automatically sets DB params for e!Genomes
        'refseq',                  # use otherfeatures RefSeq DB instead of Ensembl
        #'no_disconnect',           # disables disconnect_when_inactive
        
        # runtime options
        'most_severe',             # only return most severe consequence
        'summary',                 # only return one line per variation with all consquence types
        'per_gene',                # only return most severe per gene
        'buffer_size=i',           # number of variations to read in before analysis
        'chunk_size=s',            # size in bases of "chunks" used in internal hash structure
        'failed=i',                # include failed variations when finding existing
        'no_whole_genome',         # disables now default whole-genome mode
        'whole_genome',            # proxy for whole genome mode - now just warns user
        'gp',                      # read coords from GP part of INFO column in VCF (probably only relevant to 1KG)
        'chr=s',                   # analyse only these chromosomes, e.g. 1-5,10,MT
        'check_ref',               # check supplied reference allele against DB
        'check_existing',          # find existing co-located variations
        'check_svs',               # find overlapping structural variations
        'check_alleles',           # only attribute co-located if alleles are the same
        'check_frequency',         # enable frequency checking
        'gmaf',                    # add global MAF of existing var
        'freq_filter=s',           # exclude or include
        'freq_freq=f',             # frequency to filter on
        'freq_gt_lt=s',            # gt or lt (greater than or less than)
        'freq_pop=s',              # population to filter on
        'allow_non_variant',       # allow non-variant VCF lines through
        'individual=s',            # give results by genotype for individuals
        'phased',                  # force VCF genotypes to be interpreted as phased
        'fork=i',                  # fork into N processes
        
        # verbosity options
        'verbose|v',               # print out a bit more info while running
        'quiet',                   # print nothing to STDOUT (unless using -o stdout)
        'no_progress',             # don't display progress bars
        
        # output options
        'everything|e',            # switch on EVERYTHING :-)
        'output_file|o=s',         # output file name
        'force_overwrite',         # force overwrite of output file if already exists
        'terms|t=s',               # consequence terms to use e.g. NCBI, SO
        'coding_only',             # only return results for consequences in coding regions
        'canonical',               # indicates if transcript is canonical
        'ccds',                    # output CCDS identifer
        'xref_refseq',             # output refseq mrna xref
        'protein',                 # add e! protein ID to extra column
        'hgnc',                    # add HGNC gene ID to extra column
        'hgvs',                    # add HGVS names to extra column
        'sift=s',                  # SIFT predictions
        'polyphen=s',              # PolyPhen predictions
        'condel=s',                # Condel predictions
        'regulatory',              # enable regulatory stuff
        'cell_type=s',             # filter cell types for regfeats
        'convert=s',               # convert input to another format (doesn't run VEP)
        'filter=s',                # run in filtering mode
        'no_intergenic',           # don't print out INTERGENIC consequences
        'gvf',                     # produce gvf output
        'vcf',                     # produce vcf output
        'original',                # produce output in input format
        'no_consequences',         # don't calculate consequences
        'lrg',                     # enable LRG-based features
        'fields=s',                # define your own output fields
        'domains',                 # output overlapping protein features
        'numbers',                 # include exon and intron numbers
        
        # cache stuff
        'cache',                   # use cache
        'write_cache',             # enables writing to the cache
        'build=s',                 # builds cache from DB from scratch; arg is either all (all top-level seqs) or a list of chrs
        'no_adaptor_cache',        # don't write adaptor cache
        'prefetch',                # prefetch exons, translation, introns, codon table etc for each transcript
        'strip',                   # strips adaptors etc from objects before caching them
        'rebuild=s',               # rebuilds cache by reading in existing then redumping - probably don't need to use this any more
        'dir=s',                   # dir where cache is found (defaults to $HOME/.vep/)
        'cache_region_size=i',     # size of region in bases for each cache file
        'no_slice_cache',          # tell API not to cache features on slice
        'standalone',              # standalone renamed offline
        'offline',                 # offline mode uses minimal set of modules installed in same dir, no DB connection
        'skip_db_check',           # don't compare DB parameters with cached
        'compress=s',              # by default we use zcat to decompress; user may want to specify gzcat or "gzip -dc"
        'custom=s' => ($config->{custom} ||= []), # specify custom tabixed bgzipped file with annotation
        'tmpdir=s',                # tmp dir used for BigWig retrieval
        'plugin=s' => ($config->{plugin} ||= []), # specify a method in a module in the plugins directory
        
        # debug
        'cluck',                   # these two need some mods to Bio::EnsEMBL::DBSQL::StatementHandle to work. Clucks callback trace and SQL
        'count_queries',           # counts SQL queries executed
        'admin',                   # allows me to build off public hosts
        'debug',                   # print out debug info
        'tabix',                   # experimental use tabix cache files
    ) or die "ERROR: Failed to parse command-line flags\n";
    
    # print usage message if requested or no args supplied
    if(defined($config->{help}) || !$args) {
        &usage;
        exit(0);
    }
    
    # dir is where the cache and plugins live
    $config->{dir} ||= join '/', ($ENV{'HOME'}, '.vep');
   
    # dir gets set to the specific cache directory later on, so take a copy to use 
    # when configuring plugins

    $config->{toplevel_dir} = $config->{dir};

    # ini file?
    my $ini_file = $config->{dir}.'/vep.ini';
    
    if(-e $ini_file) {
        read_config_from_file($config, $ini_file);
    }
    
    # config file?
    if(defined $config->{config}) {
        read_config_from_file($config, $config->{config});
    }

    # can't be both quiet and verbose
    die "ERROR: Can't be both quiet and verbose!\n" if defined($config->{quiet}) && defined($config->{verbose});
    
    # check forking
    if(defined($config->{fork})) {
        die "ERROR: Fork number must be greater than 1\n" if $config->{fork} <= 1;
        
        # check we can use MIME::Base64
        eval q{ use MIME::Base64; };
        
        if($@) {
            debug("WARNING: Unable to load MIME::Base64, forking disabled") unless defined($config->{quiet});
            delete $config->{fork};
        }
        else {
            
            # try a practice fork
            my $pid = fork;
            
            if(!defined($pid)) {
                debug("WARNING: Fork test failed, forking disabled") unless defined($config->{quiet});
                delete $config->{fork};
            }
            elsif($pid) {
                waitpid($pid, 0);
            }
            elsif($pid == 0) {
                exit(0);
            }
        }
    }
    
    # check file format
    if(defined $config->{format}) {
        die "ERROR: Unrecognised input format specified \"".$config->{format}."\"\n" unless $config->{format} =~ /^(pileup|vcf|guess|hgvs|ensembl|id|vep)$/i;
    }
    
    # check convert format
    if(defined $config->{convert}) {
        die "ERROR: Unrecognised output format for conversion specified \"".$config->{convert}."\"\n" unless $config->{convert} =~ /vcf|ensembl|pileup|hgvs/i;
    }
    
    # check if user still using --standalone
    if(defined $config->{standalone}) {
        die "ERROR: --standalone replaced by --offline\n";
    }
    
    # connection settings for Ensembl Genomes
    if($config->{genomes}) {
        $config->{host} ||= 'mysql.ebi.ac.uk';
        $config->{port} ||= 4157;
    }
    
    # connection settings for main Ensembl
    else {
        $config->{species} ||= "homo_sapiens";
        $config->{host}    ||= 'ensembldb.ensembl.org';
        $config->{port}    ||= 5306;
    }
    
    # refseq or core? Core
    $config->{core_type} = 'core';
    
    # output term
    if(defined $config->{terms}) {
        die "ERROR: Unrecognised consequence term type specified \"".$config->{terms}."\" - must be one of ensembl, so, ncbi\n" unless $config->{terms} =~ /ensembl|display|so|ncbi/i;
        if($config->{terms} =~ /ensembl|display/i) {
            $config->{terms} = 'display';
        }
        else {
            $config->{terms} = uc($config->{terms});
        }
    }
    
    # everything?
    if(defined($config->{everything})) {
        my %everything = (
            sift       => 'b',
            polyphen   => 'b',
            ccds       => 1,
            hgvs       => 1,
            hgnc       => 1,
            numbers    => 1,
            domains    => 1,
            regulatory => 1,
            canonical  => 1,
            protein    => 1,
            gmaf       => 1,
        );
        
        $config->{$_} = $everything{$_} for keys %everything;
        
        # these ones won't work with offline
        delete $config->{hgvs} if defined($config->{offline});
    }
    
    # check nsSNP tools
    foreach my $tool(grep {defined $config->{lc($_)}} qw(SIFT PolyPhen Condel)) {
        die "ERROR: Unrecognised option for $tool \"", $config->{lc($tool)}, "\" - must be one of p (prediction), s (score) or b (both)\n" unless $config->{lc($tool)} =~ /^(s|p|b)/;
        
        die "ERROR: $tool not available for this species\n" unless $config->{species} =~ /human|homo/i;
        
        die "ERROR: $tool functionality is now available as a VEP Plugin - see http://www.ensembl.org/info/docs/variation/vep/vep_script.html#plugins\n" if $tool eq 'Condel';
    }
    
    # force quiet if outputting to STDOUT
    if(defined($config->{output_file}) && $config->{output_file} =~ /stdout/i) {
        delete $config->{verbose} if defined($config->{verbose});
        $config->{quiet} = 1;
    }
    
    # individual(s) specified?
    if(defined($config->{individual})) {
        $config->{individual} = [split /\,/, $config->{individual}];
        
        # force allow_non_variant
        $config->{allow_non_variant} = 1;
    }
    
    # summarise options if verbose
    if(defined $config->{verbose}) {
        my $header =<<INTRO;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version $VERSION

By Will McLaren (wm2\@ebi.ac.uk)

Configuration options:

INTRO
        print $header;
        
        my $max_length = (sort {$a <=> $b} map {length($_)} keys %$config)[-1];
        
        foreach my $key(sort keys %$config) {
            next if ref($config->{$key}) eq 'ARRAY' && scalar @{$config->{$key}} == 0;
            print $key.(' ' x (($max_length - length($key)) + 4)).(ref($config->{$key}) eq 'ARRAY' ? join "\t", @{$config->{$key}} : $config->{$key})."\n";
        }
        
        print "\n".("-" x 20)."\n\n";
    }
    
    # check custom annotations
    for my $i(0..$#{$config->{custom}}) {
        my $custom = $config->{custom}->[$i];
        
        my ($filepath, $shortname, $format, $type, $coords) = split /\,/, $custom;
        $type ||= 'exact';
        $format ||= 'bed';
        $coords ||= 0;
        
        # check type
        die "ERROR: Type $type for custom annotation file $filepath is not allowed (must be one of \"exact\", \"overlap\")\n" unless $type =~ /exact|overlap/;
        
        # check format
        die "ERROR: Format $format for custom annotation file $filepath is not allowed (must be one of \"bed\", \"vcf\", \"gtf\", \"gff\", \"bigwig\")\n" unless $format =~ /bed|vcf|gff|gtf|bigwig/;
        
        # bigwig format
        if($format eq 'bigwig') {
            # check for bigWigToWig
            die "ERROR: bigWigToWig does not seem to be in your path - this is required to use bigwig format custom annotations\n" unless `which bigWigToWig 2>&1` =~ /bigWigToWig$/;
        }
        
        else {
            # check for tabix
            die "ERROR: tabix does not seem to be in your path - this is required to use custom annotations\n" unless `which tabix 2>&1` =~ /tabix$/;
            
            # remote files?
            if($filepath =~ /tp\:\/\//) {
                my $remote_test = `tabix $filepath 1:1-1 2>&1`;
                if($remote_test =~ /fail/) {
                    die "$remote_test\nERROR: Could not find file or index file for remote annotation file $filepath\n";
                }
                elsif($remote_test =~ /get_local_version/) {
                    debug("Downloaded tabix index file for remote annotation file $filepath") unless defined($config->{quiet});
                }
            }
        
            # check files exist
            else {
                die "ERROR: Custom annotation file $filepath not found\n" unless -e $filepath;
                die "ERROR: Tabix index file $filepath\.tbi not found - perhaps you need to create it first?\n" unless -e $filepath.'.tbi';
            }
        }
        
        $config->{custom}->[$i] = {
            'file'   => $filepath,
            'name'   => $shortname || 'CUSTOM'.($i + 1),
            'type'   => $type,
            'format' => $format,
            'coords' => $coords,
        };
    }
    
    # check if using filter and original
    die "ERROR: You must also provide output filters using --filter to use --original\n" if defined($config->{original}) && !defined($config->{filter});
    
    # filter by consequence?
    if(defined($config->{filter})) {
        
        my %filters = map {$_ => 1} split /\,/, $config->{filter};
        
        # add in shortcuts
        foreach my $filter(keys %filters) {
            my $value = 1;
            if($filter =~ /^no_/) {
                delete $filters{$filter};
                $filter =~ s/^no_//g;
                $value = 0;
                $filters{$filter} = $value;
            }
            
            if(defined($FILTER_SHORTCUTS{$filter})) {
                delete $filters{$filter};
                $filters{$_} = $value for keys %{$FILTER_SHORTCUTS{$filter}};
            }
        }
        
        $config->{filter} = \%filters;
        
        $config->{filter_count} = 0;
    }
    
    # set defaults
    $config->{user}              ||= 'anonymous';
    $config->{buffer_size}       ||= 5000;
    $config->{chunk_size}        ||= '50kb';
    $config->{output_file}       ||= "variant_effect_output.txt";
    $config->{tmpdir}            ||= '/tmp';
    $config->{format}            ||= 'guess';
    $config->{terms}             ||= 'SO';
    $config->{cache_region_size} ||= 1000000;
    $config->{compress}          ||= 'zcat';
    
    # regulatory has to be on for cell_type
    if(defined($config->{cell_type})) {
        $config->{regulatory} = 1;
        $config->{cell_type} = [split /\,/, $config->{cell_type}] if defined($config->{cell_type});
    }
    
    # can't use a whole bunch of options with most_severe
    if(defined($config->{most_severe})) {
        foreach my $flag(qw(no_intergenic protein hgnc sift polyphen coding_only ccds canonical xref_refseq numbers domains summary)) {
            die "ERROR: --most_severe is not compatible with --$flag\n" if defined($config->{$flag});
        }
    }
    
    # can't use a whole bunch of options with summary
    if(defined($config->{summary})) {
        foreach my $flag(qw(no_intergenic protein hgnc sift polyphen coding_only ccds canonical xref_refseq numbers domains most_severe)) {
            die "ERROR: --summary is not compatible with --$flag\n" if defined($config->{$flag});
        }
    }
    
    # frequency filtering
    if(defined($config->{check_frequency})) {
        foreach my $flag(qw(freq_freq freq_filter freq_pop freq_gt_lt)) {
            die "ERROR: To use --check_frequency you must also specify flag --$flag\n" unless defined $config->{$flag};
        }
        
        # need to set check_existing
        $config->{check_existing} = 1;
    }
    
    $config->{check_existing} = 1 if defined $config->{check_alleles} || defined $config->{gmaf};
    
    # warn users still using whole_genome flag
    if(defined($config->{whole_genome})) {
        debug("INFO: Whole-genome mode is now the default run-mode for the script. To disable it, use --no_whole_genome") unless defined($config->{quiet});
    }
    
    $config->{whole_genome}      = 1 unless defined $config->{no_whole_genome};
    $config->{failed}            = 0 unless defined $config->{failed};
    $config->{chunk_size}        =~ s/mb?/000000/i;
    $config->{chunk_size}        =~ s/kb?/000/i;
    $config->{cache_region_size} =~ s/mb?/000000/i;
    $config->{cache_region_size} =~ s/kb?/000/i;
    
    # cluck and display executed SQL?
    $Bio::EnsEMBL::DBSQL::StatementHandle::cluck = 1 if defined($config->{cluck});
    
    # offline needs cache, can't use HGVS
    if(defined($config->{offline})) {
        $config->{cache} = 1;
        
        #die("ERROR: Cannot generate HGVS coordinates in offline mode\n") if defined($config->{hgvs});
        die("ERROR: Cannot use HGVS as input in offline mode\n") if $config->{format} eq 'hgvs';
        die("ERROR: Cannot use variant identifiers as input in offline mode\n") if $config->{format} eq 'id';
        die("ERROR: Cannot do frequency filtering in offline mode\n") if defined($config->{check_frequency});
        die("ERROR: Cannot retrieve overlapping structural variants in offline mode\n") if defined($config->{check_sv});
    }
    
    # write_cache needs cache
    $config->{cache} = 1 if defined $config->{write_cache};
    
    # no_slice_cache, prefetch and whole_genome have to be on to use cache
    if(defined($config->{cache})) {
        $config->{prefetch} = 1;
        $config->{no_slice_cache} = 1;
        $config->{whole_genome} = 1;
        $config->{strip} = 1;
    }
    
    $config->{build} = $config->{rebuild} if defined($config->{rebuild});
    
    # force options for full build
    if(defined($config->{build})) {
        $config->{prefetch} = 1;
        $config->{hgnc} = 1;
        $config->{no_slice_cache} = 1;
        $config->{cache} = 1;
        $config->{strip} = 1;
        $config->{write_cache} = 1;
        $config->{cell_type} = 1 if defined($config->{regulatory});
    }
    
    # connect to databases
    $config->{reg} = &connect_to_dbs($config);
    
    # complete dir with species name and db_version
    $config->{dir} .= '/'.(
        join '/', (
            defined($config->{offline}) ? $config->{species} : ($config->{reg}->get_alias($config->{species}) || $config->{species}),
            $config->{db_version} || $config->{reg}->software_version
        )
    );
    
    # warn user cache directory doesn't exist
    if(!-e $config->{dir}) {
        
        # if using write_cache
        if(defined($config->{write_cache})) {
            debug("INFO: Cache directory ", $config->{dir}, " not found - it will be created") unless defined($config->{quiet});
        }
        
        # want to read cache, not found
        elsif(defined($config->{cache})) {
            die("ERROR: Cache directory ", $config->{dir}, " not found");
        }
    }
    
    if(defined($config->{cache})) {
        # read cache info
        if(read_cache_info($config)) {
            debug("Read existing cache info") unless defined $config->{quiet};
        }
    }
   
    # we configure plugins here because they can sometimes switch on the 
    # regulatory config option
    configure_plugins($config);
    
    # include regulatory modules if requested
    if(defined($config->{regulatory})) {
        # do the use statements here so that users don't have to have the
        # funcgen API installed to use the rest of the script
        eval q{
            use Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryFeatureAdaptor;
            use Bio::EnsEMBL::Funcgen::DBSQL::MotifFeatureAdaptor;
            use Bio::EnsEMBL::Funcgen::MotifFeature;
            use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
            use Bio::EnsEMBL::Funcgen::BindingMatrix;
        };
        
        if($@) {
            die("ERROR: Ensembl Funcgen API must be installed to use --regulatory or plugins that deal with regulatory features\n");
        }
    }
    
    # user defined custom output fields
    if(defined($config->{fields})) {
        $config->{fields} = [split ',', $config->{fields}];
        debug("Output fields redefined (".scalar @{$config->{fields}}." defined)") unless defined($config->{quiet});
        $config->{fields_redefined} = 1;
    }
    $config->{fields} ||= \@OUTPUT_COLS;
    
    # suppress warnings that the FeatureAdpators spit if using no_slice_cache
    Bio::EnsEMBL::Utils::Exception::verbose(1999) if defined($config->{no_slice_cache});
    
    # get adaptors (don't get them in offline mode)
    unless(defined($config->{offline})) {
        
        if(defined($config->{cache}) && !defined($config->{write_cache})) {
            
            # try and load adaptors from cache
            if(!&load_dumped_adaptor_cache($config)) {
                &get_adaptors($config);
                &dump_adaptor_cache($config) if defined($config->{write_cache}) && !defined($config->{no_adaptor_cache});
            }
            
            # check cached adaptors match DB params
            else {
                my $dbc = $config->{sa}->{dbc};
            
                my $ok = 1;
                
                if($dbc->{_host} ne $config->{host}) {
                    
                    # ens-livemirror, useastdb and ensembldb should all have identical DBs
                    unless(
                        (
                            $dbc->{_host} eq 'ens-livemirror'
                            || $dbc->{_host} eq 'ensembldb.ensembl.org'
                            || $dbc->{_host} eq 'useastdb.ensembl.org'
                        ) && (
                            $config->{host} eq 'ens-livemirror'
                            || $config->{host} eq 'ensembldb.ensembl.org'
                            || $config->{host} eq 'useastdb.ensembl.org'
                        )
                    ) {
                        $ok = 0;
                    }
                    
                    unless(defined($config->{skip_db_check})) {
                        # but we still need to reconnect
                        debug("INFO: Defined host ", $config->{host}, " is different from cached ", $dbc->{_host}, " - reconnecting to host") unless defined($config->{quiet});
                        
                        &get_adaptors($config);
                    }
                }
                
                if(!$ok) {
                    if(defined($config->{skip_db_check})) {
                        debug("INFO: Defined host ", $config->{host}, " is different from cached ", $dbc->{_host}) unless defined($config->{quiet});
                    }
                    else {
                        die "ERROR: Defined host ", $config->{host}, " is different from cached ", $dbc->{_host}, ". If you are sure this is OK, rerun with -skip_db_check flag set";
                    }
                }
            }
        }
        else {
            &get_adaptors($config);
            &dump_adaptor_cache($config) if defined($config->{write_cache}) && !defined($config->{no_adaptor_cache});
        }
        
        # reg adaptors (only fetches if not retrieved from cache already)
        &get_reg_adaptors($config) if defined($config->{regulatory});
    }
    
    # check cell types
    if(defined($config->{cell_type}) && !defined($config->{build})) {
        my $cls = '';
        
        if(defined($config->{cache})) {
            $cls = $config->{cache_cell_types};
        }
        else {
            my $cta = $config->{RegulatoryFeature_adaptor}->db->get_CellTypeAdaptor();
            $cls = join ",", map {$_->name} @{$cta->fetch_all};
        }
        
        foreach my $cl(@{$config->{cell_type}}) {
            die "ERROR: cell type $cl not recognised; available cell types are:\n$cls\n" unless $cls =~ /(^|,)$cl(,|$)/;
        }
    }
    
    # get terminal width for progress bars
    unless(defined($config->{quiet})) {
        my $width;
        
        # module may not be installed
        eval q{
            use Term::ReadKey;
        };
        
        if(!$@) {
            my ($w, $h);
            
            # module may be installed, but e.g.
            eval {
                ($w, $h) = GetTerminalSize();
            };
            
            $width = $w if defined $w;
        }
        
        $width ||= 60;
        $width -= 12;
        $config->{terminal_width} = $width;
    }
    
    # jump out to build cache if requested
    if(defined($config->{build})) {
        
        if($config->{host} =~ /^(ensembl|useast)db\.ensembl\.org$/ && !defined($config->{admin})) {
            die("ERROR: Cannot build cache using public database server ", $config->{host}, "\n");
        }
        
        # build the cache
        debug("Building cache for ".$config->{species}) unless defined($config->{quiet});
        build_full_cache($config);
        
        # exit script
        debug("Finished building cache") unless defined($config->{quiet});
        exit(0);
    }
    
    
    # warn user DB will be used for SIFT/PolyPhen/HGVS/frequency/LRG
    if(defined($config->{cache})) {
        
        # these two def depend on DB
        foreach my $param(grep {defined $config->{$_}} qw(hgvs check_frequency lrg check_sv)) {
            debug("INFO: Database will be accessed when using --$param") unless defined($config->{quiet});
        }
        
        # as does using HGVS or IDs as input
        debug("INFO: Database will be accessed when using --format ", $config->{format}) if ($config->{format} eq 'id' || $config->{format} eq 'hgvs') && !defined($config->{quiet});
        
        # the rest may be in the cache
        foreach my $param(grep {defined $config->{$_}} qw(sift polyphen regulatory)) {
            next if defined($config->{'cache_'.$param});
            debug("INFO: Database will be accessed when using --$param; consider using the complete cache containing $param data (see documentation for details)") unless defined($config->{quiet});
        }
    }
    
    # get list of chrs if supplied
    if(defined($config->{chr})) {
        my %chrs;
        
        foreach my $val(split /\,/, $config->{chr}) {
            my @nnn = split /\-/, $val;
            
            foreach my $chr($nnn[0]..$nnn[-1]) {
                $chrs{$chr} = 1;
            }
        }
        
        $config->{chr} = \%chrs;
    }
    
    # get input file handle
    $config->{in_file_handle} = &get_in_file_handle($config);
    
    return $config;
}

# reads config from a file
sub read_config_from_file {
    my $config = shift;
    my $file = shift;
    
    open CONFIG, $file or die "ERROR: Could not open config file \"$file\"\n";
    
    while(<CONFIG>) {
        next if /^\#/;
        my @split = split /\s+|\=/;
        my $key = shift @split;
        $key =~ s/^\-//g;
        
        if(defined($config->{$key}) && ref($config->{$key}) eq 'ARRAY') {
            push @{$config->{$key}}, @split;
        }
        else {
            $config->{$key} ||= $split[0];
        }
    }
    
    close CONFIG;
    
    # force quiet if outputting to STDOUT
    if(defined($config->{output_file}) && $config->{output_file} =~ /stdout/i) {
        delete $config->{verbose} if defined($config->{verbose});
        $config->{quiet} = 1;
    }
    
    debug("Read configuration from $file") unless defined($config->{quiet});
}

# configures custom VEP plugins
sub configure_plugins {

    my $config = shift;
    
    $config->{plugins} = [];
    
    if (my @plugins = @{ $config->{plugin} }) {

        # add the Plugins directory onto @INC

        unshift @INC, $config->{toplevel_dir}."/Plugins";

        for my $plugin (@plugins) {

            # parse out the module name and parameters

            my ($module, @params) = split /,/, $plugin;

            # check we can use the module
            
            eval qq{
                use $module;
            };
            if ($@) {
                debug("Failed to compile plugin $module: $@") unless defined($config->{quiet});
                next;
            }
            
            # now check we can instantiate it, passing any parameters to the constructor
            
            my $instance;
            
            eval {
                $instance = $module->new($config, @params);
            };
            if ($@) {
                debug("Failed to instantiate plugin $module: $@") unless defined($config->{quiet});
                next;
            }

            # check that the versions match
            
            my $plugin_version;
            
            if ($instance->can('version')) {
                $plugin_version = $instance->version;
            }
            
            my $version_ok = 1;

            if ($plugin_version) {
                my ($plugin_major, $plugin_minor, $plugin_maintenance) = split /\./, $plugin_version;
                my ($major, $minor, $maintenance) = split /\./, $VERSION;
    
                if ($plugin_major != $major) {
                    debug("Warning: plugin $plugin version ($plugin_version) does not match the current VEP version ($VERSION)") unless defined($config->{quiet});
                    $version_ok = 0;
                }
            }
            else {
                debug("Warning: plugin $plugin does not define a version number") unless defined($config->{quiet});
                $version_ok = 0;
            }

            debug("You may experience unexpected behaviour with this plugin") unless defined($config->{quiet}) || $version_ok;

            # check that it implements all necessary methods
            
            for my $required(qw(run get_header_info check_feature_type check_variant_feature_type)) {
                unless ($instance->can($required)) {
                    debug("Plugin $module doesn't implement a required method '$required', does it inherit from BaseVepPlugin?") unless defined($config->{quiet});
                    next;
                }
            }
           
            # all's good, so save the instance in our list of plugins
            
            push @{ $config->{plugins} }, $instance;
            
            debug("Loaded plugin: $module") unless defined($config->{quiet}); 

            # for convenience, check if the plugin wants regulatory stuff and turn on the config option if so
            
            if (grep { $_ =~ /motif|regulatory/i } @{ $instance->feature_types }) {
                debug("Fetching regulatory features for plugin: $module") unless defined($config->{quiet});
                $config->{regulatory} = 1;
            }
        }
    }
} 

# connects to DBs (not done in offline mode)
sub connect_to_dbs {
    my $config = shift;
    
    # get registry
    my $reg = 'Bio::EnsEMBL::Registry';
    
    unless(defined($config->{offline})) {
        # load DB options from registry file if given
        if(defined($config->{registry})) {
            debug("Loading DB config from registry file ", $config->{registry}) unless defined($config->{quiet});
            $reg->load_all(
                $config->{registry},
                $config->{verbose},
                undef,
                $config->{no_slice_cache}
            );
        }
        
        # otherwise manually connect to DB server
        else {
            $reg->load_registry_from_db(
                -host       => $config->{host},
                -user       => $config->{user},
                -pass       => $config->{password},
                -port       => $config->{port},
                -db_version => $config->{db_version},
                -species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
                -verbose    => $config->{verbose},
                -no_cache   => $config->{no_slice_cache},
            );
        }
        
        eval { $reg->set_reconnect_when_lost() };
        
        if(defined($config->{verbose})) {
            # get a meta container adaptors to check version
            my $core_mca = $reg->get_adaptor($config->{species}, 'core', 'metacontainer');
            my $var_mca = $reg->get_adaptor($config->{species}, 'variation', 'metacontainer');
            
            if($core_mca && $var_mca) {
                debug(
                    "Connected to core version ", $core_mca->get_schema_version, " database ",
                    "and variation version ", $var_mca->get_schema_version, " database"
                );
            }
        }
    }
    
    return $reg;
}

# get adaptors from DB
sub get_adaptors {
    my $config = shift;
    
    die "ERROR: No registry" unless defined $config->{reg};
    
    $config->{vfa}   = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variationfeature');
    $config->{svfa}  = $config->{reg}->get_adaptor($config->{species}, 'variation', 'structuralvariationfeature');
    $config->{tva}   = $config->{reg}->get_adaptor($config->{species}, 'variation', 'transcriptvariation');
    $config->{pfpma} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'proteinfunctionpredictionmatrix');
    $config->{va}    = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variation');
    
    # get fake ones for species with no var DB
    if(!defined($config->{vfa})) {
        $config->{vfa}  = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($config->{species});
        $config->{svfa} = Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor->new_fake($config->{species});
        $config->{tva}  = Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor->new_fake($config->{species});
    }
    
    $config->{sa}  = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'slice');
    $config->{ga}  = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'gene');
    $config->{ta}  = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'transcript');
    $config->{mca} = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'metacontainer');
    $config->{csa} = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'coordsystem');
    
    # cache schema version
    $config->{mca}->get_schema_version if defined $config->{mca};
    
    # check we got slice adaptor - can't continue without a core DB
    die("ERROR: Could not connect to core database\n") unless defined $config->{sa};
}

# gets regulatory adaptors
sub get_reg_adaptors {
    my $config = shift;

    foreach my $type(@REG_FEAT_TYPES) {
        next if defined($config->{$type.'_adaptor'});
        
        my $adaptor = $config->{reg}->get_adaptor($config->{species}, 'funcgen', $type);
        if(defined($adaptor)) {
            $config->{$type.'_adaptor'} = $adaptor;
        }
        else {
            delete $config->{regulatory};
            last;
        }
    }
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

# gets file handle for output and adds header
sub get_out_file_handle {
    my $config = shift;
    
    # define filehandle to write to
    my $out_file_handle = new FileHandle;
    
    # check if file exists
    if(-e $config->{output_file} && !defined($config->{force_overwrite})) {
        die("ERROR: Output file ", $config->{output_file}, " already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite\n");
    }
    
    if($config->{output_file} =~ /stdout/i) {
        $out_file_handle = *STDOUT;
    }
    else {
        $out_file_handle->open(">".$config->{output_file}) or die("ERROR: Could not write to output file ", $config->{output_file}, "\n");
    }
    
    # define headers for a VCF file
    my @vcf_headers = (
        '#CHROM',
        'POS',
        'ID',
        'REF',
        'ALT',
        'QUAL',
        'FILTER',
        'INFO'
    );
    
    # file conversion, don't want to add normal headers
    if(defined($config->{convert})) {
        # header for VCF
        if($config->{convert} =~ /vcf/i) {
            print $out_file_handle "##fileformat=VCFv4.0\n";
            print $out_file_handle join "\t", @vcf_headers;
            print $out_file_handle "\n";
        }
        
        return $out_file_handle;
    }
    
    # GVF output, no header
    elsif(defined($config->{gvf}) || defined($config->{original})) {
        print $out_file_handle join "\n", @{$config->{headers}} if defined($config->{headers}) && defined($config->{original});
        return $out_file_handle;
    }
    
    elsif(defined($config->{vcf})) {
        
        # create an info string for the VCF header        
        my @new_headers;
        
        # if the user has defined the fields themselves, we don't need to worry
        if(defined $config->{fields_redefined}) {
            @new_headers = @{$config->{fields}};
        }
        else {
            @new_headers = (
                
                # get default headers, minus variation name and location (already encoded in VCF)
                grep {
                    $_ ne 'Uploaded_variation' and
                    $_ ne 'Location' and
                    $_ ne 'Extra'
                } @{$config->{fields}},
                
                # get extra headers
                map {@{$extra_headers{$_}}}
                grep {defined $config->{$_}}
                keys %extra_headers
            );
            
            # plugin headers
            foreach my $plugin_header(split /\n/, get_plugin_headers($config)) {
                $plugin_header =~ /\#\# (.+?)\t\:.+/;
                push @new_headers, $1;
            }
            
            # redefine the main headers list in config
            $config->{fields} = \@new_headers;
        }
        
        # add the newly defined headers as a header to the VCF
        my $string = join '|', @{$config->{fields}};
        my @vcf_info_strings = ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: '.$string.'">');
        
        # add custom headers
        foreach my $custom(@{$config->{custom}}) {
            push @vcf_info_strings, '##INFO=<ID='.$custom->{name}.',Number=.,Type=String,Description="'.$custom->{file}.' ('.$custom->{type}.')">';
        }
        
        # if this is already a VCF file, we need to add our new headers in the right place
        if(defined($config->{headers})) {
            
            for my $i(0..$#{$config->{headers}}) {
                if($config->{headers}->[$i] =~ /^\#CHROM\s+POS\s+ID/) {
                    splice(@{$config->{headers}}, $i, 0, @vcf_info_strings);
                }
            }
            
            print $out_file_handle join "\n", @{$config->{headers}};
            print $out_file_handle "\n";
        }
        
        else {
            print $out_file_handle "##fileformat=VCFv4.0\n";
            print $out_file_handle join "\n", @vcf_info_strings;
            print $out_file_handle "\n";
            print $out_file_handle join "\t", @vcf_headers;
            print $out_file_handle "\n";
        }
        
        return $out_file_handle;
    }
    
    # make header
    my $time = &get_time;
    my $db_string = $config->{mca}->dbc->dbname." on ".$config->{mca}->dbc->host if defined $config->{mca};
    $db_string .= "\n## Using cache in ".$config->{dir} if defined($config->{cache});
    my $version_string =
        "Using API version ".$config->{reg}->software_version.
        ", DB version ".(defined $config->{mca} && $config->{mca}->get_schema_version ? $config->{mca}->get_schema_version : '?');
    
    # add key for extra column headers based on config
    my $extra_column_keys = join "\n",
        map {'## '.$_.' : '.$extra_descs{$_}}
        sort map {@{$extra_headers{$_}}}
        grep {defined $config->{$_}}
        keys %extra_headers;
    
    my $header =<<HEAD;
## ENSEMBL VARIANT EFFECT PREDICTOR v$VERSION
## Output produced at $time
## Connected to $db_string
## $version_string
## Extra column keys:
$extra_column_keys
HEAD
   
    $header .= get_plugin_headers($config);
    
    # add headers
    print $out_file_handle $header;
    
    # add custom data defs
    if(defined($config->{custom})) {
        foreach my $custom(@{$config->{custom}}) {
            print $out_file_handle '## '.$custom->{name}."\t: ".$custom->{file}.' ('.$custom->{type}.")\n";
        }
    }
    
    # add column headers
    print $out_file_handle '#', (join "\t", @{$config->{fields}});
    print $out_file_handle "\n";
    
    return $out_file_handle;
}

sub get_plugin_headers {

    my $config = shift;

    my $header = "";

    for my $plugin (@{ $config->{plugins} }) {
        if (my $hdr = $plugin->get_header_info) {
            for my $key (keys %$hdr) {
                my $val = $hdr->{$key};
                
                $header .= "## $key\t: $val\n";
            }
        }
    }

    return $header;
}

# convert a variation feature to a line of output
sub convert_vf {
    my $config = shift;
    my $vf = shift;
    
    my $convert_method = 'convert_to_'.lc($config->{convert});
    my $method_ref   = \&$convert_method; 
    
    my $line = &$method_ref($config, $vf);
    my $handle = $config->{out_file_handle};
    
    if(scalar @$line) {
        print $handle join "\t", @$line;
        print $handle "\n";
    }
}

# converts to Ensembl format
sub convert_to_ensembl {
    my $config = shift;
    my $vf = shift;
    
    return [
        $vf->{chr} || $vf->seq_region_name,
        $vf->start,
        $vf->end,
        $vf->allele_string,
        $vf->strand,
        $vf->variation_name
    ];
}


# converts to HGVS (hackily returns many lines)
sub convert_to_hgvs {
    my $config = shift;
    my $vf = shift;
    
    # ensure we have a slice
    $vf->{slice} ||= get_slice($config, $vf->{chr});
    
    my $tvs = $vf->get_all_TranscriptVariations;
    
    my @return = values %{$vf->get_all_hgvs_notations()};
    
    if(defined($tvs)) {
        push @return, map {values %{$vf->get_all_hgvs_notations($_->transcript, 'c')}} @$tvs;
        push @return, map {values %{$vf->get_all_hgvs_notations($_->transcript, 'p')}} @$tvs;
    }
    
    return [join "\n", @return];
}

# prints a line of output from the hash
sub print_line {
    my $config = shift;
    my $line = shift;
    return unless defined($line);
    
    my $output;
    
    # normal
    if(ref($line) eq 'HASH') {
        my %extra = %{$line->{Extra}};
        
        $line->{Extra} = join ';', map { $_.'='.$line->{Extra}->{$_} } keys %{ $line->{Extra} || {} };
        
        # if the fields have been redefined we need to search through in case
        # any of the defined fields are actually part of the Extra hash
        $output = join "\t", map {
            (defined $line->{$_} ? $line->{$_} : (defined $extra{$_} ? $extra{$_} : '-'))
        } @{$config->{fields}};
    }
    
    # gvf/vcf
    else {
        $output = $$line;
    }
    
    my $fh = $config->{out_file_handle};
    print $fh "$output\n";
}

# outputs usage message
sub usage {
    my $usage =<<END;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version $VERSION

By Juan Manuel Rosa (jm.rosa\@sistemasgenomicos.com)


Usage:
perl sg_variant_effect_predictor_ultimate.pl [arguments]

Options
=======

--help                 Display this message and quit
--verbose              Display verbose output as the script runs [default: off]
--quiet                Suppress status and warning messages [default: off]
--no_progress          Suppress progress bars [default: off]
                    
--fork [num_forks]     Use forking to improve script runtime [default: off]

-i | --input_file      Input file - if not specified, reads from STDIN. Files
                       may be gzip compressed.
-o | --output_file     Output file. Write to STDOUT by specifying -o STDOUT - this
                       will force --quiet [default: "variant_effect_output.txt"]
--force_overwrite      Force overwriting of output file [default: quit if file
                       exists]
                       
--species [species]    Species to use [default: "human"]

 
--sift=[p|s|b]         Add SIFT [p]rediction, [s]core or [b]oth [default: off]
--polyphen=[p|s|b]     Add PolyPhen [p]rediction, [s]core or [b]oth [default: off]

NB: SIFT and PolyPhen predictions are currently available for human only
NB: Condel support has been moved to a VEP plugin module - see documentation

--regulatory           Look for overlaps with regulatory regions. The script can
                       also call if a variant falls in a high information position
                       within a transcription factor binding site. Output lines have
                       a Feature type of RegulatoryFeature or MotifFeature
                       [default: off]
                       
NB: Regulatory consequences are currently available for human and mouse only

--hgnc                 Add HGNC gene identifiers to output [default: off]
--hgvs                 Output HGVS identifiers (coding and protein). Requires database
                       connection [default: off]
--protein              Output Ensembl protein identifer [default: off]
--domains              Include details of any overlapping protein domains [default: off]

--check_existing       If specified, checks for existing co-located variations in the
                       Ensembl Variation database [default: off]

--gmaf                 Include global MAF of existing variant from 1000 Genomes
                       Phase 1 in output

					   
--chr [list]           Select a subset of chromosomes to analyse from your file. Any
                       data not on this chromosome in the input will be skipped. The
                       list can be comma separated, with "-" characters representing
                       a range e.g. 1-5,8,15,X [default: off]

                       
--host                 Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user            Database username [default: "anonymous"]
--port                 Database port [default: 5306]
--password             Database password [default: no password]
--genomes              Sets DB connection params for Ensembl Genomes [default: off]
--db_version=[number]  Force script to load DBs from a specific Ensembl version. Not
                       advised due to likely incompatibilities between API and DB. Default "68"
END

    print $usage;
}
