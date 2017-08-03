#!/usr/bin/perl

#=================================================#
# Getting coverage from global file to probe file #
#=================================================#


use warnings;
use strict;
use Getopt::Long;
use Parallel::ForkManager;
use Cwd;
use File::Copy;
use Statistics::Descriptive qw(:all);
use Data::Dumper;
use Tie::File;

print "Process started...".localtime()."\n";
my $threads = 8;

my ($file) = @ARGV;

if (not defined $file) { die "split_coverage.pl <file>\n";}

my $pm = new Parallel::ForkManager($threads); # Run in parallel with 8 threads max
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; } );
$pm->run_on_start( sub { my ($pid,$ident)=@_; } );
$pm->run_on_wait( sub {}, 0.5 );

foreach my $chr (1..21) {
	my $pid = $pm->start($chr) and next;
	my $command = '\'{if ($1 == "'.$chr.'" || $1 == "chr'.$chr.'") {print $0} else if ($1 == "'.($chr + 1).'" || $1 == "chr'.($chr + 1).'") {exit}}\'';
	print "awk $command $file > $file.$chr\n";
	`awk $command $file > $file.$chr`;
	$pm->finish;
}
$pm->wait_all_children;

foreach my $chr ("22","X", "Y") {
	my $pid = $pm->start($chr) and next;
	if ($chr =~ m/\d+/) {
		my $command = '\'{if ($1 == "'.$chr.'" || $1 == "chr'.$chr.'") {print $0} else if ($1 == "X" || $1 == "chrX") {exit}}\'';		
		print "awk $command $file > $file.$chr\n";
		`awk $command $file > $file.$chr`;
	}
	elsif ($chr eq "X") {
		my $command = '\'{if ($1 == "'.$chr.'" || $1 == "chr'.$chr.'") {print $0} else if ($1 == "Y" || $1 == "chrY") {exit}}\'';		
		print "awk $command $file > $file.$chr\n";
		`awk $command $file > $file.$chr`;
	}
	else {
		my $command = '\'$1 == "'.$chr.'" || $1 == "chr'.$chr.'"\'';
		print "awk $command $file > $file.$chr\n";
		`awk $command $file > $file.$chr`;
	}
	$pm->finish;
}
$pm->wait_all_children;

print "Process finished...".localtime()."\n";

exit
