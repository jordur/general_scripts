#!/usr/bin/perl -w
##################################################################
##
## Written by Jingwei Ni.
## Copyright 2011 Life Technologies Corporation. All rights reserved. 
## First created on: Wed 2011-02-09  Time: 07:44:34
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## For Research Use Only.  Not intended for any animal or human
## therapeutic or diagnostic use.
##
## The trademarks mentioned herein are the property of Life
## Technologies Corporation or their respective owners.
##
##################################################################
# 0.1 - First quick solution
# 0.2 - Added debugging and stats
# 0.25 - Full pairing stats, added unmapped bam file as input, made 
# cleaned up bam file output as an option
# 0.26 - Bug fixing per Hank's findings

print STDERR "\nCalculate pairing stats (and clean up junk entries) from 1.3.1 Bam files.\n\n";

use strict;
use Getopt::Long;
use vars qw($bamFile $memory $remove $clean $help $debug $version $stats_file $bam_unmapped);
use FileHandle;

my $VERSION = 0.26;
$debug = 0;
$clean = 1;
$remove = 1;

@ARGV == 0 and usage();
GetOptions ("mem|m=i"          => \$memory,      # numeric
            "bam|b=s"          => \$bamFile,      # string
            "u|umapped_bam=s"  => \$bam_unmapped,
            "clean_bam|clean!" => \$clean,
            "remove|r!"        => \$remove,
            "debug|d=i"        => \$debug,
            "version|v=s"      => \$version,
            "stats|s=s"        => \$stats_file,
            "help|h|?"         => \$help
    );       # flag

(!$bamFile or !$bam_unmapped or $help) and usage();

sub usage {
    print  <<"ZZ_USAGE_ZZ";
Usage:\n\t$0 [-help] -b mapped_bam_file -u unmapped_bam_file [-mem num] -s stats_file

 -help               : Print usage
 -bam|b file         : Mapped (Bopscope 1.3) Bam file, required
 -umapped_bam|u file : Umapped (Bopscope 1.3) Bam file, required
 -mem num            : Max memory size in bytes used for Bam file sorting by samtools, 
                       default 500000000 (500M) when not specified
 -clean_bam|clean    : Flag to output _Clean.bam file (after removing junk entries)
                       Default on, use -noclean_bam to turn it off
 -r                  : Flag to remove intermediate bam file sorted by bead_id
                       Default on, use -nor to keep the bead_id sorted bam file
 -stats   file       : File name for stats output file, print to STDOUT if not specified

ZZ_USAGE_ZZ
    exit;
}

print $VERSION, "\n" if $version;

if ( $memory ) { # prevents samtools sort from over-shooting mem usage
    $memory = 0.6 * $memory;
}

my $prefix = $bamFile;
$prefix =~ s/\.bam$//;
print STDERR "Sorting $bamFile by reads IDs to $prefix\_sortByName.bam\n";
my $cmd_sort = "samtools sort -n $bamFile $prefix\_sortByName\n";
$cmd_sort = "samtools sort -n -m $memory $bamFile $prefix\_sortByName\n" if $memory;
system($cmd_sort) and die "Failed to execute command $cmd_sort: $!\n";

my $cmd_bam_read = "samtools view -X $prefix\_sortByName.bam";
my $fh_input = new FileHandle " $cmd_bam_read | " or die "Failed to open input pipe: $cmd_bam_read | :$!\n";

my $fh_output;
if ( $clean ) {
    my $bamHeader = `samtools view -H $bamFile`;
    my $cmd_bam_sort = "samtools view -bSu - | samtools sort ";
    $cmd_bam_sort .= " -m $memory" if $memory;
    $cmd_bam_sort .= " - $prefix\_Clean";
    $fh_output = new FileHandle " | $cmd_bam_sort" or die "Failed to open output pipe | $cmd_bam_sort: $!\n";
    print STDERR "Cleaning up bam records and resorting by coordinates, output file to $prefix\_Clean.bam\n";
    print $fh_output $bamHeader;
}

my $counter = 0;
my ($id, @rec) = ();

my $fh_stats;
if ($stats_file) {
    $fh_stats = new FileHandle " > $stats_file " or die "Can't create file $stats_file: $!\n";
} else {
    $fh_stats = \*STDOUT;
}
my %stats = ( "p" =>	"the read is paired in sequencing",
              "P" =>	"the read is mapped in a proper pair",
              "u" =>	"the query sequence itself is unmapped",
              "U" =>	"singletons (mates unmapped)",
              "r" =>	"strand of the query is reverse",
              "R" =>	"strand of the mate is reverse",
              "1" =>	"the read is the first read in a pair",
              "2" =>	"the read is the second read in a pair",
              "s" =>	"the alignment is not primary",
              "f" =>	"the read fails platform/vendor quality checks",
              "d" =>	"the read is either a PCR or an optical duplicate"
            );
my (%count,%mapq_mapped, %mapq_unmapped);
my @flags = qw(p P u U r R 1 2 s f d);
for my $flag (@flags,qw(total_beads total_reads total_removed total_mapped mapped_1st mapped_2nd both_mapped ) ) {
    $count{$flag} = 0;
}
my $maxInsertSize = 0;
my $minInsertSize =999999999999999;
my $totalInsertSize = 0;

while (<$fh_input>) {
    my @cols = split /\t/, $_;
    
    if (defined $id and $id eq $cols[0] ) { #same id
        push @rec, \@cols;
    } else { # a new id
        @rec = &process_rec(@rec);
        &output_rec($fh_output,@rec) if $clean;
        @rec = (\@cols);
        $id = $cols[0];
    }
}

if (@rec > 0) {
    @rec = &process_rec(@rec);
    &output_rec($fh_output,@rec) if $clean;
}

close $fh_input;

print $fh_stats "\nStats without the unmapped Bam file:\n";
&output_counts(1);

my $cmd_bam_unmapped = "samtools sort -n -o";
$cmd_bam_unmapped .= " -m $memory" if $memory;
$cmd_bam_unmapped .= " $bam_unmapped - | samtools view -X -";

my $fh_input_u = new FileHandle " $cmd_bam_unmapped | " or die "Can't open pipeline | $cmd_bam_unmapped:$!\n";
while (<$fh_input_u>) {
    my @cols = split /\t/, $_;
    if (defined $id and $id eq $cols[0] ) { #same id
        push @rec, \@cols;
    } else { # a new id
        @rec = &process_rec(@rec);
        @rec = (\@cols);
        $id = $cols[0];
    }
}

if (@rec > 0) {
    @rec = &process_rec(@rec);
}
close $fh_input_u;

print $fh_stats "\nStats with the unmapped Bam file:\n";
&output_counts();

close $fh_output if $clean;
close $fh_stats;
if ($remove) {
    print STDERR "Remving intermediate file $prefix\_sortByName.bam\n";
    system("rm -rf $prefix\_sortByName.bam $prefix\_sortByName.[0-9][0-9][0-9][0-9].bam") and warn "Failed to to remove files $prefix\_sortByName.bam $prefix\_sortByName.[0-9][0-9][0-9][0-9].bam:$!\n";
    system("rm -rf $prefix\_Clean.[0-9][0-9][0-9][0-9].bam") and warn "Failed to to remove files $prefix\_Clean.[0-9][0-9][0-9][0-9].bam:$!\n";
}
print STDERR "DONE!\n\n";

sub output_rec {
    my ($fh, @rec) = @_;
    for my $rec (@rec) {
        print $fh join("\t", @$rec);
    }
}

sub compare_rec {
    my ($rec1, $rec2) = @_;
    my $status = 0;
    if ( $rec1->[2] ne "*" and $rec2->[6] ne "*" and $rec1->[2] eq ($rec2->[6]eq"="?$rec2->[2]:$rec2->[6]) and $rec1->[3] eq $rec2->[7] ) {
        $status++;
    }
    if ( $rec2->[2] ne "*" and $rec1->[6] ne "*" and $rec2->[2] eq ($rec1->[6]eq"="?$rec1->[2]:$rec1->[6]) and $rec2->[3] eq $rec1->[7] ) {
        $status++;
    }
    return $status if $status;
    if ( $rec1->[2] eq "*" and $rec1->[6] eq "*" and $rec2->[2] eq "*" and $rec2->[6] eq "*" ) {
        $status = -1;
    }
    return $status;
}

sub get_best_pair {
    my ($rec,@pairs) = @_;
    my $best = 0;
    for my $pair (@pairs) {
        $best = $pair->[2] if $pair->[2] > $best;
    }
    my $count_best = 0;
    my $best_pair;
    for my $pair (@pairs) {
        if ($pair->[2] == $best) {
            $best_pair = $pair;
            $count_best++;
        }
    }
    if ($count_best > 1 and $debug == 10) {
        print "More than one best pairs\n";
        &output_rec(\*STDOUT, @$rec);
    }
    return ($best_pair);
}

sub process_rec {
    my (@rec) = @_;
    return if @rec == 0;
    &output_rec(\*STDOUT, @rec) if @rec == $debug; # 1,2,3 or 4 
    if ( $debug == 30 ) {
        print scalar @rec, "\n";
    }
    $count{"total_beads"}++;
    if (@rec == 1) {  # singles
        &get_counts(@rec);
    } elsif (@rec == 2) {  # pairs
        my $status = compare_rec($rec[0],$rec[1]);
        if ( $status == 0 ) {
            print STDERR "Ophan pairs:\n";
            &output_rec(\*STDERR, @rec);
            exit -1;
        }
        &get_counts(@rec);
    } elsif (@rec > 2) {  # needs cleanup
        my @paired_index;
        for my $i (0..($#rec-1)) { # finding pairs
             for my $j (($i+1)..$#rec) {
                my $status = compare_rec($rec[$i],$rec[$j]);
                if ( $status > 0 ) { # paired
                    push @paired_index, [$i,$j,$status];
                }
            }
        }
        if (@paired_index == 0) { # No pairs here? Shouldn't happen for this bug
            print STDERR "Unable to find a pair:\n";
            &output_rec(\*STDERR, @rec);
            exit -2;
        } elsif (@paired_index > 1) {
            @paired_index = &get_best_pair(\@rec, @paired_index);
        }
        my @saved_index = ($paired_index[0]->[0],$paired_index[0]->[1]);
        &get_counts(@rec[@saved_index]);
        for my $i (0..$#rec) {
            next if ( $i == $saved_index[0] or $i == $saved_index[1] );
            if ( $rec[$i]->[5] =~ /I|D/ ) { # saving all gapped alignments
                push @saved_index, $i; # doesn't really happen in the PE bam file although it is allowed by bam standard
              #  print STDERR "Got one\n";
            } else {
                if ($rec[$i]->[1] !~ /u/) { # should we die here too?
                    print STDERR "Mapped record removed ($i):\n";
                    &output_rec(\*STDERR, @rec);
                    exit -3;
                }
                if ( $debug == 20 ) {  # -d output removed records to STDOUT
                    &output_rec(\*STDOUT, $rec[$i]);
                }
                $count{total_removed}++;
            }
        }
        if (@saved_index > 3) { # This shouldn't happen
            print STDERR "More than 3 records from a pair including gapped alignments:\n";
            &output_rec(\*STDERR,@rec);
            exit -4;
        } elsif (@saved_index == 0) { # This shouldn't happen either, nothing left?
            print STDERR "Nothing left after filtering:\n";
            &output_rec(\*STDERR,@rec);
            exit -5;
        }
        @rec = @rec[@saved_index];
    }
    return @rec;
}

sub get_counts {
    my @rec = @_;
    foreach my $rec (@rec) {
        $count{total_reads}++;
        my @flags = split //, $rec->[1];
        foreach my $flag (@flags) {
            $count{$flag}++;
        }
        if ($rec->[1] =~ /u/) { #unmapped
            $mapq_unmapped{$rec->[4]}++;
            if ( $rec->[1] =~ /1/ ) {
                $count{unmapped_1st}++;
            } elsif ( $rec->[1] =~ /2/ ) {
                $count{unmapped_2nd}++;
            } else {
                print STDERR "Unknown tag order:\n";
                output_rec(\*STDERR,$rec);
            }
        } else { # mapped
            $mapq_mapped{$rec->[4]}++;
            $count{total_mapped}++;
            $count{both_mapped}++ if $rec->[1] !~ /U/;
            if ( $rec->[1] =~ /1/ ) {
                $count{mapped_1st}++;
            } elsif ( $rec->[1] =~ /2/ ) {
                $count{mapped_2nd}++;
            } else {
                print STDERR "Unknown tag order:\n";
                output_rec(\*STDERR,$rec);
            }
        }
        if ( $rec->[1] =~ /P/ ) { # proper pair
            my $size = abs($rec->[8]);
            $totalInsertSize += $size;
            $maxInsertSize = $size if ($size > $maxInsertSize);
            $minInsertSize = $size if ($size < $minInsertSize);
        }
    }
}

sub output_counts {
    my ($flag) = @_;
    print $fh_stats "\nTotal $count{total_beads} beads, $count{total_reads} reads";
    print $fh_stats " (after removing $count{total_removed} entries)";
    print $fh_stats "\n\n";
    print $fh_stats "Flag stats:\n";
    
    print $fh_stats "$count{p}\t$stats{p}\n";
    printf $fh_stats "    $count{1}\tTotal first read (%.2f%% total read)\n", $count{total_reads} ? (100*$count{1} / $count{total_reads}) : 0;
    printf $fh_stats "    $count{2}\tTotal second read (%.2f%% total read)\n", $count{total_reads} ? (100*$count{2} / $count{total_reads}) : 0;;
    
    printf $fh_stats "$count{u}\t$stats{u} (%.2f%% total read)\n", $count{total_reads} ? (100*$count{u} / $count{total_reads}) : 0;
    printf $fh_stats "    $count{unmapped_1st}\tUnmapped first read (%.2f%% total first read)\n", $count{1} ? (100*$count{unmapped_1st} / $count{1}) : 0;
    printf $fh_stats "    $count{unmapped_2nd}\tUnmapped second read (%.2f%% total second read)\n", $count{2} ? (100*$count{unmapped_2nd} / $count{2}) : 0;
    
    printf $fh_stats "$count{total_mapped}\tTotal mapped reads (%.2f%% total read)\n", $count{total_reads} ? (100*$count{total_mapped} / $count{total_reads}) : 0;
    printf $fh_stats "    $count{mapped_1st}\tmapped first read (%.2f%% total mapped, %.2f%% total first read)\n", 
                                  $count{total_mapped} ? (100*$count{mapped_1st} / $count{total_mapped}) : 0,
                                  $count{1} ? (100*$count{mapped_1st} / $count{1}) : 0;
    printf $fh_stats "    $count{mapped_2nd}\tmapped second read (%.2f%% total mapped, %.2f%% total second read)\n", 
                                  $count{total_mapped} ? (100*$count{mapped_2nd} / $count{total_mapped}) : 0,
                                  $count{2} ? (100*$count{mapped_2nd} / $count{2}) : 0;
    printf $fh_stats "    $count{both_mapped}\tboth reads mapped (%.2f%% total mapped, %.2f%% total read)\n", 
                                  $count{total_mapped} ? (100*$count{both_mapped} / $count{total_mapped}) : 0,
                                  $count{total_reads} ? (100*$count{both_mapped} / $count{total_reads}) : 0;
    printf $fh_stats "    $count{P}\t$stats{P} (%.2f%% total mapped, %.2f%% total reads)\n", 
                                  $count{total_mapped} ? (100*$count{P} / $count{total_mapped}) : 0, 
                                  $count{total_reads}? (100*$count{P} / $count{total_reads}) : 0;
    printf $fh_stats "$count{U}\t$stats{U} (%.2f%%)\n", $count{total_reads} ? (100*$count{U} / $count{total_reads}) : 0;
    for my $flag qw(r R s f d) {
        print $fh_stats "$count{$flag}\t$stats{$flag}\n";
    }
    print $fh_stats "\n";
    
    return if $flag;
    
    unless ($count{"P"} == 0) {
        printf $fh_stats "%d\tMean Insert Size\n", $totalInsertSize / $count{"P"};
        printf $fh_stats "%d - %d Insert Size Range\n", $minInsertSize, $maxInsertSize;
    }
    print $fh_stats "\n";
    
    print $fh_stats "#####Mapping Quality Distribution of Mapped Reads#####\n";
    print $fh_stats "MAPQ\t#Reads\t% Total\t#Cumulative\t% Cumulative\n";
    &print_mapq(\%mapq_mapped);
    
    print $fh_stats "\n#####Mapping Quality Distribution of Unmapped Reads#####\n";
    print $fh_stats "MAPQ\t#Reads\t% Total\t#Cumulative\t% Cumulative\n";
    &print_mapq(\%mapq_unmapped);
    print $fh_stats "\n";
}


sub print_mapq {
    my ($mapq) = @_;
    my $total = 0;
    my @keys = sort {$a<=>$b } keys %$mapq;
    for (@keys) {
        $total += $mapq->{$_};
    }
    my $cum = 0;
    for my $k (@keys) {
        my $per = 100*$mapq->{$k}/$total;
        $cum +=$mapq->{$k};
        my $cum_per= 100*$cum/$total;
        printf $fh_stats "%d\t%d\t%4.2f\t%d\t%4.2f\n", $k,$mapq->{$k},$per,$cum,$cum_per;
    }
}
