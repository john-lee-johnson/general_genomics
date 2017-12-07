#!/usr/bin/env perl
#===============================================================================
#
#         FILE: set_interval_tree.pl
#
#        USAGE: perl set_interval_tree.pl
#               --query=<path/to/bed>
#               --bed=</path/to/bed>
#               --bin=<INT>
#               --flank=<INT>
#               --outfile=<path/to/outfile.txt
#               --fast=<0/1>
#               --size=<INT>

#
#  DESCRIPTION: Finds the sequencing profile for a set of regions
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: Not yet implemented for bedpe files
#        NOTES: ---
#       AUTHOR: John L Johnson (jj), johnlee.johnson@gmail.com
# ORGANIZATION: University of Pennsylvaniaperl suppo
#      VERSION: 1.0
#      CREATED: 12/06/2017 08:15:42 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use v5.26.1;
use Getopt::Long;
use Set::IntervalTree;

my ($help, $query, $bed, $bin, $flank, $outfile, $chrHash, $size, $totalReads);
my ($fast) = 0; #Set --fast=1 for loading only reads that overlap intervals
my $mid_flank = 100; #Flank size for finding intensity in the middle

usage() if ( @ARGV <1 or join("",@ARGV)!~/--/ or
    ! GetOptions( 'help' => \$help,
        'query=s' => \$query,
        'bed=s' => \$bed,
        'bin=i' => \$bin,
        'flank=i' => \$flank,
        'outfile=s' => \$outfile,
        'fast=i' => \$fast,
        'size=i' => \$size,
    ) or
    defined $help or (!defined $query or !defined $bed or !defined $bin or !defined $flank or !defined $outfile));

my $time = time;

if ($fast == 0) {
    $totalReads = `cat $bed | wc -l` / 1000000;
    $chrHash = load_intervals($bed);
}
else {
    $totalReads = $size / 1000000;
    $chrHash = load_intervals_fast($bed);
}

my $queryHash = load_query($query);

open (OUT, ">", $outfile) or die "Could not open file: $outfile $!";
foreach my $chr(keys %{$queryHash}) {
    #print "$chrHash{$chr}\n";
    unless( exists $$chrHash{$chr} ) {
        warn "\t\tSkipping $chr. It doesn't exist in tag file.\n";
        next;
    }
    warn "\t\tWorking on $chr.\n";
    for(my $i = 0; $i < @{$$queryHash{$chr}}; $i++) {
        my ($start, $name) = @{$$queryHash{$chr}[$i]};
        my $stop = $start + 2*$flank;
        my $countBins = 0;
        my @hits;
        my $mid_start = $start + $flank - $mid_flank;
        my $mid_end = $start + $flank + $mid_flank;
        for(my $j = $start; $j <= $stop; $j += $bin) {
            $countBins++;
            my $overlaps;
            if( $bin == 1 ) {
                $overlaps = $$chrHash{$chr}->fetch($j, $j);
            }
            else {
                $overlaps = $$chrHash{$chr}->fetch($j, $j+$bin-1);
            }
            push @hits, scalar(@$overlaps) / $totalReads;
        }
        my $mid_overlaps = $$chrHash{$chr}->fetch($mid_start,$mid_end);
        my $mid_counts = scalar(@$mid_overlaps)/$totalReads;
        if ($countBins % 2 == 0) {
            warn "WARNING: NOT ODD NUMBER OF BINS $countBins\n";
        }
        print OUT "$name $mid_counts @hits\n";
    }
}
warn "\t\tFinished in ".(time - $time)." seconds.\n";
print "Total time to completion: ".(time - $^T)."\n";

sub load_intervals {
    my ($bed) = @_;
    open (IN, "<:crlf", $bed) or die "Could not open file: $bed $!";
    warn "\n\t\tLoading reads from $bed\n";
    my $time = time;
    my (%hash);
    while (<IN>) {
        chomp;
        my @coords = split '\t';
        if ( !exists $hash{$coords[0]} ) {
            $hash{$coords[0]} = Set::IntervalTree->new;
        }
        my $name = join "_", @coords;
        $hash{$coords[0]}->insert($name, $coords[1], $coords[2]);
    }
    close IN;
    warn "\t\tFinished in ".(time - $time)." seconds.\n";
    return \%hash;
}

sub load_intervals_fast {
    my ($bed) = @_;
    die "ERROR: A gzipped bedfile for Tabix was not provided!\n" unless $bed =~ /[.]gz$/;
    warn "\n\t\tLoading reads from $bed\n";
    $time = time;
    chomp(my @tags = `tabix $bed -R $query`);
    my (%hash);
    foreach my $line(@tags) {
        my @coords = split '\t', $line;
        if ( !exists $hash{$coords[0]} ) {
            $hash{$coords[0]} = Set::IntervalTree->new;
        }
        my $name = join "_", @coords;
        $hash{$coords[0]}->insert($name, $coords[1], $coords[2]);
    }
    warn "\t\tFinished in ".(time - $time)." seconds.\n";
    return \%hash;
}

sub load_query {
    my ($query) = @_;
    warn "\t\tCounting reads in regions from $bed.\n";
    $time = time;
    my (%queryHash);
    open (BED, "<", $query) or die "Could not open file: $query $!";
    while (<BED>) {
        chomp;
        my @coords = split /\s+/;
        my $center = $coords[1]+int(($coords[2]-$coords[1])/2);
        my $left = $center-$flank;
        my @array = ($left, "$coords[3].$coords[1]");
        push @{$queryHash{$coords[0]}}, \@array;
    }
    close BED;
    return \%queryHash;
}

sub usage {
    print "\nUnknown option: @_\n" if ( @_ );
    print "\nUsage: perl set_interval_tree.pl [OPTIONS] --query=<path/to/bed> --bed=</path/to/bed> --bin=<INT> --flank=<INT> --outfile=<path/to/outfile.txt>".
    "\n\t[DESCRIPTION]\n".
    "\tFinds the sequencing coverage for a set of regions (i.e., the profile)\n".
    "\n\t[OPTIONS]\n".
    "\t--help\t\t\tPrint this help message.\n".
    "\t--query=<path/to/bed>\tfilepath to bed file of regions to find profile\n".
    "\t  Info: Script will take the center of the interval\n".
    "\t        Column 4 should have unique identifier\n".
    "\t--bed=<path/to/bed>\tfilepath to bedpe/bed generated from bamtobed\n".
    "\t  Info: Currently works with bed, not finished for bedpe yet\n".
    "\t--bin=<INT>\tsize of bins\n".
    "\t--flank=<INT>\tsize of flanking region\n".
    "\t--outfile=<path/to/outfile.txt>\tfilepath to write output\n".
    "\t  Info: written as space-separated file with the first column being the middle intensity\n";
    exit;
}
