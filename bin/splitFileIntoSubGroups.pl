#!/usr/bin/env perl

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper qw(Dumper);
use File::Basename;
use constant INF => 9**9**9; # an infinite number

# global variables
my ( $help, $man, $DEBUG, $noheader, $verbose );
my $delimiter='\t';
my $column;
my @range;  # Store the range information for each file


BEGIN {
	STDERR->autoflush(1);
	STDOUT->autoflush(1);
}

########## functions ##########
sub printTime {
	my $message;
	if (@_ >= 1) {
		($message) = @_;
	} else {
		$message = "Current Time";
	}
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	printf "$message: %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec; 
}


########## options ##########

GetOptions(
		'help|?' => \$help, 
		'man' => \$man,
		'debug' =>\$DEBUG,
		'verbose' => \$verbose,
		'delimiter=s' =>\$delimiter,
		'column=i' =>\$column,
		'range=s' =>\@range,
		'noheader' => \$noheader,
	) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(1) if ( !$column || !@range );
pod2usage(1) if ( !$column );
pod2usage(-verbose => 2) if $man;
# In array, the index starts with zero
$column--;


########## arguments ##########
pod2usage("***** Error: this script only accept one input argument *****") if (@ARGV != 1 );


########## program start ##########
printTime("Split results into sub groups start at") if ($verbose);

my ($filename1, $dir1, $ext1) = fileparse($ARGV[0], qr/\.[^.]*/);
print "$filename1, $dir1, $ext1\n" if $DEBUG;

print Dumper(\@ARGV) if $DEBUG;

my @filehandles;

# open the files
open IN1, "<", $ARGV[0] || die $!;

my $numRangeCount=0;
my $range_map=[];
my $description;
# Convert @range into @range_map
foreach my $aRange (@range) {
	my $singRange;
	if (($singRange,$description) = $aRange =~ m/^([0-9.inf]+),?(.*)$/ ) {
		$range_map->[$numRangeCount][0]="=";
		$range_map->[$numRangeCount][1]=$singRange;
	} else {
		my ($l, $ln, $un, $r);
		($l, $ln, $un, $r, $description) = $aRange =~ m/^([\[\(])([\d\.inf]+),([\d\.inf]+)([\)\]]),?(.*)$/;

		# ignore this range if it doesn't match the format
		if (!$r) {
			print STDERR "Range does not recognise, ignore *** $aRange ***\n";
			next;
		}

		if ($l eq "[") {
			$range_map->[$numRangeCount][0]=">=";
		} else {
			$range_map->[$numRangeCount][0]=">";
		}

		$range_map->[$numRangeCount][1]=$ln;

		if ($r eq "]") {
			$range_map->[$numRangeCount][2]="<=";
		} else {
			$range_map->[$numRangeCount][2]="<";
		}
		$range_map->[$numRangeCount][3]=$un;
	}

	# make array of $count file handles
	#localize the file glob, so FILE is unique to the inner loop
	my $FILE ;
	if ($description) {
		open($FILE, ">${dir1}${filename1}_${description}${ext1}") || die $!;
	} else {
		# Replaced symble with letter for filename
		$aRange =~ s/\[/ge/;
		$aRange =~ s/\(/gt/;
		$aRange =~ s/,/-/;
		$aRange =~ s/\]/le/;
		$aRange =~ s/\)/lt/;

		open($FILE, ">${dir1}${filename1}_${aRange}${ext1}") || die $!;
	}
	# Keep the file handle organised
	$filehandles[$numRangeCount]=$FILE;

	$numRangeCount++;
}

my $file;
my $line;


# Put the header line to each output file
if (!$noheader) {
	$line=<IN1>;
	foreach $file (@filehandles) {
		print $file $line;
	}
}


# go through each line and distribute it to the right file
while ($line=<IN1>) {
	print "line: $line    delimiter: $delimiter\n" if $DEBUG;

	if ($line!~m/^\s*$/) {
		# Get the column value based on $column and $delimiter.
		my @matches= split($delimiter,$line);
		print $matches[$column]."\n" if $DEBUG;

		my $curValue = $matches[$column];
		
		# compare with values in the $range_map and output the line to the corresponding file.
		for (my $i=0; $i < (scalar @$range_map); $i++) {

			$file = $filehandles[$i];
			my $lb = $range_map->[$i][0];
			my $min = $range_map->[$i][1];
			my $rb = $range_map->[$i][2];
			my $max = $range_map->[$i][3];

			if ($lb eq "=") {
				print $file $line if ($min == $curValue);
			} elsif ($lb eq ">") {
				if ( $curValue > $min) {
					if ($rb eq "<") {
						print $file $line if ($curValue < $max);
					} elsif ($rb eq "<=") {
						print $file $line if ($curValue <= $max);
					}
				}
			} elsif ($lb eq ">=") {
				if ($curValue >= $min) {
					if ($rb eq "<") {
						print $file $line if ($curValue < $max);
					} elsif ($rb eq "<=") {
						print $file $line if ($curValue <= $max);
					}
				}
			}
		}
	}
}


# close all files
foreach $file (@filehandles) {
  close $file;
}


printTime("Split results into sub groups ends at") if ($verbose);


__DATA__

=head1 NAME

splitFileIntoSubGroups.pl - This script will split a tabular file into sub files.


=head1 SYNOPSIS

splitFileIntoSubGroups.pl [-help] [-man] [-debug] [-verbose] [-delimiter <SYMBOL>]
[-noheader] -column <INTEGER> -range <RANGE1> [-range <RANGE2> ...] 
<INPUT FILE> 

B<Example>

splitFileIntoSubGroups.pl -column 4 -r \"0,1-Absent-in-${t1ID}\" -r \"(0,0.5),2-Highly-expressed-in-${t2ID}\" DERfile.tsv


=head1 OPTIONS

=over 4

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-debug>

Turn on the script debugging printout.

=item B<-verbose>

Prints processing time.

=item B<-delimiter>

The delimiter used in the input file.

=item B<-noheader>

Default the input file contains a header line. Turn this flag on if it is not. 

=item B<-column>

The column number used for the range checking.

=item B<-range>

The range to be split into.
The <RANGE> format: [|(LOWER_VALUE,UPPER_VALUE)|],DESCRIPTION, 
eg. '[0.5,2],Not Differentially Expressed'
'(0,0.5),Treatment 2 differentially expressed'

=item B<-numOfDEK>

Accept values:
  
	all - a read is a DER when all kmers in the read are DEK

	INTEGER - a read is a DER when the number of DEK in the read are greater than or equal to 'INTEGER'

=item B<-kmerFileFormat>

Accept values:

  singleNum - this format contains 2 columns: 1. kmer  2. a number which can be the occurrence
	meanStd - 5 columns: 1. mean value in 1st treatment (T1)  2. standard deviation in T1  3. mean value in 2nd treatment (T2)  4. std in T2  5. the p-value in t-test 
	occurNratio - 3 columns: 1. occurrence in T1  2. occurrence in T2  3. occurrence ratio (T1/T2)

=item B<-reverse>

Whether treatment the reverse compliment of a DEK as another DEK. Default off.

=item B<-outCount>

Whether outputting the 'kmer size' vs '# of kmer in a read' for deciding how many # of DEK in a read to be considered as a DER. Default off.

=item B<-delimiter>

Use a customised delimiter for the DEK files. Default delimiter is a tab '\t'

=back


=head1 DESCRIPTION

This script finds differentially expressed reads (DER) based on checking 
the number of kmer belonging to the differentially expressed kmer (DEK) input file list.


=head1 LICENSE

Copyright (c) 2012 Chon-Kit Kenneth Chan. All rights reserved.

This script is freely available for non-commercial use 
and is covered by the GNU General Public License v3 or later. 
See L<http://www.gnu.org/licenses/gpl.html>

The user of the script agrees to acknowledge the author(s) in any 
scientific publication of results based in part on the use of the script.

=cut

