#!/usr/bin/env perl

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

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
my ( $help, $man, $DEBUG, $verbose );
my $geneDERnum={};
my $keyCol=13; # default is for the standard parseBlast.pl output
my @infoCol;
my $minDERout=11;
my $outfile;
my $noHeader;
my $delimiter="\t";
GetOptions(
		'help|?' => \$help, 
		'man' => \$man,
		'keyCol=i' => \$keyCol,
		'outfile=s' => \$outfile,
		'minDERout=i' => \$minDERout,
		'infoCol=s' => \@infoCol,
		'noHeader' => \$noHeader,
		'delimiter=s' => \$delimiter,
		'debug' =>\$DEBUG,
		'verbose' => \$verbose,
	) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

########## arguments ##########
pod2usage("***** Error: At least 1 input files are required *****") if (@ARGV < 1);
pod2usage("***** Error: Output file must be set *****") if ($outfile eq "");



########## program start ##########
printTime("Generate gene-based report starts at") if ($verbose);

my $colNum=scalar @infoCol;


# array index starts with zero
$keyCol--;
for (my $i=0; $i<$colNum; $i++) {
	$infoCol[$i]--;
}

my $outFileHD;
open $outFileHD, ">$outfile" || die $!;

my $filecount;
for ( $filecount=0; $filecount<(scalar @ARGV); $filecount++) {
	my $infile=$ARGV[$filecount];
	
	printTime("Processing file: $infile\n") if ($verbose);

	my $inFileHD;
	open $inFileHD, "<$infile" || die $!; 

	my $aLine;
	my $curIndex=0;
	my $iCol;
	my @lineSplit;
	# collect the header line
	if (!$noHeader) {
		$aLine=<$inFileHD>;
		@lineSplit = split($delimiter, $aLine);
		foreach $iCol (@infoCol) {
			$geneDERnum->{'header'}[$curIndex] = $lineSplit[$iCol];
			$curIndex++;
		}
		$geneDERnum->{'header'}[$curIndex+$filecount]=basename($ARGV[$filecount]);
		$geneDERnum->{'header'}[$curIndex+$filecount+1]='Total';
	}

	while ($aLine=<$inFileHD>) {
		@lineSplit = split($delimiter, $aLine);

		$curIndex=0;
		foreach $iCol (@infoCol) {
			$geneDERnum->{$lineSplit[$keyCol]}[$curIndex] = $lineSplit[$iCol];
			$curIndex++;
		}
		$geneDERnum->{$lineSplit[$keyCol]}[$curIndex+$filecount]++;
	}

	close $inFileHD;
}


# fill in the total occurrence column
for my $key (keys %$geneDERnum) {
	if ($key ne 'header') {
		my $total=0;
		for (my $i=$colNum; $i<($filecount + $colNum); $i++) {
			print "key: $key, i: $i, number: ", $geneDERnum->{$key}[$i], " total: $total\n" if $DEBUG;
			$total += $geneDERnum->{$key}[$i] if ($geneDERnum->{$key}[$i]);
		}
		$geneDERnum->{$key}[$colNum+$filecount]=$total;
	}
}


print $outFileHD 'Accession ID';
&printGeneLine('header');

# print all info out
for my $key (keys %$geneDERnum) {
	if ($key ne 'header' && $geneDERnum->{$key}[$colNum+$filecount] >= $minDERout) {
		print $outFileHD $key;
		&printGeneLine($key);
	}
}

close $outFileHD;

printTime("Generate gene-based report ends at") if ($verbose);




sub printGeneLine {
	my ($key) = @_;
	for (my $i=0; $i<($filecount + $colNum + 1); $i++) {
		if ($geneDERnum->{$key}[$i]) {
			print $outFileHD $delimiter . $geneDERnum->{$key}[$i];
		} else {
			print $outFileHD $delimiter . "0";
		}   
	}   
	print $outFileHD "\n";
}

__DATA__

=head1 NAME

getNumOfDER4Gene.pl - This script will collect the number of lines for
  the unique <keyCol> value in each input file. 
	The output will also display the columns listing in <infoCol>. 


=head1 SYNOPSIS

getNumOfDER4Gene.pl [-help] [-man] [-debug] [-verbose] [-keyCol <INTEGER>] 
[-minDERout <INTEGER>] [-infoCol <INTEGER> [-infoCol <INTEGER> ...]] 
[-noHeader] [-delimiter <SYMBOL>] <-outfile OUT_FILENAME> <InputFile1 [InputFile2 ...]>

B<Example>

getNumOfDER4Gene.pl -infoCol 16 -outfile -outfile LG1_LG2_genes.tsv 
'LG1_LG2_0.tsv' 'LG1_LG2_(0,0.5).tsv' 'LG1_LG2_(2,inf).tsv' 'LG1_LG2_inf.tsv' 


=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-debug>

Turn on the script debugging printout.

=item B<-verbose>

Prints processing time.

=item B<-keyCol>

The column to be treated as the key. Default 13.

=item B<-outfile>

The output filename

=item B<-minDERout>

Only the key value contains at least these number of lines will be outputed. Default is 11.

=item B<-infoCol>

The column with this column number will be outputed. Multiple values are accepted.

=item B<-noHeader>

Set it when the input files do not have a header

=item B<-delimiter>

Default delimiter in input file is '\t'. 

=back


=head1 DESCRIPTION

This script will collect the number of lines for
  the unique <keyCol> value in each input file. The input file should be
	in tubular form seperated by a specific delimiter.
	The output will also display the columns listing in <infoCol>. 


=head1 LICENSE

Copyright (c) 2012 Chon-Kit Kenneth Chan. All rights reserved.

This script is freely available for non-commercial use 
and is covered by the GNU General Public License v3 or later. 
See L<http://www.gnu.org/licenses/gpl.html>

The user of the script agrees to acknowledge the author(s) in any 
scientific publication of results based in part on the use of the script.

=cut

