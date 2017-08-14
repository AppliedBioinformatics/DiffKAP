#!/usr/bin/env perl

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Copy;
use File::Path qw(make_path remove_tree);
use Data::Dumper qw(Dumper);

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

# Expand the array list by adding each element in the first list by the set of the 2nd list
sub addOneNC2List {
	my ($refNcList, $refNc) = @_;
	my $NcLength=scalar @$refNcList;
	# Handle the empty array
	if ($NcLength == 0) {
		foreach my $n (@$refNc) {
			push(@$refNcList, $n);
		}
	} else {
		# when it is not empty
		for( my $i=0; $i<$NcLength; $i++ ) {
			my $orgNc = shift(@$refNcList);
			foreach my $n (@$refNc) {
				my $newNc = $orgNc . $n;
				push(@$refNcList, $newNc);
			}
		}
	}
}


########## options ##########
my ( $help, $man, $DEBUG, $verbose, $splitDir );
my $splitNum=2;
GetOptions(
		'help|?' => \$help, 
		'man' => \$man,
		'debug' => \$DEBUG,
		'verbose' => \$verbose,
		'outDir=s' => \$splitDir,
		'splitNum=i' => \$splitNum,
	) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

########## arguments ##########
pod2usage("***** Error: Incorrect num of arguments *****") if (@ARGV != 1);


########## program start ##########
printTime("Split kmer file starts at") if ($verbose);

my ($filename1, $dir1, $ext1) = fileparse($ARGV[0], qr/\.[^.]*/);
print "$filename1, $dir1, $ext1\n" if $DEBUG;

print("Filename: $filename1\n") if ($verbose);
print("Split number: $splitNum\n") if ($verbose);

my @nc=('a','c','g','t');
my @combNC=( );
my $matchPattern='^(';

for (my $i=0; $i<$splitNum; $i++) {
	addOneNC2List(\@combNC, \@nc);
	$matchPattern.='\w';
}
$matchPattern.= ').*';

print Dumper(\@combNC) if $DEBUG;

my %ncMap;
my $count = scalar @combNC;
for (my $i=0; $i<$count; $i++) {
    $ncMap{$combNC[$i]}=$i;
}

print Dumper(\%ncMap) if $DEBUG;


my @filehandles;
# make array of $count file handles
for(my $i=0; $i<$count; $i++)
{
    #localize the file glob, so FILE is unique to the inner loop
    local *FILE;
    open(FILE, ">$dir1$filename1$ext1.$i") || die $!;
    #push the typeglobe to the end of the array
    push(@filehandles, *FILE);
}

# open the files
open IN1, "<", $ARGV[0] || die $!;

# go through each line and distribute it to the right file
my $file;
while (my $line=<IN1>) {
  my ($matches)=$line=~m/$matchPattern/;
	$matches=lc($matches);
  print "$line: $matches, $ncMap{$matches}\n" if $DEBUG;
  $file = $filehandles[$ncMap{$matches}];
  print $file "$line";
}


# close all files
for (my $i=0; $i<$count; $i++) {
  $file = $filehandles[$i];
  close $file;
}


# creat the split dir and move the splited files into the dir
$splitDir = "$dir1$filename1"."_split" unless ($splitDir);
make_path($splitDir) if (!-d $splitDir);
for (my $i=0; $i<$count; $i++) {
	move("${dir1}${filename1}${ext1}.${i}",$splitDir);
}



printTime("Split kmer file ends at") if ($verbose);



__DATA__

=head1 NAME

splitFileInMerGroup.pl - This script split a kmer file into [splitNum] of subfiles based on the first [splitNum]mer of the kmer


=head1 SYNOPSIS

splitFileInMerGroup.pl [-help] [-man] [-debug] [-verbose] [-outDir <outputDirectory>] [-splitNum <splitNum>] INFILE

B<Example>

splitFileInMerGroup.pl -splitNum 3 kmerOccurrenceFile.tsv


=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-debug>

Print out the messages for debugging.

=item B<-verbose>

Prints processing time.

=item B<-outDir>

The directory name for storing the splitted files.

=item B<-splitNum>

Control the number of nucleotide used for splitting. Default is 2.

=back


=head1 DESCRIPTION

This script split a kmer file into [splitNum] of subfiles based on the first [splitNum]mer of the kmer.
The kmer string should be in the first column of the kmer file.


=head1 LICENSE

Copyright (c) 2012 Chon-Kit Kenneth Chan. All rights reserved.

This script is freely available for non-commercial use 
and is covered by the GNU General Public License v3 or later. 
See L<http://www.gnu.org/licenses/gpl.html>

The user of the script agrees to acknowledge the author(s) in any 
scientific publication of results based in part on the use of the script.

=cut

