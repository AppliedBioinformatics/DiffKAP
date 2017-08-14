#!/usr/bin/env perl

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Bio::SeqIO;
use Data::Dumper;

use Cwd 'abs_path';
use lib dirname(abs_path($0));
use Utilities;


BEGIN {
	STDERR->autoflush(1);
	STDOUT->autoflush(1);
}

########## options ##########
my ( $help, $man, $DEBUG, $verbose );
my ( $maxNumOfLine, $outFilename1 );
my $isNoSplit=0;
my $informat="fasta";
my $numLinePerRead;
GetOptions(
		'help|?' => \$help, 
		'man' => \$man,
		'debug' =>\$DEBUG,
		'verbose' =>\$verbose,
		'n=i' => \$maxNumOfLine,
		'informat=s' => \$informat,
		'numOfLinePerRead=i' => \$numLinePerRead,
		'outFile1=s' =>\$outFilename1,
	) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;


# Check for errors
pod2usage("***** Error: At least 1 input files is required *****") if (@ARGV < 1);


########## Program starts ##########
printTime("Remove repeats and split files starts at") if ($verbose);

# if $mxNumOfLine is not set, assuming all output to one file
$isNoSplit=1 if (!$maxNumOfLine);

my %seenReads;
my $totalClone=0;
my $totalReadCount = 0;
my $totalOutCount = 0;
my $outnum = 0;

my ($outFile1, $outSeq1);

# Setup output filename for read A
if (!$outFilename1) {
	$outFilename1=$ARGV[0]; 
	$outFilename1 =~ s/\.fastq$|\.fq$|\.txt$/.fasta/;
}

print "OutFilename_1: $outFilename1\n" if ($verbose);


if ("$informat" eq "fasta") {
	my $numSeqLine=$numLinePerRead-1;
	my $numSkipLine=0;
	$totalClone = &processing($numSeqLine, $numSkipLine);
} else {
	my $numSkipLine=$numLinePerRead/2;
	my $numSeqLine=$numSkipLine-1;
	$totalClone = &processing($numSeqLine, $numSkipLine);
}

print "$totalClone possible clonal read pairs removed.\n" if ($verbose);

printTime("Remove repeats and split files ends at") if ($verbose);



sub processing {
	my ($numSeqLine, $numSkipLine) = @_;
	
	my $totalClone=0;
	
	# This flag is to prevent generating new empty file.  
	my $isNewFileOk=1;
	# Go through each input files for splitting
	while (@ARGV>0) {
		my $inFile1 = shift(@ARGV);
		
		open IN, "<$inFile1" || die $!;
		while (my $readID=<IN>) {
			$readID=~s/^@/>/;
			#print Dumper($readID);
			my $inSeq1="";
			my $tmp;
			for (my $i=0; $i<$numSeqLine; $i++) {
				$tmp=<IN>;
				chomp($tmp);
				$inSeq1="$inSeq1"."$tmp";
			}
			for (my $i=0; $i<$numSkipLine; $i++) {
				$tmp=<IN>;
			}
			# Open a new output file
			if ($totalOutCount==0 || ($isNewFileOk && !$isNoSplit && (($totalOutCount % $maxNumOfLine)==0)) ) {
				close OUT;
				$outFile1 = qq|${outFilename1}.$outnum|;
				open OUT, ">$outFile1";
				$isNewFileOk=0;
				$outnum++;
			}

			# Handle the trimming and removing repeats
			$seenReads{$inSeq1}++;
			if ($seenReads{$inSeq1} <= 1) {
				print OUT "$readID";
				print OUT "$inSeq1\n";
				$isNewFileOk=1;
				$totalOutCount++;
			} else {
				$totalClone++;
			}

			$totalReadCount++;	
		}	
	}
	
	return $totalClone;
}




__DATA__

=head1 NAME

removeRepeatSplitFile.pl - This script will remove repeat reads and splits all input fasta (or fastq) files into smaller files in the same or other format.


=head1 SYNOPSIS

removeRepeatSplitFile.pl [-help] [-man] [-debug] [-verbose] [-n <MAX_NUM_OF_ENTRIES>]
[-compareLen <INTEGER>] [-numOfRepAllow <INTEGER>] [-outFile <OUT_FILE_NAME>]
[-readType <single|pair>] 
[-outFile1 <READ_A_OUTPUT_FILE>] [-outFile2 <READ_B_OUTPUT_FILE>] 
[-informat <fasta|fastq>] [-outformat <fasta|fastq>] <IN_FILE> [<IN_FILE2> ...]

B<Example>

removeRepeatSplitFile.pl -n 1000000 -numOfRepAllow 2 -informat fastq -outformat fastq myInputFile1.fastq myInputFile2.fastq 

removeRepeatSplitFile.pl -n 1000000 -readType pair myInputFile1-A.fastq myInputFile1-B.fastq 
myInputFile2-A.fastq myInputFile2-B.fastq 

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-debug>

Turn on the script debugging printout.

=item B<-verbose>

Prints progress and # of possible clonal read removed.

=item B<-compareLen>

The nucleotide lenght which will be used for checking if a read is repeated. Default is to use the whole read.

=item B<-numOfRepAllow>

The number of repeating reads will still be outputed. Default is 1 (output only the unique read).

=item B<-readType>

If they are single reads or pair-end reads. Default is single. 
For pair-end reads, the input file list should be in a sequence of 
ReadA_seqFile ReadB_seqFile ReadA_seqFile ReadB_seqFile.

=item B<-n>

The maximum number of extries allowed in the each of the split files. Default is infinity so all will be put to a single file.

=item B<-informat>

The input file format, can be fasta or fastq format. Default is fasta.

=item B<-outformat>

The output file format. If -informat is fastq, this one can be fasta or fastq, otherwise, can only be fasta. 
Default is fasta.

=item B<-outFile1>

The output filename (for the 1st read in the pair-end read data). If not specified, the first input filename (with the extension modified) will be used.

=item B<-outFile2>

The output filename for read 2. This only has effect on pair-end read data. 


=back


=head1 DESCRIPTION

This script removes all repeating reads and splits a list of FASTA|FASTQ files 
into several multiple FASTA|FASTQ files. Each 
generated FASTA|FASTQ file will contain no more than the specified number of FASTA|FASTQ entries.


=head1 LICENSE

Copyright (c) 2012 Chon-Kit Kenneth Chan. All rights reserved.

This script is freely available for non-commercial use 
and is covered by the GNU General Public License v3 or later. 
See L<http://www.gnu.org/licenses/gpl.html>

The user of the script agrees to acknowledge the author(s) in any 
scientific publication of results based in part on the use of the script.

=cut

