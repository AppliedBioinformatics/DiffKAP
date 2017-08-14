#!/usr/bin/env perl

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper qw(Dumper);
use Bio::SeqIO;
use Bio::Perl;
use File::Basename;
use Statistics::Descriptive;
use constant PAD_N => 10; # the number of N to be pad between ReadA and ReadB
use constant PRINT_INTERVAL => 5000000; # print time after every time processed this amount
use constant INF => 9**9**9; # an infinite number

# global variables
my ( $help, $man, $DEBUG, $verbose );
my %allDEK;
my @kmerInReadCount;
my $kmerLen=0;

my ( $a, $b );

my $seqFileFormat='fasta';  # values: 'fastq'
my $numOfDEK='all';  # value: integer
my $kmerFileFormat='meanStd';  # values: 'singleNum', 'occur', 'occurNratio
my $delimiter='\t';

# whether to output counts for '# of kmer in a read'
my $outCount;
my $reverse;

BEGIN {
	STDERR->autoflush(1);
	STDOUT->autoflush(1);
}

my ($DERseqname, $DERinfoname, $DERcountname, $DERseq, $DERinfo, $DERcount);


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

# Read in DEKs and their information into a hash for all the DEK files
#   and determine kmer size based on the DEK size in file
sub readInDEK {
	  my ($delimiter) = @_;

		foreach my $kmerFile (@ARGV) {
			my $kmerIN;
			# Exit if any kmer file can not be open
			if (! open $kmerIN, "<$kmerFile") {
				print STDERR "ERROR: Can't open kmer file $kmerFile\n";
				exit;
			}

			# Go through each line to collect kmer
			while (my $kLine = <$kmerIN>) {
				my ($kmer, $kmerInfo) = $kLine =~ m/^(\w+)${delimiter}(.+)$/;

				if (!$kmer) {
					print STDERR "ERROR: Using delimiter $delimiter - could not parse this line: $kLine\n";
					next;
				}

				if (!$kmerLen) {
					$kmerLen=length($kmer);
				}
				$allDEK{$kmer}=$kmerInfo;
				$allDEK{(revcom($kmer))->seq}=$kmerInfo if ($reverse);
			}
		}
		
		if ($kmerLen == 0) {
			print STDERR "ERROR: kmer lenght equals to zero\n";
			exit;
		}
}


sub countDEKinRead {
	my ($seq_ref, $kmerLen, $totalKmer_ref, $numOfDEK_ref, $DEKinfo_ref) = @_;
	my $seqLen=length($$seq_ref);
	my $numOfKmerInCurRead = $seqLen - $kmerLen + 1;

	for (my $start=0; $start < $numOfKmerInCurRead; $start++) {
		$$totalKmer_ref++;
		my $curKmer = substr $$seq_ref, $start, $kmerLen;

		if ($allDEK{$curKmer}) {
			$$numOfDEK_ref++;
			push(@$DEKinfo_ref, $allDEK{$curKmer});
		}
	}
}


sub printDERseq {
	my ($DERseq, $seqA, $seqB, $padN) = @_;
	if ($seqB) {
		my $newSeq = Bio::Seq->new( -seq => $seqA->seq .'N'x$padN. $seqB->seq,
																-id => $seqA->id );
		$DERseq->write_seq($newSeq);
	} else {
		$DERseq->write_seq($seqA);
	}
}


sub printNumOfDERInRead {
	my ($DERcount, $countDER_ref) = @_;
	for(my $index=0; $index<(scalar @$countDER_ref); $index++) {
		print $DERcount "$index\t".$countDER_ref->[$index]."\n" if ($countDER_ref->[$index]);
	}
}


sub printDERInfo {
	my ($DERinfoFH, $type, $id_ref, $DEKinfo_ref) = @_;

	my (@occurT1, @occurT2, @occurExtra);
	if ($kmerFileFormat eq 'occurNratio') {
		if ($type eq 'header') {
			print $DERinfoFH "ReadID\tMedian-T1\tMedian-T2\tRatio of Median (RoM)\tCV-T1\tCV-T2\tMedian-Ratio\n";
		} else {
			&parse3sets($DEKinfo_ref, \@occurT1, \@occurT2, \@occurExtra, 0, 1, 2);
			&print3sets($DERinfoFH, $id_ref, \@occurT1, \@occurT2, \@occurExtra);
		}
	} elsif ($kmerFileFormat eq 'meanStd') {
		if ($type eq 'header') {
			print $DERinfoFH "ReadID\tMedian-T1\tMedian-T2\tRatio of Median (RoM)\tCV-T1\tCV-T2\tMedian-pValue\n";
		} else {
			# The meanT1, meanT2 and pvalue are is value columns 0, 2 and 4 respectively
			&parse3sets($DEKinfo_ref, \@occurT1, \@occurT2, \@occurExtra, 0, 2, 4);
			&print3sets($DERinfoFH, $id_ref, \@occurT1, \@occurT2, \@occurExtra);
		}
	} elsif ($kmerFileFormat eq 'singleNum') {
		if ($type eq 'header') {
			print $DERinfoFH "ReadID\tMedian\tMean\tCV\n";
		} else {
		}
	}
}


sub parse3sets {
	my ($DEKinfo_ref, $occurT1_ref, $occurT2_ref, $occurPvalue_ref, $one, $two, $three) = @_;

	foreach my $aInfo (@$DEKinfo_ref) {
		my @elements = split($delimiter, $aInfo);
		push(@$occurT1_ref, $elements[$one]);
		push(@$occurT2_ref, $elements[$two]);
		push(@$occurPvalue_ref, $elements[$three]);
	}
}


sub print3sets {
	my ($DERinfoFH, $id_ref, $occurT1_ref, $occurT2_ref, $occurExtra_ref) = @_;

	print "occurT1:\n" if $DEBUG;
	print Dumper($occurT1_ref) if $DEBUG;
	my ($medianT1, $meanT1, $cvT1);
	&calMedianMeanCV($occurT1_ref, \$medianT1, \$meanT1, \$cvT1);

	print "occurT2:\n" if $DEBUG;
	print Dumper($occurT2_ref) if $DEBUG;
	my ($medianT2, $meanT2, $cvT2);
	&calMedianMeanCV($occurT2_ref, \$medianT2, \$meanT2, \$cvT2);

	my $RoM; 
	if ($medianT2 == 0) {
		$RoM = INF;
	} else {
		$RoM = $medianT1 / $medianT2;
	}

	print "ocurExtra:\n" if $DEBUG;
	print Dumper($occurExtra_ref) if $DEBUG;
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@$occurExtra_ref);
	my $medianExtra = $stat->median();

	print $DERinfoFH "$$id_ref\t$medianT1\t$medianT2\t$RoM\t$cvT1\t$cvT2\t$medianExtra\n";
}


sub calMedianMeanCV {
	my ($array_ref, $median_ref, $mean_ref, $cv_ref) = @_;

	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@$array_ref);
	$$median_ref = $stat->median();
	$$mean_ref = $stat->mean();
	if ($$mean_ref == 0) {
		$$cv_ref = INF;
	} else {
		$$cv_ref = $stat->standard_deviation()/$$mean_ref*100;
	}
}



########## options ##########

GetOptions(
		'help|?' => \$help, 
		'man' => \$man,
		'debug' =>\$DEBUG,
		'verbose' =>\$verbose,
		'a=s' =>\$a,
		'b=s' =>\$b,
		'seqFileFormat=s' =>\$seqFileFormat,
		'numOfDEK=s' =>\$numOfDEK,
		'kmerFileFormat=s' =>\$kmerFileFormat,
		'reverse' =>\$reverse,
		'outCount' =>\$outCount,
		'delimiter=s' =>\$delimiter,
		'derSeqOutFile=s' =>\$DERseqname,
		'derInfoOutFile=s' =>\$DERinfoname,
	) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(1) if ( !$a );
pod2usage(1) if ($seqFileFormat ne 'fasta' && $seqFileFormat ne 'fastq' );
pod2usage(1) if ($numOfDEK ne 'all' && $numOfDEK !~ /^\d+$/ );
pod2usage(1) if ($kmerFileFormat ne 'singleNum' && $kmerFileFormat ne 'meanStd' && $kmerFileFormat ne 'occurNratio');
pod2usage(-verbose => 2) if $man;


########## arguments ##########
pod2usage("***** Error: require at least one kmer file as input argument *****") if (@ARGV < 1 );


########## program start ##########
printTime("Identify DER starts at") if ($verbose);

print Dumper(\@ARGV) if $DEBUG;


if ($outCount) {
	$DERcountname = $a;
	$DERcountname =~ s/\./_DER_count./;
	open $DERcount, ">$DERcountname" or die $!;
	print "***** -outCount is enable - only $DERcountname will be outputed. *****\n";
	print "***** Disable -outCount to output DER fasta and info files. *****\n";
	print "Output file: $DERcountname\n" if ($verbose);
} else {
	unless ($DERseqname) {
		$DERseqname = $a;
		$DERseqname =~ s/\./_DER_seq./;
	}
	$DERseq = Bio::SeqIO->new(-file => ">$DERseqname", -format => 'fasta');

	unless ($DERinfoname) {
		$DERinfoname = $a;
		$DERinfoname =~ s/\./_DER_info./;
	}
	open $DERinfo, ">$DERinfoname" or die $!;
	&printDERInfo($DERinfo, 'header');
	print "Output file1: $DERseqname\n" if ($verbose);
	print "Output file2: $DERinfoname\n" if ($verbose);
}

&readInDEK($delimiter);
print Dumper(\%allDEK) if $DEBUG;

my $seqAobj;
my $seqBobj;

# Only accept 'fasta' or 'fastq' format
# Open SeqAobj and SeqBobjbased on $seqFileFormat;
if ($seqFileFormat eq 'fasta') {
	$seqAobj = Bio::SeqIO->new(-file => $a, -format => 'fasta');
	$seqBobj = Bio::SeqIO->new(-file => $b, -format => 'fasta') if ($b);
} else {
	$seqAobj = Bio::SeqIO->new(-file => $a, -format => 'fastq');
	$seqBobj = Bio::SeqIO->new(-file => $b, -format => 'fastq') if ($b);
}


my $numOfProcessSeq=0;

# Go through SeqAObj:
while (my $seqA = $seqAobj->next_seq) {
	my $totalKmerInA=0; 
	my $numOfDEKInA=0; 
	my $totalKmerInB=0; 
	my $numOfDEKInB=0; 
	my @DEKinfo;

	&countDEKinRead(\($seqA->seq), $kmerLen, \$totalKmerInA, \$numOfDEKInA, \@DEKinfo);

	my $seqB;
	if ($b) {
		$seqB = $seqBobj->next_seq;
		&countDEKinRead(\($seqB->seq), $kmerLen, \$totalKmerInB, \$numOfDEKInB, \@DEKinfo);
	}

	my $totalKmerInRead = $totalKmerInA + $totalKmerInB;
	my $numOfDEKInRead = $numOfDEKInA + $numOfDEKInB;

	if ($outCount) {
		# Keep track of the '# of DEK in a read' data
		$kmerInReadCount[$numOfDEKInRead]++ if ($outCount);
	} else {
		# Print out the DER
		my $requireNumOfDEKInRead=$numOfDEK;
		$requireNumOfDEKInRead=$totalKmerInRead if ($numOfDEK eq 'all');

		print "totalKmerInA: $totalKmerInA, InB: $totalKmerInB; numOfDEKInA: $numOfDEKInA, InB: $numOfDEKInA; requireNumOfDEKInRead: $requireNumOfDEKInRead\n" if $DEBUG;
		print Dumper(\@DEKinfo) if $DEBUG;

		if ($numOfDEKInRead >= $requireNumOfDEKInRead) {
			&printDERseq($DERseq, $seqA, $seqB, PAD_N);
			&printDERInfo($DERinfo, 'data', \($seqA->id), \@DEKinfo);
		}
	}

	$numOfProcessSeq++;
	printTime ("Processed $numOfProcessSeq reads") if ($numOfProcessSeq % PRINT_INTERVAL == 0 && $verbose); 
}

&printNumOfDERInRead($DERcount, \@kmerInReadCount) if ($outCount);

printTime("Identify DER ends at") if ($verbose);


__DATA__

=head1 NAME

identifyDER.pl - This script finds differentially expressed reads (DER) based on checking 
									the number of kmer belonging to the differentially expressed kmer (DEK) list.


=head1 SYNOPSIS

identifyDER.pl [-help] [-man] [-debug] [-verbose] -a <FileReadA> [-b <FileReadB>] 
[-seqFileFormat <fasta|fastq>] [-numOfDEK <all|INTEGER>]
[-kmerFileFormat <singleNum|meanStd|occurNratio>] [-reverse] [-outCount] [-delimiter <DELIMITER>]
[-derSeqOutFile <FileNamee>] [-derInfoOutFile <FileName>]
<DEKfile1> [DEKfile2 DEKfile3 ...]

B<Example>

identifyDER.pl -a LG_12_1_sequence.fasta -b LG_12_2_sequence.fasta -outCount -kmerFileFormat occurNratio LG1_match_norm_LG2_DEK.tsv
identifyDER.pl -a LG_12_1_sequence.fasta -b LG_12_2_sequence.fasta -kmerFileFormat occurNratio LG1_match_norm_LG2_DEK.tsv


=head1 OPTIONS

=over 4

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-debug>

Turn on the script debugging printout.

=item B<-verbose>

Prints progress.

=item B<-seqFileFormat>

Accept values: 

  fasta - input files are fasta format

  fastq - input files are fastq format

=item B<-a>

ReadA sequence file for pair-end reads. If only -a is presents, it is a sequence file in single reads.

=item B<-b>

ReadB sequence file for pair-end reads. 

=item B<-numOfDEK>

Accept values:
  
  all - a read is a DER when all kmers in the read are DEK

  INTEGER - a read is a DER when the number of DEK in the read 
	are greater than or equal to 'INTEGER'. This number has to be greater than 0.

=item B<-kmerFileFormat>

Accept values:

  singleNum - this format contains 2 columns: 1. kmer  2. a number which can be the occurrence
  meanStd - 5 columns: 1. mean value in 1st treatment (T1)  2. standard deviation in T1  
    3. mean value in 2nd treatment (T2)  4. std in T2  5. the p-value in t-test 
  occurNratio - 3 columns: 1. occurrence in T1  2. occurrence in T2  3. occurrence ratio (T1/T2)

=item B<-reverse>

Whether treatment the reverse compliment of a DEK as another DEK. Default off.

=item B<-outCount>

Whether outputting the 'kmer size' vs '# of kmer in a read' for deciding how many # of DEK in a read to be considered as a DER. Default off.

=item B<-delimiter>

Use a customised delimiter for the DEK files. Default delimiter is a tab '\t'

=item B<derSeqOutFile>

The file name for outputting the DER sequence.

=item B<derInfoOutFile>

The file name for outputting the DER infomation. 

=back


=head1 DESCRIPTION

This script finds differentially expressed reads (DER) based on checking 
the number of kmer belonging to the differentially expressed kmer (DEK) input file list.
The DEK filename should contain the word '_DEK' before the first left-most dot.


=head1 LICENSE

Copyright (c) 2012 Chon-Kit Kenneth Chan. All rights reserved.

This script is freely available for non-commercial use 
and is covered by the GNU General Public License v3 or later. 
See L<http://www.gnu.org/licenses/gpl.html>

The user of the script agrees to acknowledge the author(s) in any 
scientific publication of results based in part on the use of the script.

=cut

