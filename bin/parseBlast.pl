#!/usr/bin/env perl

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper qw(Dumper);
use Bio::SearchIO;
use Bio::SeqIO;
use constant PRINT_INTERVAL => 5000000; # print time after every time processed this amount
use constant INF => 9**9**9; # an infinite number

# global variables
my ( $help, $man, $DEBUG, $verbose );
my $seqFileFormat='fasta';
my $inDelimiter="\t";
my $outDelimiter="\t";
my $infoNoHeader;
my ( $noSplit, $noOutHeader, $onlyHits);
my $numOfHSP=1;


my $eval=INF;
my @seqFiles;
my @infoFiles;
my $outFile;

my %allSeq;
my %allInfo;

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



# Put all sequences into a hash
sub readInSeq {
	my ($seqFiles_ref, $seqFormat, $allSeq_ref) = @_;

	foreach my $file (@$seqFiles_ref) {
		my $seqio_obj = Bio::SeqIO->new(-file => $file,
																		-format => $seqFormat );
		while (my $seq = $seqio_obj->next_seq) {
			print "seqID: $seq->id\n" if $DEBUG;
			$allSeq_ref->{$seq->id} = $seq->seq;
		}
	}
}


# Put all information into a hash
sub readInInfo {
	my ($infoFiles_ref, $delimiter, $noHeader, $allInfo_ref) = @_;

	foreach my $file (@$infoFiles_ref) {
		my $infoFH;
		if (! open $infoFH, "<$file") {
			print "Skipped file $file: $!\n";
			next;
		}

		if (! $noHeader) {
			my $aLine = <$infoFH>;
			chomp($aLine);
			$allInfo_ref->{'header'} = $aLine;
		}
		
		while (my $aLine = <$infoFH>) {
			my ($id, $rest) = $aLine =~ m/^(.+?)$delimiter(.+)$/;
			print "infoID: $id -> $rest\n" if $DEBUG;
			$allInfo_ref->{$id} = $rest;
		}
	}
}


# This subroutine can print out header (when $outputType='header') or data (when anything else)
sub printBlast {
	my ($outputType, $outFH, $dm, $hit, $hsp) = @_;

	if ($outputType eq 'header') {
		# output headers this header should match with the fields in the printFields function
		print $outFH 
			"Query String", $dm,
			"Hit String", $dm,
			"Homology String", $dm,
			"Hit Accession Num", $dm,
			"Hit Name", $dm,
			"Hit length", $dm,
			"Hit description", $dm,
			"Score (bits)", $dm,
			"Raw score", $dm,
			"E value", $dm,
			"HSP total length", $dm,
			"Num of identical", $dm,
			"% of identical", $dm,
			"Num of conserved", $dm,
			"Gaps", $dm,
			"Query frame", $dm,
			"Hit frame";
	} else {
		print $outFH
			$hsp->query_string, $dm,
			$hsp->hit_string, $dm,
			$hsp->homology_string, $dm,
		  "=hyperlink(\"http://www.ncbi.nlm.nih.gov/protein/", $hit->accession, "\",\"", $hit->accession, "\")", $dm,
			$hit->name, $dm,
			$hit->length, $dm,
			$hit->description, $dm,
			$hit->bits, $dm,
			$hit->raw_score, $dm,
			$hit->significance, $dm,
			$hsp->hsp_length, $dm,
			$hsp->num_identical, $dm,
			$hsp->percent_identity, $dm,
			$hsp->num_conserved, $dm,
			$hsp->gaps, $dm,
			$hsp->query->frame, $dm,
			$hsp->hit->frame;
	}
}




########## options ##########

GetOptions(
		'help|?' => \$help, 
		'man' => \$man,
		'debug' =>\$DEBUG,
		'verbose' =>\$verbose,
		'seqFileFormat=s' =>\$seqFileFormat,
		'inDelimiter=s' =>\$inDelimiter,
		'outDelimiter=s' =>\$outDelimiter,
		'noOutHeader' =>\$noOutHeader,
		'infoNoHeader' =>\$infoNoHeader,
		'onlyHits' =>\$onlyHits,
		'numOfHSP=i' =>\$numOfHSP,
		'eval=f' =>\$eval,
		'seq=s' =>\@seqFiles,
		'info=s' =>\@infoFiles,
		'out=s' =>\$outFile,
	) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(1) if ( !$outFile );
pod2usage(1) if ($seqFileFormat ne 'fasta' && $seqFileFormat ne 'fastq' );
pod2usage(-verbose => 2) if $man;


########## arguments ##########
pod2usage("***** Error: require at least one blast result file as input argument *****") if (@ARGV < 1 );


########## program start ##########
printTime("Parse blast results starts at") if ($verbose);

print Dumper(\@ARGV) if $DEBUG;

&readInSeq(\@seqFiles, $seqFileFormat, \%allSeq);
&readInInfo(\@infoFiles, $inDelimiter, $infoNoHeader, \%allInfo);

open my $outFH, ">$outFile" or die "Can't open for write, $!";


# output headers 
if (!$noOutHeader) {
	# print info header
	if (%allInfo && !$infoNoHeader) {
		$allInfo{'header'} =~ s/$inDelimiter/$outDelimiter/g if ($inDelimiter ne $outDelimiter);
		print $outFH $allInfo{'header'},$outDelimiter;
	} else {
		print $outFH "ReadID",$outDelimiter;
	}

	if (%allSeq) {
		print $outFH "Query Sequence",$outDelimiter;
	}
	# print blast header
	&printBlast('header', $outFH, $outDelimiter);
	print $outFH "\n";
}


foreach my $blastFile (@ARGV) {
	my $in = new Bio::SearchIO(-format => 'blast', 
														 -file   => $blastFile);
	my $flag;
	my $sequence;
	while( my $result = $in->next_result ) {

		$flag=0;
		## $result is a Bio::Search::Result::ResultI compliant object
		while( my $hit = $result->next_hit ) {
			## $hit is a Bio::Search::Hit::HitI compliant object

			print "hit significance: ". $hit->significance. "  require e-value: " .  $eval . "\n" if $DEBUG;
			# Only print value if the E-value is less than or equal to the requested one
			if ($hit->significance <= $eval) {
				print "sig less than equal to eval.\n" if $DEBUG;
				my $numOfHspPrint=0;
				while( my $hsp = $hit->next_hsp ) {
					if ($numOfHspPrint <= $numOfHSP) {

						# ReadID
						print $outFH $result->query_name . $outDelimiter;

						if (%allInfo) {
							# print info data
							$allInfo{$result->query_name} =~ s/$inDelimiter/$outDelimiter/g if ($inDelimiter ne $outDelimiter);
							print $outFH $allInfo{$result->query_name}.$outDelimiter;
						}

						if (%allSeq) {
							$sequence = $allSeq{$result->query_name};
							print $outFH $sequence,$outDelimiter;
						}

						## $hsp is a Bio::Search::HSP::HSPI compliant object
						# print blast result data
						&printBlast('data', $outFH, $outDelimiter, $hit, $hsp );
						$numOfHspPrint++;
						$flag=1;
					}
				}  
			}
		}

		# if no hits, no blast result data but to print a simple 'NoHit' word.
		if ($flag==0 && !$onlyHits && $result->query_name) {

			# ReadID
			print $outFH $result->query_name, $outDelimiter;

			if (%allInfo && $allInfo{$result->query_name}) {
				# print info data
				$allInfo{$result->query_name} =~ s/$inDelimiter/$outDelimiter/g if ($inDelimiter ne $outDelimiter);
				print $outFH $allInfo{$result->query_name}, $outDelimiter;
			}

			if (%allSeq && $allSeq{$result->query_name}) {
				$sequence = $allSeq{$result->query_name};
				print $outFH $sequence, $outDelimiter;
			}
			print $outFH "NoHit";
		}

		# print a new line for next record
		print $outFH "\n" if (!($flag==0 && $onlyHits));
	}
}

close($outFH);

printTime("Parse blast results ends at") if ($verbose);


__DATA__

=head1 NAME

parseBlast.pl - This script parse a blast result file to output into a tabular file.
										Sequence and information files can be included to be in the output file.


=head1 SYNOPSIS

parseBlast.pl [-help] [-man] [-debug] [-verbose] [-seqFileFormat <fasta|fastq>]
[-inDelimiter <Delimiter>] [-outDelimiter <Delimiter>] [-noOutHeader] [-infoNoHeader]
[-onlyHits] [-numOfHSP <INTEGER>] [-eval <E-value>] [-seq <SequenceFile1> [-seq <SequenceFile2>] ... ]
[-info <InformationFile1> [-info <InformationFile2>] ... ] -out <OutputFile>
<BlastResultFile1> [<BlastResultFile2> ... ]

B<Example>

perl parseBlast.pl -eval 1e-5 -seq C72-LTE_allSeq_uniqRead.fasta -info C72-LTE_allSeq_uniqRead_DER_info.fasta -out C72-LTE_allSeq_uniqRead_DER.tsv C72-LTE_allSeq_uniqRead_DER_seq.blastx


=head1 OPTIONS

=over 4

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-verbose>

Prints progress.

=item B<-debug>

Turn on the script debugging printout.

=item B<-seqFileFormat>

Accept values: 

  fasta - input sequence files are fasta format. (Default)

  fastq - input sequence files are fastq format

=item B<-inDelimiter>

Delimiter used in the info input file. Default is a tab '\t'.

=item B<-outDelimiter>

Delimiter used in the output file. Default is a tab '\t'.

=item B<-noOutHeader>

This flag will turn off outputing header to the output file. 

=item B<-infoNoHeader>

If no header line is in the Information files, include this flag.

=item B<-onlyHits>

Including this flag will only output reads with hits.

=item B<-numOfHSP>

Control the number of high-scoring segment pairs (HSP) to be printed 
when there are multiple HSP in the blast result. Default is 1, only print the best hit.

=item B<-eval>

Take in a numerical value for only outputting reads having the E-value less than or equal to this value. 
Default is infinity meaning outputing reads with any e-values. If -onlyHits is set, reads with e-value larger than
this e-value will be marked as 'NoHit'.

=item B<-seq>

The sequence files. It accepts multiple files.

=item B<-info>

The extra information files. It accepts multiple files.

=item B<-out>

The output file. This filename is required.

=back


=head1 DESCRIPTION

This script parse a blast result file to output into a tabular file.
Sequence and information files can be included to be in the output file.

=head1 LICENSE

Copyright (c) 2012 Chon-Kit Kenneth Chan. All rights reserved.

This script is freely available for non-commercial use 
and is covered by the GNU General Public License v3 or later. 
See L<http://www.gnu.org/licenses/gpl.html>

The user of the script agrees to acknowledge the author(s) in any 
scientific publication of results based in part on the use of the script.

=cut

