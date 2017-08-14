#!/usr/bin/env perl

package Utilities;

use Statistics::Descriptive;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(printTime getOptions loadToHash combineFields calMedianMeanCV isFaOrFq);

my $DEBUG=0;

########## functions ##########
sub printTime {
  my $message = "Current Time";
  ($message) = shift if (@_ >= 1); 

  print ("$message: " . scalar(localtime(time)) . "\n");
}


sub getOptions{
	my @std_options=("help", "man", "debug");
	my %options;
	
	GetOptions(	\%options,
		@std_options,
		@_,
	);
	
	exec("pod2usage $0") if $options{'help'};
	exec("perldoc $0") if $options{'man'};
	
	return \%options;
}


sub loadToHash{
	my ($inputFileHandle, $delimiter, $keyColumn, $keyRegExpTemplate, $printDuplicateKey) = @_;
	
	print "delimiter: $delimiter, " if $DEBUG && $delimiter;
	print "keycolumn: $keyColumn, keyRegExpTemplate: $keyRegExpTemplate\n" if $DEBUG;
	
	$delimiter="\t" if (!$delimiter);
	
	my %file;
	my $line;
	my $countDuplicate=0;
	my $countTotal=0;
	while ($line=<$inputFileHandle>) {
		print "File Line: $line\n" if $DEBUG;
		my @fields = split($delimiter, $line);
		my $curKey = parseID($fields[$keyColumn], $keyRegExpTemplate);
		print "current key: $curKey\n" if $DEBUG;
		
		splice @fields, $keyColumn, 1;
		if ($printDuplicateKey && $curKey && $file{$curKey}) {
			my @tmpFields = split($delimiter, $file{$curKey});
			$file{$curKey}=&combineFields(\@tmpFields, \@fields, $delimiter);
			$countDuplicate++;
		} else {
			$file{$curKey} = join($delimiter, @fields) if ($curKey);
		}
		$countTotal++ if ($curKey);
		
	}
	print "Number of duplicated key: $countDuplicate\n"if ($printDuplicateKey);
	print "Total number of genes: $countTotal\n";
	
	return \%file;
}



# combine two same size arrays by combining a description field and adding up all other fields.
# Input arguement: two array reference.
# return: a string of the combined fields
sub combineFields{
	my ($arrayRef1, $arrayRef2, $delimiter, $descriptionFieldCol, $descriptionCombineStr) = @_;
	
	my $desCol=0;
	my $desCombineStr="-**-";
	my $delimit="\t";
	
	$desCol=$descriptionFieldCol if ($descriptionFieldCol);
	$desCombineStr = $descriptionCombineStr if ($descriptionCombineStr);
	$delimit=$delimiter if ($delimiter);
	
	chomp $arrayRef1->[$desCol];
	chomp $arrayRef2->[$desCol];
	
	if ($arrayRef1->[$desCol] ne $arrayRef2->[$desCol]) {
		$arrayRef1->[$desCol] .= $desCombineStr.$arrayRef2->[$desCol];
	}
	
	for (my $i=1; $i<scalar(@$arrayRef1); $i++) {
		$arrayRef1->[$i] +=  $arrayRef2->[$i];
	}
	
	return join($delimit, @$arrayRef1);
}

	
## Calculate the median, mean and coefficient variance of an array of numbers
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

## input $filename
## returns ($type, $numOfLinePerRead) where $type can be 'fasta', 'fastq' or '?' for unknown.
sub isFaOrFq {                                                                                                                                               
  my ($filename) = @_; 
  my $type='';
  
  open (IN,"<$filename") || die $!; 
  
  my $charExp='';
  my $aLine;
  while ($aLine=<IN>) {
    my ($char1)=$aLine=~m/^\s*([^\s])/;
    if ($char1) {
      if ("$char1" eq "@") {
        $charExp='+';
        $type='fastq';
        last;
      } elsif ("$char1" eq ">") {
        $charExp='>';
        $type='fasta';
        last;
      } else {
        # Stop it as it's not the expected file
        return ('?', 0); 
      }   
    }   
  }
  
  my $numOfLinePerRead=0;
  my $count=0;
  while ($aLine=<IN>) {
    $count++;
    my ($char2)=$aLine=~m/^\s*([^\s])/;
    if ("$char2" eq "$charExp") {
      if ("$type" eq "fastq") {
        $numOfLinePerRead=$count*2;
      }   
      if ("$type" eq "fasta") {
        $numOfLinePerRead=$count;
      }   
      last;
    }   
  }
  
  return ($type, $numOfLinePerRead);
}



1;