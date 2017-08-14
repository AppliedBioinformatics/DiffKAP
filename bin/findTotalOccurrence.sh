#!/bin/bash
################################################################
# Description: This script outputs the total number of kmer occurrence
#							 of all input files and also the maximum total number.
# 
# Author: Chon-Kit Kenneth Chan
# Date: 2011
# 
# Copyright (c) 2011 Chon-Kit Kenneth Chan. All rights reserved.
# 
# This script is covered by the GNU General Public License v3 or later. 
# You can redistribute it and/or modify it under the terms of the GNU 
# General Public License v3 or later as published by the Free Software 
# Foundation. 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS.
# This program is distributed WITHOUT ANY WARRANTY; without even the 
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# The user of the script agrees to acknowledge the author(s) in any 
# scientific publication of results based in part on the use of the script.
#
# This copyright notice must be included in all distribution of the script.
################################################################

##### Handle arguments
usage="$0 -o <OUTPUT_FILE> [-v] IN_FILE1 [IN_FILE2 ...]";
OUTFILE="";
VERBOSE=0;
while getopts "o:v" OPTION
do
  case $OPTION in
    o) OUTFILE=${OPTARG};;
    v) VERBOSE=1;;
    *) echo "$usage"; exit 1;;
  esac
done

if test $VERBOSE -eq 1 
then
	echo "Find Total Occurrence starts at: `date`";
fi

# Make sure the output file is defined
if [ "$OUTFILE" = "" ]
then
	echo "$usage"; exit 1
fi

# shift away the -o OUTFILE
shift;shift

# Empty the output file
>$OUTFILE
for f in $@
do
	awk '{sum+=$2} END {print FILENAME "\t" sum}' FS="\t" $f >> ${OUTFILE}
done

echo "#####" >> ${OUTFILE}
awk 'max=="" || $2>max {max=$2} END {print "max\t" max}' FS="\t" $OUTFILE >> ${OUTFILE}

if test $VERBOSE -eq 1
then
	echo "Find Total Occurrence ends at: `date`";
fi
