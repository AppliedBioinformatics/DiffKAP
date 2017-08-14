#!/bin/bash
################################################################
# Description: This script outputs a summary of the DiffKAP results
# 
# Author: Chon-Kit Kenneth Chan
# Date: 2012
# 
# Copyright (c) 2012 Chon-Kit Kenneth Chan. All rights reserved.
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


dp=2; # num of decimal point
minDiff=-1;
minFoldChange=-1;
VERBOSE=0;

##### Handle arguments
usage="$0 -a <T1-ID> -b <T2-ID> -o <OUTPUT-DIR> -f <fastq|fasta> \
-k <kmer-SIZE> -x <T1-DATA-DIR> -y <T2-DATA-DIR> [-p <NUM-OF-DECIMAL-POINT-TO-DISPLAY>] \
[-d <MIN-DIFFERENCE>] [-c <MIN-FOLD-CHANGE>] [-e <EVALUE>] [-n <NUM-OF-LINE-PER-READ>] [-v] \
-m <NUM-DER-TO-REPORT>";

while getopts "a:b:o:f:k:x:y:m:d:p:c:e:n:v" OPTION
do
  case $OPTION in
    a) T1_ID=${OPTARG};;
    b) T2_ID=${OPTARG};;
    o) myOutDir=${OPTARG};;
    f) DATA_FORMAT=${OPTARG};;
    k) kmerSize=${OPTARG};;
    x) DATA1_DIR=${OPTARG};;
    y) DATA2_DIR=${OPTARG};;
    p) dp=${OPTARG};;
    m) DER2REPORT=${OPTARG};;
    n) numLinePerRead=${OPTARG};;
    d) minDiff=${OPTARG};;
    c) minFoldChange=${OPTARG};;
    e) evalue=${OPTARG};;
    v) VERBOSE=1;;
    *) echo "$usage"; exit 1;; 
  esac
done


if test $VERBOSE -eq 1
then
	echo "Generate summary report starts at: `date`";
fi


if [ -z $T1_ID ] || [ -z $T2_ID ] || [ -z $myOutDir ] || [ -z $DATA_FORMAT ] || [ -z $kmerSize ] || [ -z $DATA1_DIR ] || [ -z $DATA2_DIR ] || [ -z $DER2REPORT ]
then
	echo "USAGE - $usage"; exit 1;
fi


KMER_DIR="${myOutDir}/kmerCount";
DEK_DIR="${myOutDir}/DEK";
DER_DIR="${myOutDir}/DER";
RESULT_DIR="${myOutDir}/results";

T1_DATA_DIR_NAME=`basename $DATA1_DIR`;
outDump1="${KMER_DIR}/${T1_DATA_DIR_NAME}_${kmerSize}mer_count.tsv"
T2_DATA_DIR_NAME=`basename $DATA2_DIR`;
outDump2="${KMER_DIR}/${T2_DATA_DIR_NAME}_${kmerSize}mer_count.tsv"

UNIQREAD_FILES=`ls ${DER_DIR}/${T1_ID}-${T2_ID}_allSeq_uniqRead.fasta.*`;
DERALL_FILE="${DER_DIR}/${T1_ID}-${T2_ID}_AllDER.tsv";
DER1_FILE="${RESULT_DIR}/${T1_ID}-${T2_ID}_AllDER_4-Highly-expressed-in-${T1_ID}_neat.tsv";
DER2_FILE="${RESULT_DIR}/${T1_ID}-${T2_ID}_AllDER_2-Highly-expressed-in-${T2_ID}_neat.tsv";

ANNO_DERALL="${DER_DIR}/${T1_ID}-${T2_ID}_AnnotatedDER.tsv";
ANNO_DER1="${RESULT_DIR}/${T1_ID}-${T2_ID}_AnnotatedDER_4-Highly-expressed-in-${T1_ID}_neat.tsv";
ANNO_DER2="${RESULT_DIR}/${T1_ID}-${T2_ID}_AnnotatedDER_2-Highly-expressed-in-${T2_ID}_neat.tsv";

GENE_FILE="${RESULT_DIR}/${T1_ID}-${T2_ID}_DEG_minDER${DER2REPORT}.tsv";

### The summary output file
SUMMARY_FILE="${RESULT_DIR}/${T1_ID}-${T2_ID}_DiffKAP_summary.log";




if [ "${DATA_FORMAT}" = "fasta" ] ; then
	SYMB1='^>';
	SYMB2='^>';
else
	SYMB1='^@';
	SYMB2='^+';
fi


# if knowing the expected number of line per read, use the line count method, 
#   otherwise just count the starting symbol that maybe over counted as '@' may occur in the beginning of the quality line
if [ ${numLinePerRead} ]; then 
	numOfLineT1=`grep -vP '^\s*$' ${DATA1_DIR}/* | wc -l`;
	numOfReadT1=$(( $numOfLineT1 / $numLinePerRead ));
	numOfLineT2=`grep -vP '^\s*$' ${DATA2_DIR}/* | wc -l`;
	numOfReadT2=$(( $numOfLineT2 / $numLinePerRead ));
else
	numOfReadT1=`grep -P ${SYMB1} ${DATA1_DIR}/* | wc -l`;
	numOfReadT2=`grep -P ${SYMB1} ${DATA2_DIR}/* | wc -l`;
fi

echo -e "# of read in ${T1_ID}\t${numOfReadT1}" | tee $SUMMARY_FILE;
echo -e "# of read in ${T2_ID}\t${numOfReadT2}" | tee -a $SUMMARY_FILE;



numOfReadT1T2=`expr $numOfReadT1 + $numOfReadT2`;
echo -e "# of read in ${T1_ID} & ${T2_ID}\t$numOfReadT1T2" | tee -a $SUMMARY_FILE;
numOfUniqRead=`grep -P '^>' $UNIQREAD_FILES | wc -l`;
echo -e "# of uniq read in ${T1_ID} & ${T2_ID}\t$numOfUniqRead" | tee -a $SUMMARY_FILE;
if [ ${numOfReadT1T2} -eq 0 ]
then
	pcOfUniqRead='NAN';
else
	pcOfUniqRead=`echo "scale=$dp; ${numOfUniqRead}*100/${numOfReadT1T2}" | bc -l`
fi
echo -e "% of uniq read in ${T1_ID} & ${T2_ID}\t${pcOfUniqRead}%" | tee -a $SUMMARY_FILE;

# use sed to get the first sequnce, use grep to list all letters in lines, then count the lines 
readLenT1=`sed -ne "/${SYMB1}/,/${SYMB2}/ {/${SYMB1}/! {/${SYMB2}/! p}}" -e "/${SYMB2}/ q" ${DATA1_DIR}/* | grep -o '\w' | wc -l`
echo -e "Read length in ${T1_ID}\t${readLenT1}" | tee -a $SUMMARY_FILE;
readLenT2=`sed -ne "/${SYMB1}/,/${SYMB2}/ {/${SYMB1}/! {/${SYMB2}/! p}}" -e "/${SYMB2}/ q" ${DATA2_DIR}/* | grep -o '\w' | wc -l`
echo -e "Read length in ${T2_ID}\t${readLenT2}" | tee -a $SUMMARY_FILE;


# kmer info
echo -e "Kmer size used\t${kmerSize}" | tee -a $SUMMARY_FILE;



if [ ${minDiff} -ne -1 ]; then 
	echo -e "Min Difference used\t${minDiff}" | tee -a $SUMMARY_FILE;
fi

isSetFoldChange=`echo "${minFoldChange} != -1" | bc`;

#if [ ${minFoldChange} -ne -1 ]; then 
if [ ${isSetFoldChange} ]; then 
	echo -e "Min Differential Fold Change used\t${minFoldChange}" | tee -a $SUMMARY_FILE;
fi



numOfKmerT1=`awk '{sum+=$2} END {print sum}' ${outDump1}`;
echo -e "Total # of kmer in ${T1_ID}\t${numOfKmerT1}" | tee -a $SUMMARY_FILE;
distKmerT1=`cat ${outDump1} | wc -l`;
echo -e "# of distinct kmer in ${T1_ID}\t${distKmerT1}" | tee -a $SUMMARY_FILE;
if [ ${numOfKmerT1} -eq 0 ]
then
	pcOfDistKmerT1='NAN';
else
	pcOfDistKmerT1=`echo "scale=$dp; ${distKmerT1}*100/${numOfKmerT1}" | bc -l`;
fi
echo -e "% of distinct kmer in ${T1_ID}\t${pcOfDistKmerT1}%" | tee -a $SUMMARY_FILE;

numOfKmerT2=`cat ${outDump2} | awk '{sum+=$2} END {print sum}'`;
echo -e "Total # of kmer in ${T2_ID}\t${numOfKmerT2}" | tee -a $SUMMARY_FILE;
distKmerT2=`cat ${outDump2} | wc -l`;
echo -e "# of distinct kmer in ${T2_ID}\t${distKmerT2}" | tee -a $SUMMARY_FILE;
if [ ${numOfKmerT2} -eq 0 ]
then
	pcOfDistKmerT2='NAN';
else
	pcOfDistKmerT2=`echo "scale=$dp; ${distKmerT2}*100/${numOfKmerT2}" | bc -l`;
fi
echo -e "% of distinct kmer in ${T2_ID}\t${pcOfDistKmerT2}%" | tee -a $SUMMARY_FILE;


# DEK info
numOfDEK=`cat ${DEK_DIR}/* | wc -l`;
echo -e "# of DEK\t${numOfDEK}" | tee -a $SUMMARY_FILE;
if [ ${distKmerT1} -eq 0 ]
then
	pcOfDEKT1='NAN';
else
	pcOfDEKT1=`echo "scale=$dp; ${numOfDEK}*100/${distKmerT1}" | bc -l`;
fi
echo -e "% of DEK to distinct kmer in ${T1_ID}\t${pcOfDEKT1}%" | tee -a $SUMMARY_FILE;
if [ ${distKmerT2} -eq 0 ]
then
	pcOfDEKT2='NAN';
else
	pcOfDEKT2=`echo "scale=$dp; ${numOfDEK}*100/${distKmerT2}" | bc -l`;
fi
echo -e "% of DEK to distinct kmer in ${T2_ID}\t${pcOfDEKT2}%" | tee -a $SUMMARY_FILE;


# DER info
numOfDER=`cat ${DERALL_FILE} | awk 'BEGIN {count=0} !/^[[:space:]]*$/ {count++} END {print count-1}'`;
echo -e "# of DER\t${numOfDER}" | tee -a $SUMMARY_FILE;
if [ ${numOfUniqRead} -eq 0 ]
then
	pcOfDER='NAN';
else
	pcOfDER=`echo "scale=$dp; ${numOfDER}*100/${numOfUniqRead}" | bc -l`;
fi
echo -e "% of DER to uniq read\t${pcOfDER}%" | tee -a $SUMMARY_FILE;
numOfDERHET1=`cat ${DER1_FILE} | awk 'BEGIN {count=0} {count++} END {print count-1}'`;
echo -e "# of DER highly expressed in ${T1_ID}\t${numOfDERHET1}" | tee -a $SUMMARY_FILE;
numOfDERHET2=`cat ${DER2_FILE} | awk 'BEGIN {count=0} {count++} END {print count-1}'`;
echo -e "# of DER highly expressed in ${T2_ID}\t${numOfDERHET2}" | tee -a $SUMMARY_FILE;


if [ ${evalue} ]; then 
	echo -e "E-value used for annotation\t${evalue}" | tee -a $SUMMARY_FILE;
fi

numOfAnnoDER=`cat ${ANNO_DERALL} | awk 'BEGIN {count=0} {count++} END {print count-1}'`;
echo -e "# of annotated DER\t${numOfAnnoDER}" | tee -a $SUMMARY_FILE;
if [ ${numOfDER} -eq 0 ]
then
	pcOfAnnoDER='NAN';
else
	pcOfAnnoDER=`echo "scale=$dp; ${numOfAnnoDER}*100/${numOfDER}" | bc -l`;
fi
echo -e "% of annotated DER\t${pcOfAnnoDER}%" | tee -a $SUMMARY_FILE;
numOfAnnoDERHET1=`cat ${ANNO_DER1} | awk 'BEGIN {count=0} {count++} END {print count-1}'`;
echo -e "# of annotated DER highly expressed in ${T1_ID}\t${numOfAnnoDERHET1}" | tee -a $SUMMARY_FILE;
numOfAnnoDERHET2=`cat ${ANNO_DER2} | awk 'BEGIN {count=0} {count++} END {print count-1}'`;
echo -e "# of annotated DER highly expressed in ${T2_ID}\t${numOfAnnoDERHET2}" | tee -a $SUMMARY_FILE;


# DEG info
numOfDEG=`cat ${GENE_FILE} | awk 'BEGIN {count=0} {count++} END {print count-1}'`;
echo -e "Total # of DEG\t${numOfDEG}" | tee -a $SUMMARY_FILE;

numOfDEGlt10DER=`cat ${GENE_FILE} | awk 'BEGIN {count=0; FS="\t"} {if ($8 < 10) count++} END {print count}'`;
# The title description is counted as greater than 10 so need to minus 1
numOfDEGge10DER=`cat ${GENE_FILE} | awk 'BEGIN {count=0; FS="\t"} {if ($8 >= 10) count++} END {print count-1}'`;
if [ ${numOfDEG} -eq 0 ]
then 
	pcOfDEGlt10DER='NAN';
	pcOfDEGge10DER='NAN'; 
else 
	pcOfDEGlt10DER=`echo "scale=$dp; ${numOfDEGlt10DER}*100/${numOfDEG}" | bc -l`;
	pcOfDEGge10DER=`echo "scale=$dp; ${numOfDEGge10DER}*100/${numOfDEG}" | bc -l`; 
fi

echo -e "# of DEG with less than 10 DER\t${numOfDEGlt10DER}" | tee -a $SUMMARY_FILE;
echo -e "% of DEG with less than 10 DER\t${pcOfDEGlt10DER}%" | tee -a $SUMMARY_FILE;

echo -e "# of DEG with 10 or more DER\t${numOfDEGge10DER}" | tee -a $SUMMARY_FILE;
echo -e "% of DEG with 10 or more DER\t${pcOfDEGge10DER}%" | tee -a $SUMMARY_FILE;

if test $VERBOSE -eq 1
then
	echo "Generate summary report ends at: `date`";
fi
