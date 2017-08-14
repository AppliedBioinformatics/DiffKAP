We have developed a Differential k-mer Analysis Pipeline (DiffKAP) for the pairwise comparison 
of RNA profiles between metatranscriptomes which does not rely on mapping to reference assemblies. 
By reducing each read to component kmers and assessing the frequency of these sequences, 
we overcome statistical limitations on the lack of identical reads for pairwise comparison 
between samples and allow inference of differential gene expression for annotated reads.

The DiffKAP application consists of a series of scripts written in Perl and Linux shell scripts 
and requires Jellyfish [Marcais 2011] and BLASTx as well as access to a blast formatted protein database. 
The scripts are freely available for academic use.


===== What does DiffKAP depend on? =====

DiffKAP depends on the following things:
- Jellyfish for fast kmer counting
- blastx for sequence alignment
- Some non-standard Perl modules:
	* bioperl
	* Bio::SeqIO
	* Bio::SearchIO
	* Parallel::ForkManager
	* Statistics::Descriptive
	* Config::IniFiles
	* Log::Log4perl
	* GD::Graph::linespoints (for the script identifyKmerSize)


===== How to install? ===== 

- Download the DiffKAP package
- Uncompress it into:
	* a DiffKAP setup file
	* a README file
	* a VERSION file
	* an example data folder containing a small subset of a metatranscriptomic data
- read the README
- Install the DiffKAP setup script by running: DiffKAP_setup
*** If you like, you can add the DiffKAP path to $PATH or just use an absolute path for running DiffKAP ***


===== How to run? ===== 

1. Your project configuration file: Use the example config file in the sample data directory as a template.
2. Run the pipeline: Run DiffKAP with your config file as an input argument, for example: DiffKAP ~/sampleProj/sampleProj.cfg
* Results will be generated in the [OUT_DIR]/results where [OUT_DIR] is defined in the config file.
* The processing log is stored in /tmp/DiffKAP.log by default.


===== How to interpret the results? ===== 

- You can download the results of the sample data.
- The script "DiffKAP" generates 4 types of files in folder [OUT_DIR]/results:
	1. 5 DER files with the word 'AllDER' in the filenames. Explanation of some columns:
		* Median-T1: The median kmer occurrence represented in Treatment 1 (corresponding to T1_ID in the config file) for all kmers in the read.
		* Median-T2: Similar to Median-T1 but for Treatment 2.
		* Ratio of Median: The ratio of Median-T1 to Median-T2.
		* CV-T1: The coefficient of variation of all kmer occurrence represented in Treatment 1 for all kmers in the read. To show how confident the Median-T1 representing all kmers in the read.
		* CV-T2: Similar to CV-T1 but for Treatment 2.
	2. 5 annotated DER files with the word 'AnnotatedDER' in the filenames. These files are similar to the 5 DER above but contain only the annotated DER.
	3. A gene-centric summary with the word 'DEG' in the filename:
		* In a tabular form showing the number of DER in the specific files (in columns 3-7) annotated to the specific gene.
		* It shows the unique gene list.
		* The 'Total' column is the total number of DER annotated to such gene.
		* It is sorted by the 'Total' column in reversed order.
	4.A result summary file with 'summary.log' in the filename.


