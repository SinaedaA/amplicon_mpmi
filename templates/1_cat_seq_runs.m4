#!/bin/bash

# m4_ignore(
echo "This is just a script template, not the script (yet) - pass it to 'argbash' to fix this." >&2
exit 11  #)Created by argbash-init v2.10.0
# ARG_OPTIONAL_SINGLE([outpath], [o], [Complete path to output directory for concatenated runs (default: current ./)], [./])
# ARG_POSITIONAL_DOUBLEDASH()
# ARG_POSITIONAL_SINGLE([Metadata], [Metadata file for the samples])
# ARG_POSITIONAL_INF([RunDir], [Directories containing sequencing runs to concatenate], [2])
# ARG_DEFAULTS_POS
# ARG_HELP([<Help message>])
# ARGBASH_GO

# [ <-- needed because of Argbash

# vvv  PLACE YOUR CODE HERE  vvv
# Script that concatenates sequencing runs when there are several.
# Written by Stefania Concetta Quattro
# Modified by Sinaeda Anderssen

### Get arguments
RUNDIRS=${_arg_rundir[@]}
OUTDIR=$_arg_outpath
SAMPLE_METADATA=$_arg_metadata

### Communicate about your work
echo '#############'
echo "Run directories:"
perl -E 'say join "\n", @ARGV' ${RUNDIRS[@]}
echo 'Output directory of merged libraries is ' $OUTDIR
echo 'Working on samples included in the file ' $SAMPLE_METADATA

if [ ! -d ./$OUTDIR ]; then
  mkdir -p ./$OUTDIR;
fi
MAINDIR=$PWD
cd $OUTDIR

### Loop through samples and concatenate corresponding ones from different runs
for sample in $(cat ../$SAMPLE_METADATA | cut -f1)
do
	echo $sample
	mkdir $sample
    cd $sample
	echo "Working directory $PWD"
    cat `printf "$MAINDIR/%s/$sample-*/*_R1_001.fastq.gz " ${RUNDIRS[@]}` > ${sample}_R1_001.fastq.gz
    cat `printf "$MAINDIR/%s/$sample-*/*_R2_001.fastq.gz " ${RUNDIRS[@]}` > ${sample}_R2_001.fastq.gz
	mkdir fastqc
	fastqc --outdir fastqc --extract -t 15 *R1_001.fastq.gz *R2_001.fastq.gz
	cd ../
	echo 'Quality analysis (FASTQC) finished for sample' $sample
done

### Make symlinks to fastqc/ directory and run multiqc on fastqc/
mkdir fastqc/
cd fastqc/
for f in `cat ../sample-metadata.tsv | cut -f1`; do for qcfile in `ls ../1_cat_reads/$f/fastqc/`; do ln -s ../1_cat_reads/$f/fastqc/$qcfile ./$qcfile; done; done
cd -
multiqc fastqc/

# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
