#!/bin/bash

# m4_ignore(
echo "This is just a script template, not the script (yet) - pass it to 'argbash' to fix this." >&2
exit 11  #)Created by argbash-init v2.10.0
# ARG_HELP([<Making the directory branching for the analysis and creating symbolic links from the read files in the input directory, to the 0_raw_reads directory.>])
# ARG_POSITIONAL_SINGLE([run-dir], [Path to the directory containing the reads from the sequencing run.])
# ARGBASH_GO

# [ <-- needed because of Argbash

# vvv  PLACE YOUR CODE HERE  vvv
### Make directories
mkdir -p  0_raw_reads
mkdir -p  1_fastqc
mkdir -p  2_manifest
mkdir -p  3_analysis/3.1_cutadapt
mkdir -p  3_analysis/3.2_trimming
mkdir -p  3_analysis/3.3_taxonomy
mkdir -p  3_analysis/3.4_data_analysis

### Get arguments
RUNDIR=$_arg_run_dir

for FASTQ in $(find $(readlink -f $RUNDIR) -name "*.fastq.gz")
    do ABBREV=`echo $FASTQ | rev | cut -d"/" -f1 | rev`
    ln -s $FASTQ 0_raw_reads/$ABBREV
done

# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
