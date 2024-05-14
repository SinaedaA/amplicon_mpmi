#!/bin/bash

# m4_ignore(
echo "This is just a script template, not the script (yet) - pass it to 'argbash' to fix this." >&2
exit 11  #)Created by argbash-init v2.10.0
# ARG_OPTIONAL_SINGLE([outpath], [o], [Complete path to output directory], [1_fastqc/])
# ARG_POSITIONAL_DOUBLEDASH()
# ARG_POSITIONAL_SINGLE([metadata], [Metadata file for the samples])
# ARG_POSITIONAL_SINGLE([read-dir], [Directory containing raw read files (fastq.gz)])
# ARG_DEFAULTS_POS
# ARG_HELP([Small wrapper script for FastQC and MultiQC])
# ARGBASH_GO

# [ <-- needed because of Argbash

# vvv  PLACE YOUR CODE HERE  vvv
# Script that concatenates sequencing runs when there are several.
# Written by Stefania Concetta Quattro
# Modified by Sinaeda Anderssen

### Get arguments
READ_DIR=$_arg_read_dir
OUTDIR=$_arg_outpath
SAMPLE_METADATA=$_arg_metadata

### Communicate about your work
echo 'Output directory of FastQC analysis ' $OUTDIR
echo 'Working on samples included in the file ' $SAMPLE_METADATA

if [ ! -d ./$OUTDIR ]; then
  mkdir -p ./$OUTDIR;
fi
mkdir -p $OUTDIR/tmp/

### Loop through samples and make symlinks to RUN1 and RUN2 (to have all fastq files in same directory)
IFS=$'\n'
#for sample in $(cat $SAMPLE_METADATA | cut -f1); do
for line in `sed 1d $SAMPLE_METADATA`; do 
  sample=`echo $line | cut -f1`
  #if [[ $sample == \#* ]]; then continue; fi
  fastqc --outdir $OUTDIR --dir $OUTDIR/tmp/ --extract -t 15 $READ_DIR/$sample*R1* $READ_DIR/$sample*R2*
  echo 'Quality analysis (FASTQC) finished for sample ' $sample
done
multiqc $OUTDIR -o $OUTDIR
mv $OUTDIR/multiqc_report.html $OUTDIR/1_fastqc_multiqc_report.html


# i=1
# for READS_DIR in ${RUNDIRS[@]}; do
#     echo "Read dir is: " $READS_DIR
#     for sample in $(cat $SAMPLE_METADATA | cut -f1); do
#         if [[ $sample == \#* ]]; then continue; fi
#         R1=$(realpath $(find $READS_DIR | grep $sample"." | grep "R1" | grep "\\.fastq\\.gz"))
#         R2=$(realpath $(find $READS_DIR | grep $sample"." | grep "R2" | grep "\\.fastq\\.gz"))
#         sample_name1=`echo $R1 | rev | cut -d"/" -f1 | rev`
#         sample_name2=`echo $R2 | rev | cut -d"/" -f1 | rev`
#         ln -s $R1 $OUTDIR/RUN${i}/$sample_name1
#         ln -s $R2 $OUTDIR/RUN${i}/$sample_name2
#         fastqc --outdir $OUTDIR/RUN${i}/fastqc --dir $OUTDIR/RUN${i}/fastqc/tmp/ --extract -t 15 $OUTDIR/RUN${i}/$sample_name1 $OUTDIR/RUN${i}/$sample_name2
#         echo 'Quality analysis (FASTQC) finished for sample ' $sample ' in ' RUN${i}
#     done
#     multiqc $OUTDIR/RUN${i}/fastqc/
#     i=$(expr $i + 1)
# done

# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
