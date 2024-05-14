#!/bin/bash

# m4_ignore(
echo "This is just a script template, not the script (yet) - pass it to 'argbash' to fix this." >&2
exit 11  #)Created by argbash-init v2.10.0
# ARG_OPTIONAL_BOOLEAN([keep-untrimmed], [u], [Whether to keep untrimmed reads or not (default=off)])
# ARG_OPTIONAL_BOOLEAN([reverse-complement], [c], [Whether the reverse complement of the primers are also removed (default=off)])
# ARG_POSITIONAL_SINGLE([primer_file], [File containing the primers used for the amplification, forward and reverse separated by a space])
# ARG_DEFAULTS_POS
# ARG_HELP([<Script to trim the primers from the sequencing reads.>])
# ARGBASH_GO

# [ <-- needed because of Argbash

# vvv  PLACE YOUR CODE HERE  vvv
PRIMERS=$_arg_primer_file
KEEP_U=$_arg_keep_untrimmed
REVCOMP=$_arg_reverse_complement
### Get primers
FOR=`cat $PRIMERS | cut -d" " -f1`
REV=`cat $PRIMERS | cut -d" " -f2`
### Output directory
OUTDIR="3_analysis/TrimmedSequences"
if [ ! -d $OUTDIR ]; then
  mkdir -p $OUTDIR;
fi

RUNS=( RUN1 RUN2 )
for RUN in ${RUNS[@]}; do
    ### Actual trimming code
    ## cutadapt has "discard-untrimmed" option, so we can see activate it or not based on the input of -u parameter of script
    ## Keep untrimmed means that we want to remove option   --p-discard-untrimmed \ from cutadapt
    if [[ $KEEP_U == "off" ]]
    then
        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences 3_analysis/sequences_$RUN.qza \
            --p-cores 20 \
            --p-front-f $FOR \
            --p-front-r $REV \
            --p-discard-untrimmed \
            --verbose \
            --o-trimmed-sequences $OUTDIR/sequences_trimmed1_$RUN.qza
    elif [[ $KEEP_U == "on" ]]
    then
        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences 3_analysis/sequences_$RUN.qza \
            --p-cores 20 \
            --p-front-f $FOR \
            --p-front-r $REV \
            --verbose \
            --o-trimmed-sequences $OUTDIR/sequences_trimmed1_$RUN.qza
    fi

    ### Create summary and export
    qiime demux summarize \
        --i-data $OUTDIR/sequences_trimmed1_$RUN.qza \
        --o-visualization $OUTDIR/demux_summary_trimmed1_$RUN.qzv
    qiime tools export \
        --input-path $OUTDIR/sequences_trimmed1_$RUN.qza \
        --output-path $OUTDIR/PrimerTrimming1

    ### If we also want to remove the REVCOMP of the primers from the reads
    if [[ $REVCOMP == "on" ]]
    then
        FOR_RC=`python ../scripts/rev_complement.py $FOR`
        REV_RC=`python ../scripts/rev_complement.py $REV`
        echo $FOR_RC $REV_RC
        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $OUTDIR/sequences_trimmed1_$RUN.qza \
            --p-cores 20 \
            --p-adapter-r $FOR_RC \
            --p-adapter-f $REV_RC \
            --verbose \
            --o-trimmed-sequences $OUTDIR/sequences_trimmed2_$RUN.qza
        ### Create summary and export
        qiime demux summarize \
            --i-data $OUTDIR/sequences_trimmed2_$RUN.qza \
            --o-visualization $OUTDIR/demux_summary_trimmed2_$RUN.qzv
        qiime tools export \
            --input-path $OUTDIR/sequences_trimmed2_$RUN.qza \
            --output-path $OUTDIR/PrimerTrimming2
    fi
done
# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
