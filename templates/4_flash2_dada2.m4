#!/bin/bash

# m4_ignore(
echo "This is just a script template, not the script (yet) - pass it to 'argbash' to fix this." >&2
exit 11  #)Created by argbash-init v2.10.0
# ARG_HELP([<Uses first flash2 to join reads, followed by dada2 to denoise in paired-end mode>])
# ARG_OPTIONAL_SINGLE([])
# ARG_POSITIONAL_SINGLE([metadata])
# ARG_POSITIONAL_SINGLE([outdir])
# ARGBASH_GO

### Parsing input arguments ###

### Performing FLASH2 on data ###

echo "Using FLASH2 to join reads, followed by DADA2 to denoise in single-end mode"
mkdir -p $OUTDIR/FLASH2_DADA2/TrimmedSequences/
mkdir -p $OUTDIR/FLASH2_DADA2/JoinedReads/
mkdir -p $OUTDIR/FLASH2_DADA2/DenoisedReads/
command -v flash2 >/dev/null 2>&1 || { echo -e >&2 "I require flash2 but it's not installed in your environment. \nTry installing it with conda: 'conda install -c bioconda flash2' or from source 'https://github.com/dstreett/FLASH2/tree/master'.  Aborting."; exit 1; }
for i in `cat sample-metadata.tsv | cut -f1`; do
    #change Phred-offset to variable at some point
    flash2 \
        --phred-offset=33 \
        --min-overlap=$OVERLAP \
        --max-overlap=300 \
        --compress $i*R1_001.fastq.gz $i*R2_001.fastq.gz \
        --output-directory=$OUTDIR/FLASH2_DADA2/JoinedReads/ \
        --output-prefix=$i 2>&1 | tee flash.log
done

# [ <-- needed because of Argbash

# vvv  PLACE YOUR CODE HERE  vvv

# PUT YOUR CODE HERE

# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
