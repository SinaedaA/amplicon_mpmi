#!/bin/bash

# m4_ignore(
echo "This is just a script template, not the script (yet) - pass it to 'argbash' to fix this." >&2
exit 11  #)Created by argbash-init v2.10.0
# ARG_OPTIONAL_SINGLE([join-method], [m], [Method for joining reads (default: DADA2, others: FLASH2, vsearch, deblur)], [DADA2])
# ARG_OPTIONAL_SINGLE([trim-left-f], [], [Position at which FORWARD read sequences should be trimmed due to low quality. This trims the 5' end of the input sequences, which will be the bases that were sequenced in the first cycles.], [0])
# ARG_OPTIONAL_SINGLE([trim-left-r], [], [Position at which REVERSE read sequences should be trimmed due to low quality. This trims the 5' end of the input sequences, which will be the bases that were sequenced in the first cycles.], [0])
# ARG_OPTIONAL_SINGLE([trunc-len-f], [], [Position at which FORWARD read sequences should be truncated due to decrease in quality. This truncates the 3' end of the of the input sequences, which will be the bases that were sequenced in the last cycles.], [0])
# ARG_OPTIONAL_SINGLE([trunc-len-r], [], [Position at which REVERSE read sequences should be truncated due to decrease in quality. This truncates the 3' end of the of the input sequences, which will be the bases that were sequenced in the last cycles.], [0])
# ARG_OPTIONAL_SINGLE([overlap], [o], [Mininum overlap for joining reads], [12])
# ARG_POSITIONAL_SINGLE([qza], [Path to trimmed sequences.qza file to input for joining and denoising])
# ARG_DEFAULTS_POS
# ARG_HELP([<Joining and denoising of paired end reads. Joining method can be chosen (DADA2, FLASH2, vsearch), denoising is always done with DADA2.>])
# ARGBASH_GO

# [ <-- needed because of Argbash

# vvv  PLACE YOUR CODE HERE  vvv
### Code written by Concetta Quattro and Sinaeda Anderssen
### Contact <cquattro@uni-koeln.de> <sanders2@uni-koeln.de>

### Declaring input variables
METHOD=$_arg_join_method
OVERLAP=$_arg_overlap
SEQ=$_arg_qza
TRIM_LEFT_F=$_arg_trim_left_f
TRIM_LEFT_R=$_arg_trim_left_r
TRUNC_LEN_F=$_arg_trunc_len_f
TRUNC_LEN_R=$_arg_trunc_len_r
OUTDIR="3_analysis"

### See which join method has been chosen
if [[ $METHOD -eq "DADA2" ]]; then
    ## Do denoising with DADA2 on paired-end data
    mkdir $OUTDIR/DADA2_join_denoise/
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs $SEQ \
        --p-trim-left-f $TRIM_LEFT_F --p-trim-left-r $TRIM_LEFT_R --p-trunc-len-f $TRUNC_LEN_F --p-trunc-len-r $TRUNC_LEN_R \
        --o-table $OUTDIR/DADA2_join_denoise/table.qza --o-representative-sequences $OUTDIR/DADA2_join_denoise/rep-seqs.qza --o-denoising-stats $OUTDIR/DADA2_join_denoise/denoising-stats.qza \
        --p-n-threads 0 \
        --verbose
    qiime feature-table summarize \
        --i-table $OUTDIR/DADA2_join_denoise/table.qza \
        --o-visualization $OUTDIR/DADA2_join_denoise/table.qzv \
        --m-sample-metadata-file sample-metadata.tsv
    qiime feature-table tabulate-seqs \
        --i-data $OUTDIR/DADA2_join_denoise/rep-seqs.qza \
        --o-visualization $OUTDIR/DADA2_join_denoise/rep-seqs.qzv
    qiime metadata tabulate \
        --m-input-file $OUTDIR/DADA2_join_denoise/denoising-stats.qza \
        --o-visualization $OUTDIR/DADA2_join_denoise/denoising-stats.qzv

elif [[ $METHOD -eq "FLASH2" ]]; then
    ## Use FLASH2 to join reads, THEN DADA2 to denoise (single-end mode)
    mkdir -p $OUTDIR/FLASH2_DADA2/TrimmedSequences/
    mkdir -p $OUTDIR/FLASH2_DADA2/JoinedReads/
    mkdir -p $OUTDIR/FLASH2_DADA2/DenoisedReads/
    command -v flash2 >/dev/null 2>&1 || { echo -e >&2 "I require flash2 but it's not installed in your environment. \nTry installing it with conda: 'conda install -c bioconda flash2'.  Aborting."; exit 1; }
    ### Exporting trimmed sequences
    qiime tools export \
        --input-path 3_analysis/TrimmedSequences/sequences_trimmed1.qza \
        --output-path $OUTDIR/FLASH2_DADA2/TrimmedSequences/
    ### Loop over samples to perform flash2 on R1 and R2
    for $i in `cat sample-metadata.tsv | cut -f1`; do 
        #change Phred-offset to variable at some point
        flash2 \
            --phred-offset=33 \
            --min-overlap=$OVERLAP \
            --max-overlap=300 \
            --compress $i*R1_001.fastq.gz $i*R2_001.fastq.gz \
            --output-directory=$OUTDIR/FLASH2_DADA2/JoinedReads/ \
            --output-prefix=$i 2>&1 | tee flash.log
    done
    ### Re-do steps 2 and 3 (make manifest and import data into QIIME)
    ../scripts/2_make_manifest.sh sample-metadata.tsv $OUTDIR/FLASH2_DADA2/JoinedReads/ \
        --outdir $OUTDIR/FLASH2_DADA2/2_manifest/
    ../scripts/3_import_seq_data.sh $OUTDIR/FLASH2_DADA2/2_manifest/manifest.tsv \
        --type single-end \
        --outdir $OUTDIR/FLASH2_DADA2/3_analysis/
    ### Then DADA2 denoising, single-end
    qiime dada2 denoise-single \
        --i-demultiplexed-seqs $SEQ \
        --p-trim-left-f $TRIM_LEFT_F --p-trim-left-r $TRIM_LEFT_R --p-trunc-len-f $TRUNC_LEN_F --p-trunc-len-r $TRUNC_LEN_R \
        --o-table $OUTDIR/FLASH2_DADA2/DenoisedReads/table.qza --o-representative-sequences $OUTDIR/FLASH2_DADA2/DenoisedReads/rep-seqs.qza --o-denoising-stats $OUTDIR/FLASH2_DADA2/DenoisedReads/denoising-stats.qza \
        --p-n-threads 0 \
        --verbose
    qiime feature-table summarize \
        --i-table $OUTDIR/FLASH2_DADA2/DenoisedReads/table.qza \
        --o-visualization $OUTDIR/FLASH2_DADA2/DenoisedReads/table.qzv \
        --m-sample-metadata-file sample-metadata.tsv
    qiime feature-table tabulate-seqs \
        --i-data $OUTDIR/FLASH2_DADA2/DenoisedReads/rep-seqs.qza \
        --o-visualization $OUTDIR/FLASH2_DADA2/DenoisedReads/rep-seqs.qzv
    qiime metadata tabulate \
        --m-input-file $OUTDIR/FLASH2_DADA2/DenoisedReads/denoising-stats.qza \
        --o-visualization $OUTDIR/FLASH2_DADA2/DenoisedReads/denoising-stats.qzv

elif [[ $METHOD -eq "vsearch" ]]; then
    mkdir $OUTDIR/vsearch_DADA2/
    ### join with vsearch, then quality filtering, then denoising with DADA2
else
    echo >&2 "Joining method not recognized, please try one of the following: DADA2, deblur, FLASH2, vsearch (case-sensitive)."; exit 1;

# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
