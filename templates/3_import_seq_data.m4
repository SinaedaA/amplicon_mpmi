#!/bin/bash

# m4_ignore(
echo "This is just a script template, not the script (yet) - pass it to 'argbash' to fix this." >&2
exit 11  #)Created by argbash-init v2.10.0
# ARG_OPTIONAL_SINGLE([type], [t], [Type of data: paired-end (pe), single-end (se) (default=pe)], [pe])
# ARG_OPTIONAL_SINGLE([format], [f], [Format of import: Manifest or Casava], [Casava])
# ARG_OPTIONAL_SINGLE([outdir], [o], [Path to output directory (default: 3_analysis)], [3_analysis])
# ARG_OPTIONAL_SINGLE([quality], [q], [Type of quality encoding: phred64 or phred33 (default=phred33)], [phred33])
# ARG_POSITIONAL_SINGLE([input-path], [Path to either directory containing sequencing files, or path to manifest.tsv file])
# ARG_DEFAULTS_POS
# ARG_HELP([<This script will simply use the manifest.tsv file to import sequencing data into QIIME2.>])
# ARGBASH_GO

# [ <-- needed because of Argbash

# vvv  PLACE YOUR CODE HERE  vvv
### Define positional and optional arguments
INPATH=$_arg_input_path
FORMAT=$_arg_format
QUALITY=$_arg_quality
OUTDIR=$_arg_outdir

### Make outdir (3_analysis) if it doesn't exist yet
if [ ! -d $OUTDIR ]; then
  mkdir -p $OUTDIR;
fi

### Define the input-format to give to QIIME, based on type and quality variables
if [[ $_arg_type == "pe" || $_arg_type == "paired-end" ]]; then
    TYPE="SampleData[PairedEndSequencesWithQuality]"
    if [[ $FORMAT == "Manifest" && $QUALITY == "phred33" ]]; then
        INPUT="PairedEndFastqManifestPhred33V2"
    elif [[ $FORMAT == "Manifest" && $QUALITY == "phred64" ]]; then
        INPUT="PairedEndFastqManifestPhred64V2"
    elif [[ $FORMAT == "Casava" ]]; then
        INPUT="CasavaOneEightSingleLanePerSampleDirFmt"
    fi
	echo "Type of data: " $TYPE >> $OUTDIR/3_import_seq_data.log
	echo "Input-format: " $INPUT >> $OUTDIR/3_import_seq_data.log
elif [[ $_arg_type == "se" || $_arg_type == "single-end" ]]
then
    TYPE="SampleData[SingleEndSequencesWithQuality]"
    if [[ $FORMAT == "Manifest" && $QUALITY == "phred33" ]]; then
        INPUT="PairedEndFastqManifestPhred33V2"
    elif [[ $FORMAT == "Manifest" && $QUALITY == "phred64" ]]; then
        INPUT="PairedEndFastqManifestPhred64V2"
    elif [[ $FORMAT == "Casava" ]]; then
        INPUT="CasavaOneEightSingleLanePerSampleDirFmt"
    fi
	echo "Type of data: " $TYPE >> $OUTDIR/3_import_seq_data.log
	echo "Input-format: " $INPUT >> $OUTDIR/3_import_seq_data.log
fi

RUNS=( RUN1 RUN2 )
for RUN in ${RUNS[@]}; do
    ### Importing data
    qiime tools import \
        --type $TYPE \
        --input-path $INPATH/$RUN \
        --output-path $OUTDIR/sequences_$RUN.qza \
        --input-format $INPUT

    qiime demux summarize \
        --i-data $OUTDIR/sequences_$RUN.qza \
        --o-visualization $OUTDIR/demux_summary_$RUN.qzv
done
# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
