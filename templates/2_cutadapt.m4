#!/bin/bash

# m4_ignore(
echo "This is just a script template, not the script (yet) - pass it to 'argbash' to fix this." >&2
exit 11  #)Created by argbash-init v2.10.0
# ARG_POSITIONAL_SINGLE([primer_file], [File containing the primers used for the amplification, forward and reverse separated by a space])
# ARG_POSITIONAL_SINGLE([run_directory], [Directory in which the fastq files from sequencing can be found])
# ARG_POSITIONAL_SINGLE([path_tg], [Path to trim_galore program, or just name of program to call if installed through bioconda])
# ARG_POSITIONAL_SINGLE([metadata], [Path to metadata file])
# ARG_OPTIONAL_SINGLE([length], [l], [Minimum length of reads (default: 50)], [50])
# ARG_OPTIONAL_SINGLE([outdir], [o], [Directory for the output of cutadapt], [3_analysis/3.1_cutadapt])
# ARG_DEFAULTS_POS
# ARG_HELP([Script to trim the primers from the sequencing reads.])
# ARGBASH_GO

# [ <-- needed because of Argbash

# vvv  PLACE YOUR CODE HERE  vvv
# For example:
PRIMERS=$_arg_primer_file
RUNDIR=$_arg_run_directory
TRIM_GAL=$_arg_path_tg
LENGTH=$_arg_length
OUTDIR=$_arg_outdir
METADATA=$_arg_metadata
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

### Get primers
FWD=`cat $PRIMERS | cut -d" " -f1`
REV=`cat $PRIMERS | cut -d" " -f2`
### Get there reverse-complements
FWD_RC=`python3 $SCRIPT_DIR/rev_complement.py $FWD`
REV_RC=`python3 $SCRIPT_DIR/rev_complement.py $REV`

if [ ! -d $OUTDIR ]; then
  mkdir -p $OUTDIR;
fi
mkdir -p $OUTDIR/trim_galore/fastqc/tmp/
mkdir -p $OUTDIR/untrimmed/
mkdir -p $OUTDIR/tooshort/
mkdir -p $OUTDIR/fastqc/tmp/

echo "Forward primer: $FWD"
echo "RC of Forward: $FWD_RC"
echo "Reverse primer: $REV"
echo "RC of Reverse: $REV_RC"

## Check program is installed, and if not propose way to install it
command -v cutadapt >/dev/null 2>&1 || { echo -e >&2 "I require cutadapt but it's not installed in your environment. \nTry installing it with conda: 'conda install -c bioconda cutadapt'.  Aborting."; exit 1; }
command -v fastqc >/dev/null 2>&1 || { echo -e >&2 "I require fastqc but it's not installed in your environment. \nTry installing it with conda: 'conda install -c bioconda fastqc'.  Aborting."; exit 1; }
command -v multiqc >/dev/null 2>&1 || { echo -e >&2 "I require fastqc but it's not installed in your environment. \nTry installing it with conda: 'conda install -c bioconda multiqc'.  Aborting."; exit 1; }
command -v $TRIM_GAL >/dev/null 2>&1 || { echo -e >&2 "I require trim_galore but the path is not correct, or it is not installed in your environment. \nTry installing it with conda: 'conda install trim-galore'.  Aborting."; exit 1; }

IFS=$'\n'

## Run trim_galore first to remove the adapters
echo "Running trim_galore to remove sequencing ADAPTERS"
# if I want to run fastqc after each step
# --fastqc --fastqc_args "--outpath $OUTDIR/trim_galore/fastqc --dir $OUTDIR/trim_galore/tmp/ --extract -t 15" 
$TRIM_GAL --quality 0 --phred64 --nextera --trim-n --fastqc --fastqc_args "--outdir $OUTDIR/trim_galore/fastqc --dir $OUTDIR/trim_galore/fastqc/tmp/ --extract -t 15" --output_dir $OUTDIR/trim_galore/ --paired $RUNDIR/*fastq.gz

## Run cutadapt in a for loop
echo "Starting cutadapt on each sample to remove PRIMERS (skipping first line, as the header)"
for line in `sed 1d $METADATA`; do 
    SAMPLE=`echo $line | cut -f1`
    R1=`echo $OUTDIR/trim_galore/${SAMPLE}_R1_001_val_1.fq.gz`
    R2=`echo $OUTDIR/trim_galore/${SAMPLE}_R2_001_val_2.fq.gz`
    # cutadapt -g "Fwd_primer=^$FWD;max_error_rate=0.1...Rev_RC=$REV_RC;max_error_rate=0;rightmost" \
    #         -G "Rev_primer=^$REV;max_error_rate=0.1...Fwd_RC=$FWD_RC;max_error_rate=0;rightmost" \
    #         --report minimal #if want only minimal report
    cutadapt -g "$FWD;max_error_rate=0.1" -a "$REV_RC;max_error_rate=0" \
            -G "$REV;max_error_rate=0.1" -A "$FWD_RC;max_error_rate=0" \
            --minimum-length $LENGTH \ 
            --match-read-wildcards \
            --too-short-output $OUTDIR/tooshort/${SAMPLE}_R1_tooshort.fastq.gz \
            --too-short-paired-output $OUTDIR/tooshort/${SAMPLE}_R2_tooshort.fastq.gz \
            --untrimmed-output $OUTDIR/untrimmed/${SAMPLE}_R1_untrimmed.fastq.gz \
            --untrimmed-paired-output $OUTDIR/untrimmed/${SAMPLE}_R2_untrimmed.fastq.gz \
            -o $OUTDIR/${SAMPLE}_R1_001.fastq.gz \
            --paired-output $OUTDIR/${SAMPLE}_R2_001.fastq.gz \
            --info-file $OUTDIR/cutadapt_infofile.tsv \
            $R1 \ 
            $R2
    fastqc --outdir $OUTDIR/fastqc/ --dir $OUTDIR/fastqc/tmp/ --extract -t 15 $OUTDIR/$SAMPLE*R1* $OUTDIR/$SAMPLE*R2*
    echo 'Quality analysis (FASTQC) finished for sample ' $SAMPLE
done

multiqc $OUTDIR -o $OUTDIR
mv $OUTDIR/multiqc_report.html $OUTDIR/2_cutadapt_multiqc_report.html

# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
