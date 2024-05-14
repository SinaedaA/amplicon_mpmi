#!/bin/bash

# m4_ignore(
echo "This is just a script template, not the script (yet) - pass it to 'argbash' to fix this." >&2
exit 11  #)Created by argbash-init v2.10.0
# ARG_OPTIONAL_SINGLE([db2], [d], [Path to second database (if classifier is "both", then --db2 should be the one for IDTAXA)])
# ARG_OPTIONAL_SINGLE([outpath], [o], [Path to output directory], [3_analysis/3.3_taxonomy/])
# ARG_OPTIONAL_SINGLE([asv2seq], [a], [Path to asv2seq table (output from 4_dada2.R)], [3_analysis/3.2_trimming/dada2_flash2/asv2seq.tsv])
# ARG_POSITIONAL_SINGLE([fasta], [Path to ASV.fasta file])
# ARG_POSITIONAL_SINGLE([seqtab], [Path to seqtab_nochim.tsv as 4_dada2.R output])
# ARG_POSITIONAL_SINGLE([db], [Path to first database (IDTAXA Rdata or RDP database for vsearch)])
# ARG_POSITIONAL_SINGLE([classifier], [Which classifier to use, either "idtaxa", "vsearch", "both"])
# ARG_DEFAULTS_POS
# ARG_HELP([A bash script that will execute either IDTAXA or RDP classifier, two taxonomic classifiers for amplicon data.])
# ARGBASH_GO

# [ <-- needed because of Argbash

# vvv  PLACE YOUR CODE HERE  vvv
### Get arguments
FASTA=$_arg_fasta
OUTDIR=$_arg_outpath
ASV2SEQ=$_arg_asv2seq
DB=$_arg_db
DB2=$_arg_db2
SEQTAB=$_arg_seqtab
CLASS=$_arg_classifier
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [ ! -d ./$OUTDIR ]; then
  mkdir -p ./$OUTDIR;
fi
mkdir -p $OUTDIR/tmp/

### Check database and software presence
## check if vsearch is installed (but not if classifier == idtaxa)
if [ $CLASS != "idtaxa" ]; then
    command -v vsearch >/dev/null 2>&1 || { echo -e >&2 "I require vsearch but it's not installed. \nTry installing it with homebrew: 'brew install vsearch'.  Aborting."; exit 1; }
fi
## check whether both database exist it CLASS = "both", or whether $DB exists if CLASS = anything else
if [ $CLASS == "both" ]; then
    [ -f $DB ] || { echo "$DB doesn't exist. Aborting"; }
    [ -f $DB2 ] || { echo "$DB2 doesn't exist. Aborting"; }
else
    [ -f $DB ] || { echo "$DB doesn't exist. Aborting"; }
fi

### Execute Rscript for IDTAXA if CLASS == "idtaxa"
TIMESTAMP=$(date '+%Y%m%d-%H%M%S')
if [ $CLASS == 'idtaxa' ]; then
    Rscript $SCRIPT_DIR/5_taxonomy.R $SEQTAB $ASV2SEQ $OUTDIR $DB 2>&1 | tee logfiles/5_taxonomy_idtaxa_${TIMESTAMP}.log
elif [ $CLASS == 'vsearch' ]; then
    vsearch --sintax $FASTA --db $DB --tabbedout $OUTDIR/ASV.fa.tax --sintax_cutoff 0.6 2>&1 | tee logfiles/5_taxonomy_vsearch_${TIMESTAMP}.log
    cut -f 1,4 $OUTDIR/ASV.fa.tax | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > $OUTDIR/taxonomy_temp.txt
    awk 'BEGIN{OFS=FS="\t"} {delete a; a["k"]="unassigned";a["p"]="unassigned";a["c"]="unassigned";a["o"]="unassigned";a["f"]="unassigned";a["g"]="unassigned";a["s"]="unassigned"; split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' $OUTDIR/taxonomy_temp.txt | gsed '1 i #ASV\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies' | sed 's/#//g;s/ //g' > $OUTDIR/taxonomy_vsearch.tsv
elif [ $CLASS == 'both' ]; then
    vsearch --sintax $FASTA --db $DB --tabbedout $OUTDIR/ASV.fa.tax --sintax_cutoff 0.6 2>&1 | tee logfiles/5_taxonomy_vsearch_${TIMESTAMP}.log
    cut -f 1,4 $OUTDIR/ASV.fa.tax | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > $OUTDIR/taxonomy_temp.txt
    awk 'BEGIN{OFS=FS="\t"} {delete a; a["k"]="unassigned";a["p"]="unassigned";a["c"]="unassigned";a["o"]="unassigned";a["f"]="unassigned";a["g"]="unassigned";a["s"]="unassigned"; split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' $OUTDIR/taxonomy_temp.txt | gsed '1 i #ASV\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies' | sed 's/#//g;s/ //g' > $OUTDIR/taxonomy_vsearch.tsv
    Rscript $SCRIPT_DIR/5_taxonomy.R $SEQTAB $ASV2SEQ $OUTDIR $DB2 --taxonomy2 $OUTDIR/taxonomy_vsearch.tsv 2>&1 | tee logfiles/5_taxonomy_idtaxa_${TIMESTAMP}.log
fi

# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
