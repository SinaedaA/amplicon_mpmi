#!/bin/bash

# Created by argbash-init v2.10.0
# ARG_OPTIONAL_SINGLE([db2],[d],[Path to second database (if classifier is "both", then --db2 should be the one for IDTAXA)])
# ARG_OPTIONAL_SINGLE([outpath],[o],[Path to output directory],[3_analysis/3.3_taxonomy/])
# ARG_OPTIONAL_SINGLE([asv2seq],[a],[Path to asv2seq table (output from 4_dada2.R)],[3_analysis/3.2_trimming/dada2_flash2/asv2seq.tsv])
# ARG_OPTIONAL_SINGLE([id],[i],[ID of the column 1 of asv2seq],["ASV"])
# ARG_POSITIONAL_SINGLE([fasta],[Path to ASV.fasta file])
# ARG_POSITIONAL_SINGLE([seqtab],[Path to seqtab_nochim.tsv as 4_dada2.R output])
# ARG_POSITIONAL_SINGLE([db],[Path to first database (IDTAXA Rdata or RDP database for vsearch)])
# ARG_POSITIONAL_SINGLE([classifier],[Which classifier to use, either "idtaxa", "vsearch", "both"])
# ARG_DEFAULTS_POS()
# ARG_HELP([A bash script that will execute either IDTAXA or RDP classifier, two taxonomic classifiers for amplicon data.])
# ARGBASH_GO()
# needed because of Argbash --> m4_ignore([
### START OF CODE GENERATED BY Argbash v2.10.0 one line above ###
# Argbash is a bash code generator used to get arguments parsing right.
# Argbash is FREE SOFTWARE, see https://argbash.io for more info


die()
{
	local _ret="${2:-1}"
	test "${_PRINT_HELP:-no}" = yes && print_help >&2
	echo "$1" >&2
	exit "${_ret}"
}


begins_with_short_option()
{
	local first_option all_short_options='doaih'
	first_option="${1:0:1}"
	test "$all_short_options" = "${all_short_options/$first_option/}" && return 1 || return 0
}

# THE DEFAULTS INITIALIZATION - POSITIONALS
_positionals=()
_arg_fasta=
_arg_seqtab=
_arg_db=
_arg_classifier=
# THE DEFAULTS INITIALIZATION - OPTIONALS
_arg_db2=
_arg_outpath="3_analysis/3.3_taxonomy/"
_arg_asv2seq="3_analysis/3.2_trimming/dada2_flash2/asv2seq.tsv"
_arg_id="ASV"


print_help()
{
	printf '%s\n' "A bash script that will execute either IDTAXA or RDP classifier, two taxonomic classifiers for amplicon data."
	printf 'Usage: %s [-d|--db2 <arg>] [-o|--outpath <arg>] [-a|--asv2seq <arg>] [-i|--id <arg>] [-h|--help] <fasta> <seqtab> <db> <classifier>\n' "$0"
	printf '\t%s\n' "<fasta>: Path to ASV.fasta file"
	printf '\t%s\n' "<seqtab>: Path to seqtab_nochim.tsv as 4_dada2.R output"
	printf '\t%s\n' "<db>: Path to first database (IDTAXA Rdata or RDP database for vsearch)"
	printf '\t%s\n' "<classifier>: Which classifier to use, either \"idtaxa\", \"vsearch\", \"both\""
	printf '\t%s\n' "-d, --db2: Path to second database (if classifier is \"both\", then --db2 should be the one for IDTAXA) (no default)"
	printf '\t%s\n' "-o, --outpath: Path to output directory (default: '3_analysis/3.3_taxonomy/')"
	printf '\t%s\n' "-a, --asv2seq: Path to asv2seq table (output from 4_dada2.R) (default: '3_analysis/3.2_trimming/dada2_flash2/asv2seq.tsv')"
	printf '\t%s\n' "-i, --id: ID of the column 1 of asv2seq (default: '"ASV"')"
	printf '\t%s\n' "-h, --help: Prints help"
}


parse_commandline()
{
	_positionals_count=0
	while test $# -gt 0
	do
		_key="$1"
		case "$_key" in
			-d|--db2)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_db2="$2"
				shift
				;;
			--db2=*)
				_arg_db2="${_key##--db2=}"
				;;
			-d*)
				_arg_db2="${_key##-d}"
				;;
			-o|--outpath)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_outpath="$2"
				shift
				;;
			--outpath=*)
				_arg_outpath="${_key##--outpath=}"
				;;
			-o*)
				_arg_outpath="${_key##-o}"
				;;
			-a|--asv2seq)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_asv2seq="$2"
				shift
				;;
			--asv2seq=*)
				_arg_asv2seq="${_key##--asv2seq=}"
				;;
			-a*)
				_arg_asv2seq="${_key##-a}"
				;;
			-i|--id)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_id="$2"
				shift
				;;
			--id=*)
				_arg_id="${_key##--id=}"
				;;
			-i*)
				_arg_id="${_key##-i}"
				;;
			-h|--help)
				print_help
				exit 0
				;;
			-h*)
				print_help
				exit 0
				;;
			*)
				_last_positional="$1"
				_positionals+=("$_last_positional")
				_positionals_count=$((_positionals_count + 1))
				;;
		esac
		shift
	done
}


handle_passed_args_count()
{
	local _required_args_string="'fasta', 'seqtab', 'db' and 'classifier'"
	test "${_positionals_count}" -ge 4 || _PRINT_HELP=yes die "FATAL ERROR: Not enough positional arguments - we require exactly 4 (namely: $_required_args_string), but got only ${_positionals_count}." 1
	test "${_positionals_count}" -le 4 || _PRINT_HELP=yes die "FATAL ERROR: There were spurious positional arguments --- we expect exactly 4 (namely: $_required_args_string), but got ${_positionals_count} (the last one was: '${_last_positional}')." 1
}


assign_positional_args()
{
	local _positional_name _shift_for=$1
	_positional_names="_arg_fasta _arg_seqtab _arg_db _arg_classifier "

	shift "$_shift_for"
	for _positional_name in ${_positional_names}
	do
		test $# -gt 0 || break
		eval "$_positional_name=\${1}" || die "Error during argument parsing, possibly an Argbash bug." 1
		shift
	done
}

parse_commandline "$@"
handle_passed_args_count
assign_positional_args 1 "${_positionals[@]}"

# OTHER STUFF GENERATED BY Argbash

### END OF CODE GENERATED BY Argbash (sortof) ### ])
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
ID=$_arg_id
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
    Rscript $SCRIPT_DIR/5_taxonomy.R $SEQTAB $ASV2SEQ $OUTDIR $DB --id $ID 2>&1 | tee logfiles/5_taxonomy_idtaxa_${TIMESTAMP}.log
elif [ $CLASS == 'vsearch' ]; then
    vsearch --sintax $FASTA --db $DB --tabbedout $OUTDIR/ASV.fa.tax --sintax_cutoff 0.6 2>&1 | tee logfiles/5_taxonomy_vsearch_${TIMESTAMP}.log
    cut -f 1,4 $OUTDIR/ASV.fa.tax | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > $OUTDIR/taxonomy_temp.txt
    awk 'BEGIN{OFS=FS="\t"} {delete a; a["k"]="unassigned";a["p"]="unassigned";a["c"]="unassigned";a["o"]="unassigned";a["f"]="unassigned";a["g"]="unassigned";a["s"]="unassigned"; split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' $OUTDIR/taxonomy_temp.txt | gsed '1 i #ASV\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies' | sed 's/#//g;s/ //g' > $OUTDIR/taxonomy_vsearch.tsv
elif [ $CLASS == 'both' ]; then
    vsearch --sintax $FASTA --db $DB --tabbedout $OUTDIR/ASV.fa.tax --sintax_cutoff 0.6 2>&1 | tee logfiles/5_taxonomy_vsearch_${TIMESTAMP}.log
    cut -f 1,4 $OUTDIR/ASV.fa.tax | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > $OUTDIR/taxonomy_temp.txt
    awk 'BEGIN{OFS=FS="\t"} {delete a; a["k"]="unassigned";a["p"]="unassigned";a["c"]="unassigned";a["o"]="unassigned";a["f"]="unassigned";a["g"]="unassigned";a["s"]="unassigned"; split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' $OUTDIR/taxonomy_temp.txt | gsed '1 i #ASV\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies' | sed 's/#//g;s/ //g' > $OUTDIR/taxonomy_vsearch.tsv
    Rscript $SCRIPT_DIR/5_taxonomy.R $SEQTAB $ASV2SEQ $OUTDIR $DB2 --taxonomy2 $OUTDIR/taxonomy_vsearch.tsv --id $ID 2>&1 | tee logfiles/5_taxonomy_idtaxa_${TIMESTAMP}.log
fi

# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
