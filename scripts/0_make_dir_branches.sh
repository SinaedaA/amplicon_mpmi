#!/bin/bash

# Created by argbash-init v2.10.0
# ARG_HELP([<Making the directory branching for the analysis and creating symbolic links from the read files in the input directory, to the 0_raw_reads directory.>])
# ARG_POSITIONAL_SINGLE([run-dir],[Path to the directory containing the reads from the sequencing run.])
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
	local first_option all_short_options='h'
	first_option="${1:0:1}"
	test "$all_short_options" = "${all_short_options/$first_option/}" && return 1 || return 0
}

# THE DEFAULTS INITIALIZATION - POSITIONALS
_positionals=()
# THE DEFAULTS INITIALIZATION - OPTIONALS


print_help()
{
	printf '%s\n' "<Making the directory branching for the analysis and creating symbolic links from the read files in the input directory, to the 0_raw_reads directory.>"
	printf 'Usage: %s [-h|--help] <run-dir>\n' "$0"
	printf '\t%s\n' "<run-dir>: Path to the directory containing the reads from the sequencing run."
	printf '\t%s\n' "-h, --help: Prints help"
}


parse_commandline()
{
	_positionals_count=0
	while test $# -gt 0
	do
		_key="$1"
		case "$_key" in
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
	local _required_args_string="'run-dir'"
	test "${_positionals_count}" -ge 1 || _PRINT_HELP=yes die "FATAL ERROR: Not enough positional arguments - we require exactly 1 (namely: $_required_args_string), but got only ${_positionals_count}." 1
	test "${_positionals_count}" -le 1 || _PRINT_HELP=yes die "FATAL ERROR: There were spurious positional arguments --- we expect exactly 1 (namely: $_required_args_string), but got ${_positionals_count} (the last one was: '${_last_positional}')." 1
}


assign_positional_args()
{
	local _positional_name _shift_for=$1
	_positional_names="_arg_run_dir "

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