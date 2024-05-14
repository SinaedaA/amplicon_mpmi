#!/bin/bash
# Copyright 2021 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@univr.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# m4_ignore(
echo "This is just a script template, not the script (yet) - pass it to 'argbash' to fix this." >&2
exit 11  #)Created by argbash-init v2.10.0
# ARG_OPTIONAL_SINGLE([outdir], [o], [Path to output directory (default: 2_manifest)], [2_manifest])
# ARG_POSITIONAL_SINGLE([metadata], [Path to metadata.tsv file for the sequencing experiment])
# ARG_POSITIONAL_SINGLE([read_dir], [Directories containing reads (ex: 0_raw_reads/)])
# ARG_DEFAULTS_POS
# ARG_HELP([<The general help message of my script>])
# ARGBASH_GO

# [ <-- needed because of Argbash

# vvv  PLACE YOUR CODE HERE  vvv
### Get variables
SAMPLE_METADATA=$_arg_metadata
READ_DIR=$_arg_read_dir
OUTDIR=$_arg_outdir

### Create OUTDIR directory if it doesn't exist
if [ ! -d $OUTDIR ]; then
  mkdir -p $OUTDIR;
fi

### Loop over metadata file, to get sample-names and get complete paths to directories with reads (R1 and R2)
# Check the OSTYPE, and accordingly use either realpath or grealpath (on macOS, from coreutils, from Homebrew)
if [[ $OSTYPE == "linux-gnu"* ]]; then
	echo -e sample-id"\t"forward-absolute-filepath"\t"reverse-absolute-filepath > $OUTDIR/manifest.tsv
	for sample in $(cat $SAMPLE_METADATA | cut -f1); do
		## Skip first line, starting with #
		if [[ $sample == \#* ]]; then continue; fi
		R1=$(realpath $(find $READ_DIR | grep $sample"." | grep "R1" | grep "\\.fastq\\.gz"));
		R2=$(realpath $(find $READ_DIR | grep $sample"." | grep "R2" | grep "\\.fastq\\.gz"));
		sample_name=`echo $R1 | rev | cut -d"/" -f1 | rev | cut -d"_" -f1-2`
		echo -e $sample_name"\t"$R1"\t"$R2 >> $OUTDIR/manifest.tsv
	done
elif [[ $OSTYPE == "darwin"* ]]; then
	## Checking if grealpath exists on macOS
	command -v grealpath >/dev/null 2>&1 || { echo >&2 "I require grealpath but it's not installed. Try installing 'coreutils' with.  Aborting."; exit 1; }
	gecho -e sample-id"\t"forward-absolute-filepath"\t"reverse-absolute-filepath > $OUTDIR/manifest.tsv
	for sample in $(cat $SAMPLE_METADATA | cut -f1); do
		## Skip first line, starting with #
		if [[ $sample == \#* ]]; then continue; fi
		R1=$(grealpath $(find $READ_DIR/$RUN | grep $sample"." | grep "R1" | grep "\\.fastq\\.gz"));
		R2=$(grealpath $(find $READ_DIR/$RUN | grep $sample"." | grep "R2" | grep "\\.fastq\\.gz"));
		sample_name=`echo $R1 | rev | cut -d"/" -f1 | rev | cut -d"_" -f1-2`
		echo -e $sample_name"\t"$R1"\t"$R2 >> $OUTDIR/manifest.tsv
	done
fi



# ^^^  TERMINATE YOUR CODE BEFORE THE BOTTOM ARGBASH MARKER  ^^^

# ] <-- needed because of Argbash
