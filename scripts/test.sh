#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo $SCRIPT_DIR

TO_TEST=$1
echo $TO_TEST
command -v $TO_TEST >/dev/null 2>&1 || { echo -e >&2 "I require trim-galore but it's not installed in your environment. \nTry installing it with conda: 'conda install -c bioconda multiqc'.  Aborting."; exit 1; }