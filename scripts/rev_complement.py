#!/usr/bin/env python3
import sys

# Complementing DNA
# defining the DNA dictionary
dna = sys.argv[1]
compl = {
    "A" : "T",
    "T" : "A",
    "C" : "G",
    "G" : "C",
    "M" : "K",
    "K" : "M",
    "S" : "S",
    "W" : "W",
    "R" : "Y",
    "Y" : "R",
    "D" : "H",
    "V" : "B",
    "H" : "D",
    "B" : "V",
    "N" : "N"
}

rev = ""
for i in dna:
   rev = compl[i] + rev

print(rev)