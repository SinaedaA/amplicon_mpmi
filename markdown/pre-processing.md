# PRE-PROCESSING THE SEQUENCING DATA
## What do we get from the sequencing platform?
The files we get from the sequencing platform are in FASTQ format, which is a variant of the FASTA format commonly used for sequence data. The FASTQ format has the particularity that each base encoded in the file is associated with a quality score as well, reflecting the "certainty" with which the sequencer has called that specific base.         
![Fasta and Fastq file difference explanation.](./images/sequencing_data.png)     

As you can see on the above image, both FASTA and FASTQ files are simple text files which follow certain rules. You can open them using a simple text editor, or in the command line you can view what is inside using the `cat` or `head` command (to see the whole file, or the first 10 lines, respectively). 

```bash
head examples/example_R1.fastq.gz
```

This will return what seems like a random string of characters. That's because the file is compressed, using `gzip`. We can either unzip the file first (`gunzip`), or use the `zcat` command:

```bash
zcat examples/example_R1.fastq.gz | head
gunzip examples/example_R1.fastq.gz
## Now we can check the content of the file:
head examples/example_R1.fastq
```

To save space on the computer, we will re-zip the files. The programs that are used in bioinformatics are often able to deal with compressed files, so that will not be a problem. 

```bash
gzip examples/example_R1.fastq
```

### Exercise:
Consider a fasta file, which follows the rules:
```
>header_135
SEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCE
SEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCE
SEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCE
>header_2235
SEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCE
SEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCE
SEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCE
>header_42
SEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCE
SEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCE
SEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCE
SEQUENCESEQUENCE
...
```

How would you count the number of sequences in this file? Using the `grep -c` command, what pattern can you use to count the number of sequences?

```bash

```

What about a FASTQ file? 

```bash

```

## Create a directories
Creating a directory can easily be done in the command line using the `mkdir` command, which stands for "make directory".       
```bash
## example
mkdir 0_raw_reads
```

For the purpose of clarity, we will number our directories to reflect the order in which the analysis is done. If you follow `mkdir` with several directory names, it will create all of them at the same time.         

```bash
mkdir 0_raw_reads 1_fastqc 2_cutadapt 3_dada2 4_taxonomy
```

As you can see, the pre-processing of our amplicon sequencing data has 4 steps:         
1. **Quality check**: using the `FastQC` software, we will examine the quality of our reads.      
2. **Removing adapters**: using `cutadapt` we will remove the illumina adapters, and the primers that were used in the amplification.      
3. **Denoising**:        
4. **Assigning taxonomy**:      

## Downloading our data


## 1. Quality check
Of course, we don't want to look at the quality of the data with our naked eyes. For this, we will use the FastQC software, which will create a quality report for each of our `.fastq` files.      
Because we have to execute that on each of our files, we are going to use a `for` loop, which will automatically loop over our files and create the associated reports. 

```bash

```