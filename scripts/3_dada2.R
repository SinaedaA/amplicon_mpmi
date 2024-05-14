#!/usr/bin/env Rscript

#### ENV SET-UP AND ARGUMENTS ####
##### Loading packages silently #####
# Check if one of .libPaths() is writable for the user, and if not, create a writable package library in $HOME (~)
source("~/Rfunctions/check_libPaths.R")
writable_path <- check_libPaths(verbose = FALSE)
# Install/load pacman.
suppressPackageStartupMessages(if (!require(pacman)) {
    install.packages("pacman", lib = writable_path)
})
needed_packs <- c("BiocManager", "tidyverse", "plyr", "ggplot2", "ShortRead", "Biostrings", "cli", "argparser", "tibble", "this.path") # , "here"
github_packs <- c("benjjneb/dada2")
to_install <- needed_packs[!needed_packs %in% pacman::p_library()]
message("Installing and/or loading required packages...")
if (length(to_install) != 0) {
    pacman::p_install(to_install, character.only = TRUE, try.bioconductor = TRUE, path = writable_path)
}
pacman::p_load(needed_packs, character.only = TRUE)
pacman::p_load_gh(github_packs)

##### Getting script directory to know where we want to source our other functions #####
# get_script_dir <- function() {
#     path <- commandArgs() %>%
#         tibble::enframe(name = NULL) %>%
#         tidyr::separate(
#             col = value, into = c("key", "value"), sep = "=", fill = "right"
#         ) %>%
#         dplyr::filter(key == "--file") %>%
#             dplyr::pull(value)
#     gsub("4_dada2.R", "", path)
# }
# script_dir <- get_script_dir()
script_dir <- this.dir()
working_dir <- getinitwd()

##### Parse arguments with arg_parser #####
cli_h1("Parsing input")
parser <- arg_parser("DADA2 denoising for Amplicon Sequencing data")
parser <- add_argument(parser, "input_directory", help = "Absolute path to directory where the reads are located (output of 3_cutadapt.sh).")
parser <- add_argument(parser, "output_directory", help = "Absolute path to directory where the output should be written.")
# parser <- add_argument(parser, "metadata", help = "Metadata for the group of interest (Bacteria, Fungi or Oomycetes)")
parser <- add_argument(parser, "--join", default = "dada2", help = "Which join method to use, either dada2 or flash2. If 'flash2', reads will be merged using this software, and denoising will be done with dada2 in single-end mode.")
parser <- add_argument(parser, "--bin", default = '', help = "Path to your local install of flash2 (if conda install just activate the environment beforehand).")
parser <- add_argument(parser, "--flashout", default = "flash2_output", help = "Name of directory to put flash2 output files .")
parser <- add_argument(parser, "--overlap", default = 12, help = "Length of overlap (in bp) for read merging.")
parser <- add_argument(parser, "--primerfile", help = "Path to the primer file, in order to check if they were actually removed using cutadapt.")
parser <- add_argument(parser, "--trimLeft",
    default = 0,
    help = "The number of nucleotides to remove from the start of each read. If both truncLen and trimLeft are provided, filtered reads will have length truncLen-trimLeft. [default: no trimming]"
)
parser <- add_argument(parser, "--trimRight",
    default = 0, 
    help = "The number of nucleotides to remove from the end of each read. If both truncLen and trimRight are provided, truncation will be performed after trimRight is enforced. [default: no trimming]"
)
parser <- add_argument(parser, "--truncLen", default = 0, help = "Truncate reads after truncLen bases. Reads shorter than this are discarded. [default: no truncation]")
parser <- add_argument(parser, "--truncQ", help = "Truncate reads at the first instance of a quality score less than or equal to truncQ.")
parser <- add_argument(parser, "--minLen", default = 50, help = "Minimum length of reads for them to be considered. This is considered AFTER trimming and truncating.")
parser <- add_argument(parser, "--maxEE", default = Inf, help = "After truncation, reads with higher than maxEE 'expected errors' will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10)). [default: no EE filtering]")
argv <- parse_args(parser)

##### Sourcing functions and defining inpath #####
source(paste0(script_dir, "/dada2_functions.R"))
primers <- setNames(read.table(argv$primerfile, sep = " ", header = F), nm = c("FOR", "REV"))
cli_alert_info("FORWARD primer : {primers$FOR}, REVERSE primer : {primers$REV}")
inpath <- paste0(working_dir, "/", argv$input_directory)

##### Setting path variable, with the output directory #####
## create the directory if it doesn't exist yet
outpath <- paste0(working_dir, "/", argv$output_directory)
if (!dir.exists(outpath)) {
    dir.create(outpath, recursive = TRUE)
    cli_alert_info("Output directory {outpath} doesn't exist, creating it.")
} else {
    cli_alert_info("Output directory {outpath} already exists.")
}

##### Getting .fastq files from the inpath directory #####
## Generate list of FWD and REV filenames
cli_alert_info("Looking for fastq files")
filenames_F <- sort(list.files(inpath, pattern = "_R1_001.fastq.gz", full.names = TRUE))
filenames_R <- sort(list.files(inpath, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# i wanted to try it on the trim_galore output just to check if I actually found the primers, but I don't
# filenames_F <- sort(list.files(inpath, pattern = "_R1_001_val_1.fq.gz", full.names = TRUE))
# filenames_R <- sort(list.files(inpath, pattern = "_R2_001_val_2.fq.gz", full.names = TRUE))

##### Read primer file and extract primer sequences (plus RCs) #####
cli_h1("Checking if primers are present in reads.")
FOR_orients <- allOrients(primers$FOR)
REV_orients <- allOrients(primers$REV)
print(FOR_orients)
print(REV_orients)

##### Pre-filter reads: put the output (w/o Ns) into "filtN" sub-directory #####
cli_alert_info("Removing reads containing Ns for this purpose...")
filenames_filtN_F <- file.path(inpath, "filtN", basename(filenames_F))
filenames_filtN_R <- file.path(inpath, "filtN", basename(filenames_R))
filterAndTrim(filenames_F, filenames_filtN_F, filenames_R, filenames_filtN_R, maxN = 0, multithread = TRUE)

##### Update which filtN names actually were written, some having no reads pass the filter ######
filtN_path <- file.path(inpath, "filtN")
filenames_filtN_F <- list.files(filtN_path, full.names = TRUE, pattern = "R1")
filenames_filtN_R <- list.files(filtN_path, full.names = TRUE, pattern = "R2")

##### Count how many times our primers appear in our sequences (should be 0) #####
cli_alert_info("Making primer matching table...")
primers_after_cutadapt <- rbind(
    FOR_ForwardReads = sapply(FOR_orients, primerHits, fn = filenames_filtN_F[[1]]),
    FOR_ReverseReads = sapply(FOR_orients, primerHits, fn = filenames_filtN_R[[1]]),
    REV_ForwardReads = sapply(REV_orients, primerHits, fn = filenames_filtN_F[[1]]),
    REV_ReverseReads = sapply(REV_orients, primerHits, fn = filenames_filtN_R[[1]])
)
print(primers_after_cutadapt)
cli_alert_info("Writing this tables to {outpath}/primer_presence_after_cutadapt.csv")
write.table(primers_after_cutadapt, file = paste0(outpath, "/primer_presence_after_cutadapt.csv"), quote = FALSE, sep = "\t")

##### Perform flash2 merge before if join == flash2 #####
if (argv$join == 'flash2') {
    flash2_outdir <- argv$flashout
    flash2 <- paste0(argv$bin, "flash2")
    cli_h1("Using {flash2} to merge reads before denoising with dada2...")
    if (Sys.which(flash2) == "") {
        cli_alert_danger("I require flash2 but it's not installed in your environment or I can't find it.")
        cli_alert_danger("If it is installed, please provide the absolute path to the directory containing the executable, \nor activate the conda environment.")
        cli_abort("Try installing it with conda: 'conda install -c bioconda flash2' or from source 'https://github.com/dstreett/FLASH2/tree/master'.  Aborting.")
    }
    else {
        cli_alert_success("flash2 is installed, proceeding...")
    }
    cli_h1("Running flash2 to merge reads.")
    flash_count <- 0
    for (i in seq_along(filenames_F)) {
        base_fn <- str_split(basename(filenames_F[i]), pattern = "_R", n = 2)[[1]][1]
        print(basename(filenames_F[i]))
        print(base_fn)
        if (!file.exists(paste0(outpath, "/", base_fn, ".extendedFrags.fastq.gz"))) {
            system2(flash2, args = c(
                "--phred-offset", 33,
                "--min-overlap", argv$overlap,
                "--max-overlap", 300,
                "--compress", filenames_F[i], filenames_R[i],
                "--output-directory", outpath,
                "--output-prefix", base_fn
            ))
        } else {
            flash_count <- flash_count + 1
        }
    }
    cli_alert_warning("{flash_count} files were alread present in {flash2_outdir} directory, so skipped them.")
    cli_alert_info("If you want to re-run flash2 on all samples, please erase previous run files, or specify other output directory in --flashout option.")

    ## Then execute dada2 but single-end, and on the base_fn.extendedFrags.fastq.gz
    cli_h1("Running dada2 to denoise merged reads.")
    filenames <- sort(list.files(outpath, pattern = "extendedFrags.fastq.gz", full.names = TRUE))
    if (!file.exists(paste0(outpath, "/flash2_dada2_denoising.RDS"))){
        dada2_results <- dada2_wrap(inpath, filenames = list(filenames), argv$maxEE, argv$truncQ, argv$truncLen, argv$trimLeft, argv$trimRight, argv$minLen)
        saveRDS(dada2_results, file = paste0(outpath, "/flash2_dada2_denoising.RDS"))
    } else {
        cli_alert_warning("flash2_dada2_denoising.RDS already exists in your directory, using this data.")
        cli_alert_info("If you want to re-do the analysis and store results in another directory, change the --flashout option.")
        dada2_results <- readRDS(paste0(outpath, "/flash2_dada2_denoising.RDS"))
    }
    list2env(setNames(dada2_results, nm = c("filt_df", "filtered", "filtered_files", "dada")), envir=.GlobalEnv)
    cli_alert_success("Success, dada2 was performed on all flash2 samples")
    cli_alert_info("RDS file saved in {outpath}/flash2_dada2_denoising.RDS")
    
    ## Make sequence table and remove chimeras
    cli_h1("Removing chimeras...")
    seqtab <- makeSequenceTable(dada)
    seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
    ## Write sequence table 
    cli_h1("Writing sequence count tables (seqtab.tsv and seqtab_nochim.tsv)...")
    write.table(seqtab, file = paste0(outpath, "/seqtab.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
    write.table(seqtab_nochim, file = paste0(outpath, "/seqtab_nochim.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

    cli_alert_info("Inspecting distribution of length. Data saved to {outpath}/length_distribution.csv")
    length_distribution <- table(nchar(getSequences(seqtab_nochim)))
    write.table(length_distribution, file = paste0(outpath, "/length_distribution.csv"), quote = FALSE, row.names = TRUE, col.names = TRUE)

    ### Track reads through pipeline
    cli_alert_info("Tracking the reads through the pipeline... See 'flash2_dada2_denoising_stats.csv'")
    track <- as.data.frame(cbind(
        filt_df, sapply(dada, getN),
        rowSums(seqtab_nochim)
    )) %>% setNames(., nm = c("before_filter", "after_filter", "denoised", "non_chimeric"))
    
    stat_plot <- draw_stat_plot(denoising_stats = track)
    stat_table <- make_stat_table(track, merging = FALSE)
    ggsave(stat_plot, file = paste0(outpath, "/dada2_flash2_denoising_stats.pdf"))
    write.table(stat_table, file = paste0(outpath, "/dada2_flash2_denoising_percentages.csv"), quote = FALSE, row.names = TRUE, col.names = TRUE)
    write.table(track, file = paste0(outpath, "/flash2_dada2_denoising_stats.csv"), quote = FALSE, row.names = TRUE, col.names = TRUE)

} else if (argv$join == "dada2") {
    ##### Perform dada2 denoising first, then dada2 merging #####
    cli_h1("Performing dada2 denoising on paired-end reads...")
    if (!file.exists(paste0(outpath,"dada2_denoising.RDS"))){
        dada2_results <- dada2_wrap(inpath, filenames = list(filenames_F, filenames_R), argv$maxEE, argv$truncQ, argv$truncLen, argv$trimLeft, argv$trimRight, argv$minLen)
        saveRDS(dada2_results, file = paste0(outpath, "/dada2_denoising.RDS"))
    } else {
        cli_alert_warning("dada2_denoising.RDS already exists in your directory, using this data.")
        cli_alert_warning("If you want to re-do the analysis and store results in another directory, change the output directory argument.")
        dada2_results <- readRDS(paste0(outpath, "/dada2_denoising.RDS"))
    }
    list2env(setNames(dada2_results, nm = c("filt_df", "filtered", "filtered_F", "filtered_R", "dada_F", "dada_R")), envir = .GlobalEnv)
    
    ## Merge paired reads
    cli_h1("Merging pairs of reads with dada2 (mergePairs function)...")
    merged <- mergePairs(dada_F, filtered_F, dada_R, filtered_R, verbose = TRUE)

    ## Construct sequence table, then remove chimeras
    cli_h1("Removing chimeras...")
    seqtab <- makeSequenceTable(merged)
    seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
    write.table(seqtab, file = paste0(outpath, "/seqtab.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
    write.table(seqtab_nochim, file = paste0(outpath, "/seqtab_nochim.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
    cli_alert_info("Inspecting distribution of length. Data saved to {outpath}/length_distribution.csv")
    length_distribution <- table(nchar(getSequences(seqtab_nochim)))
    write.table(length_distribution, file = paste0(outpath, "/length_distribution.csv"), quote = FALSE, row.names = TRUE, col.names = TRUE)
    # Inspect distribution of length
    print("Inspect distribution of length")
    #print(table(nchar(getSequences(seqtab_nochim))))
    ### Track reads through pipeline
    getN <- function(x) sum(getUniques(x))
    track <- as.data.frame(cbind(
        filt_df, sapply(dada_F, getN), sapply(dada_R, getN), sapply(merged, getN),
        rowSums(seqtab_nochim)
    )) %>% setNames(., nm = c("before_filter", "after_filter", "denoised_F", "denoised_R", "merged", "non_chimeric"))
    print("Read tracking")
    #print(track)
    write.table(track, file = paste0(outpath, "/dada2_denoising_stats.csv"), quote = FALSE, row.names = TRUE, col.names = TRUE)

    stat_plot <- draw_stat_plot(denoising_stats = track)
    stat_table <- make_stat_table(track)
    ggsave(stat_plot, file = paste0(outpath, "/dada2_denoising_stats.pdf"))
    write.table(stat_table, file = paste0(outpath, "/dada2_denoising_percentages.csv"), quote = FALSE, row.names = TRUE, col.names = TRUE)
}
## write fasta file with unique ASVs
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

## make ASV table with total counts for each
# seqtab_nochim is : SEQ x FILE (SEQ = columns, FILE = sample_filename)
#seqtab_nochim <- dplyr::select(seqtab_nochim, -c("V1"))
asv_total <- as.data.frame(colSums(seqtab_nochim)) 
asv2seq <- dplyr::mutate(asv_total, ASV = paste0("ASV", seq(nrow(asv_total)))) %>%
    rownames_to_column("Seq") %>%
    dplyr::select("ASV", "Seq")

## write fasta file based on this asv_names and asv2seq
write.table(asv2seq, file = paste0(outpath, "/asv2seq.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cli_alert_success("ASV to SEQ correspondance written to {outpath}/asv2seq.tsv : ASV x SEQ")
writeFasta(setNames(as_tibble(asv2seq), nm = c("name", "seq")), paste0(outpath, "ASV.fasta"))
cli_alert_success("ASV fasta file written to {outpath}/ASV.fasta : >ASVxxx \n sequence")
