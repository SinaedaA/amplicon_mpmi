#!/usr/bin/env Rscript

##### ENV SET-UP AND ARGUMENTS #####
##### Loading packages silently #####
# Check if one of .libPaths() is writable for the user, and if not, create a writable package library in $HOME (~)
source("~/Rfunctions/check_libPaths.R")
writable_path <- check_libPaths(verbose = FALSE)
# Install/load pacman.
suppressPackageStartupMessages(if (!require(pacman)) {
    install.packages("pacman", lib = writable_path)
})
needed_packs <- c(  "tidyverse", "plyr", "cli", "argparser", "this.path", "rio") 
to_install <- needed_packs[!needed_packs %in% pacman::p_library()]
message("Installing and/or loading required packages...")
if (length(to_install) != 0) {
    pacman::p_install(to_install, character.only = TRUE, try.bioconductor = TRUE, path = writable_path)
}
pacman::p_load(needed_packs, character.only = TRUE)
#github_packs <- c("benjjneb/dada2")
#pacman::p_load_gh(github_packs)

##### Set WD, source functions and variables #####
script_dir <- this.dir() # gets dir where the script is stored
working_dir <- getinitwd() # gets dir from where the script is launched
# source(paste0(script_dir, "/data_analysis_functions.R"))
# source(paste0(script_dir, "/data_analysis_options.RData"))
columns <- c("Sample_Name", "Short_Sample", "Sample_ID", "Genotype", "Compartment")

##### Parse arguments with arg_parser #####
cli_h1("Parsing input")
parser <- arg_parser("Check the compatibility of the data (count table, metadata, taxonomy file).")
parser <- add_argument(parser, "counttable", help = "Count table with absolute abundances for each taxon (5_taxonomy.R output)")
parser <- add_argument(parser, "taxonomy", help = "Taxonomy table outputted by 5_taxonomy.R")
parser <- add_argument(parser, "metadata", help = "Metadata file, where Sample_Name column contains sample names exactly as colnames in counttable")
parser <- add_argument(parser, "--outdir", help = "Output directory for all the following steps.")
argv <- parse_args(parser)

##### Communicate arguments #####
cli_alert_success("Parsed arguments successfully: ")
cli_li("count table -> {argv$counttable}")
cli_li("taxonomy -> {argv$taxonomy}")
cli_li("metadata -> {argv$metadata}")

##### Create output directory if necessary #####
outdir <- file.path(argv$outdir)
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    cli_alert_info("Output directory {outdir} doesn't exist, creating it.")
} else {
    cli_alert_info("Output directory {outdir} already exists.")
}

#### CHECK / MODIFY COUNT TABLE FORMAT ####
##### Load count table and change column names #####
cli_h1("Checking input format, and pre-processing data")
count_table <- rio::import(argv$counttable) %>% column_to_rownames("V1") %>% t()
colnames(count_table) <- str_replace(colnames(count_table), pattern = ".extended.*", "")

##### Load taxonomy table #####
taxonomy <- rio::import(argv$taxonomy, na = "unassigned") %>% 
    dplyr::select(-Seq) %>% 
    column_to_rownames("V1")
    #dplyr::filter(!is.na(domain))

##### Load metadata table #####
metadata <- rio::import(argv$metadata)[1:length(columns)] %>% 
    setNames(nm = columns)

##### Check if colnames of count table (samples) correspond partly to Sample_Name column in metadata #####
if (length(max.col(sapply(sort(colnames(count_table)), grepl, sort(metadata$Sample_Name)))) != length(colnames(count_table))){
    cli_alert_danger("The colnames of the count table (sample names) do not correspond to the Sample_Name column in metadata")
    stop()
}

cli_alert_success("All files were compatible with eachother, you're ready to go.")
rio::export(metadata, paste0(argv$outdir, "/metadata.tsv"), row.names = TRUE)
rio::export(taxonomy, paste0(argv$outdir, "/taxonomy.tsv"), row.names = TRUE)
rio::export(count_table, paste0(argv$outdir, "/count_table.tsv"), row.names = TRUE)