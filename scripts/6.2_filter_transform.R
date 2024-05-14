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
needed_packs <- c(  "tidyverse", "plyr", "cli", "argparser", "this.path", "patchwork", "rio", "vegan") 
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
print(script_dir)
working_dir <- getinitwd() # gets dir from where the script is launched
source(paste0(script_dir, "/data_analysis_functions.R"))
#source(paste0(script_dir, "/data_analysis_options.RData"))

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
##### Load data #####
taxons <- c("phylum", "class", "order", "family", "genus", "species")
count_table <- rio::import(argv$counttable) %>% column_to_rownames("V1")
taxonomy <- rio::import(argv$taxonomy) %>% column_to_rownames("V1")
metadata <- rio::import(argv$metadata) %>% column_to_rownames("V1")

##### Create new counttable containing tax_info, and replace NA with "unassigned" #####
counttable <- rownames_to_column(taxonomy, "ASV") %>% 
    left_join(rownames_to_column(as.data.frame(count_table), "ASV"), by = "ASV") %>% 
    replace(is.na(.), "unassigned") %>% 
    column_to_rownames("ASV") %>% 
    as.data.frame()

counttable <- ReplaceUnassigned(counttable)

#### MAKE INPUT DATA IN LIST FORMAT (for each taxon) ####
##### Subset and collapse for each taxon level #####
## taxon_data_list should not put all unassigned species together for example. It should take the lowest defined taxon, instead of the species.
taxon_data_list <- MakeTaxonDataList(counttable, taxons)
cli_alert_success("Making a count table for each taxon level (order, family, genus, species)")

##### Filter data #####
cli_h1("Filtering data (removing unassigned and unidentified taxa) and transforming into relative abundances")
filtered_list <- FilterAbs(taxon_data_list, metadata, taxons, Nreads = 2000, Outdir = outdir)
cli_alert_success("Successfully filtered out unidentified and unassigned taxa")

##### Transform data #####
transformed_list <- Transform(filtered_list, taxons, Filter = "TotalAbundance")
cli_alert_success("Successfully transformed absolute to relative abundances (and removing taxa that represent < 0.01% of abundances)")

##### Save to .csv as well as .RData #####
cli_alert_info("Writing filtered and transformed tables")
l_ply(taxons, .fun = function(x) write.table(filtered_list[[x]], file = paste0(outdir, "/filtered_", x, ".csv"), quote = FALSE, sep = ";", row.names = TRUE, col.names = TRUE))
l_ply(taxons, .fun = function(x) write.table(transformed_list[[x]], file = paste0(outdir, "/transformed_", x, ".csv"), quote = FALSE, sep = ";", row.names = TRUE, col.names = TRUE))
saveRDS(filtered_list, file = paste0(outdir, "/filtered_list.RData"))
saveRDS(transformed_list, file = paste0(outdir, "/transformed_list.RData"))