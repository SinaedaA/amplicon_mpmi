#!/usr/bin/env Rscript
## add option to create a suffix (if I want to do same analysis with different count table for example)
## add rarefaction for PCoA (for beta diversity)

#### ENV SET-UP AND ARGUMENTS ####
##### Loading packages silently #####
# Check if one of .libPaths() is writable for the user, and if not, create a writable package library in $HOME (~)
source("~/Rfunctions/check_libPaths.R")
writable_path <- check_libPaths(verbose = FALSE)
# Install/load pacman.
suppressPackageStartupMessages(if (!require(pacman)) {
    install.packages("pacman", lib = writable_path)
})
# i don't think I need patchwork. gonna remove it but we'll see
# also don't think any of these packages are from github, but we'll see
needed_packs <- c(  "tidyverse", "plyr", "ggplot2", "vegan", "gridExtra", "ggplotify", "ggpubr", "ggtext",
                    "cli", "argparser", "phyloseq", "mia", "pheatmap", "spaa", "apeglm", "this.path", "ggalluvial") 
#github_packs <- c("benjjneb/dada2")
to_install <- needed_packs[!needed_packs %in% pacman::p_library()]
message("Installing and/or loading required packages...")
if (length(to_install) != 0) {
    pacman::p_install(to_install, character.only = TRUE, try.bioconductor = TRUE, path = writable_path)
}
pacman::p_load(needed_packs, character.only = TRUE)
#pacman::p_load_gh(github_packs)

##### Set WD, source functions and variables #####
script_dir <- this.dir() # gets dir where the script is stored
working_dir <- getinitwd() # gets dir from where the script is launched
source(paste0(script_dir, "/data_analysis_functions.R"))
source(paste0(script_dir, "/data_analysis_options.RData"))
##### Parse arguments with arg_parser #####
cli_h1("Parsing input")
parser <- arg_parser("Diversity analysis of Amplicon Sequencing data.")
parser <- add_argument(parser, "counttable", help = "Count table with absolute abundances for each taxon (QIIME output)")
parser <- add_argument(parser, "metadata", help = "Metadata for the group of interest (Bacteria, Fungi or Oomycetes)")
parser <- add_argument(parser, "--target", help = "Organism of interest: bacteria, fungi or oomycetes")
parser <- add_argument(parser, "--organism", help = "Which plant to focus on: Hv (for Hordeum vulgare) or At (for Arabidopsis thaliana)")
parser <- add_argument(parser, "--outdir", default = "./data_analysis/", help = "Path to output directory.")
parser <- add_argument(parser, "--normalization", default = "TSS", help = "Method to use for normalization (total sum scaling [TSS] or rarefaction [RAR])")
parser <- add_argument(parser, "--differential-abundance", default = "None", help = "Which method to use for differential abundance of taxa (supported: deseq2, None, ANCOM-BC to come)")
argv <- parse_args(parser)

##### Communicate arguments #####
cli_alert_success("Parsed arguments successfully: ")
cli_li("metadata -> {argv$metadata}")
cli_li("count table -> {argv$counttable}")
cli_li("outdir -> {argv$outdir}")
cli_li("organism -> {argv$organism}")
cli_li("target -> {argv$target}")
cli_li("normalization -> {argv$normalization}")
cli_li("differential abundance -> {argv$differential_abundance}") 

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
freqtable <- read.table(argv$counttable, sep = "\t", header = TRUE, comment.char = "")
colnames(freqtable) <- gsub("\\.", "-", colnames(freqtable))
freqtable$taxonomy <- gsub(".__", "", freqtable$taxonomy)

##### Load metadata table #####
columns <- c("Sample_Name", "Short_Sample", "Sample_ID", "Organism", "Genotype", "Compartment", "Timepoint", "Abbrev", "Summary")
metadata <- setNames(read.table(argv$metadata, header = TRUE, sep = "\t", comment.char = ""), nm = columns)
##### Check if colnames of count table (samples) correspond exactly to Sample_Name column in metadata #####
# will need to change it so it doesn't use identical maybe, if bug in the future because some samples were removed from the data in dada2
if (!identical(sort(colnames(freqtable)[-1]), sort(metadata$Sample_Name))) {
    cli_alert_danger("The colnames of the count table (sample names) do not correspond to the Sample_Name column in metadata")
    stop()
}
##### Separate taxonomy column in count table into subdivisions #####
freqtable <- freqtable %>%
    separate(col = taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
    replace(is.na(.), "unassigned") %>%
    mutate_all(list(~str_remove_all(., "D_"))) 

##### Subset tables based on chosen organism #####
cli_alert_info("You've decided to focus your analysis on {argv$organism}, subsetting metadata and count table")
metadata <- subset(metadata, metadata$Abbrev == argv$organism)
metadata$Genotype[!metadata$Genotype == "BI1"] <- "WT"
freqtable <- freqtable %>%
    dplyr::select(1:7, metadata$Sample_Name)

#### MAKE INPUT DATA IN LIST FORMAT (for each taxon) ####
##### Subset and collapse for each taxon level #####
taxon_data_list <- MakeTaxonDataList(freqtable, taxons)
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

#### ALPHA-DIVERSITY ####
##### Compute Alpha-diversity index (Shannon) #####
cli_h1("Computing alpha-diversity for each compartment")
alpha_list <- GetAlphaList(metadata, filtered_list, taxons, index = "shannon")
cli_alert_success("Successfully calculated Shannon alpha-diversity index")
cli_alert_info("Writing alpha diversity tables.")
l_ply(taxons, .fun = function(x) write.table(alpha_list[[x]], file = paste0(outdir, "/alpha-diversity_", x, ".csv"), quote = FALSE, sep = ";", row.names = TRUE, col.names = TRUE))

##### Plot Alpha-diversity (Boxplots) #####
cli_h1("Making alpha-diversity boxplots")
alpha_plots <- AlphaPlotList(alpha_list, taxons, myColors, index = "shannon")
multi_page <- marrangeGrob(grobs = alpha_plots, ncol = 2, nrow = 2)
ggsave(multi_page, filename = paste0(outdir, "/alpha_diversity.pdf"), width = 20, height = 20)
cli_alert_success("Alpha diversity plots saved to {outdir}/alpha_diversity.pdf")


### BETA-DIVERSITY ####
#### Compute Bray-Curtis pairwise distances based on TransformedList #####
cli_h1("Computing pairwise distances between all samples")
dist_list <- getPWDistList(transformed_list, taxons, "bray")
#neg_dist_list <- getPWDistList(neg_TL, taxons, "bray")

##### Perform ordination #####
cli_h1("Performing ordination (unconstrained and constrained)")
Multi_Ordination <- OrdinateFunction(transformed_list, taxons, metadata, Metric = "bray")
l_ply(taxons, .fun = function(x) {
    ggsave(Multi_Ordination[[x]], filename = paste0(outdir, "/", x, "_ordination.pdf"), width = 25, height = 25)
    cli_alert_success("Ordination plot for {x} saved to {outdir}/{x}_ordination.pdf")
})


##### BC distance and permanovas #####
## problem: betadisper -> centroidFUN -> apply -> FUN -> tapply
PermResults <- AdonisWrapper(transformed_list, taxons, metadata, permu_scheme)
saveRDS(PermResults, paste0(outdir, "/PermanovaResults.RData"))



# #### DIFFERENTIAL ABUNDANCE ####
# cli_h1("Differential abundance analysis of taxa")
# ## with ANCOM-BC2
##### Using DESeq2 #####
# if (argv$differential_abundance == "None") {
#     cli_alert_info("No differential abundance method selected, skipping. Set the --differential-abundance argument to change this behaviour.")
# }
# elseif (argv$differential_abundance == "deseq2") {
#     cli_alert_info("Performing differential abundance analysis using DESeq2 package")
#     deseq_list <- DESeqWrapper(filtered_list, taxons, metadata)
#     cli_alert_success("DESeq2 successfully completed, output files written to {outdir}/deseq_results_(taxonlevel).csv")
# }
# else {
#     cli_alert_warning("Invalid method for differential abundance testing (supported methods: deseq2). Please check case.")
# }
# hm_list <- PlotDESeqData(deseq_list, taxons, metadata, pval = 0.1, lfc = 0.58)
# multi_page <- marrangeGrob(grobs = hm_list, ncol = 2, nrow = 2)
# ggsave(multi_page, filename = paste0(outdir, "/deseq_heatmap.pdf"), width = 25, height = 25)
# cli_alert_success("DESeq2 heatmaps saved to {outdir}/deseq_heatmap.pdf")


#### RELATIVE ABUNDANCE CHANGES BETWEEN CONDITIONS FOR TAXA OF INTEREST ####
cli_h1("Plotting differences in relative abundance values for taxa of interest")
cli_alert_info("Making relative abundance plots between conditions")
phylum_boxplots <- RelAbundanceBoxPlots(transformed_list, metadata, myColors, TaxaOI = "all", TaxLevel = "phylum")
cli_alert_success("Saved relative abundance boxplots in {outdir}/Relative_abundance_differences_phylum.pdf")
class_boxplots <- RelAbundanceBoxPlots(transformed_list, metadata, myColors, TaxaOI = "all", TaxLevel = "class")
cli_alert_success("Saved relative abundance boxplots in {outdir}/Relative_abundance_differences_class.pdf")
order_boxplots <- RelAbundanceBoxPlots(transformed_list, metadata, myColors, TaxaOI = "all", TaxLevel = "order")
cli_alert_success("Saved relative abundance boxplots in {outdir}/Relative_abundance_differences_order.pdf")
family_boxplots <- RelAbundanceBoxPlots(transformed_list, metadata, myColors, TaxaOI = "all", TaxLevel = "family")
cli_alert_success("Saved relative abundance boxplots in {outdir}/Relative_abundance_differences_family.pdf")
genus_boxplots <- RelAbundanceBoxPlots(transformed_list, metadata, myColors, TaxaOI = "all", TaxLevel = "genus")
cli_alert_success("Saved relative abundance boxplots in {outdir}/Relative_abundance_differences_genus.pdf")

cli_h1("Making alluvial plot for dominant phyla (mean relative abundance > 2%).")
alluvial_plot <- AlluvialPlot(transformed_list, metadata)
ggsave(alluvial_plot, filename = paste0(outdir, "Phyla_Alluvial_Plot.pdf"), height = 15, width = 25)
cli_alert_success("Saved alluvial plot to {outdir}/Phyla_Alluvial_Plot.pdf")

#### ANALYSIS ON THE DATA CONTAINING THE NEGATIVE CONTROL ####
##### Extract negative control (Cas soil) #####
# cli_alert_info("Extracting Cas soil (negative control) before subsetting data")
# neg_metadata <- subset(metadata, metadata$Abbrev == argv$organism | metadata$Compartment == "Cas")
# neg_freqtable <- dplyr::select(freqtable, 1:7, neg_metadata$Sample_Name)

# neg_TDL <- MakeTaxonDataList(neg_freqtable, taxons)
# neg_FL <- FilterAbs(neg_TDL, neg_metadata, taxons, Nreads = 2000)
# neg_TL <- Transform(neg_FL, taxons, Filter = "TotalAbundance")
# neg_Multi_Ordination <- OrdinateFunction(neg_TL, taxons, neg_metadata, Metric = "bray", subset_compartments = "FALSE")
# l_ply(taxons, .fun = function(x) {
#     ggsave(neg_Multi_Ordination[[x]], filename = paste0(outdir, "/", x, "_ordination_NEG.pdf"), width = 25, height = 25)
#     cli_alert_success("Ordination plot for {x} saved to {outdir}/{x}_ordination_NEG.pdf")
# })