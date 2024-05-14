source("~/Documents/projects/ASP/scripts/data_analysis_options.Rdata")
# count table
#freqtable <- read.table("~/Documents/bioinfo/decrypt/ASP/fungi/count_table.tsv", sep = "\t", header = TRUE, comment.char = "")
freqtable <- read.table("~/Documents/projects/ASP/fungi/count_table.tsv", sep = "\t", header = TRUE, comment.char = "")
colnames(freqtable) <- gsub("\\.", "-", colnames(freqtable))
freqtable$taxonomy <- gsub(".__", "", freqtable$taxonomy)
# metadata
columns <- c("Sample_Name", "Short_Sample", "Sample_ID", "Organism", "Genotype", "Compartment", "Timepoint", "Abbrev", "Summary")
metadata <- setNames(read.table("~/Documents/projects/ASP/fungi/sample-metadata-withEmpty.tsv", header = TRUE, sep = "\t", comment.char = ""), nm = columns)
if (!identical(sort(colnames(freqtable)[-1]), sort(metadata$Sample_Name))) {
    cli_alert_danger("The colnames of the count table (sample names) do not correspond to the Sample_Name column in metadata")
    stop()
}

freqtable <- freqtable %>%
    separate(col = taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
    replace(is.na(.), "unassigned") %>%
    mutate_all(list(~ str_remove_all(., "D_")))

metadata <- subset(metadata, metadata$Abbrev == "Hv")
metadata$Genotype[!metadata$Genotype == "BI1"] <- "WT"
freqtable <- freqtable %>%
    dplyr::select(1:7, metadata$Sample_Name)
outdir <- "~/Documents/projects/ASP/fungi/data_analysis/"

taxon_data_list <- MakeTaxonDataList(freqtable, taxons)
filtered_list <- FilterAbs(taxon_data_list, metadata, taxons, Nreads = 2000, Outdir = outdir)
transformed_list <- Transform(filtered_list, taxons, Filter = "TotalAbundance")
alpha_list <- GetAlphaList(metadata, filtered_list, taxons, index = "shannon")
alpha_plots <- AlphaPlotList(alpha_list, taxons, myColors, index = "shannon")

## Boxplots for relative abundance differences:
transformed <- transformed_list[["genus"]]
taxa_transf <- transformed %>%
    tibble::rownames_to_column("Sample_Name") %>%
    dplyr::left_join(metadata[c("Sample_Name", "Compartment", "Genotype", "Timepoint")], by = "Sample_Name") %>%
    tidyr::pivot_longer(-c("Sample_Name", "Compartment", "Genotype", "Timepoint"), names_to = "Taxon", values_to = "Abundance") %>%
    dplyr::mutate(Compartment = factor(Compartment, levels = c("Roots", "Rhizoplane", "Rhizosphere", "Soil")))
kw_test <- taxa_transf %>%
    group_by(Compartment, Timepoint, Taxon) %>%
    group_map(~ KW_wrapper(.y, .x$Abundance, .x$Genotype)) %>%
    bind_rows()
taxa_transf <- taxa_transf %>%
    dplyr::left_join(kw_test, by = c("Compartment", "Timepoint", "Taxon"))
# skipping the looping over compartments, and only considering roots
Comp <- "Roots"
# Comp <- "Rhizosphere"
data <- dplyr::filter(taxa_transf, Compartment == Comp)
tax_data <- data %>%
    dplyr::filter(p_value <= 0.05) %>%
    dplyr::filter(Taxon != "Fusarium")
fusarium_data <- data  %>% 
    dplyr::filter(p_value <= 0.05 & Taxon == "Fusarium")
#low_tax_data <- dplyr::filter(tax_data, Taxon %in% c("Ascobolus", "Cladorrhinum", "Iodophanus", "Slopeiomyces", "Trichocladium"))
low_tax_data <- dplyr::filter(tax_data, Taxon %in% c("F_Alteromonadaceae", "F_Paenibacillaceae", "F_Parachlamydiaceae", "F_Saccharimonadaceae")) %>% 
    dplyr::filter(!is.na(Genotype))
f_low_tax_data <- dplyr::filter(tax_data, Taxon %in% c("F_Mortierellaceae", "F_unidentified_3823")) %>% 
    dplyr::filter(!is.na(Genotype))
#mid_tax_data <- dplyr::filter(tax_data, Taxon %in% c("Bipolaris", "Dactylella", "Penicillium", "Periconia", "Schizothecium"))
mid_tax_data <- dplyr::filter(tax_data, Taxon %in% c("F_Fibrobacteraceae", "F_Cellvibrionaceae", "F_Chitinophagaceae", "F_Pseudomonadaceae", "F_Streptomycetaceae", "F_Pseudonocardiaceae")) %>% 
    dplyr::filter(!is.na(Genotype))
f_mid_tax_data <- dplyr::filter(tax_data, Taxon %in% c("F_Orbiliaceae", "F_unidentified_30", "F_Olpidiaceae")) %>% 
    dplyr::filter(!is.na(Genotype))
# unidentified_30 == Helotiales (order) and Helotiales_sp
taxa <- tax_data %>%
    dplyr::select(Taxon) %>%
    unique()
Colors = myColors

## would need to order them based on whether they are up or down. Maybe separate them, then patch
## add label for timepoint at which this happens

low_plot <- f_low_tax_data %>%
    ggplot(aes(y = factor(Genotype, levels = c("Hid4", "GP")), x = Abundance, fill = Genotype)) +
    geom_boxplot(outlier.shape = NA, width = .5) + # outlier.size = 3, outlier.colour = "black", ) ++
    #coord_cartesian(xlim = quantile(f_low_tax_data$Abundance, c(0, 0.99)))+
    geom_point(aes(shape = Genotype), alpha = .7, size = 3) +
    facet_grid(Taxon ~ . , scales = "free_x") +
    scale_fill_brewer(palette = "Pastel1") +
    #scale_color_manual(values = Colors) +
    theme_bw() +
    theme(
        strip.text = element_text(size = 8), legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)
    ) +
    #guides(x = guide_axis(angle = 45)) +
    labs(y = "", x = "") 

mid_plot <- f_mid_tax_data %>%
    ggplot(aes(y = factor(Genotype, levels = c("Hid4", "GP")), x = Abundance, fill = Genotype)) +
    geom_boxplot(outlier.shape = NA, width = .5) + # outlier.size = 3, outlier.colour = "black", ) ++
    #coord_cartesian(xlim = quantile(f_low_tax_data$Abundance, c(0, 0.99)))+
    geom_point(aes(shape = Genotype), alpha = .7, size = 3) +
    facet_grid(Taxon ~ . , scales = "free_x") +
    scale_fill_brewer(palette = "Pastel1") +
    theme_bw() +
    theme(
        strip.text = element_text(size = 8), legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)
    ) +
    #guides(x = guide_axis(angle = 45)) +
    labs(y = "", x = "") 

# fusarium_plot <- fusarium_data %>%
#     ggplot(aes(y = factor(Genotype, levels = c("BI1", "WT")), x = Abundance, color = Compartment)) +
#     geom_boxplot(outlier.shape = NA, width = .5) + # outlier.size = 3, outlier.colour = "black", ) ++
#     coord_cartesian(xlim = quantile(fusarium_data$Abundance, c(0, 0.99)))+
#     geom_point(aes(shape = Genotype), alpha = .7, size = 3) +
#     #facet_grid(Taxon ~ Compartment, scales = "free_x") +
#     scale_color_manual(values = Colors) +
#     theme_bw() +
#     theme(
#         strip.text.x = element_blank(), strip.text.y = element_text(size = 18), legend.title = element_text(size = 16), legend.text = element_text(size = 14),
#         axis.title = element_text(size = 20), axis.text = element_text(size = 16)
#     ) +
#     guides(x = guide_axis(angle = 45)) +
#     labs(y = "", x = "Relative abundance") 
#ggsave(relabundance, file = "~/Documents/projects/ASP/fungi/data_analysis/root_relative_abundance_significant_differences.pdf", height = 20, width = 15, limitsize = FALSE)

alpha_genus <- alpha_plots[["genus"]]

library(patchwork)
composite_plot <- pd_plot / Plot_c | (low_plot / mid_plot + plot_layout(guides = "collect", heights = c(4,6)))
composite_plot <- (low_plot / mid_plot / fusarium_plot + plot_layout(guides = 'collect', heights = c(4,5,1))) | alpha_genus 
composite_plot
ggsave(composite_plot, file = "~/Documents/projects/ASP/fungi/data_analysis/root_relative_abundance_significant_differences.pdf", height = 20, width = 15, limitsize = FALSE)

## maybe i can add one more plot below alpha_genus...




############# Phylogenetic Alpha-diversity #############
library(abdiv)
library(picante) # contains mpd, which I can't find in abdiv
## available alpha-div indices: abdiv::faith_pd(), abdiv::shannon(), abdiv::simpson(), picante::mpd()
## available beta-div indices : abdiv::bray_curtis(), abdiv::unifrac() (unweighted AND weighted), abdiv::phylosor() 
# phylosor is basically the same thing as unifrac, but there's a factor 2 difference
# BUT, for alpha for example, I don't need to compute alpha-phylo for each taxon level, because it uses the phylogenetic tree which is created based on the ASVs.
## phylo alpha
bact_tree <- tidytree::read.tree("~/Documents/projects/hid_microbiome/bacteria/3_analysis/3.3_taxonomy/tree_unrooted.nwk")
alpha_pd <- pd(t(counttable[,-c(1:7)]), bact_tree, include.root=F) %>% 
    rownames_to_column("Sample_Name")
alpha_pd <- left_join(alpha_pd, metadata, by = "Sample_Name")
## test normality 
shapiro.test(dplyr::filter(alpha_pd, Genotype == "Hid4")$PD)$p.value # 0.26 ok
shapiro.test(dplyr::filter(alpha_pd, Genotype == "GP")$PD)$p.value #0.93 ok
## anova and tukey
anova <- aov(data = alpha_pd, PD ~ Genotype)
tukey <- TukeyHSD(anova)
cld <- as.data.frame(multcompLetters4(anova, tukey)$Genotype$Letters) %>% 
    rownames_to_column("Genotype") %>% 
    setNames(nm = c("Genotype", "Group"))
## plot
pd_plot <- ggplot(alpha_pd, aes(x = Genotype, y = PD, fill = Genotype)) + 
    geom_boxplot(outlier.shape = 9, outlier.size = 3, outlier.colour = "black", width = .5) + 
    geom_point(aes(shape = Genotype), alpha = .7, size = 3) +
    theme_bw() +
    scale_fill_brewer(palette = "Pastel1") +
    labs(x = "Genotype", title = paste0("Faith's PD (alpha phylogenetic diversity).")) +
    theme(strip.text = element_text(size = 20), legend.title = element_text(size = 16), legend.text = element_text(size = 14),
    axis.title = element_text(size = 20), axis.text = element_text(size = 16), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    guides(x = guide_axis(angle = 45))+
    geom_text(data = cld, aes(x = Genotype, y = 55, label = Group))

all_alpha <- alpha_plots$order + pd_plot + plot_layout(guides = 'collect')


#### Modifications to 6_data_analysis.R for it to work with the new output format ####
## inside 3_analysis/3.3_taxonomy/ we have : asv2seq.tsv, count_table.tsv and taxonomy.tsv
## most important is count_table, but it counts for ASV1, 2, etc.
## taxonomy.tsv takes the ASV and matches it to the taxonomy
#### CHECK / MODIFY COUNT TABLE FORMAT ####
##### Load count table and change column names #####
cli_h1("Checking input format, and pre-processing data")
freqtable <- t(read.table("~/Documents/projects/hid_microbiome/bacteria/3_analysis/3.3_taxonomy/count_table.tsv", sep = "\t", header = TRUE, comment.char = ""))
colnames(freqtable) <- str_replace(colnames(freqtable), pattern = ".extended.*", "") ## remove extendedFrags
taxonomy <- read.table("~/Documents/projects/hid_microbiome/bacteria/3_analysis/3.3_taxonomy/taxonomy.tsv", sep = "\t", header = TRUE, comment.char = "") %>% 
    dplyr::select(-Seq)
# colnames(freqtable) <- gsub("\\.", "-", colnames(freqtable))
# freqtable$taxonomy <- gsub(".__", "", freqtable$taxonomy)

##### Load metadata table #####
## PUT COLUMNS IN OPTIONS FILE SO THE USER CAN CHANGE IT THERE
columns <- c("Sample_Name", "Short_Sample", "Sample_ID", "Genotype", "Compartment")
metadata <- setNames(read.table("~/Documents/projects/hid_microbiome/bacteria/sample-metadata.tsv", header = TRUE, sep = "\t", comment.char = "")[1:length(columns)], nm = columns)
##### Check if colnames of count table (samples) correspond exactly to Sample_Name column in metadata #####
# will need to change it so it doesn't use identical maybe, if bug in the future because some samples were removed from the data in dada2
if (length(max.col(sapply(sort(colnames(freqtable)), grepl, sort(metadata$Sample_Name)))) != length(colnames(freqtable))){
    cli_alert_danger("The colnames of the count table (sample names) do not correspond to the Sample_Name column in metadata")
    stop()
}

##### Separate taxonomy column in count table into subdivisions #####
counttable <- rownames_to_column(taxonomy, "ASV") %>% 
    left_join(rownames_to_column(as.data.frame(freqtable), "ASV"), by = "ASV") %>% 
    replace(is.na(.), "unassigned") %>% 
    column_to_rownames("ASV")
# freqtable <- as.data.frame(freqtable) %>%
#     rownames_to_column("ASV") %>% 
#     left_join(rownames_to_column(taxonomy, "ASV"), by = "ASV") %>% 
#     replace(is.na(.), "unassigned") 

##### Subset tables based on chosen organism #####
cli_alert_info("You've decided to focus your analysis on {argv$organism}, subsetting metadata and count table")
metadata <- subset(metadata, metadata$Abbrev == argv$organism)
metadata$Genotype[!metadata$Genotype == "BI1"] <- "WT"
counttable <- counttable %>%
    dplyr::select(1:7, metadata$Sample_Name)



####### Growth promotion tests GP and Hid4 #######
growth_prom <- read.csv2("~/Documents/projects/hid_microbiome/metadata/Root-Shoot-FW_AMP-GP-HID4_Mock-Sv-SvR11-R11.csv")
growth_prom$Weight_mg <- as.numeric(growth_prom$Weight_mg)
growth_prom$Length_mm <- as.numeric(growth_prom$Length_mm)

### FUNCTION ###
boxplot_func <- function(data, organ, measurement, ylab, ymax = 150){
    ## sub datasets
    apm <- dplyr::filter(data, Genotype == "APM" & Organ == organ)
    gp <- dplyr::filter(data, Genotype == "GP" & Organ == organ)
    hid <- dplyr::filter(data, Genotype == "Hid4" & Organ == organ)
    measure <- as.symbol(measurement)
    ## statistics
    apm_anova <- aov(data = apm, apm[[measurement]] ~ Condition)
    apm_tukey <- TukeyHSD(apm_anova)
    apm_cld <- as.data.frame(multcompLetters4(apm_anova, apm_tukey)$Condition$Letters) %>% 
        rownames_to_column("Condition") %>% 
        setNames(nm = c("Condition", "Group"))
    gp_anova <- aov(data = gp, gp[[measurement]] ~ Condition)
    gp_tukey <- TukeyHSD(gp_anova)
    gp_cld <- as.data.frame(multcompLetters4(gp_anova, gp_tukey)$Condition$Letters) %>% 
        rownames_to_column("Condition") %>% 
        setNames(nm = c("Condition", "Group"))
    hid_anova <- aov(data = hid, hid[[measurement]] ~ Condition)
    hid_tukey <- TukeyHSD(hid_anova)
    hid_cld <- as.data.frame(multcompLetters4(hid_anova, hid_tukey)$Condition$Letters) %>% 
        rownames_to_column("Condition") %>% 
        setNames(nm = c("Condition", "Group"))
    ## plots
    apm_plot <- ggplot(apm, aes(x = Condition, y = .data[[measurement]], fill = Condition)) + 
        geom_boxplot() + 
        #geom_jitter()+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_fill_brewer(palette = "Pastel1") +
        geom_text(data = apm_cld, aes(x = Condition, y = ymax, label = Group))+
        ggtitle("APM")+
        labs(y = ylab)
    gp_plot <- ggplot(gp, aes(x = Condition, y = .data[[measurement]], fill = Condition)) + 
        geom_boxplot() + 
        #geom_jitter()+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_fill_brewer(palette = "Pastel1") +
        geom_text(data = gp_cld, aes(x = Condition, y = ymax, label = Group))+
        ggtitle("GP")+
        labs(y = ylab)
    hid_plot <- ggplot(hid, aes(x = Condition, y = .data[[measurement]], fill = Condition)) + 
        geom_boxplot() + 
        #geom_jitter()+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_fill_brewer(palette = "Pastel1") +
        geom_text(data = hid_cld, aes(x = Condition, y = ymax, label = Group))+
        ggtitle("Hid4")+
        labs(y = ylab)
    ## full plot
    full_plot <- apm_plot + gp_plot + hid_plot + plot_layout(guides = 'collect')
    return(full_plot)
}

root_weight <- boxplot_func(data = growth_prom, organ = "Root", measurement = "Weight_mg", ylab = "Root Fresh Weight (mg)", ymax = 190)
root_length <- boxplot_func(data = growth_prom, organ = "Root", measurement = "Length_mm", ylab = "Root Length (cm)", ymax = 20)

shoot_weight <- boxplot_func(data = growth_prom, organ = "Shoot", measurement = "Weight_mg", ylab = "Shoot Fresh Weight (mg)", ymax = 300)
shoot_length <- boxplot_func(data = growth_prom, organ = "Shoot", measurement = "Length_mm", ylab = "Shoot Length (cm)", ymax = 20)

### beta-diversity ###



### problem with number of ASVs in tree and table
seqtab <- read.table("~/Documents/projects/hid_microbiome/fungi/3_analysis/3.2_trimming/dada2_flash2/seqtab.tsv")
seqtab_nochim <- read.table("~/Documents/projects/hid_microbiome/fungi/3_analysis/3.2_trimming/dada2_flash2/seqtab_nochim.tsv")
#load("/Users/sinaeda/lib/IDTAXA/SILVA_SSU_r138_2019.rdata")
load("/Users/sinaeda/lib/IDTAXA/UNITE_v2021_May2021.rdata")
length(colnames(seqtab))
length(colnames(seqtab_nochim))

dna <- DNAStringSet(getSequences(colnames(seqtab_nochim)))
ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=TRUE)
ids
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") 
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(colnames(seqtab_nochim))

taxonomy <- as.data.frame(taxid) %>% 
    dplyr::mutate(ASV = paste0("ASV", seq(nrow(taxid)))) %>%
    rownames_to_column("Seq") %>%
    column_to_rownames("ASV") %>%
    relocate("domain", "phylum", "class",  "order",  "family", "genus",  "species", "Seq")

asv2seq <- setNames(as.data.frame(cbind(rownames(taxonomy), taxonomy$Seq)), nm = c("ASV", "Seq"))
## change the column names of count_table
count_table <- seqtab_nochim %>%
  rename_with(~coalesce(asv2seq$ASV[match(., asv2seq$Seq)], .))

writeFasta(setNames(as_tibble(asv2seq), nm = c("name", "seq")), "~/Documents/projects/hid_microbiome/bacteria/3_analysis/3.3_taxonomy/ASV.fasta")



### significance of vairance explained by Genotype ###
f_w_unifrac <- usedist::dist_make(t(counttable[,-c(1:7)]), weighted_unifrac, r_fung_tree)
wuf_dbrda <- capscale(f_w_unifrac ~ Genotype, data = Metadata, comm = t(counttable[,-c(1:7)]), sqrt.dist = TRUE)
Sites_c <- data.frame(scores(wuf_dbrda)$sites)
Expl_c <- data.frame(summary(eigenvals(wuf_dbrda)))[2, c(1, 2)]

## Perform betadisper tests for Genotype, Timepoint, and GenoTime
bd <- betadisper(f_w_unifrac, group = Metadata$Genotype, type = "centroid")
## Perform significance test on the dispersion analysis using permutest
# betadisper+permutest tests homogeneity of dispersion among groups (Genotype)
permu_scheme <- how(
    within = Within(type = "free"), nperm = 9999 #, blocks = Metadata$Genotype # , blocks = data$Compartment
)
disp <- permutest(bd, permutations = permu_scheme, pairwise = TRUE)
## Permanova on bc
ado <- adonis2(f_w_unifrac ~ Genotype, data = Metadata, permutations = permu_scheme, by = NULL)

core_microbiota
asv_counts <- dplyr::select(counttable[,-c(1:7)], -c("BactV5V7-BC52-Cas_soil_S151_L001"))
asv_presence <- asv_counts %>% dplyr::mutate(across(c(1:ncol(asv_counts)), ~ ifelse(. > 0, 1, 0)))

occurence <- asv_presence %>% rownames_to_column("ASV") %>% 
    pivot_longer(-ASV) %>% 
    dplyr::group_by(ASV) %>% 
    dplyr::summarise(mean(value)) %>% 
    setNames(nm = c("ASV", "Occurence"))
taxonomy <- ReplaceUnassigned(taxonomy)
core90 <- dplyr::filter(occurence, Occurence > 0.9) %>% 
    left_join(rownames_to_column(taxonomy, "ASV"))
core95 <- dplyr::filter(occurence, Occurence > 0.95)%>% 
    left_join(rownames_to_column(taxonomy, "ASV"))
core100 <- dplyr::filter(occurence, Occurence == 1) %>% 
    left_join(rownames_to_column(taxonomy, "ASV")) 

transf_counts <- decostand(asv_counts, method = "total", zap = TRUE)
test <- core100 %>% dplyr::left_join(rownames_to_column(transf_counts, "ASV")) %>% 
    dplyr::select(-c("ASV", "Occurence", "domain", "phylum", "class", "order", "genus", "species")) %>% 
    dplyr::mutate(family = ifelse(is.na(family), "unassigned", family)) %>%
    pivot_longer(cols = -c("family"), names_to = "Sample_Name", values_to = "relative") %>% 
    left_join(metadata, by = "Sample_Name") %>% 
    group_by(family, Genotype) %>%
    dplyr::mutate(fam_geno = paste0(Genotype, "_", family)) 

ggplot(test, aes(x = family, y = relative, shape = Genotype, fill = Genotype)) + 
    geom_boxplot()+
    guides(x = guide_axis(angle = 45))+
    theme_bw()

asv_pres_hid <- dplyr::select(asv_presence, contains("Hid"))
occ_hid <- asv_pres_hid %>% rownames_to_column("ASV") %>% 
    pivot_longer(-ASV) %>% 
    dplyr::group_by(ASV) %>% 
    dplyr::summarise(mean(value)) %>% 
    setNames(nm = c("ASV", "Occurence"))
hid_core90 <- dplyr::filter(occ_hid, Occurence > 0.9) %>% 
    left_join(rownames_to_column(taxonomy, "ASV"))
hid_core95 <- dplyr::filter(occ_hid, Occurence > 0.95)%>% 
    left_join(rownames_to_column(taxonomy, "ASV"))
hid_core100 <- dplyr::filter(occ_hid, Occurence == 1) %>% 
    left_join(rownames_to_column(taxonomy, "ASV"))
## core hid4 always 57 ASVs
hid_test <- hid_core100 %>% dplyr::left_join(rownames_to_column(transf_counts, "ASV")) %>% 
    dplyr::select(-c("ASV", "Occurence", "domain", "phylum", "class", "order", "genus", "species")) %>% 
    dplyr::mutate(family = ifelse(is.na(family), "unassigned", family)) %>%
    dplyr::select(c(family, contains("Hid"))) %>% 
    pivot_longer(cols = -c("family"), names_to = "Sample_Name", values_to = "relative") %>% 
    left_join(metadata, by = "Sample_Name") %>% 
    group_by(family, Genotype) %>%
    dplyr::mutate(fam_geno = paste0(Genotype, "_", family)) 
ggplot(hid_test, aes(x = family, y = relative, shape = Genotype, fill = family)) + 
    geom_boxplot()+
    guides(x = guide_axis(angle = 45))+
    theme_bw()

asv_pres_gp <- dplyr::select(asv_presence, contains("GP"))
occ_gp <- asv_pres_gp %>% rownames_to_column("ASV") %>% 
    pivot_longer(-ASV) %>% 
    dplyr::group_by(ASV) %>% 
    dplyr::summarise(mean(value)) %>% 
    setNames(nm = c("ASV", "Occurence"))
gp_core90 <- dplyr::filter(occ_gp, Occurence > 0.9) %>% 
    left_join(rownames_to_column(taxonomy, "ASV"))
gp_core95 <- dplyr::filter(occ_gp, Occurence > 0.95)%>% 
    left_join(rownames_to_column(taxonomy, "ASV"))
gp_core100 <- dplyr::filter(occ_gp, Occurence == 1) %>% 
    left_join(rownames_to_column(taxonomy, "ASV"))
## core GP always 60 strains
gp_test <- gp_core100 %>% dplyr::left_join(rownames_to_column(transf_counts, "ASV")) %>% 
    dplyr::select(-c("ASV", "Occurence", "domain", "phylum", "class", "order", "genus", "species")) %>% 
    dplyr::mutate(family = ifelse(is.na(family), "unassigned", family)) %>%
    dplyr::select(c(family, contains("GP"))) %>% 
    pivot_longer(cols = -c("family"), names_to = "Sample_Name", values_to = "relative") %>% 
    left_join(metadata, by = "Sample_Name") %>% 
    group_by(family, Genotype) %>%
    dplyr::mutate(fam_geno = paste0(Genotype, "_", family)) 
ggplot(gp_test, aes(x = family, y = relative, shape = Genotype, fill = family)) + 
    geom_boxplot()+
    guides(x = guide_axis(angle = 45))+
    theme_bw()


gp_families <- unique(gp_test$family)
hid_families <- unique(hid_test$family)
both_families <- unique(test$family)

intersect(gp_families, both_families)
intersect(hid_families, both_families)
# test %>% 
#     ungroup() %>% 
#     group_by(Sample_Name) %>% 
#     dplyr::filter(Sample_Name == "BactV5V7-BC44-Hid4_12_S143_L001") %>% 
#     mutate(TotalRelative = sum(relative))   %>%  print(n = 36)

ReplaceUnassigned <- function (FreqTable) 
{
    FreqTable <- FreqTable %>% dplyr::mutate(phylum = ifelse(phylum == 
        "unassigned", paste0("D_", domain), paste0("P_", phylum))) %>% 
        dplyr::mutate(class = ifelse(class == "unassigned", phylum, 
            paste0("C_", class))) %>% dplyr::mutate(order = ifelse(order == 
        "unassigned", class, paste0("O_", order))) %>% dplyr::mutate(family = ifelse(family == 
        "unassigned", order, paste0("F_", family))) %>% dplyr::mutate(genus = ifelse(genus == 
        "unassigned", family, paste0("G_", genus))) #%>% 
        #dplyr::mutate(species = ifelse(species == "unassigned", genus, paste0("S_", species)))
    return(FreqTable)
}
