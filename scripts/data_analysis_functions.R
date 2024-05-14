#### Set theme for plots ####
# This theme extends the 'theme_light' that comes with ggplot2.
# The "Lato" font is used as the base font. This is similar
# to the original font in Cedric's work, Avenir Next Condensed.
# theme_set(theme_light(base_family = "Lato"))

# theme_update(
#     # Remove title for both x and y axes
#     # axis.title = element_blank(),
#     # Axes labels are grey
#     axis.text = element_text(color = "grey40"),
#     # The size of the axes labels are different for x and y.
#     axis.text.x = element_text(size = 10, margin = margin(t = 5)),
#     axis.text.y = element_text(size = 10, margin = margin(r = 5)),
#     # Also, the ticks have a very light grey color
#     axis.ticks = element_line(color = "grey91", size = .5),
#     # The length of the axis ticks is increased.
#     axis.ticks.length.x = unit(.5, "lines"),
#     axis.ticks.length.y = unit(.5, "lines"),
#     # Remove the grid lines that come with ggplot2 plots by default
#     # panel.grid = element_blank(),
#     # Customize margin values (top, right, bottom, left)
#     plot.margin = margin(10, 10, 10, 10),
#     # Use a light grey color for the background of both the plot and the panel
#     plot.background = element_rect(fill = "grey98", color = "grey98"),
#     panel.background = element_rect(fill = "grey98", color = "grey98"),
#     # Customize title appearence
#     plot.title = element_text(
#         color = "grey10",
#         size = 10,
#         face = "bold",
#         margin = margin(t = 10)
#     ),
#     # Customize subtitle appearence
#     plot.subtitle = element_markdown(
#         color = "grey30",
#         size = 8,
#         lineheight = 1.35,
#         margin = margin(t = 10, b = 20)
#     ),
#     # Title and caption are going to be aligned
#     plot.title.position = "plot",
#     plot.caption.position = "plot",
#     plot.caption = element_text(
#         color = "grey30",
#         size = 8,
#         lineheight = 1.2,
#         hjust = 0,
#         margin = margin(t = 20) # Large margin on the top of the caption.
#     ),
#     # Remove legend
#     # legend.position = "none"
#     # Change the background of the legend
#     legend.background = element_rect(fill = "grey98", color = "grey98"),
#     # Add a square around the legend
#     # legend.box.background = element_rect(colour = "black")
# )

ReplaceUnassigned <- function(FreqTable){
    FreqTable <- FreqTable %>% 
        dplyr::mutate(phylum = ifelse(phylum == "unassigned", paste0("D_", domain), paste0("P_", phylum))) %>% 
        dplyr::mutate(class = ifelse(class == "unassigned", phylum, paste0("C_", class))) %>% 
        dplyr::mutate(order = ifelse(order == "unassigned", class, paste0("O_", order))) %>% 
        dplyr::mutate(family = ifelse(family == "unassigned", order, paste0("F_", family))) %>% 
        dplyr::mutate(genus = ifelse(genus == "unassigned", family, paste0("G_", genus))) %>% 
        dplyr::mutate(species = ifelse(species == "unassigned", genus, paste0("S_", species)))

    return(FreqTable)
}

#### Functions used in diversity analysis for Amplicon Sequencing data ####
MakeTaxonDataList <- function(FreqTable, TaxonLevels) {
    TaxonDataList <- list()
    for (i in 1:length(TaxonLevels)){
        TaxonLvl <- TaxonLevels[[i]]
        TaxonData <- FreqTable %>%
            select(c(!!sym(TaxonLvl), as.integer(8:ncol(FreqTable)))) %>%
            group_by(!!sym(TaxonLvl)) %>%
            mutate_at(vars(-group_cols()), ~ as.integer(.)) %>%
            summarise_all(sum) %>%
            column_to_rownames(TaxonLvl) %>%
            t()
        saveRDS(TaxonData, file="temp.RData")
        TaxonDataList[[TaxonLvl]] <- TaxonData
    }
    return(TaxonDataList)
}
FilterAbs <- function(TaxonDataList, Metadata, TaxonLevels, Nreads = 2000, Outdir) {
    # Remove unidentified and unassigned reads BEFORE filtering the rest
    # Filters out taxons where absolute abundance == 0
    # Filters out samples with less than Nreads reads
    FilteredList <- list()
    LostList <- list()
    for (i in 1:length(TaxonLevels)) {
        level <- TaxonLevels[i]
        AbsAbundance <- TaxonDataList[[level]]
        Filtered <- as_tibble(rownames_to_column(data.frame(AbsAbundance), var = "Sample_Name")) %>%
            dplyr::select_if(!names(.) %in% c("unassigned", "unidentified")) %>%
            pivot_longer(-Sample_Name) %>%
            group_by(Sample_Name) %>%
            dplyr::mutate(total = sum(value)) %>%
            filter(total > Nreads) %>%
            arrange(total) %>%
            group_by(name) %>%
            dplyr::mutate(total = sum(value)) %>%
            filter(total != 0) %>%
            ungroup() %>%
            dplyr::select(-total) %>%
            pivot_wider(id_cols = Sample_Name) %>%
            as.data.frame() %>%
            column_to_rownames("Sample_Name")
        FilteredList[[level]] <- Filtered
        lost <- Metadata %>%
            dplyr::filter(!Sample_Name %in% rownames(Filtered))
        write.table(lost,
            file = paste0(Outdir, "/lost_during_filter_", level, ".tsv"),
            append = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE
        )
    }
    return(FilteredList)
}
Transform <- function(FilteredList, TaxonLevels, Filter = "TotalAbundance") {
    TransformedList <- list()
    for (i in 1:length(TaxonLevels)) {
        level <- TaxonLevels[i]
        Filtered <- FilteredList[[level]]
        ## 1. Transform absolute abundance into relative abundance
        RelAbundance <- decostand(Filtered, method = "total", zap = TRUE)
        ## 2. Filter out taxa which represent < 0.01% of frequencies
        if (Filter == "TotalAbundance") {
            TotalAbundance <- sapply(c(1:ncol(RelAbundance)), FUN = function(x) {
                sum(RelAbundance[, x])
            }) %>% setNames(nm = colnames(RelAbundance[, 1:ncol(RelAbundance)]))
            ToKeep <- names(TotalAbundance[TotalAbundance > 0.01])
            RelAbundance <- data.frame(RelAbundance[, ToKeep])
        }
        # Duran 2021: ASVs w/ > 0.1% of relAbundance in at least 1 sample (keeps very very few orders, wouldn't do this)
        # mainly because some are very much over-represented, like Pezizales
        else if (Filter == "Duran2021") {
            ToKeep <- rownames_to_column(RelAbundance) %>%
                pivot_longer(-rowname) %>%
                group_by(name) %>%
                dplyr::filter(value >= 0.1) %>%
                dplyr::select(name) %>%
                unique()
            RelAbundance <- data.frame(RelAbundance[, ToKeep])
        }
        TransformedList[[level]] <- RelAbundance
    }
    return(TransformedList)
}
KW_wrapper <- function(keys, alpha, comparison) {
    if (length(unique(comparison)) > 1) {
        result <- kruskal.test(alpha, comparison)
        tibble(
            keys,
            p_value = result$p.value
        )
    }
}

GetAlphaList <- function(Metadata, FilteredList, TaxonLevels, index) {
    AlphaList <- list()
    for (i in 1:length(TaxonLevels)) {
        level <- TaxonLevels[i]
        data <- rownames_to_column(data.frame(FilteredList[[level]]), var = "Sample_Name")
        index_data <- setNames(data.frame(cbind(data[, 1], diversity(data[, -1], index))), nm = c("Sample_Name", index))
        index_data <- merge(Metadata, index_data, by = "Sample_Name")
        index_data[, ncol(index_data)] <- as.numeric(index_data[, ncol(index_data)])

        kw_test <- index_data %>%
            group_by(Compartment, Timepoint) %>%
            group_map(~ KW_wrapper(.y, .x$shannon, .x$Genotype)) %>%
            bind_rows()

        index_data <- index_data %>%
            dplyr::left_join(kw_test, by = c("Compartment", "Timepoint"))
        AlphaList[[level]] <- index_data
    }
    return(AlphaList)
}
AlphaPlotList <- function(AlphaList, TaxonLevels, Colors, index) {
    PlotList <- list()
    for (i in 1:length(TaxonLevels)){
        level <- TaxonLevels[i]
        data <- AlphaList[[level]] %>%
            dplyr::mutate(GenoTime = paste0(AlphaList[[level]]$Genotype, "_", gsub(" ", "", AlphaList[[level]]$Timepoint))) %>%
            dplyr::mutate(Compartment = factor(Compartment, levels = c("Roots", "Rhizoplane", "Rhizosphere", "Soil")))
        p <- data %>%
            ggplot(aes(x = factor(Genotype, levels = c("WT", "BI1")), y = !!sym(index), col = Compartment)) +
            geom_boxplot(outlier.shape = 9, outlier.size = 3, outlier.colour = "black", width = .5) + # , position = position_dodge(width = .15)) +
            geom_point(aes(shape = Genotype), alpha = .7, size = 3) + # , position = position_dodge(width = .15)) +
            scale_color_manual(values = Colors) +
            facet_grid(Timepoint ~ Compartment) +
            theme_bw() +
            theme(strip.text = element_text(size = 20), legend.title = element_text(size = 16), legend.text = element_text(size = 14),
            axis.title = element_text(size = 20), axis.text = element_text(size = 16))+
            guides(x = guide_axis(angle = 45)) +
            labs(x = "Compartment", title = paste0("Alpha diversity plot (", index, " index) at the ", level, " level.")) +
            ylim(0, 5)

        text_df <- data %>%
            dplyr::select("Compartment", "Timepoint", "p_value") %>%
            unique()
        p <- p + annotate("segment", x = .7, xend = 2.3, y = 4.5, yend = 4.5) +
            geom_text(data = text_df, mapping = aes(x = 1.5, y = 4.6, label = round(p_value, 3)), color = "black")

        PlotList[[level]] <- p
    }    
    return(PlotList)
}
getPWDistList <- function(TransformedList, TaxonLevels, Metric) {
    DistList <- list()
    for (i in 1:length(TaxonLevels)) {
        level <- TaxonLevels[i]
        RelAbundance <- TransformedList[[level]]
        Dist <- as.matrix(vegdist(RelAbundance, method = Metric))
        DistList[[level]] <- Dist
    }
    return(DistList)
}
OrdinateFunction <- function(TransformedList, TaxonLevels, Metadata, Metric, subset_compartments = "TRUE") {
    Compartments <- unique(Metadata$Compartment)
    PlotList <- list()
    for (i in 1:length(TaxonLevels)){
        level <- TaxonLevels[i]
        RelAbundance <- TransformedList[[level]]
        ## 1. Get specified distance matrix
        Dist <- as.matrix(vegdist(RelAbundance, method = Metric))
        ## 2. Select metadata corresponding to current analysis
        Metadata <- Metadata %>% dplyr::filter(Sample_Name %in% rownames(RelAbundance))
        rownames(Metadata) <- Metadata$Sample_Name
        Metadata <- Metadata[rownames(RelAbundance), ]
        if (!identical(rownames(RelAbundance), rownames(Metadata))) {
            stop("Sample Names are not in the same order in both dataframes")
        }
        Organism <- unique(Metadata$Abbrev)
        ## 3. Define PlotList and add plots
        subPlotList <- list()
        ## 4. Ordinate (PCoA and dbRDA for all compartments)
        PCoA <- capscale(Dist ~ 1, data = Metadata, comm = RelAbundance, sqrt.dist = TRUE) # , add = "lingoes")
        dbRDA <- capscale(Dist ~ Compartment + Genotype, data = Metadata, comm = RelAbundance, sqrt.dist = TRUE) # , add = "lingoes")
        ## 5. Extract data from ordinations (sites, ...) and plot
        Sites_u <- data.frame(scores(PCoA)$sites)
        Expl_u <- data.frame(summary(eigenvals(PCoA)))[2, c("MDS1", "MDS2")]
        Sites_c <- data.frame(scores(dbRDA)$sites)
        Expl_c <- data.frame(summary(eigenvals(dbRDA)))[2, c("CAP1", "CAP2")]
        Plot_u <- ggplot(Sites_u, aes(MDS1, MDS2, col = Metadata$Compartment, shape = Metadata$Genotype, label = Metadata$Sample_ID)) +
            geom_point(size = 1.5) +
            geom_hline(yintercept = 0, linetype = 3) +
            geom_vline(xintercept = 0, linetype = 3) +
            labs(x = paste0("MDS1 (", round(Expl_u$MDS1 * 100, 1), "%)"), y = paste0("MDS2 (", round(Expl_u$MDS2 * 100, 1), "%)")) +
            theme_classic() +
            geom_text(hjust = 1.5, vjust = 0) +
            ggtitle("Unconstrained PCoA (all compartments).") +
            theme(strip.text = element_text(size = 30), legend.title = element_text(size = 26), legend.text = element_text(size = 24),
            axis.title = element_text(size = 30), axis.text = element_text(size = 26))+
            stat_ellipse(geom = "polygon", alpha = .2, aes(MDS1, MDS2, col = Metadata$Compartment, fill = Metadata$Compartment), inherit.aes = FALSE)
        Plot_c <- ggplot(Sites_c, aes(CAP1, CAP2, col = Metadata$Compartment, shape = Metadata$Genotype, label = Metadata$Sample_ID)) +
            geom_point(size = 1.5) +
            geom_hline(yintercept = 0, linetype = 3) +
            geom_vline(xintercept = 0, linetype = 3) +
            theme_classic() +
            labs(x = paste0("CAP1 (", round(Expl_c$CAP1 * 100, 1), "%)"), y = paste0("CAP2 (", round(Expl_c$CAP2 * 100, 1), "%)")) +
            geom_text(hjust = 1.5, vjust = 0) +
            ggtitle("Constrained dbRDA (all compartments).") +
            theme(strip.text = element_text(size = 30), legend.title = element_text(size = 26), legend.text = element_text(size = 24),
            axis.title = element_text(size = 30), axis.text = element_text(size = 26))+
            stat_ellipse(geom = "polygon", alpha = .2, aes(CAP1, CAP2, col = Metadata$Compartment, fill = Metadata$Compartment), inherit.aes = FALSE) +
            geom_segment(aes(x = 0.0, y = 0.0, xend = scores(dbRDA)$biplot[1, 1], yend = scores(dbRDA)$biplot[1, 2]), arrow = arrow(), col = "black") +
            geom_label(data = data.frame(scores(dbRDA)$biplot), aes(CAP1, CAP2, label = rownames(scores(dbRDA)$biplot)), inherit.aes = FALSE, nudge_x = .1, nudge_y = .1) +
            geom_segment(aes(x = 0.0, y = 0.0, xend = scores(dbRDA)$biplot[2, 1], yend = scores(dbRDA)$biplot[2, 2]), arrow = arrow(), col = "black") +
            geom_segment(aes(x = 0.0, y = 0.0, xend = scores(dbRDA)$biplot[3, 1], yend = scores(dbRDA)$biplot[3, 2]), arrow = arrow(), col = "black")
        subPlotList[["PCoA"]] <- Plot_u
        subPlotList[["dbRDA"]] <- Plot_c
        ## 5. Define subsets for each compartment
        if (subset_compartments == "TRUE") {
            for (i in 1:length(Compartments)) {
                meta <- subset(Metadata, Metadata$Compartment == Compartments[i])
                CompData <- RelAbundance %>%
                    rownames_to_column() %>%
                    dplyr::filter(rowname %in% meta$Sample_Name) %>%
                    column_to_rownames("rowname")
                CompDist <- vegdist(CompData, method = Metric)
                CompPCoA <- capscale(CompDist ~ 1, data = meta, comm = CompData, sqrt.dist = TRUE, add = "lingoes")
                CompdbRDA <- capscale(CompDist ~ Genotype + Timepoint, data = meta, comm = CompData, sqrt.dist = TRUE, add = "lingoes")
                Sites_u <- data.frame(scores(CompPCoA)$sites)
                Expl_u <- data.frame(summary(eigenvals(CompPCoA)))[2, c("MDS1", "MDS2")]
                Sites_c <- data.frame(scores(CompdbRDA)$sites)
                Expl_c <- data.frame(summary(eigenvals(CompdbRDA)))[2, c("CAP1", "CAP2")]
                PlotComp_u <- ggplot(Sites_u, aes(MDS1, MDS2, col = meta$Genotype, shape = meta$Timepoint, label = meta$Sample_ID)) +
                    geom_point(size = 1.5) +
                    geom_hline(yintercept = 0, linetype = 3) +
                    geom_vline(xintercept = 0, linetype = 3) +
                    theme_classic() +
                    labs(x = paste0("MDS1 (", round(Expl_u$MDS1 * 100, 1), "%)"), y = paste0("MDS2 (", round(Expl_u$MDS2 * 100, 1), "%)")) +
                    geom_text(hjust = 1.5, vjust = 0) +
                    theme(strip.text = element_text(size = 30), legend.title = element_text(size = 26), legend.text = element_text(size = 24),
                    axis.title = element_text(size = 30), axis.text = element_text(size = 26))+
                    ggtitle(paste0("Unconstrained PCoA (", Compartments[i], ").")) +
                    stat_ellipse(geom = "polygon", alpha = .2, aes(MDS1, MDS2, col = meta$Genotype, fill = meta$Genotype), inherit.aes = FALSE)
                PlotComp_c <- ggplot(Sites_c, aes(CAP1, CAP2, col = meta$Genotype, shape = meta$Timepoint, label = meta$Sample_ID)) +
                    geom_point(size = 1.5) +
                    geom_hline(yintercept = 0, linetype = 3) +
                    geom_vline(xintercept = 0, linetype = 3) +
                    theme_classic() +
                    labs(x = paste0("CAP1 (", round(Expl_c$CAP1 * 100, 1), "%)"), y = paste0("CAP2 (", round(Expl_c$CAP2 * 100, 1), "%)")) +
                    geom_text(hjust = 1.5, vjust = 0) +
                    theme(strip.text = element_text(size = 30), legend.title = element_text(size = 26), legend.text = element_text(size = 24),
                    axis.title = element_text(size = 30), axis.text = element_text(size = 26))+
                    ggtitle(paste0("Constrained dbRDA (", Compartments[i], ").")) +
                    stat_ellipse(geom = "polygon", alpha = .2, aes(CAP1, CAP2, col = meta$Genotype, fill = meta$Genotype), inherit.aes = FALSE) +
                    geom_segment(aes(x = 0.0, y = 0.0, xend = scores(CompdbRDA)$biplot[1, 1], yend = scores(CompdbRDA)$biplot[1, 2]), arrow = arrow(), col = "black") +
                    geom_label(data = data.frame(scores(CompdbRDA)$biplot), aes(CAP1, CAP2, label = rownames(scores(CompdbRDA)$biplot)), inherit.aes = FALSE, nudge_x = .1, nudge_y = .1) +
                    geom_segment(aes(x = 0.0, y = 0.0, xend = scores(CompdbRDA)$biplot[2, 1], yend = scores(CompdbRDA)$biplot[2, 2]), arrow = arrow(), col = "black")
                PlotComp <- grid.arrange(PlotComp_u, PlotComp_c, ncol = 2, nrow = 1)
                subPlotList[[Compartments[i]]] <- PlotComp
            }
        }
        multi_page <- marrangeGrob(grobs = subPlotList, nrow = 2, ncol = 1)
        PlotList[[level]] <- multi_page
        ## 6. Print all these plots to pdf with ggarrange
        #ggsave(multi_page, filename = paste0(outdir, Organism, "_", suffix, "_ordination.pdf"), width = 25, height = 25)
    }
    ## 6. Print all these plots to pdf with ggarrange
    #ggsave(multi_page, filename = paste0("~/Documents/bioinfo/R/decrypt/diversity/", Organism, "_", suffix, "_ordination.pdf"), width = 25, height = 25)
    return(PlotList)
}
AdonisWrapper <- function(TransformedList, TaxonLevels, Metadata, permu_scheme) {
    compartments <- unique(Metadata$Compartment)
    CompList <- list()
    for (i in 1:length(TaxonLevels)) {
        level <- TaxonLevels[i]
        transformed <- TransformedList[[level]]
        ## First perform betadisper and permutest to check for dispersion of datapoints between different groups
        for (j in 1:length(compartments)) {
            comp <- compartments[j]
            print(paste0(level, "<->", comp))
            subMeta <- Metadata[Metadata$Compartment == comp, ]
            subCount <- transformed[rownames(transformed) %in% subMeta$Sample_Name, ]
            subMeta <- subMeta[subMeta$Sample_Name %in% rownames(subCount), ]
            # print(rownames(subCount))
            # print(subMeta$Sample_Name)
            subBC <- vegdist(subCount)
            ## Perform betadisper tests for Genotype, Timepoint, and GenoTime
            bd_G <- betadisper(subBC, group = subMeta$Genotype, type = "centroid")
            bd_T <- betadisper(subBC, group = subMeta$Timepoint, type = "centroid")
            #bd_GT <- betadisper(subBC, group = subMeta$GenoTime, type = "centroid")
            ## Perform significance test on the dispersion analysis using permutest
            disp_G <- permutest(bd_G, permutations = permu_scheme, pairwise = TRUE)
            disp_T <- permutest(bd_T, permutations = permu_scheme, pairwise = TRUE)
            #disp_GT <- permutest(bd_GT, permutations = permu_scheme, pairwise = TRUE)
            ## Permanova on bc
            ado_G <- adonis2(subBC ~ Genotype, data = subMeta, permutations = permu_scheme, by = NULL)
            ado_T <- adonis2(subBC ~ Timepoint, data = subMeta, permutations = permu_scheme, by = NULL)
            #ado_GT <- adonis2(subBC ~ GenoTime, data = subMeta, permutations = permu_scheme, by = NULL)
            CompList[[level]][[comp]][["Genotype"]] <- list(bd_G, disp_G, ado_G)
            CompList[[level]][[comp]][["Timepoint"]] <- list(bd_T, disp_T, ado_T)
            #CompList[[level]][[comp]][["GenoTime"]] <- list(bd_GT, disp_GT, ado_GT)
        }
    }
    return(CompList)
}
DESeqWrapper <- function(FilteredList, TaxonLevels, Metadata){
    Compartments <- unique(Metadata$Compartment)
    Metadata <- Metadata %>%
        mutate(GenoTime = paste0(Metadata$Genotype, "_", gsub(" ", "", Metadata$Timepoint)))
    DS_list <- list()
    #VDS_list <- list()
    for (i in 1:length(TaxonLevels)){
        level <- TaxonLevels[i]
        filtered <- FilteredList[[level]]
        DS_results <- llply(Compartments, .fun = function(x){
            coldata <- Metadata[Metadata$Compartment == x, ] %>%
                remove_rownames() %>%
                column_to_rownames("Sample_Name") %>% 
                mutate(GenoTime = factor(GenoTime, levels = c("WT_14dpi", "BI1_14dpi", "WT_28dpi", "BI1_28dpi")))
            countdata <- t(filtered)
            countdata <- countdata[, colnames(countdata) %in% rownames(coldata)]
            coldata <- coldata[rownames(coldata) %in% colnames(countdata), ]
            countdata <- countdata[, rownames(coldata)]
            ds2 <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, tidy = FALSE, design = ~ GenoTime)# + Timepoint + Genotype:Timepoint, tidy = FALSE)
            dds <- DESeq(ds2, fitType = "local", sfType = "poscounts", quiet = TRUE)
            ## add shrinkage step to avoid misinterpreting very low count taxa fold changes
            # vds <- varianceStabilizingTransformation(dds)
            res <- lfcShrink(dds, coef = "GenoTime_BI1_14dpi_vs_WT_14dpi") # , type = "apeglm", quiet = TRUE)
            deseq_results <- as.data.frame(res) %>%
                rownames_to_column("taxon") %>%
                arrange(pvalue, log2FoldChange)
            return(deseq_results)
        })
        names(DS_results) <- Compartments
        DS_results <- bind_rows(DS_results, .id = "Compartment")
        write.table(DS_results, file = paste0(outdir, "/deseq_results_", level, ".csv"), sep = ";")
        DS_list[[level]] <- DS_results
    }
    return(DS_list)
}
PlotDESeqData <- function(DESeqList, TaxonLevels, Metadata, pval, lfc){
    Compartments <- setNames(c(rep(NA, length(unique(Metadata$Compartment)))), nm = unique(Metadata$Compartment))
    hm_list <- list()
    for (i in 1:length(TaxonLevels)) {
        level <- TaxonLevels[i]
        deseq <- DESeqList[[level]]
        hm_input <- deseq %>%
            dplyr::filter(abs(log2FoldChange) >= lfc & padj <= pval) %>% 
            dplyr::select(taxon, Compartment, log2FoldChange) %>%
            pivot_wider(names_from = Compartment, values_from = log2FoldChange) %>%
            add_column(!!!Compartments[!names(Compartments) %in% names(.)]) %>% 
            relocate(c("soil", "Rhizosphere", "roots")) %>%
            column_to_rownames("taxon") %>%
            replace(is.na(.), 0) %>%
            setNames(nm = c("Soil", "Rhizosphere", "Roots"))
        hm <- pheatmap(hm_input,
            border_color = "white", cluster_rows = TRUE, cluster_cols = FALSE,
            display_numbers = TRUE, number_color = "black", #fontsize_number = 10,
            show_colnames = TRUE, labels_col = colnames(hm_input), angle_col = 0,
            #fontsize = 12, #cellwidth = 50, cellheight = 50,
            main = paste0("Differential abundance at the ", level, " level in BI1 mutant compared to WT.") # ,
            # filename = paste0(outdir, "/deseq_heatmap.pdf")
        )
        hm_list[[level]] <- as.grob(hm)
    }
    return(hm_list)
}
# RelAbundanceBoxPlots <- function(TransformedList, Metadata, TaxaOI = "all", Colors) {
#     transformed <- TransformedList$order
#     if (TaxaOI != "all") {
#         transformed <- transformed %>%
#             dplyr::select(any_of(TaxaOI))
#     }
#     taxa <- colnames(transformed)
#     taxa_transf <- transformed %>%
#         tibble::rownames_to_column("Sample_Name") %>%
#         dplyr::left_join(Metadata[c("Sample_Name", "Compartment", "Genotype", "Timepoint")], by = "Sample_Name") %>%
#         tidyr::pivot_longer(-c("Sample_Name", "Compartment", "Genotype", "Timepoint"), names_to = "Taxon", values_to = "Abundance") %>%
#         dplyr::mutate(Compartment = factor(Compartment, levels = c("Roots", "Rhizoplane", "Rhizosphere", "Soil")))
#     Boxplot_list <- llply(taxa, .fun = function(x) {
#         taxa_transf %>%
#             dplyr::filter(Taxon == x) %>%
#             ggplot(aes(x = factor(Genotype, levels = c("WT", "BI1")), y = Abundance, color = Compartment)) +
#             geom_boxplot(outlier.shape = 9, outlier.size = 3, outlier.colour = "black", width = .5) +
#             geom_point(aes(shape = Genotype), alpha = .7, size = 3) +
#             facet_grid(Timepoint ~ Compartment) +
#             scale_color_manual(values = Colors) +
#             theme_bw() +
#             guides(x = guide_axis(angle = 45)) +
#             labs(x = "Genotype", y = "Relative abundance") +
#             ggtitle(x) +
#             scale_x_discrete(labels = c("Rhizosphere" = "Rhizosphere", "Rhizoplane" = "Rhizoplane", "soil" = "Soil", "roots" = "Roots"))
#     })
#     names(Boxplot_list) <- taxa
#     multi_page <- marrangeGrob(grobs = Boxplot_list, ncol = 2, nrow = 3)
#     ggsave(multi_page, filename = paste0(outdir, "/Relative_abundance_differences_ALL.pdf"), height = 25, width = 25)
# }
RelAbundanceBoxPlots <- function(TransformedList, Metadata, Colors, TaxaOI = "all", TaxLevel) {
    transformed <- TransformedList[[TaxLevel]]
    taxa_transf <- transformed %>%
        tibble::rownames_to_column("Sample_Name") %>%
        dplyr::left_join(Metadata[c("Sample_Name", "Compartment", "Genotype", "Timepoint")], by = "Sample_Name") %>%
        tidyr::pivot_longer(-c("Sample_Name", "Compartment", "Genotype", "Timepoint"), names_to = "Taxon", values_to = "Abundance") %>%
        dplyr::mutate(Compartment = factor(Compartment, levels = c("Roots", "Rhizoplane", "Rhizosphere", "Soil")))
    kw_test <- taxa_transf %>%
        group_by(Compartment, Timepoint, Taxon) %>%
        group_map(~ KW_wrapper(.y, .x$Abundance, .x$Genotype)) %>%
        bind_rows()
    taxa_transf <- taxa_transf %>%
        dplyr::left_join(kw_test, by = c("Compartment", "Timepoint", "Taxon"))
    Compartments <- unique(taxa_transf$Compartment)
    for (i in 1:length(Compartments)) {
        Comp <- Compartments[i]
        data <- dplyr::filter(taxa_transf, Compartment == Comp)
        ## select taxa based on whether or not they have at least one significant value in the table (taxa_transf)
        if (TaxaOI != "all") {
            transformed <- transformed %>%
                dplyr::select(any_of(TaxaOI))
        } else {
            cli_alert_info("Filtering out taxa for which no differences in relative abundances were found...")
            tax_data <- data %>%
                dplyr::filter(p_value < 0.05) #%>%
            taxa <- tax_data %>% 
                dplyr::select(Taxon) %>%
                unique()
            # if I would literally want ALL taxa
            # taxa <- colnames(transformed)
        }
        if (nrow(taxa) == 0) {
            cli_alert_danger("No taxa with significant differences in relative abundances were found at the {TaxLevel} level in {Comp}. Skipping.")
            next
        }
        Boxplot_list <- llply(taxa$Taxon, .fun = function(x) {
            data <- dplyr::filter(data, Taxon == x)
            p <- data %>%
                ggplot(aes(x = factor(Genotype, levels = c("WT", "BI1")), y = Abundance, color = Compartment)) +
                geom_boxplot(outlier.shape = 9, outlier.size = 3, outlier.colour = "black", width = .5) +
                geom_point(aes(shape = Genotype), alpha = .7, size = 3) +
                facet_grid(Timepoint ~ Compartment) +
                scale_color_manual(values = Colors) +
                theme_bw() +
                theme(
                    strip.text = element_text(size = 20), legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                    axis.title = element_text(size = 20), axis.text = element_text(size = 16)
                ) +
                guides(x = guide_axis(angle = 45)) +
                labs(x = "Genotype", y = "Relative abundance") +
                ggtitle(x) +
                scale_x_discrete(labels = c("Rhizosphere" = "Rhizosphere", "Rhizoplane" = "Rhizoplane", "soil" = "Soil", "roots" = "Roots"))
            ymax <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2]
            text_df <- data %>%
                dplyr::select("Compartment", "Timepoint", "p_value") %>%
                unique()
            p <- p + annotate("segment", x = .7, xend = 2.3, y = ymax + ymax / 4, yend = ymax + ymax / 4) +
                geom_text(data = text_df, mapping = aes(x = 1.5, y = ymax + ymax / 3, label = round(p_value, 3)), color = "black")
            return(p)
        })
        names(Boxplot_list) <- taxa
        multi_page <- marrangeGrob(grobs = Boxplot_list, ncol = 3, nrow = 3)
        ggsave(multi_page, filename = paste0(outdir, "/Relative_abundance_differences_", TaxLevel, "_", Comp, ".pdf"), height = 25, width = 25)
    }
    return(Boxplot_list)
    # separate based on the compartment, and for each compartment show only significant taxa
    # Boxplot_list <- llply(taxa$Taxon, .fun = function(x) {
    #     data <- dplyr::filter(taxa_transf, Taxon == x)
    #     p <- data %>%
    #         ggplot(aes(x = factor(Genotype, levels = c("WT", "BI1")), y = Abundance, color = Compartment)) +
    #         geom_boxplot(outlier.shape = 9, outlier.size = 3, outlier.colour = "black", width = .5) +
    #         geom_point(aes(shape = Genotype), alpha = .7, size = 3) +
    #         facet_grid(Timepoint ~ Compartment) +
    #         scale_color_manual(values = Colors) +
    #         theme_bw() +
    #         guides(x = guide_axis(angle = 45)) +
    #         labs(x = "Genotype", y = "Relative abundance") +
    #         ggtitle(x) +
    #         scale_x_discrete(labels = c("Rhizosphere" = "Rhizosphere", "Rhizoplane" = "Rhizoplane", "soil" = "Soil", "roots" = "Roots"))
    #     ymax <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2]
    #     text_df <- data %>%
    #         dplyr::select("Compartment", "Timepoint", "p_value") %>%
    #         unique()
    #     p <- p + annotate("segment", x = .7, xend = 2.3, y = ymax + ymax / 4, yend = ymax + ymax / 4) +
    #         geom_text(data = text_df, mapping = aes(x = 1.5, y = ymax + ymax / 3, label = round(p_value, 3)), color = "black")
    #     return(p)
    # })
    # names(Boxplot_list) <- taxa
    # multi_page <- marrangeGrob(grobs = Boxplot_list, ncol = 2, nrow = 3)
    # ggsave(multi_page, filename = paste0(outdir, "/Relative_abundance_differences_", TaxLevel, "_", Comp ".pdf"), height = 25, width = 25)
}

AlluvialPlot <- function(TransformedList, Metadata) {
    transformed <- TransformedList$phylum
    #transformed <- transformed_list$family
    taxa_transf <- transformed %>%
        tibble::rownames_to_column("Sample_Name") %>%
        dplyr::left_join(metadata[c("Sample_Name", "Sample_ID", "Compartment", "Genotype", "Timepoint")], by = "Sample_Name") %>%
        tidyr::pivot_longer(-c("Sample_Name", "Sample_ID", "Compartment", "Genotype", "Timepoint"), names_to = "Taxon", values_to = "Abundance") %>%
        dplyr::mutate(Compartment = factor(Compartment, levels = c("Roots", "Rhizoplane", "Rhizosphere", "Soil"))) %>%
        dplyr::mutate(Taxon = as.factor(Taxon))
    taxa_alluvial <- dplyr::select(taxa_transf, -c(Sample_Name, Sample_ID)) %>%
        group_by(Compartment, Genotype, Timepoint, Taxon) %>%
        mutate(Abundance = mean(Abundance)) %>%
        unique() %>%
        ungroup() %>%
        group_by(Compartment, Genotype, Timepoint) %>%
        mutate(Sum = sum(Abundance)) %>%
        ungroup() %>%
        group_by(Taxon) %>%
        mutate(
            id = cur_group_id(),
            taxon_mean = mean(Abundance),
            group = paste0(Compartment, "_", Genotype, "_", str_remove(Timepoint, " dpi"))
        ) %>%
        mutate(group = factor(group, levels = c(
            "Roots_WT_14", "Roots_WT_28", "Roots_BI1_14", "Roots_BI1_28", "Rhizosphere_WT_14", "Rhizosphere_WT_28",
            "Rhizosphere_BI1_14", "Rhizosphere_BI1_28", "Soil_WT_14", "Soil_WT_28", "Soil_BI1_14", "Soil_BI1_28"
        ), labels = c("WT_14", "WT_28", "BI1_14", "BI1_28", "WT_14", "WT_28", "BI1_14", "BI1_28", "WT_14", "WT_28", "BI1_14", "BI1_28"))) %>%
        ungroup() %>%
        ## replace taxa with mean abundance < 0,02 (2%) together into "Other"
        mutate(Taxon = ifelse(taxon_mean < 0.02, "Other", as.character(Taxon)))

    alluv_plot <- ggplot(taxa_alluvial, aes(x = group, y = Abundance, stratum = Taxon, fill = Taxon, alluvium = id)) +
        facet_wrap(~Compartment, scales = "free_x") +
        geom_stratum(linetype = "blank", decreasing = FALSE) +
        geom_alluvium(decreasing = FALSE) +
        geom_flow(decreasing = FALSE) +
        theme_bw() +
        guides(x = guide_axis(angle = 45)) +
        theme(strip.text = element_text(size = 30), legend.title = element_text(size = 26), legend.text = element_text(size = 24),
        axis.title = element_text(size = 30), axis.text = element_text(size = 26))+
        scale_fill_viridis_d()
    return(alluv_plot)
}
