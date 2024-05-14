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
KW_wrapper <- function(keys, alpha, comparison) {
    if (length(unique(comparison)) > 1) {
        result <- kruskal.test(alpha, comparison)
        tibble(
            keys,
            p_value = result$p.value
        )
    }
}