### Get all orientations of primers
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer) # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(
        Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
        RevComp = Biostrings::reverseComplement(dna)
    )
    return(sapply(orients, toString)) # Convert back to character vector
}

### Count the number of times our primers appear in our fwd and rev reads, considering all primer orientations
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

#### Wrap dada2 step (without the merge)
dada2_wrap <- function(inpath, filenames = list(), maxEE, truncQ, truncLen, trimLeft, trimRight, minLen) {
    if (length(filenames) == 0) {
        cli_abort("No files given for dada2 analysis")
    }
    else if (length(filenames) == 1) {
        cli_alert_info("Performing dada2 denoising on 'single-end' reads (probably coming from Flash2 merging).")
        ## Define output filenames
        filtered_files <- file.path(inpath, "filtered", basename(filenames[[1]]))
        ## Filter the reads. Default parameters make it so that filtering is only done to remove all reads containing Ns and minLen=50 bases
        filtered <- as.data.frame(filterAndTrim(filenames[[1]], filtered_files,
            maxN = 0, maxEE = maxEE, truncQ = truncQ, truncLen = truncLen,
            trimLeft = trimLeft, trimRight = trimRight, minLen = minLen,
            rm.phix = TRUE, compress = TRUE, multithread = TRUE
        ))
        ## Filter out files that have less than 500 reads after N filtering:
        filt_df <- filtered[filtered$reads.out > 500,]
        colnames(filt_df)[length(colnames(filt_df))] <- "non_chimeric"
        print("filtered table where reads.out have > 500 reads")
        print(filt_df)
        r1s <- row.names(filt_df)
        filtered_files <- file.path(inpath, "filtered", r1s)
        ## Learning error rates
        errors <- learnErrors(filtered_files, multithread = TRUE)
        pdf(file = paste0(outpath, "/errors_estimated.pdf"))
        plotErrors(errors, nominalQ = TRUE)
        dev.off()
        ## Sample inference
        dada <- dada(filtered_files, err = errors, multithread = TRUE)
        ## Return
        return(list(filt_df, filtered, filtered_files, dada))
    }
    else if (length(filenames) == 2) {
        cli_alert_info("Performing dada2 denoising on 'paired-end' reads (not yet merged).")
        ## Define output filenames (after real filtering)
        filtered_F <- file.path(inpath, "filtered", basename(filenames[[1]]))
        filtered_R <- file.path(inpath, "filtered", basename(filenames[[2]]))
        filtered <- as.data.frame(filterAndTrim(filenames[[1]], filtered_F, filenames[[2]], filtered_R,
            maxN = 0, maxEE = maxEE, truncQ = truncQ, truncLen = truncLen,
            trimLeft = trimLeft, trimRight = trimRight, minLen = minLen,
            rm.phix = FALSE, compress = TRUE, multithread = TRUE
        ))
        ## Filter out files that have less than 500 reads after N filtering:
        filt_df <- filtered[filtered$reads.out > 500,]
        colnames(filt_df)[length(colnames(filt_df))] <- "non_chimeric"
        print("filtered table where reads.out have > 500 reads")
        print(filt_df)
        r1s <- row.names(filt_df)
        r2s <- sub("R1", "R2", row.names(filt_df))
        filtered_F <- file.path(inpath, "filtered", r1s)
        filtered_R <- file.path(inpath, "filtered", r2s)
        ## Learning error rates
        errors_F <- learnErrors(filtered_F, multithread = TRUE)
        errors_R <- learnErrors(filtered_R, multithread = TRUE)
        pdf(file = paste0(outpath, "/errors_estimated.pdf"))
        plotErrors(errors_F, nominalQ = TRUE)
        plotErrors(errors_R, nominalQ = TRUE)
        dev.off()
        ## Sample inference
        dada_F <- dada(filtered_F, err = errors_F, multithread = TRUE)
        dada_R <- dada(filtered_R, err = errors_R, multithread = TRUE)
        ## Return
        return(list(filt_df, filtered, filtered_F, filtered_R, dada_F, dada_R))
    }
    else { 
        cli_abort("Unsupported number of input files (> 2)")
    }
}

getN <- function(x) sum(getUniques(x))

draw_stat_plot <- function(denoising_stats){
    denoise_percentages <- denoising_stats %>%
        dplyr::mutate(across(everything()), (. / before_filter) * 100) %>%
        rownames_to_column("Sample") %>%
        tidyr::pivot_longer(cols = -Sample, names_to = "Time") %>%
        dplyr::mutate(Time = factor(Time, levels = c("before_filter", "after_filter", "denoised_F", "denoised_R", "merged", "non_chimeric")))
    
    denoising_plot <- ggplot(data = denoise_percentages, aes(x = Time, y = value, group = Sample)) +
        geom_line() +
        geom_point()
    return(denoising_plot)
}

make_stat_table <- function(denoising_stats, merging = TRUE){
    if (merging == TRUE){
        denoise_percentages <- denoising_stats %>%
            dplyr::mutate(perc_passed_filter = (after_filter/before_filter)*100) %>% 
            dplyr::mutate(perc_merged = (merged/before_filter)*100) %>% 
            dplyr::mutate(perc_non_chimeric = (non_chimeric/before_filter)*100)
    } else {
        denoise_percentages <- denoising_stats %>%
            dplyr::mutate(perc_passed_filter = (after_filter/before_filter)*100) %>% 
            dplyr::mutate(perc_non_chimeric = (non_chimeric/before_filter)*100)
    }
    print("Denoise percentages...")
    print(denoise_percentages)
    return(denoise_percentages)
}