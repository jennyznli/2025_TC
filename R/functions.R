### HELPER FUNCTIONS #### 

# Selects top n probes with highest standard deviation across samples
# Args: betas - numeric matrix, n - number of top variable probes to select
# Returns: subset matrix
bSubMostVariable <- function(betas, n=3000) {
    print(paste("Original dimensions:", nrow(betas), "x", ncol(betas)))
    std <- apply(betas, 1, sd, na.rm=TRUE)
    result <- betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
    print(paste("New dimensions:", nrow(result), "x", ncol(result)))
    return(result)
}

# Prepares enrichment dataframe for plotting filtered by FDR and minimum counts
# Args: df - enrichment results, n_min/n_max - min/max entries, max_fdr - FDR cutoff
# Returns: filtered and annotated dataframe ready for plotting
preparePlotFDR <- function(df, n_min, n_max, max_fdr) {
    df <- df[df$nD >0,]
    df$FDR[df$FDR==0] <- .Machine$double.xmin
    df <- df[order(df$FDR),]
    if (sum(df$FDR < max_fdr) < n_min) {
        df1 <- head(df, n=n_min)
    } else {
        df_filtered <- df[df$estimate > 0,] # enrichment only, exclude depletion
        df_filtered_sig <- df_filtered[df_filtered$FDR < max_fdr,]
        if (nrow(df_filtered_sig) < n_min) {
            df1 <- head(df_filtered, n=n_min)  # Take top n_min from filtered set
        } else {
            df1 <- head(df_filtered_sig, n=n_max)
        }
    }
    gp <- vapply(str_split(df1$group, "\\."), function(x) x[3], character(1))
    if ("Target" %in% colnames(df1)) {
        df1$db1 <- sprintf("%s: %s (%s)", gp, df1$Target, df1$dbname)
    } else {
        df1$db1 <- sprintf("%s: %s", gp, df1$dbname)
    }
    df1$db1 <- factor(df1$db1, levels=rev(df1$db1))
    df1$dbname <- factor(df1$dbname, levels=rev(df1$dbname))
    df1
}

# Plots dotplot of enrichment results showing -log10(FDR) vs term with coloring by estimate
# Args: df - enrichment results, n_min/n_max - min/max points to plot, max_fdr - cutoff
# Returns: ggplot object
plotDotFDR <- function(df, n_min = 10, n_max = 10, max_fdr = 0.05) {
    db1 <- FDR <- overlap <- estimate <- NULL
    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))
    df1 <- preparePlotFDR(df, n_min, n_max, max_fdr)
    ggplot(df1) +
        geom_point(aes(-log10(FDR), dbname, size=overlap, color=estimate)) +
        scale_color_gradientn(colors = brewer.pal(9, "Purples")[4:9]) +
        xlab("-log10(FDR)") + ylab("") +
        theme_minimal()
}
