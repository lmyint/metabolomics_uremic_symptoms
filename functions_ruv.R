################################################################
## Functions for RUV
################################################################

#' Explorations required before running RUV:
#' Picking empirical standards and looking at scree plot
#' for chosen negative controls (internal & empirical 
#' standards) to choose k
#' @param se SummarizedExperiment containing assay data
#' @return Vector of indices to metabolites serving as 
#'         negative controls
pre_ruv_explorations <- function(se) {
    ## Identify rows corresponding to internal standards
    rd <- rowData(se)
    idx_standards <- which(str_detect(tolower(rd$hmdb_id), "internal"))
    idx_non_standards <- setdiff(seq_len(nrow(se)), idx_standards)
    
    ## If there are less than 10 internal standards, identify
    ## empirical controls that are most positively correlated 
    ## with one or more of those standards
    num_empirical_standards <- 10 - length(idx_standards)
    log_abund <- t(assay(se, "log_abund"))
    cor_standards <- cor(
        log_abund[,idx_standards],
        log_abund[,idx_non_standards],
        use = "pairwise.complete.obs"
    )
    cor_maxs <- colMaxs(cor_standards, na.rm = TRUE)
    temp_idx_empirical_standards <- tail(order(cor_maxs), num_empirical_standards)
    cat("Max correlations between empirical standards and internal standards:\n", sort(cor_maxs[temp_idx_empirical_standards]), "\n")
    idx_empirical_standards <- idx_non_standards[temp_idx_empirical_standards]
    idx_all_standards <- c(idx_standards, idx_empirical_standards)
    
    ## Print out missingness information for these standards
    ## (Metabolites are in columns due to transposing)
    cat("Fraction of samples missing for the chosen negative controls:", round(colMeans(is.na(log_abund[,idx_all_standards])), 2), "\n")
    
    ## Look at PCA score plots of chosen standards to determine
    ## choice for k using "elbow" inspection
    se_standards <- se[idx_all_standards,]
    plot_pc_info(se_standards, sample_descrip = "Internal & empirical standards")
    
    idx_all_standards
}

#' Implement RUV-random for metabolomics data
#' @param se SummarizedExperiment containing assay data
#' @param idx_all_standards Indices of metabolites serving as 
#'                          negative controls for RUV (includes
#'                          internal and empirical standards)
implement_ruv <- function(se, which_log_abund = c("imp_zero", "imp_half_min", "imp_qrilc"), idx_all_standards, k) {
    ## Obtain indicated abundance info
    ## Transpose to format as samples x metabs for RUV
    which_log_abund <- match.arg(which_log_abund)
    which_log_abund <- paste0("log_abund_", which_log_abund)
    log_abund <- t(assay(se, which_log_abund))
    
    ## Run RUV
    NormalizeRUVRand(Y = log_abund, ctl = idx_all_standards, k = k)
}
