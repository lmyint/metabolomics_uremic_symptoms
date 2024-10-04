plot_qc_drift <- function(se, sample_descrip, include_subtitle_corr = FALSE) {
    if(is.unsorted(colData(se)$injection_order)) {
        o <- order(colData(se)$injection_order)
        se <- se[,o]
    }
    log_abund <- assay(se, "log_abund")
    log_abund_sample1 <- log_abund[,1]
    ## Looping so that entire correlation matrix 
    ## doesn't need to be computed
    cor_with_first <- sapply(2:ncol(log_abund), function(i) {
        cor(log_abund_sample1, log_abund[,i], use = "complete.obs")
    })
    cor_sequential <- sapply(2:ncol(log_abund), function(i) {
        cor(log_abund[,i-1], log_abund[,i], use = "complete.obs")
    })
    
    ## Plot correlations
    title1 <- sample_descrip
    title2 <- sample_descrip
    if (include_subtitle_corr) {
        title1 <- paste0(title1, ":\nCorrelation between each sample and the FIRST")
        title2 <- paste0(title2, ":\nCorrelation between SUBSEQUENT samples")
    }
    
    par(mfrow = c(1,2))
    plot(cor_with_first, xlab = "Injection order", ylab = "Correlation", main = title1)
    lines(lowess(x = seq_along(cor_with_first), y = cor_with_first), lwd = 2, col = "deeppink")
    # axis(side = 1, at = seq_along(cor_with_first), labels = tail(colnames(log_abund), -1), las = 3)
    
    plot(cor_sequential, xlab = "Injection order", ylab = "Correlation", main = title2)
    lines(lowess(x = seq_along(cor_sequential), y = cor_sequential), lwd = 2, col = "deeppink")
    # axis(side = 1, at = seq_along(cor_sequential), labels = tail(colnames(log_abund), -1), las = 3)
    
    ## MA plots
    par(mfrow = c(1,1))
    index_other <- min(c(10, ncol(log_abund)))
    plot_ma(log_abund_sample1, log_abund[,index_other], main = paste("Samples 1 and", index_other))
    
    ## Summary stats on correlations
    cat("Summary statistics: correlations with FIRST sample\n")
    print(summary(cor_with_first))
    cat("Summary statistics: correlations with SUBSEQUENT samples\n")
    print(summary(cor_sequential))
}

plot_cv_dists <- function(se_split, sample_descrip, xlim) {
    se_main <- se_split$main
    se_prefa <- se_split$qc_a
    se_prefb <- se_split$qc_b
    
    get_cv <- function(se) {
        log_met <- assay(se, "log_abund")
        log_met[is.na(log_met)] <- 0
        cvs <- rowSds(log_met, na.rm = TRUE)/rowMeans(log_met, na.rm = TRUE)
        cvs
    }
    
    drop_nan_inf <- function(x) {
        x[!is.na(x) & !is.nan(x) & !is.infinite(x)]
    }
    
    idx_internal_std <- which(str_detect(rowData(se_main)$hmdb_id %>% tolower(), "internal"))
    idx_non_std <- setdiff(seq_len(nrow(se_main)), idx_internal_std)
    
    idx_internal_std_prefa <- which(str_detect(rowData(se_prefa)$hmdb_id %>% tolower(), "internal"))
    idx_non_std_prefa <- setdiff(seq_len(nrow(se_prefa)), idx_internal_std_prefa)
    
    idx_internal_std_prefb <- which(str_detect(rowData(se_prefb)$hmdb_id %>% tolower(), "internal"))
    idx_non_std_prefb <- setdiff(seq_len(nrow(se_prefb)), idx_internal_std_prefb)
    
    cv_main <- get_cv(se_main)
    cv_prefa <- get_cv(se_prefa)
    cv_prefb <- get_cv(se_prefb)
    
    cv_main_stds <- cv_main[idx_internal_std] %>% drop_nan_inf()
    cv_main <- cv_main[idx_non_std] %>% drop_nan_inf()
    cv_prefa <- cv_prefa[idx_non_std_prefa] %>% drop_nan_inf()
    cv_prefb <- cv_prefb[idx_non_std_prefb] %>% drop_nan_inf()
    
    title <- paste0(sample_descrip, ": CV distribution")
    dens_prefb <- density(cv_prefb, from = 0)
    dens_main <- density(cv_main, from = 0)
    plot(dens_prefb, xlab = "Coefficient of variation", main = title, lwd = 3, xlim = xlim)
    lines(dens_main, lwd = 3, col = "deeppink")
    abline(v = cv_main_stds, lwd = 3, col = "dodgerblue")
    legend("topright", legend = c("PREFA", "PREFB", "LUCID study", "Internal standard"), bty = "n", lwd = 3, col = c("deeppink", "lightpink", "black", "dodgerblue"))
    
    ## Summary statistics for CVs
    cat("CVs for internal standards:", cv_main_stds, "\n")
    cat("Summary stats for CVs for study samples:\n")
    print(summary(cv_main))
    cat("Summary stats for CVs for QC samples (PREFA):\n")
    print(summary(cv_prefa))
    cat("Summary stats for CVs for QC samples (PREFB):\n")
    print(summary(cv_prefb))
}


plot_pc_info <- function(se, sample_descrip = "", which_log_abund = NULL) {
    if (is.null(which_log_abund)) {
        abund <- assay(se, "abund")
        abund[is.na(abund)] <- 0
        log_abund <- log2(abund + 1)
    } else {
        log_abund <- assay(se, which_log_abund)
    }
    
    pc_out <- prcomp(log_abund, center = TRUE, scale = TRUE)
    
    injection_order_cut <- colData(se)$injection_order %>% cut(breaks = 8)
    par(mfrow = c(1,2))
    ## Score plot
    plot(pc_out$x[,"PC1"], pc_out$x[,"PC2"], col = injection_order_cut, xlab = "PC1", ylab = "PC2", main = sample_descrip, ylim = range(pc_out$x[,"PC2"]) * c(1, 1.5))
    legend("top", legend = 1:8, col = 1:8, pch = 1, horiz = TRUE, bty = "n")
    ## Scree plot
    pc_vars <- pc_out$sdev^2
    var_explained <- pc_vars/sum(pc_vars)
    cum_var_explained <- cumsum(var_explained)
    plot(var_explained[1:20], pch = 16, ylim = c(0,1), xlab = "PC", ylab = "(Cumulative) % variance explained")
    lines(cum_var_explained[1:20], pch = 16, col = "deeppink")
    legend("right", legend = c("% var. explained", "Cumulative"), col = c("black", "deeppink"), pch = 16, bty = "n")
}

plot_technical_dups <- function(se, proton_id, ...) {
    cd <- colData(se)
    bool <- cd$proton_id==proton_id
    se_subs <- se[,bool]
    cd <- colData(se_subs)
    index_first <- which.min(cd$injection_order)
    index_second <- which.max(cd$injection_order)
    log_abund <- assay(se_subs, "log_abund")
    x1 <- log_abund[,index_first]
    x2 <- log_abund[,index_second]
    plot_ma(x1 = x1, x2 = x2, main = paste("Proton ID:", proton_id), ...)
    cat("Pearson correlation:", cor(x1, x2, method = "pearson", use = "complete.obs"), "\n")
    cat("Spearman correlation:", cor(x1, x2, method = "spearman", use = "complete.obs"), "\n")
}

