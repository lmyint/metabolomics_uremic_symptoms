############################################################
## Functions for missing data explorations
############################################################

missing_data_metab <- function(se) {
    rd <- rowData(se) %>% as.data.frame()
    log_abund <- assay(se, "log_abund")
    num_missing <- rowSums(is.na(log_abund))
    perc_missing <- rowMeans(is.na(log_abund))
    mean_abund <- rowMeans(log_abund, na.rm = TRUE)
    df_miss <- rd %>%
        mutate(
            num_missing = num_missing,
            num_missing_cat = case_when(
                num_missing==0 ~ "0",
                num_missing >= 1 & num_missing <= 10 ~ "1-10",
                num_missing >= 11 & num_missing <= 50 ~ "11-50",
                num_missing >= 50 & num_missing != ncol(se) ~ "50-~All",
                num_missing == ncol(se) ~ "All"
            ),
            perc_missing = perc_missing,
            metab_known = !is.na(metabolite),
            mean_abund = mean_abund
        )
    p1 <- ggplot(df_miss, aes(x = perc_missing)) + 
        geom_density()
    p2 <- ggplot(df_miss, aes(x = metab_known, fill = num_missing_cat)) +
        geom_bar(position = "fill") +
        labs(x = "Metabolite identity known?", y = "# samples missing") +
        theme_classic()
    p3 <- ggplot(df_miss, aes(x = mean_abund, y = perc_missing)) +
        geom_point() +
        geom_smooth() +
        labs(x = "Mean metabolite abundance", y = "% missingness") +
        theme_classic()
    print(p1)
    print(p2)
    print(p3)
}

missing_data_samples <- function(se) {
    log_abund <- assay(se, "log_abund")
    perc_missing <- colMeans(is.na(log_abund))
    print(summary(perc_missing))
    plot(density(perc_missing, from = 0, to = 1), xlab = "% metabolites with missing abundance", main = "")
}




############################################################
## Functions for investigating and implementing imputation
############################################################

#' @param dat The data (a data.frame)
#' @param cols A character vector of the columns in which data will be removed. If NULL, all columns will be used.
#' @param num Number of cases per column to make NA. Can be a single integer, in which case, the same number of cases are deleted in each column. Or can be an integer vector of the same length as `cols`, indicating the number of cases to delete in each column.
#' @return The original data with missing values introduced
create_missingness <- function(dat, cols = NULL, num, method = c("mcar", "left_cens")) {
    method <- match.arg(method)
    
    ## If cols is NULL, all columns will be used
    if (is.null(cols)) {
        cols <- colnames(dat)
    }
    
    ## Argument checking on num
    if (length(num)!=1 & length(num)!=length(cols)) {
        stop("length(num) is not equal to 1 or length(cols)")
    }
    if (length(num)==1) {
        num <- rep(num, length(cols))
    }
    
    ## Loop over desired columns
    for (i in seq_along(cols)) {
        col <- cols[i]
        n <- num[i]
        
        ## Get rows of non-NA entries
        idx_not_na <- which(!is.na(dat[[col]]))
        
        ## Update number to make missing in case of too few non-NAs
        n_actual <- min(n, length(idx_not_na))
        
        ## If all values are missing, move to next variable
        if (length(idx_not_na)==0)
            next
        
        if (method=="mcar") {
            ## Randomly select among the non-missing observations
            idx_make_na <- sample(idx_not_na, size = n_actual)
        } else if (method=="left_cens") {
            ## Select the smallest observations
            o <- order(dat[[col]])
            idx_make_na <- head(o, n_actual)
        }
        ## Make values NA
        dat[[col]][idx_make_na] <- NA
    }
    
    dat
}

get_imputation_error <- function(data_truth, data_forced_na, data_imputed, col_types) {
    ## In case some variables were dropped,
    ## obtain variables that are in both truth and imputed
    common_cols <- intersect(colnames(data_truth), colnames(data_imputed))
    
    ## All columns are of the same type
    if (length(col_types)==1) {
        col_types <- rep(col_types, length(common_cols))
    }
    
    ## Loop over variables
    errs <- sapply(seq_along(common_cols), function(i) {
        col <- common_cols[i]
        col_type <- col_types[i]
        truth <- data_truth[[col]]
        forced_na <- data_forced_na[[col]]
        guess <- data_imputed[[col]]
        
        ## Subset to only places where NAs were manually inserted
        bool <- !is.na(truth) & is.na(forced_na)
        
        ## Compute error (RMSE and misclassification rate)
        if (col_type=="q") {
            err <- mean((truth[bool]-guess[bool])^2) %>% sqrt()
        } else {
            err <- mean(truth[bool] != guess[bool])*100
        }
        err
    })
    names(errs) <- common_cols
    
    errs
}

impute_half_min <- function(data) {
    for (j in seq_len(ncol(data))) {
        bool_miss <- is.na(data[,j])
        
        if (!any(bool_miss))
            next
        
        ## Generate random values from 0.5*min - min
        min_obs <- min(data[,j], na.rm = TRUE)
        rand_vals <- runif(sum(bool_miss), min = 0.5*min_obs, max = min_obs)
        
        ## Impute with randomly-generated values
        data[bool_miss,j] <- rand_vals
    }
    
    data
}

evaluate_imp_methods_metab <- function(se) {
    ## Transpose data to be samples x metabs
    log_abund_orig <- assay(se, "log_abund") %>% t() %>% as.data.frame()
    num_metabs <- ncol(log_abund_orig)
    
    ## Create missingness in columns without too many NAs already
    bool_cols_na <- colMeans(is.na(log_abund_orig)) < 0.6
    names_cols_na <- colnames(log_abund_orig)[bool_cols_na]
    log_abund_na_lc <- create_missingness(log_abund_orig, num = 20, cols = names_cols_na, method = "left_cens")
    log_abund_na_lc_t <- t(log_abund_na_lc)
    
    ## QRILC on samples x metabs
    imp_qrilc1 <- impute.QRILC(log_abund_na_lc)
    ## QRILC on metabs x samples
    imp_qrilc2 <- impute.QRILC(log_abund_na_lc_t)
    ## KNN on metabs x samples
    imp_knn5 <- impute.knn(log_abund_na_lc_t, k = 5, rowmax = 0.5, colmax = 0.8, maxp = num_metabs)
    imp_knn10 <- impute.knn(log_abund_na_lc_t, k = 10, rowmax = 0.5, colmax = 0.8, maxp = num_metabs)
    imp_knn20 <- impute.knn(log_abund_na_lc_t, k = 20, rowmax = 0.5, colmax = 0.8, maxp = num_metabs)
    imp_knn10_rm70 <- impute.knn(log_abund_na_lc_t, k = 10, rowmax = 0.7, colmax = 0.8, maxp = num_metabs)
    imp_knn10_rm90 <- impute.knn(log_abund_na_lc_t, k = 10, rowmax = 0.9, colmax = 0.8, maxp = num_metabs)
    ## Half-min imputation on samples x metabs
    imp_half_min <- impute_half_min(log_abund_na_lc)
    
    ## Store imputed datasets in a list
    ## Ensure orientation for all inputs: samples x metabs
    imp_list <- list(
        qrilc1 = as.data.frame(imp_qrilc1[[1]]),
        qrilc2 = as.data.frame(t(imp_qrilc2[[1]])),
        knn5 = as.data.frame(t(imp_knn5[["data"]])),
        knn10 = as.data.frame(t(imp_knn10[["data"]])),
        knn20 = as.data.frame(t(imp_knn20[["data"]])),
        knn10_r70 = as.data.frame(t(imp_knn10_rm70[["data"]])),
        knn10_r90 = as.data.frame(t(imp_knn10_rm90[["data"]])),
        half_min = imp_half_min
    )
    
    ## Evaluate imputation accuracy
    imp_errors <- lapply(imp_list, function(imp_data) {
        err_vec <- get_imputation_error(
            data_truth = log_abund_orig,
            data_forced_na = log_abund_na_lc,
            data_imputed = imp_data,
            col_types = "q"
        )
        tibble(
            compound_id = names(err_vec),
            error = err_vec
        )
    })
    imp_errors <- bind_rows(imp_errors, .id = "imp_method")
    list(orig = log_abund_orig, forced_na = log_abund_na_lc, imputed = imp_list, errors = imp_errors)
}

impute_metab <- function(se) {
    ## Transpose data to be samples x metabs
    log_abund <- assay(se, "log_abund") %>% t() %>% as.data.frame()
    
    ## QRILC on samples x metabs
    imp_qrilc <- impute.QRILC(log_abund)[[1]] %>% t() %>% as.matrix()
    
    ## Half-min imputation on samples x metabs
    imp_half_min <- impute_half_min(log_abund) %>% t() %>% as.matrix()
    
    ## Zero imputation
    imp_zero <- log_abund
    imp_zero[is.na(imp_zero)] <- 0
    imp_zero <- imp_zero %>% t() %>% as.matrix()
    
    ## Store imputed datasets
    assay(se, "log_abund_imp_qrilc") <- imp_qrilc
    assay(se, "log_abund_imp_half_min") <- imp_half_min
    assay(se, "log_abund_imp_zero") <- imp_zero
    se
}





