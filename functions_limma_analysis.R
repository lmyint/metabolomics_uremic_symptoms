###################################################################
## Functions needed to conduct the main analysis: LUCID
###################################################################

#' Convert a numeric vector into percentiles
percentile <- function(x) {
    rank(x, ties.method = "min")*100/length(x)
}

#' Add QC column and its percentile version to row_data
#' @param x the QC variable
#' @param col_name the name of the QC variable
#' @return the same row_data supplemented with 2 additional columns
.add_qc_metrics <- function(row_data, x, col_name) {
    perc_col_name <- paste0("perc_", col_name)
    row_data[[col_name]] <- x
    row_data[[perc_col_name]] <- percentile(x)
    row_data
}

#' Run limma analysis for each symptom
run_analysis <- function(se, mod_form, event = c("Baseline", "Year 1", "FAIR"), coarse_symptom_levels = c(NA, 2, 3), metab_type = c("orig", "qrilc_no_ruv", "qrilc_ruv"), use_imputed_symptoms = FALSE) {
    event <- match.arg(event)
    metab_type <- match.arg(metab_type)
    
    ## Extract colData for creating design matrix later
    col_data <- colData(se)
    
    ## Subset SummExp to desired event (Baseline/Year 1)
    if (event != "FAIR") {
        bool_keep <- col_data$event_name==event
        se <- se[,bool_keep]
    }
    row_data <- rowData(se) %>% as.data.frame()
    col_data <- colData(se) %>% as.data.frame()
    
    ## Regardless of metab_type, QC info is for non-imputed data
    ## Add QC info about each metabolite to give context for top hits
    ## % missing, mean, SD, CV AND percentile versions
    log_abund <- assay(se, "log_abund") ## Non-imputed data
    mean_abund <- rowMeans(log_abund, na.rm = TRUE)
    sd_abund <- rowSds(log_abund, na.rm = TRUE)
    cv_abund <- sd_abund/mean_abund
    perc_missing <- rowMeans(is.na(log_abund))
    
    row_data <- .add_qc_metrics(row_data, mean_abund, "mean_abund")
    row_data <- .add_qc_metrics(row_data, sd_abund, "sd_abund")
    row_data <- .add_qc_metrics(row_data, cv_abund, "cv_abund")
    row_data <- .add_qc_metrics(row_data, perc_missing, "perc_missing")
    row_data$metab_id <- rownames(row_data)
    rowData(se) <- row_data
    
    ## Extract indicated metabolomics data
    if (metab_type=="orig") {
        log_abund <- assay(se, "log_abund")
    } else if (metab_type=="qrilc_no_ruv") {
        log_abund <- assay(se, "log_abund_imp_qrilc")
    } else if (metab_type=="qrilc_ruv") {
        log_abund <- assay(se, "log_abund_ruv_qrilc")
    }

    old_opt <- options()$na.action
    options(na.action = "na.pass")
    
    ## Loop over symptom numbers
    symptom_df <- get_symptom_definitions()
    
    results <- lapply(seq_len(nrow(symptom_df)), function(i) {
        ## Create `symptom` variable in col_data
        ## (a copy of the column for this symptom)
        if (use_imputed_symptoms) {
            this_symptom_name <- paste0(symptom_df$short_descrip[i], "_imp")
        } else {
            this_symptom_name <- symptom_df$short_descrip[i]
        }
        col_data$symptom <- col_data[[this_symptom_name]]
        cat("    ", this_symptom_name, "\t\t")
        
        ## Coarsen outcome variable (if specified)
        if (!is.na(coarse_symptom_levels)) {
            col_data$symptom <- coarsen_symptom(col_data$symptom, num_new_levels = coarse_symptom_levels)
        }
        col_data$symptom <- as.factor(col_data$symptom)
        
        ## Create design matrix
        design <- model.matrix(mod_form, data = col_data)
        colnames(design)[str_detect(colnames(design), "Intercept")] <- "intercept"
        
        ## Remove samples with NAs in `design`
        samples_keep <- complete.cases(design)
        design <- design[samples_keep,,drop=FALSE]
        log_abund <- log_abund[, samples_keep, drop=FALSE]
        
        cat("n =", sum(samples_keep), "\n")
        
        if (sum(samples_keep)==0) {
            stop("No samples remain")
        }
        
        ## Run limma
        fit <- lmFit(log_abund, design = design)
        fit <- eBayes(fit)
        
        ## Select symptom coefficients for F test
        cn_design <- colnames(design)
        select_coeffs <- cn_design[str_detect(cn_design, "^symptom")]
        tt <- topTable(fit, coef = select_coeffs, number = Inf)
        tt$metab_id <- rownames(tt)
        
        res <- tt %>%
            left_join(row_data, by = "metab_id") %>%
            mutate(across(contains("hmdb"), as.character))
        res_known <- res %>%
            filter(!is.na(metabolite)) %>%
            arrange(P.Value) %>%
            select(-adj.P.Val) %>%
            mutate(adj_pval = p.adjust(P.Value, method = "BH"))
        res_unknown <- res %>%
            filter(is.na(metabolite)) %>%
            arrange(P.Value) %>%
            select(-adj.P.Val) %>%
            mutate(adj_pval = p.adjust(P.Value, method = "BH"))
        res <- bind_rows(res_known, res_unknown)
        attr(res, "n") <- sum(samples_keep)
        res
    })
    names(results) <- symptom_df$short_descrip
    
    ## Reset NA options
    options(na.action = old_opt)
    
    results
}

run_all_analyses <- function(se_list, mod_form, cohort = c("LUCID", "FAIR"), coarse_symptom_levels = c(NA, 2, 3), metab_type = c("orig", "qrilc_no_ruv", "qrilc_ruv"), use_imputed_symptoms = FALSE) {
    metab_type <- match.arg(metab_type)
    cohort <- match.arg(cohort)
    
    ## Create one SE that contains all injections
    se_all <- merge_ses_across_injections(se_list)
    
    if (cohort=="LUCID") {
        event_names <- c("Baseline", "Year 1")
    } else if (cohort=="FAIR") {
        event_names <- "FAIR"
    }
    list_by_event <- lapply(event_names, function(event) {
        cat(event, ":\n")
        run_analysis(
            se_all,
            mod_form = mod_form,
            event = event,
            coarse_symptom_levels = coarse_symptom_levels,
            metab_type = metab_type,
            use_imputed_symptoms = use_imputed_symptoms
        )
    })
    names(list_by_event) <- event_names
    list_by_event
}

#' Convert output from `run_all_analyses` to data frame
results_list_to_df <- function(results) {
    ## Combine lists of results into one data frame
    res_df <- lapply(results, function(results_this_event) {
        bind_rows(results_this_event, .id = "symptom")
    }) %>% bind_rows(.id = "event")
    
    ## Add a metab_id variable that gives the metabolite name 
    ## for known metabs or the compound ID for unknown metabs
    res_df %>%
        mutate(metab_id = ifelse(is.na(metabolite), compound_id, metabolite))
}

#' Run limma analysis for mortality outcomes
run_mort_analysis <- function(se_list, mod_form, metab_type = c("orig", "qrilc_no_ruv", "qrilc_ruv"), drug_metabs) {
    metab_type <- match.arg(metab_type)
    
    ## Create one SE that contains all injections
    se <- merge_ses_across_injections(se_list)
    
    ## Subset SummExp to Baseline samples from Frenova and Canada
    col_data <- colData(se) %>% as.data.frame()
    cat("Subsetting to Baseline samples from Frenova and Canada:\n")
    bool_keep <- col_data$event_name=="Baseline" & (col_data$country=="Canada" | col_data$frenova)
    se <- se[,bool_keep]
    cat("    n =", ncol(se), "\n")
    
    ## Extract row and colData
    row_data <- rowData(se) %>% as.data.frame()
    col_data <- colData(se) %>% as.data.frame()
    row_data$metab_id <- rownames(row_data)
    rowData(se) <- row_data
    
    ## Add indoxylsulfate to colData
    index_indox <- which(row_data$metabolite=="indoxylsulfate")
    col_data$indox_log_qrilc <- as.numeric(assay(se[index_indox,], "log_abund_imp_qrilc"))
    col_data$indox_log_qrilc_ruv <- as.numeric(assay(se[index_indox,], "log_abund_ruv_qrilc"))
    
    ## Extract indicated metabolomics data
    if (metab_type=="orig") {
        log_abund <- assay(se, "log_abund")
    } else if (metab_type=="qrilc_no_ruv") {
        log_abund <- assay(se, "log_abund_imp_qrilc")
    } else if (metab_type=="qrilc_ruv") {
        log_abund <- assay(se, "log_abund_ruv_qrilc")
    }
    
    old_opt <- options()$na.action
    options(na.action = "na.pass")
    
    ## Loop over 1-year and 2-year mortality definitions
    days_mort_thresh <- 365*c(1,2)
    
    results <- lapply(days_mort_thresh, function(thresh) {
        ## Create `mort_outcome` variable in col_data
        col_data$mort_outcome <- col_data$daysdeath <= thresh & !is.na(col_data$daysdeath)
        
        ## Create design matrix
        design <- model.matrix(mod_form, data = col_data)
        colnames(design)[str_detect(colnames(design), "Intercept")] <- "intercept"
        
        ## Remove samples with NAs in `design`
        samples_keep <- complete.cases(design)
        design <- design[samples_keep,,drop=FALSE]
        log_abund <- log_abund[, samples_keep, drop=FALSE]
        
        cat("n =", sum(samples_keep), "\n")
        
        if (sum(samples_keep)==0) {
            stop("No samples remain")
        }
        
        ## Run limma
        fit <- lmFit(log_abund, design = design)
        fit <- eBayes(fit)
        
        ## Get all coefficients from 
        df_fit_coeffs <- fit$coefficients %>% 
            as.data.frame()
        colnames(df_fit_coeffs)[colnames(df_fit_coeffs)=="(Intercept)"] <- "intercept"
        colnames(df_fit_coeffs) <- paste0("mod_coeff_", colnames(df_fit_coeffs))
        df_fit_coeffs <- df_fit_coeffs %>% 
            mutate(metab_id = rownames(fit$coefficients))
        rownames(df_fit_coeffs) <- NULL
        
        ## Select mort_outcome coefficient for hypothesis testing
        tt <- topTable(fit, coef = "mort_outcomeTRUE", number = Inf)
        tt$metab_id <- rownames(tt)
        tt <- tt %>% 
            left_join(df_fit_coeffs, by = "metab_id")
        
        res <- tt %>%
            left_join(row_data, by = "metab_id") %>%
            mutate(across(contains("hmdb"), as.character))
        res <- res %>%
            filter(!(metabolite %in% drug_metabs)) # Remove drug metabs
        res_known <- res %>%
            filter(!is.na(metabolite)) %>%
            arrange(P.Value) %>%
            select(-adj.P.Val) %>%
            mutate(adj_pval = p.adjust(P.Value, method = "BH"))
        res_unknown <- res %>%
            filter(is.na(metabolite)) %>%
            arrange(P.Value) %>%
            select(-adj.P.Val) %>%
            mutate(adj_pval = p.adjust(P.Value, method = "BH"))
        res <- bind_rows(res_known, res_unknown)
        res
    })
    names(results) <- c("1_year_mort", "2_year_mort")
    
    ## Reset NA options
    options(na.action = old_opt)
    
    results
}


### Ordinary linear regression
run_linreg_analysis <- function(se, mod_form, event = c("Baseline", "Year 1", "FAIR"), metab_type = c("orig", "qrilc_no_ruv", "qrilc_ruv"), use_imputed_symptoms = FALSE) {
    event <- match.arg(event)
    metab_type <- match.arg(metab_type)
    
    ## Extract colData for creating design matrix later
    col_data <- colData(se)
    
    ## Subset SummExp to desired event (Baseline/Year 1)
    if (event != "FAIR") {
        bool_keep <- col_data$event_name==event
        se <- se[,bool_keep]
    }
    
    ## Subset to just known metabolites
    # known <- !is.na(rowData(se)$metabolite)
    # se <- se[known,]
    
    row_data <- rowData(se) %>% as.data.frame()
    col_data <- colData(se) %>% as.data.frame()
    
    ## Regardless of metab_type, QC info is for non-imputed data
    ## Add QC info about each metabolite to give context for top hits
    ## % missing, mean, SD, CV AND percentile versions
    log_abund <- assay(se, "log_abund") ## Non-imputed data
    mean_abund <- rowMeans(log_abund, na.rm = TRUE)
    sd_abund <- rowSds(log_abund, na.rm = TRUE)
    cv_abund <- sd_abund/mean_abund
    perc_missing <- rowMeans(is.na(log_abund))
    
    row_data <- .add_qc_metrics(row_data, mean_abund, "mean_abund")
    row_data <- .add_qc_metrics(row_data, sd_abund, "sd_abund")
    row_data <- .add_qc_metrics(row_data, cv_abund, "cv_abund")
    row_data <- .add_qc_metrics(row_data, perc_missing, "perc_missing")
    row_data$metab_id <- rownames(row_data)
    rowData(se) <- row_data
    
    ## Extract indicated metabolomics data
    if (metab_type=="orig") {
        log_abund <- assay(se, "log_abund")
    } else if (metab_type=="qrilc_no_ruv") {
        log_abund <- assay(se, "log_abund_imp_qrilc")
    } else if (metab_type=="qrilc_ruv") {
        log_abund <- assay(se, "log_abund_ruv_qrilc")
    }
    
    ## Loop over symptom numbers
    symptom_df <- get_symptom_definitions()
    
    results <- lapply(seq_len(nrow(symptom_df)), function(i) {
        ## Create `symptom` variable in col_data
        ## (a copy of the column for this symptom)
        if (use_imputed_symptoms) {
            this_symptom_name <- paste0(symptom_df$short_descrip[i], "_imp")
        } else {
            this_symptom_name <- symptom_df$short_descrip[i]
        }
        col_data$symptom <- col_data[[this_symptom_name]]
        col_data$symptom <- as.integer(str_remove(col_data$symptom, "grade"))
        cat("    ", this_symptom_name, "\n")
        
        ## Loop over metabolites
        mod_list <- lapply(seq_len(nrow(log_abund)), function(r) {
            col_data$metab <- log_abund[r,]
            
            mod <- try({lm(mod_form, data = col_data)}, silent = TRUE)
            if (inherits(mod, "try-error")) {
                data.frame(estimate = NA, p.value = NA)
            } else {
                broom::tidy(mod) %>%
                    filter(term=="metab") %>%
                    select(estimate, p.value)
            }
        })
        names(mod_list) <- row_data$metab_id
        tt <- bind_rows(mod_list, .id = "metab_id")
        
        res <- tt %>%
            left_join(row_data, by = "metab_id") %>%
            mutate(across(contains("hmdb"), as.character))
        res_known <- res %>%
            filter(!is.na(metabolite)) %>%
            arrange(p.value) %>%
            mutate(adj_pval = p.adjust(p.value, method = "BH"))
        res_unknown <- res %>%
            filter(is.na(metabolite)) %>%
            arrange(p.value) %>%
            mutate(adj_pval = p.adjust(p.value, method = "BH"))
        res <- bind_rows(res_known, res_unknown)
        res
    })
    names(results) <- symptom_df$short_descrip
    
    results
}

run_all_linreg_analyses <- function(se_list, mod_form, cohort = c("LUCID", "FAIR"), metab_type = c("orig", "qrilc_no_ruv", "qrilc_ruv"), use_imputed_symptoms = FALSE) {
    metab_type <- match.arg(metab_type)
    cohort <- match.arg(cohort)
    
    ## Create one SE that contains all injections
    se_all <- merge_ses_across_injections(se_list)
    
    if (cohort=="LUCID") {
        event_names <- c("Baseline", "Year 1")
    } else if (cohort=="FAIR") {
        event_names <- "FAIR"
    }
    list_by_event <- lapply(event_names, function(event) {
        cat(event, ":\n")
        run_linreg_analysis(
            se_all,
            mod_form = mod_form,
            event = event,
            metab_type = metab_type,
            use_imputed_symptoms = use_imputed_symptoms
        )
    })
    names(list_by_event) <- event_names
    list_by_event
}


