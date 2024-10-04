############################################################
## Functions for fitting individual models
############################################################

#' @param force_in A character vector of variables to force into the model, or NULL to not force any variables in
train_lasso <- function(data, force_in = NULL, lambdas = NULL) {
    if (is.null(lambdas)) {
        power_seq <- c(seq(-3, 0, length.out = 50), seq(0.01, 0.1, length.out = 50))
        lambdas <- 10^power_seq
    }
    if (!is.null(force_in)) {
        x <- model.matrix(~ .-outcome-1, data = data)
        cn <- colnames(x)
        penalty_factor <- rep(1, ncol(x))
        for (v in force_in) {
            bool_in <- str_detect(cn, v)
            penalty_factor[bool_in] <- 0
        }
        lasso_mod <- train(
            x = x,
            y = data$outcome,
            data = data,
            method = "glmnet",
            trControl = trainControl(method = "cv", number = 10),
            tuneGrid = data.frame(alpha = 1, lambda = lambdas),
            penalty.factor = penalty_factor,
            metric = "Accuracy"
        )
    } else {
        lasso_mod <- train(
            outcome ~ .,
            data = data,
            method = "glmnet",
            trControl = trainControl(method = "cv", number = 10),
            tuneGrid = data.frame(alpha = 1, lambda = lambdas),
            metric = "Accuracy"
        )
    }
    lasso_mod
}

train_rand_forest <- function(data, mtry = NULL) {
    mtry_values <- seq(2, ncol(data), length.out = 10) %>% ceiling()
    rf_mod <- train(
        outcome ~ .,
        data = data,
        method = "rf",
        metric = "Accuracy",
        trControl = trainControl(method = "oob"),
        tuneGrid = data.frame(mtry = mtry_values)
    )
    rf_mod
}

############################################################
## Functions for obtaining variable importance measures
############################################################

var_importance <- function(caret_mod, method) {
    mod <- caret_mod$finalModel
    if (method=="glmnet") {
        # best_lambda <- mod$tuneValue$lambda
        # coeffs_best <- tail(as.numeric(coef(mod, s = best_lambda)), -1)
        # first_zeros <- sapply(seq_len(nrow(coef_mat)), function(i) {
        #     is_zero <- coef_mat[i,]==0
        #     match(TRUE, is_zero)
        # })
        # first_zeros[is.na(first_zeros)] <- max(first_zeros, na.rm = TRUE)+1
        
        # Alt way based on coefficients that are nonzero in final model
        # coef(mod, caret_mod$bestTune$lambda)
        coef_mat <- mod$beta
        if (is(coef_mat, "list")) { # 3+ category outcome
            imp_list <- lapply(coef_mat, function(mat) {
                rowSums(mat != 0)
            })
            imp <- rowSums(bind_cols(imp_list))
            vars <- rownames(coef_mat[[1]])
        } else {
            imp <- rowSums(coef_mat != 0)
            vars <- rownames(coef_mat)
        }
        
        var_imp <- tibble(var = vars, importance_lasso = imp) %>%
            arrange(desc(importance_lasso))
    } else if (method=="rf") {
        var_imp <- randomForest::importance(mod, type = 2)
        var_imp <- tibble(
            var = rownames(var_imp),
            importance_rf = var_imp[,"MeanDecreaseGini"]
        ) %>%
            arrange(desc(importance_rf))
    }
    var_imp
}

############################################################
## Functions for fitting all models
############################################################

#' Fit ML models for all symptoms of interest
#' @param se SummarizedExperiment
#' @param adjustment_vars Character vector of non-metabolite
#'                        adjustment variables
#' @param only_known_metabs Only include known (HMDB-identified) 
#'                          metabolites in the models?
#' @param which_symptoms Character vector of short_descrip's of
#'                       symptoms for which to fit models
#' @param coarse_symptom_levels Integer. Number of levels for 
#'                       coarsened version of symptoms. NA for orig.
fit_ml_models <- function(se, adjustment_vars, only_known_metabs = FALSE, which_symptoms = NULL, event = c("Baseline", "Year 1", "FAIR"), coarse_symptom_levels = c(NA, 2, 3), metab_type = c("orig", "qrilc_no_ruv", "qrilc_ruv"), use_imputed_symptoms = FALSE) {
    event <- match.arg(event)
    
    col_data <- colData(se)
    
    ## Subset SummExp to desired event (Baseline/Year 1)
    if (event != "FAIR") {
        bool_keep <- col_data$event_name==event
        se <- se[,bool_keep]
    }
    row_data <- rowData(se) %>% as.data.frame()
    col_data <- colData(se) %>% as.data.frame()
    
    if (only_known_metabs) {
        ## Keep only identified peaks (which have HMDB ID)
        bool_keep <- !is.na(row_data$hmdb_id)
        se <- se[bool_keep,]
    }
    
    ## Get symptom information and col and rowData
    symptom_df <- get_symptom_definitions()
    col_data <- colData(se)
    row_data <- rowData(se)
    
    ## Remove duplicate metabolites
    bool_dup_metabs <- duplicated(row_data$metabolite) & !is.na(row_data$hmdb_id)
    se <- se[!bool_dup_metabs,]
    row_data <- rowData(se)
    
    ## Change rownames of SummExp to be metabolite names
    rownames(se) <- row_data$metabolite
    
    ## If no specific symptoms specified, choose all of them
    if (is.null(which_symptoms)) {
        which_symptoms <- symptom_df$short_descrip
    }
    
    ## Subset to chosen symptoms
    symptom_df <- symptom_df %>%
        filter(short_descrip %in% which_symptoms)
    
    ## Extract indicated metabolomics data and transpose
    if (metab_type=="orig") {
        log_abund <- assay(se, "log_abund")
    } else if (metab_type=="qrilc_no_ruv") {
        log_abund <- assay(se, "log_abund_imp_qrilc")
    } else if (metab_type=="qrilc_ruv") {
        log_abund <- assay(se, "log_abund_ruv_qrilc")
    }
    log_abund <- t(log_abund)
    
    ## Loop over symptoms
    results <- lapply(seq_len(nrow(symptom_df)), function(i) {
        if (use_imputed_symptoms) {
            this_symptom_name <- paste0(symptom_df$short_descrip[i], "_imp")
        } else {
            this_symptom_name <- symptom_df$short_descrip[i]
        }
        cat("    ", this_symptom_name, ":\n")
        ## Create dataset for fitting ML models where all 
        ## columns are used (so that y ~ . formula works)
        dat_ml <- bind_cols(
            outcome = col_data[[this_symptom_name]],
            as.data.frame(col_data[,adjustment_vars]),
            as.data.frame(log_abund)
        )
        
        ## Coarsen outcome variable (if specified)
        if (!is.na(coarse_symptom_levels)) {
            dat_ml$outcome <- coarsen_symptom(dat_ml$outcome, num_new_levels = coarse_symptom_levels)
        }
        dat_ml$outcome <- as.factor(dat_ml$outcome)
        
        ## Remove NAs
        cc <- complete.cases(dat_ml)
        cat("        n =", sum(cc), "\n")
        dat_ml <- dat_ml[cc,]
        
        ## LASSO
        lasso_mod <- train_lasso(dat_ml)
        
        ## Random forest
        rf_mod <- train_rand_forest(dat_ml)
        
        res <- list(lasso_mod, rf_mod)
        attr(res, "n") <- sum(cc)
        res
    })
    names(results) <- symptom_df$short_descrip
    
    results
}

run_all_ml_analyses <- function(se_list, adjustment_vars = NULL, cohort = c("LUCID", "FAIR"), only_known_metabs = FALSE, coarse_symptom_levels = c(NA, 2, 3), metab_type = c("orig", "qrilc_no_ruv", "qrilc_ruv"), use_imputed_symptoms = FALSE) {
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
        cat(event, "\n")
        fit_ml_models(
            se_all,
            adjustment_vars = adjustment_vars,
            only_known_metabs = only_known_metabs,
            which_symptoms = NULL,
            event = event,
            coarse_symptom_levels = coarse_symptom_levels,
            metab_type = metab_type,
            use_imputed_symptoms = use_imputed_symptoms
        )
    })
    names(list_by_event) <- event_names
    list_by_event
}

fit_ml_models_mort <- function(se_list, adjustment_vars, only_known_metabs = FALSE, metab_type = c("orig", "qrilc_no_ruv", "qrilc_ruv")) {
    metab_type <- match.arg(metab_type)
    
    ## Create one SE that contains all injections
    se <- merge_ses_across_injections(se_list)
    
    ## Subset SummExp to Baseline samples from Frenova and Canada
    cat("Subsetting to Baseline samples from Frenova and Canada:\n")
    col_data <- colData(se) %>% as.data.frame()
    bool_keep <- col_data$event_name=="Baseline" & (col_data$country=="Canada" | col_data$frenova)
    se <- se[,bool_keep]
    cat("    n =", ncol(se), "\n")
    
    if (only_known_metabs) {
        ## Keep only identified peaks (which have HMDB ID)
        bool_keep <- !is.na(rowData(se)$hmdb_id)
        se <- se[bool_keep,]
    } else {
        bool_keep <- is.na(rowData(se)$hmdb_id)
        se <- se[bool_keep,]
    }
    row_data <- rowData(se)
    
    ## Remove duplicate metabolites
    bool_dup_metabs <- duplicated(row_data$metabolite) & !is.na(row_data$hmdb_id)
    se <- se[!bool_dup_metabs,]
    row_data <- rowData(se)
    col_data <- colData(se)
    
    ## Change rownames of SummExp to be metabolite names
    rownames(se) <- row_data$metabolite
    
    ## Extract indicated metabolomics data and transpose
    if (metab_type=="orig") {
        log_abund <- assay(se, "log_abund")
    } else if (metab_type=="qrilc_no_ruv") {
        log_abund <- assay(se, "log_abund_imp_qrilc")
    } else if (metab_type=="qrilc_ruv") {
        log_abund <- assay(se, "log_abund_ruv_qrilc")
    }
    log_abund <- t(log_abund)
    
    ## Loop over 1-year and 2-year mortality definitions
    days_mort_thresh <- 365*c(1,2)
    
    results <- lapply(days_mort_thresh, function(thresh) {
        ## Create `mort_outcome` variable in col_data
        col_data$mort_outcome <- col_data$daysdeath <= thresh & !is.na(col_data$daysdeath)
        col_data$mort_outcome <- as.factor(col_data$mort_outcome)
        
        ## Create dataset for fitting ML models where all 
        ## columns are used (so that y ~ . formula works)
        dat_ml <- bind_cols(
            outcome = col_data[["mort_outcome"]],
            as.data.frame(col_data[,adjustment_vars]),
            as.data.frame(log_abund)
        )
        
        ## Remove NAs
        cc <- complete.cases(dat_ml)
        cat("        n =", sum(cc), "\n")
        dat_ml <- dat_ml[cc,]
        
        ## LASSO
        lasso_mod <- train_lasso(dat_ml)
        
        ## Random forest
        rf_mod <- train_rand_forest(dat_ml)
        
        res <- list(lasso_mod, rf_mod)
        attr(res, "n") <- sum(cc)
        res
    })
    names(results) <- c("1_year_mort", "2_year_mort")
    
    results
}

#####################################################################
## Summarizing ML results
#####################################################################
percent_rank_names <- function(x) {
    case_when(
        x > 0.9 ~ "High",
        x > 0.7 & x <= 0.9 ~ "Medium",
        x <= 0.7 ~ "Low"
    )
}

ml_results_to_df <- function(results) {
    ml_summ <- lapply(results, function(list_by_event) {
        l <- lapply(list_by_event, function(list_by_symp) {
            lasso_mod <- list_by_symp[[1]]
            rf_mod <- list_by_symp[[2]]
            vi_lasso <- var_importance(lasso_mod, method = "glmnet")
            vi_rf <- var_importance(rf_mod, method = "rf")
            full_join(vi_lasso, vi_rf, by = "var") %>%
                mutate(
                    importance_lasso = percent_rank(importance_lasso),
                    importance_rf = percent_rank(importance_rf),
                    metab_id = str_remove_all(var, "`")
                ) %>%
                mutate(across(starts_with("importance"), percent_rank_names)) %>%
                select(-var)
        })
        bind_rows(l, .id = "symptom")
    }) %>% bind_rows(.id = "event")
    ml_summ <- pivot_wider(ml_summ, names_from = event, values_from = starts_with("importance"))
    
    ml_summ
}