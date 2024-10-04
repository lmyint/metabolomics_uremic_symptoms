################################################################
## Functions for plotting
################################################################

#' Create a blank plot
blank_plot <- function() {
    plot(1, type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "", main = "")
}

#' Create a scatterplot with LOESS fit
plot_scatter_lowess <- function(x, y, ...) {
    plot(x, y, ...)
    keep <- !is.na(x) & !is.na(y)
    x_subs <- x[keep]
    y_subs <- y[keep]
    lines(lowess(x_subs, y_subs), lwd = 3, col = "deeppink")
}

#' Create an MA plot
plot_ma <- function(x1, x2, lowess = TRUE, ...) {
    A <- (x1+x2)/2
    M <- x2 - x1
    keep <- !is.na(A) & !is.na(M)
    A_subs <- A[keep]
    M_subs <- M[keep]
    plot(A_subs, M_subs, xlab = "Average", ylab = "Difference", ...)
    abline(h = 0, col = "red", lwd = 2)
    if (lowess) {
        lines(lowess(A_subs, M_subs), col = "deepskyblue", lwd = 2)
    }
}




################################################################
## Functions for reading and cleaning data
################################################################

#' Read in metabolite data and store as a SummarizedExperiment
#' @param file file path
#' @param range_abund cell range for abundance data
#' @param range_row_data cell range for row (metabolite) metadata
#' @param range_col_data cell range for col (sample) metadata
#' @param injection String description of chrom. method (C8, HILIC-pos, HILIC-neg)
read_broad_data <- function(file, range_abund, range_row_data, range_col_data, injection) {
    ## Quantitative abundance data
    abund <- read_excel(file, sheet = "Sheet1", range = range_abund, col_names = FALSE)
    abund <- as.matrix(abund)
    rownames(abund) <- NULL
    
    ## Metadata
    ## Row data on metabolites
    row_data <- read_excel(file, sheet = "Sheet1", range = range_row_data)
    colnames(row_data) <- fix_names(colnames(row_data))
    
    ## Column data on samples
    col_data <- read_excel(file, sheet = "Sheet2", range = range_col_data, col_names = TRUE, guess_max = 5000)
    colnames(col_data) <- fix_names(colnames(col_data))
    
    ## Spacer columns/rows create all NAs
    bool_remove <- rowAlls(is.na(col_data))
    
    ## Subset and update row/colnames
    col_data <- col_data[!bool_remove,,drop=FALSE]
    abund <- abund[,!bool_remove,drop=FALSE]
    colnames(abund) <- col_data$sample_id
    
    se <- SummarizedExperiment(
        assays = list(abund = abund, log_abund = log2(abund+1)),
        rowData = row_data,
        colData = col_data
    )
    rownames(se) <- paste0(injection, "_metab", seq_len(nrow(se)))
    colnames(se) <- col_data$sample_id
    se
}

#' Fix variable names in a character vector by changing
#' spaces to underscores and setting letters to lowercase
fix_names <- function(x) {
    x <- tolower(x)
    x <- str_replace_all(x, " ", "_")
    str_replace_all(x, "\\?", "")
}

#' Split a SummarizedExperiment object into 3 parts:
#' main study samples, PREFA QC samples, PREFB QC samples
split_experiment <- function(se) {
    cd <- colData(se)
    if ("sample_id" %in% colnames(cd)) {
        sids <- cd$sample_id
    } else if ("sample_id.x" %in% colnames(cd)) {
        sids <- cd$sample_id.x
    }
    idx_qc_a <- str_detect(sids, "PREFA") %>% which()
    idx_qc_b <- str_detect(sids, "PREFB") %>% which()
    idx_lucid_pool <- str_detect(sids, "Lucid-pool") %>% which()
    idx_main <- setdiff(seq_len(nrow(cd)), c(idx_qc_a, idx_qc_b, idx_lucid_pool))
    list(main = se[,idx_main], qc_a = se[,idx_qc_a], qc_b = se[,idx_qc_b], lucid_pool = se[,idx_lucid_pool])
}

check_ids <- function(ids1, ids2, ids3) {
    all_ids <- unique(c(ids1, ids2, ids3))
    
    ids_info <- tibble(
        id = all_ids,
        in1 = all_ids %in% ids1,
        in2 = all_ids %in% ids2,
        in3 = all_ids %in% ids3,
        in_all = in1 & in2 & in3
    )
    
    ids_info
}

merge_redcap_kdqol <- function(se, redcap, kdqol) {
    cd <- colData(se) %>% as.data.frame()
    
    ## Make event names uniform: use "Baseline" and "Year 1"
    cd <- cd %>%
        mutate(
            event_name = case_when(
                event_name=="Event 1" ~ "Baseline",
                event_name %in% c("Followup Year 1", "Year1") ~ "Year 1",
                TRUE ~ event_name
            )
        )
    redcap <- redcap %>%
        mutate(
            event_name = case_when(
                event_name=="Event 1" ~ "Baseline",
                event_name=="Followup Year 1" ~ "Year 1",
                TRUE ~ event_name
            )
        )
    kdqol <- kdqol %>%
        mutate(
            visit = case_when(
                visit=="1 year" ~ "Year 1",
                TRUE ~ visit
            )
        )
    
    ## Subset columns and fix column names
    redcap <- redcap %>%
        select(!starts_with("complete"), !contains("consent"))
    redcap_time_fixed <- redcap %>% select(proton_id, age, gender, race, ethnicity) %>% unique()
    redcap_time_vary <- redcap %>% select(-starts_with("age"), -c(gender, race, ethnicity))
    colnames(kdqol) <- str_remove_all(colnames(kdqol), "[\\?\\(\\)]")
    kdqol <- kdqol %>%
        select(!starts_with("complete"))
    
    ## Merge Broad metadata with Redcap and KDQOL info
    cd <- cd %>%
        left_join(redcap_time_fixed, by = "proton_id") %>%
        left_join(redcap_time_vary, by = c("proton_id", "event_name")) %>%
        left_join(kdqol, by = c("proton_id" = "proton_id", "event_name" = "visit"))
    
    ## Add helper metadata
    sids <- cd$sample_id.x
    dup_sids <- sids[duplicated(sids)]
    cd <- cd %>%
        mutate(
            country = ifelse(proton_id=="NA", "Canada", "US"),
            is_tech_dup = sample_id.x %in% dup_sids
        )
    cd <- as(cd, "DataFrame")
    
    colData(se) <- cd
    se
}

#' Check sample sizes and ID duplications
check_sample_sizes_dups <- function(se, se_split) {
    cat("Full SummarizedExperiment:", dim(se), "\n")
    cat("Main study samples:", dim(se_split$main), "\n")
    cat("PREFA QC samples:", dim(se_split$qc_a), "\n")
    cat("PREFB QC samples:", dim(se_split$qc_b), "\n")
    
    cat("Duplicated sample IDs (technical duplicates):\n")
    sids <- colData(se_split$main)$sample_id.x
    dup_sids <- sids[duplicated(sids)]
    cd <- colData(se_split$main) %>% as.data.frame()
    cd %>%
        dplyr::filter(sample_id.x %in% dup_sids) %>%
        select(starts_with("sample_id"), "proton_id", "event_name") %>%
        dplyr::arrange(proton_id, event_name) %>%
        print()
    
    cat("Sample sizes for discovery (event 1/baseline) and validation (Year 1):\n")
    table(cd$event_name, cd$country) %>% print()
    
    cat("Cases for which proton ID matches but sample ID doesn't:\n")
    colData(se_split$main) %>%
        as.data.frame() %>%
        mutate(
            sample_id_tail_x = as.integer(str_extract(sample_id.x, "[0-9]+$")),
            sample_id_tail_y = as.integer(str_extract(sample_id.y, "[0-9]+$"))
        ) %>%
        select(starts_with("sample_id"), proton_id, event_name) %>%
        dplyr::filter(sample_id_tail_x != sample_id_tail_y) %>%
        print()
}

#' Update `site` variable. Change NAs to Canada.
#' Relevant for LUCID but not for FAIR
update_site_variable <- function(se) {
    cd <- colData(se)
    cd$site <- ifelse(is.na(cd$site), "Canada", cd$site)
    colData(se) <- cd
    se
}

#' Filter metabolites
#' Remove >= 95% missingness per metab and var == 0
#' @param se SummarizedExperiment containing abundance data
filter_metabs <- function(se) {
    rd <- rowData(se)
    log_abund <- assay(se, "log_abund")
    is_high_miss <- rowMeans(is.na(log_abund)) >= 0.95
    is_low_var <- rowVars(as.matrix(log_abund), na.rm = TRUE) < 1e-6
    
    keep <- !is_high_miss & !is_low_var
    se_subs <- se[keep,]
    cat(
        "Original # metabs:", nrow(se),
        "\n# metabs with high missingness:", sum(is_high_miss),
        "\n# metabs with (near) zero variance:", sum(is_low_var, na.rm = TRUE),
        "\nNew # metabs:", nrow(se_subs), "\n"
    )
    numbers_by_metab_type <- tibble(
        known = ifelse(is.na(rd$metabolite), "unknown", "known"),
        high_miss = is_high_miss,
        low_var = is_low_var
    )
    print(dplyr::count(numbers_by_metab_type, known))
    print(dplyr::count(numbers_by_metab_type, known, high_miss, low_var))
    se_subs
}

#' Remove redundant metabolites
#' DO NOT remove internal standards yet
remove_metabolites <- function(se, df_redundancy_info, injection = c("C8-pos", "HIL-pos", "HIL-neg"), row_data_col = c("metabolite", "former_nomenclature")) {
    ## Get metabolites that should be omitted
    metabs_omit <- df_redundancy_info %>%
        filter(method==injection, can_omit=="X", !str_detect(hmdb_id, "[Ii]nternal")) %>%
        pull(metabolite)
    rd <- rowData(se) %>% as.data.frame()
    bool_omit <- rd[[row_data_col]] %in% metabs_omit
    se_subs <- se[!bool_omit,]
    cat("Orig # metabs:", nrow(se), "\nNew # metabs:", nrow(se_subs), "\n")
    se_subs
}

remove_internal_standards <- function(se) {
    rd <- rowData(se)
    bool_omit <- str_detect(tolower(rd$hmdb_id), "internal") & !is.na(rd$hmdb_id)
    se_subs <- se[!bool_omit,]
    cat("Orig # metabs:", nrow(se), "\nNew # metabs:", nrow(se_subs), "\n")
    se_subs
}

#' Combine SummarizedExperiments across all 3 injections
#' into one SummarizedExperiment
#' @param se_list A list of se_METHOD_main objects (PREFA/B removed)
merge_ses_across_injections <- function(se_list) {
    sample_id_colname <- ifelse("sample_id.x" %in% colnames(colData(se_list[[1]])), "sample_id.x", "sample_id")
    
    ## Align samples to match ordering in first injection
    sample_ids <- colData(se_list[[1]])[[sample_id_colname]]
    reorder_idx2 <- match(sample_ids, colData(se_list[[2]])[[sample_id_colname]])
    reorder_idx3 <- match(sample_ids, colData(se_list[[3]])[[sample_id_colname]])
    se_list_reordered <- list(se_list[[1]], se_list[[2]][,reorder_idx2], se_list[[3]][,reorder_idx3])
    stopifnot(identical(colData(se_list_reordered[[1]])[[sample_id_colname]], colData(se_list_reordered[[2]])[[sample_id_colname]]))
    stopifnot(identical(colData(se_list_reordered[[2]])[[sample_id_colname]], colData(se_list_reordered[[3]])[[sample_id_colname]]))
    
    ## Remove colData from 2nd and 3rd injections to allow rbind-ing
    colData(se_list_reordered[[2]]) <- NULL
    colData(se_list_reordered[[3]]) <- NULL
    se_all <- rbind(se_list_reordered[[1]], se_list_reordered[[2]], se_list_reordered[[3]])
    se_all
}



################################################################
## Functions relating to symptoms
################################################################

get_symptom_definitions <- function() {
    tibble(
        short_descrip = c("fatigue", "pruritus", "anorexia", "nausea_vomiting", "daytime_sleepiness", "difficulty_concentrating", "bodily_pain"),
        long_descrip = c("washed_out_or_drained", "itchy_skin", "lack_of_appetite", "nausea_or_upset_stomach", "have_trouble_staying_awake_during_the_day", "did_you_have_difficulty_concentrating_or_thinking", "how_much_bodily_pain_have_you_had_during_the_past_4_weeks")
    )
}

#' Coarsen a 5 or 6 level symptom variable to 
coarsen_symptom <- function(y, num_new_levels = c(2,3)) {
    num_levels <- length(setdiff(y, NA))
    if (num_new_levels==3) {
        if (num_levels==5) {
            new_y <- case_when(
                y=="grade1" ~ "level1",
                y %in% c("grade2", "grade3") ~ "level2",
                y %in% c("grade4", "grade5") ~ "level3"
            )
        } else if (num_levels==6) {
            new_y = case_when(
                y=="grade1" ~ "level1",
                y %in% c("grade2", "grade3", "grade4") ~ "level2",
                y %in% c("grade5", "grade6") ~ "level3"
            )
        }
    } else if (num_new_levels==2) {
        if (num_levels==5) {
            new_y <- case_when(
                y %in% c("grade1", "grade2") ~ "level1",
                y %in% c("grade3", "grade4", "grade5") ~ "level2"
            )
        } else if (num_levels==6) {
            new_y <- case_when(
                y %in% c("grade1", "grade2", "grade3") ~ "level1",
                y %in% c("grade4", "grade5", "grade6") ~ "level2"
            )
        }
    }
    new_y
}
