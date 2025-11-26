gene_category_pal = c(
    "Not detected" = "#bebebe",
    "Detected in all" = "#253494",
    "Detected in many" = "#2c7fb8",
    "Detected in some" = "#41b6c4",
    "Detected in single" = "#a1dab4",
    "Low tissue specificity" = "#666666",
    "Tissue enhanced" = "#984ea3",
    "Group enriched" = "#FF9D00",
    "Tissue enriched" = "#e41a1c"
)


tf_palette <- c(
    "beta-Barrel DNA-binding domains" = "#55ad89",
    "beta-Hairpin exposed by an alpha/beta-scaffold" = "#c3bc3f",
    "beta-Sheet binding to DNA" = "#8cc2ca",
    "alpha-Helices exposed by beta-structure" = "#F9d23c",
    "Helix-turn-helix domains" = "#bb7693",
    "Zinc-coordinating DNA-binding domains" = "#6388b4",
    "Basic domains" = "#ffae34",
    "Immunoglobulin fold" = "#ef6f6a",
    "Other all-alpha-helical DNA-binding domains" = "#cfcfcf",
    "Yet undefined DNA-binding domains" = "#414451"
)


calculate_tau_score <- function(X) {
    # Check if X is a tibble, where columns are samples, and rows are the genes.
    if (!tibble::is_tibble(X)) {
        stop("Input X must be a tibble.")
    }

    # Transform numeric columns using log1p transformation
    wide_data <- X %>%
        dplyr::mutate_if(is.numeric, function(x) log1p(x)) %>%
        tibble::column_to_rownames(colnames(.)[1])

    # Calculate maximum expression for each gene across all samples
    max_exp <- apply(wide_data, MARGIN = 1, function(x) max(x, na.rm = TRUE))

    # Calculate number of non-missing values for each gene
    N <- apply(wide_data, MARGIN = 1, function(x) length(which(!is.na(x))))

    # Calculate expression sum for each gene
    expression_sum <- wide_data %>%
        sweep(MARGIN = 1, STATS = max_exp, FUN = `/`) %>%
        {
            1 - .
        } %>%
        apply(MARGIN = 1, function(x) sum(x, na.rm = TRUE))

    # Calculate tau score for each gene and convert to tibble
    tau_score <- (expression_sum / (N - 1)) %>%
        tibble::enframe("gene", "tau_score")

    # Return the tau scores
    return(tau_score)
}
hpa_gene_classification <- function(
    AnnDatR,
    enr_fold = 4,
    max_group_n,
    det_lim = 1,
    inplace = TRUE
) {
    X <- AnnDatR$X
    tissue_col <- AnnDatR$obs_names_col
    gene_col <- AnnDatR$var_names_col

    data <- X %>%
        tidyr::pivot_longer(
            -1,
            names_to = 'tissue_col',
            values_to = 'expression'
        )

    # from https://github.com/maxkarlsson/Pig-Atlas/blob/91be1cb4e49aba2e4a27f313a70e887007765a51/scripts/functions_classification.R
    data_ <-
        data %>%
        dplyr::select(
            gene = gene_col,
            expression = expression,
            tissue = tissue_col
        ) %>%
        dplyr::mutate(expression = round(expression, 4))

    if (any(is.na(data_$expression))) {
        stop("NAs in expression column")
    }
    if (any(is.na(data_$gene))) {
        stop("NAs in gene column")
    }
    if (any(is.na(data_$tissue))) {
        stop("NAs in tissue column")
    }

    n_groups <- length(unique(data_$tissue))

    gene_class_info <-
        data_ %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(
            # Gene expression distribution metrics
            mean_exp = mean(expression, na.rm = T),
            min_exp = min(expression, na.rm = T),
            max_exp = max(expression, na.rm = T),
            max_2nd = sort(expression)[length(expression) - 1],

            # Expression frequency metrics
            n_exp = length(which(expression >= det_lim)),
            frac_exp = n_exp / length(expression[!is.na(expression)]) * 100,

            # Limit of enhancement metrics
            lim = max_exp / enr_fold,

            exps_over_lim = list(expression[which(
                expression >= lim & expression >= det_lim
            )]),
            n_over = length(exps_over_lim[[1]]),
            mean_over = mean(exps_over_lim[[1]]),
            min_over = ifelse(n_over == 0, NA, min(exps_over_lim[[1]])),

            max_under_lim = max(
                expression[which(expression < min_over)],
                det_lim * 0.1
            ),

            exps_enhanced = list(which(
                expression / mean_exp >= enr_fold & expression >= det_lim
            )),

            # Expression patterns
            enrichment_group = paste(
                sort(tissue[which(expression >= lim & expression >= det_lim)]),
                collapse = ";"
            ),

            n_enriched = length(tissue[which(
                expression >= lim & expression >= det_lim
            )]),
            n_enhanced = length(exps_enhanced[[1]]),
            enhanced_in = paste(
                sort(tissue[exps_enhanced[[1]]]),
                collapse = ";"
            ),
            n_na = n_groups - length(expression),
            max_2nd_or_lim = max(max_2nd, det_lim * 0.1),
            tissues_not_detected = paste(
                sort(tissue[which(expression < det_lim)]),
                collapse = ";"
            ),
            tissues_detected = paste(
                sort(tissue[which(expression >= det_lim)]),
                collapse = ";"
            )
        )

    gene_categories <-
        gene_class_info %>%

        dplyr::mutate(
            spec_category = dplyr::case_when(
                n_exp == 0 ~ "not detected",

                # Genes with expression fold times more than anything else are tissue enriched
                max_exp / max_2nd_or_lim >= enr_fold ~ "tissue enriched",

                # Genes with expression fold times more than other tissues in groups of max group_n - 1 are group enriched
                max_exp >= lim &
                    n_over <= max_group_n &
                    n_over > 1 &
                    mean_over / max_under_lim >= enr_fold ~ "group enriched",

                # Genes with expression in tissues fold times more than the mean are tissue enhance
                n_enhanced > 0 ~ "tissue enhanced",

                # Genes expressed with low tissue specificity
                T ~ "low tissue specificity"
            ),

            dist_category = dplyr::case_when(
                frac_exp == 100 ~ "detected in all",
                frac_exp >= 31 ~ "detected in many",
                n_exp > 1 ~ "detected in some",
                n_exp == 1 ~ "detected in single",
                n_exp == 0 ~ "not detected"
            ),

            spec_score = dplyr::case_when(
                spec_category == "tissue enriched" ~ max_exp / max_2nd_or_lim,
                spec_category == "group enriched" ~ mean_over / max_under_lim,
                spec_category == "tissue enhanced" ~ max_exp / mean_exp
            )
        )

    ##### Rename and format
    gene_classification <- gene_categories %>%
        dplyr::mutate(
            enriched_tissues = dplyr::case_when(
                spec_category %in%
                    c("tissue enriched", "group enriched") ~ enrichment_group,
                spec_category == "tissue enhanced" ~ enhanced_in
            ),
            n_enriched = dplyr::case_when(
                spec_category %in%
                    c("tissue enriched", "group enriched") ~ n_enriched,
                spec_category == "tissue enhanced" ~ n_enhanced
            )
        ) %>%
        dplyr::select(
            gene,
            spec_category,
            dist_category,
            spec_score,
            n_expressed = n_exp,
            fraction_expressed = frac_exp,
            max_exp = max_exp,
            enriched_tissues,
            n_enriched,
            n_na = n_na,
            tissues_not_detected,
            tissues_detected
        )

    gene_classification <- gene_classification %>%
        dplyr::left_join(calculate_tau_score(X), by = 'gene') %>%
        dplyr::mutate(
            tau_score = ifelse(spec_category == "not detected", NA, tau_score)
        )

    if (inplace) {
        AnnDatR$var <- AnnDatR$var %>%
            dplyr::left_join(
                gene_classification,
                by = stats::setNames('gene', gene_col)
            )
    }

    return(gene_classification)
}
