value_counts <- function(
    data,
    columns = NULL,
    normalize = FALSE,
    na.rm = FALSE,
    sort = TRUE,
    decreasing = TRUE
) {
    if (is.null(columns)) {
        if (!is.vector(data) && !is.factor(data)) {
            stop(
                "Input 'data' must be a vector or factor if 'columns' is NULL."
            )
        }
        counts_table <- table(data, useNA = ifelse(na.rm, "no", "ifany"))
        counts_df <- as.data.frame(counts_table)
        colnames(counts_df)[ncol(counts_df)] <- "n"
        if (length(dim(counts_table)) == 1) {
            colnames(counts_df)[1] <- names(dimnames(counts_table))
        }
    } else {
        if (!is.data.frame(data) && !is.matrix(data)) {
            stop(
                "Input 'data' must be a data.frame or matrix if 'columns' are specified."
            )
        }
        if (!all(columns %in% colnames(data))) {
            stop("One or more specified columns not found in 'data'.")
        }

        selected_data <- data[, columns, drop = FALSE]

        # Use table directly on the selected columns
        counts_table <- do.call(
            table,
            c(as.list(selected_data), useNA = ifelse(na.rm, "no", "ifany"))
        )
        counts_df <- as.data.frame.table(counts_table, responseName = "n")
    }

    # Filter out rows with zero counts
    counts_df <- counts_df[counts_df$n > 0, ]

    if (normalize) {
        counts_df$proportion <- counts_df$n / sum(counts_df$n)
    }

    if (sort) {
        counts_df <- counts_df[order(counts_df$n, decreasing = decreasing), ]
    }

    rownames(counts_df) <- NULL
    return(counts_df)
}

ensure_dir <- function(dir_path, verbose = FALSE) {
    if (!file.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
        if (verbose) message("Directory created: ", dir_path)
    } else {
        if (verbose) message("Directory already exists: ", dir_path)
    }
}


check_column_alignment <- function(adata) {
    x <- adata$X
    obs <- adata$obs
    var <- adata$var

    # Extract column names (excluding first column) and associated obs labels
    x_colnames <- colnames(x)[-1]
    obs_names <- as.character(obs[[1]])

    # Extract row identifiers from first column and var
    x_varorder <- x[[1]]
    var_order <- var[[1]]

    col_match <- identical(x_colnames, obs_names)
    row_match <- identical(x_varorder, var_order)

    if (col_match && row_match) {
        message(
            "Column names and row identifiers match and are in the same order."
        )
        return(TRUE)
    }
}


get_pal <- function(df, names, color_code) {
    df <- df %>% dplyr::select(names, color_code) %>% dplyr::distinct()
    stats::setNames(df[[color_code]], df[[names]])
}

assign_met_colors <- function(
    df,
    column,
    palette = "VanGogh2",
    reverse = FALSE,
    type = "discrete"
) {
    #assign_met_colors(adata$obs, "cell_type_name", palette = "Klimt", type = 'continuous')
    vals <- unique(df[[column]])
    n_vals <- length(vals)
    pal <- MetBrewer::met.brewer(name = palette, n = n_vals, type = type)
    if (reverse) {
        pal <- rev(pal)
    }
    color_map <- stats::setNames(pal, vals)
    color_column <- paste0(column, "_color")
    df[[color_column]] <- color_map[df[[column]]]
    attr(df, paste0(column, "_colormap")) <- color_map
    df
}
