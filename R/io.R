load_scanpy_data <- function(
    prefix_name,
    var_names = 'ensembl_ids',
    layer = NULL,
    file_dir = './'
) {
    # Check if file directory exists
    if (!dir.exists(file_dir)) {
        stop("The specified directory does not exist.")
    }

    # Initialize an empty list to store data
    data <- list()

    # Read data based on whether layer is specified
    if (is.null(layer)) {
        file_name <- paste0(prefix_name, '.tsv')
        data[['X']] <- readr::read_delim(
            file.path(file_dir, file_name),
            delim = '\t'
        ) %>%
            dplyr::arrange(!!rlang::sym(var_names))
    } else {
        file_name <- paste0(prefix_name, '_', layer, '.tsv')
        data[[layer]] <- readr::read_delim(
            file.path(file_dir, file_name),
            delim = '\t'
        ) %>%
            dplyr::arrange(!!rlang::sym(var_names))
    }

    # Read metadata files
    data$var <- readr::read_tsv(file.path(
        file_dir,
        paste0(prefix_name, '_var.tsv')
    )) %>%
        dplyr::arrange(!!rlang::sym(var_names))
    data$obs <- readr::read_tsv(file.path(
        file_dir,
        paste0(prefix_name, '_obs.tsv')
    ))

    # Return the data
    return(data)
}
