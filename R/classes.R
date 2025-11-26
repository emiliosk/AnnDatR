#' @title AnnDatR Class
#' @description A flexible, Anndata-like container for transcriptomics data.
#'
#' @importFrom R6 R6Class
#' @importFrom dplyr %>%
#'
#' @field X Primary expression data matrix (features x samples).
#' @field obs Cell/Sample metadata (data.frame).
#' @field var Feature/Gene metadata (data.frame).
#' @export

AnnDatR <- R6Class(
    "AnnDatR",
    public = list(
        X = NULL,
        var = NULL,
        obs = NULL,
        var_names_col = NULL,
        layers = list(),
        obs_names_col = NULL,
        uns = list(),
        obsm = list(),
        raw = NULL,

        initialize = function(
            prefix_name = NULL,
            var_names = "ensembl_ids",
            layer = NULL,
            file_dir = "./",
            X = NULL,
            var = NULL,
            obs = NULL
        ) {
            if (!is.null(X) && !is.null(var) && !is.null(obs)) {
                # Direct assignment
                self$X <- X
                self$var <- var
                self$obs <- obs
                self$var_names_col <- var_names
                # Attempt to guess obs_names_col if not provided or standard
                self$obs_names_col <- colnames(self$obs)[1]
            } else if (!is.null(prefix_name)) {
                # Load data
                data <- load_scanpy_data(
                    prefix_name,
                    var_names,
                    layer,
                    file_dir
                )

                if (is.null(layer)) {
                    self$X <- data[['X']]
                } else {
                    self$X <- data[[layer]]
                    self$layers[[layer]] <- data[[layer]]
                }

                self$var <- data$var
                self$obs <- data$obs
                self$var_names_col <- var_names
                self$obs_names_col <- colnames(self$obs)[1]
            } else {
                stop(
                    "Provide either 'prefix_name' to load data or 'X', 'var', 'obs' for direct assignment."
                )
            }

            # Initialize lists if null
            if (is.null(self$uns)) {
                self$uns <- list()
            }
            if (is.null(self$obsm)) {
                self$obsm <- list()
            }
            if (is.null(self$layers)) {
                self$layers <- list()
            }

            self$validate()
        },

        print = function() {
            cat(paste0(
                "AnnDatR object with n_obs x n_vars = ",
                self$n_obs,
                " x ",
                self$n_vars,
                "\n"
            ))
            cat(paste0(
                "    obs: ",
                paste(colnames(self$obs), collapse = ", "),
                "\n"
            ))
            cat(paste0(
                "    var: ",
                paste(colnames(self$var), collapse = ", "),
                "\n"
            ))

            if (length(self$uns) > 0) {
                cat(paste0(
                    "    uns: ",
                    paste(names(self$uns), collapse = ", "),
                    "\n"
                ))
            }
            if (length(self$obsm) > 0) {
                cat(paste0(
                    "    obsm: ",
                    paste(names(self$obsm), collapse = ", "),
                    "\n"
                ))
            }
            if (length(self$layers) > 0) {
                cat(paste0(
                    "    layers: ",
                    paste(names(self$layers), collapse = ", "),
                    "\n"
                ))
            }
        },

        validate = function() {
            # Basic validation can be added here
            invisible(self)
        },

        filter_obs = function(column, string, negative = FALSE) {
            # 1. Identify cells to keep
            if (negative == FALSE) {
                cells_to_keep <- self$obs[[self$obs_names_col]][
                    self$obs[[column]] %in% string
                ]
            } else {
                cells_to_keep <- self$obs[[self$obs_names_col]][
                    !self$obs[[column]] %in% string
                ]
            }

            # 2. Filter obs
            self$obs <- self$obs %>%
                dplyr::filter(!!sym(self$obs_names_col) %in% cells_to_keep)

            # 3. Filter X (Columns)
            # Assuming X has a gene column + cell columns
            self$X <- self$X %>%
                dplyr::select(all_of(c(self$var_names_col, cells_to_keep)))

            # 4. Filter layers (Columns)
            for (layer_name in names(self$layers)) {
                self$layers[[layer_name]] <- self$layers[[layer_name]] %>%
                    dplyr::select(all_of(c(self$var_names_col, cells_to_keep)))
            }

            # 5. Filter obsm (Rows)
            # obsm matrices are usually Cells x Dims. We need to subset the rows matching cells_to_keep.
            # We assume rownames of obsm matrices match obs_names.
            for (obsm_name in names(self$obsm)) {
                obj <- self$obsm[[obsm_name]]
                if (is.matrix(obj) || is.data.frame(obj)) {
                    # Check if rownames are present and match
                    if (!is.null(rownames(obj))) {
                        self$obsm[[obsm_name]] <- obj[
                            rownames(obj) %in% cells_to_keep,
                            ,
                            drop = FALSE
                        ]
                    } else {
                        warning(paste0(
                            "obsm entry '",
                            obsm_name,
                            "' has no rownames. Cannot filter safely."
                        ))
                    }
                }
            }
        },

        filter_var = function(column, string, negative = FALSE) {
            if (negative == FALSE) {
                genes_to_keep <- self$var[[self$var_names_col]][
                    self$var[[column]] %in% string
                ]
            } else {
                genes_to_keep <- self$var[[self$var_names_col]][
                    !self$var[[column]] %in% string
                ]
            }

            # Filter var
            self$var <- self$var %>%
                dplyr::filter(!!sym(self$var_names_col) %in% genes_to_keep)

            # Filter X (Rows)
            self$X <- self$X %>%
                dplyr::filter(!!sym(self$var_names_col) %in% genes_to_keep)

            # Filter layers
            for (layer_name in names(self$layers)) {
                self$layers[[layer_name]] <- self$layers[[layer_name]] %>%
                    dplyr::filter(!!sym(self$var_names_col) %in% genes_to_keep)
            }
        }
    ),

    active = list(
        n_obs = function() {
            nrow(self$obs)
        },
        n_vars = function() {
            nrow(self$var)
        },
        shape = function() {
            c(self$n_obs, self$n_vars)
        },
        obs_names = function(value) {
            if (missing(value)) {
                return(self$obs[[self$obs_names_col]])
            } else {
                stop(
                    "Setting obs_names directly is not yet supported. Modify 'obs' instead."
                )
            }
        },
        var_names = function(value) {
            if (missing(value)) {
                return(self$var[[self$var_names_col]])
            } else {
                stop(
                    "Setting var_names directly is not yet supported. Modify 'var' instead."
                )
            }
        }
    )
)
