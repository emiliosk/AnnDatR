set_PCA <- function(
    AnnDatR,
    nPcs = 40,
    transform = 'log1p',
    scale_by = 'sample',
    id = 'pca',
    layer = NULL,
    pca_method = 'svd'
) {
    if (!is.null(layer)) {
        wide_data <- AnnDatR$layers[[layer]]
    } else {
        wide_data <- AnnDatR$X
    }

    nPcs = min(nPcs, dim(wide_data)[1])

    # Apply the desired transformation
    if (transform == 'log1p') {
        transformed_data <- wide_data %>%
            dplyr::mutate_if(is.numeric, function(x) {
                log1p(x)
            }) %>%
            tibble::column_to_rownames(colnames(wide_data)[1])
    } else if (transform == 'sqrt') {
        transformed_data <- wide_data %>%
            dplyr::mutate_if(is.numeric, function(x) {
                sqrt(x)
            }) %>%
            tibble::column_to_rownames(colnames(wide_data)[1])
    } else if (transform == 'none') {
        transformed_data <- wide_data %>%
            tibble::column_to_rownames(colnames(wide_data)[1])
    } else {
        stop('Transformation not defined, only "log1p" or "sqrt" are allowed')
    }

    # Scale the data using the scaling function
    scaled_data <- scale_data(transformed_data, method = scale_by)

    # Perform PCA
    pca_results <- pcaMethods::pca(
        scaled_data,
        nPcs = nPcs,
        method = pca_method
    )

    # Store the PCA results in AnnDatR
    AnnDatR$uns[[id]] <- pca_results
    AnnDatR$obsm[[paste0('X_', id)]] <- pcaMethods::scores(AnnDatR$uns[[id]])

    summary(AnnDatR$uns[[id]])
}


set_PCA_legacy <- function(
    AnnDatR,
    nPcs = 40,
    transform = 'log1p',
    scale_by = 'sample',
    layer = NULL
) {
    if (!is.null(layer)) {
        wide_data <- AnnDatR$layers$layer
    } else {
        wide_data <- AnnDatR$X
    }

    nPcs = min(nPcs, dim(wide_data)[1])

    if (transform == 'log1p') {
        pca_results <- wide_data %>%
            dplyr::mutate_if(is.numeric, function(x) {
                log1p(x)
            }) %>%
            tibble::column_to_rownames(colnames(wide_data)[1])
    } else if (transform == 'sqrt') {
        pca_results <- wide_data %>%
            dplyr::mutate_if(is.numeric, function(x) {
                sqrt(x)
            }) %>%
            tibble::column_to_rownames(colnames(wide_data)[1])
    } else if (transform == 'none') {
        pca_results <- wide_data %>%
            tibble::column_to_rownames(colnames(wide_data)[1])
    } else {
        stop('not defined transformation, only "log1p" or "sqrt"; or "none"')
    }

    if (scale_by == 'gene') {
        pca_results <- pca_results %>%
            t() %>%
            pcaMethods::pca(nPcs = nPcs, scale = 'uv')
    } else if (scale_by == 'sample') {
        pca_results <- pca_results %>%
            pcaMethods::prep(scale = 'uv', eps = .Machine$double.eps) %>%
            t() %>%
            pcaMethods::pca(nPcs = nPcs)
    } else if (scale_by == 'pareto') {
        pca_results <- pca_results %>%
            t() %>%
            pcaMethods::pca(nPcs = nPcs, scale = 'pareto')
    } else if (scale_by == 'vector') {
        pca_results <- pca_results %>%
            t() %>%
            pcaMethods::pca(nPcs = nPcs, scale = 'vector')
    } else {
        stop(
            'Scaling method not defined, only "sample", "gene", "pareto", "vector", are allowed'
        )
    }

    AnnDatR$uns$pca <- pca_results

    AnnDatR$obsm[['X_pca']] <- pcaMethods::scores(AnnDatR$uns$pca)
    summary(AnnDatR$uns$pca)
}


kaisers_PCA_rule <- function(pca_results, with_alternative = TRUE) {
    if (inherits(pca_results, 'prcomp')) {
        return(kaisers_PCA_rule_prcomp(pca_results, with_alternative))
    }
    # Extract the squared standard deviations (eigenvalues)
    squared_devs <- pcaMethods::sDev(pca_results)^2

    # Find the first component where the squared standard deviation (eigenvalue) is less than 1
    n_comp <- which(squared_devs < 1)[1]

    # Check if n_comp is NA or exceeds the number of components
    if (is.na(n_comp) || n_comp > length(squared_devs)) {
        stop("No eigenvalue is lower than 1")
    }

    # Check cumulative R2 value at the identified component
    if (pca_results@R2cum[n_comp] < 0.8) {
        if (with_alternative) {
            print(
                'Explained variance at Kaiser rule is under 80%, suggesting at least 80% Variation'
            )
            n_comp <- which(pca_results@R2cum > 0.8)[1]
            # If no component satisfies the 0.8 threshold
            if (is.na(n_comp)) {
                stop(
                    "No principal component achieves the cumulative R2 threshold of 0.8"
                )
            }
        } else {
            print('Suggested n_comp explains less than 80% variation')
        }
    } else {
        print("Kaiser's rule is also above 80% variation")
    }

    return(n_comp)
}


kaisers_PCA_rule_prcomp <- function(pca_results, with_alternative = TRUE) {
    # Extract the squared standard deviations (eigenvalues)
    squared_devs <- (pca_results$sdev)^2

    # Find the first component where the squared standard deviation (eigenvalue) is less than 1
    n_comp <- which(squared_devs < 1)[1]

    # Check if n_comp is NA or exceeds the number of components
    if (is.na(n_comp) || n_comp > length(squared_devs)) {
        stop("No eigenvalue is lower than 1")
    }

    # Calculate cumulative variance explained (R2 equivalent)
    cumulative_variance_explained <- cumsum(squared_devs) / sum(squared_devs)

    # Check if the cumulative variance explained is below 80% at the identified component
    if (cumulative_variance_explained[n_comp] < 0.8) {
        if (with_alternative) {
            print(
                'Explained variance at Kaiser rule is under 80%, suggesting at least 80% variation'
            )
            n_comp <- which(cumulative_variance_explained > 0.8)[1]
            # If no component satisfies the 0.8 threshold
            if (is.na(n_comp)) {
                stop(
                    "No principal component achieves the cumulative variance threshold of 0.8"
                )
            }
        } else {
            print('Suggested n_comp explains less than 80% variation')
        }
    } else {
        print("Kaiser's rule is also above 80% variation")
    }

    return(n_comp)
}

set_UMAP <- function(
    AnnDatR,
    n_neighbors = 15,
    n_epochs = 1000,
    seed = 37,
    min_dist = 0.5,
    spread = 1,
    metric = "euclidean",
    pc_lim = NULL,
    pca_id = 'pca'
) {
    if (is.null(AnnDatR$uns[[pca_id]])) {
        stop(paste0(
            'AnnDatR$uns$',
            pca_id,
            ' not found. Call set_PCA function before plotting.'
        ))
    }
    pca_results <- AnnDatR$uns[[pca_id]]
    if (is.null(pc_lim)) {
        pc_lim <- which(pca_results@R2cum > 0.8)[1]
    }

    #pc_lim_sd <-
    #  rev(which(pca_results@sDev > 1))[1]

    umap_results <- pca_results@scores[, 1:pc_lim] %>%
        uwot::umap(
            n_neighbors = n_neighbors,
            n_epochs = n_epochs,
            seed = seed,
            min_dist = min_dist,
            spread = spread,
            metric = metric
        )
    AnnDatR$obsm[['X_umap']] <- umap_results
}


set_dendrogram <- function(
    AnnDatR,
    cor_method = 'spearman',
    cor_use = 'everything',
    hclust_method = 'complete',
    layer = NULL,
    prefix = NULL
) {
    if (!is.null(layer)) {
        X <- AnnDatR$layers[[layer]]
    } else {
        X <- AnnDatR$X
    }
    AnnDatR$uns[[paste0(prefix, 'hclust')]] <- X %>%
        tibble::column_to_rownames(AnnDatR$var_names_col) %>%
        stats::cor(method = cor_method, use = cor_use) %>%
        {
            1 - .
        } %>%
        stats::as.dist() %>%
        stats::hclust(method = hclust_method)
    AnnDatR$uns[[paste0(prefix, 'dendrogram')]] <- AnnDatR$uns[[paste0(
        prefix,
        'hclust'
    )]] %>%
        ggdendro::dendro_data()
}

get_distance <- function(
    AnnDatR,
    method,
    n_comp = NULL,
    id = NULL,
    pca_id = 'pca'
) {
    if (is.null(AnnDatR$uns[[pca_id]])) {
        stop(
            'AnnDatR$uns$pca not found. Call set_PCA function before plotting.'
        )
    }
    pca_results <- AnnDatR$uns[[pca_id]]
    if (is.null(n_comp)) {
        n_comp <- which(pca_results@R2cum > 0.8)[1]
    }

    if (is.null(AnnDatR$uns$distance)) {
        AnnDatR$uns$distance = list()
    }

    distance <- pca_results@scores[, 1:n_comp] %>%
        factoextra::get_dist(method = method)

    key <- if (!is.null(id)) id else method
    AnnDatR$uns$distance[[key]] <- distance
}
