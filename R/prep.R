scale_data <- function(data, method = 'sample') {
    ## like this works with the set_PCA,
    if (method == 'gene') {
        scaled_data <- data %>%
            t() %>%
            pcaMethods::prep(scale = 'uv', eps = .Machine$double.eps)
    } else if (method == 'sample') {
        scaled_data <- data %>%
            pcaMethods::prep(scale = 'uv', eps = .Machine$double.eps) %>%
            t()
    } else if (method == 'pareto') {
        scaled_data <- data %>% t() %>% pcaMethods::prep(scale = 'pareto')
    } else if (method == 'vector') {
        scaled_data <- data %>% t() %>% pcaMethods::prep(scale = 'vector')
    } else if (method == 'min-max') {
        scaled_data <- apply(t(data), 1, function(x) {
            (x - min(x)) / (max(x) - min(x))
        }) %>%
            t()
    } else if (method == 'max_sequencial') {
        scaled_data <- apply(t(data), 1, function(x) x / max(x)) %>% t()
    } else if (method == 'max') {
        #vectorized
        #scaled_data <- apply(t(data), 1, function(x) x / max(x)) %>% t()
        data <- t(data)
        scaled_data <- data / matrixStats::rowMaxs(data)
    } else {
        stop(
            'Scaling method not defined, only "sample", "gene", "pareto", "vector", "min-max", and "max" are allowed'
        )
    }

    return(scaled_data)
}

calc_tmm_normfactors <-
    function(
        object,
        method = c("TMM", "quantile"),
        refColumn = NULL,
        logratioTrim = 0.3,
        sumTrim = 0.05,
        doWeighting = TRUE,
        Acutoff = -1e10,
        quantile = 0.75
    ) {
        method <- match.arg(method)
        if (is.matrix(object)) {
            if (is.null(refColumn)) {
                refColumn <- 1
            }
            data <- object
            libsize <- colSums(data)
        } else {
            stop("calcNormFactors() only operates on 'matrix' objects")
        }

        if (refColumn == "median") {
            ref <-
                apply(data, MARGIN = 1, median)
        } else {
            ref <- data[, refColumn]
        }

        f <- switch(
            method,
            TMM = apply(
                data,
                2,
                NOISeq:::.calcFactorWeighted,
                ref = ref,
                logratioTrim = logratioTrim,
                sumTrim = sumTrim,
                doWeighting = doWeighting,
                Acutoff = Acutoff
            ),
            quantile = NOISeq:::.calcFactorQuantile(data, libsize, q = quantile)
        )
        f <- f / exp(mean(log(f)))
        return(f)
    }

calc_tmm_normfactors_NA_robust <- function(
    object,
    method = c("TMM", "quantile"),
    refColumn = NULL,
    logratioTrim = 0.3,
    sumTrim = 0.05,
    doWeighting = TRUE,
    Acutoff = -1e10,
    quantile = 0.75
) {
    method <- match.arg(method)

    # Check input
    if (!is.matrix(object)) {
        stop("calcNormFactors() only operates on 'matrix' objects")
    }

    # Clean input
    data <- object
    data[!is.finite(data)] <- NA # remove Inf/-Inf

    if (is.null(refColumn)) {
        refColumn <- 1
    }

    libsize <- colSums(data, na.rm = TRUE)

    # Get reference sample
    if (refColumn == "median") {
        ref <- apply(data, MARGIN = 1, median, na.rm = TRUE)
    } else {
        ref <- data[, refColumn]
        if (all(is.na(ref))) {
            stop("Reference column contains only NA values.")
        }
    }

    # Compute normalization factors
    f <- switch(
        method,
        TMM = apply(
            data,
            2,
            function(x) {
                tryCatch(
                    NOISeq:::.calcFactorWeighted(
                        x,
                        ref = ref,
                        logratioTrim = logratioTrim,
                        sumTrim = sumTrim,
                        doWeighting = doWeighting,
                        Acutoff = Acutoff
                    ),
                    error = function(e) NA_real_
                )
            }
        ),
        quantile = NOISeq:::.calcFactorQuantile(data, libsize, q = quantile)
    )

    # Rescale and protect against NA/Inf
    if (all(!is.finite(f))) {
        warning("All normalization factors are NA or non-finite.")
        return(rep(NA_real_, ncol(data)))
    }

    f <- f / exp(mean(log(f[f > 0 & is.finite(f)]), na.rm = TRUE))
    f[!is.finite(f)] <- NA_real_

    return(f)
}
