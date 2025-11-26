plot_PCA <- function(
    AnnDatR,
    color_by,
    color_code = NULL,
    pca_id = 'pca',
    alpha = 0.7,
    leg_ncol = 3
) {
    if (is.null(AnnDatR$uns[[pca_id]])) {
        stop(
            'AnnDatR$uns$pca not found. Call set_PCA function before plotting.'
        )
    }
    pca_results <- AnnDatR$uns[[pca_id]]
    obs <- AnnDatR$obs %>% dplyr::rename(cell_id = AnnDatR$obs_names_col)
    plot_data <- pca_results %>%
        pcaMethods::scores() %>%
        tibble::as_tibble(rownames = 'cell_id') %>%
        dplyr::left_join(obs, by = dplyr::join_by(cell_id))

    if (!is.null(color_code)) {
        plot_data <-
            plot_data %>%
            dplyr::rename(
                color_by = !!rlang::sym(color_by),
                color_code = !!rlang::sym(color_code)
            ) %>%
            dplyr::mutate(color_by = stringr::str_to_sentence(color_by))
        pal <- plot_data %>%
            dplyr::select(color_by, color_code) %>%
            dplyr::mutate(color_by = stringr::str_to_sentence(color_by)) %>%
            dplyr::arrange(color_by) %>%
            dplyr::distinct() %>%
            {
                stats::setNames(.$color_code, .$color_by)
            }

        plot <- plot_data %>%
            ggplot2::ggplot(ggplot2::aes(PC1, PC2)) +
            ggplot2::geom_point(
                ggplot2::aes(fill = color_by),
                color = "gray25",
                shape = 21,
                alpha = alpha,
                size = 2,
                stroke = 0.5
            ) +
            ggplot2::scale_fill_manual(values = pal) +
            ggplot2::xlab(paste(
                "PC1",
                pca_results@R2[1] * 100,
                "% of the variance"
            )) +
            ggplot2::ylab(paste(
                "PC2",
                pca_results@R2[2] * 100,
                "% of the variance"
            )) +
            ggplot2::theme_classic() +
            ggplot2::theme(
                legend.text = ggplot2::element_text(size = 10),
                legend.title = ggplot2::element_blank(),
                legend.position = "bottom",
                legend.spacing.y = grid::unit(-0.3, 'cm'),
                legend.spacing.x = grid::unit(-0.01, 'cm'),
                legend.key.height = grid::unit(0.01, "cm"),
                legend.key.spacing.y = grid::unit(0.05, 'cm'),
                # Adjust the height as needed
                legend.key.width = grid::unit(0.5, "cm")
            ) + # Adjust the width as needed
            ggplot2::guides(
                fill = ggplot2::guide_legend(ncol = leg_ncol, byrow = TRUE)
            ) +
            ggplot2::coord_fixed()
        # +
        # guides(shape = guide_legend(ncol = 1, byrow = TRUE))
    } else {
        plot_data <-
            plot_data %>%
            dplyr::rename(color_by = !!rlang::sym(color_by)) %>%
            dplyr::mutate(color_by = stringr::str_to_sentence(color_by))

        plot <- plot_data %>%
            ggplot2::ggplot(ggplot2::aes(PC1, PC2)) +
            ggplot2::geom_point(
                ggplot2::aes(fill = color_by),
                color = "gray25",
                shape = 21,
                alpha = alpha,
                size = 2,
                stroke = 0.5
            ) +
            ggplot2::xlab(paste(
                "PC1",
                pca_results@R2[1] * 100,
                "% of the variance"
            )) +
            ggplot2::ylab(paste(
                "PC2",
                pca_results@R2[2] * 100,
                "% of the variance"
            )) +
            ggplot2::theme_classic() +
            ggplot2::theme(
                legend.text = ggplot2::element_text(size = 10),
                legend.title = ggplot2::element_blank(),
                legend.position = "bottom",
                legend.spacing.y = grid::unit(-0.3, 'cm'),
                legend.spacing.x = grid::unit(-0.01, 'cm'),
                legend.key.spacing.y = grid::unit(0.05, 'cm'),
                legend.key.height = grid::unit(0.01, "cm"),
                # Adjust the height as needed
                legend.key.width = grid::unit(0.5, "cm")
            ) + # Adjust the width as needed) +
            ggplot2::coord_fixed() +
            ggplot2::guides(
                fill = ggplot2::guide_legend(ncol = leg_ncol, byrow = TRUE)
            ) # +
        # guides(shape = guide_legend(ncol = 1, byrow = TRUE))
    }

    return(plot)
}

plot_PCA.highlight <-
    function(
        AnnDatR,
        obs_column,
        element_to_highlight,
        plot_names = TRUE,
        pca_id = 'pca',
        leg_ncol = 3
    ) {
        if (is.null(AnnDatR$uns[[pca_id]])) {
            stop(
                'AnnDatR$uns$pca not found. Call set_PCA function before plotting.'
            )
        }
        pca_results <- AnnDatR$uns[[pca_id]]
        obs <- AnnDatR$obs %>% dplyr::rename(cell_id = AnnDatR$obs_names_col)
        plot_data <- pca_results %>%
            pcaMethods::scores() %>%
            tibble::as_tibble(rownames = 'cell_id') %>%
            dplyr::left_join(obs, by = dplyr::join_by(cell_id)) %>%
            dplyr::mutate(
                highlight_color = dplyr::case_when(
                    .[[obs_column]] == element_to_highlight ~ "#ff0000",
                    TRUE ~ "#a9a9a9"
                )
            )
        pl <- ggplot2::ggplot() +
            ggplot2::geom_point(
                data = plot_data,
                ggplot2::aes(PC1, PC2, colour = highlight_color),
                alpha = 0.6
            ) +
            ggplot2::scale_color_identity() +
            ggplot2::xlab(paste(
                "PC1",
                pca_results@R2[1] * 100,
                "% of the variance"
            )) +
            ggplot2::ylab(paste(
                "PC2",
                pca_results@R2[2] * 100,
                "% of the variance"
            )) +
            ggplot2::ggtitle(paste0(element_to_highlight)) +
            ggplot2::theme_classic()
        if (plot_names) {
            pl <-
                pl +
                ggrepel::geom_text_repel(
                    data = plot_data %>%
                        dplyr::filter(
                            !!rlang::sym(obs_column) == element_to_highlight
                        ),
                    ggplot2::aes(PC1, PC2, label = cell_id)
                )
        }
        return(pl)
    }


plot_UMAP <- function(
    AnnDatR,
    color_by,
    color_code,
    alpha = 0.7,
    size = 2,
    leg_ncol = 3
) {
    if (is.null(AnnDatR$obsm$X_umap)) {
        stop(
            'AnnDatR$obsm$X_umap not found. Call set_UMAP function before plotting.'
        )
    }

    obs <- AnnDatR$obs %>% dplyr::rename('cell_id' = AnnDatR$obs_names_col)
    plot_data <- AnnDatR$obsm$X_umap %>%
        tibble::as_tibble() %>%
        stats::setNames(paste0("UMAP", 1:ncol(.))) %>%
        dplyr::mutate(cell_id = rownames(AnnDatR$obsm$X_umap)) %>%
        dplyr::left_join(obs, by = dplyr::join_by(cell_id))

    if (!is.null(color_code)) {
        plot_data <-
            plot_data %>%
            dplyr::rename(
                color_by = !!rlang::sym(color_by),
                color_code = !!rlang::sym(color_code)
            ) %>%
            dplyr::mutate(color_by = stringr::str_to_sentence(color_by))
        pal <- plot_data %>%
            dplyr::select(color_by, color_code) %>%
            dplyr::mutate(color_by = stringr::str_to_sentence(color_by)) %>%
            dplyr::arrange(color_by) %>%
            dplyr::distinct() %>%
            {
                stats::setNames(.$color_code, .$color_by)
            }

        plot <- plot_data %>%
            ggplot2::ggplot(ggplot2::aes(UMAP1, UMAP2)) +
            ggplot2::geom_point(
                ggplot2::aes(fill = color_by),
                color = "gray25",
                shape = 21,
                alpha = alpha,
                size = size,
                stroke = 0.5
            ) +
            ggplot2::scale_fill_manual(values = pal) +
            ggplot2::theme_classic() +
            ggplot2::theme(
                legend.text = ggplot2::element_text(size = 10),
                legend.title = ggplot2::element_blank(),
                legend.position = "bottom",
                legend.spacing.y = grid::unit(-0.3, 'cm'),
                legend.spacing.x = grid::unit(-0.01, 'cm'),
                legend.key.height = grid::unit(0.01, "cm"),
                legend.key.spacing.y = grid::unit(0.05, 'cm'),
                # Adjust the height as needed
                legend.key.width = grid::unit(0.5, "cm")
            ) + # Adjust the width as needed
            ggplot2::guides(
                fill = ggplot2::guide_legend(ncol = leg_ncol, byrow = TRUE)
            ) +
            ggplot2::coord_fixed()
        # +
        # guides(shape = guide_legend(ncol = 1, byrow = TRUE))
    } else {
        plot_data <-
            plot_data %>%
            dplyr::rename(color_by = !!rlang::sym(color_by)) %>%
            dplyr::mutate(color_by = stringr::str_to_sentence(color_by))

        plot <- plot_data %>%
            ggplot2::ggplot(ggplot2::aes(UMAP1, UMAP2)) +
            ggplot2::geom_point(
                ggplot2::aes(fill = color_by),
                color = "gray25",
                shape = 21,
                alpha = 0.7,
                size = 2,
                stroke = 0.5
            ) +
            ggplot2::theme_classic() +
            ggplot2::theme(
                legend.text = ggplot2::element_text(size = 10),
                legend.title = ggplot2::element_blank(),
                legend.position = "bottom",
                legend.spacing.y = grid::unit(-0.3, 'cm'),
                legend.spacing.x = grid::unit(-0.01, 'cm'),
                legend.key.spacing.y = grid::unit(0.05, 'cm'),
                legend.key.height = grid::unit(0.01, "cm"),
                # Adjust the height as needed
                legend.key.width = grid::unit(0.5, "cm")
            ) + # Adjust the width as needed) +
            ggplot2::coord_fixed() #+
        ggplot2::guides(
            fill = ggplot2::guide_legend(ncol = leg_ncol, byrow = TRUE)
        ) # +
        # guides(shape = guide_legend(ncol = 1, byrow = TRUE))
    }

    return(plot)
}

plot_UMAP_plotly <- function(
    AnnDatR,
    color_by,
    color_code,
    hover_text_columns = NULL,
    alpha = 0.7,
    size = 2,
    leg_ncol = 3
) {
    if (is.null(AnnDatR$obsm$X_umap)) {
        stop(
            'AnnDatR$obsm$X_umap not found. Call set_UMAP function before plotting.'
        )
    }

    obs <- AnnDatR$obs %>% dplyr::rename('cell_id' = AnnDatR$obs_names_col)

    # Prepare UMAP data
    plot_data <- AnnDatR$obsm$X_umap %>%
        tibble::as_tibble() %>%
        stats::setNames(paste0("UMAP", 1:ncol(.))) %>%
        dplyr::mutate(cell_id = rownames(AnnDatR$obsm$X_umap)) %>%
        dplyr::left_join(obs, by = "cell_id")

    # Handle coloring by color_by column and color_code column
    if (!is.null(color_code)) {
        plot_data <- plot_data %>%
            dplyr::rename(
                color_by = !!rlang::sym(color_by),
                color_code = !!rlang::sym(color_code)
            ) %>%
            dplyr::mutate(color_by = stringr::str_to_sentence(color_by))

        pal <- plot_data %>%
            dplyr::select(color_by, color_code) %>%
            dplyr::distinct() %>%
            dplyr::arrange(color_by) %>%
            {
                stats::setNames(.$color_code, .$color_by)
            }
    } else {
        plot_data <- plot_data %>%
            dplyr::rename(color_by = !!rlang::sym(color_by)) %>%
            dplyr::mutate(color_by = stringr::str_to_sentence(color_by))
    }

    # Create hover text if hover_text_columns are provided
    if (!is.null(hover_text_columns)) {
        hover_text_columns <- c('cell_id', hover_text_columns, 'color_by')
        hover_text <- plot_data %>%
            dplyr::select(dplyr::all_of(hover_text_columns)) %>%
            apply(1, function(row) {
                paste(names(row), row, sep = ": ", collapse = "<br>")
            })
        plot_data <- plot_data %>% dplyr::mutate(hover_text = hover_text)
    } else {
        plot_data <- plot_data %>% dplyr::mutate(hover_text = NA)
    }

    # Build plotly plot
    plot <- plotly::plot_ly(
        data = plot_data,
        x = ~UMAP1,
        y = ~UMAP2,
        type = "scatter",
        mode = "markers",
        marker = list(color = ~color_code, size = size, opacity = alpha),
        hoverinfo = "text",
        text = ~hover_text
    ) %>%
        plotly::layout(
            legend = list(
                orientation = "h",
                x = 0.5,
                y = -0.1,
                xanchor = "center",
                ncol = leg_ncol
            ),
            xaxis = list(title = "UMAP1"),
            yaxis = list(title = "UMAP2"),
            showlegend = TRUE
        )

    return(plot)
}

plot_UMAP.highlight <- function(
    AnnDatR,
    obs_column,
    element_to_highlight,
    plot_names = TRUE,
    leg_ncol = 3
) {
    if (is.null(AnnDatR$obsm$X_umap)) {
        stop(
            'AnnDatR$obsm$X_umap not found. Call set_PCA function before plotting.'
        )
    }

    obs <- AnnDatR$obs %>% dplyr::rename('cell_id' = AnnDatR$obs_names_col)
    plot_data <- AnnDatR$obsm$X_umap %>%
        tibble::as_tibble() %>%
        stats::setNames(paste0("UMAP", 1:ncol(.))) %>%
        dplyr::mutate(cell_id = rownames(AnnDatR$obsm$X_umap)) %>%
        dplyr::left_join(obs, by = dplyr::join_by(cell_id)) %>%
        dplyr::mutate(
            highlight_color = dplyr::case_when(
                .[[obs_column]] == element_to_highlight ~ "#ff0000",
                TRUE ~ "#a9a9a9"
            )
        )
    pl <- ggplot2::ggplot() +
        ggplot2::geom_point(
            data = plot_data,
            ggplot2::aes(UMAP1, UMAP2, colour = highlight_color),
            alpha = 0.6
        ) +
        ggplot2::scale_color_identity() +
        ggplot2::ggtitle(paste0(element_to_highlight)) +
        ggplot2::theme_classic()
    if (plot_names) {
        pl <- pl +
            ggrepel::geom_text_repel(
                data = plot_data %>%
                    dplyr::filter(
                        !!rlang::sym(obs_column) == element_to_highlight
                    ),
                ggplot2::aes(UMAP1, UMAP2, label = cell_id)
            )
    }
    return(pl)
}


plot_dendrogram = function(
    AnnDatR,
    color_code,
    x_expansion = 2,
    y_expansion = 0.4,
    cor_method = 'spearman'
) {
    if (is.null(AnnDatR$uns$dendrogram)) {
        stop(
            'AnnDatR$uns$dendrogram not found. Call set_dendrogram function before plotting.'
        )
    }
    plot_dendro_sc <- AnnDatR$uns$dendrogram
    plot_dendro_sc$labels$label <- stringr::str_to_sentence(
        plot_dendro_sc$labels$label
    )
    plot_order_cell_type_name <- plot_dendro_sc$labels$label

    dendro_plot_data <-
        dplyr::left_join(
            plot_dendro_sc$segments,
            plot_dendro_sc$labels,
            by = c("x" = "x", "yend" = "y")
        )

    pal <- stats::setNames(
        AnnDatR$obs[[color_code]],
        stringr::str_to_sentence(AnnDatR$obs_names())
    )

    left_plot <-
        dendro_plot_data %>%
        ggplot2::ggplot() +
        ggplot2::geom_segment(ggplot2::aes(
            x = y,
            y = x,
            xend = yend,
            yend = xend
        )) +
        # geom_rect(aes(xmin=0, ymin=x + 0.5,
        #               xmax=-0.02, ymax=xend - 0.5,
        #               fill = tissue),
        #           show.legend = F) +
        ggplot2::geom_point(
            data = dendro_plot_data %>% tidyr::drop_na(),
            ggplot2::aes(
                x = 0,
                y = x,
                #  shape = species,
                fill = label,
                color = label
            ),
            #color = "gray25",
            alpha = 1,
            size = 3,
            stroke = 0.7
        ) +
        ggplot2::geom_text(
            data = dendro_plot_data %>% tidyr::drop_na(),
            ggplot2::aes(x = -0.02, y = x, label = label),
            hjust = 0,
            show.legend = F
        ) +
        # scale_shape_manual(values = shape_def ) +
        ggplot2::scale_color_manual(values = pal, guide = 'none') +
        ggplot2::scale_fill_manual(values = pal, guide = "none") +
        ggplot2::scale_x_reverse(
            expand = ggplot2::expansion(x_expansion),
            position = "top"
        ) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(y_expansion)) +
        ggplot2::xlab(paste0("1 - ", cor_method)) +

        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            legend.title = ggplot2::element_blank(),
            plot.margin = grid::unit(c(1, 1, 1, 1), units = "mm"),
            panel.background = ggplot2::element_blank()
        )

    return(left_plot)
}
