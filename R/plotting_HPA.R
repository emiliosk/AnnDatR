plot_specificity_barplot = function(
    gene_classification = gene_classification,
    title = "Specificity category per cell type",
    return_data = FALSE
) {
    ordered_names_sp <-
        gene_classification %>%
        dplyr::filter(
            spec_category %in%
                c("group enriched", "tissue enriched", "tissue enhanced")
        ) %>%
        tidyr::separate_rows(enriched_tissues, sep = ";") %>%
        dplyr::group_by(spec_category, enriched_tissues) %>%
        dplyr::count() %>%
        dplyr::ungroup() %>%
        dplyr::group_by(enriched_tissues) %>%
        dplyr::summarise(sum = sum(n)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(sum) %>%
        .$enriched_tissues %>%
        stringr::str_to_sentence()

    specificity_palette <- rev(gene_category_pal)

    plot_data <- gene_classification %>%
        dplyr::filter(
            spec_category %in%
                c("group enriched", "tissue enriched", "tissue enhanced")
        ) %>%
        tidyr::separate_rows(enriched_tissues, sep = ";") %>%
        dplyr::group_by(spec_category, enriched_tissues) %>%
        dplyr::count() %>%
        #mutate(enriched_tissues = factor(str_to_sentence(enriched_tissues), ordered_names_sp)) %>%
        dplyr::mutate(
            spec_category = factor(
                stringr::str_to_sentence(spec_category),
                c("Tissue enriched", "Group enriched", "Tissue enhanced")
            )
        )

    if (return_data) {
        return(plot_data)
    } else {
        plot <- plot_data %>%
            dplyr::mutate(
                enriched_tissues = factor(
                    stringr::str_to_sentence(enriched_tissues),
                    ordered_names_sp
                )
            ) %>%
            ggplot2::ggplot(ggplot2::aes(
                x = n,
                y = enriched_tissues,
                fill = spec_category
            )) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_manual(values = specificity_palette) +
            ggplot2::xlab("Number of genes") +
            ggplot2::theme_classic() +
            ggplot2::theme(
                axis.title.y = ggplot2::element_blank(),
                #  axis.title.x = element_blank(),
                axis.line.y = ggplot2::element_blank(),
                legend.position = "bottom",
                legend.title = ggplot2::element_blank()
            ) +
            ggplot2::scale_x_continuous(
                expand = ggplot2::expansion(mult = 0, add = 0)
            ) +
            ggplot2::ggtitle(title)
        return(plot)
    }
}

plot_distribution_barplot = function(
    gene_classification = gene_classification,
    title = "Distribution category per cell type",
    return_data = FALSE
) {
    ordered_names_sp <-
        gene_classification %>%
        dplyr::filter(
            dist_category %in%
                c(
                    "detected in all",
                    "detected in many",
                    "detected in some",
                    "detected in single"
                )
        ) %>%
        tidyr::separate_rows(tissues_detected, sep = ";") %>%
        dplyr::group_by(dist_category, tissues_detected) %>%
        dplyr::count() %>%
        dplyr::ungroup() %>%
        dplyr::group_by(tissues_detected) %>%
        dplyr::summarise(sum = sum(n)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(sum) %>%
        .$tissues_detected %>%
        stringr::str_to_sentence()

    specificity_palette <- rev(gene_category_pal)

    plot_data <- gene_classification %>%
        dplyr::filter(
            dist_category %in%
                c(
                    "detected in all",
                    "detected in many",
                    "detected in some",
                    "detected in single"
                )
        ) %>%
        tidyr::separate_rows(tissues_detected, sep = ";") %>%
        dplyr::group_by(dist_category, tissues_detected) %>%
        dplyr::count() %>%
        dplyr::mutate(
            tissues_detected = factor(
                stringr::str_to_sentence(tissues_detected),
                ordered_names_sp
            )
        ) %>%
        dplyr::mutate(
            dist_category = factor(
                stringr::str_to_sentence(dist_category),
                rev(c(
                    "Detected in all",
                    "Detected in many",
                    "Detected in some",
                    "Detected in single"
                ))
            )
        )
    if (return_data) {
        return(plot_data)
    } else {
        plot <- plot_data %>%
            ggplot2::ggplot(ggplot2::aes(
                x = n,
                y = tissues_detected,
                fill = dist_category
            )) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_manual(values = specificity_palette) +
            ggplot2::xlab("Number of genes") +
            ggplot2::theme_classic() +
            ggplot2::theme(
                axis.title.y = ggplot2::element_blank(),
                #  axis.title.x = element_blank(),
                axis.line.y = ggplot2::element_blank(),
                legend.position = "bottom",
                legend.title = ggplot2::element_blank()
            ) +
            ggplot2::scale_x_continuous(
                expand = ggplot2::expansion(mult = 0, add = 0)
            ) +
            ggplot2::ggtitle(title)

        return(plot)
    }
}


plot_specificity_tau = function(
    gene_classification = gene_classification,
    title = NULL
) {
    specificity_palette <- rev(gene_category_pal)
    p1 <- gene_classification %>%
        dplyr::mutate(
            spec_category = factor(
                stringr::str_to_sentence(spec_category),
                levels = names(specificity_palette)
            ),
            enriched_tissues = stringr::str_to_sentence(enriched_tissues)
        ) %>%
        dplyr::filter(spec_category != "Not detected") %>%
        ggplot2::ggplot(ggplot2::aes(
            x = tau_score,
            y = spec_category,
            fill = spec_category
        )) +
        ggplot2::geom_violin() +
        ggplot2::scale_fill_manual(
            values = gene_category_pal,
            name = "Specificity"
        ) +
        ggplot2::xlab("Tau score") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        )

    #ggsave("./final_plots/classification/tau_to_specificity.pdf",width = 5.5, height = 4)

    p2 <- gene_classification %>%
        ggplot2::ggplot(ggplot2::aes(tau_score)) +
        ggplot2::geom_histogram(bins = 100) +
        ggplot2::theme_classic() +
        ggplot2::ylab("Count") +
        ggplot2::theme(
            panel.background = ggplot2::element_rect("gray90"),
            #axis.title.y = element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank()
        ) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = 0, add = 0)
        )

    plot <- patchwork::wrap_plots(
        p2 + ggplot2::ggtitle(title),
        p1,
        ncol = 1,
        heights = c(1, 3)
    )
    return(plot)
    #ggsave("./main_plots/classification/cell_type_tau_specificity.pdf",height = 4, width = 6)
}

plot_alluvial_categories = function(
    gene_classification
) {
    width_alluvial = 0.1

    alluv_1 <-
        gene_classification %>%
        dplyr::mutate(
            Specificity = stringr::str_to_sentence(spec_category),
            Distribution = stringr::str_to_sentence(dist_category)
        ) %>%
        dplyr::select(Specificity, Distribution) %>%
        dplyr::mutate(row_n = dplyr::row_number()) %>%
        tidyr::gather(bar, chunk, -row_n) %>%
        dplyr::mutate(color_vars = 1) %>%
        dplyr::group_by(row_n) %>%
        dplyr::mutate(
            chunk_color = chunk[match(
                c("Specificity", "Distribution")[color_vars],
                bar
            )]
        ) %>%
        dplyr::ungroup() %>%

        dplyr::mutate(
            chunk = factor(
                chunk,
                levels = c(
                    'Tissue enriched',
                    'Group enriched',
                    'Tissue enhanced',
                    'Low tissue specificity',
                    'Detected in single',
                    'Detected in some',
                    'Detected in many',
                    'Detected in all',
                    'Not detected'
                )
            ),
            bar = factor(bar, levels = c("Specificity", "Distribution"))
        ) %>%

        ggplot2::ggplot(ggplot2::aes(
            x = bar,
            stratum = chunk,
            alluvium = row_n,
            y = 1
        )) +

        ggalluvial::geom_flow(
            ggplot2::aes(fill = chunk_color),
            show.legend = F,
            width = width_alluvial,
            knot.pos = 1 / 6
        ) +
        ggalluvial::geom_stratum(
            ggplot2::aes(fill = chunk),
            show.legend = F,
            color = NA,
            width = width_alluvial
        ) +

        ggplot2::scale_x_discrete(expand = c(.1, .1), position = "top") +
        ggplot2::scale_fill_manual(values = c(gene_category_pal)) +

        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 18, face = "bold"),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank()
        ) +
        ggplot2::coord_flip()

    # alluv_1

    flow_data <-
        ggplot2::ggplot_build(alluv_1)$data[[1]] %>%
        tibble::as_tibble() %>%
        {
            if ("side" %in% names(.)) {
                .
            } else {
                dplyr::mutate(
                    .,
                    side = dplyr::case_when(
                        flow == "from" ~ "start",
                        flow == "to" ~ "end"
                    )
                )
            }
        }

    stratum_data <-
        ggplot2::ggplot_build(alluv_1)$data[[2]]

    flow_data_labels <-
        flow_data %>%
        tibble::as_tibble() %>%
        dplyr::select(x, stratum, group, side, ymin, ymax) %>%
        tidyr::pivot_wider(
            names_from = side,
            values_from = c(x, stratum, ymin, ymax)
        ) %>%
        dplyr::mutate_at(
            c(
                "x_end",
                "ymax_end",
                "ymin_end",
                "x_start",
                "ymax_start",
                "ymin_start"
            ),
            as.numeric
        ) %>%
        dplyr::group_by(stratum_start, stratum_end, x_start, x_end) %>%
        dplyr::summarise(
            y_end = (min(ymin_end) + max(ymax_end)) / 2,
            y_start = (min(ymin_start) + max(ymax_start)) / 2,
            size = max(ymax_start) - min(ymin_start)
        )

    alluv_2 <-
        alluv_1 +
        ggplot2::geom_text(
            data = flow_data_labels,
            ggplot2::aes(
                x = x_start + width_alluvial / 2,
                y = y_start,
                label = size
            ),
            inherit.aes = F,
            size = 3,
            angle = -90,
            hjust = 1,
            vjust = 0.5
        ) +
        ggplot2::geom_text(
            data = flow_data_labels,
            ggplot2::aes(
                x = x_end - width_alluvial / 2,
                y = y_end,
                label = size
            ),
            inherit.aes = F,
            size = 3,
            angle = -90,
            hjust = 0,
            vjust = 0.5
        ) +

        # Stratum label

        ggplot2::geom_text(
            data = stratum_data %>%
                dplyr::filter(x == 1),
            ggplot2::aes(x = x - width_alluvial / 2, y = y, label = stratum),
            size = 4,
            vjust = 1.5,
            inherit.aes = F
        ) +
        ggplot2::geom_text(
            data = stratum_data %>%
                dplyr::filter(x == 2),
            ggplot2::aes(x = x + width_alluvial / 2, y = y, label = stratum),
            size = 4,
            vjust = -0.5,
            inherit.aes = F
        ) +

        ggplot2::geom_text(
            data = stratum_data,
            ggplot2::aes(x = x, y = y, label = ymax - ymin),
            size = 4,
            fontface = "bold",
            color = "white",
            inherit.aes = F
        )

    return(alluv_2)
}

plot_alluvial <- function(alluv_1_data, width = 0.1) {
    alluv_1 <- alluv_1_data %>%
        ggplot2::ggplot(ggplot2::aes(
            x = bar,
            stratum = chunk,
            alluvium = row_n,
            y = 1
        )) +

        ggalluvial::geom_flow(
            ggplot2::aes(fill = chunk_color),
            show.legend = F,
            width = width,
            knot.pos = 1 / 6
        ) +
        ggalluvial::geom_stratum(
            ggplot2::aes(fill = chunk),
            show.legend = F,
            color = NA,
            width = width
        ) +

        ggplot2::scale_x_discrete(expand = c(.1, .1), position = "top") +
        ggplot2::scale_fill_manual(values = c(gene_category_pal)) +

        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 18, face = "bold"),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank()
        ) +
        ggplot2::coord_flip()

    # alluv_1

    flow_data <-
        ggplot2::ggplot_build(alluv_1)$data[[1]] %>%
        tibble::as_tibble() %>%
        {
            if ("side" %in% names(.)) {
                .
            } else {
                dplyr::mutate(
                    .,
                    side = dplyr::case_when(
                        flow == "from" ~ "start",
                        flow == "to" ~ "end"
                    )
                )
            }
        }

    stratum_data <-
        ggplot2::ggplot_build(alluv_1)$data[[2]]

    flow_data_labels <-
        flow_data %>%
        tibble::as_tibble() %>%
        dplyr::select(x, stratum, group, side, ymin, ymax) %>%
        tidyr::pivot_wider(
            names_from = side,
            values_from = c(x, stratum, ymin, ymax)
        ) %>%
        dplyr::mutate_at(
            c(
                "x_end",
                "ymax_end",
                "ymin_end",
                "x_start",
                "ymax_start",
                "ymin_start"
            ),
            as.numeric
        ) %>%
        dplyr::group_by(stratum_start, stratum_end, x_start, x_end) %>%
        dplyr::summarise(
            y_end = (min(ymin_end) + max(ymax_end)) / 2,
            y_start = (min(ymin_start) + max(ymax_start)) / 2,
            size = max(ymax_start) - min(ymin_start)
        )
    alluv_2 <-
        alluv_1 +
        ggplot2::geom_text(
            data = flow_data_labels,
            ggplot2::aes(x = x_start + width / 2, y = y_start, label = size),
            inherit.aes = F,
            size = 3,
            angle = -90,
            hjust = 1,
            vjust = 0.5
        ) +
        ggplot2::geom_text(
            data = flow_data_labels,
            ggplot2::aes(x = x_end - width / 2, y = y_end, label = size),
            inherit.aes = F,
            size = 3,
            angle = -90,
            hjust = 0,
            vjust = 0.5
        ) +

        # Stratum label

        ggplot2::geom_text(
            data = stratum_data %>%
                dplyr::filter(x == 1),
            ggplot2::aes(x = x - width / 2, y = y, label = stratum),
            size = 4,
            vjust = 1.5,
            inherit.aes = F
        ) +
        ggplot2::geom_text(
            data = stratum_data %>%
                dplyr::filter(x == 2),
            ggplot2::aes(x = x + width / 2, y = y, label = stratum),
            size = 4,
            vjust = -0.5,
            inherit.aes = F
        ) +

        ggplot2::geom_text(
            data = stratum_data,
            ggplot2::aes(x = x, y = y, label = ymax - ymin),
            size = 4,
            fontface = "bold",
            color = "white",
            inherit.aes = F
        )

    return(alluv_2)
}
