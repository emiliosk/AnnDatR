get_possible_ensembl_versions <- function() {
    return(biomaRt::listEnsemblArchives()[['version']])
}

get_possible_datasets <- function(ensembl_version) {
    ensembl <- biomaRt::useEnsembl(biomart = "genes", version = ensembl_version)
    return(biomaRt::listDatasets(ensembl)[["dataset"]])
}

get_mito_genes <- function(dataset, ensembl_version) {
    ensembl <- biomaRt::useEnsembl(
        biomart = 'genes',
        dataset = dataset,
        version = ensembl_version
    )

    mito_genes <- biomaRt::getBM(
        attributes = c('ensembl_gene_id', 'external_gene_name'),
        filters = 'chromosome_name',
        values = 'MT',
        mart = ensembl
    )
    return(mito_genes)
}


get_basic_gene_info <- function(dataset, ensembl_version) {
    ensembl <- biomaRt::useEnsembl(
        biomart = 'genes',
        dataset = dataset,
        version = ensembl_version
    )

    gene_data <- biomaRt::getBM(
        attributes = c(
            'ensembl_gene_id',
            'external_gene_name',
            'chromosome_name',
            'gene_biotype',
            'start_position',
            'end_position'
        ),
        mart = ensembl
    ) %>%
        dplyr::mutate('length' = end_position - start_position)
    return(gene_data)
}
