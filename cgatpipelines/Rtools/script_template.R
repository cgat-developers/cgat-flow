#' filtering single cell data based on QC metrics
#'
#' WARNING: This script is work-in-progress
#' 
#' Example usage:
#' 
#' cgat sc-counts2counts --counts-filename=featurecounts.tsv --phenotypes-filename=phenodata.tsv --factor=group,mouse_id,collection_date,slice_depth,slice_number,pipette_visual,timepoint > filtered_counts.tsv
#'
#' `feature_counts.tsv` is a table (tab-separated) of ngenes x ncells,
#' that is the genes are in rows and the columns are cells.
#'
#' `phenodata.tsv` is a table (tab-separated) of ncells x nfeatures,
#' that is rows are cells and features are in columns. The table should contain
#' a column called `sample_id` that will match the columns in the table
#' `feature_counts.tsv`.
#'
#' Features can then be selected in the `--factor` option to be
#' plotted.
#'
#' -> todo: parameterize detection of ERCC (pattern?)
#' -> todo: parameterize definition of mitochondrial genes - currently hardcoded for mouse.

## conda dependencies: bioconductor-scater r-cairo

suppressMessages(library(futile.logger))
suppressMessages(library(getopt))

source(file.path(Sys.getenv("R_ROOT"), "experiment.R"))


run <- function(opt) {

    ## do something
    
}

main <- function() {

    option_list <- list(
        make_option(
            "--counts-filename",
            dest = "counts_filename",
            type = "character",
            default = "featurecounts.tsv",
            help = paste("filename with input data of counts")
        ),
        make_option(
            "--factor",
            dest = "factors",
            type = "character",
            default = "group,collection_date",
            help = paste("factors to colour QC plots by.")
        )
    )
    opt <- experiment_start(option_list = option_list,
                            description = description)

    if (!is.null(opt$factors)) {
        opt$factors = unlist(strsplit(opt$factors, ","))
    }
    run(opt)
    
    experiment_stop()
}

main()
