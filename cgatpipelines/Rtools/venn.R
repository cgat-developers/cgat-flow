#' Wrapper for VennDiagram with basic customization and command line functionality
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
suppressMessages(library(VennDiagram))


source(file.path(Sys.getenv("R_ROOT"), "experiment.R"))


run <- function(opt) {
  
  dataframes <- lapply(opt$filenames, read.table)
  names(dataframes) <- opt$names
  x <- lapply(dataframes, function(x) subset(x, padj < opt$pvalue))
  if(opt$sign == "positive"){
    x <- lapply(x, function(y) subset(y, log2FoldChange >0))}
  if(opt$sign == "negative"){
    x <- lapply(x, function(y) subset(y, log2FoldChange <0))}
  x <- lapply(x, rownames)
  flog.info(str(x))
  lapply(calculate.overlap(x), write, "venn.txt", append=TRUE, ncolumns=1000)


  venn.diagram(x,
               filename = 'venn_diagramm.png',
               output = TRUE ,
               imagetype="png" ,
               height = 480 ,
               width = 480 ,
               resolution = 300,
               compression = "lzw",
               lwd = 2,
               lty = 'blank',
               fill = opt$colours,
               cex = 0.8,
               fontface = "bold",
               fontfamily = "sans",
               cat.cex = 0.4,
               cat.fontface = "bold",
               cat.default.pos = "outer",
               cat.fontfamily = "sans"
  )
}

main <- function() {
  
  option_list <- list(
    make_option(
      "--filenames",
      dest = "filenames",
      type = "character",
      default = "file1.tsv,file2.tsv",
      help = paste("filenames with results data")
    ),
    make_option(
      "--colours",
      dest = "colours",
      type = "character",
      default = "blue,red",
      help = paste("Venn Diagram Colours")
    ),
    make_option(
      "--names",
      dest = "names",
      type = "character",
      default = "A,B",
      help = paste("Names")
    ),
    make_option(
      "--pvalue",
      dest = "pvalue",
      type = "numeric",
      default = 0.05,
      help = paste("p value threshold")
    ),
    make_option(
      "--sign",
      dest = "sign",
      type = "character",
      default = "",
      help = paste("p value threshold")
    )
  )
  opt <- experiment_start(option_list = option_list,
                          description = description)
  
  if (!is.null(opt$filenames)) {
    opt$filenames = unlist(strsplit(opt$filenames, ","))
  }
  if (!is.null(opt$colours)) {
    opt$colours = unlist(strsplit(opt$colours, ","))
  }
  if (!is.null(opt$names)) {
    opt$names = unlist(strsplit(opt$names, ","))
  }
  run(opt)
  
  experiment_stop()
}

main()
