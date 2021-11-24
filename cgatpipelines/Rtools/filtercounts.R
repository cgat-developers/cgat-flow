#' Reading in Data and Basic filtering analysis Script
#'
#' 
#' Example usage:
#' 
#' Rscript PATH/TO/filtercounts.R 
#' 
#' input: directory containing read count files or tsv file containing reads
#' additional input variables: method used to generate file, model
#' output: `experiment_out.rds` is an experiment object after filtering
#' 

suppressMessages(library(futile.logger))
suppressMessages(library(getopt))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(csaw))
suppressMessages(library(rhdf5))
suppressMessages(library(tximport))
suppressMessages(library(DEXSeq))

source(file.path(Sys.getenv("R_ROOT"), "experiment.R"))

run <- function(opt) {
  
  ### READING DATA ###
  # Read in sampleData Table
  futile.logger::flog.info(paste("reading sampleData table from", normalizePath(opt$sampleData)))
  sampleData <- read.table(opt$sampleData, header = TRUE)
  sampleData <-sampleData[sampleData$include ==1, ]
  futile.logger::flog.info(paste("read sampleData ", paste(dim(sampleData), collapse = ",")))
  rownames(sampleData) <- sampleData$track

  futile.logger::flog.info(paste("reading in data from ", opt$source))
  if(opt$source %in% c("salmon", "kallisto")){
    # Read in Transcript to Gene Map
    tx2gene <- read_tsv(opt$tx2gene)
    colnames(tx2gene) <- c("ensembl_transcript_id", "ensembl_gene_id")
    if(opt$tx2gene_regex != "None"){
      tx2gene <- filter(tx2gene, !grepl(opt$tx2gene_regex,ensembl_transcript_id))
    }
    # Read in Data
    futile.logger::flog.info(opt$counts_dir)
    futile.logger::flog.info(sampleData$track)
    files <- file.path(opt$counts_dir, sampleData$track, "quant.sf")
    names(files) <- sampleData$track
    txi <- tximport(files, type = opt$source, tx2gene = tx2gene)
    if(opt$method == "deseq2"){
      dataset <- DESeqDataSetFromTximport(txi, colData = sampleData, design = formula(opt$model))
    } else if(opt$method == "edger"){
      cts <- txi$counts
      normMat <- txi$length
      normMat <- normMat/exp(rowMeans(log(normMat)))
      normCts <- cts/normMat
      normMat <- sweep(normMat, 2, eff.lib, "*")
      normMat <- log(normMat)
      dataset <- DGEList(counts=cts, 
                   samples=sampleData)
      dataset <- scaleOffset(dataset,normMat)
    } else if(opt$method == "Sleuth"){
      stop("Sleuth method not yet implemented. Sorry.")
      
    } else{
      stop("Method not defined. Allowable methods are \"DESeq2\", \"EdgeR\" or \"Sleuth\"")
    }
  }
  
  if(opt$source == "dexseq"){
    # Read in Data
    files <- file.path(opt$counts_dir, paste0(sampleData$track, ".txt"))
    names(files) <- sampleData$track
    if(opt$method != "dexseq"){
      stop("DEXSeq input is handled by diffexonexpression. Please correct the method argument.")
    }
    dataset = DEXSeqDataSetFromHTSeq(
      files,
      sampleData=sampleData,
      design=formula(opt$model),
      flattenedfile=opt$flattenedFile)
  } else if(opt$source == "counts_table"){
    # Read in Data
    raw <- read.table(file = gzfile(opt$counts_tsv), header=TRUE, row.name=1)
    experiment_tsv <- raw[,sampleData$track,drop=FALSE]
    if(opt$method == "deseq2"){
      dataset = DESeqDataSetFromMatrix(experiment_tsv, sampleData, design = formula(opt$model))
    } else if(opt$method == "edger"){
      dataset <- DGEList(counts=experiment_tsv, samples=sampleData)
    } else if(opt$method == "sleuth"){
      stop("Sleuth method not yet implemented. Sorry.")
    } else{
      stop("Method not defined. Allowable methods are \"DESeq2\", \"EdgeR\" or \"Sleuth\"")
    }
  }
  
  ### FILTERING ###
  if(opt$filter == TRUE) {
    futile.logger::flog.info(paste("filtering data ", opt$source))
    if(opt$method == "edger"){
      futile.logger::flog.info(paste("Counts before filtering ", paste(dim(dataset$counts), collapse = ",")))
      keep <- filterByExpr(dataset)
      dataset <- dataset[keep, , keep.lib.sizes=FALSE]
      counts_table <- dataset$counts
      futile.logger::flog.info(paste("Counts after filtering ", paste(dim(dataset$counts), collapse = ",")))
    } else if(opt$method == "deseq2"){
      futile.logger::flog.info(paste("Counts before filtering ", paste(dim(counts(dataset)), collapse = ",")))
      keep <- rowSums(counts(dataset)) >= 10
      dataset <- dataset[keep,]
      counts_table <- counts(dataset)
      futile.logger::flog.info(paste("Counts after filtering ", paste(dim(counts(dataset)), collapse = ",")))
    } else if(opt$method == "dexseq"){
      futile.logger::flog.info(paste("Filtering for DEXSeq not implemented "))
      counts_table <- counts(dataset)
    }
  }
  else {
    futile.logger::flog.info(paste("No filtering on dataset performed.", opt$source))
  }

  ### SAVING DATA ###
  file = get_output_filename(paste0(opt$outdir,"/experiment_out.rds"))
  flog.info(paste("saving experiment data to", file))
  saveRDS(dataset, file = file)
  
  ## Set up gene lengths for RPKM
  flog.info("outputting counts data")
  write.table(counts(dataset),
              file = paste0(opt$outdir,"/Counts_Data.tsv"),
              sep = "\t",
              quote = FALSE,
              row.names = TRUE,
              col.names = NA)
  if(opt$source %in% c("salmon", "kallisto")){
    write.table(txi$abundance,
                file = paste0(opt$outdir,"/tpm.tsv"),
                sep = "\t",
                quote = FALSE,
                row.names = TRUE,
                col.names = NA)
  }
}

  
main <- function() {
  option_list <- list(
    make_option(
      "--counts-dir",
      dest="counts_dir",
      type="character",
      help=paste("directory containing expression estimates",
                 "from salmon/kallisto/DEXSeq.")
    ),
    make_option(
      "--counts-tsv",
      dest="counts_tsv",
      type="character",
      help=paste("file containing counts generated",
                 "by e.g. featurecounts.")
    ),
    make_option(
      c("-d", "--sampleData"),
      dest="sampleData",
      type="character",
      default = "",
      help=paste("input file with experimental design/sample info")
    ),
    make_option(
      c("--outdir"),
      dest="outdir",
      type="character",
      default = "results.dir",
      help=paste("output directory")
    ),
    make_option(
      "--model",
      dest = "model",
      type = "character",
      default = "~group",
      help = paste("formula for multivariate model")
    ),
    make_option(
      c("-s", "--source"),
      dest="source", 
      type="character",
      default="salmon",
      help=paste("Source of data. Possible options are ",
                 "\"salmon\", \"kallisto\", \"counts_table\", \"dexseq\"")
    ),
    make_option(
      c("--tx2gene"),
      dest="tx2gene", 
      type="character",
      default="transcript2geneMap.tsv",
      help=paste("Path to transcript to gene tsv.")
    ),
    make_option(
      c("--tx2gene_regex"),
      dest="tx2gene_regex", 
      type="character",
      default="None",
      help=paste("Regex/Prefix for removal of certain features from ",
                 "experiment (e.g. removal of spike-ins)")
    ),  
    make_option(
      "--method",
      dest="method", 
      type="character",
      default="deseq2",
      help=paste("differential expression method to apply ")
    ),
    make_option(
      "--filter",
      dest="filter",
      type="logical",
      default=TRUE,
      help=paste("adopt filtering strategy. ",
                 "For EDGER, the default strategy is applied. ",
                 "For DESeq2 basic rowsum filtering < 10 is applied.")
    ),
    make_option(
      "--dexseq-flattened-file",
      dest="flattenedFile",
      type="character",
      help=paste("directory containing flat gtf for dexseq. DEXSeq ",
                 "expects this to be generated by the",
                 "DEXSeq_prepare_annotations.py script")
    )
  )
   
  opt <- experiment_start(option_list = option_list,
                          description = description)
  if (!is.null(opt$method)) {
    opt$method = str_to_lower(opt$method)
  }
  
  if (!is.null(opt$source)) {
    opt$source = str_to_lower(opt$source)
  }
  run(opt)
  
  experiment_stop()
}

main()
