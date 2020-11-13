#' Running DEXSeq
#'
#' WARNING: This script is work-in-progress
#' 
#' Example usage:
#' 
#' Rscript PATH/TO/diffexonexpression.R --rds-filename=experiment.rds --model=~group --factor=group,mouse_id,collection_date,slice_depth,slice_number,pipette_visual,timepoint
#'
#' `experiment.rds` is a DEXSeq experiment object after filtering
#'
#'
#' Features can then be selected in the `--factor` option to be
#' plotted.
#'
#' -> todo: ideas

## conda dependencies: r-cairo

suppressMessages(library(futile.logger))
suppressMessages(library(getopt))
suppressMessages(library(fgsea))
suppressMessages(library(data.table))
suppressMessages(library(sva))	
suppressMessages(library(biomaRt))
suppressMessages(library(tidyverse))
suppressMessages(library(RUVSeq))
suppressMessages(library(Cairo))
suppressMessages(library(goseq))
suppressMessages(library(DEXSeq))

source(file.path(Sys.getenv("R_ROOT"), "experiment.R"))


mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "jul2018.archive.ensembl.org")
getmart <- function(values){
  data<- getBM(
    filters= "ensembl_gene_id", 
    attributes= c("ensembl_gene_id", "external_gene_name", "description","entrezgene"),
    values= values,
    mart= mart)
  data$description <- gsub("\t", "", data$description)
  return(data)
}

start_plot <- function(section, outdir = "", height = 6, width = 6, type = "png") {
  file = get_output_filename(paste0(outdir,"/",section, ".", type))
  Cairo(file = file,
        type = type,
        width = width,
        height = height,
        units="in",
        dpi = 300,
        bg = "white")
  #opar <- par(lwd=0.5)
}

end_plot <- function() {
  dev.off()
}


# plotTPMs function
# Function: plots a table of TPMs from a data frame, ensuring appropriate formating for ggplot
# Input: data frame
# Output: plot object
plotTPMs <- function(dftemp, contrast_name){
  dftemp %>% 
    gather(key = "var", value="value", -contrast, -track) %>% 
    mutate(var = factor(var, levels=unique(var))) %>%
    ggplot(aes(x = contrast, y = value, color = contrast)) +
    geom_point(position = position_jitter(w = 0.15, h = 0)) +
    facet_wrap(~ var, scales = "free") + theme_bw() +
    ylab("normalised counts") + xlab (contrast_name) + guides(color = "none")
}


run <- function(opt) {
  
  futile.logger::flog.info(paste("reading Experiment object from", normalizePath(opt$rds_filename)))
  experiment <- readRDS(opt$rds_filename)
  if (class(experiment) != "DEXSeqDataSet"){
    stop("This script supports DEXSeqDataSet objects only.")
  }
  flog.info("Running DEXSeq")
  dxd = estimateSizeFactors(experiment)
  dxd = estimateDispersions(dxd)
  dxd = testForDEU(dxd)
  dxd = estimateExonFoldChanges( dxd, fitExpToVar=opt$contrast)
  res = DEXSeqResults(dxd)
  
  flog.info("... plotting dispersion estimates")
  ## Plot dispersion estimates
  start_plot("Dispersion", opt$outdir)
  plotDispEsts(dxd)
  end_plot()
  
  flog.info("... plotting MA")
  ## MA Plot
  start_plot("MAPlot", opt$outdir)
  plotMA(res, ylim = c(-3,3))
  end_plot()
  
  flog.info("... saving DE data")
  ## Save DE data
  #resSig <- subset(res, padj < opt$alpha)
  #data <- getmart(as.character(map(res$transcripts, 1)))
  #res$symbol<-data$external_gene_name[match(data$ensembl_transcript_id, as.character(map(res$transcripts, 1)))]
  #res$desc<-data$description[match(data$ensembl_transcript_id, as.character(map(res$transcripts, 1)))]
  resSig <- subset(res, padj < opt$alpha)
  write.table(resSig, paste0(opt$outdir,"/","results.tsv"), sep = "\t")
  write.table(res, paste0(opt$outdir,"/","results_full.tsv"), sep = "\t")
  
  flog.info("... plotting P histogram")
  ## Plot P value Histogram
  start_plot("PHistogram", opt$outdir)
  hist(res$pvalue,breaks=50, col='skyblue', xlab="p-value", main="P-value Histogram")
  end_plot()
  
  flog.info("... plotting downregulated genes")    
  ## Plot Top Downregulated Genes
  genelist <- rownames(res[ order( res[,grepl("log2",colnames(res))] ), ][0:9,])
  for(gene in unique(substr(genelist,0,15))){
    start_plot(paste0("Downregulated_",gene), opt$outdir)
      plotDEXSeq(res, gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 ,fitExpToVar = opt$contrast)
    end_plot()
  }

  
  flog.info("... plotting upregulated genes")
  genelist <- rownames(res[ order( -res[,grepl("log2",colnames(res))] ), ][0:9,])
  for(gene in unique(substr(genelist,0,15))){
    start_plot(paste0("Upregulated_",gene), opt$outdir)
    plotDEXSeq(res, gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 ,fitExpToVar = opt$contrast)
    end_plot()
  }

  flog.info("... plotting significant genes")
  genelist <- rownames(resSig[order(resSig$padj), ])
  for(gene in unique(substr(genelist,0,15))){
    start_plot(paste0("Significant_",gene), opt$outdir)
    plotDEXSeq(res, gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 ,fitExpToVar = opt$contrast)
    end_plot()
  }
}

main <- function() {
  
  option_list <- list(
    make_option(
      "--rds-filename",
      dest = "rds_filename",
      type = "character",
      default = "experiment.rds",
      help = paste("filename with input data of counts")
    ),
    make_option(
      "--model",
      dest = "model",
      type = "character",
      default = "~ group",
      help = paste("model formula for GLM")
    ),
    make_option(
      "--contrast",
      dest = "contrast",
      type = "character",
      default = "group",
      help = paste("contrast/factor to use for comparison")
    ),
    make_option(
      "--refgroup",
      dest = "refgroup",
      type = "character",
      default = "CTR",
      help = paste("reference group, e.g. CTR")
    ),
    make_option(
      "--alpha",
      dest = "alpha",
      type = "numeric",
      default = 0.05,
      help = paste("Adjusted P value threshold")
    ),
   make_option(
        "--outdir",
        dest = "outdir",
        type = "character",
        default = "experiment",
        help = paste("Libraries for gsea")
    ),
    make_option(
      "--permute",
      dest = "perm",
      type = "logical",
      default = FALSE,
      help = paste("number of permutations")
    ),
    make_option(
      "--pathways",
      dest = "pathways",
      type = "character",
      default = "",
      help = paste("Libraries for gsea")
    )
  )
  opt <- experiment_start(option_list = option_list,
                          description = description)
  
  if (!is.null(opt$pathways)) {
    opt$pathways = unlist(strsplit(opt$pathways, ","))
  }
  run(opt)
  
  experiment_stop()
}

main()
