#' Differential expression analysis Script
#'
#' 
#' Example usage:
#' 
#' cgatflow R exploratory --model=~ group --contrast=group --factor=mouse_id,collection_date,slice_depth,slice_number,pipette_visual,timepoint
#'
#' Primary Input: `filename.rds` - output from the readandfiltercounts.R file or DESeqExpermient object or DEGList object
#' Additional inputs: model and contrast as well as factors of interest
#'
#' Output: png images of: PCA, clustering, Surrogate Variable Analysis, Heatmaps
#'
suppressMessages(library(getopt))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(rhdf5))
suppressMessages(library(sva))
suppressMessages(library(ggplot2))
suppressMessages(library(biomaRt))
suppressMessages(library(Cairo))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggforce))


source(file.path(Sys.getenv("R_ROOT"), "io.R"))
source(file.path(Sys.getenv("R_ROOT"), "experiment.R"))


# getmart
# Function: wrapper for getBM for ENSEMBL gene list
# Dependencies: biomaRt package
# Input: vector of ensembl gene ids
# Output: data frame with description, symbol and entrez gene
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


start_plot <- function(section, height = 6, width = 6, type = "png", outdir="") {
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


run <- function(opt) {

  ### READING DATA ###
  futile.logger::flog.info(paste("reading Experiment object from", normalizePath(opt$rds_filename)))
  experiment <- readRDS(opt$rds_filename)
  if (class(experiment) == "DESeqDataSet"){
    dds <- experiment
  } else if (class(experiment) == "DGEList"){
    dds = DESeqDataSetFromMatrix(experiment$counts, experiment$sample, design = formula(opt$model))
  }
  futile.logger::flog.info(paste("reading  Experiment object", paste(dim(counts(experiment)), collapse = ",")))
  
  ### SVA - ANALYSIS OF BIASES ###
  futile.logger::flog.info(paste("Performing Surrogate Variable Analysis"))
  dds <- estimateSizeFactors(dds)
  dat  <- counts(dds, normalized = TRUE)
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  mod  <- model.matrix(formula(opt$model), colData(dds))
  mod0 <- model.matrix(~ 1, colData(dds))
  svseq <- svaseq(dat, mod, mod0, n.sv = 2)
  for(factor in opt$factors){
    start_plot(paste0('SVA for ', factor), outdir=opt$outdir)
    par(mfrow = c(2, 1), mar = c(3,5,3,1))
    for (i in 1:2) {
      stripchart(svseq$sv[, i] ~ colData(dds)[, factor], vertical = TRUE, main = paste0("SV", i))
      abline(h = 0)
    }
    end_plot()
  }
  
  ### TRANSFORMATION OF DATA ###
  futile.logger::flog.info(paste("Transforming data"))
  rld<- rlog(dds)
  vsd<- vst(dds)
  df <- bind_rows(
    as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
      mutate(transformation = "log2(x + 1)"),
    as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
    as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
  colnames(df)[1:2] <- c("x", "y")
  start_plot('Variance_Transformations', outdir=opt$outdir)
    print(ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
      coord_fixed() + facet_grid( . ~ transformation))
  end_plot()
  
  ### PRINCIPAL COMPONENT ANALYSIS ###
  futile.logger::flog.info(paste("Performing Principal Component Analysis"))
  pca = prcomp(t(assay(vsd)))
  for(factor in opt$factors){
    variable.group <- colData(dds)[, factor]
    names(variable.group) <- factor
    percentVar <- round(100 * summary(pca)$importance[2,])
    futile.logger::flog.info(paste("PCA", dim(pca$x)))
    scores <- data.frame(variable.group, pca$x[,1:8])
    start_plot(paste0('PCA_', factor), outdir=opt$outdir)
      print(qplot(x=PC1, y=PC2, data=scores, colour=factor(variable.group)) +
        theme(legend.position="right") +  
        labs(colour=factor, x=paste0("PC1 (", percentVar[1],"% of variance)"),
             y=paste0("PC2 (", percentVar[2],"% of variance)")) + 
        ggtitle("Principal Component Analysis") + theme_grey(base_size = 15) +
        theme(plot.title = element_text(lineheight=1, face="bold"))  + geom_point(size=2) +
        theme(text=element_text(family='serif')))
    end_plot()
    start_plot(paste0('PCA_13_', factor), outdir=opt$outdir)
          print(qplot(x=PC1, y=PC3, data=scores, colour=factor(variable.group)) +
        theme(legend.position="right") +  
        labs(colour=factor, x=paste0("PC1 (", percentVar[1],"% of variance)"),
             y=paste0("PC3 (", percentVar[3],"% of variance)")) + 
        ggtitle("Principal Component Analysis") + theme_grey(base_size = 15) +
        theme(plot.title = element_text(lineheight=1, face="bold"))  + geom_point(size=2) +
        theme(text=element_text(family='serif')))
    end_plot()
  }
  variable.group <- colData(dds)[, opt$contrast]
  names(variable.group) <- opt$contrast
  scores <- data.frame(variable.group, pca$x[,1:8])
  start_plot(paste0('PCA_grid'), outdir=opt$outdir)
  print(ggplot(scores, aes(x = .panel_x, y = .panel_y, fill = variable.group, colour = variable.group)) + 
    geom_point(shape = 16, size = 0.5, position = 'auto') + 
    geom_autodensity(alpha = 0.3, colour = NA, position = 'identity') + 
    facet_matrix(vars(PC1:PC8), layer.diag = 2))
  end_plot()
  
  loadings <- pca$rotation[,1:8]
  loadings <- data.frame(loadings[order(loadings[,1]), ])
  data <- getmart(rownames(loadings))
  loadings$symbol<-data$mgi_symbol[match(rownames(loadings), data$ensembl_gene_id)]
  loadings$description <- data$description[match(rownames(loadings), data$ensembl_gene_id)]
  write.table(loadings, file=paste0(opt$outdir,"/",'PCA_loadings.tsv'), quote=FALSE, sep='\t', row.names = TRUE, col.names = NA) 
  
  ### HEATMAPS ###
  futile.logger::flog.info(paste("Performing Heatmap"))
  df <- as.data.frame(colData(dds)[,opt$factors])
  rownames(df) <- colData(dds)$track
  # Heatmap of Top 20 Expressed Genes
  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
  start_plot('Heatmap_topExpressed', outdir=opt$outdir)
    pheatmap(assay(vsd)[select,], cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, annotation_col=df)
  end_plot()
  # Heatmap of Top 20 Variable Genes
  topVarGenes <- head(order(rowVars(assay(vsd)),decreasing=TRUE),20)
  mat <- assay(vsd)[ topVarGenes, ]
  temp <- getmart(rownames(mat))
  row.names(temp) <- temp$ensembl_gene_id
  rownames(mat) <- temp[rownames(mat),"external_gene_name"]
  start_plot('Heatmap_topVariable', outdir=opt$outdir)
    pheatmap(mat, annotation_col=df, cluster_rows=FALSE,fontsize_row = 6)
  end_plot()
  # Heatmap of Genes of interest
  if (!is.null(opt$genes_of_interest)) {
    mat <- assay(vsd)[opt$genes_of_interest, ]
    mat <- mat - rowMeans(mat)
    temp <- getmart(rownames(mat))
    row.names(temp) <- temp$ensembl_gene_id
    rownames(mat) <- temp[rownames(mat),"external_gene_name"]
    start_plot('Heatmap_ofInterest', outdir=opt$outdir)
      pheatmap(mat, annotation_col=df, scale = "row", cluster_cols = FALSE)
    end_plot()
  }

  ### CLUSTERING ###
  futile.logger::flog.info(paste("Clustering"))
  # for all genes
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste(colData(vsd)[,opt$contrast], vsd$track, sep = " - " )
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  start_plot('Heatmap_all', outdir=opt$outdir)
    pheatmap(sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors)
  end_plot()
  # for top 500 genes
  top500 <- head(assay(vsd)[ order(rowMeans(assay(vsd)), decreasing=TRUE ),  ], n=500)
  distsRL500 <- dist(t(top500))
  sampleDistMatrix500 <- as.matrix( distsRL500 )
  rownames(sampleDistMatrix500) <- paste(colData(vsd)[,opt$contrast], vsd$track, sep = " - " )
  colnames(sampleDistMatrix500) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  start_plot('Heatmap_top500', outdir=opt$outdir)
    pheatmap(sampleDistMatrix500,
    clustering_distance_rows = distsRL500,
    clustering_distance_cols = distsRL500,
    col = colors)
  end_plot()
  
  ### EXPLORE BATCH EFFECTS ###
  futile.logger::flog.info(paste("Exploring Batch Effects"))
  for(factor in opt$factors){
    factor_transformed <- vsd
    assay(factor_transformed) <- limma::removeBatchEffect(assay(factor_transformed), colData(factor_transformed)[,factor])
    pca = prcomp(t(assay(factor_transformed)))
    sample.group <- as_factor(colData(dds)[, factor])
    variable.group <- colData(dds)[, opt$contrast]
    percentVar <- round(100 * summary(pca)$importance[2,])
    scores <- data.frame(variable.group, sample.group, pca$x[,1:8])
    start_plot(paste0('PCA_', factor, '_removed'), outdir=opt$outdir)
    print(qplot(x=PC1, y=PC2, data=scores, colour=factor(variable.group), shape=factor(sample.group)) +
            theme(legend.position="right") +  
            labs(colour=opt$contrast, shape=factor, x=paste0("PC1 (", percentVar[1],"% of variance)"),
                 y=paste0("PC2 (", percentVar[2],"% of variance)")) + 
            ggtitle(paste0("Principal Component Analysis\n after batch correction for ", factor)) + 
            theme_grey(base_size = 15) +
            theme(plot.title = element_text(lineheight=1, face="bold"))  + geom_point(size=2) +
            theme(text=element_text(family='serif')))
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
      help = paste("formula for multivariate model")
    ),
    make_option(
      "--contrast",
      dest = "contrast",
      type = "character",
      default = "group",
      help = paste("formula for multivariate model")
    ),
    make_option(
      "--factors",
      dest = "factors",
      type = "character",
      default = "",
      help = paste("formula for multivariate model")
    ),
    make_option(
      "--genes_of_interest",
      dest = "genes_of_interest",
      type = "character",
      default = NULL,
      help = paste("genes of interest for plotting")
    ),
    make_option(
      c("-o","--outdir"),
      dest = "outdir",
      type = "character",
      default = "",
      help = paste("genes of interest for plotting")
    )
  )
  opt <- experiment_start(option_list = option_list,
                          description = description)
  
  if (!is.null(opt$factors)) {
    opt$factors = unlist(strsplit(opt$factors, ","))
  }
  if (!is.null(opt$genes_of_interest)) {
    opt$genes_of_interest = unlist(strsplit(opt$genes_of_interest, ","))
  }
  run(opt)
  
  experiment_stop()
}

main()
