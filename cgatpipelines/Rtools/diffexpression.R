#' filtering single cell data based on QC metrics
#'
#' WARNING: This script is work-in-progress
#' 
#' Example usage:
#' 
#' cgat-singlecell sc_diffexpression --rds-filename=sce.rds --phenotypes-filename=phenodata.tsv --factor=group,mouse_id,collection_date,slice_depth,slice_number,pipette_visual,timepoint > filtered_counts.tsv
#'
#' `sce.rds` is a single cell experiment object after filtering
#'
#'
#' Features can then be selected in the `--factor` option to be
#' plotted.
#'
#' -> todo: parameterize detection of ERCC (pattern?)
#' -> todo: parameterize definition of mitochondrial genes - currently hardcoded for mouse.

## conda dependencies: bioconductor-scater r-cairo

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

source(file.path(Sys.getenv("R_ROOT"), "experiment.R"))


mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "jul2018.archive.ensembl.org")
getmart <- function(values){
  data<- getBM(
    filters= "ensembl_gene_id", 
    attributes= c("ensembl_gene_id", "external_gene_name", "description","entrezgene", 'chromosome_name',
                   'start_position', 'end_position'),
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

# makeTPMtable
# Function: combines data from a genelist with abundance data to create a 
# Input: 1. genelist vector in ENSEMBL format, 2. abundance matrix from tximport
# Output: data frame with gene symbol names
makeTPMtable  <- function(genelist, abundance, design, contrast){
  genelist.df <- getmart(genelist)
  genelist.df <- genelist.df[!duplicated(genelist.df[,1]),]
  rownames(genelist.df) <- genelist.df$ensembl_gene_id
  genelist.names <- genelist.df[genelist,]$external_gene_name
  dftemp <- as_tibble(t(abundance[genelist,]), rownames = "track")
  dftemp <- dftemp %>% rename_at(vars(genelist), ~ genelist.names)
  dftemp$contrast <- design[dftemp$track,][,contrast]
  return(dftemp)
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
  futile.logger::flog.info(paste("Dataset", class(experiment)))
  if (class(experiment) == "DESeqDataSet"){
    flog.info("Running DESeq2")
    suppressMessages(library(DESeq2))
    colData(experiment)[,opt$contrast] <- relevel(colData(experiment)[,opt$contrast], ref = opt$refgroup)
    design(experiment) <- formula(opt$model)
    dds <- DESeq(experiment, betaPrior=FALSE)
    futile.logger::flog.info(paste("coef", opt$coef))
    futile.logger::flog.info(paste("coef", resultsNames(dds)))
    res <- results(dds, name=opt$coef)
    resLFC <- lfcShrink(dds, coef=opt$coef, type=opt$shrinkage)
  } 
  if (class(experiment) == "DGEList"){
    dds = DESeqDataSetFromMatrix(experiment$counts, experiment$sample, design = formula(opt$model))
  } 
  if (class(experiment) == "DEXSeqDataSet"){
    stop("please use diffexonexpression R script")
  }

  flog.info("... plotting dispersion estimates")
  ## Plot dispersion estimates
  start_plot("Dispersion", outdir=opt$outdir)
    plotDispEsts(dds)
  end_plot()
  
  flog.info("... plotting MA")
  ## MA Plot
  start_plot("MAPlot", outdir=opt$outdir)
    DESeq2::plotMA(resLFC, ylim = c(-3,3))
  end_plot()
  
  flog.info("... saving DE data")
  ## Save DE data
  resSig <- subset(resLFC, padj < opt$alpha)
  data <- getmart(rownames(resLFC))
  resLFC$symbol<-data$external_gene_name[match(rownames(resLFC), data$ensembl_gene_id)]
  resLFC$desc<-data$description[match(rownames(resLFC), data$ensembl_gene_id)]
  resLFC$chromosome<-data$chromosome_name[match(rownames(resLFC), data$ensembl_gene_id)]
  resLFC$start<-data$start_position[match(rownames(resLFC), data$ensembl_gene_id)]
  resSig <- subset(resLFC, padj < opt$alpha)
  write.table(resSig, paste0(opt$outdir,"/","results.tsv"), sep = "\t")
  write.table(resLFC, paste0(opt$outdir,"/","results_full.tsv"), sep = "\t")
  resSig2 <- subset(res, padj < opt$alpha)
  write.table(resSig2, paste0(opt$outdir,"/","results_noLFCshrinkage.tsv"), sep = "\t")
  resdf <- data.frame(geneid=rownames(resLFC), pvalue=resLFC$pvalue*sign(resLFC$log2FoldChange))
  write.table(resdf, paste0(opt$outdir,"/","px_results_pvalue.gene.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  resdf <- data.frame(geneid=rownames(resLFC), l2fc=resLFC$log2FoldChange)
  write.table(resdf, paste0(opt$outdir,"/","px_results_l2fc.gene.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  resdf <- data.frame(geneid=rownames(resLFC), l2fc=resLFC$log2FoldChange)
  write.table(resdf, paste0(opt$outdir,"/","px_results_l2fc.gene.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  genelengths <- rowMeans(assays(experiment)$avgTxLength)
  write.table(genelengths, paste0(opt$outdir,"/","gene_lengths.tsv"), sep = "\t", row.names = TRUE, quote = FALSE)
  file = get_output_filename(paste0(opt$outdir,"/","results_table.rds"))
  flog.info(paste("saving results data to", file))
  saveRDS(resLFC, file = file)


  flog.info("... plotting P histogram")
  ## Plot P value Histogram
  start_plot("PHistogram", outdir=opt$outdir)
    hist(res$pvalue,breaks=50, col='skyblue', xlab="p-value", main="P-value Histogram")
  end_plot()

  flog.info("... plotting downregulated genes")    
  ## Plot Top Downregulated Genes
  genelist <- rownames(resLFC[ order( resLFC[,grepl("log2",colnames(resLFC))] ), ][0:9,])
  dftemp <- makeTPMtable(genelist, counts(dds, normalized=TRUE), colData(dds), opt$contrast)
  start_plot("Downregulated", outdir=opt$outdir)
  print(plotTPMs(dftemp, opt$contrast))
  end_plot()
  
  flog.info("... plotting upregulated genes")
  genelist <- rownames(resLFC[ order( -resLFC[,grepl("log2",colnames(resLFC))] ), ][0:9,])
  dftemp <- makeTPMtable(genelist, counts(dds, normalized=TRUE), colData(dds), opt$contrast)
  start_plot("Upregulated", outdir=opt$outdir)
    print(plotTPMs(dftemp, opt$contrast))
  dev.off()
  
  flog.info("... plotting significant genes")
  resSig <- subset(resLFC, padj < 0.05)
  genelist <- rownames(resSig[order(resSig$padj), ])
  if(length(genelist) > 9){
      genelist <- genelist[0:9]}
  dftemp <- makeTPMtable(genelist, counts(dds, normalized=TRUE), colData(dds), opt$contrast)
  start_plot("significant", outdir=opt$outdir)
    print(plotTPMs(dftemp, opt$contrast))
  dev.off()

  flog.info("... plotting user-defined genes")
  genelist <- unlist(strsplit(opt$userlist, ","))
  dftemp <- makeTPMtable(genelist, counts(dds, normalized=TRUE), colData(dds), opt$contrast)
  start_plot("Userdefined", outdir=opt$outdir)
    print(plotTPMs(dftemp, opt$contrast))
  dev.off()
  
  flog.info("... performing gene ontology (GO) enrichment analysis")
  res.nona <- subset(resLFC, (!is.na(padj)))
  de.genes <- as.integer(res.nona$padj < 0.05)
  names(de.genes) <- rownames(res.nona)
  temp <- rowMeans(assays(experiment)$avgTxLength[match(names(de.genes), rownames(assays(experiment)$avgTxLength)),])
  pwf=nullp(de.genes,bias.data=temp)
  all = goseq(pwf,"hg38","ensGene", method="Hypergeometric")
  sigall <- all
  names(sigall) <- c("category","pvalue","underrepresented_pvalue","numberDE", "numberTOT", "term", "ontology")
  sigall$pvalue <- p.adjust(sigall$pvalue, method="BH")
  sigall$percent <- sigall$numberDE/sigall$numberTOT
  sigall <- sigall[which(sigall[,2] < 0.05),]
  cats <- sigall$category
  write.table(sigall[,c("category", "pvalue")], file=paste0(opt$outdir,"/",'GO.tsv'), quote=FALSE, sep='\t', row.names = FALSE)
  write.table(sigall, file=paste0(opt$outdir,"/",'GOall.tsv'), quote=FALSE, sep='\t', row.names = FALSE)
  
  
  
  flog.info("... performing gene set enrichment analysis (GSEA)")
  resrnk<-resLFC
  data <- getmart(rownames(resrnk))
  resrnk$symbol <- data$external_gene_name[match(rownames(resrnk), data$ensembl_gene_id)]
  resrnk$desc <-data$description[match(rownames(resrnk), data$ensembl_gene_id)]
  resrnk$entrezgene <-data$entrezgene[match(rownames(resrnk), data$ensembl_gene_id)]
  rnk.df <- as_tibble(resrnk[,c("entrezgene","log2FoldChange")]) %>% na.omit()
  rnk <- rnk.df$log2FoldChange
  names(rnk) <- rnk.df$entrezgene
  rnk <- rnk[isUnique(names(rnk))]
  
  for(pathway in opt$pathways){
      pathways <- gmtPathways(pathway)
      fgseaRes <- fgsea(pathways = pathways, 
                        stats = rnk,
                        minSize=15,
                        maxSize=500,
                        nperm=10000)
      fwrite(fgseaRes, file=paste0(opt$outdir,"/",'GSEA_',sub("([^.]+)\\.[[:alnum:]]+$", "\\1", (basename(pathway))),'.tsv'), sep="\t", sep2=c("", " ", ""))
      topPathwaysUp <- fgseaRes[ES > 0,][head(order(pval), n=10),]$pathway
      png(paste0(opt$outdir,"/",'GSEA_GO_up_',sub("([^.]+)\\.[[:alnum:]]+$", "\\1", (basename(pathway))),'.png'),
          width =15, height = 3, units = 'in', res = 600)
      plotGseaTable(pathways[topPathwaysUp], rnk, fgseaRes, 
                    gseaParam = 0.5, colwidths = c(10,2,1,1,1))
      dev.off()
  
  
      topPathwaysDown <- fgseaRes[ES < 0,][head(order(pval), n=10),]$pathway
      plotGseaTable(pathways[topPathwaysDown], rnk, fgseaRes, 
      gseaParam = 0.5)
      png(paste0(opt$outdir,"/",'GSEA_GO_down_',sub("([^.]+)\\.[[:alnum:]]+$", "\\1", (basename(pathway))),'.png'),
          width =15, height = 3, units = 'in', res = 600)
      plotGseaTable(pathways[topPathwaysDown], rnk, fgseaRes, 
                    gseaParam = 0.5, colwidths = c(10,2,1,1,1))
      dev.off()
  }
  
  if(opt$perm >0){
    flog.info("... performing permutations")
    designperm <- colData(dds)
    ddsperm = dds
    x = 0
    i = 1
    y= list()
    while(i <= opt$perm) {
      designperm[, opt$contrast] = as.factor(sample(levels(designperm[, opt$contrast]), length(designperm[, opt$contrast]), replace=TRUE, prob=c(0.5, 0.5)))
      if(sum(designperm$group == levels(designperm[, opt$contrast])[1]) != table(colData(dds)[,opt$contrast])[1]) next
      colData(ddsperm) <- designperm
      ddsperm <- DESeq(ddsperm)
      resperm <- results(ddsperm)
      x[i] = length(subset(resperm, padj < 0.05)$padj)
      y[[i]] = ddsperm$group
      i = i+1
    }
    start_plot("Simulations", outdir=opt$outdir)
    theme_set(theme_gray(base_size = 18))
    sims = qplot(x,
                 geom="histogram",
                 breaks=seq(0, 20, by = 1),fill=I("grey"), col=I("black"),
                 main = "Histogram of DE experiments\n with random group labels", 
                 xlab = "Number of differentially expressed genes",
                 ylab = "Number of simulations") +
      geom_vline(xintercept = length(rownames(subset(res, padj < 0.05)))) +
      theme_classic() + theme(plot.title = element_text(hjust = 0.5, size=22))
    print(sims)
    end_plot()
    flog.info(paste0("... Permutation p value: ",
                     length(x[x < length(rownames(subset(res, padj < 0.05)))])/length(x)))
    z <- list()
    for(i in 0:length(x)){
      z[i] <- paste( unlist(y[i]), collapse=' ')
    }
    z <- unlist(z)
    df <- data.frame(number=x, combination=z)
    df
    write_tsv(df,paste0(opt$outdir,"/","Simulations.tsv"))
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
            default = "~group",
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
          "--coef",
          dest = "coef",
          type = "character",
          default = "",
          help = paste("Comparison of interest for DESeq2: e.g. group_MUT_vs_CTR, where MUT and CTR are in group")
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
            "--shrinkage",
            dest = "shrinkage",
            type = "character",
            default = "apeglm",
            help = paste("Method for LFC shrinkage. Options are apeglm (default), ashr, normal")
        ),
         make_option(
            "--userlist",
            dest = "userlist",
            type = "character",
            default = "ENSG00000205927,ENSG00000130675,ENSG00000016082,ENSG00000070748,ENSG00000064300,ENSG00000078018",
            help = paste("User defined genes to plot")
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
            type = "integer",
            default = 0,
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
