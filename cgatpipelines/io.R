#' read counts and phenotype data from tab-separated files
#' and return a SingleCellExperiment object.
#' 


read_single_cell_experiment_from_tables <- function(counts_filename,
                                                    phenotypes_filename,
                                                    remove_all_zero = TRUE) {
   
    futile.logger::flog.info(paste("reading counts data from", normalizePath(counts_filename)))

    counts_data <-read.table(counts_filename, header = TRUE, row.names = 1, sep = "\t")
    futile.logger::flog.info(paste("read counts data", paste(dim(counts_data), collapse = ",")))

    futile.logger::flog.info(paste("reading phenotype data from", normalizePath(phenotypes_filename)))
    annotation_data <- read.table(phenotypes_filename, sep = "\t", header = TRUE)
    futile.logger::flog.info(paste("read phenotype data", paste(dim(annotation_data), collapse = ",")))
    
    annotation_data$sample_id <- gsub("-", ".", annotation_data$sample_id)
    row.names(annotation_data) <- annotation_data$sample_id

    futile.logger::flog.info("building SingleCellExperiment data set")
    all_sceset <- SingleCellExperiment::SingleCellExperiment(
                                            assays = list(counts = as.matrix(counts_data)),
                                            colData = annotation_data)
    futile.logger::flog.info(paste("built SingleCellExperiment data set", paste(dim(all_sceset), collapse=",")))

    if (remove_all_zero) {
        futile.logger::flog.info("removing genes not expressed in any cell")
        keep_feature <- rowSums(counts(all_sceset) > 0) > 0
        futile.logger::flog.info(paste("keeping", sum(keep_feature), "genes"))
        all_sceset <- all_sceset[keep_feature, ]
    }
    
    ercc <- rownames(all_sceset)[grep("ERCC", rownames(all_sceset))]
    mt <- c("ENSMUSG00000064336","ENSMUSG00000064337","ENSMUSG00000064338",
            "ENSMUSG00000064339","ENSMUSG00000064340","ENSMUSG00000064341",
            "ENSMUSG00000064342","ENSMUSG00000064343","ENSMUSG00000064344",
            "ENSMUSG00000064345","ENSMUSG00000064346","ENSMUSG00000064347",
            "ENSMUSG00000064348","ENSMUSG00000064349","ENSMUSG00000064350",
            "ENSMUSG00000064351","ENSMUSG00000064352","ENSMUSG00000064353",
            "ENSMUSG00000064354","ENSMUSG00000064355","ENSMUSG00000064356",
            "ENSMUSG00000064357","ENSMUSG00000064358","ENSMUSG00000064359",
            "ENSMUSG00000064360","ENSMUSG00000064361","ENSMUSG00000064363",
            "ENSMUSG00000064364","ENSMUSG00000064365","ENSMUSG00000064366",
            "ENSMUSG00000064367","ENSMUSG00000064368","ENSMUSG00000064369",
            "ENSMUSG00000064370","ENSMUSG00000064371","ENSMUSG00000064372",
            "ENSMUSG00000065947")

    is.spike <- (rownames(all_sceset) %in% ercc)
    isSpike(all_sceset, type="ERCC") <- is.spike        

    is.mito <- (rownames(all_sceset) %in% mt)
    isSpike(all_sceset, type="Mt") <- is.mito
    
    futile.logger::flog.info(paste("marking", sum(is.spike), "spike-in genes"))
    futile.logger::flog.info(paste("marking", sum(is.mito), "mitochondrial genes"))

    assay(all_sceset, "logcounts_raw") <- log2(counts(all_sceset) + 1)
    return(all_sceset)
}

read_single_cell_experiment_from_rds <- function(rds_filename,
                                                 remove_all_zero = TRUE) {
   
    futile.logger::flog.info(paste("reading single cell experiment object from", normalizePath(rds_filename)))
    sce <- readRDS(rds_filename)
    futile.logger::flog.info(paste("read single cell experiment object", paste(dim(counts(sce)), collapse = ",")))

    if (remove_all_zero) {
        futile.logger::flog.info("removing genes not expressed in any cell")
        keep_feature <- rowSums(counts(sce) > 0) > 0
        futile.logger::flog.info(paste("keeping", sum(keep_feature), "genes"))
        sce <- sce[keep_feature, ]
    }
    return(sce)
}
