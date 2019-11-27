#' read counts and design data from tab-separated files
#' and returns an Experiment object.
#' 


read_single_cell_experiment_from_tables <- function(counts_filename,
                                                    design_filename,
                                                    remove_all_zero = TRUE) {
   
    futile.logger::flog.info(paste("reading counts data from", normalizePath(counts_filename)))

    counts_data <-read.table(counts_filename, header = TRUE, row.names = 1, sep = "\t")
    futile.logger::flog.info(paste("read counts data", paste(dim(counts_data), collapse = ",")))

    futile.logger::flog.info(paste("reading design data from", normalizePath(design_filename)))
    annotation_data <- read.table(design_filename, sep = "\t", header = TRUE)
    futile.logger::flog.info(paste("read design data", paste(dim(annotation_data), collapse = ",")))
    
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
    is.spike <- (rownames(all_sceset) %in% ercc)
    isSpike(all_sceset, type="ERCC") <- is.spike        

    is.mito <- (rownames(all_sceset) %in% mt)
    isSpike(all_sceset, type="Mt") <- is.mito
    
    futile.logger::flog.info(paste("marking", sum(is.spike), "spike-in genes"))
    futile.logger::flog.info(paste("marking", sum(is.mito), "mitochondrial genes"))

    assay(all_sceset, "logcounts_raw") <- log2(counts(all_sceset) + 1)
    return(all_sceset)
}


read_experiment_from_rds <- function(rds_filename,
                                     remove_all_zero = FALSE) {
   
    futile.logger::flog.info(paste("reading Experiment object from", normalizePath(rds_filename)))
    experiment <- readRDS(rds_filename)
    futile.logger::flog.info(paste("read Experiment object", paste(dim(counts(experiment)), collapse = ",")))

    if (remove_all_zero) {
        futile.logger::flog.info("removing genes not expressed in any cell")
        keep_feature <- rowSums(counts(experiment) > 0) > 0
        futile.logger::flog.info(paste("keeping", sum(keep_feature), "genes"))
        experiment <- experiment[keep_feature, ]
    }
    return(experiment)
}
