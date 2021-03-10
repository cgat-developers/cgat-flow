#' differential expression of IRFinder 
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


source(file.path("/ifs/projects/jakubs/cgat-developers-v2/cgat-flow/cgatpipelines" , "experiment.R"))


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


DESeqDataSetFromIRFinder = function(filePaths,designMatrix,designFormula){
    res=c()
    libsz=c()
    spl=c()
    irtest=read.table(filePaths[1])
    if (irtest[1,1]=="Chr"){irtest=irtest[-1,]}
    irnames=unname(apply(as.matrix(irtest),1,FUN=function(x){return(paste0(x[4],"/",x[1],":",x[2],"-",x[3],":",x[6]))}))
    n=1
    for (i in filePaths){
        print(paste0("processing file ",n," at ",i))
        irtab=read.table(i)
        if (irtab[1,1]=="Chr"){irtab=irtab[-1,]}
        #rn=unname(apply(irtab,1,FUN=function(x){return(paste0(x[4],"/",x[1],":",x[2],"-",x[3],":",x[6]))}))
        #row.names(irtab)=rn
        #tmp1=round(as.numeric(as.vector(irtab[irnames,9])))
        #tmp2=as.numeric(as.vector(irtab[irnames,19]))
        tmp1=as.numeric(as.vector(irtab[,9]))
        tmp2=as.numeric(as.vector(irtab[,19]))
        tmp3=tmp1+tmp2
        tmp4=as.numeric(as.vector(irtab[,17]))
        tmp5=as.numeric(as.vector(irtab[,18]))
        tmp6=pmax(tmp4,tmp5, na.rm=T)
        res=cbind(res,tmp1)
        libsz=cbind(libsz,tmp2)
        spl=cbind(spl,tmp6)
        n=n+1
    }
    res.rd=round(res)
    libsz.rd=round(libsz)
    spl.rd=round(spl)
    colnames(res.rd)=paste("intronDepth",as.vector(designMatrix[,1]),sep=".")
    rownames(res.rd)=irnames
    colnames(libsz.rd)=paste("totalSplice",as.vector(designMatrix[,1]),sep=".")
    rownames(libsz.rd)=irnames
    colnames(spl.rd)=paste("maxSplice",as.vector(designMatrix[,1]),sep=".")
    rownames(spl.rd)=irnames
    
    ir=c(rep("IR",dim(designMatrix)[1]),rep("Splice",dim(designMatrix)[1]))
    group=rbind(designMatrix,designMatrix)
    group$IRFinder=ir
    group$IRFinder=factor(group$IRFinder,levels=c("Splice","IR"))
    
    #counts.IRFinder=cbind(res.rd,libsz.rd)
    counts.IRFinder=cbind(res.rd,spl.rd)
    
    dd = DESeqDataSetFromMatrix(countData = counts.IRFinder, colData = group, design = designFormula)
    sizeFactors(dd)=rep(1,dim(group)[1])
    rownames(dd)=irnames
    final=list(dd,res,libsz,spl)
    names(final)=c("DESeq2Object","IntronDepth","SpliceDepth","MaxSplice")
    return(final)
}


run <- function(opt) {
    library(DESeq2)
    futile.logger::flog.info(paste("reading Experiment object from", normalizePath(opt$counts_file)))
    results = read.table(opt$counts_file)
    paths = as.vector(results$V1)                                            # File names must be saved in a vector
    sampleData = read.table(opt$design,header=T)                       
    sampleData$group=factor(sampleData$group)    # Set WT as the baseline in the analysis
    rownames(sampleData)=NULL 
    metaList=DESeqDataSetFromIRFinder(filePaths=paths, designMatrix=sampleData, designFormula=~1)

    experiment = metaList$DESeq2Object                       # Extract DESeq2 Object with normalization factors ready
    colData(experiment)[,opt$contrast] <- relevel(colData(experiment)[,opt$contrast], ref = opt$refgroup)
    design(experiment) <- formula(opt$model)     # Build a formula of GLM. Read below for more details. 
    dds = DESeq(experiment)                                  # Estimate parameters and fit to model
    futile.logger::flog.info(resultsNames(dds))                                 # Check the actual variable name assigned by DESeq2

    REF = paste0(opt$contrast,opt$refgroup,".IRFinderIR")
    res.REF = results(dds, name = REF)
    REF.IR_vs_Splice=2^res.REF$log2FoldChange
    IRratio.REF = REF.IR_vs_Splice/(1+REF.IR_vs_Splice)

    COMP = paste0(opt$contrast,opt$compgroup,".IRFinderIR")
    res.COMP = results(dds, name = COMP)
    COMP.IR_vs_Splice=2^res.COMP$log2FoldChange
    IRratio.COMP = COMP.IR_vs_Splice/(1+COMP.IR_vs_Splice)

    # Finally we can test the difference of (intronic reads/normal spliced reads) ratio between WT and KO
    res.diff = results(dds, contrast=list(COMP,REF))
    resSig <- subset(res.diff, padj < opt$alpha)
    write.table(resSig, paste0(opt$outdir,"/","results.tsv"), sep = "\t")
    write.table(res.diff, paste0(opt$outdir,"/","results_full.tsv"), sep = "\t")
    
    # We can plot the changes of IR ratio with p values
    # In this example we defined significant IR changes as
    # 1) IR changes no less than 10% (both direction) and 
    # 2) with adjusted p values less than 0.05

    IR.change = IRratio.COMP - IRratio.REF
    write.table(IR.change, paste0(opt$outdir,"/","IR_change.tsv"), sep = "\t")


    start_plot("Dispersion", outdir=opt$outdir)
        plot(IR.change,col=ifelse(res.diff$padj < 0.05 & abs(IR.change)>=0.1, "red", "black"))
    end_plot()


    if(opt$perm >0){
        flog.info("... performing permutations")
        designperm <- colData(experiment)
        experiment_perm = experiment
        x = 0
        i = 1
        y= list()
        while(i <= opt$perm) {
            designperm[, opt$contrast] = as.factor(sample(levels(designperm[, opt$contrast]), length(designperm[, opt$contrast]), replace=TRUE))
            if(sum(designperm$group == levels(designperm[, opt$contrast])[1]) != table(colData(dds)[,opt$contrast])[1]) next
            colData(experiment_perm) <- designperm
            ddsperm <- DESeq(experiment_perm)
            REF = paste0(opt$contrast,opt$refgroup,".IRFinderIR")
            res.REF = results(ddsperm, name = REF)
            REF.IR_vs_Splice=2^res.REF$log2FoldChange
            IRratio.REF = REF.IR_vs_Splice/(1+REF.IR_vs_Splice)

            COMP = paste0(opt$contrast,opt$compgroup,".IRFinderIR")
            res.COMP = results(ddsperm, name = COMP)
            COMP.IR_vs_Splice=2^res.COMP$log2FoldChange
            IRratio.COMP = COMP.IR_vs_Splice/(1+COMP.IR_vs_Splice)

            # Finally we can test the difference of (intronic reads/normal spliced reads) ratio between WT and KO
            res.diff_tmp = results(ddsperm, contrast=list(COMP,REF))
            x[i] = length(subset(res.diff_tmp, padj < 0.05)$padj)
            y[[i]] = ddsperm$group
            i = i+1
        }
        start_plot("Simulations", outdir=opt$outdir)
        theme_set(theme_gray(base_size = 18))
        sims = qplot(x,
            geom="histogram",
            breaks=seq(0, length(rownames(subset(res.diff, padj < 0.05))), by = signif(length(rownames(subset(res.diff, padj < 0.05)))/10,digits=1),fill=I("grey"), col=I("black"),
            main = "Histogram of DE experiments\n with random group labels", 
            xlab = "Number of differentially expressed genes",
            ylab = "Number of simulations") +
            geom_vline(xintercept = length(rownames(subset(res.diff, padj < 0.05)))) +
            theme_classic() + theme(plot.title = element_text(hjust = 0.5, size=22))
        print(sims)
        end_plot()
        flog.info(paste0("... Permutation p value: ",
                         length(x[x < length(rownames(subset(res.diff, padj < 0.05)))])/length(x)))
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
            "--counts",
            dest = "counts_file",
            type = "character",
            default = "",
            help = paste("path to file with paths to IRFinder outputs")
        ),
        make_option(
            "--design",
            dest = "design",
            type = "character",
            default = "",
            help = paste("path to designfile")
        ),
        make_option(
            "--model",
            dest = "model",
            type = "character",
            default = "~group+group:IRfinder",
            help = paste("model formula for GLM, must include IRFinder")
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
          "--compgroup",
          dest = "compgroup",
          type = "character",
          default = "MUT",
          help = paste("comparison group, e.g. MUT")
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
            type = "integer",
            default = 0,
            help = paste("number of permutations")
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
