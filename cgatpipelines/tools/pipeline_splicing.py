"""===========================
Pipeline Splicing
===========================


Overview
========

This pipeline enables differential exon usage testing through the
implementation of
* rMATS-turbo
* DEXSeq

rMATS is a computational tool to detect differential alternative splicing
events from RNA-Seq data. The statistical model of MATS calculates the P-value
and false discovery rate that the difference in the isoform ratio of a gene
between two conditions exceeds a given user-defined threshold. From the
RNA-Seq data, MATS can automatically detect and analyze alternative splicing
events corresponding to all major types of alternative splicing patterns.
MATS handles replicate RNA-Seq data from both paired and unpaired study design.

DEXSeq is a bioconductor R package to detect differential exon usage between
conditions in RNA-Seq experiments. Relative exon usage is defined as
(number of reads from exon)/(number of reads from its gene). It uses a similar
model to DESeq to estimate dispersion parameters prior to differential
testing.

Principal targets
-----------------

full
    compute all functions

Optional targets
----------------

permute
    repeat rMATS after permuting sample group labels


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use cgat pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.
cgatReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.yml` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_splicing.py config

Input files
-----------

".bam" files generated using STAR or Tophat2. Other mappers
may also work. Bam indexation is not required.

Design_files ("*.design.tsv") are used to specify sample variates. The
minimal design file is shown below, where include specifies if the
sample should be included in the analysis, group specifies the sample
group and pair specifies whether the sample is paired. Note, multiple
design files may be included, for example so that multiple models can
be fitted to different subsets of the data

(tab-seperated values)

sample    include    group    pair
WT-1-1    1    WT    0
WT-1-2    1    WT    0
Mutant-1-1    1    Mutant    0
Mutant-1-2    1    Mutant    0

The pipeline can only handle comparisons between two conditions with
replicates. If further comparisons are needed, further design files
should be used.

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

Requirements:

* samtools
* DEXSeq
* rMATS-turbo
* pysam
* HTSeqCounts

Pipeline output
===============

For each experiment, the output from rMATS is placed in the results.dir
folder. Each experiment is found in a subdirectory named designfilename.dir

rMATS output is further described here:
http://rnaseq-mats.sourceforge.net/user_guide.htm



Glossary
========

.. glossary::


Code
====

"""
from ruffus import *
import sys
import os
import glob
import shutil
import sqlite3
import pandas as pd
from rpy2.robjects import r as R
import cgat.BamTools.bamtools as BamTools
import cgatcore.iotools as iotools
import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatpipelines.tasks.tracks as tracks
import cgatpipelines.tasks.splicing as splicing
import cgatpipelines

###################################################################
###################################################################
###################################################################
# Load options and annotations
###################################################################

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# add configuration values from associated pipelines
PARAMS = P.PARAMS
PARAMS.update(P.peek_parameters(
    PARAMS["annotations_dir"],
    "genesets",
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))  # add config values from associated pipelines

# The DEXSeq R directory contains important python helper functions
PYTHONSCRIPTSDIR = R('''system.file("python_scripts", package="DEXSeq")''')[0]


###################################################################
###################################################################
###################################################################
# Utility functions
###################################################################

def connect():
    '''Connect to database (sqlite by default)

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


class MySample(tracks.Sample):
    attributes = tuple(PARAMS["attributes"].split(","))


TRACKS = tracks.Tracks(MySample).loadFromDirectory(
    glob.glob("*.bam"), "(\S+).bam")

Sample = tracks.AutoSample
DESIGNS = tracks.Tracks(Sample).loadFromDirectory(
    glob.glob("*.design.tsv"), "(\S+).design.tsv")


###################################################################
###################################################################
###################################################################
# DEXSeq workflow
###################################################################

@mkdir("results.dir")
@files(PARAMS["annotations_interface_geneset_all_gtf"],
       "geneset_flat.gff")
def buildGff(infile, outfile):
    '''Creates a gff for DEXSeq

    This takes the gtf and flattens it to an exon based input
    required by DEXSeq. The required python script is provided by DEXSeq
    and uses HTSeqCounts.

    Parameters
    ----------

    infile : string
       Input filename in :term:`gtf` format

    outfile : string
        A :term:`gff` file for use in DEXSeq

    annotations_interface_geneset_all_gtf : string
       :term:`PARAMS`. Filename of :term:`gtf` file containing
       all ensembl annotations
    '''

    tmpgff = P.get_temp_filename(".")
    statement = "gunzip -c %(infile)s > %(tmpgff)s"
    P.run(statement)

    ps = PYTHONSCRIPTSDIR
    statement = '''python %(ps)s/dexseq_prepare_annotation.py
                %(tmpgff)s %(outfile)s'''
    P.run(statement)

    os.unlink(tmpgff)


@mkdir("counts.dir")
@transform(glob.glob("*.bam"),
           regex("(\S+).bam"),
           add_inputs(buildGff),
           r"counts.dir/\1.txt")
def countDEXSeq(infiles, outfile):
    '''create counts for DEXSeq

    Counts bam reads agains exon features in flattened gtf.
    The required python script is provided by DEXSeq
    and uses HTSeqCounts.

    Parameters
    ----------

    infile[0]: string
        :term:`bam` file input

    infile[1]: string
        :term:`gff` output from buildGff function

    outfile : string
        A :term:`txt` file containing results

    DEXSeq_strandedness : string
       :term:`PARAMS`. Specifies strandedness, options
       are 'yes', 'no' and 'reverse'

    '''

    infile, gfffile = infiles
    ps = PYTHONSCRIPTSDIR
    if BamTools.is_paired(infile):
        paired = "yes"
    else:
        paired = "no"
    strandedness = PARAMS["DEXSeq_strandedness"]

    statement = '''python %(ps)s/dexseq_count.py
    -p %(paired)s
    -s %(strandedness)s
    -r pos
    -f bam  %(gfffile)s %(infile)s %(outfile)s'''
    P.run(statement)


@collate(countDEXSeq,
         regex("counts.dir/([^.]+)\.txt"),
         r"summarycounts.tsv")
def aggregateExonCounts(infiles, outfile):
    ''' Build a matrix of counts with exons and tracks dimensions.

    Uses `combine_tables.py` to combine all the `txt` files output from
    countDEXSeq into a single :term:`tsv` file named
    "summarycounts.tsv". A `.log` file is also produced.

    Parameters
    ---------
    infiles : list
        a list of `tsv.gz` files from the counts.dir that were the
        output from dexseq_count.py

    outfile : string
        a filename denoting the file containing a matrix of counts with genes
        as rows and tracks as the columns - this is a `tsv.gz` file      '''

    infiles = " ".join(infiles)
    statement = '''cgat combine_tables
    --columns=1
    --take=2
    --use-file-prefix
    --regex-filename='([^.]+)\.txt'
    --no-titles
    --log=%(outfile)s.log
    %(infiles)s
    > %(outfile)s '''

    P.run(statement)


@follows(aggregateExonCounts)
@mkdir("results.dir/DEXSeq")
@subdivide(["%s.design.tsv" % x.asFile().lower() for x in DESIGNS],
           regex("(\S+).design.tsv"),
           add_inputs(buildGff),
           r"results.dir/DEXSeq/\1/experiment_out.rds")
def filterDEXSeq(infiles, outfile):
    ''' Load counts into RDS object and filter'''

    design, gfffile = infiles
    countsdir = "counts.dir/"
    designname = design.split(".")[0]
    model = PARAMS["DEXSeq_model_%s" % designname]

    
    outdir = os.path.dirname(outfile)
    r_root = os.path.abspath(os.path.dirname(cgatpipelines.__file__))
    scriptpath = os.path.join(r_root, "Rtools/filtercounts.R")
    

    statement = '''
    export R_ROOT=%(r_root)s &&
    Rscript %(scriptpath)s
    --counts-dir %(countsdir)s
    --source dexseq
    --method dexseq
    --dexseq-flattened-file %(gfffile)s
    --sampleData %(design)s
    --outdir %(outdir)s
    --model %(model)s
    > %(outdir)s/filter.log;
    '''

    P.run(statement) 


@transform(filterDEXSeq,
           regex("(.+)\/(?P<DESIGN>\S+)\/experiment_out.rds"),
           r"\1/\2/results.tsv",
           r"\2")
def runDEXSeq(infile, outfile, design):
    ''' DEXSeq is run using the R scripts from the
    cgat code collection. Output is standardised to
    correspond to differential gene expression output
    from DESeq2 or Sleuth.

    Will currently only test 2 groups.

    Parameters
    ---------
    infiles : string
        filename and path of design file

    outfile : string
        a filename denoting the file containing a standard results
        output with full results of all tested exons

    DEXSeq_model_% : string
    DEXSeq_contrast_% : string
    DEXSeq_refgroup_% : string
       :term:`PARAMS`. Specifies model, contrast and reference
       group for DEXSeq analysis
    '''


    outdir = os.path.dirname(outfile)
    DEXSeq_fdr = 0.05
    DEXSeq_permutations = 0
    
    designname = design.split(".")[0]
    model = PARAMS["DEXSeq_model_%s" % designname]
    reducedmodel = PARAMS["DEXSeq_reducedmodel_%s" % designname]
    contrast = PARAMS["DEXSeq_contrast_%s" % designname]
    refgroup = PARAMS["DEXSeq_refgroup_%s" % designname]
    r_root = os.path.abspath(os.path.dirname(cgatpipelines.__file__))
    scriptpath = os.path.join(os.path.abspath(os.path.dirname(cgatpipelines.__file__)), "Rtools/diffexonexpression.R")

    statement = '''
    export R_ROOT=%(r_root)s &&
    Rscript %(scriptpath)s
    --rds-filename %(infile)s   
    --model %(model)s
    --reducedmodel %(reducedmodel)s
    --contrast %(contrast)s
    --refgroup %(refgroup)s
    --alpha %(DEXSeq_fdr)s
    --outdir %(outdir)s
    --permute %(DEXSeq_permutations)s

    > %(outdir)s/dexseq.log;
    '''

    P.run(statement)



###################################################################
###################################################################
###################################################################
# IRFinder workflow
###################################################################

@mkdir("IRFinder.dir")
@originate("IRFinder.dir/REF.log")
def buildIRReference(outfile):
    '''Creates a mapping reference for IRFinder using STAR

    This downloads a gtf from ensembl and the fitting FASTA and
    first runs STAR and then checks mapability

    Parameters
    ----------

    infile : string
       FTP path to correct ensembl annotation

    outfile : string
        One of the IRFinder reference files indicating
        completion of the run

    annotations_interface_geneset_all_gtf : string
       :term:`PARAMS`. Filename of :term:`gtf` file containing
       all ensembl annotations
    '''

    extra = PARAMS['IRFinder_extra']
    bedfile = PARAMS['IRFinder_bed']
    gtf = PARAMS["annotations_interface_geneset_all_gtf"]
    star = PARAMS['IRFinder_ensembl_star']
    genome = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))
    irfinder = PARAMS["IRFinder_singularity"]

    job_threads = PARAMS["IRFinder_threads"]
    job_memory = PARAMS["IRFinder_memory"]

    tmpgtf = P.get_temp_filename(".")
    statement = "gunzip -c %(gtf)s > %(tmpgtf)s"
    P.run(statement)

    statement = '''singularity run -H $PWD:/home
                   -B %(star)s,%(gtf)s,%(genome)s'''
    if extra is not None:
        statement += ",%(extra)s"
    if bedfile is not None:
        statement += ",%(bedfile)s "
    statement += '''%(irfinder)s BuildRefFromSTARRef 
                    -r IRFinder.dir/REF
                    -x %(star)s
                    -g %(tmpgtf)s
                    -f %(genome)s
                    -t %(job_threads)s '''
    if extra is not None:
        statement += "-e %(extra)s "
    if bedfile is not None:
        statement += "-b %(bedfile)s "

    statement +=  " > IRFinder.dir/REF.log"

    P.run(statement)

    os.unlink(tmpgtf)

@transform(glob.glob("*.bam"),
           regex("(\S+).bam"),
           add_inputs(buildIRReference),           
           r"IRFinder.dir/\1/IRFinder-IR-nondir.txt")
def runIRFinder(infiles, outfile):
    '''
    Maps reads using STAR and creates Intron quantification tables.

    Parameters
    ----------
    infile: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    IRFinder_singularity: str
        :term:`PARAMS`
        path to IRFinder singularity executable

    outfile: str
        :term:`txt` filename to write .
    '''

    infile, reference = infiles
    ref_dir = P.snip(reference)
    singularity_dir = os.path.dirname(os.getcwd())


    job_threads = PARAMS["IRFinder_threads"]
    job_memory = PARAMS["IRFinder_memory"]
    outdir = os.path.dirname(outfile)
    executable = PARAMS["IRFinder_singularity"]

    statement = '''
    singularity run -H $PWD:/home
    -B %(singularity_dir)s
    %(executable)s BAM
    -r %(ref_dir)s 
    -d %(outdir)s
    -t %(IRFinder_threads)s
    %(infile)s;'''

    P.run(statement)


@collate(runIRFinder,
         regex("IRFinder.dir/(.*)/IRFinder-IR-nondir.txt"),
         r"IRFinder.dir/filelist.tsv")
def aggregateIRFinder(infiles, outfile):
    ''' Build a matrix of counts with exons and tracks dimensions.

    Uses `combine_tables.py` to combine all the `txt` files output from
    countDEXSeq into a single :term:`tsv` file named
    "summarycounts.tsv". A `.log` file is also produced.

    Parameters
    ---------
    infiles : list
        a list of `tsv.gz` files from the counts.dir that were the
        output from dexseq_count.py

    outfile : string
        a filename denoting the file containing a matrix of counts with genes
        as rows and tracks as the columns - this is a `tsv.gz` file      '''

    # if stranded/directional output exists use that instead
    infiles_dir = [infile.replace("nondir","dir") for infile in infiles]
    if os.path.exists(infiles_dir[1]):
        infiles = infiles_dir

    with open(outfile, 'a') as outputfile:
        for infile in infiles:
            outputfile.write(os.path.abspath(infile+"\n"))


@mkdir("results.dir/IRFinder")
@subdivide(["%s.design.tsv" % x.asFile().lower() for x in DESIGNS],
           regex("(\S+).design.tsv"),
           add_inputs(aggregateIRFinder),
           r"results.dir/IRFinder/\1/results.tsv")
def diffIRFinder(infiles, outfile):

    design, counts = infiles
    designname = design.split(".")[0]
    model = PARAMS["IRFinder_model_%s" % designname]
    contrast = PARAMS["IRFinder_contrast_%s" % designname]
    REF = PARAMS["IRFinder_reference_%s" % designname]
    COMP = PARAMS['IRFinder_comparator_%s' % designname]
    permute = PARAMS["IRFinder_permute"]
    pvalue = PARAMS["IRFinder_pvalue"]
    r_root = os.path.abspath(os.path.dirname(cgatpipelines.__file__))
    scriptpath = os.path.join(r_root, "Rtools/IRFinder.R")

    
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    statement = '''
    export R_ROOT=%(r_root)s &&
    Rscript %(scriptpath)s
    --counts %(counts)s
    --design %(design)s
    --model %(model)s
    --contrast %(contrast)s
    --refgroup %(REF)s
    --compgroup %(COMP)s
    --permute %(permute)s
    --alpha %(pvalue)s
    --outdir %(outdir)s > %(outdir)s/IRFinder.log
    '''

    P.run(statement)
    return



###################################################################
###################################################################
###################################################################
# rMATS workflow
###################################################################

@mkdir("results.dir/rMATS")
@subdivide(["%s.design.tsv" % x.asFile().lower() for x in DESIGNS],
           regex("(\S+).design.tsv"),
           add_inputs(PARAMS["annotations_interface_geneset_all_gtf"]),
           [r"results.dir/rMATS/\1.dir/%s.MATS.JC.txt" % x for x in ["SE", "A5SS", "A3SS", "MXE", "RI"]])
def runMATS(infile, outfiles):
    '''run rMATS-turbo

    Runs rMATS command.

    Parameters
    ---------
    infiles[0] : string
        filename and path of design file

    infiles[1] : string
        filename and path of :term:`gtf` file

    outfile : list
        a list of filenames denoting the file containing a standard results
        output with full results for all five tested differential exon
        usage conditions.

    MATS_libtype : string
       :term:`PARAMS`. Specifies library type. Can be "fr-firstrand",
       "fr-secondstrand" or "fr-unstranded"
    '''

    design, gtffile = infile
    strand = PARAMS["MATS_libtype"]
    cutoff = PARAMS["MATS_cutoff"]
    outdir = os.path.dirname(outfiles[0])
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    splicing.runRMATS(gtffile=gtffile, designfile=design,
                              pvalue=cutoff,
                              strand=strand, outdir=outdir)


@follows(runMATS)
@transform(runMATS,
           regex("results.dir/rMATS/(\S+).dir/(\S+).MATS.JC.txt"),
           r"results.dir/rMATS/rMATS_\1_.dir/\2_JC.load")
def loadMATS(infile, outfile):
    '''load RMATS results into relational database

    Loads rMATS results into relational database.
    Continues if table empty.

    Parameters
    ----------
    infile: term:`tsv` file containing one type of rMATS results.
    outfile: .load file
    '''
    
    try:
        P.load(infile, outfile)
    except:
        iotools.touch_file(outfile)


@follows(runMATS)
@transform(runMATS,
           regex("results.dir/rMATS/(\S+).dir/(\S+).MATS.JC.txt"),
           r"results.dir/rMATS/\1.dir/\2.novelJunc.JC.tsv")
def novelJuncMATS(infile, outfile):
    '''Combine novel Junction result with MATS output.

    Parameters
    ----------
    infile: term:`tsv` file containing one type of rMATS results.
    outfile: .load file
    '''

    indir = os.path.dirname(infile)
    event = os.path.basename(infile).split(".")[0]
    
    statement = "cut -f 4- %(indir)s/fromGTF.novelJunction.%(event)s.txt |grep -Fwf - %(infile)s > %(outfile)s"
    P.run(statement)

    

@follows(novelJuncMATS)
@collate(runMATS,
         regex("results.dir/rMATS/(\S+).dir/\S+.MATS.JC.txt"),
         r"results.dir/rMATS/rMATS_\1_results.summary")
def collateMATS(infiles, outfile):
    '''collates summary results from all events

    Collates number of events below FDR threshold from all
    five events into simple table, split for up and 
    downregulated events

    Parameters
    ----------
    infiles: list
        list of results files from rMATS

    MATS_fdr : string
       :term:`PARAMS`. User specified threshold for result counting

    outfile: string
        summary file containing number of results below FDR threshold
    '''

    indir = os.path.dirname(infiles[1])
    design = P.snip(os.path.basename(os.path.normpath(indir)))


    total = [design, "all", "Total"]
    up = [design, "all", "Sample1HigherInclusion"]
    down = [design, "all", "Sample2HigherInclusion"]
    total_newJunc = [design, "novelJunc", "Total"]
    up_newJunc = [design, "novelJunc", "Sample1HigherInclusion"]
    down_newJunc = [design, "novelJunc", "Sample2HigherInclusion"]
    #total_newSS = [design, "novelSS", "Total"]
    #up_newSS = [design, "novelSS", "Sample1HigherInclusion"]
    #down_newSS = [design, "novelSS", "Sample2HigherInclusion"]


    for event in ["SE", "A5SS", "A3SS", "MXE", "RI"]:
        temp = pd.read_csv("%s/%s.MATS.JC.txt" %
                           (indir, event), sep='\t')
        total.append(int(len(temp[(temp['FDR'] <
                                float(PARAMS['MATS_fdr'])) & (abs(temp['IncLevelDifference']) > 0.1)])))
        up.append(int(len(temp[(temp['FDR'] <
                                float(PARAMS['MATS_fdr'])) & (temp['IncLevelDifference'] > 0.1)])))
        down.append(int(len(temp[(temp['FDR'] <
                                  float(PARAMS['MATS_fdr'])) & (temp['IncLevelDifference'] < -0.1)])))
        temp = pd.read_csv("%s/%s.novelJunc.JC.tsv" %
                           (indir, event), sep='\t')
        total_newJunc.append(int(len(temp[(temp['FDR'] <
                                float(PARAMS['MATS_fdr'])) & (abs(temp['IncLevelDifference']) > 0.1)])))
        up_newJunc.append(int(len(temp[(temp['FDR'] <
                                float(PARAMS['MATS_fdr'])) & (temp['IncLevelDifference'] > 0.1)])))
        down_newJunc.append(int(len(temp[(temp['FDR'] <
                                  float(PARAMS['MATS_fdr'])) & (temp['IncLevelDifference'] < -0.1)])))
        #experimental feature - deactivated       
        #temp = pd.read_csv("%s/%s.novelSS.JC.tsv" %
        #                   (indir, event), sep='\t')
        #total_newSS.append(int(len(temp[(temp['FDR'] <
        #                        float(PARAMS['MATS_fdr'])) & (abs(temp['IncLevelDifference']) > 0.1)])))
        #up_newSS.append(int(len(temp[(temp['FDR'] <
        #                        float(PARAMS['MATS_fdr'])) & (temp['IncLevelDifference'] > 0.1)])))
        #down_newSS.append(int(len(temp[(temp['FDR'] <
        #                          float(PARAMS['MATS_fdr'])) & (temp['IncLevelDifference'] < -0.1)])))

    eventdf = pd.DataFrame([total,up,down,total_newJunc,up_newJunc,down_newJunc], columns=["Design","Type","Subset","SE","A5SS","A3SS","MXE","RI"])
    eventdf.loc[:,"TOTAL"] = eventdf.sum(numeric_only=True, axis=1)

    eventdf.to_csv(outfile, sep="\t", index=False)


@transform(collateMATS,
           suffix(".summary"),
           ".load")
def loadCollateMATS(infile, outfile):
    '''load rMATS summary into relational database

    Loads rMATS summary results into relational database.

    Parameters
    ----------
    infile: file containing summary table of rMATS results
    outfile: .load file
    '''

    P.load(infile, outfile)


@active_if(PARAMS["permute"] == 1)
@subdivide(["%s.design.tsv" % x.asFile().lower() for x in DESIGNS],
           regex("(\S+).design.tsv"),
           r"results.dir/rMATS/\1.dir/permutations/run*.dir/init",
           r"results.dir/rMATS/\1.dir/permutations")
def permuteMATS(infile, outfiles, outdir):
    '''creates directories for permutation testing

    Creates directories for permutation testing and leaves dummy
    init file in directory (for timestamping)
    Only becomes active if :term:`PARAMS` permute is set to 1

    Parameters
    ----------
    infile: string
        name and path to design

    outfile: list
        list of unknown length, capturing all permutations
        retrospectively

    outdir: string
        directory to generate permutations in

    permutations : string
       :term:`PARAMS`. number of directories to be generated
    '''

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for i in range(0, PARAMS["permutations"]):
        if not os.path.exists("%s/run%i.dir" % (outdir, i)):
            os.makedirs("%s/run%i.dir" % (outdir, i))
        iotools.touch_file("%s/run%i.dir/init" % (outdir, i))


@transform(permuteMATS,
           regex("results.dir/rMATS/(\S+).dir/permutations/(\S+).dir/init"),
           add_inputs(PARAMS["annotations_interface_geneset_all_gtf"]),
           r"results.dir/rMATS/\1.dir/permutations/\2.dir/result.tsv",
           r"\1.design.tsv")
def runPermuteMATS(infiles, outfile, design):
    '''run rMATS-turbo permutation testing

    Runs rMATS command on permutations and then collates results into
    small summary table for each permutation
 
    Parameters
    ---------
    infiles[0] : string
        filename and path of design file

    infiles[1] : string
        filename and path of :term:`gtf` file

    outfile : :term:`tsv` file
        file containing summary results meeting the user-specified FDR
        threshold

    design : string
        name and path of design file

    MATS_libtype : string
       :term:`PARAMS`. Specifies library type. Can be "fr-firstrand",
       "fr-secondstrand" or "fr-unstranded"

    MATS_fdr : string
       :term:`PARAMS`. User specified threshold for result counting

    '''

    init, gtffile = infiles
    directory = os.path.dirname(init)
    strand = PARAMS["MATS_libtype"]
    cutoff = PARAMS["MATS_cutoff"]
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    splicing.runRMATS(gtffile=gtffile, designfile=design,
                              pvalue=cutoff,
                              strand=strand, outdir=directory, permute=1)

    collate = []
    with open(os.path.dirname(init) + "/b1.txt", "r") as f:
        collate.append(f.readline())
    with open(os.path.dirname(init) + "/b2.txt", "r") as f:
        collate.append(f.readline())
    for event in ["SE", "A5SS", "A3SS", "MXE", "RI"]:
        temp = pd.read_csv("%s/%s.MATS.JC.txt" %
                           (os.path.dirname(outfile), event), sep='\t')
        collate.append(str(len(temp[(temp['FDR'] <
                                float(PARAMS['MATS_fdr'])) & (abs(temp['IncLevelDifference']) > 0.1)])))
    with open(outfile, "w") as f:
        f.write("Group1\tGroup2\tSE\tA5SS\tA3SS\tMXE\tRI\n")
        f.write('\t'.join(collate))


@collate(runPermuteMATS,
         regex("results.dir/rMATS/(\S+).dir/permutations/\S+.dir/result.tsv"),
         r"results.dir/rMATS/rMATS_\1_permutations.summary")
def collatePermuteMATS(infiles, outfile):
    '''collates summary table of all permutations

    Collates number of events below FDR threshold from all
    permutation runs.

    Parameters
    ----------
    infiles: list
        list of rMATS result summaries from all permutation runs

    outfile: string
        summary file containing a table with all permutation run
        results
    '''

    collate = []
    for infile in infiles:
        collate.append(pd.read_csv(infile, sep='\t'))
    pd.concat(collate).to_csv(outfile, sep='\t', index=0)


@transform(collatePermuteMATS,
           suffix(".summary"),
           ".load")
def loadPermuteMATS(infile, outfile):
    '''load rMATS permutation results

    Loads rMATS permutation summary results into relational database.

    Parameters
    ----------
    infile: file containing summary table of rMATS permutation results
    outfile: .load file
    '''

    P.load(infile, outfile)


@mkdir("results.dir/sashimi")
@transform(runMATS,
           regex("results.dir/rMATS/(\S+).dir/(\S+).MATS.JC.txt"),
           add_inputs(r"\1.design.tsv"),
           r"results.dir/sashimi/\1.dir/\2")
def runSashimi(infiles, outfile):
    '''draws sashimi plots

    Draws Sashimi plots (pdf files) for all results below FDR threshold
    from all five rMATS events


    Parameters
    ----------
    infiles: list
        list of results files from rMATS

    MATS_fdr : string
       :term:`PARAMS`. User specified threshold for result drawing

    outfile: string
        summary file containing number of results below FDR threshold
    '''

    infile, design = infiles
    fdr = PARAMS["MATS_fdr"]
    plotmax = int(PARAMS["MATS_plotmax"])
    if not os.path.exists(outfile):
        os.makedirs(outfile)

    splicing.rmats2sashimi(infile, design, fdr, outfile, plotmax)


###################################################################
###################################################################
###################################################################
# Pipeline management
###################################################################

@follows(loadMATS,
         loadCollateMATS,
         loadPermuteMATS,
         runSashimi,
         runDEXSeq,
         diffIRFinder)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
