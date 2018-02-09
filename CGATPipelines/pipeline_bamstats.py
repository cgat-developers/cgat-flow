"""
===========================
Pipeline bamstats
===========================

:Author: Adam Cribbs
:Release: $Id$
:Date: |today|
:Tags: Python

The intention of this pipeline is to perform QC statistics on
   `.bam` files that are produced following mapping of fastq
   files.

The pipeline requires a `.bam` file as an input.

Overview
========

The pipeline perform the following stats in each folder:
    * IdxStats     -Samtools idxstats is ran and this calculates
                    the number of mapped and unmapped reads per contig.
    * BamStats     -This is a CGAT script (bam2stats) that performs stats
                    on a bam file and outputs alignment statistics.
    * PicardStats  -this runs to CollectRnaSeqMetrics picard tools.
    * StrandSpec   -Gives a measure of the proportion of reads that map to
                    each strand. Is used to work out strandness of library
                    if unknown.
    * nreads       -Calculates the number of reads in the bam file.
    * Paired_QC    -This contains metrics that are only required for paired
                    end. Most of the statistics collate metrics regarding
                    splicing.
                    Transcript profile is across the upstream,exons and
                    downstream because this is usually specific to rna seq
                    analysis. ### May need to remove this to make it single ended...........



This pipeline computes the word frequencies in the configuration
files :file:``pipeline.ini` and :file:`conf.py`.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.


Configuration
-------------

This pipeline requires the user to run pipeline_gtf_subset.py. The
location of the database then needs to be set in the pipeline.ini
file.

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_bamstats.py config

Input files
-----------

The pipeline configuration files need to be generated by running:

   python <srcdir>/pipeline_bamstats.py config

Once the config file  (pipeline.ini) is generated this should be modified
before running the pipeline.


The pipeline requires `.bam` files to be loacted within the same directory
that the piepline is ran.

Requirements
------------

The pipeline requires the gtf file produced from
:doc:`pipeline_gtf_subset`. Set the configuration variable
:py:data:`gtf_database`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:


+--------------+----------+------------------------------------+
|*Program*     |*Version* |*Purpose*                           |
+--------------+----------+------------------------------------+
|samtools      |>=0.1.16  |bam/sam files                       |
+--------------+----------+------------------------------------+
|cgat tools    |          |bam2stats script                    |
+--------------+----------+------------------------------------+
|picard        |>=1.42    |bam/sam files. The .jar files need  |
|              |          |to be in your CLASSPATH environment |
|              |          |variable.                           |
+--------------+----------+------------------------------------+
|bamstats_     |>=1.22    |from CGR, Liverpool                 |
+--------------+----------+------------------------------------+




Pipeline output
===============

The major output of the pipeline is the database file :file:`csvdb`.

SQL query of this database forms the basis of the final reports.

The following reports are generated as part of running:

    python <srcdir>/pipeline_bamstats.py make build_report

    * Jupyter notebook - a python implimentation. The output files
                         are located in Jupyter_report.dir. To view
                         the report open the _site/CGAT_FULL_BAM_STATS_REPORT.html.
                         You can navigate throught the various report
                         pages through here.

    * Rmarkdown        - an R markdown report implimentation.The output
                         report os located in the R_report.dir/_site
                         directory and can be accessed by opening any of
                         the html files.

    * multiQC          - this builds a basic report using the multiqc -
                         http://multiqc.info/ external tool. There is the
                         potential for customising multiQC so it can be used
                         to generate reports from CGAT tools, however at presnt this
                         is not possible because of development stage of multiQC.

Example
=======

Example data is available at:
..........Add data...............

python <srcdir>/pipeline_bamstats.py config
python <srcdir>/pipeline_bamstats.py make full


Glossary
========

.. glossary::

.. _bamstats: http://www.agf.liv.ac.uk/454/sabkea/samStats_13-01-2011


Code
====

"""

# load modules for use in the pipeline

import sys
import os
import sqlite3
import CGATCore.IOTools as IOTools

from ruffus import *


import CGATCore.Pipeline as P
import CGATPipelines.PipelineBamStats as PipelineBamStats


# load options from the config file
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS

# Add parameters from the gtf_subset pipeline, but
# only the interface section. All PARAMS options
# will have the prefix `annotations_`

PARAMS.update(P.peekParameters(
    PARAMS["gtf_dir"],
    "pipeline_genesets.py",
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))

# -----------------------------------------------
# Utility functions


def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])

    if not os.path.exists(PARAMS["gtf_database"]):
        raise ValueError(
            "can't find database '%s'" %
            PARAMS["gtf_database"])

    statement = '''ATTACH DATABASE '%s' as annotations''' % \
                (PARAMS["gtf_database"])

    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

# Determine whether the gemone is paired

SPLICED_MAPPING = PARAMS["bam_paired_end"]


#########################################################################
# Count reads as some QC targets require it
#########################################################################


@follows(mkdir("nreads.dir"))
@transform("*.bam",
           suffix(".bam"),
           r"nreads.dir/\1.nreads")
def countReads(infile, outfile):
    '''Count number of reads in input files.'''

    statement = '''printf "nreads \\t" >> %(outfile)s'''

    P.run()

    statement = '''samtools view %(infile)s | wc -l | xargs printf >> %(outfile)s'''

    P.run()

#########################################################################
# QC tasks start here
#########################################################################


@follows(mkdir("StrandSpec.dir"))
@transform("*.bam",
           suffix(".bam"),
           r"StrandSpec.dir/\1.strand")
def strandSpecificity(infile, outfile):
    '''This function will determine the strand specificity of your library
    from the bam file'''

    iterations = "1000000"

    PipelineBamStats.getStrandSpecificity(infile,
                                          outfile,
                                          iterations)


@follows(mkdir("BamFiles.dir"))
@transform("*.bam",
           regex("(.*).bam$"),
           r"BamFiles.dir/\1.bam")
def intBam(infile, outfile):
    '''make an intermediate bam file if there is no sequence infomation.
    If there is no sequence quality then make a softlink. Picard tools
    has an issue when quality score infomation is missing'''

    if PARAMS["bam_sequence_stipped"] is True:
        PipelineBamStats.addPseudoSequenceQuality(infile,
                                                  outfile)
    else:
        PipelineBamStats.copyBamFile(infile,
                                     outfile)


@follows(mkdir("Picard_stats.dir"))
@P.add_doc(PipelineBamStats.buildPicardAlignmentStats)
@transform(intBam,
           regex("BamFiles.dir/(.*).bam$"),
           add_inputs(os.path.join(PARAMS["genome_dir"],
                                   PARAMS["genome"] + ".fa")),
           r"Picard_stats.dir/\1.picard_stats")
def buildPicardStats(infiles, outfile):
    ''' build Picard alignment stats '''
    infile, reffile = infiles

    # patch for mapping against transcriptome - switch genomic reference
    # to transcriptomic sequences
    if "transcriptome.dir" in infile:
        reffile = "refcoding.fa"

    PipelineBamStats.buildPicardAlignmentStats(infile,
                                               outfile,
                                               reffile)


@P.add_doc(PipelineBamStats.buildPicardDuplicationStats)
@transform(intBam,
           regex("BamFiles.dir/(.*).bam$"),
           r"Picard_stats.dir/\1.picard_duplication_metrics")
def buildPicardDuplicationStats(infile, outfile):
    '''Get duplicate stats from picard MarkDuplicates '''
    PipelineBamStats.buildPicardDuplicationStats(infile, outfile)


@follows(mkdir("BamStats.dir"))
@follows(countReads)
@transform(intBam,
           regex("BamFiles.dir/(.*).bam$"),
           add_inputs(r"nreads.dir/\1.nreads"),
           r"BamStats.dir/\1.readstats")
def buildBAMStats(infiles, outfile):
    '''count number of reads mapped, duplicates, etc.

    Excludes regions overlapping repetitive RNA sequences

    Parameters
    ----------
    infiles : list
    infiles[0] : str
       Input filename in :term:`bam` format
    infiles[1] : str
       Input filename with number of reads per sample

    outfile : str
       Output filename with read stats

    annotations_interface_rna_gtf : str
        :term:`PARMS`. :term:`gtf` format file with repetitive rna
    '''

    rna_file = PARAMS["annotations_interface_rna_gff"]

    job_memory = "32G"

    bamfile, readsfile = infiles

    nreads = PipelineBamStats.getNumReadsFromReadsFile(readsfile)
    track = P.snip(os.path.basename(readsfile),
                   ".nreads")

    # if a fastq file exists, submit for counting
    if os.path.exists(track + ".fastq.gz"):
        fastqfile = track + ".fastq.gz"
    elif os.path.exists(track + ".fastq.1.gz"):
        fastqfile = track + ".fastq.1.gz"
    else:
        fastqfile = None

    if fastqfile is not None:
        fastq_option = "--fastq-file=%s" % fastqfile
    else:
        fastq_option = ""

    statement = '''
    cgat bam2stats
         %(fastq_option)s
         --force-output
         --mask-bed-file=%(rna_file)s
         --ignore-masked-reads
         --num-reads=%(nreads)i
         --output-filename-pattern=%(outfile)s.%%s
    < %(bamfile)s
    > %(outfile)s
    '''

    P.run()


@follows(intBam)
@transform(PARAMS["annotations_interface_genomic_context_bed"],
           regex("^\/(.+\/)*(.+).bed.gz"),
           r"BamStats.dir/\2.bed.gz")
def processGenomicContext(infile, outfile):
    '''
    This module process genomic context file.
    It assigns each and every features of context
    file to a specific catagory. It helps us to
    understand heiarchical classification
    of features.
    '''
    PipelineBamStats.defineBedFeatures(infile, outfile)


@follows(processGenomicContext)
@P.add_doc(PipelineBamStats.summarizeTagsWithinContext)
@transform(intBam,
           regex("BamFiles.dir/(.*).bam$"),
           add_inputs(processGenomicContext),
           r"BamStats.dir/\1.contextstats.tsv.gz")
def buildContextStats(infiles, outfile):
    ''' build mapping context stats '''
    PipelineBamStats.summarizeTagsWithinContext(
        infiles[0], infiles[1], outfile)


@follows(mkdir("IdxStats.dir"))
@transform(intBam,
           regex("BamFiles.dir/(.*).bam$"),
           r"IdxStats.dir/\1.idxstats")
def buildIdxStats(infile, outfile):
    '''gets idxstats for bam file so number of reads per chromosome can
    be plotted later'''

    statement = '''samtools idxstats %(infile)s > %(outfile)s'''

    P.run()

# ------------------------------------------------------------------
# QC specific to spliced mapping
# ------------------------------------------------------------------


@follows(mkdir("Paired_QC.dir"))
@active_if(SPLICED_MAPPING)
@transform(intBam,
           regex("BamFiles.dir/(.*).bam$"),
           add_inputs(PARAMS["annotations_interface_geneset_coding_exons_gtf"]),
           r"Paired_QC.dir/\1.exon.validation.tsv.gz")
def buildExonValidation(infiles, outfile):
    '''Compare the alignments to the exon models to quantify exon
    overrun/underrun

    Expectation is that reads should not extend beyond known exons.

    Parameters
    ----------
    infiles : list
    infiles[0] : str
       Input filename in :term:`bam` format
    infiles[1] : str
       Input filename in :term:`gtf` format

    outfile : str
       Output filename in :term:`gtf` format with exon validation stats
    '''

    infile, exons = infiles
    statement = '''cat %(infile)s
    | cgat bam_vs_gtf
         --exons-file=%(exons)s
         --force-output
         --log=%(outfile)s.log
         --output-filename-pattern="%(outfile)s.%%s.gz"
    | gzip
    > %(outfile)s
    '''

    P.run()


@active_if(SPLICED_MAPPING)
@transform(intBam,
           regex("BamFiles.dir/(.*).bam$"),
           add_inputs(PARAMS["annotations_interface_geneset_coding_exons_gtf"]),
           r"Paired_QC.dir/\1.transcript_counts.tsv.gz")
def buildTranscriptLevelReadCounts(infiles, outfile):
    '''count reads in gene models

    Count the reads from a :term:`bam` file which overlap the
    positions of protein coding transcripts in a :term:`gtf` format
    transcripts file.

    Parameters
    ----------
    infiles : list of str
    infiles[0] : str
       Input filename in :term:`bam` format
    infiles[1] : str
       Input filename in :term:`gtf` format

    outfile : str
       Output filename in :term:`tsv` format


    .. note::
       In paired-end data sets each mate will be counted. Thus
       the actual read counts are approximately twice the fragment
       counts.

    '''
    infile, geneset = infiles

    job_memory = "8G"

    statement = '''
    zcat %(geneset)s
    | cgat gtf2table
    --reporter=transcripts
    --bam-file=%(infile)s
    --counter=length
    --column-prefix="exons_"
    --counter=read-counts
    --column-prefix=""
    --counter=read-coverage
    --column-prefix=coverage_
    -v 0
    | gzip
    > %(outfile)s
    ''' % locals()

    P.run()


@active_if(SPLICED_MAPPING)
@transform(intBam,
           regex("BamFiles.dir/(.*).bam$"),
           add_inputs(PARAMS["annotations_interface_geneset_intron_gtf"]),
           r"Paired_QC.dir/\1.intron_counts.tsv.gz")
def buildIntronLevelReadCounts(infiles, outfile):
    '''count reads in gene models
    Count the reads from a :term:`bam` file which overlap the
    positions of introns in a :term:`gtf` format transcripts file.
    Parameters
    ----------
    infiles : list of str
       infile :term:`str`
          Input filename in :term:`bam` format
       geneset :term:`str`
          Input filename in :term:`gtf` format
    outfile : str
       Output filename in :term:`tsv` format
    .. note::
       In paired-end data sets each mate will be counted. Thus
       the actual read counts are approximately twice the fragment
       counts.
    '''
    infile, exons = infiles

    job_memory = "4G"

    if "transcriptome.dir" in infile:
        P.touch(outfile)
        return

    statement = '''
    zcat %(exons)s
    | awk -v OFS="\\t" -v FS="\\t" '{$3="exon"; print}'
    | cgat gtf2table
          --reporter=genes
          --bam-file=%(infile)s
          --counter=length
          --column-prefix="introns_"
          --counter=read-counts
          --column-prefix=""
          --counter=read-coverage
          --column-prefix=coverage_
    | gzip
    > %(outfile)s
    '''

    P.run()


@active_if(SPLICED_MAPPING)
@transform(intBam,
           regex("BamFiles.dir/(.*).bam$"),
           add_inputs(PARAMS["annotations_interface_geneset_coding_exons_gtf"]),
           r"Paired_QC.dir/\1.transcriptprofile.gz")
def buildTranscriptProfiles(infiles, outfile):
    '''build gene coverage profiles

    PolyA-RNA-Seq is expected to show a bias towards the 3' end of
    transcripts. Here we generate a meta-profile for each sample for
    the read depth from the :term:`bam` file across the gene models
    defined in the :term:`gtf` gene set

    In addition to the outfile specified by the task, plots will be
    saved with full and focus views of the meta-profile

    Parameters
    ----------
    infiles : list of str
    infiles[0] : str
       Input filename in :term:`bam` format
    infiles[1] : str`
       Input filename in :term:`gtf` format

    outfile : str
       Output filename in :term:`tsv` format
    '''

    bamfile, gtffile = infiles

    job_memory = "8G"

    statement = '''cgat bam2geneprofile
    --output-filename-pattern="%(outfile)s.%%s"
    --force-output
    --reporter=transcript
    --use-base-accuracy
    --method=geneprofileabsolutedistancefromthreeprimeend
    --normalize-profile=all
    %(bamfile)s %(gtffile)s
    | gzip
    > %(outfile)s
    '''

    P.run()


@active_if(SPLICED_MAPPING)
@P.add_doc(PipelineBamStats.buildPicardRnaSeqMetrics)
@transform(intBam,
           regex("BamFiles.dir/(.*).bam$"),
           add_inputs(PARAMS["annotations_interface_ref_flat"]),
           r"Picard_stats.dir/\1.picard_rna_metrics")
def buildPicardRnaSeqMetrics(infiles, outfile):
    '''Get duplicate stats from picard RNASeqMetrics '''
    # convert strandness to tophat-style library type
    if PARAMS["strandness"] == ("RF" or "R"):
        strand = "SECOND_READ_TRANSCRIPTION_STRAND"
    elif PARAMS["strandness"] == ("FR" or "F"):
        strand = "FIRST_READ_TRANSCRIPTION_STRAND"
    else:
        strand = "NONE"
    PipelineBamStats.buildPicardRnaSeqMetrics(infiles, strand, outfile)


##########################################################################
# Database loading statements
##########################################################################


@P.add_doc(PipelineBamStats.loadPicardAlignmentStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildPicardStats, "Picard_stats.dir/picard_stats.load")
def loadPicardStats(infiles, outfile):
    '''merge alignment stats into single tables.'''
    PipelineBamStats.loadPicardAlignmentStats(infiles, outfile)


@P.add_doc(PipelineBamStats.loadPicardDuplicationStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildPicardDuplicationStats, ["picard_duplication_stats.load",
                                     "picard_duplication_histogram.load"])
def loadPicardDuplicationStats(infiles, outfiles):
    '''merge alignment stats into single tables.'''

    PipelineBamStats.loadPicardDuplicationStats(infiles, outfiles)


@P.add_doc(PipelineBamStats.loadBAMStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildBAMStats, "bam_stats.load")
def loadBAMStats(infiles, outfile):
    ''' load bam statistics into bam_stats table '''
    PipelineBamStats.loadBAMStats(infiles, outfile)


@P.add_doc(PipelineBamStats.loadSummarizedContextStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@follows(loadBAMStats)
@merge(buildContextStats, "context_stats.load")
def loadContextStats(infiles, outfile):
    ''' load context mapping statistics into context_stats table '''
    PipelineBamStats.loadSummarizedContextStats(infiles, outfile)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildIdxStats, "idxstats_reads_per_chromosome.load")
def loadIdxStats(infiles, outfile):
    '''merge idxstats files into single dataframe and load
    to database

    Loads tables into the database
       * mapped_reads_per_chromosome

    Arguments
    ---------
    infiles : list
        list where each element is a string of the filename containing samtools
        idxstats output. Filename format is expected to be 'sample.idxstats'
    outfile : string
        Logfile. The table name will be derived from `outfile`.'''

    PipelineBamStats.loadIdxstats(infiles, outfile)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@active_if(SPLICED_MAPPING)
@merge(buildExonValidation, "exon_validation.load")
def loadExonValidation(infiles, outfile):
    ''' load individual and merged exon validation stats

    For each sample, the exon validation stats are loaded into a table
    named by sample and mapper
    [sample]_[mapper]_overrun

    The merge alignment stats for all samples are merged and loaded
    into single table called exon_validation

    Parameters
    ----------
    infiles : list
       Input filenames with exon validation stats
    outfile : str
       Output filename
    '''

    suffix = ".exon.validation.tsv.gz"

    P.mergeAndLoad(infiles, outfile, suffix=suffix)
    for infile in infiles:
        track = P.snip(infile, suffix)
        o = "%s_overrun.load" % track
        P.load(infile + ".overrun.gz", o)


@P.add_doc(PipelineBamStats.loadPicardRnaSeqMetrics)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildPicardRnaSeqMetrics, ["picard_rna_metrics.load",
                                  "picard_rna_histogram.load"])
def loadPicardRnaSeqMetrics(infiles, outfiles):
    '''merge alignment stats into single tables.'''
    PipelineBamStats.loadPicardRnaSeqMetrics(infiles, outfiles)


@P.add_doc(PipelineBamStats.loadCountReads)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@follows(loadPicardRnaSeqMetrics)
@merge(countReads, "count_reads.load")
def loadCountReads(infiles, outfile):
    ''' load read counts count_reads table '''
    PipelineBamStats.loadCountReads(infiles, outfile)


@active_if(SPLICED_MAPPING)
@P.add_doc(PipelineBamStats.loadTranscriptProfile)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@follows(loadCountReads)
@merge(buildTranscriptProfiles, "transcript_profile.load")
def loadTranscriptProfile(infiles, outfile):
    ''' merge transcript profiles into a single table'''
    PipelineBamStats.loadTranscriptProfile(infiles, outfile)


@P.add_doc(PipelineBamStats.loadStrandSpecificity)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@follows(loadTranscriptProfile)
@merge(strandSpecificity, "strand_spec.load")
def loadStrandSpecificity(infiles, outfile):
    ''' merge strand specificity data into a single table'''
    PipelineBamStats.loadStrandSpecificity(infiles, outfile)


# ---------------------------------------------------
# Generic pipeline tasks
# These tasks allow ruffus to pipeline tasks together


@follows(loadPicardStats,
         loadPicardDuplicationStats,
         loadBAMStats,
         loadContextStats,
         buildIntronLevelReadCounts,
         loadIdxStats,
         loadExonValidation,
         loadPicardRnaSeqMetrics,
         loadTranscriptProfile,
         loadStrandSpecificity)
def full():
    '''a dummy task to run all tasks in the pipeline'''
    pass


# --------------------------------------------------
# Reporting tasks
# --------------------------------------------------
@follows(mkdir("R_report.dir"))
def renderRreport():
    '''build R markdown report '''

    report_path = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               'pipeline_docs',
                                               'pipeline_bamstats',
                                               'R_report'))

    statement = '''cp %(report_path)s/* R_report.dir ; cd R_report.dir ; R -e "rmarkdown::render_site()"'''

    P.run()


@follows(mkdir("Jupyter_report.dir"))
def renderJupyterReport():
    '''build Jupyter notebook report'''

    report_path = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               'pipeline_docs',
                                               'pipeline_bamstats',
                                               'Jupyter_report'))

    statement = ''' cp %(report_path)s/* Jupyter_report.dir/ ; cd Jupyter_report.dir/;
                    jupyter nbconvert --ExecutePreprocessor.timeout=None --allow-errors --to html --execute *.ipynb;
                    mkdir _site;
                    mv -t _site *.html cgat_logo.jpeg oxford.png'''

    P.run()


@follows(mkdir("MultiQC_report.dir"))
@originate("MultiQC_report.dir/multiqc_report.html")
def renderMultiqc(infile):
    '''build mulitqc report'''

    statement = '''LANG=en_GB.UTF-8 multiqc . -f;
                   mv multiqc_report.html MultiQC_report.dir/'''

    P.run()


@follows(renderRreport,
         renderJupyterReport,
         renderMultiqc)
def build_report():
    '''report dummy task to build reports'''
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
