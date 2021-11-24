"""========================================
RNA-Seq Differential expression pipeline
========================================


The RNA-Seq differential expression pipeline performs differential
expression analysis. It requires three inputs:

   1. A geneset in :term:`gtf` formatted file
   2. Mapped reads in :term:`bam` formatted files and/or unaligned reads in
      :term:`fastq` formatted files
   3. Design files as :term:`tsv`-separated format

This pipeline works on a single genome.

Overview
========

The pipeline performs the following:

   * Compute tag (counts) at the transcript and gene level
     The following counting methods are implemented:
      * featureCounts_
      * gtf2table

    and/or

   * Gene expression estimates (TPM and counts) at the transcript and
     gene level. The following alignment-free expression estimation
     methods are implemented:
      * kallisto_
      * salmon_

   * Perform differential expression analysis. The methods currently
     implemented are:

      * deseq2_
      * edger_
      * sleuth_

Background
============

Quantification:

Transcripts are the natural choice to measure expression. However
other quantities might be of interest. Some quantities are biological
meaningful, for example differential expression from a promotor shared
by several trancripts. Other quantities might no biologically
meaningful but are necessary as a technical comprise.  For example,
the overlapping transcripts might be hard to resolve and thus might
need to be aggregated per gene. Furthermore, functional annotation is
primarily associated with genes and not individual transcripts.

This pipeline estimates transcript and gene-level expression and
performs differential expression analysis on both levels.

The quantification tools fall into two categories:
   * Alignment-free
      Quantification is performed directly from the raw reads against
      a reference transcriptome using "pseduo-alignment". In essence,
      the tool attempts to identify the compatible transcripts for
      each read without identifying the exact alignment position of
      the read on the transcript or genome. Following this the set of
      expression estimates which best explain the observed reads are
      obtained using an Expectation Maximisation approach

      The available tools are:
      * kallisto_
      * salmon_

   * Alignment-based
      Quantification is performed on the aligned reads using the
      position of features described in the reference geneset
      gtf. Reads are discretely assigned to one feature (may be
      performed at the transcript or gene-level).  It should be noted
      that transcript-level quantification with tag counting methods
      is inherrently inaccurate since a read which aligns to an exon
      present in multiple isoforms of the same gene can not be naively
      assigned to a single transcript.

      The available tools are:
      * featurecounts_
      * gtf2table (in-house script)

The alignment-free methods should be preffered over featureCounts and
gtf2table in almost all instances. However, many analyses still use
tag counting so it may be neccessary to repeat other groups
analyses. In addition gtf2table provides more detailed tag counting
which may be useful when exploring problematic RNA-Seq
data. Alignment-free methods also provide estimated counts per
transcript which can be rounded to integer counts.


Differential expression:

Differential expression can be performed on (non normalised) counts per
transcript/gene or on expression estimates (Transripts Per Million =
TPM) per transcript/gene.

Count-based expression estimates (alignment-free & alignment-based)
are well modelled by a negative binomial distribution and differential
expression can therefore be performed with a negative binomial
Generalised Linear Model (GLM). This is the approach taken by DESeq2
and edgeR which are both used here. Most simply, a Wald test can be
performed to identify genes where the log2fold change between two
levels of a factor (e.g factor = genotype, levels = WT, KO) is
significantly different from 0. Where the factor has more than one
level (e.g factor = genotype, levels = WT, KO1, KO2), a Likelihood
Ratio Test (LRT) can be performed to identify genes where a full model
including the (e.g genotype) factor is a signficantly better fit than
a reduced model not including the said factor. This pipeline performs Wald
test only. Please see the deseq2_/edger_ vingettes for LRT.

Log TPM (alignment-free only) are well modelled by a gaussian
distribution and differential expression can therefore be performed
with a linear model. This is the approach taken by sleuth which is
used here. In addition, slueth_ uses the bootstrap estimates from
kallisto_/salmon_ to estimate the proportion of the variance
which is technical and therefore the proportion which is biological.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use cgat pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.

The sphinxreport report requires a :file:`conf.py` and
:file:`sphinxreport.yml` file (see :ref:`PipelineReporting`). To start
with, use the files supplied with the Example_ data.

Input
-----

Reads
+++++

Reads are imported by placing :term:`bam` formatted files are linking
to files in the :term:`working directory`.

The default file format assumes the following convention::

   <samplename>.bam (aligned reads)
   <samplename>.fastq.gz (fastq.1.gz and .sra are also accepted for raw reads)

To compare alignment-free and alignment-based methods, the raw reads
and aligned reads must both be supplied

Geneset
++++++++

The Geneset is specified by the "geneset" parameter

Design matrices
+++++++++++++++

Design matrices are imported by placing :term:`tsv` formatted files
into the :term:`working directory`. A design matrix describes the
experimental design to test. The design files should be named
design*.tsv.

Each design file has at leasr four columns but may contain any number
of columns after the 'pair' column:

      track   include group   pair
      CW-CD14-R1      0       CD14    1
      CW-CD14-R2      0       CD14    1
      CW-CD14-R3      1       CD14    1
      CW-CD4-R1       1       CD4     1
      FM-CD14-R1      1       CD14    2
      FM-CD4-R2       0       CD4     2
      FM-CD4-R3       0       CD4     2
      FM-CD4-R4       0       CD4     2

track
     name of track - should correspond to a sample name.
include
     flag to indicate whether or not to include this data
group
     group indicator - experimental group
pair
     pair that sample belongs to (for paired tests) - set to 0 if the
     design is not paired.


Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default cgat setup, the pipeline requires the following
software to be in the path:

+--------------+----------+------------------------------------+
|*Program*     |*Version* |*Purpose*                           |
+--------------+----------+------------------------------------+
|samtools      |>=0.1.16  |bam/sam files                       |
+--------------+----------+------------------------------------+
|bedtools      |          |working with intervals              |
+--------------+----------+------------------------------------+
|deseq2_/edgeR_|          |count-based differential expression |
+--------------+----------+------------------------------------+
|sleuth_       |          |TPM-based differential expression   |
+--------------+----------+------------------------------------+
|samtools      |>=0.1.16  |bam/sam files                       |
+--------------+----------+------------------------------------+
|featureCounts_|>=1.4.6   |alignment-based quantification      |
+--------------+----------+------------------------------------+
|gtf2table     |          |alignment-based quantification      |
+--------------+----------+------------------------------------+
|kallisto_     |>=0.43.0  |alignment-free quantification       |
+--------------+----------+------------------------------------+
|salmon_       |>=0.7.2   |alignment-free quantification       |
+--------------+----------+------------------------------------+



Pipeline output
===============

Quantification
--------------

The quantification estimates from each method are outputted to:
[method].dir/[sample]/[level].tsv.gz,
where [method] is the quantification method, [sample] is the sample
name and [level] is the feature level (transcript or gene)

Each tool also generates specific log files etc which are outputted,
along with the raw quantification outfile in the directory:
[method].dir/[sample]

For each method, the merged counts are outputted to:
[method].dir/[level].tsv.gz


Differential gene expression results
-------------------------------------

Results are stored per method in subdirectories
such as :file:`deseq.dir`, :file:`edger.dir` or :file:`sleuth.dir`

Plots from the differential expression analyses are also contained
within the directories.


Glossary
========

.. glossary::

   kallisto
      kallisto_ - alignment-free quantification
   salmon
      salmon_ - alignment-free quantification
   featureCounts
      featurecounts_ - alignment-free quantification
   deseq
      deseq_ - differential expression analysis
   edger
      edger_ - differential expression analysis
   sleuth
      sleuth_ - differential expression analysis

.. _featurecounts: http://bioinf.wehi.edu.au/featureCounts/
.. _kallisto: https://pachterlab.github.io/kallisto/
.. _salmon: https://combine-lab.github.io/salmon/
.. _deseq: http://www-huber.embl.de/users/anders/DESeq/
.. _edger: http://bioconductor.org/packages/release/bioc/html/edgeR.html
.. _sleuth: https://github.com/pachterlab/sleuth


ChangeLog
=========

28.03.2014  Andreas Heger
            added automated selection of paired counting to featureCounts
            and gtf2table counting.

11.4.2014   Andreas Heger
            changed workflow. Multiple counters are applied and
            differential expression is computed on all.

15.10.2015  Charlotte George, Sebastian Luna-Valero
            SCRUM Oct 2015. Updating documentation.

10.10.2016  cgat Fellows. SCRUM Oct 2016. Complete re-write of
            pipeline and modules to simplify workflow and add alignment-free
            methods


###########################################################################
Possible improvements:
###########################################################################
Tom Smith 17 OCT 16:
Add Stringtie + Cufflinks2? - Need to work out how to extract counts
A python script (prepDE.py) is available from:
https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#deseq but
this doesn't appear to be under version control and does not work
directly from the stringtie/cufflinks2 output

###########################################################################

Code
====

"""

# load modules
from ruffus import *
from ruffus.combinatorics import *

import cgatcore.experiment as E
# import cgat.scrum_expression as SE

import sys
import os
import re
import glob
import pandas as pd
import sqlite3
import cgat.GTF as GTF
import cgatcore.iotools as iotools
from cgatcore import pipeline as P

import cgatpipelines.tasks.geneset as geneset
import cgatpipelines.tasks.rnaseq as rnaseq
import cgatpipelines.tasks.tracks as tracks
import cgatpipelines

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

PARAMS = P.PARAMS
PARAMS.update(P.peek_parameters(
    PARAMS["annotations_dir"],
    "genesets",
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))

geneset.PARAMS = PARAMS

# Helper functions mapping tracks to conditions, etc
# determine the location of the input files (reads).
try:
    PARAMS["input"]
except KeyError:
    DATADIR = "."
else:
    if PARAMS["input"] == 0:
        DATADIR = "."
    elif PARAMS["input"] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS["input"]  # not recommended practise.

Sample = tracks.AutoSample

# collect sra and fastq.gz tracks
BAM_TRACKS = tracks.Tracks(Sample).loadFromDirectory(
    glob.glob("*.bam"), "(\S+).bam")

DESIGNS = tracks.Tracks(Sample).loadFromDirectory(
    glob.glob("design*.tsv"), "design(\S+).tsv")


# do not use - legacy methods
# here only to stop ruffus erroring. Remove once pipleine scrum is
# complete and old code has been removed
GENESETS = tracks.Tracks(Sample).loadFromDirectory(
    glob.glob("*.gtf.gz"), "(\S+).gtf.gz")
TRACKS = tracks.Tracks(Sample).loadFromDirectory(
    glob.glob("*.bam"), "(\S+).bam")
# group by experiment (assume that last field is a replicate identifier)
EXPERIMENTS = tracks.Aggregate(
    BAM_TRACKS, labels=("condition", "tissue"))


###############################################################################
# Utility function
###############################################################################

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


###############################################################################
# build indexes
###############################################################################

@mkdir('geneset.dir')
@transform(PARAMS['geneset'],
           regex("(\S+).gtf.gz"),
           r"geneset.dir/\1.fa")
def buildReferenceTranscriptome(infile, outfile):
    '''
    Builds a reference transcriptome from the provided GTF geneset - generates
    a fasta file containing the sequence of each feature labelled as
    "exon" in the GTF.
    --fold-at specifies the line length in the output fasta file

    Parameters
    ----------
    infile: str
        path to the GTF file containing transcript and gene level annotations
    genome_dir: str
        :term: `PARAMS` the directory of the reference genome
    genome: str
        :term: `PARAMS` the filename of the reference genome (without .fa)
    outfile: str
        path to output file
    '''

    genome_file = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))

    statement = '''
    zcat %(infile)s |
    awk '$3=="exon"'|
    cgat gff2fasta
    --is-gtf --genome-file=%(genome_file)s --fold-at=60 -v 0
    --log=%(outfile)s.log > %(outfile)s &&
    samtools faidx %(outfile)s
    '''

    P.run(statement)


@transform(buildReferenceTranscriptome,
           suffix(".fa"),
           ".kallisto.index")
def buildKallistoIndex(infile, outfile):
    '''
    Builds a kallisto index for the reference transcriptome

    Parameters
    ----------
    infile: str
       path to reference transcriptome - fasta file containing transcript
       sequences
    kallisto_kmer: int
       :term: `PARAMS` kmer size for Kallisto.  Default is 31.
       Kallisto will ignores transcripts shorter than this.
    outfile: str
       path to output file

    '''

    job_memory = "12G"

    statement = '''
    kallisto index -i %(outfile)s -k %(kallisto_kmer)s %(infile)s
    '''

    P.run(statement)

@transform(buildReferenceTranscriptome,
           suffix(".fa"),
           ".gentrome.fa")
def buildGentrome(infile, outfile):
    '''
    Builds a "Gentrome", a concatenation of the reference transcriptome and genome
    FASTA file, required by Salmon indexing

    Parameters
    ----------
    infile: str
        path to the reference transcriptome generated in the previous function
    genome_dir: str
        :term: `PARAMS` the directory of the reference genome
    genome: str
        :term: `PARAMS` the filename of the reference genome (without .fa)
    outfile: str
        path to output file
    '''

    genome_file = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))

    outdir = os.path.dirname(outfile)

    statement = '''
    grep "^>" %(genome_file)s |
    cut -d " " -f 1 > %(outdir)s/decoys.txt &&
    sed -i.bak -e 's/>//g' %(outdir)s/decoys.txt &&
    cat %(infile)s %(genome_file)s > %(outfile)s
    '''

    P.run(statement)
    

@transform(buildGentrome,
           suffix(".gentrome.fa"),
           ".salmon.index")
def buildSalmonIndex(infile, outfile):
    '''
    Builds a salmon index for the reference transriptome
    Parameters
    ----------
    infile: str
       path to reference transcriptome - fasta file containing transcript
       sequences
    salmon_kmer: int
       :term: `PARAMS` kmer size for salmon.  Default is 31.
       Salmon will ignores transcripts shorter than this.
    salmon_index_options: str
       :term: `PARAMS` string to append to the salmon index command to
       provide specific options e.g. --force --threads N
    outfile: str
       path to output file
    '''

    job_memory = "unlimited"
    threads = 12
    if PARAMS["salmon_threads"] is not None:
        threads = PARAMS["salmon_threads"]
    outdir = os.path.dirname(outfile)
    
    # need to remove the index directory (if it exists) as ruffus uses
    # the directory timestamp which wont change even when re-creating
    # the index files
    statement = '''
    rm -rf %(outfile)s;
    salmon index -p %(threads)s -k %(salmon_kmer)i -t %(infile)s -d %(outdir)s/decoys.txt -i %(outfile)s
    '''

    P.run(statement)


@originate("transcript2geneMap.tsv")
def getTranscript2GeneMap(outfile):
    ''' Extract a 1:1 map of transcript_id to gene_id from the geneset '''

    iterator = GTF.iterator(iotools.open_file(PARAMS['geneset']))
    transcript2gene_dict = {}

    for entry in iterator:

        # Check the same transcript_id is not mapped to multiple gene_ids!
        if entry.transcript_id in transcript2gene_dict:
            if not entry.gene_id == transcript2gene_dict[entry.transcript_id]:
                raise ValueError('''multipe gene_ids associated with
                the same transcript_id %s %s''' % (
                    entry.gene_id,
                    transcript2gene_dict[entry.transcript_id]))
        else:
            transcript2gene_dict[entry.transcript_id] = entry.gene_id

    with iotools.open_file(outfile, "w") as outf:
        outf.write("transcript_id\tgene_id\n")
        for key, value in sorted(transcript2gene_dict.items()):
            outf.write("%s\t%s\n" % (key, value))


###################################################
# count-based quantifiers
###################################################

@active_if("featurecounts" in P.as_list(PARAMS["quantifiers"]))
@follows(mkdir("featurecounts.dir"))
@transform(["%s.bam" % x.asFile() for x in BAM_TRACKS],
           regex("(\S+).bam"),
           add_inputs(PARAMS['geneset']),
           [r"featurecounts.dir/\1/transcripts.tsv.gz",
            r"featurecounts.dir/\1/genes.tsv.gz"])
def runFeatureCounts(infiles, outfiles):
    '''
    Counts reads falling into "features" - in each transcript and
    each gene.

    A read is counted as overlapping with a feature if at least one bp
    overlaps.

    Pairs and strandedness can be used to resolve reads falling into
    more than one feature. Reads that cannot be resolved to a single
    feature are ignored.

    The raw output of featureCounts is sent to
    featurecounts.dir/SAMPLEID/transcripts.tsv.raw.gz and
    featurecounts.dir/SAMPLEID/genes.tsv.raw.gz

    Parsed output is sent to featurecounts.dir/SAMPLEID/transcripts.tsv.gz
    and featurecounts.dir/SAMPLEID/genes.tsv.gz

    Other default output files are stored in the featurecounts.dir/SAMPLEID
    directory.

    See feature counts manual http://bioinf.wehi.edu.au/featureCounts/
    for information about :term:`PARAMS` options or look at
    featureCounts --help

    Parameters
    ----------
    infiles : list
        List with two components:
        0 - list of file names of the bam formatted files containing the
        aligned reads
        1 - file name of the GTF file containing the features over which
        to count.

    featurecounts_threads : int
        :term:`PARAMS` - number of threads to run feature counts. This is
        specified in pipeline.yml

    featurecounts_strand : int
        :term:`PARAMS`
        0: unstranded
        1: the first read of the pair is on sense relative to transcript
        2: the first read of the pair is on antisense relative to the
        transcript

    featurecounts_options : string
        :term:`PARAMS` - options string for running feature counts e.g.
         -Q flag specifies minimum mapping quality. Set to 10 to be
         -M will allow multi mapping reads
         -O will allow reads to overlap more than one feature
        more in featureCounts --help

    transcript_outfile/gene_outfile : string used to denote output
        files from feature counts using transcript_ids or gene_ids.

    Eight output files are produced for each input :term:`bam`.
        SAMPLENAME/transcripts.tsv.raw.gz
        SAMPLENAME/genes.tsv.raw.gz
        These are the raw count tables produced by featureCounts

        SAMPLENAME/transcripts.tsv.gz
        SAMPLENAME/genes.tsv.gz
        These are standardised count tables with two columns showing the
        feature ID and the number of reads overlapping with this feature

        SAMPLENAME/transcripts.tsv.raw.summary
        SAMPLENAME/genes.tsv.raw.summary
        Summary of the reads which have been counted specifiying the number
        of multimapping reads, duplicates etc.

        SAMPLENAME/transcripts.tsv.gz.log
        SAMPLENAME/genes.tsv.gz.log
        Log file for each analysis
    '''
    bamfile, annotations = infiles
    transcript_outfile, gene_outfile = outfiles
    Quantifier = rnaseq.FeatureCountsQuantifier(
        infile=bamfile,
        transcript_outfile=transcript_outfile,
        gene_outfile=gene_outfile,
        job_threads=PARAMS['featurecounts_threads'],
        strand=PARAMS['featurecounts_strand'],
        options=PARAMS['featurecounts_options'],
        annotations=annotations)

    Quantifier.run_all()


@follows(mkdir("gtf2table.dir"))
@transform(["%s.bam" % x.asFile() for x in BAM_TRACKS],
           regex("(\S+).bam"),
           add_inputs(PARAMS['geneset']),
           [r"gtf2table.dir/\1/transcripts.tsv.gz",
            r"gtf2table.dir/\1/genes.tsv.gz"])
def runGTF2Table(infiles, outfiles):
    '''
    Compute read counts and coverage of transcripts and genes using the
    cgat gtf2table tools.

    Takes a list of :term:`bam` files defined in "BAM_TRACKS" and a
    :term:`gtf` file containing transcript and gene level annotations
    and produces `.tsv.gz`
    file using gtf2table.py detailing coverage of genes and transcripts by
    reads for each bam.

    Appropriate lines in the GTF are labelled as "exon".

    The following gtf2table "counters" provide columns in the output
    read-coverage - outputs the number of reads overlapping with the
    transcript / gene model, the number of bases in the overlap and summary
    statistics of the coverage per base
    length - outputs the number of exons in each transcript / gene and exon
    length summary statistics
    read-counts or readpair-counts - outputs the number of reads overlapping
    with a gene or transcript.

    .. note::
        This tool ignores multimapping reads

    Parameters
    ----------
    infiles : list
        List with two components:
        0: list of bam file names
        1: path to gtf file over which to count

    transcript_outfile/gene_outfile : string used to denote output
        from feature counts using transcript_ids or gene_ids.

    Six output files are produced for each input :term:`bam`:
        SAMPLENAME/transcripts.tsv.raw.gz
        SAMPLENAME/genes.tsv.raw.gz
        These are the raw count tables produced by gtf2table

        SAMPLENAME/transcripts.tsv.gz
        SAMPLENAME/genes.tsv.gz
        These are standardised count tables with two columns showing the
        feature ID and the number of reads overlapping with this feature

        SAMPLENAME/transcripts.tsv.gz.log
        SAMPLENAME/genes.tsv.gz.log
        Log file for each analysis

    '''
    bamfile, annotations = infiles
    transcript_outfile, gene_outfile = outfiles
    Quantifier = rnaseq.Gtf2tableQuantifier(
        infile=bamfile,
        transcript_outfile=transcript_outfile,
        gene_outfile=gene_outfile,
        annotations=annotations)

    Quantifier.run_all()


###################################################
###################################################
# alignment-free quantifiers
###################################################
###################################################

###################################################
# Define quantification regex and output
###################################################

SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.sra")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

# enable multiple fastqs from the same sample to be analysed together
if "merge_pattern_input" in PARAMS and PARAMS["merge_pattern_input"]:
    SEQUENCEFILES_REGEX = regex(
        r"%s.(fastq.1.gz|fastq.gz|sra)" % (
            PARAMS["merge_pattern_input"].strip()))

    # the last expression counts number of groups in pattern_input
    SEQUENCEFILES_KALLISTO_OUTPUT = [
        r"kallisto.dir/%s/transcripts.tsv.gz" % (
            PARAMS["merge_pattern_output"].strip()),
        r"kallisto.dir/%s/genes.tsv.gz" % (
            PARAMS["merge_pattern_output"].strip())]

    SEQUENCEFILES_SALMON_OUTPUT = [
        r"salmon.dir/%s/transcripts.tsv.gz" % (
            PARAMS["merge_pattern_output"].strip()),
        r"salmon.dir/%s/genes.tsv.gz" % (
            PARAMS["merge_pattern_output"].strip())]

else:
    SEQUENCEFILES_REGEX = regex(
        "(\S+).(fastq.1.gz|fastq.gz|sra)")

    SEQUENCEFILES_KALLISTO_OUTPUT = [
        r"kallisto.dir/\1/transcripts.tsv.gz",
        r"kallisto.dir/\1/genes.tsv.gz"]

    SEQUENCEFILES_SALMON_OUTPUT = [
        r"salmon.dir/\1/transcripts.tsv.gz",
        r"salmon.dir/\1/genes.tsv.gz"]

###################################################


@follows(mkdir("kallisto.dir"))
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(buildKallistoIndex, getTranscript2GeneMap),
         SEQUENCEFILES_KALLISTO_OUTPUT)
def runKallisto(infiles, outfiles):
    '''
    Computes read counts across transcripts and genes based on a fastq
    file and an indexed transcriptome using Kallisto.

    Runs the kallisto "quant" function across transcripts with the specified
    options.  Read counts across genes are counted as the total in all
    transcripts of that gene (based on the getTranscript2GeneMap table)

    Parameters
    ----------
    infiles: list
        list with three components
        0 - list of strings - paths to fastq files to merge then quantify
        across using Kallisto
        1 - string - path to Kallisto index file
        2 - string - path totable mapping transcripts to genes

    kallisto_threads: int
       :term: `PARAMS` the number of threads for Kallisto
    kallisto_memory: str
       :term: `PARAMS` the job memory for Kallisto
    kallisto_options: str
       :term: `PARAMS` string to append to the Kallisto quant command to
       provide specific
       options, see https://pachterlab.github.io/kallisto/manual
    kallisto_bootstrap: int
       :term: `PARAMS` number of bootstrap samples to run.
       Note, you need to bootstrap for differential expression with sleuth
       if there are no technical replicates. If you only need point estimates,
       set to 1.  Note that bootstrap must be set to at least 1
    kallisto_fragment_length: int
       :term: `PARAMS` Fragment length for Kallisto, required for single end
       reads only
    kallisto_fragment_sd: int
       :term: `PARAMS` Fragment length standard deviation for Kallisto,
       required for single end reads only.
    outfiles: list
       paths to output files for transcripts and genes
    '''

    # TS more elegant way to parse infiles and index?
    fastqfile = [x[0] for x in infiles]
    index = infiles[0][1]
    transcript2geneMap = infiles[0][2]

    transcript_outfile, gene_outfile = outfiles
    Quantifier = rnaseq.KallistoQuantifier(
        infile=fastqfile,
        transcript_outfile=transcript_outfile,
        gene_outfile=gene_outfile,
        annotations=index,
        job_threads=PARAMS["kallisto_threads"],
        job_memory=PARAMS["kallisto_memory"],
        options=PARAMS["kallisto_options"],
        bootstrap=PARAMS["kallisto_bootstrap"],
        fragment_length=PARAMS["kallisto_fragment_length"],
        fragment_sd=PARAMS["kallisto_fragment_sd"],
        transcript2geneMap=transcript2geneMap)

    Quantifier.run_all()


@follows(mkdir("salmon.dir"))
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(buildSalmonIndex, getTranscript2GeneMap),
         SEQUENCEFILES_SALMON_OUTPUT)
def runSalmon(infiles, outfiles):
    '''
    Computes read counts across transcripts and genes based on a fastq
    file and an indexed transcriptome using Salmon.

    Runs the salmon "quant" function across transcripts with the specified
    options.  Read counts across genes are counted as the total in all
    transcripts of that gene (based on the getTranscript2GeneMap table)

    Parameters
    ----------
    infiles: list
        list with three components
        0 - list of strings - paths to fastq files to merge then quantify
        across using salmon
        1 - string - path to salmon index file
        2 - string - path to table mapping transcripts to genes

    salmon_threads: int
       :term: `PARAMS` the number of threads for salmon
    salmon_memory: str
       :term: `PARAMS` the job memory for salmon
    salmon_options: str
       :term: `PARAMS` string to append to the salmon quant command to
       provide specific
       options, see https://salmon.readthedocs.io/en/latest/salmon.html#description-of-important-optionsç∂
    salmon_bootstrap: int
       :term: `PARAMS` number of bootstrap samples to run.
       Note, you need to bootstrap for differential expression with sleuth
       if there are no technical replicates. If you only need point estimates,
       set to 1.
    salmon_libtype: str
       :term: `PARAMS` salmon library type, https://salmon.readthedocs.io/en/latest/library_type.html
    outfiles: list
       paths to output files for transcripts and genes
    '''
    fastqfile = [x[0] for x in infiles]
    index = infiles[0][1]
    transcript2geneMap = infiles[0][2]

    transcript_outfile, gene_outfile = outfiles
    Quantifier = rnaseq.SalmonQuantifier(
        infile=fastqfile,
        transcript_outfile=transcript_outfile,
        gene_outfile=gene_outfile,
        annotations=index,
        job_threads=PARAMS["salmon_threads"],
        job_memory=PARAMS["salmon_memory"],
        options=PARAMS["salmon_options"],
        bootstrap=PARAMS["salmon_bootstrap"],
        libtype=PARAMS['salmon_libtype'],
        kmer=PARAMS['salmon_kmer'],
        transcript2geneMap=transcript2geneMap)

    Quantifier.run_all()


###################################################
###################################################
# Create quantification targets
###################################################

QUANTTARGETS = []
mapToQuantTargets = {'kallisto': (runKallisto,),
                     'salmon': (runSalmon,),
                     'featurecounts': (runFeatureCounts,),
                     'gtf2table': (runGTF2Table,)}

for x in P.as_list(PARAMS["quantifiers"]):
    QUANTTARGETS.extend(mapToQuantTargets[x])


###################################################


@collate(QUANTTARGETS,
         regex("(\S+).dir/(\S+)/transcripts.tsv.gz"),
         [r"\1.dir/transcripts.tsv.gz",
          r"\1.dir/genes.tsv.gz"])
def mergeCounts(infiles, outfiles):
    ''' merge counts for alignment-based methods'''

    transcript_infiles = [x[0] for x in infiles]
    gene_infiles = [x[1] for x in infiles]

    transcript_outfile, gene_outfile = outfiles

    def mergeinfiles(infiles, outfile):
        final_df = pd.DataFrame()

        for infile in infiles:
            tmp_df = pd.read_table(infile, sep="\t", index_col=0)
            final_df = final_df.merge(
                tmp_df, how="outer",  left_index=True, right_index=True)

        final_df = final_df.round()
        final_df.sort_index(inplace=True)
        final_df.to_csv(outfile, sep="\t", compression="gzip")

    mergeinfiles(transcript_infiles, transcript_outfile)
    mergeinfiles(gene_infiles, gene_outfile)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(mergeCounts, regex("(\S+).dir/transcripts.tsv.gz"),
           [r"\1.dir/\1_transcripts.load",
            r"\1.dir/\1_genes.load"])
def loadMergedCounts(infiles, outfiles):
    P.load(infiles[0], outfiles[0])
    P.load(infiles[1], outfiles[1])


@active_if("featurecounts" in P.as_list(PARAMS["quantifiers"]))
@collate(runFeatureCounts,
         regex("featurecounts.dir/([^.]+)/([^.]+).tsv.gz"),
         r"featurecounts.dir/genelength.tsv.gz")
def mergeLengths(infiles, outfile):
    ''' build a matrix of "genelengths" derived from FeatureCounts
    with genes and tracks dimensions.
    '''
    raw_infiles = [x[1].replace("gz", "raw.gz") for x in infiles]
    final_df = pd.DataFrame()

    for infile in raw_infiles:
        tmp_df = pd.read_table(infile, sep="\t", index_col=0,
                               comment='#', usecols=["Geneid", "Length"])
        m = re.search('featurecounts.dir\/(.+?)\/genes.tsv.raw.gz', infile)
        if m:
            tmp_df.columns = [m.group(1)]
        final_df = final_df.merge(
            tmp_df, how="outer",  left_index=True, right_index=True)

    final_df.sort_index(inplace=True)
    final_df.to_csv(outfile, sep="\t", compression="gzip")


@active_if("featurecounts" in P.as_list(PARAMS["quantifiers"]))
@transform(mergeLengths,
           suffix(".tsv.gz"),
           ".load")
def loadMergedLengths(infile, outfile):
    '''load genength table into database'''
    P.load(infile, outfile, "--add-index=gene_id")


@follows(*QUANTTARGETS)
@follows(loadMergedCounts,
         loadMergedLengths)
def count():
    ''' dummy task to define upstream quantification tasks'''

###################################################
# Loading and filtering counts table
###################################################

@mkdir("DEresults.dir/deseq2")
@product(mergeCounts,
         formatter(".*/(?P<QUANTIFIER>\S+).dir/transcripts.tsv.gz"),
         ["design%s.tsv" % x.asFile() for x in DESIGNS],
         formatter(".*/design(?P<DESIGN>\S+).tsv$"),
         "DEresults.dir/deseq2/{QUANTIFIER[0][0]}_{DESIGN[1][0]}/experiment_out.rds",
         "{DESIGN[1][0]}",
         "{QUANTIFIER[0][0]}")
def filterDESeq2(infiles, outfile, design_name, quantifier_name):
    ''' Load counts into RDS object and filter'''

    counts, design = infiles
    transcripts, genes = counts
    quantifier_name = quantifier_name.lower()
    counts = "--counts-dir %s.dir" % quantifier_name     
    model = PARAMS.get('deseq2_model%s' % design_name, None)
        
    if not quantifier_name in ("salmon", "kallisto"):
        counts = "--counts-tsv %s" % genes  
        quantifier_name = "counts_table"    
    
    outdir = os.path.dirname(outfile)
    r_root = os.path.abspath(os.path.dirname(cgatpipelines.__file__))
    scriptpath = os.path.join(r_root, "Rtools/filtercounts.R")
    

    statement = '''
    export R_ROOT=%(r_root)s &&
    Rscript %(scriptpath)s %(counts)s
    --sampleData %(design)s
    --outdir %(outdir)s
    --model %(model)s
    --method deseq2
    --filter %(filter_deseq2)s
    --source %(quantifier_name)s
    --tx2gene_regex %(filter_regex)s
    > %(outdir)s/filter.log;
    '''
    P.run(statement) 


@mkdir("DEresults.dir/edger")
@product(mergeCounts,
         formatter(".*/(?P<QUANTIFIER>\S+).dir/transcripts.tsv.gz"),
         ["design%s.tsv" % x.asFile() for x in DESIGNS],
         formatter(".*/design(?P<DESIGN>\S+).tsv$"),
         "DEresults.dir/edger/{QUANTIFIER[0][0]}_{DESIGN[1][0]}/experiment_out.rds",
         "{DESIGN[1][0]}",
         "{QUANTIFIER[0][0]}")
def filterEdgeR(infiles, outfile, design_name, quantifier_name):
    ''' Load counts into RDS object and filter'''

    counts, design = infiles
    transcripts, genes = counts
    quantifier_name = quantifier_name.lower()
    counts = "--counts-dir %s.dir" % quantifier_name     
    model = PARAMS.get('edger_model%s' % design_name, None)
        
    if not quantifier_name in ("salmon", "kallisto"):
        counts = "--counts-tsv %s" % genes  
        quantifier_name = "counts_table"    
    
    outdir = os.path.dirname(outfile)
    r_root = os.path.abspath(os.path.dirname(cgatpipelines.__file__))
    scriptpath = os.path.join(r_root, "Rtools/filtercounts.R")
    

    statement = '''
    export R_ROOT=%(r_root)s &&
    Rscript %(scriptpath)s %(counts)s
    --sampleData %(design)s
    --outdir %(outdir)s
    --model %(model)s
    --method edger
    --filter %(filter_edger)s
    --source %(quantifier_name)s
    --tx2gene_regex %(filter_regex)s
    > %(outdir)s/filter.log;
    '''
    P.run(statement) 


###################################################
# Exploratory Analysis
###################################################

@transform(filterDESeq2,
           formatter("DEresults.dir/(?P<DETOOL>\S+)/(?P<QUANTIFIER>\S+)_(?P<DESIGN>\S+)/experiment_out.rds"),
           "DEresults.dir/{DETOOL[0]}/{QUANTIFIER[0]}_{DESIGN[0]}/Heatmap_top500.png",
           "{DESIGN[0]}",
           "{QUANTIFIER[0]}",
           "{DETOOL[0]}")
def exploratoryAnalysis(infile, outfile, design_name, quantifier_name, detool_name):


    design = "design" + design_name + ".tsv"
    outdir = os.path.dirname(outfile)
    r_root = os.path.abspath(os.path.dirname(cgatpipelines.__file__))
    scriptpath = os.path.join(r_root, "Rtools/exploratory.R")

    model = PARAMS.get('%s_model%s' % (detool_name, design_name), None)
    contrast = PARAMS.get('%s_contrast%s' % (detool_name, design_name), None)

    if model is None:
        raise ValueError("{}_model{} is not specified".format(
            (detool_name, design_name)))
    if contrast is None:
        raise ValueError("{}_contrast{} is not specified".format(
            (detool_name, design_name)))

    statement = '''
    export R_ROOT=%(r_root)s &&
    Rscript %(scriptpath)s
    --rds-filename %(infile)s
    --model %(model)s
    --contrast %(contrast)s
    --factors %(exploratory_factors)s
    --genes_of_interest %(exploratory_goi)s
    --outdir %(outdir)s
    > %(outdir)s/exploratory.log;
    '''
    P.run(statement)


###################################################
# Differential Expression
###################################################


@transform(filterDESeq2,
           formatter("DEresults.dir/deseq2/(?P<QUANTIFIER>\S+)_(?P<DESIGN>\S+)/experiment_out.rds"),
           "DEresults.dir/deseq2/{QUANTIFIER[0]}_{DESIGN[0]}/results_full.tsv",
           "{DESIGN[0]}",
           "{QUANTIFIER[0]}")
def runDESeq2(infile, outfile, design_name, quantifier_name):
    ''' run DESeq2 to identify differentially expression transcripts/genes'''

    design = "design" + design_name + ".tsv"
    outdir = os.path.dirname(outfile)
    r_root = os.path.abspath(os.path.dirname(cgatpipelines.__file__))
    scriptpath = os.path.join(r_root, "Rtools/diffexpression.R")

    model = PARAMS.get('deseq2_model%s' % design_name, None)
    contrast = PARAMS.get('deseq2_contrast%s' % design_name, None)
    refgroup = PARAMS.get('deseq2_refgroup%s' % design_name, None)
    coef = PARAMS.get('deseq2_coef%s' % design_name, None)

    if model is None:
        raise ValueError("deseq2_model{} is not specified".format(
            design_name))
    if contrast is None:
        raise ValueError("deseq2_contrast{} is not specified".format(
            design_name))
    if refgroup is None:
        raise ValueError("deseq2_refgroup{} is not specified".format(
            design_name))
    if coef is None:
        raise ValueError("deseq2_coef{} is not specified".format(
            design_name))

    statement = '''
    export R_ROOT=%(r_root)s &&
    Rscript %(scriptpath)s
    --rds-filename %(infile)s
    --model %(model)s
    --contrast %(contrast)s
    --refgroup %(refgroup)s
    --coef %(coef)s
    --alpha %(deseq2_fdr)s
    --outdir %(outdir)s
    --shrinkage %(deseq2_shrinkage)s
    --userlist %(pathways_usergenes)s
    --permute %(deseq2_permutations)s
    --pathways %(pathways_GSEA_datasets)s
    > %(outdir)s/deseq2.log;
    '''
    P.run(statement)

@transform(filterEdgeR,
           formatter("DEresults.dir/edger/(?P<QUANTIFIER>\S+)_(?P<DESIGN>\S+)/experiment_out.rds"),
           "DEresults.dir/edger/{QUANTIFIER[0]}_{DESIGN[0]}/results_table.rds",
           "{DESIGN[0]}",
           "{QUANTIFIER[0]}")
def runEdgeR(infile, outfile, design_name, quantifier_name):
    ''' run EdgeR to identify differentially expression transcripts/genes'''

    design = "design" + design_name + ".tsv"
    outdir = os.path.dirname(outfile)
    r_root = os.path.abspath(os.path.dirname(cgatpipelines.__file__))
    scriptpath = os.path.join(r_root, "Rtools/diffexpression.R")

    model = PARAMS.get('edger_model%s' % design_name, None)
    contrast = PARAMS.get('edger_contrast%s' % design_name, None)
    refgroup = PARAMS.get('edger_refgroup%s' % design_name, None)
    coef = PARAMS.get('edger_coef%s' % design_name, None)

    if model is None:
        raise ValueError("edger_model{} is not specified".format(
            design_name))
    if contrast is None:
        raise ValueError("edger_contrast{} is not specified".format(
            design_name))
    if refgroup is None:
        raise ValueError("edger_refgroup{} is not specified".format(
            design_name))
    if coef is None:
        raise ValueError("edger_coef{} is not specified".format(
            design_name))

    statement = '''
    export R_ROOT=%(r_root)s &&
    Rscript %(scriptpath)s
    --rds-filename %(infile)s
    --model %(model)s
    --contrast %(contrast)s
    --refgroup %(refgroup)s
    --coef %(coef)s
    --alpha %(edger_fdr)s
    --outdir %(outdir)s
    --userlist %(pathways_usergenes)s
    --permute %(edger_permutations)s
    --pathways %(pathways_GSEA_datasets)s
    > %(outdir)s/edger.log;
    '''
    P.run(statement)


@mkdir("DEresults.dir/sleuth")
@product(mergeCounts,
         formatter(
             ".*/(?P<QUANTIFIER>(kallisto|salmon)).dir/transcripts.tsv.gz"),
         ["design%s.tsv" % x.asFile() for x in DESIGNS],
         formatter(".*/design(?P<DESIGN>\S+).tsv$"),
         ["DEresults.dir/sleuth/{QUANTIFIER[0][0]}_{DESIGN[1][0]}_transcripts_results.tsv",
          "DEresults.dir/sleuth/{QUANTIFIER[0][0]}_{DESIGN[1][0]}_genes_results.tsv"],
         "{DESIGN[1][0]}",
         "{QUANTIFIER[0][0]}")
def runSleuth(infiles, outfiles, design_name, quantifier):
    ''' run sleuth to identify differentially expression transcripts/genes'''
    ''' NOT FUNCTIONAL'''


    # to estimate sleuth memory, we need to know the number of
    # samples, transcripts and boostraps
    number_transcripts = 0
    with iotools.open_file(transcripts, "r") as inf:
        for line in inf:
            if line.startswith(">"):
                number_transcripts += 1
    
    number_samples = CALCULATENUMBEROFSAMPLESHERE

    job_memory = rnaseq.estimateSleuthMemory(
        PARAMS["%(quantifier)s_bootstrap" % locals()],
        number_samples, number_transcripts)

    statement = '''
    python -m cgatpipelines.tasks.counts2table
    --design-tsv-file=%(design)s
    --output-filename-pattern=%(transcript_prefix)s
    --log=%(transcript_log)s
    --method=sleuth
    --fdr=%(sleuth_fdr)s
    --model=%(model)s
    --contrast=%(contrast)s
    --sleuth-counts-dir=%(quantifier)s.dir
    --reference-group=%(refgroup)s
    --de-test=%(detest)s
    '''
    if detest == "lrt":
        statement += '''
        --reduced-model=%(reduced_model)s
        '''
    statement += '''
    -v 0
    >%(transcript_out)s
    '''

    P.run(statement)



# Define the task for differential expression and normalisation
DETARGETS = []
NORMTARGETS = []
mapToDETargets = {'edger': (runEdgeR, ),
                  'deseq2': (runDESeq2,),
                  'sleuth': (runSleuth,)}

for x in P.as_list(PARAMS["de_tools"]):
    DETARGETS.extend(mapToDETargets[x])


@follows(*DETARGETS)
def differentialExpression():
    ''' dummy task to define upstream differential expression tasks'''


# AH: see below
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(DETARGETS, "differential_expression.load")
def loadDifferentialExpression(infiles, outfiles):
    for infile in iotools.flatten(infiles):
        outfile = P.snip(infile, ".tsv") + ".load"
        P.load(infile, outfile)


###################################################
# target functions for code execution             #
###################################################


@follows(exploratoryAnalysis,
         loadDifferentialExpression)
def full():
    ''' collects DE tasks'''


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
