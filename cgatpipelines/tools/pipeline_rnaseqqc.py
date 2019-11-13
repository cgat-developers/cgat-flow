
"""
====================
RNASeqQC pipeline
====================



Overview
========

This pipeline should be run as the first step in your RNA seq analysis
work flow. It will help detect error and biases within your raw
data. The output of the pipeline can be used to filter out problematic
cells in a standard RNA seq experiment. For single cell RNA seq the
pipeline_rnaseqqc.py should be run instead.

Salmon is used to perform rapid alignment-free transcript
quantification and hisat is used to align a subset of reads to the
reference genome provided by the genesets pipeline.

From the salmon and hisat output, a number of analyses are
performed, either within the pipeline or during the reporting:

- Proportion of reads aligned to annotated features
    (rRNA, protein coding, lincRNA etc)
- Sequencing depth saturation curves Per Sample
- Per-sample expression distributions
- Strandedness assesment
- Assessment of sequence biases
- Expression of top genes and expression of genes of interest

Most of the above analysis will group samples by the sample factors
(see Important configuration options below for details on how factors
are identified)


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning`
on general information how to use cgat pipelines.


Input
-----

Reads are imported by placing files or linking to files in the :term:
`working directory`.

The following suffixes/file types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the :file:
   `fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq.2.gz
   Paired-end reads in fastq format.
   The two fastq files must be sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input files.
   Thus it might be difficult to mix different formats.

Important configuration options
===============================

To determine the experimental factors in your experiment, name files
with factors separated by ``-``, for example::

   sample1-mRNA-10k-R1-L01.fastq.1.gz
   sample1-mRNA-10k-R1-L01.fastq.2.gz
   sample1-mRNA-10k-R1-L02.fastq.1.gz
   sample1-mRNA-10k-R1-L02.fastq.2.gz
   sample1-mRNA-150k-R1-L01.fastq.1.gz
   sample1-mRNA-150k-R1-L01.fastq.2.gz
   sample1-mRNA-150k-R1-L02.fastq.1.gz
   sample1-mRNA-150k-R1-L02.fastq.2.gz

and then set the ``factors`` variable in :file:`pipeline.yml` to::

   factors=experiment-source-replicate-lane

If you want to include additional factors which are not identifiable
from the sample names you can specfify these in an optional file
"additional_factors.tsv". This file must contain the sample names in
the first columns and then an additional column for each factor (tab
separated). See below for an example to include the additional factors
"preparation_date" and "rna_quality":

sample    preparation_date    rna_quality
sample1-mRNA-10k-R1-L01    01-01-2016    poor
sample1-mRNA-10k-R1-L01    01-01-2016    poor
sample1-mRNA-10k-R1-L02    04-01-2016    good
sample1-mRNA-10k-R1-L02    04-01-2016    good


Pipeline output
===============

The major output is a set of HTML pages and plots reporting on the
apparent biases in transcript abudance within the sequence archive
The pipeline also produces output in the database file:`csvdb`.

Example
=======

Example data is available at
http://www.cgat.org/~andreas/sample_data/pipeline_rnaseqqc.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_readqc.tgz
   tar -xvzf pipeline_readqc.tgz
   cd pipeline_readqc
   python <srcdir>/pipeline_readqc.py make full

Requirements:

+---------+------------+------------------------------------------------+
|*Program*|*Version*   |*Purpose*                                       |
+---------+------------+------------------------------------------------+
|salmon   |>=0.9.0     |pseudo alignment                               |
+---------+------------+------------------------------------------------+
|hisat    |>=0.1.6     |read mapping                                   |
+---------+------------+------------------------------------------------+
|samtools |>=0.1.16    |bam/sam files
+---------+------------+------------------------------------------------+
|bedtools |            |work with intervals
+---------+------------+------------------------------------------------+
|picard   |>=1.42      |bam/sam files
+---------+------------+------------------------------------------------+
|bamstats |>=1.22      |from CGR, liverpool
+---------+------------+------------------------------------------------+
|sra-tools|            |extracting sra files
+---------+------------+------------------------------------------------+


Glossary
========

.. glossary::

   hisat
      hisat_- a read mapper used in the pipeline because it is
              relatively quick to run
   salmon
      salmon_-a pseudoaligner that is used for quantifying the
                abundance transcripts
.._hisat: http://ccb.jhu.edu/software/hisat/manual.shtml
.. salmon: https://github.com/

Code
====

"""

###################################################
###################################################
###################################################
# load modules
###################################################

# import ruffus
from ruffus import transform, suffix, regex, merge, \
    follows, mkdir, originate, add_inputs, jobs_limit, split, \
    subdivide, formatter, collate

# import useful standard python modules
import sys
import os
import sqlite3
import re
import pandas as pd
import numpy as np
import itertools
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter

import cgatcore.experiment as E
import cgat.GTF as GTF
import cgatcore.iotools as iotools
import cgatpipelines.tasks.mapping as mapping
import cgatpipelines.tasks.windows as windows
import cgatpipelines.tasks.mappingqc as mappingqc
import cgatpipelines.tasks.rnaseq as rnaseq
from cgatcore import pipeline as P

import json

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


#########################################################################
#########################################################################
#########################################################################
# define input files
SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.fa.gz",
                    "*.sra",
                    "*.export.txt.gz",
                    "*.csfasta.gz",
                    "*.csfasta.F3.gz",
                    )

SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(
    r"(.*\/)*(\S+).(fastq.1.gz|fastq.gz|fa.gz|sra|"
    "csfasta.gz|csfasta.F3.gz|export.txt.gz)")

###################################################################
# Pipeline Utility functions
###################################################################


def findSuffixedFile(prefix, suffixes):
    for check_suffix in suffixes:
        check_infile = prefix + check_suffix
        if os.path.exists(check_infile):
            return (check_infile, check_suffix)

###################################################################
# count number of reads
###################################################################


@follows(mkdir("nreads.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"nreads.dir/\2.nreads")
def countReads(infile, outfile):
    '''Count number of reads in input files.'''
    m = mapping.Counter()
    statement = m.build((infile,), outfile)
    P.run(statement)

###################################################################
# build geneset
###################################################################


@follows(mkdir("geneset.dir"))
@transform(PARAMS["annotations_interface_geneset_coding_exons_gtf"],
           regex("(\S+).gtf.gz"), "geneset.dir/refcoding.junctions")
def buildJunctions(infile, outfile):
    '''build file with splice junctions from gtf file.

    Identify the splice junctions from a gene set :term:`gtf`
    file. A junctions file is a better option than supplying a GTF
    file, as parsing the latter often fails. See:

    http://seqanswers.com/forums/showthread.php?t=7563

    Parameters
    ----------
    infile : str
       Input filename in :term:`gtf` format
    outfile: str
       Output filename

    '''

    outf = iotools.open_file(outfile, "w")
    njunctions = 0
    for gffs in GTF.transcript_iterator(
            GTF.iterator(iotools.open_file(infile, "r"))):

        gffs.sort(key=lambda x: x.start)
        end = gffs[0].end
        for gff in gffs[1:]:
            # subtract one: these are not open/closed coordinates but
            # the 0-based coordinates
            # of first and last residue that are to be kept (i.e., within the
            # exon).
            outf.write("%s\t%i\t%i\t%s\n" %
                       (gff.contig, end - 1, gff.start, gff.strand))
            end = gff.end
            njunctions += 1

    outf.close()

    if njunctions == 0:
        E.warn('no junctions found in gene set')
        return
    else:
        E.info('found %i junctions before removing duplicates' % njunctions)

    # make unique
    statement = '''mv %(outfile)s %(outfile)s.tmp &&
                   cat < %(outfile)s.tmp | sort | uniq > %(outfile)s &&
                   rm -f %(outfile)s.tmp; '''
    P.run(statement)


@transform(PARAMS["annotations_interface_geneset_coding_exons_gtf"],
           regex("(\S+)"),           
           "geneset.dir/refcoding.fasta")
def buildTranscriptFasta(infile, outfile):
    """build geneset where all exons within a gene
    are merged.
    """
    dbname = outfile[:-len(".fasta")]

    statement = '''zcat %(infile)s
    | cgat gff2fasta
    --is-gtf
    --genome=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log
    | cgat index_fasta
    %(dbname)s --force-output -
    > %(dbname)s.log
    '''
    P.run(statement)


@transform(PARAMS["annotations_interface_geneset_coding_exons_gtf"],
            regex("(\S+).gtf.gz"),
           "geneset.dir/transcript2gene.tsv")
def buildTranscriptGeneMap(infile, outfile):
    """build a map of transcript ids to gene ids."""

    statement = """
    zcat %(infile)s
    | cgat gtf2tsv
    --attributes-as-columns
    --output-only-attributes
    | cgat csv_cut transcript_id gene_id
    > %(outfile)s"""
    P.run(statement, job_memory="16G")

###################################################################
# subset fastqs
###################################################################


@follows(mkdir("fastq.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"fastq.dir/\2.subset")
def subsetSequenceData(infile, outfile):
    """subset fastq files"""
    m = mapping.SubsetHead(limit=PARAMS["sample_size"])
    statement = m.build((infile,), outfile)
    P.run(statement,
          ignore_pipe_errors=True,
          ignore_errors=True)

    iotools.touch_file(outfile)


@follows(mkdir("fastq.dir"))
@merge(countReads,
       "fastq.dir/highest_depth_sample.sentinel")
def identifyHighestDepth(infiles, outfile):
    ''' identify the sample with the highest depth'''

    highest_depth = 0
    for count_inf in infiles:
        for line in iotools.open_file(count_inf, "r"):
            if not line.startswith("nreads"):
                continue
            nreads = int(line[:-1].split("\t")[1])
            if nreads > highest_depth:
                highest_depth = nreads
                highest_depth_sample = os.path.basename(
                    P.snip(count_inf, ".nreads"))

    assert highest_depth_sample, ("unable to identify the sample "
                                  "with the highest depth")

    infile, inf_suffix = findSuffixedFile(highest_depth_sample,
                                          [x[1:] for x in SEQUENCESUFFIXES])
    infile = os.path.abspath(infile)

    assert infile, ("unable to find the raw data for the "
                    "sample with the highest depth")

    dst = os.path.abspath(P.snip(outfile, ".sentinel") + inf_suffix)

    def forcesymlink(src, dst):
        try:
            os.symlink(src, dst)
        except:
            os.remove(dst)
            os.symlink(src, dst)

    forcesymlink(infile, dst)

    # if paired end fastq, need to link the paired end too!
    if inf_suffix == ".fastq.1.gz":
        dst2 = P.snip(outfile, ".sentinel") + ".fastq.2.gz"
        forcesymlink(infile.replace(".fastq.1.gz", ".fastq.2.gz"), dst2)

    forcesymlink("%s.nreads" % highest_depth_sample,
                 "nreads.dir/highest_depth_sample.nreads")

    iotools.touch_file(outfile)


@split(identifyHighestDepth,
       "fastq.dir/highest_counts_subset_*")
def subsetRange(infile, outfiles):
    '''subset highest depth sample to 10%-100% depth '''

    outfile = "fastq.dir/highest_counts_subset.sentinel"
    infile_prefix = P.snip(os.path.basename(infile), ".sentinel")
    nreads_inf = "nreads.dir/%s.nreads" % infile_prefix

    for line in iotools.open_file(nreads_inf, "r"):
        if not line.startswith("nreads"):
            continue
        nreads = int(line[:-1].split("\t")[1])

    infile, inf_suffix = findSuffixedFile(P.snip(infile, ".sentinel"),
                                          [x[1:] for x in SEQUENCESUFFIXES])

    # mapping.Counter double counts for paired end
    # Note: this wont handle sra. Need to add a call to Sra.peak to check for
    # paired end files in SRA
    if inf_suffix == ".fastq.1.gz":
        nreads = nreads / 2

    subset_depths = list(range(10, 110, 10))
    limits = [int(nreads / (100.0 / int(depth)))
              for depth in subset_depths]

    m = mapping.SubsetHeads(limits=limits)
    statement = m.build((infile,), outfile)

    P.run(statement,
          ignore_pipe_errors=True,
          ignore_errors=True)

    iotools.touch_file(outfile)


@follows(subsetSequenceData)
def subset():
    pass

###################################################################
# map reads
###################################################################


@follows(mkdir("hisat.dir"))
@transform(subsetSequenceData,
           regex("fastq.dir/(.*).subset"),
           add_inputs(buildJunctions),
           r"hisat.dir/\1.hisat.bam")
def mapReadsWithHisat(infiles, outfile):
    '''
    Map reads using Hisat  (spliced reads).

    Parameters
    ----------
    infiles: list
        contains two filenames -

    infiles[0]: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    infiles[1]: str
        filename with suffix .junctions containing a list of known
        splice junctions.

    hisat_threads: int
        :term:`PARAMS`
        number of threads with which to run hisat

    hisat_memory: str
        :term:`PARAMS`
        memory required for hisat job

    hisat_executable: str
        :term:`PARAMS`
        path to hisat executable

    hisat_library_type: str
        :term:`PARAMS`
        hisat rna-strandess parameter, see
        https://ccb.jhu.edu/software/hisat/manual.shtml#command-line

    hisat_options: str
        options string for hisat, see
        https://ccb.jhu.edu/software/hisat/manual.shtml#command-line

    hisat_index_dir: str
        path to directory containing hisat indices

    strip_sequence: bool
        :term:`PARAMS`
        if set, strip read sequence and quality information

    outfile: str
        :term:`bam` filename to write the mapped reads in bam format.

    .. note::
    If hisat fails with an error such as::

       Error: segment-based junction search failed with err =-6
       what():  std::bad_alloc

    it means that it ran out of memory.

    '''

    job_threads = PARAMS["hisat_threads"]
    job_memory = PARAMS["hisat_memory"]

    m = mapping.Hisat(
        executable=P.substitute_parameters(
            **locals())["hisat_executable"],
        strip_sequence=PARAMS["strip_sequence"])

    infile, junctions = infiles
    infile = P.snip(infile, ".subset") + ".fastq.gz"
    if not os.path.exists(infile):
        infile = P.snip(infile, ".fastq.gz") + ".fastq.1.gz"

    statement = m.build((infile,), outfile)

    P.run(statement)


###################################################################
# build mapping stats
###################################################################


@transform(mapReadsWithHisat,
           regex("(.*)/(.*)\.(.*).bam"),
           r"\1/\2.\3.readstats")
def buildBAMStats(infile, outfile):
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

    job_memory = "16G"

    track = P.snip(os.path.basename(infile), ".hisat.bam")

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
         --num-reads=%(sample_size)i
         --output-filename-pattern=%(outfile)s.%%s
    < %(infile)s
    > %(outfile)s
    '''

    P.run(statement)


@P.add_doc(mappingqc.loadBAMStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildBAMStats, "bam_stats.load")
def loadBAMStats(infiles, outfile):
    ''' load bam statistics into bam_stats table '''
    mappingqc.loadBAMStats(infiles, outfile)


@P.add_doc(windows.summarizeTagsWithinContext)
@transform(mapReadsWithHisat,
           suffix(".bam"),
           add_inputs(
               PARAMS["annotations_interface_genomic_context_bed"]),
           ".contextstats.tsv.gz")
def buildContextStats(infiles, outfile):
    ''' build mapping context stats '''
    windows.summarizeTagsWithinContext(
        infiles[0], infiles[1], outfile)


@P.add_doc(windows.loadSummarizedContextStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@follows(loadBAMStats)
@merge(buildContextStats, "context_stats.load")
def loadContextStats(infiles, outfile):
    ''' load context mapping statistics into context_stats table '''
    windows.loadSummarizedContextStats(infiles, outfile)


@originate("geneset.dir/altcontext.bed.gz")
def buildBedContext(outfile):
    ''' Generate a bed file that can be passed into buildAltContextStats '''

    dbh = P.connect()

    tmp_bed_sorted_filename = P.get_temp_filename(shared=True)

    sql_statements = [
        '''SELECT DISTINCT GTF.contig, GTF.start, GTF.end, "lincRNA"
        FROM gene_info GI
        JOIN geneset_lincrna_exons_gtf GTF
        ON GI.gene_id=GTF.gene_id
        WHERE GI.gene_biotype == "lincRNA"''',
        '''SELECT DISTINCT GTF.contig, GTF.start, GTF.end, "snoRNA"
        FROM gene_info GI
        JOIN geneset_noncoding_exons_gtf GTF
        ON GI.gene_id=GTF.gene_id
        WHERE GI.gene_biotype == "snoRNA"''',
        '''SELECT DISTINCT GTF.contig, GTF.start, GTF.end, "miRNA"
        FROM gene_info GI
        JOIN geneset_noncoding_exons_gtf GTF
        ON GI.gene_id=GTF.gene_id
        WHERE GI.gene_biotype == "miRNA"''',
        '''SELECT DISTINCT GTF.contig, GTF.start, GTF.end, "protein_coding"
        FROM gene_info GI
        JOIN geneset_coding_exons_gtf GTF
        ON GI.gene_id=GTF.gene_id
        WHERE GI.gene_biotype == "protein_coding"''']

    with iotools.open_file(tmp_bed_sorted_filename, "w") as tmp_bed_sorted:
        for sql_statement in sql_statements:
            state = dbh.execute(sql_statement)
            for line in state:
                tmp_bed_sorted.write(("%s\n") % "\t".join(map(str, line)))

    statement = '''
    sort -k1,1 -k2,2n -k3,3n
    < %(tmp_bed_sorted_filename)s
    | bgzip
    > %(outfile)s'''

    P.run(statement)

    os.unlink(tmp_bed_sorted_filename)


@P.add_doc(windows.summarizeTagsWithinContext)
@follows(buildBedContext)
@transform(mapReadsWithHisat,
           suffix(".bam"),
           add_inputs(buildBedContext),
           ".altcontextstats.tsv.gz")
def buildAltContextStats(infiles, outfile):
    ''' build mapping context stats of snoRNA, miRNA,
        lincRNA, protein coding '''

    infile, bed = infiles

    windows.summarizeTagsWithinContext(
        infile, bed,  outfile)


@P.add_doc(windows.loadSummarizedContextStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@follows(loadContextStats)
@merge(buildAltContextStats, "altcontext_stats.load")
def loadAltContextStats(infiles, outfile):
    ''' load context mapping statistics into context_stats table '''
    windows.loadSummarizedContextStats(infiles,
                                               outfile,
                                               suffix=".altcontextstats.tsv.gz")


###################################################################
# alignment-free quantification
###################################################################

@follows(mkdir("salmon.dir"))
@transform(buildTranscriptFasta,
           regex("(\S+)"),
           "salmon.dir/transcripts.salmon.index")
def indexForSalmon(infile, outfile):
    '''create a salmon index'''

    statement = '''
    salmon index -k %(salmon_kmer)i -t %(infile)s -i %(outfile)s >& %(outfile)s.log'''
    P.run(statement, job_memory="8G", job_condaenv="salmon")


@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(indexForSalmon,
                    buildTranscriptGeneMap),
           [r"salmon.dir/\2/lib_format_counts.json",
            r"salmon.dir/\2/transcripts.tsv.gz",
            r"salmon.dir/\2/genes.tsv.gz"])
def runSalmon(infiles, outfiles):
    '''quantify abundance using Salmon'''

    fastqfile = [x[0] for x in infiles]
    index = infiles[0][1]
    transcript2geneMap = infiles[0][2]
    libformat, transcript_outfile, gene_outfile = outfiles

    Quantifier = rnaseq.SalmonQuantifier(
        infile=fastqfile,
        transcript_outfile=transcript_outfile,
        gene_outfile=gene_outfile,
        annotations=index,
        job_threads=PARAMS["salmon_threads"],
        job_memory=PARAMS["salmon_memory"],
        options=PARAMS["salmon_options"],
        bootstrap=1,
        libtype='A',
        kmer=PARAMS['salmon_kmer'],
        transcript2geneMap=transcript2geneMap)

    Quantifier.run_all()


@subdivide(runSalmon,
         regex(".*"),
         ["salmon.dir/transcripts.tsv.gz",
          "salmon.dir/genes.tsv.gz"])
def mergeSalmonResults(infiles, outfiles):
    ''' merge counts for alignment-based methods'''

    transcript_infiles = [x[1] for x in infiles]
    gene_infiles = [x[2] for x in infiles]

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
@collate(mergeSalmonResults,
        regex(r"(.*).tsv.gz"),
        r"\1.load")
def loadSalmonResults(infile, outfile):
    P.load(infile[0], outfile)

###################################################################
# strand bias
###################################################################


@follows(mkdir("geneset.dir"))
@merge(PARAMS["annotations_interface_geneset_all_gtf"],
       "geneset.dir/refflat.txt")
def buildRefFlat(infile, outfile):
    '''build flat geneset for Picard RnaSeqMetrics.'''

    tmpflat = P.get_temp_filename(".")

    statement = '''
    gtfToGenePred -genePredExt -geneNameAsName2 %(infile)s %(tmpflat)s &&
    paste <(cut -f 12 %(tmpflat)s) <(cut -f 1-10 %(tmpflat)s)
    > %(outfile)s
    '''
    P.run(statement)
    os.unlink(tmpflat)


@P.add_doc(mappingqc.buildPicardRnaSeqMetrics)
@transform(mapReadsWithHisat,
           suffix(".bam"),
           add_inputs(buildRefFlat),
           ".picard_rna_metrics")
def buildPicardRnaSeqMetrics(infiles, outfile):
    '''Get duplicate stats from picard RNASeqMetrics '''
    # convert strandness to tophat-style library type
    if PARAMS["hisat_library_type"] == ("RF" or "R"):
        strand = "SECOND_READ_TRANSCRIPTION_STRAND"
    elif PARAMS["hisat_library_type"] == ("FR" or "F"):
        strand = "FIRST_READ_TRANSCRIPTION_STRAND"
    else:
        strand = "NONE"

    mappingqc.buildPicardRnaSeqMetrics(infiles, strand, outfile)


@P.add_doc(mappingqc.loadPicardRnaSeqMetrics)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildPicardRnaSeqMetrics, ["picard_rna_metrics.load",
                                  "picard_rna_histogram.load"])
def loadPicardRnaSeqMetrics(infiles, outfiles):
    '''merge alignment stats into single tables.'''
    mappingqc.loadPicardRnaSeqMetrics(infiles, outfiles)


###################################################################
# saturation analysis
###################################################################

@transform(subsetRange,
           regex("fastq.dir/highest_counts_subset_(\d+)."
                 "(fastq.1.gz|fastq.gz|fa.gz|sra|"
                 "csfasta.gz|csfasta.F3.gz|export.txt.gz)"),
           add_inputs(indexForSalmon),
           r"salmon.dir/highest_counts_subset_\1/quant.sf")
def runSalmonSaturation(infiles, outfile):
    '''quantify abundance of transcripts with increasing subsets of the data'''
    
    infile, index = infiles

    job_threads = PARAMS["salmon_threads"]
    job_memory = PARAMS["salmon_memory"]

    salmon_bootstrap = 20
    salmon_libtype = 'A'
    salmon_options = PARAMS["salmon_options"]

    m = mapping.Salmon()

    statement = m.build((infile,), outfile)

    P.run(statement)


@jobs_limit(1, "R")
@mkdir("salmon.dir/plots.dir")
@merge(runSalmonSaturation,
       "salmon.dir/plots.dir/saturation_plots.sentinel")
def plotSalmonSaturation(infiles, outfile):
    ''' Plot the relationship between sample sequencing depth and
    quantification accuracy'''

    plotfile_base = P.snip(outfile, ".sentinel")
    bootstrap_sat_plotfile = plotfile_base + "_boostrap_cv.png"
    accuracy_sat_plotfile = plotfile_base + "_accuracy.png"

    quant_dir = os.path.dirname(os.path.dirname(infiles[0]))

    # This is currently hardcoded to expect 10 infiles named:
    # (quant.dir)/highest_counts_subset_(n)/quant.sf,
    # where (n) is the subset index (0-9)

    R('''
    library(reshape2)
    library(ggplot2)
    library(Hmisc)

    Path = "%(quant_dir)s"

    # following code to read Salmon binary files borrows from Rob
    # Patro's Wasabi R package for making salmon output
    # compatable with sleuth
    minfo <- rjson::fromJSON(file=file.path(
      Path, 'highest_counts_subset_9', "aux", "meta_info.json"))

    numBoot <- minfo$num_bootstraps

    point_df = read.table(file.path(Path, 'highest_counts_subset_9', "quant.sf"),
                          sep="\t", header=T, row.names=1)

    final_cols = NULL

    for (ix in seq(0,9,1)){
      bootCon <- gzcon(file(file.path(
        Path, paste0('highest_counts_subset_', ix), 'aux',
                     'bootstrap', 'bootstraps.gz'), "rb"))

      # read in binary data
      boots <- readBin(bootCon, "double",
                       n = minfo$num_targets * minfo$num_bootstraps)
      close(bootCon)

      # turn data into dataframe
      boots_df = t(data.frame(matrix(unlist(boots),
                              nrow=minfo$num_bootstraps, byrow=T)))

      # add rownames
      rownames(boots_df) = rownames(point_df)

      final_cols[[paste0("sample_", ix)]] = apply(boots_df, 1,
                                                  function(x) sd(x)/mean(x))
    }

    # make final dataframe with boostrap CVs
    final_df = data.frame(do.call("cbind", final_cols))

    # add expression values, subset to transcripts with >1 read and bin exp
    final_df$max_exp = point_df$NumReads
    final_df = final_df[final_df$max_exp>1,]
    final_df$max_exp = as.numeric(cut2(final_df$max_exp, g=10))

    # melt and aggregate
    melted_df = melt(final_df, id="max_exp")
    melted_df = melted_df[is.finite(melted_df$value),]
    aggdata <-aggregate(melted_df$value,
                        by=list(melted_df$max_exp, melted_df$variable),
                        FUN=mean)
    aggdata$Group.1 = as.factor(aggdata$Group.1)

    m_txt = element_text(size=20)
    my_theme = theme(
    axis.text=m_txt,
    axis.title=m_txt,
    legend.text=m_txt,
    legend.title=m_txt,
    aspect.ratio=1)

    p = ggplot(aggdata, aes(10*as.numeric(Group.2), x,
                            colour=Group.1, group=Group.1)) +
    geom_line() +
    theme_bw() +
    xlab("Sampling depth (%%)") +
    ylab("Average Coefficient of variance") +
    scale_colour_manual(name="Exp. Decile",
                        values=colorRampPalette(c("yellow","purple"))(10)) +
    scale_x_continuous(breaks=seq(10,100,10), limits=c(10,100)) +
    my_theme

    ggsave("%(bootstrap_sat_plotfile)s")

    # read in the point estimate data

    tpm_est = NULL

    ref_point_df = read.table(
      file.path(Path, 'highest_counts_subset_9', "quant.sf"),
      sep="\t", header=T, row.names=1)

    for (ix in seq(0,9,1)){

    point_df = read.table(
      file.path(Path, paste0('highest_counts_subset_', ix), "quant.sf"),
      sep="\t", header=T, row.names=1)
`
    tpm_est[[paste0("sample_", ix)]] = (
      abs(point_df$TPM - ref_point_df$TPM) / ref_point_df$TPM)
    }

    tpm_est_df = data.frame(do.call("cbind", tpm_est))

    # add expression values, subset to transcripts with >1 read and bin exp.
    tpm_est_df$max_exp = point_df$NumReads
    tpm_est_df = tpm_est_df[point_df$NumReads>1,]
    tpm_est_df$max_exp = as.numeric(cut2(tpm_est_df$max_exp, g=10))

    # melt and aggregate
    melted_df = melt(tpm_est_df, id="max_exp")
    melted_df = melted_df[is.finite(melted_df$value),]
    aggdata <-aggregate(melted_df$value,
                        by=list(melted_df$max_exp, melted_df$variable),
                        FUN=mean)
    aggdata$Group.1 = as.factor(aggdata$Group.1)

    p = ggplot(aggdata, aes(10*as.numeric(Group.2), x,
                            colour=Group.1, group=Group.1)) +
    geom_line() +
    theme_bw() +
    xlab("Sampling depth (%%)") +
    ylab("Abs. difference in exp. estimate (normalised)") +
    scale_colour_manual(name="Exp. Decile",
                        values=colorRampPalette(c("yellow","purple"))(10)) +
    scale_x_continuous(breaks=seq(10,90,10), limits=c(10,90)) +
    my_theme

    ggsave("%(accuracy_sat_plotfile)s")

    ''' % locals())

    iotools.touch_file(outfile)


###################################################################
# gene coverage profiles
###################################################################


@follows(mkdir("transcriptprofiles.dir"))
@transform(mapReadsWithHisat,
           regex("hisat.dir/(\S+).hisat.bam"),
           add_inputs(PARAMS["annotations_interface_geneset_coding_exons_gtf"]),
           r"transcriptprofiles.dir/\1.transcriptprofile.gz")
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
    --method=geneprofile
    --method=geneprofileabsolutedistancefromthreeprimeend
    --normalize-profile=all
    %(bamfile)s %(gtffile)s
    | gzip
    > %(outfile)s
    '''

    P.run(statement)


@merge(buildTranscriptProfiles,
       "transcriptprofiles.dir/threeprimebiasprofiles.load")
def loadTranscriptProfiles(infiles, outfile):
    ''' concatenate and load the transcript profiles
    Retain sample name as column = "track'''

    regex = ("transcriptprofiles.dir/(\S+).transcriptprofile.gz."
             "geneprofileabsolutedistancefromthreeprimeend.matrix.tsv.gz")

    infiles = [
        x + ".geneprofileabsolutedistancefromthreeprimeend.matrix.tsv.gz" for x in infiles]

    P.concatenate_and_load(infiles, outfile, regex_filename=regex)


@merge(SEQUENCEFILES,
       "experiment.tsv")
def buildExperimentTable(infiles, outfile):

    d = os.getcwd()
    # TODO: read from config file
    project_id = "unknown"
    with iotools.open_file(outfile, "w") as outf:
        outf.write("id\tname\tproject_id\tdirectory\ttitle\n")
        outf.write("\t".join(
            ("1",
             "unknown",
             project_id,
             d,
             PARAMS.get("title", ""))) + "\n")


@merge(SEQUENCEFILES,
       "samples.tsv")
def buildSamplesTable(infiles, outfile):

    with iotools.open_file(outfile, "w") as outf:
        outf.write("id\texperiment_id\tsample_name\n")

        for sample_id, filename in enumerate(sorted(infiles)):
            sample_name, suffix = os.path.basename(filename).split(".", 1)
            outf.write("\t".join(
                (str(sample_id + 1), "1", sample_name)) + "\n")


@merge(SEQUENCEFILES,
       "factors.tsv")
def buildFactorTable(infiles, outfile):

    if "factors" not in PARAMS:
        raise ValueError("factors not defined in config file")

    factor_names = PARAMS.get("factors")
    if factor_names is None or factor_names == "!?":
        raise ValueError("factors not defined in config file")
    factor_names = factor_names.split("-")

    sampleID2sampleName = {}

    with iotools.open_file(outfile, "w") as outf:
        outf.write("sample_id\tfactor\tfactor_value\n")

        for sample_id, filename in enumerate(sorted(infiles)):

            sample_name, suffix = os.path.basename(filename).split(".", 1)
            sampleID2sampleName[sample_name] = sample_id + 1

            parts = sample_name.split("-")

            if len(parts) != len(factor_names):
                raise ValueError(
                    "unexpected number of factors in sample {}: "
                    "expected={}, got={}".format(
                        filename, factor_names, parts))

            for factor, factor_value in zip(factor_names, parts):
                if factor == "_":
                    continue
                outf.write("\t".join((str(sample_id + 1),
                                      factor, factor_value)) + "\n")
            outf.write("\t".join((str(sample_id + 1), "genome",
                                  PARAMS["genome"])) + "\n")

        if os.path.exists("additional_factors.tsv"):
            with iotools.open_file("additional_factors.tsv", "r") as inf:
                header = next(inf)
                header = header.strip().split("\t")
                additional_factors = header[1:]
                for line in inf:
                    line = line.strip().split("\t")
                    sample_name = line[0]
                    factors_values = line[1:]
                    for factor_ix in range(0, len(additional_factors)):
                        try:
                            outf.write("\t".join((
                                str(sampleID2sampleName[sample_name]),
                                additional_factors[factor_ix],
                                factors_values[factor_ix])) + "\n")
                        except KeyError as ke:
                            sample_names = [os.path.basename(x).split(".")[0]
                                            for x in infiles]
                            raise KeyError(
                                "Sample name in additional_factors table does "
                                " not match up with sample names from raw "
                                "infiles: %s not in %s" % (
                                    ke, ",".join(sample_names)))


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform((buildExperimentTable, buildSamplesTable, buildFactorTable),
           suffix(".tsv"),
           ".load")
def loadMetaInformation(infile, outfile):
    P.load(infile, outfile,
           options="--map=id:int "
           "--map=sample_id:int "
           "--map=experiment_id:int "
           "--add-index=id "
           "--add-index=experiment_id "
           "--add-index=sample_id ")


@transform(buildTranscriptFasta,
           suffix("refcoding.fasta"),
           "transcripts_attributes.tsv.gz")
def characteriseTranscripts(infile, outfile):
    ''' obtain attributes for transcripts '''

    statement = '''
    cat %(infile)s | cgat fasta2table
    --split-fasta-identifier --section=na,dn,length -v 0
    | gzip > %(outfile)s'''

    P.run(statement)


@transform(characteriseTranscripts,
           regex("transcripts_attributes.tsv.gz"),
           add_inputs(mergeSalmonResults),
           "bias_binned_means.tsv")
def summariseBias(infiles, outfile):

    def percentage(x):
        return float(x[0]) / float(x[1])

    def lin_reg_grad(x, y):
        slope, intercept, r, p, stderr = linregress(x, y)
        return slope

    attributes, genes, transcripts = infiles

    atr = pd.read_table(attributes, sep='\t', index_col="id")
    atr = atr.rename(columns={'pGC': 'GC_Content'})

    for di in itertools.product("ATCG", repeat=2):
        di = di[0] + di[1]
        temp_df = atr.loc[:, [di, "length"]]
        atr[di] = temp_df.apply(percentage, axis=1)

    drop_cols = (["nAT", "nGC", "pAT", "pA", "pG", "pC", "pT", "nA",
                  "nG", "nC", "nT", "ncodons",
                  "mCountsOthers", "nUnk", "nN", "pN"])
    atr = atr.drop(drop_cols, axis=1)

    atr["length"] = np.log2(atr["length"])

    E.info("loading transcripts from {}".format(transcripts))
    exp = pd.read_csv(transcripts, sep='\t', index_col="transcript_id")
    exp['LogTPM'] = np.log2(exp['TPM'] + 0.1)

    merged = atr.join(exp[['sample_id', 'LogTPM']])

    def norm(array):
        array_min = array.min()
        array_max = array.max()
        return pd.Series([(x - array_min) / (array_max - array_min) for x in array])

    def bin2floats(qcut_bin):
        return [qcut_bin.left, qcut_bin.right]

    def aggregate_by_factor(df, attribute, sample_names, bins, function):

        temp_dict = dict.fromkeys(sample_names, function)

        temp_dict[attribute] = function
        means_df = df[["LogTPM", "sample_id"]].groupby(
            ["sample_id", pd.qcut(df.ix[:, attribute], bins)])

        means_df = pd.DataFrame(means_df.agg(function))
        means_df.reset_index(inplace=True)

        atr_values = means_df[attribute]
        means_df.drop(attribute, axis=1, inplace=True)

        means_df["LogTPM_norm"] = list(
            means_df.groupby("sample_id")["LogTPM"].apply(norm))

        means_df[attribute] = [np.mean(bin2floats(x)) for x in atr_values]
        means_df = pd.melt(means_df, id_vars=[attribute, "sample_id"])
        means_df.columns = ["bin", "sample_id", "variable", "value"]
        means_df["bias_factor"] = [attribute, ] * len(means_df)

        return means_df

    means_binned_df = pd.DataFrame()

    samples = set(exp.index)
    factors = atr.columns.tolist()

    for factor in factors:
        tmp_df = aggregate_by_factor(
            merged, factor, samples, PARAMS["bias_bin"], np.mean)

        means_binned_df = pd.concat([means_binned_df, tmp_df], axis=0)

    means_binned_df.to_csv(outfile, sep="\t",
                           index=False, float_format='%.6f')


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(summariseBias,
           suffix(".tsv"),
           ".load")
def loadBias(infile, outfile):
    P.load(infile, outfile, options="--add-index=sample_id")


###################################################################
# top genes
###################################################################


@mkdir("salmon.dir/plots.dir")
@follows(loadSalmonResults, loadMetaInformation)
@originate("salmon.dir/plots.dir/top_expressed.sentinel")
def plotTopGenesHeatmap(outfile):
    '''extract the top 1000 genes (by expression) for each sample and
    plot a heatmap of the intersection'''

    # if someone can find a nice heatmap plotter from a dissimilarity
    # matrix which is compatable with cgatReport, the sqlite and
    # pandas code should be changed into a tracker

    exp_select_cmd = '''
    SELECT * FROM genes
    '''

    dbh = P.connect()

    exp_df = pd.read_sql(exp_select_cmd, dbh)
    exp_df = exp_df.set_index("id")

    factors_select_cmd = '''
    SELECT factor, factor_value, sample_name
    FROM samples AS A
    JOIN factors AS B
    ON A.id = B.sample_id
    '''

    top_n = 1000

    factors_df = pd.read_sql(factors_select_cmd, dbh)
   
    # extract the top genes per sample
    top_genes = {}
    for col in exp_df.columns:
        top_genes[col] = exp_df[col].sort_values(
            ascending=False)[0:top_n].index

    # set up the empty df
    intersection_df = pd.DataFrame(
        index=range(0, len(exp_df.columns) **
                    2 - len(exp_df.columns)),
        columns=["sample1", "sample2", "intersection", "fraction"])

    # populate the df
    n = 0
    for col1, col2 in itertools.combinations_with_replacement(exp_df.columns, 2):
        s1_genes = top_genes[col1]
        s2_genes = top_genes[col2]
        intersection = set(s1_genes).intersection(set(s2_genes))
        fraction = len(intersection) / float(top_n)
        intersection_df.ix[n] = [col1, col2, len(intersection), fraction]
        n += 1

        # if the samples are different, calculate the reverse intersection too
        if col1 != col2:
            intersection_df.ix[n] = [col2, col1,
                                     len(intersection), fraction]
            n += 1


    # pivot to format for heatmap.3 plotting
    intersection_df['fraction'] = intersection_df['fraction'].astype('float')
    intersection_pivot = pd.pivot_table(
        intersection_df, index="sample1", columns="sample2", values="fraction")

    for factor in set(factors_df['factor'].tolist()):
        # don't want to plot coloured by genome
        if factor == "genome":
            continue

        plotfile = "%s_%s.png" % (P.snip(outfile, ".sentinel"), factor)

        plotHeatmap = R('''
        function(int_df, fact_df){

        library(GMD)
        library(RColorBrewer)

        # subset fact_df to required factor and
        # refactor to remove unwanted levels
        fact_df = fact_df[fact_df$factor=="%(factor)s",]
        rownames(fact_df) = fact_df$sample_name
        fact_df$factor_value = factor(fact_df$factor_value)

        # set up heatmap side colours
        colours = colorRampPalette(
          brewer.pal(length(levels(fact_df$factor_value)),"Dark2"))(
            length(levels(fact_df$factor_value)))
        side_colours = colours[as.numeric((fact_df$factor_value))]

        # plot
        png("%(plotfile)s", width=1000, heigh=1000)
        heatmap.3(as.dist(1- as.matrix(int_df)),
                  Rowv=FALSE, Colv=FALSE,
                  ColIndividualColors = side_colours,
                  RowIndividualColors = side_colours,
                  breaks=100, main="%(factor)s")
        dev.off()
        }
        ''' % locals())

        with localconverter(ro.default_converter + pandas2ri.converter):
            r_intersection_pivot = ro.conversion.py2ri(intersection_pivot)
            r_factors_df = ro.conversion.py2ri(factors_df)
        
        plotHeatmap(r_intersection_pivot,
                    r_factors_df)

    iotools.touch_file(outfile)


###################################################################
# Plot expression distribution
###################################################################

@mkdir("salmon.dir/plots.dir")
@follows(loadMetaInformation,
         loadSalmonResults)
@originate("salmon.dir/plots.dir/expression_distribution.sentinel")
def plotExpression(outfile):
    "Plot the per sample expression distibutions coloured by factor levels"

    # Note: this was being done within the pipeline but the size of
    # the dataframe seemed to be causing errors:
    # "Data values must be of type string or None."
    # See RnaseqqcReport.ExpressionDistribution tracker

    dbh = P.connect()

    statement = """
    SELECT * FROM transcripts"""

    df = pd.read_sql(statement, dbh)
    df = pd.melt(df, id_vars="id", value_name="TPM", var_name="sample_id")
    
    df['logTPM'] = df['TPM'].apply(lambda x: np.log2(x + 1))
    
    factors = dbh.execute("SELECT DISTINCT factor FROM factors")
    factors = [x[0] for x in factors if x[0] != "genome"]

    for factor in factors:

        plotfile = P.snip(outfile, ".sentinel") + "_%s.png" % factor

        factor_statement = '''
        select *
        FROM factors
        JOIN samples
        ON factors.sample_id = samples.id
        WHERE factor = "%(factor)s"''' % locals()

        factor_df = pd.read_sql(factor_statement, dbh)
        full_df = pd.merge(df, factor_df, left_on="sample_id",
                           right_on="sample_name")

        plotDistribution = R('''
        function(df){

        library(ggplot2)
        p = ggplot(df, aes(x=logTPM, group=sample_name,
                           colour=as.factor(factor_value))) +
        geom_density() +
        xlab("Log2(TPM)") + ylab("Density") +
        scale_colour_discrete(name="Factor Level") +
        theme_bw() +
        ggtitle("%(factor)s")

        ggsave("%(plotfile)s")
        }
        ''' % locals())
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_full_df = ro.conversion.py2ri(full_df)
        
        plotDistribution(r_full_df)

    iotools.touch_file(outfile)

###################################################################
# Run Salmon To Autodetect Strandedness
###################################################################

@merge(runSalmon, "strandedness.tsv")
def checkStrandednessSalmon(infiles, outfile):
    '''
    Read the output from salmon used to determine strandedness
    and write a table containing the number of alignments
    consistent with each type of library.
    The possible types are described here:
    http://salmon.readthedocs.io/en/latest/library_type.html
    '''
    results = pd.DataFrame()
    strandfiles = [x[0] for x in infiles]
    for strandfile in strandfiles:
        j = json.load(open(strandfile, "r"))
        vals = list(j.values())
        cols = list(j.keys())
        D = pd.DataFrame(vals, index=cols).T
        D['sample'] = strandfile.split("/")[-2]
        results = results.append(D)
    results = results[["sample", "expected_format",
                       "compatible_fragment_ratio",
                       "num_compatible_fragments",
                       "num_assigned_fragments",
                       "num_frags_with_consistent_mappings",
                       "num_frags_with_inconsistent_or_orphan_mappings",
                       "MSF", "OSF", "ISF", "MSR",
                       "OSR", "ISR", "SF", "SR",
                       "MU", "OU", "IU", "U"]]

    results.to_csv(outfile, sep="\t", index=None)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(checkStrandednessSalmon,
           suffix(".tsv"),
           ".load")
def loadStrandednessSalmon(infile, outfile):
    P.load(infile, outfile)


@transform(checkStrandednessSalmon, suffix(".tsv"), ".png")
def plotStrandednessSalmon(infile, outfile):
    '''
    Plots a bar plot of the salmon strandness estimates
    as counts per sample.
    '''
    sns.set_style('ticks')
    tab = pd.read_csv(infile, sep="\t")
    counttab = tab[tab.columns[7:]]
    f = plt.figure(figsize=(10, 7))
    a = f.add_axes([0.1, 0.1, 0.6, 0.75])
    x = 0
    colors = sns.color_palette("Dark2", 10)
    a.set_ylim(0, max(counttab.values[0]) + max(counttab.values[0]) * 0.1)
    for item in counttab.columns:
        a.bar(range(x, x + len(tab)), tab[item], color=colors)
        x += len(tab)
    a.ticklabel_format(style='plain')
    a.vlines(np.arange(-0.4, a.get_xlim()[1], len(tab)),
             a.get_ylim()[0], a.get_ylim()[1], lw=0.5)
    a.set_xticks(np.arange(0 + len(tab) / 2, a.get_xlim()[1],
                           len(tab)))
    a.set_xticklabels(counttab.columns)
    sns.despine()
    patches = []
    for c in colors[0:len(tab)]:
        patches.append(mpatches.Patch(color=c))
    l = f.legend(labels=list(tab['sample']), handles=patches, loc=1)
    f.suptitle('Strandedness Estimates')
    f.savefig(outfile)

###################################################################
# Main pipeline tasks
###################################################################


@follows(loadContextStats,
         loadBAMStats,
         loadTranscriptProfiles,
         loadSalmonResults,
         loadMetaInformation,
         loadBias,
         loadPicardRnaSeqMetrics,
         loadAltContextStats,
         plotSalmonSaturation,
         plotTopGenesHeatmap,
         plotExpression,
         plotStrandednessSalmon,
         loadStrandednessSalmon)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
