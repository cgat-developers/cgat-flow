##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""===========================
Pipeline transcriptdiffexpression
===========================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python


Overview
========
RNA-Seq differential expression analysis can, broadly speaking, be
performed at two levels. Gene-level and transcript-level.

As transcripts are the true unit of expression, differential
expression at the transcript-level is more ideal. However,
quantification of transcript-level expression is complicated by reads
which align with multiple transcripts from the same gene, especially
with short read techonologies.  In addition transcript-level
quantification may be hindered by inadequate genome annotation.

The Kallisto and Sleuth tools from Lior Pachter's group are designed
for rapid transcript differential expression analysis

Kallisto estimates transcript expression using a "psuedoalignment"
proceedure and bootstraps the transcript expression estimates. Sleuth
uses these bootstrap estiamtes to estimate the transcript-wise
techincal variance which is subtracted from the total variance, thus
leaving an estimate of the remaining biological variance. Sleuth then
allows the user to fit a transcript-wise general linear model to the
expression data to identify transcripts which are signficantly
differentially expressed.

Kallisto requires that the user has prior knowledge of the transcripts
which may be present in the sample. Here, we use the results of the
annotation pipeline to generate the transcripts input file.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_transcriptdiffexpression.py config

Input files
-----------

Sequence read input. Can be fastq or sra, single or paired end.

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

If further variates need to be given, e.g the General linear model is
defined as ~ group + replicate, these can be specified in further columns:

sample    include    group    pair    replicate
WT-1-1    1    WT    0    1
WT-1-2    1    WT    0    2
Mutant-1-1    1    Mutant    0    1
Mutant-1-2    1    Mutant    0    2

For each design file, the pipeline.ini must specify a model and contrasts


Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

.. Add any additional external requirements such as 3rd party software
   or R modules below:

Requirements:

* kallisto >= 0.42.1
* sleuth >= 0.27.1

Pipeline output
===============

The main outputs of the pipeline are results tables and plots from the
differential expression analysis. Outputs are generated for each
*.design.tsv file and each contrast specfied and placed in DEresults.dir

  `results_[design]_counts.tsv"`
    counts table for all samples within the design

  `results_[design]_tpm.tsv"`
    Transcripts Per Million (tpm) table for all samples within the design

  `results_[design]_[contrast]_sleuth_ma.png`
    MA plot using sleuth function

DEresults.dir contains further plots summarising the sleuth analysis

  `results_[design]_[contrast]_sleuth_vars.png`
    technical vs. observed variance plot from sleuth

  `results_[design]_[contrast]_sleuth_mean_var.png`
    mean-variance plot from sleuth

  `results_[design]_[contrast]_MA_plot.png`
    MA plot from sleuth results table
    (for direction comparison with other methods)

  `results_[design]_[contrast]_volcano_plot.png`
    volcano plot from sleuth results table
    (for direction comparison with other methods)

The summary_plots directory contains further plots summarising the
expression estimates across the samples

Glossary
========

.. glossary::


Code
====

"""

# To do:

# Once Kallisto is upgraded > 0.42.2, include alignment stats from parsing sam

# Once sleuth is capable of performing within gene comparisons of
# transcript expression, add this analysis here too

# Add power test using counts2counts.py?

# Allow user to supply transcripts as a gtf/multifasta rather than
# using annotations pipeline results

from ruffus import *

import sys
import os
import sqlite3
import glob
import subprocess
import pandas as pd
import numpy as np

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Counts as Counts
import CGAT.Expression as Expression

import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineTracks as PipelineTracks

import PipelineTranscriptDiffExpression as TranscriptDiffExpression

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

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

# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.
PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))


# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.

# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.sra")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                      for suffix_name in SEQUENCESUFFIXES])

# enable multiple fastqs from the same sample to be analysed together
if "merge_pattern_input" in PARAMS and PARAMS["merge_pattern_input"]:
    SEQUENCEFILES_REGEX = regex(
        r"%s/%s\.(fastq.1.gz|fastq.gz|sra)" % (
            DATADIR, PARAMS["merge_pattern_input"].strip()))
    # the last expression counts number of groups in pattern_input
    SEQUENCEFILES_OUTPUT = r"quant.dir/%s/abundance.h5" % (
        PARAMS["merge_pattern_output"].strip())

else:
    SEQUENCEFILES_REGEX = regex(
        r".*/(\S+).(fastq.1.gz|fastq.gz|sra)")
    SEQUENCEFILES_OUTPUT = r"quant.dir/\1/abundance.h5"

Sample = PipelineTracks.AutoSample
DESIGNS = PipelineTracks.Tracks(Sample).loadFromDirectory(
    glob.glob("*.design.tsv"), "(\S+).design.tsv")

GENESET = glob.glob("*.gtf.gz")

###############################################################################
# Create kallisto index
###############################################################################


if PARAMS["geneset_auto_generate"]:

    # TS:
    # to do: should we enable filtering by the transcript_biotype rather
    # than gene_biotype?
    @originate("index.dir/gene_ids.tsv")
    def identifyGenes(outfile):
        '''output a list of gene identifiers where biotype matches filter'''

        dbh = connect()

        table = os.path.basename(PARAMS["annotations_interface_table_gene_info"])

        where_cmd = " OR ".join(["WHERE gene_biotype = '%s'" % x for x in
                                 PARAMS["geneset_biotypes"].split(",")])

        select = dbh.execute("""
        SELECT DISTINCT gene_id
        FROM annotations.%(table)s
        %(where_cmd)s""" % locals())

        with IOTools.openFile(outfile, "w") as outf:
            outf.write("gene_id\n")
            outf.write("\n".join((x[0] for x in select)) + "\n")

    @transform(identifyGenes,
               regex("index.dir/gene_ids.tsv"),
               "index.dir/transcripts.gtf.gz")
    def buildGeneSet(mapfile, outfile):
        ''' build a gene set with only transcripts from genes which
        pass filter '''

        geneset = PARAMS['annotations_interface_geneset_all_gtf']

        statement = '''
        zcat %(geneset)s
        | python %(scriptsdir)s/gtf2gtf.py
        --method=filter
        --filter-method=gene
        --map-tsv-file=%(mapfile)s
        --log=%(outfile)s.log
        | gzip
        > %(outfile)s
        '''
        P.run()

else:
    # if a reference gtf is provided, just soft link to this
    assert len(GENESET) > 0, ("if not auto generating a geneset, must"
                              "provide a geneset in a *.gtf.gz file")

    @files(GENESET[0], "index.dir/transcripts.gtf.gz")
    def buildGeneSet(infile, outfile):
        ''' link to the geneset provided'''
        P.clone(os.path.abspath(infile), os.path.abspath(outfile))


@transform(buildGeneSet,
           suffix(".gtf.gz"),
           ".fa")
def buildReferenceTranscriptome(infile, outfile):
    ''' build reference transcriptome from geneset'''

    gtf_file = P.snip(infile, ".gz")

    genome_file = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))

    # sed statement replaces e.g ">1 transcript_id" with ">transcript_id"
    statement = '''
    zcat %(infile)s
    | awk '$3 ~ /exon|UTR/' > %(gtf_file)s;
    gtf_to_fasta %(gtf_file)s %(genome_file)s %(outfile)s;
    sed -i 's/>[0-9]\+ />/g' %(outfile)s;
    checkpoint;
    samtools faidx %(outfile)s
    '''
    P.run()


@transform(buildReferenceTranscriptome,
           suffix(".fa"),
           ".kallisto.index")
def buildKallistoIndex(infile, outfile):
    ''' build a kallisto index'''

    statement = '''
    kallisto index -i %(outfile)s -k %(kallisto_kmer)s %(infile)s
    '''

    P.run()


@follows(mkdir("index.dir"),
         buildKallistoIndex)
def index():
    pass


###############################################################################
# Simulation
###############################################################################

@transform(buildReferenceTranscriptome,
           suffix(".fa"),
           "_kmers.tsv",
           output_dir="simulation.dir")
def countKmers(infile, outfile):
    ''' count the number of unique and non-unique kmers per transcript '''

    job_memory = PARAMS["simulation_kmer_memory"]

    statement = '''
    python %(scriptsdir)s/fasta2unique_kmers.py --input-fasta=%(infile)s
    --kmer-size=%(kallisto_kmer)s -L %(outfile)s.log > %(outfile)s '''

    P.run()


@follows(buildReferenceTranscriptome)
@files([("index.dir/transcripts.fa",
         ("simulation.dir/simulated_reads_%i.fastq.1.gz" % x,
          "simulation.dir/simulated_read_counts_%i.tsv" % x))
        for x in range(0, PARAMS["simulation_iterations"])])
def simulateRNASeqReads(infile, outfiles):
    ''' simulate RNA-Seq reads from the transcripts fasta file '''

    outfile, outfile_counts = outfiles
    outfile2 = outfile.replace(".1.gz", ".2.gz")

    if PARAMS["simulation_paired"]:
        options = "--output-paired-end"
    else:
        options = ""

    statement = '''
    cat %(infile)s |
    python %(scriptsdir)s/fasta2fastq.py
    --output-read-length=%(simulation_read_length)s
    --insert-length-mean=%(simulation_insert_mean)s
    --insert-length-sd=%(simulation_insert_sd)s
    --reads-per-entry-min=%(simulation_min_reads_per_transcript)s
    --reads-per-entry-max=%(simulation_max_reads_per_transcript)s
    --sequence-error-phred=%(simulation_phred)s
    --output-fastq2=%(outfile2)s
    --output-counts=%(outfile_counts)s
    --output-quality-format=33 -L %(outfile)s.log
    %(options)s |
    gzip > %(outfile)s'''

    job_memory = "1G"

    P.run()


@mkdir("simulation.dir/quant.dir")
@transform(simulateRNASeqReads,
           regex("simulation.dir/simulated_reads_(\d+).fastq.1.gz"),
           add_inputs(buildKallistoIndex),
           r"simulation.dir/quant.dir/simulated_reads_\1/abundance.h5")
def runKallistoSimulation(infiles, outfile):
    ''' quantify trancript abundance from simulated reads with kallisto'''

    # TS more elegant way to parse infiles and index?
    infiles, index = infiles
    infile, counts = infiles

    # multithreading not supported until > v0.42.1
    # job_threads = PARAMS["kallisto_threads"]
    job_threads = 1
    job_memory = "6G"

    kallisto_options = PARAMS["kallisto_options"]

    # single bootstrap should be fine for our purposes
    bootstrap = 1

    m = PipelineMapping.Kallisto()
    statement = m.build((infile,), outfile)

    P.run()


@transform(runKallistoSimulation,
           suffix(".h5"),
           ".txt")
def extractKallistoCounts(infile, outfile):
    ''' run kalliso h5dump to extract txt file'''

    outfile_dir = os.path.dirname(os.path.abspath(infile))

    statement = '''kallisto h5dump -o %(outfile_dir)s %(infile)s'''

    P.run()


@transform(extractKallistoCounts,
           regex("simulation.dir/quant.dir/simulated_reads_(\d+)/abundance.txt"),
           r"simulation.dir/quant.dir/simulated_reads_\1/results.tsv",
           r"simulation.dir/simulated_read_counts_\1.tsv")
def mergeAbundanceCounts(infile, outfile, counts):
    ''' merge the kallisto abundance and simulation counts files for
    each simulation '''

    df_abund = pd.read_table(infile, sep="\t", index_col=0)
    df_counts = pd.read_table(counts, sep="\t", index_col=0)

    df_merge = pd.merge(df_abund, df_counts, left_index=True, right_index=True)
    df_merge.index.name = "id"
    df_merge.to_csv(outfile, sep="\t")


@merge(mergeAbundanceCounts,
       "simulation.dir/simulation_results.tsv")
def concatSimulationResults(infiles, outfile):
    ''' concatenate all simulation results '''

    df = pd.DataFrame()

    for inf in infiles:
        df_tmp = pd.read_table(inf, sep="\t")
        df = pd.concat([df, df_tmp], ignore_index=True)

    df.to_csv(outfile, sep="\t", index=False)


@merge((concatSimulationResults,
        countKmers),
       "simulation.dir/simulation_correlations.tsv")
def calculateCorrelations(infiles, outfile):
    ''' calculate correlation across simulation iterations per transcript'''

    abund, kmers = infiles

    df_abund = pd.read_table(abund, sep="\t", index_col=0)
    df_kmer = pd.read_table(kmers, sep="\t", index_col=0)

    # this is hacky, it's doing all against all correlations for the
    # two columns and subsetting
    df_agg = df_abund.groupby(level=0)[[
        "est_counts", "read_count"]].corr().ix[0::2, 'read_count']

    # drop the "read_count" level, make into dataframe and rename column
    df_agg.index = df_agg.index.droplevel(1)
    df_agg = pd.DataFrame(df_agg)
    df_agg.columns = ["cor"]

    # merge and bin the unique fraction values
    df_final = pd.merge(df_kmer, df_agg, left_index=True, right_index=True)
    df_final['fraction_bin'] = (
        np.digitize(df_final["fraction_unique"]*100, bins=range(0, 100, 1),
                    right=True))/100.0

    df_abund_sum = df_abund.groupby(level=0)["est_counts", "read_count"].sum()
    df_final = pd.merge(df_final, df_abund_sum,
                        left_index=True, right_index=True)
    df_final['log2diff'] = np.log2(df_final['est_counts'] /
                                   df_final['read_count'])

    df_final['log2diff_thres'] = [x if abs(x) < 1 else x/abs(x)
                                  for x in df_final['log2diff']]

    df_final.to_csv(outfile, sep="\t", index=True)


@transform(calculateCorrelations,
           suffix(".tsv"),
           ".load")
def loadCorrelation(infile, outfile):
    ''' load the correlations data table'''
    P.load(infile, outfile)


@transform(calculateCorrelations,
           regex("simulation.dir/simulation_correlations.tsv"),
           "simulation.dir/flagged_transcripts.tsv")
def identifyLowConfidenceTranscript(infile, outfile):
    '''
    identify the transcripts which cannot be confidently quantified
    these fall into two categories:

    1. Transcripts whose with poor accuracy of estimated counts

       - transcripts with >2 absolute fold difference between the
         sum of ground truths and the sum of estimated counts are
         flagged

    2. Transcripts with poor correlation between estimated counts

       - spline fitted to relationship between correlation and kmer fraction.
         cut-off of 0.9 used to define minimum kmer fraction threshold.
         transcripts below threshold are flagged

    2. is not yet implemented. Currently the minimum kmer fraction is
    hardcoded as 0.03. Need to implement automatic threshold
    generation from data
    '''

    job_memory = "2G"

    TranscriptDiffExpression.identifyLowConfidenceTranscripts(
        infile, outfile, submit=True)


@transform(identifyLowConfidenceTranscript,
           suffix(".tsv"),
           ".load")
def loadLowConfidenceTranscripts(infile, outfile):
    ''' load the low confidence transcripts '''
    P.load(infile, outfile)


@mkdir("simulation.dir")
@follows(loadCorrelation,
         loadLowConfidenceTranscripts)
def simulate():
    pass


###############################################################################
# Remove flagged transcripts
###############################################################################

# Add task to optionally remove flagged transcripts

###############################################################################
# Estimate transcript abundance
###############################################################################

@mkdir("quant.dir")
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(buildKallistoIndex),
         SEQUENCEFILES_OUTPUT)
def runKallisto(infiles, outfile):
    ''' quantify trancript abundance with kallisto'''

    # TS more elegant way to parse infiles and index?
    infile = [x[0] for x in infiles]
    index = infiles[0][1]

    # multithreading not supported until > v0.42.1
    # job_threads = PARAMS["kallisto_threads"]
    job_threads = 1
    job_memory = "6G"

    kallisto_options = PARAMS["kallisto_options"]
    bootstrap = PARAMS["kallisto_bootstrap"]

    m = PipelineMapping.Kallisto()
    statement = m.build(infile, outfile)

    P.run()


@follows(runKallisto)
def quantify():
    pass


###############################################################################
# Differential isoform expression analysis
###############################################################################


@follows(runKallisto)
@mkdir("DEresults.dir")
@subdivide(["%s.design.tsv" % x.asFile().lower() for x in DESIGNS],
           regex("(\S+).design.tsv"),
           add_inputs(buildReferenceTranscriptome),
           [r"DEresults.dir/\1_results.tsv",
            r"DEresults.dir/\1_counts.tsv",
            r"DEresults.dir/\1_tpm.tsv"])
def runSleuth(infiles, outfiles):
    ''' run Sleuth to perform differential testing '''

    design, transcripts = infiles
    outfile, counts, tpm = outfiles

    Design = Expression.ExperimentalDesign(design)
    number_samples = sum(Design.table['include'])

    number_transcripts = 0
    with IOTools.openFile(transcripts, "r") as inf:
        for line in inf:
            if line.startswith(">"):
                number_transcripts += 1

    # TS: rough estimate is 24 bytes * bootstraps * samples * transcripts
    # (https://groups.google.com/forum/#!topic/kallisto-sleuth-users/mp064J-DRfI)
    # I've found this to be a serious underestimate so this is a more
    # conservative estimate
    memory_estimate = (48 * PARAMS["kallisto_bootstrap"] * number_samples *
                       number_transcripts)
    job_memory = "%fG" % ((memory_estimate / 1073741824))

    design_id = P.snip(design, ".design.tsv")
    model = PARAMS["sleuth_model_%s" % design_id]

    contrasts = PARAMS["sleuth_contrasts_%s" % design_id].split(",")

    for contrast in contrasts:

        TranscriptDiffExpression.runSleuth(
            design, "quant.dir", model, contrast,
            outfile, counts, tpm, PARAMS["sleuth_fdr"],
            submit=True, job_memory=job_memory)


@transform(runSleuth,
           regex("(\S+)_(counts|tpm).tsv"),
           r"\1_\2.load")
def loadSleuthTables(infile, outfile):
    ''' load tables from Sleuth '''
    P.load(infile, outfile)


@transform(runSleuth,
           suffix("_results.tsv"),
           "_withBiotypes.tsv")
def addTranscriptBiotypes(infile, outfile):
    ''' add the transcript biotypes to the results outfile'''
    # TS: This could be done when the report is built but saves time just
    # to just do it once here

    df = pd.read_table(infile, sep="\t", index_col=0)
    df.set_index('test_id', inplace=True)
    df.index.names = ["transcript_id"]

    dbh = connect()

    table = os.path.basename(
        PARAMS["annotations_interface_table_transcript_info"])

    where_cmd = " OR ".join(["WHERE gene_biotype = '%s'" % x
                             for x in PARAMS["geneset_biotypes"].split(",")])

    select = dbh.execute("""
    SELECT DISTINCT
    transcript_id, transcript_biotype, gene_id, gene_name
    FROM annotations.%(table)s
    %(where_cmd)s""" % locals())

    df_annotations = pd.DataFrame.from_records(
        select, index="transcript_id",
        columns=("transcript_id", "transcript_biotype",
                 "gene_id", "gene_name"))

    df = df.join(df_annotations, sort=False)

    df.to_csv(outfile, index=True, sep="\t")


@transform(addTranscriptBiotypes,
           suffix("_withBiotypes.tsv"),
           "_DEresults.load")
def loadSleuthResults(infile, outfile):
    ''' load Sleuth results '''
    P.load(infile, outfile)


@follows(loadSleuthTables,
         loadSleuthResults)
def differentialExpression():
    pass


###############################################################################
# Expression summary plots
###############################################################################

@mkdir("summary_plots")
@transform(runSleuth,
           regex("DEresults.dir/(\S+)_tpm.tsv"),
           add_inputs(r"\1.design.tsv"),
           r"summary_plots/\1_plots.log")
def expressionSummaryPlots(infiles, logfile):
    ''' make summary plots for expression values for each design file'''

    counts_inf, design_inf = infiles

    job_memory = "4G"

    TranscriptDiffExpression.makeExpressionSummaryPlots(
        counts_inf, design_inf, logfile, submit=True, job_memory=job_memory)


###############################################################################
# Generic pipeline tasks
###############################################################################


@follows(index,
         quantify,
         differentialExpression,
         expressionSummaryPlots,
         simulate)
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
