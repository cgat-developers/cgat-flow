"""rnaseq.py - Tasks for RNAseq analysis
==============================================

This module provides tasks related to RNAseq analysis.


Requirements:

* HiddenMarkov >= 1.8.0
* MASS >= 7.3.34
* RColorBrewer >= 1.0.5
* featureCounts >= 1.4.3
* samtools >= 1.1

Reference
----------

"""

import cgatcore.experiment as E
import cgatcore.csvutils as CSV
import cgat.Sra as Sra

import collections
import glob
import itertools
import math
import numpy as np
import os
import pandas as pd
import re

import cgat.BamTools.bamtools as BamTools
import cgatpipelines.tasks.counts as Counts
import cgatcore.database as Database
import cgatpipelines.tasks.expression as Expression
import cgat.GTF as GTF
import cgatcore.iotools as iotools
from cgatcore import pipeline as P
import cgatpipelines.tasks.mapping as mapping

from cgatcore.pipeline import cluster_runnable




def normaliseCounts(counts_inf, outfile, method):
    ''' normalise a counts table'''

    counts = Counts.Counts(counts_inf)
    counts.normalise(method=method)
    counts.table.index.name = "id"

    if outfile.endswith(".gz"):
        counts.table.to_csv(outfile, sep="\t", compression="gzip")
    else:
        counts.table.to_csv(outfile, sep="\t")


def parse_table(sample, outfile_raw, outfile, columnname):
    '''
    parse the output of featurecounts or alignment free qauntifiers
    and extract number of reads for downstream quantification
    '''

    column_ix = findColumnPosition(outfile_raw, columnname)
    sample = sample

    if outfile_raw.endswith("gz"):
        grep = "zgrep"
    else:
        grep = "grep"

    statement = '''
    echo -e "id\\t%(sample)s" | gzip > %(outfile)s;
    %(grep)s -v '^#' %(outfile_raw)s |
    cut -f1,%(column_ix)s | awk 'NR>1' | gzip >> %(outfile)s;
    ''' % locals()
    P.run(statement)


def estimateSleuthMemory(bootstraps, samples, transcripts):
    ''' The memory usage of Sleuth is dependent upon the number of
    samples, transcripts and bootsraps.

    A rough estimate is:
    24 bytes * bootstraps * samples * transcripts
    (https://groups.google.com/forum/#!topic/kallisto-sleuth-users/mp064J-DRfI)

    TS: I've found this to be a serious underestimate so we use a
    more conservative estimate here (48 bytes * ... ) with a default of 2G
    '''

    estimate = (48 * bootstraps * samples * transcripts)

    job_memory = "%fG" % (max(2.0, (estimate / 1073741824.0)))

    return job_memory


def findColumnPosition(infile, column):
    ''' find the position in the header of the specified column
    The returned value is one-based (e.g for bash cut)'''
    with iotools.open_file(infile, "r") as inf:
        while True:
            header = inf.readline()
            if not header.startswith("#"):
                head = header.strip().split("\t")
                j = 0
                for h in head:
                    if column in h:
                        column_ix = j
                    j += 1
                # column_ix = header.strip().split("\t").index(column)
                break
        if column_ix:
            return column_ix + 1
        else:
            raise ValueError("could not find %s in file header" % column)


class Quantifier(object):
    ''' base class for transcript and gene-level quantification from a
    BAM or fastq

    Note: All quantifier classes are designed to perform both
    transcript-level and gene-level analyses in a single run. The
    runGene method
    '''

    def __init__(self, infile, transcript_outfile, gene_outfile, annotations,
                 job_threads=None, job_memory=None, strand=None,
                 options=None, bootstrap=None,
                 fragment_length=None, fragment_sd=None,
                 transcript2geneMap=None, libtype=None, kmer=None,
                 biascorrect=None):
        '''
        Attributes
        ----------
        infile: string
           Input  filename
        transcript_outfile: string
           Outfile of transcript quantifications in :term: `gz.raw` format
        gene_outfile: string
           Outfile of gene quantifications in :term: `gz.raw` format
        job_threads: string
           Number of threads per job
        strand: int
           For FeatureCounts the strand is specified as either 0, 1, 2
        options: string
           Options specified as a string
        annotations: string
           Filename with gene set in :term:`gtf` format.
        bootstrap: int
           Number of boostrap values for alignment free quantifiers
        job_memory: str
           Amount of memory available for job
        frangment_length: int
           Must specify the expected fragment length for single-end reads
           This is specified in pipeline_ini.
           :term:`PARAMS` - fragment_length option.
        frangment_sd: int
           Must specify the expected fragment length sd for single-end reads
           This is specified in pipeline_ini.
           :term:`PARAMS` - fragment_sd option.
        libtype: string
           This is specified in pipeline_ini
           :term:`PARAMS` - library type option.
        kmer: int
           This is specified in the pipeline.yml
           :term:`PARAMS` - kmer size for aligment free quant.
        '''

        self.infile = infile
        self.transcript_outfile = transcript_outfile
        self.gene_outfile = gene_outfile
        self.job_threads = job_threads
        self.strand = strand
        self.options = options
        self.annotations = annotations
        self.bootstrap = bootstrap
        self.job_memory = job_memory
        self.fragment_length = fragment_length
        self.fragment_sd = fragment_sd
        self.t2gMap = transcript2geneMap
        self.libtype = libtype
        self.kmer = kmer
        self.biascorrect = biascorrect

        # TS: assume sample name is directory for outfile which is
        # should be for pipeline_rnaseqdiffexpression. This would be
        # better handled in pipeline though
        self.sample = os.path.basename(os.path.dirname(self.gene_outfile))

    def run_transcript(self):
        ''' generate transcript-level quantification estimates'''

    def run_gene(self):
        ''' generate gene-level quantification estimates'''

    def run_all(self):
        ''' '''
        self.run_transcript()
        self.run_gene()


class FeatureCountsQuantifier(Quantifier):
    ''' quantifier class to run featureCounts '''

    def run_featurecounts(self, level="gene_id"):
        ''' function to run featureCounts at the transcript-level or gene-level

        If `bamfile` is paired, paired-end counting is enabled and the bam
        file automatically sorted.

        '''

        bamfile = self.infile
        job_threads = self.job_threads
        strand = self.strand
        options = self.options
        annotations = self.annotations
        sample = self.sample

        if level == "gene_id":
            outfile = self.gene_outfile
        elif level == "transcript_id":
            outfile = self.transcript_outfile
        else:
            raise ValueError("level must be gene_id or transcript_id!")

        # -p -B specifies count fragments rather than reads, and both
        # reads must map to the feature
        # for legacy reasons look at feature_counts_paired
        if BamTools.is_paired(bamfile):
            # sort bamfile
            bam_tmp = '${TMPDIR}/' + os.path.basename(bamfile)
            # select paired end mode, additional options
            paired_options = "-p -B"
            # sort by read name
            paired_processing = \
                """samtools
                sort -@ %(job_threads)i -n -o %(bam_tmp)s %(bamfile)s;
                """ % locals()
            bamfile = bam_tmp
        else:
            paired_options = ""
            paired_processing = ""

        # raw featureCounts output saved to ".raw" file
        outfile_raw = P.snip(outfile, ".gz") + ".raw"
        outfile_dir = os.path.dirname(outfile)
        if not os.path.exists(outfile_dir):
            os.makedirs(outfile_dir)

        statement = '''zcat %(annotations)s > ${TMPDIR}/geneset.gtf;
                       %(paired_processing)s
                       featureCounts %(options)s
                                     -T %(job_threads)i
                                     -s %(strand)s
                                     -a ${TMPDIR}/geneset.gtf
                                     %(paired_options)s
                                     -o %(outfile_raw)s -g %(level)s
                                     %(bamfile)s
                        >& %(outfile)s.log;
                        gzip -f %(outfile_raw)s
        '''
        P.run(statement)

        # parse output to extract counts
        parse_table(self.sample, outfile_raw + ".gz",
                    outfile, "%s.bam" % self.sample)

    def run_transcript(self):
        ''' generate transcript-level quantification estimates'''
        self.run_featurecounts(level="transcript_id")

    def run_gene(self):
        ''' generate gene-level quantification estimates'''
        self.run_featurecounts(level="gene_id")


class Gtf2tableQuantifier(Quantifier):
    ''' quantifier class to run gtf2table'''

    def run_gtf2table(self, level="gene_id"):
        ''' function to run gtf2table script at the transcript-level
        or gene-level'''

        bamfile = self.infile
        annotations = self.annotations
        sample = self.sample

        # define the quantification level
        if level == "gene_id":
            outfile = self.gene_outfile
            reporter = "genes"
        elif level == "transcript_id":
            outfile = self.transcript_outfile
            reporter = "transcripts"
        else:
            raise ValueError("level must be gene_id or transcript_id!")

        if BamTools.is_paired(bamfile):
            counter = 'readpair-counts'
        else:
            counter = 'read-counts'

        outfile_raw = P.snip(outfile, ".gz") + ".raw"

        outfile_dir = os.path.dirname(outfile)
        if not os.path.exists(outfile_dir):
            os.makedirs(outfile_dir)

        # ignore multi-mapping reads ("--multi-mapping-method=ignore")
        statement = '''
        zcat %(annotations)s
        | cgat gtf2table
              --reporter=%(reporter)s
              --bam-file=%(bamfile)s
              --counter=length
              --column-prefix="exons_"
              --counter=%(counter)s
              --column-prefix=""
              --counter=read-coverage
              --column-prefix=coverage_
              --min-mapping-quality=%(counting_min_mapping_quality)i
              --multi-mapping-method=ignore
              --log=%(outfile_raw)s.log
        > %(outfile_raw)s;
        gzip -f %(outfile_raw)s
        '''

        P.run(statement)
        # parse output to extract counts
        parse_table(self.sample, outfile_raw + ".gz", outfile, 'counted_all')

    def run_transcript(self):
        ''' generate transcript-level quantification estimates'''
        self.run_gtf2table(level="transcript_id")

    def run_gene(self):
        ''' generate gene-level quantification estimates'''
        self.run_gtf2table(level="gene_id")


class AF_Quantifier(Quantifier):
    ''' Parent class for all alignment-free quantification methods'''

    def run_gene(self):
        ''' Aggregate transcript counts to generate gene-level counts
        using a map of transript_id to gene_id '''

        transcript_df = pd.read_table(self.transcript_outfile,
                                      sep="\t", index_col=0)
        transcript2gene_df = pd.read_table(self.t2gMap, sep="\t", index_col=0)
        transcript_df = pd.merge(transcript_df, transcript2gene_df,
                                 left_index=True, right_index=True,
                                 how="inner")
        gene_df = pd.DataFrame(transcript_df.groupby(
             'gene_id')[self.sample].sum())
        gene_df.index.name = 'id'

        gene_df.to_csv(
            self.gene_outfile, compression="gzip", sep="\t")


class KallistoQuantifier(AF_Quantifier):
    ''' quantifier class to run kallisto'''

    def run_transcript(self):
        ''' '''
        fastqfile = self.infile
        index = self.annotations
        job_threads = self.job_threads
        job_memory = self.job_memory
        kallisto_options = self.options
        kallisto_bootstrap = self.bootstrap
        kallisto_fragment_length = self.fragment_length
        kallisto_fragment_sd = self.fragment_sd
        outfile = os.path.join(
            os.path.dirname(self.transcript_outfile), "abundance.h5")
        sample = self.sample

        # kallisto output is in binary (".h5") format
        # Supplying a "readable_suffix" to the mapping.Kallisto
        # ensures an additional human readable file is also generated
        readable_suffix = ".tsv"
        m = mapping.Kallisto(readable_suffix=readable_suffix)

        statement = m.build((fastqfile), outfile)

        P.run(statement)

        outfile_readable = outfile + readable_suffix

        # parse the output to extract the counts
        parse_table(self.sample, outfile_readable,
                    self.transcript_outfile, 'est_counts')


class SalmonQuantifier(AF_Quantifier):
    '''quantifier class to run salmon'''
    def run_transcript(self):
        fastqfile = self.infile
        index = self.annotations
        job_threads = self.job_threads
        job_memory = self.job_memory
        biascorrect = self.biascorrect

        salmon_options = self.options
        salmon_bootstrap = self.bootstrap
        salmon_libtype = self.libtype
        outfile = os.path.join(
            os.path.dirname(self.transcript_outfile), "quant.sf")
        sample = self.sample

        m = mapping.Salmon(bias_correct=biascorrect)

        statement = m.build((fastqfile), outfile)

        P.run(statement)

        # parse the output to extract the counts
        parse_table(self.sample, outfile,
                    self.transcript_outfile, 'NumReads')

