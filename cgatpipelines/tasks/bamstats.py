"""
================================================
bamstats.py - Tasks for QC'ing Bam Files
================================================

:Author: Adam Cribbs
:Release: $Id$
:Date: |today|
:Tags: Python

Perfroming stats on a bam file is an important task in quality control
of `.bam` files. The majority of tools focus on quality checking fastq files,
however mapping QC should be performed to assess the quality of the alignment
accross the whole genome. This module provides utility functions to abstract
some of these variations.



Requirements:
* picardtools >= 1.106
* cgat tools
* bamstats >=1.22
* pandas >= 0.18.0

Reference
---------
"""

# load modules
import os
import re
import gzip
import pysam
import pandas as pd

# load cgat specific modules
import cgatcore.Experiment as E
import cgatcore.IOTools as IOTools
import cgat.BamTools.bamtools as BamTools
import cgatcore.Pipeline as P

PICARD_MEMORY = "12G"


def getNumReadsFromReadsFile(infile):
    '''get number of reads from a .nreads file.'''
    with IOTools.open_file(infile) as inf:
        line = inf.readline()
        if not line.startswith("nreads"):
            raise ValueError(
                "parsing error in file '%s': "
                "expected first line to start with 'nreads'")
        nreads = line.split("\t")[1]
        nreads = int(nreads)
    return nreads


def buildPicardInsertSizeStats(infile, outfile, genome_file):
    '''run Picard:CollectInsertSizeMetrics
    Collect insert size statistics.
    Arguments
    ---------
    infile : string
        Input filename in :term:`BAM` format.
    outfile : string
        Output filename with picard output.
    genome_file : string
        Filename with genomic sequence.
    '''
    job_memory = PICARD_MEMORY
    picard_opts = '-Xmx%(job_memory)s -XX:+UseParNewGC -XX:+UseConcMarkSweepGC' % locals()
    job_threads = 3

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        IOTools.touch_file(outfile)
        return

    statement = '''picard %(picard_opts)s CollectInsertSizeMetrics
    INPUT=%(infile)s
    REFERENCE_SEQUENCE=%(genome_file)s
    ASSUME_SORTED=true
    OUTPUT=%(outfile)s
    VALIDATION_STRINGENCY=SILENT
    >& %(outfile)s'''

    P.run(statement, job_memory=PICARD_MEMORY)


def addPseudoSequenceQuality(infile, outfile):
    '''to allow multiQC to pick up the picard metric files
    sequence quality needs to be added if it is stripped.
    Therefore an intermediate sam file needs to be generated
    Arguments
    ---------
    infile : string
        Input file in :term:`BAM` format.
    outfile : string
        Output file in :term: `BAM` format.
    '''

    statement = '''cat %(infile)s
    | cgat bam2bam -v 0
    --method=set-sequence > %(outfile)s'''

    P.run(statement)

    statement = '''samtools index %(outfile)s
    '''

    P.run(statement)


def copyBamFile(infile, outfile):
    '''Make softlinks of the bam files

    Arguments
    ---------
    infile : string
        Input file in :term:`BAM` format.
    outfile : string
        Output file in :term: `BAM` format.
    '''

    statement = '''ln -s ../%(infile)s
    %(outfile)s'''

    P.run(statement)

    statement = '''samtools index %(outfile)s
    '''

    P.run(statement)


def buildPicardAlignmentStats(infile, outfile, genome_file):
    '''run picard:CollectMultipleMetrics
    Arguments
    ---------
    infile : string
        Input filename in :term:`BAM` format.
    outfile : string
        Output filename with picard output.
    genome_file : string
        Filename with genomic sequence.
    '''

    job_memory = PICARD_MEMORY
    picard_opts = '-Xmx%(job_memory)s -XX:+UseParNewGC -XX:+UseConcMarkSweepGC' % locals()
    job_threads = 3

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        IOTools.touch_file(outfile)
        return

    statement = '''picard %(picard_opts)s CollectMultipleMetrics
    INPUT=%(infile)s
    REFERENCE_SEQUENCE=%(genome_file)s
    ASSUME_SORTED=true
    OUTPUT=%(outfile)s
    VALIDATION_STRINGENCY=SILENT
    >& %(outfile)s'''

    P.run(statement)


def buildPicardDuplicationStats(infile, outfile):
    '''run picard:MarkDuplicates
    Record duplicate metrics using Picard, the marked records
    are discarded.
    Arguments
    ---------
    infile : string
        Input filename in :term:`BAM` format.
    outfile : string
        Output filename with picard output.
    '''

    job_memory = PICARD_MEMORY
    picard_opts = '-Xmx%(job_memory)s -XX:+UseParNewGC -XX:+UseConcMarkSweepGC' % locals()
    job_threads = 3

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        IOTools.touch_file(outfile)
        return

    # currently, MarkDuplicates cannot handle split alignments from gsnap
    # these can be identified by the custom XT tag.
    if ".gsnap.bam" in infile:
        tmpf = P.get_temp_file(".")
        tmpfile_name = tmpf.name
        statement = '''samtools view -h %(infile)s
        | awk "!/\\tXT:/"
        | samtools view /dev/stdin -S -b > %(tmpfile_name)s;
        ''' % locals()
        data_source = tmpfile_name
    else:
        statement = ""
        data_source = infile

    statement += ''' picard %(picard_opts)s MarkDuplicates
    INPUT=%(data_source)s
    ASSUME_SORTED=true
    METRICS_FILE=%(outfile)s
    OUTPUT=/dev/null
    VALIDATION_STRINGENCY=SILENT
    >& %(outfile)s.log
    '''
    P.run(statement)

    if ".gsnap.bam" in infile:
        os.unlink(tmpfile_name)


def buildPicardDuplicateStats(infile, outfile):
    '''run picard:MarkDuplicates
    Record duplicate metrics using Picard and keep the dedupped .bam
    file.
    Pair duplication is properly handled, including inter-chromosomal
    cases. SE data is also handled.  These stats also contain a
    histogram that estimates the return from additional sequecing.  No
    marked bam files are retained (/dev/null...)  Note that picards
    counts reads but they are in fact alignments.
    Arguments
    ---------
    infile : string
        Input filename in :term:`BAM` format.
    outfile : string
        Output filename with picard output.
    '''
    job_memory = PICARD_MEMORY
    picard_opts = '-Xmx%(job_memory)s -XX:+UseParNewGC -XX:+UseConcMarkSweepGC' % locals()
    job_threads = 3

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        IOTools.touch_file(outfile)
        return

    statement = '''picard %(picard_opts)s MarkDuplicates
    INPUT=%(infile)s
    ASSUME_SORTED=true
    METRICS_FILE=%(outfile)s.duplicate_metrics
    OUTPUT=%(outfile)s
    VALIDATION_STRINGENCY=SILENT
    >& %(outfile)s.log &&
    samtools index %(outfile)s'''
    P.run(statement)


def buildPicardCoverageStats(infile, outfile, baits, regions):
    '''run picard:CollectHSMetrics
    Generate coverage statistics for regions of interest from a bed
    file using Picard.
    Arguments
    ---------
    infile : string
        Input filename in :term:`BAM` format.
    outfile : string
        Output filename with picard output.
    baits : :term:`bed` formatted file of bait regions
    regions : :term:`bed` formatted file of target regions
    '''

    job_memory = PICARD_MEMORY
    picard_opts = '-Xmx%(job_memory)s -XX:+UseParNewGC -XX:+UseConcMarkSweepGC' % locals()
    job_threads = 3

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        IOTools.touch_file(outfile)
        return

    statement = '''picard %(picard_opts)s CollectHsMetrics
    BAIT_INTERVALS=%(baits)s
    TARGET_INTERVALS=%(regions)s
    INPUT=%(infile)s
    OUTPUT=%(outfile)s
    VALIDATION_STRINGENCY=LENIENT''' % locals()
    P.run(statement)


def buildPicardGCStats(infile, outfile, genome_file):
    """picard:CollectGCBiasMetrics
    Collect GC bias metrics.
    Arguments
    ---------
    infile : string
        Input filename in :term:`BAM` format.
    outfile : string
        Output filename with picard output.
    genome_file : string
        Filename with genomic sequence.
    """

    job_memory = PICARD_MEMORY
    picard_opts = '-Xmx%(job_memory)s -XX:+UseParNewGC -XX:+UseConcMarkSweepGC' % locals()
    job_threads = 3

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        IOTools.touch_file(outfile)
        return

    statement = '''picard %(picard_opts)s CollectGcBiasMetrics
    INPUT=%(infile)s
    REFERENCE_SEQUENCE=%(genome_file)s
    OUTPUT=%(outfile)s
    VALIDATION_STRINGENCY=SILENT
    CHART_OUTPUT=%(outfile)s.pdf
    SUMMARY_OUTPUT=%(outfile)s.summary
    >& %(outfile)s'''

    P.run(statement)


def defineBedFeatures(infile, outfile):
    '''
    This module process genomic context file.
    It assigns each and every features of context
    file to a specific catagory. It helps us to
    understand heiarchical classification
    of features.
    '''
    category_feature = {'IG_': 'IG_gene', 'TR_': 'TR_gene',
                        'antisense': 'long_noncoding_RNA',
                        'overlapping': 'long_noncoding_RNA',
                        'processed_transcript': 'unclassified_noncoding_RNA',
                        'non_coding': 'long_noncoding_RNA',
                        'TEC': 'to_be_experimentally_confirmed',
                        'protein': 'protein_coding_gene',
                        'intergenic': 'INTERGENIC',
                        'prime_utr': 'UTR'
                        }
    category_RNA = {'ncRNA': 'long_noncoding_RNA',
                    '-RTE': 'repeats',
                    '-Deu': 'repeats'}

    '''Read genomic context file'''
    infile = IOTools.open_file(infile)
    '''Create processed genomic context file'''
    outfile = IOTools.open_file(outfile, "w")
    '''Start processing of genomic context file'''
    for line in infile:
        line_copy = line.split("\t")
        name = ""
        if 'RNA' in line and 'pseudogene' not in line:
            flag = 0
            for i in category_RNA.keys():
                if i in line:
                    name = category_RNA[i]
                    flag = 1
                    break
            if(flag == 0):
                if(line_copy[3] == "RNA\n"):
                    name = "unclassified_RNA"
                else:
                    name = "short_noncoding_RNA"
        else:
            flag = 0
            for i in category_feature.keys():
                if i in line and 'pseudogene' not in line:
                    name = category_feature[i]
                    flag = 1
                    break
            if(flag == 0):
                if('pseudogene' in line):
                    name = "pseudogene"
                elif(line_copy[3] == "intron\n"):
                    name = "introns"
                elif('intron' in line):
                    name = "long_noncoding_RNA"
                else:
                    name = "repeats"
        outfile.write(line)
        outfile.write(("%s\t%s\t%s\t%s\n") % (line_copy[0], str(line_copy[1]), str(line_copy[2]), name))
    outfile.close()
    infile.close()


def summarizeTagsWithinContext(tagfile,
                               contextfile,
                               outfile,
                               min_overlap=0.5,
                               job_memory="15G"):
    '''count occurances of tags in genomic context.

    Examines the genomic context to where tags align.

    A tag is assigned to the genomic context that it
    overlaps by at least 50%. Thus some reads mapping
    several contexts might be dropped.

    Arguments
    ---------
    tagfile : string
        Filename with tags. The file can be :term:`bam` or :term:`bed` format.
    contextfile : string
        Filename of :term:`bed` formatted files with named intervals (BED4).
    outfile : string
        Output in :term:`tsv` format.
    min_overlap : float
        Minimum overlap (fraction) to count features as overlapping.
    job_memory : string
        Memory to reserve.
    '''

    tmpfile = P.get_temp_filename(shared=True)
    tmpfiles = ["%s_%i" % (tmpfile, x) for x in range(2)]
    statement = '''
    cgat bam_vs_bed
    --min-overlap=%(min_overlap)f
    --log=%(outfile)s.log
    %(tagfile)s %(contextfile)s
    > %(tmpfile)s_0
    '''

    P.run(statement)

    statement = '''
    printf "intergenic\\t" >> %(tmpfile)s_1'''

    P.run(statement)

    statement = '''
    bedtools intersect -a %(tagfile)s
    -b %(contextfile)s
    -bed -v | wc -l
    | xargs printf
    >> %(tmpfile)s_1
    '''
    P.run(statement)

    files = " ".join(tmpfiles)
    statement = '''
    sort --merge  %(files)s
    | gzip > %(outfile)s
    '''
    P.run(statement)

    for x in tmpfiles:
        os.unlink(x)


def loadTranscriptProfile(infiles, outfile,
                          suffix="transcript_profile",
                          tablename=None):
    '''load transcript profiles into one table.
    Arguments
    ---------
    infiles : string
        Filenames of files with matrix from bam2geneprofile. Each file
        corresponds to a different track.
    outfile : string
        Logfile.
    suffix : string
        Suffix to append to table name.
    pipeline_suffix : string
        Suffix to remove from track name.
    tablename : string
        Tablename to use. If unset, the table name will be derived
        from `outfile` and suffix as ``to_table(outfile) + "_" +
        suffix``.
    '''

    if not tablename:
        tablename = "%s" % (suffix)

    outf = P.get_temp_file(".")

    table_count = 0
    table_join = None

    for infile in infiles:

        matrix_file = str(infile) + ".geneprofileabsolutedistancefromthreeprimeend.matrix.tsv.gz"
        name = P.snip(os.path.basename(infile), ".transcriptprofile.gz")

        table = pd.read_csv(matrix_file, sep="\t")
        table.rename(columns={'none': name}, inplace=True)
        table.drop(["area", "counts", "background"], axis=1, inplace=True)

        if table_count == 0:
            table_join = table
            table_count += 1
        else:
            table_join = table.merge(table_join,
                                     on=["bin", "region", "region_bin"],
                                     how="left")
    table_join.to_csv(outf, sep="\t", index=False)

    outf.close()

    P.load(infile=outf.name,
           outfile=outfile,
           tablename=tablename,
           options="--add-index=bin")

    os.unlink(outf.name)


def loadStrandSpecificity(infiles, outfile,
                          suffix="strand",
                          tablename=None):
    '''
    '''

    if not tablename:
        tablename = "%s_%s" % (P.to_table(outfile), suffix)

    outf = P.get_temp_file(".")

    table_count = 0
    table_join = None

    for infile in infiles:
        name = P.snip(os.path.basename(infile), ".strand")

        table = pd.read_csv(infile, sep="\t", comment="#")
        table["track"] = name

        if table_count == 0:
            table_join = table
            table_count += 1
        else:
            table_join = table.merge(table_join,
                                     on=["MSR", "ISR", "OSR", "ISF", "MSF", "OSF", "SF", "SR", "track"],
                                     how="outer")

    table_join.to_csv(outf, sep="\t", index=False)

    outf.close()

    P.load(infile=outf.name,
           outfile=outfile,
           tablename=tablename,
           options="--add-index=track")

    os.unlink(outf.name)


def loadCountReads(infiles, outfile,
                   suffix="nreads",
                   pipeline_suffix=".nreads",
                   tablename=None):
    '''load read counts.
    Arguments
    ---------
    infiles : string
        Filenames of files with number of reads per sample. Each file
        corresponds to a different track.
    outfile : string
        Logfile.
    suffix : string
        Suffix to append to table name.
    pipeline_suffix : string
        Suffix to remove from track name.
    tablename : string
        Tablename to use. If unset, the table name will be derived
        from `outfile` and suffix as ``to_table(outfile) + "_" +
        suffix``.
    '''

    if not tablename:
        tablename = "%s_%s" % (P.to_table(outfile), suffix)

    outf = P.get_temp_file(".")

    outf.write("%s\t%s\n" % ("track", "nreads"))

    for filename in infiles:
        track = P.snip(os.path.basename(filename), pipeline_suffix)

        if not os.path.exists(filename):
            E.warn("File %s missing" % filename)
            continue

        lines = IOTools.open_file(filename, "r").readlines()

        for line in lines:
            count = line.split("\t")[1]
            outf.write("%s\t%s\n" % (track, count))

    outf.close()

    P.load(infile=outf.name,
           outfile=outfile,
           tablename=tablename,
           options="--add-index=track")

    os.unlink(outf.name)


def loadPicardMetrics(infiles, outfile, suffix,
                      pipeline_suffix=".picard_stats",
                      tablename=None):
    '''load picard metrics.
    Arguments
    ---------
    infiles : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfile : string
        Logfile.
    suffix : string
        Suffix to append to table name.
    pipeline_suffix : string
        Suffix to remove from track name.
    tablename : string
        Tablename to use. If unset, the table name will be derived
        from `outfile` and suffix as ``to_table(outfile) + "_" +
        suffix``.
    '''

    if not tablename:
        tablename = "%s_%s" % (P.to_table(outfile), suffix)

    outf = P.get_temp_file(".")

    filenames = ["%s.%s" % (x, suffix) for x in infiles]

    first = True
    for filename in filenames:
        track = P.snip(os.path.basename(filename), "%s.%s" %
                       (pipeline_suffix, suffix))

        if not os.path.exists(filename):
            E.warn("File %s missing" % filename)
            continue

        lines = IOTools.open_file(filename, "r").readlines()

        # extract metrics part
        rx_start = re.compile("## METRICS CLASS")
        for n, line in enumerate(lines):
            if rx_start.search(line):
                lines = lines[n + 1:]
                break

        for n, line in enumerate(lines):
            if not line.strip():
                lines = lines[:n]
                break

        if len(lines) == 0:
            E.warn("no lines in %s: %s" % (track, filename))
            continue

        if first:
            outf.write("%s\t%s" % ("track", lines[0]))
            fields = lines[0][:-1].split("\t")
        else:
            f = lines[0][:-1].split("\t")
            if f != fields:
                raise ValueError(
                    "file %s has different fields: expected %s, got %s" %
                    (filename, fields, f))

        first = False
        for i in range(1, len(lines)):
            outf.write("%s\t%s" % (track, lines[i]))

    outf.close()

    P.load(outf.name,
           outfile,
           tablename=tablename,
           options="--add-index=track --allow-empty-file")

    os.unlink(outf.name)


def loadPicardHistogram(infiles, outfile, suffix, column,
                        pipeline_suffix=".picard_stats", tablename=False):
    '''extract a histogram from a picard output file and load
    it into database.
    Arguments
    ---------
    infiles : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfile : string
        Logfile.
    suffix : string
        Suffix to append to table name.
    column : string
        Column name to take from the histogram.
    pipeline_suffix : string
        Suffix to remove from track name.
    tablename : string
        Tablename to use. If unset, the table name will be derived
        from `outfile` and suffix as ``to_table(outfile) + "_" +
        suffix``.
    '''

    if not tablename:
        tablename = "%s_%s" % (P.to_table(outfile), suffix)
        tablename = tablename.replace("_metrics", "_histogram")

    # some files might be missing
    xfiles = [x for x in infiles if os.path.exists("%s.%s" % (x, suffix))]

    if len(xfiles) == 0:
        E.warn("no files for %s" % tablename)
        return

    header = ",".join([P.snip(os.path.basename(x), pipeline_suffix)
                       for x in xfiles])
    filenames = " ".join(["%s.%s" % (x, suffix) for x in xfiles])

    # there might be a variable number of columns in the tables
    # only take the first ignoring the rest

    load_statement = P.build_load_statement(
        tablename,
        options="--add-index=track "
        " --header-names=%s,%s"
        " --allow-empty-file"
        " --replace-header" % (column, header))

    statement = """cgat combine_tables
    --regex-start="## HISTOGRAM"
    --missing-value=0
    --take=2
    %(filenames)s
    | %(load_statement)s
    >> %(outfile)s
    """

    P.run(statement)


def loadPicardAlignmentStats(infiles, outfile):
    '''load all output from Picard's CollectMultipleMetrics into database.
    Loads tables into database with prefix derived from outfile:
       * [outfile]_alignment_summary_metric
       * [outfile]_insert_size_metrics
       * [outfile]_quality_by_cycle_metrics
       * [outfile]_quality_distribution_metrics
       * [outfile]_insert_size_metrics
    Arguments
    ---------
    infiles : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfile : string
        Logfile. The table name will be derived from `outfile`.
    '''

    loadPicardMetrics(infiles, outfile, "alignment_summary_metrics")

    # insert size metrics only available for paired-ended data
    loadPicardMetrics(infiles, outfile, "insert_size_metrics")

    histograms = (("quality_by_cycle_metrics", "cycle"),
                  ("quality_distribution_metrics", "quality"),
                  ("insert_size_metrics", "insert_size"))

    for suffix, column in histograms:
        loadPicardHistogram(infiles, outfile, suffix, column)


def loadPicardDuplicationStats(infiles, outfiles):
    '''load picard duplicate filtering stats into database.
    Loads two tables into the database
       * picard_duplication_metrics
       * picard_complexity_histogram
    Arguments
    ---------
    infiles : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfile : string
        Logfile. The table name will be derived from `outfile`.
    '''
    # SNS: added to enable naming consistency

    outfile_metrics, outfile_histogram = outfiles

    suffix = "picard_duplication_metrics"

    # the loading functions expect "infile_name.pipeline_suffix" as the infile
    # names.
    infile_names = [x[:-len("." + suffix)] for x in infiles]

    loadPicardMetrics(infile_names, outfile_metrics, suffix, "",
                      tablename="picard_duplication_metrics")

    infiles_with_histograms = []

    # The complexity histogram is only present for PE data, so we must check
    # because by design the pipeline does not track endedness
    for infile in infile_names:
        with_hist = False
        with open(".".join([infile, suffix]), "r") as open_infile:
            for line in open_infile:
                if line.startswith("## HISTOGRAM"):
                    infiles_with_histograms.append(infile)
                    break

    if len(infiles_with_histograms) > 0:
        loadPicardHistogram(infiles_with_histograms,
                            outfile_histogram,
                            suffix,
                            "coverage_multiple",
                            "",
                            tablename="picard_complexity_histogram")
    else:
        with open(outfile_histogram, "w") as ofh:
            ofh.write("No histograms detected, no data loaded.")


def loadPicardDuplicateStats(infiles, outfile, pipeline_suffix=".bam"):
    '''load picard duplicate filtering stats.
    Arguments
    ---------
    infiles : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfile : string
        Logfile. The table name will be derived from `outfile`.
    pipeline_suffix : string
        Suffix appended to pipeline output file, will be removed to
        define track.
    '''

    loadPicardMetrics(
        infiles, outfile, "duplicate_metrics",
        pipeline_suffix=pipeline_suffix)
    loadPicardHistogram(infiles,
                        outfile,
                        "duplicate_metrics",
                        "duplicates",
                        pipeline_suffix=pipeline_suffix)


def loadPicardCoverageStats(infiles, outfile):
    '''import coverage statistics into database.
    Arguments
    ---------
    infiles : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfile : string
        Logfile. The table name will be derived from `outfile`.
    '''

    outf = P.get_temp_file(".")
    first = True
    for f in infiles:
        track = P.snip(os.path.basename(f), ".cov")
        lines = [x for x in open(f, "r").readlines()
                 if not x.startswith("#") and x.strip()]
        if first:
            outf.write("%s\t%s" % ("track", lines[0]))
        first = False
        outf.write("%s\t%s" % (track, lines[1]))
    outf.close()
    P.load(outf.name,
           outfile,
           options="--ignore-empty --add-index=track")
    os.unlink(outf.name)


def buildBAMStats(infile, outfile):
    '''Count number of reads mapped, duplicates, etc. using
    bam2stats.py
    Arguments
    ---------
    infile : string
        Input file in :term:`BAM` format
    outfile : string
        Output file in :term:`tsv` format.
    '''

    statement = '''cgat bam2stats
    --force-output
    --output-filename-pattern=%(outfile)s.%%s
    < %(infile)s
    > %(outfile)s'''
    P.run(statement)


def loadBAMStats(infiles, outfile):
    '''load output of :func:`buildBAMStats` into database.
    Arguments
    ---------
    infiles : string
        Input files, output from :func:`buildBAMStats`.
    outfile : string
        Logfile. The table name will be derived from `outfile`.
    '''

    header = ",".join([P.snip(os.path.basename(x), ".readstats")
                       for x in infiles])
    filenames = " ".join(["<( cut -f 1,2 < %s)" % x for x in infiles])
    tablename = P.to_table(outfile)

    load_statement = P.build_load_statement(
        tablename,
        options="--add-index=track "
        " --allow-empty-file")

    E.info("loading bam stats - summary")
    statement = """cgat combine_tables
    --header-names=%(header)s
    --missing-value=0
    --ignore-empty
    %(filenames)s
    | perl -p -e "s/bin/track/"
    | cgat table2table --transpose
    | %(load_statement)s
    > %(outfile)s"""
    P.run(statement)

    for suffix in ("nm", "nh"):
        E.info("loading bam stats - %s" % suffix)
        filenames = " ".join(["%s.%s" % (x, suffix) for x in infiles])

        load_statement = P.build_load_statement(
            "%s_%s" % (tablename, suffix),
            options="--allow-empty-file")

        statement = """cgat combine_tables
        --header-names=%(header)s
        --skip-titles
        --missing-value=0
        --ignore-empty
        %(filenames)s
        | perl -p -e "s/bin/%(suffix)s/"
        | %(load_statement)s
        >> %(outfile)s """
        P.run(statement)

    # load mapping qualities, there are two columns per row
    # 'all_reads' and 'filtered_reads'
    # Here, only filtered_reads are used (--take=3)
    for suffix in ("mapq",):
        E.info("loading bam stats - %s" % suffix)
        filenames = " ".join(["%s.%s" % (x, suffix) for x in infiles])

        load_statement = P.build_load_statement(
            "%s_%s" % (tablename, suffix),
            options=" --allow-empty-file")

        statement = """cgat combine_tables
        --header-names=%(header)s
        --skip-titles
        --missing-value=0
        --ignore-empty
        --take=3
        %(filenames)s
        | perl -p -e "s/bin/%(suffix)s/"
        | %(load_statement)s
        >> %(outfile)s """
        P.run(statement)


def loadPicardRnaSeqMetrics(infiles, outfiles):
    '''load picard rna stats into database.
    Loads tables into the database
       * picard_rna_metrics
       * picard_rna_histogram
    Arguments
    ---------
    infile : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfiles : string
        Logfile. The table names will be derived from `outfile`.
    '''

    outfile_metrics, outfile_histogram = outfiles

    suffix = "picard_rna_metrics"

    # the loading functions expect "infile_name.pipeline_suffix" as the infile
    # names.
    infile_names = [x[:-len("." + suffix)] for x in infiles]

    loadPicardMetrics(infile_names, outfile_metrics, suffix, "",
                      tablename="picard_rna_metrics")

    infiles_with_histograms = []

    # Checking if histogram is present (?is this necessary)
    for infile in infile_names:
        with_hist = False
        with open(".".join([infile, suffix]), "r") as open_infile:
            for line in open_infile:
                if line.startswith("## HISTOGRAM"):
                    infiles_with_histograms.append(infile)
                    break

    if len(infiles_with_histograms) > 0:
        loadPicardHistogram(infiles_with_histograms,
                            outfile_histogram,
                            suffix,
                            "coverage_multiple",
                            "",
                            tablename="picard_rna_histogram")
    else:
        with open(outfile_histogram, "w") as ofh:
            ofh.write("No histograms detected, no data loaded.")


def loadIdxstats(infiles, outfile):
    '''take list of file paths to samtools idxstats output files
    and merge to create single dataframe containing mapped reads per
    contig for each track. This dataframe is then loaded into
    database.
    Loads tables into the database
        * idxstats_reads_per_chromosome
    Arguments
    ---------
    infiles : list
        list where each element is a string of the filename containing samtools
        idxstats output. Filename format is expected to be 'sample.idxstats'
    outfile : string
        Logfile. The table name will be derived from `outfile`.
    '''

    outf = P.get_temp_file(".")
    dfs = []
    for f in infiles:
        track = P.snip(f, ".idxstats").split('/')[-1]

        if not os.path.exists(f):
            E.warn("File %s missing" % f)
            continue

        # reformat idx stats
        df = pd.read_csv(f, sep='\t', header=None)
        df.columns = ['region', 'length', 'mapped', 'unmapped']

        # calc total reads mapped & unmappedpep
        total_reads = df.unmapped.sum() + df.mapped.sum()
        total_mapped_reads = df.mapped.sum()

        reformatted_df = pd.DataFrame([['total_mapped_reads',
                                        total_mapped_reads],
                                       ['total_reads', total_reads],
                                       ['track', track]],
                                      columns=(['region', 'mapped']))

        # reformat the df
        df = df.append(reformatted_df, ignore_index=True)
        df.set_index('region', inplace=True)
        df1 = df[['mapped']].T
        dfs.append(df1)

    # merge dataframes into single table
    master_df = pd.concat(dfs)
    master_df.drop('*', axis=1, inplace=True)
    master_df.to_csv(outf, sep='\t', index=False)
    outf.close()

    P.load(outf.name,
           outfile,
           options="--ignore-empty --add-index=track")
    os.unlink(outf.name)


def loadSummarizedContextStats(infiles,
                               outfile,
                               suffix=".contextstats.tsv.gz"):
    """merge output from :func:`summarizeTagsWithinContex` and load into database.

    Arguments
    ---------
    infiles : list
        List of filenames in :term:`tsv` format. The files should end
        in suffix.
    outfile : string
        Output filename, the table name is derived from `outfile`.
    suffix : string
        Suffix to remove from filename for track name.

    """

    header = ",".join([P.snip(os.path.basename(x), suffix)
                       for x in infiles])
    filenames = " ".join(infiles)

    load_statement = P.build_load_statement(
        P.to_table(outfile),
        options="--add-index=track")

    statement = """cgat combine_tables
    --header-names=%(header)s
    --missing-value=0
    --skip-titles
    %(filenames)s
    | perl -p -e "s/bin/track/; s/\?/Q/g"
    | cgat table2table --transpose
    | %(load_statement)s
    > %(outfile)s
    """
    P.run(statement)

# dont know if this is used anymore


def mergeAndFilterGTF(infile, outfile, logfile,
                      genome,
                      max_intron_size=None,
                      remove_contigs=None,
                      rna_file=None):

    '''sanitize transcripts file for cufflinks analysis.
    Merge exons separated by small introns (< 5bp).
    Transcripts will be ignored that
       * have very long introns (`max_intron_size`) (otherwise,
         cufflinks complains)
       * are located on contigs to be ignored (usually: chrM, _random, ...)
    Optionally remove transcripts based on repetitive sequences by
       supplying a repetitve `rna_file`.
    This method preserves all features in a gtf file (exon, CDS, ...).
    Arguments
    ---------
    infile : string
       Input filename in :term:`gtf` format
    outfile : string
       Output filename in :term:`gtf` format
    logfile : string
       Output filename for logging information.
    genome : string
       Filename (without extension) of indexed genome file
       in :term:`fasta` format.
    max_intron_size : int
       Remove transripts with introns larger than this value.
    remove_contigs : string
       Remove transcripts on contigs matching this regular
       expression string.
    rna_file : string
       Filename of :term:`gff` formatted file with repetetive
       sequences. If given, all transcripts overlapping any regions
       in this file will be removed.
    Returns
    -------
    kept_genes : dict
        a dictionary of all gene_ids that have been kept.
        '''

    c = E.Counter()

    outf = gzip.open(outfile, "w")

    E.info("filtering by contig and removing long introns")
    contigs = set(IndexedFasta.IndexedFasta(genome).getContigs())

    rx_contigs = None
    #
    if remove_contigs is not None:
        rx_contigs = re.compile(remove_contigs)
        E.info("removing contigs %s" % remove_contigs)

    rna_index = None
    if rna_file is not None:
        if not os.path.exists(rna_file):
            E.warn("file '%s' to remove repetetive rna does not exist" %
                   rna_file)
        else:
            rna_index = GTF.readAndIndex(
                GTF.iterator(IOTools.open_file(rna_file, "r")))
            E.info("removing ribosomal RNA in %s" % rna_file)

    gene_ids = {}

    logf = IOTools.open_file(logfile, "w")
    logf.write("gene_id\ttranscript_id\treason\n")

    for all_exons in GTF.transcript_iterator(
            GTF.iterator(IOTools.open_file(infile))):

        c.input += 1

        e = all_exons[0]
        # filtering
        if e.contig not in contigs:
            c.missing_contig += 1
            logf.write(
                "\t".join((e.gene_id, e.transcript_id,
                           "missing_contig")) + "\n")
            continue

        if rx_contigs and rx_contigs.search(e.contig):
            c.remove_contig += 1
            logf.write(
                "\t".join((e.gene_id, e.transcript_id,
                           "remove_contig")) + "\n")
            continue

        if rna_index and all_exons[0].source != 'protein_coding':
            found = False
            for exon in all_exons:
                if rna_index.contains(e.contig, e.start, e.end):
                    found = True
                    break
            if found:
                logf.write(
                    "\t".join((e.gene_id, e.transcript_id,
                               "overlap_rna")) + "\n")
                c.overlap_rna += 1
                continue

        is_ok = True

        # keep exons and cds separate by grouping by feature
        all_exons.sort(key=lambda x: x.feature)
        new_exons = []

        for feature, exons in itertools.groupby(
                all_exons, lambda x: x.feature):

            tmp = sorted(list(exons), key=lambda x: x.start)

            gene_ids[tmp[0].transcript_id] = tmp[0].gene_id

            l, n = tmp[0], []

            for e in tmp[1:]:
                d = e.start - l.end
                if max_intron_size and d > max_intron_size:
                    is_ok = False
                    break
                elif d < 5:
                    l.end = max(e.end, l.end)
                    c.merged += 1
                    continue

                n.append(l)
                l = e

            n.append(l)
            new_exons.extend(n)

            if not is_ok:
                break

        if not is_ok:
            logf.write(
                "\t".join((e.gene_id, e.transcript_id,
                           "bad_transcript")) + "\n")
            c.skipped += 1
            continue

        new_exons.sort(key=lambda x: (x.start, x.gene_id, x.transcript_id))

        for e in new_exons:
            outf.write("%s\n" % str(e))
            c.exons += 1

        c.output += 1

    outf.close()
    L.info("%s" % str(c))

    return gene_ids


def resetGTFAttributes(infile, genome, gene_ids, outfile):
    """set GTF attributes in :term:`gtf` formatted file so that they are
    compatible with cufflinks.
    This method runs cuffcompare with `infile` against itself to add
    attributes such as p_id and tss_id.
    Arguments
    ---------
    infile : string
        Filename of :term:`gtf`-formatted input file
    genome : string
       Filename (without extension) of indexed genome file
       in :term:`fasta` format.
    gene_ids : dict
       Dictionary mapping transcript ids to gene ids.
    outfile : string
       Output filename in :term:`gtf` format
    """
    tmpfile1 = P.get_temp_filename(".")
    tmpfile2 = P.get_temp_filename(".")

    #################################################
    E.info("adding tss_id and p_id")

    # The p_id attribute is set if the fasta sequence is given.
    # However, there might be some errors in cuffdiff downstream:
    #
    # cuffdiff: bundles.cpp:479: static void HitBundle::combine(const std::
    # vector<HitBundle*, std::allocator<HitBundle*> >&, HitBundle&): Assertion
    # `in_bundles[i]->ref_id() == in_bundles[i-1]->ref_id()' failed.
    #
    # I was not able to resolve this, it was a complex
    # bug dependent on both the read libraries and the input reference gtf
    # files
    job_memory = "5G"

    statement = '''
    cuffcompare -r <( gunzip < %(infile)s )
         -T
         -s %(genome)s.fa
         -o %(tmpfile1)s
         <( gunzip < %(infile)s )
         <( gunzip < %(infile)s )
    > %(outfile)s.log
    '''
    P.run(statement)

    #################################################
    E.info("resetting gene_id and transcript_id")

    # reset gene_id and transcript_id to ENSEMBL ids
    # cufflinks patch:
    # make tss_id and p_id unique for each gene id
    outf = IOTools.open_file(tmpfile2, "w")
    map_tss2gene, map_pid2gene = {}, {}
    inf = IOTools.open_file(tmpfile1 + ".combined.gtf")

    def _map(gtf, key, val, m):
        if val in m:
            while gene_id != m[val]:
                val += "a"
                if val not in m:
                    break
        m[val] = gene_id

        gtf.setAttribute(key, val)

    for gtf in GTF.iterator(inf):
        transcript_id = gtf.oId
        gene_id = gene_ids[transcript_id]
        gtf.setAttribute("transcript_id", transcript_id)
        gtf.setAttribute("gene_id", gene_id)

        # set tss_id
        try:
            tss_id = gtf.tss_id
        except AttributeError:
            tss_id = None
        try:
            p_id = gtf.p_id
        except AttributeError:
            p_id = None

        if tss_id:
            _map(gtf, "tss_id", tss_id, map_tss2gene)
        if p_id:
            _map(gtf, "p_id", p_id, map_pid2gene)

        outf.write(str(gtf) + "\n")

    outf.close()

    # sort gtf file
    geneset.sortGTF(tmpfile2, outfile)

    # make sure tmpfile1 is NEVER empty
    assert tmpfile1
    for x in glob.glob(tmpfile1 + "*"):
        os.unlink(x)
    os.unlink(tmpfile2)


def buildPicardRnaSeqMetrics(infiles, strand, outfile):
    '''run picard:RNASeqMetrics



    Arguments
    ---------
    infiles : string
        Input filename in :term:`BAM` format.
        Genome file in refflat format
            (http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat)
    outfile : string
        Output filename with picard output.

    '''
    job_memory = PICARD_MEMORY
    picard_opts = '-Xmx%(job_memory)s -XX:+UseParNewGC -XX:+UseConcMarkSweepGC' % locals()
    job_threads = 3
    infile, genome = infiles

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        IOTools.touch_file(outfile)
        return

    statement = '''picard %(picard_opts)s CollectRnaSeqMetrics
    REF_FLAT=%(genome)s
    INPUT=%(infile)s
    ASSUME_SORTED=true
    OUTPUT=%(outfile)s
    STRAND=%(strand)s
    VALIDATION_STRINGENCY=SILENT
    '''
    P.run(statement)
