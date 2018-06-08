"""PipelineReadqc.py - Tasks for QC'ing short read data sets
============================================================

The majority of the functions in this module are for running
and processing the information from the fastqc_ tool.

Reference
---------

"""

import os
import re
import glob
import collections
from io import StringIO
import pandas as pd
from CGATCore import Pipeline as P
import CGATCore.IOTools as IOTools
import CGATCore.CSV2DB as CSV2DB


def FastqcSectionIterator(infile):
    """iterate over FASTQC output file and yield each section.

    Sections in FASTQC output start with `>>`` and end with
    ``>>END_MODULE``.

    Yields
    ------
    name : string
        Section name
    status : string
        Section status
    header : string
        Section header
    data : list
        Lines within section

    Arguments
    ---------
    infile : iterator
        Iterator over contents of Fastqc output.

    """
    data = []
    name, status, header, data = None, None, None, None
    for line in infile:
        if line.startswith(">>END_MODULE"):
            yield name, status, header, data
        elif line.startswith(">>"):
            name, status = line[2:-1].split("\t")
            data = []
        elif line.startswith("#"):
            header = "\t".join([x for x in line[1:-1].split("\t") if x != ""])
        else:
            data.append(
                "\t".join([x for x in line[:-1].split("\t") if x != ""]))


def collectFastQCSections(infiles, section, datadir):
    '''iterate over all fastqc files and extract a particular section.

    Arguments
    ---------
    infiles : list
        List of filenames with fastqc output (logging information). The
        track name is derived from that.
    section : string
        Section name to extract
    datadir : string
        Location of actual Fastqc output to be parsed.

    Returns
    -------
    results : list
        List of tuples, one tuple per input file. Each tuple contains
        track, status, header and data of `section`.

    '''
    results = []
    for infile in infiles:
        track = P.snip(os.path.basename(infile), ".fastqc")
        filename = os.path.join(datadir, track + "*_fastqc", "fastqc_data.txt")
        for fn in glob.glob(filename):
            for name, status, header, data in FastqcSectionIterator(
                    IOTools.open_file(fn)):
                if name == section:
                    results.append((track, status, header, data))
    return results


def loadFastqc(filename,
               database_url):
    '''load FASTQC statistics into database.

    Each section will be uploaded to its own table.

    Arguments
    ----------
    filename : string
        Filename with FASTQC data
    database_url : string
        Database backend.
    '''

    parser = CSV2DB.buildParser()
    (options, args) = parser.parse_args([])

    options.database_url = database_url
    options.database_schema = None
    options.allow_empty = True

    for fn in glob.glob(filename):
        prefix = os.path.basename(os.path.dirname(fn))
        results = []

        for name, status, header, data in FastqcSectionIterator(
                IOTools.open_file(fn)):
            # do not collect basic stats, see loadFastQCSummary
            if name == "Basic Statistics":
                continue
            options.tablename = prefix + "_" + re.sub(" ", "_", name)

            inf = StringIO("\n".join([header] + data) + "\n")
            CSV2DB.run(inf, options)
            results.append((name, status))

        # load status table
        options.tablename = prefix + "_status"

        inf = StringIO(
            "\n".join(["name\tstatus"] +
                      ["\t".join(x) for x in results]) + "\n")
        CSV2DB.run(inf, options)


def buildFastQCSummaryStatus(infiles, outfile, datadir):
    '''collect fastqc status results from multiple runs into a single table.

    Arguments
    ---------
    infiles : list
        List of filenames with fastqc output (logging information). The
        track name is derived from that.
    outfile : list
        Output filename in :term:`tsv` format.
    datadir : string
        Location of actual Fastqc output to be parsed.

    '''

    outf = IOTools.open_file(outfile, "w")
    names = set()
    results = []
    for infile in infiles:
        track = P.snip(os.path.basename(infile), ".fastqc")
        filename = os.path.join(datadir,
                                track + "*_fastqc",
                                "fastqc_data.txt")
        # there can be missing sections
        for fn in glob.glob(filename):
            stats = collections.defaultdict(str)
            for name, status, header, data in FastqcSectionIterator(
                    IOTools.open_file(fn)):
                stats[name] = status

            results.append((track, fn, stats))
            names.update(list(stats.keys()))

    names = sorted(names)
    outf.write("track\tfilename\t%s\n" % "\t".join(names))
    for track, fn, stats in results:
        outf.write("%s\t%s\t%s\n" %
                   (track, os.path.dirname(fn),
                    "\t".join(stats[x] for x in names)))
    outf.close()


def buildFastQCSummaryBasicStatistics(infiles, outfile, datadir):
    '''collect fastqc summary results from multiple runs into a single table.

    Arguments
    ---------
    infiles : list
        List of filenames with fastqc output (logging information). The
        track name is derived from that.
    outfile : list
        Output filename in :term:`tsv` format.
    datadir : string
        Location of actual Fastqc output to be parsed.

    '''

    data = collectFastQCSections(infiles, "Basic Statistics", datadir)

    outf = IOTools.open_file(outfile, "w")
    first = True
    for track, status, header, rows in data:
        rows = [x.split("\t") for x in rows]
        if first:
            headers = [row[0] for row in rows]
            outf.write("track\t%s\n" % "\t".join(headers))
            first = False
        outf.write("%s\t%s\n" % (track, "\t".join([row[1] for row in rows])))
    outf.close()


def buildExperimentReadQuality(infiles, outfile, datadir):
    """build per-experiment read quality summary.

    Arguments
    ---------
    infiles : list
        List of filenames with fastqc output (logging information). The
        track name is derived from that.
    outfile : list
        Output filename in :term:`tsv` format.
    datadir : string
        Location of actual Fastqc output to be parsed.

    """
    data = collectFastQCSections(infiles,
                                 "Per sequence quality scores",
                                 datadir)
    first = True

    if len(data) == 0:
        raise ValueError("received no data")

    for track, status, header, rows in data:
        T = track
        rows = [list(map(float, x.split("\t"))) for x in rows]
        header = header.split("\t")
        if first:
            first = False
            df_out = pd.DataFrame(rows)
            df_out.columns = header
            df_out.rename(columns={"Count": track}, inplace=True)
        else:
            df = pd.DataFrame(rows)
            df.columns = header
            df.rename(columns={"Count": track}, inplace=True)
            df_out = df_out.merge(df, how="outer", on="Quality", sort=True)

    df_out.set_index("Quality", inplace=True)
    df_out = pd.DataFrame(df_out.sum(axis=1))
    df_out.columns = ["_".join(T.split("-")[:-1]), ]

    df_out.to_csv(IOTools.open_file(outfile, "w"), sep="\t")


def read_fastqc(infiles, track_regex, sep="-"):
    """merge multiple fastq output into multiple dataframes.

    Arguments
    ---------
    infiles : string
        Input filename with fastqscreen output.
    regex_track: string
        Regular expression to extract track name from filename.
    sep: char
        Separator for merging multiple capture groups in regex.

    Returns
    -------
    dataframes
    """

    dfs, tracks = collections.defaultdict(list), []
    for infile in infiles:
        try:
            track = sep.join(re.search(track_regex, infile).groups())
        except AttributeError:
            raise ValueError("regex {} did not match file {}".format(
                track_regex, infile))
        tracks.append(track)
        with IOTools.open_file(infile) as inf:
            for name, status, header, data in FastqcSectionIterator(inf):
                records = [x.split("\t") for x in data]
                df = pd.DataFrame.from_records(records, columns=header.split("\t"))
                dfs[name].append(df)

    result = {}
    for key, dd in dfs.items():
        df = pd.concat(dd, keys=tracks, names=["track"])
        df.index = df.index.droplevel(1)
        key = re.sub(" ", "_", key.lower())
        result[key] = df
    return result


def read_fastq_screen(infiles, track_regex, sep="-"):
    """merge fastqscreen output into dataframes.

    Arguments
    ---------
    infiles : string
        Input filename with fastqscreen output.
    regex_track: string
        Regular expression to extract track name from filename.
    sep: char
        Separator for merging multiple capture groups in regex.

    Returns
    -------
    multiple dataframes
    """

    dfs, tracks, summaries = [], [], []
    for infile in infiles:

        try:
            track = sep.join(re.search(track_regex, infile).groups())
        except AttributeError:
            raise ValueError("regex {} did not match file {}".format(
                track_regex, infile))

        with IOTools.open_file(infile) as inf:
            lines = inf.readlines()
        version, aligner, reads = re.search(
            "#Fastq_screen version: (\S+)\t#Aligner: (\S+)\t#Reads in subset: (\d+)\n",
            lines.pop(0)).groups()
        percent_no_hit = re.search(
            "%Hit_no_genomes: (\S+)\n", lines.pop(-1)).groups()[0]

        summaries.append((version, aligner, reads, percent_no_hit))

        records = [x[:-1].split("\t") for x in lines if x.strip()]
        df = pd.DataFrame.from_records(records[1:], columns=records[0])
        df = df.rename(columns={
            'Genome': "genome",
            '#Reads_processed': "reads_processed",
            '#Unmapped': "reads_unmapped",
            '%Unmapped': "reads_unmapped_percent",
            '#One_hit_one_genome': "one_hit_one_genome",
            '%One_hit_one_genome': "one_hit_one_genome_percent",
            '#Multiple_hits_one_genome': "multiple_hits_one_genome",
            '%Multiple_hits_one_genome': "multiple_hits_one_genome_percent",
            '#One_hit_multiple_genomes': "one_hit_multiple_genomes",
            '%One_hit_multiple_genomes': "one_hit_multiple_genomes_percent",
            'Multiple_hits_multiple_genomes': "multiple_hits_multiple_genomes",
            '%Multiple_hits_multiple_genomes': "multiple_hits_multiple_genomes"})
        dfs.append(df)
        tracks.append(track)
    df_details = pd.concat(dfs, keys=tracks, names=["track"])
    df_details.index = df_details.index.droplevel(1)
    df_summary = pd.DataFrame.from_records(
        summaries, columns=["version", "aligner", "nreads", "nohit_percent"],
        index=tracks)
    df_summary.index.name = "track"
    return df_summary, df_details
