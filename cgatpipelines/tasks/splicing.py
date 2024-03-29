'''
splicing.py - wrap various differential expression tools
===========================================================

Purpose
-------

This module provides tools for differential splicing analysis, which
are abstracted modules that can ease the use of the tools and provide
additional functionality. They generate a command line statement.

Methods implemented are:

   rMATS

DEXSeq is implemented elsewhere (cgat code collection: counts2table),
as it is closely related to counts-based differential expression tools.


Usage
-----

The basic usage inside a pipeline task is as such::

    @transform()
    def runMATS(infile, outfile):

        # build statment and run statement
        statement = splicing.runRMATS(
             gtf_string, designfile_string, PARAMS['pvalue_threshold'],
             PARAMS['strandedness'], outdir_string, permute=1)


When implementing a tool, avoid specifying algorithmic options as
class variables. Instead use an option string that can be set in
:file:`pipeline.yml`. The only arguments to a tool constructor should
pertain to pipeline integration, such as filenames, index locations,
threading and in general processing options that change the tools
input/output, as these need to be tracked by the pipeline.

The benefit of this approach is to provide compleeat control to the
user and is likely to work with different versions of a tool, i.e., if
a command line option to a tool changes, just the configuration file
needs to be changed, but the code remains the same.

Requirements
------------

* rMATS-turbo



'''

import os
import random
import shutil
import pandas as pd
import cgat.BamTools.bamtools as BamTools
import cgatcore.experiment as E
import cgatpipelines.tasks.expression as Expression
from cgatcore import pipeline as P
from cgatpipelines.tasks.mapping import SequenceCollectionProcessor


def runRMATS(gtffile, designfile, pvalue, strand, outdir, permute=0):
    '''Module to generate rMATS statment

    Module offers the option to permute group name labels and
    calculates readlength, which must be identical in all reads.

    Arguments
    ---------
    gtffile: string
        path to :term:`gtf` file
    designfile: string
        path to design file
    pvalue: string
        threshold for FDR testing
    strand: string
        strandedness option: can be 'fr-unstranded', 'fr-firststrand',
        or 'fr-secondstrand'
    outdir: string
        directory path for rMATS results
    permute : 1 or 0
        option to activate random shuffling of sample groups
    '''

    design = Expression.ExperimentalDesign(designfile)
    if permute == 1:
        permutelist = design.table.group.tolist()
        random.shuffle(permutelist)
        design.table.group = permutelist
    group1 = ",".join(
        ["%s.bam" % x for x in design.getSamplesInGroup(design.groups[0])])
    with open(outdir + "/b1.txt", "w") as f:
        f.write(group1)
    group2 = ",".join(
        ["%s.bam" % x for x in design.getSamplesInGroup(design.groups[1])])
    with open(outdir + "/b2.txt", "w") as f:
        f.write(group2)
    readlength = BamTools.estimateTagSize(design.samples[0]+".bam")
    tmpdir = P.get_temp_dir()

    statement = '''rmats.py
    --b1 %(outdir)s/b1.txt
    --b2 %(outdir)s/b2.txt
    --gtf <(gunzip -c %(gtffile)s)
    --od %(outdir)s
    --readLength %(readlength)s
    --cstat %(pvalue)s
    --libType %(strand)s
    --tmp %(tmpdir)s
    ''' % locals()

    # if Paired End Reads
    if BamTools.is_paired(design.samples[0]+".bam"):
        statement += '''-t paired''' % locals()

    statement += '''
    > %(outdir)s/%(designfile)s.log
    '''

    P.run(statement)


def rmats2sashimi(infile, designfile, FDR, outfile, plotmax=20):
    '''Module to generate sashimi plots from rMATS output

    Module generates a statement to call rmats2sashimiplot and provides
    it with correct arguments. Only results containing no NA in results
    and below FDR threshold are drawn to prevent unneccassary compute and
    memory use.

    Arguments
    ---------
    infile: string
        path to rMATS results file (can be one of five types)
    designfile: string
        path to design file
    FDR: string
        FDR threshold for drawing plots'
    outfile: string
        directory path for sashimiplot output
    '''

    Design = Expression.ExperimentalDesign(designfile)
    if len(Design.groups) != 2:
        raise ValueError("Please specify exactly 2 groups per experiment.")

    g1 = Design.getSamplesInGroup(Design.groups[0])
    g2 = Design.getSamplesInGroup(Design.groups[1])

    if len(g1) != len(g2):
        g1 = g1[:min(len(g1), len(g2))]
        g2 = g2[:min(len(g1), len(g2))]
        E.info("The two groups compared were of unequal size. For  " +
               "visual display using sashimi they have been truncated " +
               "to the same length")

    group1 = ",".join(["%s.bam" % x for x in g1])
    group2 = ",".join(["%s.bam" % x for x in g2])
    group1name = Design.groups[0]
    group2name = Design.groups[1]
    event = os.path.basename(os.path.normpath(outfile))
    if "MXE" in infile:
        column = "22"
        sortby = "25"
    else:
        column = "20"
        sortby = "23"

    temp = pd.read_csv(infile, sep='\t')
    temp = temp.dropna()
    temp.sort_values(by=['FDR'], inplace=True)
    temp = temp[abs(temp['IncLevelDifference']) > 0.2]
    plotnum = min(int(len(temp[temp['FDR'] < float(FDR)])), plotmax)
    infile = P.snip(infile)
    temp.iloc[1:plotnum,:].to_csv("%s_plot.txt" % infile, sep='\t', index=False)

    statement = '''
    rmats2sashimiplot
    --b1 %(group1)s
    --b2 %(group2)s
    -t %(event)s
    -e %(infile)s_plot.txt
    --l1 %(group1name)s
    --l2 %(group2name)s
    -o %(outfile)s
    > %(outfile)s/%(event)s.log
    ''' % locals()

    P.run(statement)


def diffIR(singularity, designfile, outdir, irratio="0.05", permute=0):
    '''Module to generate IRFinder Diff statment

    Module offers the option to permute group name labels.

    Arguments
    ---------
    singularity: string
        path to IRFinder :term:`singularity` file
    designfile: string
        path to design file
    irratio: string
        minimum IR ratio met in at least one sample
    outdir: string
        directory path for rMATS results
    permute : 1 or 0
        option to activate random shuffling of sample groups
    '''

    design = Expression.ExperimentalDesign(designfile)
    if permute == 1:
        permutelist = design.table.group.tolist()
        random.shuffle(permutelist)
        design.table.group = permutelist

    group1_list = ["IRFinder.dir/%s/IRFinder-IR-nondir.txt" % x for x in design.getSamplesInGroup(design.groups[0])]
    group1_list_dir = [filename.replace("nondir","dir") for filename in group1_list]
    if os.path.exists(group1_list_dir[1]):
        group1_list = group1_list_dir
    groupname1 = design.groups[0]
    group1 = " ".join(["".join(["-g:",groupname1," ",item]) for item in group1_list])  

    group2_list = ["IRFinder.dir/%s/IRFinder-IR-nondir.txt" % x for x in design.getSamplesInGroup(design.groups[1])]
    group2_list_dir = [filename.replace("nondir","dir") for filename in group2_list]
    if os.path.exists(group2_list_dir[1]):
        group2_list = group2_list_dir
    groupname2 = design.groups[1]  
    group2 = " ".join(["".join(["-g:",groupname2," ",item]) for item in group2_list])

    statement = '''
    singularity run -H $PWD:/home
    %(singularity)s Diff
    %(group1)s
    %(group2)s
    -ir %(irratio)s
    -o %(outdir)s
    ''' % locals()

    statement += '''
    > %(outdir)s/%(designfile)s.log
    '''

    P.run(statement)
