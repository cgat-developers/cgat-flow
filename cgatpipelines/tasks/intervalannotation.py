'''
intervalannotation.py - Tasks associated with annotation of genomic intervals
=====================================================================================
'''

import shutil
import os
import sqlite3
import CGATCore.IOTools as IOTools
import CGAT.IndexedGenome as IndexedGenome
import CGAT.Bed as Bed
from CGATCore import Pipeline as P

############################################################
############################################################
############################################################
# Pipeline configuration
P.get_parameters(["%s.yml" % __file__[:-len(".py")], "pipeline.yml"])
PARAMS = P.PARAMS

############################################################
############################################################
############################################################


def exportIntervalsAsBed(database, query, outfile):
    '''export intervals from SQlite database as bed files. '''

    dbhandle = sqlite3.connect(database)
    cc = dbhandle.cursor()
    cc.execute(query)

    outs = IOTools.open_file(outfile, "w")
    for result in cc:
        contig, start, end, interval_id, score = result
        outs.write("%s\t%i\t%i\t%s\t%s\n" %
                   (contig, start, end, str(interval_id), str(score)))
    cc.close()
    outs.close()

############################################################
############################################################
############################################################


def BedFileVenn(infiles, outfile):
    '''merge :term:`bed` formatted *infiles* by intersection
    and write to *outfile*.

    Only intervals that overlap in all files are retained.
    Interval coordinates are given by the first file in *infiles*.

    Bed files are normalized (overlapping intervals within 
    a file are merged) before intersection. 

    Intervals are renumbered starting from 1.
    '''
    bed1, bed2 = infiles
    liver_name = P.snip(os.path.basename(liver), ".replicated.bed")
    testes_name = P.snip(os.path.basename(testes), ".replicated.bed")
    to_cluster = True

    statement = '''cat %(liver)s %(testes)s | mergeBed -i stdin | awk 'OFS="\\t" {print $1,$2,$3,"CAPseq"NR}' > replicated_intervals/liver.testes.merge.bed;
                   echo "Total merged intervals" > %(outfile)s; cat replicated_intervals/liver.testes.merge.bed | wc -l >> %(outfile)s; 
                   echo "Liver & testes" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.merge.bed -b %(liver)s -u | intersectBed -a stdin -b %(testes)s -u > replicated_intervals/liver.testes.shared.bed; cat replicated_intervals/liver.testes.shared.bed | wc -l >> %(outfile)s; 
                   echo "Testes only" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.merge.bed -b %(liver)s -v > replicated_intervals/%(testes_name)s.liver.testes.unique.bed; cat replicated_intervals/%(testes_name)s.liver.testes.unique.bed | wc -l >> %(outfile)s; 
                   echo "Liver only" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.merge.bed -b %(testes)s -v > replicated_intervals/%(liver_name)s.liver.testes.unique.bed; cat replicated_intervals/%(liver_name)s.liver.testes.unique.bed | wc -l >> %(outfile)s;                   
                   sed -i '{N;s/\\n/\\t/g}' %(outfile)s; '''

    if len(infiles) == 1:
        shutil.copyfile(infiles[0], outfile)

    elif len(infiles) == 2:

        if IOTools.is_empty(infiles[0]) or IOTools.isEmpty(infiles[1]):
            IOTools.touch_file(outfile)
        else:
            statement = '''
        intersectBed -u -a %s -b %s 
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        > %%(outfile)s 
        ''' % (infiles[0], infiles[1])
            P.run(statement)

    else:

        tmpfile = P.get_temp_filename(".")

        # need to merge incrementally
        fn = infiles[0]
        if IOTools.is_empty(infiles[0]):
            IOTools.touch_file(outfile)
            return

        statement = '''mergeBed -i %(fn)s > %(tmpfile)s'''
        P.run(statement)

        for fn in infiles[1:]:
            if IOTools.is_empty(infiles[0]):
                IOTools.touch_file(outfile)
                os.unlink(tmpfile)
                return

            statement = '''mergeBed -i %(fn)s | intersectBed -u -a %(tmpfile)s -b stdin > %(tmpfile)s.tmp; mv %(tmpfile)s.tmp %(tmpfile)s'''
            P.run(statement)

        statement = '''cat %(tmpfile)s
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        > %(outfile)s '''
        P.run(statement)

        os.unlink(tmpfile)


############################################################
############################################################
############################################################
def makeIntervalCorrelation(infiles, outfile, field, reference):
    '''compute correlation of interval properties between sets
    '''

    dbhandle = sqlite3.connect(PARAMS["database_name"])

    tracks, idx = [], []
    for infile in infiles:
        track = P.snip(infile, ".bed")
        tablename = "%s_intervals" % P.tablequote(track)
        cc = dbhandle.cursor()
        statement = "SELECT contig, start, end, %(field)s FROM %(tablename)s" % locals(
        )
        cc.execute(statement)
        ix = IndexedGenome.IndexedGenome()
        for contig, start, end, peakval in cc:
            ix.add(contig, start, end, peakval)
        idx.append(ix)
        tracks.append(track)
    outs = IOTools.open_file(outfile, "w")
    outs.write("contig\tstart\tend\tid\t" + "\t".join(tracks) + "\n")

    for bed in Bed.iterator(infile=open(reference, "r")):

        row = []
        for ix in idx:
            try:
                intervals = list(ix.get(bed.contig, bed.start, bed.end))
            except KeyError:
                row.append("")
                continue

            if len(intervals) == 0:
                peakval = ""
            else:
                peakval = str((max([x[2] for x in intervals])))
            row.append(peakval)

        outs.write(str(bed) + "\t" + "\t".join(row) + "\n")

    outs.close()
