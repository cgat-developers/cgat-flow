"""
====================
Exome pipeline
====================


The exome pipeline imports unmapped reads from one or more fastq or
sra files and aligns them to the genome using BWA.  Post alignment
quality control is performed using Picard.  The pipeline then performs
local realignment around indels and base quality score recalibration
using GATK.  Next variants (SNVs and indels) are called, annotated,
and filtered according to various inheritance models (de novo,
dominant, and recessive).


   1. Align to genome using gapped alignment (BWA-MEM)
   2. Check alignment quality and target region coverage (Picard)
   3. Local realignment and BQSR in GATK and deduplication in Picard
   4. Calculate the ratio of reads on the X and Y chromosomes to assert sex
   5. Variant calling in families (SNVs & indels) using GATK HaplotypeCaller
   6. Comparison of join VCF genotypes to assess relatedness (somalier)
   7. Variant annotation using SNPeff, GATK VariantAnnotator, and SnpSift
   8. Variant quality score recalibration (GATK)
   9. Flags variants within genes of interest (such as known disease genes)
      (GATK) (optional)
   10. Filters potential de novo variants
   11. Filters potential de novo variants using lower stringency criteria
   12. Filters potential dominant mutations
   13. Filters potential homozygous recessive mutations
   14. Filters potential compound heterozygous mutations
   15. Generates summary statistics for unfiltered vcf file
   16. Generates report

.. note::

   1. Great care should be taken in interpreting lower stringency de
      novo variants as it is expected almost all will be false
      positives. Users should examine them manually in the csvdb
      database and html report.

   2. Great care should be taken when interpreting compound
      heterozygous changes.  Gemini is very permissive and users
      should examine the genotypes in the csvdb table to make sure
      they are consistent with a recessive inheritance pattern.

To do:
   1. Allow users to add other training sets for variant quality score
      recalibration
   2. Allow users to add annotations using SnpSift

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use cgat pipelines.

Configuration
-------------

Input
-----

Reads are imported by placing files or linking to files in the
:term:`working directory`.

The default file format assumes the following convention:

   <family>_<sample>-<condition>-<replicate>.<suffix>

``family`` = "Single", "Trio", or "Multiplex" followed by numerical
identifier.  ``sample`` and ``condition`` make up an
:term:`experiment`, while ``replicate`` denotes the :term:`replicate`
within an :term:`experiment`.  The ``suffix`` determines the file
type. The following suffixes/file types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the
   :file:`fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq.2.gz
   Paired-end reads in fastq format. The two fastq files must be sorted
   by read-pair.

.. note::

   Quality scores need to be of the same scale for all input
   files. Thus it might be difficult to mix different formats.

If you are submitting families then a .ped file for each family and an
all_samples ped file (cat *ped > all_samples.ped) must be supplied
within your working directory.  This is a tab-delimited file named
<family>.ped (where <family> is the family ID in the title of the
corresponding fastq files) with no header and one individual per line
according to the following pattern:

family_id sample_id father_id mother_id sex phenotype

family_id and sample_id should correspond to <family> and
<family>-<sample> in the sra/fastq filenames, father_id and mother_id
should be '0' if unknown, sex should be '1' if male, '2' if female and
'0' if unknown, and phenotype should be '1' if unaffected, '2' if
affected and '0' if unknown.

If you are running the functions to look for compound heterozygotes in
Multiplex families then there is a further requirement for the .ped
files.  The phasing tools expect a trio and therefore any other family
members (other than parents and one child) must be labelled as
unrelated.  That is, the first additional family member could be
labelled "family0" in the family_id column, and subsequent additional
family members could be "family1", "family2" and so on.  For example,
a multiplex family called Multiplex1 may have two parents and two
affected children.  The .ped file would look like this:

Multiplex1 ID1 0 0 1 1
Multiplex1 ID2 0 0 2 1
Multiplex1 ID3 ID1 ID2 1 2
Family0 ID4 ID1 ID2 2 2

Documentation
-------------

If you would like the genes of interest to be flagged in your vcf,
make add_genes_of_interest=1 (default=0) and provide a list of comma
separated genes (without spaces) in the ini file.

Pipeline output
===============

The major output is a csvdb containing quality control information by
sample and variant information by family and an html report with
similar information.

Example
=======

ToDo: make exome sequencing example

Requirements
------------

On top of the default cgat setup, the pipeline requires the following
software to be in the path:

Requirements:

* BWA >= 0.7.8
* picardtools >= 1.106
* samtools >= 1.1
* GATK >= 2.7
* snpEff >= 4.0
* somalier
* Gemini >= ?
* VCFtools >= 0.1.8a

Code
====

"""

# load modules
from ruffus import transform, mkdir, follows, merge, regex, suffix, \
    jobs_limit, files, collate, add_inputs, formatter, \
    active_if, originate
from ruffus.combinatorics import permutations
import sys
import os
import csv
import glob
import re
import shutil
import decimal
import pandas as pd
import cgatcore.experiment as E
import cgatcore.iotools as iotools
from cgatcore import pipeline as P
import cgatpipelines.tasks.mapping as mapping
import cgatpipelines.tasks.bamstats as bamstats
import cgatpipelines.tasks.exome as exome

###############################################################################
###############################################################################
###############################################################################
# load options from the config file


PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0], "pipeline.yml"])

INPUT_FORMATS = ("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz")
REGEX_FORMATS = regex(r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz)")

matches = glob.glob("*.fastq.1.gz") + glob.glob(
    "*.fastq.gz") + glob.glob("*.sra") + glob.glob("*.csfasta.gz")

PICARD_MEMORY = PARAMS["picard_memory"]

GATK_MEMORY = PARAMS["gatk_memory"]


###############################################################################
###############################################################################
###############################################################################
# Load target and sample data into database
# The following functions are designed to upload meta-data to the csvdb
# These haven't been fully implemented yet

@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@files(PARAMS["roi_bed"], "roi.load")
def loadROI(infile, outfile):
    '''Import regions of interest bed file into SQLite.'''
    header = "chr,start,stop,feature"
    tablename = P.to_table(outfile)
    statement = '''cat %(infile)s
            | cgat csv2db %(csv2db_options)s
              --ignore-empty
              --retry
              --header-names=%(header)s
              --table=%(tablename)s
            > %(outfile)s  '''
    P.run(statement)

@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@files(PARAMS["roi_to_gene"], "roi2gene.load")
def loadROI2Gene(infile, outfile):
    '''Import genes mapping to regions of interest bed file into SQLite.'''
    tablename = P.to_table(outfile)
    statement = '''cat %(infile)s
            | cgat csv2db %(csv2db_options)s
              --ignore-empty
              --retry
              --table=%(tablename)s
            > %(outfile)s  '''
    P.run(statement)

@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@files(PARAMS["samples"], "samples.load")
def loadSamples(infile, outfile):
    '''Import sample information into SQLite.'''
    tablename = P.to_table(outfile)
    statement = '''cat %(infile)s
            | cgat csv2db %(csv2db_options)s
              --ignore-empty
              --retry
              --table=%(tablename)s
            > %(outfile)s  '''
    P.run(statement)

@originate("genome.dict")
def createSequenceDictionary(outfile):
    '''creates a sequence dictionary required for running
    gatk/picard tools'''
    genome_file = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))

    statement = '''picard CreateSequenceDictionary 
                        -R %(genome_file)s -O %(outfile)s'''
    P.run(statement)



###############################################################################
###############################################################################
###############################################################################
# Pre-Processing for Variant Discovery

###############################################################################
# Alignment & post-alignment QC & GATK preparation

@follows(mkdir("bam"))
@transform(INPUT_FORMATS, REGEX_FORMATS, r"bam/\1.bam")
def mapReads(infiles, outfile):
    '''Map reads to the genome using BWA-MEM (output=SAM), convert to BAM,
    sort and index BAM file'''
    track = P.snip(os.path.basename(outfile), ".bam")
    m = mapping.BWAMEM(remove_unique=PARAMS["bwa_remove_non_unique"])
    statement = m.build((infiles,), outfile)
    P.run(statement, job_memory="8G", job_threads=PARAMS["bwa_threads"])


@transform(mapReads, regex(r"bam/(\S+).bam"), r"bam/\1.picard_stats")
def PicardAlignStats(infile, outfile):
    '''Run Picard CollectMultipleMetrics on each BAM file'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    bamstats.buildPicardAlignmentStats(infile, outfile, genome, PICARD_MEMORY)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(PicardAlignStats, "picard_stats.load")
def loadPicardAlignStats(infiles, outfile):
    '''Merge Picard alignment stats into single table and load into SQLite.'''
    bamstats.loadPicardAlignmentStats(infiles, outfile)


@follows(loadPicardAlignStats, mkdir("gatk"))
@transform(mapReads, regex(r"bam/(\S+).bam"),
           add_inputs(createSequenceDictionary),
           r"gatk/\1.readgroups.bam")
def GATKReadGroups(infiles, outfile):
    '''Reorders BAM according to reference fasta and adds read groups using
    GATK'''

    infile, dictionary = infiles
    track = re.sub(r'-\w+-\w+\.bam', '', os.path.basename(infile))
    job_threads = PARAMS["gatk_threads"]
    library = PARAMS["readgroup_library"]
    platform = PARAMS["readgroup_platform"]
    platform_unit = PARAMS["readgroup_platform_unit"]
    exome.GATKReadGroups(infile, outfile, dictionary,
                                 library, platform,
                                 platform_unit, track,
                                 GATK_MEMORY)
    iotools.zap_file(infile)


###############################################################################
# Remove duplicates, realign and recalibrate lane-by-lane

@transform(GATKReadGroups,
           regex(r"gatk/(\S+).readgroups.bam"),
           r"gatk/\1.dedup.bam")
def RemoveDuplicatesLane(infile, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    bamstats.buildPicardDuplicateStats(infile, outfile, PICARD_MEMORY)
    iotools.zap_file(infile)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(RemoveDuplicatesLane, "picard_duplicate_stats_lane.load")
def loadPicardDuplicateStatsLane(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    bamstats.loadPicardDuplicateStats(infiles, outfile)


@transform(RemoveDuplicatesLane,
           regex(r"gatk/(\S+).dedup.bam"),
           r"gatk/\1.bqsr.bam")
def GATKBaseRecal(infile, outfile):
    '''recalibrates base quality scores using GATK'''
    intrack = P.snip(os.path.basename(infile), ".bam")
    outtrack = P.snip(os.path.basename(outfile), ".bam")
    dbsnp = PARAMS["gatk_dbsnp"]
    solid_options = PARAMS["gatk_solid_options"]
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    intervals = PARAMS["roi_intervals"]
    padding = PARAMS["roi_padding"]
    if PARAMS["targetted"]:
        shutil.copyfile(infile, outfile)
        shutil.copyfile(intrack + ".bai", outtrack + ".bai")
    else:
        exome.GATKBaseRecal(infile, outfile, genome, intervals,
                                    padding, dbsnp, solid_options,
                                    GATK_MEMORY)
    iotools.zap_file(infile)


###############################################################################
# Merge BAMs across different lanes for the same sample

@collate(GATKBaseRecal,
         regex(r"gatk/(.*).bqsr.bam"),
         r"gatk/\1.merged.bam")
def mergeBAMs(infiles, outfile):
    '''merges BAMs for a single sample over multiple lanes'''
    inputfiles = " INPUT=".join(infiles)
    outf = iotools.open_file(outfile + ".count", "w")
    outf.write(str(len(infiles)))
    outf.close()
    statement = ("picard MergeSamFiles "
                 "INPUT=%(inputfiles)s "
                 "OUTPUT=%(outfile)s "
                 "ASSUME_SORTED=true "
                 ">& %(outfile)s.log && "
                 "samtools index %(outfile)s " % locals())
    P.run(statement, job_memory="8G")

    for inputfile in infiles:
        iotools.zap_file(inputfile)


###############################################################################
# Remove duplicates sample-by-sample

@transform(mergeBAMs,
           regex(r"gatk/(\S+).merged.bam"),
           add_inputs(r"gatk/\1.merged.bam.count"),
           r"gatk/\1.dedup2.bam")
def RemoveDuplicatesSample(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    infile, countfile = infiles
    countf = open(countfile, "r")
    if countf.read() != '1':
        bamstats.buildPicardDuplicateStats(infile, outfile, PICARD_MEMORY)
    else:
        shutil.copyfile(infile, outfile)
        shutil.copyfile(infile + ".bai", outfile + ".bai")
        iotools.zap_file(infile)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(RemoveDuplicatesSample, "picard_duplicate_stats_sample.load")
def loadPicardDuplicateStatsSample(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    bamstats.loadPicardDuplicateStats(infiles, outfile)


###############################################################################
# Coverage of targetted area


@transform(RemoveDuplicatesSample,
           regex(r"gatk/(\S+).dedup2.bam"),
           r"gatk/\1.cov")
def buildCoverageStats(infile, outfile):
    '''Generate coverage statistics for regions of interest from a bed
    file using Picard'''
    baits = PARAMS["roi_baits"]
    regions = PARAMS["roi_regions"]
    bamstats.buildPicardCoverageStats(infile, outfile,
                                               baits, regions, PICARD_MEMORY)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildCoverageStats, "coverage_stats.load")
def loadCoverageStats(infiles, outfile):
    '''Import coverage statistics into SQLite'''
    bamstats.loadPicardCoverageStats(infiles, outfile)



###############################################################################
###############################################################################
###############################################################################
# GATK germline cohort calling


@follows(mkdir("variants"))
@transform(RemoveDuplicatesSample, regex(r"gatk/(\S+).dedup2.bam"),
           r"variants/\1.haplotypeCaller.g.vcf")
def haplotypeCaller(infile, outfile):
    '''Call SNVs and indels using GATK HaplotypeCaller in individuals'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    job_threads = PARAMS["gatk_threads"]
    dbsnp = PARAMS["gatk_dbsnp"]
    intervals = PARAMS["roi_intervals"]
    padding = PARAMS["roi_padding"]
    options = PARAMS["gatk_hc_options"]
    exome.haplotypeCaller(infile, outfile, genome, dbsnp,
                                  intervals, padding, options, GATK_MEMORY)


@merge(haplotypeCaller, "variants/genomicsdb.log")
def consolidateGVCFs(infiles, outfile):
    '''generates a GenomicsDB workspace from all GVCF files, this is an
       easy-access database of all samples developed by the Intel-Broad team.'''
    inputlen = len(infiles)
    inputfiles = " -V ".join(infiles)
    
    
    DB_MEMORY  = PARAMS["gatk_dbmem"]

    exome.consolidateGVCFs(inputfiles, outfile, inputlen, DB_MEMORY)


@transform(consolidateGVCFs, regex(r"variants/genomicsdb.log"),
          "variants/all_samples.vcf")
def genotypeGVCFs(infile, outfile):
    '''Joint genotyping of all samples together'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    options = PARAMS["gatk_genotype_options"]
    
    exome.genotypeGVCFs(infile, outfile, genome, options)


@transform(genotypeGVCFs,
           regex(r"variants/all_samples.vcf"),
           r"variants/all_samples_excesshet.vcf.gz")
def filterExcessHets(infile, outfile):
    '''
    Hard filter on Excess Heterozygotes
    Should not make a difference on small cohorts unless related
    54.69 corresponds to a Z score of -4.5
    '''
    job_memory = GATK_MEMORY
    statement = '''
    gatk --java-options "-Xmx%(job_memory)s -Xms%(job_memory)s" 
    VariantFiltration
    -V %(infile)s
    --filter-expression "ExcessHet > 54.69"
    --filter-name ExcessHet
    -O %(outfile)s
    > %(outfile)s.log 2>&1
    '''
    P.run(statement)


@transform(filterExcessHets,
           regex(r"variants/all_samples_excesshet.vcf.gz"),
           r"variants/all_samples_sitesonly.vcf.gz")
def sitesOnlyVcf(infile,outfile):
    '''
    Create sites-only VCF with MakeSitesOnlyVcf
    '''
    job_memory = GATK_MEMORY
    statement = '''gatk MakeSitesOnlyVcf
    -I %(infile)s
    -O %(outfile)s
    > %(outfile)s.log 2>&1 '''
    P.run(statement)



###############################################################################
# SNP Recalibration

@transform(sitesOnlyVcf,
           regex(r"variants/all_samples_sitesonly.vcf.gz"),
           r"variants/all_samples.snp_vqsr.recal")
def variantRecalibratorSnps(infile, outfile):
    '''Create variant recalibration file'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    dbsnp = PARAMS["gatk_dbsnp"]
    job_threads = PARAMS["gatk_threads"]
    kgenomes = PARAMS["gatk_kgenomes"]
    hapmap = PARAMS["gatk_hapmap"]
    omni = PARAMS["gatk_omni"]
    mode = 'SNP'
    exome.variantRecalibrator(infile, outfile, genome, mode, dbsnp,
                                      kgenomes, hapmap, omni, GATK_MEMORY)


@follows(variantRecalibratorSnps)
@transform(filterExcessHets,
           regex(r"variants/all_samples_excesshet.vcf.gz"),
           add_inputs(r"variants/all_samples.snp_vqsr.recal",
                      r"variants/all_samples.snp_vqsr.tranches"),
           r"variants/all_samples.snp_vqsr.vcf.gz")
def applyVQSRSnps(infiles, outfile):
    '''Perform variant quality score recalibration using GATK '''
    vcf, recal, tranches = infiles
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    mode = 'SNP'
    exome.applyVQSR(vcf, recal, tranches,
                                            outfile, genome, mode, GATK_MEMORY)

###############################################################################
# Indel recalibration


@transform(sitesOnlyVcf,
           regex(r"variants/all_samples_sitesonly.vcf.gz"),
           r"variants/all_samples.indel_vqsr.recal")
def variantRecalibratorIndels(infile, outfile):
    '''Create variant recalibration file'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    job_threads = PARAMS["gatk_threads"]
    mills = PARAMS["gatk_mills"]
    dbsnp = PARAMS["gatk_dbsnp"]
    axiom = PARAMS["gatk_axiom"]
    mode = 'INDEL'
    exome.variantRecalibrator(infile, outfile, genome, mode,
                                      mills=mills, axiom=axiom, dbsnp=dbsnp, gatkmem=GATK_MEMORY)


@follows(variantRecalibratorIndels)
@transform(applyVQSRSnps,
           regex(r"variants/all_samples.snp_vqsr.vcf.gz"),
           add_inputs(r"variants/all_samples.indel_vqsr.recal",
                      r"variants/all_samples.indel_vqsr.tranches"),
           r"variants/all_samples.vqsr.vcf.gz")
def applyVQSRIndels(infiles, outfile):
    '''Perform variant quality score recalibration using GATK '''
    vcf, recal, tranches = infiles
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    mode = 'INDEL'
    exome.applyVQSR(vcf, recal, tranches,
                                            outfile, genome, mode, GATK_MEMORY)



###############################################################################
###############################################################################
###############################################################################
# Variant Annotation


###############################################################################
# Somalier - assessing relatedness of samples

@transform(applyVQSRIndels,
           regex(r"variants/all_samples.vqsr.vcf.gz"),
           r"variants/somalier/somalier.log")
def extractSomalier(infile,outfile):
    '''extract list of polymorphic sites using Somalier'''

    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    somalier = PARAMS['annotation_somalier_path']
    sites = PARAMS['annotation_somalier_sites']
    outdir = os.path.dirname(outfile)
    statement = '''%(somalier)s extract 
    -d %(outdir)s
    --sites %(sites)s
    -f %(genome)s 
    %(infile)s
    > %(outfile)s 2>&1
    '''
    P.run(statement)


@transform(extractSomalier,
           regex(r"variants/somalier/somalier.log"),
           r"variants/somalier/relate.log")
def relateSomalier(infile,outfile):
    '''Calculate Ancestry based on 1000G project'''

    somalier = PARAMS['annotation_somalier_path']
    outdir = os.path.dirname(outfile)
    statement ='''%(somalier)s relate
    %(outdir)s/*.somalier
    > %(outfile)s 2>&1'''
    P.run(statement)


@transform(extractSomalier,
           regex(r"variants/somalier/somalier.log"),
           r"variants/somalier/ancestry.log")
def ancestrySomalier(infile,outfile):
    '''Calculate Ancestry based on 1000G project'''

    somalier = PARAMS['annotation_somalier_path']
    labels = PARAMS['annotation_somalier_labels']
    reference = PARAMS['annotation_somalier_1kg']
    outdir = os.path.dirname(outfile)
    statement ='''%(somalier)s ancestry 
    --labels %(labels)s %(reference)s/*.somalier ++ 
    %(outdir)s/*.somalier
    > %(outfile)s 2>&1'''
    P.run(statement)


###############################################################################
# SnpSift

@transform(applyVQSRIndels,
           regex(r"variants/all_samples.vqsr.vcf.gz"),
           r"variants/all_samples.snpeff.vcf.gz")
def annotateVariantsSNPeff(infile, outfile):
    '''Annotate variants using SNPeff'''
    job_memory = PARAMS["annotation_memory"]
    job_threads = PARAMS["annotation_threads"]
    snpeff_genome = PARAMS["annotation_snpeff_genome"]
    config = PARAMS["annotation_snpeff_config"]
    outfile = P.snip(outfile,".gz")
    statement = ''' snpEff -Xmx%(job_memory)s 
    -c %(config)s 
    -v %(snpeff_genome)s 
    %(infile)s > %(outfile)s 2> %(outfile)s.log;
    bgzip %(outfile)s;
    tabix -p vcf %(outfile)s.gz'''  % locals()
    P.run(statement, job_memory=PARAMS["annotation_memory"])


@transform(annotateVariantsSNPeff,
           regex(r"variants/all_samples.snpeff.vcf.gz"),
           r"variants/all_samples.snpeff.table")
def vcfToTableSnpEff(infile, outfile):
    '''Converts vcf to tab-delimited file'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    columns = PARAMS["annotation_snpeff_to_table"]
    exome.vcfToTable(infile, outfile, genome, columns, GATK_MEMORY)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(vcfToTableSnpEff, regex(r"variants/(\S+).table"),
           r"variants/\1.table.load")
def loadTableSnpEff(infile, outfile):
    '''Load VCF annotations into database'''
    P.load(infile, outfile, options="--retry --ignore-empty")

###############################################################################
# SnpSift

@transform(annotateVariantsSNPeff,
           regex(r"variants/all_samples.snpeff.vcf.gz"),
           r"variants/all_samples.snpsift_1.vcf.gz")
def annotateVariantsDBNSFP(infile, outfile):
    '''Add annotations using SNPsift'''
    job_memory = "6G"
    job_threads = PARAMS["annotation_threads"]
    dbNSFP = PARAMS["annotation_dbnsfp"]
    dbNSFP_columns = PARAMS["annotation_dbnsfp_columns"]
    outfile = P.snip(outfile,".gz")
    if not dbNSFP_columns:
        annostring = ""
    else:
        annostring = "-f %s" % dbNSFP_columns

    statement = '''SnpSift dbnsfp 
    -db %(dbNSFP)s -v %(infile)s
    %(annostring)s >
    %(outfile)s 2> %(outfile)s.log;
    bgzip %(outfile)s;
    tabix -p vcf %(outfile)s.gz'''
    P.run(statement)


@transform(annotateVariantsDBNSFP,
           regex(r"variants/all_samples.snpsift_1.vcf.gz"),
           r"variants/all_samples.snpsift_2.vcf.gz")
def annotateVariantsClinvar(infile, outfile):
    '''Add annotations using SNPsift'''
    job_threads = PARAMS["annotation_threads"]
    clinvar = PARAMS["annotation_clinvar"]
    outfile = P.snip(outfile,".gz")
    exome.snpSift(infile,outfile,clinvar)


@transform(annotateVariantsClinvar,
           regex(r"variants/all_samples.snpsift_2.vcf.gz"),
           r"variants/all_samples.snpsift_3.vcf.gz")
def annotateVariantsGWASC(infile, outfile):
    '''Add annotations using SNPsift'''
    job_memory = "6G"
    job_threads = PARAMS["annotation_threads"]
    gwas_catalog = PARAMS["annotation_gwas_catalog"]
    outfile = P.snip(outfile,".gz")

    statement = """SnpSift gwasCat -db %(gwas_catalog)s
    %(infile)s > %(outfile)s  2> %(outfile)s.log;
    bgzip %(outfile)s;
    tabix -p vcf %(outfile)s.gz;"""
    P.run(statement)


@transform(annotateVariantsGWASC,
           regex(r"variants/all_samples.snpsift_3.vcf.gz"),
           r"variants/all_samples.snpsift_4.vcf.gz")
def annotateVariantsPhastcons(infile, outfile):
    '''Add annotations using SNPsift'''
    job_memory = "6G"
    job_threads = PARAMS["annotation_threads"]
    genome_index = "%s/%s.fa.fai" % (
        PARAMS['genome_dir'],
        PARAMS['genome'])
    phastcons = PARAMS["annotation_phastcons"]
    outfile = P.snip(outfile,".gz")

    statement = """ln -sf %(genome_index)s %(phastcons)s/genome.fai &&
    SnpSift phastCons %(phastcons)s %(infile)s >
    %(outfile)s  2> %(outfile)s.log;
    bgzip %(outfile)s;
    tabix -p vcf %(outfile)s.gz;"""
    P.run(statement)


@transform(annotateVariantsPhastcons,
           regex(r"variants/all_samples.snpsift_4.vcf.gz"),
           r"variants/all_samples.snpsift_5.vcf.gz")
def annotateVariants1000G(infile, outfile):
    '''Add annotations using SNPsift'''
    job_memory = "6G"
    job_threads = PARAMS["annotation_threads"]
    outfile = P.snip(outfile,".gz")

    vcfs = []
    for f in os.listdir(PARAMS["annotation_tgdir"]):
        if f.endswith(".vcf.gz"):
            vcfs.append("%s/%s" % (PARAMS['annotation_tgdir'], f))

    T = P.get_temp_filename(".")
    shutil.copy(infile, T)
    tempin = T
    tempout = P.get_temp_filename(".")

    for vcf in vcfs:
        statement = """SnpSift annotate
                       %(vcf)s
                       %(tempin)s > %(tempout)s &&
                       mv %(tempout)s %(tempin)s"""
        P.run(statement)

    shutil.move(tempin, outfile)
    statement = '''  bgzip %(outfile)s;
    tabix -p vcf %(outfile)s.gz''';
    P.run(statement)
 

@transform(annotateVariants1000G,
           regex(r"variants/all_samples.snpsift_5.vcf.gz"),
           r"variants/all_samples_dbsnp.snpsift_final.vcf.gz")
def annotateVariantsDBSNP(infile, outfile):
    '''Add annotations using SNPsift'''
    job_threads = PARAMS["annotation_threads"]
    dbsnp = PARAMS["annotation_dbsnp"]
    exome.snpSift(infile,outfile,dbsnp)
    P.run(statement)


@follows(annotateVariantsDBSNP)
def annotateVariantsSNPsift():
    pass


@transform(annotateVariantsDBSNP,
           regex(r"variants/all_samples.snpsift_final.vcf.gz"),
           r"variants/all_samples.vep.vcf")
def annotateVariantsVEP(infile, outfile):
    '''
    Adds annotations as specified in the pipeline.yml using Ensembl
    variant effect predictor (VEP).
    '''
    # infile - VCF
    # outfile - VCF with vep annotations
    job_memory = "6G"
    job_threads = 4

    VEP = PARAMS["annotation_vepannotators"].split(",")
    vep_annotators = PARAMS["annotation_vepannotators"]
    vep_path = PARAMS["annotation_veppath"]
    vep_cache = PARAMS["annotation_vepcache"]
    vep_species = PARAMS["annotation_vepspecies"]
    vep_assembly = PARAMS["annotation_vepassembly"]
    if len(vep_annotators) != 0:
        annostring = vep_annotators
        statement = '''perl %(vep_path)s/variant_effect_predictor.pl
                       --cache --dir %(vep_cache)s --vcf
                       --species %(vep_species)s
                       --fork 2
                       --assembly %(vep_assembly)s --input_file %(infile)s
                       --output_file %(outfile)s --force_overwrite
                       %(annostring)s --offline;'''
        P.run(statement)


@follows(mkdir("variant_tables"))
@transform(RemoveDuplicatesSample, regex(r"gatk/(.*).dedup2.bam"),
           add_inputs(annotateVariantsVEP), r"variant_tables/\1.tsv")
def makeAnnotationsTables(infiles, outfile):
    '''
    Converts the multi sample vcf generated with Haplotype caller into
    a single table for each sample contain only positions called as variants
    in that sample and all columns from the vcf.
    '''
    bamname = infiles[0]
    inputvcf = infiles[1]
    TF = P.get_temp_filename(".")
    samplename = bamname.replace(".dedup2.bam",
                                 ".bam").replace("gatk/", "")

    statement = '''bcftools view -h %(inputvcf)s |
                   awk -F '=|,' '$1=="##INFO" || $1=="##FORMAT"
                   {printf("%%s\\t%%s\\n", $3, $9)}'
                   | sed 's/>//g' > %(TF)s'''
    P.run(statement)
    cols = []
    colds = []
    cols2 = []
    for line in iotools.open_file(TF).readlines():
        if line.split("\t")[0] != "Samples":
            cols.append("[%%%s]" % line.split("\t")[0])
            colds.append(line.split("\t")[1].strip().replace(" ", "_"))
            cols2.append(line.split("\t")[0])
    l1 = "\t".join(cols2)
    l2 = "\t".join(colds)
    out = open(outfile, "w")
    out.write('''CHROM\tPOS\tQUAL\tID\tFILTER\tREF1\tALT\tGT\t%s\n\
                 chromosome\tposition\tquality\tid\tfilter\tref\talt\t\
                 genotype\t%s\n''' % (l1, l2))
    out.close()
    cstring = "\\t".join(cols)
    cstring = "%CHROM\\t%POS\\t%QUAL\\t%ID\\t%FILTER\\t%REF\\t\
               %ALT\\t[%GT]\\t" + cstring
    if PARAMS['test'] == 1:
        statement = '''bcftools query
                   -f '%(cstring)s\\n'
                   -i 'FILTER=="PASS" && GT!="0/0" && GT!="./."'
                   %(inputvcf)s >> %(outfile)s'''
    else:
        statement = '''bcftools query -s %(samplename)s
                   -f '%(cstring)s\\n'
                   -i 'FILTER=="PASS" && GT!="0/0" && GT!="./."'
                   %(inputvcf)s >> %(outfile)s'''
    P.run(statement)
        

###############################################################################
# Quality Filtering #

@follows(mkdir("variant_tables_highqual"))
@transform(makeAnnotationsTables, regex("variant_tables/(.*).tsv"),
           [r'variant_tables_highqual/\1.tsv',
            r'variant_tables_highqual/\1_failed.tsv'])
def qualityFilterVariants(infile, outfiles):
    '''
    Filter variants based on quality.  Columns to filter on and
    how they should be filtered can be specified in the pipeline.yml.
    Currently only implemented to filter numeric columns.  "." is assumed
    to mean pass.
    '''
    qualstring = PARAMS['filtering_quality']
    qualfilter = PARAMS['filtering_quality_ft']
    exome.filterQuality(infile, qualstring, qualfilter,
                                outfiles, submit=True)


@follows(mkdir("variant_tables_rare"))
@transform(qualityFilterVariants, regex("variant_tables_highqual/(.*).tsv"),
           [r'variant_tables_rare/\1.tsv',
            r'variant_tables_rare/\1_failed.tsv'])
def rarityFilterVariants(infiles, outfiles):
    '''
    Filter out variants which are common in any of the exac or other
    population datasets as specified in the pipeline.yml.
    '''
    infile = infiles[0]
    thresh = PARAMS['filtering_rarethresh']
    freqs = PARAMS['filtering_freqs']
    exac = PARAMS['filtering_exac']
    exome.filterRarity(infile, exac, freqs, thresh, outfiles,
                               submit=True)


@follows(mkdir("variant_tables_damaging"))
@transform(rarityFilterVariants, regex("variant_tables_rare/(.*).tsv"),
           [r'variant_tables_damaging/\1.tsv',
            r'variant_tables_damaging/\1_failed.tsv'])
def damageFilterVariants(infiles, outfiles):
    '''
    Filter variants which have not been assessed as damaging by any
    of the specified tools.
    Tools and thresholds can be specified in the pipeline.yml.

    Does not account for multiple alt alleles - if any ALT allele has
    been assessed as damaging with any tool the variant is kept,
    regardless of if this is the allele called in the sample.

    '''
    infile = infiles[0]
    exome.filterDamage(infile, PARAMS['filtering_damage'], outfiles,
                               submit=True)


@follows(mkdir("variant_tables_family"))
@transform(damageFilterVariants, regex("variant_tables_damaging/(.*).tsv"),
           [r'variant_tables_family/\1.tsv',
            r'variant_tables_family/\1_failed.tsv'])
def familyFilterVariants(infiles, outfiles):
    '''
    Filter variants according to the output of calculateFamily -
    only variants shared by both members of a family will be kept.
    BROKEN - NEED TO ADD FAMILY INPUTS
    '''
    if len(matches) > 1:
        infile = infiles[0][0]

        infilenam = infile.split("/")[-1]
        infilestem = "/".join(infile.split("/")[:-1])

        # figure out who is related to who
        families = [line.strip().split("\t")[:2]
                    for line in iotools.open_file(infiles[1]).readlines()]
        infam = [line[0] for line in families] + [line[1] for line in families]

        # no relatives - copy the input file to the output file and generate
        # a blank "failed" file
        if infilenam not in infam or PARAMS['filtering_family'] == 0:
            shutil.copy(infile, outfiles[0])
            o = iotools.open_file(outfiles[1], "w")
            o.close()
        else:
            for line in families:
                if infilenam in line:
                    i = line.index(infilenam)
                    if i == 0:
                        infile2 = "%s/%s" % (infilestem, line[1])
                    else:
                        infile2 = "%s/%s" % (infilestem, line[0])
            exome.filterFamily(infile, infile2, outfiles)
    else:
        infile = infiles[0][0]
        shutil.copy(infile, outfiles[0])
        out = iotools.open_file(outfiles[1], "w")
        out.close()


@follows(familyFilterVariants)
def filter():
    pass


@follows(mkdir("genelists"))
@transform(familyFilterVariants, regex("variant_tables_family/(.*).tsv"),
           [r'genelists/\1.tsv',
            r'genelists/\1.bed'])
def makeGeneLists(infiles, outfiles):
    infile = infiles[0]
    outfile = outfiles[1]
    genes = PARAMS['general_geneset']

    # four %s because they need to be escaped in generating the statement
    # then again when submitting the P.run(statement)
    statement = '''awk 'NR > 2 {printf("%%%%s\\t%%%%s\\t%%%%s\\n",\
    $1, $2, $2 + 1)}'\
    %(infile)s |\
    bedtools intersect -wo -a stdin -b %(genes)s > %(outfile)s''' % locals()
    P.run(statement)

    geneids = set()
    with iotools.open_file(outfile) as inp:
        for line in inp:
            line = line.strip().split("\t")
            details = line[11].split(";")
            for detail in details:
                r = re.search('gene_id', detail)
                if r:
                    geneid = detail.split(" ")[-1]
                    geneids.add(geneid.replace("\"", ""))
    out = iotools.open_file(outfiles[0], "w")
    for geneid in geneids:
        out.write("%s\n" % geneid)
    out.close()


@follows(mkdir('final_variant_tables'))
@collate((makeGeneLists, familyFilterVariants), regex('(.*)/(.*).tsv'),
         r'final_variant_tables/\2.tsv')
def finalVariantTables(infiles, outfile):
    genes = infiles[0]
    variants = infiles[1][0]
    cols = PARAMS['filtering_columns'].split(",")
    exome.CleanVariantTables(genes, variants, cols, outfile,
                                     submit=True)

###############################################################################
###############################################################################
###############################################################################
# Genes of interest


@active_if(PARAMS["annotation_add_genes_of_interest"] == 1)
@transform((annotateVariantsVEP),
           regex(r"variants/all_samples.snpsift.vcf"),
           r"variants/all_samples.genes.vcf")
def findGenes(infile, outfile):
    '''Adds expression "GENE_OF_INTEREST" to the FILTER column of the vcf
    if variant is within a gene of interest as defined in the ini
    file'''
    geneList = P.as_list(PARAMS["annotation_genes_of_interest"])
    expression = '\'||SNPEFF_GENE_NAME==\''.join(geneList)
    statement = '''GenomeAnalysisTK -T VariantFiltration
                   -R %%(genome_dir)s/%%(genome)s.fa
                   --variant %(infile)s
                   --filterExpression "SNPEFF_GENE_NAME=='%(expression)s'"
                   --filterName "GENE_OF_INTEREST" -o %(outfile)s''' % locals()
    P.run(statement)

###############################################################################
###############################################################################
###############################################################################
# Tabulation


TABULATION_INPUT = {'': annotateVariantsVEP,
                    0: annotateVariantsVEP,
                    1: findGenes}


@transform(TABULATION_INPUT[PARAMS.get("annotation_add_genes_of_interest", 0)],
           regex(r"variants/all_samples.(snpsift|genes).vcf"),
           r"variants/all_samples.snpsift.table")
def vcfToTable(infile, outfile):
    '''Converts vcf to tab-delimited file'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    exome.vcfToTable(infile, outfile, genome, columns, GATK_MEMORY)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(vcfToTable, regex(r"variants/all_samples.snpsift.table"),
           r"variants/all_samples.snpsift.table.load")
def loadVariantAnnotation(infile, outfile):
    '''Load VCF annotations into database'''
    P.load(infile, outfile, options="--retry --ignore-empty")

##############################################################################
###############################################################################
###############################################################################
# Confirm parentage (do novo trios only)


@follows(loadVariantAnnotation)
@transform("*Trio*.ped",
           regex(r"(\S*Trio\S+).ped"),
           add_inputs(r"variants/all_samples.snpsift.vcf"),
           r"variants/\1.parentage")
def confirmParentage(infiles, outfile):
    '''Filter variants according to autosomal recessive disease model'''
    pedfile, infile = infiles
    pedigree = csv.DictReader(open(pedfile), delimiter='\t', fieldnames=[
        'family', 'sample', 'father', 'mother', 'sex', 'status'])
    trio = P.snip(os.path.basename(pedfile), ".ped")
    trio = trio.replace(".", "_").replace("-", "_")
    database = PARAMS["database_name"]
    proband = None
    mother = None
    father = None
    for row in pedigree:
        if row['status'] == '2':
            proband = row['sample'].replace(".", "_").replace("-", "_")
            mother = row['mother'].replace(".", "_").replace("-", "_")
            father = row['father'].replace(".", "_").replace("-", "_")
    E.info("proband:" + proband)
    E.info("mother:" + mother)
    E.info("father:" + father)

    query = '''SELECT '%(trio)s' as family,
               hom_nonref_trio/(hom_nonref_parents+0.0) as hom_nonref_conc,
               child_mat_het_nonref/(maternal_hom_nonref+0.0) as maternal_conc,
               child_pat_het_nonref/(paternal_hom_nonref+0.0) as paternal_conc,
                *
               FROM
               (SELECT count(*) as hom_nonref_trio
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND %(mother)s_PL LIKE '%%%%,0'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND  %(father)s_PL LIKE '%%%%,0'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND %(proband)s_PL LIKE '%%%%,0'
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as hom_nonref_parents
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND %(mother)s_PL LIKE '%%%%,0'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND  %(father)s_PL LIKE '%%%%,0'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as maternal_hom_nonref
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND %(mother)s_PL LIKE '%%%%,0'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND  %(father)s_PL LIKE '0,%%%%'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as child_mat_het_nonref
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND %(mother)s_PL LIKE '%%%%,0'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND  %(father)s_PL LIKE '0,%%%%'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND %(proband)s_PL LIKE '%%%%,0,%%%%'
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as paternal_hom_nonref
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND  %(father)s_PL LIKE '%%%%,0'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND  %(mother)s_PL LIKE '0,%%%%'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as child_pat_het_nonref
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND  %(father)s_PL LIKE '%%%%,0'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND  %(mother)s_PL LIKE '0,%%%%'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND %(proband)s_PL LIKE '%%%%,0,%%%%'
               AND cast(%(proband)s_DP as INTEGER) > 10)''' % locals()
    statement = '''sqlite3 %(database)s "%(query)s" > %(outfile)s
                   2> %(outfile)s.log''' % locals()
    P.run(statement)

###############################################################################
###############################################################################
###############################################################################
# De novos


@follows(loadVariantAnnotation)
@transform("*Trio*.ped",
           regex(r"(\S*Trio\S+).ped"),
           add_inputs(r"variants/all_samples.snpsift.vcf"),
           r"variants/\1.filtered.vcf")
def deNovoVariants(infiles, outfile):
    '''Filter de novo variants based on provided jexl expression'''
    job_memory = GATK_MEMORY
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    pedfile, infile = infiles
    pedigree = csv.DictReader(
        iotools.open_file(pedfile), delimiter='\t', fieldnames=[
            'family', 'sample', 'father', 'mother', 'sex', 'status'])
    for row in pedigree:
        if row['status'] == '2':
            father = row['father']
            mother = row['mother']
            child = row['sample']
    select = '''vc.getGenotype("%(father)s").getDP()>=10&&vc.getGenotype("%(mother)s").getDP()>=10&&vc.getGenotype("%(child)s").getPL().0>20&&vc.getGenotype("%(child)s").getPL().1==0&&vc.getGenotype("%(child)s").getPL().2>0&&vc.getGenotype("%(father)s").getPL().0==0&&vc.getGenotype("%(father)s").getPL().1>20&&vc.getGenotype("%(father)s").getPL().2>20&&vc.getGenotype("%(mother)s").getPL().0==0&&vc.getGenotype("%(mother)s").getPL().1>20&&vc.getGenotype("%(mother)s").getPL().2>20&&vc.getGenotype("%(child)s").getAD().1>=3&&((vc.getGenotype("%(child)s").getAD().1)/(vc.getGenotype("%(child)s").getDP().floatValue()))>=0.25&&(vc.getGenotype("%(father)s").getAD().1==0||(vc.getGenotype("%(father)s").getAD().1>0&&((vc.getGenotype("%(father)s").getAD().1)/(vc.getGenotype("%(father)s").getDP().floatValue()))<0.05))&&(vc.getGenotype("%(mother)s").getAD().1==0||(vc.getGenotype("%(mother)s").getAD().1>0&&((vc.getGenotype("%(mother)s").getAD().1)/(vc.getGenotype("%(mother)s").getDP().floatValue()))<0.05))&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals()
    exome.selectVariants(infile, outfile, genome, select)


@transform(deNovoVariants,
           regex(r"variants/(\S+).filtered.vcf"),
           r"variants/\1.filtered.table")
def tabulateDeNovos(infile, outfile):
    '''Tabulate de novo variants'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    exome.vcfToTable(infile, outfile, genome, columns, GATK_MEMORY)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(tabulateDeNovos,
           regex(r"variants/(\S+).filtered.table"),
           r"variants/\1.filtered.table.load")
def loadDeNovos(infile, outfile):
    '''Load de novos into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################


@follows(loadVariantAnnotation)
@transform("*Trio*.ped",
           regex(r"(\S*Trio\S+).ped"),
           add_inputs(r"variants/all_samples.snpsift.vcf"),
           r"variants/\1.denovos.vcf")
def lowerStringencyDeNovos(infiles, outfile):
    '''Filter lower stringency de novo variants based on provided jexl
    expression'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    pedfile, infile = infiles
    pedigree = csv.DictReader(
        iotools.open_file(pedfile), delimiter='\t', fieldnames=[
            'family', 'sample', 'father', 'mother', 'sex', 'status'])
    for row in pedigree:
        if row['status'] == '2':
            father = row['father']
            mother = row['mother']
            child = row['sample']
    select = '''vc.getGenotype("%(child)s").getPL().1==0&&vc.getGenotype("%(father)s").getPL().0==0&&vc.getGenotype("%(mother)s").getPL().0==0&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals(
    )
    exome.selectVariants(infile, outfile, genome, select)


@transform(lowerStringencyDeNovos,
           regex(r"variants/(\S+).denovos.vcf"),
           r"variants/\1.denovos.table")
def tabulateLowerStringencyDeNovos(infile, outfile):
    '''Tabulate lower stringency de novo variants'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    exome.vcfToTable(infile, outfile, genome, columns)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(tabulateLowerStringencyDeNovos,
           regex(r"variants/(\S+).denovos.table"),
           r"variants/\1.denovos.table.load")
def loadLowerStringencyDeNovos(infile, outfile):
    '''Load lower stringency de novos into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# Dominant


@follows(loadVariantAnnotation)
@transform("*Multiplex*.ped",
           regex(r"(\S*Multiplex\S+).ped"),
           add_inputs(r"variants/all_samples.snpsift.vcf"),
           r"variants/\1.dominant.vcf")
def dominantVariants(infiles, outfile):
    '''Filter variants according to autosomal dominant disease model'''
    pedfile, infile = infiles
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    pedigree = csv.DictReader(open(pedfile), delimiter='\t',
                              fieldnames=['family', 'sample', 'father',
                                          'mother', 'sex', 'status'])
    affecteds = []
    unaffecteds = []
    for row in pedigree:
        if row['status'] == '2':
            affecteds += [row['sample']]
        if row['status'] == '1':
            unaffecteds += [row['sample']]
    affecteds_exp = '").getPL().1==0&&vc.getGenotype("'.join(affecteds)
    affecteds_exp2 = '").getAD().1>=1&&vc.getGenotype("'.join(affecteds)
    if len(unaffecteds) == 0:
        unaffecteds_exp = ''
    else:
        unaffecteds_exp = '&&vc.getGenotype("' + \
            ('").isHomRef()&&vc.getGenotype("'.join(unaffecteds)) + \
            '").isHomRef()'
    # for some weird reason the 1000G filter doesn't work on these particular
    # files - will add later when I've figured out what's wrong
    # currently 1000G filter is performed at the report stage (not in csvdb)
    select = '''vc.getGenotype("%(affecteds_exp)s").getPL().1==0&&vc.getGenotype("%(affecteds_exp2)s").getAD().1>=1%(unaffecteds_exp)s&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals()
    exome.selectVariants(infile, outfile, genome, select)

###############################################################################


@transform(dominantVariants,
           regex(r"variants/(\S+).dominant.vcf"),
           r"variants/\1.dominant.table")
def tabulateDoms(infile, outfile):
    '''Tabulate dominant disease candidate variants'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    exome.vcfToTable(infile, outfile, genome, columns)

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(tabulateDoms, regex(r"variants/(\S+).dominant.table"),
           r"variants/\1.dominant.table.load")
def loadDoms(infile, outfile):
    '''Load dominant disease candidates into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# Recessive


@follows(loadVariantAnnotation)
@transform("*.ped",
           regex(
               r"(\S*Trio\S+|\S*Multiplex\S+).ped"),
           add_inputs(r"variants/all_samples.snpsift.vcf"),
           r"variants/\1.recessive.vcf")
def recessiveVariants(infiles, outfile):
    '''Filter variants according to autosomal recessive disease model'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    pedfile, infile = infiles
    pedigree = csv.DictReader(open(pedfile), delimiter='\t',
                              fieldnames=['family', 'sample', 'father',
                                          'mother', 'sex', 'status'])
    affecteds = []
    parents = []
    unaffecteds = []
    for row in pedigree:
        if row['status'] == '2':
            affecteds += [row['sample']]
            if row['father'] != '0':
                parents += [row['father']]
            if row['mother'] != '0':
                parents += [row['mother']]
        elif row['status'] == '1' and row['sample'] not in parents:
            unaffecteds += [row['sample']]
    affecteds_exp = '").getPL().2==0&&vc.getGenotype("'.join(affecteds)
    if len(unaffecteds) == 0:
        unaffecteds_exp = ''
    else:
        unaffecteds_exp = '&&vc.getGenotype("' + \
            ('").getPL().2!=0&&vc.getGenotype("'.join(unaffecteds)) + \
            '").getPL().2!=0'
    if len(parents) == 0:
        parents_exp = ''
    else:
        parents_exp = '&&vc.getGenotype("' + \
            ('").getPL().1==0&&vc.getGenotype("'.join(parents)) + \
            '").getPL().1==0'
    # need a way of specifying that other unaffecteds eg. sibs can't be
    # homozygous for alt allele
    select = '''vc.getGenotype("%(affecteds_exp)s").getPL().2==0%(unaffecteds_exp)s%(parents_exp)s&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals(
    )
    exome.selectVariants(infile, outfile, genome, select)

###############################################################################


@transform(recessiveVariants,
           regex(r"variants/(\S+).recessive.vcf"),
           r"variants/\1.recessive.table")
def tabulateRecs(infile, outfile):
    '''Tabulate potential homozygous recessive disease variants'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    exome.vcfToTable(infile, outfile, genome, columns)

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(tabulateRecs, regex(r"variants/(\S+).recessive.table"),
           r"variants/\1.recessive.table.load")
def loadRecs(infile, outfile):
    '''Load homozygous recessive disease candidates into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# X-linked


@follows(loadVariantAnnotation)
@transform("*.ped", regex(r"(\S*Trio\S+|\S*Multiplex\S+).ped"),
           add_inputs(r"variants/all_samples.snpsift.vcf"),
           r"variants/\1.xlinked.vcf")
def xlinkedVariants(infiles, outfile):
    '''Find maternally inherited X chromosome variants in male patients'''
    track = P.snip(os.path.basename(outfile), ".vcf")
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    pedfile, infile = infiles
    pedigree = csv.DictReader(open(pedfile), delimiter='\t',
                              fieldnames=['family', 'sample', 'father',
                                          'mother', 'sex', 'status'])
    affecteds = []
    mothers = []
    male_unaffecteds = []
    female_unaffecteds = []
    for row in pedigree:
        if row['status'] == '2' and row['mother'] != '0':
            if row['mother'] not in mothers:
                mothers += [row['mother']]
            affecteds += [row['sample']]
        elif row['status'] == '1' and row['sex'] == '1':
            male_unaffecteds += [row['sample']]
        elif row['status'] == '1' and row['sex'] != '1' and row['sample'] not in mothers:
            female_unaffecteds += [row['sample']]
    if len(affecteds) == 0:
        affecteds_exp = ''
    else:
        affecteds_exp = '").isHomRef()&&!vc.getGenotype("'.join(affecteds)
#        affecteds_exp = '&&!vc.getGenotype("' + \
#        ('").isHomRef()&&!vc.getGenotype("'.join(affecteds)) + \
#        '").isHomRef()'
    if len(mothers) == 0:
        mothers_exp = ''
    else:
        mothers_exp = '&&vc.getGenotype("' + \
                      ('").getPL().1==0&&vc.getGenotype("'.join(mothers)) + \
                      '").getPL().1==0'
    if len(male_unaffecteds) == 0:
        male_unaffecteds_exp = ''
    else:
        male_unaffecteds_exp = '&&vc.getGenotype("' + \
                               ('").isHomRef()&&vc.getGenotype("'.join(
                                   male_unaffecteds)) + \
                               '").isHomRef()'
    if len(female_unaffecteds) == 0:
        female_unaffecteds_exp = ''
    else:
        female_unaffecteds_exp = '&&vc.getGenotype("' + \
                                 ('").getPL().2!=0&&vc.getGenotype("'.join(
                                     female_unaffecteds)) + \
            '").getPL().2!=0'
    select = '''CHROM=="X"&&!vc.getGenotype("%(affecteds_exp)s").isHomRef()%(mothers_exp)s%(male_unaffecteds_exp)s%(female_unaffecteds_exp)s&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals()
    exome.selectVariants(infile, outfile, genome, select)

###############################################################################


@transform(xlinkedVariants,
           regex(r"variants/(\S+).xlinked.vcf"),
           r"variants/\1.xlinked.table")
def tabulateXs(infile, outfile):
    '''Tabulate potential X-linked disease variants'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    exome.vcfToTable(infile, outfile, genome, columns)

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(tabulateXs, regex(r"variants/(\S+).xlinked.table"),
           r"variants/\1.xlinked.table.load")
def loadXs(infile, outfile):
    '''Load X-linked disease candidates into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# Compound hets


# why does this not work from snpsift VCF file? DS
# @transform(annotateVariantsSNPsift,
#    regex(r"variants/(\S*Trio\S+|\S*Multiplex\S+).haplotypeCaller.snpsift.vcf"),
@transform(annotateVariantsSNPeff,
           regex(
               r"variants/all_samples.snpeff.vcf"),
           add_inputs(r"all_samples.ped", r"gatk/all_samples.list"),
           r"variants/all_samples.phased.vcf")
def phasing(infiles, outfile):
    '''phase variants with GATK'''
    infile, pedfile, bamlist = infiles
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    statement = '''GenomeAnalysisTK -T PhaseByTransmission
                   -R %(genome)s
                   -V %(infile)s
                   -ped %(pedfile)s
                   -mvf %(infile)s.mvf
                   -o %(outfile)s ;''' % locals()
    P.run(statement)

###############################################################################


@transform(phasing, regex(r"variants/all_samples.phased.vcf"),
           add_inputs(r"gatk/all_samples.list"),
           r"variants/all_samples.rbp.vcf")
def readbackedphasing(infiles, outfile):
    '''phase variants with ReadBackedPhasing'''
    job_memory = "32G"
    infile, bamlist = infiles
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    statement = '''GenomeAnalysisTK -T ReadBackedPhasing -nt 4
                   -R %(genome)s
                   -I %(bamlist)s
                   -V %(infile)s
                   -o %(outfile)s; ''' % locals()
    P.run(statement)

###############################################################################


@follows(readbackedphasing)
@transform("*.ped",
           regex(r"(\S*Multiplex\S+|\S*Trio\S+).ped"),
           add_inputs(r"no_multiallelic_all_samples.rbp.vcf",
                      r"all_samples.ped"),
           r"variants/\1.compound_hets.table")
def compoundHets(infiles, outfile):
    '''Identify potentially pathogenic compound heterozygous variants
    (phasing with GATK followed by compound het calling using Gemini'''
    family, infile, pedfile = infiles
    family_id = P.snip(os.path.basename(family), ".ped")
    statement = '''gemini load -v %(infile)s
                   -p %(pedfile)s -t snpEff %(family_id)s.db &&
                   gemini comp_hets
                   --families %(family_id)s
                   --columns "chrom, start, end, ref, alt, codon_change, gene, qual, depth"
                   --filter
                   "(impact_severity = 'HIGH' OR impact_severity = 'MED')
                   AND (in_esp = 0 OR aaf_esp_all < 0.01)
                   AND (in_1kg = 0 OR aaf_1kg_all < 0.01)"
                   %(family_id)s.db > %(outfile)s;'''
    # rm -f %(family_id)s.db'''
    P.run(statement)

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(compoundHets, regex(r"variants/(\S+).compound_hets.table"),
           r"variants/\1.compound_hets.table.load")
def loadCompoundHets(infile, outfile):
    '''Load compound heterozygous variants into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# coverage over candidate genes

@merge(RemoveDuplicatesSample, "gatk/all_samples.list")
def listOfBAMs(infiles, outfile):
    '''generates a file containing a list of BAMs for use in VQSR'''
    with iotools.open_file(outfile, "w") as outf:
        for infile in infiles:
            outf.write(infile + '\n')

@active_if(PARAMS["coverage_calculate"] == 1)
@transform(listOfBAMs, regex(r"gatk/all_samples.list"), r"candidate.sample_interval_summary")
def candidateCoverage(infile, outfile):
    '''Calculate coverage over exons of candidate genes'''
    all_exons = PARAMS["coverage_all_exons"]
    candidates = PARAMS["coverage_candidates"]
    candidates = candidates.replace(",", " -e ")
    genome = PARAMS["genome_dir"] + "/" + PARAMS["genome"] + ".fa"
    threshold = PARAMS["coverage_threshold"]
    statement = '''zcat %(all_exons)s | grep -e %(candidates)s
                   | awk '{print $1 ":" $4 "-" $5}' - | sed 's/chr//' - >
                   candidate.interval_list ; ''' % locals()
    statement += '''zcat %(all_exons)s | grep -e %(candidates)s
                    | awk '{print $16}' - | sed 's/"//g;s/;//g' - >
                    candidate_gene_names.txt ;''' % locals()
    statement += '''GenomeAnalysisTK -T DepthOfCoverage -R %(genome)s
                    -o candidate -I %(infile)s -ct %(threshold)s -L candidate.interval_list ;''' % locals()
    P.run(statement)

###############################################################################


@active_if(PARAMS["coverage_calculate"] == 1)
@transform(candidateCoverage, regex(r"candidate.sample_interval_summary"), r"candidate_coverage_plot.pdf")
def candidateCoveragePlots(infile, outfile):
    '''Produce plots of coverage'''
    rscript = PARAMS["coverage_rscript"]
    threshold = PARAMS["coverage_threshold"]
    statement = '''Rscript %(rscript)s %(infile)s candidate_gene_names.txt %(threshold)s %(outfile)s ;'''
    P.run(statement)

###############################################################################
###############################################################################
###############################################################################
# vcf statistics


@transform((annotateVariantsVEP), regex(
    r"variants/all_samples.snpsift.vcf"), r"variants/all_samples.vcfstats")
def buildVCFstats(infile, outfile):
    '''Calculate statistics on VCF file'''
    statement = '''vcf-stats %(infile)s > %(outfile)s
                   2>>%(outfile)s.log;''' % locals()
    P.run(statement)

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildVCFstats, "vcf_stats.load")
def loadVCFstats(infiles, outfile):
    '''Import variant statistics into SQLite'''
    filenames = " ".join(infiles)
    tablename = P.to_table(outfile)
    E.info("Loading vcf stats...")
    statement = '''cgat vcfstats2db %(filenames)s >>
                   %(outfile)s && '''
    statement += '''cat vcfstats.txt | cgat csv2db
                    %(csv2db_options)s --allow-empty-file --add-index=track
                    --table=vcf_stats >> %(outfile)s &&'''
    statement += '''cat sharedstats.txt | cgat csv2db
                    %(csv2db_options)s --allow-empty-file --add-index=track
                    --table=vcf_shared_stats >> %(outfile)s &&'''
    statement += '''cat indelstats.txt | cgat csv2db
                    %(csv2db_options)s --allow-empty-file --add-index=track
                    --table=indel_stats >> %(outfile)s &&'''
    statement += '''cat snpstats.txt | cgat csv2db
                    %(csv2db_options)s --allow-empty-file --add-index=track
                    --table=snp_stats >> %(outfile)s &&'''
    P.run(statement)

###############################################################################
###############################################################################
###############################################################################
# Targets


@follows(loadVariantAnnotation,
         finalVariantTables)
def testFromVariantRecal():
    pass


@follows(loadROI,
         loadROI2Gene,
         loadSamples)
def loadMetadata():
    pass


@follows(mapReads,
         loadPicardAlignStats)
def mapping_tasks():
    pass


@follows(GATKBaseRecal,
         loadPicardDuplicateStatsLane,
         RemoveDuplicatesSample,
         loadCoverageStats)
def gatk():
    pass


@follows(haplotypeCaller,
         genotypeGVCFs)
def callVariants():
    pass


@follows(loadTableSnpEff,
         listOfBAMs,
         loadVariantAnnotation,
         finalVariantTables)
def annotation():
    pass


@follows(loadDeNovos)
def denovo():
    pass


@follows(loadLowerStringencyDeNovos)
def denovo2():
    pass


@follows(loadDoms)
def dominant():
    pass


@follows(loadRecs)
def recessive():
    pass


@follows(loadXs)
def xlinked():
    pass


@follows(loadCompoundHets)
def compoundHet():
    pass


@follows(denovo,
         dominant,
         recessive,
         compoundHet)
def filtering():
    pass


@follows(buildVCFstats,
         loadVCFstats)
def vcfstats():
    pass


@follows(mapping_tasks,
         gatk,
         callVariants,
         annotation,
         filtering,
         makeAnnotationsTables)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
