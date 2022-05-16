
'''
PEDIGREE CODE FROM PIPELINE EXOME
REMOVED FOR CLARITY
'''


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
