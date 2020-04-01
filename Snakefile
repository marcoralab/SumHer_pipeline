'''Snakefile for sex stratified GWAS
   Version 0.2'''

from scripts.parse_config import parser, getcom, getloads, subs
from scripts.parse_config import proc_cohorts, proc_covars
from getpass import getuser
import glob
import os
import socket

isMinerva = "hpc.mssm.edu" in socket.getfqdn()

configfile: "config.yaml"
shell.executable("/bin/bash")

BPLINK = ["bed", "bim", "fam"]
CHROM = parser(config)


def gsc(wildcards, k):
    return study_cohorts[wildcards.sample][k]


study_cohorts = {k: proc_cohorts(config['studies'][k], config)
                 for k in config['studies'].keys()}
study_covars = {k: proc_covars(k, config, study_cohorts)
                for k in config['studies'].keys()}

if isMinerva:
    anacondapath = sys.exec_prefix + "/bin"
    shell.prefix(". ~/.bashrc; PATH={}:$PATH; ".format(anacondapath))
    tempdir = "/sc/orga/scratch/{}/temp/".format(getuser())
else:
    import tempfile
    tempdir = tempfile.gettempdir() + '/'

com = getcom({'plink': 'plink --keep-allele-order', 'plink2': 'plink',
              'plink3': 'plink --keep-allele-order',
              'bcftools': 'bcftools', 'R': 'Rscript', 'R2': 'R', 'METAL': 'metal'},
             {'plink': '--memory 3800', 'plink3': '--memory 15800'})
loads = getloads(config,
                 {'plink': 'plink', 'bcftools': 'bcftools', 'METAL': 'metal',
                  'R': ['R', 'pandoc', 'udunits']},
                 lsfextra={'R': 'RSTUDIO_PANDOC=$(which pandoc)'})


def filldir(string):
    if '{basedir}' in string:
        return string.format(basedir=config['base_dir'])
    return string


infiles = {k: {kk: filldir(vv) for kk, vv in v.items() if isinstance(vv, str)}
           for k, v in config['studies'].items()}

maf = config['maf']

COV = list(list(study_covars.values())[0].keys())
MAF = list(maf.keys())

MANEXT = ["manhattan.png", "qq.png"]

GWASMAN = ["empP", "filtered"] if config['perm'] else ["filtered"]


def maybetemp(x):
    return x if config['keepgenos'] else temp(x)


OUTFILES = expand(
    "GWAS/cov-{cov}.maf-{maf}/{sample}.{allsex}_{Ptype}.assoc.{tt}.{ext}",
    cov=COV, maf=MAF, Ptype=GWASMAN, ext=MANEXT, tt=config['test'],
    sample='ADGC', allsex=['male', 'female', 'interaction'])

if config['keepgenos']:
    OUTFILES += expand(
        "filtered/cov-{cov}.maf-{maf}/{sample}.{allsex}.chr{chrom}.{ext}",
        cov=COV, maf=MAF, sample='ADGC', chrom=CHROM, ext=BPLINK,
        allsex=['male', 'female', 'interaction'])

rule all:
    input: OUTFILES

rule make_samplist:
    input: lambda wildcards: infiles[wildcards.sample]["pheno"]
    output:
        pheno = "phenotypes/{sample}.pheno",
        ikeep = "{sample}.chosen.ikeep"
    params:
        default_cht = lambda wildcards: gsc(wildcards, 'DEFAULT_COHORT'),
        filts = lambda wildcards: subs(config[wildcards.sample]['filter']),
        ec = lambda wildcards: gsc(wildcards, 'EXCLUDECOHORTS')
    shell:
        """
{loads[R]}
{com[R]} scripts/make_plink_phenos.R {input} {output.pheno} {output.ikeep} \
 {params.filts} {params.default_cht} {params.ec}
"""

rule filter_adgc:
    input:
        plink = lambda wildcards: expand(infiles[wildcards.sample]["geno"] + ".{ext}", ext=BPLINK),
        keep = "{sample}.chosen.ikeep"
    output: temp(expand("genotypes/{{sample}}.{{sex}}.maf-{{maf}}.chr{{chrom}}.{ext}", ext=BPLINK))
    params:
        i = lambda wildcards: infiles[wildcards.sample]["geno"],
        o = "genotypes/{sample}.{sex}.maf-{maf}.chr{chrom}",
        maf = lambda wildcards: maf[wildcards.maf],
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.i} --keep {input.keep} --filter-{wildcards.sex}s \
--maf {params.maf} --mac 10 --chr {wildcards.chrom} --geno 0.05 \
--hwe 0.000001 midp --hardy midp gz --make-bed --out {params.o}
"""

#look for significant nonrandom missingness
#make sure case-control is in file

rule test_miss:
    input: rules.filter_adgc.params.o + '.bim'
    output: "miss/{sample}.{sex}.maf-{maf}.chr{chrom}.missing",
            "miss/{sample}.{sex}.maf-{maf}.chr{chrom}.exclude",
            "miss/{sample}.{sex}.maf-{maf}.chr{chrom}.include"
    params:
        i = rules.filter_adgc.params.o,
        o = "miss/{sample}.{sex}.maf-{maf}.chr{chrom}",
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.i} --test-missing midp \
--out {params.o}

sed -r 's/[[:blank:]]+/ /g;s/^\s|\s$//g' {output[0]} | \
awk 'NR > 1 && $5 < 0.000001 {{print $2}}' > {output[1]}

awk 'NR == FNR {{a[$1]; next}} !($2 in a) {{print $2}}' \
<(cat <(echo exclude) {output[1]}) {input} > {output[2]}
"""

rule combine_miss:
    input:
        expand("miss/{{sample}}.{sex}.maf-{{maf}}.chr{{chrom}}.include",
               sex=['male', 'female'])
    output:
        "miss/allsex/{sample}.maf-{maf}.chr{chrom}.include"
    shell:
        """
awk 'NR == FNR {{a[$1]; next}} $1 in a {{print}}' {input} > {output}
"""

rule prep_GWAS:
    input:
        geno = rules.filter_adgc.output,
        keep = "miss/{sample}.{sex}.maf-{maf}.chr{chrom}.exclude",
    output: maybetemp(multiext("filtered/cov-{cov}.maf-{maf}/{sample}.{sex,male|female}.chr{chrom}", ".bed", ".bim", ".fam"))
    params:
        i = rules.filter_adgc.params.o,
        o = "filtered/cov-{cov}.maf-{maf}/{sample}.{sex}.chr{chrom}",
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.i} --exclude {input.keep} --make-bed --out {params.o}
"""

rule prep_GWAS_interact:
    input:
        geno = lambda wildcards: expand(infiles[wildcards.sample]["geno"] + ".{ext}", ext=BPLINK),
        ikeep = "{sample}.chosen.ikeep",
        keep = "miss/allsex/{sample}.maf-{maf}.chr{chrom}.include",
    output: maybetemp(multiext("filtered/cov-{cov}.maf-{maf}/{sample}.interaction.chr{chrom}", ".bed", ".bim", ".fam"))
    params:
        i = rules.filter_adgc.params.i,
        o = "filtered/cov-{cov}.maf-{maf}/{sample}.interaction.chr{chrom}",
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.i} --extract {input.keep} --keep {input.ikeep} --make-bed --out {params.o}
"""

rule do_GWAS:
    input:
        geno = rules.prep_GWAS.output,
        phen = "phenotypes/{sample}.pheno"
    output:
        "GWAS/cov-{{cov}}.maf-{{maf}}/{{sample}}.{{sex,male|female}}.chr{{chrom}}.assoc.{tt}".format(tt=config['test']),
        "GWAS/cov-{{cov}}.maf-{{maf}}/{{sample}}.{{sex,male|female}}.chr{{chrom}}.assoc.{tt}.perm".format(tt=config['test']) if config['perm'] else []
    params:
        i = rules.prep_GWAS.params.o,
        o = "GWAS/cov-{cov}.maf-{maf}/{sample}.{sex,male|female}.chr{chrom}",
        cov = lambda wildcards: study_covars[wildcards.study][wildcards.cov],
        perm = 'perm ' if config['perm'] else '',
        pname = lambda wildcards: config['studies'][wildcards.sample]['phenoname']
    shell:
        """
{loads[plink]}
{com[plink3]} --bfile {params.i} \
--pheno {input.phen} --pheno-name {params.pname} \
--covar {input.phen} --covar-name {params.cov} \
--{config[test]} {params.perm} genotypic beta --ci 0.99 --out {params.o}
"""


def gettests(wildcards):
    ncov = len(covariates[wildcards.cov].split(', '))
    return '1, 3-{}, {}-{}'.format(
        ncov + 2, 3 * ncov + 3, 3 * ncov + 4)


rule do_GWAS_interact:
    input:
        geno = rules.prep_GWAS_interact.output,
        phen = "phenotypes/{sample}.pheno"
    output:
        "GWAS/cov-{{cov}}.maf-{{maf}}/{{sample}}.interaction.chr{{chrom}}.assoc.{tt}".format(tt=config['test']),
        "GWAS/cov-{{cov}}.maf-{{maf}}/{{sample}}.interaction.chr{{chrom}}.assoc.{tt}.perm".format(tt=config['test']) if config['perm'] else []
    params:
        i = rules.prep_GWAS_interact.params.o,
        o = "GWAS/cov-{cov}.maf-{maf}/{sample}.interaction.chr{chrom}",
        cov = lambda wildcards: study_covars[wildcards.study][wildcards.cov],
        perm = 'perm ' if config['perm'] else '',
        pname = lambda wildcards: config['studies'][wildcards.sample]['phenoname'],
        tests = gettests
    shell:
        """
{loads[plink]}
{com[plink3]} --bfile {params.i} --parameters {params.tests} \
--pheno {input.phen} --pheno-name {params.pname} \
--covar {input.phen} --covar-name {params.cov} \
--{config[test]} {params.perm} genotypic sex interaction beta \
--ci 0.99 --out {params.o}
"""

#--parameters 1, 4, 6-7
rule fix_gwas:
    input: expand("GWAS/cov-{{cov}}.maf-{{maf}}/{{sample}}.{{allsex}}.chr{chrom}.assoc.{tt}", chrom=CHROM, tt=config['test'])
    output: "GWAS/cov-{{cov}}.maf-{{maf}}/{{sample}}.{{allsex}}_filtered.assoc.{tt}".format(tt=config['test'])
    shell:
        r"""
sed -r 's/[[:blank:]]+/ /g;s/^\s|\s$//g' {input} | \
awk 'NR == 1 || ($5 == "ADD" && $7 != "NA")' | \
awk 'BEGIN {{FS=" |:"}} NR == 1 {{print $0, "A2"}} NR != 1 {{print $0, $4}}' > \
{output}
"""

rule gwas_manhattan:
    input: rules.fix_gwas.output
    output: [rules.fix_gwas.output[0] + '.' + x for x in MANEXT]
    log: "GWAS/cov-{cov}.maf-{maf}/{sample}.{allsex}_filtered.plots.log"
    shell:
        """
{loads[R]}
scripts/manhattan.R {input} &> {log}
"""

rule add_emp:
    input:
        sstats = rules.fix_gwas.output,
        emp = expand("GWAS/cov-{{cov}}.maf-{{maf}}/{{sample}}.{{allsex}}.chr{chrom}.assoc.{tt}.perm", chrom=CHROM, tt=config['test'])
    output: "GWAS/cov-{{cov}}.maf-{{maf}}/{{sample}}.{{allsex}}_empP.assoc.{tt}".format(tt=config['test'])
    shell:
        r"""
awk 'NR == FNR {{emp[$2] = $3}} NR != FNR {{print $0, emp[$2]}}' \
<(cat {input.emp} | sed -r 's/[[:blank:]]+/ /g;s/^\s|\s$//g') \
{input.sstats} > {output}
"""

rule emp_manhattan:
    input: rules.add_emp.output
    output: [rules.add_emp.output[0] + '.' + x for x in MANEXT]
    log: "GWAS/cov-{cov}.maf-{maf}/{sample}.{allsex}_empP.plots.log"
    shell:
        """
{loads[R]}
scripts/manhattan.emp.R {input} &> {log}
"""
