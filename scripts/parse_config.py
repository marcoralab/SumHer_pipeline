import re
import os
from collections import defaultdict
import socket


def parse_chrom(chrs):
    clist = [x.split(":") for x in chrs.split(",")]
    parsed = []
    for chrs in clist:
        if len(chrs) == 2:
            chrs = [str(c) for c in range(int(chrs[0]), int(chrs[1]) + 1)]
        elif len(chrs) != 1:
            raise ValueError("Invalid chromosome list.")
        parsed += chrs
    return parsed


def fix_path(path):
    path = os.path.abspath(path) + '/'
    return path


def parser(config):
    # Construct chromosome list
    chrom = parse_chrom(config["chroms"])

    return chrom


def subs(st):
    st = st.replace(' ', 'space')
    st = st.replace('/', 'slash')
    st = st.replace('(', 'openper')
    st = st.replace(')', 'closper')
    st = st.replace('"', 'quo')
    st = st.replace('&', 'andsymbol')
    st = st.replace('|', 'orsymbol')
    st = st.replace('>', '_gt_')
    st = st.replace('<', '_lt_')
    st = st.replace('!', 'excl')
    return(st)


def getcom(both, lsf, local={}):
    isMinerva = "hpc.mssm.edu" in socket.getfqdn()
    addons = defaultdict(lambda: '', lsf if isMinerva else local)
    return {k: ' '.join([v, addons[k]]) for k, v in both.items()}

def getloads(config, lsf, local=None, lsfextra=None):
    def buildlsf(packages):
        pk = packages if isinstance(packages, list) else [packages]
        return ' '.join([config['modules'][i] for i in pk])
    if "hpc.mssm.edu" in socket.getfqdn():
        ml = {k: 'module load {}'.format(buildlsf(v)) for k, v in lsf.items()}
        if isinstance(lsfextra, dict):
            addons = defaultdict(lambda: '', lsfextra)
            return {k: '; '.join([v, addons[k]]) for k, v in ml.items()}
        else:
            return ml
    else:
        rr = {k: 'echo running {}'.format(k) for k in lsf.keys()}
        if isinstance(local, dict):
            addons = defaultdict(lambda: '', local)
            return {k: '; '.join([v, addons[k]]) for k, v in rr.items()}
        else:
            return rr


def proc_cohorts(study, config):
    if "cohorts" not in study:
        return []
    if "cohorts_all" in study:
        ALLCOHORTS = set(study["cohorts_all"])
        COHORTS = set(study["cohorts"])
        if 'cohorts_fam' in study:
            FAMCOHORTS = study['cohorts_fam']
        else:
            FAMCOHORTS = []
        if not COHORTS <= ALLCOHORTS:
            print(COHORTS - ALLCOHORTS)
        assert COHORTS <= ALLCOHORTS, "Invalid cohort name in config."
        if study['exclude_cohorts']:
            if study['exclude_fams']:
                EXCLUDECOHORTS = COHORTS | FAMCOHORTS
            else:
                EXCLUDECOHORTS = COHORTS
            COHORTS = ALLCOHORTS - EXCLUDECOHORTS
        else:
            EXCLUDECOHORTS = ALLCOHORTS - COHORTS
    else:
        if study['exclude_cohorts']:
            raise NameError('you must define cohorts_all to exclude cohorts')
        COHORTS = study["cohorts"]
    if "default_cohort" in study:
        DEFAULT_COHORT = study["default_cohort"]
    else:
        DEFAULT_COHORT = COHORTS[0]
    assert DEFAULT_COHORT in COHORTS
    EXCLUDECOHORTS = list(EXCLUDECOHORTS)
    COHORTS.remove(DEFAULT_COHORT)
    DUMMY_COHORTS = ['cohort_' + x for x in COHORTS]
    COHORTS = [DEFAULT_COHORT] + list(COHORTS)
    return {'COHORTS': COHORTS,
            'DUMMY_COHORTS': ', '.join(DUMMY_COHORTS),
            'DEFAULT_COHORT': 'cohort_' + DEFAULT_COHORT,
            'EXCLUDECOHORTS': EXCLUDECOHORTS}


def proc_covars(studyname, config, study_cohorts):
    sc = config['studies'][studyname]
    covs = sc['covariates']
    dc = study_cohorts[studyname]['DUMMY_COHORTS']
    pclist = [sc['PC_format'].format(x) for x in range(1, config['PC_N'] + 1)]
    PCS = ", ".join(pclist)
    return {k: v.format(pcs=PCS, cht=dc) for k, v in covs.items()}
