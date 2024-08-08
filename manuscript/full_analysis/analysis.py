"""
runs the kanpig manuscript analysis scripts

writes output dict to a joblib file 

Structure:
    Level1 Keys:
        - 'base_gt_dist': dict of truvari.GT and counts from the truth vcf
        - 'dist_gt_dist': dict of truvari.GT and counts from the discovery vcf
        - 's8': section 8 results. first key is program name
        - 's5': section 5 results. first key is program name
    s5 keys: (after you access by program)
        - 'table': dataframe of overall stats
        - 'gt_dist': dict of truvari.GT and counts from the genotyped vcf
        - 'type': table of stats of only the discovery vcf by svtype
        - 'neigh': table of stats of only the discovery vcf by number of neighbors
    s8 keys: (after you access by program)
        - 'type': table of stats by svtype
        - 'neigh': table of stats by number of neighbors
        - 'gt_dist': dict of truvari.GT and counts from the genotyped vcf
"""
import os
import ast
import json
import argparse
import logging
import itertools

import pysam
import joblib
import shutil
import truvari

import numpy as np
import pandas as pd

from io import StringIO
from pysam import bcftools
from collections import Counter, defaultdict

SIZEMIN = 50
SIZEMAX = 10000
NEIGHLIMIT = 11
HAS_SVJEDI = ["HG002", "HG00621", "HG00673", "HG00735", "HG01928",
              "HG02572", "HG03098", "HG03492", "HG03516"]


def make_tables():
    gtmatrix = {}
    for key in truvari.GT:
        gtmatrix[key.name] = 0

    m_cnt_base = {'Concordant': 0,
                  'Discordant': 0,
                  'Missing': 0,
                  'Filtered': 0,
                  'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}

    # by neighbors
    m_neigh_cnt = {}
    for i in range(0, NEIGHLIMIT):
        m_neigh_cnt[i] = dict(m_cnt_base)

    # by svtype
    m_cnt = {'DEL': dict(m_cnt_base),
             'INS': dict(m_cnt_base),
             'TOT': dict(m_cnt_base)}
    return gtmatrix, m_cnt, m_neigh_cnt


def truth_single_compare(truth_fn, bed_fn, vcf_fn):
    """
    Merges a vcf file to the baseline
    Returns a pandas Series with Program Technology, and all the metrics we collect
    """
    gtmatrix, m_cnt, m_neigh_cnt = make_tables()

    # merge... if I can abstract this out, I can reuse the rest for grabbing stats
    out_fn1 = f"temp/merge.vcf.gz"
    bcftools.merge("-m", "none", "--force-samples", "-O", "z", "-o", out_fn1,
                   truth_fn, vcf_fn, catch_stdout=False)
    pysam.tabix_index(out_fn1, force=True, preset="vcf")

    vcf = pysam.VariantFile(out_fn1)
    bed = truvari.build_region_tree(vcfA=vcf, includebed=bed_fn)
    vcf_i = truvari.region_filter(vcf, bed)

    for entry in vcf_i:
        sz = truvari.entry_size(entry)
        if sz < SIZEMIN or sz > SIZEMAX:
            continue
        gt_base = truvari.get_gt(entry.samples[0]["GT"])
        gt_comp = truvari.get_gt(entry.samples[1]["GT"])
        gtmatrix[gt_comp.name] += 1
        svtype = truvari.entry_variant_type(entry).name

        if None in entry.samples[0]["GT"]:
            m_cnt[svtype]["Filtered"] += 1
            m_cnt['TOT']["Filtered"] += 1
            continue
        if None in entry.samples[1]["GT"]:
            m_cnt[svtype]['Missing'] += 1
            m_cnt['TOT']['Missing'] += 1
            continue
        c = sorted(list(entry.samples[0]['GT']))
        b = sorted(list(entry.samples[1]['GT']))
        state = "Concordant" if b == c else "Discordant"
        m_cnt[svtype][state] += 1
        m_cnt['TOT'][state] += 1

        numneigh = min(entry.info['NumNeighbors'], NEIGHLIMIT-1)
        m_neigh_cnt[numneigh][state] += 1

        for i, j in zip(b, c):
            state = "T" if i == j else "F"
            condition = "P" if j == 1 else "N"
            m_cnt[svtype][state+condition] += 1
            m_cnt['TOT'][state+condition] += 1
            m_neigh_cnt[numneigh][state+condition] += 1

    ret = []
    for k, v in m_cnt.items():
        v['svtype'] = k
        ret.append(v)
    bytype = pd.DataFrame(ret)

    ret = []
    for k, v in m_neigh_cnt.items():
        v['num_neigh'] = k
        ret.append(v)
    byneigh = pd.DataFrame(ret)

    return {'type': bytype, 'neigh': byneigh, 'gt_dist': gtmatrix}


def parse_summary_json(path):
    with open(path + '/summary.json', 'r') as fh:
        d = json.load(fh)

    d['gt_matrix']['(0, 0)'] = Counter()
    d['Concordant'] = d['TP-comp_TP-gt']
    d['Discordant'] = d['TP-comp_FP-gt']
    d['Missing'] = 0
    d['Present FP'] = 0
    d['Present TP'] = 0
    v = pysam.VariantFile(path + "/fp.vcf.gz")

    for entry in v:
        cgt = entry.samples[0]['GT']
        d['gt_matrix']['(0, 0)'][str(cgt)] += 1
        if None in cgt:
            d['Missing'] += 1
        elif 1 in cgt:  # Only penalize for non-refhom
            d['Discordant'] += 1
            d['Present FP'] += 1

    # reclassify discordant which should be considered missing
    v = pysam.VariantFile(path + '/tp-comp.vcf.gz')
    for entry in v:
        cgt = entry.samples[0]['GT']
        if None in cgt:
            d['Missing'] += 1
            d['Discordant'] -= 1
        elif 1 in cgt:
            d['Present TP'] += 1

    d['gt_accuracy'] = d['gt_concordance']
    d['gt_concordance'] = d['Concordant'] / (d['Concordant'] + d['Discordant'])

    astate = []
    d['as_TP'] = 0
    d['as_TN'] = 0
    d['as_FP'] = 0
    d['as_FN'] = 0
    d['as_filtered'] = 0
    d['as_missing'] = 0
    for b_key, b_data in d['gt_matrix'].items():
        b = list(ast.literal_eval(b_key))
        b = sorted([_ if _ is not None else -1 for _ in b])
        for c_key, c_data in b_data.items():
            c = list(ast.literal_eval(c_key))
            c = sorted([_ if _ is not None else -1 for _ in c])
            for i, j in zip(b, c):
                if i == -1:
                    d['as_filtered'] += c_data
                elif j == -1:
                    d['as_missing'] += c_data
                else:
                    state = "T" if i == j else "F"
                    condition = "P" if j == 1 else "N"
                    d[f'as_{state + condition}'] += c_data
    # Cleaning
    del (d['gt_matrix'])
    dkeys = [_ for _ in d.keys() if _.endswith('-gt')]
    for i in dkeys:
        del (d[i])
    return pd.DataFrame([d])


def run_truvari_single(base_vcf, base_bed, comp_vcf):
    # need a temp dir to put things in
    logging.info("Running truvari bench")
    bench_dir = truvari.make_temp_filename()
    r = truvari.cmd_exe(f"""
        truvari bench -b {base_vcf} --includebed {base_bed} \
            --sizemin {SIZEMIN} --sizefilt {SIZEMIN} --sizemax {SIZEMAX} -c {comp_vcf} --short -o {bench_dir}
    """)
    return bench_dir


def run_truvari_multi(base_vcf, base_bed, comp_vcf):
    # need a temp dir to put things in
    logging.info("Running truvari bench")
    bench_dir = truvari.make_temp_filename()
    r = truvari.cmd_exe(f"""
        truvari bench -b {base_vcf} --includebed {base_bed} \
            --no-ref a --pick ac --sizemin {SIZEMIN} --sizefilt {SIZEMIN} \
            --sizemax {SIZEMAX} -c {comp_vcf} --short -o {bench_dir}
    """)
    return bench_dir


def pull_gt_dist(vcf_fn, bed_fn, sample=0):
    """
    """
    logging.info("Pulling %s gt counts", vcf_fn)
    matrix = {}
    for key in truvari.GT:
        matrix[key.name] = 0

    vcf = pysam.VariantFile(vcf_fn)
    bed = truvari.build_region_tree(vcfA=vcf, includebed=bed_fn)
    vcf_i = truvari.region_filter(vcf, bed)
    for entry in vcf_i:
        gt = truvari.get_gt(entry.samples[sample]['GT'])
        matrix[gt.name] += 1

    return matrix


def dir_count(vcf_fn, is_tp, m_cnt, m_neigh_cnt, sample=0):
    """
    Updates in-place
    """
    vcf = pysam.VariantFile(vcf_fn)
    for entry in vcf:
        if not is_tp:
            state = "Discordant"
        elif entry.info['GTMatch'] == 0:
            state = 'Concordant'
        else:
            state = 'Discordant'

        svtype = truvari.entry_variant_type(entry).name
        m_cnt[svtype][state] += 1
        m_cnt['TOT'][state] += 1
        numneigh = min(entry.info['NumNeighbors'], NEIGHLIMIT-1)
        m_neigh_cnt[numneigh][state] += 1

        # I question this
        gt = truvari.get_gt(entry.samples[sample]['GT'])

        # We can only do allele state on matched variants.
        # actually, all the FPs can be parsed...
        if gt == truvari.GT.HET:
            b = (0, 1)
            if state == 'Discordant':
                c = (1, 1)
            else:
                c = (0, 1)
        elif gt == truvari.GT.HOM:
            b = (1, 1)
            if state == 'Discordant':
                c = (0, 1)
            else:
                b = (1, 1)

        if not is_tp:
            c = (0, 0)
            continue

        for i, j in zip(b, c):
            state = "T" if i == j else "F"
            condition = "P" if j == 1 else "N"
            m_cnt[svtype][state+condition] += 1
            m_cnt['TOT'][state+condition] += 1
            m_neigh_cnt[numneigh][state+condition] += 1


def parse_bench_dir(bdir, sample=0):
    """
    Make the bytype and byneigh dataframes
    This is only considering the baseline variants
    """
    _, m_cnt, m_neigh_cnt = make_tables()
    dir_count(os.path.join(bdir, 'tp-base.vcf.gz'),
              True, m_cnt, m_neigh_cnt, sample)
    dir_count(os.path.join(bdir, 'fn.vcf.gz'),
              False, m_cnt, m_neigh_cnt, sample)

    ret = []
    for k, v in m_cnt.items():
        v['svtype'] = k
        ret.append(v)
    bytype = pd.DataFrame(ret)

    ret = []
    for k, v in m_neigh_cnt.items():
        v['num_neigh'] = k
        ret.append(v)
    byneigh = pd.DataFrame(ret)

    return bytype, byneigh


def bench_single_compare(base_fn, bed_fn, vcf_fn):
    """
    Run truvari bench and then calculate the statistics
    """
    gtmatrix = pull_gt_dist(vcf_fn, bed_fn)
    bdir = run_truvari_single(base_fn, bed_fn, vcf_fn)
    ty, neigh = parse_bench_dir(bdir)
    return {'table': parse_summary_json(bdir), 'type': ty, 'neigh': neigh, 'gt_dist': gtmatrix}


def bench_multi_compare(base_fn, bed_fn, vcf_fn, sample):
    """
    Run truvari bench and then calculate the statistics
    """
    gtmatrix = pull_gt_dist(vcf_fn, bed_fn)
    bdir = run_truvari_multi(base_fn, bed_fn, vcf_fn)
    ty, neigh = parse_bench_dir(bdir, sample)
    return {'table': parse_summary_json(bdir), 'type': ty, 'neigh': neigh, 'gt_dist': gtmatrix}


def add_summary_08(data):
    """
    Add summary stats to a dataframe
    I might not be using this anymore?
    """
    data['Total Calls'] = data[['Concordant',
                                'Discordant', 'Missing', 'Filtered']].sum(axis=1)
    data['Total Genotyped'] = data['Total Calls'] - data['Missing']
    data['Missing Rate'] = data['Missing'] / data['Total Genotyped']
    data['GT Concordance'] = data['Concordant'] / data['Total Genotyped']
    data['compP'] = data['TP'] + data['FP']
    data['baseP'] = data['TP'] + data['FN']
    data['compN'] = data['TN'] + data['FN']
    data['baseN'] = data['TN'] + data['FP']
    data['ppv'] = data['TP'] / data['compP']
    data['tpr'] = data['TP'] / data['baseP']
    data['tnr'] = data['TN'] / data['baseN']
    data['npv'] = data['TN'] / data['compN']
    data['acc'] = (data['TP'] + data['TN']) / (data['baseP'] + data['baseN'])
    data['ba'] = (data['tpr'] + data['tnr']) / 2
    data['f1'] = 2 * ((data['ppv'] * data['tpr']) /
                      (data['ppv'] + data['tpr']))


def prepare_truth(truth, bed):
    """
    no chrX / chrY
    add num neigh
    returns new truth/bed file names
    """
    logging.info("Preparing truth")
    with open(bed, 'r') as fh, open('temp/new.bed', 'w') as fout:
        for line in fh:
            d = line.strip().split('\t')
            if not ('_' in d[0] or d[0] in ['chrX', 'chrY', 'chrM', 'chrEBV']):
                fout.write(line)

    # Do I need to remove large SVs before counting neighbors?
    anno_neigh(truth, "temp/baseline.neigh.vcf")
    return "temp/baseline.neigh.vcf.gz", "temp/new.bed"

def anno_neigh(in_vcf, out_fn):
    """
    annotate numneigh
    """
    r = truvari.cmd_exe(
        f"truvari anno numneigh -o {out_fn} {in_vcf}")
    truvari.compress_index_vcf(out_fn)


def clean_svjedi_8(orig_vcf):
    """
    The ad format is messed up in their output
    """
    logging.info("Cleaning svjedi 8")
    out_vcf = "temp/svjedi08.vcf.gz"
    bcftools.annotate("-x", "FORMAT/AD", "-O", "z", "-o",
                      out_vcf, orig_vcf, catch_stdout=False)
    pysam.tabix_index(out_vcf, force=True, preset="vcf")
    return "temp/svjedi08.vcf.gz"


if __name__ == '__main__':
    programs = ['kanpig', 'svjedi', 'sniffles', 'cutesv']
    parser = argparse.ArgumentParser(prog="analysis", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=str, required=True,
                        help="Output joblib file")
    parser.add_argument("--sample", type=str, required=True,
                        help="Sample name (for auto skipping some svjedi")
    parser.add_argument("--truth", type=str, required=True,
                        help="Dipcall truth VCF")
    parser.add_argument("--bed", type=str, required=True,
                        help="Dipcall high confidence bed")
    parser.add_argument("--discovery", metavar="DISC", type=str, required=True,
                        help="Single sample discovery vcf")
    for i in programs:
        # 7, 9... which will just use pretty much the same thing as 5...
        for j in ['5', '8']:
            # Which means we need to parse the tp-base fn to get neigh and maybe type...
            # HERE
            parser.add_argument(f"--{i}-{j}", metavar=f"{i[0].upper()}{j}", type=str, required=True,
                                help=f"{i} section {j} vcf")
    args = parser.parse_args()
    truvari.setup_logging(show_version=False)

    logging.info("Running analysis")
    if not os.path.exists('temp'):
        os.mkdir('temp')
    # Should I clean temp? nah.

    args.truth, args.bed = prepare_truth(args.truth, args.bed)

    anno_neigh(args.discovery, "temp/discovery.neigh.vcf")
    args.discovery = "temp/discovery.neigh.vcf.gz"

    base_gt_dist = pull_gt_dist(args.truth, args.bed)

    # Easier to index in a loop
    d_args = dict(args._get_kwargs())

    # Can't process some svjedi
    if args.sample not in HAS_SVJEDI:
        programs.remove('svjedi')
    else:
        d_args['svjedi_8'] = clean_svjedi_8(d_args['svjedi_8'])

    disc_gt_dist = pull_gt_dist(args.discovery, args.bed)

    t1_parts = {}
    d1_parts = {}
    for p in programs:
        logging.info("Processing %s", p)
        t1_parts[p] = truth_single_compare(
            args.truth, args.bed, d_args[f'{p}_8'])
        d1_parts[p] = bench_single_compare(
            args.discovery, args.bed, d_args[f'{p}_5'])
        # tm_parts p_9
        # dm_parts p_7

    out = {'base_gt_dist': base_gt_dist,
           'dist_gt_dist': disc_gt_dist,
           's8': t1_parts,
           's5': d1_parts}
    joblib.dump(out, args.output)
