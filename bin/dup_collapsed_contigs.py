#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
duplicate collapsed contigs by copy number
"""

import argparse
import logging
import os
import os.path as op
import sys
import re

import numpy as np
import pandas as pd

from collections import OrderedDict
from copy import deepcopy
from rich.logging import Console, RichHandler

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=Console(stderr=True))]
)
logger = logging.getLogger(__name__)
logging.getLogger('numexpr').setLevel(logging.ERROR)

class clusterGroup:
    def __init__(self, line):
        line_list = line.strip().split("\t")
        self.line = line.strip()
        self.group = line_list[0]
        self.nContigs = line_list[1]
        self.Contigs = line_list[2].split()
    
    def __repr__(self):
        return self.line 

    __str__ = __repr__
    
class clusterTable:
    def __init__(self, filename):
        self.filename = filename
        logger.info(f'Loading `{self.filename}` ...')
        self.get_data()

    def parse(self):
        with open(self.filename, 'r') as fp:
            for line in fp:
                if line[0] == "#":
                    continue
                yield clusterGroup(line)

    def get_data(self):
        self.data = OrderedDict()
        for i in self.parse():
            self.data[i.group] = i.Contigs

    @property
    def contig_groups(self):
        _res = OrderedDict()
        for group in self.data:
            for contig in self.data[group]:
                _res[contig] = group
        
        return _res
        
    @property
    def groups(self):
        _group = [i.group for i in self.parse()]
        return _group

    @property
    def ncontigs(self):
        _ncontigs = OrderedDict()
        for i in self.parse():
            _ncontigs[i.group] = int(i.nContigs)
        return _ncontigs

    @property
    def contigs(self):
        _contigs = []
        for i in self.parse():
            _contigs.extend(i.Contigs)

        return _contigs
    
    def write(self, output):
        for group in self.data:
            print(group, len(self.data[group]), 
                    " ".join(self.data[group]), 
                    sep='\t', file=output)

    def get_max_group(self):
        res = 0
        for i in self.ncontigs:
            if res < self.ncontigs[i]:
                res = self.ncontigs[i]
        
        return res 

def import_allele_table(table, ploidy):
    logger.info(f'Loading `{table}` ...')
    names = ['chrom', 'site'] + list(range(0, ploidy))
    df = pd.read_csv(table, sep='\t', header=None, names=names)
    df = df.drop_duplicates(subset=range(ploidy))

    return df

def import_cn(infile):
    logger.info(f'Loading `{infile}` ...')
    df = pd.read_csv(infile, sep='\t', header=0, index_col=0)
    df.columns = [1]
    df[df[1] == '-'] = np.nan
    df = df[1].astype('float32')
    
    return df

def import_cluster(cluster_table):
    cluster_data = clusterTable(cluster_table)
    max_group_length = cluster_data.get_max_group()
    
    df = pd.DataFrame({key: pd.Series(value) for key, value 
                            in cluster_data.data.items()})
    
    df = df.T

    return df 
    
def find_contig(df, contig):
    cond = df == contig
    res = df.where(cond)
    res = res.dropna(how='all')

    return res.index

def find_cluster_group(x, cluster_data):
    if x is np.nan:
        return np.nan
    try:
        res = cluster_data.contig_groups[x]
        return res
    except:
        return np.nan

def get_chrom_group(chrom):
    regex = re.compile(r'(Chr\d+)\w+.*')
    group = regex.match(chrom).groups()[0]
    
    return group

def get_groups(group, chrom_list):
    regex = re.compile(r"{}\w+.*".format(group))
    groups = list(filter(regex.match, chrom_list))
    
    return groups

def split_dataframe_by_group(df):
    data = OrderedDict()
    for chrom, item in df.groupby('chrom'):
        item = item.drop(['chrom', 'site'], axis=1)
        data[chrom] = item

    return data

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('allele', 
            help='allele table from ALLHiC')
    pReq.add_argument('cn', 
            help='`06.genes.round.cn` from popCNV')
    pReq.add_argument('cluster',
            help='cluster file from ALLHiC')
    pReq.add_argument('ploidy',  type=int,
            help='ploidy')
    pOpt.add_argument('-o', '--outdir', default='./',
             help='output directory [default: .]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    allele_table = args.allele
    cn_table = args.cn
    cluster_table = args.cluster
    ploidy = args.ploidy
    outdir = args.outdir
    if not op.exists(outdir):
        os.makedirs(outdir)
    else:
        logger.warning(f'Output directory `{outdir}` existed, previous results will be covered.')

    ## out file
    out_cluster = f'{outdir}/collapsed.rescued.cluster.txt'
    out_rescued = f'{outdir}/collapsed.rescued.txt'
    out_unrescued = f'{outdir}/collapsed.unrescued.txt'

    cn_df = import_cn(cn_table)
    allele_df = import_allele_table(allele_table, ploidy)
    allele_split = split_dataframe_by_group(allele_df)
    cluster_data = clusterTable(cluster_table)
    ## filter collapsed contigs
    collapsed_cn_df = cn_df[cn_df >= 2]
    collapsed_cn_df = collapsed_cn_df.astype('int8')
    logger.info(f'Get {len(collapsed_cn_df)} collapsed contigs.')
    
    new_cluster_data = deepcopy(cluster_data)
    rescued_res = []
    unrescued_res = []

    logger.info(f'Starting to rescue collapsed contigs ...')
    for contig, value in collapsed_cn_df.items():
        try:
            real_chrom = cluster_data.contig_groups[contig]
        except KeyError:
            for i, j in enumerate(range(2, value+1)):
                new_contig_name = f'{contig}_d{j}'
                unrescued_res.append((new_contig_name, 'unanchor'))
            continue
        real_group = get_chrom_group(real_chrom)
        groups = get_groups(real_group, cluster_data.groups)
        groups_count = pd.Series(dict(zip(groups, [0]*len(groups))))
        allele_chrom_df = allele_split[real_group]
        tmp_df = allele_chrom_df.loc[find_contig(allele_chrom_df, contig)]
        
        if tmp_df.empty is True: 
            for i, j in enumerate(range(2, value+1)):
                new_contig_name = f'{contig}_d{j}'
                unrescued_res.append((new_contig_name, 'unallelic'))
            continue

        tmp_df = tmp_df.loc[tmp_df.apply(lambda x: x.dropna().count(), 
                    axis=1).sort_values(ascending=False).index]

        tmp_df.reset_index(inplace=True, drop=True)
        tmp_group_df = tmp_df.applymap(lambda x: find_cluster_group(x, cluster_data))
        groups_count += tmp_group_df.apply(pd.value_counts).count(axis=1)
        groups_count = groups_count.fillna(0).sort_values(ascending=True)
        if (groups_count == len(tmp_df)).all():
            for i, j in enumerate(range(2, value+1)):
                new_contig_name = f'{contig}_d{j}'
                unrescued_res.append((new_contig_name, 'redundancy'))
            continue
        for i, j in enumerate(range(2, value+1)):
            new_contig_name = f'{contig}_d{j}'
            new_group = groups_count.index[i]
            if groups_count.loc[new_group] == len(tmp_df):
                unrescued_res.append((new_contig_name, 'redundancy'))
                continue
            rescued_res.append((new_contig_name, new_group))
            new_cluster_data.data[new_group].append(new_contig_name)
    
    
    with open(out_cluster, 'w') as out:
        new_cluster_data.write(out)

    with open(out_rescued, 'w') as out:
        for contig, chrom in rescued_res:
            out.write(f'{contig}\t{chrom}\n')
    
    with open(out_unrescued, 'w') as out:
        for contig, res in unrescued_res:
            out.write(f'{contig}\t{res}\n')

    logger.info(f'Result: [rescued collapsed contigs: {len(rescued_res)}]')
    logger.info(f'\t[unrescued collapsed contigs: {len(unrescued_res)}]')
    logger.info('Done.')
    
if __name__ == "__main__":
    main(sys.argv[1:])