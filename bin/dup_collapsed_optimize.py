#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
optimization of rescued collapsed contigs
"""

import argparse
import logging
import os
import os.path as op
import sys

import glob
import re
import numpy as np
import pandas as pd 

from collections import OrderedDict
from joblib import Memory
from pytools import natsorted

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


class clmLine():
    def __init__(self, line):
        self.line = line.strip()
        self.line_list = self.line.split("\t")
        self.ctg1, self.ctg2 = self.line_list[0].split()
        self.strand1 = self.ctg1[-1]
        self.strand2 = self.ctg2[-1]
        self.ctg1 = self.ctg1[:-1]
        self.ctg2 = self.ctg2[:-1]
        self.count = int(self.line_list[1])
        self.distances = list(map(int, self.line_list[2].split()))
    
    @property
    def dk(self):
        
        reciprocal_values = list(map(lambda x: 1/(x + 1), self.distances))
        return sum(reciprocal_values)
    
    def __str__(self):
        return self.line
 

class clm(object):
    def __init__(self, clmfile, mem_cache='.'):
        
        self.file = clmfile
        self.parse()
        # self.mem_cache = mem_cache
        # self.memory = Memory(mem_cache, verbose=0)
        # self.dk_df = self.memory.cache(self._dk_df)

    def parse(self):
        self.data = {}
        with open(self.file, 'r') as fp:
            for line in fp:
                cl = clmLine(line)
                if (cl.ctg1, cl.ctg2) not in self.data:
                    self.data[(cl.ctg1, cl.ctg2)] = []
                
                self.data[(cl.ctg1, cl.ctg2)].append(cl.dk)
    
    @property
    def dk_df(self):
        res_dk_df = pd.DataFrame(self.data)
        res_dk_df = res_dk_df.T 
        res_dk_df.columns = self.strands
        
        return res_dk_df

    @property
    def strands(self):
        return [('+', '+'), ('+', '-'),
                ('-', '+'), ('-', '-')]
      

def complementary(strand):
    s1, s2 = strand
    if s1 == '+':
        new_s1 = '-'
    else:
        new_s1 = '+'
    
    if s2 == '+':
        new_s2 = '-'
    else:
        new_s2 = '+'
    
    return (new_s1, new_s2)

def import_countRE(countRE_table):
    logger.info('Loadding `{}`'.format(countRE_table))
    df = pd.read_csv(countRE_table, header=0, index_col=0,
                    sep='\t')
    
    return df


def symmetric_pairs(pairs_df):
    """
    symmetric pairs dataframe
    """
    pairs_df2 = pairs_df.copy()
    pairs_df2['Contig1'], pairs_df2['Contig2'] = pairs_df['Contig2'], pairs_df['Contig1']
    pairs_df2['#X'], pairs_df2['Y'] = pairs_df['Y'], pairs_df['#X']
    #pairs_df2['Group1'], pairs_df2['Group2'] = pairs_df['Group2'], pairs_df['Group1']
    pairs_df2['RE1'], pairs_df2['RE2'] = pairs_df['RE2'], pairs_df['RE1']
    pairs_df = pd.concat([pairs_df, pairs_df2])

    return pairs_df 

def import_pairs(pairs_table):

    logger.info('Loading `{}`'.format(pairs_table))
    df = pd.read_csv(pairs_table, header=0, index_col=None,
                    sep='\t')
    df = df.astype(
        {'Contig1': 'category', 'Contig2': 'category',
        '#X': 'int32', 'Y': 'int32', 'RE1': 'int64', 
        'RE2': 'int64'}
    )

    return df

def tail(infile, n, offset=0):
    from subprocess import Popen, PIPE
    p = Popen(['tail', '-n', f'{n + offset}', infile], stdout=PIPE)
    lines = p.stdout.readlines()

    return lines

def rm_orientation(contig):
    return re.sub('[+-]$', '', contig)

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('rescue', 
            help='rescued collapsed contig table.')
    pReq.add_argument('pairs',
            help='pairs file from ALLHiC.')
    pReq.add_argument('clm', 
            help='clm file from ALLHiC.')
    pReq.add_argument('tour', 
            help='tour')
    pOpt.add_argument('-o', '--outdir', default='./',
             help='output directory [default: .]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    outdir = args.outdir
    if not op.exists(outdir):
        os.makedirs(outdir)
    else:
        logger.warning(f'Output directory `{outdir}` existed, '
                            'previous results will be covered.')

    if op.exists('clm.dk.df'):
        clm_dk_df = pd.read_csv('clm.dk.df', sep='\t', 
                    header=0, index_col=[0, 1])
        clm_dk_df.columns = [('+', '+'), ('+', '-'),
                ('-', '+'), ('-', '-')]
    else:
        clm_dk_df = clm(args.clm).dk_df
        clm_dk_df.to_csv('clm.dk.df', sep='\t', 
                            index=True, header=True)
    
    rescue_table = dict(i.strip().split() for i in open(args.rescue))

    pairs_df = import_pairs(args.pairs)
    pairs_df = symmetric_pairs(pairs_df)
    pairs_df = pairs_df[pairs_df['Label'] == 'ok']
    pairs_df = pairs_df.set_index(['Contig1', 'Contig2'])
    pairs_df['score'] = pairs_df['ObservedLinks'] / (pairs_df['RE1'] + pairs_df['RE2'])
    
    tour_file_list = natsorted(glob.glob(f"{args.tour}/*tour"))
    tour_name_list = natsorted(list(map(lambda x: op.basename(x).replace('.tour', ''), 
                            tour_file_list)))
    tour_file_db = dict(zip(tour_name_list, tour_file_list))
    tour_contents = list(map(lambda x: tail(x, 1), tour_file_list))
    tour_contents = list(map(lambda x: list(map(lambda x: x.decode(), x[0].split())), 
                                tour_contents))
    tour_contents_rm_orientation = list(map(lambda x: list(map(rm_orientation, x)), 
                                tour_contents))
    tour_only_orientation_list = list(map(lambda x: list(map(lambda x: x[-1], x)), tour_contents))
    tour_only_orientation_db = {}
    for i, name in enumerate(tour_name_list):
        tour_only_orientation_db[name] = dict(zip(tour_contents_rm_orientation[i], 
                                                tour_only_orientation_list[i]))

    tour_db = dict(zip(tour_name_list, tour_contents))
    tour_rm_orientation_db = dict(zip(tour_name_list, tour_contents_rm_orientation))

    
    for contig, group in rescue_table.items():
        chrom_tour_rm_orientation = tour_rm_orientation_db[group]
        contig_raw_name = re.sub(r'_d[0-9]+$', '', contig)
        try:
            contig_pairs = pairs_df.loc[contig_raw_name]
        except KeyError:
            # logger.warning(f'No contacts found for {contig}')
            continue
        contig_pairs = contig_pairs[contig_pairs.index.isin(chrom_tour_rm_orientation)]
        if contig_pairs.empty is True:
            # logger.warning(f'No contact of {contig} in {group}')
            continue
        contig_pairs = contig_pairs.sort_values(by='score', ascending=False)
        target_contig_pairs = contig_pairs.iloc[0]

        target_contig = target_contig_pairs.name 
        target_real_orientation = tour_only_orientation_db[group][target_contig]
        target_contig_idx = chrom_tour_rm_orientation.index(target_contig)
        
        flag = 0
        try: 
            tmp_dk = clm_dk_df.loc[contig_raw_name, target_contig]
            
        except KeyError:
            flag = 1
            try:
                tmp_dk = clm_dk_df.loc[target_contig, contig_raw_name]
            except KeyError:
                continue
        
        if tmp_dk.empty is True:
            continue

        tmp_dk = tmp_dk.sort_values(ascending=False)
        tmp_dk = tmp_dk.to_frame().reset_index()
        correct_orientation = tmp_dk.iloc[0]['index']
        correct_orientation = correct_orientation.tolist()[0]
        if flag == 1:
            source_orientation, target_orientation = correct_orientation
        
        else:
            source_orientation, target_orientation = complementary(correct_orientation)
            
        if target_orientation == '+':
            if target_real_orientation == '+':
                tour_db[group].insert(target_contig_idx+1, f'{contig}{source_orientation}')
                print(contig, target_contig, 'right', source_orientation)
            else:
                if source_orientation == '+':
                    source_orientation = '-'
                else: 
                    source_orientation = '+'
                tour_db[group].insert(target_contig_idx, f'{contig}{source_orientation}')
                print(contig, target_contig, 'left', source_orientation)
        else:
            if target_real_orientation == '+':
                tour_db[group].insert(target_contig_idx, f'{contig}{source_orientation}')
                print(contig, target_contig, 'left', source_orientation)
            else:
                if source_orientation == '+':
                    source_orientation = '-'
                else: 
                    source_orientation = '+'
                tour_db[group].insert(target_contig_idx+1, f'{contig}{source_orientation}')
                print(contig, target_contig, 'right', source_orientation)

    
    for group in tour_db:
        with open(f'{outdir}/{group}.tour', 'w') as out:
            print(' '.join(tour_db[group]), file=out)

        # if target_contig_idx == 0:
        #     pass
        # elif target_contig_idx == len(chrom_tour_rm_orientation) - 1:
        #     pass
        # else:
        #     up_idx = target_contig_idx - 10
        #     down_idx = target_contig_idx + 10
        #     up_contigs = chrom_tour_rm_orientation[up_idx: target_contig_idx]
        #     down_contigs = chrom_tour_rm_orientation[target_contig_idx + 1 : down_idx - 1]
            
        #     up_score = contig_pairs[contig_pairs.index.isin(up_contigs)]['score'].mean()
        #     up_score = 0 if up_score is np.nan else up_score

        #     down_score = contig_pairs[contig_pairs.index.isin(down_contigs)]['score'].mean()
        #     down_score = 0 if up_score is np.nan else down_score
            
        #     if up_score >= down_score:
        #         print(contig, target_contig, 'up')
        #     elif up_score < down_score:
        #         print(contig, target_contig, 'down')
            

if __name__ == "__main__":
    main(sys.argv[1:])
