#!/usr/bin/env python

import argparse
import os
import os.path as op
import sys
import shutil
import pandas as pd
import re

from pyfaidx import Fasta

def import_table(table):
    df = pd.read_csv(table, sep=' ', header=None, index_col=None)

    return df

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('table', 
            help='table from `dup_collapsed_optimize.py`')
    pReq.add_argument('fasta', 
            help='fasta of draft assembly')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    table = args.table
    fasta = args.fasta
    
    df = import_table(table)
    fasta = Fasta(fasta)
    with open(args.fasta, 'r') as fp:
        shutil.copyfileobj(fp, args.output)
    
    for i, contig in df[0].iteritems():
        contig_raw_name = re.sub(r'_d[0-9]+$', '', contig)
        seq = fasta[contig_raw_name]
        print(f'>{contig}\n{seq}', file=args.output)


if __name__ == "__main__":
    main(sys.argv[1:])
