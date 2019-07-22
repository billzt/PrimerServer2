#!/usr/bin/env python3

'''PrimerServer2: a high-throughput primer design and specificity-checking platform
Github: https://github.com/billzt/PrimerServer2
'''

__author__ = 'Tao Zhu'
__copyright__ = 'Copyright 2019'
__license__ = 'GPL'
__version__ = '0.1'
__email__ = 'taozhu@mail.bnu.edu.cn'
__status__ = 'Development'

import argparse
import re
import os
import sys
import json
import shutil

from distutils.version import LooseVersion

from primerserver2.core import make_sites, make_primers, design_primer, run_blast, sort_primers, output


def error(msg, judge):
    if judge is True:
        print(json.dumps({'ERROR': msg}))
    raise Exception(msg)

def make_args():
    parser = argparse.ArgumentParser(description='PrimerServer2: a high-throughput primer \
        design and specificity-checking platform', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('query', help='query file. (STDIN is acceptable)', type=argparse.FileType('r'))
    parser.add_argument('template', help='template file in FASTA format')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--no-specificity', help="Don't check specificity; Only design primers", action='store_true')
    group.add_argument('--only-specificity', help="Only check specificity and skip primer designs", action='store_true')
    parser.add_argument('--type', choices=['SEQUENCE_TARGET', 'SEQUENCE_INCLUDED_REGION', 'FORCE_END'],\
        help='designing primer types', default='SEQUENCE_TARGET')
    parser.add_argument('-p', '--cpu', type=int, help='Used CPU number.', default=2)
    parser.add_argument('-d', '--Tm-diff', type=int, \
        help='The mininum melting temperature (in ℃) suggested to produce off-target amplicon. Suggest >10', default=20)
    parser.add_argument('-3', '--use-3-end', help='If turned on, primer pairs having at least one mismatch at the 3 end\
        position with templates would not be considered to produce off-target amplicon, even if their melting temperatures \
            are high. Turn on this would find more candidate primers, but might also have more false positives\
                ', action='store_true')
    parser.add_argument('-l', '--checking-size-min', type=int, help='Lower limit of the checking amplicon size range (bp).', \
        default=70)
    parser.add_argument('-u', '--checking-size-max', type=int, help='Upper limit of the checking amplicon size range (bp).', \
        default=1000)
    parser.add_argument('-r', '--primer-num-retain', type=int, help='The maximum number of primers for each site to return.', \
        default=10)
    parser.add_argument('-a', '--report-amplicon-seqs', help="Get amplicon seqs (might be slow)", action='store_true')
    parser.add_argument('--json-debug', help="Output debug information in JSON mode", action='store_true')
    parser.add_argument('-o', '--out', help="Output primers in JSON format. (Default is STDIN)", type=argparse.FileType('w'))
    parser.add_argument('-t', '--tsv', help="Output primers in TSV format", type=argparse.FileType('w'))
    args = parser.parse_args()
    return args

def check_environments(args):
    if shutil.which('samtools') is None:
        error('No samtools detected in your system', args.json_debug)

    samtools_version = os.popen('samtools --version').readlines()[0].strip().split(' ')[1]
    if LooseVersion(samtools_version) < LooseVersion('1.9'):
        error(f'Your samtools version is v{samtools_version}, but >=v1.9 is required', args.json_debug)

    if shutil.which('blastn') is None:
        error('No NCBI-BLAST+ (blastn) detected in your system', args.json_debug)

    if shutil.which('makeblastdb') is None:
        error('No NCBI-BLAST+ (makeblastdb) detected in your system', args.json_debug)

def check_template(args):
    if os.path.isfile(args.template) is False:
        error(f'File not found: {args.template}', args.json_debug)
    if os.path.isfile(args.template+'.fai') is False:
        code = os.system(f'samtools faidx {args.template} 2>/dev/null')
        if code != 0:
            error(f'File {args.template} cannot be indexed by samtools faidx. Perhaps it is not in FASTA format', args.json_debug)

def run(args):
    ###################  Design primers ###################
    query_string = args.query.read()
    if args.only_specificity is True:
        primers = make_primers.make_primers(query=query_string)
    else:
        sites = make_sites.build(query=query_string, template_file=args.template, primer_type=args.type)
        primers = design_primer.multiple(sites, cpu=args.cpu)

    ###################  Checking specificity  #############
    dbs = [args.template]
    if args.no_specificity is False:
        primers = run_blast.run_blast_parallel(primers=primers, dbs=dbs, cpu=args.cpu,\
            checking_size_max=args.checking_size_max, checking_size_min=args.checking_size_min, \
                report_amplicon_seq=args.report_amplicon_seqs, Tm_diff=args.Tm_diff, use_3_end=args.use_3_end)
        primers = sort_primers.sort_rank(primers=primers, dbs=dbs, max_num_return=args.primer_num_retain)

    ###################  Output  ###########################
    if args.out is not None:
        print(json.dumps(primers, indent=4), file=args.out)
    else:
        print(json.dumps(primers, indent=4))

    if args.tsv is not None:
        if args.no_specificity is False:
            output.after_check(primers, dbs, file=args.tsv)

def main():
    args = make_args()
    check_environments(args)
    check_template(args)
    run(args)

if __name__ == "__main__":
    main()

