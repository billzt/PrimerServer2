#!/usr/bin/env python3

'''PrimerServer2: a high-throughput primer design and specificity-checking platform
Github: https://github.com/billzt/PrimerServer2
This is a CLI script for PrimerServer2
'''

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
    parser = argparse.ArgumentParser(description='primertool: the command-line version of PrimerServer2', \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parent arguments
    # These arguments are used by all the three modes
    parent_parser_all = argparse.ArgumentParser(add_help=False)
    group_io = parent_parser_all.add_argument_group('Input (REQUIRED)')
    group_io.add_argument('query', help='query file. (STDIN is acceptable)', \
        type=argparse.FileType('r'))
    group_io.add_argument('templates', help='template file in FASTA format. \
        Allowing multiple files (separated by comma), where the first one is used to design primers and/or order the primer \
            specificity')

    group_all = parent_parser_all.add_argument_group('Overall Setttings')
    group_all.add_argument('--primer-num-retain', type=int, help='The maximum number of primers to retain in each \
        site in the final report.', default=10)
    group_all.add_argument('-p', '--cpu', type=int, help='Used CPU number.', default=2)
    group_all.add_argument('--json-debug', help="Output debug information in JSON mode", action='store_true')
    group_all.add_argument('-o', '--out', help="Output primers in JSON format. (Default is STDIN)", type=argparse.FileType('w'))
    group_all.add_argument('-t', '--tsv', help="Output primers in TSV format", type=argparse.FileType('w'))

    # These arguments are used by design and full
    parent_parser_design = argparse.ArgumentParser(add_help=False)
    group_design = parent_parser_design.add_argument_group('Design Primers')
    group_design.add_argument('--type', choices=['SEQUENCE_TARGET', 'SEQUENCE_INCLUDED_REGION', 'FORCE_END'],\
        help='designing primer types', default='SEQUENCE_TARGET')
    group_design.add_argument('--product-size-min', type=int, help='Lower limit of the product amplicon size range (bp).', \
        default=70)
    group_design.add_argument('--product-size-max', type=int, help='Upper limit of the product amplicon size range (bp).', \
        default=1000)
    group_design.add_argument('--primer-num-return', type=int, help='The maximum number of primers to return in Primer3 \
        designing results.', default=30)
    
    # These arguments are used by check and full
    parent_parser_check = argparse.ArgumentParser(add_help=False)
    group_check = parent_parser_check.add_argument_group('Check Specificity')
    group_check.add_argument('--Tm-diff', type=int, \
        help='The mininum melting temperature (in ℃) suggested to produce off-target amplicon. Suggest >10', default=20)
    group_check.add_argument('-3', '--use-3-end', help='If turned on, primer pairs having at least one mismatch at the 3 end\
        position with templates would not be considered to produce off-target amplicon, even if their melting temperatures \
            are high. Turn on this would find more candidate primers, but might also have more false positives\
                ', action='store_true')
    group_check.add_argument('--checking-size-min', type=int, help='Lower limit of the checking amplicon size range (bp).', \
        default=50)
    group_check.add_argument('--checking-size-max', type=int, help='Upper limit of the checking amplicon size range (bp).', \
        default=2000)
    group_check.add_argument('-a', '--report-amplicon-seqs', help="Get amplicon seqs (might be slow)", action='store_true')

    # sub commands
    subparsers = parser.add_subparsers(help='Sub commands: running mode')
    subparsers.required = True
    subparsers.dest = 'run_mode'
    # sub commands: full (design primers and check specificity)
    parser_full = subparsers.add_parser('full', help='design primers, check their specificity and output specific ones', \
        parents=[parent_parser_all, parent_parser_design, parent_parser_check], \
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # sub commands: design (design primers only, don't check specificity)
    parser_design = subparsers.add_parser('design', help="design primers only, don't check specificity", \
        parents=[parent_parser_all, parent_parser_design], formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # sub commands: check (check specificity only, skip the design step)
    parser_check = subparsers.add_parser('check', help='check specificity only, skip the design step', \
        parents=[parent_parser_all, parent_parser_check], formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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

def check_templates(args):
    for template in args.templates.split(','):
        if os.path.isfile(template) is False:
            error(f'File not found: {template}', args.json_debug)
        if os.path.isfile(template+'.fai') is False:
            code = os.system(f'samtools faidx {template} 2>/dev/null')
            if code != 0:
                error(f'File {template} cannot be indexed by samtools faidx. Perhaps it is not in FASTA format', args.json_debug)

def run(args):
    ###################  Design primers ###################
    query_string = args.query.read()
    dbs = args.templates.split(',')
    if args.run_mode=='check':
        primers = make_primers.make_primers(query=query_string)
        if 'error' in primers:
            error(primers['error'], args.json_debug)
    else:
        if args.run_mode=='design':
            sites = make_sites.build(query=query_string, template_file=dbs[0], primer_type=args.type, \
                primer_num_return=args.primer_num_retain, size_min=args.product_size_min, size_max=args.product_size_max)
        else:
            sites = make_sites.build(query=query_string, template_file=dbs[0], primer_type=args.type, \
                primer_num_return=args.primer_num_return, size_min=args.product_size_min, size_max=args.product_size_max)
            if 'error' in sites:
                error(sites['error'], args.json_debug)
        primers = design_primer.multiple(sites, cpu=args.cpu)

    ###################  Checking specificity  #############
    if args.run_mode!='design':
        primers = run_blast.run_blast_parallel(primers=primers, dbs=dbs, cpu=args.cpu,\
            checking_size_max=args.checking_size_max, checking_size_min=args.checking_size_min, \
                report_amplicon_seq=args.report_amplicon_seqs, Tm_diff=args.Tm_diff, use_3_end=args.use_3_end)
        primers = sort_primers.sort_rank(primers=primers, dbs=dbs, max_num_return=args.primer_num_retain)

    ###################  Output  ###########################
    if args.out is not None:
        print(json.dumps({'meta':{'mode':args.run_mode, 'dbs':dbs}, 'primers':primers}, indent=4), file=args.out)
    else:
        print(json.dumps({'meta':{'mode':args.run_mode, 'dbs':dbs}, 'primers':primers}, indent=4))

    if args.tsv is not None:
        print(output.tsv(primers, dbs), file=args.tsv)

def main():
    args = make_args()
    check_environments(args)
    check_templates(args)
    run(args)

if __name__ == "__main__":
    main()

