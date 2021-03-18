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

from primerserver2.core import make_sites, make_primers, design_primer, run_blast, sort_primers, output, analysis_blast, global_var
from primerserver2.core import multiplex
from primerserver2.core import version

def make_args():
    parser = argparse.ArgumentParser(description='primertool: the command-line version of PrimerServer2', \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s '+version.get())

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
    group_all.add_argument('--check-multiplex', help='Checking dimers between primers in different sites, which is useful in \
        multiplex PCR.', action='store_true')
    group_all.add_argument('--Tm-diff', type=int, help='The mininum difference of melting temperature (℃) \
        suggested to produce off-target amplicon or primer dimers. Suggest >10', default=20)
    group_all.add_argument('-p', '--cpu', type=int, help='Used CPU number.', default=2)
    group_all.add_argument('-o', '--out', help="Output primers in JSON format. default: {query}.out", \
        type=argparse.FileType('w'))
    group_all.add_argument('-t', '--tsv', help="Output primers in TSV format. default: {query}.tsv", \
        type=argparse.FileType('w'))

    # These arguments are used by design and full
    parent_parser_design = argparse.ArgumentParser(add_help=False)
    group_design = parent_parser_design.add_argument_group('Design Primers')
    group_design.add_argument('--type', choices=['SEQUENCE_TARGET', 'SEQUENCE_INCLUDED_REGION', 'FORCE_END'],\
        help='designing primer types', default='SEQUENCE_TARGET')
    group_design.add_argument('--pick-oligo', action='store_true', help='Pick internal Oligos (Probes) for qRT-PCR')
    group_design.add_argument('--product-size-min', type=int, help='Lower limit of the product amplicon size range (bp).', \
        default=70)
    group_design.add_argument('--product-size-max', type=int, help='Upper limit of the product amplicon size range (bp).', \
        default=1000)
    group_design.add_argument('--Tm-opt', type=int, help='Optmized melting temperature for primers (℃).', \
        default=60)
    group_design.add_argument('--primer-num-return', type=int, help='The maximum number of primers to return in Primer3 \
        designing results.', default=30)
    # qRT-PCR specific
    group_design.add_argument('--junction', help='Primer pair must be separated by at least one intron on the \
        corresponding genomic DNA; or primers must span an exon-exon junction. Junction data in JSON format should be prepared by \
            the command-line primertool-junctions', action='store_true')
    
    # These arguments are used by check and full
    parent_parser_check = argparse.ArgumentParser(add_help=False)
    group_check = parent_parser_check.add_argument_group('Check Specificity')
    group_check.add_argument('-3', '--use-3-end', help='If turned on, primer pairs having at least one mismatch at the 3 end\
        position with templates would not be considered to produce off-target amplicon, even if their melting temperatures \
            are high. Turn on this would find more candidate primers, but might also have more false positives\
                ', action='store_true')
    group_check.add_argument('--checking-size-min', type=int, help='Lower limit of the checking amplicon size range (bp).', \
        default=50)
    group_check.add_argument('--checking-size-max', type=int, help='Upper limit of the checking amplicon size range (bp).', \
        default=2000)
    group_check.add_argument('--amplicon-num-max', type=int, help='The maximum number of amplicons for checking.', \
        default=10)
    group_check.add_argument('-a', '--report-amplicon-seqs', help="Get amplicon seqs (might be slow)", action='store_true')
    # qRT-PCR specific
    group_check.add_argument('--isoform', help="Allow primers targeting on alternative isoforms and still regard them \
        as specific ones. Isoform data in JSON format should be prepared by the command-line primertool-isoforms", action='store_true')

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
        raise Exception('No samtools detected in your system')

    samtools_version = os.popen('samtools --version').readlines()[0].strip().split(' ')[1]
    if LooseVersion(samtools_version) < LooseVersion('1.9'):
        raise Exception(f'Your samtools version is v{samtools_version}, but >=v1.9 is required')

    if shutil.which('blastn') is None:
        raise Exception('No NCBI-BLAST+ (blastn) detected in your system')

    if shutil.which('makeblastdb') is None:
        raise Exception('No NCBI-BLAST+ (makeblastdb) detected in your system')

def check_templates(args):
    for (i, template) in enumerate(args.templates.split(',')):
        if os.path.isfile(template) is False:
            raise Exception(f'File not found: {template}')
        if os.path.isfile(template+'.fai') is False:
            raise Exception(f'File {template} has not been indexed by samtools faidx.')
        if os.path.isfile(re.sub('[.]gz$', '', template)+'.nhr') is False and os.path.isfile(re.sub('\.gz$', '', template)+'.nal') is False:
            raise Exception(f'File {template} has not been indexed by makeblastdb.')
        
        # qPCR specific
        if i==0:
            if os.path.isfile(template+'.junctions.json') is False and 'junction' in args and args.junction is True:
                raise Exception(f'The junction file for {template} is not ready. Parameter --junction is not allowed')
            if os.path.isfile(template+'.isoforms.json') is False and 'isoform' in args and args.isoform is True:
                raise Exception(f'The isoform file for {template} is not ready. Parameter --isoform is not allowed')

def check_qPCR(args):
    if args.run_mode!='check' and args.junction is True and args.type!='SEQUENCE_INCLUDED_REGION':
        raise Exception('Parameter --junction must be used with --type=SEQUENCE_INCLUDED_REGION')

def run(args):
    ###################  Design primers ###################
    query_string = args.query.read()
    dbs = args.templates.split(',')
    if args.run_mode=='check':
        primers = make_primers.make_primers(query=query_string)
        if 'error' in primers:
            raise Exception(primers['error'])
    else:
        if args.run_mode=='design':
            primer_num_return = args.primer_num_retain
        else:
            primer_num_return = args.primer_num_return
        if make_sites.judge_input_type(query_string)=='pos':
            sites = make_sites.build_by_pos(query=query_string, template_file=dbs[0], primer_type=args.type, \
                primer_num_return=primer_num_return, size_min=args.product_size_min, size_max=args.product_size_max, \
                    pick_internal=args.pick_oligo, use_junction=args.junction, Tm_opt=args.Tm_opt)
        else:
            sites = make_sites.build_by_seq(query=query_string, primer_type=args.type, \
                primer_num_return=primer_num_return, size_min=args.product_size_min, size_max=args.product_size_max, \
                    pick_internal=args.pick_oligo, Tm_opt=args.Tm_opt)
        if 'error' in sites:
            raise Exception(sites['error'])
        primers = design_primer.multiple(sites, cpu=args.cpu)

    ###################  Checking specificity  #############
    if args.run_mode!='design':
        primers = run_blast.run_blast_parallel(primers=primers, dbs=dbs, cpu=args.cpu,\
            checking_size_max=args.checking_size_max, checking_size_min=args.checking_size_min, \
                report_amplicon_seq=args.report_amplicon_seqs, Tm_diff=args.Tm_diff, use_3_end=args.use_3_end, \
                    max_amplicon=args.amplicon_num_max)
        primers = sort_primers.sort_rank(primers=primers, dbs=dbs, max_num_return=args.primer_num_retain, use_isoforms=args.isoform)

    ###################  Checking multiplex  ###############
    dimers = {}
    if args.check_multiplex is True:
        dimers = multiplex.extract_fake_pair(primers, Tm_diff=args.Tm_diff, cpu=args.cpu)

    ###################  Output  ###########################
    region_type = args.type if 'type' in args else 'NA'
    if args.out is None:
        output_fh = open(f'{args.query.name}.json', mode='w')
    else:
        output_fh = args.out
    print(json.dumps({'meta':{'mode':args.run_mode, 'dbs':dbs, 'region_type': region_type, 'check_multiplex':args.check_multiplex}, \
        'primers':primers, 'dimers':dimers}, indent=4), file=output_fh)

    if args.tsv is None:
        tsv_fh = open(f'{args.query.name}.tsv', mode='w')
    else:
        tsv_fh = args.tsv
    print(output.tsv(primers, dbs, args.run_mode), file=tsv_fh)
    if args.check_multiplex is True:
        print(output.dimer_list(dimers), file=tsv_fh)
    output_fh.close()
    tsv_fh.close()


def main():
    global_var.init()
    args = make_args()
    check_environments(args)
    check_templates(args)
    check_qPCR(args)
    run(args)

if __name__ == "__main__":
    main()

