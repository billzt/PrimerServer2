#!/usr/bin/env python3

'''PrimerServer2: a high-throughput primer design and specificity-checking platform
Github: https://github.com/billzt/PrimerServer2
This is a CLI script to configure PrimerServer2:
    (1) pre-load databases (by makeblastdb in NCBI BLAST+)
    (2) define numbers of CPUs to use
'''

import argparse
import os
import sys
import re
import json

from primerserver2.cmd.primertool import check_environments

def make_args():
    parser = argparse.ArgumentParser(description='primerserver-config: configure PrimerServer2', \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--template', help='a template file in FASTA format.')
    parser.add_argument('-d', '--description', help='a human readable string describe your template FASTA file. \
        default: The same as -t / --template')
    parser.add_argument('-g', '--group', help='Assigning this template to a group.', default='group')
    parser.add_argument('-p', '--cpu', type=int, help='Used CPU number. default: 2')
    args = parser.parse_args()
    return args

def load():
    web_config = json.load(open(os.path.join(os.path.dirname(__file__), '../data/web_config.json')))
    return web_config

def main():
    # init conf
    if os.path.isfile(os.path.join(os.path.dirname(__file__), '../data/web_config.json')) is False:
        os.system(f"cp {os.path.join(os.path.dirname(__file__), '../data/web_config_sample.json')} \
            {os.path.join(os.path.dirname(__file__), '../data/web_config.json')}")

    web_config = load()

    # args
    args = make_args()
    check_environments(args)
    db_dir = os.path.join(os.path.dirname(__file__), '../templates/')

    # template, junction, isoform, description, group
    if args.template is not None:
        db = args.template
        print(f'Begin configuring template... {db}', file=sys.stderr)

        # check templates
        if os.path.isfile(db) is False:
            raise Exception(f'File not found: {db}')
        if os.path.isfile(db+'.fai') is False:
            code = os.system(f'samtools faidx {db} 2>/dev/null')
            if code != 0:
                raise Exception(f'File {db} cannot be indexed by samtools faidx. Perhaps it is not in FASTA format')
        if os.path.isfile(db+'.nhr') is False:
            code = os.system(f'makeblastdb -dbtype nucl -in {db} 2>/dev/null')
            if code != 0:
                raise Exception(f'File {db} cannot be indexed by makeblastdb.')

        # db load
        db_basename = os.path.basename(db)

        # Name can't be example.fa
        if db_basename=='example.fa' or db_basename=='example2.fa':
            raise Exception('name: example.fa is not allowed')
        
        # Index and Copy
        if os.path.isfile(db_dir+'/'+db_basename):
            print(f'{db} has already existed in PrimerServer. Skip copying BLAST databases', file=sys.stderr)
        else:
            # Index
            if os.path.isfile(db+'.nhr') is False:
                print(f'{db} has not been indexed. Start indexing it...', file=sys.stderr)
                code = os.system(f'makeblastdb -dbtype nucl -in {db} 2>/dev/null')
                if code != 0:
                    raise Exception(f'Error in indexing {db}. Perhaps it is not a valid DNA file in FASTA format')
                print(f'Finish indexing it {db}', file=sys.stderr)
            
            # Copy
            print(f'Start copying {db} to PrimerServer...', file=sys.stderr)
            code = os.system(f'cp {db} {db}.fai {db}.nhr {db}.nin {db}.nsq {db_dir} 2>/dev/null')
            if code != 0:
                raise Exception(f'Error in copying {db}. Make sure you have write permisson to {db_dir} \
                    and the space is enough', file=sys.stderr)
        
        # Add to config
        if db_basename not in web_config['templates']:
            web_config['templates'][db_basename] = {}
        web_config['templates'][db_basename] = {'description': db_basename, 'group': args.group, 'junction': False, 'isoform': False}

        # junction and isoform
        if os.path.isfile(db+'.isoforms.json') is True and os.path.isfile(db_dir+'/'+db+'.isoforms.json') is not True:
            code = os.system(f'cp {db}.isoforms.json {db_dir} 2>/dev/null')
            web_config['templates'][db_basename]['junction'] = True
        if os.path.isfile(db+'.junctions.json') is True and os.path.isfile(db_dir+'/'+db+'.junctions.json') is not True:
            code = os.system(f'cp {db}.junctions.json {db_dir} 2>/dev/null')
            web_config['templates'][db_basename]['isoform'] = True
        
        # description and group
        if args.description is not None:
            web_config['templates'][db_basename]['description'] = args.description
        if args.group is not None:
            web_config['templates'][db_basename]['group'] = args.group
        

    # CPUs
    if args.cpu is not None:
        print('Begin configuring CPUs...', file=sys.stderr)
        if args.cpu<0:
            raise Exception('CPU numbers must be larger than 0')
        if args.cpu>os.cpu_count()-1:
            raise Exception(f'CPU numbers cannot be larger than your system resources: {os.cpu_count()}')
        web_config['cpu'] = args.cpu

    # IDs
    for db_file in web_config['templates']:
        if 'IDs' not in web_config['templates'][db_file]:
            ids = os.popen(f'sort -k2,2nr {db_dir}/{db_file}.fai | cut -f 1 | head -n 3').read().splitlines()
            web_config['templates'][db_file]['IDs'] = '; '.join(ids)
    
    # write configure file
    with open(os.path.join(os.path.dirname(__file__), '../data/web_config.json'), 'w') as f:
        json.dump(web_config, f, indent=4)
        print('Configure finished', file=sys.stderr)


if __name__ == "__main__":
    main()
    


