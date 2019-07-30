#!/usr/bin/env python3

'''PrimerServer2: a high-throughput primer design and specificity-checking platform
Github: https://github.com/billzt/PrimerServer2
This is a CLI script to configure PrimerServer2:
    (1) pre-load databases (by makeblastdb in NCBI BLAST+)
    (2) define numbers of CPUs to use

    The --description-file (-d) file format (TSV):
    template1.fa(dbname, without path)      Name1         (Optional)GroupA
    template2.fa(dbname, without path)      Name2         (Optional)GroupA
    template3.fa(dbname, without path)      Name3         (Optional)GroupB
'''

import argparse
import os
import sys
import re
import json

from primerserver2.cmd.primertool import check_environments, check_templates

def make_args():
    parser = argparse.ArgumentParser(description='primerserver-config: configure PrimerServer2', \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--json-debug', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('-t', '--templates', help='template file in FASTA format. \
        Allowing multiple files (separated by comma)')
    parser.add_argument('-d', '--description-file', help='a TSV file including descriptions for the templates', \
        type=argparse.FileType('r'))
    parser.add_argument('-p', '--cpu', type=int, help='Used CPU number.')
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

    # templates
    if args.templates is not None:
        print('Begin configuring templates...', file=sys.stderr)

        check_templates(args)

        # db load
        for db in args.templates.split(','):
            db_basename = os.path.basename(db)

            # Name can't be example.fa
            if db_basename=='example.fa':
                raise Exception('name: example.fa is not allowed')
            
            # Index and Copy
            if os.path.isfile(db_dir+'/'+db_basename):
                print(f'{db} has already existed in PrimerServer. Skip it', file=sys.stderr)
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
            web_config['templates'][db_basename] = {'description': db_basename, 'group': 'group'}

        
    # db descriptions
    if args.description_file is not None:
        print('Begin configuring template descriptions...', file=sys.stderr)
        for line in args.description_file:
            (db_file, *db_desc) = line.strip().split('\t')
            if os.path.isfile(db_dir+'/'+db_file) is False:
                print(f'template {db_file} has not been loaded into PrimerServer. Skip this', file=sys.stderr)
            if db_file not in web_config['templates']:
                web_config['templates'][db_file] = {}
            if len(db_desc)==2:
                (name, group) = db_desc
                web_config['templates'][db_file] = {'description': name, 'group': group}
            elif len(db_desc)==1:
                web_config['templates'][db_file] = {'description': db_desc[0], 'group': 'group'}
            else:
                raise Exception(f'the file {args.description_file.name} does not have two or three columns seprated by tabs')
    
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
    


