#!/usr/bin/env python3

'''PrimerServer2: a high-throughput primer design and specificity-checking platform
Github: https://github.com/billzt/PrimerServer2
This is a CLI script to configure PrimerServer2
'''

import argparse
import os
import sys
import json
import re

def load():
    home_dir = os.environ['HOME']
    web_config = json.load(open(f'{home_dir}/.primerserver.json'))
    return web_config

def prepare():
    home_dir = os.environ['HOME']
    if os.path.isfile(f'{home_dir}/.primerserver.json') is False:
        data_path = os.path.join(os.path.dirname(__file__), '../data')
        os.system(f"cp {data_path}/web_config_sample.json {home_dir}/.primerserver.json")
    print(f'Configure file {home_dir}/.primerserver.json is ready. Please modify it')

def check():
    home_dir = os.environ['HOME']
    status = 'error'
    msg = ''
    if os.path.isfile(f'{home_dir}/.primerserver.json') is False:
        msg = f'The configure file {home_dir}/.primerserver.json is not found. Perhaps you are in a fresh \
            installation? Please run "primerserver-config" to load and edit your configure files'
        return {'status': status, 'msg': msg}
    web_config = json.load(open(f'{home_dir}/.primerserver.json'))
    for key in ['cpu', 'templates_directory', 'templates']:
        if key not in web_config:
            msg = f'The configure file {home_dir}/.primerserver.json is broken. It lacks the "{key}" key.'
            return {'status': status, 'msg': msg}
    if isinstance(web_config['cpu'], int) is False:
        msg = f'There is some error in the configure file {home_dir}/.primerserver.json: the value for cpu must \
            be int.'
        return {'status': status, 'msg': msg}
    db_dir = web_config['templates_directory']
    if os.path.isdir(db_dir) is False:
        msg = f'There is some error in the configure file {home_dir}/.primerserver.json: the path for \
            templates_directory ({db_dir}) does not exist'
        return {'status': status, 'msg': msg}
    for (dbname, dbinfo) in web_config['templates'].items():
        if os.path.isfile(f'{db_dir}/{dbname}') is False:
            if dbname=='example.fa' or dbname=='example2.fa':
                data_path = os.path.join(os.path.dirname(__file__), '../data')
                os.system(f"cp {data_path}/{dbname}* {db_dir}/")
            else:
                msg = f'File {db_dir}/{dbname} does not exist'
                return {'status': status, 'msg': msg}
        if os.path.isfile(re.sub('[.]gz$', '', f'{db_dir}/{dbname}')+'.nhr') is False and os.path.isfile(re.sub('[.]gz$', '', f'{db_dir}/{dbname}')+'.nal') is False:
            msg = f'File {db_dir}/{dbname} has no index. Use "makeblastdb" to index {db_dir}/{dbname}'
            return {'status': status, 'msg': msg}
        if os.path.isfile(f'{db_dir}/{dbname}.fai') is False:
            msg = f'File {db_dir}/{dbname}.fai does not exist. Use "samtools faidx" to index {db_dir}/{dbname}'
            return {'status': status, 'msg': msg}
        for key in ['IDs', 'description', 'group']:
            if key not in dbinfo:
                msg = f'The configure file {home_dir}/.primerserver.json is broken. It lacks the "{key}" key.'
                return {'status': status, 'msg': msg}
            if isinstance(dbinfo[key], str) is False:
                msg = f'There is some error in the configure file {home_dir}/.primerserver.json: the value for {key} must \
                    be strings'
                return {'status': status, 'msg': msg}
    status = 'success'
    return {'status': status, 'msg': msg}
    
if __name__ == "__main__":
    print(json.dumps(check()))
