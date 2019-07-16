import re
import json
import os
import multiprocessing as mp

import primer3

def calculate_GC(seq):
    gc = len([x for x in seq.upper() if x=='G' or x=='C'])
    return gc/len(seq)*100

def make_primers(query_file):
    '''
        In case of user's input, make list of primer dicts similar as the result of design_primer module
    '''
    primers = []
    with open(query_file) as f:
        for line in f:
            (id, seq_F, seq_R) = re.split(r'\s+', line.strip())
            primers.append({'id': id, 'primers': [{
                'PRIMER_PAIR_NUM_RETURNED': 1,
                'PRIMER_PAIR_0_PENALTY': 0,
                'PRIMER_LEFT_0_SEQUENCE': seq_F,
                'PRIMER_RIGHT_0_SEQUENCE': seq_R,
                'PRIMER_LEFT_0': [-1, len(seq_F)],
                'PRIMER_RIGHT_0': [-1, len(seq_R)],
                'PRIMER_LEFT_0_TM': primer3.calcTm(seq_F),
                'PRIMER_RIGHT_0_TM': primer3.calcTm(seq_R),
                'PRIMER_LEFT_0_GC_PERCENT': calculate_GC(seq_F),
                'PRIMER_RIGHT_0_GC_PERCENT': calculate_GC(seq_R),
                'PRIMER_PAIR_0_PRODUCT_SIZE': -1
            }]})
    
    return primers

def run_blast(query_and_db):
    primer = query_and_db['primer']
    db = query_and_db['db']
    query_file = query_and_db['query_file']
    if os.path.isfile(db+'.nhr')==False:
        raise Exception(f'The database file is not complete: file {db}.nhr is not found')
    with open(query_file, 'w') as f:
        primer_id = primer['id']
        for i in primer['primers']: # usually there's only one item in primer['primers'], except for FORCE_END
            primer_num = i['PRIMER_PAIR_NUM_RETURNED']
            if primer_num>0:
                for j in range(0, primer_num-1):
                    print(f'>{primer_id}.{j}.Primer0\n'+i[f'PRIMER_LEFT_{j}_SEQUENCE']+'\n', file=f)
                    print(f'>{primer_id}.{j}.Primer1\n'+i[f'PRIMER_RIGHT_{j}_SEQUENCE']+'\n', file=f)
    
    cmd = f'blastn -task blastn-short -query {query_file} -db {db} -evalue 30000 -word_size 7' \
        + '-perc_identity 60 -dust no -reward 1 -penalty -1 -max_hsps 500 -ungapped' \
        + ' -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send sstrand"'
    os.popen(cmd)

def run_blast_parallel(primers, dbs, cpu=2):
    pool = mp.Pool(processes=cpu)
    multi_res = []
    for primer in primers:
        for db in dbs:
            query_and_db = {
                'primer': primer,
                'db': db,
                'query_file': primer+db
            }
            multi_res.append(pool.apply_async(run_blast, (query_and_db,)))
    return [x.get() for x in multi_res]

if __name__ == "__main__":
    print('GC Content: AAGGCAGGCTAAGGAGGAGA:', calculate_GC('AAGGCAGGCTAAGGAGGAGA'))
    print('make_primers:')
    print(json.dumps(make_primers('tests/query_2')))