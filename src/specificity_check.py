import re
import json
import os
import multiprocessing as mp

import primer3

from analysis_blast import filter_len, filter_Tm
from make_sites import faidx
from make_primers import make_primers

def run_blast(query_and_db):
    primer = query_and_db['primer']
    db = query_and_db['db']
    query_file = query_and_db['query_file']
    if os.path.isfile(db+'.nhr')==False:
        raise Exception(f'The database file is not complete: file {db}.nhr is not found')
    with open(query_file, 'w') as f:
        for (primer_id, primer_data) in primer.items(): # only one
            primer_num = primer_data['PRIMER_PAIR_NUM_RETURNED']
            if primer_num>0:
                for j in range(0, primer_num):
                    print(f'>{primer_id}.{j}.Primer.LEFT\n'+primer_data[f'PRIMER_LEFT_{j}_SEQUENCE']+'\n', file=f)
                    print(f'>{primer_id}.{j}.Primer.RIGHT\n'+primer_data[f'PRIMER_RIGHT_{j}_SEQUENCE']+'\n', file=f)
    
    cmd = f'blastn -task blastn-short -query {query_file} -db {db} -evalue 30000 -word_size 7' \
        + '-perc_identity 60 -dust no -reward 1 -penalty -1 -max_hsps 500 -ungapped' \
        + ' -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send sstrand"'
    blast_out = os.popen(cmd).read()
    amplicons = filter_len(blast_out=blast_out, len_min=query_and_db['checking_size_min'], \
        len_max=query_and_db['checking_size_max'], region_file=f'{query_file}'+'.regions')
    hits_seqs = faidx(template_file=db, region_file=f'{query_file}'+'.regions')
    report_primers = filter_Tm(amplicons, query_primers=primer, hits_seqs=hits_seqs)
    return report_primers

def run_blast_parallel(primers, dbs, cpu=2):
    pool = mp.Pool(processes=cpu)
    multi_res = []
    for primer in primers:
        for db in dbs:
            query_and_db = {
                'primer': primer,
                'db': 'tests/'+db,
                'query_file': 'tests/'+primer+db,
                'checking_size_min': 70,
                'checking_size_max': 1000
            }
            multi_res.append(pool.apply_async(run_blast, (query_and_db,)))
    report_primers = {}
    for result in multi_res:
        report_primers.update(result.get())
    return report_primers

if __name__ == "__main__":
    primers = make_primers('tests/query_check')
    dbs = ['example.fa']
    #report_primers = run_blast_parallel(primers, dbs, cpu=10)
    
    primer = primers[0]
    db = 'example.fa'
    report_primers = run_blast({
        'primer': primer,
        'db': 'tests/'+db,
        'query_file': 'tests/'+primer+db,
        'checking_size_min': 70,
        'checking_size_max': 1000
    })
    print(json.dumps(report_primers))

# analysis_blast 的对象可以不要绑定id，因为只有一个组
