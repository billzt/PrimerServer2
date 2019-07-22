import re
import json
import os
import subprocess
import multiprocessing as mp

import primer3

from primerserver2.core.analysis_blast import filter_len, filter_Tm, add_amplicon_seq
from primerserver2.core.make_sites import faidx
from primerserver2.core.make_primers import make_primers

def run_blast(p3_input):
    '''
        run blastn for a single primer pair, then use analysis_blast.py to obtain amplicons within Tm thredsholds
    '''
    seq_L = p3_input['seq_L']
    seq_R = p3_input['seq_R']
    db = p3_input['db']
    if os.path.isfile(db+'.nhr')==False:
        raise Exception(f'The database file is not complete: file {db}.nhr is not found')
    blast_query = f'>LEFT\n'+seq_L+'\n' + f'>RIGHT\n'+seq_R+'\n'
    cmd = f'blastn -task blastn-short -db {db} -evalue 30000 -word_size 7 ' \
        + '-perc_identity 60 -dust no -reward 1 -penalty -1 -max_hsps 500 -ungapped ' \
        + ' -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send sstrand"'
    blast_out = subprocess.run(cmd, input=blast_query, stdout=subprocess.PIPE, shell=True, encoding='ascii').stdout
    amplicons = filter_len(blast_out=blast_out, len_min=p3_input['checking_size_min'], len_max=p3_input['checking_size_max'])
    hits_seqs = faidx(template_file=db, region_string=amplicons['regions_primer'])
    report_amplicons = filter_Tm(amplicons['amplicons'], query_primer_seq={'LEFT': seq_L, 'RIGHT': seq_R}, hits_seqs=hits_seqs)
    if p3_input['report_amplicon_seq']==True:
        report_amplicons = add_amplicon_seq(amplicons=report_amplicons, template_file=db)
    return {'id': p3_input['id'], 'rank': p3_input['rank'], 'db': os.path.basename(db), 'amplicons': report_amplicons}

def run_blast_parallel(primers, dbs, cpu=2, checking_size_min=70, checking_size_max=1000, report_amplicon_seq=False):
    pool = mp.Pool(processes=cpu)
    multi_res = []
    for (id, primer) in primers.items():
        primer_num = primer['PRIMER_PAIR_NUM_RETURNED']
        for rank in range(0, primer_num):
            for db in dbs:
                p3_input = {
                    'id': id,
                    'rank': rank,
                    'db': db,
                    'seq_L': primer[f'PRIMER_LEFT_{rank}_SEQUENCE'],
                    'seq_R': primer[f'PRIMER_RIGHT_{rank}_SEQUENCE'],
                    'checking_size_min': checking_size_min,
                    'checking_size_max': checking_size_max,
                    'report_amplicon_seq': report_amplicon_seq
                }
                multi_res.append(pool.apply_async(run_blast, (p3_input,)))
    
    for result in multi_res:
        result_data = result.get()
        id = result_data['id']
        rank = result_data['rank']
        db = result_data['db']

        if db not in primers[id]:
            primers[id][db] = {}
        primers[id][db][f'PRIMER_PAIR_{rank}_AMPLICONS'] = result_data['amplicons']

    return primers

if __name__ == "__main__":
    with open('tests/query_check_multiple') as f:
        primers = make_primers(query=f.read())
        dbs = ['tests/example.fa']
        report_primers = run_blast_parallel(primers, dbs, cpu=10, report_amplicon_seq=True)
        print(json.dumps(report_primers, indent=4))
