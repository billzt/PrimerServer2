import re
import json
import os
import sys
import subprocess
import multiprocessing as mp
import time

import primer3
import progressbar

from primerserver2.core.analysis_blast import filter_len, filter_Tm, add_amplicon_seq, add_isoform_annotation
from primerserver2.core.make_sites import faidx
from primerserver2.core.make_primers import make_primers
from primerserver2.core import global_var

def run_blast(p3_inputs):
    '''
        run blastn for a single primer pair, then use analysis_blast.py to obtain amplicons within Tm thredsholds
    '''
    if global_var.stop_run is True:
        return {'db': '', 'amplicons': {}}

    if len(p3_inputs)==0 :  # unfortunately, these sites have no primers
        return {'db': '', 'amplicons': {}}
    db = p3_inputs[0]['db']  # dbs are the same
    blast_query = ''
    query_primer_seq_dict = {}  # for the func: filter_Tm
    for x in p3_inputs:
        id = x['id']
        rank = x['rank']
        seq_L = x['seq_L']
        seq_R = x['seq_R']
        blast_query += f'>{id}.{rank}.L\n{seq_L}\n>{id}.{rank}.R\n{seq_R}\n'
        query_primer_seq_dict[f'{id}.{rank}.L'] = seq_L
        query_primer_seq_dict[f'{id}.{rank}.R'] = seq_R

    blast_db = re.sub('[.]gz$', '', db)
    if os.path.isfile(blast_db+'.nhr')==False and os.path.isfile(blast_db+'.nal')==False:
        raise Exception(f'The database file is not complete: file {blast_db}.nhr or {blast_db}.nal is not found')
    cmd = f'blastn -task blastn-short -db {blast_db} -evalue 30000 -word_size 7 ' \
        + '-perc_identity 60 -dust no -reward 1 -penalty -1 -max_hsps 500 -ungapped ' \
        + ' -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send sstrand"'
    blast_out = subprocess.run(cmd, input=blast_query, stdout=subprocess.PIPE, shell=True, encoding='ascii').stdout
    amplicons = filter_len(blast_out=blast_out, len_min=p3_inputs[0]['checking_size_min'], \
        len_max=p3_inputs[0]['checking_size_max'])
    hits_seqs = faidx(template_file=db, region_string=amplicons['regions_primer'])
    report_amplicons = filter_Tm(amplicons['amplicons'], query_primer_seq=query_primer_seq_dict, \
        hits_seqs=hits_seqs, Tm_diff=p3_inputs[0]['Tm_diff'], use_3_end=p3_inputs[0]['use_3_end'], \
            max_amplicons=p3_inputs[0]['max_amplicon'])
    if p3_inputs[0]['report_amplicon_seq'] is True:
        report_amplicons = add_amplicon_seq(amplicons=report_amplicons, template_file=db)
    return {'db': os.path.basename(db), 'amplicons': report_amplicons}

def run_blast_parallel(primers, dbs, cpu=2, checking_size_min=70, checking_size_max=1000, \
    report_amplicon_seq=False, Tm_diff=20, use_3_end=False, monitor=True, max_amplicon=10):
    if global_var.stop_run is True:
        return {'error': 'Stop running'}

    # distribute tasks
    pool = mp.Pool(processes=cpu)
    multi_res = []
    global_var.all_tasks_num = 0
    global_var.current_task = 'blast'
    for db in dbs:
        p3_inputs = []
        for (id, primer) in primers.items():
            primer_num = primer['PRIMER_PAIR_NUM_RETURNED']
            for rank in range(0, primer_num):
                p3_inputs.append({
                    'id': id,
                    'rank': rank,
                    'db': db,
                    'seq_L': primer[f'PRIMER_LEFT_{rank}_SEQUENCE'],
                    'seq_R': primer[f'PRIMER_RIGHT_{rank}_SEQUENCE'],
                    'checking_size_min': checking_size_min,
                    'checking_size_max': checking_size_max,
                    'report_amplicon_seq': report_amplicon_seq,
                    'Tm_diff': Tm_diff,
                    'use_3_end': use_3_end,
                    'max_amplicon': max_amplicon
                })
        for i in range(0, int(len(p3_inputs)/5)+1):
            sub_start = i*5
            sub_end = min(len(p3_inputs), sub_start+5)
            global_var.all_tasks_num += 5
            multi_res.append(pool.apply_async(run_blast, (p3_inputs[sub_start:sub_end],)))

    # monitor
    if monitor is True:
        widgets = ['Checking specificity: ', progressbar.Counter(),\
            ' Finished', ' (', progressbar.Percentage(), ')', \
                progressbar.Bar(), progressbar.ETA()]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=global_var.all_tasks_num).start()

    while True:
        global_var.complete_count = sum([1 for x in multi_res if x.ready()])*5
        if global_var.complete_count == global_var.all_tasks_num:
            if monitor is True:
                bar.finish()
            break
        if monitor is True:
            bar.update(global_var.complete_count)
        time.sleep(1)
        
    # Results
    for result in multi_res:
        result_data = result.get()  # db and amplicons
        db = result_data['db']
        if db== '': # unfortunately, these sites have no primers
            continue
        amplicons = result_data['amplicons']
        for id in amplicons.keys():
            if db not in primers[id]:
                primers[id][db] = {}
            for rank in amplicons[id].keys():
                primers[id][db][f'PRIMER_PAIR_{rank}_AMPLICONS'] = amplicons[id][rank]

    return primers

if __name__ == "__main__":
    global_var.init()
    with open('tests/query_check_multiple') as f:
        primers = make_primers(query=f.read())
        dbs = ['tests/example.fa']
        report_primers = run_blast_parallel(primers, dbs, cpu=10, report_amplicon_seq=True)
        print(json.dumps(report_primers, indent=4))
