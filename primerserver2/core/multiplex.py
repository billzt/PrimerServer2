import json
import multiprocessing as mp
import time

from itertools import product
from collections import defaultdict

import primer3
import progressbar

from primerserver2.core import global_var

def judge_two_site(site_pair):

    if global_var.stop_run is True:
        return []

    dimers = []
    for (rank_1, primers_1) in site_pair['site_1'].items():
        for (rank_2, primers_2) in site_pair['site_2'].items():
            for (seq_1, seq_2) in product([primers_1['seq_L'], primers_1['seq_R']], [primers_2['seq_L'], primers_2['seq_R']]):
                min_Tm = min(primers_1['min_Tm'], primers_2['min_Tm'])-site_pair['Tm_diff']

                Heterodimer = primer3.calcHeterodimer(seq_1, seq_2, output_structure=True)
                if Heterodimer.tm > min_Tm:
                    dimers.append([site_pair['id_1'], site_pair['id_2'], rank_1, rank_2, seq_1, seq_2, Heterodimer.tm, Heterodimer.ascii_structure])
                    break
                EndStability = primer3.bindings.calcEndStability(seq_1, seq_2)
                if EndStability.tm > min_Tm:
                    dimers.append([site_pair['id_1'], site_pair['id_2'], rank_1, rank_2, seq_1, seq_2, EndStability.tm, ''])
                    break
    return dimers


def extract_fake_pair(primers_dict, Tm_diff=10, cpu=2, monitor=True):

    if global_var.stop_run is True:
        return {'error': 'Stop running'}

    # make sure more than two sites has primers
    site_num = 0
    for (id, primers) in primers_dict.items():
        if primers['PRIMER_PAIR_NUM_RETURNED']>0:
            site_num += 1
    if site_num < 2:
        return {'error': f'Only {site_num} site has primers. Skip the multiplex module'}
    
    # data group
    primer_sites = defaultdict(dict)
    for (id, primers) in primers_dict.items():
        if primers['PRIMER_PAIR_NUM_RETURNED']==0:
            continue
        if 'PRIMER_PAIR_NUM_RETURNED_FINAL' in primers:
            primer_num = primers['PRIMER_PAIR_NUM_RETURNED_FINAL']
        else:
            primer_num = primers['PRIMER_PAIR_NUM_RETURNED']
        for amplicon_rank in range(0, primer_num):
            if 'PRIMER_PAIR_NUM_RETURNED_FINAL' in primers:
                raw_rank = primers[f'PRIMER_PAIR_AMPLICON_NUM_RANK_{amplicon_rank}']
            else:
                raw_rank = amplicon_rank
            seq_L = primers[f'PRIMER_LEFT_{raw_rank}_SEQUENCE']
            seq_R = primers[f'PRIMER_RIGHT_{raw_rank}_SEQUENCE']
            min_Tm = min(primers[f'PRIMER_LEFT_{raw_rank}_TM'], primers[f'PRIMER_RIGHT_{raw_rank}_TM'])
            primer_sites[id][amplicon_rank] = {"seq_L": seq_L, "seq_R": seq_R, "min_Tm": min_Tm}

    # distribute tasks
    pool = mp.Pool(processes=cpu)
    multi_res = []
    global_var.all_tasks_num = 0
    global_var.current_task = 'multiplex'
    primer_ids = list(primer_sites.keys())
    for i in range(0, len(primer_ids)):
        id_1 = primer_ids[i]
        for j in range(i+1, len(primer_ids)):
            id_2 = primer_ids[j]
            site_pair = {'id_1':id_1, 'site_1':primer_sites[id_1], 'id_2':id_2, 'site_2':primer_sites[id_2], 'Tm_diff':Tm_diff}
            global_var.all_tasks_num += 1
            multi_res.append(pool.apply_async(judge_two_site, (site_pair,)))

    # monitor
    if monitor is True:
        widgets = ['Checking Multiplex: ', progressbar.Counter(),\
            ' Finished', ' (', progressbar.Percentage(), ')', \
                progressbar.Bar(), progressbar.ETA()]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=global_var.all_tasks_num).start()

    while True:
        global_var.complete_count = sum([1 for x in multi_res if x.ready()])
        if global_var.complete_count == global_var.all_tasks_num:
            if monitor is True:
                bar.finish()
            break
        if monitor is True:
            bar.update(global_var.complete_count)
        time.sleep(0.1)

    # results
    dimers_dict = defaultdict(dict)
    for result in multi_res:
        for dimer in result.get():
            (id_1, id_2, rank_1, rank_2, seq_1, seq_2, tm, ascii_structure) = dimer
            dimers_dict[id_1][rank_1] = defaultdict(dict)
            dimers_dict[id_1][rank_1][id_2] = defaultdict(dict)
            dimers_dict[id_1][rank_1][id_2][rank_2] = {'seq_1':seq_1, 'seq_2':seq_2, 'Tm':round(tm, 2), 'ascii_structure':ascii_structure}
    return dimers_dict

if __name__ == "__main__":
    global_var.init()
    primers_dict = json.load(open('tests/_internal_/multisites.json'))
    dimers_dict = extract_fake_pair(primers_dict['primers'])
    print(json.dumps(dimers_dict, indent=4))