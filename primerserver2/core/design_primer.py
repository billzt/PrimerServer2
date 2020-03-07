
import json
import os
import time
import multiprocessing as mp
import sys

import primer3
import progressbar

from primerserver2.core import global_var

def single(site):
    '''
        Using the primer3-module to design primers for a single site, return a list (one or two items)
        in DICT format

        Site format: 
        site['id'], site['template'], site['type'], site['pos'], site['length'], 
        site['size_max'], site['size_min'], site['primer_num_return'], site['pick_internal']
    '''
    if global_var.stop_run is True:
        return []

    id = site['id']
    template = site['template']
    type = site['type']
    pos = site['pos']
    length = site['length']
    size_max = site['size_max']
    size_min = site['size_min']
    primer_num_return = site['primer_num_return']
    Tm_opt = site['Tm_opt'] if 'Tm_opt' in site else 60

    p3_json = os.path.join(os.path.dirname(__file__), '../data/p3_settings.json')
    p3_settings = dict(json.load(open(p3_json)))
    p3_settings['PRIMER_PRODUCT_SIZE_RANGE'] = [[size_min, size_max]]
    p3_settings['PRIMER_NUM_RETURN'] = primer_num_return
    p3_settings['PRIMER_OPT_TM'] = Tm_opt
    p3_settings['PRIMER_MIN_TM'] = Tm_opt-3
    p3_settings['PRIMER_MAX_TM'] = Tm_opt+3
    p3_settings['PRIMER_INTERNAL_OPT_TM'] = Tm_opt+10
    p3_settings['PRIMER_INTERNAL_MIN_TM'] = p3_settings['PRIMER_INTERNAL_OPT_TM']-3
    p3_settings['PRIMER_INTERNAL_MAX_TM'] = p3_settings['PRIMER_INTERNAL_OPT_TM']+3
    p3_settings['PRIMER_MAX_SELF_ANY_TH'] = Tm_opt-15
    p3_settings['PRIMER_PAIR_MAX_COMPL_ANY_TH'] = Tm_opt-15
    p3_settings['PRIMER_MAX_SELF_END_TH'] = Tm_opt-15
    p3_settings['PRIMER_PAIR_MAX_COMPL_END_TH'] = Tm_opt-15
    p3_settings['PRIMER_MAX_HAIRPIN_TH'] = Tm_opt-15
    if site['pick_internal'] is True:
        p3_settings['PRIMER_PICK_INTERNAL_OLIGO'] = 1

    #  mode: target or include
    if type=='SEQUENCE_TARGET' or type=='SEQUENCE_INCLUDED_REGION':
        primer3.bindings.setP3Globals(p3_settings)
        if len(site['junction'])==0:
            p3 = primer3.bindings.designPrimers(
                {
                    'SEQUENCE_ID': id,
                    'SEQUENCE_TEMPLATE': template,
                    type: [pos,length]
                })
            p3['SEQUENCE_RELATIVE_TARGET_START'] = pos
            return [{'id': id, 'primers': p3}]
        else:
            # qPCR: primer pairs separated by at least one intron on the corresponding genomic DNA
            p3_A = primer3.bindings.designPrimers(
                {
                    'SEQUENCE_ID': id+'-PARTA',
                    'SEQUENCE_TEMPLATE': template,
                    'SEQUENCE_INCLUDED_REGION': [pos,length],
                    'SEQUENCE_TARGET': [[x,1] for x in site['junction']]
                })
            p3_A['SEQUENCE_RELATIVE_TARGET_START'] = pos

            # qPCR: primers must span an exon-exon junction
            p3_B = primer3.bindings.designPrimers(
                {
                    'SEQUENCE_ID': id+'-PARTB',
                    'SEQUENCE_TEMPLATE': template,
                    'SEQUENCE_INCLUDED_REGION': [pos,length],
                    'SEQUENCE_OVERLAP_JUNCTION_LIST': site['junction']
                })
            p3_B['SEQUENCE_RELATIVE_TARGET_START'] = pos
        
            return([{'id': id+'-PARTA', 'primers': p3_A},
                    {'id': id+'-PARTB', 'primers': p3_B}
            ])

    #  mode: force
    if type=='FORCE_END':
        p3_settings['PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE'] = -1
        p3_settings['PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE'] = 3
        primer3.bindings.setP3Globals(p3_settings)
        p3_L = primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': id+'-PARTA',
                'SEQUENCE_TEMPLATE': template,
                'SEQUENCE_FORCE_LEFT_END': pos
            })
        p3_L['SEQUENCE_RELATIVE_TARGET_START'] = pos

        p3_settings['PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE'] = 3
        p3_settings['PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE'] = -1
        primer3.bindings.setP3Globals(p3_settings)
        p3_R = primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': id+'-PARTB',
                'SEQUENCE_TEMPLATE': template,
                'SEQUENCE_FORCE_RIGHT_END': pos
            })
        p3_R['SEQUENCE_RELATIVE_TARGET_START'] = pos

        return([{'id': id+'-PARTA', 'primers': p3_L},
                {'id': id+'-PARTB', 'primers': p3_R}
        ])


def multiple(sites, cpu=2, monitor=True):
    if global_var.stop_run is True:
        return {'error': 'Stop running'}
    
    # distribute
    pool = mp.Pool(processes=cpu)
    multi_res = []
    global_var.all_tasks_num = 0
    global_var.current_task = 'design'
    for site in sites:
        global_var.all_tasks_num += 1
        multi_res.append(pool.apply_async(single, (site,)))

    # monitor
    if monitor is True:
        widgets = ['Designing Primers: ', progressbar.Counter(), ' Finished', ' (', progressbar.Percentage(), ')', \
            progressbar.Bar(), progressbar.ETA()]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=global_var.all_tasks_num).start()

        while True:
            global_var.complete_count = sum([1 for x in multi_res if x.ready()])
            if global_var.complete_count == global_var.all_tasks_num:
                bar.finish()
                break
            bar.update(global_var.complete_count)
            time.sleep(0.1)

    # results
    primers = {}
    for result in multi_res:
        primers_in_each_site = result.get()
        for primer in primers_in_each_site:
            primers[primer['id']] = primer['primers']
        
    return primers



if __name__ == "__main__":
    global_var.init()
    template = ''
    with open('tests/example.fa') as f:
        data = f.read().split('>')
        seq1 = data[1].splitlines()
        del seq1[0]
        template = ''.join([x.strip() for x in seq1])
        print(json.dumps(multiple([
                {
                    'id':'A-SEQUENCE_INCLUDED_REGION-1-600', 
                    'template':template,
                    'type':'SEQUENCE_INCLUDED_REGION',
                    'pos':1,
                    'length':600,
                    'size_min':75,
                    'size_max':1000,
                    'primer_num_return':10,
                    'pick_internal': False,
                    'junction': []
                }, 
                {
                    'id':'B-SEQUENCE_TARGET-200-10', 
                    'template':template,
                    'type':'SEQUENCE_TARGET',
                    'pos':200,
                    'length':10,
                    'size_min':75,
                    'size_max':1000,
                    'primer_num_return':10,
                    'pick_internal': False,
                    'junction': []
                },
                {
                    'id':'C-FORCE_END-319-1', 
                    'template':template,
                    'type':'FORCE_END',
                    'pos':319,
                    'length':1,
                    'size_min':75,
                    'size_max':1000,
                    'primer_num_return':10,
                    'pick_internal': False,
                    'junction': []
                }
                ]), indent=4))
