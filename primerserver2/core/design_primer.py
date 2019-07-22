
import json
import os
import multiprocessing as mp

import primer3

p3_settings = dict(json.load(open('primerserver2/data/p3_settings.json')))

def single(site):
    '''
        Using the primer3-module to design primers for a single site, return a list (one or two items)
        in DICT format
        site: site['id'], site['template'], site['type'], site['pos'], site['length'], 
        site['size_max'], site['size_min']
    '''
    id = site['id']
    template = site['template']
    type = site['type']
    pos = site['pos']
    length = site['length']
    size_max = site['size_max']
    size_min = site['size_min']
    p3_settings['PRIMER_PRODUCT_SIZE_RANGE'] = [[size_min, size_max]]

    if type=='SEQUENCE_TARGET' or type=='SEQUENCE_INCLUDED_REGION':
        primer3.bindings.setP3Globals(p3_settings)
        p3 = primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': id,
                'SEQUENCE_TEMPLATE': template,
                type: [pos,length]
            })
        return([{'id': id, 'primers': p3}])
    if type=='FORCE_END':
        p3_settings['PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE'] = -1
        p3_settings['PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE'] = 3
        primer3.bindings.setP3Globals(p3_settings)
        p3_L = primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': id+'-LEFT',
                'SEQUENCE_TEMPLATE': template,
                'SEQUENCE_FORCE_LEFT_END': pos
            })

        p3_settings['PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE'] = 3
        p3_settings['PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE'] = -1
        primer3.bindings.setP3Globals(p3_settings)
        p3_R = primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': id+'-RIGHT',
                'SEQUENCE_TEMPLATE': template,
                'SEQUENCE_FORCE_RIGHT_END': pos
            })

        return([{'id': id+'-LEFT', 'primers': p3_L},
                {'id': id+'-RIGHT', 'primers': p3_R}
        ])


def multiple(sites, cpu=2):
    pool = mp.Pool(processes=cpu)
    multi_res = []
    for site in sites:
        multi_res.append(pool.apply_async(single, (site,)))

    primers = {}
    for result in multi_res:
        primers_in_each_site = result.get()
        for primer in primers_in_each_site:
            primers[primer['id']] = primer['primers']
        
    return primers



if __name__ == "__main__":
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
                'size_max':1000
            }, 
            {
                'id':'B-SEQUENCE_TARGET-200-10', 
                'template':template,
                'type':'SEQUENCE_TARGET',
                'pos':200,
                'length':10,
                'size_min':75,
                'size_max':1000
            },
            {
                'id':'C-FORCE_END-319-1', 
                'template':template,
                'type':'FORCE_END',
                'pos':319,
                'length':1,
                'size_min':75,
                'size_max':1000
            }
            ]), indent=4))
