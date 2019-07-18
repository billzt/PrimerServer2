
import re
import os
import json

import pysam

def faidx(template_file, region_file):
    '''
        return a list of sequence dict using samtools faidx based on sequence file (fa)
        and region file
    '''
    if os.path.isfile(template_file)==False:
        raise Exception(f'File {template_file} is not found')

    seqs = pysam.faidx('-r', region_file, template_file).split('>')
    del seqs[0]
    result_seqs = {}
    for seq in seqs:
        seq_split = seq.splitlines()
        result_seqs[seq_split[0]] = ''.join(seq_split[1:])
    return result_seqs

def build(query_file, template_file, region_file, primer_type):
    '''
        return a list of primer sites dict that can be passed to the design_primer module
    '''
    primer_sites = []
    retrieve_region2raw_region = {}
    with open(region_file, 'w') as fw:
        with open(query_file) as f:
            for line in f:
                query_data = re.split(r'\s+', line.strip())
                (chr, pos, length, size_min, size_max) = (query_data[0], 1, 1, 70, 1000)
                if len(query_data)>1 :     # seq ID and pos
                    pos = int(query_data[1].replace(',', ''))
                if len(query_data)>2 :
                    length = int(query_data[2])
                if len(query_data)>3 :
                    size_min = int(query_data[3])
                    size_max = int(query_data[4])

                # retrieve_start and retrieve_end are used by samtools to extract template sequences
                retrieve_start = max(pos-size_max, 1)
                retrieve_end = pos+length+size_max
                print(f'{chr}:{retrieve_start}-{retrieve_end}', file=fw)
                retrieve_region2raw_region[f'{chr}:{retrieve_start}-{retrieve_end}'] = [chr, pos, length, size_min, size_max, retrieve_start]
    
    result_seqs = faidx(template_file, region_file)
    for (result_seq_id, result_seq) in result_seqs.items():
        (chr, pos, length, size_min, size_max, retrieve_start) = retrieve_region2raw_region[result_seq_id]
        primer_sites.append({
            'id': chr+'-'+str(pos)+'-'+str(length), 
            'template': result_seq,
            'type': primer_type,
            'pos': pos-retrieve_start,
            'length':length,
            'size_min':size_min,
            'size_max':size_max
        })
    
    return primer_sites

if __name__ == "__main__":
    primer_sites = build(query_file='tests/query_design', template_file='tests/example.fa', 
        region_file='tests/region.tmp', primer_type='SEQUENCE_TARGET')
    print(json.dumps(primer_sites))
