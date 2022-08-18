import re
import json
import sys
from itertools import product

import primer3

from primerserver2.core import global_var
from primerserver2.core.Santalucia_NN_Tm import transform_degenerate

def calculate_GC(seq):
    gc = len([x for x in seq.upper() if x=='G' or x=='C'])
    return gc/len(seq)*100

def check_primer_seq(seq):
    if len(seq)<10 or len(seq)>40:
        return False
    if re.search('[^ATGCNRYMKSWHBIVD]', seq, re.RegexFlag.IGNORECASE) is not None:
        return False

d = {
'A': ['A'],
'C': ['C'],
'G': ['G'],
'T': ['T'],
'R': ['A', 'G'],
'Y': ['C', 'T'],
'S': ['G', 'C'],
'W': ['A', 'T'],
'K': ['G', 'T'],
'M': ['A', 'C'],
'B': ['C', 'G', 'T'],
'D': ['A', 'G', 'T'],
'H': ['A', 'C', 'T'],
'V': ['A', 'C', 'G'],
'N': ['A', 'C', 'G', 'T'],
'I': ['A', 'C', 'G', 'T'],}

def expand_degenerate_bases(seq):
   """return a list of all possible k-mers given a degenerate base"""
   return list(map("".join, product(*map(d.get, seq))))


def make_primers(query):
    '''
        In case of user's input, make list of primer dicts similar as the result of design_primer module
        Input:
            query: a string in multi-lines
    '''
    if global_var.stop_run is True:
        return {'error': 'Stop running'}

    primers = {}

    for line in query.splitlines():
        if line.strip() == '':
            continue
        line_data = re.split(r'\s+', line.strip())
        if len(line_data)==4:
            (id, rank, seq_F, seq_R) = line_data
            rank = int(rank)
        elif len(line_data)==3:
            (id, seq_F, seq_R) = line_data
            rank = 0
        else:
            return {'error': f'Your input: {line} does not have three or four columns'}
        
        # whether seq_F and seq_R are valid primer seqs
        if check_primer_seq(seq_F) is False:
            return {'error': f'Your input: {seq_F} does not seem like a vaild primer seq'}
        if check_primer_seq(seq_R) is False:
            return {'error': f'Your input: {seq_R} does not seem like a vaild primer seq'}
        
        # generate primers
        seq_F = seq_F.upper()
        seq_R = seq_R.upper()
        seq_F_expands = expand_degenerate_bases(seq_F)
        seq_R_expands = expand_degenerate_bases(seq_R)
        id_sub_num = 0
        for seq_F_expand in seq_F_expands:
            for seq_R_expand in seq_R_expands:
                id_sub_num += 1
                if f'{id}|||{id_sub_num}' not in primers:
                    primers[f'{id}|||{id_sub_num}'] = {}
                primers[f'{id}|||{id_sub_num}']['PRIMER_PAIR_NUM_RETURNED'] = rank+1
                primers[f'{id}|||{id_sub_num}'][f'PRIMER_PAIR_{rank}_PENALTY'] = 0
                primers[f'{id}|||{id_sub_num}'][f'PRIMER_LEFT_{rank}_SEQUENCE'] = seq_F_expand
                primers[f'{id}|||{id_sub_num}'][f'PRIMER_RIGHT_{rank}_SEQUENCE'] =  seq_R_expand
                primers[f'{id}|||{id_sub_num}'][f'PRIMER_LEFT_{rank}'] = [-1, len(seq_F_expand)]
                primers[f'{id}|||{id_sub_num}'][f'PRIMER_RIGHT_{rank}'] = [-1, len(seq_R_expand)]
                primers[f'{id}|||{id_sub_num}'][f'PRIMER_LEFT_{rank}_TM'] = primer3.calcTm(seq_F_expand)
                primers[f'{id}|||{id_sub_num}'][f'PRIMER_RIGHT_{rank}_TM'] = primer3.calcTm(seq_R_expand)
                primers[f'{id}|||{id_sub_num}'][f'PRIMER_LEFT_{rank}_GC_PERCENT'] = calculate_GC(seq_F_expand)
                primers[f'{id}|||{id_sub_num}'][f'PRIMER_RIGHT_{rank}_GC_PERCENT'] = calculate_GC(seq_R_expand)
                primers[f'{id}|||{id_sub_num}'][f'PRIMER_PAIR_{rank}_PRODUCT_SIZE'] = -1
                primers[f'{id}|||{id_sub_num}']['SEQUENCE_RELATIVE_TARGET_START'] = 0
    
    return primers

if __name__ == "__main__":
    global_var.init()
    with open('tests/_internal_/query_check') as f:
        primers = make_primers(f.read())
        print(json.dumps(primers, indent=4))
