import re
import json

import primer3

from primerserver2.core import global_var
from primerserver2.core.Santalucia_NN_Tm import transform_degenerate

def calculate_GC(seq):
    gc = len([x for x in seq.upper() if x=='G' or x=='C'])
    return gc/len(seq)*100

def check_primer_seq(seq):
    if len(seq)<10 or len(seq)>40:
        return False
    if re.search('[^ATGCNRYMKSWHBVD]', seq, re.RegexFlag.IGNORECASE) is not None:
        return False

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
        if id not in primers:
            primers[id] = {}
        primers[id]['PRIMER_PAIR_NUM_RETURNED'] = rank+1
        primers[id][f'PRIMER_PAIR_{rank}_PENALTY'] = 0
        primers[id][f'PRIMER_LEFT_{rank}_SEQUENCE'] = seq_F
        primers[id][f'PRIMER_RIGHT_{rank}_SEQUENCE'] =  seq_R
        primers[id][f'PRIMER_LEFT_{rank}'] = [-1, len(seq_F)]
        primers[id][f'PRIMER_RIGHT_{rank}'] = [-1, len(seq_R)]
        primers[id][f'PRIMER_LEFT_{rank}_TM'] = primer3.calcTm(transform_degenerate(seq_F))
        primers[id][f'PRIMER_RIGHT_{rank}_TM'] = primer3.calcTm(transform_degenerate(seq_R))
        primers[id][f'PRIMER_LEFT_{rank}_GC_PERCENT'] = calculate_GC(transform_degenerate(seq_F))
        primers[id][f'PRIMER_RIGHT_{rank}_GC_PERCENT'] = calculate_GC(transform_degenerate(seq_R))
        primers[id][f'PRIMER_PAIR_{rank}_PRODUCT_SIZE'] = -1
        primers[id]['SEQUENCE_RELATIVE_TARGET_START'] = 0
    
    return primers

if __name__ == "__main__":
    global_var.init()
    with open('tests/_internal_/query_check') as f:
        primers = make_primers(f.read())
        print(json.dumps(primers, indent=4))
