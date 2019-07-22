import re
import json

import primer3

def calculate_GC(seq):
    gc = len([x for x in seq.upper() if x=='G' or x=='C'])
    return gc/len(seq)*100

def make_primers(query):
    '''
        In case of user's input, make list of primer dicts similar as the result of design_primer module
        Input:
            query: a string in multi-lines
    '''
    primers = {}

    for line in query.splitlines():
        line_data = re.split(r'\s+', line.strip())
        if len(line_data)==4:
            (id, rank, seq_F, seq_R) = line_data
            rank = int(rank)
        else:
            (id, seq_F, seq_R) = line_data
            rank = 0
        
        if id not in primers:
            primers[id] = {}
        primers[id]['PRIMER_PAIR_NUM_RETURNED'] = rank+1
        primers[id][f'PRIMER_PAIR_{rank}_PENALTY'] = 0
        primers[id][f'PRIMER_LEFT_{rank}_SEQUENCE'] = seq_F
        primers[id][f'PRIMER_RIGHT_{rank}_SEQUENCE'] =  seq_R
        primers[id][f'PRIMER_LEFT_{rank}'] = [-1, len(seq_F)]
        primers[id][f'PRIMER_RIGHT_{rank}'] = [-1, len(seq_R)]
        primers[id][f'PRIMER_LEFT_{rank}_TM'] = primer3.calcTm(seq_F)
        primers[id][f'PRIMER_RIGHT_{rank}_TM'] = primer3.calcTm(seq_R)
        primers[id][f'PRIMER_LEFT_{rank}_GC_PERCENT'] = calculate_GC(seq_F)
        primers[id][f'PRIMER_RIGHT_{rank}_GC_PERCENT'] = calculate_GC(seq_R)
        primers[id][f'PRIMER_PAIR_{rank}_PRODUCT_SIZE'] = -1
    
    return primers

if __name__ == "__main__":
    with open('tests/_internal_/query_check') as f:
        primers = make_primers(f.read())
        print(json.dumps(primers, indent=4))
