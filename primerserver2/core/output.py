import json
import sys
import os

def tsv(primers_dict, dbs, file):
    print('#Site_ID', 'Primer_Rank', 'Primer_Seq_Left', 'Primer_Seq_Right', 'Target_Amplicon_Size', \
        'Primer_Pair_Penalty_Score', 'Database', 'Possible_Amplicon_Number', \
            'Primer_Rank_in_Primer3_output', 'Tm_left', 'Tm_Right', sep='\t', file=file)
    for db in dbs:
        for (id, primers) in primers_dict.items():
            # No primers
            if primers['PRIMER_PAIR_NUM_RETURNED']==0:
                print(id, 'No_Primer\t', sep='\t', end='', file=file)
                if 'PRIMER_LEFT_EXPLAIN' in primers:
                    print(primers['PRIMER_LEFT_EXPLAIN'].replace(' ', '_'), \
                        primers['PRIMER_RIGHT_EXPLAIN'].replace(' ', '_'), \
                            primers['PRIMER_PAIR_EXPLAIN'].replace(' ', '_'), \
                                sep='\t', file=file)
                print('###', file=file)
                continue

            if 'PRIMER_PAIR_NUM_RETURNED_FINAL' in primers:
                primer_num = primers['PRIMER_PAIR_NUM_RETURNED_FINAL']
            else:
                primer_num = primers['PRIMER_PAIR_NUM_RETURNED']
            for amplicon_rank in range(0, primer_num):
                db_base = os.path.basename(db)
                print_data = []
                print_data.append(id)
                print_data.append(str(amplicon_rank))
                if 'PRIMER_PAIR_NUM_RETURNED_FINAL' in primers:
                    raw_rank = primers[f'PRIMER_PAIR_AMPLICON_NUM_RANK_{amplicon_rank}']
                else:
                    raw_rank = amplicon_rank
                print_data.append(primers[f'PRIMER_LEFT_{raw_rank}_SEQUENCE'])
                print_data.append(primers[f'PRIMER_RIGHT_{raw_rank}_SEQUENCE'])
                print_data.append(str(primers[f'PRIMER_PAIR_{raw_rank}_PRODUCT_SIZE']))
                print_data.append(str(round(primers[f'PRIMER_PAIR_{raw_rank}_PENALTY'],2)))
                print_data.append(db_base)
                if db_base in primers:
                    print_data.append(str(len(primers[db_base][f'PRIMER_PAIR_{raw_rank}_AMPLICONS'])))
                else:
                    print_data.append('NA')
                print_data.append(str(raw_rank))
                print_data.append(str(round(primers[f'PRIMER_LEFT_{raw_rank}_TM'], 2)))
                print_data.append(str(round(primers[f'PRIMER_RIGHT_{raw_rank}_TM'], 2)))
                print('\t'.join(print_data), file=file)
            print('###', file=file)


if __name__ == "__main__":
    primers = json.load(open('tests/_internal_/sort_primers.json'))
    dbs = ['example.fa']
    tsv(primers, dbs, file=sys.stdout)
