import json

def text(primers, dbs):
    print('#Site_ID', 'Primer_Rank', 'Primer_Seq_Left', 'Primer_Seq_Right', 'Target_Amplicon_Size', \
        'Primer_Pair_Penalty_Score', 'Database', 'Possible_Amplicon_Number', \
            'Primer_Rank_in_Primer3_output', sep='\t')
    for db in dbs:
        for (id, primers) in primers.items():
            if 'PRIMER_PAIR_NUM_RETURNED_FINAL' not in primers:
                print(id, 'No_Primer\t', sep='\t', end='')
                if 'PRIMER_LEFT_EXPLAIN' in primers:
                    print(primers['PRIMER_LEFT_EXPLAIN'].replace(' ', '_'), \
                        primers['PRIMER_RIGHT_EXPLAIN'].replace(' ', '_'), \
                            primers['PRIMER_PAIR_EXPLAIN'].replace(' ', '_'), \
                                sep='\t')
                print('###')
                continue

            for amplicon_rank in range(0, primers['PRIMER_PAIR_NUM_RETURNED_FINAL']):
                print_data = []
                print_data.append(id)
                print_data.append(str(amplicon_rank))
                raw_rank = primers[f'PRIMER_PAIR_AMPLICON_NUM_RANK_{amplicon_rank}']
                print_data.append(primers[f'PRIMER_LEFT_{raw_rank}_SEQUENCE'])
                print_data.append(primers[f'PRIMER_RIGHT_{raw_rank}_SEQUENCE'])
                print_data.append(str(primers[f'PRIMER_PAIR_{raw_rank}_PRODUCT_SIZE']))
                print_data.append(str(primers[f'PRIMER_PAIR_{raw_rank}_PENALTY']))
                print_data.append(db)
                print_data.append(str(len(primers[db][f'PRIMER_PAIR_{raw_rank}_AMPLICONS'])))
                print_data.append(str(raw_rank))
                print('\t'.join(print_data))
            print('###')

if __name__ == "__main__":
    primers = json.load(open('tests/sort_primers.json'))
    dbs = ['example.fa']
    text(primers, dbs)
