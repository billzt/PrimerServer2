import json
import sys
import os

def tsv(primers_dict, dbs):
    print_line = []
    header_line = '\t'.join(['#Site_ID', 'Primer_Rank', 'Primer_Seq_Left', 'Primer_Seq_Right', 'Target_Amplicon_Size', \
        'Primer_Pair_Penalty_Score', 'Primer_Rank_in_Primer3_output', 'Tm_left', 'Tm_Right', \
            'Database\tPossible_Amplicon_Number\tAmplicon_Size_in_DB\t'*len(dbs)])
    print_line.append(header_line)
    for (id, primers) in primers_dict.items():
        # No primers
        if primers['PRIMER_PAIR_NUM_RETURNED']==0:
            print_data = []
            print_data.append(id)
            print_data.append('No_Primer')
            if 'PRIMER_LEFT_EXPLAIN' in primers:
                print_data.append(primers['PRIMER_LEFT_EXPLAIN'].replace(' ', '_'))
                print_data.append(primers['PRIMER_RIGHT_EXPLAIN'].replace(' ', '_'))
                print_data.append(primers['PRIMER_PAIR_EXPLAIN'].replace(' ', '_'))
            print_line.append('\t'.join(print_data))
            print_line.append('###')
            continue

        if 'PRIMER_PAIR_NUM_RETURNED_FINAL' in primers:
            primer_num = primers['PRIMER_PAIR_NUM_RETURNED_FINAL']
        else:
            primer_num = primers['PRIMER_PAIR_NUM_RETURNED']
        for amplicon_rank in range(0, primer_num):
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
            print_data.append(str(raw_rank))
            print_data.append(str(round(primers[f'PRIMER_LEFT_{raw_rank}_TM'], 2)))
            print_data.append(str(round(primers[f'PRIMER_RIGHT_{raw_rank}_TM'], 2)))
            for db in dbs:
                db_base = os.path.basename(db)
                print_data.append(db_base)
                if db_base in primers:
                    if f'PRIMER_PAIR_{raw_rank}_AMPLICONS' in primers[db_base]:
                        amplicon_num_in_db = len(primers[db_base][f'PRIMER_PAIR_{raw_rank}_AMPLICONS'])
                        print_data.append(str(amplicon_num_in_db))
                        if amplicon_num_in_db==1:
                            print_data.append(str(primers[db_base][f'PRIMER_PAIR_{raw_rank}_AMPLICONS'][0]['product_size']))
                        else:
                            print_data.append('NA')
                    else:
                        print_data.append('NA')
                        print_data.append('NA')
                else:
                    print_data.append('NA')
                    print_data.append('NA')
            

            print_line.append('\t'.join(print_data))
        print_line.append('###')

    return '\n'.join(print_line)


if __name__ == "__main__":
    primers = json.load(open('tests/_internal_/sort_primers2.json'))
    dbs = ["Ghir.JGI.genomic", "Gbar.IGDB.genomic"]
    print(tsv(primers, dbs))
