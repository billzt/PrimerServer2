import json
import sys
import os

def tsv(primers_dict, dbs):
    print_line = []
    header_line = '\t'.join(['#Site_ID', 'Primer_Rank', 'Primer_Seq_Left', 'Primer_Seq_Right', 'Primer_Seq_Oligo', 'Target_Amplicon_Size', \
        'Primer_Pair_Penalty_Score', 'Primer_Rank_in_Primer3_output', 'Tm_left', 'Tm_Right', 'Tm_Oligo', \
            'Database\tPossible_Amplicon_Number\tAmplicon_Size_in_DB\tAmplicon_Region_in_DB\t'*len(dbs)])
    print_line.append(header_line)
    for (id, primers) in primers_dict.items():
        template_id = id.split('-')[0]

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
            if f'PRIMER_INTERNAL_{raw_rank}_SEQUENCE' in primers:
                print_data.append(primers[f'PRIMER_INTERNAL_{raw_rank}_SEQUENCE'])
            else:
                print_data.append('NA')
            print_data.append(str(primers[f'PRIMER_PAIR_{raw_rank}_PRODUCT_SIZE']))
            print_data.append(str(round(primers[f'PRIMER_PAIR_{raw_rank}_PENALTY'],2)))
            print_data.append(str(raw_rank))
            print_data.append(str(round(primers[f'PRIMER_LEFT_{raw_rank}_TM'], 2)))
            print_data.append(str(round(primers[f'PRIMER_RIGHT_{raw_rank}_TM'], 2)))
            if f'PRIMER_INTERNAL_{raw_rank}_TM' in primers:
                print_data.append(str(round(primers[f'PRIMER_INTERNAL_{raw_rank}_TM'], 2)))
            else:
                print_data.append('NA')
            for db in dbs:
                db_base = os.path.basename(db)
                print_data.append(db_base)
                if db_base in primers:
                    if f'PRIMER_PAIR_{raw_rank}_AMPLICONS' in primers[db_base]:
                        amplicon_num_in_db = 0
                        for amplicon in primers[db_base][f'PRIMER_PAIR_{raw_rank}_AMPLICONS']:
                            if 'isoform' not in amplicon:
                                amplicon_num_in_db+=1
                            elif amplicon['isoform'] is False:
                                amplicon_num_in_db+=1
                            elif amplicon['plus']['sseqid'] == template_id:
                                amplicon_num_in_db+=1

                        print_data.append(str(amplicon_num_in_db))
                        if amplicon_num_in_db==1:
                            print_data.append(str(primers[db_base][f'PRIMER_PAIR_{raw_rank}_AMPLICONS'][0]['product_size']))
                            print_data.append(str(primers[db_base][f'PRIMER_PAIR_{raw_rank}_AMPLICONS'][0]['region']))
                        else:
                            print_data.extend(('NA', 'NA'))
                    else:
                        print_data.extend(('NA', 'NA', 'NA'))
                else:
                    print_data.extend(('NA', 'NA', 'NA'))
            

            print_line.append('\t'.join(print_data))
        print_line.append('###')

    return '\n'.join(print_line)

def dimer_list(dimers):
    print_line = ['################### Dimers #################']
    if 'error' in dimers:
        print_line.append('ERROR: '+dimers['error'])
    else:
        header_line = '\t'.join(['#Site_1', 'Primer_Rank_1', 'Seq_1', 'Site_2', 'Primer_Rank_2', 'Seq_2', 'Tm'])
        print_line.append(header_line)
        for site_1 in dimers.keys():
            for primer_1 in dimers[site_1].keys():
                for site_2 in dimers[site_1][primer_1].keys():
                    for primer_2 in dimers[site_1][primer_1][site_2].keys():
                        dimer = dimers[site_1][primer_1][site_2][primer_2]
                        print_line.append('\t'.join([site_1, str(primer_1), dimer['seq_1'], \
                            site_2, str(primer_2), dimer['seq_2'], str(dimer['Tm'])]))
            print_line.append('###')
    return '\n'.join(print_line)


if __name__ == "__main__":
    primers = json.load(open('tests/_internal_/sort_primers2.json'))
    dbs = ["Ghir.JGI.genomic", "Gbar.IGDB.genomic"]
    print(tsv(primers, dbs))
