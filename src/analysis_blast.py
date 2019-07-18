import re
import json
import sys

from make_sites import faidx
from Santalucia_NN_Tm import complement, rev_complement, NN_Tm
from make_primers import make_primers

def filter_len(blast_out, len_min, len_max, region_file):
    '''
        return hit-pairs; make a region file
    '''

    # collect hits data
    hits = {}
    for line in blast_out.splitlines():
        (qseqid, qlen, qstart, qend, sseqid, slen, sstart, send, sstrand) = line.strip().split('\t')
        (primer_id, primer_rank, primer_sub) = re.match(r'(.*)\.(\d+)\.Primer\.(\w+)$', qseqid).group(1,2,3)
        if primer_id not in hits:
            hits[primer_id] = {}
        if primer_rank not in hits[primer_id]:
            hits[primer_id][primer_rank] = []
        hits[primer_id][primer_rank].append({
            'sub': primer_sub,  # LEFT or RIGHT
            'qseqid': qseqid,
            'qlen': int(qlen),
            'qstart': int(qstart),
            'qend': int(qend),
            'sseqid': sseqid,
            'slen': int(slen),
            'sstart': int(sstart),
            'send': int(send),
            'sstrand': sstrand
        })
    
    # filter length between plus and minus
    amplicons = {}
    hits_regions = {}
    for primer_id in hits.keys():
        for primer_rank in hits[primer_id].keys():
            primer_hits = hits[primer_id][primer_rank]
            if len(primer_hits)==1:
                continue

            primer_hits = sorted(primer_hits, key=lambda i: i['sstart'])
            for i in range(0, len(primer_hits)):
                if primer_hits[i]['sstrand'] != 'plus':
                    continue
                for j in range(i+1, len(primer_hits)):
                    if primer_hits[j]['sstrand'] != 'minus':
                        continue
                    size = primer_hits[j]['sstart']-primer_hits[i]['sstart']
                    if size>len_min and size<len_max and primer_hits[i]['sseqid']==primer_hits[j]['sseqid']:
                        # extend this pair:
                        # plus
                        start_plus = primer_hits[i]['sstart']-(primer_hits[i]['qstart']-1)
                        if start_plus<1:    # extension failed
                            continue
                        end_plus = primer_hits[i]['send']+(primer_hits[i]['qlen']-primer_hits[i]['qend'])
                        if end_plus>primer_hits[i]['slen']: # extension failed
                            continue
                        
                        # minus
                        start_minus = primer_hits[j]['sstart']+(primer_hits[j]['qstart']-1)
                        if start_minus>primer_hits[j]['slen']: # extension failed
                            continue
                        end_minus = primer_hits[j]['send']-(primer_hits[j]['qlen']-primer_hits[j]['qend'])
                        if end_minus<1:  # extension failed
                            continue

                        # store this pair
                        if primer_id not in amplicons:
                            amplicons[primer_id] = {}
                        if primer_rank not in amplicons[primer_id]:
                            amplicons[primer_id][primer_rank] = []
                        amplicons[primer_id][primer_rank].append({
                            'plus': {
                                'qseqid': primer_hits[i]['qseqid'],
                                'sub': primer_hits[i]['sub'],
                                'sseqid': primer_hits[i]['sseqid'],
                                'sstart': start_plus,
                                'send': end_plus
                            },
                            'minus': {
                                'qseqid': primer_hits[j]['qseqid'],
                                'sub': primer_hits[j]['sub'],
                                'sseqid': primer_hits[j]['sseqid'],
                                'sstart': end_minus,    # flip
                                'send': start_minus
                            }
                        })

                        # store regions
                        hits_regions[primer_hits[i]['sseqid']+':'+str(start_plus)+'-'+str(end_plus)] = 1
                        hits_regions[primer_hits[j]['sseqid']+':'+str(end_minus)+'-'+str(start_minus)] = 1 # flip
    
    # make regions file
    with open(region_file, 'w') as f:
        for region in hits_regions.keys():
            print(region, file=f)

    return amplicons

def filter_Tm(amplicons, query_primers, hits_seqs, Tm_diff=20, max_amplicons=10):
    for primer_id in amplicons.keys():
        # calculate Tm for each primer within this site
        rank_to_amplicon_valid = {}
        for primer_rank in amplicons[primer_id].keys():
            # for each primer, there are multiple amplicons
            for i in range(0, len(amplicons[primer_id][primer_rank])):
                pair = amplicons[primer_id][primer_rank][i]

                # hit seqs (template)
                region_plus = pair['plus']['sseqid']+':'+str(pair['plus']['sstart'])+'-'+str(pair['plus']['send'])
                region_minus = pair['minus']['sseqid']+':'+str(pair['minus']['sstart'])+'-'+str(pair['minus']['send'])
                hit_seq_plus = hits_seqs[region_plus]
                hit_seq_minus = rev_complement(hits_seqs[region_minus])  # samtools faidx hasn't do rev_complement yet

                # query primers (query)
                sub_plus = pair['plus']['sub']  # LEFT or RIGHT
                sub_minus = pair['minus']['sub']  # LEFT or RIGHT
                query_primer_plus = query_primers[primer_id][f'PRIMER_{sub_plus}_{primer_rank}_SEQUENCE']
                query_primer_minus = query_primers[primer_id][f'PRIMER_{sub_minus}_{primer_rank}_SEQUENCE']

                # Tm
                if len(query_primer_plus) != len(hit_seq_plus):
                    raise Exception(f'Not the same length in {primer_id}-{primer_rank}-{sub_plus}: \n \
                        Primer: {query_primer_plus} \n \
                        Hit: {region_plus} {hit_seq_plus}')
                Tm_plus = float(NN_Tm(seq=query_primer_plus, compl_seq=complement(hit_seq_plus)))
                amplicons[primer_id][primer_rank][i]['plus']['Tm'] = Tm_plus

                if len(query_primer_minus) != len(hit_seq_minus):
                    raise Exception(f'Not the same length in {primer_id}-{primer_rank}-{sub_minus}: \n \
                        Primer: {query_primer_minus} \n \
                        Hit: {region_minus} {hit_seq_minus}')
                Tm_minus = float(NN_Tm(seq=query_primer_minus, compl_seq=complement(hit_seq_minus)))
                amplicons[primer_id][primer_rank][i]['minus']['Tm'] = Tm_minus

                amplicons[primer_id][primer_rank][i]['Tm'] = min(Tm_plus, Tm_minus)
            
            # sort amplicons by Tm (reverse) to get perfect Tm
            amplicons_this_primer = sorted(amplicons[primer_id][primer_rank], key=lambda item: item['Tm'], reverse=True)
            Tm_perfect = amplicons_this_primer[0]['Tm']

            # add amplicons to the primer dict
            amplicon_num_filter = 0
            for amplicon in amplicons_this_primer:
                if amplicon['Tm']<Tm_perfect-Tm_diff:  # ignore unstable amplicons
                    continue
                if f'PRIMER_PAIR_{primer_rank}_AMPLICONS' not in query_primers[primer_id]:
                    query_primers[primer_id][f'PRIMER_PAIR_{primer_rank}_AMPLICONS'] = []
                query_primers[primer_id][f'PRIMER_PAIR_{primer_rank}_AMPLICONS'].append(amplicon)
                amplicon_num_filter += 1
                if amplicon_num_filter==max_amplicons:
                    break

            query_primers[primer_id][f'PRIMER_PAIR_{primer_rank}_AMPLICONS_NUM'] = amplicon_num_filter
            rank_to_amplicon_valid[primer_rank] = amplicon_num_filter
        
        # sort primers within this site
        for (i, primer_rank) in enumerate(sorted(rank_to_amplicon_valid, key=rank_to_amplicon_valid.__getitem__)):
            query_primers[primer_id][f'PRIMER_PAIR_TM_RANK_{i}'] = primer_rank


    return query_primers



if __name__ == "__main__":
    blast_out = open('tests/query_blast.fa.out').read()
    amplicons = filter_len(blast_out=blast_out, len_min=75, len_max=1000, region_file='tests/query_blast.fa.out.regions')
    hits_seqs = faidx(template_file='tests/example.fa', region_file='tests/query_blast.fa.out.regions')
    report_primers = filter_Tm(amplicons, query_primers=make_primers('tests/query_check'), hits_seqs=hits_seqs)
    print(json.dumps(report_primers))
