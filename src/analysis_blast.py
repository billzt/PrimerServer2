import re
import json

def filter_len(blast_out, len_min, len_max, region_file):
    '''
        return hit-pairs; make a region file
    '''

    # collect hits data
    hits = {}
    for line in blast_out.splitlines():
        (qseqid, qlen, qstart, qend, sseqid, slen, sstart, send, sstrand) = line.strip().split('\t')
        (primer_id, primer_rank, primer_sub) = re.match(r'(.*)\.(\d+)\.Primer(\d+)', qseqid).group(1,2,3)
        if primer_id not in hits:
            hits[primer_id] = {}
        if primer_rank not in hits[primer_id]:
            hits[primer_id][primer_rank] = []
        hits[primer_id][primer_rank].append({
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
    hits_pair = {}
    for primer_id in hits.keys():
        for primer_rank in hits[primer_id].keys():
            primer_hits = hits[primer_id][primer_rank]
            if len(primer_hits)==1:
                continue

            primer_hits = sorted(primer_hits, key=lambda i: i['sstart'])
            for i in range(0, len(primer_hits)-1):
                if primer_hits[i]['sstrand'] != 'plus':
                    continue
                for j in range(i+1, len(primer_hits)-1):
                    if primer_hits[j]['sstrand'] != 'minus':
                        continue
                    size = primer_hits[j]['sstart']-primer_hits[i]['sstart']
                    if size>len_min and size<len_max:
                        # extend this pair:
                        # plus
                        start_plus = primer_hits[i]['sstart']-(primer_hits[i]['qstart']-1)
                        if start_plus<1:    # extension failed
                            continue
                        end_plus = primer_hits[i]['send']+(primer_hits[i]['qlen']-primer_hits[i]['qend'])
                        if end_plus>primer_hits[i]['slen']: # extension failed
                            continue
                        
                        # minus
                        start_minus = primer_hits[j]['sstart']+(primer_hits[j]['qlen']-primer_hits[j]['qend'])
                        if start_minus>primer_hits[j]['slen']: # extension failed
                            continue
                        end_minus = primer_hits[j]['send']-(primer_hits[i]['qstart']-1)
                        if end_minus<1:  # extension failed
                            continue

                        # store this pair
                        if primer_id not in hits_pair:
                            hits_pair[primer_id] = {}
                        if primer_rank not in hits_pair[primer_id]:
                            hits_pair[primer_id][primer_rank] = []
                        hits_pair[primer_id][primer_rank].append({
                            'plus': {
                                'qseqid': primer_hits[i]['qseqid'],
                                'sseqid': primer_hits[i]['sseqid'],
                                'sstart': start_plus,
                                'send': end_plus
                            },
                            'minus': {
                                'qseqid': primer_hits[j]['qseqid'],
                                'sseqid': primer_hits[j]['sseqid'],
                                'sstart': end_minus,    # flip
                                'send': start_minus
                            }
                        })

                        # store regions
                        hits_regions[primer_hits[i]['sseqid']+':'+start_plus+'-'+end_plus] = 1
                        hits_regions[primer_hits[j]['sseqid']+':'+start_minus+'-'+end_minus] = 1
    
    # make regions file
    with open(region_file, 'w') as f:


    return hits_pair



if __name__ == "__main__":
    blast_out = open('tests/query_blast.fa.out').read()
    print(json.dumps(filter_len(blast_out=blast_out, len_min=75, len_max=1000)))
        