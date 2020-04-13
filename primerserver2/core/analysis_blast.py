import re
import json
import sys

from primerserver2.core.make_sites import faidx
from primerserver2.core.Santalucia_NN_Tm import complement, rev_complement, NN_Tm, transform_degenerate
from primerserver2.core.make_primers import make_primers

def filter_len(blast_out, len_min, len_max, chr_main_only=True):    # 1s
    '''
        return 
            'amplicons': 
                primer_id=>
                    primer_rank=> [hit-pairs]
            'regions_primer': a region string for each primer used by samtools faidx;
    '''

    # collect hits data
    primer_hits_all_sites = {}
    for line in blast_out.splitlines():
        (qseqid, qlen, qstart, qend, sseqid, slen, sstart, send, sstrand) = line.strip().split('\t')
        
        # include main chrs only
        # if 'random' in sseqid or 'alt' in sseqid or 'Un' in sseqid or 'scaffold' in sseqid:
        #     continue

        # collect
        # qseqid: {id}.{rank}.L or {id}.{rank}.R
        qseqid_data = qseqid.split('.')
        primer_rank = qseqid_data[-2]
        primer_id = '.'.join(qseqid_data[0:-2])
        if primer_id not in primer_hits_all_sites:
            primer_hits_all_sites[primer_id] = {}
        if primer_rank not in primer_hits_all_sites[primer_id]:
            primer_hits_all_sites[primer_id][primer_rank] = {}
        if sseqid not in primer_hits_all_sites[primer_id][primer_rank]:
            primer_hits_all_sites[primer_id][primer_rank][sseqid] = []
        primer_hits_all_sites[primer_id][primer_rank][sseqid].append({
            'qseqid': qseqid,
            'qlen': int(qlen),
            'qstart': int(qstart),
            'qend': int(qend),
            'slen': int(slen),
            'sstart': int(sstart),
            'send': int(send),
            'sstrand': sstrand
        })
    
    # # filter length between plus and minus
    amplicons = {}
    hits_regions = {}   # make region file (FASTA)

    for primer_id in primer_hits_all_sites.keys():
        amplicons[primer_id] = {}
        for primer_rank in primer_hits_all_sites[primer_id].keys():
            if primer_rank not in amplicons[primer_id]:
                amplicons[primer_id][primer_rank] = []
            for sseqid in primer_hits_all_sites[primer_id][primer_rank].keys():
                primer_hits = sorted(primer_hits_all_sites[primer_id][primer_rank][sseqid], key=lambda i: i['sstart'])
                primer_hits_plus = [x for x in primer_hits if x['sstrand']=='plus']
                primer_hits_minus = [x for x in primer_hits if x['sstrand']=='minus']
                for hit_plus in primer_hits_plus: 
                    for hit_minus in primer_hits_minus:
                        size = hit_minus['sstart']-hit_plus['sstart']
                        if size>len_min and size<len_max:
                            # extend this pair:
                            # plus
                            start_plus = hit_plus['sstart']-(hit_plus['qstart']-1)
                            if start_plus<1:    # extension failed
                                continue
                            end_plus = hit_plus['send']+(hit_plus['qlen']-hit_plus['qend'])
                            if end_plus>hit_plus['slen']: # extension failed
                                continue
                            
                            # minus
                            start_minus = hit_minus['sstart']+(hit_minus['qstart']-1)
                            if start_minus>hit_minus['slen']: # extension failed
                                continue
                            end_minus = hit_minus['send']-(hit_minus['qlen']-hit_minus['qend'])
                            if end_minus<1:  # extension failed
                                continue

                            # store this pair
                            amplicons[primer_id][primer_rank].append({
                                'plus': {
                                    'qseqid': hit_plus['qseqid'], # LEFT or RIGHT
                                    'sseqid': sseqid,
                                    'sstart': start_plus,
                                    'send': end_plus
                                },
                                'minus': {
                                    'qseqid': hit_minus['qseqid'], # LEFT or RIGHT
                                    'sseqid': sseqid,
                                    'sstart': end_minus,    # flip
                                    'send': start_minus
                                }
                            })

                            # store regions
                            hits_regions[sseqid+':'+str(start_plus)+'-'+str(end_plus)] = 1
                            hits_regions[sseqid+':'+str(end_minus)+'-'+str(start_minus)] = 1 # flip

    return {'amplicons':amplicons, 'regions_primer': '\n'.join(hits_regions.keys()) }

def filter_Tm(amplicons, query_primer_seq, hits_seqs, Tm_diff=20, max_amplicons=10, use_3_end=False):
    '''
        Input:
            query_primer_seq:
                "LEFT": XXXXXXXXXXXXXX
                "RIGHT": XXXXXXXXXXXX
    '''
    amplicons_filter = {}
    for primer_id in amplicons.keys():
        amplicons_filter[primer_id] = {}
        for primer_rank in amplicons[primer_id]:
            if primer_rank not in amplicons_filter[primer_id]:
                amplicons_filter[primer_id][primer_rank] = []
            for amplicon in amplicons[primer_id][primer_rank]:
                # hit seqs (template)
                region_plus = amplicon['plus']['sseqid']+':'+str(amplicon['plus']['sstart'])+'-'+str(amplicon['plus']['send'])
                region_minus = amplicon['minus']['sseqid']+':'+str(amplicon['minus']['sstart'])+'-'+str(amplicon['minus']['send'])
                hit_seq_plus = hits_seqs[region_plus]
                hit_seq_minus = rev_complement(hits_seqs[region_minus])  # samtools faidx hasn't do rev_complement yet

                # query primers (query)
                query_primer_plus = query_primer_seq[amplicon['plus']['qseqid']]  # LEFT or RIGHT
                query_primer_minus = query_primer_seq[amplicon['minus']['qseqid']]  # LEFT or RIGHT

                # perfect Tm for this primer pair
                Tm_perfect = min(float(NN_Tm(seq=transform_degenerate(query_primer_plus), \
                    compl_seq=complement(transform_degenerate(query_primer_plus)))), \
                    float(NN_Tm(seq=transform_degenerate(query_primer_minus), \
                        compl_seq=complement(transform_degenerate(query_primer_minus)))))

                # 3' end in this amplicon
                if use_3_end is True:
                    # plus
                    if query_primer_plus[-1].upper() != hit_seq_plus[-1].upper():
                        continue
                    # minus
                    if query_primer_minus[-1].upper() != hit_seq_minus[-1].upper():
                        continue

                # Tm in this amplicon
                # plus
                if len(query_primer_plus) != len(hit_seq_plus):
                    raise Exception(f'Not the same length: \n \
                        Primer: {query_primer_plus} \n \
                        Hit: {region_plus} {hit_seq_plus}')
                Tm_plus = float(NN_Tm(seq=transform_degenerate(query_primer_plus), \
                    compl_seq=complement(transform_degenerate(hit_seq_plus))))
                if Tm_plus<Tm_perfect-Tm_diff:
                    continue

                # minus
                if len(query_primer_minus) != len(hit_seq_minus):
                    raise Exception(f'Not the same length: \n \
                        Primer: {query_primer_minus} \n \
                        Hit: {region_minus} {hit_seq_minus}')
                Tm_minus = float(NN_Tm(seq=transform_degenerate(query_primer_minus), \
                    compl_seq=complement(transform_degenerate(hit_seq_minus))))
                if Tm_minus<Tm_perfect-Tm_diff:
                    continue

                # store them
                amplicon['plus']['Tm'] = Tm_plus
                amplicon['plus']['qseq'] = query_primer_plus
                amplicon['plus']['sseq'] = hit_seq_plus
                amplicon['minus']['Tm'] = Tm_minus
                amplicon['minus']['qseq'] = query_primer_minus
                amplicon['minus']['sseq'] = hit_seq_minus
                amplicon['Tm'] = min(Tm_plus, Tm_minus)
                amplicon['region'] = amplicon['plus']['sseqid']+':'+str(amplicon['plus']['sstart'])\
                    +'-'+str(amplicon['minus']['send'])
                amplicon['product_size'] = amplicon['minus']['send']-amplicon['plus']['sstart']+1
                amplicon['product_seq'] = ''
                amplicons_filter[primer_id][primer_rank].append(amplicon)

                if len(amplicons_filter[primer_id][primer_rank])==max_amplicons:
                    break
    return amplicons_filter

def add_amplicon_seq(amplicons, template_file):
    amplicon_regions = ''
    for primer_id in amplicons.keys():
        for primer_rank in amplicons[primer_id]:
            for amplicon in amplicons[primer_id][primer_rank]:
                amplicon_regions += amplicon['region']+'\n'
    amplicon_seqs = faidx(template_file=template_file, region_string=amplicon_regions)
    for primer_id in amplicons.keys():
        for primer_rank in amplicons[primer_id]:
            for (i, amplicon) in enumerate(amplicons[primer_id][primer_rank]):
                amplicons[primer_id][primer_rank][i]['product_seq'] = amplicon_seqs[amplicon['region']]
    return amplicons

def add_isoform_annotation(amplicons, dbs):
    isoform_data = json.load(open(dbs[0]+'.isoforms.json'))
    for primer_id in amplicons.keys():
        for primer_rank in amplicons[primer_id]:
            for (i, amplicon) in enumerate(amplicons[primer_id][primer_rank]):
                sseqid = amplicon['plus']['sseqid']
                qseqid = amplicon['plus']['qseqid'].split('-')[0]
                if qseqid in isoform_data and sseqid in isoform_data[qseqid]:
                    amplicons[primer_id][primer_rank][i]['isoform'] = True
    return amplicons


if __name__ == "__main__":
    blast_out = open('tests/_internal_/query_blast.fa.out.large').read()
    amplicons = filter_len(blast_out=blast_out, len_min=75, len_max=1000)
    # filter_len(blast_out=blast_out, len_min=75, len_max=1000)
    # hits_seqs = faidx(template_file='tests/example.fa', region_string=amplicons['regions_primer'])
    hits_seqs = faidx(template_file='tests/_internal_/Ghir.JGI.genomic', region_string=amplicons['regions_primer'])
    # query_primer_seq={'P1.0.L':'CTTCTGCAATGCCAAGTCCAG',\
    #         'P1.0.R': 'GTGGTGAAGGGTCGGTTGAA'}
    query_primer_seq_lines = open('tests/_internal_/query_blast.fa.out.large.primer.seq').read().splitlines()
    i = iter(query_primer_seq_lines)
    query_primer_seq = dict(zip(i, i))
    report_amplicons = filter_Tm(amplicons['amplicons'], query_primer_seq=query_primer_seq, hits_seqs=hits_seqs)
    #report_amplicons = add_amplicon_seq(amplicons=report_amplicons, template_file='tests/example.fa')
    report_amplicons = add_amplicon_seq(amplicons=report_amplicons, template_file='tests/_internal_/Ghir.JGI.genomic')
    print(json.dumps(report_amplicons, indent=4))
