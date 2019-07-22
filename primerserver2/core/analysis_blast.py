import re
import json
import sys

from primerserver2.core.make_sites import faidx
from primerserver2.core.Santalucia_NN_Tm import complement, rev_complement, NN_Tm
from primerserver2.core.make_primers import make_primers

def filter_len(blast_out, len_min, len_max):
    '''
        return 
            'amplicons': hit-pairs 
            'regions_primer': a region string for each primer used by samtools faidx;
    '''

    # collect hits data
    primer_hits = []
    for line in blast_out.splitlines():
        (qseqid, qlen, qstart, qend, sseqid, slen, sstart, send, sstrand) = line.strip().split('\t')
        primer_hits.append({
            'qseqid': qseqid,   # LEFT or RIGHT
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
    amplicons = []
    hits_regions = {}   # make region file (FASTA)

    if len(primer_hits)==1:
        return []

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
                amplicons.append({
                    'plus': {
                        'qseqid': primer_hits[i]['qseqid'], # LEFT or RIGHT
                        'sseqid': primer_hits[i]['sseqid'],
                        'sstart': start_plus,
                        'send': end_plus
                    },
                    'minus': {
                        'qseqid': primer_hits[j]['qseqid'], # LEFT or RIGHT
                        'sseqid': primer_hits[j]['sseqid'],
                        'sstart': end_minus,    # flip
                        'send': start_minus
                    }
                })

                # store regions
                hits_regions[primer_hits[i]['sseqid']+':'+str(start_plus)+'-'+str(end_plus)] = 1
                hits_regions[primer_hits[j]['sseqid']+':'+str(end_minus)+'-'+str(start_minus)] = 1 # flip

    return {'amplicons':amplicons, 'regions_primer': '\n'.join(hits_regions.keys()) }

def filter_Tm(amplicons, query_primer_seq, hits_seqs, Tm_diff=20, max_amplicons=10, use_3_end=False):
    '''
        Input:
            query_primer_seq:
                "LEFT": XXXXXXXXXXXXXX
                "RIGHT": XXXXXXXXXXXX
    '''
    amplicons_filter = []
    for amplicon in amplicons:
        # hit seqs (template)
        region_plus = amplicon['plus']['sseqid']+':'+str(amplicon['plus']['sstart'])+'-'+str(amplicon['plus']['send'])
        region_minus = amplicon['minus']['sseqid']+':'+str(amplicon['minus']['sstart'])+'-'+str(amplicon['minus']['send'])
        hit_seq_plus = hits_seqs[region_plus]
        hit_seq_minus = rev_complement(hits_seqs[region_minus])  # samtools faidx hasn't do rev_complement yet

        # query primers (query)
        query_primer_plus = query_primer_seq[amplicon['plus']['qseqid']]  # LEFT or RIGHT
        query_primer_minus = query_primer_seq[amplicon['minus']['qseqid']]  # LEFT or RIGHT

        # perfect Tm for this primer pair
        Tm_perfect = min(float(NN_Tm(seq=query_primer_plus, compl_seq=complement(query_primer_plus))), \
            float(NN_Tm(seq=query_primer_minus, compl_seq=complement(query_primer_minus))))

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
        Tm_plus = float(NN_Tm(seq=query_primer_plus, compl_seq=complement(hit_seq_plus)))
        if Tm_plus<Tm_perfect-Tm_diff:
            continue

        # minus
        if len(query_primer_minus) != len(hit_seq_minus):
            raise Exception(f'Not the same length: \n \
                Primer: {query_primer_minus} \n \
                Hit: {region_minus} {hit_seq_minus}')
        Tm_minus = float(NN_Tm(seq=query_primer_minus, compl_seq=complement(hit_seq_minus)))
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
        amplicons_filter.append(amplicon)

        if len(amplicons_filter)==max_amplicons:
            break

    return amplicons_filter

def add_amplicon_seq(amplicons, template_file):
    amplicon_regions = '\n'.join([x['region'] for x in amplicons])
    amplicon_seqs = faidx(template_file=template_file, region_string=amplicon_regions)
    for (i, amplicon) in enumerate(amplicons):
        amplicons[i]['product_seq'] = amplicon_seqs[amplicon['region']]
    return amplicons


if __name__ == "__main__":
    blast_out = open('tests/_internal_/query_blast.fa.out').read()
    amplicons = filter_len(blast_out=blast_out, len_min=75, len_max=1000)
    hits_seqs = faidx(template_file='tests/example.fa', region_string=amplicons['regions_primer'])
    report_amplicons = filter_Tm(amplicons['amplicons'], query_primer_seq={'LEFT':'CTTCTGCAATGCCAAGTCCAG',\
        'RIGHT': 'GTGGTGAAGGGTCGGTTGAA'}, hits_seqs=hits_seqs)
    report_amplicons = add_amplicon_seq(amplicons=report_amplicons, template_file='tests/example.fa')
    print(json.dumps(report_amplicons, indent=4))
