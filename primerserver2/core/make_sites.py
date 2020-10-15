
import re
import os
import json
import subprocess

from primerserver2.core import global_var

def get_template_len(template_file):
    '''
        return a dict of each seq's len based on the faidx
    '''
    template_len_for = {}
    with open(template_file+'.fai') as f:
        for line in f:
            (id, length) = line.split('\t')[0:2]
            template_len_for[id] = int(length)
    return template_len_for

def faidx(template_file, region_string):
    '''
        return a list of sequence dict using samtools faidx based on sequence file (fa)
        and region_string: 
            seq1:1-100
            seq2:2-200
            ...
    '''
    if os.path.isfile(template_file)==False:
        raise Exception(f'File {template_file} is not found')

    cmd = f'samtools faidx -r - {template_file}'
    p = subprocess.run(cmd, input=region_string, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True, encoding='ascii')
    seqs = p.stdout.split('>')
    del seqs[0]
    result_seqs = {}
    for seq in seqs:
        seq_split = seq.splitlines()
        result_seqs[seq_split[0]] = ''.join(seq_split[1:])
    return result_seqs

def read_junction(template_file):
    if os.path.isfile(f'{template_file}.junctions.json') is True:
        return json.load(open(f'{template_file}.junctions.json'))
    else:
        return {}

def get_junction_poses(junction_data, id):
    if id not in junction_data:
        return []
    else:
        return [x['pos'] for x in junction_data[id]['junctions'] if x['intron_len']>100]


def judge_input_type(query):
    '''
        whether input by pos or by seq?
        return: pos or seq
    '''
    if '>' in query:
        return 'seq'
    else:
        return 'pos'

def build_by_seq(query, primer_type, primer_num_return=30, size_min=70, \
    size_max=1000, pick_internal=False, Tm_opt=60):
    '''
        Input:
            query: a string in FASTA format
        Return:
            a list of primer sites dict that can be passed to the design_primer module
            [{
                'id': , 
                'template': ,
                'type': ,
                'pos': ,
                'length':,
                'size_min':,
                'size_max':,
                'primer_num_return':,
            }]
    '''
    if global_var.stop_run is True:
        return {'error': 'Stop running'}

    primer_sites = []
    primer_site_rank = 0
    for site in query.split('>'):
        # site name and templates
        if site.strip()=='':
            continue
        (name_data, *seq_data) = site.strip().split('\n')
        (name, *dummy) = name_data.strip().split()
        seq = ''.join(seq_data)
        if re.search('[^ATGCNRYMKSWHBVD]', seq.replace('[','').replace(']',''), re.RegexFlag.IGNORECASE) is not None:
            continue
        primer_site_rank += 1

        # SEQUENCE_TARGET, FORCE_END
        if primer_type=='SEQUENCE_TARGET' or primer_type=='FORCE_END':
            if '[' not in seq or ']' not in seq:
                return {'error': f'seq {seq} has no [target regions]'}
            if seq.count('[')>1 or seq.count(']')>1:
                return {'error': f'multiple [target regions] in {seq}'}
            pos = seq.find('[')
            length = seq.find(']')-seq.find('[')
            primer_sites.append({
                'id': name, 
                'template': seq.replace('[','').replace(']',''),
                'type': primer_type,
                'pos': pos,
                'length':length,
                'size_min':size_min,
                'size_max':size_max,
                'primer_num_return':primer_num_return,
                'pick_internal': pick_internal,
                'junction': [],
                'Tm_opt': Tm_opt
            })
        # SEQUENCE_INCLUDED_REGION
        if primer_type=='SEQUENCE_INCLUDED_REGION':
            if '[' not in seq and ']' not in seq:
                pos = 1
                length = len(seq)-1
            else:
                pos = seq.find('[')
                length = seq.find(']')-seq.find('[')
            primer_sites.append({
                'id': 'S'+str(primer_site_rank)+'-'+str(pos)+'-'+str(length), 
                'template': seq.replace('[','').replace(']',''),
                'type': primer_type,
                'pos': pos,
                'length':length,
                'size_min':size_min,
                'size_max':size_max,
                'primer_num_return':primer_num_return,
                'pick_internal': pick_internal,
                'junction': [],
                'Tm_opt': Tm_opt
            })
    return primer_sites



def build_by_pos(query, template_file, primer_type, primer_num_return=30, size_min=70, \
    size_max=1000, pick_internal=False, use_junction=False, Tm_opt=60):
    '''
        Input:
            query: a string in multi-lines
        Return:
            a list of primer sites dict that can be passed to the design_primer module
            [{
                'id': , 
                'template': ,
                'type': ,
                'pos': ,
                'length':,
                'size_min':,
                'size_max':,
                'primer_num_return':,
            }]
    '''
    if global_var.stop_run is True:
        return {'error': 'Stop running'}
        
    template_len_for = get_template_len(template_file)
    primer_sites = []
    retrieve_region2raw_region = {}
    for line in query.splitlines():
        if line.strip()=='':
            continue
        query_data = re.split(r'\s+', line.strip())
        if query_data[0]=='':
            return {'error': f'Your input {line.strip()} has no template seqs'}
        (chr, pos, length) = (query_data[0], -1, -1)
        if chr not in template_len_for:
            return {'error': f'Your input {line.strip()} has no template seqs'}
        chr_len = template_len_for[chr]
        if len(query_data)>1 :     # pos is provided
            pos = int(query_data[1].replace(',', ''))
            if pos<0 or pos>chr_len:
                return {'error': f'Your input {line.strip()} has wrong pos: {pos}: it must be between 1 and {chr_len}'}
        if len(query_data)>2 :      # len is provided
            length = int(query_data[2])
            if pos+length > chr_len:
                return {'error': f'Your input {line.strip()} has wrong region: {pos}+{length}: it cannot be exceed {chr_len}'}
        if len(query_data)>3 :
            size_min = int(query_data[3])
            size_max = int(query_data[4])
        
        # some adjustion
        if primer_type=='SEQUENCE_TARGET' or primer_type=='FORCE_END':
            if pos==-1:
                return {'error': f'Your input {line.strip()} lack position, which is required by {primer_type}'}
            if length==-1:
                length=1
        if primer_type=='SEQUENCE_INCLUDED_REGION':
            if pos==-1:
                pos=1
            if length==-1:
                length=chr_len-pos

        # retrieve_start and retrieve_end are used by samtools to extract template sequences
        retrieve_start = max(pos-size_max, 1)
        retrieve_end = pos+length+size_max
        retrieve_region2raw_region[f'{chr}:{retrieve_start}-{retrieve_end}'] = [chr, pos, length, size_min, size_max, retrieve_start]
    
    if (len(retrieve_region2raw_region)==0):
        return {'error': f'Your input has no template seqs'}
    retrieve_region_string = '\n'.join(retrieve_region2raw_region.keys())
    result_seqs = faidx(template_file, retrieve_region_string)
    if use_junction is True:
        junction_data = read_junction(template_file)
    for (result_seq_id, result_seq) in result_seqs.items():
        (chr, pos, length, size_min, size_max, retrieve_start) = retrieve_region2raw_region[result_seq_id]
        if result_seq=='':
            return {'error': f'Your input: seq={chr} pos={pos} length={length} has no template seqs'}
        junctions = []
        if use_junction is True:
            junctions = get_junction_poses(junction_data, chr)
        primer_sites.append({
            'id': chr+'-'+str(pos)+'-'+str(length), 
            'template': result_seq,
            'type': primer_type,
            'pos': pos-retrieve_start,
            'length':length,
            'size_min':size_min,
            'size_max':size_max,
            'primer_num_return':primer_num_return,
            'pick_internal': pick_internal,
            'junction': junctions,
            'Tm_opt': Tm_opt
        })
    
    return primer_sites

if __name__ == "__main__":
    global_var.init()
    import argparse
    parser = argparse.ArgumentParser(description='Your title here')
    parser.add_argument('query_design', help='input file', type=argparse.FileType('r'))
    parser.add_argument('type', help='SEQUENCE_TARGET, SEQUENCE_INCLUDED_REGION, FORCE_END')
    args = parser.parse_args()
    query = args.query_design.read()
    input_type = judge_input_type(query)
    if input_type=='pos':
        primer_sites = build_by_pos(query=query, template_file='tests/example.fa', primer_type=args.type, primer_num_return=30)
    else:
        primer_sites = build_by_seq(query=query, primer_type=args.type, primer_num_return=30)
    print(json.dumps(primer_sites, indent=4))
