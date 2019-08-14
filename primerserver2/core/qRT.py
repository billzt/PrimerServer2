import json
import re

from operator import itemgetter

def get_junction(gff_file, features=['exon'], rna_ID_key='Parent', rna_print_ID_key='Parent', remove_version=True):
    with open(gff_file) as f:
        id2exons = {}
        id2strand = {}
        id2print = {}
        for line in f:
            if line.startswith('#'):
                continue
            (feature, start, end, strand, note) = itemgetter(2,3,4,6,-1)(line.strip().split('\t'))
            if feature not in features:
                continue
            notes = note.split(';')
            for x in notes:
                if f'{rna_ID_key}=' in x:
                    rna_ID = x.replace(f'{rna_ID_key}=', '')
                    if rna_ID not in id2exons:
                        id2exons[rna_ID] = []
                if f'{rna_print_ID_key}=' in x:
                    rna_print_ID = x.replace(f'{rna_print_ID_key}=', '')
                    if remove_version is True:
                        rna_print_ID = re.sub(r'\.\d+$', '', rna_print_ID)
            id2exons[rna_ID].append([int(start), int(end)])
            id2strand[rna_ID] = strand
            id2print[rna_ID] = rna_print_ID
        
        junction_data = {}
        for rna_ID in id2strand.keys():
            exons = id2exons[rna_ID]
            if len(exons)==1:
                continue
            strand = id2strand[rna_ID]
            rna_print_ID = id2print[rna_ID]
            
            # sort exons according to strand
            if strand=='-':
                exons = sorted(exons, key=lambda x:x[0], reverse=True)
            
            # junctions
            junction_poses = []
            junction_pos = 0
            for i in range(0, len(exons)-1):
                (start1, end1) = exons[i]
                (start2, end2) = exons[i+1]
                intron_length = start2-end1-1 if strand=='+' else start1-end2-1
                junction_pos += end1-start1+1
                junction_poses.append({'pos': junction_pos, 'intron_len': intron_length})
            junction_data[rna_print_ID] = {'strand': strand, 'junctions':junction_poses}
        return junction_data
            

def get_AS(gff_file, features=['mRNA'], rna_print_ID_key='ID', remove_version=True):
    with open(gff_file) as f:
        gene2rnas = {}
        rna2gene = {}
        for line in f:
            if line.startswith('#'):
                continue
            (feature, note) = itemgetter(2,-1)(line.strip().split('\t'))
            if feature not in features:
                continue
            notes = note.split(';')
            for x in notes:
                if 'Parent=' in x:
                    gene_ID = x.replace('Parent=', '')
                    if gene_ID not in gene2rnas:
                        gene2rnas[gene_ID] = []
                if f'{rna_print_ID_key}=' in x:
                    rna_print_ID = x.replace(f'{rna_print_ID_key}=', '')
                    if remove_version is True:
                        rna_print_ID = re.sub(r'\.\d+$', '', rna_print_ID)
                        rna2gene[rna_print_ID] = gene_ID
            gene2rnas[gene_ID].append(rna_print_ID)
        
        rna2partners = {}
        for (rna_print_ID, gene_ID) in rna2gene.items():
            partners = gene2rnas[gene_ID]
            rna2partners[rna_print_ID] = partners
        return rna2partners



if __name__ == "__main__":
    junction_data = get_junction(gff_file='tests/_internal_/ref_GRCh38.p12_top_level.gff3', rna_print_ID_key='transcript_id')
    as_data = get_AS(gff_file='tests/_internal_/ref_GRCh38.p12_top_level.gff3', rna_print_ID_key='transcript_id')
    #print(json.dumps(junction_data, indent=4))
    print(json.dumps(as_data, indent=4))