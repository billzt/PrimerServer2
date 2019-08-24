#!/usr/bin/env python3

'''PrimerServer2: a high-throughput primer design and specificity-checking platform
Github: https://github.com/billzt/PrimerServer2
This is a CLI script for preparing the isoform database (in JSON format) for PrimerServer2
'''

import json
import argparse

from primerserver2.core.qRT import get_AS
from primerserver2.core import version

def main():
    parser = argparse.ArgumentParser(description='primertool-isoforms: preparing the isoform database for PrimerServer2', \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s '+version.get())
    parser.add_argument('template', help='the template file in FASTA format')
    parser.add_argument('gff', help='the annotation file in GFF3 format')
    parser.add_argument('-f', '--features', help='features in GFF3 (the 3rd column) that should be parsed as isoforms. \
        Multiple features separated by comma are allowed', default='mRNA')
    parser.add_argument('-i', '--id', help='the key in GFF3 (in the last column) that should be treated as the isoform IDs. \
        ', default='ID')
    parser.add_argument('-r', '--remove-version', help='remove the version code (such as .1, .2, ...) from isoform IDs. \
        ', action='store_true')
    args = parser.parse_args()

    isoform_data = get_AS(gff_file=args.gff, features=args.features.split(','), rna_print_ID_key=args.id, remove_version=args.remove_version)

    with open(args.template+'.isoforms.json', 'w') as f:
        json.dump(isoform_data, f, indent=4)

if __name__ == "__main__":
    main()