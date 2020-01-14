#!/usr/bin/env python
"""
Author: Rens Holmer
"""

__author__ = 'rensholmer'

import sys
import gff_toolkit as gt

def main(gff_file, fasta_file):
    """Converts a ratt formatted gff3 file (which is not properly formatted according to the specs)
    to a proper gff3 file. Conversion is done by the gff_toolkit
    Prints gff3 lines to stdout

    Args:
        gff_file (str): ratt formatted gff3 file
        fasta_file (str): fasta file with the genome sequence to which the annotation belongs
    Returns:
        None
    """
    gff = gt.parser(gff_file, fasta_file = fasta_file, filetype='ratt')
    for gene in gff.getitems(featuretype='gene'):
        for sub in gff.get_children(gene):
            print '\t'.join(sub.gff_fields)

if __name__ == '__main__':
    if len(sys.argv) == 3:
        main(*sys.argv[1:])
    else:
        print 'Usage: python {0} <ratt formatted gff file> <genome fasta file> '.format(sys.argv[0])
