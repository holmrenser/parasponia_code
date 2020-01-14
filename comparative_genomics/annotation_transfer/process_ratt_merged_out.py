#!/usr/bin/env python
"""
Author: Rens Holmer
"""

__author__ = 'rensholmer'

import sys
import gff_toolkit as gt

START_CODONS = ('ATG','CTG')
STOP_CODONS = ('TAA','TGA','TAG')

def get_first_cds(transcript):
    """
    Returns the first CDS feature of a given mRNA
    """
    subfeatures = transcript.container.get_children(transcript, featuretype = 'CDS')
    sorted_subfeatures = sorted(subfeatures, key = lambda x: x.get_start(), reverse = transcript.reverse)
    return sorted_subfeatures[0]

def main(gff_file, fasta_file):
    gff = gt.parser(gff_file, fasta_file = fasta_file)
    for gene in gff.getitems(featuretype = 'gene'):
        correct = True
        for transcript in gff.get_children(gene, featuretype = 'mRNA'):
            first_cds = get_first_cds(transcript)
            phase = first_cds.phase
            #internal stops are not allowed
            if '*' in transcript.pep[:-1]:
                correct = False
            #more than 50% X amino acids is not allowed
            if transcript.pep.count('X') > len(transcript.pep) / 2:
                correct = False
            #very short annotations are discarded
            if len(transcript.pep) < 25:
                correct = False
            #leading X is not allowed
            if transcript.pep[0] == 'X':
                correct = False
            #remove single exon predictions without start and stop
            if transcript.seq[phase:][0:3] not in START_CODONS and transcript.seq[-3:] not in STOP_CODONS and len(transcript.children) == 1:
                correct = False
        if not correct:
            continue
        for sub in gff.get_children(gene):
            print '\t'.join(sub.gff_fields)

if __name__ == '__main__':
    main(*sys.argv[1:])
