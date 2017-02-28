#!/usr/bin/env python
"""
Author: Rens Holmer
"""

__author__ = 'rensholmer'

import sys
import gff_toolkit as gt

def get_partials(partial_id_file):
    """
    Returns a list with gene IDs that are partially transferred
    """
    partials = set()
    with open(partial_id_file,'rU') as fh:
        for line in fh:
            parts = line.strip()
            if not parts:
                continue
            ID = parts.lstrip('/locus_tag="').split('.')[0]
            partials.add(ID)
    return partials

def main(ratt_file,fasta_file,partial_id_file):
    partials = get_partials(partial_id_file)
    gff = gt.parser(ratt_file,filetype='ratt',fasta_file=fasta_file)
    for gene in gff.getitems(featuretype='gene'):
        partial = False
        print gene.ID
        if gene.ID in partials:
            partial = True
        transcripts = [transcript for transcript in gff.get_children(gene,featuretype='mRNA')]
        if len(transcripts) > 1:
            sorted_transcripts = sorted(transcripts,key=lambda x: len(x.seq))
            longest = sorted_transcripts.pop()
        else:
            longest = transcripts[0]
        if partial:
            transcript.set_attribute('Partial','True')
        else:
            transcript.set_attribute('Partial','False')
        #make sure the gene line is correct 
        new_gene = gt.GffSubPart(*longest.gff_fields)
        new_gene.ID = gene.ID
        new_gene.featuretype = 'gene'
        new_gene.parents = []
        print '\t'.join(new_gene.gff_fields)

        print '\t'.join(longest.gff_fields)
        for cds in gff.get_children(longest,featuretype='CDS'):
            if partial:
                cds.set_attribute('Partial','True')
            else:
                cds.set_attribute('Partial','False')
            print '\t'.join(cds.gff_fields)

if __name__ == '__main__':
    if len(sys.argv) == 4:
        main(*sys.argv[1:])
    else:
        print 'Usage: python {0} <ratt_transfer.gff> <genome.fasta> <ratt_not_transferred.txt>'.format(sys.argv[0])
        print 'Notes:\t- Format the output of ratt'
        print '\t- if a gene ID is in both <ratt_transfer.gff> and <ratt_not_transferred.txt> it will be flagged as partial transfer'
