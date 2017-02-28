#!/usr/bin/env python
"""
Author: Rens Holmer
External dependencies: mafft
All other dependencies are installable with pip
"""

__author__ = 'rensholmer'

import sys
import gff_toolkit as gt
import numpy as np
from subprocess import Popen,PIPE
import json
import pysam
from itertools import groupby
import networkx as nx

class AnnotationCluster(object):
    """
    Holds overlapping gene annotations on two genomes and methods to choose between them
    Keeping track of original annotations and RATT transferred annotations is done with the indices 1,2,-1,-2
    In the case of Parasponia/Trema this ammounts to the following:
        1 : Parasponia
        2 : Trema
        -1 : Trema transferred to Parasponia
        -2 : Parasponia transferred to trema
    """
    #these are the only allowed indices
    _indices = (1,2,-1,-2)
    #these are all possible cominations of the indices
    _pairs = ((1,2),(1,-2),(-1,2),(-1,-2))
    def __init__(self, anndic, bamdic, outdic = None, *args, **kwargs):
        #combinations of overlapping annotations
        self.combinations = ((0,0),(0,0))
        #which genes are partial (no start/ no stop) (sets of IDs)
        self.partials = {i:set() for i in self._indices}
        #pysam objects of bamfiles
        self.bamdic = bamdic
        #gff_toolkit objects of annotations
        self.anndic = anndic
        #filehandles for output
        self.outdic = outdic
        #relative intronscore
        self.intronscore_percentage = {p:1 for p in self._pairs} 
        #absolute intronscore
        self.intronscore_absolute = {p:0 for p in self._pairs}
        #protein score 
        self.proteinscore = {p:0 for p in self._pairs}
        #genes in cluster 
        self.cluster = {i:set() for i in self._indices}
        #introns per ann
        self.introns = {i:{0:0} for i in self._indices}
        #number of introns with support
        self.intron_match = {i:0 for i in self._indices}
        #intron percentage identity: number of supported introns / total introns
        self.intron_identity = {i:None for i in self._indices}
        #all potential splicesites per ann
        self.splicesites = {i:{} for i in self._indices}
        #coverage per ann
        self.coverage = {i:[] for i in self._indices}
        #peptide sequence per ann
        self.pep = {i:'' for i in self._indices}
        #final outcome
        self.solution =    None
    def __getitem__(self, value):
        """
        Extracts the annotations based on one of the allowed indices (1,2,-1,-2)
        Returns a list
        """
        if value not in self._indices:
            e = '<{0}> is not a valid index for object of type AnnotationCluster, valid options are [1,2,-1,-2]'.format(value)
            raise IndexError(e)
        return list(self.cluster[value])
    def __str__(self):
        return str(self.cluster)
    def add(self, genename, index):
        """
        add geneID to cluster for the given index
        """
        if index not in self._indices:
            e = '{0} is not a valid index for object of type AnnotationCluster'.format(index)
            raise ValueError(e)
        self.cluster[index].add(genename)
    def classify(self):
        """
        Determine type of overlap based on number of annotations per track
        Sets self.combinations to ((ann01,ratt01),(ann02,ratt02))
        """
        self.combinations = ((len(self.cluster[1]),len(self.cluster[-1])),(len(self.cluster[2]),len(self.cluster[-2])))
        for n in self.cluster:
            gff = self.anndic[n]
            for gene_name in self.cluster[n]:
                for gene in gff[gene_name]:
                    if gene.attributes.get('Partial',False) == ['True']:
                        self.partials[n].add(gene.ID)

    def size(self):
        """
        return total number of annotations found in all four possible tracks
        """
        return len([gene for genes in self.cluster.values() for gene in genes])
    def _get_introns(self, index):
        """
        this is called by self.score_introns()
        sets self.introns
        """
        gff = self.anndic[index]
        for gene_name in self.cluster[index]:
            #select transcripts
            for transcript in gff.get_children(gene_name, featuretype = 'mRNA'):
                #transferred genes that are partial (no start/no stop) are allowed if the original annotation is partial
                partial = False
                if index < 0:
                    if transcript.seq[transcript.phase:][0:3] not in ('CTG','ATG') or transcript.pep[-1] != '*':
                        if transcript.ID not in self.partials[-index]:
                            partial = True
                #premature stops are never allowed
                if '*' in transcript.pep[:-1] or partial:
                    self.introns[index] = {0:0}
                    return

                cds_list = 
                cds_list = sorted(gff.get_children(transcript, featuretype = 'CDS'), key = lambda x: x.get_start(), reverse = transcript.reverse)
                x = max(self.introns[index])
                for i,c in enumerate(cds_list):
                    if i != len(cds_list) - 1:
                        self.introns[index][i + x] = c.get_end()
                    if i != 0:
                        ii = (self.introns[index][i + x - 1],c.get_start())
                        self.introns[index][i + x - 1] = (min(ii),max(ii) - 1)

    def _get_intron_score(self, index):
        """
        this is called by self.score_introns()
        sets self.coverage
        """
        gff = self.anndic[index]
        for gene_name in self.cluster[index]:
            #for gene in gff[gene_name]:
            for transcript in gff.get_children(gene_name, featuretype = 'mRNA'):
                #iterate over indexed bamfiles
                cov = 0
                for bam in self.bamdic[index]:
                    #get reads mapping within annotation
                    for read in bam.fetch(transcript.seqid,transcript.start,transcript.end):
                        cov += read.infer_query_length()
                        pos = read.reference_start
                        #parse cigarstring for all possible splicesites 'cigar[0] == 3' equals 'cigar N'
                        for cigar in read.cigartuples:
                            if cigar[0] != 3:
                                pos += cigar[1]
                            else:
                                s1 = pos
                                pos += cigar[1]
                                s2 = pos
                                if (s1,s2) in self.splicesites:
                                    self.splicesites[index][(s1,s2)] += 1
                                else:
                                    self.splicesites[index][(s1,s2)] = 1
                    mean_cov = float(cov) / len(transcript.seq)
                    self.coverage[index].append(mean_cov)

    def score_introns(self):
        """
        Identifies all introns in all annotations in the cluster
        Calculates coverage of the predicted introns in the RNA-seq data
        If the percentage coveraged introns is equal we calculate the absolute amount of covered introns
        """
        #get introns
        for index in self._indices:
            self._get_introns(index)
            
        #check if introns are all equal (automatically true for single exon genes)
        if self.introns[1] != self.introns [-1] or self.introns[2] != self.introns[-2]:
            self._get_intron_score(index)
            #check if all annotations have sufficient coverage
            if not [cov for cov_list in self.coverage.values() for cov in cov_list if cov <= 1]:
                #iterate over annotations in introns dictionary
                for n in self.introns:
                    #iterate over introns in annotation
                    for intron in self.introns[n]:
                        #check if predicted intron in potential splicesites
                        if self.introns[n][intron] in self.splicesites[n]:
                            self.intron_match[n] += 1
                    for key in self.intronscore_percentage:
                        if n in key:
                            if self.introns[n]:
                                self.intronscore_percentage[key] *= float(self.intron_match[n]) / len(self.introns[n])
                                self.intronscore_absolute[key] += self.intron_match[n]
        else:
            #all introns are equal, return original annotations
            return {(1,2):100.0}
        #calculate percentage of predicted introns covered by RNA-seq
        percentagewin = {key:value for key,value in self.intronscore_percentage.iteritems() if value == max(self.intronscore_percentage.values()) and self.check_solution(key)}
        if len(percentagewin) == 1:
            return percentagewin
        #calculate absolute amount of predicted introns covered by RNA-seq
        absolutewin = {key:value for key,value in self.intronscore_absolute.iteritems() if value == max(self.intronscore_absolute.values()) and self.check_solution(key)}
        return absolutewin
    
    def _get_protein(self, index):
        """
        This is called by self._score_protein()
        """
        gff = self.anndic[index]
        for gene_name in self.cluster[index]:
            for transcript in gff.get_children(gene_name, featuretype = 'mRNA'):
                #annotations with internal stop codons are not taken into account
                if '*' in transcript.pep[:-1]:
                    self.pep[index] = ''
                    return
                else:
                    self.pep[index] += transcript.pep

    def _score_protein(self):
        """
        Make a multiple sequence alignment (MSA) of the protein translation of two annotations with mafft
        """
        for index in self._indices:
            self._get_protein(index)

        for ann1,ann2 in self.proteinscore:
            if not self.pep[ann1] or not self.pep[ann2]:
                self.proteinscore[(ann1,ann2)] = 0
                continue
            alignment = []
            match = 0
            mafft_string = '>{0} {1}\n{2}\n'.format(ann1, self.cluster[ann1], self.pep[ann1])
            mafft_string += '>{0} {1}\n{2}\n'.format(ann2, self.cluster[ann2], self.pep[ann2])
            mafft_command = 'mafft --anysymbol -'
            p = Popen(mafft_command.split(), stdin = PIPE, stdout = PIPE, stderr = PIPE)
            out,err = p.communicate(input = mafft_string)
            out_lines = (x for x in out.split('\n') if x.strip())
            faiter = (x[1] for x in groupby(out_lines, lambda line: line[0] == '>'))
            for header in faiter:
                header = header.next()[1:].strip()
                seq = ''.join(s.strip() for s in faiter.next())
                alignment.append(seq)
            for pair in zip(alignment[0], alignment[1]):
                if pair[0] == pair[1]:
                    match += 1
            self.proteinscore[(ann1, ann2)] = float(match) / len(alignment[0])

        return {(ann1,ann2):value for (ann1,ann2),value in self.proteinscore.iteritems() if value == max(self.proteinscore.values())}
        
        #if percentage identity is below 90%, keep the original annatations
        if win.values()[0] < 0.90:
            return {(1,2):False}
        else:
            return win
    def _score(self):
        intronwin = self.score_introns()
        if len(intronwin) == 1:
            self.solution = intronwin.keys()[0]
        else:
            proteinwin = self._score_protein()
            if len(proteinwin) == 1:
                self.solution = proteinwin.keys()[0]
            else:
                self.solution = (1,2)
    def choose(self):
        """
        Determine best pair based on transferability or transcriptome and/or protein alignments
        combinations: ((ann01,ratt01),(ann02,ratt02))
        """
        self.classify()
        if self.combinations == ((1,1),(1,1)) and self._exact_match():
            self.solution = (1,2)
        elif self.combinations == ((1,1),(1,0)):
            self._score()
        elif self.combinations == ((1,1),(0,1)):
            self._score()
        elif self.combinations == ((1,0),(0,1)):
            self.solution = (1,-2)
        elif self.combinations == ((1,0),(0,0)):
            self.solution = (1,0)
        else:
            self._score()

        #check if solution has protein coding genes
        if not self.check_solution(self.solution):
            if self.solution == (1,0):
                return False
            self.solution = (1,2)
            if not self.check_solution(self.solution):
                raise Exception('ORIGINAL ANNOTATIONS ARE WRONG, THIS SHOULD NOT HAPPEN')
        #set putative ortholog attribute
        self.set_ortholog_attribute()

        return True 

    def _exact_match(self):
        """
        Determine if all CDS exons are identical in both annotations and transfers
        """
        #first get transcript names from both species
        gene_name1 = list(self.cluster[1])[0]
        gene_name2 = list(self.cluster[2])[0]

        for n in (1,-2):
            #create two sorted lists with CDS features
            c1 = sorted(self.anndic[n].get_children(gene_name1, featuretype = 'CDS') , key = lambda x:x.get_start())
            c2 = sorted(self.anndic[-n].get_children(gene_name2, featuretype = 'CDS'), key = lambda x:x.get_start())
            
            for x1,x2 in zip(c1,c2):
                if x1.start != x2.start or x1.end != x2.end:
                    return False
        return True 

    def check_solution(self,solution):
        for s in solution:
            if s not in self._indices + (0,):
                e = '{0} is not a valid solution!'.format(solution)
                raise ValueError()
            #original solutions is always valid
            if s >= 0:
                continue
            #partial results only valid if the original annotation is partial
            for n in self.cluster[s]:
                gff = self.anndic[s]
                for gene in gff[n]:
                    if gene.attributes.get('partial',False) == ['True']:
                        if gene.ID not in self.partial[-s]:
                            return False
            for i in self.introns[s]:
                if self.introns[s][i] == 0:
                    continue
                if max(self.introns[s][i]) - min(self.introns[s][i]) < 10:
                    return False
        return True

    def set_ortholog_attribute(self):
        """
        Add an attribute to the gff annotation referencing the putative orthologous annotation that was used in this method
        """
        s1 = self.solution[0]
        s2 = self.solution[1]
        if s2 == 0:
            for transcript in [y for x in self.cluster[s1] for y in self.anndic[s1][x]]:
                transcript.attributes.setdefault('putative orthologs',[]).append('None')
            return True
        #if annotations are transferred, orthologs are set based on ID
        #negative index is transferred gene, positive index is original
        #if multiplied indices are negative, one is transferred and the other is original, thus orthologs can be set based on ID
        if s1 * s2 < 0:
            for s in s1,s2:
                for ID in self.cluster[s]:
                    for transcript in self.anndic[s][ID]:
                        transcript.set_attribute('putative orthologs',transcript.ID)
        else:
            lists = {}
            lists[1] = [y for x in self.cluster[s1] for y in self.anndic[s1][x]]
            lists[2] = [y for x in self.cluster[s2] for y in self.anndic[s2][x]]
            toggle = {1:2,2:1}
            for t in toggle:
                for transcript in lists[t]:
                    for other in lists[toggle[t]]:
                        transcript.set_attribute('putative orthologs',other.ID)
        return True

    def save(self):
        """
        Save to filehandles in self.outdic
        """
        if not self.outdic:
            e = 'Must set Cluster.outdic value first!'
            raise Exception(e)
        if self.solution == None:
            e = 'I don\'t know what to save! First run cluster.score()'
            raise Exception(e)
        #determine if singleton
        if self.solution[1] == 0:
            singleton = 1
        else:
            singleton = 0
        for solution in self.solution:
            #solution == 0 means singleton gene in the other annotation
            if solution == 0:
                continue
            gff = self.anndic[solution]
            outfile = self.outdic[solution][singleton]
            for gene_name in self.cluster[solution]:
                for gene in gff[gene_name]:
                    for gffsubpart in gff.get_children(gene_name):
                        gffsubpart.attributes.pop('Partial',None)
                        line = '\t'.join(gffsubpart.gff_fields) + '\n'
                        outfile.write(line)

def _overlap(ann,ratt):
    """
    This uses pybedtools, not advised
    Input two gff annotations on the same genome (ann & ratt)
    Return networkx Graph object with overlapping genemodels
    """
    print 'Determining {0} vs {1} overlap'.format(ann,ratt)
    G = nx.Graph()
    ann_bed = pbt.BedTool(ann.stringify(),from_string=True)
    ratt_bed = pbt.BedTool(ratt.stringify(),from_string=True)
    intersect = ann_bed.intersect(ratt_bed,wao=True)
    for parts in intersect:
        #add annotation node
        if parts[2] == 'CDS':
            ann_attributes = { a.split('=')[0]: a.split('=')[1].split(',') for a in parts[8].split(';')}
            ann_ID = ann_attributes['ID'][0]
            for x in (parent for feature in ann[ann_ID] for parent in feature.parents):
                G.add_node(x,source='ann')
        #add transferred node
        if parts[11] == 'CDS':
            ratt_attributes = { a.split('=')[0]: a.split('=')[1].split(',') for a in parts[17].split(';')}
            ratt_ID = ratt_attributes['ID'][0]
            for x in (parent for feature in ratt[ratt_ID] for parent in feature.parents):
                G.add_node(x,source='ratt')
        #if CDS and on the same strand, add edge
        if parts[2] == 'CDS' and parts[11] == 'CDS' and parts[6] == parts[15]:
            for x in (parent for feature in ann[ann_ID] for parent in feature.parents):
                for y in (parent for feature in ratt[ratt_ID] for parent in feature.parents):
                    G.add_edge(x,y)
    return G

def get_parents(gff,name):
    """
    Recursive generator that returns all nested parent names
    """
    for sub in gff[name.ID]:
        yield sub
        for p in (y for x in sub.parents for y in gff[x]):
            for q in get_parents(gff,p):
                yield q

def get_gene(gff,name):
    """
    Returns the corresponding gene GffSubPart object for any feature
    """
    for sub in get_parents(gff,name):
        if sub.featuretype == 'gene':
            yield sub

def overlap(ann, ratt):
    """
    This uses gff_toolkit, advised (less dependencies, no external processes)
    
    Args:
    	ann (Gff):
    	ratt (Gff):
    Returns:
    	networkX graph object with all overlapping genemodels
    """
    print 'Determining {0} vs {1} overlap'.format(ann,ratt)
    overlaps = nx.Graph()

    for ann_gene in ann.getitems(featuretype='gene'):
        for ann_transcript in ann.get_children(ann_gene,featuretype='mRNA'):
            overlaps.add_node(ann_gene.ID,source='ann')
            for ann_cds in ann.get_children(ann_transcript,featuretype='CDS'):
                for ratt_cds in ratt.getitems(seqid=ann_cds.seqid,start=ann_cds.start,end=ann_cds.end,strand=ann_cds.strand,featuretype='CDS'):
                    for ratt_gene in get_gene(ratt,ratt_cds):
                        overlaps.add_node(ratt_gene.ID,source='ratt')
                        overlaps.add_edge(ann_gene.ID,ratt_gene.ID)
    for ratt_gene in ratt.getitems(featuretype='gene'):
        if ratt_gene.ID not in overlaps.nodes():
            overlaps.add_node(ratt_gene.ID,source='ratt')
    return overlaps

def main(config_file):
    with open(config_file,'rU') as fh:
        config = json.load(fh)

    ratt_file1 = config['ratt1']
    ratt_file2 = config['ratt2']
    ann_file1 = config['gff1']
    out_file1 = 'transfer.out.{0}'.format(ann_file1)
    singleton_file1 = 'transfer.singleton.{0}'.format(ann_file1)
    ann_file2 = config['gff2']
    out_file2 = 'transfer.out.{0}'.format(ann_file2)
    singleton_file2 = 'transfer.singleton.{0}'.format(ann_file2)
    fasta_file1 = config['fasta1']
    fasta_file2 = config['fasta2']
    bam_files1 = config['bam1']
    bam_files2 = config['bam2']

    print 'parsing',ratt_file1
    ratt1 = gt.parser(ratt_file1, fasta_file = fasta_file1)
    print 'parsing',ratt_file2
    ratt2 = gt.parser(ratt_file2, fasta_file = fasta_file2)
    print 'parsing',ann_file1
    ann1 = gt.parser(ann_file1, fasta_file = fasta_file1)
    print 'parsing',ann_file2
    ann2 = gt.parser(ann_file2, fasta_file = fasta_file2)

    bam1 = []
    for bam_file1 in bam_files1:
        print bam_file1
        bam1.append(pysam.AlignmentFile(bam_file1, 'rb'))
    bam2 = []
    for bam_file2 in bam_files2:
        print bam_file2
        bam2.append(pysam.AlignmentFile(bam_file2, 'rb'))

    anndic = {1: ann1, -1: ratt1, 2: ann2, -2: ratt2}
    bamdic = {1: bam1, -1: bam1, 2: bam2, -2: bam2}

    #this toggle is vital to switching between the two genomes and keeping track of all the right files
    toggle = {1: 2, 2: 1}

    seen = set()

    correct = {1: [], 2: []}

    total = {1: 0, 2: 0}

    G1 = overlap(anndic[1], anndic[-1])
    G2 = overlap(anndic[2], anndic[-2])

    graphdic = {1: G1, 2: G2}

    counter = 0

    clusters = []

    print 'Finding clusters'
    with open(out_file1,'w') as out1, open(out_file2,'w') as out2, open(singleton_file1,'w') as singleton1, open(singleton_file2,'w') as singleton2:
        outdic = {1: (out1, singleton1), -1: (out1, singleton1), 2: (out2, singleton2), -2: (out2, singleton2)}
        for t in toggle:
            anndic = {1: anndic[t], -1: anndic[-t], 2: anndic[toggle[t]], -2: anndic[-toggle[t]]}
            bamdic = {1: bamdic[t], -1: bamdic[-t], 2: bamdic[toggle[t]], -2: bamdic[-toggle[t]]}
            outdic = {1: outdic[t], -1: outdic[-t], 2: outdic[toggle[t]], -2: outdic[-toggle[t]]}
            print 'toggle',t,anndic
            for component in nx.connected_components(graphdic[t]):
                print '---'
                print 'component',component
                print [x for x in component if x in seen]
                if [x for x in component if x in seen]:
                    print 'continue'
                    continue
                counter += 1
                if counter % 1000 == 0:
                    print 'Processed {0} clusters'.format(counter)
                names = set()
                cluster = AnnotationCluster(anndic,bamdic,outdic)
                for gene in component:
                    if graphdic[t].node[gene]['source'] == 'ann':
                        names.add(gene)
                        cluster.add(gene, 1)
                        if gene in anndic[-2]:
                            cluster.add(gene, -2)
                print 'starting',cluster
                if cluster.size() == 0:
                    continue
                clustersize = 0
                #while loop to find all overlapping genes in both genomes
                while clustersize < cluster.size():
                    clustersize = cluster.size()
                    newnames = set()
                    for n1 in names:
                        if n1 in seen:
                            continue
                        seen.add(n1)
                        for n2 in graphdic[t][n1]:
                            cluster.add(n2, -1)
                            cluster.add(n2, 2)
                            seen.add(n2)
                            for newname in graphdic[toggle[t]][n2]:
                                cluster.add(newname, -2)
                                cluster.add(newname, 1)
                                newnames.add(newname)
                                seen.add(newname)
                    names = newnames
                print 'finished',cluster
                clusters.append(cluster)
                #score the overlapping genes based on transcriptomes and/or protein similarity
                if cluster.choose():
                    #write the best models to their respective files (singletons saved separately)
                    print 'SAVE',cluster.solution
                    cluster.save()

if __name__ == '__main__':
    if len(sys.argv) == 2:
        main(*sys.argv[1:])
    else:
        print 'Usage: python {0} <config.json>'.format(sys.argv[0])







