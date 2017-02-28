#!/usr/bin/env python
"""
Author: Rens Holmer
"""

__author__ = 'rensholmer'

import re
import os
import shutil
from subprocess import Popen,PIPE,STDOUT
import sys
import multiprocessing as mp

EvmUtils = 'EVM_r2012-06-25/EvmUtils/'


def format_evidence(line):
    """
    Changes column three of a gff3 line (annotation type) into gene/mRNA/exon/CDS because EVM only works with these types
    
    Args:
        line (str):
    Returns:
        String with reformatted gff3 line
    """
    parts = line.strip().split('\t')
    name = parts[2]

    illegal_types = ('match','expressed_sequence_match',
    	'protein_match','translated_nucleotide_match')

    if name in illegal_types:
        gene = '\t'.join(parts[0:2] + ['gene'] + parts[3:])
        mrna = '\t'.join(parts[0:2] + ['mRNA'] + parts[3:])
        new_line = [gene,mrna]
    elif name == 'match_part':
        exon = '\t'.join(parts[0:2] + ['exon'] + parts[3:])
        cds = '\t'.join(parts[0:2] + ['CDS'] + parts[3:])
        new_line = [exon,cds]
    else:
        raise Exception
    return '\n'.join(new_line) + '\n'

def split_files(maker_gff, prefix):
    """
    Args:
        maker_gff (str):
        prefix (str):
    Returns:
        Dictionary with filenames per annotation type
    """
    file_dic = {}
    file_name_dic = {}
    folder = '{0}.EVIDENCEMODELER.INFILES'.format(prefix)
    type_dic = {'blastn': 'transcript',
                'cdna2genome': 'transcript',
                'est2genome': 'transcript',
                'transdecoder': 'transcript',
                'maker': 'ab_initio',
                'genemark': 'ab_initio',
                'repeatmasker': 'repeat',
                'repeatrunner': 'repeat',
                'snap_masked': 'ab_initio',
                'augustus': 'ab_initio',
                'tblastx': 'transcript',
                'blastx': 'protein',
                'protein2genome': 'protein'}
    try:
        shutil.rmtree(folder)
        print '{0} exists, replacing'.format(folder)
        os.mkdir(folder)
    except:
        os.mkdir(folder)

    for suffix in ('transcript', 'protein', 'ab_initio', 'repeat'): 
        splitfile = '{0}/{1}.{2}.gff'.format(folder, maker_gff.split('.gff')[0], suffix)
        try:
            os.remove(splitfile)
            print '{0} exists, overwriting'.format(splitfile)
        except:
            pass
        file_dic[suffix] = open(splitfile,'a')
        file_name_dic[suffix] = os.path.abspath(splitfile)

    with open(maker_gff,'rU') as inhandle:
        for line in inhandle:
            if not line.strip():
                continue
            elif line[0] == '#':
                continue
            parts = line.strip().split()
            if parts[1] == '.':
                continue
            tool = parts[1]
            if 'blast' in tool:
                continue
            evidence_type = type_dic[tool]
            if evidence_type in file_dic:
                outfile = file_dic[evidence_type]
            else:
                continue
            if tool == 'maker':
                if parts[2] in ('gene','mRNA','exon','CDS'):
                    outfile.write(line)
            elif evidence_type == 'repeat':
                if 'Simple_repeat' not in parts[8] and 'Low_complexity' not in parts[8]:
                    outfile.write(line)
            else:
                outfile.write(format_evidence(line))
    for splitfile in file_dic.keys():
        file_dic[splitfile].close()
    return file_name_dic

def prep_evm(genome, file_name_dic, weights, prefix):
    """
    Args:
        genome (str):
        file_name_dic (dict):
        weights (str):
        prefix (str):
    Returns:
        List of EVM commands and list of partitions 
    """
    folder = '{0}.EVIDENCEMODELER.PARTITIONS'.format(prefix)
    partition_list = os.path.abspath('{0}/{1}.partitions.list'.format(folder, prefix))
    partition_log = os.path.abspath('{0}/{1}.partitions.log'.format(folder, prefix))
    evm_commands = []
    try:
        shutil.rmtree(folder)
        print '{0} exists, overwriting'.format(folder)
        os.mkdir(folder)
    except:
        os.mkdir(folder)
    cwd = os.getcwd()
    os.chdir(folder)
    

    partition = ('{0}partition_EVM_inputs.pl --genome {1} '
    			'--gene_predictions {2} --protein {3} '
    			'--transcript_alignments {4}  --repeats {5} '
    			'--segmentSize 150000 --overlapSize 80000 '
    			'--partition_listing {6}')
    
    partition.format(EvmUtils, genome, file_name_dic['ab_initio'], file_name_dic['protein'], 
    	file_name_dic['transcript'], file_name_dic['repeat'], partition_list)
   
    print partition
    with open(partition_log,'w') as log:
        p1 = Popen(partition.split(),stdout=log,stderr=log)
    stderr = p1.communicate()

    write_commands = ('{0}write_EVM_commands.pl --genome {1} '
    				'--weights {2} --gene_predictions {3} --protein {4} '
    				'--transcript_alignments {5} --repeats {6} '
    				'--output_file_name {7}.evm.out --partitions {8}')
    
    write_commands.format(EvmUtils,genome,weights,file_name_dic['ab_initio'],
    	file_name_dic['protein'],file_name_dic['transcript'],file_name_dic['repeat'],
    	prefix,partition_list)
    

    print write_commands
    p2 = Popen(write_commands.split(), stdout = PIPE, stderr = PIPE)
    stdout2,stderr2 = p2.communicate()
    if stderr2:
        print stderr2
    for line in stdout2.splitlines():
        if not line.strip():
            continue
        evm_commands.append(line.strip().split())

    os.chdir(cwd)
    return evm_commands,partition_list

def run_evm(command):
    """
    Args: 
        command (list):
    Returns:
        None
    """
    outfile = command.pop()
    command.pop()
    with open(outfile,'w') as outhandle:
        p = Popen(command, stdout = outhandle, stderr = PIPE)
        stderr = p.communicate()
        if stderr:
            return ['ERROR -- ',stderr]
        else:
            return ['SUCCES -- ',command]

def merge_evm_output(genome, prefix, partition_list):
    """
    Args:
        genome (str):
        prefix (str):
        partition_list (str):
    Returns:
        None
    """
    gff_files = []
    evm_output = '{0}.EVIDENCEMODELER.gff3'.format(prefix)
   
    recombine = ('{0}recombine_EVM_partial_outputs.pl --partitions {1} '
    	'--output_file_name {2}.evm.out')
    
    recombine.format(EvmUtils, partition_list, prefix)
    
    convert = ('{0}convert_EVM_outputs_to_GFF3.pl --partitions {1} '
    '--output_file_name {2}.evm.out --genome {3}')
    
    convert.format(EvmUtils, partition_list, prefix, genome)    
    
    print recombine
    p1 = Popen(recombine.split(), stdout = PIPE, stderr = STDOUT)
    stdout1,stderr1 = p1.communicate()
    for line in stdout1.splitlines():
        if not line.strip():
            continue
        parts = line.split()
        if parts[0] == 'writing':
            gff = parts[-1]+'.gff3'
            gff_files.append(gff)
    print convert
    p2 = Popen(convert.split(), stdout = PIPE, stderr = PIPE)
    stdout2,stderr2 = p2.communicate()

    try:
        os.remove(evm_output)
        print '{0} exists, overwriting'.format(evm_output)
    except:
        pass
    with open(evm_output,'w') as outhandle:
        for gff_file in gff_files:
            with open(gff_file,'rU') as inhandle:
                for line in inhandle:
                    if line.strip():
                        outhandle.write(line)

def main(genome, maker_gff, prefix, weights, threads):
    genome = os.path.abspath(genome)
    weights = os.path.abspath(weights)
    file_name_dic = split_files(maker_gff, prefix)
    commands,partition_list = prep_evm(genome, file_name_dic, weights, prefix)
    pool = mp.Pool(processes = int(threads))
    out = pool.map(run_evm, commands)
    merge_evm_output(genome, prefix, partition_list)


if __name__ == '__main__':
    if len(sys.argv) == 6:
        main(*sys.argv[1:])
    else:
        print 'Usage: python {0} <genome fasta> <maker gff3> <prefix> <weightsfile> <threads>'.format(sys.argv[0])
