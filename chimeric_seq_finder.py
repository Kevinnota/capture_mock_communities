#!/usr/bin/env python

import argparse
from tqdm import tqdm
import pysam

parser = argparse.ArgumentParser(prog='mapping summary', description='This tool will look for chimeric sequences in a quick and dirty way. It examines local alignments to determine if a read is mapping partially to two or more reference sequences, with one part in the reverse direction and another part in the forward direction. It then prints a summary with three identity scores over the whole sample. Don't use this for detecting chimerics sequences in complex samples')
parser.add_argument('--bam', "-b",  help='input bam from Bowtie2 local alignment with -k>2')
parser.add_argument('--print', "-p", action='store_true', help='print a nice output for the terminal, keep false for loops and making tables', default=True)
args=parser.parse_args()

def reverse_cigar(cigar_string):
    cigar_parts = []
    current_part = ''
    for char in cigar_string:
        if char.isalpha():
            cigar_parts.append(current_part + char)
            current_part = ''
        else:
            current_part += char
    
    reversed_cigar_parts = cigar_parts[::-1]
    reversed_cigar = ''.join(reversed_cigar_parts)
    return reversed_cigar

def calculate_total_mapped_length(cigar_list, read_length):
    mapped_length = 0
    count=0
    mapping_dict = {}
    for pos in range(read_length):
        mapping_dict[pos]=0

    for cigar_string in cigar_list:
        if cigar_string is None:
            continue
        if "I" in cigar_string:
            continue
        if "D" in cigar_string:
            continue
        translation_list = translate_cigar_to_list(cigar_string)
        
        for pos in range(len(translation_list)):
            mapping_dict[pos]+=translation_list[pos]
    for pos in mapping_dict.keys():
        if mapping_dict[pos]!=0:
            count+=1

    return count/read_length

def translate_cigar_to_list(cigar_string):
    cigar_parts = []
    current_part = ''
    for char in cigar_string:
        if char.isalpha():
            cigar_parts.append(current_part + char)
            current_part = ''
        else:
            current_part += char
    
    translation = []
    for part in cigar_parts:
        if part.endswith('M'):
            translation.extend([1] * int(part[:-1]))  # Add '1's for mapped base pairs
        else:
            translation.extend([0] * int(part[:-1]))  # Add '0's for soft-clipped or unmapped base pairs

    return translation

def reading_in_bam():
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(args.bam, "rb", check_sq=False)
    pysam.set_verbosity(save)
    previous_read_name = str()
    first_read=True
    #seq_dict={}
    cigars = []
    coverage_list=[]
    i=0
    for read in samfile:
        if first_read == True:
            first_read=False
            previous_read_name=read.query_name
            length = read.query_length

        if previous_read_name != read.query_name:
            #print(cigars)
            coverage_list.append(calculate_total_mapped_length(cigars, length))
            previous_read_name=read.query_name
            length = read.query_length
            cigars = []

            i+=1
        
        if read.is_reverse:
            cigars.append(reverse_cigar(read.cigarstring))
        else:
            cigars.append(read.cigarstring)

    #print('proportion at 95%', sum(1 for value in coverage_list if value > 0.95)/len(coverage_list))
    #print('proportion at 85%', sum(1 for value in coverage_list if value > 0.85)/len(coverage_list))
    #print('proportion at 75%', sum(1 for value in coverage_list if value > 0.75)/len(coverage_list))
    print(str(sum(1 for value in coverage_list if value > 0.95)/len(coverage_list))+"\t"+str(sum(1 for value in coverage_list if value > 0.85)/len(coverage_list))+"\t"+str(sum(1 for value in coverage_list if value > 0.75)/len(coverage_list)))

            
        #seq_dict[read.reference_name] = [read.is_forward, read.reference_start, read.reference_end]


if __name__ == '__main__':
    #print('start main')
    reading_in_bam()
