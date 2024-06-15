#!/usr/bin/env python

import pysam
import argparse
from tqdm import tqdm
import re
import statistics
from collections import Counter
import math
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

# this pyhthon script takes a bam produced by bowtie2 where reads and mapped to baits with --local(-sentisive) flag. It will find the largest ungabt part of the allignemnt to a bait and 

parser = argparse.ArgumentParser(prog='parse mapping data to baits', description='This Python script processes a BAM file produced by Bowtie2, where reads are mapped to baits using the --local (--sensitive) flag. It finds the largest ungapped part of the alignment to a bait and calculates the GC-content, ignoring the mismatches. Use the --summarizing_output flag if you want to calculate a median over multiple mappings when Bowtie2 is used with -k >1. The script will produce a median value for each read.')
parser.add_argument('--bam', "-b",  help='input bamfile', required=False)
parser.add_argument('--output', "-o",  help='file name', required=False)
parser.add_argument('--summarizing_output', "-so",  help='Do you want to summarize by calculating a median over multiple mappings (bowtie2 -k >1)', action='store_true', default=False)

args=parser.parse_args()

def shannon_entropy(seq):
    # Count occurrences of each nucleotide
    counts = Counter(seq)
    
    # Calculate probability for each nucleotide
    probabilities = [count / len(seq) for count in counts.values()]
    
    # Calculate Shannon entropy
    entropy = -sum(p * math.log2(p) for p in probabilities)
    
    return entropy


def read_bamfile():
    progress=tqdm(total=None, desc='total number of reads processed')
    write_file = open(args.output, "w")

    write_file.write("query_name"+"\t"\
              +"cigarstring"+"\t"\
              +"MD"+"\t"\
              +"NM"+"\t"\
              +"mismatch_count"+"\t"\
              +'query_sequence'+"\t"\
              +'query_sequence_mismatch'+"\t"\
              +'longest_match_sequence'+'\t' \
              +'longest_match_sequence_ref'+'\t' \
              +"gc_content_longest_match"+"\t"\
              +"gc_content_mismatch"+"\t"\
              +"gc_content_mapped"+"\n")
    
    lo=0
    LL=0
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(args.bam, "rb", check_sq=False, require_index=False)
    pysam.set_verbosity(save)
    for read in samfile:
        
        mismatch_count=0
        if read.is_unmapped == True:
            write_file.write(read.query_name+"\t"\
              +str("NA")+"\t"\
              +str("NA")+"\t"\
              +str("NA")+"\t"\
              +str("NA")+"\t"\
             +str("NA")+"\t"\
              +str("NA")+"\t"\
              +"NA"+"\t"\
              +"NA"+"\t"\
              +"NA"+"\t"\
              +"NA"+"\t"\
              +"NA"+"\n")

            continue
        md_tag = read.get_tag("MD")
        longest_match_length = 0
        longest_match_start = 0
        current_match_length = 0
        current_match_start = -1
        # Parse the CIGAR string to extract matching segments
        start_read=True
        start_read_ofset=0
        for op, length in read.cigartuples:
            if op == 4 and start_read == True:
                start_read=False
                start_read_ofset = length-1
        
            if op == 0:  # M (matching bases)
                if current_match_start == -1:
                    current_match_start = read.query_alignment_start  # Set the start position for the current match
                    first_match=current_match_start
                current_match_length += length

                if current_match_length > longest_match_length:
                    longest_match_length = current_match_length
                    longest_match_start = current_match_start
            
            else:  # Reset the current match when encountering a different operation
                current_match_start = -1
                current_match_length = 0     
            start_read=False
            
        sum_indels=0
        for op, length in read.cigartuples:
            if op == 1:
                sum_indels+=length
            if op == 0:  #M (matching bases)         
                if length != longest_match_length:
                    longest_match_start = longest_match_start+length

                if length == longest_match_length:
                    break
            

            #print(sum_indels)
            # Extract the longest matching segment for the current read
        #print(longest_match_sequence)
        if sum_indels == 0 and "I" not in read.cigarstring:
            longest_match_sequence = read.query_sequence[longest_match_start : longest_match_start+longest_match_length]
            longest_match_sequence_ref = read.get_reference_sequence()
        elif sum_indels == 0 and "I" in read.cigarstring:
            #print('longest_match_length', longest_match_length)
            longest_match_sequence = read.query_sequence[longest_match_start : longest_match_start+longest_match_length]
            
            if start_read_ofset == 0:
                longest_match_sequence_ref = read.get_reference_sequence()[0:longest_match_length]
            else:
                longest_match_sequence_ref = read.get_reference_sequence()[(longest_match_start-start_read_ofset-1):longest_match_length]
            #print("True")
        else :
            #print('test', 0-0)
            longest_match_sequence = read.query_sequence[longest_match_start + sum_indels : longest_match_start + longest_match_length + sum_indels]
            longest_match_sequence_ref = read.get_reference_sequence()[(longest_match_start-start_read_ofset-sum_indels):longest_match_start + longest_match_length - start_read_ofset-sum_indels]
            
        if not any(char.isalpha() for char in read.get_tag("MD")):
            mismatch_count=read.get_tag("NM")
        elif not any(letter in read.cigarstring for letter in ['I', 'D', 'S']):
            mismatch_count=read.get_tag("NM")
        elif read.cigarstring.count('M') == 1:
            mismatch_count=read.get_tag("NM")
        else:
            #print(read.cigarstring)
            #print(read.get_tag("MD"))
            #matches = re.findall(r"(\d+)([A-Z]*)", read.get_tag("MD"))
            matches = re.findall(r"(\d+)(\^[A-Z]+|[A-Z])?",  read.get_tag("MD"))
            # Convert the extracted matches into the desired format
            result_list = [(int(count), base) if base else (int(count), '') for count, base in matches]
            output_list = []
            cumulative_sum = first_match
            for i in range(len(result_list)):
                cumulative_sum += result_list[i][0]
                result_list[i] = (cumulative_sum, result_list[i][1])

            for item in result_list:
                if "^" in item:
                    continue
                elif item[1] == "" :
                    continue
                elif item[0] > longest_match_start and item[0] < longest_match_start+longest_match_length:
                    mismatch_count+=1
        
        regexp = re.compile(r'[a-z]')
        ref_seq=longest_match_sequence_ref
        read_seq=longest_match_sequence

        if regexp.search(ref_seq):
            for i in range(len(read_seq)-1):
                if regexp.search(ref_seq[i:i+1]):
                    ref_seq = ref_seq[:i] + "x" + ref_seq[i+1:]
                    read_seq = read_seq[:i] + "x" + read_seq[i+1:]

        ref_seq_2=read.get_reference_sequence()
        read_seq_2=read.query_alignment_sequence
        
        if regexp.search(read.get_reference_sequence()):
            for i in range(len(read.get_reference_sequence())-1):
                if regexp.search(read.get_reference_sequence()[i:i+1]):
                    ref_seq_2 = ref_seq_2[:i] + "x" + ref_seq_2[i+1:]
                    read_seq_2 = read_seq_2[:i] + "x" + read_seq_2[i+1:]
    
        sequence = Seq(read_seq)
        
        gc_count = longest_match_sequence.count('G') + longest_match_sequence.count('C')
        gc_content = gc_count / len(longest_match_sequence) * 100
        gc_count = read_seq.count('G') + read_seq.count('C')
        binding_potential = (read_seq.count('G') + read_seq.count('C')*3) + (read_seq.count('A') + read_seq.count('T')*2)
        gc_content_2 = gc_count / len(read_seq) * 100
        gc_count = read_seq_2.count('G') + read_seq_2.count('C')
        gc_content_3 = gc_count / len(read_seq_2) * 100
        write_file.write(read.query_name+"\t"\
              +str(read.cigarstring)+"\t"\
              +str(read.get_tag("MD"))+"\t"\
              +str(read.get_tag("NM"))+"\t"\
              +str(mismatch_count)+"\t"\
              +read.query_sequence+"\t"\
              +read_seq+"\t"\
              +longest_match_sequence+"\t"\
              +longest_match_sequence_ref+"\t"\
              +str(gc_content)+'\t'\
              +str(gc_content_2)+'\t'\
              +str(gc_content_3)+'\n')
        lo+=1
        progress.update()

def write_summery_table():

    in_file = open(args.output)
    
    summary_table = open(re.sub(".tsv", "_median_table.tsv", args.output), 'w')
    summary_table.write("\t".join(['query_name',
                "mismatch_count", 
                "longest_match_sequence",
                "NM", 
                "gc_content_longest_match",
                "gc_content_mismatch",
                "gc_content_mapped",
                "number_of_probes"])+"\n")

    merged_data = {}
    read_name=None
    i=0
    first=True
    for row in in_file:
        if 'query_name' in row:
            col_names=row.strip().split('\t')
            continue
        #print(read_name, row.split('\t')[0])
        if 'NA' in row:
            continue
        if read_name != row.split('\t')[0]:
            if first == True:
                first=False
                data={}
            elif first == False:    

                summary_table.write("\t".join([data['query_name'][0],
                str(1-(statistics.median([float(x) for x in data["mismatch_count"]]))), 
                str(statistics.median([float(x) for x in data["longest_match_sequence"]])), 
                str(statistics.median([float(x) for x in data["NM"]])),
                str(statistics.median([float(x) for x in data["gc_content_longest_match"]])),
                str(statistics.median([float(x) for x in data["gc_content_mismatch"]])),
                str(statistics.median([float(x) for x in data["gc_content_mapped"]])),
                str(len(data["gc_content_longest_match"]))])+'\n')
                
            #if i == 10:
            #    quit()
            data={}
            for item in col_names:
                data[item]=[] 
            read_name = row.split('\t')[0]
            i+=1
            
        col_data = dict(zip(col_names, row.strip().split('\t')))
        #print(col_data)
        #print(col_data.keys())
        for item in col_data.keys():
            if item == 'longest_match_sequence':
                data[item].append(len(col_data[item]))
            
            elif item == 'mismatch_count':
                data[item].append(len(col_data[item])/len(col_data['longest_match_sequence']))
            else:
                data[item].append(col_data[item])
                read_name = row.split('\t')[0]
        


if __name__ == '__main__' :
    #print("start main")
    if args.output != None:
        if not args.output.endswith('.tsv'):
            args.output=re.sub("\\..*", ".tsv", args.output)
 
    if args.output != None and args.summarizing_output== True:
        read_bamfile()
        write_summery_table()

    if args.output != None and args.summarizing_output== False:

        read_bamfile()

    if args.output == None and args.summarizing_output== True:
        args.output = args.summarizing_in
        write_summery_table()
