import os
import sys
import pandas as pd
import gzip
import re
import gzip
import argparse
from collections import defaultdict
from collections import Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import Levenshtein
import difflib

import math
import numpy as np
import time

from Bio import SeqIO
import gzip
from time import clock
import gc

'''The code is used to seperate reads by different cell barcode. 
Update in 2019.11.01 by Bailey. \nCopyright (c) Bailey.
What's new:
    The pooled reads with some sample without barcode,so the Tn5 is started in 0base.
    mark these sample with "SpikeIn"
'''


parser = argparse.ArgumentParser(description="The code is used to seperate reads by different cell barcode. Update in 2019.07.01 by Bailey. \nCopyright (c) Bailey")

parser.add_argument("-R1", "--Read1_fastq", dest = 'Read1_fastq',
                    help="Read1 fastq file path.")
parser.add_argument("-R2", "--Read2_fastq",  dest = 'Read2_fastq',
                    help="Read2 fastq file path.")
parser.add_argument("-b","--barcode_seq_file",  dest = 'barcode_seq_file',
                    help="barcode seq ref file: barcode_id\tbarcode_seq.")

parser.add_argument("-o", "--outdir", dest = 'outdir',
                    help="The outdir of output files")

parser.add_argument("--min_dis",
                    type = int,
                    default=2,
                    help="min dist of (index+barcode seq) and seq (Default 3).")
parser.add_argument("--Tn5",
                    type = str,
                    default='AGATGTGTATAAGAGACAG',
                    help="The Tn5 seq after barcode(Default AGATGTGTATAAGAGACAG)..")

parser.add_argument("--Tn5_min_dis",
                    type = int,
                    default=5,
                    help="min dist of Tn5.")

parser.add_argument("--prefix",
                    type = str,
                    help="The out name ")


args = parser.parse_args()

Read1_fastq = args.Read1_fastq
Read2_fastq = args.Read2_fastq
barcode_seq_file = args.barcode_seq_file
outdir = args.outdir

Tn5_minSta=20 #0+14+7-1
Tn5_maxSta=27 #4+15+7+1

min_dis = args.min_dis
Tn5 = args.Tn5
bridge1 = "TCGTCGGCAGCGTC" #14bp
bridge2 = "GTCTCGTGGGCTCGG"#15bp

basename = args.prefix
if not args.prefix:
     basename =re.split(r'[._]R1', os.path.basename(Read1_fastq))[0]



def get_meta_barcode_seq(file_in):
    '''barcode file:
        "bar_1S  GATATG"
    or meta tags file:
        "META16-1        GGCACCGAAAA"
    '''
    barcode_seq = {}
    bac = pd.read_table(file_in, header=None)
    barcode_seq = bac.set_index(0).to_dict()[1]
    return barcode_seq



def choose_barcode(seq_record, barcode_seq):
    read_seq = seq_record.seq

    tn5_sta = 0
    tn5_mindis = 9999
    #for ss in range(Tn5_minSta, len(read_seq)-len(Tn5)+1):
    for ss in range(Tn5_minSta, Tn5_maxSta):
        cut_seq = read_seq[ss:ss+len(Tn5)]
        tn5_dis = Levenshtein.distance(str(cut_seq), Tn5)
        if tn5_dis < tn5_mindis:
            tn5_mindis = tn5_dis
            tn5_sta = ss
    print "Tn5:", tn5_mindis
    if tn5_mindis > 10:
        bar_final = "WrongTn5"
        return bar_final, seq_record

    bar_tmp = 'None'
    bar_tmp_score = 9999

    for bar in barcode_seq.keys():
        bar_seq = barcode_seq[bar]
        
        query_cut_seq = str(read_seq[tn5_sta-len(bar_seq):tn5_sta])
        bar_dis = Levenshtein.distance(query_cut_seq, bar_seq)
        if bar_dis < bar_tmp_score:
            bar_tmp_score = bar_dis
            bar_tmp = bar
    
        
    if bar_tmp_score <= min_dis+1:
        bar_final = bar_tmp
        new_seq = SeqRecord(read_seq[tn5_sta+len(Tn5):], id=seq_record.id, description =seq_record.description, name=seq_record.name)
        new_seq.letter_annotations["phred_quality"] = seq_record.letter_annotations['phred_quality'][tn5_sta+len(Tn5):]
        return bar_final, new_seq
    else:
        bar_final = "PuzzBar"
        return bar_final, seq_record
    
def main ():
 #   if not os.path.exists(outdir):
 #       os.makedirs(outdir)

    barcode_seq = get_meta_barcode_seq(barcode_seq_file)

    out_basename = os.path.abspath(outdir) + '/' + basename
    print out_basename

    global max_barcode
    global mean_barcode
    max_barcode = 0
    all_l = 0
    for kk,vv in barcode_seq.items():
        if len(vv)>max_barcode:
            max_barcode = len(vv)
        all_l += len(vv)
    mean_barcode = all_l / len(barcode_seq)
    
    t0=clock()
    
   
    
    handle1 = gzip.open(Read1_fastq, "r")
    
    handle2 = gzip.open(Read2_fastq, "r")
    records2 = SeqIO.parse(handle2, "fastq")
    
    #aa=0
    for seq_record1 in SeqIO.parse(handle1, "fastq"):
        seq_record2 = records2.next()
        old_id = seq_record1.id
     #   if aa<5:
     #       print old_id
     #   aa += 1

        
        meta_tag1, new_seq1 = choose_barcode(seq_record1, barcode_seq)
        meta_tag2, new_seq2 = choose_barcode(seq_record2, barcode_seq)
        merged_id = "%s:%s_%s" % (seq_record1.id, meta_tag1, meta_tag2)
        new_seq1.id = merged_id
        new_seq2.id = merged_id
        
        seq1_quality = ''.join([chr(cc+33) for cc in new_seq1.letter_annotations['phred_quality']])
        seq2_quality = ''.join([chr(cc+33) for cc in new_seq2.letter_annotations['phred_quality']])
        
        out1 = gzip.open('%s.%s.%s.R1.fastq.gz' % (out_basename, meta_tag1, meta_tag2), "a")
        out2 = gzip.open('%s.%s.%s.R2.fastq.gz' % (out_basename, meta_tag1, meta_tag2), "a")
        
        out1.write('@%s\n'%merged_id)
        out1.write(str(new_seq1.seq)+"\n")
        out1.write('+\n')
        out1.write(seq1_quality+"\n")
        
        out2.write('@%s\n'%merged_id)
        out2.write(str(new_seq2.seq)+"\n")
        out2.write('+\n')
        out2.write(seq2_quality+"\n")

        out1.close()
        out2.close()
    
    handle1.close()
    handle2.close()
        

    
    t1 = clock()
    print 'Time spend all: %s' % str(t1-t0)
    

if __name__ == '__main__':
    main()








    


