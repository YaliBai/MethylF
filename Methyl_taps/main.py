#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
#from mismatch_record import main_run
from mismatch_record_i import main_run
from get_args import get_args
from CpG_ratio import Get_mC_ratio
from CpG_ratio import get_sta

#from  multiprocessing import Pool
import multiprocess
from get_args import get_args
import pysam
import signal
import pandas as pd
import re
import time
import subprocess
from time import clock

Read_record_header = [ "read_id", "chrom", "Flagment_sta", "Flagment_end",
"Flagment_length","r1_flag",
"CpG_all_num",
"CT_CpG_num", "GA_CpG_num", "All_mCpG_num",
"CpG_methy_ratio", "Pair_type", "read_strand",
"CpG_methy_info", "CXX"]
Read_record_header = '\t'.join(map(lambda x:str(x),Read_record_header))


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def sta_f(table_in, cpg_file_out_sta, sampleid):
    cpg_dict_tmp = pd.read_table(table_in, header=None, index_col=None,sep ="\t")
    kk = ['SampleID', 'Sites_total_num', 'methy_sites_num', 'methy_sites_ratio',
      'CHG_sites_num', 'CHG_sites_ratio',
      'CHH_sites_num', 'CHH_sites_ratio',
      'all_depth', 'methy_depth', 'methy_depth_ratio'
     ]
    #print cpg_dict_tmp.head(n=2)
    sta_all = {}
    sta_all['SampleID'] = sampleid 
    sta_all['Sites_total_num'] = cpg_dict_tmp.shape[0]
    sta_all['methy_sites_num'] = cpg_dict_tmp.loc[cpg_dict_tmp[9]>0.0,:].shape[0]
    sta_all['methy_sites_ratio'] = 1.0*sta_all['methy_sites_num'] / sta_all['Sites_total_num']
    sta_all['CHG_sites_num'] = cpg_dict_tmp.loc[(cpg_dict_tmp[9]>0.0)
                           &((cpg_dict_tmp[5]=='CAG') | (cpg_dict_tmp[5]=='CTG') |(cpg_dict_tmp[5]=='CCG')),:].shape[0]
    sta_all['CHG_sites_ratio'] = 1.0*sta_all['CHG_sites_num']/sta_all['Sites_total_num']
    sta_all['CHH_sites_num'] = cpg_dict_tmp.loc[(cpg_dict_tmp[9]>0.0)
                           &(cpg_dict_tmp[5]!='CAG')
                           & (cpg_dict_tmp[5]!='CTG')
                           &(cpg_dict_tmp[5]!='CCG'),:].shape[0]
    sta_all['CHH_sites_ratio'] = 1.0*sta_all['CHH_sites_num']/sta_all['Sites_total_num']

    sta_all['all_depth'] = sum(cpg_dict_tmp[8])
    sta_all['methy_depth'] = sum(cpg_dict_tmp[6])
    sta_all['methy_depth_ratio'] = 1.0 * sta_all['methy_depth'] / sta_all['all_depth']


    ttt = pd.DataFrame.from_dict(sta_all,orient='index').T
    #print ttt
    ttt = ttt.loc[:,kk]
    ttt.to_csv(cpg_file_out_sta, header=True,index=None,sep='\t')


def run_all(args):
    global lock
    t0=clock()

    bam_file = args['input_bam']
    file_out = args['output_dir'] +"/"+ args['prefix']

    samfile = pysam.AlignmentFile(bam_file, "rb")
    assert samfile.check_index(),  "ErrorType: %s file does not have index file." % bam_file

    header = pd.DataFrame(samfile.header['SQ'])
    Chrs = header['SN']
    samfile.close()

    threads = args['threads']
    lock = multiprocess.Lock()
    pool = multiprocess.Pool(threads,init_worker)
    a_args=multiprocess.Manager().dict()
    a_args = args
    results_q = []
    tmp_file_out = []
    tmp_ReadOut = []
    for chrom in Chrs:
        if re.search(r'_',chrom):
            continue
        if chrom == "chrM":
            continue
	print chrom
        file_out_tmp = file_out + '.' + chrom + ".ReadOut.tmp"
        baseCH_file_out = file_out + '.' + chrom + ".CH.tmp"
        baseCpG_file_out = file_out + '.' + chrom + ".CpG.tmp"
	vars_a = [chrom, file_out_tmp, baseCpG_file_out, baseCH_file_out, args]
	#main_run(vars_a)
	result = pool.apply_async(main_run, (vars_a, ))
        results_q.append(result)
	tmp_file_out.append(baseCH_file_out)
	tmp_file_out.append(baseCpG_file_out)
	tmp_ReadOut.append(file_out_tmp)

    try:
        print "Waiting 1 seconds"
        time.sleep(1)
    except KeyboardInterrupt:
        print "Caught KeyboardInterrupt, terminating workers"
        pool.terminate()
        pool.join()
    else:
        print "Quitting normally"
        pool.close()
        pool.join()
        print 'Pool finished'

    for result in results_q:
        print(result.get())

    t1=clock()
    print "Call methy spend %s" % str(t1-t0)

    ReadOut = file_out + ".ReadOut.txt"
    Args_m = ['cat %s' % (" ".join(tmp_ReadOut))]
    subprocess.check_call(Args_m, stdout = open(ReadOut, 'w'),shell=True)
 
    run_r = []
    tmp_CpG_out = []
    tmp_CH_out = []
    base_out = args['output_dir'] +"/"+ args['prefix'] + '.base'
    pool = multiprocess.Pool(args['threads'],init_worker)
    for chrom in Chrs:
        if re.search(r'_',chrom):
            continue
        if chrom == "chrM":
            continue
        baseCH_file_out = file_out + '.' + chrom + ".CH.tmp"
        baseCpG_file_out = file_out + '.' + chrom + ".CpG.tmp"
        baseCH_file_out2 = file_out + '.' + chrom + ".CH.txt"
        baseCpG_file_out2 = file_out + '.' + chrom + ".CpG.txt"


        r1_t = pool.apply_async(Get_mC_ratio, args=(baseCH_file_out, baseCH_file_out2, ))
        r2_t = pool.apply_async(Get_mC_ratio, args=(baseCpG_file_out, baseCpG_file_out2, ))	
	run_r.extend([r1_t, r2_t])

	tmp_CpG_out.append(baseCpG_file_out2)
	tmp_CH_out.append(baseCH_file_out2)

	tmp_file_out.append(baseCpG_file_out2)
	tmp_file_out.append(baseCH_file_out2)
    try:
        print "Waiting 1 seconds"
        time.sleep(1)
    except KeyboardInterrupt:
        print "Caught KeyboardInterrupt, terminating workers"
        pool.terminate()
        pool.join()
    else:
        print "Quitting normally"
        pool.close()
        pool.join()
        print 'Pool finished'

    for result in run_r:
        print(result.get())

    t2=clock()
    print "merge and get ratio spend %s" % str(t2-t1)

    baseCpG_file_out = file_out + ".CpG.txt"
    baseCH_file_out = file_out + ".CH.txt"
    baseCpG_file_sta = file_out + ".CpG.sta.txt"
    baseCH_file_sta = file_out + ".CH.sta.txt"
    
    print tmp_CpG_out
    Args_m = ['cat %s' % (" ".join(tmp_CpG_out))]
    subprocess.check_call(Args_m, stdout = open(baseCpG_file_out, 'w'),shell=True)

    print tmp_CH_out
    Args_m = ['cat %s' % (" ".join(tmp_CH_out))]
    subprocess.check_call(Args_m, stdout = open(baseCH_file_out, 'w'),shell=True) 

    get_sta(baseCpG_file_out, baseCpG_file_sta, args['prefix'],"CpG")
    get_sta(baseCH_file_out, baseCH_file_sta, args['prefix'],"CH")

    for my_file in tmp_file_out:
	if os.path.exists(my_file):
    	    os.remove(my_file)    

    for my_file in tmp_ReadOut:
	if os.path.exists(my_file):
            os.remove(my_file)
    t3=clock()
    print 'Time spend all: %s' % str(t3-t0)

if __name__ == "__main__":
    ''' main to run
    '''
    #fasta_file = '/share/newdata4/baiyl/basic/database/human/hg19/Homo_sapiens_seq/WholeGenomeFasta/genome.fa'
    args = get_args(sys.argv[1:])
    fastafile = pysam.Fastafile(args['reference'])
    run_all(args)
   


