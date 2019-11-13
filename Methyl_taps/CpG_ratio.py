import pandas as pd
import numpy as np
import os
import sys

import collections



def Get_mC_ratio(cpg_file, cpg_file_out):
    if os.stat(cpg_file).st_size == 0:
	print "empty"
	with open(cpg_file_out, 'w') as filout:
                filout.write("")
	return
    cpg_f = pd.read_table(cpg_file, header=None, index_col=None, sep='_')
    cpg_f['index'] = ["%s_%s" % (cpg_f.loc[mm,0], cpg_f.loc[mm,1]) for mm in cpg_f.index]
    cpg_dict = collections.defaultdict(list)

    for row in cpg_f.itertuples(index=True, name='Pandas'):
        site = row[-1]
        if cpg_dict[site] == []:
	    cpg_dict[site].extend(row[1:7])
	    cpg_dict[site].extend(row[8:-1])
        else:
	    cpg_dict[site][6] += row[8]
	    cpg_dict[site][7] += row[9]
	    cpg_dict[site][8] += row[10]
#    print len(cpg_dict)
 #   if len(cpg_dict) == 2:
#	print cpg_dict
    #print cpg_dict
    if cpg_dict:
    	cpg_dict_tmp = pd.DataFrame.from_dict(cpg_dict,orient='index')
    	cpg_dict_tmp['ratio'] = [1.0*cpg_dict_tmp.loc[mm,6] / cpg_dict_tmp.loc[mm,8] for mm in cpg_dict_tmp.index]
    	cpg_dict_tmp.to_csv(cpg_file_out, header=None,index=None,sep='\t')
    else:
	#cpg_dict_tmp = pd.DataFrame([])
	#cpg_dict_tmp.to_csv(cpg_file_out, header=None,index=None,sep='\t')
	with open(cpg_file_out, 'w') as filout:
		filout.write("")

	#print cpg_dict[site]
#chr1_133390_C_T_+_CGC_Both_2_0_2

#    hhh = ["chrom", "Pos", "Ref", "Strand", "Content", "Mod", "Unmod", "Depth", "Mod_ratio"]
#"chr17 82790490        C       +       CGT     1       0       1       1.0"
#		   0         1  2  3  4    5  6  7  8
#chr1_31373256  chr1  31373256  C  T  +  CGT  2  0  2


def get_sta(cpg_file_out, cpg_file_out_sta, sampleid, content):
    if os.stat(cpg_file_out).st_size == 0:
        print "empty"
        with open(cpg_file_out_sta, 'w') as filout:
                filout.write("")
        return
    cpg_dict_tmp = pd.read_table(cpg_file_out, header=None, index_col=None,sep ="\t")

    sta_all = {}

    sta_all['SampleID'] = sampleid
    sta_all['Sites_total_num'] = cpg_dict_tmp.shape[0]
    sta_all['methy_sites_num'] = cpg_dict_tmp.loc[cpg_dict_tmp[9]>0.0,:].shape[0]
    sta_all['methy_sites_ratio'] = 1.0*sta_all['methy_sites_num'] / sta_all['Sites_total_num']

    sta_all['all_depth'] = sum(cpg_dict_tmp[8])
    sta_all['methy_depth'] = sum(cpg_dict_tmp[6])
    sta_all['methy_depth_ratio'] = 1.0 * sta_all['methy_depth'] / sta_all['all_depth']


    if content == "CpG":
	kk = ['SampleID','Sites_total_num', 'methy_sites_num', 'methy_sites_ratio',
          'all_depth', 'methy_depth', 'methy_depth_ratio']
    else:
   	kk = ['SampleID','Sites_total_num', 'methy_sites_num', 'methy_sites_ratio',
   	  'CHG_sites_num', 'CHG_sites_ratio',
   	  'CHH_sites_num', 'CHH_sites_ratio',
   	  'all_depth', 'methy_depth', 'methy_depth_ratio'
   	 ]
	
        sta_all['CHG_sites_num'] = cpg_dict_tmp.loc[(cpg_dict_tmp[9]>0.0)
                               &((cpg_dict_tmp[5]=='CAG') | (cpg_dict_tmp[5]=='CTG') |(cpg_dict_tmp[5]=='CCG')),:].shape[0]
        sta_all['CHG_sites_ratio'] = 1.0*sta_all['CHG_sites_num']/sta_all['Sites_total_num']
        sta_all['CHH_sites_num'] = cpg_dict_tmp.loc[(cpg_dict_tmp[9]>0.0)
                               &(cpg_dict_tmp[5]!='CAG')
                               & (cpg_dict_tmp[5]!='CTG')
                               &(cpg_dict_tmp[5]!='CCG'),:].shape[0]
        sta_all['CHH_sites_ratio'] = 1.0*sta_all['CHH_sites_num']/sta_all['Sites_total_num']



    ttt = pd.DataFrame.from_dict(sta_all,orient='index').T
    ttt = ttt.loc[:,kk]
    ttt.to_csv(cpg_file_out_sta, header=True,index=None,sep='\t')


if __name__ == "__main__":
    cpg_file = sys.argv[1]
    cpg_file_out = sys.argv[2]
    cpg_file_out_sta = sys.argv[3]
    cpg_dict = Get_mC_ratio(cpg_file, cpg_file_out)

    base = os.path.basename(cpg_file_out_sta)
    name = os.path.splitext(base)[0]
 
    get_sta(cpg_file_out, cpg_file_out_sta, name,"CpG")
