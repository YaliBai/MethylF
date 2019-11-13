#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import os
import pysam
import Bio
import re
from multiprocess import Pool
import time
import sys, getopt
import argparse
import subprocess
import glob
from reverse_complement import reverse_complement
from soft_process import *
from collections import defaultdict
import gzip
import multiprocess

'''
The code to Call the modifications in the Bam/Sam files.
Written by Bailey.
The last update is 2019.06.11.
'''


### Sample info and input files
#bam_in = args.input_bam
#fasta_file = args.reference
#outdir = args.outdir
#sample_id = args.sample_id


#library = args.library
#context = args.context
#modified_method = args.modified_method
#skip_clip_overlap = args.skip_clip_overlap
##minMQ = args.minMQ
#minBQ = args.minBQ
#ignore_orphans = args.ignore_orphans
#threads = args.threads
### Output Arguments
#compress = args.compress

#print
#'''
#%s input file is %s, Reference file is %s, Outdir is %s.
#Library is %s, call the modification %s in %s context.
#With the Arguments:
#skip_clip_overlap: %s;
#minMQ: %s; minBQ: %s;
#ignore_orphans: %s;
#mark_ends: %s
#Start running with the threads of %s ...
#''' % (sample_id, bam_in, fasta_file, outdir,
#library, modified_method, context,
#skip_clip_overlap, minMQ, minBQ, ignore_orphans, mark_ends, threads)
#singC_f="$indir/record/${ii}.mismatch.singleRead.singleC.txt"

Read_record_header = [ "read_id", "chrom", "Flagment_sta", "Flagment_end",
"Flagment_length","r1_flag","r2_flag",
"CpG_all_num",
"CT_CpG_num", "GA_CpG_num", "All_mCpG_num",
"CpG_methy_ratio", "Pair_type", "read_strand",
"CpG_methy_info", "CXX"]
Read_record_header = '\t'.join(map(lambda x:str(x),Read_record_header))

#global args
#global baseCpG_file_out
#global baseCH_file_out
#global fastafile
#global base_need
#global base_noneed
#global final_record
#
global lock

lock = multiprocess.Lock()

def get_short_seq (chrom, pos, strand, is_CpG=True):
	''' Get the short seq of target pos.
	strand = '+' / '-'
	short_seq: XXX
	is_CpG = True, The short seq will get from positive strand. which represent the CpG contenct in the ref.
	'''
	if strand == '+':
		short_seq = fastafile.fetch(chrom, pos, pos+3).upper()
	elif strand == '-' and not is_CpG:
		if pos > 1:
			short_seq = fastafile.fetch(chrom, pos-2, pos+1).upper()
			short_seq = reverse_complement(short_seq)
		else:
			short_seq = "Na"
	elif strand == '-' and is_CpG:
		short_seq = fastafile.fetch(chrom, pos-1, pos+2).upper()
		#short_seq = reverse_complement(short_seq)
	else:
		short_seq = 'Na'
	return short_seq


def get_mismatch_V(chrom, mismatch, read_aligned_pairs, base_qualities, formatSeq, read_type):
	'''
	mismatch = [12,29]
	read_aligned_pairs = [(0, 1052981, 'C'), (1, 1052982, 'T'), (2, 1052983, 'C')...]  get alt
	reads1.query_alignment_qualities
	base_qualities = array('B', [28, 26, 29, 34, 33, 28, 32, 32, 29, 35, 34])
	formatSeq = 'ATCTTTAACCGGG'
	read_type = 'read1' or 'read2'
	mismatch_sites = [[pos, ref, alt, read_type]]
	such as:    [[12,'G','A','read1'], [29,'G','A','read1']]
	'''
	mismatch_sites = []
	for mm in mismatch:
		if read_aligned_pairs[mm][2].upper() == base_noneed[0] or read_aligned_pairs[mm][2].upper() == base_noneed[1]:
			continue
		if int(base_qualities[mm]) < minBQ:
		#	print base_qualities[mm], minBQ
			continue
		m_tmp = [read_aligned_pairs[mm][1], read_aligned_pairs[mm][2].upper()]
		m_tmp.append(formatSeq[mm].upper())
		m_tmp.append(read_type)
		mismatch_sites.append(m_tmp)
	return mismatch_sites


def get_sites(read, read_type, chrom):
	''' Merge functions and get mismatch_sites
	Final: read_mismatchs = [[12,'G','A','read1'], [29,'G','A','read1']]
	read_aligned_pairs = [(0, 1052981, 'C'), (1, 1052982, 'T'), (2, 1052983, 'C')...]  get alt
	read_mismatchs = [[12,'G','A','read1'], [29,'G','A','read1']]
	'''
	read_MDstring = read.get_tag('MD')
	read_MD_list = get_MD_list(read_MDstring)
	read_mismatchs = get_mis_site(read_MD_list, read.reference_length)

	read_cigarstring = read.cigarstring
	read_query_s = read.query_sequence
	read_aligned_pairs = read.get_aligned_pairs(with_seq=True,matches_only=True)
	formatSeq = formatSeqByCigar_I(read_query_s, read_cigarstring, read_MD_list)
	# change the formatSeq by read.query_alignment_sequence from cut soft clipped in orginal seq
#	formatSeq = read.query_alignment_sequence
	base_qualities = read.query_alignment_qualities
	read_mismatchs = get_mismatch_V(chrom, read_mismatchs, read_aligned_pairs, base_qualities, formatSeq, read_type)
	return read_mismatchs


def get_all_CXX(chrom, start , end):
	'''Get all C, G, CpG position
	'''
	short_seq = fastafile.fetch(chrom, start , end).upper()
	G_start = start
	C_end = end
	if start != 0:
		G_start = start - 1
	if end != fastafile.get_reference_length(chrom)-1:
		C_end = end + 1
	CpG_seq_C = fastafile.fetch(chrom, start , C_end).upper()
	CpG_seq_G = fastafile.fetch(chrom, G_start, end).upper()
	all_C = find_dinucleotide(short_seq, 'C')
	all_G = find_dinucleotide(short_seq, 'G')
	all_CpG_C = find_dinucleotide(CpG_seq_C, 'CG')
	all_CpG_G = find_dinucleotide(CpG_seq_G, 'CG')

	all_C_noCpG = []
	all_G_noCpG = []

	CpGX = []
	if all_C:
		all_C_noCpG = list(set(all_C) - set(all_CpG_C))
		all_C_noCpG = map(lambda x: x + start, all_C_noCpG)
	if all_G:
		if start != 0:
			all_G_noCpG = list(set(all_G) - (set(all_CpG_G)))
		else:
			all_CpG_G = map(lambda x: x + 1, all_CpG_G)
			all_G_noCpG = list(set(all_G) - (set(all_CpG_G)))
			#print "get c,cg:", all_G, all_CpG_G, all_CpG_G
		all_G_noCpG = map(lambda x: x + start, all_G_noCpG)

	if all_CpG_C:
		all_CpG_C = map(lambda x: x + start, all_CpG_C)
	else:
		all_CpG_C = []

	if all_CpG_G:
		all_CpG_G = map(lambda x: x + start, all_CpG_G)
	else:
		all_CpG_G = []
	return all_C_noCpG, all_G_noCpG, all_CpG_C, all_CpG_G

def Get_CXX_methy (chrom, read_strand, CXX, GXX, All_mCXX_mis):
	'''Get all non CpG <C or G> mismatch list
	CXX_all_short = [chrom_pos_ref_alt_strand_shortseq_paired, ]
	'''
	CXX_all_short = []
	if read_strand == 'C_strand':
		if len(CXX) >0 :
			for ii in CXX:
				#short_seq = fastafile.fetch(chrom, ii, ii+2).upper()
				#short_seq = get_CXX(chrom, ii, 'C')
				ref = '*'
				alt = '*'
				strand = '+'
				pairend = '*'
				short_seq = get_short_seq (chrom, ii, strand, is_CpG=False)
				if len(All_mCXX_mis) > 0:
					if ii in list(All_mCXX_mis[0]):
						ref = 'C'
						alt = All_mCXX_mis.loc[All_mCXX_mis[0] == ii, 2].values[0]
						strand = '+'
						pairend = All_mCXX_mis.loc[All_mCXX_mis[0] == ii, 4].values[0]
				tmp_CpG = '_'.join(map(lambda x:str(x),[chrom,ii,ref, alt, strand, short_seq,pairend]))
				CXX_all_short.append(tmp_CpG)
	elif read_strand == 'G_strand':
		if len(GXX) >0 :
			for ii in GXX:
				#short_seq = fastafile.fetch(chrom, ii-1, ii+1).upper()
				#get_CXX (chrom, ii, 'G')
				ref = '*'
				alt = '*'
				strand = '*'
				pairend = '*'
				short_seq = get_short_seq (chrom, ii, strand, is_CpG=False)
				if len(All_mCXX_mis) > 0:
					if ii in list(All_mCXX_mis[0]):
						ref = 'G'
						alt = All_mCXX_mis.loc[All_mCXX_mis[0] == ii, 2].values[0]
						strand = '-'
						pairend = All_mCXX_mis.loc[All_mCXX_mis[0] == ii, 4].values[0]
				tmp_CpG = '_'.join(map(lambda x:str(x),[chrom,ii,ref, alt, strand, short_seq,pairend]))
				CXX_all_short.append(tmp_CpG)
	else:
		ref = '*'
		alt = '*'
		pairend = '*'
		if len(CXX) >0 :
			strand = '+'
			for ii in CXX:
				#short_seq = get_CXX (chrom, ii, 'C')
				short_seq = get_short_seq (chrom, ii, strand)
				tmp_CpG = '_'.join(map(lambda x:str(x),[chrom,ii,ref, alt, strand, short_seq,pairend]))
				CXX_all_short.append(tmp_CpG)
		if len(GXX) >0 :
			strand = '-'
			for ii in GXX:
				#short_seq = get_CXX (chrom, ii, 'G')
				short_seq = get_short_seq (chrom, ii, strand)
				tmp_CpG = '_'.join(map(lambda x:str(x),[chrom,ii,ref, alt, strand, short_seq,pairend]))
				CXX_all_short.append(tmp_CpG)
	CXX_all_short = ','.join(CXX_all_short)
	if CXX_all_short == '':
		CXX_all_short = 'Na'
	return CXX_all_short

def Get_all_sites_info (chrom, strand, sites_list, All_mis, depth, is_CpG=True):
	'''Get all non CpG <C or G> mismatch list
	strand = '+/-'
	sites_list = [pos1, pos2]
	All_mis = [[101987, 'G', 'A', 'read2'], [101831, 'G', 'A', 'read1']]
	all_sites_info = [chrom_pos_ref_alt_strand_shortseq_paired, ]
	is_CpG = True, The G_pos will minus 1 to be equal to C pos in this CpG site. G_pos = G_pos - 1
	'''
	all_sites_info = []
	All_mis = pd.DataFrame(All_mis)
	if strand == '+':
		ref = 'C'
		alt = "T"
	elif strand == '-':
		ref = 'G'
		alt = "A"
	if len(sites_list) > 0:
		#if All_mis.shape[0] > 0:
			#print "  "All_mis[0].tolist()
		for site in sites_list:
			pairend = '*'
			mod = 0
			unmod = depth
			short_seq = get_short_seq(chrom, site, strand, is_CpG)


			if All_mis.shape[0] > 0:
				if site in All_mis[0].tolist():
					pairend = All_mis.loc[All_mis[0] == site, 3].values[0]
					if pairend == "dread1" or pairend == "dread2":
						mod = 1
						unmod = 1
					else:
	                                        mod = depth
						unmod = 0
					
			if is_CpG and strand == '-' and ref == "G":
				site = str(int(site) - 1)
			tmp_CpG = '_'.join(map(lambda x:str(x),[chrom, site, ref, alt, strand, short_seq, pairend, mod, unmod, depth]))
			all_sites_info.append(tmp_CpG)
	return all_sites_info


def write_out_file(file_out, results):
	#lock.acquire()
	fo = open(file_out, "a+")
	fo.write(results+'\n')
	fo.close()
	#lock.release()

def Write_out_new(file_out, final_record):
	#lock.acquire()
	handle_out = gzip.open(file_out, 'wt')
	handle_out.write(Read_record_header)
	for record in final_record.itervalues():
		handle_out.write(record)
	handle_out.close()
	#lock.release()

def base_write(base_file_out,base_record):
	#lock.acquire()
	handle_out = open(base_file_out, 'wt')
	handle_out.write(base_record)
	for record in final_record.itervalues():
		handle_out.write(record)
	handle_out.close()
	#lock.release()

def Single_read(read, chrom):
	global lock
	r1_sta = read.reference_start
	r1_end = read.reference_end
	r1_seq = read.get_reference_sequence
	r1_map_len = read.reference_length
	r1_query_len = read.query_length
	read_id = read.query_name
	pair_type = 'single_end'

	if args['clip_5_end']!=0 or args['clip_3_end']!=0:
		r1_end_gap = r1_query_len - read.query_alignment_end
		if read.is_reverse:
			if read.query_alignment_start< args['clip_3_end']:
				r1_sta = r1_sta + (args['clip_3_end'] - read.query_alignment_start)
			if r1_end_gap < args['clip_5_end']:
				r1_end = r1_end - (args['clip_5_end'] - r1_end_gap)

		else:
			if read.query_alignment_start < args['clip_5_end']:
				r1_sta = r1_sta + (args['clip_5_end'] - read.query_alignment_start)
			if r1_end_gap < args['clip_3_end']:
				r1_end = r1_end - (args['clip_3_end'] - r1_end_gap)

	if r1_end < r1_sta:
		return

	r_mismatchs = get_sites(read, 'read', chrom) # [12,'G','A','read1', mod, depth] get mismatchs in each read

	all_noCpG_C, all_noCpG_G, all_CpG_C, all_CpG_G = get_all_CXX(chrom, r1_sta, r1_end) # get ref C/G positions in each read
	CpG_all_num = 0.5 * (len(all_CpG_C) + len(all_CpG_G))

	'''get the mismatch type.  And get strand type.
	'''
	mis_CpG_C = []
	mis_CpG_G = []
	mis_noCpG_C = []
	mis_noCpG_G = []
	CT_al = 0
	GA_al = 0
	CT_CpG_num = 0
	GA_CpG_num = 0
	All_mCpG_num = 0
	depth = 1

	## Match mismatchs in the C/G positions
	if r_mismatchs != []:
		All_arrange_mis = pd.DataFrame(r_mismatchs)

		mis_CpG_C = [All_arrange_mis.loc[pos,:] for pos in range(All_arrange_mis.shape[0]) if All_arrange_mis.loc[pos,0] in all_CpG_C and All_arrange_mis.loc[pos,2] == 'T']
		mis_CpG_G = [All_arrange_mis.loc[pos,:] for pos in range(All_arrange_mis.shape[0]) if All_arrange_mis.loc[pos,0] in all_CpG_G and All_arrange_mis.loc[pos,2] == 'A']
		mis_noCpG_C = [All_arrange_mis.loc[pos,:] for pos in range(All_arrange_mis.shape[0]) if All_arrange_mis.loc[pos,0] in all_noCpG_C and All_arrange_mis.loc[pos,2] == 'T']
		mis_noCpG_G = [All_arrange_mis.loc[pos,:] for pos in range(All_arrange_mis.shape[0]) if All_arrange_mis.loc[pos,0] in all_noCpG_G and All_arrange_mis.loc[pos,2] == 'A']

		CT_CpG_num = len(mis_CpG_C)
		GA_CpG_num = len(mis_CpG_G)
		CT_al = CT_CpG_num + len(mis_noCpG_C)
		GA_al = GA_CpG_num + len(mis_noCpG_G)
		All_mCpG_num = CT_CpG_num + GA_CpG_num

		strand_t, CpG_m_ratio = strand_type(CpG_all_num, CT_al, GA_al, CT_CpG_num, GA_CpG_num)
	else:
		strand_t, CpG_m_ratio = ['N_strand', 0.0]
	if strand_t == 'P_strand':
		return


	All_CpG_m = []
	All_noCpG_m = []
	if strand_t == 'C_strand':
		all_CpG_C_info = Get_all_sites_info(chrom, '+', all_CpG_C, mis_CpG_C, depth, is_CpG=True)
		all_noCpG_C_info = Get_all_sites_info(chrom, '+', all_noCpG_C, mis_noCpG_C, depth, is_CpG=False)
		All_CpG_m.extend(all_CpG_C_info)
		All_noCpG_m.extend(all_noCpG_C_info)
	elif strand_t == 'G_strand':
		all_CpG_G_info = Get_all_sites_info(chrom, '-', all_CpG_G, mis_CpG_G, depth, is_CpG=True)
		all_noCpG_G_info = Get_all_sites_info(chrom, '-', all_noCpG_G, mis_noCpG_G, depth, is_CpG=False)
		All_CpG_m.extend(all_CpG_G_info)
		All_noCpG_m.extend(all_noCpG_G_info)
	elif strand_t == 'N_strand' or strand_t == 'P_strand':
		all_CpG_C_info = Get_all_sites_info(chrom, '+', all_CpG_C, mis_CpG_C, depth, is_CpG=True)
		all_noCpG_C_info = Get_all_sites_info(chrom, '+', all_noCpG_C, mis_noCpG_C, depth, is_CpG=False)
		All_CpG_m.extend(all_CpG_C_info)
		All_noCpG_m.extend(all_noCpG_C_info)

	if All_CpG_m != []:
		lock.acquire()
		handle_out = open(baseCpG_file_out, 'a+')
		for ii,rr in enumerate(All_CpG_m):
			handle_out.write(rr+'\n')
		handle_out.close()
		lock.release()
	if All_noCpG_m != []:
		lock.acquire()
		handle_out = open(baseCH_file_out, 'a+')
		for ii,rr in enumerate(All_noCpG_m):
			handle_out.write(rr+'\n')
		handle_out.close()
		lock.release()

	All_CpG_m = ','.join(All_CpG_m)
	All_noCpG_m = ','.join(All_noCpG_m)
	if All_CpG_m == '':
		All_CpG_m = 'Na'
	if All_noCpG_m == '':
		All_noCpG_m = 'Na'

	results = [read_id, chrom, r1_sta, r1_end, r1_end - r1_sta+1, read.flag, "Na",
			   CpG_all_num, CT_al, GA_al, CT_CpG_num, GA_CpG_num, All_mCpG_num,
			   CpG_m_ratio, pair_type, strand_t, All_CpG_m, All_noCpG_m]
	final_record[read_id] = results

def Pair_read(read1, read2, chrom):
	global lock
	''' Estimate the converion type of pair end reads, and get the converated C and non converated C from CpGs and nonCpGs.
	'''
	read_id = read1.query_name

	r1_sta = read1.reference_start
	r2_sta = read2.reference_start
	r1_end = read1.reference_end
	r2_end = read2.reference_end
	r1_seq = read1.get_reference_sequence().upper()
	r2_seq = read2.get_reference_sequence().upper()
	r1_map_len = read1.reference_length
	r2_map_len = read2.reference_length
	r1_query_len = read1.query_length
	r2_query_len = read2.query_length

#	print read_id, r1_sta, r1_end, r2_sta, r2_end
	# Get overlap region (over_sta, over_end), and fragment region (new_sta, new_end)
        if r1_sta <= r2_sta:
                new_sta = r1_sta
        elif r1_sta > r2_sta:
                new_sta = r2_sta
        if r1_end <= r2_end:
                new_end = r2_end
        elif r1_end > r2_end:
                new_end = r1_end

	if new_end - new_sta <= read_len:
		clip_3_end = 10
	elif new_end - new_sta <= read_len + 10:
		clip_3_end = new_end - new_sta - read_len
	else:
		clip_3_end = 0

	if args['clip_5_end']!=0 or clip_3_end!=0:
		r1_end_gap = r1_query_len - read1.query_alignment_end
		r2_end_gap = r2_query_len - read2.query_alignment_end
		if read1.is_reverse:
			if read1.query_alignment_start < clip_3_end:
				r1_sta = r1_sta + (clip_3_end - read1.query_alignment_start)
			if r1_end_gap < args['clip_5_end']:
				r1_end = r1_end - (args['clip_5_end'] - r1_end_gap)

			if read2.query_alignment_start < args['clip_5_end']:
				r2_sta = r2_sta + (args['clip_5_end'] - read2.query_alignment_start)
			if r2_end_gap < clip_3_end:
				r2_end = r2_end - (clip_3_end - r2_end_gap)
		else:
			if read1.query_alignment_start < args['clip_5_end']:
				r1_sta = r1_sta + (args['clip_5_end'] - read1.query_alignment_start)
			if r1_end_gap < clip_3_end:
				r1_end = r1_end - (clip_3_end - r1_end_gap)

			if read2.query_alignment_start < clip_3_end:
				r2_sta = r2_sta + (clip_3_end - read2.query_alignment_start)
			if r2_end_gap < args['clip_5_end']:
				r2_end = r2_end - (args['clip_5_end'] - r2_end_gap)
	#print read_id, r1_sta, r1_end, r2_sta, r2_end
	if r2_end < r2_sta and r1_end < r1_sta:
		return
	elif r1_end < r1_sta:
	#	print read_id, "single"
		Single_read(read2, chrom)
		return
	elif r2_end < r2_sta:
	#	print read_id, "single"
		Single_read(read1, chrom)
		return

	# Get overlap region (over_sta, over_end), and fragment region (new_sta, new_end)
	if r1_sta <= r2_sta:
		new_sta = r1_sta
		over_sta = r2_sta
	elif r1_sta > r2_sta:
		new_sta = r2_sta
		over_sta = r1_sta
	if r1_end <= r2_end:
		new_end = r2_end
		over_end = r1_end
	elif r1_end > r2_end:
		new_end = r1_end
		over_end = r2_end

	## Get mismatch sites in the reads. [[12,'G','A','read1'], [29,'G','A','read1']]
	r1_mismatchs = get_sites(read1, 'read1', chrom)
	r2_mismatchs = get_sites(read2, 'read2', chrom)

	all_noCpG_C = []
	all_noCpG_G = []
	all_CpG_C = []
	all_CpG_G = []


	if over_end - over_sta > 0:
		'''This is Overlap type, if skip_clip_overlap is true, will not consider overlap region'''
		pair_type = 'overlap_pair'
		All_arrange_mis = []
		r1_mismatchs = pd.DataFrame(r1_mismatchs)
		r2_mismatchs = pd.DataFrame(r2_mismatchs)

		if r1_mismatchs.shape[0]>0 and r2_mismatchs.shape[0]>0:
			'''Find the intersect site between r1 r2. Look if they are the same.
			if same merge them, else label it with dread1 or dread2'''
			inter_sites = list(set(list(r1_mismatchs[0])).intersection(set(list(r2_mismatchs[0]))))
			if not args['skip_clip_overlap']:
				for sites in inter_sites:
					r1_sites = r1_mismatchs.loc[r1_mismatchs[0]==sites,:].values[0]
					r2_sites = r2_mismatchs.loc[r2_mismatchs[0]==sites,:].values[0]
					if r1_sites[2] == r2_sites[2]:
						tmp_d = r1_sites.tolist()
						tmp_d[3] = 'Both'
						All_arrange_mis.append(tmp_d)
					else:
						r1_sites[3] = 'dread1'
						r2_sites[3] = 'dread2'
						All_arrange_mis.append(r1_sites.tolist())
						All_arrange_mis.append(r2_sites.tolist())

			diff_sites = list(set(list(set(list(r1_mismatchs[0])).union(set(list(r2_mismatchs[0]))))).difference(set(inter_sites)))
			for sites in diff_sites:
				if sites in list(r1_mismatchs[0]):
					r1_sites = r1_mismatchs.loc[r1_mismatchs[0]==sites,:].values[0].tolist()
					All_arrange_mis.append(r1_sites)
				elif sites in list(r2_mismatchs[0]):
					r2_sites = r2_mismatchs.loc[r2_mismatchs[0]==sites,:].values[0].tolist()
					All_arrange_mis.append(r2_sites)
		elif r1_mismatchs.shape[0]>0:
			All_arrange_mis = np.array(r1_mismatchs).tolist()
		elif r2_mismatchs.shape[0]>0:
			All_arrange_mis = np.array(r2_mismatchs).tolist()


		all_noCpG_C, all_noCpG_G, all_CpG_C, all_CpG_G = get_all_CXX(chrom, new_sta, new_end)
		CpG_all_num = 0.5 * (len(all_CpG_C) + len(all_CpG_G))

		'''get the mismatch type.  And get strand type.
		'''
		mis_CpG_C = []
		mis_CpG_G = []
		mis_noCpG_C = []
		mis_noCpG_G = []
		CT_al = 0
		GA_al = 0
		CT_CpG_num = 0
		GA_CpG_num = 0
		All_mCpG_num = 0
		if All_arrange_mis != []:
			All_arrange_mis = pd.DataFrame(All_arrange_mis)

			mis_CpG_C = [All_arrange_mis.loc[pos,:] for pos in range(All_arrange_mis.shape[0]) if All_arrange_mis.loc[pos,0] in all_CpG_C and All_arrange_mis.loc[pos,2] == 'T']
			mis_CpG_G = [All_arrange_mis.loc[pos,:] for pos in range(All_arrange_mis.shape[0]) if All_arrange_mis.loc[pos,0] in all_CpG_G and All_arrange_mis.loc[pos,2] == 'A']
			mis_noCpG_C = [All_arrange_mis.loc[pos,:] for pos in range(All_arrange_mis.shape[0]) if All_arrange_mis.loc[pos,0] in all_noCpG_C and All_arrange_mis.loc[pos,2] == 'T']
			mis_noCpG_G = [All_arrange_mis.loc[pos,:] for pos in range(All_arrange_mis.shape[0]) if All_arrange_mis.loc[pos,0] in all_noCpG_G and All_arrange_mis.loc[pos,2] == 'A']

			CT_CpG_num = len(mis_CpG_C)
			GA_CpG_num = len(mis_CpG_G)
			CT_al = CT_CpG_num + len(mis_noCpG_C)
			GA_al = GA_CpG_num + len(mis_noCpG_G)
			All_mCpG_num = CT_CpG_num + GA_CpG_num

			strand_t, CpG_m_ratio = strand_type(CpG_all_num, CT_al, GA_al, CT_CpG_num, GA_CpG_num)
		else:
			strand_t, CpG_m_ratio = ['N_strand', 0.0]

		if strand_t == 'P_strand':
			return
	
		suppl_noCpG_C = []
        	suppl_noCpG_G = []
        	suppl_CpG_C = []
	        suppl_CpG_G = []

		
		suppl_noCpG_C1, suppl_noCpG_G1, suppl_CpG_C1, suppl_CpG_G1 = get_all_CXX(chrom, new_sta, over_sta)
                suppl_noCpG_C2, suppl_noCpG_G2, suppl_CpG_C2, suppl_CpG_G2 = get_all_CXX(chrom, over_end, new_end)
                suppl_noCpG_C.extend(suppl_noCpG_C1)
                suppl_noCpG_C.extend(suppl_noCpG_C2)
                suppl_noCpG_G.extend(suppl_noCpG_G1)
                suppl_noCpG_G.extend(suppl_noCpG_G2)
                suppl_CpG_C.extend(suppl_CpG_C1)
                suppl_CpG_C.extend(suppl_CpG_C2)
                suppl_CpG_G.extend(suppl_CpG_G1)
                suppl_CpG_G.extend(suppl_CpG_G2)


		All_CpG_m = []
		All_noCpG_m = []
		depth = 1
		if strand_t == 'C_strand':
			all_CpG_C_info = Get_all_sites_info(chrom, '+', suppl_CpG_C, mis_CpG_C, depth, is_CpG=True)
			all_noCpG_C_info = Get_all_sites_info(chrom, '+', suppl_noCpG_C, mis_noCpG_C, depth, is_CpG=False)
			All_CpG_m.extend(all_CpG_C_info)
			All_noCpG_m.extend(all_noCpG_C_info)
		elif strand_t == 'G_strand':
			all_CpG_G_info = Get_all_sites_info(chrom, '-', suppl_CpG_G, mis_CpG_G, depth, is_CpG=True)
			all_noCpG_G_info = Get_all_sites_info(chrom, '-', suppl_noCpG_G, mis_noCpG_G, depth, is_CpG=False)
			All_CpG_m.extend(all_CpG_G_info)
			All_noCpG_m.extend(all_noCpG_G_info)
		elif strand_t == 'N_strand':
			all_CpG_C_info = Get_all_sites_info(chrom, '+', suppl_CpG_C, mis_CpG_C, depth, is_CpG=True)
			all_noCpG_C_info = Get_all_sites_info(chrom, '+', suppl_noCpG_C, mis_noCpG_C, depth, is_CpG=False)
			All_CpG_m.extend(all_CpG_C_info)
			All_noCpG_m.extend(all_noCpG_C_info)
			

		if not args['skip_clip_overlap']:
                        over_noCpG_C, over_noCpG_G, over_CpG_C, over_CpG_G = get_all_CXX(chrom, over_sta, over_end)
			depth = 1
			if strand_t == 'C_strand':
	                        over_CpG_C_info = Get_all_sites_info(chrom, '+', over_CpG_C, mis_CpG_C, depth, is_CpG=True)
	                        over_noCpG_C_info = Get_all_sites_info(chrom, '+', over_noCpG_C, mis_noCpG_C, depth, is_CpG=False)
	                        All_CpG_m.extend(over_CpG_C_info)
	                        All_noCpG_m.extend(over_noCpG_C_info)
	                elif strand_t == 'G_strand':
	                        over_CpG_G_info = Get_all_sites_info(chrom, '-', over_CpG_G, mis_CpG_G, depth, is_CpG=True)
	                        over_noCpG_G_info = Get_all_sites_info(chrom, '-', over_noCpG_G, mis_noCpG_G, depth, is_CpG=False)
	                        All_CpG_m.extend(over_CpG_G_info)
	                        All_noCpG_m.extend(over_noCpG_G_info)
	                elif strand_t == 'N_strand' or strand_t == 'P_strand':
				depth = 1
	                        over_CpG_C_info = Get_all_sites_info(chrom, '+', over_CpG_C, mis_CpG_C, depth, is_CpG=True)
	                        over_noCpG_C_info = Get_all_sites_info(chrom, '+', over_noCpG_C, mis_noCpG_C, depth, is_CpG=False)
	                        All_CpG_m.extend(over_CpG_C_info)
	                        All_noCpG_m.extend(over_noCpG_C_info)


		if All_CpG_m != []:
			lock.acquire()
			handle_out = open(baseCpG_file_out, 'a+')
			for ii,rr in enumerate(All_CpG_m):
				handle_out.write(rr+'\n')
			handle_out.close()
			lock.release()
		if All_noCpG_m != []:
			lock.acquire()
			handle_out = open(baseCH_file_out, 'a+')
			for ii,rr in enumerate(All_noCpG_m):
				handle_out.write(rr+'\n')
			handle_out.close()
			lock.release()

		All_CpG_m = ','.join(All_CpG_m)
		All_noCpG_m = ','.join(All_noCpG_m)
		if All_CpG_m == '':
			All_CpG_m = 'Na'
		if All_noCpG_m == '':
			All_noCpG_m = 'Na'

		results = [read_id, chrom, new_sta, new_end, new_end - new_sta+1, read1.flag, read2.flag,
				   CpG_all_num, CT_al, GA_al, CT_CpG_num, GA_CpG_num, All_mCpG_num,
				   CpG_m_ratio, pair_type, strand_t, All_CpG_m, All_noCpG_m]
		final_record[read_id] = results

	else:
		'''This is Gap type'''
		'''Merge the mismatchs from read1 and read2'''
		pair_type = 'gap_pair'
		All_arrange_mis = []
		All_arrange_mis.extend(r1_mismatchs)
		All_arrange_mis.extend(r2_mismatchs)

		All_arrange_mis = pd.DataFrame(All_arrange_mis)

		r1_noCpG_C, r1_noCpG_G, r1_CpG_C, r1_CpG_G = get_all_CXX(chrom, r1_sta, r1_end)
		r2_noCpG_C, r2_noCpG_G, r2_CpG_C, r2_CpG_G = get_all_CXX(chrom, r2_sta, r2_end)

		all_noCpG_C.extend(r1_noCpG_C)
		all_noCpG_C.extend(r2_noCpG_C)
		all_noCpG_G.extend(r1_noCpG_G)
		all_noCpG_G.extend(r2_noCpG_G)

		all_CpG_C.extend(r1_CpG_C)
		all_CpG_C.extend(r2_CpG_C)
		all_CpG_G.extend(r1_CpG_G)
		all_CpG_G.extend(r2_CpG_G)

		CpG_all_num = 0.5 * (len(all_CpG_C) + len(all_CpG_G))

		'''get the mismatch type.  And get strand type.
		'''
		mis_CpG_C = []
		mis_CpG_G = []
		mis_noCpG_C = []
		mis_noCpG_G = []
		CT_al = 0
		GA_al = 0
		CT_CpG_num = 0
		GA_CpG_num = 0
		All_mCpG_num = 0
		if All_arrange_mis.shape[0]>0:
			mis_CpG_C = [All_arrange_mis.loc[pos,:] for pos in range(All_arrange_mis.shape[0]) if All_arrange_mis.loc[pos,0] in all_CpG_C and All_arrange_mis.loc[pos,2] == 'T']
			mis_CpG_G = [All_arrange_mis.loc[pos,:] for pos in range(All_arrange_mis.shape[0]) if All_arrange_mis.loc[pos,0] in all_CpG_G and All_arrange_mis.loc[pos,2] == 'A']
			mis_noCpG_C = [All_arrange_mis.loc[pos,:] for pos in range(All_arrange_mis.shape[0]) if All_arrange_mis.loc[pos,0] in all_noCpG_C and All_arrange_mis.loc[pos,2] == 'T']
			mis_noCpG_G = [All_arrange_mis.loc[pos,:] for pos in range(All_arrange_mis.shape[0]) if All_arrange_mis.loc[pos,0] in all_noCpG_G and All_arrange_mis.loc[pos,2] == 'A']

			CT_CpG_num = len(mis_CpG_C)
			GA_CpG_num = len(mis_CpG_G)
			CT_al = CT_CpG_num + len(mis_noCpG_C)
			GA_al = GA_CpG_num + len(mis_noCpG_G)
			All_mCpG_num = CT_CpG_num + GA_CpG_num
			strand_t, CpG_m_ratio = strand_type(CpG_all_num, CT_al, GA_al, CT_CpG_num, GA_CpG_num)
		else:
			strand_t, CpG_m_ratio = ['N_strand', 0.0]
		if strand_t == 'P_strand':
			return

		All_CpG_m = []
		All_noCpG_m = []
		depth = 1
		if strand_t == 'C_strand':
			all_CpG_C_info = Get_all_sites_info(chrom, '+', all_CpG_C, mis_CpG_C, depth, is_CpG=True)
			all_noCpG_C_info = Get_all_sites_info(chrom, '+', all_noCpG_C, mis_noCpG_C, depth, is_CpG=False)
			All_CpG_m.extend(all_CpG_C_info)
			All_noCpG_m.extend(all_noCpG_C_info)
		elif strand_t == 'G_strand':
			all_CpG_G_info = Get_all_sites_info(chrom, '-', all_CpG_G, mis_CpG_G, depth, is_CpG=True)
			all_noCpG_G_info = Get_all_sites_info(chrom, '-', all_noCpG_G, mis_noCpG_G, depth, is_CpG=False)
			All_CpG_m.extend(all_CpG_G_info)
			All_noCpG_m.extend(all_noCpG_G_info)
		elif strand_t == 'N_strand' or strand_t == 'P_strand':
			all_CpG_C_info = Get_all_sites_info(chrom, '+', all_CpG_C, mis_CpG_C, depth, is_CpG=True)
			all_noCpG_C_info = Get_all_sites_info(chrom, '+', all_noCpG_C, mis_noCpG_C, depth, is_CpG=False)
			All_CpG_m.extend(all_CpG_C_info)
			All_noCpG_m.extend(all_noCpG_C_info)

		if All_CpG_m != []:
			lock.acquire()
			handle_out = open(baseCpG_file_out, 'a+')
			for ii,rr in enumerate(All_CpG_m):
				handle_out.write(rr+'\n')
			handle_out.close()
			lock.release()
		if All_noCpG_m != []:
			lock.acquire()
			handle_out = open(baseCH_file_out, 'a+')
			for ii,rr in enumerate(All_noCpG_m):
				handle_out.write(rr+'\n')
			handle_out.close()
			lock.release()

		All_CpG_m = ','.join(All_CpG_m)
		All_noCpG_m = ','.join(All_noCpG_m)
		if All_CpG_m == '':
			All_CpG_m = 'Na'
		if All_noCpG_m == '':
			All_noCpG_m = 'Na'

		results = [read_id, chrom, new_sta, new_end, new_end - new_sta+1,
				read1.flag, read2.flag,
				CpG_all_num, CT_al, GA_al, CT_CpG_num, GA_CpG_num, All_mCpG_num,
				CpG_m_ratio, pair_type, strand_t, All_CpG_m, All_noCpG_m]

		final_record[read_id] = results
#        print results

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]



def methy_call(bam_file, chrom):
	samfile = pysam.AlignmentFile(bam_file, "rb")
	read_name = set()
	for read in samfile.fetch(chrom, multiple_iterators=True):
		readID = read.query_name
		Flag = read.flag
		proper_pair = [99, 147, 83, 163]
		#if readID in read_name:
		#	continue
		#else:
		#	read_name.add(read.query_name)
		#print readID, Flag

		if (Flag & 4 > 0) & (Flag & 8 > 0):  #print Flag, 'both unmap'
			continue
		if (Flag&512 > 0) | (Flag&1024 > 0):  #print Flag, 'Not uniquely'
			continue
		if (Flag&2>0): # not proper map
			if ((Flag & 16 > 0) & (Flag & 32 > 0)) | ((Flag & 16 == 0) & (Flag & 32 == 0)): 	#print Flag, 'wrong orientation'
				if args['Dip_C']:
					if not read.mapping_quality < minMQ:
						Single_read(read, chrom)
				else:
					continue
			elif (Flag & 256 > 0) | (Flag & 2048 > 0): # supplmentary align + second align
				if args['Dip_C']:
					if not read.mapping_quality < minMQ:
						Single_read(read, chrom)
				else:
					continue
		elif not args['ignore_orphans']: #proper map
			if (Flag & 4 == 0) & (Flag & 8 > 0):  #print Flag, 'single map'
				if not read.mapping_quality < minMQ:
					Single_read(read, chrom)
					#print Flag, read
			elif ((Flag & 16 > 0) & (Flag & 32 > 0)) | ((Flag & 16 == 0) & (Flag & 32 == 0)):   #print Flag, 'wrong orientation'
				if args['Dip_C']:
					if not read.mapping_quality < minMQ:
						Single_read(read, chrom)
				else:
					continue
			elif (Flag & 256 > 0) | (Flag & 2048 > 0):
				if args['Dip_C']:
					if not read.mapping_quality < minMQ:
						Single_read(read, chrom)
				else:
					continue
			else:     #print ff, 'Mapped uniquely, but with wrong insert size'
				if args['Dip_C']:
					if not read.mapping_quality < minMQ:
						Single_read(read, chrom)
				elif args['ignore_wrong_insert']:
		#			print "ignore_wrong_insert"
					if not read.mapping_quality < minMQ:
                                                Single_read(read, chrom)
	
	for read1, read2 in read_pair_generator(samfile,chrom):
		#print "PE", read1.query_name
		#print type(args['minMQ']), type(read1.mapping_quality)
		if read1.mapping_quality >= minMQ:
			if read2.mapping_quality >= minMQ:
#				print read2.query_name, "pair"
				Pair_read(read1, read2, chrom)
		#		print "PE", read1.query_name
			else:
#				print read1.query_name, "single"
				Single_read(read1, chrom)
		elif read2.mapping_quality >= minMQ:
#			print read2.query_name, "single"
			Single_read(read2, chrom)

global final_record
final_record = defaultdict(list)
def main_run(vars_a):
    global lock
    global args
    global baseCpG_file_out
    global baseCH_file_out
    global fastafile
    global base_need
    global base_noneed
    chrom, file_out, baseCpG_file_out, baseCH_file_out, args = vars_a

    print chrom, file_out, baseCpG_file_out, baseCH_file_out, args
    global minMQ
    global minBQ
    global read_len
    minMQ = int(args['minMQ'])
    minBQ = int(args['minBQ'])
    read_len = int(args['read_len'])

    handle_out = open(baseCpG_file_out, 'wt')
    handle_out.close()
    handle_out = open(baseCH_file_out, 'wt')
    handle_out.close()


    reference = args['reference']
    fastafile = pysam.Fastafile(reference)
    if args['modified_method'] == 'mCtoT':
        base_need = ['C', 'G']
        base_noneed = ['A', 'T']

    methy_call(args['input_bam'], chrom)
    #Write_out_new(file_out, final_record)
    f_w = pd.DataFrame.from_dict(final_record, orient='index')
    f_w.to_csv(file_out,header=None, index=None, sep='\t')


if __name__ == "__main__":
	''' main to run
	'''
	#args = get_args(sys.argv[1:])
	chrom, file_out, final_record, args = get_args(sys.argv[1:])
	fastafile = pysam.Fastafile(args['reference'])
	main_run(chrom, file_out, final_record, args)
	#fasta_file = '/share/newdata4/baiyl/basic/database/human/hg19/Homo_sapiens_seq/WholeGenomeFasta/genome.fa'
#	samfile = pysam.AlignmentFile(bam_file, "rb")
#	main_run(chrom, file_out, args)
