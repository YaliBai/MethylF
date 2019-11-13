#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import sys
import os



def get_cigar_list (cigarstring):
    '''Get the cigar list.
    cigarstring:'67M', '52M1I13M'. -> [(67, M)], [(52,M),(1,I),(13,M)]
    '''
    values = re.split(r'[SMDIH]', cigarstring)[:-1]
    types = re.split(r'\d+', cigarstring)[1:]
    cigar = [(v,d) for v,d in zip(values, types)]
    return cigar

def get_MD_list (MDstring):
    '''
    Get the MD tag list. MD_tag: MD:Z:12G29G24
    MDstring: 12G29G24
    MD list: [[12,'G'], [29,'G']]
    '''
    pat = "[0-9]+[GCAT]+"
    MD_s_list = re.findall(pat,MDstring)
    MD_list = [re.split(r'(\d+)',aa)[1:] for aa in MD_s_list]
    return MD_list

def formatSeqByCigar_I(seq, cigarString, MD_list):
    ''' Get cigarString.
    seq = query seq
    cigarString = '52M1I13M'
    cigar = [(52,M),(1,I),(13,M)]
    MD_list = [[12,'G'], [29,'G']]
    formatSeq = 'ATTTCCC' <which remove insert, add deletion, get mapped seq>
    '''
    formatSeq = ''
    pointer = 0; qstart = 0; qend = -1; origin_seq_len = 0
    cigar = get_cigar_list (cigarString)

    if cigar[0][1] == "S":
        qstart = int(cigar[0][0])
    if cigar[-1][1] == "S":
        qend = - int(cigar[-1][0]) - 1

    pointer = int(qstart)
    for pair in cigar:
        operation = pair[1]
        cigar_len = pair[0]
        if operation == "M":
            formatSeq += seq[pointer:(pointer+int(cigar_len))]
            tmp_length = 0
            flag = 0
            for pp, ss in MD_list:
                if (int(pp) + tmp_length) < int(cigar_len):
                    tmp_length = int(pp) + tmp_length + 1
                    flag += 1
                else:
                    break
            while flag > 0:
                MD_list.pop(0)
                flag -= 1

            pointer += int(cigar_len)
            origin_seq_len += int(cigar_len)

        elif operation == "I":
            pointer += int(cigar_len)
            origin_seq_len += int(cigar_len)

        elif operation == "D":
            formatSeq += 'D'*int(cigar_len)
#        elif operation == "N":
       #     formatSeq += 'N'*int(cigar_len)
#            print formatSeq,cigarString
#            breakpoint.append(len(formatSeq))
#            breakpoint.append(len(formatSeq)+int(cigar_len))
        elif operation == "H" or operation == "S":
            origin_seq_len += int(cigar_len)
            continue
        else:
            raise TypeError("There are cigar besides S/M/D/I/H\n")
    return formatSeq


def get_mis_site(MD_list, ref_length):
    ''' Get the mismatch sites position.
    MD_list: [[12,'G'], [29,'G']]
    ref_length: read reference seq length
    mismatch: [12, 42] <pos start from 0>
    '''
    ll = 0
    mismatch = []
    for pair in MD_list:
        ll += int(pair[0]) + 1
        mismatch.append(ll-1)
    return mismatch

def find_dinucleotide (seqin, di_n):
    '''
    find first pos of dinucleotides (such as CpG or GpC) in the seqin.
    di_n = "CG" or di_n = "CG" and so on.
    di_pos = [3, 13]
    '''
    sta = 0
    di_pos = []
    while sta < len(seqin):
        aa = seqin.find(di_n,sta)
        if aa == -1:
            return di_pos
            break
        di_pos.append(aa)
        sta=aa+1
    return di_pos


def get_mismatch_V(chrom, mismatch, read_aligned_pairs, formatSeq, read_type):
    '''
    mismatch = [12,29]
    read_aligned_pairs = [(0, 1052981, 'C'), (1, 1052982, 'T'), (2, 1052983, 'C')...]  get alt
    formatSeq = 'ATCTTTAACCGGG'
    read_type = 'read1' or 'read2'
    mismatch_sites = [[pos, ref, alt, read_type]]
    such as:    [[12,'G','A','read1'], [29,'G','A','read1']]
    '''
    mismatch_sites = []
    for mm in mismatch:
        if read_aligned_pairs[mm][2].upper() == 'A' or read_aligned_pairs[mm][2].upper() == 'T':
            next

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
    read_mismatchs = get_mismatch_V(chrom, read_mismatchs, read_aligned_pairs, formatSeq, read_type)
    return read_mismatchs

def strand_type(CpG_all_num, CT_al, GA_al, CT_mCpG_num, GA_mCpG_num):
    ''' Calculate the CT GA ratio to define strand type.
    Strand type includ:
        'C_strand', mainly happen C->T conversion
        'G_strand', mainly happen G->A conversion
        'P_strand', puzzled strand, can not define which conversion
        'N_strand', No conversion happed
     CpG_all_num = CpG_number in a seqin
     CpG_methy_ratio = mCpGs / all CpGs
    '''
    if CpG_all_num!=0:
        all_methy = CT_mCpG_num +GA_mCpG_num
        if all_methy > 0:
            CpG_methy_ratio = 1.0 * all_methy / CpG_all_num
            if 1.0 * CT_mCpG_num / all_methy > 0.6:
                read_strand = 'C_strand'
            elif 1.0 * CT_mCpG_num / all_methy < 0.4:
                read_strand = 'G_strand'
            else:
                read_strand = 'P_strand'
        else:
            all_change = CT_al + GA_al
            if all_change >0:
                if 1.0 * CT_al / all_change > 0.6:
                    read_strand = 'C_strand'
                elif 1.0 * CT_al / all_change < 0.4:
                    read_strand = 'G_strand'
                else:
                    read_strand = 'P_strand'
            else:
                read_strand = 'N_strand'
            CpG_methy_ratio = 0
    else:
        all_change = CT_al + GA_al
        CpG_methy_ratio = 0
        if all_change >0:
            if 1.0 * CT_al / all_change > 0.6:
                read_strand = 'C_strand'
            elif 1.0 * CT_al / all_change < 0.4:
                read_strand = 'G_strand'
            else:
                read_strand = 'P_strand'
        else:
            read_strand = 'N_strand'
    return read_strand, CpG_methy_ratio
