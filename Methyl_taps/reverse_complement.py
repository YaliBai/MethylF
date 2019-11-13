#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import os
import sys

def reverse_complement(seqin):
    '''Get the reverse and complement seq of a seqin.
    AACTG --> CAGTT
    '''
    seqout = seqin.lower()
    seqout = [re.sub(r't',"A",tt) for tt in [re.sub(r'a',"T",aa) for aa in [re.sub(r'g',"C",gg) for gg in [re.sub(r'c',"G",cc) for cc in seqout]]]]
    seqout = reverse1(seqout)
    return seqout

def reverse1(s):
    '''Get reverse seq'''
    s2 = s[::-1]
    s3 = ''.join(s2)
    return s3

def reverse2(s):
    s2 = list(s)
    s2.reverse()
    s3 = ''.join(s2)
    return s3

if __name__ == "__main__":
    args = sys.argv
    seqin = args[1]
    reversed_seq = reverse_complement(seqin)
    print reversed_seq
