#!/usr/bin/env python

"""
File:           Algorithm3_Dynamic_Programming.py
Author:         Shengyuan Wang
Date:           Apr 1, 2020

Description:    Coursera_Genomic_Algorithm_quiz3
"""

from Algorithm1_Naive_Exact import *


def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(y)+1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    fewest_edit = 999999
    for i in range(len(D)):
        if D[i][-1] < fewest_edit:
            fewest_edit = D[i][-1]

    return fewest_edit


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match


def smart_overlap_map(reads, k):
    olaps = {}
    result = {}
    for read in reads:
        for i in range(len(read)-k+1):
            if read[i:i+k] not in olaps:
                olaps[read[i:i+k]] = [read]
            else:
                olaps[read[i:i+k]].append(read)

    count = 0
    for read in reads:
        read_suffix = read[-k:]
        for possible_read in olaps[read_suffix]:
            if possible_read != read:
                olen = overlap(read, possible_read, k)
                if olen > 0:
                    count += 1
                    result[(read, possible_read)] = olen

    return result, count


if __name__ == '__main__':
    filename = 'chr1.GRCh38.excerpt.fasta'

    t = readGenome(filename)

    # Question1
    p = 'GCTGATCGATCGTACG'
    print(editDistance(t, p))

    # Question2
    p = 'GATTTACCAGATTGAG'
    print(editDistance(t, p))

    # Question3, 4
    # Test
    reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
    print(len(smart_overlap_map(reads, 5)))

    reads_filename = 'ERR266411_1.for_asm.fastq'
    reads, _ = readFastq(reads_filename)
    pairs_olen, pairs_count = smart_overlap_map(reads, 30)
    print("How many pairs of reads overlap:", pairs_count)

    reads_involved = []
    for key, value in pairs_olen:
        reads_involved.append(key)
    print("How many reads have a suffix involved in an overlap:", len(set(reads_involved)))



