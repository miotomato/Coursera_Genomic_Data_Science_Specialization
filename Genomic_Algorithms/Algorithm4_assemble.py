#!/usr/bin/env python

"""
File:           Algorithm4_Assemble.py
Author:         Shengyuan Wang
Date:           Apr 2, 2020

Description:    Coursera_Genomic_Algorithm_quiz4
"""

import itertools
import operator
from Algorithm3_Dynamic_programming import readFastq, smart_overlap_map


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match


def revised_scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = []
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        shortest_sup.append(sup)  # found shorter superstring
    shortest_len = len(ss) * len(ss[0])
    for sup in shortest_sup:
        if len(sup) <= shortest_len:
            shortest_len = len(sup)
    shortest_sup = [sup for sup in shortest_sup if len(sup) == shortest_len]
    return list(set(shortest_sup))  # return shortest


def scs(ss):
    """ Returns shortest common superstring of given strings,
        assuming no string is a strict substring of another """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            #sup += ssperm[i+1][-(len(ssperm[i+1])-olen):]
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest


def pick_maximal_overlap(reads, k):
    """ Return a pair of reads from the list with a
        maximal suffix/prefix overlap >= k.  Returns
        overlap length 0 if there are no such overlaps."""
    reada, readb = None, None
    best_olen = 0
    for a, b in itertools.permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > best_olen:
            reada, readb = a, b
            best_olen = olen
    return reada, readb, best_olen


def smart_greedy_scs(reads, k):
    """ Greedy shortest-common-superstring merge.
        Repeat until no edges (overlaps of length >= k)
        remain. """
    pairs_olen, pairs_count = smart_overlap_map(reads, k)
    sorted_pairs_olen = sorted(pairs_olen.items(), key=operator.itemgetter(1), reverse=True)
    read_a, read_b, olen = sorted_pairs_olen[0][0][0], sorted_pairs_olen[0][0][1], sorted_pairs_olen[0][1]
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        pairs_olen, pairs_count = smart_overlap_map(reads, k)
        if pairs_olen != {}:
            sorted_pairs_olen = sorted(pairs_olen.items(), key=operator.itemgetter(1), reverse=True)
            read_a, read_b, olen = sorted_pairs_olen[0][0][0], sorted_pairs_olen[0][0][1], sorted_pairs_olen[0][1]
        else:
            read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return ''.join(reads)


if __name__ == '__main__':
    # Question 1, 2
    shortest_sup_list = revised_scs(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'])
    print("Length of the shortest common superstring:", len(shortest_sup_list[0]))
    print("How many different shortest common superstrings:", len(shortest_sup_list))

    # Question 3, 4
    reads_filename = 'ads1_week4_reads.fastq'
    fastq_reads, _ = readFastq(reads_filename)

    genome = smart_greedy_scs(fastq_reads, 10)
    print(len(genome))
    print(genome.count('A'))
    print(genome.count('T'))
