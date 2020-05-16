#!/usr/bin/env python

"""
File:           Algorithm1_naive_exact.py
Author:         Shengyuan Wang
Date:           Mar 30, 2020

Description:    Coursera_Genomic_Algorithm_quiz1
"""

import matplotlib.pyplot as plt


def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


def naive_with_rc(p, t):
    p_rev = reverseComplement(p)
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if not match:
            match = True
            for j in range(len(p)):  # loop over characters
                if t[i + j] != p_rev[j]:  # compare characters
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences


def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        count_mismatch = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                count_mismatch += 1
        if count_mismatch <= 2:
            occurrences.append(i)  # all chars matched; record
    return occurrences


def phred33ToQ(qual):
    return ord(qual) - 33


def createHist(qualities):
    # Create a histogram of quality scores
    hist = [0]*len(qualities[0])
    for qual in qualities:
        for i in range(len(qual)):
            q = phred33ToQ(qual[i])
            hist[i] += q
    return hist


if __name__ == '__main__':
    genome_file = 'lambda_virus.fa'

    genome = readGenome(genome_file)
    # Question1
    print(len(naive_with_rc('AGGT', genome)))

    # Question2
    print(len(naive_with_rc('TTAA', genome)))

    # Question3
    print(min(naive_with_rc('ACTAAGT', genome)))

    # Question4
    print(min(naive_with_rc('AGTCGA', genome)))

    # Question5
    print(len(naive_2mm('TTCAAGCC', genome)))

    # Question6
    print(min(naive_2mm('AGGAGGTT', genome)))

    # Question7
    reads_file = 'ERR037900_1.first1000.fastq'
    _, quals = readFastq(path, reads_file)
    h = createHist(quals)
    plt.plot(range(len(h)), h)
    plt.show()
