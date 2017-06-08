#!/usr/bin/env python

import sys
import subprocess

k = 6
fastaFile = 'kmer.fa'
result = {}
kmerCmd = '../kmer-counter --fasta --k=%d %s' % (k, fastaFile)
try:
    output = subprocess.check_output(kmerCmd, shell=True)
    for line in output.splitlines():
        (header, counts) = line.strip().split('\t')
        header = header[1:]
        kmers = dict((k,int(v)) for (k,v) in [d.split(':') for d in counts.split(' ')])
        result[header] = kmers
except subprocess.CalledProcessError as error:
    output = error

print result
