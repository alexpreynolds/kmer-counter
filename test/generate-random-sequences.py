#!/usr/bin/env python

import sys
import random

def random_dna_sequence(length, seed, i):
    if seed:
        random.seed(seed + i)
    return ''.join(random.choice('ACTG') for _ in range(length))

n=int(sys.argv[1])
length=int(sys.argv[2])

try:
    seed=int(sys.argv[3])
except IndexError as err:
    seed=None

for i in range(n):
    sys.stdout.write(">%d\n%s\n" % (i, random_dna_sequence(length, seed, i)))
