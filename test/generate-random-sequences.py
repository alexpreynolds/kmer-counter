#!/usr/bin/env python

import sys
import random

def random_dna_sequence(length):
    return ''.join(random.choice('ACTG') for _ in range(length))

n=int(sys.argv[1])
length=int(sys.argv[2])

for i in range(n):
    sys.stdout.write(">%d\n%s\n" % (i, random_dna_sequence(length)))
