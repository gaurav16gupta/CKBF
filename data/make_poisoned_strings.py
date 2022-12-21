import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--numQ", default=500, type=int)
# parser.add_argument("--poison", default=1, type=int)
parser.add_argument("--RandomQlen", default=True, type=bool)
parser.add_argument("--in_filename", default='input.fastq', type=str)
# parser.add_argument("--out_filename", default='query.txt', type=str)
args = parser.parse_args()

def switch(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'G'
    if c == 'G':
        return 'C'
    if c == 'C':
        return 'A'
    if c == 'N':
        return 'X'

def genQueries(poison):
    n = args.numQ
    qlen = args.RandomQlen
    in_filename = 'fastqFiles/'+args.in_filename
    out_filename = 'queries/'+args.in_filename.split('.')[0] + 'query.p' +str(poison)

    lines = None
    with open(in_filename, "r") as f:
        lines = f.readlines()

    lines = lines[:n]
    newlines = []
    i =0   
    for l in lines:
        if ((i+3)%4==0):
            if args.RandomQlen:
                qst = np.random.randint(0, len(l) - 1-32)
                l = l[qst:]
            if poison==0:
                newlines.append(l)
            else:
                idxs = np.random.randint(0, len(l) - 1, (poison,))
                larr = list(l)
                for idx in idxs:
                    larr[idx] = switch(l[idx])
                newlines.append(''.join(larr))
        i+=1

    with open(out_filename, "w") as f:
        f.writelines(newlines)

for poison in [0, 1, 2, 10, 20]:
    genQueries(poison)

