import numpy as np
import sys

n = int(sys.argv[1])
poison = int(sys.argv[2])
in_filename = str(sys.argv[3])
out_filename = str(sys.argv[4])

lines = None
with open(in_filename, "r") as f:
    lines = f.readlines()

lines = lines[:n]
newlines = []
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

for l in lines:
    idxs = np.random.randint(0, len(l) - 1, (poison,))
    larr = list(l)
    for idx in idxs:
        larr[idx] = switch(l[idx])
    newlines.append(''.join(larr))

with open(out_filename, "w") as f:
    f.writelines(newlines)
