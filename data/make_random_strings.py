import numpy as np
import sys

n = int(sys.argv[1])
l = int(sys.argv[2])
filename = str(sys.argv[3])

idx = np.random.randint(0, 4, (n,l))
base = np.array(['A', 'T' , 'C' , 'G'])
strings = base[idx]
lines = [''.join(strings[i])+'\n' for i in range(n)]
with open(filename, "w") as f:
    f.writelines(lines)


