import numpy as np
import sys

l = int(sys.argv[1])
filename = str(sys.argv[2])

idx = np.random.randint(0, 4, l*200)
base = np.array(['A', 'T' , 'C' , 'G'])
strings = base[idx]
randamarkers = np.cumsum(np.random.randint(35, 200, (l-1,)))
strings[randamarkers] = '\n'
with open(filename, "w") as f:
    f.writelines(strings)
