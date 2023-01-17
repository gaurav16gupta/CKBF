```
./debug/program SRR649954 1
range: 2147483648 bits, or 256 Mbs
```
Murmur hash

|              | Insert |       |  Query |       |
|:------------:|:------:|:-----:|:------:|:-----:|
|              | Disk | RAM | Disk | RAM |
|     hash     |   52 sec  |  27 |   1711 us |  1654 us |
| insert/check |   added  |  added  |   3354 us|  2431 us |
|  # False Pos |    2   |  -    |    2   |   - |

```
./debug/program SRR649954 1
range: 3147483648 bits, or 375 Mbs
```
Murmur hash

|              | Insert |       |  Query |       |
|:------------:|:------:|:-----:|:------:|:-----:|
|              | Disk | RAM | Disk | RAM |
|     hash     |   37 sec  |  27 |   2482 us |  1933 us |
| insert/check |   added  |  added  |   4219 us|  2610 us |
|  # False Pos |    2   |  -    |    2   |   - |





# CKBF
cache efficient kmer Bloom Filter

1) Download the data files using SRA toolkit. For installing SRA toolkit use this guide https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
2) Download any sequence file in data/ folder: 
fastq-dump SRR649944  -O .
fastq-dump SRR649954  -O .
3) Parameter: set Bloom filter size and k in main.cpp
4) compile : make
5) run: ./build/program 0


