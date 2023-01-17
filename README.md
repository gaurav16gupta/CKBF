#01/16/2023
```
taskset -c 2 ./debug/program ./configs/fuzzy/fuzzy.268435456.4.12.1024.cfg
taskset -c 2 ./debug/program ./configs/murmur/murmur.268435456.4.cfg
```
|              | Insert |       |  Query |       |
|:------------:|:------:|:-----:|:------:|:-----:|
|              | murmur | fuzzy | murmur | fuzzy |
|     hash     |   45K  |  330K |   40K  |  338K |
| insert/check |   70K  |  52K  |   16K  |  12K  |
|      FPR     |    0   |  8e-6 |    0   |  8e-6 |


# CKBF
cache efficient kmer Bloom Filter

1) Download the data files using SRA toolkit. For installing SRA toolkit use this guide https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
2) Download any sequence file in data/ folder: 
fastq-dump SRR649944  -O .
fastq-dump SRR649954  -O .
3) Parameter: set Bloom filter size and k in main.cpp
4) compile : make
5) run: ./build/program 0


