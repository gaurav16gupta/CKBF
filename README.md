# CKBF
cache efficient kmer Bloom Filter

1) Download the data files using SRA toolkit. For installing SRA toolkit use this guide https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
2) Download any sequence file in data/ folder: 
fastq-dump SRR649944  -O .
fastq-dump SRR649954  -O .
3) Parameter: set Bloom filter size and k in main.cpp
4) compile : make
5) run: ./build/program 0


