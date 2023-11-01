# IDentity with Locality: An ideal hash for efficient gene sequence search --- the official code repository

Please follow the instructions below for reproducing the experiments. The code is tested on a machine running Ubuntu 20.04 with CMake 3.16.3 and GNU Make 4.2.1.

1) Download the data files using SRA toolkit. For installing SRA toolkit use this guide https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
2) Download the sequence files to `data/` folder: 
```
cd data
fastq-dump ERR337856 -O .
cd ..
```
3) Build the program:
```
mkdir build && cd build
cmake .. && make -j
cd ..
```
4) Run all experiments for vanilla Bloom filters vs IDL-Bloom filters, with the results written into `results` folder
```
mkdir results
bash run.sh
```
5) To run all experiments for vanilla COBS vs IDL-COBS, switch to the `abf` branch by using `git checkout abf` and follow the instructions in README.md to run the experiments.
