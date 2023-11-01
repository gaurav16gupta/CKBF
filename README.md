# IDentity with Locality: An ideal hash for efficient gene sequence search --- the official code repository

Please follow the instructions below for reproducing the COBS experiments. The code is tested on a machine running Ubuntu 20.04 with CMake 3.16.3 and GNU Make 4.2.1.

1) Download the data files using SRA toolkit. For installing SRA toolkit use this guide https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
2) Download the sequence files to `data/` folder: 
```
cd data
./download.sh
cd ..
```
3) Build the program:
```
mkdir build && cd build
cmake .. && make -j
cd ..
```
4) Run all experiments (on RAM and on disk) for vanilla COBS vs IDL-COBS, with the results written into `results` folder
```
mkdir results
bash run_cobs.sh
```