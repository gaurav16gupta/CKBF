
1) Download the data files using SRA toolkit. For installing SRA toolkit use this guide https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
2) Download any sequence file in data/ folder: 
```
mkdir data && cd data
fastq-dump ERR337856  -O .
```
run data/download.sh for IDL-COBS and IDL-RAMBO experiments. The corresponding gene sequence file names are also mentioned in data/fileList.txt.

3) Build the program:
```
mkdir debug
bash build.sh
```
4) Run the program:
```
mkdir results
bash run.sh
```
