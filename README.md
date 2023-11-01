# IDL-RAMBO

The code is tested on a machine running Ubuntu 20.04 with CMake 3.16.3 and GNU Make 4.2.1.

Compile the IDL-RAMBO with
```
make index
make query
```

Download the data files using SRA toolkit. For installing SRA toolkit use this guide https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump

Download sequence files in data/ folder with script /data/download.sh
   
Construct index on 64 CPUs
```
taskset -c 0-63 ./index <hash> <L> <m>
```

Run query on single CPU
```
taskset -c 0 ./query <hash> <L> <m>
```

For example:
```
taskset -c 0-63 ./index IDL 2048 200000000
taskset -c 0 ./query IDL 2048 200000000
taskset -c 0 ./query murmur 2048 2000000000
```

Run script_index.sh and then script_query.sh for replicating the IDL-RAMBO experiments in the paper.

