# IDL-RAMBO

compile the IDL-RAMBO with
```
make index
make query
```

Download the data files using SRA toolkit. For installing SRA toolkit use this guide https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump

Download sequence files in data/ folder with script /data/download.sh
   
Construct index on 64 CPUs
```
taskset -c 0-63 ./index IDL 2048 200000000
```

Run query on single CPU
```
taskset -c 0 ./query <hash> <L> <m>
```

For example:
```
./query IDL 2048 200000000
./query murmur 2048 2000000000
```


