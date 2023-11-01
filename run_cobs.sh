for disk in 0 1
do
    for k in 3 4
    do
        for rangefactor in 0.05 0.1 0.2 0.4 0.8 1.6
        do
            ./build/cobs ./configs/ArBF_10.cfg hash=murmur rangefactor=${rangefactor} k=${k} exp_dir=results/murmur_${rangefactor}_${k}_${disk} query_only=0 result_file=result.txt disk=${disk}
            ./build/cobs ./configs/ArBF_10.cfg hash=fuzzyexp rangefactor=${rangefactor} k=${k} universal_hash_range=8192 exp_dir=results/idl_${rangefactor}_${k}_${disk} query_only=0 result_file=result.txt disk=${disk}
        done
    done
done