out_path=results

for range in 268435456 536870912 1073741824 2147483648 4294967296 8589934592 17179869184 34359738368 68719476736
do
    ./build/single_bf ./configs/default.cfg hash=murmur range=${range} > ${out_path}/murmur_5_${range}
    ./build/single_bf ./configs/default.cfg hash=fuzzyexp range=${range} > ${out_path}/idl_5_${range}
    ./build/single_bf ./configs/default.cfg hash=murmur k=6 range=${range} > ${out_path}/murmur_6_${range}
    ./build/single_bf ./configs/default.cfg hash=fuzzyexp k=6 range=${range} universal_hash_range=1024 > ${out_path}/idl_6_${range}
done
