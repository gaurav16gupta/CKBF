for range in 268435456 536870912 1073741824 2147483648 4294967296 8589934592 17179869184 34359738368 68719476736
do
    taskset -c 0 ./debug/program ./configs/default.cfg hash=murmur range=${range} > results/murmur_5_${range}
    taskset -c 0 ./debug/program ./configs/default.cfg hash=fuzzyexp range=${range} > results/idl_5_${range}
    taskset -c 0 ./debug/program ./configs/default.cfg hash=murmur k=6 range=${range} > results/murmur_6_${range}
    taskset -c 0 ./debug/program ./configs/default.cfg hash=fuzzyexp k=6 range=${range} > results/idl_6_${range}
done
