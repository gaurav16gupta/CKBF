taskset -c 128-191 ./index IDL 2048 200000000
taskset -c 128-191 ./index IDL 2048 500000000
taskset -c 128-191 ./index IDL 2048 1000000000
taskset -c 128-191 ./index IDL 2048 1500000000
taskset -c 128-191 ./index IDL 2048 2000000000
taskset -c 128-191 ./index IDL 2048 2500000000

taskset -c 64-127 ./index IDL 4096 200000000
taskset -c 64-127 ./index IDL 4096 500000000
taskset -c 64-127 ./index IDL 4096 1000000000
taskset -c 64-127 ./index IDL 4096 1500000000
taskset -c 64-127 ./index IDL 4096 2000000000
taskset -c 64-127 ./index IDL 4096 2500000000

taskset -c 0-63 ./index murmur 2048 200000000
taskset -c 0-63 ./index murmur 2048 500000000
taskset -c 0-63 ./index murmur 2048 1000000000
taskset -c 0-63 ./index murmur 2048 1500000000
taskset -c 0-63 ./index murmur 2048 2000000000
taskset -c 0-63 ./index murmur 2048 2500000000




