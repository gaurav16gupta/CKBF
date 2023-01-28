from argparse import ArgumentParser
from os import path
import random

parser = ArgumentParser('create queries')
parser.add_argument('--data-dir', type=str, default='/scratch1/gg29/CKBF/data/fastqFiles')
parser.add_argument('--file-name', type=str, default='SRR649954.fastq')
parser.add_argument('--len-range', type=int, nargs='+', default=[64, 128])
parser.add_argument('--poisoning', type=int, default=2)
parser.add_argument('--num-queries', type=int, default=10000)
args = parser.parse_args()

cnt = 0

with open(path.join(args.data_dir, args.file_name), 'r') as fastq_file:
    for line_index, line in enumerate(fastq_file):
        if (line_index - 1) % 4 == 0:
            seq = line.strip()
            query_len = random.randint(args.len_range[0], args.len_range[1])
            query_start = random.randint(0, len(seq) - query_len)
            query = seq[query_start:query_start+query_len]
            poisoning_pos = random.sample(range(len(query)), args.poisoning)
            for p in poisoning_pos:
                query = query[:p] + '*' + query[(p+1):]
            
            print(query)
            cnt += 1
            if cnt >= args.num_queries:
                break
