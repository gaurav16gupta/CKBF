from argparse import ArgumentParser
from os import path
import random
from tqdm import tqdm

parser = ArgumentParser('create queries')
parser.add_argument('--data-dir', type=str, default='/scratch/gg29/CKBF/data/fastqFiles')
parser.add_argument('--file-name', type=str, default='data/fileList.txt')
parser.add_argument('--num-per-file', type=int, default=1000)
parser.add_argument('--len-range', type=int, nargs='+', default=[64, 128])
parser.add_argument('--poisoning', type=int, default=2)
parser.add_argument('--num-queries', type=int, default=10000)
parser.add_argument('--blacklist', type=str, nargs='+', default=['SRR1524101.fastq', 'SRR1723794.fastq', 'SRR1148176.fastq'])
args = parser.parse_args()

total_cnt = 0

with open(args.file_name, 'r') as file_list:
    for l in file_list:
        file_name = l.strip()
        cnt = 0
        if len(file_name) > 0 and file_name not in args.blacklist:
            with open(path.join(args.data_dir, file_name), 'r') as fastq_file:
                for line_index, line in enumerate(fastq_file):
                    seq = line.strip()
                    if (line_index - 1) % 4 == 0 and len(seq) >= args.len_range[0]:
                        query_len = random.randint(args.len_range[0], min(args.len_range[1], len(seq)))
                        query_start = random.randint(0, len(seq) - query_len)
                        query = seq[query_start:query_start+query_len]
                        poisoning_pos = random.sample(range(len(query)), args.poisoning)
                        for p in poisoning_pos:
                            query = query[:p] + '*' + query[(p+1):]
                        
                        print(query)
                        cnt += 1
                        if cnt >= args.num_per_file:
                            break

                        total_cnt += 1
                        if total_cnt >= args.num_queries:
                            exit(0)
