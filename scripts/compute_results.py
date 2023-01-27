from argparse import ArgumentParser
from os import path
from tqdm import tqdm

parser = ArgumentParser('create queries')
parser.add_argument('--data-dir', type=str, default='/scratch1/gg29/CKBF/data/fastqFiles')
parser.add_argument('--file-list', type=str, default='data/fileList.txt')
parser.add_argument('--result-file', type=str, default='results/queryResults.txt')
parser.add_argument('--first-n-files', type=int, default=None)
args = parser.parse_args()

file_names = []

with open(args.file_list, 'r') as file_list:
    for file_name in file_list:
        file_names.append(path.join(args.data_dir, file_name.strip()))

file_name_to_queries = {}
query_count = 0

with open(args.result_file, 'r') as result_file:
    for result_line in result_file:
        query_and_file_indices = result_line.strip().split(',')
        if len(query_and_file_indices) > 0:
            query_count += args.first_n_files if isinstance(args.first_n_files, int) else len(file_names)

        for i in range(1, len(query_and_file_indices)):
            if file_names[int(query_and_file_indices[i])] not in file_name_to_queries:
                file_name_to_queries[file_names[int(query_and_file_indices[i])]] = []
            
            file_name_to_queries[file_names[int(query_and_file_indices[i])]].append(query_and_file_indices[0])

fpc = 0
for file_name in file_name_to_queries:
    with open(file_name, 'r') as fastq_file:
        contents = fastq_file.read()
        for query in tqdm(file_name_to_queries[file_name]):
            if query not in contents:
                fpc += 1

print(f'False positive rate {fpc / query_count:.4e}')
