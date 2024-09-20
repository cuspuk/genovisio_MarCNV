import argparse
import gzip
import json
import time

import pandas as pd

from classify_cnv import find_intersections, allowed_chromosomes
from src.acmg.acmg_classify import evaluate_from_dict
from src.acmg.helpers import is_duplication
from src.mongo import get_mongo_database

# Global variable to store the start time
_start_time = None


def progress_bar(iteration: int, total: int, prefix: str = 'Annotating: ', suffix: str = 'Complete',
                 decimals: int = 2, length: int = 50, fill: str = '█', print_end: str = '\r') -> None:
    """
    Display a progress bar in the terminal with an estimated time until completion.

    :param iteration: int - current iteration
    :param total: int - total iterations
    :param prefix: str - prefix string (default: 'Annotating: ')
    :param suffix: str - suffix string (default: 'Complete')
    :param decimals: int - number of decimal places to show in percent complete (default: 1)
    :param length: int - character length of the bar (default: 100)
    :param fill: str - bar fill character (default: '█')
    :param print_end: str - end character (e.g., '\r', '\r\n') (default: '\r')
    :return: None
    """
    global _start_time

    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)

    if iteration > 0:
        elapsed_time = time.time() - _start_time
        estimated_total_time = elapsed_time / iteration * total
        remaining_time = estimated_total_time - elapsed_time
        remaining_time_str = time.strftime('%H:%M:%S', time.gmtime(remaining_time))
    else:
        _start_time = time.time()
        remaining_time_str = '...'

    print(f'\r{prefix} |{bar}| {percent}% {suffix} | Time Remaining: {remaining_time_str}', end=print_end)

    # Print New Line on Complete
    if iteration == total:
        print()


if __name__ == '__main__':
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Classify CNV in a batch mode.')
    parser.add_argument('input', help='A .tsv file with input CNVs.')
    parser.add_argument('output', help='Filename, where to store the resulting classified CNV data.')
    parser.add_argument('--generate-jsons', '-j', action='store_true',
                        help='Generate gzipped jsons for each CNV with name "data_{chrom}_{start}_{end}_{cnv_type}.json.gz".')
    parser.add_argument('--json-directory', '-d', help='Directory, where to store jsons. Ignored if flag "--generate-jsons" not specified',
                        default='json_outputs')
    parser.add_argument('--verbose', '-v', action='store_true', help='Be verbose.')
    parser.add_argument('--skip', type=int, help='Skip first lines from file. Time estimation is not right when this is used.', default=0)
    parser.add_argument('--limit', type=int, help='Limit to maximal number of lines.', default=10 ** 100)
    parser.add_argument('--mongodb_uri', help='MongoDB full URI', default='mongodb://localhost:27017/')
    parser.add_argument('--db_name', help='MongoDB database name', default='genovisio')

    args = parser.parse_args()

    # Load file
    input_table = pd.read_csv(args.input, sep='\t')

    # Adjust columns
    chromosome_column = 'chromosome' if 'chromosome' in input_table.columns else 'chrom'
    assert chromosome_column in input_table.columns, '"chromosome" or "chrom" not in table columns!'
    assert 'start' in input_table.columns, '"start" not in table columns!'
    assert 'end' in input_table.columns, '"end" not in table columns!'
    assert 'cnv_type' in input_table.columns, '"cnv_type" not in table columns!'
    input_table['chromosome'] = input_table[chromosome_column]

    # MongoDB connection
    db = get_mongo_database(args.mongodb_uri, args.db_name)

    # List of collections to traverse
    collection_names = ['Benign_CNV', 'Benign_CNV_GS_inner', 'Benign_CNV_GS_outer', 'Regulatory', 'GnomAD', 'HI_gene', 'HI_region', 'Genes']

    # Go through the file
    if not args.verbose:
        progress_bar(0, min(len(input_table) - args.skip, args.limit))
    for i, row in input_table.iloc[args.skip:args.skip + args.limit].iterrows():
        # Get CNV params
        cnv_params = {k: v for k, v in row.items() if k in ['chromosome', 'start', 'end', 'cnv_type']}

        # Print what is being evaluated
        if args.verbose:
            print(f'Classifying ({i + 1:>6}/{len(input_table):<6}): '
                  f'{cnv_params["chromosome"]}:{cnv_params["start"]}-{cnv_params["end"]}/{cnv_params["cnv_type"]}')

        # Check if the chromosome is good
        if not cnv_params['chromosome'] in allowed_chromosomes:
            if not f'chr{cnv_params["chromosome"]}' in allowed_chromosomes:
                print(f'Not allowed chromosome: "{cnv_params["chromosome"]}"!')
                continue
            else:
                cnv_params['chromosome'] = f'chr{cnv_params["chromosome"]}'

        # Check if the cnv_type is good
        if is_duplication(cnv_params['cnv_type']) is None:
            print(f'Not allowed cnv_type: "{cnv_params["cnv_type"]}"!')
            continue
        else:
            cnv_params['cnv_type'] = 'gain' if is_duplication(cnv_params['cnv_type']) else 'loss'

        # Convert start/end to integers
        try:
            cnv_params['start'] = int(cnv_params['start'])
            cnv_params['end'] = int(cnv_params['end'])
        except ValueError:
            print(f'Cannot convert to int: {cnv_params["start"]} or {cnv_params["start"]}')
            continue

        # Traverse collections and find intersections
        results = {collection_name: find_intersections(db[collection_name], cnv_params) for collection_name in collection_names}
        data_dict = {'cnv': cnv_params, **results}

        # Store jsons?
        if args.generate_jsons:
            with gzip.open(f'{args.json_directory}/data_{cnv_params[chromosome_column]}_{cnv_params["start"]}_{cnv_params["end"]}'
                           f'_{cnv_params["cnv_type"]}.json.gz', 'wt', encoding='utf-8') as f:
                f.write(json.dumps(data_dict, default=str, indent=4))

        # Classify and store
        classification = evaluate_from_dict(data_dict)
        input_table.at[i, 'Score'] = classification['score']
        input_table.at[i, 'Severity'] = classification['severity']
        for section in classification['criteria']:
            num = section['section']
            input_table.at[i, f'Section{num}_option'] = section['option']
            input_table.at[i, f'Section{num}_score'] = section['score']
            input_table.at[i, f'Section{num}_reason'] = section['reason'].replace('\n', '|')

        # Print progress
        if not args.verbose:
            progress_bar(i + 1, min(len(input_table) - args.skip, args.limit))

    # Finally output the updated table
    input_table.to_csv(args.output, sep='\t', index=False)
