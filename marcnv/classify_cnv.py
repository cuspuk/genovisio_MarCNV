import argparse
import gzip
import json
import re
import typing

import pymongo

from marcnv.src.acmg.acmg_classify import evaluate_from_dict
from marcnv.src.mongo import get_mongo_database

# define allowed chromosomes
allowed_chromosomes = [f'chr{chr_id}' for chr_id in list(range(1, 23)) + ['X', 'Y']]


def parse_input(input_str: str) -> dict[str, str | int]:
    """
    Parse the input string of the form 'chr1:10000-20000/DEL' and return a dictionary.
    """
    match = re.match(r'(chr[\dXY]+):(\d+)-(\d+)/(\w+)', input_str)
    cnv_type_map = {'del': 'loss', 'dup': 'gain', 'gain': 'gain', 'loss': 'loss'}
    if not match or match.group(4).lower() not in cnv_type_map.keys() or not match.group(1) in allowed_chromosomes:
        raise ValueError(f'Input format must be "chr1:10000-20000/del". CNV type should be {"/".join(cnv_type_map.keys())}. '
                         f'Chromosome should be {"/".join(allowed_chromosomes)}')

    return {
        'chromosome': match.group(1),
        'start': int(match.group(2)),
        'end': int(match.group(3)),
        'cnv_type': cnv_type_map[match.group(4).lower()]
    }


def find_intersections(collection: pymongo.collection.Collection, search_params: dict[str, str | int]) -> list[dict[str, typing.Any]]:
    """
    Find intersecting items within a single collection using indexed queries.
    """
    # Use indexed query to find potential intersections
    query = {
        'chromosome': search_params['chromosome'],
        'start': {'$lte': search_params['end']},  # search_params.start <= other.end
        'end': {'$gte': search_params['start']}  # search_params.end >= other.start
    }

    return list(collection.find(query))


def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Classify CNV and/or find intersecting items in MongoDB collections.')
    parser.add_argument('input', help='Input string in the form "chr1:10000-20000/del". CNV type should be del/dup/loss/gain.')
    parser.add_argument('--json-output', '-j', help='Filename, where to store the resulting json. (add .gz for gzipped output)', default=None)
    parser.add_argument('--skip-classification', '-s', action='store_true', help='Skips the automatic classification')
    parser.add_argument('--mongodb_uri', help='MongoDB full URI', default='mongodb://localhost:27017/')
    parser.add_argument('--db_name', help='MongoDB database name', default='genovisio')

    args = parser.parse_args()

    # Parse the input string
    search_params = parse_input(args.input)

    # MongoDB connection
    db = get_mongo_database(args.mongodb_uri, args.db_name)

    # List of collections to traverse
    collection_names = ['Benign_CNV', 'Benign_CNV_GS_inner', 'Benign_CNV_GS_outer', 'Regulatory', 'GnomAD', 'HI_gene', 'HI_region', 'Genes']

    # Traverse collections and find intersections
    results = {collection_name: find_intersections(db[collection_name], search_params) for collection_name in collection_names}

    # Convert the results to JSON format (add cnv info to beginning)
    data_dict = {'cnv': search_params, **results}
    results_json = json.dumps(data_dict, default=str, indent=4)

    # Write into a file?
    if args.json_output is not None:
        with (gzip.open(args.json_output, 'wt', encoding='utf-8') if args.json_output.endswith('.gz') else
              open(args.json_output, 'w', encoding='utf-8') as f):
            f.write(results_json)

    # Classify?
    if not args.skip_classification:
        print(json.dumps(evaluate_from_dict(data_dict), default=str, indent=4))


if __name__ == '__main__':
    main()
