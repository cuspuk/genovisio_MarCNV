import functools
import os

import pandas as pd


@functools.lru_cache()
def cytobands_list(chromosome: str, start: int, end: int) -> list[str]:
    """
    Get all cytobands in the genomic location
    :param chromosome: str - full chromosome identifier ('chr12')
    :param start: int - location start
    :param end: int - location end
    :return: list - all cytobands in the location
    """
    # Load cytoband file
    cytoband_path = os.path.normpath(os.path.join(os.path.dirname(__file__), 'cytoBand_hg38.txt'))
    cytobands = pd.read_csv(cytoband_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'color'])

    # Get overlap with the required location
    cyto_overlap = (cytobands['chrom'] == chromosome) & (cytobands['end'] >= start) & (cytobands['start'] <= end)

    # Return list of cytoband names
    return list(cytobands[cyto_overlap]['name'])


def cytobands_desc(chromosome: str, start: int, end: int) -> str:
    """
    Get cytoband description/range as a human-readable string.
    :param chromosome: str - full chromosome identifier ('chr12')
    :param start: int - location start
    :param end: int - location end
    :return: str - cytobands description of the range
    """
    # Get all cytobands
    all_cytobands = cytobands_list(chromosome, start, end)

    # Return error if none:
    if len(all_cytobands) == 0:
        return 'unknown'

    # Return one if one:
    if len(all_cytobands) == 1:
        return all_cytobands[0]

    # Return range
    return f'{all_cytobands[0]}-{all_cytobands[-1]}'
