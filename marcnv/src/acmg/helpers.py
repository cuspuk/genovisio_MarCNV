import bisect
import enum
import functools
import math
import os

import pandas as pd

# Gain and loss additional info
SECTION_STRING_LOSS = {
    1: {
        'name': 'Initial Assessment of Genomic Content',
        'info': ''
    },
    2: {
        'name': 'Overlap with Established/Predicted HI or Established Benign Genes/Genomic Regions',
        'info': 'Skip to Section 3 if your copy number loss DOES NOT overlap these types of genes/regions'
    },
    3: {
        'name': 'Evaluation of Gene Number',
        'info': ''
    },
    4: {
        'name': 'Detailed Evaluation of Genomic Content Using Published Literature, Public Databases, and/or Internal Lab Data',
        'info': 'Skip to Section 5 if either your CNV overlapped with an established HI gene/region in Section 2, OR there have been no reports '
                'associating either the CNV or any genes within the CNV with human phenotypes caused by loss of function (LOF) or copy number loss'
    },
    5: {
        'name': 'Evaluation of Inheritance Pattern/Family History for Patient Being Studied',
        'info': ''
    }
}

SECTION_STRING_GAIN = {
    1: {
        'name': 'Initial Assessment of Genomic Content',
        'info': ''
    },
    2: {
        'name': 'Overlap with Established Triplosensitive (TS), Haploinsufficient (HI), or Benign Genes or Genomic Regions',
        'info': 'Skip to Section 3 if the copy number gain does not overlap these types of genes/regions'
    },
    3: {
        'name': 'Evaluation of Gene Number', 'info': ''
    },
    4: {
        'name': 'Detailed Evaluation of Genomic Content Using Published Literature, Public Databases, and/or Internal Lab Data',
        'info': 'If there have been no reports associating either the copy number gain or any of the genes therein with '
                'human phenotypes caused by triplosensitivity, skip to Section 5'
    },
    5: {
        'name': 'Evaluation of Inheritance Pattern/Family History for Patient Being Studied',
        'info': ''
    }
}

# sufficient HI/TS scores for dosage sensitivity
# sufficient_HI_TS_scores = ['1', '2', '3', 1, 2, 3]  # -> hi_all == 1
sufficient_HI_TS_scores = ['3', 3]  # -> hi_all == 0


class RiskEnum(enum.Enum):
    LOW = 'low'
    MEDIUM = 'medium'
    HIGH = 'high'


class OverlapTypeEnum(enum.Enum):
    all = 0  # return all types of overlap
    span_whole_only = 1  # return only ranges, that spans the full query range
    partial_start = 2  # return only ranges, that overlap start of the range (and not end)
    partial_end = 3  # return only ranges, that overlap end of the range (and not start)
    partial_both = 4  # return only ranges, that overlap end xor start (and not the second)
    inside_only = 5  # return only ranges, that are completely inside the full query range


def is_overlapping(start: int, end: int, start2: int, end2: int, overlap_type: OverlapTypeEnum) -> bool:
    """
    Return True if the ranges overlap under the criterion specified by OverlapTypeEnum.

    :param start: int - start of the first range (query range)
    :param end: int - end of the first range (query range)
    :param start2: int - start of the second range (target range)
    :param end2: int - end of the second range (target range)
    :param overlap_type: OverlapTypeEnum - the criterion for overlap
    :return: bool - True if the ranges overlap under the specified criterion, otherwise False
    """

    if overlap_type == OverlapTypeEnum.all:
        # Any overlap: (start2 < end and end2 > start) means ranges overlap
        return start2 < end and end2 > start

    elif overlap_type == OverlapTypeEnum.span_whole_only:
        # Target range must span the entire query range
        return start2 <= start and end2 >= end

    elif overlap_type == OverlapTypeEnum.partial_start:
        # Target range overlaps the start of the query range but not the end
        return start2 <= start < end2 <= end

    elif overlap_type == OverlapTypeEnum.partial_end:
        # Target range overlaps the end of the query range but not the start
        return start < start2 <= end <= end2

    elif overlap_type == OverlapTypeEnum.partial_both:
        # Target range overlaps either start or end of the query range but not both
        return (start2 <= start < end2 < end) or (start < start2 < end <= end2)

    elif overlap_type == OverlapTypeEnum.inside_only:
        # Target range is completely inside the query range
        return start <= start2 and end2 <= end

    else:
        raise ValueError('Invalid overlap type')


# Severity values
class SeverityValuesEnum(enum.Enum):
    BENIGN = 'Benign'
    LBENIGN = 'Likely benign'
    VOUS = 'Uncertain'
    LPATHOGENIC = 'Likely pathogenic'
    PATHOGENIC = 'Pathogenic'

    @classmethod
    def values(cls):
        """Returns a list of all enum values."""
        return [member.value for member in cls]

    @classmethod
    def from_score(cls, score: float):
        """Returns the enum value from score."""
        epsilon = 0.00000001
        severity_thresholds = [-1.0 + epsilon, -0.9 + epsilon, 0.9, 1.0]
        severity_index = bisect.bisect(severity_thresholds, score)
        return cls.values()[severity_index]


@functools.lru_cache(maxsize=2)
def load_acmg(duplication: bool) -> pd.DataFrame:
    """
    Prepare acmg criteria table.
    :param duplication: bool - load duplication acmg texts?
    :return: pd.DataFrame - acmg criteria table as a dataframe
    """
    # load table
    acmg_path = os.path.normpath(os.path.join(os.path.dirname(__file__), 'data/acmg_gain.tsv' if duplication else 'data/acmg_loss.tsv'))
    acmg_table = pd.read_csv(acmg_path, sep='\t', index_col='Index',
                             dtype={'Min Score': float, 'Max Score': float, 'Suggested points': float, 'Pretext': object})
    acmg_table = acmg_table.fillna({'Evidence Type': '', 'Evidence': '', 'Pretext': '', 'Tooltip': ''})

    return acmg_table


@functools.lru_cache(maxsize=2)
def get_acmg_criteria(duplication: bool) -> dict:
    """
    Endpoint for ACMG texts.
    :param duplication: bool - duplication? else deletion
    :return: dict - ACMG texts and ranges as json/dict
    """
    # Load data
    acmg_table = load_acmg(duplication)
    final_dict = SECTION_STRING_GAIN if duplication else SECTION_STRING_LOSS

    # define columns
    columns_names = ['Evidence Type', 'Evidence', 'Suggested points', 'Min Score', 'Max Score', 'Pretext', 'Tooltip']

    # Fill in the texts
    for i, row in acmg_table.iterrows():
        # fix types
        row['Min Score'] = None if math.isnan(row['Min Score']) else float(row['Min Score'])
        row['Max Score'] = None if math.isnan(row['Max Score']) else float(row['Max Score'])
        row['Suggested points'] = None if math.isnan(row['Suggested points']) else float(row['Suggested points'])
        # Parse values:
        line_dict = {column: row[column] for column in columns_names}
        line_dict['Pretext'] = acmg_table.at[row['Pretext'], 'Evidence'] if row['Pretext'] != '' else ''
        # Add to final dict
        index = int(i[:1])
        final_dict[index][i] = line_dict

    # Return final dict
    return final_dict


def get_overlap(a_beg, a_end, b_beg, b_end):
    """
    Gets overlap between two intervals.
    :param a_beg: int - beginning of the first interval
    :param a_end: int - end of the first interval
    :param b_beg: int - beginning of the second interval
    :param b_end: int - end of the second interval
    :return: int - length of the overlap
    """
    return max(0, min(a_end, b_end) - max(a_beg, b_beg))


def normalize_type(cnv_type: str) -> str | None:
    """
    Normalize type of CNV
    :param cnv_type: str - CNV type
    :return: str - CNV type in normalized form
    """
    if str(cnv_type).upper() in ['DUP', 'DUPLICATION', 'GAIN', 'INSERTION']:
        return 'gain'
    if str(cnv_type).upper() in ['DEL', 'DELETION', 'LOSS']:
        return 'loss'
    if str(cnv_type).upper() in ['GAIN+LOSS', 'LOSS+GAIN']:
        return 'gain+loss'
    return None


def is_duplication(cnv_type: str) -> bool | None:
    """
    Is the string duplication of deletion?
    :param cnv_type: str - CNV type
    :return: bool - duplication/deletion/unknown
    """
    normalized_cnv_type = normalize_type(cnv_type)
    if normalized_cnv_type == 'gain': return True
    if normalized_cnv_type == 'loss': return False

    return None


def both_types_cnv(cnv_type: str) -> bool:
    """
    Is the CNV of both types?
    :param cnv_type: str - CNV type
    :return: bool - True if the CNV is of both types
    """
    return normalize_type(cnv_type) == 'gain+loss'


def return_dict(section: str, option: str, reason: str, duplication: bool) -> dict:
    """
    Return dictionary with final evaluation of a section.
    :param section: str - section number as str
    :param option: str - final chosen option
    :param reason: str - reason for choosing
    :param duplication: bool - duplication or deletion
    :return: dict - ACMG option with reason and score
    """
    # Get info about criteria
    score = 0
    evidence = ''
    if option != '':
        # get proper option dict
        texts = get_acmg_criteria(duplication)
        acmg_option = texts[int(option[:1])][option]

        # extract score, evidence name
        score = None if acmg_option['Suggested points'] is None or math.isnan(acmg_option['Suggested points']) else acmg_option['Suggested points']
        evidence = acmg_option['Evidence']
    # return dictionary
    return {
        'section': section,
        'option': option,
        'score': score,
        'reason': reason,
        'evidence': evidence,
    }
