import bisect
import json
import math
import typing

from src.acmg.helpers import (both_types_cnv, get_overlap, is_duplication, is_overlapping, OverlapTypeEnum,
                              return_dict, RiskEnum, SeverityValuesEnum, sufficient_HI_TS_scores)


def evaluate_transcript(transcript: dict, start: int, end: int, forward_gene: bool) -> (str, str):
    """
    Evaluate transcript for overlap with the CNV (Section 2 in ACMG loss).
    :param transcript: dict - transcript dictionary
    :param start: int - start of the CNV
    :param end: int - end of the CNV
    :param forward_gene: bool - if the gene is on forward strand
    :return: str, str - option and reason for ACMG criteria evaluation
    """
    # Get overlapping CDS, 5', and 3'
    overlapping_cds = []
    overlap_cds = []
    if 'CDS' in transcript:
        overlapping_cds = list(filter(lambda x: get_overlap(x['start'], x['end'], start, end) > 0, transcript['CDS']))
        overlap_cds = list(map(lambda x: get_overlap(x['start'], x['end'], start, end), overlapping_cds))
    overlapping_5p = []
    overlap_5p = []
    if 'five_prime_UTR' in transcript:
        overlapping_5p = list(filter(lambda x: get_overlap(x['start'], x['end'], start, end) > 0, transcript['five_prime_UTR']))
        overlap_5p = list(map(lambda x: get_overlap(x['start'], x['end'], start, end), overlapping_5p))
    overlapping_3p = []
    overlap_3p = []
    if 'three_prime_UTR' in transcript:
        overlapping_3p = list(filter(lambda x: get_overlap(x['start'], x['end'], start, end) > 0, transcript['three_prime_UTR']))
        overlap_3p = list(map(lambda x: get_overlap(x['start'], x['end'], start, end), overlapping_3p))

    # If transcript's 5' is in CNV
    if (forward_gene and start <= transcript['start'] <= end) or (not forward_gene and start <= transcript['end'] <= end):
        if len(overlapping_5p) > 0:
            text_5p = f'Overlaps {len(overlapping_5p)} start (5\') regions with a total length of {sum(overlap_5p)}bp.'
        else:
            text_5p = f'Overlaps no start (5\') regions.'

        # Count the number overlapping CDS
        if len(overlapping_cds) == 0:
            return '2C-2', f'{text_5p}\nOverlaps no CDS regions.'
        else:
            return '2C-1', f'{text_5p}\nOverlaps {len(overlapping_cds)} CDS regions with a total length of {sum(overlap_cds)}bp.'

    # If transcript's 3' is in CNV
    if (forward_gene and start <= transcript['end'] <= end) or (not forward_gene and start <= transcript['start'] <= end):
        if len(overlapping_3p) > 0:
            text_3p = f'Overlaps {len(overlapping_3p)} end (3\') regions with total length of {sum(overlap_3p)}bp.'
        else:
            text_3p = f'Overlaps no end (3\') regions.'
        # Count the number overlapping CDS
        if len(overlapping_cds) == 0:
            return '2D-1', f'{text_3p}\nOverlaps no CDS regions.'
        elif len(overlapping_cds) == 1:
            return '2D-3', f'{text_3p}\nOnly last CDS is involved with overlap of {overlap_cds[0]}bp.'
        else:
            return '2D-4', f'{text_3p}\nOverlaps {len(overlapping_cds)} CDS regions with a total length of {sum(overlap_cds)}bp.'

    # If CNV is in the gene (not very likely)
    if transcript['start'] <= end <= transcript['end'] and transcript['start'] <= start <= transcript['end']:
        return '2E', f'A transcript completely contains the target CNV.\n' \
                     f'Overlaps {len(overlapping_cds)} CDS regions with a total length of {sum(overlap_cds)}bp.'

    # print(f'WARNING: transcript {transcript["ID"]} is not evaluated - it is wholly inside the CNV, although the gene is not.') TODO do we need this?

    return '', ''


def evaluate_gene(gene_info: dict, start: int, end: int) -> (str, str):
    """
    Evaluate gene for overlap with the CNV (Section 2 in ACMG loss).
    :param gene_info: dict - gene dictionary
    :param start: int - start of the CNV
    :param end: int - end of the CNV
    :return: str, str - option and reason for ACMG criteria evaluation (most severe form all transcripts)
    """
    # Check if "gene_type": "protein_coding"
    if gene_info['gene_type'] != 'protein_coding':
        print(f'WARNING: evaluated GENE TYPE is {gene_info["gene_type"]}')

    # Mark the direction (strand)
    forward_gene = gene_info['strand'] == '+'

    # Evaluate all transcripts and pick the most severe (evaluate from the longest)
    most_severe_transcript_name = 'UNKNOWN'
    most_severe_transcript_reason = ''
    most_severe_transcript_option = ''
    severity_seq = ['2C-1', '2D-4', '2D-3', '2C-2', '2E', '2D-2', '2D-1', '']
    for transcript in sorted(gene_info['transcript'], key=lambda x: x['end'] - x['start'], reverse=True):
        option, reason = evaluate_transcript(transcript, start, end, forward_gene)
        if option != '':
            rev_severity = severity_seq.index(option)
            if rev_severity < severity_seq.index(most_severe_transcript_option):
                most_severe_transcript_name = transcript['ID']
                most_severe_transcript_reason = reason
                most_severe_transcript_option = option

    assert most_severe_transcript_option != '', 'Transcripts were evaluated incorrectly.'

    return most_severe_transcript_option, (f'{most_severe_transcript_reason} (Gene name: {gene_info["gene_name"]}, '
                                           f'Transcript ID: {most_severe_transcript_name})')


def get_reason_from_regions(start: int, end: int, regions: list[dict], region_name: str) -> str:
    """
    Get reason from HI/TS regions.
    :param start: int - start of CNV
    :param end: int - end of CNV
    :param regions: list[dict] - info about HI/TS regions
    :param region_name: str - TS/HI
    :return: str - string reason for this selection
    """
    names = [region['ISCA Region Name'] for region in regions]
    overlaps = [get_overlap(start, end, region['start'], region['end']) for region in regions]
    names_str = (', '.join(names[:10]) + ', ...') if len(names) > 10 else ', '.join(names)
    return f'Partially overlaps established {region_name} region(s) ({len(regions)} - {names_str}, max. overlap {max(overlaps)}bp).'


def get_reason_from_benigncnvs(start: int, end: int, benign_cnvs: list[dict]) -> str:
    """
    Get reason from Benign CNVs.
    :param start: int - start of CNV
    :param end: int - end of CNV
    :param benign_cnvs: list[dict] - info about benign CNVs
    :return: str - string reason for this selection
    """
    reasons = []
    for benign_cnv in benign_cnvs:
        overlap = get_overlap(benign_cnv['start'], benign_cnv['end'], start, end)
        reasons.append(f'Accession number: {benign_cnv["variantaccession"]}, pubmedid: {benign_cnv.get("pubmedid", "UNKNOWN")}. '
                       f'Overlap: {overlap}bp ({overlap / float(end - start) * 100.0:.1f}%).')
    reason = f'Overlaps {len(benign_cnvs)} established benign CNV, but contains additional protein coding genes.'
    reason += ' Benign CNVs: (displaying only first 5)\n' if len(reasons) > 5 else ' Benign CNVs:\n'
    reason += '\n'.join(reasons[:5])
    return reason


def section2_classify(start: int, end: int, duplication: bool, hi_ranges: list, hi_genes: list, benign_cnvs: list,
                      genes_getter: typing.Callable[..., list[dict]],
                      find_gene: typing.Callable[[str], typing.Optional[dict]]) -> dict:
    """
    Evaluate Section 2 from ACMG criteria.

    :param start: int - start of the interval
    :param end: int - end of the interval
    :param duplication: bool - True if the CNV is a duplication, False if it is a deletion
    :param hi_ranges: list - List of high-impact ranges on the chromosome
    :param hi_genes: list - List of high-impact genes
    :param benign_cnvs: list - List of known benign CNVs
    :param genes_getter: Callable - Function to get genes in the specified range;
                             accepts start, end, optional gene_type, and overlap_type; returns a list of dictionaries
    :param find_gene: Callable - Function to find a gene by name; accepts gene_name and returns a dictionary or None
    :return: dict - Dictionary containing Section 2 evaluation and additional information
    """
    protein_genes = genes_getter(start, end, gene_type='protein_coding')

    # cont_evaluation clauses:
    cont_eval = {}

    if duplication:
        # Complete containment of a ts_range/ts_gene
        for hi_range in hi_ranges:
            if hi_range['start'] >= start and hi_range['end'] <= end and str(hi_range['Triplosensitivity Score']) in sufficient_HI_TS_scores:
                reason = (f'Completely contains an established TS region {hi_range["ISCA Region Name"]} with TS score '
                          f'{hi_range["Triplosensitivity Score"]}.')
                return return_dict('2', '2A', reason, duplication)
        for hi_gene in hi_genes:
            if hi_gene['start'] >= start and hi_gene['end'] <= end and str(hi_gene['Triplosensitivity Score']) in sufficient_HI_TS_scores:
                reason = f'Completely contains an established TS gene {hi_gene["Gene Symbol"]} with TS score {hi_gene["Triplosensitivity Score"]}.'
                return return_dict('2', '2A', reason, duplication)

        # Partial overlap with ts_range
        ts_ranges = [hi_range for hi_range in hi_ranges if str(hi_range['Triplosensitivity Score']) in sufficient_HI_TS_scores]
        if len(ts_ranges) > 0:
            reason = get_reason_from_regions(start, end, ts_ranges, 'TS')
            cont_eval = return_dict('2', '2B', reason, duplication)

        # Load (protein) coding genes and genes on breakpoints
        genes = sorted([gene['gene_name'] for gene in genes_getter(start, end, overlap_type=OverlapTypeEnum.inside_only)])
        protein_genes_on_breakpoints = genes_getter(start, end, gene_type='protein_coding', overlap_type=OverlapTypeEnum.partial_both)

        # Compare protein coding genes in benign CNV - searching for identical gene content
        for benign_cnv in benign_cnvs:
            genes_cnv = sorted([gene['gene_name'] for gene in benign_cnv['genes']])
            # here we assume that benign CNV start and end in between genes
            if genes_cnv == genes and len(protein_genes_on_breakpoints) == 0:
                reason = f'Identical in gene content ({len(genes)} genes) to a benign CNV gain (variant_accesion={benign_cnv["variantaccession"]}).'
                return return_dict('2', '2C', reason, duplication)

        # Smaller than established benign CNV, breakpoints are OK
        for benign_cnv in benign_cnvs:
            if benign_cnv['start'] <= start and benign_cnv['end'] >= end and len(protein_genes_on_breakpoints) == 0:
                reason = (f'Smaller than an established benign CNV gain (variant_accesion={benign_cnv["variantaccession"]}), breakpoints do not '
                          f'interrupt protein-coding genes.')
                return return_dict('2', '2D', reason, duplication)

        # Smaller than established benign CNV, breakpoints are NOT OK
        for benign_cnv in benign_cnvs:
            if benign_cnv['start'] <= start and benign_cnv['end'] >= end and len(protein_genes_on_breakpoints) > 0:
                gene_names = [g['gene_name'] for g in protein_genes_on_breakpoints]
                reason = (f'Smaller than an established benign CNV gain (variant_accesion={benign_cnv["variantaccession"]}), but breakpoints '
                          f'potentially interrupt protein-coding gene(s) ({", ".join(gene_names)}).')
                if cont_eval == {}:
                    cont_eval = return_dict('2', '2E', reason, duplication)

        # Larger than established benign CNV, identical protein coding genes
        for benign_cnv in benign_cnvs:
            if benign_cnv['start'] > start and benign_cnv['end'] < end and genes_getter(benign_cnv['start'], benign_cnv['end'],
                                                                                        gene_type='protein_coding') == protein_genes:
                reason = (f'Larger than an established benign CNV gain (variant_accesion={benign_cnv["variantaccession"]}), does not include '
                          f'additional protein-coding genes.')
                return return_dict('2', '2F', reason, duplication)

        # Overlapping a benign CNV
        if len(benign_cnvs) > 0 and cont_eval == {}:
            reason = get_reason_from_benigncnvs(start, end, benign_cnvs)
            cont_eval = return_dict('2', '2G', reason, duplication)

        # Complete containment of an HI gene
        for hi_gene in hi_genes:
            if hi_gene['start'] >= start and hi_gene['end'] <= end and hi_gene['Haploinsufficiency Score'] in sufficient_HI_TS_scores:
                reason = f'Completely contains an established HI gene {hi_gene["Gene Symbol"]} with HI score {hi_gene["Haploinsufficiency Score"]}.'
                if cont_eval == {}:
                    cont_eval = return_dict('2', '2H', reason, duplication)

        # Breakpoints with a hi_gene
        for hi_gene in hi_genes:
            start_in_gene = hi_gene['start'] <= start <= hi_gene['end']
            end_in_gene = hi_gene['start'] <= end <= hi_gene['end']

            if start_in_gene and end_in_gene and hi_gene['Haploinsufficiency Score'] in sufficient_HI_TS_scores:
                reason = (f'Both breakpoints are within the same HI gene {hi_gene["Gene Symbol"]} - gene-level sequence variant, possibly resulting '
                          f'in loss of function (LOF).')
                return return_dict('2', '2I', reason, duplication)

            if (start_in_gene or end_in_gene) and hi_gene['Haploinsufficiency Score'] in sufficient_HI_TS_scores:
                reason = f'One breakpoint is within an established HI gene {hi_gene["Gene Symbol"]}, the patient’s phenotype is unknown.'
                if cont_eval == {}:
                    cont_eval = return_dict('2', '2J', reason, duplication)

        # Breakpoint within other gene
        genes_on_breakpoints = genes_getter(start, end, overlap_type=OverlapTypeEnum.partial_both)
        if len(genes_on_breakpoints) > 0:
            gene_names = [g['gene_name'] for g in genes_on_breakpoints]
            reason = f'One or both breakpoints are within gene(s) of no established clinical significance ({gene_names}).'
            if cont_eval == {}:
                cont_eval = return_dict('2', '2L', reason, duplication)

        # Skip section if there are no supporting data
        if cont_eval != {}:
            return cont_eval
        else:
            reason = 'The section is skipped due to lack of supporting data (no TS/HI regions/genes and benign CNVs).'
            return return_dict('2', '2Skip', reason, duplication)

    else:  # deletion
        # Complete containment of a hi_range/hi_gene
        for hi_range in hi_ranges:
            if hi_range['start'] >= start and hi_range['end'] <= end and 'Haploinsufficiency Score' in hi_range.keys() and hi_range[
                'Haploinsufficiency Score'] in sufficient_HI_TS_scores:
                reason = (f'Completely contains an established HI region {hi_range["ISCA Region Name"]} with HI score '
                          f'{int(hi_range["Haploinsufficiency Score"])}.')
                return return_dict('2', '2A', reason, duplication)
        for hi_gene in hi_genes:
            if hi_gene['start'] >= start and hi_gene['end'] <= end and 'Haploinsufficiency Score' in hi_gene.keys() and hi_gene[
                'Haploinsufficiency Score'] in sufficient_HI_TS_scores:
                reason = (f'Completely contains an established HI gene {hi_gene["Gene Symbol"]} with HI score '
                          f'{int(hi_gene["Haploinsufficiency Score"])}.')
                return return_dict('2', '2A', reason, duplication)

        # Evaluation of every single HI gene:
        for hi_gene in hi_genes:
            if 'Haploinsufficiency Score' in hi_gene.keys() and hi_gene['Haploinsufficiency Score'] in sufficient_HI_TS_scores:
                # Get all gene data from gencode
                gene_info = find_gene(hi_gene['Gene Symbol'])
                if gene_info is None:
                    print(f'WARNING: gene {hi_gene["Gene Symbol"]} NOT FOUND in GenCode, probably mismatch in coordinates GenCode/Clingen.')
                    continue
                # Evaluate the gene and all its transcripts, return the worst value
                option, reason = evaluate_gene(gene_info, start, end)
                if option != '' and option != '2D-1':  # 2D-1 is "Continue evaluation"
                    return return_dict('2', option, reason, duplication)
                if option == '2D-1':
                    cont_eval = return_dict('2', option, reason, duplication)

        # Partial overlap with hi_range
        hi_ranges = [hi_range for hi_range in hi_ranges if
                     'Haploinsufficiency Score' in hi_range.keys() and hi_range['Haploinsufficiency Score'] in sufficient_HI_TS_scores]
        if len(hi_ranges) > 0:
            reason = get_reason_from_regions(start, end, hi_ranges, 'HI')
            return return_dict('2', '2B', reason, duplication)

        # Completely contained within a benign CNV
        for benign_cnv in benign_cnvs:
            if benign_cnv['start'] <= start and benign_cnv['end'] >= end:
                reason = (f'An established benign CNV {benign_cnv["variantaccession"]} with pubmedid {benign_cnv.get("pubmedid", "UNKNOWN")} '
                          f'completely contains the target CNV.')
                return return_dict('2', '2F', reason, duplication)

        # Overlapping a benign CNV
        if len(benign_cnvs) > 0 and cont_eval == {}:
            reason = get_reason_from_benigncnvs(start, end, benign_cnvs)
            cont_eval = return_dict('2', '2G', reason, duplication)

        # check HI predictors
        for gene in genes_getter(start, end):
            hi_predictors = ['HIPred', 'Huang', 'GHIS', 'gnomeAD', 'ExAC']
            minimum_predicted = 2
            high_risk = []
            for hi_predictor in hi_predictors:
                try:
                    if gene['conservancy'][hi_predictor]['risk']['loss'] == RiskEnum.HIGH.value:
                        high_risk.append(hi_predictor)
                except KeyError:
                    pass

            if len(high_risk) >= minimum_predicted:
                reason = (f'Overlaps a gene with a predicted high risk of haplo-insufficiency (Gene name: {gene["gene_name"]}, Predictors: '
                          f'{", ".join(high_risk)})')
                return return_dict('2', '2H', reason, duplication)

        # Skip section if no data available
        if cont_eval != {}:
            return cont_eval
        else:
            reason = 'The section is skipped due to lack of supporting data (no HI regions/genes and benign CNVs).'
            return return_dict('2', '2Skip', reason, duplication)


def section1_from_dict(input_dict: dict) -> dict:
    """
    Evaluate Section 1 from ACMG criteria.
    :param input_dict: dict - dict with all annotations of yaml
    :return: dict - section 1 evaluation and info
    """
    # get attributes
    duplication = is_duplication(input_dict['cnv']['cnv_type'])

    # Find number of genes and regulatory elements
    gene_count = len([gene for gene in input_dict['Genes'] if gene['gene_type'] == 'protein_coding'])
    enhancers_count = len([region for region in input_dict['Regulatory'] if region['type'] == 'enhancer'])

    # Pick final option and assign reason
    if gene_count + enhancers_count == 0:
        option = '1B'
        reason = 'The number of overlapping protein-coding genes and regulatory elements is zero.'
    else:
        option = '1A'
        reason = f'The number of overlapping protein-coding genes ({gene_count}) or enhancers ({enhancers_count}) is more than zero.'

    # get info about criteria
    return return_dict('1', option, reason, duplication)


def section2_from_dict(input_dict: dict) -> dict:
    """
    Evaluate Section 2 from ACMG criteria.
    :param input_dict: dict - dict with all annotations of yaml
    :return: dict - section 2 evaluation and info
    """
    # get attributes
    duplication = is_duplication(input_dict['cnv']['cnv_type'])
    start = int(input_dict['cnv']['start'])
    end = int(input_dict['cnv']['end'])

    # Load data
    hi_ranges = input_dict['HI_region']
    hi_genes = input_dict['HI_gene']
    benign_cnvs = input_dict['Benign_CNV_GS_outer']  # allowed options: "Benign_CNV", "Benign_CNV_GS_outer", "Benign_CNV_GS_inner"

    # filter for minimal frequency
    min_frequency = 0.005
    benign_cnvs = [cnv for cnv in benign_cnvs if
                   cnv.get('frequency', '') == '' or math.isnan(cnv['frequency']) or float(cnv['frequency']) >= float(min_frequency)]

    if duplication:  # loss_benign_cnvs_with_gains == 0, for 1 remove this line and de-indent next one
        benign_cnvs = [cnv for cnv in benign_cnvs if both_types_cnv(cnv['cnv_type']) or is_duplication(cnv['cnv_type']) == duplication]

    def genes_getter(start_int: int, end_int: int, gene_type: str = None, overlap_type: OverlapTypeEnum = OverlapTypeEnum.all) -> list[dict]:
        """
        Genes getter function.
        :param start_int: int - start of the interval
        :param end_int: int - end of the interval
        :param gene_type: str - type of gene, None if all
        :param overlap_type: OverlapTypeEnum - type of overlap
        :return: list - list of genes
        """
        # get genes
        if gene_type is None:
            gene_candidates = input_dict['Genes']
        else:
            gene_candidates = [gene for gene in input_dict['Genes'] if gene['gene_type'] == gene_type]

        # filter according to region (the gene is the 'target' region - if inside_only, gene is completely inside the CNV)
        return [gene for gene in gene_candidates if is_overlapping(start_int, end_int, gene['start'], gene['end'], overlap_type)]

    def find_gene(gene_name: str) -> dict | None:
        """
        Find a gene from the input dictionary's 'Genes' list where the gene's 'gene_name' matches the provided gene_name.
        :param gene_name: str - The name of the gene to find.
        :return: dict - The first gene dictionary that matches the gene_name, or an empty dictionary if not found.
        """
        for gene in input_dict['Genes']:
            if gene.get('gene_name') == gene_name:
                return gene
        return None

    return section2_classify(start, end, duplication, hi_ranges, hi_genes, benign_cnvs, genes_getter, find_gene)


def section3_from_dict(input_dict: dict) -> dict:
    """
    Evaluate Section 3 from ACMG criteria.
    :param input_dict: dict - dict with all annotations of yaml
    :return: dict - section 3 evaluation and info
    """
    # get attributes
    duplication = is_duplication(input_dict['cnv']['cnv_type'])

    # Get protein coding genes
    protein_genes = [gene for gene in input_dict['Genes'] if gene['gene_type'] == 'protein_coding']

    # Define threshold
    threshold = [35, 50] if duplication else [25, 35]

    # Evaluate their number and return
    index = bisect.bisect(threshold, len(protein_genes))
    option = ['3A', '3B', '3C'][index]
    reason = f'Overlaps {len(protein_genes)} protein-coding genes.'

    return return_dict('3', option, reason, duplication)


def section4_from_dict(input_dict: dict) -> dict:
    """
    Evaluate Section 4 from ACMG criteria.
    :param input_dict: dict - dict with all annotations of yaml
    :return: dict - section 4q evaluation and info
    """
    # get attributes
    duplication = is_duplication(input_dict['cnv']['cnv_type'])
    chrom = input_dict['cnv']['chromosome']
    start = int(input_dict['cnv']['start'])
    end = int(input_dict['cnv']['end'])

    # check 4O (the only option that is available to automatic interpretation) - consider nfe population (non-finish european)
    common_variability = [variability for variability in input_dict['GnomAD']
                          if variability['start'] <= start and variability['end'] >= end and
                          variability['population'] == 'nfe' and is_duplication(variability['svtype']) == duplication]

    # check frequency >= 1%
    for region in common_variability:
        homref = region['frequencies']['all']['HOMREF']['count']
        homalt = region['frequencies']['all']['HOMALT']['count']
        het = region['frequencies']['all']['HET']['count']

        # skip those with 0 count
        if homalt + homref + het == 0:
            continue

        frequency = (homalt + het) / float(homalt + homref + het)

        if frequency >= 0.01:
            reason = f'Common population variation {chrom}:{start}-{end} for population {region["population"]} has frequency of {frequency * 100.0}%.'
            return return_dict('4', '4O', reason, duplication)

    # Pick final option and assign reason
    option = '4Skip'
    reason = 'Manual decision needed.'

    return return_dict('4', option, reason, duplication)


def section5_from_dict(input_dict: dict) -> dict:
    """
    Evaluate Section 5 from ACMG criteria.
    :param input_dict: dict - dict with all annotations of yaml
    :return: dict - section 5 evaluation and info
    """
    # get attributes
    duplication = is_duplication(input_dict['cnv']['cnv_type'])

    # Pick final option and assign reason
    option = '5F'
    reason = 'No family history is available.'

    return return_dict('5', option, reason, duplication)


def evaluate_from_dict(data_dict: dict) -> dict:
    """
    Evaluate all the ACMG criteria from dictionary of all data.
    :param data_dict: dict - data dictionary with all data needed for evaluation
    :return: dict - evaluated criteria
    """
    # evaluate all criteria
    criteria = [
        section1_from_dict(data_dict),
        section2_from_dict(data_dict),
        section3_from_dict(data_dict),
        section4_from_dict(data_dict),
        section5_from_dict(data_dict)
    ]
    # Get final score and evaluate pathogenic status
    final_score = round(sum([value['score'] for value in criteria if value['score'] is not None]), 8)

    return {
        'criteria': criteria,
        'score': final_score,
        'severity': SeverityValuesEnum.from_score(final_score)
    }


def evaluate_from_json(json_file: str,
                       chromosome: typing.Union[str, None] = None,
                       start: typing.Union[int, None] = None,
                       end: typing.Union[int, None] = None,
                       duplication: typing.Union[bool, None] = None) -> dict:
    """
    Evaluate all the ACMG criteria from JSON file.
    :param json_file: str - JSON file with all annotations
    :param chromosome: str - chromosome
    :param start: int - start of the interval
    :param end: int - end of the interval
    :param duplication: bool/None - cnv type, estimate if None
    :return: dict - evaluated criteria
    """
    # load the JSON file
    with open(json_file) as f:
        json_dict = json.load(f)

    # estimate the CNV type if not set
    if duplication is None:
        try:
            duplication = is_duplication(json_dict['cnv']['cnv_type'])
        except KeyError:
            # estimate according to number of occurrences in JSON
            cnt_dup = json.dumps(json_dict).count('gain') + json.dumps(json_dict).count('Duplication')
            cnt_del = json.dumps(json_dict).count('loss') + json.dumps(json_dict).count('Deletion')
            duplication = cnt_dup >= cnt_del

    # write the data to proper format and add if missing
    json_dict['cnv']['cnv_type'] = 'gain' if duplication else 'loss'
    if 'start' not in json_dict['cnv']:
        json_dict['cnv']['start'] = int(start)
    if 'end' not in json_dict['cnv']:
        json_dict['cnv']['end'] = int(end)
    if 'chromosome' not in json_dict['cnv']:
        json_dict['cnv']['chromosome'] = chromosome if chromosome.startswith('chr') else 'chr' + chromosome
    elif not json_dict['cnv']['chromosome'].startswith('chr'):
        json_dict['cnv']['chromosome'] = 'chr' + json_dict['cnv']['chromosome']

    return evaluate_from_dict(json_dict)
