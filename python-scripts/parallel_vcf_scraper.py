#!/usr/bin/env python3
__description__ =\
"""
Purpose: From a vcf, generate a set of potential SNPs.
"""
__author__ = "Erick Samera"
__version__ = "0.0.1"
__comments__ = "stable;"
# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    RawTextHelpFormatter)
from pathlib import Path
# --------------------------------------------------
import pandas as pd
from datetime import datetime
from multiprocessing import Pool
from collections import Counter
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description=__description__,
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        'input_path',
        type=Path,
        help="path (.vcf)")
    parser.add_argument(
        '--range',
        metavar='INT',
        type=int,
        default="350",
        help='container range (bp) to consider for finding SNPs (default=350)')
    parser.add_argument(
        '--region_min',
        metavar='INT',
        type=int,
        default="100",
        help='minimUM container range (bp) to consider for finding SNPs (default=100)')
    parser.add_argument(
        '--snp_min',
        metavar='INT',
        type=int,
        default=2,
        help='minumum SNP count to be considered in the summary (default=2)')
    parser.add_argument(
        '--density_min',
        metavar='INT',
        type=int,
        default=0,
        help='minumum density to be considered in the summary (default=0)')
    parser.add_argument(
        '--snp_rep_min',
        metavar='FLOAT',
        type=float,
        default=0.9,
        help='minumum sample representation of SNP to be counted')

    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------

    return args
# --------------------------------------------------
def _find_SNP_regions(args, _chromosome_dict: dict, _chromosome_name: str, _anchor: tuple):
    i, (key, value) = _anchor

    snp_regions: list = []
    genotypes: list = []

    for ii, (key_ii, value_ii) in enumerate(list(_chromosome_dict.items())[i:]):
        region_range: int = (int(key_ii) - int(key))
        if args.region_min < region_range < args.range:

            samples: list = [sample.strip() for sample in list(value_ii.keys()) if sample not in ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')]
            if not genotypes:
                genotypes = list(['']*len(samples))
            snp_representation: list = [value_ii[sample] for sample in samples if value_ii[sample] != '.:.:.:.:.:.:.:.']
            genotypes_at_this_position = [value_ii[sample].split(':')[0] for sample in samples]
            if len(snp_representation)/len(samples) < args.snp_rep_min: continue
            for i, genotype in enumerate(genotypes_at_this_position):
                genotypes[i] += f"|{genotype}"
            

            variability = len(Counter(genotypes))/ii

            number_of_genotypes = len(Counter([genotype for genotype in genotypes if genotype.count('.')/ii < 0.5]))
            adjusted_variability = len(Counter([genotype for genotype in genotypes if genotype.count('.')/ii < 0.5]))/ii


            region_range_str: str = f"{_chromosome_name}:{key}-{key_ii}"
            
            snp_count = int(ii+1)
            if snp_count < args.snp_min: continue
            snp_density = snp_count/region_range
            snp_regions.append((region_range_str, region_range, snp_count, snp_density, len(snp_representation)/len(samples), variability, adjusted_variability, number_of_genotypes))
        elif region_range > args.range: break
    return snp_regions


def _mp_find_SNP_regions(args, _chromosome_dict: dict, _chromosome_name: str):
    """
    """
    snp_regions: list = []

    map_args = []
    for i, (key, value) in enumerate(_chromosome_dict.items()):
        map_args.append(tuple([args, _chromosome_dict, _chromosome_name, (i, (key, value))]))

    map_args = tuple(map_args)
    with Pool() as pool:
        snp_regions_raw = pool.starmap(_find_SNP_regions, map_args)

    for snp_region in snp_regions_raw:
        if snp_region: snp_regions += snp_region

    return snp_regions

def _parse_vcf(args, _vcf_path: Path) -> dict:
    """
    """

    vcf_dict: dict = {}
    with open(_vcf_path) as vcf_file:
        for line in [line.strip() for line in vcf_file.readlines()]:
            if line.startswith('#CHROM\t'): headers = [line_info.strip() for line_info in line[1:].split('\t')]
            if line.startswith('#'): continue
            else:
                line_info = [line_info.strip() for line_info in line.split('\t')]
                line_dict: dict = {key: value for key, value in zip(headers, line_info)}
                if line_dict['CHROM'] not in vcf_dict: vcf_dict[line_dict['CHROM']] = {}
                vcf_dict[line_dict['CHROM']].update({line_dict['POS']: line_dict})
    return vcf_dict

def main() -> None:
    """ Insert docstring here """

    args = get_args()

    print(f'{" ".join(datetime.now().isoformat(timespec="seconds").split("T"))} Parsing the VCF ...')
    vcf_dict = _parse_vcf(args, args.input_path)
    print(f'{" ".join(datetime.now().isoformat(timespec="seconds").split("T"))} Parsed!')

    print(f'{" ".join(datetime.now().isoformat(timespec="seconds").split("T"))} Iterating across the SNPs ...')
    snp_region_list: list = []
    for chrom in vcf_dict:
        snp_region_list += _mp_find_SNP_regions(args, vcf_dict[chrom], chrom)

    print(f'{" ".join(datetime.now().isoformat(timespec="seconds").split("T"))} Found {len(snp_region_list)} SNP region(s) across {len(vcf_dict)} chromosome(s)/scaffold(s).')

    snp_region_dataframe = pd.DataFrame(
        sorted(snp_region_list, reverse=True, key=lambda x: x[3]),
        columns=[
            'chrom:start-end', 
            'range length (bp)', 
            'SNPs contained', 
            'SNP density (0-1)', 
            '%% rep. in samples (0-1)', 
            'variability (# of genotypes/# of SNPs)', 
            'adjusted variability',
            'number of genotypes'])
    
    output_file = Path(__file__).parent.joinpath('SNP_regions.csv')
    snp_region_dataframe.to_csv(output_file, index=False)
    print(f'{" ".join(datetime.now().isoformat(timespec="seconds").split("T"))} Exported to {output_file.name}')

    return None
# --------------------------------------------------
if __name__ == '__main__':
    main()