from pandas_schema import Column
import numpy as np
from ss_validate.helpers import InInclusiveRangeValidation, match_regex, in_list, in_range, is_dtype, p_value_validation

#=====================================#
#Summary Statistics Validation Schema #
#=====================================#


SCHEMA = {
    'fields': {
        'VAR_ID': {
            'label': 'variant_id',
            'dtype': str,
            'mandatory': False,
            'description': 'Variant identifier in the form of \
                          <chromosome>_<base_pair_location>_<other_allele>_<effect_allele>',
            'validation': [is_dtype(str),
                           match_regex('^[0-9]+_[0-9]+_[ACTG]+_[ACTG]+$')]
        },
        'RSID': {
            'label': 'rsid',
            'dtype': str,
            'mandatory': False,
            'description': 'rsID',
            'validation': [is_dtype(str),
                           match_regex('^rs[0-9]+$')]
        },
        'PVAL': {
            'label': 'p_value',
            'column_index': 7,
            'dtype': float,
            'mandatory': True,
            'column_type': 'p_value',
            'description': 'P-value of the association statistic',
            'validation': [is_dtype(float),
                           p_value_validation]
        },
        'NEG_LOG_PVAL': {
            'label': 'neg_log_10_p_value',
            'column_index': 7,
            'dtype': float,
            'mandatory': True,
            'column_type': 'p_value',
            'description': 'Negative log10 P-value of the association statistic',
            'validation': [is_dtype(float),
                           in_range(lower=0)]
        },
        'CHR': {
            'label': 'chromosome',
            'column_index': 0,
            'dtype': int,
            'mandatory': True,
            'description': 'Chromosome where the variant is located (X=23, Y=24, MT=25)',
            'validation': [in_list([str(c) for c in range(1,26)])]
        },
        'BP': {
            'label': 'base_pair_location',
            'column_index': 1,
            'dtype': int,
            'mandatory': True,
            'description': 'The first position of the variant in the reference, counting on the bases, from 1 (1-based)',
            'validation': [is_dtype(int),
                           in_range(1, 999999999)]
        },
        'OR': {
            'label': 'odds_ratio',
            'column_index': 4,
            'dtype': float,
            'mandatory': True,
            'field_type': 'effect_size',
            'description': 'Odds ratio',
            'validation': [is_dtype(float),
                           in_range(lower=0)]
        },
        'BETA': {
            'label': 'beta',
            'column_index': 4,
            'dtype': float,
            'mandatory': True,
            'field_type': 'effect_size',
            'description': 'Beta',
            'validation': [is_dtype(float),
                           in_range()]
        },
        'HR': {
            'label': 'hazard_ratio',
            'column_index': 4,
            'dtype': float,
            'mandatory': True,
            'field_type': 'effect_size',
            'description': 'Hazard ratio',
            'validation': [is_dtype(float),
                           in_range(lower=0)]
        },
        'RANGE_U': {
            'label': 'ci_upper',
            'dtype': float,
            'mandatory': False,
            'description': 'Upper confidence interval',
            'validation': [is_dtype(float),
                           in_range()]
        },
        'RANGE_L': {
            'label': 'ci_lower',
            'dtype': float,
            'mandatory': False,
            'description': 'Lower confidence interval',
            'validation': [is_dtype(float),
                           in_range()]
        },
        'SE': {
            'label': 'standard_error',
            'column_index': 5,
            'dtype': float,
            'mandatory': True,
            'description': 'Standard error',
            'validation': [is_dtype(float),
                           in_range()]
        },
        'EFFECT': {
            'label': 'effect_allele',
            'column_index': 2,
            'dtype': str,
            'mandatory': True,
            'description': 'Allele associated with the effect',
            'validation': [match_regex('^[ACTGactg]+$')]
        },
        'OTHER':{
            'label': 'other_allele',
            'column_index': 3,
            'dtype': str,
            'mandatory': True,
            'description': 'The non-effect allele',
            'validation': [match_regex('^[ACTGactg]+$')]
        },
        'EAF': {
            'label': 'effect_allele_frequency',
            'column_index': 6,
            'dtype': float,
            'mandatory': True,
            'description': 'Frequency of the effect allele',
            'validation': [is_dtype(float),
                           in_range(0, 1)]
        },
        'INFO': {
            'label': 'info',
            'dtype': float,
            'mandatory': False,
            'description': 'Imputation information metric',
            'validation': [is_dtype(float),
                           in_range()]
        },
        'HM_CODE': {
            'label': 'hm_code',
            'dtype': int,
            'mandatory': False,
            'description': 'Harmonisation code, which can be looked up in the metadata to determine the transformation',
            'validation': [is_dtype(int)]
        },
        'SAMPLE_SIZE': {
            'label': 'n',
            'dtype': int,
            'mandatory': False,
            'description': 'Sample size',
            'validation': [is_dtype(int),
                           in_range(lower=0)]
        },
        'REF': {
            'label': 'ref_allele',
            'dtype': str,
            'mandatory': False,
            'description': 'Denote whether the effect or the other allele is the reference allele',
            'validation': [in_list(["ea","oa","NA"])]
        }
    },
    'minimum_rows': 100000,
    'valid_file_extensions': [
        ".tsv",
        ".tsv.gz"
    ]
}
