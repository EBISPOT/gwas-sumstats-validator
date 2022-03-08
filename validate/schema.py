from pandas_schema import Column
import numpy as np
from helpers import InInclusiveRangeValidation, match_regex, in_list, in_range, is_dtype, p_value_validation


#==============================#
#Summary Statistics Validation #
#==============================#

SCHEMA = {
    'fields': {
        'VAR_ID': {
            'label': 'variant_id',
            'dtype': str,
            'mandatory': False,
            'description': 'Variant identifier in the form of \
                          <chromosome>_<base_pair_location>_<other_allele>_<effect_allele>',
            'validation': [is_dtype(str),
                           match_regex('[1-25]_[0-9]+_[ACTG]+_[ACTG]+')]
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
            'dtype': float,
            'mandatory': True,
            'dependency': 'NEG_LOG_PVAL',
            'description': 'P-value of the association statistic',
            'validation': [is_dtype(float),
                           p_value_validation]
        },
        'NEG_LOG_PVAL': {
            'label': 'neg_log_10_p_value',
            'dtype': float,
            'mandatory': True,
            'dependency': 'PVAL',
            'description': 'Negative log10 P-value of the association statistic',
            'validation': [is_dtype(float)]
        },
        'CHR': {
            'label': 'chromosome',
            'dtype': int,
            'mandatory': True,
            'description': 'Chromosome where the variant is located (X=23, Y=24, MT=25)',
            'validation': [in_list([str(c) for c in range(1,26)])]
        },
        'BP': {
            'label': 'base_pair_location',
            'dtype': int,
            'mandatory': True,
            'description': 'The first position of the variant in the reference, counting on the bases, from 1 (1-based)',
            'validation': [is_dtype(int),
                           in_range(1, 999999999)]
        },
        'OR': {
            'label': 'odds_ratio',
            'dtype': float,
            'mandatory': True,
            'dependency': 'BETA',
            'description': 'Odds ratio',
            'validation': [is_dtype(float)]
        },
        'BETA': {
            'label': 'beta',
            'dtype': float,
            'mandatory': True,
            'dependency': 'OR',
            'description': 'Beta',
            'validation': [is_dtype(float)]
        },
        'RANGE_U': {
            'label': 'ci_upper',
            'dtype': float,
            'mandatory': False,
            'description': 'Upper confidence interval',
            'validation': [is_dtype(float)]
        },
        'RANGE_L': {
            'label': 'ci_lower',
            'dtype': float,
            'mandatory': False,
            'description': 'Lower confidence interval',
            'validation': [is_dtype(float)]
        },
        'SE': {
            'label': 'standard_error',
            'dtype': float,
            'mandatory': True,
            'description': 'Standard error',
            'validation': [is_dtype(float)]
        },
        'EFFECT': {
            'label': 'effect_allele',
            'dtype': str,
            'mandatory': True,
            'description': 'Allele associated with the effect',
            'validation': [match_regex('^[ACTGactg]+$')]
        },
        'OTHER':{
            'label': 'other_allele',
            'dtype': str,
            'mandatory': True,
            'description': 'The non-effect allele',
            'validation': [match_regex('^[ACTGactg]+$')]
        },
        'EAF': {
            'label': 'effect_allele_frequency',
            'dtype': float,
            'mandatory': True,
            'description': 'Frequency of the effect allele',
            'validation': [is_dtype(float)]
        },
        'INFO': {
            'label': 'info',
            'dtype': float,
            'mandatory': False,
            'description': 'Imputation information metric',
            'validation': [is_dtype(float)]
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
            'validation': [is_dtype(int)]
        }
    },
    'minimum_rows': 100000,
    'valid_file_extensions': [
        ".tsv",
        ".tsv.gz"
    ]
}


#=============================================#
# Field 'validation' for pandas schema to apply #
#=============================================#

VALIDATORS = {
    SCHEMA[fields][VAR_ID_DSET][''label'']: Column(SCHEMA[fields][VAR_ID_DSET][''label''], Column(SCHEMA[fields][VAR_ID_DSET][''validation''], allow_empty = not Column(SCHEMA[fields][VAR_ID_DSET][''mandatory''])
    SNP_DSET: Column(SNP_DSET, [CanConvertValidation(DSET_TYPES[SNP_DSET]), MatchesPatternValidation(r'^rs[0-9]+$')], allow_empty=False),
    CHR_DSET: Column(CHR_DSET, [InListValidation(VALID_CHROMOSOMES)], allow_empty=True),
    BP_DSET: Column(BP_DSET, [CanConvertValidation(DSET_TYPES[BP_DSET]), InExclusiveRangeValidation(0, 999999999)], allow_empty=True),
    PVAL_DSET: Column(PVAL_DSET, [CanConvertValidation(DSET_TYPES[PVAL_DSET]),
                                  InRangeValidationUpperInclusive(0, 1) |
                                  (
                                          CustomSeriesValidation(
                                              lambda x: pd.to_numeric(x.str.split('e|E', expand=True)[1].fillna(value=np.nan)
                                                                      , errors='coerce') < -1,
                                              'Numbers should be between 0 and 1') &
                                          CustomSeriesValidation(
                                              lambda x: pd.to_numeric(x.str.split('e|E', expand=True)[0].fillna(value=np.nan)
                                                                      , errors='coerce') > 0,
                                              'Numbers should be between 0 and 1')
                                  )
                                  ], allow_empty=False),
    OR_DSET: Column(OR_DSET, [CanConvertValidation(DSET_TYPES[OR_DSET])], allow_empty=True),
    RANGE_U_DSET: Column(RANGE_U_DSET, [CanConvertValidation(float)], allow_empty=True),
    RANGE_L_DSET: Column(RANGE_L_DSET, [CanConvertValidation(float)], allow_empty=True),
    BETA_DSET: Column(BETA_DSET, [CanConvertValidation(float)], allow_empty=True),
    SE_DSET: Column(SE_DSET, [CanConvertValidation(float)], allow_empty=True),
    EFFECT_DSET: Column(EFFECT_DSET, [MatchesPatternValidation(r'^[ACTGNactgn]+$')], allow_empty=True),
    OTHER_DSET: Column(OTHER_DSET, [MatchesPatternValidation(r'^[ACTGNactgn]+$')], allow_empty=True),
    FREQ_DSET: Column(FREQ_DSET, [CanConvertValidation(float)], allow_empty=True)
}
