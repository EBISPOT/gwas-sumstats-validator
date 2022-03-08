from pandas_schema import Column
import numpy as np
from helpers import InInclusiveRangeValidation, match_regex, in_list, in_range, is_dtype, p_value_validation


#==============================#
#Summary Statistics Validation #
#==============================#

VALIDATION_SCHEMA = {
    fields: {
        VAR_ID_DSET: {
            label: 'variant_id',
            dtype: str,
            mandatory: False,
            description: 'Variant identifier in the form of \
                          <chromosome>_<base_pair_location>_<other_allele>_<effect_allele>',
            validation: [is_dtype(str),
                         match_regex('[1-25]_[0-9]+_[ACTG]+_[ACTG]+')]
        },
        RSID_DSET: {
            label: 'rsid',
            dtype: str,
            mandatory: False,
            description: 'rsID',
            validation: [is_dtype(str),
                         match_regex('^rs[0-9]+$')]
        },
        PVAL_DSET: {
            label: 'p_value',
            dtype: float,
            mandatory: True,
            dependency: NEG_LOG_PVAL_DSET,
            description: 'P-value of the association statistic',
            validation: [is_dtype(float),
                         p_value_validation]
        },
        NEG_LOG_PVAL_DSET: {
            label: 'neg_log_10_p_value',
            dtype: float,
            mandatory: True,
            dependency: PVAL_DSET,
            description: 'Negative log10 P-value of the association statistic',
            validation: is_dtype(float)
        },
        CHR_DSET: {
            label: 'chromosome',
            dtype: int,
            mandatory: True,
            description: None,
            valiation: in_list([str(c) for c in range(1,26)])
        },
        BP_DSET: {
            label: 'base_pair_location',
            dtype: int,
            mandatory: True
        },
        OR_DSET: {
            label: 'odds_ratio',
            dtype: float,
            mandatory: True,
            dependency: BETA_DSET
        },
        BETA_DSET: {
            label: 'beta',
            dtype: float,
            mandatory: True,
            dependency: OR_DSET
        },
        RANGE_U_DSET: {
            label: 'ci_upper',
            dtype: float,
            mandatory: False,
            description: None
        },
        RANGE_L_DSET: {
            label: 'ci_lower',
            dtype: float,
            mandatory: False
        },
        SE_DSET: {
            label: 'standard_error',
            dtype: float,
            mandatory: True,
            description: None
        },
        EFFECT_DSET: {
            label: 'effect_allele',
            dtype: str,
            mandatory: True
        },
        OTHER_DSET:{
            label: 'other_allele',
            dtype: str,
            mandatory: True
        },
        FREQ_DSET: {
            label: 'effect_allele_frequency',
            dtype: float,
            mandatory: True
        },
        INFO_DSET: {
            label: 'info',
            dtype: float,
            mandatory: False
        },
        HM_CODE_DSET: {
            label: 'hm_code',
            dtype: int,
            mandatory: False
        },
        SAMPLE_SIZE_DSET: {
            label: 'n',
            dtype: int,
            mandatory: False
        }
    },
    minimum_rows: 100000,
    valid_file_extensions: [
        ".tsv",
        ".tsv.gz"
    ]
}


#=============================================#
# Field validation for pandas schema to apply #
#=============================================#

VALIDATORS = {
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
