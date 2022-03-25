import math
from pandas_schema import Column
import pandas as pd
import numpy as np
from pandas_schema.validation import MatchesPatternValidation,  InListValidation,  CanConvertValidation, CustomSeriesValidation, _SeriesValidation


SNP_DSET = 'variant_id'
MANTISSA_DSET = 'mantissa'
EXP_DSET = 'exponent'
PVAL_DSET = 'p_value'
TRAIT_DSET = 'trait'
STUDY_DSET = 'study_accession'
CHR_DSET = 'chromosome'
BP_DSET = 'base_pair_location'
OR_DSET = 'odds_ratio'
RANGE_U_DSET = 'ci_upper'
RANGE_L_DSET = 'ci_lower'
BETA_DSET = 'beta'
SE_DSET = 'standard_error'
EFFECT_DSET = 'effect_allele'
OTHER_DSET = 'other_allele'
FREQ_DSET = 'effect_allele_frequency'
HM_OR_DSET = 'hm_odds_ratio'
HM_RANGE_U_DSET = 'hm_ci_upper'
HM_RANGE_L_DSET = 'hm_ci_lower'
HM_BETA_DSET = 'hm_beta'
HM_EFFECT_DSET = 'hm_effect_allele'
HM_OTHER_DSET = 'hm_other_allele'
HM_FREQ_DSET = 'hm_effect_allele_frequency'
HM_VAR_ID = 'hm_variant_id'
HM_CODE = 'hm_code'


DSET_TYPES = {SNP_DSET: str, STUDY_DSET: int, PVAL_DSET: float, MANTISSA_DSET: float, EXP_DSET: int, TRAIT_DSET: str,
              CHR_DSET: int, BP_DSET: int, OR_DSET: float, RANGE_U_DSET: float, RANGE_L_DSET: float, BETA_DSET: float, SE_DSET: float,
              EFFECT_DSET: str, OTHER_DSET: str, FREQ_DSET: float, HM_EFFECT_DSET: str,
              HM_OTHER_DSET: str, HM_BETA_DSET: float, HM_OR_DSET: float, HM_FREQ_DSET: float, HM_CODE: int,
              HM_VAR_ID: str, HM_RANGE_L_DSET: float, HM_RANGE_U_DSET: float}

MININMUM_ROWS = 100000

STD_COLS = (PVAL_DSET, CHR_DSET, BP_DSET, SNP_DSET)
STD_COLS_NO_SNP = (PVAL_DSET, CHR_DSET, BP_DSET)
STD_COLS_NO_POS = (PVAL_DSET, SNP_DSET)
VALID_COLS = (PVAL_DSET, OR_DSET, RANGE_L_DSET, RANGE_U_DSET, BETA_DSET, SE_DSET, FREQ_DSET , EFFECT_DSET, OTHER_DSET, CHR_DSET, BP_DSET, SNP_DSET)

CURATOR_STD_MAP = {

    # variant id
    'snp': SNP_DSET,
    # p-value
    'pval': PVAL_DSET,
    # chromosome
    'chr': CHR_DSET,
    # base pair location
    'bp': BP_DSET,
    # odds ratio
    'or': OR_DSET,
    # ci lower
    'ci_lower': RANGE_L_DSET,
    # ci upper
    'ci_upper': RANGE_U_DSET,
    # beta
    'beta': BETA_DSET,
    # standard error
    'se': SE_DSET,
    # effect allele
    'effect_allele': EFFECT_DSET,
    # other allele
    'other_allele': OTHER_DSET,
    # effect allele frequency
    'eaf': FREQ_DSET
}

VALID_CHROMOSOMES = ['1', '2', '3', '4', '5', '6', '7', '8', 
                     '9', '10', '11', '12', '13', '14', '15', '16', 
                     '17', '18', '19', '20', '21', '22', '23', '24', 
                     '25']

BUILD_MAP = {'28': 'NCBI28',
             '29': 'NCBI29',
             '30': 'NCBI30',
             '31': 'NCBI31',
             '33': 'NCBI33',
             '34': 'NCBI34',
             '35': 'NCBI35',
             '36': 'NCBI36',
             '37': 'GRCh37',
             '38': 'GRCh38'}

VALID_FILE_EXTENSIONS = [".tsv", ".csv", ".tsv.gz", ".csv.gz", ".gz", ".gzip", ".tsv.gzip", ".csv.gzip"]


class InInclusiveRangeValidation(_SeriesValidation):
    """
    Checks that each element in the series is within a given inclusive
    numerical range like so x <= n < y.
    Doesn't care if the values are not numeric - it will try anyway.
    """
    def __init__(self, min: float = -math.inf, max: float = math.inf, **kwargs):
        """
        :param min: The minimum (inclusive) value to accept
        :param max: The maximum (inclusive) value to accept
        """
        self.min = min
        self.max = max
        super().__init__(**kwargs)

    @property
    def default_message(self):
        return 'was not >= {} and <= {})'.format(self.min, self.max)

    def validate(self, series: pd.Series) -> pd.Series:
        series = pd.to_numeric(series, errors='coerce')
        return (series >= self.min) & (series <= self.max)


class InExclusiveRangeValidation(_SeriesValidation):
    """
    Checks that each element in the series is within a given exclusive numerical range.
    """
    def __init__(self, min: float = -math.inf, max: float = math.inf, **kwargs):
        """
        :param min: The minimum (exclusive) value to accept
        :param max: The maximum (exclusive) value to accept
        """
        self.min = min
        self.max = max
        super().__init__(**kwargs)

    @property
    def default_message(self):
        return 'was not > {} and < {})'.format(self.min, self.max)

    def validate(self, series: pd.Series) -> pd.Series:
        series = pd.to_numeric(series, errors='coerce')
        return (series > self.min) & (series < self.max)


class InRangeValidationUpperInclusive(_SeriesValidation):
    """
    Checks that each element in the series is within a given numerical range.
    """
    def __init__(self, min: float = -math.inf, max: float = math.inf, **kwargs):
        """
        :param min: The minimum (exclusive) value to accept
        :param max: The maximum (inclusive) value to accept
        """
        self.min = min
        self.max = max
        super().__init__(**kwargs)

    @property
    def default_message(self):
        return 'was not > {} and <= {})'.format(self.min, self.max)

    def validate(self, series: pd.Series) -> pd.Series:
        series = pd.to_numeric(series, errors='coerce')
        return (series > self.min) & (series <= self.max)


SNP_VALIDATORS = {
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

POS_VALIDATORS = {
    SNP_DSET: Column(SNP_DSET, [CanConvertValidation(DSET_TYPES[SNP_DSET]), MatchesPatternValidation(r'^rs[0-9]+$')], allow_empty=True),
    CHR_DSET: Column(CHR_DSET, [InListValidation(VALID_CHROMOSOMES)], allow_empty=False),
    BP_DSET: Column(BP_DSET, [CanConvertValidation(DSET_TYPES[BP_DSET]), InExclusiveRangeValidation(0, 999999999)], allow_empty=False),
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