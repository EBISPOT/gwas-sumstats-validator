from pandas_schema import Column
import numpy as np
from pandas_schema.validation import MatchesPatternValidation, InListValidation,  CanConvertValidation, CustomSeriesValidation
from helpers import InInclusiveRangeValidation

#===============================#
#Summary Statistic Field Labels #
#===============================#

VAR_ID_DSET = 'variant_id'
RSID_DSET = 'rsid'
MANTISSA_DSET = 'mantissa'
EXP_DSET = 'exponent'
PVAL_DSET = 'p_value'
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
HM_CODE = 'hm_code'

#================================#
# Summary Statistic Field Dtypes #
#================================#

DSET_TYPES = {VAR_ID_DSET: str,
              STUDY_DSET: int,
              PVAL_DSET: float,
              MANTISSA_DSET: float,
              EXP_DSET: int,
              CHR_DSET: int,
              BP_DSET: int,
              OR_DSET: float,
              RANGE_U_DSET: float,
              RANGE_L_DSET: float,
              BETA_DSET: float,
              SE_DSET: float,
              EFFECT_DSET: str,
              OTHER_DSET: str,
              FREQ_DSET: float,
              HM_CODE: int
              }

#================================#
# Minimum number of rows in file #
#================================#

MININMUM_ROWS = 100000

#====================#
# Fields to validate #
#====================#

VALID_COLS = (PVAL_DSET,
              OR_DSET,
              RANGE_L_DSET,
              RANGE_U_DSET,
              BETA_DSET,
              SE_DSET,
              FREQ_DSET,
              EFFECT_DSET,
              OTHER_DSET,
              CHR_DSET,
              BP_DSET,
              VAR_ID_DSET
              )

#==================#
# Mandatory fields #
#==================#


#==================================#
# Accepted chromosome values:      #
# 1..25 where 23=X, 24=Y and 25=MT #
#==================================#

VALID_CHROMOSOMES = [str(c) for c in range(1,26)]

#=================#
# File extensions #
#=================#

VALID_FILE_EXTENSIONS = [".tsv", ".tsv.gz"]


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
