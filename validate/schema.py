import sys
import numpy as np
from pandas_schema import Column
from pandas_schema.validation import MatchesPatternValidation, InRangeValidation, InListValidation, CustomSeriesValidation, CustomElementValidation, CanConvertValidation, IsDtypeValidation, CanCallValidation
from validate.helpers import InInclusiveRangeValidation

from validate.common_constants import *


MININMUM_ROWS = 100000

STD_COLS = (PVAL_DSET, CHR_DSET, BP_DSET, SNP_DSET) #OR_DSET, RANGE_L_DSET, RANGE_U_DSET, BETA_DSET, SE_DSET, FREQ_DSET , EFFECT_DSET, OTHER_DSET)
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

VALID_FILE_EXTENSIONS = [".tsv", ".csv", ".tsv.gz", ".csv.gz", "gz", "gzip", ".tsv.gzip", ".csv.gzip"]


#VALIDATORS = {
#    PVAL_DSET: Column(PVAL_DSET, [CanConvertValidation(DSET_TYPES[PVAL_DSET]), InInclusiveRangeValidation(0, 1)], allow_empty=False),
#    OR_DSET: Column(OR_DSET, [CanConvertValidation(DSET_TYPES[OR_DSET])], allow_empty=True),
#    RANGE_U_DSET: Column(RANGE_U_DSET, [CanConvertValidation(float)], allow_empty=True),
#    RANGE_L_DSET: Column(RANGE_L_DSET, [CanConvertValidation(float)], allow_empty=True),
#    BETA_DSET: Column(BETA_DSET, [CanConvertValidation(float)], allow_empty=True),
#    SE_DSET: Column(SE_DSET, [CanConvertValidation(float)], allow_empty=True),
#    EFFECT_DSET: Column(EFFECT_DSET, [MatchesPatternValidation(r'^[ACTGNactgn]+$')], allow_empty=True),
#    OTHER_DSET: Column(OTHER_DSET, [MatchesPatternValidation(r'^[ACTGNactgn]+$')], allow_empty=True),
#    FREQ_DSET: Column(FREQ_DSET, [CanConvertValidation(float)], allow_empty=True)
#}

SNP_VALIDATORS = {
    SNP_DSET: Column(SNP_DSET, [CanConvertValidation(DSET_TYPES[SNP_DSET]), MatchesPatternValidation(r'^rs[0-9]+$')], allow_empty=False),
    CHR_DSET: Column(CHR_DSET, [InListValidation(VALID_CHROMOSOMES)], allow_empty=True),
    BP_DSET: Column(BP_DSET, [CanConvertValidation(DSET_TYPES[BP_DSET]), InInclusiveRangeValidation(1, 999999999)], allow_empty=True),
    PVAL_DSET: Column(PVAL_DSET, [CanConvertValidation(DSET_TYPES[PVAL_DSET]), InInclusiveRangeValidation(0, 1)], allow_empty=False),
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
    BP_DSET: Column(BP_DSET, [CanConvertValidation(DSET_TYPES[BP_DSET]), InInclusiveRangeValidation(1, 999999999)], allow_empty=False),
    PVAL_DSET: Column(PVAL_DSET, [CanConvertValidation(DSET_TYPES[PVAL_DSET]), InInclusiveRangeValidation(0, 1)], allow_empty=False),
    OR_DSET: Column(OR_DSET, [CanConvertValidation(DSET_TYPES[OR_DSET])], allow_empty=True),
    RANGE_U_DSET: Column(RANGE_U_DSET, [CanConvertValidation(float)], allow_empty=True),
    RANGE_L_DSET: Column(RANGE_L_DSET, [CanConvertValidation(float)], allow_empty=True),
    BETA_DSET: Column(BETA_DSET, [CanConvertValidation(float)], allow_empty=True),
    SE_DSET: Column(SE_DSET, [CanConvertValidation(float)], allow_empty=True),
    EFFECT_DSET: Column(EFFECT_DSET, [MatchesPatternValidation(r'^[ACTGNactgn]+$')], allow_empty=True),
    OTHER_DSET: Column(OTHER_DSET, [MatchesPatternValidation(r'^[ACTGNactgn]+$')], allow_empty=True),
    FREQ_DSET: Column(FREQ_DSET, [CanConvertValidation(float)], allow_empty=True)
}

