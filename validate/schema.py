import sys
import numpy as np
from pandas_schema import Column
from pandas_schema.validation import MatchesPatternValidation, InRangeValidation, InListValidation, CustomSeriesValidation, CustomElementValidation, CanConvertValidation, IsDtypeValidation, CanCallValidation
from validate.helpers import InInclusiveRangeValidation

from validate.common_constants import *

VALID_COLS = TO_LOAD_DSET_HEADERS_DEFAULT

VALID_CHROMOSOMES = ['1', '2', '3', '4', '5', '6', '7', '8', 
                     '9', '10', '11', '12', '13', '14', '15', '16', 
                     '17', '18', '19', '20', '21', '22', '23', '24', 
                     '25']


VALIDATORS = {
    PVAL_DSET: Column(PVAL_DSET, [CanConvertValidation(DSET_TYPES[PVAL_DSET]), InInclusiveRangeValidation(0, 1)], allow_empty=True),
    BETA_DSET: Column(BETA_DSET, [CanConvertValidation(float)], allow_empty=True),
    SNP_DSET: Column(SNP_DSET, [CanConvertValidation(DSET_TYPES[SNP_DSET]), MatchesPatternValidation(r'^chr[0-9]+_[0-9]+_[ACTGNactgn]+_[ACTGNactgn]+|LONG_STRING$')], allow_empty=True),
    CHR_DSET: Column(CHR_DSET, [InListValidation(VALID_CHROMOSOMES)], allow_empty=True),
    BP_DSET: Column(BP_DSET, [CanConvertValidation(DSET_TYPES[BP_DSET]), InInclusiveRangeValidation(1, 999999999)], allow_empty=True),
    EFFECT_DSET: Column(EFFECT_DSET, [MatchesPatternValidation(r'^[ACTGNactgn]+|LONG_STRING$')], allow_empty=True),
    OTHER_DSET: Column(OTHER_DSET, [MatchesPatternValidation(r'^[ACTGNactgn]+|LONG_STRING$')], allow_empty=True),
    FREQ_DSET: Column(FREQ_DSET, [CanConvertValidation(float)], allow_empty=True),
    RSID_DSET: Column(RSID_DSET, [CanConvertValidation(DSET_TYPES[RSID_DSET]), MatchesPatternValidation(r'^rs[0-9]+$')], allow_empty=True),
    MUTATION_DSET: Column(MUTATION_DSET, [CanConvertValidation(DSET_TYPES[MUTATION_DSET])], allow_empty=True),
    PHEN_DSET: Column(PHEN_DSET, [CanConvertValidation(DSET_TYPES[PHEN_DSET])], allow_empty=True),
    AC_DSET: Column(AC_DSET, [CanConvertValidation(DSET_TYPES[AC_DSET])], allow_empty=True),
    AN_DSET: Column(AN_DSET, [CanConvertValidation(DSET_TYPES[AN_DSET])], allow_empty=True),
    R2_DSET: Column(R2_DSET, [CanConvertValidation(DSET_TYPES[R2_DSET])], allow_empty=True),
    GENE_DSET: Column(GENE_DSET, [CanConvertValidation(DSET_TYPES[GENE_DSET])], allow_empty=True),
    MTO_DSET: Column(MTO_DSET, [CanConvertValidation(DSET_TYPES[MTO_DSET])], allow_empty=True),
    EXPR_DSET: Column(EXPR_DSET, [CanConvertValidation(DSET_TYPES[EXPR_DSET])], allow_empty=True)
}

