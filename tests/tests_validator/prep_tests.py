from tests_validator.test_values import *
import pandas as pd


class SSTestFile:
    def __init__(self, filename="test_file.tsv", sep="\t"):
        self.filename = filename
        self.sep = sep
        self.test_data_dict = {}

    def set_test_data_dict(self):
        self.test_data_dict = self.prepare_dictionary()

    def prep_test_file(self):
        if not self.test_data_dict:
            self.set_test_data_dict()
        df = pd.DataFrame.from_dict(self.test_data_dict)
        df.to_csv("./tests/data/{}".format(self.filename), sep=self.sep, index=False, mode='w')

    def prepare_dictionary(self):
        return {SNP_DSET: snpsarray, PVAL_DSET: pvalsarray, CHR_DSET: chrarray, OR_DSET: orarray, BP_DSET: bparray,
                    EFFECT_DSET: effectarray, OTHER_DSET: otherarray, FREQ_DSET: frequencyarray, SE_DSET: searray, BETA_DSET: betaarray,
                    RANGE_L_DSET: rangearray, RANGE_U_DSET: rangearray}

