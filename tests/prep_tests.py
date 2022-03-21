from tests.test_values import *
from ss_validate.schema import SCHEMA
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
        return {
                SCHEMA['fields']['RSID']['label']: snpsarray,
                SCHEMA['fields']['PVAL']['label']: pvalsarray,
                SCHEMA['fields']['CHR']['label']: chrarray,
                SCHEMA['fields']['OR']['label']: orarray,
                SCHEMA['fields']['BP']['label']: bparray,
                SCHEMA['fields']['EFFECT']['label']: effectarray,
                SCHEMA['fields']['OTHER']['label']: otherarray,
                SCHEMA['fields']['EAF']['label']: frequencyarray,
                SCHEMA['fields']['SE']['label']: searray,
                SCHEMA['fields']['BETA']['label']: betaarray,
                SCHEMA['fields']['RANGE_U']['label']: rangearray,
                SCHEMA['fields']['RANGE_L']['label']: rangearray
                }

