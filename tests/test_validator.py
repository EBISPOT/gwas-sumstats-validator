import unittest
import shutil
import os
import tests.prep_tests as prep
import validate.validator as v
from validate.common_constants import *
import tests.test_values as test_arrays


class BasicTestCase(unittest.TestCase):
    def setUp(self):
        self.test_storepath = "./tests/data"
        os.makedirs(self.test_storepath, exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self.test_storepath)        

    def test_validate_good_file_extension(self):
        test_filepath = os.path.join(self.test_storepath, "test_file.tsv.gz")
        validator = v.Validator(test_filepath, "gwas-upload", logfile=test_filepath + ".LOG")
        valid_ext = validator.validate_file_extenstion()
        self.assertTrue(valid_ext)
        # alternative
        test_filepath = os.path.join(self.test_storepath, "test_file.csv.gz")
        validator = v.Validator(test_filepath, "gwas-upload", logfile=test_filepath + ".LOG")
        valid_ext = validator.validate_file_extenstion()
        self.assertTrue(valid_ext)

    def test_validate_bad_file_extension(self):
        test_filepath = os.path.join(self.test_storepath, "test_file.zip")
        validator = v.Validator(test_filepath, "gwas-upload", logfile=test_filepath + ".LOG")
        valid_ext = validator.validate_file_extenstion()
        self.assertFalse(valid_ext)

    def test_validate_good_file_headers(self):
        test_filepath = os.path.join(self.test_storepath, "test_file.tsv")
        setup_file = prep.SSTestFile()
        setup_file.prep_test_file()
        validator = v.Validator(test_filepath, "gwas-upload", logfile=test_filepath + ".LOG")
        valid_headers = validator.validate_headers()
        self.assertTrue(valid_headers)

    def test_validate_bad_file_headers(self):
        test_filepath = os.path.join(self.test_storepath, "test_file.tsv")
        setup_file = prep.SSTestFile()
        setup_file.set_test_data_dict()
        setup_file.test_data_dict.pop(PVAL_DSET) # remove a mandatory field
        setup_file.prep_test_file()
        validator = v.Validator(test_filepath, "gwas-upload", logfile=test_filepath + ".LOG")
        valid_headers = validator.validate_headers()
        self.assertFalse(valid_headers)

    def test_validate_good_file_data(self):
        test_filepath = os.path.join(self.test_storepath, "test_file.tsv")
        logfile=test_filepath.replace('tsv', 'LOG')
        setup_file = prep.SSTestFile()
        setup_file.prep_test_file()
        validator = v.Validator(test_filepath, "gwas-upload", logfile=logfile)
        valid_data = validator.validate_data()

        self.assertTrue(valid_data)

    def test_validate_bad_madatory_file_data(self):
        test_filename = "bad_pval.tsv"
        test_filepath = os.path.join(self.test_storepath, test_filename)
        logfile=test_filepath.replace('tsv', 'LOG')
        setup_file = prep.SSTestFile(filename=test_filename)
        setup_file.set_test_data_dict()
        setup_file.test_data_dict[PVAL_DSET] = ["invalid", -123, "another string", 1.5] # set bad pvalue
        setup_file.prep_test_file()
        validator = v.Validator(file=test_filepath, filetype="gwas-upload", logfile=logfile)
        valid_data = validator.validate_data()
        self.assertEqual(len(validator.bad_rows), 4)
        self.assertFalse(valid_data)


if __name__ == '__main__':
    unittest.main()
