import sys
import gzip
import csv
import os
import argparse
import pathlib
import logging
from tqdm import tqdm
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
from pandas_schema import Schema, Column

from ss_validate.schema import SCHEMA
from ss_validate.helpers import get_version, p_value_validation_allow_zero, is_dtype

"""
GWAS Summary statistics file validator using pandas_schema https://github.com/TMiguelT/PandasSchema
"""


# set field size limit but catch overflow errors
max_int = sys.maxsize
while True:
    try:
        csv.field_size_limit(max_int)
        break
    except OverflowError:
        max_int = int(max_int/10)


logging.basicConfig(level=logging.INFO, format='(%(levelname)s): %(message)s')
logger = logging.getLogger(__name__)


class Validator:
    def __init__(self,
                 file,
                 logfile='VALIDATE.log',
                 schema=SCHEMA,
                 error_limit=1000,
                 minrows=SCHEMA['minimum_rows'],
                 dropbad=False,
                 zero_pvalues=False,
                 chunksize=100000):
        self.file = file
        self.schema = schema
        self.header = []
        self.conditional_fields = []
        self.rows_to_drop = []
        self.cols_to_validate = []
        self.sep = get_seperator(self.file)
        self.errors = []
        self.valid_extensions = SCHEMA['valid_file_extensions']
        self.error_limit = int(error_limit) if dropbad is False else None
        self.minrows = int(minrows)
        self.nrows = None
        self.psplit_row_index = -99
        self.chunksize = chunksize
        self.zero_pvalues = zero_pvalues
        if self.zero_pvalues is True:
            self.allow_zero_pvalues()

        handler = logging.FileHandler(logfile)
        handler.setLevel(logging.INFO)
        logger.addHandler(handler)

    def setup_field_validation(self):
        fields = [f['label'] for f in SCHEMA['fields'].values()]
        self.header = self.get_header()
        self.cols_to_validate = [h for h in self.header if h in fields]

    def get_header(self):
        first_row = pd.read_csv(self.file, sep=self.sep, comment='#', nrows=1, index_col=False)
        return first_row.columns.values

    def validate_file_squareness(self):
        self.setup_field_validation()
        square_file = self.open_file_and_check_for_squareness()
        if square_file is False:
            logger.error("Please fix the table. Some rows have different numbers of columns to the header")
            logger.info("Rows with different numbers of columns to the header are not validated")
            logger.info("File is invalid")
            return False
        return True

    def validate_rows(self):
        if self.nrows < self.minrows:
            logger.error("There are only {} rows detected in the file, but the minimum requirement is {}".format(str(self.nrows), str(self.minrows)))
            logger.info("File is invalid")
            return False
        return True

    def allow_zero_pvalues(self):
        self.schema['fields']['PVAL']['validation'] = [is_dtype(float), p_value_validation_allow_zero]

    def validate_data(self):
        self.setup_field_validation()
        with tqdm(total=self.nrows) as pbar:
            for chunk in self.df_iterator():
                to_validate = self.setup_df_for_validation(chunk)
                pd_schema = self.construct_validator(self.cols_to_validate)
                errors = pd_schema.validate(to_validate)
                self.store_errors(errors, self.errors)
                stop = self.check_if_exceeding_line_limit()
                self.evaluate_errors()
                pbar.update(self.chunksize)
                if stop:
                    break
            if self.rows_to_drop:
                logger.info("File is invalid - {} rows with errors, limit set to {}".format(len(self.rows_to_drop), self.error_limit))
                return False
            logger.info("File is valid")
            return True

    def setup_df_for_validation(self, df):
        """
        A dummy row with a scientific notation style pvalue is added.
        If we don't do this we can't apply the pvalue validation to every chunk
        because they may not have any pvalues in scientific notation.
        set the psplit row index to -99 so that we can filter it out.
        """
        to_validate = df[self.cols_to_validate]
        p_val_label = self.schema['fields']['PVAL']['label']
        psplit_row = pd.Series({p_val_label: '1000e1000'}, name=self.psplit_row_index)
        to_validate = to_validate.append(psplit_row, ignore_index=False)
        return to_validate

    def store_errors(self, errors, store):
        for error in errors:
            if error.row != self.psplit_row_index:
                store.append(error)

    def check_if_exceeding_line_limit(self):
        if self.error_limit and len(self.errors) >= self.error_limit:
            logger.error("Reached limit of {} errors. Stopping validation process now.".format(self.error_limit))
            return True
        return False

    def get_schema_columns_order(self):
        cols_to_check = {}
        for f in self.schema['fields'].values():
            if 'column_index' in f:
                if f['column_index'] in cols_to_check:
                    cols_to_check[f['column_index']].append(f['label'])
                else:
                    cols_to_check[f['column_index']] = [f['label']]
        return cols_to_check

    def find_column_dependency(self, column_to_check):
        dependent_column = None
        for fields in self.conditional_fields:
            if column_to_check in fields:
                dependent_column = list(fields - {column_to_check})[0]
        return dependent_column

    def evaluate_errors(self):
        if len(self.errors) > 0:
            for error in self.errors:
                if error.row not in self.rows_to_drop:
                    self.rows_to_drop.append(error.row)
                    if self.error_limit:
                        if len(self.rows_to_drop) <= self.error_limit:
                            logger.error(error)
                        else:
                            break

    def construct_validator(self, column_list):
        validator_list = []
        for column_label in column_list:
            field_id = self.field_id_from_column_label(column_label)
            if field_id:
                validator_list.append(Column(column_label,
                                             self.prop_from_field(field_id, 'validation'),
                                             allow_empty = not self.prop_from_field(field_id, 'mandatory')
                                             )
                                      )
        validator = Schema(validator_list) if validator_list else None
        return validator

    def prop_from_field(self, field_id, property):
        return self.schema['fields'][field_id][property]

    def field_id_from_column_label(self, column_label):
        column_id = None
        for field, props in self.schema['fields'].items():
            if props['label'] == column_label:
                column_id = field
        return column_id

    def write_valid_lines_to_file(self):
        newfile = self.file + ".valid"
        first_chunk = True
        with tqdm(total=self.nrows) as pbar:
            for chunk in self.df_iterator():
                chunk.drop(self.rows_to_drop, inplace=True, errors='ignore')
                if first_chunk:
                    chunk.to_csv(newfile, mode='w', sep='\t', index=False, na_rep='NA')
                    first_chunk = False
                else:
                    chunk.to_csv(newfile, mode='a', header=False, sep='\t', index=False, na_rep='NA')
                pbar.update(self.chunksize)

    def validate_file_extension(self):
        check_exts = [check_ext(self.file, ext) for ext in self.valid_extensions]
        if not any(check_exts):
            logger.error("File extension should be in {}".format(self.valid_extensions))
            return False
        return True


    def df_iterator(self):
        df = pd.read_csv(self.file, 
                         sep=self.sep, 
                         dtype=str, 
                         error_bad_lines=False,
                         warn_bad_lines=False,
                         comment='#', 
                         chunksize=self.chunksize)
        return df

    def check_rows(self, csv_file):
        square = True
        dialect = csv.Sniffer().sniff(csv_file.readline())
        csv_file.seek(0)
        reader = csv.reader(csv_file, dialect)
        self.nrows = 0
        try:
            for row in reader:
                if len(row) != len(self.header):
                    logger.error("Length of row {c} is: {l} instead of {h}".format(c=self.nrows,
                                                                                   l=str(len(row)),
                                                                                   h=str(len(self.header))))
                    square = False
                self.nrows += 1
        except csv.Error as e:
            logger.error("There was the following error when checking the squareness of the csv: {}".format(e))
            square = False
        return square

    def open_file_and_check_for_squareness(self):
        if pathlib.Path(self.file).suffix in [".gz", ".gzip"]:
             with gzip.open(self.file, 'rt') as f:
                 return self.check_rows(f)
        else: 
            with open(self.file, 'r') as f:
                 return self.check_rows(f)

    def validate_headers(self):
        """
        Assumes that the fields in the schema with a 'column_index' are the mandatory fields.
        Checks that these fields are at their correct index. If more than one field label is
        permitted e.g. 'beta' or 'odds_ratio', one of these must match.
        """
        missing = []
        self.setup_field_validation()
        orders_fields = self.get_schema_columns_order()
        for index, fields in orders_fields.items():
            if len(self.header) > index:
                if self.header[index] not in fields:
                    missing.append(fields)
            else:
                missing.append(fields)
        if len(missing) > 0:
            logger.error("The following fields where either missing or in the wrong order: {}".format(missing))
            return False
        return True


def check_ext(filename, ext):
    if filename.endswith(ext):
        return True
    return False

def get_seperator(file):
    filename, file_extension = os.path.splitext(file)
    sep = '\t'
    if '.csv' in file_extension:
        sep = ','
    return sep


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-f", "--file",
                           help='The path to the summary statistics file to be validated')
    argparser.add_argument("-l", "--logfile",
                           help='Provide the filename for the logs',
                           default='VALIDATE.log')
    argparser.add_argument("-e", "--linelimit",
                           help='Stop when this number of rows with errors has been found',
                           default=1000)
    argparser.add_argument("-m", "--minrows",
                           help='Minimum number of rows acceptable for the file',
                           default=SCHEMA['minimum_rows'])
    argparser.add_argument("-d", "--drop-bad-rows",
                           help='Store the good lines from the file in a file named <summary-stats-file>.valid. \
                                 If this option is used, --linelimit will be set to None',
                           action='store_true',
                           dest='dropbad')
    argparser.add_argument("-v", "--version",
                           help='Just return the version of the validator',
                           action='store_true')
    argparser.add_argument("-z", "--zero_pvalues",
                           help="Use if you want allow p-values of zero",
                           action='store_true')
    args = argparser.parse_args()

    file_to_validate = args.file
    error_limit = args.linelimit
    minrows = args.minrows
    drop_bad = args.dropbad
    logfile = args.logfile
    print_version = args.version
    zero_pvalues = args.zero_pvalues

    if print_version:
        print(get_version())
        sys.exit(0)
    else:
        if not file_to_validate:
            logger.error("the following arguments are required: -f/--file")
            sys.exit()

    validator = Validator(file=file_to_validate,
                          logfile=logfile,
                          error_limit=error_limit,
                          minrows=minrows,
                          dropbad=drop_bad,
                          zero_pvalues=zero_pvalues)

    logger.info("Validating file extension...")
    if not validator.validate_file_extension():
        logger.info("Invalid file extesion: {}".format(file_to_validate))
        logger.info("Exiting before any further checks")
        sys.exit()
    else:
        logger.info("ok")

    logger.info("Validating headers...")
    if not validator.validate_headers():
        logger.info("Invalid headers...exiting before any further checks")            
        sys.exit()
    else:
        logger.info("ok")

    logger.info("Validating file for squareness...")
    if not validator.validate_file_squareness():
        logger.info("Rows are malformed..exiting before any further checks")
        sys.exit()
    else:
        logger.info("ok")

    logger.info("Validating rows...")
    if not validator.validate_rows():
        logger.info("File contains too few rows..exiting before any further checks")
        sys.exit()
    else:
        logger.info("ok")

    logger.info("Validating data...")
    validator.validate_data()
    if drop_bad:
        logger.info("Writing good lines to {}.valid".format(file_to_validate))
        validator.write_valid_lines_to_file()


if __name__ == '__main__':
    main()
