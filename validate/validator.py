import sys
import gzip
import csv
import os
import argparse
import pathlib
import logging
import pandas as pd
from pandas_schema import Schema
from validate.schema import *

"""
GWAS Summary statistics file validator 
- using pandas_schema https://github.com/TMiguelT/PandasSchema
It can be run on files pre- ('standard') or post- harmonisation ('harmonised') 
and also in the 'curator' format. 
File names and numbers of fields on each row are checked. 
Fields validated in the standard and harmonised stage are all the required fields 
for HDF5 convertion. The curator format validation only checks the file name, 
the table shape and the pvalue.
"""


csv.field_size_limit(sys.maxsize)

logging.basicConfig(level=logging.INFO, format='(%(levelname)s): %(message)s')
logger = logging.getLogger(__name__)


class Validator:
    def __init__(self, file, logfile="VALIDATE.log", error_limit=1000):
        self.file = file
        self.schema = None
        self.header = []
        self.cols_to_validate = []
        self.cols_to_read = []
        self.sep = get_seperator(self.file) 
        self.errors = []
        self.bad_rows = []
        self.logfile = logfile
        self.error_limit = int(error_limit)

        handler = logging.FileHandler(self.logfile)
        handler.setLevel(logging.ERROR)
        logger.addHandler(handler)


    def setup_field_validation(self):
        self.header = self.get_header()
        self.cols_to_validate = [h for h in self.header if h in VALID_COLS]
        
    def setup_schema(self):
        self.setup_field_validation()
        self.schema = Schema([VALIDATORS[h] for h in self.cols_to_validate])

    def get_header(self):
        first_row = pd.read_csv(self.file, sep=self.sep, comment='#', nrows=1, index_col=False)
        return first_row.columns.values

    def validate_data(self):
        self.setup_field_validation()
        if not self.open_file_and_check_for_squareness():
            logger.error("Please fix the table. Some rows have different numbers of columns to the header")
            logger.info("Rows with different numbers of columns to the header are not validated")
        for chunk in self.df_iterator():
            to_validate = chunk[self.cols_to_validate]
            self.setup_schema()
            errors = self.schema.validate(to_validate)
            self.store_errors(errors, self.errors)
            self.process_errors()
            if len(self.bad_rows) >= self.error_limit:
                break
        if not self.bad_rows:
            logger.info("File is valid")
            return True
        else:
            logger.info("File is invalid - {} bad rows, limit set to {}".format(len(self.bad_rows), self.error_limit))
            return False

    def process_errors(self):
        for error in self.errors:
            if len(self.bad_rows) < self.error_limit:
                logger.error(error)
                if error.row not in self.bad_rows:
                    self.bad_rows.append(error.row)
        self.errors = []


    @staticmethod
    def store_errors(errors, store):
        for error in errors:
            store.append(error)

    def write_valid_lines_to_file(self):
        newfile = self.file + ".valid"
        first_chunk = True
        for chunk in self.df_iterator():
            chunk.drop(self.bad_rows, inplace=True, errors='ignore')
            if first_chunk:
                chunk.to_csv(newfile, mode='w', sep='\t', index=False, na_rep='NA')
                first_chunk = False
            else:
                chunk.to_csv(newfile, mode='a', header=False, sep='\t', index=False, na_rep='NA')

    def df_iterator(self):
        df = pd.read_csv(self.file, 
                         sep=self.sep, 
                         dtype=str, 
                         error_bad_lines=False,
                         warn_bad_lines=False,
                         comment='#', 
                         chunksize=1000000)
        return df

    def check_file_is_square(self, csv_file):
        square = True
        dialect = csv.Sniffer().sniff(csv_file.readline())
        csv_file.seek(0)
        reader = csv.reader(csv_file, dialect)
        count = 1
        for row in reader:
            if (len(row) != len(self.header)):
                logger.error("Length of row {c} is: {l} instead of {h}".format(c=count, l=str(len(row)), h=str(len(self.header))))
                square = False
            count += 1
        return square

    def open_file_and_check_for_squareness(self):
        if pathlib.Path(self.file).suffix in [".gz", ".gzip"]:
             with gzip.open(self.file, 'rt') as f:
                 return self.check_file_is_square(f)
        else: 
            with open(self.file) as f:
                 return self.check_file_is_square(f)

    def validate_headers(self):
        self.setup_field_validation()
        required_is_subset = set(VALID_COLS).issubset(self.header)
        if not required_is_subset:
            logger.error("Required headers: {} are not in the file header: {}".format(VALID_COLS, self.header))
        return required_is_subset 
        
def get_seperator(file):
    filename, file_extension = os.path.splitext(file)
    sep =  '\t'
    if '.csv' in file_extension:
        sep = ','
    return sep


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-f", help='The path to the summary statistics file to be validated', required=True)
    argparser.add_argument("--logfile", help='Provide the filename for the logs', default='VALIDATE.log')
    argparser.add_argument("--linelimit", help='Stop when this number of bad rows has been found', default=1000)
    argparser.add_argument("--drop-bad-lines", help='Store the good lines from the file in a file named <summary-stats-file>.valid', action='store_true', dest='dropbad')
    
    args = argparser.parse_args()

    logfile = args.logfile
    linelimit = args.linelimit
    
    validator = Validator(file=args.f, logfile=args.logfile, error_limit=linelimit)
    
    logger.info("Validating headers...")
    if not validator.validate_headers():
        logger.info("Invalid headers...exiting before any further checks")            
        sys.exit()

    logger.info("Validating data...")
    validator.validate_data()
    if args.dropbad:
        logger.info("Writing good lines to {}.valid".format(args.f))
        validator.write_valid_lines_to_file()


if __name__ == '__main__':
    main()
