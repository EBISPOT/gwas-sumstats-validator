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
    def __init__(self, file, filetype, logfile="VALIDATE.log", error_limit=1000, minrows=MININMUM_ROWS):
        self.file = file
        self.filetype = filetype
        self.schema = None
        self.header = []
        self.cols_to_validate = []
        self.cols_to_read = []
        self.sep = get_seperator(self.file) 
        self.errors = []
        self.bad_rows = []
        self.snp_errors = []
        self.pos_errors = []
        #self.required_fields = STD_COLS
        self.valid_extensions = VALID_FILE_EXTENSIONS
        self.logfile = logfile
        self.error_limit = int(error_limit)
        self.minrows = int(minrows)
        self.nrows = None

        if self.filetype == 'curated' or self.filetype == 'gwas-upload':
            # if curator format allow for more chromosome values
            VALID_CHROMOSOMES.extend(['X', 'x', 'Y', 'y', 'MT', 'Mt', 'mt'])
    
        handler = logging.FileHandler(self.logfile)
        handler.setLevel(logging.ERROR)
        logger.addHandler(handler)


    def setup_field_validation(self):
        self.header = self.get_header()
        if self.filetype == 'curated':
            self.required_fields = [key for key, value in CURATOR_STD_MAP.items() if value == PVAL_DSET]
            self.cols_to_validate = [CURATOR_STD_MAP[h] for h in self.header if h in self.required_fields]
        else:
            self.cols_to_validate = [h for h in self.header if h in VALID_COLS]
        self.cols_to_read = [h for h in self.header if h in VALID_COLS]
        
    def setup_schema(self):
        self.setup_field_validation()
        self.schema = Schema([VALIDATORS[h] for h in self.cols_to_validate])

    def get_header(self):
        first_row = pd.read_csv(self.file, sep=self.sep, comment='#', nrows=1, index_col=False)
        return first_row.columns.values

    def validate_data(self):
        self.setup_field_validation()
        square_file = self.open_file_and_check_for_squareness()
        enough_rows = True
        if self.nrows < self.minrows:
            logger.error("There are only {} rows detected in the file, but the minimum requirement is {}".format(str(self.nrows), str(self.minrows)))
            enough_rows = False
        if square_file is False:
            logger.error("Please fix the table. Some rows have different numbers of columns to the header")
            logger.info("Rows with different numbers of columns to the header are not validated")
        for chunk in self.df_iterator():
            to_validate = chunk[self.cols_to_read] 
            to_validate.columns = self.cols_to_validate # sets the headers to standard format if neeeded
            # validate the snp column if present
            if SNP_DSET in self.header:
                self.schema = Schema([SNP_VALIDATORS[h] for h in self.cols_to_validate])
                errors = self.schema.validate(to_validate)
                self.store_errors(errors, self.snp_errors)
            if CHR_DSET and BP_DSET in self.header:
                self.schema = Schema([POS_VALIDATORS[h] for h in self.cols_to_validate])
                errors = self.schema.validate(to_validate)
                self.store_errors(errors, self.pos_errors)
            self.process_errors()
            if len(self.bad_rows) >= self.error_limit:
                break
        if enough_rows is False or square_file is False:
            logger.info("File is invalid")
        elif not self.bad_rows:
            logger.info("File is valid")
            return True
        else:
            logger.info("File is invalid - {} bad rows, limit set to {}".format(len(self.bad_rows), self.error_limit))
            return False

    def process_errors(self):
        snp_rows = [error.row for error in self.snp_errors]
        pos_rows = [error.row for error in self.pos_errors]
        intersect_errors = [error for error in self.pos_errors if error.row in snp_rows] # error in both the snp and pos (one or the other is fine)
        for error in intersect_errors:
            if len(self.bad_rows) < self.error_limit:
                logger.error(error)
                if error.row not in self.bad_rows:
                    self.bad_rows.append(error.row)
        self.snp_errors = []
        self.pos_errors = []


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

    def validate_file_extension(self):
        check_exts = [check_ext(self.file, ext) for ext in self.valid_extensions]
        if not any(check_exts):
            self.valid_ext = False
            logger.error("File extension should be in {}".format(self.valid_extensions))
            return False
        else:
            self.valid_ext = True
        return True

    def validate_filename(self):
        self.validate_file_extension()
        pmid, study, trait, build = None, None, None, None
        filename = self.file.split('/')[-1].split('.')[0]
        filename_parts = filename.split('-')
        if len(filename_parts) != 4:
            logger.error("Filename: {} should follow the pattern <pmid>-<study>-<trait>-build>.tsv".format(filename))
            return False
        else:
            pmid, study, trait, build = filename_parts
        if not check_build_is_legit(build):
            logger.error("Build: {} is not an accepted build value".format(build))
            return False
        logger.info("Filename looks good!")
        return True

    def df_iterator(self):
        df = pd.read_csv(self.file, 
                         sep=self.sep, 
                         dtype=str, 
                         error_bad_lines=False,
                         warn_bad_lines=False,
                         comment='#', 
                         chunksize=1000000)
        return df

    def check_rows(self, csv_file):
        square = True
        dialect = csv.Sniffer().sniff(csv_file.readline())
        csv_file.seek(0)
        reader = csv.reader(csv_file, dialect)
        self.nrows = 0
        for row in reader:
            if (len(row) != len(self.header)):
                logger.error("Length of row {c} is: {l} instead of {h}".format(c=self.nrows, l=str(len(row)), h=str(len(self.header))))
                square = False
            self.nrows += 1
        return square

    def open_file_and_check_for_squareness(self):
        if pathlib.Path(self.file).suffix in [".gz", ".gzip"]:
             with gzip.open(self.file, 'rt') as f:
                 return self.check_rows(f)
        else: 
            with open(self.file) as f:
                 return self.check_rows(f)

    def validate_headers(self):
        self.setup_field_validation()
        required_is_subset = set(STD_COLS).issubset(self.header)
        if not required_is_subset:
            # check if everything but snp:
            required_is_subset = set(STD_COLS_NO_SNP).issubset(self.header)
            if not required_is_subset:
                required_is_subset = set(STD_COLS_NO_POS).issubset(self.header)
                logger.error("Required headers: {} are not in the file header: {}".format(STD_COLS, self.header))
        return required_is_subset 
        

def check_ext(filename, ext):
    if filename.endswith(ext):
        return True
    return False

def check_build_is_legit(build):
    build_string = build.lower()
    build_number = build_string.replace('build', '')
    if build_number in BUILD_MAP.keys():
        return True
    return False

def get_seperator(file):
    filename, file_extension = os.path.splitext(file)
    sep =  '\t'
    if '.csv' in file_extension:
        sep = ','
    return sep


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-f", help='The path to the summary statistics file to be validated', required=True)
    argparser.add_argument("--filetype", help='The type of file/stage in the process the file to be validated is in. Recommended to leave as default if unknown.', default='gwas-upload', choices=['gwas-upload','curated','standard','harmonised'])
    argparser.add_argument("--logfile", help='Provide the filename for the logs', default='VALIDATE.log')
    argparser.add_argument("--linelimit", help='Stop when this number of bad rows has been found', default=1000)
    argparser.add_argument("--minrows", help='Minimum number of rows acceptable for the file', default=MININMUM_ROWS)
    argparser.add_argument("--drop-bad-lines", help='Store the good lines from the file in a file named <summary-stats-file>.valid', action='store_true', dest='dropbad')
    
    args = argparser.parse_args()

    logfile = args.logfile
    linelimit = args.linelimit
    minrows = args.minrows
    
    validator = Validator(file=args.f, filetype=args.filetype, logfile=args.logfile, error_limit=linelimit, minrows=minrows)
    
    if args.filetype == "curated":
        logger.info("Validating filename...")
        if not validator.validate_filename():
            logger.info("Invalid filename: {}".format(args.f)) 
            logger.info("Exiting before any further checks")
            sys.exit()
    else:
        logger.info("Validating file extension...")
        if not validator.validate_file_extension():
            logger.info("Invalid file extesion: {}".format(args.f)) 
            logger.info("Exiting before any further checks")
            sys.exit()
    
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
