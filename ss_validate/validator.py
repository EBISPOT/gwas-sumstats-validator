import sys
import gzip
import csv
import os
import argparse
import pathlib
import logging
from tqdm import tqdm
from tabulate import tabulate
from pandas_schema import Schema, Column
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
from ss_validate.schema import SCHEMA


"""
GWAS Summary statistics file validator using pandas_schema https://github.com/TMiguelT/PandasSchema
- check filetype
- read in header - check for mandatory
- for each column:
        - read in whole column 
        - sani check row number (all columns need the same numnber else not square)
        - row number needs to exceed min rows
        - if in schema: 
                ss_validate - store error row number in column object?
                - if not conditional column - stop at error limit
- calculate errors:
        - some cols are conditional on others
"""


# set field size limit but catch overflow errors
max_int = sys.maxsize
while True:
    try:
        csv.field_size_limit(max_int)
        break
    except OverflowError:
        max_int = int(max_int/10)

CHUNKSIZE = 100000

logging.basicConfig(level=logging.INFO, format='(%(levelname)s): %(message)s')
logger = logging.getLogger(__name__)


class Validator:
    def __init__(self, file, logfile='VALIDATE.log', schema=SCHEMA, error_limit=1000, minrows=SCHEMA['minimum_rows'], dropbad=False):
        self.file = file
        self.file_label = self.get_file_label()
        self.error_outfile = '.temp_validator_{fl}.parquet'.format(fl=self.file_label)
        self.errors_df = pd.DataFrame(columns=['row'])
        self.schema = schema
        self.pd_schema = None
        self.header = []
        self.conditional_fields = []
        self.rows_to_drop = []
        self.cols_to_validate = []
        self.columns_with_errors = []
        self.sep = get_seperator(self.file) 
        self.errors = []
        self.valid_extensions = SCHEMA['valid_file_extensions']
        self.error_limit = int(error_limit) if dropbad is False else None
        self.minrows = int(minrows)
        self.nrows = None

        handler = logging.FileHandler(logfile)
        handler.setLevel(logging.INFO)
        logger.addHandler(handler)

    def setup_field_validation(self):
        fields = [f['label'] for f in SCHEMA['fields'].values()]
        self.header = self.get_header()
        self.cols_to_validate = [h for h in self.header if h in fields]

    def setup_schema(self):
        self.setup_field_validation()
        self.pd_schema = Schema([VALIDATORS[h] for h in self.cols_to_validate])

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

    def validate_data(self):
        self.setup_field_validation()
        self.write_temp_error_file()
        with tqdm(total=len(self.header)) as pbar:
            for column in self.header:
                self.validate_column(column)
                pbar.update(1)
        status = self.evaluate_errors()
        return status

    def write_temp_error_file(self):
        self.errors_df.to_parquet(self.error_outfile, index=True)

    def temp_error_file_to_df(self):
        df = pd.read_parquet(self.error_outfile)
        return df

    def remove_temp_error_file(self):
        os.remove(self.error_outfile)

    def validate_column(self, column_label):
        logger.info("Validating column: {}".format(column_label))
        column_df = self.col_to_df(column_label)
        pd_schema = self.construct_validator(column_label)
        col_errors_df = pd.DataFrame()
        if pd_schema:
            errors = pd_schema.validate(column_df)
            error_data = [(error.row, ' '.join([str(error.value), error.message])) for error in errors]
            if error_data:
                self.columns_with_errors.append(column_label)
                col_errors_df = col_errors_df.from_records(error_data, columns=['row', column_label])
                self.errors_df = self.merge_with_existing_errors(col_errors_df)
                self.write_temp_error_file()
                self.check_if_exceeding_line_limit()

    def check_if_exceeding_line_limit(self):
        if self.error_limit and len(self.errors_df) >= self.error_limit:
            logger.error("Reached limit of {} errors. Stopping validation process now.".format(self.error_limit))
            self.evaluate_errors()
            sys.exit()

    def find_column_dependency(self, column_to_check):
        dependent_column = None
        for pair in self.conditional_fields:
            if column_to_check in pair:
                dependent_column = list(pair - {column_to_check})[0]
        return dependent_column

    def merge_with_existing_errors(self, col_errors_df):
        existing_errors = self.temp_error_file_to_df()
        merged_df = pd.merge(existing_errors, col_errors_df, on='row', how='outer')
        merged_df.reindex(['row'])
        return merged_df

    def evaluate_errors(self):
        error_df = pd.read_parquet(self.error_outfile)
        if len(error_df) > 0:
            error_df.sort_values(by=['row'], inplace=True)
            self.rows_to_drop = pd.read_parquet(self.error_outfile, columns=['row'])['row'].unique().tolist()
            error_count = error_df.drop(columns=['row']).count()
            logger.error("Some rows have errors. Summary is as follows: \n{}".format(error_count))
            logger.error("A full list of errors is as follows:\n{}".format(tabulate(error_df,
                                                                    headers='keys',
                                                                    tablefmt='psql',
                                                                    showindex=False)))
            return False
        return True

    def col_to_df(self, column_label):
        df = pd.read_table(self.file,
                           sep=self.sep,
                           dtype=str,
                           usecols=[column_label])
        # Below we check if p_value column and add a valid scientific notation value
        # which gives the custom validator an 'e' to split on.
        p_val_label = self.schema['fields']['PVAL']['label']
        if column_label == p_val_label:
            psplit_row = pd.Series({p_val_label:'1e-1000'})
            psplit_row.name = -99
            df = df.append(psplit_row)
        return df

    def construct_validator(self, column_label):
        validator = None
        field_id = self.field_id_from_column_label(column_label)
        if field_id:
            validator = Schema([Column(column_label,
                               self.prop_from_field(field_id, 'validation'),
                               allow_empty = not self.prop_from_field(field_id, 'mandatory')
                               )])
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
                pbar.update(CHUNKSIZE)

    def validate_file_extension(self):
        check_exts = [check_ext(self.file, ext) for ext in self.valid_extensions]
        if not any(check_exts):
            logger.error("File extension should be in {}".format(self.valid_extensions))
            return False
        return True

    def validate_filename(self):
        self.validate_file_extension()
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

    def get_file_label(self):
        return pathlib.Path(self.file).stem

    def df_iterator(self):
        df = pd.read_csv(self.file, 
                         sep=self.sep, 
                         dtype=str, 
                         error_bad_lines=False,
                         warn_bad_lines=False,
                         comment='#', 
                         chunksize=CHUNKSIZE)
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
        self.setup_field_validation()
        mandatory_fields = self.get_mandatory_fields()
        self.conditional_fields = self.get_conditional_fields()
        required_is_subset = set(mandatory_fields).issubset(self.header)
        if not required_is_subset:
            missing_fields = list(set(mandatory_fields) - set(self.header))
            report_missing = []
            for field in missing_fields:
                if self.find_column_dependency(field) not in self.header:
                    report_missing.append(field)
            if len(report_missing) > 0:
                logger.error("Mandatory field(s) missing: {}".format(report_missing))
                return False
        return True

    def get_mandatory_fields(self):
        required = [f['label'] for f in self.schema['fields'].values() if f['mandatory']]
        return required

    def get_conditional_fields(self):
        conditional = [frozenset([f['label'], self.schema['fields'][f['dependency']]['label']]) for f in self.schema['fields'].values() if 'dependency' in f]
        conditional_list = [set(i) for i in set(conditional)]
        return conditional_list


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
    sep = '\t'
    if '.csv' in file_extension:
        sep = ','
    return sep


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-f",
                           help='The path to the summary statistics file to be validated',
                           required=True)
    argparser.add_argument("--logfile",
                           help='Provide the filename for the logs',
                           default='VALIDATE.log')
    argparser.add_argument("--errorlimit",
                           help='Stop when this number of bad rows has been found',
                           default=1000)
    argparser.add_argument("--minrows",
                           help='Minimum number of rows acceptable for the file',
                           default=SCHEMA['minimum_rows'])
    argparser.add_argument("--drop-bad-lines",
                           help='Store the good lines from the file in a file named <summary-stats-file>.valid. \
                                 If this option is used, --errorlimit will be set to None',
                           action='store_true',
                           dest='dropbad')
    args = argparser.parse_args()

    errorlimit = args.errorlimit
    minrows = args.minrows
    drop_bad = args.dropbad
    logfile = args.logfile


    validator = Validator(file=args.f,
                          logfile=logfile,
                          error_limit=errorlimit,
                          minrows=minrows,
                          dropbad=args.dropbad)

    logger.info("Validating file extension...")
    if not validator.validate_file_extension():
        logger.info("Invalid file extesion: {}".format(args.f))
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
        logger.info("Writing good lines to {}.valid".format(args.f))
        validator.write_valid_lines_to_file()

    validator.remove_temp_error_file()


if __name__ == '__main__':
    main()
