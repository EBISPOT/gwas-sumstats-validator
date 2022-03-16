import sys
import gzip
import csv
import os
import argparse
import pathlib
import logging
from tqdm import tqdm
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

log_levels = {
    'critical': logging.CRITICAL,
    'error': logging.ERROR,
    'warn': logging.WARNING,
    'warning': logging.WARNING,
    'info': logging.INFO,
    'debug': logging.DEBUG
}

logging.basicConfig(level=logging.DEBUG, format='(%(levelname)s): %(message)s')
logger = logging.getLogger(__name__)


class Validator:
    def __init__(self, file, logfile='VALIDATE.log', loglevel=logging.INFO, schema=SCHEMA, error_limit=1000, minrows=SCHEMA['minimum_rows'], dropbad=False):
        self.file = file
        self.file_label = self.get_file_label()
        self.error_outfile = '.temp_validator_{fl}.parquet'.format(fl=self.file_label)
        self.errors_df = pd.DataFrame(columns=['row'])
        self.schema = schema
        self.pd_schema = None
        self.header = []
        self.row_count_dict = {}
        self.cols_to_validate = []
        self.cols_to_read = []
        self.sep = get_seperator(self.file) 
        self.errors = []
        self.bad_rows = []
        self.snp_errors = []
        self.pos_errors = []
        self.valid_extensions = SCHEMA['valid_file_extensions']
        self.error_limit = int(error_limit) if dropbad is False else None
        self.minrows = int(minrows)
        self.nrows = None

        handler = logging.FileHandler(logfile)
        handler.setLevel(loglevel)
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
        self.square_file = self.open_file_and_check_for_squareness()
        if self.square_file is False:
            logger.error("Please fix the table. Some rows have different numbers of columns to the header")
            logger.info("Rows with different numbers of columns to the header are not validated")
            logger.info("File is invalid")
            return False
        return True

    def validate_rows(self):
        self.enough_rows = True
        if self.nrows < self.minrows:
            logger.error("There are only {} rows detected in the file, but the minimum requirement is {}".format(str(self.nrows), str(self.minrows)))
            self.enough_rows = False
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
        self.evaluate_errors()

    def write_temp_error_file(self):
        self.errors_df.to_parquet(self.error_outfile, index=True)

    def temp_error_file_to_df(self):
        df = pd.read_parquet(self.error_outfile)
        return df

    def remove_temp_error_file(self):
        os.remove(self.error_outfile)

    def validate_column(self, column_label):
        logger.info("Validating column: {}".format(column_label))
        column_dependency = self.find_column_dependency(column_label)
        column_df = self.col_to_df(column_label)
        row_count = len(column_df)
        logger.debug("Row count: {}".format(row_count))
        self.row_count_dict[column_label] = row_count
        pd_schema = self.construct_validator(column_label)
        col_errors_df = pd.DataFrame()
        if pd_schema:
            errors = pd_schema.validate(column_df)
            error_data = [(error.row, error.value) for error in errors]
            if error_data:
                col_errors_df = col_errors_df.from_records(error_data, columns=['row', column_label])
                self.errors_df = self.merge_with_existing_errors(col_errors_df)
                self.write_temp_error_file()


    def find_column_dependency(self, column_to_check):
        dependent_column = None
        for pair in self.conditional_fields:
            if column_to_check in pair:
                dependent_column = pair - {column_to_check}
        if dependent_column and dependent_column in self.header:
            return dependent_column
        else:
            return None

    def merge_with_existing_errors(self, col_errors_df):
        existing_errors = self.temp_error_file_to_df()
        merged_df = pd.merge(existing_errors, col_errors_df, on='row', how='outer')
        merged_df.reindex(['row'])
        return merged_df

    def evaluate_errors(self):
        pass

    """
    TODO: store validation details:
    dict: 
        {column1: row_count, column2: row_count...}
    dataframe:
        index=[row_nums],
        data = { 'column1': [errors], 'row_nums': [row_nums]}
    store the df on disk (in case of a large number of errors) perhaps as parquet, indexing on row number.
    conditions:
        if not_equal(row counts):
            report not square and fail
        for column not in dependency_list:
            read errors from file
            if any(error):
                report error and fail
        for each column pair in dependency_list:
            read errors from file
            if errors.row in both columns:
                report error and fail   
    remove error files
    """

    def col_to_df(self, column_label):
        df = pd.read_table(self.file,
                           sep=self.sep,
                           dtype=str,
                           usecols=[column_label])
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

    def validate_data_bak(self):
        self.setup_field_validation()
        with tqdm(total=self.nrows) as pbar:
            for chunk in self.df_iterator():
                to_validate = chunk[self.cols_to_read]
                to_validate.columns = self.cols_to_validate # sets the headers to standard format if needed
                # ss_validate by snp only and then position only
                # first we need to add a dummy row with a scientific notation style pvalue
                # if we don't do this we can't apply the pvalue validation to every chunk
                # because they may not have any pvalues in scientific notation.
                # set the psplit row index to -99 so that we can filter it out.
                self.psplit_row_index = -99
                psplit_row = pd.Series({PVAL_DSET:'1000e1000'}, name=self.psplit_row_index)
                to_validate = to_validate.append(psplit_row, ignore_index=False)
                if SNP_DSET in self.header:
                    self.pd_schema = Schema([SNP_VALIDATORS[h] for h in self.cols_to_validate])
                    errors = self.pd_schema.validate(to_validate)
                    self.store_errors(errors, self.snp_errors)
                if CHR_DSET and BP_DSET in self.header:
                    self.pd_schema = Schema([POS_VALIDATORS[h] for h in self.cols_to_validate])
                    errors = self.pd_schema.validate(to_validate)
                    self.store_errors(errors, self.pos_errors)
                self.process_errors()
                pbar.update(CHUNKSIZE)
                if self.error_limit:
                    if len(self.bad_rows) >= self.error_limit:
                        break
            if not self.bad_rows:
                logger.info("File is valid")
                return True
            else:
                logger.info("File is invalid - {} bad rows, limit set to {}".format(len(self.bad_rows), self.error_limit))
                return False

    def process_errors(self):
        snp_rows = [error.row for error in self.snp_errors]
        errors = []
        if SNP_DSET in self.header and CHR_DSET in self.header:
            errors = [error for error in self.pos_errors if error.row in snp_rows] # error in both the snp and pos (one or the other is fine)
        elif SNP_DSET not in self.header:
            errors = [error for error in self.pos_errors]
        elif CHR_DSET not in self.header:
            errors = [error for error in self.snp_errors]
        for error in errors:
            if self.error_limit:
                if len(self.bad_rows) < self.error_limit:
                    logger.error(error)
            else:
                logger.error(error)
            if error.row not in self.bad_rows:
                    self.bad_rows.append(error.row)
        self.snp_errors = []
        self.pos_errors = []

    def store_errors(self, errors, store):
        for error in errors:
            if error.row != self.psplit_row_index:
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
                if (len(row) != len(self.header)):
                    logger.error("Length of row {c} is: {l} instead of {h}".format(c=self.nrows, l=str(len(row)), h=str(len(self.header))))
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
        self.mandatory_fields = self.get_mandatory_fields()
        self.conditional_fields = self.get_conditional_fields()
        required_is_subset = set(self.mandatory_fields).issubset(self.header)
        if not required_is_subset:
            missing_fields = list(set(self.mandatory_fields) - set(self.header))
            for field in missing_fields:
                for condition in self.conditional_fields:
                    if not (field in condition and set(condition - {field}).issubset(self.header)):
                        logger.error("Required headers: {} are not in the file header: {}".format(self.mandatory_fields, self.header))
        return True
        
    def get_mandatory_fields(self):
        required = [f['label'] for f in SCHEMA['fields'].values() if f['mandatory']]
        return required

    def get_conditional_fields(self):
        conditional = [frozenset([f['label'], SCHEMA['fields'][f['dependency']]['label']]) for f in SCHEMA['fields'].values() if 'dependency' in f]
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
    argparser.add_argument("-f", help='The path to the summary statistics file to be validated', required=True)
    argparser.add_argument("--logfile", help='Provide the filename for the logs', default='VALIDATE.log')
    argparser.add_argument("--loglevel", help='Log level', default='info', choices=[l for l in log_levels.keys()])
    argparser.add_argument("--linelimit", help='Stop when this number of bad rows has been found', default=1000)
    argparser.add_argument("--minrows", help='Minimum number of rows acceptable for the file', default=SCHEMA['minimum_rows'])
    argparser.add_argument("--drop-bad-lines", help='Store the good lines from the file in a file named <summary-stats-file>.valid', action='store_true', dest='dropbad')
    args = argparser.parse_args()

    linelimit = args.linelimit
    minrows = args.minrows

    logfile = args.logfile
    loglevel = log_levels[args.loglevel]

    validator = Validator(file=args.f, logfile=logfile, loglevel=loglevel, error_limit=linelimit, minrows=minrows, dropbad=args.dropbad)

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

    #logger.info("Validating file for squareness...")
    #if not validator.validate_file_squareness():
    #    logger.info("Rows are malformed..exiting before any further checks")
    #    sys.exit()
    #else:
    #    logger.info("ok")
    #
    #logger.info("Validating rows...")
    #if not validator.validate_rows():
    #    logger.info("File contains too few rows..exiting before any further checks")
    #    sys.exit()
    #else:
    #    logger.info("ok")

    #logger.info("Validating data...")
    validator.validate_data()
    if args.dropbad:
        logger.info("Writing good lines to {}.valid".format(args.f))
        validator.write_valid_lines_to_file()


if __name__ == '__main__':
    main()
