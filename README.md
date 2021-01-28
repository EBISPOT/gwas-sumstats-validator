# Summary Statistics TSV file Validator

A file validator for validating GWAS summary statistics TSV files prior to and post [harmonisation](https://github.com/EBISPOT/sum-stats-formatter/tree/master/harmonisation) using [pandas_schema](https://tmiguelt.github.io/PandasSchema/). The purpose is to validate files before their [conversion to HDF5](https://github.com/EBISPOT/SumStats/). 

## Requirements
- python3

## Installation
- `pip install ss-validate`

## Running the validator
To run the validator on a file:
- `ss-validate -f <file_to_validate.tsv> --logfile <logfile_name>`

Information and errors are logged to the console and errors logged to the file specified. A console output might look like:
```
(INFO): Filename is good!
(INFO): Validating file...
(ERROR): Length of row 7 is: 16 instead of 15
(ERROR): Please fix the table. Some rows have different numbers of columns to the header
(INFO): Rows with different numbers of columns to the header are not validated
(ERROR): {row: 1, column: "p_value"}: "-99" was not in the range [0, 1)
```
The errors from the output tell us that row seven has too many columns and row one does not have a valid pvalue. If these rows are not fixed, they will later be dropped and not converted to HDF5. 

### Addional options
- `--linelimit` : _int, default 1000_

   Once this number of erroneous rows has been reached, stop looking for more.
- `--minrows` : _int, default 100000_

   The minimum number of rows the file is required to have in order to validate sucZZcessfully.
- `--drop-bad-lines` : _bool, default False_

   Drops the the lines with errors from the file and writes it to a new file called <file_to_validate.tsv.valid>
- `--stage` : _{'standard', 'harmonised', 'curated'}, default 'standard'_

   The stage the file is in. It is either standard format ('standard'), harmonised ('harmonised') or pre-standard in the custom curated format ('curated'). Recommended to leave as default.
