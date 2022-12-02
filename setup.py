from setuptools import setup

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='ss-validate',
    description='GWAS summary statistics file validator',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version='1.0.0-dev1',
    packages=['ss_validate'],
    license='Apache License, Version 2.0',
    entry_points={
        "console_scripts": ['ss-validate = ss_validate.validator:main']
    },
    url='https://github.com/EBISPOT/gwas-sumstats-validator',
    author='EBI SPOT',
    author_email='gwas-info@ebi.ac.uk',
    install_requires=['pandas_schema>=0.3.4',
                  'tqdm>=4.48.2'
                 ]
)
