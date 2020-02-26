from setuptools import setup

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='eqss-validate',
    description='Summary statistics file validator',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version='0.0.1',
    packages=['validate'],
    license='Apache License, Version 2.0',
    entry_points={
        "console_scripts": ['eqss-validate = validate.validator:main']
    },
    author='EBI SPOT',
    author_email='gwas-info@ebi.ac.uk',
    install_requires=['pandas_schema']
)
