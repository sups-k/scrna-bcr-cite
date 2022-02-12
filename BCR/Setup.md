# BCR Lineage Tracing Using Immcantation

## Setting up virtual python3 environment on O2

Be present in the home directory.

`module load gcc/9.2.0`
`module load python/3.8.12`

`which python3`
> /n/app/python/3.8.12/bin/python3

`which virtualenv`
> /n/app/python/3.8.12/bin/virtualenv

`mkdir immcantation-python3.8.12/ && cd immcantation-python3.8.12/`

`virtualenv imm-env --system-site-packages`

To activate the environment (same command to activate later on):
`source imm-env/bin/activate`

To deactivate the environment:
`deactivate`


## Installing required modules for Change-O, and Change-O

Be present in `~/immcantation-python3.8.12/`

Within the virtual environment, run:

`pip3 install setuptools`

`pip3 install scipy`

`pip3 install pandas`

`pip3 install biopython`

`pip3 install presto`

`pip3 install airr`

`pip3 install changeo --user`


