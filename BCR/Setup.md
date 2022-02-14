# BCR Lineage Tracing Using Immcantation (v.4.3.0) - Setting Up Tools

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


## To install other required non-python modules (Be present in `~/immcantation-python3.8.12/imm-env/bin`):

### For MUSCLE 3.8.31

`wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz`

`tar -zxvf muscle3.8.31_i86linux64.tar.gz`

`rm -r *.tar.gz`

`chmod +x muscle3.8.31_i86linux64`

### For tbl2asn

`wget https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz`

`gunzip linux64.tbl2asn.gz`

`mv linux64.tbl2asn tbl2asn`

`chmod +x tbl2asn`

### For IgBLAST

Go to https://ftp.ncbi.nih.gov/blast/executables/igblast/release/
This opens in Finder on MacBook
Copy the folder `1.17.1/` to `Downloads/`
Using FileZilla, transfer `1.17.1/` to `~/immcantation-python3.8.12/imm-env/bin/1.17.1/`

In the home directory on O2:

`cd immcantation-python3.8.12/imm-env/bin/1.17.1/`

`tar -zxvf ncbi-igblast-1.17.1-x64-linux.tar.gz`

Remove all other files using:
`rm -r ncbi-igblast-1.17.1-*`

`rm ChangeLog`

Place the folder in imm-env/ using:
`mv ncbi-igblast-1.17.1/ ../ncbi-igblast-1.17.1`

`cd ..`

`rm -r 1.17.1/`

### For IgPhyML

`module load git/2.9.5`

`git clone https://bitbucket.org/kleinstein/igphyml`

`cd igphyml/`

`./configure --prefix=/home/suk571/immcantation-python3.8.12/imm-env/bin/`

`make`

`make check`

`make install`

`make installcheck`

`./make_phyml_omp`

`export PATH=/home/suk571/immcantation-python3.8.12/imm-env/bin/igphyml/src:$PATH`

### Installing other Immcantation packages

Alakazam, SHazaM, TIgGER, SCOPer, prestoR, RDI are on CRAN.
Except prestoR. Installation: https://immcantation.readthedocs.io/en/stable/packages/prestor.html
RAbHIT is only available on R-4.1+. I’m using 4.0.1

### Installing PHYLIP

`wget http://evolution.gs.washington.edu/phylip/download/phylip-3.697.tar.gz`

`tar -zxvf phylip-3.697.tar.gz`

`rm -r *.gz`

`cd phylip-3.697/`

`cd src`

`make -f Makefile.unx install`

`cd`

`nano ~/.bashrc`

Add the following lines in bashrc:
`# Run path for using PHYLIP`
`export LD_LIBRARY_PATH=/home/suk571/immcantation-python3.8.12/imm-env/bin/phylip-3.697/exe:$LOAD_LIBRARY_PATH`

Exit `nano` and run:
`source ~/.bashrc`

### Adding Immcantation Docker pipeline scripts to `~/immcantation-python3.8.12/imm-env/bin`

Taken from https://bitbucket.org/kleinstein/immcantation/src/c79aba1f73470e3acbefde8cdb19f27d28fc9752/pipelines/

`cd immcantation-python3.8.12/imm-env/bin/`

Using `nano`, copy all scripts



## Setting Up Databases for IgBLAST

Use the `fetch_imgtdb.sh` script to get the human database. Delete all other databases that get downloaded.
OR Manually download IMGT reference sequences from http://www.imgt.org/vquest/refseqh.html and move via FileZilla.
Saved in folder `~/immcantation-python3.8.12/imm-env/imgt/human/vdj/`

Run the following in Linux:

`cd ~/immcantation-python3.8.12/imm-env/imgt/human/vdj`

Group heavy and light chain files together:
`cat imgt_human_IGHV.fasta imgt_human_IGKV.fasta imgt_human_IGLV.fasta > imgt_human_IGV.fasta`

`cat imgt_human_IGHJ.fasta imgt_human_IGKJ.fasta imgt_human_IGLJ.fasta > imgt_human_IGJ.fasta`

`cd ../../../`

**Ensure that Python environment is activated before running these lines to setup the database for Ig chains:**

Make directories to store the output databases (while in `imm-env/`):

`mkdir db_output`

`mkdir db_output/imgt_human_ig_v`

`mkdir db_output/imgt_human_ig_d`

`mkdir db_output/imgt_human_ig_j`

For V genes:

`bin/clean_imgtdb.py imgt/human/vdj/imgt_human_IGV.fasta db_output/imgt_human_ig_v/imgt_human_IGV.fasta`

**This step is important – it removes gaps from IMGT reference sequences to build the IgBLAST database.**

`bin/ncbi-igblast-1.17.1/bin/makeblastdb -parse_seqids -dbtype nucl -in db_output/imgt_human_ig_v/imgt_human_IGV.fasta`

For D genes:

`bin/clean_imgtdb.py imgt/human/vdj/imgt_human_IGHD.fasta db_output/imgt_human_ig_d/imgt_human_IGD.fasta`

`bin/ncbi-igblast-1.17.1/bin/makeblastdb -parse_seqids -dbtype nucl -in db_output/imgt_human_ig_d/imgt_human_IGD.fasta`

For J genes:

`bin/clean_imgtdb.py imgt/human/vdj/imgt_human_IGJ.fasta db_output/imgt_human_ig_j/imgt_human_IGJ.fasta`

`bin/ncbi-igblast-1.17.1/bin/makeblastdb -parse_seqids -dbtype nucl -in db_output/imgt_human_ig_j/imgt_human_IGJ.fasta`

**Still within imm-env/:**

`mv db_output/ bin/ncbi-igblast-1.17.1/bin/`

`cd bin/ncbi-igblast-1.17.1/bin/`

`mv db_output/ database/`

Also move all folders to bin so that the directory structure looks like (showing all folders and select files):

bin

├── blastdbcmd

├── database

│   ├── imgt_human_ig_d

│   │   ├── imgt_human_IGD.fasta

│   │   ├── imgt_human_IGD.fasta.ndb

│   │   ├── imgt_human_IGD.fasta.nhr

│   │   ├── imgt_human_IGD.fasta.nin

│   │   ├── imgt_human_IGD.fasta.nog

│   │   ├── imgt_human_IGD.fasta.nos

│   │   ├── imgt_human_IGD.fasta.not

│   │   ├── imgt_human_IGD.fasta.nsq

│   │   ├── imgt_human_IGD.fasta.ntf

│   │   └── imgt_human_IGD.fasta.nto

│   ├── imgt_human_ig_j

│   └── imgt_human_ig_v

│
├── igblastn

├── igblastp

├── internal_data

│
├── makeblastdb

└── optional_file


**Before running IgBLAST, export IGDATA to path to bin folder in bashrc.**

`nano ~/.bashrc`

`# Path to IgBLAST`
`export IGDATA=/home/suk571/immcantation-python3.8.12/imm-env/bin/ncbi-igblast-1.17.1/bin`

Exit and run:
`source ~/.bashrc`
