# BCR Lineage Tracing Using Immcantation

Before you begin, follow setup instructions in `Setup.md`. Once complete, follow the steps outlined below.

### Using the virtual python3 environment on O2

Be present in the home directory in a compute node.

`module load gcc/9.2.0`

`module load python/3.8.12`

`source immcantation-python3.8.12/imm-env/bin/activate`

To deactivate the environment:
`deactivate`


## Step 1: Running IgBLAST

Running IgBLAST to convert 10X VDJ FASTA files into AIRR Rearrangement TSV format.

`AssignGenes.py igblast` runs IgBLAST but it is recommended to use IgBLAST directly.
In my data, I have 4 patients and 4 healthy controls, and 4 replicates of each sample. CITE-seq hashing labels are as follows:
Healthy - 1, 3, 5, 6
Patient - 2, 4, 7, 8

Before running IgBLAST,change all duplicate contig names in the FASTA and annotation files:
Barcodes should match the first part of the contig name. After the `barcodes` column comes the `is_cell` column and then the `contigs` column.
Changes must be made within the R program that generates these files, `lineage_tracing_bcr_annotations.R` for the annotations file followed by `lineage_tracing_bcr_contigs.R` for the FASTA files.

Now run IgBLAST:

```
/home/suk571/immcantation-python3.8.12/imm-env/bin/ncbi-igblast-1.17.1/bin/igblastn \
-query /home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/B_VDJ/patient_v2/8_filtered_contig.fasta \
-out /home/suk571/immcantation-python3.8.12/imm-env/igblast_output/patient/8_filtered_contig_igblast.fmt7 \
-num_threads 4 \
-ig_seqtype Ig \
-organism human \
-auxiliary_data /home/suk571/immcantation-python3.8.12/imm-env/bin/ncbi-igblast-1.17.1/bin/optional_file/human_gl.aux \
-germline_db_V /home/suk571/immcantation-python3.8.12/imm-env/bin/ncbi-igblast-1.17.1/bin/database/imgt_human_ig_v/imgt_human_IGV.fasta \
-germline_db_D /home/suk571/immcantation-python3.8.12/imm-env/bin/ncbi-igblast-1.17.1/bin/database/imgt_human_ig_d/imgt_human_IGD.fasta \
-germline_db_J /home/suk571/immcantation-python3.8.12/imm-env/bin/ncbi-igblast-1.17.1/bin/database/imgt_human_ig_j/imgt_human_IGJ.fasta \
-outfmt "7 std qseq sseq btop" \
-domain_system imgt
```
>Note: In the code, I have named the folder `healthy` and `patient` but on the console I'm using folder names `healthy_v2` and `patient_v2` since I created new files using the code mentioned above and put them there. The `v1` folders have files generated using an old code that is not present on this git.


