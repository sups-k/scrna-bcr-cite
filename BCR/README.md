# BCR Lineage Tracing Using Immcantation (v.4.3.0)

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


## Step 2: Create tab-delimited database file to store sequence alignment information

For this step, use the original IMGT reference files (before gaps were removed using `clean_imgtdb.py` (refer `Setup.md`)) present in `imgt/human/vdj/`. Create a new folder `makedb_output` within `imm-env` to store the results.
Within `imm-env/`:

```
MakeDb.py igblast \
-i igblast_output/patient/8_filtered_contig_igblast.fmt7 \
-s /home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/B_VDJ/patient_v2/8_filtered_contig.fasta \
-r imgt/human/vdj/imgt_human_IGV.fasta imgt/human/vdj/imgt_human_IGHD.fasta imgt/human/vdj/imgt_human_IGJ.fasta \
-o makedb_output/8_filtered_contig_igblast_db-pass.tab \
--10x /home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/B_VDJ/patient_v2/8_filtered_contig_annotations.csv \
--format changeo \
--extended
```


## Step 3: Filter database to only keep functional sequences

Withim `imm-env/`:

```
ParseDb.py select -d makedb_output/8_filtered_contig_igblast_db-pass.tab \
-f FUNCTIONAL \
-u T \
--outname 8_functional
```

Move all output files to `makedb_output/functional`

`mv *parse-select.tab functional/`

## Step 4: Parse database output into heavy and light chain ChangeO files

Parse output into heavy chain ChangeO files:
```
ParseDb.py select -d makedb_output/functional/8_functional_parse-select.tab \
-f LOCUS \
-u "IGH" \
--logic all \
--regex \
--outname 8_heavy
```
Parse output into light chain ChangeO files:
```
ParseDb.py select -d makedb_output/functional/8_functional_parse-select.tab \
-f LOCUS \
-u "IG[LK]" \
--logic all \
--regex \
--outname 8_light
```

Move light chain files into `makedb_output/functional/light` and heavy chain files into `makedb_output/functional/heavy`.
```
mv *heavy*.tab heavy/
mv *light*.tab light/
```

## Step 5: Obtain Predicted Thresholds for Step 6

Run `01_VDJ_genotype_and_threshold-bimodal.R` which results in a `predicted_thresholds.csv` file if the distribution of thresholds is bimodal. If it is bimodal, **do not** run SCOPer. If the distribution is unimodal, run `01_VDJ_genotype_and_threshold-unimodal.R`.

## Step 6: Split BCR sequences by clone

### Did not run SCOPer
Get threshold values from `~/immcantation-python3.8.12/imm-env/makedb_output/functional/predicted_thresholds.csv` if the distribution of thresholds is bimodal and you didn’t have to run SCOPer.

Run the following block of code for each sample, i.e., IGHV-genotyped_Hashtag1.tab to IGHV-genotyped_Hashtag8.tab. Enter the corresponding threshold values from `predicted_threshold.csv` in the `dist` argument. For every sample, change these arguments - `d`, `dist`, `outname`, `log`.
For example, for sample 1, the `outname` is `IGHV-genotyped_healthy_1` and `log` is `healthy1_DefineClones.log`.

```
DefineClones.py -d ~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/IGHV-genotyped_Hashtag8.tab \
--act set \
--model ham \
--norm len \
--dist 0.195959 \
--format changeo \
--outname IGHV-genotyped_patient_8 \
--log ~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/logs/patient8_DefineClones.log \
--nproc 1
```
Split heavy chain clones with different light chains. Run this block for each sample. For example, for sample 8, the `d` is `IGHV-genotyped_patient_8_clone-pass.tab` (the output from `DefineClones.py`), the `e` is `8_light_parse-select.tab`, and the `o` is `8_10X_clone-pass.tab`.
```
light_cluster.py -d ~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/IGHV-genotyped_healthy_1_clone-pass.tab \
-e ~/immcantation-python3.8.12/imm-env/makedb_output/functional/light/1_light_parse-select.tab \
-o 1_10X_clone-pass.tab \
--format changeo
```
Reconstruct heavy chains V and J germline genes. Run this block for each sample. For example, for sample 8, the `d` is `8_10X_clone-pass.tab` (the output from `light_cluster.py`), the `outname` is `8_ph_genotyped` and the `log` is `8_CreateGermlines.log`.
```
CreateGermlines.py -d ~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/1_10X_clone-pass.tab \
-g dmask \
--cloned \
-r ~/immcantation-python3.8.12/imm-env/imgt/human/vdj/imgt_human_IGHD.fasta \
~/immcantation-python3.8.12/imm-env/imgt/human/vdj/imgt_human_IGHJ.fasta \
~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/fasta/IGHV_genotype_Hashtag8.fasta \
--vf V_CALL_GENOTYPED \
--format changeo \
--outname 1_ph_genotyped \
--log ~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/1_CreateGermlines.log
```

### Ran SCOPer

If the distribution of thresholds is not bimodal, use single cell mode in SCOPer (`01_VDJ_genotype_and_threshold-unimodal.R`).

First, run `DefineClones.py` for each sample but in the `dist` argument, enter the threshold value obtained after running SCOPer. Use the same threshold for all the samples.
Second, skip the `light_cluster.py` function above.
Third, directly reconstruct all chains V and J germline genes using the output of SCOPer:
```
CreateGermlines.py -d ~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/Hashtag8_IGV-genotyped.tab \
-g dmask \
--cloned \
-r ~/immcantation-python3.8.12/imm-env/imgt/human/vdj/imgt_human_IGHD.fasta \
~/immcantation-python3.8.12/imm-env/imgt/human/vdj/imgt_human_IGJ.fasta \
~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/fasta/IGV_genotype_Hashtag8.fasta \
--vf V_CALL_GENOTYPED \
--format changeo \
--outname 8_ph_genotyped \
--log ~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/8_CreateGermlines.log
```
Run the above block for all samples. For example, for sample 1, `d` is `Hashtag1_IGV-genotyped.tab` (output of `01_VDJ_genotype_and_threshold-unimodal.R`), `outname` is `1_ph_genotyped`, and `log` is `1_CreateGermlines.log`.

## Step 7: Run `02_VDJ_mutation_quant.R`
This step identifies mutational load using SHazaM. This results in the file `VDJseq_mutation_quant.tab` present in `~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/create_germlines/`.

You can now use the `VDJseq_mutation_quant.tab` file to perform analysis of clonal frequencies (`clone_freq.R` in the `code` folder of this git) or lineage tracing. You can also try running the `plot_mutational_freq.R` file from the `code` folder, which I got from https://github.com/angelettilab/scMouseBcellFlu. From the same website, I got the `lineage_tracing.R` code to build lineage trees using Alakazam (it works but need to figure out how to make the image better).
I'm using `plot_isotypes_bcr.R` to find the isotype frequencies in each cluster.

I have used the output file of this table in conjuction with the Seurat metadata to obtain clustering information, resulting in the file `bcr_annot.csv`. Using `bcr_annot.csv`, I used the code `plot_clonal_freq.R` (which gives the same output as the function `occupiedscRepertoire()` from the `scRepertoire` package) to show clonal expansion in each cluster in each patient.
