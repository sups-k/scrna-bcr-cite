# Workflow for scRNA-seq

1. Demultiplex data using CITE-seq hashtags. Save information in metadata and combine all technical replicates into one Seurat object.
Output file: `merged_seurat.RData`

2. Run QC, filter bad samples, run QC again. Remove all BCR genes except heavy chain genes - Ighd, Ighm, Ighg1, Ighg2c, Ighg2b, Ighg3, Igha, Ighe.
Output file: `sauron_filtered.rds`

3. Split Seurat object by sample of origin. Run SCTransform to normalize, scale, and regress all objects.
Output file: `split_sauron_regressed.rds`

4. Integrate all Seurat objects so that the same cell type, regardless of sample of origin, belongs to the same cluster.
Output file: `sauron_integrated_Jan23.rds`

5. Perform PCA again based on the number of significant PCs. Check at what resolution over-clustering happens.
Output file: `sauron_integrated_Jan23_pc43.rds`

6. Once you identify any contamination, use `remove_cells.R` to specify which clusters you want to remove from the raw Seurat object, then re-run the entire analysis. Here, I've removed clusters 19 (T/NK contamination) and 20 (monocyte/macrophage contamination).
