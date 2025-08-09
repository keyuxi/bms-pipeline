---
output:
  pdf_document: default
  html_document: default
---
# bms-pipeline

## Starting point

**Frozen data:**

Treat these data as is.

- ArchR project: `/oak/stanford/groups/wjg/kyx/data/bms/archr/BMS_diff_m`

- Seurat object: `/oak/stanford/groups/wjg/kyx/data/bms/seurat/gex_gr.rds`

- Motif data: `data/atac/motif/motif_metadata.csv`, `data/atac/motif/pwm_list.rds`

**Other read-only data:**

Data generated from the above frozen data using scripts in `workflow`,
but are read-only to downstream analysis in R Markdown notebooks.

- BPCells fragments file: `/oak/stanford/groups/wjg/kyx/data/bms/bpcells/bms_diff.bpfrags`

- Biological pseudobulks
    - PCA of these pseudobulks, including PC reduction and rotation matrices