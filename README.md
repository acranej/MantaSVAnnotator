# MantaSVAnnotator
Takes Manta VCFs or BEDPE file formats and annotates the structural variants

use `Manta_SV_Annotator_2`

## Example usage

```bash
Rscript /.../Manta_SV_Annotation/Manta_SV_Annotator_2.R 
-i "/.../Manta_Bedpe_SVtools/Filtered/CDS-0b4jFH.somaticSV_filtered.bedpe" 
-r "/.../Manta_SV_Annotation/inputfiles/gencode_hg38_annotations_table.txt" 
-c "/.../Manta_SV_Annotation/inputfiles/Census_cancer_genes_allTue_May_2018.tsv" 
-t "/.../Manta_SV_Annotation/inputfiles/imr90fibroblast_tad_boundaries.txt"
```
