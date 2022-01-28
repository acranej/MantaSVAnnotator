# VCF_to Bedpe Function

The vcf to bedpe workflow is desgined to prepare the Manta calls for SV annotation. During the preperation a few predetermined filtering steps are applied.
	1. Manta calls that are less than 1kb are removed. The removed calls are mostly indels and potential germline calls. 
	2. Manta calls labeled as 'IMPRECISE' by Manta are removed. These lack additional metadata such as homology that impact downstream analysis. These calls are also ones that may lack a precise breakpoint location.
	3. Manta calls that align to any chromosome other than 1-22 and X are removed. 

These Manta calls are written out in the same output directory in a minimally processed file for further investigation if desired. This file is labeled as the sample name with '.removed_calls' ending

## Example Usage
```bash
Rscript MANTA_vcf2bedpe.R -i <path to vcf file> -o <output directory path>
```


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
