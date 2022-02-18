# MANTA_vcf2bedpe

The vcf to bedpe workflow is desgined to prepare the Manta calls for SV annotation. During the preperation a few predetermined filtering steps are applied.
	
	1. Manta calls that are less than 50bp are removed. 
	
	2. Manta calls labeled as 'IMPRECISE' by Manta are removed. These lack additional metadata such as homology that impact downstream analysis. These calls are also ones that may lack a precise breakpoint location.
	
	3. Manta calls that align to any chromosome other than 1-22, X, and Y are removed. 

These Manta calls are written out in the same output directory in a minimally processed file for further investigation if desired. This file is labeled as the sample name with {%.removed_calls} ending.

This code is an R implementation of the svtools vcf2bedpe function which can also be used. The filters applied differ between the two functions.

## Example Usage
```bash
Rscript MANTA_vcf2bedpe.R -i <path to vcf file> -o <output_directory_path/>
```


# MantaSVAnnotator
Takes Manta bedpe from either the MANTA_vcf2bedpe function or from the svtools. Annotates each breakpoint to determine if it is in a gencode identified region. 
This function also outputs genes that are present in the TAD the SV occurs in.

Uses fuzzy filtering based on gnomad germline SVs to determine somatic events. {%.sv.annotated.bedpe} contains both germline and somatic annotated events.  {%.somatic_only_sv.annotated.bedpe} contains all SV annotations that were not within 100bp of a perfect match in the gnomad germline SV reference or have a span greater than 1000bp.

use `Manta_SV_Annotator_2` (REQUIRED inputs are input fiile path, output directory, gene annotation file, and germline annotation file)

use `hg38_ensembl_genelocations_formatted.txt` as the gene annotation file

use `gnomad_germline_hg38all.txt` as the germline annotation file

## Example usage

```bash
Rscript Manta_SV_Annotator_2.R 
-i <Required:input input_filepath> 
-o <Required:output directory_path/> 
-r <Required:gene annotation_filepath> 
-g <Required:germline reference_filepath> 
-c <Optional:cores (default = 1)>
```

## Reference files 
Ensembl gene locations were downloaded from biomart. Formatting restricted the table to gene boundaries, chromosome, and gene ID to reduce file size. 

Gnomad germline SV reference files were downloaded from the gnomad project public database (controls only). They were then restricted to those that 'PASSED' gnomad's filtering process. The rtracklayer implimentation of liftOver was used to translate the hg19 positions to hg38. Insertions had endpoints adjusted to reflect the size of the insertion in relation to the reference genome.
