#### import libraries and packages
cat("Loading packages...\n")
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(GenomicRanges))
suppressPackageStartupMessages(require(gUtils))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(library(parallel))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
cat("...done.\n")

#' 
#' Breakend annotation
#' @name annotate_sv
#' 
#' @param i: passed from lapply to iterate
#' @param bedpe_inp: Manta Bedpe returned from SVtools
#' @return SV data table with columns added indicating if the SV affects a gene, the name of the gene, and the type of the gene
#' @description Determines if each breakend of a SV occurs within a gene's boundaries
#' @export
#'
annotate_sv = function(i, bedpe_inp) {
  
  ### build granges for reference
  hg38_annotgenes <- hg38_gene_gr_reference
  
  gene_locations <- GRanges(hg38_annotgenes$chrom, IRanges(hg38_annotgenes$start_pos, 
                                                    hg38_annotgenes$end_pos))
  mcols(gene_locations) <- hg38_annotgenes[,c(2,4:10)]
  
  
  bedpe <- bedpe_inp[i,]
  
  ### Annotate break point A by giving genes that exist where the breakpoint occurs
  temp_A_gr <- GRanges(paste0("chr",bedpe$CHROM_A[1]), 
                         IRanges(bedpe$START_A[1], bedpe$END_A[1]))
  temp_A_gr_ant <-  suppressWarnings(gene_locations %&% temp_A_gr)
  temp_A_dt_ant <- as.data.table(temp_A_gr_ant)
    
  ### decision table for 
  if(!nrow(temp_A_dt_ant) == 0) {
    if(nrow(temp_A_dt_ant) == 1) {
        
      bedpe$BREAK_A_GENE[1] <- temp_A_dt_ant$GENE_NAME[1]
      bedpe$BREAK_A_GENE_TYPE[1] <- temp_A_dt_ant$GENE_TYPE[1]
        
    } else {
        
      bedpe$BREAK_A_GENE[1] <- paste0(temp_A_dt_ant$GENE_NAME, collapse = ";")
      bedpe$BREAK_A_GENE_TYPE[1] <- paste0(temp_A_dt_ant$GENE_TYPE, collapse = ";")
        
    }
  } else {
      
    bedpe$BREAK_A_GENE[1] <- ""
    bedpe$BREAK_A_GENE_TYPE[1] <- ""
      
  }
    
  ### Annotate break point B
  temp_B_gr <- GRanges(paste0("chr",bedpe$CHROM_B[1]), 
                         IRanges(bedpe$START_B[1], bedpe$END_B[1]))
  temp_B_gr_ant <-  suppressWarnings(gene_locations %&% temp_B_gr)
  temp_B_dt_ant <- as.data.table(temp_B_gr_ant)
    
  if(!nrow(temp_B_dt_ant) == 0) {
    if(nrow(temp_B_dt_ant) == 1) {
        
      bedpe$BREAK_B_GENE[1] <- temp_B_dt_ant$GENE_NAME[1]
      bedpe$BREAK_B_GENE_TYPE[1] <- temp_B_dt_ant$GENE_TYPE[1]
        
    } else {
        
      bedpe$BREAK_B_GENE[1] <- paste0(temp_B_dt_ant$GENE_NAME, collapse = ";")
      bedpe$BREAK_B_GENE_TYPE[1] <- paste0(temp_B_dt_ant$GENE_TYPE, collapse = ";")
        
    }
  } else {
      
    bedpe$BREAK_B_GENE[1] <- ""
    bedpe$BREAK_B_GENE_TYPE[1] <- ""
      
  }
return(bedpe)
}

#' 
#' TAD annotation
#' @name annotate_tad
#' 
#' @param i: passed from lapply to iterate
#' @param bedpe: Manta Bedpe returned from annotate_sv
#' @return SV data table with columns added indicating if the SV occurs in a tad, oncogenes in the TAD, and non-oncogenes in the TAD
#' @description Determines if each breakend of a SV occurs within a TAD's boundaries
#' @export
#'

annotate_tad = function(i, bedpe) {
  tad_boundaries <- hg38_tad_gr_reference
  
  tad_gr <- GRanges(tad_boundaries$chrom, IRanges(tad_boundaries$sttart_pos, 
                                                    tad_boundaries$end_pos))
  
  hg38_annotgenes <- hg38_gene_gr_reference
  
  gene_locations <- GRanges(hg38_annotgenes$chrom, IRanges(hg38_annotgenes$start_pos, 
                                                           hg38_annotgenes$end_pos))
  mcols(gene_locations) <- hg38_annotgenes[,c(2,4:10)]
  
  cosmic.genes <- cancer_genes
  
  rearrangement <- bedpe[i,]
    
   # intra SVs can be treated as their own range
  if(rearrangement$CHROM_A == rearrangement$CHROM_B) {
      
    rearrangement_gr <- GRanges(paste0("chr", rearrangement$CHROM_A), 
                                  IRanges(as.integer(mean(c(rearrangement$START_A, rearrangement$END_A))), 
                                          as.integer(mean(c(rearrangement$START_B, rearrangement$END_B)))))
      
    ### gets TADs affected by the rearrangement
    tad_gr_in <- suppressWarnings(tad_gr %&% rearrangement_gr)
      
    ### some SV locations may not be in TADs (maybe)
    if(length(tad_gr_in) > 0) {
        
      ### gets genes in tad
      genes_in_tad <- suppressWarnings(gene_locations %&% tad_gr_in)
      genes_in_tad_dt <- gr2dt(genes_in_tad)
      genes_in_tad_dt_pc <- genes_in_tad_dt[gsub(" ", "", as.character(GENE_TYPE)) == "protein_coding"]### limit to protein coding, could expand to lnRNA, miRNA, etc.
        
      ### identifies known oncogenic and non oncogenic genes in tad
      onco.tmp <- gsub(" ","", as.character(genes_in_tad_dt_pc[gsub(" ", "",as.character(GENE_NAME)) %in% cosmic.genes]$GENE_NAME))
      non_onco.tmp <- gsub(" ","", as.character(genes_in_tad_dt_pc[!gsub(" ", "",as.character(GENE_NAME)) %in% cosmic.genes]$GENE_NAME))
        
      # give TAD locations
      tad_gr_in_dt <- gr2dt(tad_gr_in)
      tad_gr_in_dt[, TAD_NAME := paste0(seqnames, "[", start, "-", end, "]")]
        
      ### adds to bedpe
      rearrangement$TAD_CORD[1] <- paste0(tad_gr_in_dt$TAD_NAME, collapse = ";")
      rearrangement$ONCO_TAD_GENE[1] <- paste0(onco.tmp, collapse = ";")
      rearrangement$TAD_GENE[1] <- paste0(non_onco.tmp, collapse = ";")
        
        
        
    } else {
        
      rearrangement$TAD_CORD[1] <- ""
      rearrangement$ONCO_TAD_GENE[1] <- ""
      rearrangement$TAD_GENE[1] <- ""
        
    }
  } else {
      
    ### if the SV is a BND, each point must be treated individually
    rearrangement_gr_A <- GRanges(paste0("chr", rearrangement$CHROM_A), 
                                    IRanges(rearrangement$START_A, rearrangement$END_A))
    rearrangement_gr_B <- GRanges(paste0("chr", rearrangement$CHROM_B), 
                                    IRanges(rearrangement$START_B, rearrangement$END_B))
      
    ### gets TADs affected by the rearrangement
    tad_gr_in_A <- suppressWarnings(tad_gr %&% rearrangement_gr_A) 
    tad_gr_in_B <- suppressWarnings(tad_gr %&% rearrangement_gr_B) 
      
    ### this should only return 2 Tads, for now do not annotate if not 2
    if(!(length(tad_gr_in_A) + length(tad_gr_in_B)) == 2) {
        
      rearrangement$TAD_CORD[1] <- ""
      rearrangement$ONCO_TAD_GENE[1] <- ""
      rearrangement$TAD_GENE[1] <- ""
        
    } else {
        
      genes_in_tad_A <- suppressWarnings(gene_locations %&% tad_gr_in_A)
      genes_in_tad_dt_A <- gr2dt(genes_in_tad_A)
      genes_in_tad_dt_pc_A <- genes_in_tad_dt_A[gsub(" ", "",as.character(GENE_TYPE)) == "protein_coding"]### limit to protein coding, could expand to lnRNA, miRNA, etc.
        
      genes_in_tad_B <- suppressWarnings(gene_locations %&% tad_gr_in_B)
      genes_in_tad_dt_B <- gr2dt(genes_in_tad_B)
      genes_in_tad_dt_pc_B <- genes_in_tad_dt_B[gsub(" ", "",as.character(GENE_TYPE)) == "protein_coding"]### limit to protein coding, could expand to lnRNA, miRNA, etc.
        
      onco.tmp_A <- gsub(" ","", as.character(genes_in_tad_dt_pc_A[gsub(" ", "",as.character(GENE_NAME)) %in% cosmic.genes]$GENE_NAME))
      non_onco.tmp_A <- gsub(" ","", as.character(genes_in_tad_dt_pc_A[!gsub(" ", "",as.character(GENE_NAME)) %in% cosmic.genes]$GENE_NAME))
        
      onco.tmp_B <- gsub(" ","", as.character(genes_in_tad_dt_pc_B[gsub(" ", "",as.character(GENE_NAME)) %in% cosmic.genes]$GENE_NAME))
      non_onco.tmp_B <- gsub(" ","", as.character(genes_in_tad_dt_pc_B[!gsub(" ", "",as.character(GENE_NAME)) %in% cosmic.genes]$GENE_NAME))
        
      tad_gr_in_dt_A <- gr2dt(tad_gr_in_A)
      tad_gr_in_dt_A[, TAD_NAME := paste0("A","|",seqnames, "[", start, "-", end, "]")]
        
      tad_gr_in_dt_B <- gr2dt(tad_gr_in_B)
      tad_gr_in_dt_B[, TAD_NAME := paste0("B","|",seqnames, "[", start, "-", end, "]")]
        
      t_name <- paste0(tad_gr_in_dt_A$TAD_NAME, ";", tad_gr_in_dt_B$TAD_NAME)
      rearrangement$TAD_CORD[1] <- t_name
        
      if((length(onco.tmp_A) + length(onco.tmp_B)) > 0) {
          
        rearrangement$ONCO_TAD_GENE[1] <- paste0("A","{", paste0(onco.tmp_A, collapse = ";"), "}", "|","B","{", paste0(onco.tmp_B, collapse = ";"),"}")
          
      } else {
          
        rearrangement$ONCO_TAD_GENE[1] <- ""
          
      }
        
      if((length(non_onco.tmp_A) + length(non_onco.tmp_B)) > 0) {
          
        rearrangement$TAD_GENE[1] <- paste0("A","{", paste0(non_onco.tmp_A, collapse = ";"), "}", "|","B","{", paste0(non_onco.tmp_B, collapse = ";"),"}")
          
      } else {
          
        rearrangement$TAD_GENE[1] <- ""
          
      }
    }
  }
return(rearrangement)
}




#### MAIN FUNCTION #####

option_list <- list(
  
  make_option(c("-i", "--input"),  type = "character", default = NULL,  help = "Input bedpe directory path"),
  
  make_option(c("-r", "--genomeant"),  type = "character", default = NULL,  help = "Gencode Annotations: path to gencode annotatioon file"),
  
  make_option(c("-c", "--cancergenes"),  type = "character", default = NULL,  help = "Cancer Genes: path to cosmic cancer genes"),
  
  make_option(c("-t", "--tadboundaries"),  type = "character", default = NULL,  help = "TAD Boundaries: path to TAD boundaries")
  
)

# inputs
parseobj = OptionParser(option_list = option_list)
opt = parse_args(parseobj)

input_pth = paste0(opt$input)
output_pth = sprintf("%s.sv_annotation", input_pth)


# get reference gene annotations and tad locations as granges
hg38_gene_gr_reference = fread(paste0(opt$genomeant))
hg38_tad_gr_reference = fread(paste0(opt$tadboundaries))
cancer_genes <- fread(paste0(opt$cancergenes))$'Gene Symbol'

if (!file.exists(input_pth)) {
  
  stop(sprintf("Input file '%s' does not exist!", input_pth))
  
} else if (is.null(opt$genomeant) | is.null(opt$tadboundaries) | is.null(opt$cancergenes)) {
  
  stop(sprintf("Please include reference files"))
  
} else {
  
  cat("Reading file...\n")
  bedpe_inp <- fread(paste0(input_pth))
  
  cat("Adding gene annotations...\n")
  bedpe_gene_annotation <- rbindlist(lapply(1:nrow(bedpe_inp), annotate_sv, bedpe_inp))
  cat("Adding TAD annotations...\n")
  bedpe_gene_tad_annotation <- rbindlist(lapply(1:nrow(bedpe_gene_annotation), annotate_tad, bedpe_gene_annotation))
  
  write.table(bedpe_gene_tad_annotation, output_pth, row.names = F, col.names = T, sep = "\t", quote = F)
}
  
