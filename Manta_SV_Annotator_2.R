#### import libraries and packages
cat("Loading packages...\n")
time_begin <- Sys.time()
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
#' @param bedpe_inp: Manta Bedpe returned from SVtools or MANTA_vcf2bedpe
#' @return SV data table with columns added indicating if the SV affects a gene, the name of the gene, and the type of the gene
#' @description Determines if each breakend of a SV occurs within a gene's boundaries
#' @export
#'

annotate_sv = function(i, bedpe_inp) {
  
  bedpe <- bedpe_inp[i,]
  
  ### Annotate break point A by giving genes that exist where the breakpoint occurs
  temp_A_gr <- GRanges(bedpe$CHROM_A[1], 
                         IRanges(bedpe$START_A[1], bedpe$END_A[1]))
  temp_A_gr_ant <-  suppressWarnings(gene_locations %&% temp_A_gr)
  temp_A_dt_ant <- as.data.table(temp_A_gr_ant)
    
  ### collapse all genes and gene types that occur in breakpoint A
  if(!nrow(temp_A_dt_ant) == 0) {
    if(nrow(temp_A_dt_ant) == 1) {
        
      bedpe$BREAK_A_Ensembl_Gene[1] <- temp_A_dt_ant$EnsemblID[1]
      bedpe$BREAK_A_Gene_Name[1] <- temp_A_dt_ant$GeneName[1]
        
    } else {
        
      bedpe$BREAK_A_Ensembl_Gene[1] <- paste0(temp_A_dt_ant$EnsemblID, collapse = ";")
      bedpe$BREAK_A_Gene_Name[1] <- paste0(temp_A_dt_ant$GeneName, collapse = ";")
        
    }
  } else {
      
    bedpe$BREAK_A_Ensembl_Gene[1] <- ""
    bedpe$BREAK_A_Gene_Name[1] <- ""
      
  }
    
  ### Annotate break point B
  temp_B_gr <- GRanges(bedpe$CHROM_B[1], 
                         IRanges(bedpe$START_B[1], bedpe$END_B[1]))
  temp_B_gr_ant <-  suppressWarnings(gene_locations %&% temp_B_gr)
  temp_B_dt_ant <- as.data.table(temp_B_gr_ant)
  
  ### collapse all genes and gene types that occur in breakpoint B 
  if(!nrow(temp_B_dt_ant) == 0) {
    if(nrow(temp_B_dt_ant) == 1) {
        
      bedpe$BREAK_B_Ensembl_Gene[1] <- temp_B_dt_ant$EnsemblID[1]
      bedpe$BREAK_B_Gene_Name[1] <- temp_B_dt_ant$GeneName[1]
        
    } else {
        
      bedpe$BREAK_B_Ensembl_Gene[1] <- paste0(temp_B_dt_ant$EnsemblID, collapse = ";")
      bedpe$BREAK_B_Gene_Name[1] <- paste0(temp_B_dt_ant$GeneName, collapse = ";")
        
    }
  } else {
      
    bedpe$BREAK_B_Ensembl_Gene[1] <- ""
    bedpe$BREAK_B_Gene_Name[1] <- ""
      
  }
return(bedpe)
}

#' 
#' Germline filtering
#' @name fuzzy_filter_germline
#' 
#' @param i: passed from lapply to iterate
#' @param bedpe: Manta Bedpe returned from annotate_sv function
#' @return SV data table with columns added indicating germline or somatic, germline is defined as <=200bp away from agnostic perfect match in reference
#' @description Determines if each SV should be considered germline by hard filtering
#' @export
#'
#'
fuzzy_filter_germline <- function(i, bed) {
 
   sub <- bed[i,]
  
  ## reorder for filtering
  if(sub$CHROM_A > sub$CHROM_B) {
    sub_ord <- cbind(CHROM_A=sub$CHROM_B, START_A=sub$START_B, END_A=sub$END_B, CHROM_B=sub$CHROM_A, START_B=sub$START_A, END_B=sub$END_A, sub[,7:ncol(bed)])
  } else if (sub$CHROM_A == sub$CHROM_B & sub$START_A > sub$START_B) {
    sub_ord <- cbind(CHROM_A=sub$CHROM_B, START_A=sub$START_B, END_A=sub$END_B, CHROM_B=sub$CHROM_A, START_B=sub$START_A, END_B=sub$END_A, sub[,7:ncol(bed)])
  } else {
    sub_ord <- sub
  }
  
  ### change to integers to match reference germline
  sub_ord[CHROM_A == "X", CHROM_A := 23]
  sub_ord[CHROM_B == "X", CHROM_B := 23]
  sub_ord[CHROM_A == "Y", CHROM_A := 24]
  sub_ord[CHROM_B == "Y", CHROM_B := 24]
  
  ### subset reference to matching chromosome
  ref_sub <- hg38_germline_gnomad[chrom1 == sub_ord$CHROM_A & chrom2 == sub_ord$CHROM_B]
  
  ### calculate distances
  ref_sub[,str_dist := abs(start - sub_ord$START_A)]
  ref_sub[,end_dist := abs(end - sub_ord$START_B)]
  ref_sub[, tot_dist := (str_dist + end_dist)]
  
  ### choose closest match
  ref_min <- ref_sub[which.min(ref_sub$tot_dist)]
  
  if (nrow(ref_min) < 1) {
    sub$Filter[1] <- "Somatic"
  
  } else if (ref_min$tot_dist > 1000 & sub$SPAN > 0 & sub$SPAN < 1000) {
    sub$Filter[1] <- "Germline"
  }
  else if (ref_min$tot_dist > 1000) {
    sub$Filter[1] <-paste0("Somatic","(", ref_min$tot_dist,")")
    
  } else {
    sub$Filter[1] <-paste0("Germline","(", ref_min$tot_dist,")")
    
  }
  return(sub)
}

#' 
#' TAD annotation
#' @name annotate_tad
#' 
#' @param i: passed from lapply to iterate
#' @param bedpe: Manta Bedpe returned from annotate_sv function
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
    ### some breakpoints ooccur outside a TAD boundary
    if(!(length(tad_gr_in_A) + length(tad_gr_in_B)) == 2) {
        
      rearrangement$TAD_CORD[1] <- ""
      rearrangement$ONCO_TAD_GENE[1] <- ""
      rearrangement$TAD_GENE[1] <- ""
        
    } else {
      ### genes in the TAD of breakpoint A
      genes_in_tad_A <- suppressWarnings(gene_locations %&% tad_gr_in_A)
      genes_in_tad_dt_A <- gr2dt(genes_in_tad_A)
      genes_in_tad_dt_pc_A <- genes_in_tad_dt_A[gsub(" ", "",as.character(GENE_TYPE)) == "protein_coding"]### limit to protein coding, could expand to lnRNA, miRNA, etc.
      
      ### genes in the TAD of breakpoint B
      genes_in_tad_B <- suppressWarnings(gene_locations %&% tad_gr_in_B)
      genes_in_tad_dt_B <- gr2dt(genes_in_tad_B)
      genes_in_tad_dt_pc_B <- genes_in_tad_dt_B[gsub(" ", "",as.character(GENE_TYPE)) == "protein_coding"]### limit to protein coding, could expand to lnRNA, miRNA, etc.
      
      ### distinguish onco genes from non onco genes --> lots of genes in a TAD, onco genes are more interesting
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
      
      ### only annotate oncogenes if one of the breakpoints occurs in a TAD with an onco gene
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
  
  make_option(c("-r", "--genomeant"),  type = "character", default = NULL,  help = "Path to gene location file"),
  
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output directory path"),
  
  make_option(c("-g","--germline"), type = "character", default = NULL, help = "Germline reference path"),
  
  make_option(c("-c", "--cores"), type = "integer", default = 1, help = "Number of cores to run on")
  
  ### to be included in V2
  #make_option(c("-c", "--cancergenes"),  type = "character", default = NULL,  help = "Cancer Genes: path to cosmic cancer genes"),
  
  #make_option(c("-t", "--tadboundaries"),  type = "character", default = NULL,  help = "TAD Boundaries: path to TAD boundaries")
  
)

# inputs
parseobj = OptionParser(option_list = option_list)
opt = parse_args(parseobj)

if (!file.exists(opt$input)) {
  
  stop(sprintf("Input file '%s' does not exist!", input_pth))
  
} else if (is.null(opt$genomeant) | is.null(opt$output) | is.null(opt$germline)) {
  
  stop(sprintf("Please include reference files and output directory"))
  
} else {
  input_pth = paste0(opt$input)
  sample <- unlist(strsplit(input_pth, "/"))[length(unlist(strsplit(input_pth, "/")))]
  sample <- gsub(".vcf.bedpe","", sample)
  out_pth = paste0(opt$output)
  
  cat('Building references... \n')
  
  # get reference gene annotations and germline
  hg38_annotgenes = fread(paste0(opt$genomeant))
  hg38_germline_gnomad <<- fread(paste0(opt$germline))
  
  ### build granges for reference
  gene_locations <<- GRanges(hg38_annotgenes$chromosome, IRanges(hg38_annotgenes$start, 
                                                           hg38_annotgenes$end), EnsemblID = hg38_annotgenes$Ensembl_ID, GeneName = hg38_annotgenes$gene_name)
  
  cat("Reading file...\n")
  bedpe_inp <- fread(paste0(input_pth))
  
  cat("Adding gene annotations...\n")
  bedpe_gene_annotation <- rbindlist(mclapply(1:nrow(bedpe_inp), annotate_sv, bedpe_inp, mc.cores = opt$cores))
  
  cat("Fuzzy filtering germline...\n")
  bedpe_fuzzy_filtered <- rbindlist(mclapply(1:nrow(bedpe_gene_annotation), fuzzy_filter_germline, bedpe_gene_annotation, mc.cores = opt$cores))
  
  bedpe_somatic_only <- bedpe_fuzzy_filtered[grep("Somatic", Filter)]
  output_somatic_only <- paste0(out_pth, sample, ".somatic_only_sv.annotated.bedpe")
  output_all <- paste0(out_pth, sample, ".sv.annotated.bedpe")
  
  write.table(bedpe_fuzzy_filtered, output_all, sep = '\t', row.names = F, col.names = T, quote = F)
  write.table(bedpe_somatic_only, output_somatic_only, sep = '\t', row.names = F, col.names = T, quote = F)
  time_end <- Sys.time()
  cat(paste0("Began at ", time_begin,"\n"))
  cat(paste0("Ended at ", time_end,"\n"))
  
}
