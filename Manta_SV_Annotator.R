# Usage Rscript Manta_SV_Annotator.R 
## takes bedpe file format
LOCAL = F
#### import libraries and packages
cat("Loading packages...\n")
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(GenomicRanges))
suppressPackageStartupMessages(require(gUtils))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(parallel))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
cat("...done.\n")

########################################
# Functions to build reference locations
########################################

#### Generates hg38 GRanges list with gene as metadata --> can include hg19
build_hg38_geneloci_gr = function() {
  
  cat("Gathering hg38 gene locations...")
    if(LOCAL) {
      
      hg38_annotgenes <- readRDS('/Volumes/xchip_beroukhimlab/Alex/Manta_SV_Annotation/inputfiles/gencode_hg38_annotations_table.rds')
      
    } else {
      
      hg38_annotgenes <- readRDS('/xchip/beroukhimlab/Alex/Manta_SV_Annotation/inputfiles/gencode_hg38_annotations_table.rds')
      
    }
    

  hg38_gr <- GRanges(hg38_annotgenes$chrom, IRanges(hg38_annotgenes$start_pos, 
                                                    hg38_annotgenes$end_pos))
  mcols(hg38_gr) <- hg38_annotgenes[,c(2,4:10)]
  
  cat("done.\n")
  return(hg38_gr)
}

build_hg38_tad_gr = function() {
  cat("Gathering fibroblast TAD reference...")
  if(LOCAL) {
    
    tad_boundaries <- readRDS("/Volumes/xchip_beroukhimlab/Alex/Manta_SV_Annotation/inputfiles/imr90fibroblast_tad_boundaries.rds")
    
  } else {
    
    tad_boundaries <- readRDS("/xchip/beroukhimlab/Alex/Manta_SV_Annotation/inputfiles/imr90fibroblast_tad_boundaries.rds")
    
  }
  
  hg38_tad <- GRanges(tad_boundaries$chrom, IRanges(tad_boundaries$sttart_pos, 
                                                    tad_boundaries$end_pos))
  cat("done.\n")
  
  return(hg38_tad)
}

###################################
# Functions to convert vcf to bedpe 
###################################

vcf_to_bedpe <- function(path) {
  cat("Converting VCF...\n")
  
  vcf.input <- data.table::fread(cmd=paste("grep -v '^#'", path),sep='\t')
  
  if (nrow(vcf.input) == 0) {return (vcf_dt) }
  
  if (ncol(vcf.input)==10) {
    setnames(vcf.input, paste0("V",seq(1:10)), c("seqnames","start","ID","REF",
                                                 "ALT","QUAL","FILTER","INFO",
                                                 "FORMAT","OCILY12"))
  }
  
  if ("INFO" %in% colnames(vcf.input)) {
    # Extract info from columns...
    
    ### SV TYPE
    vcf.input[, TYPE := gsub(".*?SVTYPE=([A-Z]+).*", "\\1", INFO)]
    
    ### SPAN
    vcf.input[, SPAN := as.numeric(gsub(".*?SVLEN=([-0-9]+).*","\\1",INFO))]
    vcf.input$SPAN <-  gsub("-","", vcf.input$SPAN)
    vcf.input[TYPE == "BND", SPAN := -1] 
    
    ### remove those that are not over 1000 bp
    ### manta identifies these differently
    
    vcf.input_s <- vcf.input[as.numeric(SPAN) >= 1000 | as.numeric(SPAN) == -1]
    
    ### remove imprecise calls, these do not have homology 
    
    vcf.input_ss <- vcf.input_s[!(grepl("IMPRECISE", INFO))]
    
    ### HOMSEQ
    vcf.input_ss[grepl("HOMSEQ", INFO), HOMSEQ := gsub(".*?HOMSEQ=([A-Z]+).*","\\1",INFO)]
    
    ### HOMLEN
    vcf.input_ss[, HOMLEN := as.numeric(gsub(".*?HOMLEN=([-0-9]+).*","\\1",INFO))]
    vcf.input_ss[is.na(HOMLEN), HOMSEQ := ""]
    
    ### remove non standard seqs (chromosomes), also for alts
    vcf.input_ss$seqnames <- gsub("chr", "",vcf.input_ss$seqnames)
    vcf.input_com_chrom <- vcf.input_ss[seqnames %in% c(1:22,"X")] 
    vcf.input_com_chrom_s <- vcf.input_com_chrom[!(grepl("chrUn",ALT))]
    vcf.input_com_chrom_ss <- vcf.input_com_chrom_s[!(grepl("random",ALT))]
    vcf.input_com_chrom_sss <- vcf.input_com_chrom_ss[!(grepl("alt",ALT))]
    
    #### build bedpe from non BND calls
    cat("Building bedpe...")
    non_bnd <- vcf.input_com_chrom_sss[TYPE %in% c("DEL","DUP","INV")]
    
    non_bnd[,CHROM_A := seqnames ]
    non_bnd[,CHROM_B := seqnames ]
    
    non_bnd[, START_A := start]
    non_bnd[, END_A := start + HOMLEN]
    non_bnd[is.na(END_A), END_A := START_A]
    
    non_bnd[, START_B := as.numeric(gsub(".*?END=([-0-9]+).*","\\1",INFO))]
    non_bnd[, END_B := START_B + HOMLEN]
    non_bnd[is.na(END_B), END_B := START_B]
    
    non_bnd[, STRAND_A := "+"]
    non_bnd[, STRAND_B := "-"]
    
    non_bnd[, NAME_A := ID]
    non_bnd[, NAME_B := "."]
    
    non_bnd[, REF_A := REF]
    non_bnd[, ALT_A := ALT]
    
    non_bnd[, REF_B := "."]
    non_bnd[, ALT_B := "."]
    
    non_bnd[, INFO_A := INFO]
    non_bnd[, INFO_B := "."]
    
    non_bnd[, uuid := paste0(unlist(strsplit(path,".somaticSV"))[1], ":",CHROM_A,
                             ":", START_A,":", END_A,":", CHROM_B, ":", START_B,
                             ":", END_B)]
    
    non_bnd_save <- non_bnd[,c("CHROM_A","START_A","END_A","CHROM_B","START_B",
                               "END_B", "ID", "QUAL", "STRAND_A","STRAND_B","TYPE",
                               "FILTER","NAME_A","REF_A", "ALT_A","NAME_B","REF_B",
                               "ALT_B","INFO_A","INFO_B", "FORMAT", "OCILY12","SPAN",
                               "HOMSEQ", "HOMLEN","uuid")]
    
    
    #### build bedpe from BND, really could be improved
    bnd_ <- vcf.input_com_chrom_sss[TYPE == "BND"]
    bnd_[,MATE_ID := unlist(strsplit(unlist(strsplit(INFO, "MATEID="))[2],"[;]"))[1], by = "ID"]
    bnd_bed <- data.table()
    
    mates_used <- NULL
    for(i in 1:length(bnd_$MATE_ID)) {
      which_mate_A <- bnd_[i,]
      which_mate_B <- bnd_[which_mate_A$ID == MATE_ID]
      
      if(nrow(which_mate_A) == 1 && nrow(which_mate_B) == 1) {
        if(!(which_mate_A$ID %in% mates_used) & !(which_mate_B$ID %in% mates_used)) {
          
          bed_temp <- as.data.table(cbind(which_mate_A$seqnames, which_mate_A$start, which_mate_A$start + which_mate_A$HOMLEN, 
                            which_mate_B$seqnames, which_mate_B$start, which_mate_B$start + which_mate_B$HOMLEN))
          colnames(bed_temp) <- c("CHROM_A","START_A","END_A","CHROM_B","START_B",
                                  "END_B")
          bed_temp[is.na(END_A), END_A := START_A]
          bed_temp[is.na(END_B), END_B := START_B]
          
          #### determine orientation
          strA <- strsplit(which_mate_A$ALT, "")
          str_m <- grep("[[]",strA)
          str_p <- grep("[]]", strA)
          
          strandB <- ifelse((length(str_m) > 0), "-", ifelse((length(str_p) > 0), "+", "*"))
          
          strB <- strsplit(which_mate_B$ALT, "")
          str_m <- grep("[[]",strB)
          str_p <- grep("[]]", strB)
        
          strandA <- ifelse((length(str_m) > 0), "-", ifelse((length(str_p) > 0), "+", "*"))
          
          bed_temp <- cbind(bed_temp, which_mate_A$ID, which_mate_A$QUAL, strandA,strandB, which_mate_A$TYPE,
                            which_mate_A$FILTER, which_mate_A$ID, which_mate_A$REF, 
                            which_mate_A$ALT, which_mate_B$ID, which_mate_B$REF, which_mate_B$ALT,
                            which_mate_A$INFO, which_mate_B$INFO, which_mate_A$FORMAT, which_mate_A$OCILY12, which_mate_A$SPAN) #### neeed to use the alt seq to adjust start positons of these (nchar(gsub[A-Z]))
          
          colnames(bed_temp)[7:23] <- c("ID", "QUAL", "STRAND_A","STRAND_B","TYPE",
                                        "FILTER","NAME_A","REF_A", "ALT_A","NAME_B","REF_B",
                                        "ALT_B","INFO_A","INFO_B", "FORMAT", "OCILY12","SPAN")
          
          bed_temp <- cbind(bed_temp, which_mate_A$HOMSEQ, which_mate_A$HOMLEN)
          colnames(bed_temp)[24:25] <- c("HOMSEQ", "HOMLEN")
          bed_temp[, uuid := paste0(unlist(strsplit(path,".somaticSV"))[1], ":",CHROM_A,
                                    ":", START_A,":", END_A,":", CHROM_B, ":", START_B,
                                    ":", END_B)]
          
          mates_used <- c(mates_used, which_mate_A$ID, which_mate_B$ID)
          
          bnd_bed <- rbind(bnd_bed, bed_temp)
        }
      } 
    }
    
    bedpe <- rbind(non_bnd_save, bnd_bed)
    bedpe$START_A <- as.numeric(as.character(bedpe$START_A))
    bedpe$START_B <- as.numeric(as.character(bedpe$START_B))
    bedpe$END_A <- as.numeric(as.character(bedpe$END_A))
    bedpe$END_B <- as.numeric(as.character(bedpe$END_B))
  }
  return(bedpe)
}

##########################
# Functions to annotate SVs
##########################

#### Determines if each breakpoint occurs in an encode defined gene annotated region
annotate_sv = function(bedpe, gene_locations) {
  cat("Adding gene annotations...\n")
  
  for(i in 1:nrow(bedpe)) {
    if((i %% 50) == 0) {cat(i,"/", nrow(bedpe), " SVs annotated for transcription disruption\n")} ### gives idea of progress
    
    ### Annotate break point A by giving genes that exist where the breakpoint occurs
    temp_A_gr <- GRanges(paste0("chr",bedpe$CHROM_A[i]), 
                         IRanges(bedpe$START_A[i], bedpe$END_A[i]))
    temp_A_gr_ant <-  suppressWarnings(gene_locations %&% temp_A_gr)
    temp_A_dt_ant <- as.data.table(temp_A_gr_ant)
    
    ### decision table for 
    if(!nrow(temp_A_dt_ant) == 0) {
      if(nrow(temp_A_dt_ant) == 1) {
        
        bedpe$BREAK_A_GENE[i] <- temp_A_dt_ant$GENE_NAME[1]
        bedpe$BREAK_A_GENE_TYPE[i] <- temp_A_dt_ant$GENE_TYPE[1]
        
      } else {
        
        bedpe$BREAK_A_GENE[i] <- paste0(temp_A_dt_ant$GENE_NAME, collapse = ";")
        bedpe$BREAK_A_GENE_TYPE[i] <- paste0(temp_A_dt_ant$GENE_TYPE, collapse = ";")
        
      }
    } else {
      
      bedpe$BREAK_A_GENE[i] <- ""
      bedpe$BREAK_A_GENE_TYPE[i] <- ""
      
    }
    
    ### Annotate break point B
    temp_B_gr <- GRanges(paste0("chr",bedpe$CHROM_B[i]), 
                         IRanges(bedpe$START_B[i], bedpe$END_B[i]))
    temp_B_gr_ant <-  suppressWarnings(gene_locations %&% temp_B_gr)
    temp_B_dt_ant <- as.data.table(temp_B_gr_ant)
    
    if(!nrow(temp_B_dt_ant) == 0) {
      if(nrow(temp_B_dt_ant) == 1) {
        
        bedpe$BREAK_B_GENE[i] <- temp_B_dt_ant$GENE_NAME[1]
        bedpe$BREAK_B_GENE_TYPE[i] <- temp_B_dt_ant$GENE_TYPE[1]
        
      } else {
        
        bedpe$BREAK_B_GENE[i] <- paste0(temp_B_dt_ant$GENE_NAME, collapse = ";")
        bedpe$BREAK_B_GENE_TYPE[i] <- paste0(temp_B_dt_ant$GENE_TYPE, collapse = ";")
        
      }
    } else {
      
      bedpe$BREAK_B_GENE[i] <- ""
      bedpe$BREAK_B_GENE_TYPE[i] <- ""
      
    }
  }
  
  cat("...done.\n")
  return(bedpe)
}

### gives list of protein coding genes in the tad
### makes 2 columns of those that are known oncogenic and those that are not
annotate_tad = function(bedpe, tad_gr, gene_locations) {
  cat("Adding TAD annotations...\n")
  
  ### get oncogene list
  if(LOCAL) {
    
    cosmic.genes <- fread("/Volumes/xchip_beroukhimlab/Alex/Manta_SV_Annotation/inputfiles/Census_cancer_genes_allTue_May_2018.tsv")$'Gene Symbol'
    
  } else {
    
    cosmic.genes <- fread("/xchip/beroukhimlab/Alex/Manta_SV_Annotation/inputfiles/Census_cancer_genes_allTue_May_2018.tsv")$'Gene Symbol' 
  
  }
  for(i in 1:nrow(bedpe)) {
    
    if((i %% 50) == 0) {cat(i,"/", nrow(bedpe), "SVs annotated in TADs\n")} ### gives idea of annotation progress
    rearrangement <- bedpe[i]
    
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
        bedpe$TAD_CORD[i] <- paste0(tad_gr_in_dt$TAD_NAME, collapse = ";")
        bedpe$ONCO_TAD_GENE[i] <- paste0(onco.tmp, collapse = ";")
        bedpe$TAD_GENE[i] <- paste0(non_onco.tmp, collapse = ";")

        
        
      } else {
        
        bedpe$TAD_CORD[i] <- ""
        bedpe$ONCO_TAD_GENE[i] <- ""
        bedpe$TAD_GENE[i] <- ""
        
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
        
        bedpe$TAD_CORD[i] <- ""
        bedpe$ONCO_TAD_GENE[i] <- ""
        bedpe$TAD_GENE[i] <- ""
        
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
        bedpe$TAD_CORD[i] <- t_name
        
        if((length(onco.tmp_A) + length(onco.tmp_B)) > 0) {
          
          bedpe$ONCO_TAD_GENE[i] <- paste0("A","{", paste0(onco.tmp_A, collapse = ";"), "}", "|","B","{", paste0(onco.tmp_B, collapse = ";"),"}")
          
        } else {
          
          bedpe$ONCO_TAD_GENE[i] <- ""
          
        }
        
        if((length(non_onco.tmp_A) + length(non_onco.tmp_B)) > 0) {
          
          bedpe$TAD_GENE[i] <- paste0("A","{", paste0(non_onco.tmp_A, collapse = ";"), "}", "|","B","{", paste0(non_onco.tmp_B, collapse = ";"),"}")
          
        } else {
          
          bedpe$TAD_GENE[i] <- ""
          
        }
      }
    }
  }
  
  cat("...done.\n")
  return(bedpe)
}



#### MAIN FUNCTION #####

option_list <- list(
  
  make_option(c("-i", "--input"),  type = "character", default = NULL,  help = "Input bedpe directory path"),
  
  make_option(c("-f", "--format"), type = "character", default = "bedpe", help = "Format: bedpe or vcf"),
  
  make_option(c("-c", "--cores"),  type = "numeric", default = 1,  help = "Number of cores"),
  
  make_option(c("-r", "--refgenome"),  type = "character", default = "hg38",  help = "Reference Genome: hg38")
  
)

# inputs
parseobj = OptionParser(option_list = option_list)
opt = parse_args(parseobj)
  
input_pth = paste0(opt$input)
output_file = sprintf("%s.sv_annotation", input_pth)
num_cores = as.numeric(opt$cores)

# get reference gene annotations and tad locations as granges
hg38_gene_gr_reference = build_hg38_geneloci_gr()
hg38_tad_gr_reference = build_hg38_tad_gr()

if (!file.exists(input_pth)) {
    
  stop(sprintf("Input file '%s' does not exist!", input_pth))
    
} else {
  
  cat("Reading file...\n")
  if (as.character(opt$format) == "vcf") {
    
    bedpe.input <- vcf_to_bedpe(input_pth)
    cat("done.\n")
    
    bedpe.annotated = annotate_sv(bedpe.input, hg38_gene_gr_reference)
    bedpe.annotated_tad = annotate_tad (bedpe.annotated, hg38_tad_gr_reference, hg38_gene_gr_reference)
    
    write.table(bedpe.annotated_tad, paste0(output_file), row.names = F, col.names = T, sep = "\t")
  } else {
    
    bedpe.input = suppressWarnings(fread(paste0(input_pth)))
    cat("done.\n")
    
    bedpe.annotated = annotate_sv(bedpe.input, hg38_gene_gr_reference)
    bedpe.annotated_tad = annotate_tad (bedpe.annotated, hg38_tad_gr_reference, hg38_gene_gr_reference)
  
    write.table(bedpe.annotated_tad, paste0(output_file), row.names = F, col.names = T, sep = "\t")
  }
}
  



