suppressPackageStartupMessages(require(data.table))

###################################
# Functions to convert vcf to bedpe 
###################################

# Function description ?????
bnd_matching <- function(id) {
  which_mate_A <- bnd_[grepl(id, ID)]
  which_mate_B <- bnd_[which_mate_A$ID == MATE_ID]
  
  # ???
  if(nrow(which_mate_A) == 1 && nrow(which_mate_B) == 1) {
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
      # ????
      strandA <- ifelse((length(str_m) > 0), "-", ifelse((length(str_p) > 0), "+", "*"))
      
      bed_temp <- cbind(bed_temp, which_mate_A$ID, which_mate_A$QUAL, strandA,strandB, which_mate_A$TYPE,
                        which_mate_A$FILTER, which_mate_A$ID, which_mate_A$REF, 
                        which_mate_A$ALT, which_mate_B$ID, which_mate_B$REF, which_mate_B$ALT,
                        #### need to use the alt seq to adjust start positons of these (nchar(gsub[A-Z]))
                        which_mate_A$INFO, which_mate_B$INFO, which_mate_A$FORMAT, which_mate_A$OCILY12, which_mate_A$SPAN)
      
      colnames(bed_temp)[7:23] <- c("ID", "QUAL", "STRAND_A","STRAND_B","TYPE",
                                    "FILTER","NAME_A","REF_A", "ALT_A","NAME_B","REF_B",
                                    "ALT_B","INFO_A","INFO_B", "FORMAT", "OCILY12","SPAN")
      
      bed_temp <- cbind(bed_temp, which_mate_A$HOMSEQ, which_mate_A$HOMLEN)
      colnames(bed_temp)[24:25] <- c("HOMSEQ", "HOMLEN")
      bed_temp[, uuid := paste0(file_name, ":",CHROM_A,
                                ":", START_A,":", END_A,":", CHROM_B, ":", START_B,
                                ":", END_B)]

     return(bed_temp)
  }
}


# ?????
vcf_to_bedpe <- function(path) {
  cat("Converting VCF...\n")
  vcf.input <- data.table::fread(cmd=paste("grep -v '^#'", path), sep='\t')
  
  if (nrow(vcf.input) == 0) {stop(sprintf('This file is empty!'))}
  
  if (ncol(vcf.input)==10) {
    setnames(vcf.input, paste0("V",seq(1:10)), c("seqnames","start","ID","REF",
                                                 "ALT","QUAL","FILTER","INFO",
                                                 "FORMAT","OCILY12"))
  }
  # Extract info from columns...
  if ("INFO" %in% colnames(vcf.input)) {
    ### SV TYPE
    vcf.input[, TYPE := gsub(".*?SVTYPE=([A-Z]+).*", "\\1", INFO)]
    
    ### SPAN
    vcf.input[, SPAN := as.numeric(gsub(".*?SVLEN=([-0-9]+).*","\\1",INFO))]
    vcf.input$SPAN <-  gsub("-","", vcf.input$SPAN)
    vcf.input[TYPE == "BND", SPAN := -1] 
    vcf.input$seqnames <- gsub("chr", "",vcf.input$seqnames)
    
    ### remove those that are not over 1000 bp
    vcf.input_s <- vcf.input[as.numeric(SPAN) >= 1000 | as.numeric(SPAN) == -1]
    
    ### remove imprecise calls, these do not have homology 
    vcf.input_ss <- vcf.input_s[!(grepl("IMPRECISE", INFO))]
    
    ### HOMSEQ
    vcf.input_ss[grepl("HOMSEQ", INFO), HOMSEQ := gsub(".*?HOMSEQ=([A-Z]+).*","\\1",INFO)]
    
    ### HOMLEN
    vcf.input_ss[, HOMLEN := as.numeric(gsub(".*?HOMLEN=([-0-9]+).*","\\1",INFO))]
    vcf.input_ss[is.na(HOMLEN), HOMSEQ := ""]
    
    ### remove non standard seqs (chromosomes), also for alts
    vcf.input_com_chrom <- vcf.input_ss[seqnames %in% c(1:22,"X")] 
    vcf.input_com_chrom_s <- vcf.input_com_chrom[!(grepl("chrUn",ALT))]
    vcf.input_com_chrom_ss <- vcf.input_com_chrom_s[!(grepl("random",ALT))]
    vcf.input_com_chrom_sss <- vcf.input_com_chrom_ss[!(grepl("alt",ALT))]
    
    ### make removed calls file and write out
    removed <- vcf.input[!(ID %in% vcf.input_com_chrom_sss$ID)]
    removed_name <- sprintf("%s.removed_calls", file_name)
    write.table(removed, paste0(output_dir,'/',removed_name), sep = '\t', col.names = T, row.names = F, quote = F)
    
    #### build bedpe from non BND calls
    cat("Building bedpe...\n")
    non_bnd <- vcf.input_com_chrom_sss[TYPE %in% c("DEL","DUP","INV")]
    
    non_bnd[,CHROM_A := seqnames ]
    non_bnd[,CHROM_B := seqnames ]
    # ??? what happens exactly?
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
    
    non_bnd[, uuid := paste0(file_name, ":",CHROM_A,
                             ":", START_A,":", END_A,":", CHROM_B, ":", START_B,
                             ":", END_B)]
    
    non_bnd_save <- non_bnd[,c("CHROM_A","START_A","END_A","CHROM_B","START_B",
                               "END_B", "ID", "QUAL", "STRAND_A","STRAND_B","TYPE",
                               "FILTER","NAME_A","REF_A", "ALT_A","NAME_B","REF_B",
                               "ALT_B","INFO_A","INFO_B", "FORMAT", "OCILY12","SPAN",
                               "HOMSEQ", "HOMLEN","uuid")]
    
    
    #### build bedpe from BND
    # ??? what happens exactly
    bnd_ = vcf.input_com_chrom_sss[TYPE == "BND"]
    bnd_[,MATE_ID := unlist(strsplit(unlist(strsplit(INFO, "MATEID="))[2],"[;]"))[1], by = "ID"]
    
    cat('Matching breakends...\n')
    bnd_bed <- lapply(bnd_$ID, bnd_matching)
    bnd_bed <- rbindlist(bnd_bed)
    
    which_dup <- NULL
    for(i in 1:nrow(bnd_bed)) {
      tmp_id <- bnd_bed$NAME_A[i]
      id_mate <- which(bnd_bed$NAME_B == tmp_id)
      if (i < id_mate) {
        which_dup <- rbind(which_dup, cbind(i, id_mate))
      } else {
        which_dup <- rbind(which_dup, cbind(id_mate, i))
      }
    }
    which_dup_dedup <- as.data.frame(which_dup[!(duplicated(which_dup)),])
    
    bnd_bed_dedup <- bnd_bed[which_dup_dedup$i,]
    
    bedpe <- rbind(non_bnd_save, bnd_bed_dedup)
    bedpe$START_A <- as.numeric(as.character(bedpe$START_A))
    bedpe$START_B <- as.numeric(as.character(bedpe$START_B))
    bedpe$END_A <- as.numeric(as.character(bedpe$END_A))
    bedpe$END_B <- as.numeric(as.character(bedpe$END_B))
  }
  return(bedpe)
}

option_list <- list(
  make_option(c("-i", "--input"),  type = "character", default = NULL,  help = "Input vcf path"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output directory")
)

parseobj = OptionParser(option_list = option_list)
opt = parse_args(parseobj)

input_pth = paste0(opt$input)
file_name = unlist(strsplit(input_pth, '/'))[[length(unlist(strsplit(input_pth, '/')))]]
output_dir = paste0(opt$output)
output_file = sprintf("%s.bedpe", file_name)

out_bed <- vcf_to_bedpe(input_pth)
write.table(out_bed, paste0(output_dir, '/',output_file), sep = '\t', row.names = F, col.names = T, quote = F)
