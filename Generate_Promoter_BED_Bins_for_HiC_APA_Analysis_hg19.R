# ==============================================================================
# Script: Generate Promoter BED Bins for Hi-C APA Analysis (hg19)
# Purpose: Retrieve TSS for all isoforms and bin them at 5Kb resolution
# ==============================================================================

# 1. Load Libraries
library(biomaRt)

# 2. Settings & Paths
res <- 5000  # 5Kb Resolution
input_path  <- "/Users/80030577/Desktop/HiC_analysis/HiC_icenorm_AggregatePeakAnalysis/Target_Gene_List_for_LCL.txt"
output_path <- "/Users/80030577/Desktop/HiC_analysis/HiC_icenorm_AggregatePeakAnalysis/Promoter_Bins_5Kb_hg19.bed"

# 3. Connect to Ensembl hg19 (GRCh37) Archive
# Using the specific host is critical for coordinate consistency in hg19 projects
ensembl <- useMart("ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   host = "https://grch37.ensembl.org")

# 4. Load Gene List
my_genes <- readLines(input_path)

# 5. Retrieve Transcription Start Sites (TSS) for all Isoforms
# 'transcription_start_site' returns the 5' end for every transcript
my_attributes <- c("external_gene_name", 
                   "ensembl_transcript_id", 
                   "chromosome_name", 
                   "transcription_start_site", 
                   "strand")

gene_coords <- getBM(attributes = my_attributes,
                     filters = "external_gene_name",
                     values = my_genes,
                     mart = ensembl,
                     useCache = FALSE)

# 6. Data Cleaning & Manual Overrides
# Keep only primary chromosomes (1-22, X, Y) to match Hi-C matrices
gene_coords <- gene_coords[gene_coords$chromosome_name %in% c(1:22, "X", "Y"), ]

# Manually add AK6 if it was not returned by the biomaRt query
if (!"AK6" %in% gene_coords$external_gene_name) {
  ak6_row <- data.frame(
    external_gene_name = "AK6",
    ensembl_transcript_id = "MANUAL",
    chromosome_name = "5",
    transcription_start_site = 68665840, # Official hg19 TSS for AK6
    strand = -1
  )
  gene_coords <- rbind(gene_coords, ak6_row)
}

# 7. Formatting for Hi-C (Binning Logic)
# 
# Snaps exact TSS coordinates to the start of the 5Kb window
gene_coords$bin_start <- floor(gene_coords$transcription_start_site / res) * res
gene_coords$bin_end   <- gene_coords$bin_start + res

# 8. Create Final BED Output
# Use unique() to ensure multiple isoforms in the same bin are only counted once
duplicated(gene_coords[, c("chromosome_name", "bin_start", "bin_end", "external_gene_name")])
bed_output <- unique(gene_coords[, c("chromosome_name", "bin_start", "bin_end", "external_gene_name")])

# Format Chromosomes (add 'chr') and Sort (essential for Hi-C tools)
bed_output$chromosome_name <- paste0("chr", bed_output$chromosome_name)
bed_output <- bed_output[order(bed_output$chromosome_name, bed_output$bin_start), ]

# 9. Save and Report
write.table(bed_output, output_path, sep="\t", quote=F, row.names=F, col.names=F)

cat("--- Processing Complete ---\n")
cat("Total Genes in Input:", length(unique(my_genes)), "\n")
cat("Total Unique Bins (including isoforms):", nrow(bed_output), "\n")
cat("BED file saved to:", output_path, "\n")

# Preview results
head(bed_output, 10)

row_number(gene_coords)
