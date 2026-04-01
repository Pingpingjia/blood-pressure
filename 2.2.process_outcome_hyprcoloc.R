# Title: Process Outcome HyPrColoc Results
# Author: Analysis Script
# Date: 2025-11-11
# Purpose: Filter and combine hyprcoloc results from multiple folders
# ===========================================================
library(data.table)
library(dplyr)
library(stringr)

# Set paths
input_dir <- "E:/results/analysis_results/outcome_hyprcoloc"
gene_file <- "D:/OneDrive/工作/1.工作/comorbidity/results/hg38_gene_list.txt"
output_dir <- "D:/OneDrive/工作/1.工作/comorbidity/results/combined_results_11/outcome_hyprcoloc_processed"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# =========================================================
# 1. Read gene annotation file
# =========================================================
cat("Reading gene annotation file...\n")
gene_dt <- fread(gene_file, header = FALSE, stringsAsFactors = FALSE)
colnames(gene_dt) <- c("chr", "start_pos", "end_pos", "gene", "ensembl_id")

# Ensure chromosome format is consistent (remove 'chr' prefix if present)
gene_dt[, chr := gsub("^chr", "", chr)]
gene_dt[, chr := as.character(chr)]

cat("Gene annotation loaded:\n")
print(head(gene_dt))
cat("Total genes:", nrow(gene_dt), "\n\n")

# =========================================================
# 2. Set key for gene_dt for efficient joining
# =========================================================
setkey(gene_dt, chr, start_pos, end_pos)

# =========================================================
# 3. Process each folder
# =========================================================
cat("Scanning folders in:", input_dir, "\n")

# Get all subdirectories
folders <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)

# Filter out folders containing "glgc" in their name
folders <- folders[!grepl("glgc", basename(folders), ignore.case = TRUE)]

cat("Found", length(folders), "folders to process (excluding glgc folders)\n\n")

# List to store results from all folders
all_results <- list()

for (folder in folders) {
    folder_name <- basename(folder)
    cat("Processing folder:", folder_name, "\n")
    
    # Find all .txt files in the folder
    txt_files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
    
    if (length(txt_files) == 0) {
        cat("  No .txt files found, skipping...\n\n")
        next
    }
    
    # Process each file in the folder
    for (file_path in txt_files) {
        file_name <- basename(file_path)
        cat("  Reading file:", file_name, "\n")
        
        # Read the file
        tryCatch({
            dt <- fread(file_path, header = TRUE, stringsAsFactors = FALSE)
            
            # Check if file has the expected columns
            if (ncol(dt) < 14) {
                cat("    Warning: File has fewer than 14 columns, skipping...\n")
                next
            }
            
            # Keep only first 14 columns
            dt <- dt[, 1:14]
            
            # Expected column names (first 14)
            expected_cols <- c("iteration", "traits", "posterior_prob", "regional_prob", 
                             "candidate_snp", "posterior_explained_by_snp", "dropped_trait", 
                             "ref_trait", "chr", "start", "stop", "SNP_num", "CHR", "BP")
            
            # If column names match, use them; otherwise assign expected names
            if (ncol(dt) == 14) {
                colnames(dt) <- expected_cols
            }
            
            # Ensure CHR and BP are properly formatted
            dt[, CHR := gsub("^chr", "", as.character(CHR))]
            dt[, BP := as.numeric(BP)]
            
            # Map SNPs to gene regions (using the same logic as 6.1.hyprcoloc.R)
            cat("    Mapping SNPs to genes...\n")
            dt[, region := {
                if (!is.na(CHR) && !is.na(BP)) {
                    matched <- gene_dt[chr == CHR & start_pos <= BP & end_pos >= BP, gene]
                    if (length(matched) > 0) {
                        paste(unique(matched), collapse = ";")
                    } else {
                        "intergenic"
                    }
                } else {
                    "intergenic"
                }
            }, by = seq_len(nrow(dt))]
            
            # For each unique combination of traits and region, keep only the row with highest posterior_prob
            cat("    Filtering for maximum posterior probability...\n")
            dt_filtered <- dt[order(-posterior_prob), .SD[1], by = .(traits, region)]
            
            # Add folder name as source
            dt_filtered[, source_folder := folder_name]
            
            cat("    Kept", nrow(dt_filtered), "out of", nrow(dt), "rows\n")
            
            # Add to results list
            all_results[[paste(folder_name, file_name, sep = "_")]] <- dt_filtered
            
        }, error = function(e) {
            cat("    Error processing file:", conditionMessage(e), "\n")
        })
    }
    cat("\n")
}

# =========================================================
# 4. Combine all results
# =========================================================
if (length(all_results) > 0) {
    cat("Combining all results...\n")
    combined_dt <- rbindlist(all_results, fill = TRUE)
    
    cat("Total rows before final deduplication:", nrow(combined_dt), "\n")
    
    # Final deduplication: for each unique combination of traits and region,
    # keep only the row with highest posterior_prob across ALL files
    cat("Performing final deduplication by traits and region...\n")
    combined_dt <- combined_dt[order(-posterior_prob), .SD[1], by = .(traits, region)]
    
    cat("Total rows after final deduplication:", nrow(combined_dt), "\n")
    
    # Reorder columns
    col_order <- c("source_folder", "iteration", "traits", "posterior_prob", "regional_prob", 
                   "candidate_snp", "posterior_explained_by_snp", "dropped_trait", 
                   "ref_trait", "chr", "start", "stop", "SNP_num", "CHR", "BP", "region")
    
    combined_dt <- combined_dt[, ..col_order]
    
    # Sort by posterior_prob (descending)
    setorder(combined_dt, -posterior_prob)
    
    # =========================================================
    # 5. Save results
    # =========================================================
    cat("\nSaving results...\n")
    
    output_file <- file.path(output_dir, "combined_hyprcoloc_results.txt")
    fwrite(combined_dt, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat("\n=================================================\n")
    cat("Processing complete!\n")
    cat("=================================================\n")
    cat("Total rows processed:", nrow(combined_dt), "\n")
    cat("Unique trait combinations:", length(unique(combined_dt$traits)), "\n")
    cat("Unique regions:", length(unique(combined_dt$region)), "\n")
    cat("Results saved to:", output_file, "\n")
    cat("=================================================\n")
    
    # Preview results
    cat("\nPreview of results (first 20 rows):\n")
    print(head(combined_dt, 20))
    
    # =========================================================
    # 6. Create summary statistics
    # =========================================================
    cat("\nCreating summary statistics...\n")
    
    # Summary by traits
    traits_summary <- combined_dt[, .(
        n_regions = length(unique(region)),
        mean_posterior = mean(posterior_prob, na.rm = TRUE),
        max_posterior = max(posterior_prob, na.rm = TRUE),
        n_records = .N
    ), by = traits]
    
    setorder(traits_summary, -n_records)
    
    summary_file <- file.path(output_dir, "traits_summary.txt")
    fwrite(traits_summary, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Traits summary saved to:", summary_file, "\n")
    
    # Summary by region
    region_summary <- combined_dt[region != "intergenic", .(
        n_trait_combos = length(unique(traits)),
        mean_posterior = mean(posterior_prob, na.rm = TRUE),
        max_posterior = max(posterior_prob, na.rm = TRUE),
        n_records = .N
    ), by = region]
    
    setorder(region_summary, -n_records)
    
    region_summary_file <- file.path(output_dir, "region_summary.txt")
    fwrite(region_summary, region_summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Region summary saved to:", region_summary_file, "\n")
    
    # =========================================================
    # 7. Additional processing: Filter for traits with at least 3 elements
    # =========================================================
    cat("\n=================================================\n")
    cat("Processing results for traits with >= 3 elements...\n")
    cat("=================================================\n")
    
    # Define trait groups
    lipid_group <- c("ldlc", "hdlc", "tg")
    glucose_group <- c("glu", "diabetes", "fbg")
    blood_pressure_group <- c("sbp", "dbp", "hypertension")
    obesity_group <- c("bmi", "wc", "whr")
    
    # Count the number of traits in each combination (separated by commas)
    combined_dt[, n_traits := sapply(strsplit(traits, ","), length)]
    
    # Filter for traits with at least 3 elements
    combined_dt_3 <- combined_dt[n_traits >= 3]
    
    cat("Rows with >= 3 traits:", nrow(combined_dt_3), "out of", nrow(combined_dt), "\n")
    
    # Add group combination column
    cat("Adding group combination column...\n")
    combined_dt_3[, group_combination := {
        # Split traits by comma and remove whitespace
        trait_list <- strsplit(traits, ",")[[1]]
        trait_list <- trimws(trait_list)
        
        # Extract trait names (remove prefixes like bbj_, ieu_, ckb_)
        trait_names <- gsub("^(bbj_|ieu_|ckb_|glgc_)", "", tolower(trait_list))
        
        # Determine which groups are present
        has_lipid <- any(trait_names %in% lipid_group)
        has_glucose <- any(trait_names %in% glucose_group)
        has_bp <- any(trait_names %in% blood_pressure_group)
        has_obesity <- any(trait_names %in% obesity_group)
        
        # Build group combination string
        groups <- c()
        if (has_lipid) groups <- c(groups, "lipid")
        if (has_glucose) groups <- c(groups, "glucose")
        if (has_bp) groups <- c(groups, "blood_pressure")
        if (has_obesity) groups <- c(groups, "obesity")
        
        paste(groups, collapse = "_")
    }, by = seq_len(nrow(combined_dt_3))]
    
    # Remove the n_traits column before saving
    combined_dt_3[, n_traits := NULL]
    
    # Save filtered results
    output_file_3 <- file.path(output_dir, "combined_hyprcoloc_results_3.txt")
    fwrite(combined_dt_3, output_file_3, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Results with >= 3 traits saved to:", output_file_3, "\n")
    
    # Summary by traits (for >= 3 traits), with all regions listed
    traits_summary_3 <- combined_dt_3[, .(
        n_regions = length(unique(region)),
        mean_posterior = mean(posterior_prob, na.rm = TRUE),
        max_posterior = max(posterior_prob, na.rm = TRUE),
        n_records = .N,
        all_regions = paste(unique(region), collapse = ";"),
        group_combination = unique(group_combination)[1]
    ), by = traits]
    
    setorder(traits_summary_3, -n_records)
    
    summary_file_3 <- file.path(output_dir, "traits_summary_3.txt")
    fwrite(traits_summary_3, summary_file_3, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Traits summary (>= 3 traits) saved to:", summary_file_3, "\n")
    
    # Summary by region (for >= 3 traits)
    region_summary_3 <- combined_dt_3[region != "intergenic", .(
        n_trait_combos = length(unique(traits)),
        mean_posterior = mean(posterior_prob, na.rm = TRUE),
        max_posterior = max(posterior_prob, na.rm = TRUE),
        n_records = .N
    ), by = region]
    
    setorder(region_summary_3, -n_records)
    
    region_summary_file_3 <- file.path(output_dir, "region_summary_3.txt")
    fwrite(region_summary_3, region_summary_file_3, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Region summary (>= 3 traits) saved to:", region_summary_file_3, "\n")
    
    # =========================================================
    # 8. Update original traits_summary to include all_regions column
    # =========================================================
    cat("\nUpdating original traits_summary with all regions...\n")
    
    # Recreate traits_summary with all_regions column
    traits_summary <- combined_dt[, .(
        n_regions = length(unique(region)),
        mean_posterior = mean(posterior_prob, na.rm = TRUE),
        max_posterior = max(posterior_prob, na.rm = TRUE),
        n_records = .N,
        all_regions = paste(unique(region), collapse = ";")
    ), by = traits]
    
    setorder(traits_summary, -n_records)
    
    # Overwrite the original traits_summary file
    fwrite(traits_summary, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Updated traits_summary.txt with all_regions column\n")
    
    cat("\nAll files created:\n")
    cat("  1. combined_hyprcoloc_results.txt (main results)\n")
    cat("  2. traits_summary.txt (summary by trait combination with all regions)\n")
    cat("  3. region_summary.txt (summary by gene region)\n")
    cat("  4. combined_hyprcoloc_results_3.txt (results with >= 3 traits)\n")
    cat("  5. traits_summary_3.txt (summary for >= 3 traits with all regions)\n")
    cat("  6. region_summary_3.txt (region summary for >= 3 traits)\n")
    
    # =========================================================
    # 9. Extract PDF files containing _3_ or _4_ from plot folders
    # =========================================================
    cat("\n=================================================\n")
    cat("Extracting PDF plots with _3_ or _4_ patterns...\n")
    cat("=================================================\n")
    
    # Create output plot directory
    plot_output_dir <- file.path(output_dir, "plot")
    if (!dir.exists(plot_output_dir)) {
        dir.create(plot_output_dir, recursive = TRUE)
    }
    
    # Counter for copied files
    pdf_count <- 0
    
    # Scan all folders for plot subdirectories
    for (folder in folders) {
        folder_name <- basename(folder)
        plot_folder <- file.path(folder, "plot")
        
        if (dir.exists(plot_folder)) {
            # Find all PDF files in the plot folder
            pdf_files <- list.files(plot_folder, pattern = "\\.pdf$", full.names = TRUE)
            
            if (length(pdf_files) > 0) {
                # Filter PDF files containing _3_ or _4_
                filtered_pdfs <- pdf_files[grepl("_3_|_4_", basename(pdf_files))]
                
                if (length(filtered_pdfs) > 0) {
                    cat(sprintf("  Found %d matching PDFs in %s\n", length(filtered_pdfs), folder_name))
                    
                    # Copy each matching PDF to output directory
                    for (pdf_file in filtered_pdfs) {
                        pdf_name <- basename(pdf_file)
                        dest_path <- file.path(plot_output_dir, pdf_name)
                        
                        # Copy file
                        tryCatch({
                            file.copy(pdf_file, dest_path, overwrite = TRUE)
                            pdf_count <- pdf_count + 1
                        }, error = function(e) {
                            cat(sprintf("    Error copying %s: %s\n", pdf_name, conditionMessage(e)))
                        })
                    }
                }
            }
        }
    }
    
    cat("\n=================================================\n")
    cat("PDF extraction complete!\n")
    cat("Total PDFs copied:", pdf_count, "\n")
    cat("PDFs saved to:", plot_output_dir, "\n")
    cat("=================================================\n")
    
} else {
    cat("No results to combine. Please check the input directory and files.\n")
}

cat("\nScript completed!\n")
