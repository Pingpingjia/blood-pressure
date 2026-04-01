# Title: Extract BP-Proteins HyPrColoc Results
# Author: Analysis Script
# Date: 2026-01-23
# Purpose: Extract and combine hyprcoloc results containing blood pressure traits and proteins
# ===========================================================
library(data.table)
library(dplyr)
library(stringr)

# Set paths
input_dir <- "E:/results/analysis_results/ckb_ukb_bp_proteins"
output_dir <- "D:/OneDrive/工作/1.工作/blood pressure/ckb_ukb_bp_proteins"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# =========================================================
# 1. Get all target folders
# =========================================================
cat("Scanning folders in:", input_dir, "\n")

# Get all subdirectories matching the pattern
folders <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)

# Filter for folders matching the pattern: ckb_sbp_ckb_dbp_ckb_pp_ckb_map_ukb_ht_seq.*
folders <- folders[grepl("ckb_sbp_ckb_dbp_ckb_pp_ckb_map_ukb_ht_seq", basename(folders))]

cat("Found", length(folders), "folders matching the pattern\n\n")

if (length(folders) == 0) {
    stop("No matching folders found!")
}

# =========================================================
# 2. Define blood pressure traits to check
# =========================================================
bp_traits <- c("ckb_sbp", "ckb_dbp", "ckb_pp", "ckb_map", "ukb_ht")

# =========================================================
# 3. Process each folder
# =========================================================
all_results <- list()

for (folder in folders) {
    folder_name <- basename(folder)
    cat("Processing folder:", folder_name, "\n")
    
    # Look for all_significant_region_results.txt
    file_path <- file.path(folder, "all_significant_region_results.txt")
    
    if (!file.exists(file_path)) {
        cat("  Warning: all_significant_region_results.txt not found, skipping...\n\n")
        next
    }
    
    # Read the file
    tryCatch({
        dt <- fread(file_path, header = TRUE, stringsAsFactors = FALSE)
        
        cat("  Read", nrow(dt), "rows with", ncol(dt), "columns\n")
        
        if (nrow(dt) == 0) {
            cat("  Warning: File is empty, skipping...\n\n")
            next
        }
        
        # Check if 'traits' column exists
        if (!"traits" %in% colnames(dt)) {
            cat("  Warning: 'traits' column not found, skipping...\n\n")
            next
        }
        
        # Extract protein names from column names
        # Proteins have columns like: PROTEIN_ES, PROTEIN_SE, PROTEIN_P
        col_names <- colnames(dt)
        protein_cols <- grep("_ES$", col_names, value = TRUE)
        
        # Remove BP trait columns to get protein columns only
        bp_pattern <- paste0("^(", paste(bp_traits, collapse = "|"), ")_ES$")
        protein_cols <- protein_cols[!grepl(bp_pattern, protein_cols)]
        
        # Extract protein names (remove _ES suffix)
        protein_names <- gsub("_ES$", "", protein_cols)
        
        if (length(protein_names) == 0) {
            cat("  Warning: No protein columns found, skipping...\n\n")
            next
        }
        
        cat("  Found", length(protein_names), "proteins:", paste(head(protein_names, 5), collapse = ", "))
        if (length(protein_names) > 5) cat(", ...")
        cat("\n")
        
        # Filter rows where traits contain at least one protein name
        # Split traits by comma and check if any protein is present
        dt[, has_protein := {
            trait_list <- strsplit(traits, ",")[[1]]
            trait_list <- trimws(trait_list)
            # Check if any protein name is in the trait list
            any(protein_names %in% trait_list)
        }, by = seq_len(nrow(dt))]
        
        # Keep only rows with proteins
        dt_filtered <- dt[has_protein == TRUE]
        dt_filtered[, has_protein := NULL]
        
        cat("  Filtered to", nrow(dt_filtered), "rows containing proteins\n")
        
        if (nrow(dt_filtered) > 0) {
            # Add source folder
            dt_filtered[, source_folder := folder_name]
            
            # Add to results list
            all_results[[folder_name]] <- dt_filtered
        }
        
    }, error = function(e) {
        cat("  Error processing file:", conditionMessage(e), "\n")
    })
    
    cat("\n")
}

# =========================================================
# 4. Combine all results
# =========================================================
if (length(all_results) > 0) {
    cat("=================================================\n")
    cat("Combining all results...\n")
    cat("=================================================\n")
    
    combined_dt <- rbindlist(all_results, fill = TRUE)
    
    cat("Total rows combined:", nrow(combined_dt), "\n")
    
    # Move source_folder to first column
    if ("source_folder" %in% colnames(combined_dt)) {
        col_order <- c("source_folder", setdiff(colnames(combined_dt), "source_folder"))
        combined_dt <- combined_dt[, ..col_order]
    }
    
    # Sort by posterior_prob (descending) if the column exists
    if ("posterior_prob" %in% colnames(combined_dt)) {
        setorder(combined_dt, -posterior_prob)
    }
    
    # =========================================================
    # 5. Save results
    # =========================================================
    output_file <- file.path(output_dir, "all_significant_bp_proteins_hyprcoloc.txt")
    fwrite(combined_dt, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat("\n=================================================\n")
    cat("Processing complete!\n")
    cat("=================================================\n")
    cat("Total rows:", nrow(combined_dt), "\n")
    cat("Total columns:", ncol(combined_dt), "\n")
    cat("Unique trait combinations:", length(unique(combined_dt$traits)), "\n")
    cat("Results saved to:", output_file, "\n")
    cat("=================================================\n")
    
    # =========================================================
    # 6. Create summary statistics
    # =========================================================
    cat("\nCreating summary statistics...\n")
    
    # Extract protein from each trait combination
    combined_dt[, protein_in_traits := {
        trait_list <- strsplit(traits, ",")[[1]]
        trait_list <- trimws(trait_list)
        # Find traits that are not BP traits
        proteins <- setdiff(trait_list, bp_traits)
        paste(proteins, collapse = ";")
    }, by = seq_len(nrow(combined_dt))]
    
    # Summary by protein
    protein_summary <- combined_dt[, .(
        n_records = .N,
        n_trait_combos = length(unique(traits)),
        mean_posterior = mean(posterior_prob, na.rm = TRUE),
        max_posterior = max(posterior_prob, na.rm = TRUE)
    ), by = protein_in_traits]
    
    setorder(protein_summary, -n_records)
    
    protein_summary_file <- file.path(output_dir, "protein_summary.txt")
    fwrite(protein_summary, protein_summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Protein summary saved to:", protein_summary_file, "\n")
    
    # Preview results
    cat("\nPreview of results (first 10 rows):\n")
    print(head(combined_dt[, .(source_folder, traits, posterior_prob, CHR, BP, protein_in_traits)], 10))
    
    cat("\nTop 10 proteins by number of records:\n")
    print(head(protein_summary, 10))
    
} else {
    cat("No results to combine. Please check the input directory and files.\n")
}

cat("\nScript completed!\n")
