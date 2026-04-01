# Title: Process HyPrColoc Results (simplified)
# Purpose: Organize colocalization results into summary tables
# ===========================================================
library(data.table)
library(dplyr)
library(stringr)

input_file <- "D:/OneDrive/工作/1.工作/blood pressure/hyprcoloc/hyprcoloc_all_significant_results.txt"
gene_file <- "D:/OneDrive/工作/1.工作/comorbidity/results/hg38_gene_list.txt"
output_dir <- "D:/OneDrive/工作/1.工作/blood pressure/hyprcoloc"

# Read data
coloc_dt <- fread(input_file, header = TRUE, stringsAsFactors = FALSE)
gene_dt <- fread(gene_file, header = FALSE, stringsAsFactors = FALSE)
colnames(gene_dt) <- c("chr", "start_pos", "end_pos", "gene", "ensembl_id")

# Normalize chromosome columns
coloc_dt[, CHR := gsub("^chr", "", CHR)]
coloc_dt[, CHR := as.character(CHR)]
coloc_dt[, BP := as.numeric(BP)]

gene_dt[, chr := gsub("^chr", "", chr)]
gene_dt[, chr := as.character(chr)]
gene_dt[, start_pos := as.numeric(start_pos)]
gene_dt[, end_pos := as.numeric(end_pos)]

# Remove rows with missing coordinates to avoid foverlaps errors
gene_dt <- gene_dt[!is.na(chr) & !is.na(start_pos) & !is.na(end_pos)]

# Trait groups (only blood pressure vs other)
bp_traits <- c("ckb_sbp", "ckb_dbp", "ckb_pp", "ckb_map", "ukb_hypertension")

# Map SNPs to genes using foverlaps
snp_dt <- coloc_dt[, .(CHR, BP, row_id = .I)]
snp_dt[, `:=`(start_pos = BP, end_pos = BP)]
 # Drop SNP rows with missing positions
snp_dt <- snp_dt[!is.na(CHR) & !is.na(start_pos) & !is.na(end_pos)]
setkey(gene_dt, chr, start_pos, end_pos)
setkey(snp_dt, CHR, start_pos, end_pos)

hits_dt <- foverlaps(
    snp_dt,
    gene_dt,
    by.x = c("CHR", "start_pos", "end_pos"),
    by.y = c("chr", "start_pos", "end_pos"),
    type = "within",
    nomatch = NA
)

snp_gene_dt <- hits_dt[, .(
    snp_gene = if (all(is.na(gene))) NA_character_ else paste(unique(gene), collapse = ";")
), by = row_id]

coloc_dt[, snp_gene := snp_gene_dt$snp_gene[match(.I, snp_gene_dt$row_id)]]
coloc_dt[, region_genes := snp_gene]

# Build summary table (same columns as 6.1.hyprcoloc.R)
coloc_expanded <- coloc_dt[, .(
    protein = Entrez.Gene.Name,
    protein_seq = protein_seq_name,
    trait = trimws(outcome),
    ref_trait = ref_trait,
    region = region_genes,
    region_range = paste0("chr", chr, ":", start, "-", stop),
    P1 = posterior_prob,
    rsID = candidate_snp,
    chr = CHR,
    gene = snp_gene,
    P2 = posterior_explained_by_snp,
    iteration = iteration
)]

coloc_expanded[, group := ifelse(trait %in% bp_traits, "blood_pressure", "other")]
coloc_expanded[, `:=`(
    region = ifelse(is.na(region) | region == "", "intergenic", region),
    gene = ifelse(is.na(gene) | gene == "", "intergenic", gene)
)]

result_dt <- coloc_expanded[order(-P1), .SD[1], by = .(protein, trait, gene)]
result_dt <- result_dt[, .(
    protein = protein,
    protein_seq = protein_seq,
    trait = trait,
    group = group,
    region = region,
    region_range = region_range,
    ref_trait = ref_trait,
    P1 = round(P1, 4),
    rsID = rsID,
    chr = chr,
    gene = gene,
    P2 = round(P2, 4)
)]
setorder(result_dt, protein, trait)

# Main summary table
fwrite(result_dt, file.path(output_dir, "hyprcoloc_summary_table.txt"),
             sep = "\t", quote = FALSE, row.names = FALSE)

# Additional summaries (keep existing outputs)
protein_group_summary <- result_dt[, .(
    n_coloc = .N,
    traits = paste(unique(trait), collapse = ";"),
    regions = length(unique(region))
), by = .(protein, group)]

fwrite(protein_group_summary,
             file.path(output_dir, "hyprcoloc_protein_group_summary.txt"),
             sep = "\t", quote = FALSE, row.names = FALSE)

protein_summary <- result_dt[, .(
    n_coloc = .N,
    groups = paste(unique(group), collapse = ";"),
    traits = paste(unique(trait), collapse = ";"),
    n_regions = length(unique(region))
), by = protein] 

setorder(protein_summary, -n_coloc)
fwrite(protein_summary,
             file.path(output_dir, "hyprcoloc_protein_summary.txt"),
             sep = "\t", quote = FALSE, row.names = FALSE)

# New table: max posterior per protein-trait-region
max_prob_region <- coloc_expanded[, .SD[which.max(P1)], by = .(protein_seq, trait, region)]
max_prob_region[, max_posterior_prob := P1]
max_prob_region[, P1 := NULL]
max_prob_region <- max_prob_region[!is.na(protein_seq)]
fwrite(max_prob_region,
       file.path(output_dir, "hyprcoloc_max_posterior_prob_by_region.txt"),
       sep = "\t", quote = FALSE, row.names = FALSE)

# New table: max posterior per protein-trait
max_prob_1 <- coloc_expanded[, .SD[which.max(P1)], by = .(protein_seq, trait)]
max_prob_1[, max_posterior_prob := P1]
max_prob_1[, P1 := NULL]
max_prob_2 <- max_prob_1[!is.na(protein_seq)]

fwrite(max_prob_2,
    file.path(output_dir, "hyprcoloc_max_posterior_prob.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE)

