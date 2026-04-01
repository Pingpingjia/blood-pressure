#!/usr/bin/env Rscript
packages <- c("ggplot2","RColorBrewer","dplyr","readr","pheatmap","tibble","rio")
for(p in packages) if(!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cran.rstudio.com/")
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(readr)
library(pheatmap)
library(tibble)
library(rio)

base_dir <- "D:/OneDrive/工作/1.工作/blood pressure"
models <- c("model_1","model_2","model_3")

# GTEx tissue expression data path
gtex_path <- "D:/OneDrive/工作/1.工作/comorbidity/results/protein_organization_all.csv"

if(!file.exists(gtex_path)){
  stop("GTEx data not found at: ", gtex_path, "\nPlease update the path in the script.")
}

# Load GTEx tissue expression data
tissue_info <- import(gtex_path)

# Load GTEx tissue expression data
tissue_info <- import(gtex_path)

# Function to generate tissue expression heatmap
generate_tissue_heatmap <- function(genes, output_path, title){
  exp <- tissue_info %>% 
    filter(external_gene_name %in% genes) %>%
    rename(Gene = external_gene_name) %>%
    select(Gene, Adipose_Subcutaneous:Whole_Blood)
  
  if(nrow(exp) == 0){
    message("  No genes found in GTEx data")
    return(NULL)
  }
  
  # Convert to matrix
  exp_mat <- exp %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  # Define tissue groups
  tissue_group <- data.frame(
    Tissue = colnames(exp_mat)
  ) %>%
    mutate(Group = case_when(
      grepl("Adipose", Tissue) ~ "Adipose",
      grepl("Muscle", Tissue) ~ "Muscle",
      grepl("Brain", Tissue) ~ "Brain",
      grepl("Liver", Tissue) ~ "Liver",
      grepl("Blood", Tissue) ~ "Blood",
      grepl("Heart", Tissue) ~ "Heart",
      grepl("Kidney", Tissue) ~ "Kidney",
      TRUE ~ "Other"
    ))
  
  # Column annotation for pheatmap
  annotation_col <- data.frame(Group = tissue_group$Group)
  rownames(annotation_col) <- tissue_group$Tissue
  
  # Color palette for groups
  group_colors <- brewer.pal(length(unique(annotation_col$Group)), "Set2")
  names(group_colors) <- unique(annotation_col$Group)
  ann_colors <- list(Group = group_colors)
  
  # Color scheme matching the original heatmap
  my_colors <- colorRampPalette(c("#337CB2", "#FFFFFF", "#D64643"))(100) 
  # Generate heatmap
  p <- pheatmap(
    mat = exp_mat,
    color = my_colors,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "complete",
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 9,
    fontsize_col = 9,
    border_color = "white",
    legend = TRUE,
    cellheight = 8,
    main = title,
    filename = output_path,
    width = 12,
    height = max(6, 0.3 * nrow(exp_mat) + 2)
  )
  
  message("  Generated heatmap with ", nrow(exp_mat), " genes")
  return(p)
}

detect_protein_col <- function(df){
  cols <- names(df)
  hits <- grep("(?i)^(protein|gene|symbol|name)", cols, perl=TRUE, value=TRUE)
  if(length(hits)) return(hits[1])
  nonnum <- cols[!sapply(df, is.numeric)]
  if(length(nonnum)) return(nonnum[1])
  return(cols[1])
}

detect_expression_cols <- function(df, protein_col){
  cols <- setdiff(names(df), protein_col)
  exprs <- cols[tolower(cols) %in% c("expression","expr","value","value_mean","mean_expression","mean")]
  if(length(exprs)) return(exprs[1])
  numcols <- cols[sapply(df[cols], is.numeric)]
  if(length(numcols)) return(numcols[1])
  return(NULL)
}

for(model in models){
  message("Processing ", model)
  csv_path <- file.path(base_dir, model, "heatmap", "Level3_Regression_MR_Coloc","circular_heatmap_data.csv")
  out_dir <- file.path(base_dir, model, "expression")
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  
  if(!file.exists(csv_path)){ 
    warning("CSV not found: ", csv_path)
    next 
  }
  
  df <- read.csv(csv_path, stringsAsFactors=FALSE, check.names=FALSE)
  
  # Check if gene name column exists
  gene_col <- NULL
  if("Entrez.Gene.Name" %in% names(df)){
    gene_col <- "Entrez.Gene.Name"
  } else {
    gene_candidates <- grep("(?i)(gene.*name|entrez|symbol)", names(df), perl=TRUE, value=TRUE)
    if(length(gene_candidates) > 0){
      gene_col <- gene_candidates[1]
    }
  }
  
  if(is.null(gene_col)){
    warning("No gene name column found in ", csv_path)
    next
  }
  
  # Get unique genes
  all_genes <- unique(df[[gene_col]])
  all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
  
  message("  Found ", length(all_genes), " unique genes")
  
  # Generate tissue expression heatmap for all genes
  out_file <- file.path(out_dir, paste0("tissue_expression_", model, ".pdf"))
  generate_tissue_heatmap(
    genes = all_genes,
    output_path = out_file,
    title = paste0("Tissue expression - ", model)
  )
  
  # Process validated proteins
  # Process validated proteins
  val_path <- file.path(base_dir, model, "validation", "validated_proteins_list.txt")
  if(file.exists(val_path)){
    # Read validated list as a table
    val_df <- tryCatch({
      read.table(val_path, header=TRUE, sep="\t", stringsAsFactors=FALSE, fill=TRUE, comment.char="")
    }, error=function(e){
      tryCatch({
        read.table(val_path, header=TRUE, sep="", stringsAsFactors=FALSE, fill=TRUE, comment.char="")
      }, error=function(e2){
        data.frame(gene=readLines(val_path), stringsAsFactors=FALSE)
      })
    })
    
    # Find gene name column
    gene_col_candidates <- c("Entrez.Gene.Name", "Gene.Name", "gene", "Gene", "symbol", "Symbol")
    val_gene_col <- NULL
    for(cand in gene_col_candidates){
      if(cand %in% names(val_df)){
        val_gene_col <- cand
        break
      }
    }
    if(is.null(val_gene_col)){
      val_gene_col <- names(val_df)[1]
    }
    
    vals <- val_df[[val_gene_col]]
    vals <- trimws(vals)
    vals <- vals[vals!="" & !is.na(vals)]
    
    message("  Found ", length(vals), " validated proteins")
    
    # Generate tissue expression heatmap for validated genes
    out_file2 <- file.path(out_dir, paste0("tissue_expression_validated_", model, ".pdf"))
    generate_tissue_heatmap(
      genes = vals,
      output_path = out_file2,
      title = paste0("Tissue expression (validated) - ", model)
    )
  } else {
    message("  Validated list not found for ", model)
  }
}

message("Done")
