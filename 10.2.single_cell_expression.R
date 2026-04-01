# Load required packages
library("ggplot2")
library("dplyr")
library("tidyr")
library("tibble")
library("pheatmap")
library("viridis")
library("RColorBrewer")
library("patchwork")  # For combining plots

#********************************************************
#* 1. download data from https://cellxgene.cziscience.com/
#* 2. delete the first 9 rows and make sure the variable names in the first row
#* 3. run this code, please check the line 47 and 51
#********************************************************
# Set working directory
setwd("D:/OneDrive/工作/1.工作/blood pressure")
print(paste("Current working directory:", getwd()))

# Read data
data_path <- "D:/OneDrive/工作/1.工作/blood pressure"
base_dir <- "D:/OneDrive/工作/1.工作/blood pressure/CELLxGENE_gene_expression"
data <- read.csv(data_path, stringsAsFactors = FALSE)
print("Data preview:")
print(head(data))
print(paste("Data dimensions:", nrow(data), "rows,", ncol(data), "columns"))

# Set Cell Count filter threshold (adjust as needed)
min_cell_count <- 100  # Keep rows with cell count >= 100

print(paste("Rows before filtering:", nrow(data)))
print(paste("Cell Count range:", min(data$Cell.Count, na.rm = TRUE), "-", max(data$Cell.Count, na.rm = TRUE)))

# Filter by Cell Count
data <- data %>% filter(Cell.Count >= min_cell_count)

# delete aggregated value
data <- data[-which(data$Cell.Type == "aggregated"),]

print(paste("Rows after filtering (Cell.Count >=", min_cell_count, "):", nrow(data)))

# Show all available tissues
print("Available tissues:")
print(unique(data$Tissue))

# Define target tissues (match names from data; may need adjustment)
target_tissues <- c("heart", "liver", "spleen", "kidney", "brain","pancreas","adipose tissue")

# Try fuzzy matching tissue names
tissues_in_data <- unique(data$Tissue)
matched_tissues <- tissues_in_data[grepl("heart|liver|spleen|kidney|brain|pancreas|adipose tissue", 
                                         tissues_in_data, ignore.case = TRUE)]


print("Matched tissues:")
print(matched_tissues)

# Filter tissue data
if(length(matched_tissues) > 0) {
  filtered_data <- data %>% filter(Tissue %in% matched_tissues)
} else {
  # If no match, use user-specified tissue names
  filtered_data <- data %>% filter(Tissue %in% target_tissues)
}

print(paste("Rows after tissue filter:", nrow(filtered_data)))

# Show gene list
print("Genes in data:")
print(unique(filtered_data$Gene.Symbol))

# Method 1: gene x tissue heatmap (mean expression)
print("Generating gene-tissue heatmap...")
tissue_expression <- filtered_data %>%
  group_by(Gene.Symbol, Tissue) %>%
  summarise(Mean_Expression = mean(Expression, na.rm = TRUE),
            Mean_Scaled = mean(Expression..Scaled, na.rm = TRUE),
            Cells_Expressing = sum(Number.of.Cells.Expressing.Genes, na.rm = TRUE),
            Total_Cells = sum(Cell.Count, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Percent_Expressing = Cells_Expressing / Total_Cells * 100)

# Draw dot heatmap
tryCatch({
  pdf(paste0(base_dir, "/single_cell_heatmap_gene_by_tissue.pdf"), width = 12, height = 10)
}, error = function(e) {
  print(paste("Error: unable to create PDF file. Please ensure:"))
  print("1. All open PDF files are closed")
  print("2. You have write permission in the working directory")
  print(paste("Current working directory:", getwd()))
  stop(e)
})

p1 <- ggplot(tissue_expression, aes(x = Tissue, y = Gene.Symbol)) +
  geom_point(aes(size = Percent_Expressing, color = Mean_Scaled)) +
  scale_color_gradientn(colors = c("#F7F8B3", "#DC2F24", "#440154"), 
                        name = "Scaled Gene\nExpression") +
  scale_size_continuous(range = c(1, 10), name = "% Cells\nExpressing") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()) +
  labs(title = "Gene expression across tissues",
       x = "Tissue", y = "Gene")
print(p1)
dev.off()
print(paste("Saved:", paste0(base_dir, "/single_cell_heatmap_gene_by_tissue.pdf")))

# Method 2: gene x cell type heatmap (per tissue)
print("Generating gene-cell type heatmaps...")
cell_expression <- filtered_data %>%
  group_by(Gene.Symbol, Tissue, Cell.Type) %>%
  summarise(Mean_Expression = mean(Expression, na.rm = TRUE),
            Mean_Scaled = mean(Expression..Scaled, na.rm = TRUE),
            Cell_Count = sum(Cell.Count, na.rm = TRUE),
            Cells_Expressing = sum(Number.of.Cells.Expressing.Genes, na.rm = TRUE),
            .groups = "drop") %>%
  filter(Cell_Count > 0) %>%  # Keep combinations with cells
  mutate(Percent_Expressing = Cells_Expressing / Cell_Count * 100)

# Store plots for later combination
tissue_plots <- list()

# Create a separate dot heatmap for each tissue
for(tissue in unique(cell_expression$Tissue)) {
  tissue_data <- cell_expression %>% 
    filter(Tissue == tissue) %>%
    filter(tolower(Cell.Type) != "cell")  # Exclude Cell.Type == "cell"
  
  if(nrow(tissue_data) > 0) {
    # Adjust figure size to cell type count
    n_celltypes <- n_distinct(tissue_data$Cell.Type)
    n_genes <- n_distinct(tissue_data$Gene.Symbol)
    width <- max(10, n_celltypes * 0.3 + 3)
    height <- max(8, n_genes * 0.3 + 2)
    
    filename <- paste0("single_cell_heatmap_", gsub(" ", "_", tissue), "_by_celltype.pdf")
    
    # Create the plot object first
    p <- ggplot(tissue_data, aes(x = Cell.Type, y = Gene.Symbol)) +
      geom_point(aes(size = Percent_Expressing, color = Mean_Scaled)) +
      scale_color_gradientn(colors = c("#F7F8B3", "#DC2F24", "#440154"), 
                            name = "Scaled Gene\nExpression") +
      scale_size_continuous(range = c(1, 10), name = "% Cells\nExpressing") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(size = 8),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank()) +
      labs(title = tissue,
       x = "Cell type", y = "Gene")
    
    # Store plot for combination
    tissue_plots[[tissue]] <- p
    
    # Save individual plot
    tryCatch({
      pdf(paste0(base_dir, "/", filename), width = width, height = height)
      print(p)
      dev.off()
      print(paste("Saved:", filename))
    }, error = function(e) {
      print(paste("Warning: unable to create", filename, "- skipping"))
    })
  }
}

# Method 3: combined heatmap (gene x tissue-cell type) - SKIPPED (not needed)
# Summary by cell type
summary_by_celltype <- filtered_data %>%
  group_by(Gene.Symbol, Tissue, Cell.Type) %>%
  summarise(
    Cell_Count = sum(Cell.Count, na.rm = TRUE),
    Expression = mean(Expression, na.rm = TRUE),
    Scaled_Expression = mean(Expression..Scaled, na.rm = TRUE),
    Cells_Expressing = sum(Number.of.Cells.Expressing.Genes, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Cell_Count > 0)

write.csv(summary_by_celltype, paste0(base_dir, "/single_cell_expression_summary_by_celltype.csv"), row.names = FALSE)
print("Saved: single_cell_expression_summary_by_celltype.csv")

# ========== Filter low-expression plots ==========
print("Generating filtered plots (remove low expression across all genes)...")

# Define low-expression threshold
low_expression_threshold <- 10  # Percent expressing < 10% is considered low

# Method 1 (filtered): drop tissues with low expression across all genes
print("Generating filtered gene-tissue heatmap...")
# Identify tissues with max expression below threshold across genes
tissues_to_keep <- tissue_expression %>%
  group_by(Tissue) %>%
  summarise(max_percent = max(Percent_Expressing, na.rm = TRUE)) %>%
  filter(max_percent >= low_expression_threshold) %>%
  pull(Tissue)

tissue_expression_filtered <- tissue_expression %>%
  filter(Tissue %in% tissues_to_keep)

if(nrow(tissue_expression_filtered) > 0) {
  tryCatch({
    pdf(paste0(base_dir, "/single_cell_heatmap_gene_by_tissue_filtered.pdf"), width = 12, height = 10)
    
    p1_filtered <- ggplot(tissue_expression_filtered, aes(x = Tissue, y = Gene.Symbol)) +
      geom_point(aes(size = Percent_Expressing, color = Mean_Scaled)) +
      scale_color_gradientn(colors = c("#F7F8B3", "#DC2F24", "#440154"), 
                            name = "Scaled Gene\nExpression") +
      scale_size_continuous(range = c(1, 10), name = "% Cells\nExpressing") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank()) +
      labs(title = "Gene expression across tissues (filtered)",
           x = "Tissue", y = "Gene")
    print(p1_filtered)
    dev.off()
    print(paste("Saved:", paste0(base_dir, "/single_cell_heatmap_gene_by_tissue_filtered.pdf")))
  }, error = function(e) {
    print(paste("Warning: unable to create filtered tissue heatmap"))
  })
}

# Method 2 (filtered): per-tissue cell type heatmaps
print("Generating filtered gene-cell type heatmaps...")
# Store filtered plots for combination
tissue_plots_filtered <- list()

for(tissue in unique(cell_expression$Tissue)) {
  tissue_data <- cell_expression %>% 
    filter(Tissue == tissue) %>%
    filter(tolower(Cell.Type) != "cell")  # Exclude Cell.Type == "cell"
  
  if(nrow(tissue_data) > 0) {
    # Identify cell types with max expression below threshold across genes
    celltypes_to_keep <- tissue_data %>%
      group_by(Cell.Type) %>%
      summarise(max_percent = max(Percent_Expressing, na.rm = TRUE)) %>%
      filter(max_percent >= low_expression_threshold) %>%
      pull(Cell.Type)
    
    tissue_data_filtered <- tissue_data %>%
      filter(Cell.Type %in% celltypes_to_keep)
    
    if(nrow(tissue_data_filtered) > 0 && length(celltypes_to_keep) > 0) {
      n_celltypes <- n_distinct(tissue_data_filtered$Cell.Type)
      n_genes <- n_distinct(tissue_data_filtered$Gene.Symbol)
      width <- max(10, n_celltypes * 0.3 + 3)
      height <- max(8, n_genes * 0.3 + 2)
      
      filename <- paste0("single_cell_heatmap_", gsub(" ", "_", tissue), "_by_celltype_filtered.pdf")
      
      # Create the plot object first
      p <- ggplot(tissue_data_filtered, aes(x = Cell.Type, y = Gene.Symbol)) +
        geom_point(aes(size = Percent_Expressing, color = Mean_Scaled)) +
        scale_color_gradientn(colors = c("#F7F8B3", "#DC2F24", "#440154"), 
                              name = "Scaled Gene\nExpression") +
        scale_size_continuous(range = c(1, 10), name = "% Cells\nExpressing") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              axis.text.y = element_text(size = 8),
              panel.grid.major = element_line(color = "grey90"),
              panel.grid.minor = element_blank()) +
        labs(title = tissue,
         x = "Cell type", y = "Gene")
      
      # Store filtered plot for combination
      tissue_plots_filtered[[tissue]] <- p
      print(paste("Stored filtered plot for tissue:", tissue))
      
      # Save individual plot
      tryCatch({
        pdf(paste0(base_dir, "/", filename), width = width, height = height)
        print(p)
        dev.off()
        print(paste("Saved:", paste0(base_dir, "/", filename)))
      }, error = function(e) {
        print(paste("Warning: unable to create", paste0(base_dir, "/", filename), "- skipping"))
      })
    } else {
      print(paste("Skipping tissue", tissue, "- no cell types passed filter threshold"))
    }
  } else {
    print(paste("Skipping tissue", tissue, "- no data after initial filter"))
  }
}

# Method 3 (filtered): combined heatmap - SKIPPED (not needed)

# ========== Combine all tissue plots into one large figure ==========
print("Combining all filtered tissue plots into one figure...")
print(paste("Number of filtered tissue plots stored:", length(tissue_plots_filtered)))
print(paste("Filtered tissue plot names:", paste(names(tissue_plots_filtered), collapse = ", ")))

if(length(tissue_plots_filtered) > 0) {
  tryCatch({
    # Arrange all filtered tissue plots vertically (one column)
    n_plots <- length(tissue_plots_filtered)
    
    # Combine plots using patchwork - arrange in 2 columns
    # Don't collect guides so each subplot keeps its own legend on the right
    combined_plot <- wrap_plots(tissue_plots_filtered, ncol = 2) +
      plot_annotation(
        title = "Gene Expression Across Tissues and Cell Types (Filtered)",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
    
    # Save combined plot - adjust dimensions for 2-column layout
    combined_width <- 24  # Wider for 2 columns
    combined_height <- ceiling(n_plots / 2) * 8  # Height based on number of rows
    
    pdf(paste0(base_dir, "/single_cell_heatmap_all_tissues_combined_filtered.pdf"), 
        width = combined_width, height = combined_height)
    print(combined_plot)
    dev.off()
    
    print(paste("Saved combined filtered plot:", paste0(base_dir, "/single_cell_heatmap_all_tissues_combined_filtered.pdf")))
    print(paste("Combined plot dimensions:", combined_width, "x", combined_height))
  }, error = function(e) {
    print(paste("Warning: unable to create combined filtered plot:", e$message))
    print(paste("Error details:", conditionMessage(e)))
  })
} else {
  print("No filtered tissue plots available to combine")
}
