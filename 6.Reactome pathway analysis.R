library(ggplot2)
library(ggrepel)
library(scales)
library(dplyr)
library(tidyr)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(tidyverse)
library(stringr)
library(openxlsx)
library(ReactomePA)
library(gridExtra)
library(cowplot)
library(patchwork)
# Suppress R CMD check warnings
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("GeneRatio", "GeneCount", "TotalGenes", "GeneRatio_percent",
                           "Description", "Count", "trait", "p.adjust", "intersection"))
}

# 避免命名空间冲突
select <- dplyr::select
filter <- dplyr::filter
n_top <- 5

#**************************************     通路分析 - 仅针对血压    ************************
#     直接使用BP_three_levels_protein_lists.xlsx中的三个层级数据
#     Level 1: Regression
#     Level 2: Regression+MR
#     Level 3: Regression+MR+Colocalization
#*************************************************************************************

# 通路分析函数
perform_pathway_analysis <- function(entrez_ids, trait_name) {
  results <- list()
  
  if(length(entrez_ids) < 3) {
    cat("Skipping", trait_name, "- insufficient genes (<3)\n")
    return(results)
  }
  
  cat("Analyzing", trait_name, "with", length(entrez_ids), "genes\n")
  
  # GO BP
  tryCatch({
    go_bp <- enrichGO(gene = entrez_ids, OrgDb = "org.Hs.eg.db", ont = "BP", 
                      pvalueCutoff = 0.05, readable = TRUE)
    if(!is.null(go_bp) && nrow(go_bp) > 0) {
      go_bp_df <- as.data.frame(go_bp) %>% arrange(`p.adjust`)
      go_bp_df$trait <- trait_name
      results$go_bp <- go_bp_df
      results$go_bp_top <- go_bp_df %>% slice_head(n=n_top)
    }
  }, error = function(e) { cat("Error in GO BP for", trait_name, "\n") })
  
  # GO CC
  tryCatch({
    go_cc <- enrichGO(gene = entrez_ids, OrgDb = "org.Hs.eg.db", ont = "CC", 
                      pvalueCutoff = 0.05, readable = TRUE)
    if(!is.null(go_cc) && nrow(go_cc) > 0) {
      go_cc_df <- as.data.frame(go_cc) %>% arrange(`p.adjust`)
      go_cc_df$trait <- trait_name
      results$go_cc <- go_cc_df
      results$go_cc_top <- go_cc_df %>% slice_head(n=n_top)
    }
  }, error = function(e) { cat("Error in GO CC for", trait_name, "\n") })
  
  # GO MF
  tryCatch({
    go_mf <- enrichGO(gene = entrez_ids, OrgDb = "org.Hs.eg.db", ont = "MF", 
                      pvalueCutoff = 0.05, readable = TRUE)
    if(!is.null(go_mf) && nrow(go_mf) > 0) {
      go_mf_df <- as.data.frame(go_mf) %>% arrange(`p.adjust`)
      go_mf_df$trait <- trait_name
      results$go_mf <- go_mf_df
      results$go_mf_top <- go_mf_df %>% slice_head(n=n_top)
    }
  }, error = function(e) { cat("Error in GO MF for", trait_name, "\n") })
  
  # KEGG
  tryCatch({
    kegg <- enrichKEGG(gene = entrez_ids, organism = "hsa", pvalueCutoff = 0.05)
    if(!is.null(kegg) && nrow(kegg) > 0) {
      kegg_df <- as.data.frame(kegg) %>% arrange(`p.adjust`)
      kegg_df$trait <- trait_name
      results$kegg <- kegg_df
      results$kegg_top <- kegg_df %>% slice_head(n=n_top)
    }
  }, error = function(e) { cat("Error in KEGG for", trait_name, "\n") })
  
  # Reactome Pathway
  tryCatch({
    reactome <- enrichPathway(gene = entrez_ids, organism = "human", pvalueCutoff = 0.05, readable = TRUE)
    if(!is.null(reactome) && nrow(reactome) > 0) {
      reactome_df <- as.data.frame(reactome) %>% arrange(`p.adjust`)
      reactome_df$trait <- trait_name
      results$reactome <- reactome_df
      results$reactome_top <- reactome_df %>% slice_head(n=n_top)
    }
  }, error = function(e) { cat("Error in Reactome for", trait_name, "\n") })
  
  return(results)
}

# 绘制气泡图函数
create_bubble_plot <- function(plotdf, title_text) {
  if(is.null(plotdf) || nrow(plotdf) == 0) return(NULL)
  
  plotdf <- plotdf %>%
    separate(GeneRatio, into = c("GeneCount", "TotalGenes"), sep = "/", convert = TRUE, remove = FALSE) %>%
    mutate(GeneRatio_percent = GeneCount / TotalGenes * 100) %>%
    mutate(Description = stringr::str_wrap(Description, width = 60))
  
  count_range <- range(plotdf$Count, na.rm = TRUE)
  
  if(count_range[2] <= 3) {
    size_range <- c(6, 16)
  } else if(count_range[2] <= 5) {
    size_range <- c(5, 14)
  } else if(count_range[2] <= 10) {
    size_range <- c(4, 12)
  } else {
    size_range <- c(3, 10)
  }
  
  color_var <- if("pathway_type" %in% colnames(plotdf)) "pathway_type" else "trait"
  
  p <- ggplot(plotdf, aes(x = GeneRatio_percent, y = reorder(Description, GeneRatio_percent), 
                          size = Count, color = !!sym(color_var))) +
    geom_point(alpha = 0.8) +
    scale_color_brewer(name = "Group", type = "qual", palette = "Set2") +
    scale_size_continuous(name = "Gene Count", range = size_range) +
    labs(x = "Gene Ratio (%)", y = "Pathway", title = title_text) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 12, hjust = 1),
          axis.text.x = element_text(size = 11),
          axis.title.x = element_text(size = 13, face = "bold"),
          axis.title.y = element_text(size = 13, face = "bold"),
          plot.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12, face = "bold"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  
  return(p)
}

# 绘制热图函数（参考4. Reactome pathway analysis_v2.R）
create_heatmap_plot <- function(plotdf, title_text, pathway_type = "GO-BP") {
  if(is.null(plotdf) || nrow(plotdf) == 0) return(NULL)
  
  plotdf <- plotdf %>%
    separate(GeneRatio, into = c("GeneCount", "TotalGenes"), sep = "/", convert = TRUE, remove = FALSE) %>%
    mutate(GeneRatio_percent = GeneCount / TotalGenes * 100) %>%
    mutate(Description = stringr::str_wrap(Description, width = 60))
  
  pathway_long <- plotdf %>%
    dplyr::select(Description, trait, GeneRatio_percent, Count, p.adjust) %>%
    complete(Description, trait, fill = list(GeneRatio_percent = 0, Count = 0, p.adjust = 1))
  
  count_range <- range(pathway_long$Count[pathway_long$Count > 0], na.rm = TRUE)
  
  if(length(count_range) == 0 || is.infinite(count_range[1])) {
    size_range <- c(4, 12)
  } else if(count_range[2] <= 3) {
    size_range <- c(5, 14)
  } else if(count_range[2] <= 5) {
    size_range <- c(4, 12)
  } else if(count_range[2] <= 10) {
    size_range <- c(3, 10)
  } else {
    size_range <- c(2, 8)
  }
  
  p <- ggplot(pathway_long, aes(x = trait, y = reorder(Description, GeneRatio_percent), 
                                fill = GeneRatio_percent, size = Count)) +
    geom_point(shape = 21, alpha = 0.8) +
    scale_fill_gradient2(name = "Gene Ratio (%)", low = "white", high = "#d73027") +
    scale_size_continuous(name = "Gene Count", range = size_range) +
    labs(x = "Phenotype", y = "Pathway", title = title_text) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 12, hjust = 1),
          axis.title.x = element_text(size = 13, face = "bold"),
          axis.title.y = element_text(size = 13, face = "bold"),
          plot.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12, face = "bold"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  
  return(p)
}

# 创建横向合并的热图函数
create_combined_heatmap <- function(go_bp_df, kegg_df, reactome_df, title_prefix, 
                                   suffixes = c("GO-BP", "KEGG", "REACTOME")) {
  heatmap_list <- list()
  
  # 创建GO-BP热图
  if(!is.null(go_bp_df) && nrow(go_bp_df) > 0) {
    h1 <- create_heatmap_plot(go_bp_df, paste(title_prefix, "-", suffixes[1], "Heatmap"), suffixes[1])
    if(!is.null(h1)) {
      h1 <- h1 + theme(legend.position = "none", axis.text.y = element_text(size = 10))
      heatmap_list[["GO_BP"]] <- h1
    }
  }
  
  # 创建KEGG热图
  if(!is.null(kegg_df) && nrow(kegg_df) > 0) {
    h2 <- create_heatmap_plot(kegg_df, paste(title_prefix, "-", suffixes[2], "Heatmap"), suffixes[2])
    if(!is.null(h2)) {
      h2 <- h2 + theme(legend.position = "none", axis.text.y = element_text(size = 10), 
                       axis.title.y = element_blank())
      heatmap_list[["KEGG"]] <- h2
    }
  }
  
  # 创建REACTOME热图
  if(!is.null(reactome_df) && nrow(reactome_df) > 0) {
    h3 <- create_heatmap_plot(reactome_df, paste(title_prefix, "-", suffixes[3], "Heatmap"), suffixes[3])
    if(!is.null(h3)) {
      h3 <- h3 + theme(axis.text.y = element_text(size = 10), axis.title.y = element_blank())
      heatmap_list[["REACTOME"]] <- h3
    }
  }
  
  # 使用+操作符横向合并
  if(length(heatmap_list) == 1) {
    return(heatmap_list[[1]])
  } else if(length(heatmap_list) == 2) {
    names_list <- names(heatmap_list)
    return(heatmap_list[[names_list[1]]] + heatmap_list[[names_list[2]]])
  } else if(length(heatmap_list) == 3) {
    return(heatmap_list[["GO_BP"]] + heatmap_list[["KEGG"]] + heatmap_list[["REACTOME"]])
  } else {
    return(NULL)
  }
}

# 创建气泡+柱状图组合（参考4.2.bar_bubble_plot_kegg.r）
create_bubble_bar_combined <- function(pathway_data, n_top = 10, title_text = "Pathway Enrichment", 
                                       bubble_color = "#C496CD", bar_color = "#D5C7E1") {
  if(is.null(pathway_data) || nrow(pathway_data) == 0) return(NULL)
  
  # 数据预处理
  pathway_data_df <- as.data.frame(pathway_data)
  pathway_data_df$log10_p_value <- -log10(pathway_data_df$p.adjust)
  
  pathway_data_df <- pathway_data_df %>%
    separate(GeneRatio, into = c("GeneCount", "TotalGenes"), sep = "/", convert = TRUE, remove = FALSE) %>%
    mutate(GeneRatio_percent = GeneCount / TotalGenes * 100) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_top) %>%
    mutate(
      Description = forcats::fct_reorder(Description, Count),
      Description_short = ifelse(nchar(as.character(Description)) > 60, 
                                 paste0(substr(as.character(Description), 1, 60), "..."), 
                                 as.character(Description))
    )
  
  # 创建左侧气泡图
  bubble_plot <- ggplot(pathway_data_df, 
                        aes(x = GeneRatio_percent, y = Description)) +
    geom_point(aes(size = Count), color = bubble_color, alpha = 0.8) +
    scale_size_continuous(
      name = "Gene Count",
      range = c(4, 12),
      breaks = c(5, 10, 15, 20),
      guide = guide_legend(override.aes = list(size = 8, alpha = 1))
    ) +
    guides(
      size = guide_legend(override.aes = list(size = 8))
    ) +
    scale_x_continuous(
      name = "Gene Ratio (%)",
      limits = c(0, max(pathway_data_df$GeneRatio_percent) * 1.1),
      expand = expansion(mult = c(0.02, 0.1))
    ) +
    scale_y_discrete(name = "") +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.size = grid::unit(1.0, "cm"),
      legend.spacing = grid::unit(0.5, "cm"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      plot.margin = margin(1, 0, 1, 1, "cm")
    )
  
  # 创建右侧柱状图
  bar_plot <- ggplot(pathway_data_df, 
                     aes(x = log10_p_value, y = Description)) +
    geom_col(fill = bar_color, alpha = 0.7, width = 0.4) +
    geom_text(aes(label = Description_short, x = 0.1), 
              color = "black", size = 4.5, hjust = 0, fontface = "bold") +
    scale_x_continuous(
      name = "-Log10(P.adjust)",
      expand = expansion(mult = c(0, 0.1))
    ) +
    scale_y_discrete(name = "") +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      plot.margin = margin(1, 1, 1, 0, "cm")
    )
  
  # 使用patchwork组合两个图 - 气泡图占1份，柱状图占4份
  combined_plot <- bubble_plot + bar_plot + plot_layout(widths = c(1, 4))
  
  # 添加标题
  final_plot <- combined_plot + 
    plot_annotation(
      title = title_text,
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    )
  
  return(final_plot)
}

#**************************************     主循环    ************************
# 针对每个模型，分析三个层级的血压蛋白
#*******************************************************************************

for (model in c("model_1","model_2","model_3")){
  cat("\n==============================================\n")
  cat("Processing model:", model, "\n")
  cat("==============================================\n\n")
  
  # 创建pathway文件夹
  pathway_dir <- paste0("D:/OneDrive/工作/1.工作/blood pressure/", model, "/pathway")
  dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 直接读取BP_three_levels_protein_lists.xlsx文件
  bp_levels_file <- paste0("D:/OneDrive/工作/1.工作/blood pressure/", model, "/upset_venn/BP_three_levels_protein_lists.xlsx")
  
  if(!file.exists(bp_levels_file)) {
    cat("Warning: BP_three_levels_protein_lists.xlsx not found, skipping", model, "\n")
    next
  }
  
  cat("Reading BP protein lists from:", bp_levels_file, "\n\n")
  
  # 存储三个层级的Entrez ID列表（只分析Level1和Level2）
  level_genes <- list()
  level_names <- c("Level1_Regression", "Level2_Regression_MR")
  sheet_names <- c("1_Regression", "2_Regression_MR")
  
  # 读取两个层级的数据并转换为Entrez ID
  for(i in 1:2) {
    tryCatch({
      level_data <- openxlsx::read.xlsx(bp_levels_file, sheet = sheet_names[i])
      if(nrow(level_data) > 0) {
        gene_symbols <- level_data$Entrez.Gene.Name
        gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
        if(length(gene_symbols) > 0) {
          gene_entrez_ids <- bitr(geneID = gene_symbols, fromType = "SYMBOL",
                                  toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
          level_genes[[level_names[i]]] <- gene_entrez_ids$ENTREZID
          cat(sprintf("%s: %d genes\n", level_names[i], length(gene_entrez_ids$ENTREZID)))
        }
      }
    }, error = function(e) {
      cat(sprintf("Warning: Could not read %s data\n", level_names[i]))
    })
  }
  
  cat("\n")
  
  if(length(level_genes) == 0) {
    cat("No gene lists loaded, skipping", model, "\n")
    next
  }
  
  # 保存Entrez ID列表
  gene_list_df <- data.frame()
  for(name in names(level_genes)) {
    if(length(level_genes[[name]]) > 0) {
      temp_df <- data.frame(
        Level = name,
        Entrez_ID = level_genes[[name]],
        stringsAsFactors = FALSE
      )
      gene_list_df <- rbind(gene_list_df, temp_df)
    }
  }
  
  if(nrow(gene_list_df) > 0) {
    wb_gene_list <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_gene_list, "Gene_Lists")
    openxlsx::writeData(wb_gene_list, "Gene_Lists", gene_list_df)
    openxlsx::saveWorkbook(wb_gene_list, file.path(pathway_dir, "gene_lists_bp_entrez.xlsx"), overwrite = TRUE)
    cat("Entrez gene lists saved to:", file.path(pathway_dir, "gene_lists_bp_entrez.xlsx"), "\n\n")
  }
  
  # 对每个层级进行通路分析
  all_results <- list()
  all_bubble_bar_plots <- list()
  figure_counter <- 1  # 初始化图号计数器
  
  for(level_name in names(level_genes)) {
    cat(sprintf("\n========== Analyzing %s ==========\n", level_name))
    
    results <- perform_pathway_analysis(level_genes[[level_name]], level_name)
    all_results[[level_name]] <- results
    
    # 创建气泡+柱状图组合（每个数据库单独保存）
    bubble_bar_list <- list()
    
    # Reactome 气泡+柱状图
    if(!is.null(results$reactome) && nrow(results$reactome) > 0) {
      reactome_bb <- create_bubble_bar_combined(
        results$reactome, 
        n_top = 10,
        title_text = sprintf("%s - Reactome Pathways", level_name),
        bubble_color = "#FFB366",
        bar_color = "#FFD4A3"
      )
      if(!is.null(reactome_bb)) {
        bubble_bar_list[["Reactome"]] <- reactome_bb
        ggsave(file.path(pathway_dir, sprintf("%d_%s_Reactome_bubble_bar.pdf", figure_counter, level_name)),
               reactome_bb, width = 14, height = 8, device = "pdf")
        figure_counter <- figure_counter + 1
      }
    }
    
    # KEGG 气泡+柱状图
    if(!is.null(results$kegg) && nrow(results$kegg) > 0) {
      kegg_bb <- create_bubble_bar_combined(
        results$kegg, 
        n_top = 10,
        title_text = sprintf("%s - KEGG Pathways", level_name),
        bubble_color = "#C496CD",
        bar_color = "#D5C7E1"
      )
      if(!is.null(kegg_bb)) {
        bubble_bar_list[["KEGG"]] <- kegg_bb
        ggsave(file.path(pathway_dir, sprintf("%d_%s_KEGG_bubble_bar.pdf", figure_counter, level_name)),
               kegg_bb, width = 14, height = 8, device = "pdf")
        figure_counter <- figure_counter + 1
      }
    }
    
    # GO-BP 气泡+柱状图
    if(!is.null(results$go_bp) && nrow(results$go_bp) > 0) {
      go_bp_bb <- create_bubble_bar_combined(
        results$go_bp, 
        n_top = 10,
        title_text = sprintf("%s - GO Biological Process", level_name),
        bubble_color = "#F09898",
        bar_color = "#FFB6B9"
      )
      if(!is.null(go_bp_bb)) {
        bubble_bar_list[["GO_BP"]] <- go_bp_bb
        ggsave(file.path(pathway_dir, sprintf("%d_%s_GO_BP_bubble_bar.pdf", figure_counter, level_name)),
               go_bp_bb, width = 14, height = 8, device = "pdf")
        figure_counter <- figure_counter + 1
      }
    }
    
    # GO-CC 气泡+柱状图
    if(!is.null(results$go_cc) && nrow(results$go_cc) > 0) {
      go_cc_bb <- create_bubble_bar_combined(
        results$go_cc, 
        n_top = 10,
        title_text = sprintf("%s - GO Cellular Component", level_name),
        bubble_color = "#7AA8CC",
        bar_color = "#A8C5E0"
      )
      if(!is.null(go_cc_bb)) {
        bubble_bar_list[["GO_CC"]] <- go_cc_bb
        ggsave(file.path(pathway_dir, sprintf("%d_%s_GO_CC_bubble_bar.pdf", figure_counter, level_name)),
               go_cc_bb, width = 14, height = 8, device = "pdf")
        figure_counter <- figure_counter + 1
      }
    }
    
    # GO-MF 气泡+柱状图
    if(!is.null(results$go_mf) && nrow(results$go_mf) > 0) {
      go_mf_bb <- create_bubble_bar_combined(
        results$go_mf, 
        n_top = 10,
        title_text = sprintf("%s - GO Molecular Function", level_name),
        bubble_color = "#9ACC7A",
        bar_color = "#C5E0A8"
      )
      if(!is.null(go_mf_bb)) {
        bubble_bar_list[["GO_MF"]] <- go_mf_bb
        ggsave(file.path(pathway_dir, sprintf("%d_%s_GO_MF_bubble_bar.pdf", figure_counter, level_name)),
               go_mf_bb, width = 14, height = 8, device = "pdf")
        figure_counter <- figure_counter + 1
      }
    }
    
    # 保存该层级的所有气泡+柱状图
    all_bubble_bar_plots[[level_name]] <- bubble_bar_list
    
    cat(sprintf("Bubble+bar combined plots saved for %s\n", level_name))
  }
  
  # 创建跨层级的组合图
  cat("\n========== Creating cross-level combined plots ==========\n")
  
  # 为每个数据库创建三个层级的垂直组合图
  for(db_type in c("Reactome", "KEGG", "GO_BP", "GO_CC", "GO_MF")) {
    plots_to_combine <- list()
    
    for(level_name in names(all_bubble_bar_plots)) {
      if(!is.null(all_bubble_bar_plots[[level_name]][[db_type]])) {
        plots_to_combine[[level_name]] <- all_bubble_bar_plots[[level_name]][[db_type]]
      }
    }
    
    if(length(plots_to_combine) > 0) {
      # 垂直组合
      if(length(plots_to_combine) == 2) {
        combined_plot <- plots_to_combine[[1]] / plots_to_combine[[2]] +
          plot_annotation(
            title = sprintf("BP Two Levels - %s Comparison", gsub("_", "-", db_type)),
            subtitle = "Level 1: Regression | Level 2: Regression+MR",
            theme = theme(
              plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40")
            )
          )
      } else {
        combined_plot <- plots_to_combine[[1]]
      }
      
      ggsave(file.path(pathway_dir, sprintf("%d_All_Levels_%s_bubble_bar_combined.pdf", figure_counter, db_type)),
             combined_plot, width = 16, height = 8 * length(plots_to_combine), device = "pdf")
      figure_counter <- figure_counter + 1
      
      cat(sprintf("  Combined %s plot saved\n", db_type))
    }
  }
  
  #**************************************     保存Excel结果    ************************
  cat("\n========== Saving Excel Results ==========\n")
  
  # 创建Excel工作簿
  wb <- openxlsx::createWorkbook()
  
  for(level_name in names(all_results)) {
    results <- all_results[[level_name]]
    
    # 保存每个数据库的结果
    if(!is.null(results$go_bp) && nrow(results$go_bp) > 0) {
      sheet_name <- paste0(gsub("Level[0-9]_", "", level_name), "_GO_BP")
      if(nchar(sheet_name) > 31) sheet_name <- substr(sheet_name, 1, 31)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, results$go_bp)
    }
    
    if(!is.null(results$go_cc) && nrow(results$go_cc) > 0) {
      sheet_name <- paste0(gsub("Level[0-9]_", "", level_name), "_GO_CC")
      if(nchar(sheet_name) > 31) sheet_name <- substr(sheet_name, 1, 31)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, results$go_cc)
    }
    
    if(!is.null(results$go_mf) && nrow(results$go_mf) > 0) {
      sheet_name <- paste0(gsub("Level[0-9]_", "", level_name), "_GO_MF")
      if(nchar(sheet_name) > 31) sheet_name <- substr(sheet_name, 1, 31)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, results$go_mf)
    }
    
    if(!is.null(results$kegg) && nrow(results$kegg) > 0) {
      sheet_name <- paste0(gsub("Level[0-9]_", "", level_name), "_KEGG")
      if(nchar(sheet_name) > 31) sheet_name <- substr(sheet_name, 1, 31)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, results$kegg)
    }
    
    if(!is.null(results$reactome) && nrow(results$reactome) > 0) {
      sheet_name <- paste0(gsub("Level[0-9]_", "", level_name), "_Reactome")
      if(nchar(sheet_name) > 31) sheet_name <- substr(sheet_name, 1, 31)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, results$reactome)
    }
  }
  
  # 保存Excel文件（只在有数据时保存）
  if(length(names(wb)) > 0) {
    excel_file <- file.path(pathway_dir, "BP_pathway_results_all.xlsx")
    openxlsx::saveWorkbook(wb, excel_file, overwrite = TRUE)
    cat("Excel results saved:", excel_file, "\n\n")
  } else {
    cat("No pathway results to save\n\n")
  }
  
  cat("==============================================\n")
  cat("Model", model, "processing completed!\n")
  cat("==============================================\n\n")
}

cat("\n\n")
cat("##########################################################\n")
cat("# All BP pathway analysis completed!                     #\n")
cat("##########################################################\n")
