library("dplyr")
library("ggplot2")
library("ggrepel")
library("gridExtra")
library("grid")
library("stringr")

select <- dplyr::select

# 设置路径
base_dir <- "D:/OneDrive/工作/1.工作/blood pressure"

# Blood pressure 相关表型（包括pp和map）
bp_traits <- c("sbp", "dbp", "hypertension", "pp", "map")

# 创建输出文件夹
for (model in c("model_1", "model_2", "model_3")) {
  # 设置模型特定的路径
  bp_volcano_dir <- file.path(base_dir, model, "volcano_plots")
  if(!dir.exists(bp_volcano_dir)) dir.create(bp_volcano_dir, recursive = TRUE)
  
  cat(paste("\n", rep("=", 70), collapse = ""), "\n")
  cat(paste("Processing Blood Pressure Volcano Plots for", model, "\n"))
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # 定义读取结果文件的函数
  read_result <- function(trait, model){
    # 从obs_all_results.txt文件中读取数据
    file_path <- file.path(base_dir, model, "obs_all_results.txt")
    if(file.exists(file_path)){
      tryCatch({
        # 读取整个文件
        all_res <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        # 筛选特定trait的数据
        res <- all_res %>% filter(traits == trait)
        
        if(nrow(res) == 0) {
          warning(paste("No data found for trait:", trait))
          return(NULL)
        }
        
        cat(paste("  Loaded", nrow(res), "records for", trait, "\n"))
        return(res)
      }, error = function(e){
        warning(paste("Failed to read file for trait", trait, ":", e$message))
        return(NULL)
      })
    } else {
      warning(paste("File not found at path:", file_path))
      return(NULL)
    }
  }
  
  # ========== 第一个图：按照原来的逻辑绘制组火山图（SBP-based） ==========
  cat("Generating BP group volcano plot (SBP-based)...\n")
  
  plot_bp_group <- function(traits, out_dir, model){
    primary_trait <- "sbp"  # 以SBP为主要表型
    
    # 读取每个trait的结果
    res_list <- lapply(traits, read_result, model = model)
    names(res_list) <- traits
    
    # 检查主要表型是否有数据
    primary_data <- res_list[[primary_trait]]
    if(is.null(primary_data) || nrow(primary_data) == 0) {
      warning(paste("No data available for primary trait", primary_trait))
      return(NULL)
    }
    
    # 构建主表，以主要表型为基础
    if("Entrez.Gene.Name" %in% colnames(primary_data)) {
      master <- primary_data %>% 
        dplyr::select(seqName, Entrez.Gene.Name, Estimate, adjust_P)
    } else {
      master <- primary_data %>% 
        dplyr::select(seqName, Estimate, adjust_P) %>%
        mutate(Entrez.Gene.Name = NA)
    }
    
    # 计算主要表型的显著性和效应方向
    master$negLog10AdjP <- -log10(master$adjust_P + 1e-300)
    master$is_primary_signif <- master$adjust_P < 0.05
    master$effect_direction <- ifelse(master$Estimate > 0, "Positive", "Negative")
    
    # 计算其他表型的显著性支持度
    master$other_signif_count <- 0
    
    for(tr in traits) {
      if(tr == primary_trait) next
      
      other_data <- res_list[[tr]]
      if(is.null(other_data) || nrow(other_data) == 0) {
        cat(paste("  - Skipping trait", tr, "(no data)\n"))
        next
      }
      
      # 合并其他表型的显著性信息
      other_signif <- other_data %>% 
        dplyr::select(seqName, !!paste0("adjustP_", tr) := adjust_P) %>%
        mutate(!!paste0("signif_", tr) := !!sym(paste0("adjustP_", tr)) < 0.05)
      
      master <- master %>% 
        left_join(other_signif, by = "seqName") %>%
        mutate(other_signif_count = other_signif_count + 
                 ifelse(is.na(!!sym(paste0("signif_", tr))), 0, as.numeric(!!sym(paste0("signif_", tr)))))
      
      cat(paste("  - Added support from trait", tr, "\n"))
    }
    
    # 计算总显著次数（包括主要表型）
    master$total_signif_count <- ifelse(master$is_primary_signif, 1, 0) + master$other_signif_count
    
    # 创建综合的颜色分组
    max_total_count <- max(master$total_signif_count, na.rm = TRUE)
    n_total_traits <- length(traits)
    
    # 设置颜色方案
    if(max_total_count == 0) {
      master$color_group <- "NS"
      palette_colors <- c("NS" = "grey70")
      palette_labels <- c("NS" = "Not Significant")
    } else {
      master$color_group <- ifelse(master$total_signif_count == 0, "NS",
                                    paste0(ifelse(master$effect_direction == "Positive", "Pos", "Neg"), 
                                           "_", master$total_signif_count))
      
      # 红色系：正相关
      pos_colors <- colorRampPalette(c("#FAE3D9", "#F8D4CC", "#F6C5BF", "#F4B6B2", "#F2A7A5", "#F09898", "#EE898B", "#EC7A7E", "#EA6B71", "#E85C64"))(max_total_count)
      # 蓝色系：负相关
      neg_colors <- colorRampPalette(c("#EFF6FC", "#E0EDFA", "#C8E0F7", "#A8D1F2", "#88C2ED", "#6AB3E8", "#4FA4E3", "#3495DE", "#1A86D9", "#0077D4"))(max_total_count)
      if(max_total_count > 0) neg_colors[1] <- "#E0EDFA"
      
      palette_colors <- c("NS" = "grey70")
      palette_labels <- c("NS" = "Not Significant")
      
      for(i in 1:max_total_count) {
        pos_key <- paste0("Pos_", i)
        neg_key <- paste0("Neg_", i)
        palette_colors[pos_key] <- pos_colors[i]
        palette_colors[neg_key] <- neg_colors[i]
        palette_labels[pos_key] <- paste0("Positive (", i, "/", n_total_traits, ")")
        palette_labels[neg_key] <- paste0("Negative (", i, "/", n_total_traits, ")")
      }
    }
    
    # 创建图形
    p <- ggplot(master, aes(x = Estimate, y = negLog10AdjP, color = color_group)) +
      geom_point(alpha = 0.7, size = 2.5) +
      scale_color_manual(values = palette_colors, 
                         labels = palette_labels,
                         name = "Significance") +
      theme_bw() +
      labs(title = "Blood Pressure Group (SBP-based)", 
           x = "Effect Size (SBP)", 
           y = "-log10(Adjusted P-value, SBP)") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 9),
            legend.position = "right",
            panel.border = element_blank(),
            axis.line = element_line(color = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    # 添加显著性阈线
    p <- p + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7)
    
    # 标注显著蛋白
    if(any(master$total_signif_count > 0)) {
      max_support <- max(master$total_signif_count, na.rm = TRUE)
      
      for(support_level in max_support:max(1, max_support-1)) {
        top_hits <- master %>% 
          filter(total_signif_count == support_level) %>% 
          arrange(-negLog10AdjP) %>% 
          head(ifelse(support_level == max_support, 15, 10))
        
        if(nrow(top_hits) > 0) {
          label_col <- ifelse(!is.na(top_hits$Entrez.Gene.Name) & 
                                top_hits$Entrez.Gene.Name != "" & 
                                top_hits$Entrez.Gene.Name != "NA", 
                              top_hits$Entrez.Gene.Name, 
                              top_hits$seqName)
          
          p <- p + geom_text_repel(data = top_hits, 
                                   aes(label = label_col), 
                                   size = 3.2, 
                                   box.padding = 0.3,
                                   point.padding = 0.2,
                                   max.overlaps = 25,
                                   force = 2,
                                   min.segment.length = 0.1)
          break
        }
      }
    }
    
    # 保存图形
    out_file_png <- file.path(out_dir, "BP_group_sbp_based_volcano.png")
    ggsave(out_file_png, p, width = 10, height = 7, dpi = 300, bg = "white")
    
    out_file_pdf <- file.path(out_dir, "BP_group_sbp_based_volcano.pdf")
    ggsave(out_file_pdf, p, width = 10, height = 7, device = "pdf", bg = "white")
    
    cat(paste("\nBP group volcano plot saved:\n"))
    cat(paste("  PNG:", out_file_png, "\n"))
    cat(paste("  PDF:", out_file_pdf, "\n"))
    cat(paste("  Total proteins:", nrow(master), "\n"))
    cat(paste("  Proteins significant in SBP:", sum(master$is_primary_signif, na.rm = TRUE), "\n"))
    cat(paste("  Max total trait consistency:", max(master$total_signif_count, na.rm = TRUE), "out of", n_total_traits, "\n\n"))
    
    return(list(png = out_file_png, pdf = out_file_pdf, plot = p))
  }
  
  bp_group_result <- plot_bp_group(bp_traits, bp_volcano_dir, model)
  
  
  # ========== 第二个图：分不同子表型绘制火山图 ==========
  cat("Generating individual BP sub-phenotype volcano plots...\n")
  
  plot_volcano_single <- function(res, trait, out_dir){
    if(is.null(res) || nrow(res) == 0) {
      warning(paste("No data available for trait:", trait))
      return(NULL)
    }
    
    # 检查必要的列
    required_cols <- c("Estimate", "adjust_P", "seqName")
    missing_cols <- setdiff(required_cols, colnames(res))
    if(length(missing_cols) > 0){
      warning(paste("Missing required columns for trait", trait, ":", paste(missing_cols, collapse = ", ")))
      return(NULL)
    }
    
    res$negLog10AdjP <- -log10(res$adjust_P + 1e-300)
    
    # 根据显著性和效应方向设置颜色分组
    res$color_group <- ifelse(res$adjust_P >= 0.05, "NS", 
                              ifelse(res$Estimate > 0, "Positive", "Negative"))
    
    # 使用BP组火山图上最深的红色和蓝色
    p <- ggplot(res, aes(x = Estimate, y = negLog10AdjP, color = color_group)) +
      geom_point(alpha = 0.7, size = 2.5) +
      scale_color_manual(values = c("NS" = "grey40", "Positive" = "#E85C64", "Negative" = "#0077D4"),
                         labels = c("NS" = "Not Significant", "Positive" = "Positive", "Negative" = "Negative"),
                         name = "Association") +
      theme_bw() +
      labs(title = paste0(toupper(trait)), 
           x = "Effect Size", 
           y = "-log10(Adjusted P-value)") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10),
            panel.border = element_blank(),
            axis.line = element_line(color = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    # 标注显著的top 15
    significant_hits <- res %>% filter(adjust_P < 0.05) %>% arrange(adjust_P) %>% head(15)
    if(nrow(significant_hits) > 0){
      label_col <- if("Entrez.Gene.Name" %in% colnames(significant_hits)) {
        ifelse(!is.na(significant_hits$Entrez.Gene.Name) & 
                 significant_hits$Entrez.Gene.Name != "" & 
                 significant_hits$Entrez.Gene.Name != "NA", 
               significant_hits$Entrez.Gene.Name, 
               significant_hits$seqName)
      } else {
        significant_hits$seqName
      }
      
      p <- p + geom_text_repel(data = significant_hits, 
                               aes(label = label_col), 
                               size = 3.0, 
                               box.padding = 0.3,
                               point.padding = 0.3,
                               max.overlaps = 20,
                               force = 1.5)
    }
    
    # 添加显著性阈线
    p <- p + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7)
    
    return(p)
  }
  
  # 为每个子表型绘图并单独保存PDF
  bp_subplots <- list()
  for(trait in bp_traits){
    cat(paste("  Processing trait:", trait, "\n"))
    res <- read_result(trait, model)
    p <- plot_volcano_single(res, trait, bp_volcano_dir)
    
    # 保存单独的PDF和PNG
    if(!is.null(p)) {
      single_png <- file.path(bp_volcano_dir, paste0(trait, "_volcano.png"))
      single_pdf <- file.path(bp_volcano_dir, paste0(trait, "_volcano.pdf"))
      
      ggsave(single_png, p, width = 8, height = 6, dpi = 300, bg = "white")
      ggsave(single_pdf, p, width = 8, height = 6, device = "pdf", bg = "white")
      
      cat(paste("    Saved:", single_png, "\n"))
      cat(paste("    Saved:", single_pdf, "\n"))
    }
    
    bp_subplots[[trait]] <- p
  }
  
  # 创建组合图
  cat("\nCreating combined BP sub-phenotype plot...\n")
  
  # 定义更清晰的子图标题
  trait_titles <- list(
    "sbp" = "SBP",
    "dbp" = "DBP",
    "hypertension" = "Hypertension",
    "pp" = "PP",
    "map" = "MAP"
  )
  
  # 修改每个子图的标题
  for(trait in names(bp_subplots)) {
    if(!is.null(bp_subplots[[trait]]) && trait %in% names(trait_titles)) {
      bp_subplots[[trait]] <- bp_subplots[[trait]] + 
        labs(title = trait_titles[[trait]]) +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    }
  }
  
  # 移除NULL值
  valid_subplots <- bp_subplots[!sapply(bp_subplots, is.null)]
  
  if(length(valid_subplots) >= 2) {
    # 根据子图数量调整布局（5个子图用2行布局）
    n_plots <- length(valid_subplots)
    if(n_plots <= 3) {
      ncol_val <- n_plots
      nrow_val <- 1
      width_val <- 6 * n_plots
      height_val <- 6
    } else {
      ncol_val <- 3
      nrow_val <- ceiling(n_plots / 3)
      width_val <- 18
      height_val <- 6 * nrow_val
    }
    
    combined_bp_plot <- grid.arrange(
      grobs = valid_subplots,
      ncol = ncol_val, nrow = nrow_val,
      top = textGrob(paste0("Blood Pressure Sub-phenotypes (", model, ")"), 
                     gp = gpar(fontsize = 16, fontface = "bold"))
    )
    
    # 保存组合图
    combined_png <- file.path(bp_volcano_dir, "BP_subphenotypes_combined.png")
    combined_pdf <- file.path(bp_volcano_dir, "BP_subphenotypes_combined.pdf")
    
    ggsave(combined_png, combined_bp_plot, width = width_val, height = height_val, dpi = 300, bg = "white")
    ggsave(combined_pdf, combined_bp_plot, width = width_val, height = height_val, device = "pdf", bg = "white")
    
    cat(paste("\nBP sub-phenotype combined plot saved:\n"))
    cat(paste("  PNG:", combined_png, "\n"))
    cat(paste("  PDF:", combined_pdf, "\n"))
  } else {
    warning("Not enough valid plots to create combined figure")
  }
  
  # 输出总结
  cat(paste("\n", rep("=", 70), collapse = ""), "\n")
  cat(paste("Blood Pressure Volcano Plots Completed for", model, "\n"))
  cat(paste("All plots saved to:", bp_volcano_dir, "\n"))
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
}

cat("\n*** All Blood Pressure volcano plots generated successfully! ***\n")
