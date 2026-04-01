# =============================================================================
# 血压相关蛋白质-表型热图分析 (Blood Pressure Specific)
# 功能：分析OBS和MR结果的一致性，生成热图可视化
# 表型：SBP, DBP, PP, MAP, Hypertension
# 逻辑：1.读入数据 -> 2.筛选蛋白 -> 3.绘图
# =============================================================================

library(dplyr)
library(tidyr)
library(xlsx)
library(ggplot2)
library(reshape2)
library(openxlsx)
library(tibble)
library(RColorBrewer)
library(circlize)  # 用于环状热图

# 设置血压专用路径
base_dir <- "D:/OneDrive/工作/1.工作/blood pressure"

# 定义血压相关trait和MR映射
trait_mapping <- data.frame(
  obs_trait = c("sbp", "dbp", "pp", "map", "hypertension"),
  mr_trait = c("ckb_sbp", "ckb_dbp", "ckb_pp", "ckb_map", "ukb_ht"),
  group = c("blood_pressure", "blood_pressure", "blood_pressure", "blood_pressure", "blood_pressure"),
  stringsAsFactors = FALSE
)

# 添加获取显著性信息的函数
get_significance_info <- function(merged_data, protein_id, obs_trait, mr_trait) {
  # 从merged数据中获取该蛋白质的记录
  protein_row <- merged_data[merged_data$protein_id == protein_id, ]
  
  if (nrow(protein_row) == 0) {
    return(list(obs_significant = FALSE, union_significant = FALSE))
  }
  
  # 构建obs和union列名
  obs_col <- paste0(obs_trait, ".obs")
  union_col <- paste0(mr_trait, ".union")
  
  # 获取显著性信息
  obs_significant <- FALSE
  union_significant <- FALSE
  
  if (obs_col %in% colnames(protein_row)) {
    obs_significant <- as.logical(protein_row[[obs_col]])
    if (is.na(obs_significant)) obs_significant <- FALSE
  }
  
  if (!is.na(mr_trait) && union_col %in% colnames(protein_row)) {
    union_significant <- as.logical(protein_row[[union_col]])
    if (is.na(union_significant)) union_significant <- FALSE
  }
  
  return(list(obs_significant = obs_significant, union_significant = union_significant))
}

# 创建热图绘制函数
create_heatmap_plot <- function(heatmap_data, model_name, level_name) {
  cat("创建血压专用热图\n")
  
  if (is.null(heatmap_data) || nrow(heatmap_data) == 0) {
    cat("无数据可用于绘图\n")
    return(NULL)
  }
  
  # 筛选有数据的蛋白质-trait组合
  plot_data <- heatmap_data %>%
    filter(!is.na(heatmap_value) & heatmap_value != 0) %>%
    dplyr::select(Entrez.Gene.Name, trait, group, heatmap_value, consistency, obs_direction, mr_direction)
  
  if (nrow(plot_data) == 0) {
    cat("无有效数据用于热图\n")
    return(NULL)
  }
  
  # 定义血压性状顺序
  trait_order <- data.frame(
    trait = c("sbp", "dbp", "pp", "map", "hypertension"),
    group = c("blood_pressure", "blood_pressure", "blood_pressure", "blood_pressure", "blood_pressure"),
    stringsAsFactors = FALSE
  )
  
  # 为trait创建排序因子
  trait_levels <- trait_order$trait
  plot_data$trait <- factor(plot_data$trait, levels = trait_levels)
  
  # 按蛋白质-trait组合数排序
  protein_importance <- plot_data %>%
    group_by(Entrez.Gene.Name) %>%
    summarise(
      trait_count = n(),
      consistency_score = sum(consistency == "both_significant_consistent", na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(trait_count), desc(consistency_score))
  
  plot_data$Entrez.Gene.Name <- factor(plot_data$Entrez.Gene.Name, 
                                       levels = protein_importance$Entrez.Gene.Name)
  
  # 定义离散颜色值和标签
  discrete_values <- c(-1.0, -0.6, -0.3, -0.2, 0, 0.2, 0.3, 0.6, 1.0)
  discrete_colors <- c(
    "#3182BD",  # -1.0 (深蓝)
    "#6BAED6",  # -0.6 (中蓝)
    "#9ECAE1",  # -0.3 (浅中蓝)
    "#DEEBF7",  # -0.2 (极浅蓝)
    "#F7F7F7",  # 0 (白色)
    "#FCBCBB",  # 0.2 (极浅红)
    "#FA9C9A",  # 0.3 (浅红)
    "#F26D6C",  # 0.6 (中红)
    "#E34544"   # 1.0 (深红)
  )
  
  categorical_labels <- c(
    "both significant -",
    "both significant, regression -, MR +",
    "only regression significant, -",
    "only MR significant, -",
    "not significant",
    "only MR significant, +",
    "only regression significant, +",
    "both significant, regression +, MR -",
    "both significant +"
  )
  
  # 将数值映射为分类标签
  plot_data <- plot_data %>%
    mutate(heatmap_category = factor(heatmap_value, 
                                     levels = discrete_values, 
                                     labels = categorical_labels))
  
  plot_data$heatmap_category <- factor(plot_data$heatmap_category, levels = categorical_labels)
  
  # 创建热图（横向布局）
  p <- ggplot(plot_data, aes(x = Entrez.Gene.Name, y = trait, fill = heatmap_category)) +
    geom_tile(color = "white", linewidth = 0.5, width = 0.98, height = 0.98) +
    scale_fill_manual(
      values = setNames(discrete_colors, categorical_labels),
      na.value = "#f7f7f7",
      name = "Significance Pattern",
      drop = FALSE,
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        ncol = 1,
        byrow = FALSE,
        override.aes = list(size = 0.5),
        keywidth = unit(0.3, "cm"),
        keyheight = unit(0.3, "cm"),
        label.theme = element_text(size = 6),
        title.theme = element_text(size = 7)
      )
    ) +
    labs(
      title = paste("Blood Pressure Protein-Trait Associations -", model_name, "-", level_name),
      x = "Protein (Gene Symbol)",
      y = "Blood Pressure Traits"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10, face = "bold", hjust = 1),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.x = unit(0.3, "cm"),
      legend.spacing.y = unit(0.05, "cm"),
      legend.margin = margin(t = 2, b = 2),
      legend.box.margin = margin(t = 0, b = 0),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7, face = "bold"),
      legend.key.spacing.y = unit(0.05, "cm"),
      legend.key.size = unit(0.3, "cm"),
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold", margin = margin(b = 3)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = 5, r = 20, b = 70, l = 20, unit = "pt"),
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 15))
    ) +
    scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
    coord_cartesian(xlim = c(0.5, length(unique(plot_data$Entrez.Gene.Name)) + 0.5))
  
  return(p)
}

# 创建环状热图函数（带染色体信息）
create_circular_heatmap <- function(heatmap_data, model_name, level_name, gene_pos = NULL) {
  cat("创建血压专用环状热图\n")
  
  if (is.null(heatmap_data) || nrow(heatmap_data) == 0) {
    cat("无数据可用于环状热图\n")
    return(NULL)
  }
  
  # 处理数据格式
  if ("heatmap_value" %in% colnames(heatmap_data)) {
    # 长格式数据处理
    plot_data <- heatmap_data %>%
      filter(!is.na(heatmap_value) & heatmap_value != 0)
    if (nrow(plot_data) == 0) {
      cat("无有效数据用于环状热图\n")
      return(NULL)
    }
    mat <- plot_data %>%
      dplyr::select(Entrez.Gene.Name, trait, heatmap_value) %>%
      mutate(
        Entrez.Gene.Name = as.character(Entrez.Gene.Name),
        trait = as.character(trait),
        heatmap_value = as.numeric(heatmap_value)
      ) %>%
      distinct(Entrez.Gene.Name, trait, .keep_all = TRUE) %>%
      pivot_wider(names_from = trait, values_from = heatmap_value, values_fill = 0) %>%
      column_to_rownames("Entrez.Gene.Name") %>%
      as.matrix()
  } else {
    # 宽格式数据处理
    mat <- heatmap_data
    if ("Entrez.Gene.Name" %in% colnames(mat)) {
      mat <- mat %>% dplyr::distinct(Entrez.Gene.Name, .keep_all = TRUE) %>% as.data.frame()
      rownames(mat) <- mat$Entrez.Gene.Name
      mat <- as.matrix(mat[, !(colnames(mat) %in% c("Entrez.Gene.Name", "protein_id")), drop = FALSE])
    } else {
      mat <- as.matrix(mat)
    }
  }
  
  # 确保矩阵有数据
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    cat("转换后的矩阵为空\n")
    return(NULL)
  }
  
  # 按染色体对蛋白质排序
  chr_info <- NULL
  if (!is.null(gene_pos) && nrow(gene_pos) > 0) {
    cat("按染色体对蛋白质排序\n")
    
    # 创建蛋白质-染色体映射
    protein_chr <- data.frame(
      protein = rownames(mat),
      stringsAsFactors = FALSE
    ) %>%
      left_join(gene_pos[, c("gene", "chr", "chr_num")], 
                by = c("protein" = "gene"))
    
    # 分类：有染色体信息 vs 无染色体信息
    proteins_with_chr <- protein_chr %>% filter(!is.na(chr_num))
    proteins_without_chr <- protein_chr %>% filter(is.na(chr_num))
    
    if (nrow(proteins_with_chr) > 0) {
      # 只按染色体排序，不需要精确位置
      proteins_with_chr <- proteins_with_chr %>%
        arrange(chr_num)
      
      # 合并排序：先有染色体信息的，再无染色体信息的
      protein_order <- c(proteins_with_chr$protein, proteins_without_chr$protein)
      
      # 重新排列矩阵
      mat <- mat[protein_order, ]
      
      # 保存染色体信息用于绘制
      chr_info <- proteins_with_chr
      
      cat(sprintf("按染色体排序完成：%d个蛋白有染色体信息\n", nrow(proteins_with_chr)))
    }
  }
  
  # 调试：检查矩阵内容
  cat("\n=== 调试：环状热图矩阵信息 ===\n")
  cat(sprintf("矩阵维度: %d行 x %d列\n", nrow(mat), ncol(mat)))
  cat(sprintf("列名: %s\n", paste(colnames(mat), collapse = ", ")))
  cat(sprintf("非零值数量: %d\n", sum(mat != 0, na.rm = TRUE)))
  cat("矩阵值统计:\n")
  print(summary(as.vector(mat)))
  cat("前5行前5列:\n")
  print(mat[1:min(5, nrow(mat)), 1:min(5, ncol(mat))])
  cat("\n")
  
  # 绘制环状热图
  tryCatch({
    circos.clear()
    
    circos.par(
      start.degree = 0,
      gap.after = 30,
      track.margin = c(0.01, 0.01),
      cell.padding = c(0, 0, 0, 0),
      canvas.xlim = c(-1.0, 1.0),
      canvas.ylim = c(-1.0, 1.0),
      points.overflow.warning = FALSE
    )
    
    discrete_values <- c(-1.0, -0.6, -0.3, -0.2, 0, 0.2, 0.3, 0.6, 1.0)
    discrete_colors <- c(
      "#3182BD", "#6BAED6", "#9ECAE1", "#DEEBF7", "#F7F7F7",
      "#FCBCBB", "#FA9C9A", "#F26D6C", "#E34544"
    )
    
    categorical_labels <- c(
      "both significant -",
      "both significant, regression -, MR +",
      "only regression significant, -",
      "only MR significant, -",
      "not significant",
      "only MR significant, +",
      "only regression significant, +",
      "both significant, regression +, MR -",
      "both significant +"
    )
    
    col_fun <- colorRamp2(discrete_values, discrete_colors)
    
    n <- nrow(mat)
    name_cex <- if (n > 200) 0.7 else if (n > 110) 0.9 else 1.1
    
    split_all <- rep("all", n)
    sector_name <- "all"
    
    circos.heatmap(
      mat, 
      col = col_fun,
      split = split_all,           
      track.height = 0.25,
      rownames.side = "outside",
      rownames.col = "#2C3E50",
      rownames.cex = name_cex,
      cluster = FALSE,
      dend.side = "none"
    )
    
    track_idx <- get.current.track.index()
    
    # 添加白色边框
    n_traits <- ncol(mat)
    for (i in seq_len(n)) {
      for (j in seq_len(n_traits)) {
        circos.rect(
          xleft = i - 1.0, ybottom = j - 1,
          xright = i + 0.0, ytop = j,
          sector.index = sector_name, track.index = track_idx,
          col = NA, border = "white", lwd = 0.5
        )
      }
    }
    
    # 添加trait标签
    trait_names <- colnames(mat)
    for (j in seq_len(n_traits)) {
      trait_idx <- n_traits - j + 1
      circos.text(
        x = n + 0.5,
        y = j - 0.5,
        labels = trait_names[trait_idx],
        sector.index = sector_name,
        track.index = track_idx,
        facing = "inside",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 0.7,
        col = "#2C3E50",
        font = 2
      )
    }
    
    # 在最内圈添加染色体信息轨道
    if (!is.null(chr_info) && nrow(chr_info) > 0) {
      cat("添加染色体信息轨道\n")
      cat(sprintf("chr_info 行数: %d\n", nrow(chr_info)))
      cat(sprintf("sector_name: %s\n", sector_name))
      
      chr_color_group1 <- "#B2182B"
      chr_color_group2 <- "#2166AC"
      default_chr_color <- "#999999"

      chr_order <- unique(chr_info$chr[!is.na(chr_info$chr) & chr_info$chr != ""])
      chr_color_map <- rep(chr_color_group1, length(chr_order))
      chr_color_map[seq_along(chr_order) %% 2 == 0] <- chr_color_group2
      names(chr_color_map) <- chr_order

      # 保存 sector_name 到局部变量
      current_sector <- sector_name
      
      circos.track(ylim = c(0, 1), track.height = 0.06, bg.border = NA,
                   panel.fun = function(x, y) {
                     cat(sprintf("在 panel.fun 中，chr_info 行数: %d\n", nrow(chr_info)))
                     
                     # 为每个蛋白质绘制染色体标记（无分隔）
                     for (i in seq_len(nrow(chr_info))) {
                        chr <- chr_info$chr[i]
                        chr_key <- chr
                        if (is.na(chr_key) || chr_key == "" || !(chr_key %in% names(chr_color_map))) {
                          chr_color <- default_chr_color
                        } else {
                          chr_color <- chr_color_map[[chr_key]]
                        }
                        if (is.na(chr_color)) {
                          chr_color <- default_chr_color
                        }
                       
                       # 绘制染色体标记（参考group_list逻辑）
                       circos.rect(
                         xleft = i - 1.0, ybottom = 0.1,
                         xright = i + 0.0, ytop = 0.9,
                         sector.index = current_sector,
                         col = chr_color,
                         border = NA, lwd = 0.1
                       )
                     }
                     cat("染色体轨道绘制完成\n")
                   })
      
      # 添加染色体编号标签
      cat("添加染色体编号标签\n")
      
      # 识别染色体边界和中心位置
      chr_changes <- data.frame(
        chr = character(),
        start_idx = numeric(),
        end_idx = numeric(),
        center_idx = numeric(),
        stringsAsFactors = FALSE
      )
      
      current_chr <- chr_info$chr[1]
      start_idx <- 1
      
      for (i in seq_len(nrow(chr_info))) {
        if (i == nrow(chr_info) || chr_info$chr[i + 1] != current_chr) {
          # 染色体结束
          end_idx <- i
          center_idx <- (start_idx + end_idx) / 2
          
          chr_changes <- rbind(chr_changes, data.frame(
            chr = current_chr,
            start_idx = start_idx,
            end_idx = end_idx,
            center_idx = center_idx,
            stringsAsFactors = FALSE
          ))
          
          if (i < nrow(chr_info)) {
            current_chr <- chr_info$chr[i + 1]
            start_idx <- i + 1
          }
        }
      } 
      
      chr_changes$chr <- as.character(chr_changes$chr)
      cat(sprintf("识别到 %d 个染色体\n", nrow(chr_changes)))
      
      # 在染色体轨道上添加文本标签
      track_idx_chr <- get.current.track.index()
      
      for (i in seq_len(nrow(chr_changes))) {
        chr_label <- chr_changes$chr[i]
        if (is.na(chr_label) || chr_label == "") {
          chr_label <- chr_info$chr[chr_changes$start_idx[i]]
        }
        chr_label <- as.character(chr_label)
        center_pos <- chr_changes$center_idx[i] - 0.5  # 调整到中心位置
        
        # 只标注数字染色体和X、Y染色体（按范围调整字号以减少重叠）
        chr_width <- chr_changes$end_idx[i] - chr_changes$start_idx[i] + 1
        label_cex <- ifelse(chr_width >= 3, 0.6, 0.45)

        circos.text(
          x = center_pos,
          y = 0.5,
          labels = chr_label,
          sector.index = sector_name,
          track.index = track_idx_chr,
          facing = "inside",
          niceFacing = TRUE,
          adj = c(0.5, 0.5),
          cex = label_cex,
          col = "white",
          font = 2
        )
      }
      
      cat("染色体编号标签添加完成\n")
      
      # 添加染色体图例
      legend("bottomright",
         legend = c("Chromosome Group 1", "Chromosome Group 2"),
         fill = c(chr_color_group1, chr_color_group2),
         border = NA,
         cex = 0.6,
         pt.cex = 0.8,
         ncol = 1,
         title = "Chromosomes",
         title.adj = 0,
         bty = "n",
         x.intersp = 0.5,
         y.intersp = 0.8)
    } else {
      cat("警告：chr_info 为空或没有数据，跳过染色体轨道\n")
    }
    
    title(paste("Circular Heatmap - BP Associations (", model_name, "-", level_name, ")", sep = ""),
          cex.main = 1.4, font.main = 2, col.main = "#2C3E50")
    
    legend("topright", 
           legend = categorical_labels,
           fill = discrete_colors,
           border = NA,
           cex = 0.5,
           pt.cex = 0.7, 
           ncol = 1,
           title = "Significance Pattern",
           title.adj = 0.5,
           bty = "n",
           x.intersp = 0.5,
           y.intersp = 0.8)
    
    cat("环状热图创建成功\n")
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf("创建环状热图时出错: %s\n", e$message))
    circos.clear()
    return(NULL)
  })
}

# =============================================================================
# 主分析开始 - 循环处理每个模型
# =============================================================================

models <- c("model_1", "model_2", "model_3")
all_results <- list()

for (model_name in models) {
  cat(sprintf("\n=== 处理 %s ===\n", model_name))
  
  # =============================================================================
  # 步骤1：读入数据
  # =============================================================================
  
  # 读取血压专用merged文件
  merged_file <- file.path(base_dir, model_name, "merged_obs_allmr_hyprcoloc.txt")
  if (!file.exists(merged_file)) {
    warning(sprintf("Merged file not found: %s", merged_file))
    next
  }
  merged_data <- read.table(merged_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat(sprintf("读取merged数据完成：%d个蛋白质\n", nrow(merged_data)))
  
  # 读取obs原始结果文件（包含Estimate）
  obs_file <- file.path(base_dir, model_name, "obs_all_results.txt")
  if (!file.exists(obs_file)) {
    warning(sprintf("OBS file not found: %s", obs_file))
    next
  }
  obs_data <- read.table(obs_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat(sprintf("读取OBS数据完成：%d行\n", nrow(obs_data)))
  
  # 读取MR原始结果文件
  mr_file <- file.path(base_dir, "mr", "mr_all_results.txt")
  if (file.exists(mr_file)) {
    mr_data <- read.table(mr_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # 筛选IVW_Egger方法
    mr_data <- mr_data[mr_data$method == "IVW_Egger", ]
    cat(sprintf("读取MR数据完成：%d行\n", nrow(mr_data)))
  } else {
    mr_data <- data.frame()
    cat("MR文件不存在\n")
  }
  
  # 读取GSMR原始结果文件
  gsmr_file <- file.path(base_dir, "gsmr", "gsmr_all_results.txt")
  if (file.exists(gsmr_file)) {
    gsmr_data <- read.table(gsmr_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    cat(sprintf("读取GSMR数据完成：%d行\n", nrow(gsmr_data)))
  } else {
    gsmr_data <- data.frame()
    cat("GSMR文件不存在\n")
  }
  
  # 读取注释信息
  anno_info <- read.csv("E:/results/anno.info_updat.csv", header = TRUE, stringsAsFactors = FALSE)
  
  # 读取基因位置信息（用于染色体排序）
  gene_pos_file <- "D:/OneDrive/工作/1.工作/blood pressure/hg38_gene_list.txt"
  if (file.exists(gene_pos_file)) {
    gene_pos <- read.table(
      gene_pos_file,
      col.names = c("chr", "start", "end", "gene", "ensembl"),
      stringsAsFactors = FALSE,
      #fill = TRUE,  # 允许列数不一致
      na.strings = c("NA", "")  # 处理空值
    )
    # 只需要染色体信息，不需要精确位置
    gene_pos$chr_num <- as.numeric(ifelse(
      gene_pos$chr %in% c("X", "Y", "M", "MT"),
      ifelse(gene_pos$chr == "X", 23,
             ifelse(gene_pos$chr == "Y", 24, 25)),
      gene_pos$chr
    ))
    cat(sprintf("读取了%d个基因的染色体位置信息\n", nrow(gene_pos)))
  } else {
    gene_pos <- data.frame()
    cat(sprintf("警告：未找到基因位置文件: %s\n", gene_pos_file))
  }
  
  # 读取BP三个层级的蛋白质列表
  protein_list_file <- file.path(base_dir, model_name, "upset_venn", "BP_three_levels_protein_lists.xlsx")
  if (!file.exists(protein_list_file)) {
    warning(sprintf("Protein list file not found: %s", protein_list_file))
    next
  }
  
  sheet_names <- getSheetNames(protein_list_file)
  cat(sprintf("Available sheets: %s\n", paste(sheet_names, collapse = ", ")))
  
  # =============================================================================
  # 处理两个层级的热图
  # =============================================================================
  
  level_configs <- list(
    list(sheet_name = "2_Regression_MR", level_name = "Level2_Regression_MR", description = "Regression+MR"),
    list(sheet_name = "3_Regression_MR_Coloc", level_name = "Level3_Regression_MR_Coloc", description = "Regression+MR+Colocalization")
  )
  
  for (level_config in level_configs) {
    sheet_name <- level_config$sheet_name
    level_name <- level_config$level_name
    level_description <- level_config$description
    
    cat(sprintf("\n=== 处理 %s: %s ===\n", level_name, level_description))
    
    if (!(sheet_name %in% sheet_names)) {
      warning(sprintf("Sheet %s not found in protein list file", sheet_name))
      next
    }
    
    # 读取该层级的蛋白质列表
    protein_list <- openxlsx::read.xlsx(protein_list_file, sheet = sheet_name)
    if (nrow(protein_list) == 0) {
      warning(sprintf("No proteins found in sheet %s", sheet_name))
      next
    }
    
    cat(sprintf("读取到 %d 个蛋白质\n", nrow(protein_list)))
    
    # 确保有protein_id和Entrez.Gene.Name列
    if (!all(c("protein_id", "Entrez.Gene.Name") %in% colnames(protein_list))) {
      warning(sprintf("Required columns not found in sheet %s", sheet_name))
      next
    }
    
    # 只保留有Entrez.Gene.Name的蛋白质
    protein_list <- protein_list[!is.na(protein_list$Entrez.Gene.Name) & protein_list$Entrez.Gene.Name != "", ]
    multi_group_proteins <- protein_list$protein_id
    
    cat(sprintf("有效蛋白质数：%d\n", length(multi_group_proteins)))
    
    # =============================================================================
    # 步骤2：创建热图数据
    # =============================================================================
    
    cat("\n=== 创建热图数据 ===\n")
    
    heatmap_data_list <- list()
    
    # 为每个trait处理数据
    for (i in seq_len(nrow(trait_mapping))) {
      obs_trait <- trait_mapping$obs_trait[i]
      mr_trait <- trait_mapping$mr_trait[i]
      group <- trait_mapping$group[i]
      
      # 为每个筛选的蛋白质创建记录
      for (protein in multi_group_proteins) {
        # 获取protein_id对应的Entrez.Gene.Name
        protein_info <- protein_list[protein_list$protein_id == protein, ][1, ]
        entrez_gene_name <- protein_info$Entrez.Gene.Name
        
        # 获取显著性信息
        sig_info <- get_significance_info(merged_data, protein, obs_trait, mr_trait)
        
        # 从obs_data获取效应估计值和方向
        obs_estimate <- NA
        obs_direction <- NA
        protein_obs_data <- obs_data[obs_data$traits == obs_trait & obs_data$seqName == protein, ]
        if (nrow(protein_obs_data) > 0) {
          obs_estimate <- protein_obs_data$Estimate[1]
          if (!is.na(obs_estimate)) {
            obs_direction <- ifelse(obs_estimate > 0, "positive", "negative")
          }
        }
        
        # 从MR/GSMR数据获取效应估计值和方向
        mr_estimate <- NA
        mr_direction <- NA
        
        # 先尝试从MR数据获取
        protein_mr_data <- data.frame()
        if (nrow(mr_data) > 0 && !is.na(mr_trait)) {
          protein_mr_data <- mr_data[mr_data$traits == mr_trait & mr_data$protein_id == protein, ]
        }
        
        # 如果MR数据为空，尝试从GSMR数据获取
        if (nrow(protein_mr_data) == 0 && nrow(gsmr_data) > 0 && !is.na(mr_trait)) {
          protein_mr_data <- gsmr_data[gsmr_data$traits == mr_trait & gsmr_data$protein_id == protein, ]
          if (nrow(protein_mr_data) > 0) {
            # GSMR使用bxy作为效应估计
            if ("bxy" %in% colnames(protein_mr_data)) {
              colnames(protein_mr_data)[which(colnames(protein_mr_data) == "bxy")] <- "b"
            }
          }
        }
        
        if (nrow(protein_mr_data) > 0) {
          mr_estimate <- protein_mr_data$b[1]
          if (!is.na(mr_estimate)) {
            mr_direction <- ifelse(mr_estimate > 0, "positive", "negative")
          }
        }
        
        # 创建记录
        protein_record <- data.frame(
          protein_id = protein,
          Entrez.Gene.Name = entrez_gene_name,
          trait = obs_trait,
          group = group,
          obs_estimate = obs_estimate,
          mr_estimate = mr_estimate,
          obs_direction = obs_direction,
          mr_direction = mr_direction,
          obs_significant = sig_info$obs_significant,
          union_significant = sig_info$union_significant,
          stringsAsFactors = FALSE
        )
        
        heatmap_data_list[[paste0(protein, "_", obs_trait)]] <- protein_record
      }
    }
    
    # 合并所有数据
    heatmap_data <- do.call(rbind, heatmap_data_list)
    
    if (is.null(heatmap_data) || nrow(heatmap_data) == 0) {
      cat("没有生成热图数据\n")
      next
    }
    
    # 计算一致性和热图值
    heatmap_data <- heatmap_data %>%
      mutate(
        consistency = case_when(
          obs_significant & union_significant & !is.na(obs_direction) & !is.na(mr_direction) & obs_direction == mr_direction ~ "both_significant_consistent",
          obs_significant & union_significant & !is.na(obs_direction) & !is.na(mr_direction) & obs_direction != mr_direction ~ "both_significant_inconsistent",
          obs_significant & !union_significant ~ "only_obs_significant",
          !obs_significant & union_significant ~ "only_union_significant",
          !obs_significant & !union_significant ~ "neither_significant",
          TRUE ~ "unknown"
        ),
        heatmap_value = case_when(
          consistency == "both_significant_consistent" & obs_direction == "positive" ~ 1,
          consistency == "both_significant_consistent" & obs_direction == "negative" ~ -1,
          consistency == "both_significant_inconsistent" & obs_direction == "positive" ~ 0.6,
          consistency == "both_significant_inconsistent" & obs_direction == "negative" ~ -0.6,
          consistency == "only_obs_significant" & !is.na(obs_direction) & obs_direction == "positive" ~ 0.3,
          consistency == "only_obs_significant" & !is.na(obs_direction) & obs_direction == "negative" ~ -0.3,
          consistency == "only_union_significant" & !is.na(mr_direction) & mr_direction == "positive" ~ 0.2,
          consistency == "only_union_significant" & !is.na(mr_direction) & mr_direction == "negative" ~ -0.2,
          TRUE ~ 0
        )
      )
    
    cat(sprintf("热图数据创建完成：%d个蛋白质-trait组合\n", nrow(heatmap_data)))
    
    # 调试：检查heatmap_value分布
    cat("\n=== 调试：Heatmap值分布 ===\n")
    cat(sprintf("非零值数量: %d\n", sum(heatmap_data$heatmap_value != 0, na.rm = TRUE)))
    cat(sprintf("NA值数量: %d\n", sum(is.na(heatmap_data$heatmap_value))))
    cat("Heatmap值统计:\n")
    print(table(heatmap_data$heatmap_value))
    cat("\n一致性统计:\n")
    print(table(heatmap_data$consistency))
    cat("\n显著性统计:\n")
    cat(sprintf("obs显著: %d\n", sum(heatmap_data$obs_significant, na.rm = TRUE)))
    cat(sprintf("union显著: %d\n", sum(heatmap_data$union_significant, na.rm = TRUE)))
    cat("\n前10行数据:\n")
    print(head(heatmap_data[, c("Entrez.Gene.Name", "trait", "obs_significant", "union_significant", "obs_direction", "mr_direction", "consistency", "heatmap_value")], 10))
    cat("\n")
    
    # =============================================================================
    # 步骤3：创建输出文件夹并绘图
    # =============================================================================
     
    heatmap_path <- file.path(base_dir, model_name, "heatmap", level_name)
    if (!dir.exists(heatmap_path)) {
      dir.create(heatmap_path, recursive = TRUE)
    }
    
    # 保存数据
    write.csv(heatmap_data, file.path(heatmap_path, "heatmap_data.csv"), row.names = FALSE)
    
    # 创建宽格式数据用于环状热图
    circ_heatmap_mat <- heatmap_data %>%
      dplyr::select(protein_id, Entrez.Gene.Name, trait, heatmap_value) %>%
      pivot_wider(names_from = trait, values_from = heatmap_value, values_fill = 0) %>%
      distinct(Entrez.Gene.Name, .keep_all = TRUE)
    
    write.csv(circ_heatmap_mat, file.path(heatmap_path, "circular_heatmap_data.csv"), row.names = FALSE)
    
    # 1. 绘制横向热图
    cat("\n=== 绘制横向热图 ===\n")
    bp_heatmap <- create_heatmap_plot(heatmap_data, model_name, level_name)
    
    if (!is.null(bp_heatmap)) {
      ggsave(file.path(heatmap_path, "BP_heatmap.pdf"), bp_heatmap, width = 18, height = 6)
      ggsave(file.path(heatmap_path, "BP_heatmap.png"), bp_heatmap, width = 18, height = 6, dpi = 300)
      cat(sprintf("横向热图已保存到 %s\n", heatmap_path))
    }
    
    # 2. 绘制环状热图
    cat("\n=== 绘制环状热图 ===\n")
    pdf(file.path(heatmap_path, "BP_circular_heatmap.pdf"), width = 12, height = 12)
    create_circular_heatmap(circ_heatmap_mat, model_name, level_name, gene_pos)
    dev.off()
    
    png(file.path(heatmap_path, "BP_circular_heatmap.png"), width = 12, height = 12, units = "in", res = 300)
    create_circular_heatmap(circ_heatmap_mat, model_name, level_name, gene_pos)
    dev.off()
    
    cat(sprintf("环状热图已保存到 %s\n", heatmap_path))
    
    # 输出统计信息
    cat("\n=== 统计信息 ===\n")
    cat(sprintf("蛋白质-表型组合总数: %d\n", nrow(heatmap_data)))
    cat(sprintf("唯一蛋白质数: %d\n", length(unique(heatmap_data$Entrez.Gene.Name))))
    cat(sprintf("唯一表型数: %d\n", length(unique(heatmap_data$trait))))
    cat(sprintf("双重显著且一致: %d\n", sum(heatmap_data$consistency == "both_significant_consistent", na.rm = TRUE)))
    cat(sprintf("双重显著但不一致: %d\n", sum(heatmap_data$consistency == "both_significant_inconsistent", na.rm = TRUE)))
    cat(sprintf("仅obs显著: %d\n", sum(heatmap_data$consistency == "only_obs_significant", na.rm = TRUE)))
    cat(sprintf("仅union显著: %d\n", sum(heatmap_data$consistency == "only_union_significant", na.rm = TRUE)))
    
    # 保存结果
    all_results[[paste0(model_name, "_", level_name)]] <- heatmap_data
  }
}

cat("\n========== 血压热图分析完成 ==========\n")
cat("分析逻辑：1.读入数据 -> 2.创建热图数据 -> 3.绘图\n")
cat("热图分析已完成！\n")
