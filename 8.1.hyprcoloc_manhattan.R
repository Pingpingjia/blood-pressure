# Title: HyPrColoc Manhattan Plot with MR Validation
# Purpose: Create Manhattan plots for colocalization results with MR validation highlights
# ===========================================================
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(readr)
library(data.table)
if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)

# =================== 1. 定义颜色方案 ===================
original_colors <- c(
  "#918BD9",  # 颜色1
  "#B3B4DF",  # 颜色2
  "#CFCCE3",  # 颜色3
  "#F6DFD6",  # 颜色4
  "#F8B2A2",  # 颜色5
  "#F1B6BB",  # 颜色6
  "#E9687A"   # 颜色7
)

# 生成22条染色体颜色向量
chr_colors <- rep(original_colors, length.out = 22)
# =================== 2. 读取数据 ===================
# 共定位数据
coloc_file <- "D:/OneDrive/工作/1.工作/blood pressure/hyprcoloc/hyprcoloc_max_posterior_prob.txt"
coloc_data_raw <- fread(coloc_file, header = TRUE, stringsAsFactors = FALSE)
coloc_data <- coloc_data_raw[coloc_data_raw$protein!="",]
coloc_data <- coloc_data%>%rename(region_chr = chr, region_gene = gene)

# 定义模型列表
models <- c("model_1", "model_2", "model_3")

# 为每个模型处理数据和生成图表
for (model in models) {
  cat(sprintf("\n=== 处理 %s ===\n", model))
  
  # MR验证数据
  mr_file <- file.path("D:/OneDrive/工作/1.工作/blood pressure", model, "merged_obs_allmr_hyprcoloc.txt")
  mr_data <- fread(mr_file, header = TRUE, stringsAsFactors = FALSE)

  # 基因注释数据（需要获取染色体位置信息）
  gene_file <- "D:/OneDrive/工作/1.工作/blood pressure/hg38_gene_list.txt"
  gene_dt <- fread(gene_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(gene_dt) <- c("chr", "start_pos", "end_pos", "gene", "ensembl_id")

 # =================== 3. 数据预处理 ===================
  # 处理MR数据，创建显著蛋白列表
  mr_significant <- mr_data %>%
    mutate(
      sbp_sig = ifelse(sbp.obs == 1 & ckb_sbp.union == 1, TRUE, FALSE),
      dbp_sig = ifelse(dbp.obs == 1 & ckb_dbp.union == 1, TRUE, FALSE),
      pp_sig = ifelse(pp.obs == 1 & ckb_pp.union == 1, TRUE, FALSE),
      map_sig = ifelse(map.obs == 1 & ckb_map.union == 1 , TRUE, FALSE),
      hypertension_sig = ifelse(hypertension.obs == 1 & ukb_hypertension.union == 1, TRUE, FALSE)
    )
  # 创建表型映射
  trait_mapping <- c(
    "ckb_sbp" = "sbp_sig",
    "ckb_dbp" = "dbp_sig", 
    "ckb_pp" = "pp_sig",
    "ckb_map" = "map_sig",
    "ukb_hypertension" = "hypertension_sig"
  )
  # 处理共定位数据
  coloc_processed <- coloc_data %>%
    filter(!is.na(protein_seq) & !is.na(trait) & !is.na(max_posterior_prob)) %>%
    # 添加基因注释
    mutate(
      gene_name = ifelse(!is.na(protein), protein, "Unknown")
    ) %>%
    # 关联真实的基因位置数据
    left_join(
      gene_dt %>% 
        select(gene, chr, start_pos, end_pos) %>%
        # 处理重复的基因，只保留第一个匹配
        distinct(gene, .keep_all = TRUE) %>%
        rename(gene_chr = chr, gene_start = start_pos, gene_end = end_pos),
      by = c("gene_name" = "gene")
    ) %>%
    mutate(
      # 处理染色体数据（可能是字符型）
      gene_chr_clean = case_when(
        is.na(gene_chr) ~ "1",
        gene_chr == "X" ~ "23",
        gene_chr == "Y" ~ "24",
        gene_chr == "MT" ~ "25",
        TRUE ~ as.character(gene_chr)
      ),
      # 如果没有找到基因位置，使用默认值
      chr = as.integer(ifelse(is.na(gene_chr), 25, 
                             ifelse(gene_chr_clean %in% c("23", "24", "25"), gene_chr_clean, gene_chr_clean)))
    ) %>%
    # 清理临时列（但保留 gene_name）
    select(-c(gene_chr, gene_start, gene_end, gene_chr_clean)) %>%
    # 标记MR验证的蛋白
    left_join(
      mr_significant %>% select(protein_id, Entrez.Gene.Name, sbp_sig, dbp_sig, pp_sig, map_sig, hypertension_sig),
      by = c("protein_seq" = "protein_id")
    ) %>%
    mutate(
      mr_validated = case_when(
        trait == "ckb_sbp" & sbp_sig == TRUE ~ TRUE,
        trait == "ckb_dbp" & dbp_sig == TRUE ~ TRUE,
        trait == "ckb_pp" & pp_sig == TRUE ~ TRUE,
        trait == "ckb_map" & map_sig == TRUE ~ TRUE,
        trait == "ukb_hypertension" & hypertension_sig == TRUE ~ TRUE,
        TRUE ~ FALSE
      ),
      # 重命名后验概率列
      posterior_prob = max_posterior_prob
    ) %>%
    ungroup()
test <- coloc_processed[coloc_processed$mr_validated == TRUE,]
# =================== 4. 绘图函数 ===================
plot_coloc_manhattan <- function(data, main_title, highlight_mr = TRUE) {
  
  # 根据是否突出显示MR验证点调整数据分层
  if (highlight_mr) {
    p <- ggplot(data, aes(x = chr, y = posterior_prob)) +
      # 1. 普通点（非MR验证）
      geom_point(
        data = ~filter(., !mr_validated),
        aes(color = as.factor(chr)),
        alpha = 0.8,
        size = 2.5,
        stroke = 0.2
      ) +
      # 2. MR验证点（突出显示）
      geom_point(
        data = ~filter(., mr_validated),
        aes(fill = as.factor(chr)),
        color = "black",
        size = 2.0,
        alpha = 1.0,
        stroke = 1.0,
        shape = 21
      ) +
      # 3. MR验证点标签
      ggrepel::geom_text_repel(
        data = ~filter(., mr_validated),
        aes(label = gene_name),
        color = "black",
        size = 3.5,
        fontface = "bold",
        box.padding = 0.3,
        point.padding = 0.2,
        segment.alpha = 0.6,
        segment.size = 0.3,
        max.overlaps = 50,
        min.segment.length = 0.2,
        force = 2,
        force_pull = 1,
        direction = "both"
      )
  } else {
    p <- ggplot(data, aes(x = chr, y = posterior_prob)) +
      # 所有点使用相同样式
      geom_point(
        aes(color = as.factor(chr)),
        alpha = 0.8,
        size = 2.5,
        stroke = 0.2
      )
  }
  
  p <- p +
    # 4. 阈值线 - 0.5
    geom_hline(
      yintercept = 0.5,
      color = "red",
      linetype = "dashed",
      size = 0.5,
      alpha = 0.8
    ) +
    # 5. 阈值线 - 0.8
    geom_hline(
      yintercept = 0.8,
      color = "darkred",
      linetype = "dashed",
      size = 0.5,
      alpha = 0.8
    ) +
    # 6. 阈值线标注
    annotate(
      "text",
      x = 22.2, y = 0.53,
      label = "",
      color = "red",
      size = 3.0,
      hjust = 1
    ) +
    annotate(
      "text",
      x = 22.2, y = 0.83,
      label = "",
      color = "darkred",
      size = 3.0,
      hjust = 1
    ) +
    
    # 7. 颜色设置
    scale_color_manual(
      values = setNames(chr_colors, 1:22),
      guide = "none"
    ) +
    scale_fill_manual(
      values = setNames(chr_colors, 1:22),
      guide = "none"
    ) +
    
    # 8. 坐标轴设置
    scale_x_continuous(
      breaks = 1:22,
      labels = 1:22,
      limits = c(0.2, 22.8),
      name = "Chromosome",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, max(data$posterior_prob, na.rm = TRUE) * 1.1),
      breaks = seq(0, 1, 0.2),  # 添加更多Y轴刻度
      name = "Posterior Probability",
      expand = c(0, 0)
    ) +
    
    # 9. 主题设置
    theme_bw() +
    theme(
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.line.y = element_line(color = "black", size = 0.5),
      axis.title.x = element_text(
        color = "black", size = 12, face = "bold", margin = margin(t = 10)
      ),
      axis.title.y = element_text(
        color = "black", size = 12, face = "bold", margin = margin(r = 10)
      ),
      axis.text.x = element_text(color = "black", size = 10, hjust = 0.5),
      axis.text.y = element_text(color = "black", size = 10),
      plot.title = element_text(
        color = "black", size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)
      ),
      panel.grid.major = element_line(color = "#ECF0F1", size = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = 10, r = 40, b = 10, l = 40, unit = "mm")
    ) +
    labs(title = main_title)
  
  return(p)
}

# 综合曼哈顿图绘制函数（显示所有表型结果）
plot_combined_traits_manhattan <- function(data, main_title, highlight_mr = TRUE) {
  
  # 根据是否突出显示MR验证点调整数据分层
  if (highlight_mr) {
    p <- ggplot(data, aes(x = chr, y = posterior_prob)) +
      # 1. 普通点（非MR验证）
      geom_point(
        data = ~filter(., !mr_validated),
        aes(color = as.factor(chr)),
        alpha = 0.8,
        size = 2.5,
        stroke = 0.2
      ) +
      # 2. MR验证点（突出显示）
      geom_point(
        data = ~filter(., mr_validated),
        aes(fill = as.factor(chr)),
        color = "black",
        size = 2.0,
        alpha = 1.0,
        stroke = 1.0,
        shape = 21
      ) +
      # 3. MR验证点标签（包含表型信息）
      ggrepel::geom_text_repel(
        data = ~filter(., mr_validated),
        aes(label = protein_with_traits),
        color = "black",
        size = 3.2,
        fontface = "bold",
        box.padding = 0.4,
        point.padding = 0.3,
        segment.alpha = 0.6,
        segment.size = 0.3,
        max.overlaps = 50,
        min.segment.length = 0.2,
        force = 2,
        force_pull = 1,
        direction = "both",
        hjust = 0.5,  # 水平居中
        vjust = 0.5   # 垂直居中
      )
  } else {
    p <- ggplot(data, aes(x = chr, y = posterior_prob)) +
      # 所有点使用染色体颜色
      geom_point(
        aes(color = as.factor(chr)),
        alpha = 0.8,
        size = 2.5,
        stroke = 0.2
      )
  }
  
  p <- p +
    # 4. 阈值线 - 0.5
    geom_hline(
      yintercept = 0.5,
      color = "red",
      linetype = "dashed",
      size = 0.5,
      alpha = 0.8
    ) +
    # 5. 阈值线 - 0.8
    geom_hline(
      yintercept = 0.8,
      color = "darkred",
      linetype = "dashed",
      size = 0.5,
      alpha = 0.8
    ) +
    # 6. 阈值线标注
    annotate(
      "text",
      x = 22.2, y = 0.53,
      label = "",
      color = "red",
      size = 3.0,
      hjust = 1
    ) +
    annotate(
      "text",
      x = 22.2, y = 0.83,
      label = "",
      color = "darkred",
      size = 3.0,
      hjust = 1
    ) +
    
    # 7. 颜色设置
    scale_color_manual(
      values = setNames(chr_colors, 1:22),
      guide = "none"
    ) +
    scale_fill_manual(
      values = setNames(chr_colors, 1:22),
      guide = "none"
    ) +
    
    # 8. 坐标轴设置
    scale_x_continuous(
      breaks = 1:22,
      labels = 1:22,
      limits = c(0.2, 22.8),
      name = "Chromosome",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, max(data$posterior_prob, na.rm = TRUE) * 1.1),
      breaks = seq(0, 1, 0.2),  # 添加更多Y轴刻度
      name = "Posterior Probability",
      expand = c(0, 0)
    ) +
    
    # 9. 主题设置
    theme_bw() +
    theme(
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.line.y = element_line(color = "black", size = 0.5),
      axis.title.x = element_text(
        color = "black", size = 12, face = "bold", margin = margin(t = 10)
      ),
      axis.title.y = element_text(
        color = "black", size = 12, face = "bold", margin = margin(r = 10)
      ),
      axis.text.x = element_text(color = "black", size = 10, hjust = 0.5),
      axis.text.y = element_text(color = "black", size = 10),
      plot.title = element_text(
        color = "black", size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)
      ),
      panel.grid.major = element_line(color = "#ECF0F1", size = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = 10, r = 40, b = 10, l = 40, unit = "mm")
    ) +
    labs(title = main_title)
  
  return(p)
}

  # =================== 5. 生成图表 ===================
  output_dir <- file.path("D:/OneDrive/工作/1.工作/blood pressure", model, "hyprcoloc_manhattan")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 5.1 汇总曼哈顿图（所有表型合并）
  summary_data <- coloc_processed %>%
    # 对于同一个蛋白，如果在任何表型中都有MR验证，就标记为验证
    group_by(protein_seq) %>%
    mutate(
      protein_mr_validated = any(mr_validated, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    # 每个蛋白只保留一行（选择最高后验概率的）
    group_by(protein_seq) %>%
    slice_max(order_by = posterior_prob, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(mr_validated = protein_mr_validated)

  # 5.1.2 创建综合表型数据（去重并汇总所有表型信息）
  combined_traits_data <- coloc_processed %>%
    # 为每个蛋白汇总所有表型信息
    group_by(protein_seq, gene_name, chr) %>%
    summarise(
      # 选择最高的后验概率作为代表
      max_posterior_prob = max(posterior_prob, na.rm = TRUE),
      # 汇总所有显著表型（后验概率>0.25）
      significant_traits = {
        sig_traits <- trait[posterior_prob > 0.25]
        if (length(sig_traits) > 0) {
          trait_labels <- case_when(
            sig_traits == "ckb_sbp" ~ "SBP",
            sig_traits == "ckb_dbp" ~ "DBP",
            sig_traits == "ckb_pp" ~ "PP",
            sig_traits == "ckb_map" ~ "MAP",
            sig_traits == "ukb_hypertension" ~ "HT",
            TRUE ~ sig_traits
          )
          paste0("(", paste(unique(sort(trait_labels)), collapse = ","), ")")
        } else {
          ""
        }
      },
      # 检查是否有任何MR验证
      has_mr_validation = any(mr_validated == TRUE, na.rm = TRUE),
      # 检查是否有任何显著结果
      has_significant = any(posterior_prob > 0.25, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # 重命名列以保持一致性
      posterior_prob = max_posterior_prob,
      mr_validated = has_mr_validation,
      # 创建带表型信息的蛋白标签
      protein_with_traits = ifelse(
        significant_traits != "" & mr_validated & has_significant,
        paste0(gene_name, "\n", significant_traits),
        gene_name
      )
    )

  p_summary <- plot_coloc_manhattan(
    data = summary_data,
    main_title = sprintf("HyPrColoc Results Summary - %s (MR Validated Proteins Highlighted)", model),
    highlight_mr = TRUE
  )

  # 5.2 各表型单独的曼哈顿图
  trait_names <- c("ckb_sbp", "ckb_dbp", "ckb_pp", "ckb_map", "ukb_hypertension")
  trait_titles <- c(
    "ckb_sbp" = sprintf("SBP Colocalization - %s (MR Validated)", model),
    "ckb_dbp" = sprintf("DBP Colocalization - %s (MR Validated)", model),
    "ckb_pp" = sprintf("PP Colocalization - %s (MR Validated)", model),
    "ckb_map" = sprintf("MAP Colocalization - %s (MR Validated)", model),
    "ukb_hypertension" = sprintf("Hypertension Colocalization - %s (MR Validated)", model)
  )

  # 创建各表型图表
  trait_plots <- list()
  for (trait_name in trait_names) {
    trait_data <- coloc_processed %>%
      filter(trait == trait_name) %>%
      arrange(desc(posterior_prob))
    
    if (nrow(trait_data) > 0) {
      trait_plots[[trait_name]] <- plot_coloc_manhattan(
        data = trait_data,
        main_title = trait_titles[trait_name],
        highlight_mr = TRUE
      )
    }
  }

  # 5.3 生成综合表型图（所有表型在一个图中）
  p_combined_traits <- plot_combined_traits_manhattan(
    data = combined_traits_data,
    main_title = sprintf("HyPrColoc All Traits Combined - %s (MR Validated with Trait Labels)", model),
    highlight_mr = TRUE
  )

  # =================== 6. 输出图表 ===================
  # 6.1 保存汇总图
  ggsave(
    filename = file.path(output_dir, sprintf("HyPrColoc_Manhattan_Summary_%s.pdf", model)),
    plot = p_summary,
    width = 16,
    height = 10,
    device = "pdf"
  )

  # 6.1.2 保存综合表型图
  ggsave(
    filename = file.path(output_dir, sprintf("HyPrColoc_Manhattan_Combined_Traits_%s.pdf", model)),
    plot = p_combined_traits,
    width = 16,
    height = 10,
    device = "pdf"
  )

  # 6.2 保存各表型图
  for (trait_name in names(trait_plots)) {
    ggsave(
      filename = file.path(output_dir, sprintf("HyPrColoc_Manhattan_%s_%s.pdf", trait_name, model)),
      plot = trait_plots[[trait_name]],
      width = 16,
      height = 8,
      device = "pdf"
    )
  }

  # 6.3 保存组合图（所有表型）
  if (length(trait_plots) > 0) {
    p_combined <- wrap_plots(trait_plots, ncol = 1)
    
    ggsave(
      filename = file.path(output_dir, sprintf("HyPrColoc_Manhattan_All_Traits_%s.pdf", model)),
      plot = p_combined,
      width = 16,
      height = 40,
      device = "pdf"
    )
  }

  # =================== 7. 输出统计信息 ===================
  cat(sprintf("\n=== %s Manhattan Plot Generation Complete ===\n", model))
  cat("总共定位结果数量:", nrow(coloc_processed), "\n")
  cat("MR验证的蛋白数量:", sum(coloc_processed$any_mr_validated, na.rm = TRUE), "\n")
  
  # 各表型统计
  for (trait_name in trait_names) {
    trait_count <- sum(coloc_processed$trait == trait_name, na.rm = TRUE)
    mr_count <- sum(coloc_processed$trait == trait_name & coloc_processed$mr_validated == TRUE, na.rm = TRUE)
    cat(sprintf("%s: 总数=%d, MR验证=%d\n", trait_name, trait_count, mr_count))
  }
  
  cat("\n输出文件保存在:", output_dir, "\n")
  cat(sprintf("- HyPrColoc_Manhattan_Summary_%s.pdf (汇总图)\n", model))
  cat(sprintf("- HyPrColoc_Manhattan_Combined_Traits_%s.pdf (综合表型图)\n", model))
  for (trait_name in trait_names) {
    cat(sprintf("- HyPrColoc_Manhattan_%s_%s.pdf\n", trait_name, model))
  }
  cat(sprintf("- HyPrColoc_Manhattan_All_Traits_%s.pdf (组合图)\n", model))

} # 结束模型循环

cat("\n=== 所有模型处理完成 ===")