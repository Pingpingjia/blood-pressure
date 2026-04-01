library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(UpSetR)
library(ggvenn)
library(ggVennDiagram)
library(RColorBrewer)
library(patchwork)

select <- dplyr::select

# Matplotlib Set2 palette (from seaborn/matplotlib)
set2_colors <- c(
  "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
  "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"
)

# Blood Pressure 专用Venn图和UpSet图
# 设置路径
base_dir <- "D:/OneDrive/工作/1.工作/blood pressure"

# Blood pressure 相关表型（包括pp和map）
obs_traits <- c("sbp", "dbp", "hypertension", "pp", "map")

# 创建obs和mr的对应关系映射（只保留BP相关）
obs_mr_mapping <- data.frame(
  obs = c("sbp", "dbp", "hypertension", "pp", "map"),
  mr = c("ckb_sbp", "ckb_dbp", "ukb_hypertension", "ckb_pp", "ckb_map"),
  stringsAsFactors = FALSE
)

methods_causal <- c("mr", "gsmr", "union")

for (model in c("model_1","model_2","model_3")) {
  cat(paste("\n", rep("=", 70), collapse = ""), "\n")
  cat(paste("Processing Blood Pressure UpSet/Venn plots for", model, "\n"))
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # 创建upset_venn子文件夹
  upset_venn_path <- file.path(base_dir, model, "upset_venn")
  if (!dir.exists(upset_venn_path)) {
    dir.create(upset_venn_path, recursive = TRUE)
  }
  
  # 读取BP相关的merged数据
  merged_file <- file.path(base_dir, model, "merged_obs_allmr_hyprcoloc.txt")
  if(!file.exists(merged_file)) {
    warning(paste("Merged file not found:", merged_file))
    next
  }
  
  merged <- read.table(merged_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  merged[is.na(merged)] <- 0
  n_rows <- nrow(merged)
  

  # ========== 1. BP所有子表型的UpSet图 (obs) ==========
  cat("Generating BP obs traits UpSet plot...\n")
  obs_sets <- paste0(obs_traits, ".obs")
  obs_sets <- obs_sets[obs_sets %in% colnames(merged)]
  
  # 创建名称映射
  trait_name_mapping <- c(
    "sbp.obs" = "SBP",
    "dbp.obs" = "DBP",
    "hypertension.obs" = "Hypertension",
    "pp.obs" = "PP",
    "map.obs" = "MAP"
  )
  
  if (length(obs_sets) > 1) {
    # 过滤掉所有trait都为0的行
    merged_filter_obs <- merged[obs_sets]
    merged_filter_obs <- merged_filter_obs[rowSums(merged_filter_obs, na.rm = TRUE) > 0, ]

    # 按照每个 trait 集合大小排序并移除为0的集合
    set_sizes <- colSums(merged_filter_obs, na.rm = TRUE)
    nonzero_sets <- set_sizes > 0
    if (any(nonzero_sets)) {
      set_sizes <- set_sizes[nonzero_sets]
      merged_filter_obs <- merged_filter_obs[, names(set_sizes), drop = FALSE]
    } else {
      next
    }

    set_order <- order(set_sizes, decreasing = TRUE)
    merged_filter_obs <- merged_filter_obs[, set_order, drop = FALSE]
    ordered_obs_sets <- names(set_sizes)[set_order]
    ordered_trait_names <- trait_name_mapping[ordered_obs_sets]

    # 重命名列名为友好的名称，并保持顺序
    colnames(merged_filter_obs) <- ordered_trait_names
    
    upset_all_obs <- upset(
      merged_filter_obs,
      sets = ordered_trait_names,
      keep.order = TRUE,
      order.by = "freq",
      decreasing = TRUE,
      main.bar.color = "#efa3a3",
      sets.bar.color = "#9BC5E4",
      empty.intersections = "off",
      sets.x.label = "Set Size",
      text.scale = c(2.0, 1.8, 1.8, 1.8, 2.8, 2.5),
      mb.ratio = c(0.6, 0.4),
      number.angles = 0,
      point.size = 3.5,
      shade.color = "gray88",
      shade.alpha = 0.25
    )
    
    png(file.path(upset_venn_path, "1.BP_obs_traits_upset.png"), width = 18, height = 10, units = "in", res = 400)
    print(upset_all_obs)
    dev.off()

    pdf(file.path(upset_venn_path, "1.BP_obs_traits_upset.pdf"), width = 18, height = 10)
    print(upset_all_obs)
    dev.off()
    
    cat("  Saved: 1.BP_obs_traits_upset.png/pdf\n")
  }
  
  
  # ========== 2. BP子表型的OBS Venn图 ==========
  cat("Generating BP OBS sub-phenotypes Venn diagram...\n")
  bp_venn_list <- list()
  trait_labels <- c(
    "sbp" = "SBP",
    "dbp" = "DBP", 
    "hypertension" = "Hypertension",
    "pp" = "PP",
    "map" = "MAP"
  )
  
  for (trait in obs_traits) {
    col <- paste0(trait, ".obs")
    if(col %in% colnames(merged)) {
      label <- trait_labels[[trait]]
      if(is.null(label)) label <- toupper(trait)
      bp_venn_list[[label]] <- merged$protein_id[merged[[col]] == 1]
    }
  }

  bp_venn_list <- bp_venn_list[sapply(bp_venn_list, length) > 0]
  
    if (length(bp_venn_list) >= 2 && length(bp_venn_list) <= 5) {
    venn_plot_obs <- ggVennDiagram(
      bp_venn_list,
      set_color = set2_colors[1:length(bp_venn_list)],
      set_size = 6,
      label_alpha = 0,
      label = "count",
      label_color = "black",
      label_size = 5,
      edge_size = 0.5
    ) + 
      #scale_fill_gradientn(colors = c('#DEEBF7', '#4E79A7', '#4292C6')) +
     # scale_fill_gradientn(colors = c("#FF6B6B", "#FFA726", "#FFD166", "#8AC926", "#EF476F"))+
      scale_fill_gradientn(colors = c('#FFE5D9', '#FFB4A2', '#E76F51'))+
      scale_x_continuous(expand = expansion(mult = 0.15)) +
      scale_y_continuous(expand = expansion(mult = 0.15)) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      ) +
      labs(title = paste0("Blood Pressure OBS Sub-phenotypes Overlap (", model, ")"))
    
    png(file.path(upset_venn_path,"2.BP_obs_subphenotypes_venn.png"), width = 14, height = 12, units = "in", res = 300)
    print(venn_plot_obs)
    dev.off()
    
    pdf(file.path(upset_venn_path,"2.BP_obs_subphenotypes_venn.pdf"), width = 14, height = 12)
    print(venn_plot_obs)
    dev.off()
    
    cat("  Saved: 2.BP_obs_subphenotypes_venn.png/pdf\n")
  }
  
  # ========== 3. BP子表型的OBS+UNION Venn图 ==========
  cat("Generating BP OBS+UNION sub-phenotypes Venn diagram...\n")
  bp_union_venn_list <- list()
  
  for (trait in obs_traits) {
    obs_col <- paste0(trait, ".obs")
    
    # 对于有union数据的表型
    mr_trait <- obs_mr_mapping$mr[obs_mr_mapping$obs == trait]
    if (!is.na(mr_trait) && length(mr_trait) > 0) {
      union_col <- paste0(mr_trait, ".union")
      if (obs_col %in% colnames(merged) && union_col %in% colnames(merged)) {
        label <- trait_labels[[trait]]
        if(is.null(label)) label <- toupper(trait)
        # OBS显著且UNION也显著的蛋白
        bp_union_venn_list[[label]] <- merged$protein_id[merged[[obs_col]] == 1 & merged[[union_col]] == 1]
      }
    } else {
      # 没有MR数据的表型（如hypertension, pp, map），只用OBS
      if (obs_col %in% colnames(merged)) {
        label <- trait_labels[[trait]]
        if(is.null(label)) label <- toupper(trait)
        bp_union_venn_list[[label]] <- merged$protein_id[merged[[obs_col]] == 1]
      }
    }
  }
  bp_union_venn_list <- bp_union_venn_list[sapply(bp_union_venn_list, length) > 0]
  
  if (length(bp_union_venn_list) >= 2) {
    venn_plot_union <- ggVennDiagram(
      bp_union_venn_list,
      set_color = set2_colors[1:length(bp_union_venn_list)],
      set_size = 6,
      label_alpha = 0,
      label = "count",
      label_color = "black",
      label_size = 5,
      edge_size = 0.5
    ) + 
      scale_fill_gradientn(colors = c('#FFE5D9', '#FFB4A2', '#E76F51'))+
      scale_x_continuous(expand = expansion(mult = 0.15)) +
      scale_y_continuous(expand = expansion(mult = 0.15)) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      ) +
      labs(title = paste0("Blood Pressure OBS+UNION Sub-phenotypes Overlap (", model, ")"))
    
    png(file.path(upset_venn_path,"3.BP_obs_union_subphenotypes_venn.png"), width = 14, height = 12, units = "in", res = 300)
    print(venn_plot_union)
    dev.off()
    
    pdf(file.path(upset_venn_path,"3.BP_obs_union_subphenotypes_venn.pdf"), width = 14, height = 12)
    print(venn_plot_union)
    dev.off()
    
    cat("  Saved: 3.BP_obs_union_subphenotypes_venn.png/pdf\n")
  }
  
  
  # ========== 4. OBS和UNION两类堆叠条形图 ==========
  cat("Generating BP OBS vs UNION stacked bar plot...\n")
  
  # 创建OBS和UNION两类数据框
  obs_union_counts <- data.frame(
    trait=character(), 
    obs_only=numeric(),
    obs_union=numeric(),
    stringsAsFactors=FALSE
  )
  
  for (trait in obs_traits) {
    obs_col <- paste0(trait, ".obs")
    if (!(obs_col %in% colnames(merged))) next

    obs_sig <- merged[[obs_col]] == 1
    mr_trait <- obs_mr_mapping$mr[obs_mr_mapping$obs == trait]

    # 检查是否有union列
    union_col <- paste0(mr_trait, ".union")
    union_exists <- !is.na(mr_trait) && union_col %in% colnames(merged)
    
    if (union_exists) {
      # 有union数据的表型
      union_sig <- merged[[union_col]] == 1
      union_sig[is.na(union_sig)] <- FALSE
      
      # OBS显著但UNION不显著
      obs_only_count <- sum(obs_sig & !union_sig)
      # OBS和UNION都显著
      obs_union_count <- sum(obs_sig & union_sig)
      
      # 特殊处理hypertension的显示名称
      trait_display <- ifelse(trait == "hypertension", "Hypertension", toupper(trait))
      
      obs_union_counts <- rbind(obs_union_counts, data.frame(
        trait = trait_display,
        obs_only = obs_only_count,
        obs_union = obs_union_count,
        stringsAsFactors = FALSE
      ))
    } else {
      # 没有union数据的表型（如hypertension），全部算作obs_only
      obs_count <- sum(obs_sig)
      trait_display <- ifelse(trait == "hypertension", "Hypertension", toupper(trait))
      
      obs_union_counts <- rbind(obs_union_counts, data.frame(
        trait = trait_display,
        obs_only = obs_count,
        obs_union = 0,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if(nrow(obs_union_counts) > 0) {
    # 按总数排序
    obs_union_counts$total <- obs_union_counts$obs_only + obs_union_counts$obs_union
    obs_union_counts <- obs_union_counts[order(obs_union_counts$total, decreasing = TRUE), ]
    trait_order_obs_union <- obs_union_counts$trait
    obs_union_counts$trait <- factor(obs_union_counts$trait, levels = trait_order_obs_union)
    
    # 转换为长格式
    obs_union_long <- reshape2::melt(
      obs_union_counts[,c("trait", "obs_only", "obs_union")], 
      id.vars="trait", variable.name="type", value.name="count"
    )
    
    # 设置顺序：从底到顶为 obs_only -> obs_union
    obs_union_long$type <- factor(obs_union_long$type, 
                                   levels=c("obs_union", "obs_only"))
    obs_union_long$trait <- factor(obs_union_long$trait, levels = trait_order_obs_union)
    
    # 绘制OBS vs UNION堆叠条形图
    plot_obs_union <- ggplot(obs_union_long, aes(x=trait, y=count, fill=type)) +
      geom_bar(stat="identity", position="stack", width=0.7) +
      scale_fill_manual(
        values=c(
          "obs_union"="#9BC5E4",
          "obs_only"="#efa3a3"
        ),
        labels=c(
          "OBS+UNION",
          "OBS only"
        )
      ) +
      geom_text(
        data = subset(obs_union_long, count > 0),
        aes(label = count),
        color = "black", size = 5.2, fontface = "bold",
        position = position_stack(vjust = 0.5)
      ) +
      # 在右侧添加总数标签
      geom_text(
        data = obs_union_counts %>% filter(total > 0),
        aes(x = trait, y = total + max(total)*0.03, label = paste0("n=", total)),
        inherit.aes = FALSE,
        hjust = 0, size = 5, fontface = "bold", color = "black"
      ) +
      labs(
        x="Blood Pressure Trait",
        y="Number of Significant Proteins",
        fill="Significant in"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x=element_text(angle=0, hjust=0.5, size=14, color="black"),
        axis.text.y=element_text(size=14, color="black"),
        axis.title=element_text(size=14),
        plot.title=element_text(face="bold", hjust=0.5, size=16),
        legend.position="top",
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank()
      ) +
      coord_flip()
  }
  # ========== 5. OBS+UNION内部MR/GSMR组成比例图 ==========
  cat("Generating BP OBS+UNION composition proportion plot...\n")
  
  # 创建OBS+UNION内部组成数据框
  union_composition <- data.frame(
    trait=character(), 
    mr_only=numeric(),
    gsmr_only=numeric(),
    both=numeric(),
    stringsAsFactors=FALSE
  )
  
  for (trait in obs_traits) {
    obs_col <- paste0(trait, ".obs")
    if (!(obs_col %in% colnames(merged))) next

    obs_sig <- merged[[obs_col]] == 1
    mr_trait <- obs_mr_mapping$mr[obs_mr_mapping$obs == trait]

    # 检查是否有union列
    union_col <- paste0(mr_trait, ".union")
    union_exists <- !is.na(mr_trait) && union_col %in% colnames(merged)
    
    if (union_exists) {
      # 有union数据的表型
      union_sig <- merged[[union_col]] == 1
      union_sig[is.na(union_sig)] <- FALSE
      
      # OBS和UNION都显著的
      obs_union_sig <- obs_sig & union_sig
      
      # 获取MR和GSMR列
      mr_col <- paste0(mr_trait, ".mr")
      gsmr_col <- paste0(mr_trait, ".gsmr")
      mr_exists <- mr_col %in% colnames(merged)
      gsmr_exists <- gsmr_col %in% colnames(merged)
      
      mr_sig <- if (mr_exists) {
        tmp <- merged[[mr_col]] == 1
        tmp[is.na(tmp)] <- FALSE
        tmp
      } else {
        rep(FALSE, n_rows)
      }
      
      gsmr_sig <- if (gsmr_exists) {
        tmp <- merged[[gsmr_col]] == 1
        tmp[is.na(tmp)] <- FALSE
        tmp
      } else {
        rep(FALSE, n_rows)
      }
      
      # 只统计OBS+UNION内部的组成
      both_count <- sum(obs_union_sig & mr_sig & gsmr_sig)
      mr_only_count <- sum(obs_union_sig & mr_sig & !gsmr_sig)
      gsmr_only_count <- sum(obs_union_sig & gsmr_sig & !mr_sig)
      
      # 只有当OBS+UNION有数据时才添加
      if ((both_count + mr_only_count + gsmr_only_count) > 0) {
        trait_display <- ifelse(trait == "hypertension", "Hypertension", toupper(trait))
        
        union_composition <- rbind(union_composition, data.frame(
          trait = trait_display,
          mr_only = mr_only_count,
          gsmr_only = gsmr_only_count,
          both = both_count,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  if(nrow(union_composition) > 0) {
    # 按总数排序
    union_composition$total <- union_composition$mr_only + union_composition$gsmr_only + union_composition$both
    union_composition <- union_composition[order(union_composition$total, decreasing = TRUE), ]
    trait_order_comp <- union_composition$trait
    union_composition$trait <- factor(union_composition$trait, levels = trait_order_comp)
    
    # 转换为长格式
    union_comp_long <- reshape2::melt(
      union_composition[,c("trait", "gsmr_only", "mr_only", "both")], 
      id.vars="trait", variable.name="type", value.name="count"
    )
    
    # 设置顺序：从底到顶为 gsmr_only -> mr_only -> both
    union_comp_long$type <- factor(union_comp_long$type, 
                                    levels=c("both", "mr_only", "gsmr_only"))
    union_comp_long$trait <- factor(union_comp_long$trait, levels = trait_order_comp)
    
    # 绘制比例堆叠条形图
    plot_composition <- ggplot(union_comp_long, aes(x=trait, y=count, fill=type)) +
      geom_bar(stat="identity", position="fill", width=0.7) +
      scale_fill_manual(
        values=c(
          "both"="#4A8BC5",        # Both MR and GSMR 
          "mr_only"="#9BC5E4",     # MR only
          "gsmr_only"="#D0E6F7"    # GSMR only
        ),
        labels=c(
          "Both MR and GSMR",
          "MR only",
          "GSMR only"
        )
      ) +
      geom_text(
        data = subset(union_comp_long, count > 0),
        aes(label = count),
        color = "black", size = 5.2, fontface = "bold",
        position = position_fill(vjust = 0.5)
      ) +
      # 在右侧添加总数标签
      geom_text(
        data = union_composition %>% filter(total > 0),
        aes(x = trait, y = 1.02, label = paste0("n=", total)),
        inherit.aes = FALSE,
        hjust = 0, size = 5, fontface = "bold", color = "black"
      ) +
      labs(
        x="Blood Pressure Trait",
        y="Proportion",
        fill="MR composition"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08)), labels = scales::percent) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x=element_text(angle=0, hjust=0.5, size=14, color="black"),
        axis.text.y=element_text(size=14, color="black"),
        axis.title=element_text(size=14),
        plot.title=element_text(face="bold", hjust=0.5, size=16),
        legend.position="top",
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank()
      ) +
      coord_flip()
    
  }
  
  # ========== 6. 总体Regression vs MR饼图 ==========
  cat("Generating BP overall Regression vs MR pie chart...\n")
  
  # 为每个蛋白标记其在所有表型中的状态
  merged$any_obs <- 0
  merged$any_obs_union <- 0
  
  for (trait in obs_traits) {
    obs_col <- paste0(trait, ".obs")
    if (obs_col %in% colnames(merged)) {
      merged$any_obs <- pmax(merged$any_obs, merged[[obs_col]], na.rm = TRUE)
    }
    
    mr_trait <- obs_mr_mapping$mr[obs_mr_mapping$obs == trait]
    union_col <- paste0(mr_trait, ".union")
    union_exists <- !is.na(mr_trait) && union_col %in% colnames(merged)
    
    if (union_exists && obs_col %in% colnames(merged)) {
      # OBS和UNION都显著才算
      obs_union_both <- merged[[obs_col]] == 1 & merged[[union_col]] == 1
      obs_union_both[is.na(obs_union_both)] <- FALSE
      merged$any_obs_union <- pmax(merged$any_obs_union, as.numeric(obs_union_both), na.rm = TRUE)
    } else if (obs_col %in% colnames(merged)) {
      # 没有union数据的表型，用obs
      merged$any_obs_union <- pmax(merged$any_obs_union, merged[[obs_col]], na.rm = TRUE)
    }
  }
  
  # 统计总数
  regression_only_count <- sum(merged$any_obs == 1 & merged$any_obs_union == 0)
  regression_mr_count <- sum(merged$any_obs == 1 & merged$any_obs_union == 1)
  total_regression <- regression_only_count + regression_mr_count
  
  if (total_regression > 0) {
    # 创建饼图数据框
    pie_data <- data.frame(
      type = c("Regression only", "Regression + MR"),
      count = c(regression_only_count, regression_mr_count),
      stringsAsFactors = FALSE
    )
    
    # 计算百分比和标签位置
    pie_data$percentage <- pie_data$count / total_regression * 100
    pie_data$label <- sprintf("%d (%.1f%%)", 
                              pie_data$count, 
                              pie_data$percentage)
    
    # 设置因子顺序
    pie_data$type <- factor(pie_data$type, levels = c("Regression + MR", "Regression only"))
    
    # 计算标签位置
    pie_data$pos <- cumsum(pie_data$count) - pie_data$count / 2
    
    # 绘制饼图
    plot_pie <- ggplot(pie_data, aes(x = "", y = count, fill = type)) +
      geom_bar(stat = "identity", width = 1, color = "white", size = 2) +
      coord_polar("y", start = 0) +
      scale_fill_manual(
        values = c(
          "Regression + MR" = "#9BC5E4",
          "Regression only" = "#efa3a3"
        )
      ) +
      geom_text(
        aes(y = pos, label = label),
        color = "black",
        size = 6,
        fontface = "bold"
      ) +
      labs(
        title = sprintf("Blood Pressure: Regression vs MR\n(Total: n=%d)", total_regression),
        fill = "Significant in"
      ) +
      theme_void(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18, margin = margin(b = 20)),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"),
        legend.key.size = unit(1, "cm")
      )
    
    ggsave(
      filename = file.path(upset_venn_path, "6.BP_overall_regression_mr_pie.png"),
      plot = plot_pie,
      width = 12,
      height = 10,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(upset_venn_path, "6.BP_overall_regression_mr_pie.pdf"),
      plot = plot_pie,
      width = 12,
      height = 10
    )
    
    cat("  Saved: 6.BP_overall_regression_mr_pie.png/pdf\n")
    cat(sprintf("  Total proteins: %d (Regression only: %d, Regression+MR: %d, MR validation rate: %.1f%%)\n",
               total_regression, regression_only_count, regression_mr_count,
               100 * regression_mr_count / total_regression))
  }
  
  # ========== 7. OBS、OBS+UNION和OBS+UNION+HYPRCOLOC三类堆叠条形图 ==========
  cat("Generating BP OBS vs UNION vs HYPRCOLOC stacked bar plot...\n")
  
  # 创建数据框
  obs_union_hyprcoloc_counts <- data.frame(
    trait=character(), 
    obs_only=numeric(),
    union_only=numeric(),
    union_hyprcoloc=numeric(),
    stringsAsFactors=FALSE
  )
  
  for (trait in obs_traits) {
    obs_col <- paste0(trait, ".obs")
    if (!(obs_col %in% colnames(merged))) next
    
    obs_sig <- merged[[obs_col]] == 1
    mr_trait <- obs_mr_mapping$mr[obs_mr_mapping$obs == trait]
    
    # 检查是否有union和hyprcoloc列
    union_col <- paste0(mr_trait, ".union")
    hyprcoloc_col <- paste0(mr_trait, ".hyprcoloc")
    union_exists <- !is.na(mr_trait) && union_col %in% colnames(merged)
    hyprcoloc_exists <- !is.na(mr_trait) && hyprcoloc_col %in% colnames(merged)
    
    if (union_exists && hyprcoloc_exists) {
      union_sig <- merged[[union_col]] == 1
      union_sig[is.na(union_sig)] <- FALSE
      hyprcoloc_sig <- merged[[hyprcoloc_col]] == 1
      hyprcoloc_sig[is.na(hyprcoloc_sig)] <- FALSE
      
      # OBS显著但UNION不显著
      obs_only_count <- sum(obs_sig & !union_sig)
      # OBS+UNION显著但HYPRCOLOC不显著
      union_only_count <- sum(obs_sig & union_sig & !hyprcoloc_sig)
      # OBS+UNION+HYPRCOLOC都显著
      union_hyprcoloc_count <- sum(obs_sig & union_sig & hyprcoloc_sig)
      
      trait_display <- ifelse(trait == "hypertension", "Hypertension", toupper(trait))
      obs_union_hyprcoloc_counts <- rbind(obs_union_hyprcoloc_counts, data.frame(
        trait = trait_display,
        obs_only = obs_only_count,
        union_only = union_only_count,
        union_hyprcoloc = union_hyprcoloc_count,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if(nrow(obs_union_hyprcoloc_counts) > 0) {
    # 按总数排序
    obs_union_hyprcoloc_counts$total <- obs_union_hyprcoloc_counts$obs_only + 
                                         obs_union_hyprcoloc_counts$union_only + 
                                         obs_union_hyprcoloc_counts$union_hyprcoloc
    obs_union_hyprcoloc_counts <- obs_union_hyprcoloc_counts[order(obs_union_hyprcoloc_counts$total, decreasing = TRUE), ]
    trait_order_ouh <- obs_union_hyprcoloc_counts$trait
    obs_union_hyprcoloc_counts$trait <- factor(obs_union_hyprcoloc_counts$trait, levels = trait_order_ouh)
    
    # 转换为长格式
    obs_union_hyprcoloc_long <- reshape2::melt(
      obs_union_hyprcoloc_counts[,c("trait", "obs_only", "union_only", "union_hyprcoloc")], 
      id.vars="trait", variable.name="type", value.name="count"
    )
    
    # 设置顺序：从底到顶为 obs_only -> union_only -> union_hyprcoloc
    obs_union_hyprcoloc_long$type <- factor(obs_union_hyprcoloc_long$type, 
                                            levels=c("union_hyprcoloc", "union_only", "obs_only"))
    obs_union_hyprcoloc_long$trait <- factor(obs_union_hyprcoloc_long$trait, levels = trait_order_ouh)
    
    # 绘制堆叠条形图
    plot_obs_union_hyprcoloc <- ggplot(obs_union_hyprcoloc_long, aes(x=trait, y=count, fill=type)) +
      geom_bar(stat="identity", position="stack", width=0.7) +
      scale_fill_manual(
        values=c(
          "union_hyprcoloc"="#B58CD1",
          "union_only"="#9BC5E4",
          "obs_only"="#efa3a3"
        ),
        labels=c(
          "Regression+MR+Colocalization",
          "Regression+MR",
          "Regression only"
        )
      ) +
      geom_text(
        data = subset(obs_union_hyprcoloc_long, count > 0),
        aes(label = count),
        color = "black", size = 5.2, fontface = "bold",
        position = position_stack(vjust = 0.5)
      ) +
      geom_text(
        data = obs_union_hyprcoloc_counts %>% filter(total > 0),
        aes(x = trait, y = total + max(total)*0.03, label = paste0("n=", total)),
        inherit.aes = FALSE,
        hjust = 0, size = 5, fontface = "bold", color = "black"
      ) +
      labs(
        x="Blood Pressure Trait",
        y="Number of Significant Proteins",
        fill="Significant in"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x=element_text(angle=0, hjust=0.5, size=14, color="black"),
        axis.text.y=element_text(size=14, color="black"),
        axis.title=element_text(size=14),
        plot.title=element_text(face="bold", hjust=0.5, size=16),
        legend.position="top",
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank()
      ) +
      coord_flip()
    
    ggsave(filename=file.path(upset_venn_path, "7.BP_obs_union_hyprcoloc_stacked.png"), 
           plot=plot_obs_union_hyprcoloc, width=14, height=10, dpi=300)
    ggsave(filename=file.path(upset_venn_path, "7.BP_obs_union_hyprcoloc_stacked.pdf"), 
           plot=plot_obs_union_hyprcoloc, width=14, height=10)
    
    cat("  Saved: 7.BP_obs_union_hyprcoloc_stacked.png/pdf\n")
  }
  
  
  # ========== 8. BP子表型的OBS+UNION+HYPRCOLOC Venn图 ==========
  cat("Generating BP OBS+UNION+HYPRCOLOC sub-phenotypes Venn diagram...\n")
  bp_hyprcoloc_venn_list <- list()
  
  for (trait in obs_traits) {
    obs_col <- paste0(trait, ".obs")
    if (!(obs_col %in% colnames(merged))) next
    
    mr_trait <- obs_mr_mapping$mr[obs_mr_mapping$obs == trait]
    union_col <- paste0(mr_trait, ".union")
    hyprcoloc_col <- paste0(mr_trait, ".hyprcoloc")
    
    if (!is.na(mr_trait) && union_col %in% colnames(merged) && hyprcoloc_col %in% colnames(merged)) {
      label <- trait_labels[[trait]]
      if(is.null(label)) label <- toupper(trait)
      # OBS+UNION+HYPRCOLOC都显著的蛋白
      bp_hyprcoloc_venn_list[[label]] <- merged$protein_id[
        merged[[obs_col]] == 1 & merged[[union_col]] == 1 & merged[[hyprcoloc_col]] == 1
      ]
    }
  }
  
  bp_hyprcoloc_venn_list <- bp_hyprcoloc_venn_list[sapply(bp_hyprcoloc_venn_list, length) > 0]
  
  if (length(bp_hyprcoloc_venn_list) >= 2) {
    venn_plot_hyprcoloc <- ggVennDiagram(
      bp_hyprcoloc_venn_list,
      set_color = set2_colors[1:length(bp_hyprcoloc_venn_list)],
      set_size = 6,
      label_alpha = 0,
      label = "count",
      label_color = "black",
      label_size = 5,
      edge_size = 0.5
    ) + 
      #scale_fill_gradientn(colors = c('#DEEBF7', '#4E79A7', '#4292C6')) +
      scale_fill_gradientn(colors = c('#FFE5D9', '#FFB4A2', '#E76F51'))+
      scale_x_continuous(expand = expansion(mult = 0.15)) +
      scale_y_continuous(expand = expansion(mult = 0.15)) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      ) +
      labs(title = paste0("Blood Pressure MR+Colocalization Sub-phenotypes Overlap (", model, ")"))
    
    png(file.path(upset_venn_path,"8.BP_union_hyprcoloc_venn.png"), width = 14, height = 12, units = "in", res = 300)
    print(venn_plot_hyprcoloc)
    dev.off()
    
    pdf(file.path(upset_venn_path,"8.BP_union_hyprcoloc_venn.pdf"), width = 14, height = 12)
    print(venn_plot_hyprcoloc)
    dev.off()
    
    cat("  Saved: 8.BP_union_hyprcoloc_venn.png/pdf\n")
  }
  
  
  # ========== 9. 总体Regression vs MR vs Colocalization三类饼图 ==========
  cat("Generating BP overall Regression vs MR vs Colocalization pie chart...\n")
  
  # 为每个蛋白标记其在任何trait中的最高验证级别
  merged$has_any_obs <- 0
  merged$has_any_union <- 0
  merged$has_any_hyprcoloc <- 0
  
  for (trait in obs_traits) {
    obs_col <- paste0(trait, ".obs")
    if (!(obs_col %in% colnames(merged))) next
    
    # 标记是否有OBS显著
    obs_sig <- merged[[obs_col]] == 1
    merged$has_any_obs <- pmax(merged$has_any_obs, as.numeric(obs_sig), na.rm = TRUE)
    
    mr_trait <- obs_mr_mapping$mr[obs_mr_mapping$obs == trait]
    union_col <- paste0(mr_trait, ".union")
    hyprcoloc_col <- paste0(mr_trait, ".hyprcoloc")
    
    union_exists <- !is.na(mr_trait) && union_col %in% colnames(merged)
    hyprcoloc_exists <- !is.na(mr_trait) && hyprcoloc_col %in% colnames(merged)
    
    if (union_exists) {
      union_sig <- merged[[union_col]] == 1
      union_sig[is.na(union_sig)] <- FALSE
      # 只有当OBS和UNION都显著时才标记has_any_union
      obs_union_both <- obs_sig & union_sig
      merged$has_any_union <- pmax(merged$has_any_union, as.numeric(obs_union_both), na.rm = TRUE)
    }
    
    if (hyprcoloc_exists && union_exists) {
      union_sig <- merged[[union_col]] == 1
      union_sig[is.na(union_sig)] <- FALSE
      hyprcoloc_sig <- merged[[hyprcoloc_col]] == 1
      hyprcoloc_sig[is.na(hyprcoloc_sig)] <- FALSE
      # 只有当OBS、UNION和HYPRCOLOC都显著时才标记has_any_hyprcoloc
      obs_union_hyprcoloc_all <- obs_sig & union_sig & hyprcoloc_sig
      merged$has_any_hyprcoloc <- pmax(merged$has_any_hyprcoloc, as.numeric(obs_union_hyprcoloc_all), na.rm = TRUE)
    }
  }
  
  # 根据最高验证级别进行互斥分类
  obs_only_count_pie <- sum(merged$has_any_obs == 1 & merged$has_any_union == 0)
  union_only_count_pie <- sum(merged$has_any_obs == 1 & merged$has_any_union == 1 & merged$has_any_hyprcoloc == 0)
  union_hyprcoloc_count_pie <- sum(merged$has_any_obs == 1 & merged$has_any_union == 1 & merged$has_any_hyprcoloc == 1)
  total_all <- obs_only_count_pie + union_only_count_pie + union_hyprcoloc_count_pie
  
  if (total_all > 0) {
    # 创建饼图数据框
    pie_data_all <- data.frame(
      type = c("Regression only", "Regression + MR", "Regression + MR + Colocalization"),
      count = c(obs_only_count_pie, union_only_count_pie, union_hyprcoloc_count_pie),
      stringsAsFactors = FALSE
    )
    
    # 计算百分比和标签位置
    pie_data_all$percentage <- pie_data_all$count / total_all * 100
    pie_data_all$label <- sprintf("%d (%.1f%%)", 
                                  pie_data_all$count, 
                                  pie_data_all$percentage)
    
    # 设置因子顺序
    pie_data_all$type <- factor(pie_data_all$type, 
                                 levels = c("Regression + MR + Colocalization", 
                                           "Regression + MR", 
                                           "Regression only"))
    
    # 计算标签位置
    pie_data_all$pos <- cumsum(pie_data_all$count) - pie_data_all$count / 2
    
    # 绘制饼图
    plot_pie_all <- ggplot(pie_data_all, aes(x = "", y = count, fill = type)) +
      geom_bar(stat = "identity", width = 1, color = "white", size = 2) +
      coord_polar("y", start = 0) +
      scale_fill_manual(
        values = c(
          "Regression + MR + Colocalization" = "#B58CD1",
          "Regression + MR" = "#9BC5E4",
          "Regression only" = "#efa3a3"
        )
      ) +
      geom_text(
        aes(y = pos, label = label),
        color = "black",
        size = 6,
        fontface = "bold"
      ) +
      labs(
        title = sprintf("Blood Pressure: Multi-method Validation\n(Total: n=%d)", total_all),
        fill = "Significant in"
      ) +
      theme_void(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18, margin = margin(b = 20)),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"),
        legend.key.size = unit(1, "cm")
      )
    
    ggsave(
      filename = file.path(upset_venn_path, "9.BP_regression_mr_colocalization_pie.png"),
      plot = plot_pie_all,
      width = 12,
      height = 10,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(upset_venn_path, "9.BP_regression_mr_colocalization_pie.pdf"),
      plot = plot_pie_all,
      width = 12,
      height = 10
    )
    
    cat("  Saved: 9.BP_regression_mr_colocalization_pie.png/pdf\n")
    cat(sprintf("  Total proteins: %d (Regression only: %d, Regression+MR: %d, Regression+MR+Colocalization: %d)\n",
               total_all, obs_only_count_pie, union_only_count_pie, union_hyprcoloc_count_pie))
    cat(sprintf("  MR validation rate: %.1f%%, Colocalization validation rate: %.1f%%\n",
               100 * (union_only_count_pie + union_hyprcoloc_count_pie) / total_all,
               100 * union_hyprcoloc_count_pie / total_all))
  }
  
  
  # ========== 组合图4和图5 ==========
  if(exists("plot_obs_union") && exists("plot_composition")) {
    cat("Combining plots 4 and 5...\n")
    
    # 使用patchwork横向组合两个图
    combined_plot <- plot_obs_union + plot_composition + 
      plot_layout(widths = c(1, 1)) +
      plot_annotation(
        title = 'Number of Proteins Identified Regression and Mendelian Randomization',
        theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
      )
    
    ggsave(filename=file.path(upset_venn_path, "4_5.BP_obs_union_combined.png"), 
          plot=combined_plot, width=28, height=10, dpi=300)
    ggsave(filename=file.path(upset_venn_path, "4_5.BP_obs_union_combined.pdf"), 
          plot=combined_plot, width=28, height=10)
    
    cat("  Saved: 4_5.BP_obs_union_combined.png/pdf\n")
  }
  
  # ========== 10. 输出三个层级的血压蛋白列表 ==========
  cat("Generating three-level BP protein lists...\n")
  
  # 读取注释信息
  anno_info <- read.csv("E:/results/anno.info_updat.csv", header = TRUE, stringsAsFactors = FALSE)
  
  # 创建Excel工作簿
  library(openxlsx)
  wb_levels <- createWorkbook()
  
  # Level 1: Regression (obs) - 任何一个BP表型的obs显著
  level1_proteins <- merged$protein_id[merged$has_any_obs == 1]
  if (length(level1_proteins) > 0) {
    level1_df <- data.frame(
      protein_id = level1_proteins,
      stringsAsFactors = FALSE
    )
    level1_df <- left_join(level1_df, 
                          anno_info[, c("AptName", "Entrez.Gene.Name")], 
                          by = c("protein_id" = "AptName"))
    level1_df <- level1_df[!is.na(level1_df$Entrez.Gene.Name), ]
    level1_df$level <- "Regression"
    level1_df <- level1_df[, c("level", "protein_id", "Entrez.Gene.Name")]
    
    addWorksheet(wb_levels, "1_Regression")
    writeData(wb_levels, "1_Regression", level1_df)
    cat(sprintf("  Level 1 (Regression): %d proteins\n", nrow(level1_df)))
  }
  
  # Level 2: Regression+MR (obs+union) - 任何一个BP表型的obs+union都显著
  level2_proteins <- merged$protein_id[merged$has_any_union == 1]
  if (length(level2_proteins) > 0) {
    level2_df <- data.frame(
      protein_id = level2_proteins,
      stringsAsFactors = FALSE
    )
    level2_df <- left_join(level2_df, 
                          anno_info[, c("AptName", "Entrez.Gene.Name")], 
                          by = c("protein_id" = "AptName"))
    level2_df <- level2_df[!is.na(level2_df$Entrez.Gene.Name), ]
    level2_df$level <- "Regression+MR"
    level2_df <- level2_df[, c("level", "protein_id", "Entrez.Gene.Name")]
    
    addWorksheet(wb_levels, "2_Regression_MR")
    writeData(wb_levels, "2_Regression_MR", level2_df)
    cat(sprintf("  Level 2 (Regression+MR): %d proteins\n", nrow(level2_df)))
  }
  
  # Level 3: Regression+MR+Colocalization (obs+union+hyprcoloc) - 任何一个BP表型的三者都显著
  level3_proteins <- merged$protein_id[merged$has_any_hyprcoloc == 1]
  if (length(level3_proteins) > 0) {
    level3_df <- data.frame(
      protein_id = level3_proteins,
      stringsAsFactors = FALSE
    )
    level3_df <- left_join(level3_df, 
                          anno_info[, c("AptName", "Entrez.Gene.Name")], 
                          by = c("protein_id" = "AptName"))
    level3_df <- level3_df[!is.na(level3_df$Entrez.Gene.Name), ]
    level3_df$level <- "Regression+MR+Colocalization"
    level3_df <- level3_df[, c("level", "protein_id", "Entrez.Gene.Name")]
    
    addWorksheet(wb_levels, "3_Regression_MR_Coloc")
    writeData(wb_levels, "3_Regression_MR_Coloc", level3_df)
    cat(sprintf("  Level 3 (Regression+MR+Colocalization): %d proteins\n", nrow(level3_df)))
  }
  
  # 保存Excel文件
  if (length(names(wb_levels)) > 0) {
    saveWorkbook(wb_levels, 
                file.path(upset_venn_path, "BP_three_levels_protein_lists.xlsx"), 
                overwrite = TRUE)
    cat("  Saved: BP_three_levels_protein_lists.xlsx\n\n")
  }
  
  cat(paste("\n", rep("=", 70), collapse = ""), "\n")
  cat(paste("Blood Pressure UpSet/Venn plots completed for", model, "\n"))
  cat(paste("All plots saved to:", upset_venn_path, "\n"))
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
}

cat("\n*** All Blood Pressure UpSet/Venn plots generated successfully! ***\n")

