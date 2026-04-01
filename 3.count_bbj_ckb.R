library(dplyr)
library(tidyr)
library(xlsx)
library(openxlsx)
#设置路径
out_dir <- "D:/OneDrive/工作/1.工作/blood pressure"
anno_info <- read.csv("E:/results/anno.info_updat.csv", header = TRUE, stringsAsFactors = FALSE)

# 读取数据
for (model in c("model_1","model_2","model_3")){
  # 为每个model创建输出目录 combined_results/model1 等
  out_path <- file.path(out_dir, model)
  obs_wide <- read.table(file.path(out_path, "obs_protein_trait_group_status.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(obs_wide)[2:7] <- paste0(colnames(obs_wide)[2:7], ".obs")
  mr_wide <- read.table(file.path(out_dir, "mr", "mr_protein_trait_group_status.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  #mr_wide <- select(mr_wide,-c(bbj_sbp,ckb_sbp,bbj_dbp,ckb_dbp))
  mr_wide <- mr_wide[mr_wide$method %in% c("IVW_Egger"), ]
  colnames(mr_wide)[3:8] <- paste0(colnames(mr_wide)[3:8], ".mr")
  hyprcoloc_wide <- read.table(file.path(out_dir,"hyprcoloc", "hyprcoloc_protein_trait_group_status.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  hyprcoloc_wide <- hyprcoloc_wide %>%dplyr::select(-Entrez.Gene.Name)%>% dplyr::rename(protein_id = protein_seq_name)
  colnames(hyprcoloc_wide)[2:7] <- paste0(colnames(hyprcoloc_wide)[2:7], ".hyprcoloc")
  gsmr <- read.table(file.path(out_dir, "gsmr", "gsmr_significant_protein_status.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  #gsmr <- select(gsmr,-c(bbj_sbp,ckb_sbp,bbj_dbp,ckb_dbp))
  colnames(gsmr)[2:7] <- paste0(colnames(gsmr)[2:7], ".gsmr")
  mr_raw_file <- file.path(out_dir, "mr", "mr_all_results.txt")
  gsmr_raw_file <- file.path(out_dir, "gsmr", "gsmr_all_results.txt")
  obs_raw_file <- file.path(out_path, "obs_all_results.txt")
  
  mr_raw_data <- read.table(mr_raw_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  mr_raw_data <- mr_raw_data[mr_raw_data$method == "IVW_Egger", ]
  gsmr_raw_data <- read.table(gsmr_raw_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  obs_raw_data <- read.table(obs_raw_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  #*************************************** 继续处理     ******************************
  #***********************************************************************************
  # MR宽表 - 保留所有列
  mr_wide_all <- mr_wide
  mr_wide_all[is.na(mr_wide_all)] <- 0
  
  gsmr_wide_all <- gsmr
  gsmr_wide_all[is.na(gsmr_wide_all)] <- 0
  
  mr_gsmr_union <- data.frame(
    protein_id = unique(c(mr_wide_all$protein_id, gsmr_wide_all$protein_id)),
    stringsAsFactors = FALSE
  )

  allowed_combined_traits <- c("sbp", "dbp", "pp", "map")
  source_prefixes <- c("bbj_", "ckb_", "combined_")

  for (prefix in source_prefixes) {
    prefix_mr_cols <- grep(paste0("^", prefix), colnames(mr_wide_all), value = TRUE)
    for (mr_col in prefix_mr_cols) {
      trait <- sub("\\.mr$", "", sub(paste0("^", prefix), "", mr_col))
      if (prefix == "combined_" && !(trait %in% allowed_combined_traits)) {
        next
      }
      gsmr_col <- paste0(prefix, trait, ".gsmr")
      if (!(gsmr_col %in% colnames(gsmr_wide_all))) {
        next
      }
      temp <- full_join(
        mr_wide_all[, c("protein_id", mr_col)],
        gsmr_wide_all[, c("protein_id", gsmr_col)],
        by = "protein_id"
      )
      temp[[mr_col]][is.na(temp[[mr_col]])] <- 0
      temp[[gsmr_col]][is.na(temp[[gsmr_col]])] <- 0
      union_values <- pmax(temp[[mr_col]], temp[[gsmr_col]])
      union_col_name <- paste0(prefix, trait, ".union")
      temp_df <- data.frame(protein_id = temp$protein_id, value = union_values)
      colnames(temp_df)[2] <- union_col_name
      mr_gsmr_union <- left_join(mr_gsmr_union, temp_df, by = "protein_id")
    }
  }

  mr_gsmr_union[is.na(mr_gsmr_union)] <- 0

  # 合并所有结果
  merged_all <- obs_wide %>%
    left_join(mr_wide_all, by = "protein_id") %>%
    left_join(gsmr_wide_all, by = "protein_id") %>%
    left_join(mr_gsmr_union, by = "protein_id") %>%
    left_join(hyprcoloc_wide, by = "protein_id")
  merged_all <- merged_all %>% dplyr::select(protein_id, Entrez.Gene.Name, everything() )

  write.table(merged_all, file = file.path(out_path,"merged_obs_allmr_hyprcoloc.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  # 从merged_all中提取所有带.obs后缀的列名（这些是具体的表型）
  obs_cols <- grep("\\.obs$", colnames(merged_all), value = TRUE)
  trait_names <- gsub("\\.obs$", "", obs_cols)  # 去掉.obs后缀得到表型名（如sbp, dbp, hypertension等）
  
  # 求交集的函数：对于指定的trait和数据源，求指定维度的交集
  get_trait_intersections <- function(df, trait, source_prefix, suffixes, seq_col = "protein_id") {
    # 构建obs列名（obs列没有前缀）
    obs_col <- paste0(trait, ".obs")
    
    # 构建其他列名（mr/gsmr/hyprcoloc有前缀）
    other_cols <- paste0(source_prefix, trait, suffixes[suffixes != ".obs"])
    
    # 合并所有需要的列
    status_cols <- c(obs_col, other_cols)
    
    # 检查所有列是否存在
    existing_cols <- intersect(status_cols, colnames(df))
    if (length(existing_cols) != length(status_cols)) {
      return(data.frame())
    }
    
    # 选择所有指定列都为1的行
    sel <- df %>% 
      filter(rowSums(dplyr::select(., all_of(existing_cols)) == 1) == length(existing_cols))
    
    if (nrow(sel) > 0) {
      sel %>% dplyr::select(all_of(seq_col))
    } else {
      data.frame()
    }
  }

  # 构建结果列表：按trait和数据源组织
  intersections_list <- list()
  
  # 对每个trait和每个数据源进行处理
  for (trait in trait_names) {
    for (prefix in source_prefixes) {
      # 构建标识符（如bbj_sbp, ckb_dbp等）
      trait_id <- paste0(prefix, trait)
      
      # 仅obs
      intersections_list[[paste0(trait_id, "_obs")]] <- 
        get_trait_intersections(merged_all, trait, prefix, c(".obs"))
      
      # obs + mr
      intersections_list[[paste0(trait_id, "_obs_mr")]] <- 
        get_trait_intersections(merged_all, trait, prefix, c(".obs", ".mr"))
      
      # obs + gsmr
      intersections_list[[paste0(trait_id, "_obs_gsmr")]] <- 
        get_trait_intersections(merged_all, trait, prefix, c(".obs", ".gsmr"))
      
      # obs + (mr 或 gsmr)
      intersections_list[[paste0(trait_id, "_obs_mr_or_gsmr")]] <- 
        rbind(
          get_trait_intersections(merged_all, trait, prefix, c(".obs", ".mr")),
          get_trait_intersections(merged_all, trait, prefix, c(".obs", ".gsmr"))
        ) %>% distinct()
      
      # obs + mr + hyprcoloc
      intersections_list[[paste0(trait_id, "_obs_mr_hyp")]] <- 
        get_trait_intersections(merged_all, trait, prefix, c(".obs", ".mr", ".hyprcoloc"))
      
      # obs + gsmr + hyprcoloc
      intersections_list[[paste0(trait_id, "_obs_gsmr_hyp")]] <- 
        get_trait_intersections(merged_all, trait, prefix, c(".obs", ".gsmr", ".hyprcoloc"))
      
      # obs + (mr 或 gsmr) + hyprcoloc
      intersections_list[[paste0(trait_id, "_obs_mr_or_gsmr_hyp")]] <- 
        rbind(
          get_trait_intersections(merged_all, trait, prefix, c(".obs", ".mr", ".hyprcoloc")),
          get_trait_intersections(merged_all, trait, prefix, c(".obs", ".gsmr", ".hyprcoloc"))
        ) %>% distinct()
    }
  }
  
  # 移除空的结果
  intersections_list <- intersections_list[sapply(intersections_list, nrow) > 0]
  
  # 计算blood_pressure作为所有trait的并集
  for (prefix in source_prefixes) {
    # 收集该数据源下所有trait的结果
    trait_id <- paste0(prefix, "blood_pressure")
    
    # 对每种组合类型都计算并集
    combination_types <- c("_obs", "_obs_mr", "_obs_gsmr", "_obs_mr_or_gsmr", 
                          "_obs_mr_hyp", "_obs_gsmr_hyp", "_obs_mr_or_gsmr_hyp")
    
    for (comb_type in combination_types) {
      combined_df <- data.frame()
      for (trait in trait_names) {
        key <- paste0(prefix, trait, comb_type)
        if (key %in% names(intersections_list)) {
          combined_df <- rbind(combined_df, intersections_list[[key]])
        }
      }
      if (nrow(combined_df) > 0) {
        intersections_list[[paste0(trait_id, comb_type)]] <- combined_df %>% distinct()
      }
    }
  }

  # 计算各类交集的蛋白质数量
  intersection_counts <- data.frame()
  for (nm in names(intersections_list)) {
    df <- intersections_list[[nm]]
    if (nrow(df) > 0) {
      cnt <- n_distinct(df$protein_id)
      # 统一命名规则
      display_name <- nm
      display_name <- gsub("obs_mr_or_gsmr_hyp", "obs_union_hypr", display_name)
      display_name <- gsub("obs_mr_or_gsmr", "obs_union", display_name)
      display_name <- gsub("obs_mr_hyp", "obs_mr_hypr", display_name)
      display_name <- gsub("obs_gsmr_hyp", "obs_gsmr_hypr", display_name)
      display_name <- gsub("blood_pressure", "bp", display_name)
      display_name <- gsub("hypertension", "ht", display_name)
      
      intersection_counts <- rbind(intersection_counts, 
                                   data.frame(intersection = display_name, count = cnt))
    }
  }
  
  write.csv(intersection_counts, file = file.path(out_path,"intersection_counts_combined.csv"), row.names = FALSE)
  
  # 创建按数据源分组的汇总统计表
  source_summary <- data.frame()
  
  for (prefix in source_prefixes) {
    source_name <- gsub("_$", "", prefix)  # 去掉末尾的下划线
    
    # 对每个trait统计
    for (trait in trait_names) {
      trait_id <- paste0(prefix, trait)
      
      # obs + mr
      obs_mr_key <- paste0(trait_id, "_obs_mr")
      obs_mr_count <- if (obs_mr_key %in% names(intersections_list)) {
        n_distinct(intersections_list[[obs_mr_key]]$protein_id)
      } else { 0 }
      
      # obs + gsmr
      obs_gsmr_key <- paste0(trait_id, "_obs_gsmr")
      obs_gsmr_count <- if (obs_gsmr_key %in% names(intersections_list)) {
        n_distinct(intersections_list[[obs_gsmr_key]]$protein_id)
      } else { 0 }
      
      # obs_union (obs + mr或gsmr)
      obs_union_key <- paste0(trait_id, "_obs_mr_or_gsmr")
      obs_union_count <- if (obs_union_key %in% names(intersections_list)) {
        n_distinct(intersections_list[[obs_union_key]]$protein_id)
      } else { 0 }
      
      # obs + mr + hyprcoloc
      obs_mr_hypr_key <- paste0(trait_id, "_obs_mr_hyp")
      obs_mr_hypr_count <- if (obs_mr_hypr_key %in% names(intersections_list)) {
        n_distinct(intersections_list[[obs_mr_hypr_key]]$protein_id)
      } else { 0 }
      
      # obs + gsmr + hyprcoloc
      obs_gsmr_hypr_key <- paste0(trait_id, "_obs_gsmr_hyp")
      obs_gsmr_hypr_count <- if (obs_gsmr_hypr_key %in% names(intersections_list)) {
        n_distinct(intersections_list[[obs_gsmr_hypr_key]]$protein_id)
      } else { 0 }
      
      # obs_union_hypr (obs + mr或gsmr + hyprcoloc)
      obs_union_hypr_key <- paste0(trait_id, "_obs_mr_or_gsmr_hyp")
      obs_union_hypr_count <- if (obs_union_hypr_key %in% names(intersections_list)) {
        n_distinct(intersections_list[[obs_union_hypr_key]]$protein_id)
      } else { 0 }
      
      if (obs_mr_count > 0 || obs_gsmr_count > 0 || obs_union_count > 0 || 
          obs_mr_hypr_count > 0 || obs_gsmr_hypr_count > 0 || obs_union_hypr_count > 0) {
        source_summary <- rbind(source_summary,
                               data.frame(
                                 data_source = source_name,
                                 trait = trait,
                                 obs_mr = obs_mr_count,
                                 obs_gsmr = obs_gsmr_count,
                                 obs_union = obs_union_count,
                                 obs_mr_hypr = obs_mr_hypr_count,
                                 obs_gsmr_hypr = obs_gsmr_hypr_count,
                                 obs_union_hypr = obs_union_hypr_count,
                                 stringsAsFactors = FALSE
                               ))
      }
    }
    
    # 添加blood_pressure汇总行（所有trait的并集）
    bp_mr_key <- paste0(prefix, "blood_pressure_obs_mr")
    bp_mr_count <- if (bp_mr_key %in% names(intersections_list)) {
      n_distinct(intersections_list[[bp_mr_key]]$protein_id)
    } else { 0 }
    
    bp_gsmr_key <- paste0(prefix, "blood_pressure_obs_gsmr")
    bp_gsmr_count <- if (bp_gsmr_key %in% names(intersections_list)) {
      n_distinct(intersections_list[[bp_gsmr_key]]$protein_id)
    } else { 0 }
    
    bp_union_key <- paste0(prefix, "blood_pressure_obs_mr_or_gsmr")
    bp_union_count <- if (bp_union_key %in% names(intersections_list)) {
      n_distinct(intersections_list[[bp_union_key]]$protein_id)
    } else { 0 }
    
    bp_mr_hypr_key <- paste0(prefix, "blood_pressure_obs_mr_hyp")
    bp_mr_hypr_count <- if (bp_mr_hypr_key %in% names(intersections_list)) {
      n_distinct(intersections_list[[bp_mr_hypr_key]]$protein_id)
    } else { 0 }
    
    bp_gsmr_hypr_key <- paste0(prefix, "blood_pressure_obs_gsmr_hyp")
    bp_gsmr_hypr_count <- if (bp_gsmr_hypr_key %in% names(intersections_list)) {
      n_distinct(intersections_list[[bp_gsmr_hypr_key]]$protein_id)
    } else { 0 }
    
    bp_union_hypr_key <- paste0(prefix, "blood_pressure_obs_mr_or_gsmr_hyp")
    bp_union_hypr_count <- if (bp_union_hypr_key %in% names(intersections_list)) {
      n_distinct(intersections_list[[bp_union_hypr_key]]$protein_id)
    } else { 0 }
    
    if (bp_mr_count > 0 || bp_gsmr_count > 0 || bp_union_count > 0 || 
        bp_mr_hypr_count > 0 || bp_gsmr_hypr_count > 0 || bp_union_hypr_count > 0) {
      source_summary <- rbind(source_summary,
                             data.frame(
                               data_source = source_name,
                               trait = "blood_pressure_all",
                               obs_mr = bp_mr_count,
                               obs_gsmr = bp_gsmr_count,
                               obs_union = bp_union_count,
                               obs_mr_hypr = bp_mr_hypr_count,
                               obs_gsmr_hypr = bp_gsmr_hypr_count,
                               obs_union_hypr = bp_union_hypr_count,
                               stringsAsFactors = FALSE
                             ))
    }
  }
  
  write.csv(source_summary, file = file.path(out_path,"source_trait_summary.csv"), row.names = FALSE)

  # 为交集蛋白质添加Entrez.Gene.Name和生成Excel文件
  wb2_id <- createWorkbook()
  all_intersections_data <- data.frame()
  
  for (nm in names(intersections_list)) {
    df <- intersections_list[[nm]]
    if (nrow(df) > 0) {
      # 关联基因名
      df <- left_join(df, anno_info[, c("AptName", "Entrez.Gene.Name")], 
                      by = c("protein_id" = "AptName"))
      df <- df[!is.na(df$Entrez.Gene.Name), ]
      
      if (nrow(df) > 0) {
        # 不重复的基因列表（带ID）
        seq_list_with_id <- df %>% 
          dplyr::select(protein_id, Entrez.Gene.Name) %>% 
          distinct()
        
        # 统一命名规则
        display_name <- nm
        display_name <- gsub("obs_mr_or_gsmr_hyp", "obs_union_hypr", display_name)
        display_name <- gsub("obs_mr_or_gsmr", "obs_union", display_name)
        display_name <- gsub("obs_mr_hyp", "obs_mr_hypr", display_name)
        display_name <- gsub("obs_gsmr_hyp", "obs_gsmr_hypr", display_name)
        display_name <- gsub("blood_pressure", "bp", display_name)
        display_name <- gsub("hypertension", "ht", display_name)
        
        # 添加intersection列
        seq_list_with_id$intersection <- display_name
        
        # 合并到总表
        all_intersections_data <- rbind(all_intersections_data, seq_list_with_id)
      }
    }
  }
  
  # 调整列顺序：intersection, protein_id, Entrez.Gene.Name
  if (nrow(all_intersections_data) > 0) {
    all_intersections_data <- all_intersections_data %>%
      dplyr::select(intersection, protein_id, Entrez.Gene.Name)
    
    addWorksheet(wb2_id, "All_Intersections")
    writeData(wb2_id, "All_Intersections", all_intersections_data)
  }
  
  saveWorkbook(wb2_id, file.path(out_path,"intersection_protein_id.xlsx"), overwrite = TRUE)
}
