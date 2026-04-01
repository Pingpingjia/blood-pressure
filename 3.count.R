library(dplyr)
library(tidyr)
library(openxlsx)

# 设置路径
out_dir <- "D:/OneDrive/工作/1.工作/blood pressure"
anno_info <- read.csv("E:/results/anno.info_updat.csv", header = TRUE, stringsAsFactors = FALSE)

# 读取数据
for (model in c("model_1", "model_2", "model_3")) {
  out_path <- file.path(out_dir, model)
  
  # ===== 步骤1: 生成merged_obs_allmr_hyprcoloc.txt文件 =====
  cat("\n========== 模型:", model, "- 生成merged文件 ==========\n")
  
  # 读取各个数据源
  obs_wide <- read.table(file.path(out_path, "obs_protein_trait_group_status.txt"), 
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(obs_wide)[2:7] <- paste0(colnames(obs_wide)[2:7], ".obs")
  
  mr_wide <- read.table(file.path(out_dir, "mr", "mr_protein_trait_group_status.txt"), 
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  mr_wide <- mr_wide[mr_wide$method %in% c("IVW_Egger"), ]
  colnames(mr_wide)[3:8] <- paste0(colnames(mr_wide)[3:8], ".mr")
  
  gsmr <- read.table(file.path(out_dir, "gsmr", "gsmr_significant_protein_status.txt"), 
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gsmr)[2:7] <- paste0(colnames(gsmr)[2:7], ".gsmr")
  
  # read hyprcoloc max posterior probability and convert to 0/1 based on threshold
  hyprcoloc_prob_threshold <- 0.1
  hyprcoloc_max_prob <- read.table(file.path(out_dir, "hyprcoloc", "hyprcoloc_max_posterior_prob.txt"), 
                                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Convert to wide format with 0/1 based on threshold
  hyprcoloc_wide <- hyprcoloc_max_prob %>%
    mutate(value = as.integer(max_posterior_prob >= hyprcoloc_prob_threshold)) %>%
    dplyr::select(protein_seq_name, outcome, value) %>%
    tidyr::pivot_wider(names_from = outcome, values_from = value, values_fill = 0) %>%
    dplyr::rename(protein_id = protein_seq_name)
  
  # Add blood_pressure column
  bp_traits_cols <- intersect(c("ckb_sbp", "ckb_dbp", "ckb_pp", "ckb_map", "ukb_hypertension"), colnames(hyprcoloc_wide))
  if (length(bp_traits_cols) > 0) {
    hyprcoloc_wide$blood_pressure <- as.integer(rowSums(hyprcoloc_wide[, bp_traits_cols, drop = FALSE]) > 0)
  }
  
  colnames(hyprcoloc_wide)[2:ncol(hyprcoloc_wide)] <- paste0(colnames(hyprcoloc_wide)[2:ncol(hyprcoloc_wide)], ".hyprcoloc")
  
  # 将NA替换为0
  mr_wide[is.na(mr_wide)] <- 0
  gsmr[is.na(gsmr)] <- 0
  
  # 创建union列（mr和gsmr的合并，有一个显著即为显著）
  mr_gsmr_union <- data.frame(
    protein_id = unique(c(mr_wide$protein_id, gsmr$protein_id)),
    stringsAsFactors = FALSE
  )
  
  # 对ckb_sbp, ckb_dbp, ckb_pp, ckb_map, ukb_hypertension生成union列
  traits_to_union <- c("ckb_sbp", "ckb_dbp", "ckb_pp", "ckb_map", "ukb_hypertension")
  
  for (trait in traits_to_union) {
    mr_col <- paste0(trait, ".mr")
    gsmr_col <- paste0(trait, ".gsmr")
    
    if (mr_col %in% colnames(mr_wide) && gsmr_col %in% colnames(gsmr)) {
      temp <- full_join(
        mr_wide[, c("protein_id", mr_col)],
        gsmr[, c("protein_id", gsmr_col)],
        by = "protein_id"
      )
      temp[[mr_col]][is.na(temp[[mr_col]])] <- 0
      temp[[gsmr_col]][is.na(temp[[gsmr_col]])] <- 0
      # union: 有一个为1即为1（使用pmax）
      union_values <- pmax(temp[[mr_col]], temp[[gsmr_col]])
      union_col_name <- paste0(trait, ".union")
      temp_df <- data.frame(protein_id = temp$protein_id, value = union_values)
      colnames(temp_df)[2] <- union_col_name
      mr_gsmr_union <- left_join(mr_gsmr_union, temp_df, by = "protein_id")
    }
  }
  
  mr_gsmr_union[is.na(mr_gsmr_union)] <- 0
  
  # 合并所有结果
  merged_all <- obs_wide %>%
    left_join(mr_wide, by = "protein_id") %>%
    left_join(gsmr, by = "protein_id") %>%
    left_join(mr_gsmr_union, by = "protein_id") %>%
    left_join(hyprcoloc_wide, by = "protein_id")
  merged_all <- merged_all %>% dplyr::select(protein_id, Entrez.Gene.Name, everything())
  
  # 保存merged文件
  write.table(merged_all, file = file.path(out_path, "merged_obs_allmr_hyprcoloc.txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("merged_obs_allmr_hyprcoloc.txt 已生成\n")
  
  # ===== 步骤2: 统计5个表型在各层级的蛋白数量 =====
  cat("\n========== 模型:", model, "- 统计蛋白数量 ==========\n")
  
  # 将NA替换为0
  char_cols <- c("protein_id", "Entrez.Gene.Name", "method")
  num_cols <- setdiff(colnames(merged_all), char_cols)
  for (col in num_cols) {
    if (col %in% colnames(merged_all)) {
      merged_all[[col]][is.na(merged_all[[col]])] <- 0
      merged_all[[col]] <- as.numeric(merged_all[[col]])
    }
  }
  
  # 打印调试信息
  cat("总行数:", nrow(merged_all), "\n")
  cat("各表型obs=1的行数:\n")
  cat("  sbp.obs:", sum(merged_all$sbp.obs == 1, na.rm=TRUE), "\n")
  cat("  dbp.obs:", sum(merged_all$dbp.obs == 1, na.rm=TRUE), "\n")
  cat("  pp.obs:", sum(merged_all$pp.obs == 1, na.rm=TRUE), "\n")
  cat("  map.obs:", sum(merged_all$map.obs == 1, na.rm=TRUE), "\n")
  cat("  hypertension.obs:", sum(merged_all$hypertension.obs == 1, na.rm=TRUE), "\n\n")
  
  # 定义5个表型及其对应的数据源列
  # obs列没有前缀，mr/gsmr/union/hyprcoloc列有前缀
  # sbp, dbp, pp, map使用CKB数据，hypertension使用UKB数据
  phenotypes <- list(
    sbp = list(
      obs_col = "sbp.obs",
      mr_col = "ckb_sbp.mr",
      gsmr_col = "ckb_sbp.gsmr",
      union_col = "ckb_sbp.union",
      hyprcoloc_col = "ckb_sbp.hyprcoloc"
    ),
    dbp = list(
      obs_col = "dbp.obs",
      mr_col = "ckb_dbp.mr",
      gsmr_col = "ckb_dbp.gsmr",
      union_col = "ckb_dbp.union",
      hyprcoloc_col = "ckb_dbp.hyprcoloc"
    ),
    pp = list(
      obs_col = "pp.obs",
      mr_col = "ckb_pp.mr",
      gsmr_col = "ckb_pp.gsmr",
      union_col = "ckb_pp.union",
      hyprcoloc_col = "ckb_pp.hyprcoloc"
    ),
    map = list(
      obs_col = "map.obs",
      mr_col = "ckb_map.mr",
      gsmr_col = "ckb_map.gsmr",
      union_col = "ckb_map.union",
      hyprcoloc_col = "ckb_map.hyprcoloc"
    ),
    hypertension = list(
      obs_col = "hypertension.obs",
      mr_col = "ukb_hypertension.mr",
      gsmr_col = "ukb_hypertension.gsmr",
      union_col = "ukb_hypertension.union",
      hyprcoloc_col = "ukb_hypertension.hyprcoloc"
    )
  )
  
  # 用于存储各层级的蛋白质列表和统计结果
  results_list <- list()
  summary_table <- data.frame()
  
  # 对每个表型计算7个层级
  for (pheno_name in names(phenotypes)) {
    pheno <- phenotypes[[pheno_name]]
    
    # 确保列存在
    if (!all(c(pheno$obs_col, pheno$mr_col, pheno$gsmr_col, pheno$union_col, pheno$hyprcoloc_col) %in% colnames(merged_all))) {
      cat("警告：表型", pheno_name, "的某些列不存在，跳过\n")
      next
    }
    
    # 1. obs
    obs_proteins <- merged_all %>%
      filter(!!sym(pheno$obs_col) == 1) %>%
      dplyr::select(protein_id, Entrez.Gene.Name)
    
    # 2. obs + mr
    obs_mr_proteins <- merged_all %>%
      filter(!!sym(pheno$obs_col) == 1 & !!sym(pheno$mr_col) == 1) %>%
      dplyr::select(protein_id, Entrez.Gene.Name)
    
    # 3. obs + gsmr
    obs_gsmr_proteins <- merged_all %>%
      filter(!!sym(pheno$obs_col) == 1 & !!sym(pheno$gsmr_col) == 1) %>%
      dplyr::select(protein_id, Entrez.Gene.Name)
    
    # 4. obs + union (obs + (mr或gsmr)的union列)
    obs_union_proteins <- merged_all %>%
      filter(!!sym(pheno$obs_col) == 1 & !!sym(pheno$union_col) == 1) %>%
      dplyr::select(protein_id, Entrez.Gene.Name)
    
    # 5. obs + mr + hyprcoloc
    obs_mr_hyprcoloc_proteins <- merged_all %>%
      filter(!!sym(pheno$obs_col) == 1 & !!sym(pheno$mr_col) == 1 & !!sym(pheno$hyprcoloc_col) == 1) %>%
      dplyr::select(protein_id, Entrez.Gene.Name)
    
    # 6. obs + gsmr + hyprcoloc
    obs_gsmr_hyprcoloc_proteins <- merged_all %>%
      filter(!!sym(pheno$obs_col) == 1 & !!sym(pheno$gsmr_col) == 1 & !!sym(pheno$hyprcoloc_col) == 1) %>%
      dplyr::select(protein_id, Entrez.Gene.Name)
    
    # 7. obs + union + hyprcoloc (obs + (mr或gsmr)的union列 + hyprcoloc)
    obs_union_hyprcoloc_proteins <- merged_all %>%
      filter(!!sym(pheno$obs_col) == 1 & !!sym(pheno$union_col) == 1 & !!sym(pheno$hyprcoloc_col) == 1) %>%
      dplyr::select(protein_id, Entrez.Gene.Name)
    
    # 存储结果
    results_list[[paste0(pheno_name, "_obs")]] <- obs_proteins
    results_list[[paste0(pheno_name, "_obs_mr")]] <- obs_mr_proteins
    results_list[[paste0(pheno_name, "_obs_gsmr")]] <- obs_gsmr_proteins
    results_list[[paste0(pheno_name, "_obs_union")]] <- obs_union_proteins
    results_list[[paste0(pheno_name, "_obs_mr_hyprcoloc")]] <- obs_mr_hyprcoloc_proteins
    results_list[[paste0(pheno_name, "_obs_gsmr_hyprcoloc")]] <- obs_gsmr_hyprcoloc_proteins
    results_list[[paste0(pheno_name, "_obs_union_hyprcoloc")]] <- obs_union_hyprcoloc_proteins
    
    # 统计数量
    summary_table <- rbind(summary_table, data.frame(
      phenotype = pheno_name,
      obs = nrow(obs_proteins),
      obs_mr = nrow(obs_mr_proteins),
      obs_gsmr = nrow(obs_gsmr_proteins),
      obs_union = nrow(obs_union_proteins),
      obs_mr_hyprcoloc = nrow(obs_mr_hyprcoloc_proteins),
      obs_gsmr_hyprcoloc = nrow(obs_gsmr_hyprcoloc_proteins),
      obs_union_hyprcoloc = nrow(obs_union_hyprcoloc_proteins),
      stringsAsFactors = FALSE
    ))
  }
  
  # 计算五个表型合并去重后的数量
  combined_proteins <- list()
  for (level in c("obs", "obs_mr", "obs_gsmr", "obs_union", "obs_mr_hyprcoloc", "obs_gsmr_hyprcoloc", "obs_union_hyprcoloc")) {
    combined_df <- data.frame()
    for (pheno_name in names(phenotypes)) {
      key <- paste0(pheno_name, "_", level)
      if (key %in% names(results_list) && nrow(results_list[[key]]) > 0) {
        combined_df <- rbind(combined_df, results_list[[key]])
      }
    }
    if (nrow(combined_df) > 0) {
      combined_proteins[[level]] <- combined_df %>% distinct(protein_id, .keep_all = TRUE)
    }
  }
  
  # 添加合并去重后的统计
  summary_table <- rbind(summary_table, data.frame(
    phenotype = "all_combined",
    obs = if ("obs" %in% names(combined_proteins)) nrow(combined_proteins$obs) else 0,
    obs_mr = if ("obs_mr" %in% names(combined_proteins)) nrow(combined_proteins$obs_mr) else 0,
    obs_gsmr = if ("obs_gsmr" %in% names(combined_proteins)) nrow(combined_proteins$obs_gsmr) else 0,
    obs_union = if ("obs_union" %in% names(combined_proteins)) nrow(combined_proteins$obs_union) else 0,
    obs_mr_hyprcoloc = if ("obs_mr_hyprcoloc" %in% names(combined_proteins)) nrow(combined_proteins$obs_mr_hyprcoloc) else 0,
    obs_gsmr_hyprcoloc = if ("obs_gsmr_hyprcoloc" %in% names(combined_proteins)) nrow(combined_proteins$obs_gsmr_hyprcoloc) else 0,
    obs_union_hyprcoloc = if ("obs_union_hyprcoloc" %in% names(combined_proteins)) nrow(combined_proteins$obs_union_hyprcoloc) else 0,
    stringsAsFactors = FALSE
  ))
  
  # 保存统计结果
  write.csv(summary_table, file = file.path(out_path, "phenotype_summary.csv"), row.names = FALSE)
  
  # 创建Excel文件保存详细的蛋白质列表
  # 格式：每个表型一个sheet，不同层级作为列（1/0表示）
  wb <- createWorkbook()
  
  # 为每个表型创建一个sheet
  for (pheno_name in names(phenotypes)) {
    # 收集该表型在所有层级出现的蛋白
    all_proteins_for_pheno <- data.frame()
    for (level in c("obs", "obs_mr", "obs_gsmr", "obs_union", "obs_mr_hyprcoloc", "obs_gsmr_hyprcoloc", "obs_union_hyprcoloc")) {
      key <- paste0(pheno_name, "_", level)
      if (key %in% names(results_list) && nrow(results_list[[key]]) > 0) {
        all_proteins_for_pheno <- rbind(all_proteins_for_pheno, results_list[[key]])
      }
    }
    
    if (nrow(all_proteins_for_pheno) > 0) {
      # 去重，获取所有蛋白
      all_proteins_for_pheno <- all_proteins_for_pheno %>%
        filter(!is.na(Entrez.Gene.Name) & Entrez.Gene.Name != "" & Entrez.Gene.Name != "0") %>%
        distinct(protein_id, .keep_all = TRUE)
      
      # 为每个层级添加1/0列
      for (level in c("obs", "obs_mr", "obs_gsmr", "obs_union", "obs_mr_hyprcoloc", "obs_gsmr_hyprcoloc", "obs_union_hyprcoloc")) {
        key <- paste0(pheno_name, "_", level)
        if (key %in% names(results_list)) {
          level_proteins <- results_list[[key]]$protein_id
          all_proteins_for_pheno[[level]] <- ifelse(all_proteins_for_pheno$protein_id %in% level_proteins, 1, 0)
        } else {
          all_proteins_for_pheno[[level]] <- 0
        }
      }
      
      # 添加sheet
      addWorksheet(wb, pheno_name)
      writeData(wb, pheno_name, all_proteins_for_pheno)
    }
  }
  
  # 创建all_combined sheet（所有表型合并去重）
  all_combined_df <- data.frame()
  for (level in c("obs", "obs_mr", "obs_gsmr", "obs_union", "obs_mr_hyprcoloc", "obs_gsmr_hyprcoloc", "obs_union_hyprcoloc")) {
    if (level %in% names(combined_proteins) && nrow(combined_proteins[[level]]) > 0) {
      all_combined_df <- rbind(all_combined_df, combined_proteins[[level]])
    }
  }
  
  if (nrow(all_combined_df) > 0) {
    all_combined_df <- all_combined_df %>%
      filter(!is.na(Entrez.Gene.Name) & Entrez.Gene.Name != "" & Entrez.Gene.Name != "0") %>%
      distinct(protein_id, .keep_all = TRUE)
    
    # 为每个层级添加1/0列
    for (level in c("obs", "obs_mr", "obs_gsmr", "obs_union", "obs_mr_hyprcoloc", "obs_gsmr_hyprcoloc", "obs_union_hyprcoloc")) {
      if (level %in% names(combined_proteins)) {
        level_proteins <- combined_proteins[[level]]$protein_id
        all_combined_df[[level]] <- ifelse(all_combined_df$protein_id %in% level_proteins, 1, 0)
      } else {
        all_combined_df[[level]] <- 0
      }
    }
    
    addWorksheet(wb, "all_combined")
    writeData(wb, "all_combined", all_combined_df)
  }
  
  saveWorkbook(wb, file.path(out_path, "phenotype_proteins_detailed.xlsx"), overwrite = TRUE)
  
  cat("模型", model, "处理完成\n")
}
