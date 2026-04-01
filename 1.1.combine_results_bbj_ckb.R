# =========================================================
# 批量安装和加载R包
#packages <- c("tidyr", "stringr","dplyr","stringr")
#installed <- packages %in% rownames(installed.packages())
#if (any(!installed)) {
#    install.packages(packages[!installed])
#}
#lapply(packages, library, character.only = TRUE)
library(dplyr)
library(stringr)
library(purrr)
library(UpSetR)
#===================================================================================================
#                        1. 合并所有观察性研究结果
#===================================================================================================
data_dir <- "E:/results/analysis_results"
out_dir <- "D:/OneDrive/工作/1.工作/blood pressure"
obs_data_dir <- "D:/OneDrive/工作/1.工作/comorbidity/1.Proteomics/results"
if (!dir.exists(out_dir)) {
    dir.create(out_dir)
}
anno_info <- read.csv("E:/results/anno.info_updat.csv", header = TRUE, stringsAsFactors = FALSE)
for (model in c("model_1", "model_2", "model_3")) {
    model_path <- file.path(obs_data_dir, model)
    out_path <- file.path(out_dir, model)
    if (!dir.exists(out_path)) {
        dir.create(out_path)
        }
    #folders <- folders[grepl("_results$", basename(folders))]
    # 读取并处理每个文件夹下的观察性研究结果
    obs_files <- c("sbp_result.csv", "dbp_result.csv","pp_result.csv","map_result.csv","hypertension_result.csv")
    obs_traits <- c("sbp", "dbp","pp","map","hypertension")
    blood_pressure_group = c("sbp", "dbp","pp","map","hypertension")
    # 读取所有结果（不过滤显著性）
    obs_all_results <- lapply(seq_along(obs_files), function(i) {
    file_path <- file.path(model_path, obs_files[i])
    if (file.exists(file_path)) {
        df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
        df <- df[, c("seqName", "Entrez.Gene.Name", "Estimate", "se", "z_value", "P_value", "adjust_P")]
        df$traits <- obs_traits[i]
        return(df)
    } else {
        return(NULL)
    }
    })
  
  obs_all_results <- do.call(rbind, obs_all_results)
  obs_all_results <- na.omit(obs_all_results)
  
  # 只保留显著结果
  obs_results <- obs_all_results[obs_all_results$adjust_P <= 0.05, ]
  
  # 按照traits去重
  obs_results_unique <- obs_results %>%
      distinct(traits, seqName, .keep_all = FALSE)
  
  # 输出所有结果
  write.table(obs_all_results, file = file.path(out_path, "obs_all_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # 输出显著结果
  write.table(obs_results, file = file.path(out_path, "obs_all_significant_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  # 构建宽格式：每个protein_id，按trait和group分别标记
  obs_wide <- obs_results_unique %>%
    mutate(
      sbp = as.integer(traits == "sbp"),
      dbp = as.integer(traits == "dbp"),
      pp = as.integer(traits == "pp"),
      map = as.integer(traits == "map"),
      hypertension = as.integer(traits == "hypertension"),
      blood_pressure = as.integer(traits %in% blood_pressure_group)
    ) %>%
    group_by(seqName) %>%
    summarise(
      sbp = max(sbp),
      dbp = max(dbp),
      pp = max(pp),
      map = max(map),
      hypertension = max(hypertension),
      blood_pressure = max(blood_pressure)
    ) %>%
    rename(protein_id = seqName)
  obs_wide <- left_join(obs_wide,anno_info[,c("AptName","Entrez.Gene.Name")],by = c("protein_id" = "AptName"))
  write.table(obs_wide, file = file.path(out_path, "obs_protein_trait_group_status.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  # 按照traits统计显著protein_id数量（不重复计数）
  obs_trait_counts <- obs_results_unique %>%
      group_by(traits) %>%
      summarise(significant_protein_count = n_distinct(seqName), .groups = "drop")
  
  # 按照group统计显著protein_id数量（不重复计数）
  obs_group_counts <- obs_results_unique %>%
      mutate(group = "blood_pressure") %>%
      group_by(group) %>%
      summarise(significant_protein_count = n_distinct(seqName), .groups = "drop")

  # 合并traits和group统计结果
  obs_trait_counts <- obs_trait_counts %>% mutate(type = "trait", name = traits) %>% select(type, name, significant_protein_count)
  obs_group_counts <- obs_group_counts %>% mutate(type = "group", name = group) %>% select(type, name, significant_protein_count)
  obs_counts_combined <- bind_rows(obs_trait_counts, obs_group_counts)
  
  
  # 输出合并统计结果
  write.table(obs_counts_combined, file = file.path(out_path, "obs_significant_protein_counts_by_trait_and_group.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}
#===================================================================================================
#                       2. 合并所有mr结果
#===================================================================================================
folders <- list.dirs(path = file.path(data_dir, "mr_results"), full.names = TRUE, recursive = FALSE)
folders <- folders[grepl("_results$", basename(folders))]
# 只筛选血压相关的结果
folders <- folders[grepl("(bbj|ckb|ukb)_(sbp|dbp|pp|map|hypertension)", basename(folders))]
out_path <- file.path(out_dir, "mr")
if (!dir.exists(out_path)) {
    dir.create(out_path)
}   
# trait分组定义
bp_traits <- c("bbj_sbp", "bbj_dbp","bbj_pp","bbj_map","ckb_sbp", "ckb_dbp","ckb_pp","ckb_map","ukb_hypertension","combined_sbp","combined_dbp","combined_pp","combined_map")
# 过滤掉不需要的文件夹
# 定义需要的列,包括mr和多效性检验的结果
cols_needed <- c("id.exposure", "id.outcome", "outcome", "exposure", "method", "nsnp", "b", "se", "pval", "pleiotropy_se","pleiotropy_pval","protein_id")

# 读取并处理每个文件夹下的all_mr_results.tsv
# 先检查并合并每个folder下的mr_results.txt为all_mr_results.tsv, 同时合并pleiotropy.txt
for (folder in folders) {
    outcome_name <- str_extract(basename(folder), "bbj_[^_]+|ckb_[^_]+")
    subfolders <- list.dirs(folder, full.names = TRUE, recursive = FALSE)

    # 合并 mr_results.txt
    mr_outfile <- file.path(folder, "all_mr_results.tsv")
    if (!file.exists(mr_outfile)) {
        mr_files <- file.path(subfolders, "mr_results.txt")
        mr_existing <- mr_files[file.exists(mr_files)]
        if (length(mr_existing) > 0) {
            mr_dfs <- lapply(mr_existing, function(f) read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
            mr_combined <- dplyr::bind_rows(mr_dfs)
            write.table(mr_combined, file = mr_outfile, sep = "\t", row.names = FALSE, quote = FALSE)
            cat(sprintf("Outcome: %s, merged %d mr_results.txt files\n", outcome_name, length(mr_existing)))
        } else {
            cat(sprintf("Outcome: %s, no mr_results.txt files found\n", outcome_name))
        }
    }

    # 合并 pleiotropy.txt
    pleio_outfile <- file.path(folder, "all_pleiotropy_results.tsv")
    if (!file.exists(pleio_outfile)) {
        pleio_files <- file.path(subfolders, "pleiotropy.txt")
        if (length(pleio_files) > 0) {
            pleio_dfs <- lapply(pleio_files, function(f) {
                df <- tryCatch(read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
                                error = function(e) NULL)
            })
            pleio_dfs <- pleio_dfs[!vapply(pleio_dfs, is.null, logical(1))]
            if (length(pleio_dfs) > 0) {
                pleio_combined <- dplyr::bind_rows(pleio_dfs) %>% dplyr::distinct()
                # 只保留需要的列（若存在）
                needed_cols <- c("id.exposure","id.outcome","outcome","exposure","egger_intercept","se","pval")
                existing_needed <- intersect(needed_cols, colnames(pleio_combined))
                pleio_combined <- pleio_combined[, existing_needed, drop = FALSE]
                write.table(pleio_combined, file = pleio_outfile, sep = "\t", row.names = FALSE, quote = FALSE)
                cat(sprintf("Outcome: %s, merged %d pleiotropy files\n", outcome_name, length(pleio_files)))
            } else {
                cat(sprintf("Outcome: %s, pleiotropy files unreadable\n", outcome_name))
            }
        } else {
            cat(sprintf("Outcome: %s, no pleiotropy files found\n", outcome_name))
        }
    }
}

# 再读取all_mr_results.tsv并处理
mr_all_results <- lapply(folders, function(folder) {
    file_path <- file.path(folder, "all_mr_results.tsv")
    file_path_pleio <- file.path(folder, "all_pleiotropy_results.tsv") 
    outcome_name <- str_extract(basename(folder), "(bbj|ckb|ukb)_(sbp|dbp|map|pp|hypertension)")
    if (file.exists(file_path)) {
        df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        pleio_df <- read.table(file_path_pleio, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        pleio_df <- pleio_df[, intersect(c("protein_id", "id.exposure", "id.outcome", "egger_intercept", "se", "pval"), colnames(pleio_df)), drop = FALSE]
        pleio_df <- pleio_df %>% rename(pleiotropy_se = se, pleiotropy_pval = pval)
        # left_join 按 protein_id
        df_all <- left_join(df, pleio_df, by = c("id.exposure", "id.outcome"))
        #df$adjust_P <- p.adjust(df$pval,method = "BH")
        # 合并 IVW (pleiotropy_pval >= 0.05) 和 Egger (pleiotropy_pval < 0.05) 结果为 IVW_Egger
        ivw<- df_all[df_all$method == "Inverse variance weighted", ]
        ivw_df <- df_all[df_all$method == "Inverse variance weighted" & df_all$pleiotropy_pval >= 0.05, ]
        ivw_df <- ivw_df[!is.na(ivw_df$pleiotropy_pval),]
        egger_df <- df_all[df_all$method == "MR Egger" & df_all$pleiotropy_pval < 0.05, ]
        egger_df <- egger_df[!is.na(egger_df$pleiotropy_pval),]
        ivw_egger_df <- bind_rows(ivw_df, egger_df) %>% distinct()
        ivw_egger_df$method <- "IVW_Egger"
        # 合并所有结果并去重
        df_all2 <- bind_rows(df_all, ivw_egger_df) %>% distinct()
        # 保留所有结果和需要的列
        df_all_complete <- df_all2 %>% select(all_of(cols_needed))
        # 只保留显著且需要的列
        df_all3 <- df_all2 %>% filter(pval < 0.05) %>% select(all_of(cols_needed))
        df_all_complete$traits <- outcome_name
        df_all3$traits <- outcome_name
        return(list(all_results = df_all_complete, significant_results = df_all3))
    } else {
        return(list(all_results = NULL, significant_results = NULL))
    }
})

# 合并所有结果和显著结果
mr_all_complete_results <- lapply(mr_all_results, function(x) x$all_results)
mr_all_complete_results <- do.call(rbind, mr_all_complete_results)
mr_all_complete_results <- na.omit(mr_all_complete_results)
#=================================================
#     sbp,dbp,pp,map来自于bbj和ckb,需要合并
#=================================================
# 从MR数据中提取bbj和ckb数据
combine_bp_data_mr <- function(mr_data, trait_pair) {
  # trait_pair: list(bbj="bbj_sbp", ckb="ckb_sbp", combined="combined_sbp")
  cat(sprintf("\n合并 %s 和 %s MR数据为 %s\n", trait_pair$bbj, trait_pair$ckb, trait_pair$combined))
  #mr_data = mr_raw_data
  #trait_pair <- list(bbj="bbj_dbp", ckb="ckb_dbp", combined="combined_dbp")
  #protein = all_proteins[1]
  # 从MR数据中提取bbj和ckb数据
  bbj_data <- mr_data[mr_data$traits == trait_pair$bbj, ]
  ckb_data <- mr_data[mr_data$traits == trait_pair$ckb, ]
  
  # 获取所有蛋白质ID
  all_proteins <- unique(c(bbj_data$protein_id, ckb_data$protein_id))
  
  combined_results <- list()
  
  for (protein in all_proteins) {
    bbj_row <- bbj_data[bbj_data$protein_id == protein, ]
    ckb_row <- ckb_data[ckb_data$protein_id == protein, ]
    
    # 判断MR显著性（p < 0.05）
    bbj_sig <- nrow(bbj_row) > 0 && !is.na(bbj_row$pval[1]) && bbj_row$pval[1] < 0.05
    ckb_sig <- nrow(ckb_row) > 0 && !is.na(ckb_row$pval[1]) && ckb_row$pval[1] < 0.05
    
    # 选择规则
    selected_row <- NULL
    selection_reason <- ""
    
    if (bbj_sig && nrow(bbj_row) > 0) {
      # BBJ显著且有数据，直接使用BBJ
      selected_row <- bbj_row
      selection_reason <- ifelse(ckb_sig, "both_sig_prioritize_bbj", "bbj_significant")
    } else if (nrow(ckb_row) > 0) {
      # BBJ不显著或没有数据，使用CKB
      selected_row <- ckb_row
      selection_reason <- ifelse(ckb_sig, "ckb_significant_only", "neither_sig_use_ckb")
    } else if (nrow(bbj_row) > 0) {
      # CKB没有数据，回退到BBJ
      selected_row <- bbj_row
      selection_reason <- "ckb_no_data_fallback_to_bbj"
    }
    
    if (exists("selected_row") && nrow(selected_row) > 0) {
      selected_row$traits <- trait_pair$combined
      selected_row$selection_reason <- selection_reason
      combined_results[[protein]] <- selected_row
    }
  }
  
  combined_df <- do.call(rbind, combined_results)
  cat(sprintf("MR合并完成：%d个蛋白质\n", length(combined_results)))
  
  # 显示选择统计
  if (!is.null(combined_df) && nrow(combined_df) > 0) {
    selection_summary <- table(combined_df$selection_reason)
    cat("MR选择原因统计：\n")
    print(selection_summary)
  }
  
  return(combined_df)
}
mr_raw_data <- mr_all_complete_results
mr_raw_data <- mr_raw_data[mr_raw_data$method == "IVW_Egger", ]
# 合并sbp和dbp数据
combined_sbp_mr <- combine_bp_data_mr(
  mr_raw_data,list(bbj = "bbj_sbp", ckb = "ckb_sbp", combined = "combined_sbp")
)
combined_dbp_mr <- combine_bp_data_mr(
  mr_raw_data,list(bbj = "bbj_dbp", ckb = "ckb_dbp", combined = "combined_dbp")
)
combined_pp_mr <- combine_bp_data_mr(
  mr_raw_data,list(bbj = "bbj_pp", ckb = "ckb_pp", combined = "combined_pp")
)
combined_map_mr <- combine_bp_data_mr(
  mr_raw_data,list(bbj = "bbj_map", ckb = "ckb_map", combined = "combined_map")
)
mr_all_complete_results2 <- rbind(mr_all_complete_results,combined_sbp_mr[,-c(14)],combined_dbp_mr[,-c(14)],combined_pp_mr[,-c(14)],combined_map_mr[,-c(14)])

#=================================================
#     显著的结果
#=================================================
mr_significant_results <- mr_all_complete_results2[mr_all_complete_results2$pval< 0.05,]
mr_all_results <- na.omit(mr_significant_results)
# 按照method, protein_id, traits去重
mr_all_results_unique <- mr_all_results %>% distinct(method, protein_id, traits, .keep_all = FALSE)

# 按照不同的method，统计每个traits显著的protein_id数量（不重复计数）
protein_counts <- mr_all_results_unique %>%
  group_by(method, traits) %>%
  summarise(significant_protein_count = n_distinct(protein_id), .groups = "drop") %>%
  mutate(type = "trait", name = traits) %>%
  select(type, method, name, significant_protein_count)

group_counts <- mr_all_results_unique %>%
  mutate(group = "blood_pressure") %>%
  group_by(method, group) %>%
  summarise(significant_protein_count = n_distinct(protein_id), .groups = "drop") %>%
  mutate(type = "group", name = group) %>%
  select(type, method, name, significant_protein_count)

# 合并trait和group统计结果
mr_counts_combined <- bind_rows(protein_counts, group_counts)

# 输出合并统计结果
write.table(mr_counts_combined, file = file.path(out_path, "mr_significant_protein_counts_by_method_trait_and_group.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# 输出所有MR结果
write.table(mr_all_complete_results2, file = file.path(out_path, "mr_all_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# 输出所有显著结果
write.table(mr_all_results, file = file.path(out_path, "mr_all_significant_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
#===================================================================================================
#                        分组查看结果
#===================================================================================================
# 构建宽格式：每个protein_id、method，不同trait和group的0/1结果
all_trait_names <- unique(bp_traits)
group_names <- c("blood_pressure")

mr_protein_group_status <- mr_all_results_unique %>%
  mutate(value = 1) %>%
  tidyr::pivot_wider(
    id_cols = c(protein_id, method),
    names_from = traits,
    values_from = value,
    values_fill = 0
  ) %>%
  mutate(
    blood_pressure = as.integer(rowSums(select(., all_of(bp_traits))) > 0)
  )

# 输出seqID, method, 不同trait和group的0/1结果
write.table(mr_protein_group_status, file = file.path(out_path, "mr_protein_trait_group_status.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
#===================================================================================================
#                       2.2 合并所有gsmr结果
#===================================================================================================
gsmr_folders <- list.dirs(path = file.path(data_dir, "gsmr_0905"), full.names = TRUE, recursive = FALSE)
# 检查并合并每个folder下的gsmr_main_result.txt为all_gsmr_main_result.txt
out_path <- file.path(out_dir, "gsmr")
if (!dir.exists(out_path)) {
    dir.create(out_path)
}   
for (folder in gsmr_folders) {
    file_path <- file.path(folder, "all_gsmr_main_result.txt")
    if (!file.exists(file_path)) {
        print(folder)
        subfolders <- list.dirs(folder, full.names = TRUE, recursive = FALSE)
        gsmr_files <- file.path(subfolders, "gsmr_main_result.txt")
        existing_files <- gsmr_files[file.exists(gsmr_files)]
        if (length(existing_files) > 0) {
            dfs <- lapply(existing_files, function(f) read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
            combined_df <- do.call(rbind, dfs)
            write.table(combined_df, file = file.path(folder, "all_gsmr_main_result.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
            cat(sprintf("Merged %d gsmr_main_result.txt files in %s\n", length(existing_files), basename(folder)))
        } else {
            cat(sprintf("No gsmr_main_result.txt files found in %s\n", basename(folder)))
        }
    }
}
# 读取并处理每个文件夹下的all_gsmr_main_result.txt
gsmr_results <- lapply(gsmr_folders, function(folder) {
    file_path <- file.path(folder, "all_gsmr_main_result.txt")
    outcome_name <- str_extract(basename(folder), "(bbj|ckb|ukb)_(sbp|dbp|map|pp|hypertension)")
    if (file.exists(file_path)) {
        df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        # 保存所有结果
        df_all <- df[, c("protein_id", "bxy", "bxy_se", "bxy_pval")]
        df_all$method <- "gsmr"
        df_all$traits <- outcome_name
        # 保存显著结果
        df_sig <- df[df$bxy_pval < 0.05, c("protein_id", "bxy", "bxy_se", "bxy_pval")]
        df_sig$method <- "gsmr"
        df_sig$traits <- outcome_name
        return(list(all_results = df_all, significant_results = df_sig))
    } else {
        return(list(all_results = NULL, significant_results = NULL))
    }
})

# 合并所有结果和显著结果
gsmr_all_complete_results <- lapply(gsmr_results, function(x) x$all_results)
gsmr_all_complete_results <- do.call(rbind, gsmr_all_complete_results)
gsmr_all_complete_results <- na.omit(gsmr_all_complete_results)
#=================================================
#     sbp和dbp来自于bbj和ckb,需要合并
#=================================================
# 从GSMR数据中提取bbj和ckb数据
combine_bp_data_gsmr <- function(gsmr_data, trait_pair) {
  # trait_pair: list(obs="sbp", bbj="bbj_sbp", ckb="ckb_sbp", combined="combined_sbp")
  cat(sprintf("\n合并 %s 和 %s GSMR数据为 %s\n", trait_pair$bbj, trait_pair$ckb, trait_pair$combined))
  # seq.10036.201
  # 从GSMR数据中提取bbj和ckb数据
  bbj_data <- gsmr_data[gsmr_data$traits == trait_pair$bbj, ]
  ckb_data <- gsmr_data[gsmr_data$traits == trait_pair$ckb, ]
  
  # 获取所有蛋白质ID
  all_proteins <- unique(c(bbj_data$protein_id, ckb_data$protein_id))
  
  combined_results <- list()
  
  for (protein in all_proteins) {
    bbj_row <- bbj_data[bbj_data$protein_id == protein, ]
    ckb_row <- ckb_data[ckb_data$protein_id == protein, ]
    
    # 判断GSMR显著性（p < 0.05），处理NA值
    bbj_sig <- nrow(bbj_row) > 0 && !is.na(bbj_row$bxy_pval[1]) && bbj_row$bxy_pval[1] < 0.05
    ckb_sig <- nrow(ckb_row) > 0 && !is.na(ckb_row$bxy_pval[1]) && ckb_row$bxy_pval[1] < 0.05
    
    # 将NA转换为FALSE
    bbj_sig <- ifelse(is.na(bbj_sig), FALSE, bbj_sig)
    ckb_sig <- ifelse(is.na(ckb_sig), FALSE, ckb_sig)
    
    # 选择规则
    selected_row <- NULL
    selection_reason <- ""
    if (bbj_sig && nrow(bbj_row) > 0) {
      # BBJ显著且有数据，直接使用BBJ
      selected_row <- bbj_row
      selection_reason <- ifelse(ckb_sig, "both_sig_prioritize_bbj", "bbj_significant")
    } else if (nrow(ckb_row) > 0) {
      # BBJ不显著或没有数据，使用CKB
      selected_row <- ckb_row
      selection_reason <- ifelse(ckb_sig, "ckb_significant_only", "neither_sig_use_ckb")
    } else if (nrow(bbj_row) > 0) {
      # CKB没有数据，回退到BBJ
      selected_row <- bbj_row
      selection_reason <- "ckb_no_data_fallback_to_bbj"
    }
    if (exists("selected_row") && nrow(selected_row) > 0) {
      selected_row$traits <- trait_pair$combined
      selected_row$selection_reason <- selection_reason
      combined_results[[protein]] <- selected_row
    }
  }
  
  combined_df <- do.call(rbind, combined_results)
  cat(sprintf("GSMR合并完成：%d个蛋白质\n", length(combined_results)))
  
  # 显示选择统计
  if (!is.null(combined_df) && nrow(combined_df) > 0) {
    selection_summary <- table(combined_df$selection_reason)
    cat("GSMR选择原因统计：\n")
    print(selection_summary)
  }
  
  return(combined_df)
}
gsmr_raw_data <- gsmr_all_complete_results
# 合并sbp和dbp数据
combined_sbp_gsmr <- combine_bp_data_gsmr(
  gsmr_raw_data, list(bbj = "bbj_sbp", ckb = "ckb_sbp", combined = "combined_sbp")
)
combined_dbp_gsmr <- combine_bp_data_gsmr(
  gsmr_raw_data,list(bbj = "bbj_dbp", ckb = "ckb_dbp", combined = "combined_dbp")
)
combined_pp_gsmr <- combine_bp_data_gsmr(
  gsmr_raw_data,list(bbj = "bbj_pp", ckb = "ckb_pp", combined = "combined_pp")
)
combined_map_gsmr <- combine_bp_data_gsmr(
  gsmr_raw_data,list(bbj = "bbj_map", ckb = "ckb_map", combined = "combined_map")
) 

# 更新raw_gsmr_data,增加combined的行
gsmr_all_complete_results2 <- rbind(gsmr_all_complete_results,combined_sbp_gsmr[,-c(7)],combined_dbp_gsmr[,-c(7)],combined_pp_gsmr[,-c(7)],combined_map_gsmr[,-c(7)])
#=================================================
#     显著的结果
#=================================================
gsmr_significant_results <- gsmr_all_complete_results2[gsmr_all_complete_results2$bxy_pval < 0.05,]
gsmr_all_results <- na.omit(gsmr_significant_results)

# 按照protein_id, traits去重
gsmr_all_results_unique <- gsmr_all_results %>% distinct(protein_id, traits, .keep_all = FALSE)
gsmr_protein_group_status <- gsmr_all_results_unique %>%
  mutate(value = 1) %>%
  tidyr::pivot_wider(
    id_cols = protein_id,
    names_from = traits,
    values_from = value,
    values_fill = 0
  ) %>%
  mutate(
    blood_pressure = as.integer(rowSums(select(., all_of(bp_traits))) > 0)
  )

# trait统计
gsmr_trait_counts <- gsmr_all_results_unique %>%
  group_by(traits) %>%
  summarise(significant_protein_count = n_distinct(protein_id), .groups = "drop") %>%
  mutate(type = "trait", method = "gsmr", name = traits) %>%
  select(type, method, name, significant_protein_count)

# group统计
gsmr_group_counts <- gsmr_all_results_unique %>%
  mutate(group = "blood_pressure") %>%
  group_by(group) %>%
  summarise(significant_protein_count = n_distinct(protein_id), .groups = "drop") %>%
  mutate(type = "group", method = "gsmr", name = group) %>%
  select(type, method, name, significant_protein_count)

gsmr_counts_combined <- bind_rows(gsmr_trait_counts, gsmr_group_counts)

# 合并到mr_counts_combined
mr_counts_combined <- bind_rows(mr_counts_combined, gsmr_counts_combined)

# 输出所有gsmr结果
write.table(gsmr_all_complete_results2, file = file.path(out_path, "gsmr_all_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# 输出gsmr显著结果
write.table(gsmr_all_results, file = file.path(out_path, "gsmr_all_significant_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


# 构建宽格式：每个protein_id，不同trait和group的0/1结果（gsmr）
gsmr_trait_names <- unique(bp_traits)

gsmr_protein_group_status <- gsmr_all_results_unique %>%
  mutate(value = 1) %>%
  tidyr::pivot_wider(
    id_cols = protein_id,
    names_from = traits,
    values_from = value,
    values_fill = 0
  ) %>%
  mutate(
    blood_pressure = as.integer(rowSums(select(., all_of(bp_traits))) > 0)
  )

write.table(gsmr_protein_group_status, file = file.path(out_path, "gsmr_significant_protein_status.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# 输出gsmr统计结果
write.table(gsmr_counts_combined, file = file.path(out_path, "gsmr_significant_protein_counts_by_trait_and_group.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
#===================================================================================================
#                       3. 合并所有hyprcoloc结果
#===================================================================================================
# 合并所有hyprcoloc结果
# trait分组定义（仅保留血压相关）
bp_traits <- c("ckb_sbp", "ckb_dbp", "ckb_pp", "ckb_map","ukb_hypertension")
all_traits <- bp_traits
# hyprcoloc结果主目录
outcome_folders <- list.dirs(file.path(data_dir, "hyprcoloc_0905"), full.names = TRUE, recursive = FALSE)
outcome_folders <- outcome_folders[grepl("(ckb|ukb)_(sbp|dbp|pp|map|hypertension)", basename(outcome_folders))]
out_path <- file.path(out_dir, "hyprcoloc")

if (!dir.exists(out_path)) {
    dir.create(out_path)
}   

# 合并所有结果
hyprcoloc_results_raw <- map_dfr(outcome_folders, function(outcome_folder) {
    outcome_name <- str_extract(basename(outcome_folder), "(ckb|ukb)_(sbp|dbp|map|pp|hypertension)")
    subfolders <- list.dirs(outcome_folder, full.names = TRUE, recursive = FALSE)
    map_dfr(subfolders, function(subfolder) {
        file_path <- file.path(subfolder, "all_significant_region_results.txt")
        if (file.exists(file_path)) {
            tmp_df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
            # 只保留前14列
            df <- tmp_df %>% select(1:14)
            protein_seq_name <- str_extract(basename(subfolder), "seq.*")
            df$outcome <- outcome_name
            df$protein_seq_name <- protein_seq_name
            return(df)
        } else {
            return(NULL)
        } 
    })
})
# 去掉glgc的结果
hyprcoloc_results <- hyprcoloc_results_raw
# 输出所有hyprcoloc结果（hyprcoloc结果本身就是显著的共定位结果）
hyprcoloc_results2 <- hyprcoloc_results %>% left_join(anno_info[,c("AptName","Entrez.Gene.Name")], by = c("protein_seq_name" = "AptName")) %>% distinct()

write.table(hyprcoloc_results2, file = file.path(out_path, "hyprcoloc_all_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(hyprcoloc_results2, file = file.path(out_path, "hyprcoloc_all_significant_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# 只保留蛋白和trait信息
hyprcoloc_simple <- hyprcoloc_results2 %>%
    select(protein_seq_name, outcome,Entrez.Gene.Name) %>%
    distinct()
#去除protein_seq_name空白的
hyprcoloc_simple <- hyprcoloc_simple[!is.na(hyprcoloc_simple$protein_seq_name),]
# 构建宽格式矩阵
hyprcoloc_wide <- hyprcoloc_simple %>%
    mutate(value = 1) %>%
    tidyr::pivot_wider(names_from = outcome, values_from = value, values_fill = 0)
# 增加分组标记
hyprcoloc_wide <- hyprcoloc_wide %>%
    mutate(
    blood_pressure = as.integer(rowSums(select(., all_of(bp_traits))) > 0)
    )

# 输出宽格式结果
write.table(hyprcoloc_wide, file = file.path(out_path, "hyprcoloc_protein_trait_group_status.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# trait统计
hyprcoloc_trait_counts <- hyprcoloc_wide %>%
    summarise(across(all_of(all_traits), ~sum(.)), .groups = "drop") %>%
    tidyr::pivot_longer(cols = everything(), names_to = "name", values_to = "colocalized_protein_count") %>%
    mutate(type = "trait") %>%
    select(type, name, colocalized_protein_count)

# group统计
hyprcoloc_group_counts <- hyprcoloc_wide %>%
  summarise(blood_pressure = sum(blood_pressure)) %>%
    tidyr::pivot_longer(cols = everything(), names_to = "name", values_to = "colocalized_protein_count") %>%
    mutate(type = "group") %>%
    select(type, name, colocalized_protein_count)

# 合并统计结果
hyprcoloc_counts_combined <- bind_rows(
  hyprcoloc_trait_counts,
  hyprcoloc_group_counts
)

# 输出统计结果
write.table(hyprcoloc_counts_combined, file = file.path(out_path, "hyprcoloc_colocalized_protein_counts_by_trait_and_group.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

