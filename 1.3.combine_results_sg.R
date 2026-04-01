# Title: colocalization
# Author: Pingping Jia
# Date: 2025-06-05
# part1: define regions for clump
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
library(xlsx)
#===================================================================================================
#                        1. 合并所有观察性研究结果
#===================================================================================================
#mac
#data_dir <- "/Users/pingpingjia/Library/CloudStorage/OneDrive-TheChineseUniversityofHongKong/工作/1.工作/comorbidity/1.Proteomics/results"
#out_dir <- "/Users/pingpingjia/Library/CloudStorage/OneDrive-TheChineseUniversityofHongKong/工作/1.工作/comorbidity/results/combined_results"
out_dir <- "D:/OneDrive/工作/1.工作/blood pressure"
obs_data_dir <- "D:/OneDrive/工作/1.工作/comorbidity/1.Proteomics/results"
setwd("E:/results")
anno_info <- read.csv("E:/results/anno.info_updat.csv", header = TRUE, stringsAsFactors = FALSE)
for (model in c("sg_model_1", "sg_model_2", "sg_model_3")) {
    model_path <- file.path(data_dir, model)
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
      dplyr::distinct(traits, seqName, .keep_all = FALSE)
  
  # 输出所有结果
  write.table(obs_all_results, file = file.path(out_path, "obs_all_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # 输出显著结果
  write.table(obs_results, file = file.path(out_path, "obs_all_significant_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  # 构建宽格式：每个protein_id，按trait和group分别标记
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
