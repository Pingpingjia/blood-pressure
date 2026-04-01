library(dplyr)

# 读取数据
out_dir <- "D:/OneDrive/工作/1.工作/blood pressure"
out_path <- file.path(out_dir, "model_1")

merged_all <- read.table(file.path(out_path, "merged_obs_allmr_hyprcoloc.txt"), 
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("总行数:", nrow(merged_all), "\n")
cat("列名:\n")
print(colnames(merged_all))

cat("\n列的类型:\n")
print(sapply(merged_all, class))

# 将NA替换为0
merged_all[is.na(merged_all)] <- 0

# 检查sbp相关列
cat("\nsbp.obs列中值为1的行数:", sum(merged_all$sbp.obs == 1, na.rm=TRUE), "\n")
cat("ckb_sbp.mr列中值为1的行数:", sum(merged_all$ckb_sbp.mr == 1, na.rm=TRUE), "\n")
cat("ckb_sbp.gsmr列中值为1的行数:", sum(merged_all$ckb_sbp.gsmr == 1, na.rm=TRUE), "\n")
cat("ckb_sbp.hyprcoloc列中值为1的行数:", sum(merged_all$ckb_sbp.hyprcoloc == 1, na.rm=TRUE), "\n")

cat("\nsbp.obs=1 & ckb_sbp.mr=1的行数:", sum(merged_all$sbp.obs==1 & merged_all$ckb_sbp.mr==1, na.rm=TRUE), "\n")

# 查看前几行
cat("\n前5行的sbp相关列:\n")
print(merged_all[1:5, c("protein_id", "sbp.obs", "ckb_sbp.mr", "ckb_sbp.gsmr", "ckb_sbp.hyprcoloc")])
