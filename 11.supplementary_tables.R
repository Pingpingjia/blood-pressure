#library(rio)
library(openxlsx)
library(dplyr)
# 路径前缀（请根据实际情况修改）
prefix <- "D:/OneDrive/工作/1.工作/blood pressure/"

## Keep only observed MR-style sheets (obs_mr / obs_gsmr / obs_union) and exclude hypr entries
sheet_content <- list(
  data.frame(),
  function() read.table(paste0(prefix, "model_1/obs_all_results.txt"), header = TRUE, sep = "\t", check.names = FALSE),
  function() read.table(paste0(prefix, "model_2/obs_all_results.txt"), header = TRUE, sep = "\t", check.names = FALSE),
  function() read.table(paste0(prefix, "model_3/obs_all_results.txt"), header = TRUE, sep = "\t", check.names = FALSE),
  data.frame(),
  function() read.table(paste0(prefix, "mr/mr_all_results.txt"), header = TRUE, sep = "\t", check.names = FALSE),
  function() read.table(paste0(prefix, "gsmr/gsmr_all_results.txt"), header = TRUE, sep = "\t", check.names = FALSE),
  function() read.table(paste0(prefix,  "hyprcoloc/hyprcoloc_all_significant_results.txt"), header = TRUE, sep = "\t", check.names = FALSE),
  function() read.table(paste0(prefix, "model_1/merged_obs_allmr_hyprcoloc.txt"), header = TRUE, sep = "\t", check.names = FALSE),
  function() read.table(paste0(prefix, "model_2/merged_obs_allmr_hyprcoloc.txt"), header = TRUE, sep = "\t", check.names = FALSE),
  function() read.table(paste0(prefix, "model_3/merged_obs_allmr_hyprcoloc.txt"), header = TRUE, sep = "\t", check.names = FALSE),
  data.frame(),
  function() read.table(paste0(prefix, "sg_model_1/obs_all_results.txt"), header = TRUE, sep = "\t", check.names = FALSE),
  function() read.table(paste0(prefix, "sg_model_2/obs_all_results.txt"), header = TRUE, sep = "\t", check.names = FALSE),
  function() read.table(paste0(prefix, "sg_model_3/obs_all_results.txt"), header = TRUE, sep = "\t", check.names = FALSE)
)


# 创建 workbook
wb <- createWorkbook()

# 每个sheet的表头
sheet_headers <- c(
  "Characteristics of the Study Population in Discovery and Validation Cohorts",
  "Regression Results of Model 1 for Blood Pressure Phenotypes and All Proteins",
  "Regression Results of Model 2 for Blood Pressure Phenotypes and All Proteins",
  "Regression Results of Model 3 for Blood Pressure Phenotypes and All Proteins",
  "GWAS Summary Data Included for Mendelian Randomization (MR) Analysis",
  "Two-Sample Mendelian Randomization(Two-sample MR) Results for Blood Pressure Phenotypes",
  "Generalized Summary Mendelian Randomization (GSMR) Results for Blood Pressure Phenotypes",
  "Colocalization Results between Proteins and blood pressure Phenotypes",
  "Integrated Results of Model 1, Two-sample MR, GSMR, and Colocalization for Blood Pressure Phenotypes",
  "Integrated Results of Model 2, Two-sample MR, GSMR, and Colocalization for Blood Pressure Phenotypes",
  "Integrated Results of Model 3, Two-sample MR, GSMR, and Colocalization for Blood Pressure Phenotypes",
  "The Detailed Information (location and drug) of the Identified Core Proteins",
  "Regression results of model 1 for the Blood Pressure Phenotypes and all the proteins in the validation cohort",
  "Regression results of model 2 for the Blood Pressure Phenotypes and all the proteins in the validation cohort",
  "Regression results of model 3 for the Blood Pressure Phenotypes and all the proteins in the validation cohort"
)

# ============================================================
# ✅ 先创建目录页 (Table of Contents)
# ============================================================
toc_df <- data.frame(
  Table = paste0("Supplementary Table S", seq_along(sheet_headers)),
  Description = sheet_headers,
  stringsAsFactors = FALSE
)

addWorksheet(wb, "Contents")
writeData(wb, sheet = "Contents", x = "List of Supplementary Tables", startRow = 1, colNames = FALSE)
writeData(wb, sheet = "Contents", x = toc_df, startRow = 3, colNames = TRUE)

# 样式设置
titleStyle <- createStyle(
  fontName = "Times New Roman",
  fontSize = 14,
  textDecoration = "bold",
  halign = "left"
)
headerStyle <- createStyle(
  fontName = "Times New Roman",
  fontSize = 12,
  textDecoration = "bold",
  border = "Bottom",
  halign = "center"
)
bodyStyle <- createStyle(fontName = "Times New Roman", fontSize = 12)

# 应用样式
addStyle(wb, "Contents", titleStyle, rows = 1, cols = 1, gridExpand = TRUE)
addStyle(wb, "Contents", headerStyle, rows = 3, cols = 1:2, gridExpand = TRUE)
addStyle(wb, "Contents", bodyStyle, rows = 4:(nrow(toc_df) + 3), cols = 1:2, gridExpand = TRUE)

# 自动调整列宽
setColWidths(wb, "Contents", cols = 1:2, widths = "auto")

# ============================================================
# ✅ 循环生成各 Supplementary Sheets
# ============================================================
for (i in seq_along(sheet_content)) {
  sheet_name <- paste0("Table S", i)
  addWorksheet(wb, sheet_name)
  
  dat <- if (is.function(sheet_content[[i]])) sheet_content[[i]]() else sheet_content[[i]]
  
  # 表头文本
  header_text <- paste0("Supplementary Table S", i, ". ", sheet_headers[i])
  writeData(wb, sheet = sheet_name, x = data.frame(header_text), startRow = 1, colNames = FALSE)
  
  # 表头样式
  titleStyle <- createStyle(
    fontName = "Times New Roman",
    fontSize = 12,
    textDecoration = "bold",
    halign = "left",
    valign = "center"
  )
  addStyle(wb, sheet = sheet_name, style = titleStyle, rows = 1, cols = 1, gridExpand = TRUE)
  
  if (nrow(dat) > 0 || ncol(dat) > 0) {
    writeData(wb, sheet = sheet_name, x = dat, startRow = 2, colNames = TRUE)
    
    headerStyle <- createStyle(
      fontName = "Times New Roman",
      fontSize = 12,
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "TopBottom"
    )
    addStyle(wb, sheet = sheet_name, style = headerStyle,
             rows = 2, cols = seq_len(ncol(dat)), gridExpand = TRUE)
    
    allStyle <- createStyle(fontName = "Times New Roman", fontSize = 12)
    addStyle(wb, sheet = sheet_name, style = allStyle,
             rows = 1:(nrow(dat) + 2), cols = seq_len(ncol(dat)), gridExpand = TRUE, stack = TRUE)
    
    # ✅ 底部边框
    bottomBorderStyle <- createStyle(border = "Bottom", borderColour = "black")
    last_row <- nrow(dat) + 2
    addStyle(
      wb, sheet = sheet_name, style = bottomBorderStyle,
      rows = last_row, cols = seq_len(ncol(dat)),
      gridExpand = TRUE, stack = TRUE
    )
    
  } else {
    allStyle <- createStyle(fontName = "Times New Roman", fontSize = 12)
    addStyle(wb, sheet = sheet_name, style = allStyle, rows = 1, cols = 1, gridExpand = TRUE)
  }
}

# ============================================================
# ✅ 保存文件
# ============================================================
# 保存
saveWorkbook(wb, paste0(prefix, "supplementary.xlsx"), overwrite = TRUE)

# ============================================================
# ✅ 读取并输出单个Excel文件
# ============================================================
cat("\n正在生成单独的补充表文件...\n")

# 读取刚保存的supplementary.xlsx
wb_full <- loadWorkbook(paste0(prefix, "supplementary.xlsx"))
sheet_names <- names(wb_full)

# 排除Contents页，只处理Table S开头的sheet
table_sheets <- sheet_names[grepl("^Table S", sheet_names)]

for (sheet_name in table_sheets) {
  # 创建新的workbook
  wb_single <- createWorkbook()
  
  # 添加sheet（使用相同的名称）
  addWorksheet(wb_single, sheet_name)
  
  # 读取原始数据
  sheet_data <- readWorkbook(paste0(prefix, "supplementary.xlsx"), 
                             sheet = sheet_name, 
                             colNames = FALSE, 
                             rowNames = FALSE)
  
  # 写入数据
  writeData(wb_single, sheet = sheet_name, x = sheet_data, 
            startRow = 1, colNames = FALSE)
  
  # 应用样式（复制原始样式）
  # 第一行：标题样式
  titleStyle <- createStyle(
    fontName = "Times New Roman",
    fontSize = 12,
    textDecoration = "bold",
    halign = "left",
    valign = "center"
  )
  addStyle(wb_single, sheet = sheet_name, style = titleStyle, 
           rows = 1, cols = 1, gridExpand = TRUE)
  
  # 如果有数据行（第2行是表头）
  if (nrow(sheet_data) > 1) {
    # 表头样式
    headerStyle <- createStyle(
      fontName = "Times New Roman",
      fontSize = 12,
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "TopBottom"
    )
    addStyle(wb_single, sheet = sheet_name, style = headerStyle,
             rows = 2, cols = seq_len(ncol(sheet_data)), gridExpand = TRUE)
    
    # 数据区域样式
    allStyle <- createStyle(fontName = "Times New Roman", fontSize = 12)
    addStyle(wb_single, sheet = sheet_name, style = allStyle,
             rows = 1:nrow(sheet_data), cols = seq_len(ncol(sheet_data)), 
             gridExpand = TRUE, stack = TRUE)
    
    # 底部边框
    bottomBorderStyle <- createStyle(border = "Bottom", borderColour = "black")
    last_row <- nrow(sheet_data)
    addStyle(wb_single, sheet = sheet_name, style = bottomBorderStyle,
             rows = last_row, cols = seq_len(ncol(sheet_data)),
             gridExpand = TRUE, stack = TRUE)
  }
  
  # 保存单独的文件，文件名格式：Supplementary Table S1.xlsx
  output_filename <- paste0(prefix, "Supplementary ", sheet_name, ".xlsx")
  saveWorkbook(wb_single, output_filename, overwrite = TRUE)
  
  cat(paste0("已生成: ", output_filename, "\n"))
}

cat("\n所有单独的补充表文件已生成完成！\n")

