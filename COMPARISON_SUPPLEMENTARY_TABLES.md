# 补充表脚本对比详情

## 文件对比

### 关键代码差异

#### 1. 数据路径 (Prefix)

**10.supplementary_tables.R:**
```r
prefix <- "D:/OneDrive/工作/1.工作/comorbidity/results/combined_results_11/"
```

**11.supplementary_tables.R:**
```r
prefix <- "D:/OneDrive/工作/1.工作/blood pressure/"
```

---

#### 2. 表型范围说明

**10.supplementary_tables.R:**
```r
## 包含所有表型的观察性研究和MR分析结果
## 表型组包括：
## - 血压: sbp, dbp, pp, map, hypertension
## - 血脂: ldlc, hdlc, tg
## - 血糖: glu, diabetes, fbg
## - 肥胖: bmi, wc, whr
```

**11.supplementary_tables.R:**
```r
## 包含血压相关表型的观察性研究和MR分析结果
## 血压表型：sbp, dbp, pp, map, hypertension
```

---

#### 3. Hyprcoloc数据文件

**10.supplementary_tables.R:**
```r
function() read.table(paste0(prefix, "outcome_hyprcoloc_processed/combined_hyprcoloc_results_3.txt"), 
                      header = TRUE, sep = "\t", check.names = FALSE)
```

**11.supplementary_tables.R:**
```r
function() read.table(paste0(prefix, "hyprcoloc/hyprcoloc_all_significant_results.txt"), 
                      header = TRUE, sep = "\t", check.names = FALSE)
```

---

#### 4. 表头描述

**10.supplementary_tables.R:**
```r
sheet_headers <- c(
  "Characteristics of the Study Population in Discovery and Validation Cohorts",
  "Regression Results of Model 1 for All Phenotypes (Blood Pressure, Lipids, Glucose, Obesity) and All Proteins",
  "Regression Results of Model 2 for All Phenotypes (Blood Pressure, Lipids, Glucose, Obesity) and All Proteins",
  ...
)
```

**11.supplementary_tables.R:**
```r
sheet_headers <- c(
  "Characteristics of the Study Population in Discovery and Validation Cohorts",
  "Regression Results of Model 1 for Blood Pressure Phenotypes and All Proteins",
  "Regression Results of Model 2 for Blood Pressure Phenotypes and All Proteins",
  ...
)
```

---

#### 5. 输出文件名

**10.supplementary_tables.R:**
```r
output_file <- paste0(prefix, "supplementary_all_phenotypes.xlsx")
```

**11.supplementary_tables.R:**
```r
output_file <- paste0(prefix, "supplementary_blood_pressure.xlsx")
```

---

#### 6. 脚本结束信息

**10.supplementary_tables.R:**
```r
cat("\n脚本执行完成！\n")
cat("=================================================\n")
cat("此脚本处理所有表型组：\n")
cat("  - 血压 (Blood Pressure): sbp, dbp, pp, map, hypertension\n")
cat("  - 血脂 (Lipids): ldlc, hdlc, tg\n")
cat("  - 血糖 (Glucose): glu, diabetes, fbg\n")
cat("  - 肥胖 (Obesity): bmi, wc, whr\n")
cat("=================================================\n")
```

**11.supplementary_tables.R:**
```r
cat("\n脚本执行完成！\n")
cat("=================================================\n")
cat("此脚本专门处理血压相关表型：\n")
cat("  - SBP (收缩压)\n")
cat("  - DBP (舒张压)\n")
cat("  - PP (脉压)\n")
cat("  - MAP (平均动脉压)\n")
cat("  - Hypertension (高血压)\n")
cat("=================================================\n")
```

---

## 代码结构对比

### 共同特性

两个脚本共享以下结构和特性：

1. **相同的库依赖**
   ```r
   library(openxlsx)
   library(dplyr)
   ```

2. **相同的数据读取结构**
   - 15个sheet内容
   - 相同的model_1/2/3结构
   - 相同的MR和GSMR数据结构
   - 相同的验证队列数据（sg_model_1/2/3）

3. **相同的Excel生成逻辑**
   - 创建目录页（Contents）
   - 15个数据表（Table S1-S15）
   - 统一的样式设置（Times New Roman, 12pt）
   - 生成单独的Excel文件

4. **相同的错误处理机制**
   - tryCatch包装数据读取
   - 数据文件缺失时显示警告
   - 继续执行不中断

### 关键差异

| 特性 | 10.supplementary_tables.R | 11.supplementary_tables.R |
|------|---------------------------|---------------------------|
| **研究范围** | 多病共病（Multi-morbidity） | 血压专项（BP-specific） |
| **表型数量** | ~15个（4组） | 5个（仅血压） |
| **数据目录** | comorbidity/results/combined_results_11/ | blood pressure/ |
| **Hyprcoloc路径** | outcome_hyprcoloc_processed/combined_hyprcoloc_results_3.txt | hyprcoloc/hyprcoloc_all_significant_results.txt |
| **Hyprcoloc内容** | ≥3个trait的组合 | 血压相关所有显著结果 |
| **输出文件** | supplementary_all_phenotypes.xlsx | supplementary_blood_pressure.xlsx |
| **表格标题** | "All Phenotypes (Blood Pressure, Lipids, Glucose, Obesity)" | "Blood Pressure Phenotypes" |

---

## 使用场景

### 10.supplementary_tables.R 适用于：
- 多病共病研究项目
- 需要分析血压、血脂、血糖、肥胖等多种表型的关联
- Hyprcoloc分析包含多种表型组合（≥3个trait）
- 更全面的表型覆盖

### 11.supplementary_tables.R 适用于：
- 血压专项研究
- 只关注血压相关的5个表型
- 血压专门的Hyprcoloc分析
- 更聚焦的血压表型研究

---

## 数据依赖关系

### 10.supplementary_tables.R 依赖：
```
comorbidity/results/combined_results_11/
├── model_1/
│   ├── obs_all_results.txt
│   └── merged_obs_allmr_hyprcoloc.txt
├── model_2/
│   ├── obs_all_results.txt
│   └── merged_obs_allmr_hyprcoloc.txt
├── model_3/
│   ├── obs_all_results.txt
│   └── merged_obs_allmr_hyprcoloc.txt
├── mr/
│   └── mr_all_results.txt
├── gsmr/
│   └── gsmr_all_results.txt
├── outcome_hyprcoloc_processed/
│   └── combined_hyprcoloc_results_3.txt
├── sg_model_1/
│   └── obs_all_results.txt
├── sg_model_2/
│   └── obs_all_results.txt
└── sg_model_3/
    └── obs_all_results.txt
```

### 11.supplementary_tables.R 依赖：
```
blood pressure/
├── model_1/
│   ├── obs_all_results.txt
│   └── merged_obs_allmr_hyprcoloc.txt
├── model_2/
│   ├── obs_all_results.txt
│   └── merged_obs_allmr_hyprcoloc.txt
├── model_3/
│   ├── obs_all_results.txt
│   └── merged_obs_allmr_hyprcoloc.txt
├── mr/
│   └── mr_all_results.txt
├── gsmr/
│   └── gsmr_all_results.txt
├── hyprcoloc/
│   └── hyprcoloc_all_significant_results.txt
├── sg_model_1/
│   └── obs_all_results.txt
├── sg_model_2/
│   └── obs_all_results.txt
└── sg_model_3/
    └── obs_all_results.txt
```

---

## 输出对比

### 10.supplementary_tables.R 输出：
```
comorbidity/results/combined_results_11/
├── supplementary_all_phenotypes.xlsx (主文件)
├── Supplementary Table S1.xlsx
├── Supplementary Table S2.xlsx
├── ...
└── Supplementary Table S15.xlsx
```

### 11.supplementary_tables.R 输出：
```
blood pressure/
├── supplementary_blood_pressure.xlsx (主文件)
├── Supplementary Table S1.xlsx
├── Supplementary Table S2.xlsx
├── ...
└── Supplementary Table S15.xlsx
```

---

## 修改建议

如果需要在两个脚本之间切换或修改，只需要注意以下几点：

1. **修改prefix路径**到正确的数据目录
2. **确认hyprcoloc文件路径**是否正确
3. **检查输出文件名**避免覆盖
4. **更新表头描述**以反映实际的表型范围
5. **确认数据文件存在**，特别是hyprcoloc文件

---

## 版本信息

- **创建日期**: 2026-05-09
- **版本**: 1.0
- **维护**: 请根据实际数据路径和需求修改
