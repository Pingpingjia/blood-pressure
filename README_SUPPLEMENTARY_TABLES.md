# Supplementary Tables Scripts - 使用说明

## 概述

本项目包含两个补充表生成脚本，用于不同范围的表型分析：

## 脚本说明

### 10.supplementary_tables.R - 全表型补充表
**用途**: 生成包含所有表型的补充表（血压、血脂、肥胖、血糖）

**表型范围**:
- **血压 (Blood Pressure)**: sbp, dbp, pp, map, hypertension
- **血脂 (Lipids)**: ldlc, hdlc, tg
- **血糖 (Glucose)**: glu, diabetes, fbg
- **肥胖 (Obesity)**: bmi, wc, whr

**数据路径**: `D:/OneDrive/工作/1.工作/comorbidity/results/combined_results_11/`

**输出文件**: 
- `supplementary_all_phenotypes.xlsx` - 主文件（包含所有表格）
- `Supplementary Table S1.xlsx` - 单独的表格文件（逐个）

**关键特性**:
- 处理多种表型组合的hyprcoloc结果
- 使用`combined_hyprcoloc_results_3.txt`（≥3个trait的组合）
- 适合多病共病研究

---

### 11.supplementary_tables.R - 血压专项补充表
**用途**: 专门生成血压相关表型的补充表

**表型范围**:
- **SBP (收缩压)**: Systolic Blood Pressure
- **DBP (舒张压)**: Diastolic Blood Pressure  
- **PP (脉压)**: Pulse Pressure
- **MAP (平均动脉压)**: Mean Arterial Pressure
- **Hypertension (高血压)**: Hypertension status

**数据路径**: `D:/OneDrive/工作/1.工作/blood pressure/`

**输出文件**:
- `supplementary_blood_pressure.xlsx` - 主文件（包含所有表格）
- `Supplementary Table S1.xlsx` - 单独的表格文件（逐个）

**关键特性**:
- 专注于血压相关的5个表型
- 使用血压专项的hyprcoloc结果
- 适合血压单一研究

---

## 两个脚本的主要区别

| 特性 | 10.supplementary_tables.R | 11.supplementary_tables.R |
|------|---------------------------|---------------------------|
| **表型数量** | ~15个（4组表型） | 5个（仅血压） |
| **数据路径** | comorbidity/results/combined_results_11/ | blood pressure/ |
| **Hyprcoloc文件** | combined_hyprcoloc_results_3.txt | hyprcoloc_all_significant_results.txt |
| **研究范围** | 多病共病研究 | 血压单一研究 |
| **输出文件名** | supplementary_all_phenotypes.xlsx | supplementary_blood_pressure.xlsx |

---

## 补充表结构

两个脚本都生成15张补充表，结构相同：

1. **Table S1**: 研究人群特征
2. **Table S2-S4**: Model 1-3 的观察性研究回归结果
3. **Table S5**: GWAS汇总数据
4. **Table S6**: Two-sample MR结果
5. **Table S7**: GSMR结果
6. **Table S8**: Hyprcoloc共定位结果
7. **Table S9-S11**: Model 1-3 的整合结果
8. **Table S12**: 核心蛋白信息
9. **Table S13-S15**: 验证队列Model 1-3结果

---

## 使用方法

### 运行脚本

```r
# 对于全表型分析
source("10.supplementary_tables.R")

# 对于血压专项分析
source("11.supplementary_tables.R")
```

### 前提条件

1. 确保已安装必需的R包：
```r
install.packages(c("openxlsx", "dplyr"))
```

2. 确保数据文件路径正确：
   - 修改脚本开头的 `prefix` 变量
   - 确认所有数据文件存在

### 错误处理

两个脚本都包含错误处理机制：
- 如果数据文件不存在，会显示警告但继续执行
- 会生成主Excel文件和单独的表格文件
- 所有错误信息会输出到控制台

---

## 注意事项

1. **数据路径**: 根据实际情况修改 `prefix` 变量
2. **内存使用**: 处理大量数据时可能需要较大内存
3. **文件覆盖**: 脚本会覆盖已存在的输出文件
4. **编码**: 使用UTF-8编码保存脚本文件

---

## 常见问题

**Q: 如何选择使用哪个脚本？**
- 如果研究涉及多种表型（血压、血脂、肥胖、血糖），使用 `10.supplementary_tables.R`
- 如果只研究血压相关表型，使用 `11.supplementary_tables.R`

**Q: 数据文件缺失怎么办？**
- 脚本会继续执行，缺失的表格会显示错误信息
- 检查数据路径和文件名是否正确

**Q: 输出文件在哪里？**
- 输出文件保存在 `prefix` 指定的路径下
- 包括一个主Excel文件和多个单独的表格文件

---

## 更新日志

**2026-05-09**
- 创建 `10.supplementary_tables.R` 用于全表型分析
- 更新 `11.supplementary_tables.R` 增强错误处理
- 添加详细的注释和文档
- 统一两个脚本的代码结构

---

## 联系方式

如有问题或建议，请联系项目维护者。
