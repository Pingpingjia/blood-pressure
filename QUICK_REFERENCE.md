# 快速参考：补充表脚本选择指南

## 🔍 我应该使用哪个脚本？

### 使用 `10.supplementary_tables.R` 如果：
✅ 你的研究涉及**多种表型**（不仅仅是血压）  
✅ 你需要分析**血压、血脂、血糖、肥胖**的关联  
✅ 你的数据在 `comorbidity/results/combined_results_11/` 目录  
✅ 你需要处理包含**≥3个trait的hyprcoloc结果**  

### 使用 `11.supplementary_tables.R` 如果：
✅ 你的研究**只关注血压**  
✅ 你只需要分析 **SBP, DBP, PP, MAP, Hypertension** 这5个表型  
✅ 你的数据在 `blood pressure/` 目录  
✅ 你需要处理**血压专项的hyprcoloc结果**  

---

## 📊 表型对比

| 表型组 | 10.supplementary_tables.R | 11.supplementary_tables.R |
|--------|---------------------------|---------------------------|
| **血压** | sbp, dbp, pp, map, hypertension | sbp, dbp, pp, map, hypertension |
| **血脂** | ldlc, hdlc, tg | ❌ |
| **血糖** | glu, diabetes, fbg | ❌ |
| **肥胖** | bmi, wc, whr | ❌ |

---

## 🚀 快速使用

### 步骤 1: 确认数据路径
```r
# 对于 10.supplementary_tables.R
prefix <- "D:/OneDrive/工作/1.工作/comorbidity/results/combined_results_11/"

# 对于 11.supplementary_tables.R
prefix <- "D:/OneDrive/工作/1.工作/blood pressure/"
```

### 步骤 2: 运行脚本
```r
# 安装依赖包（如果还没安装）
install.packages(c("openxlsx", "dplyr"))

# 运行脚本
source("10.supplementary_tables.R")  # 或
source("11.supplementary_tables.R")
```

### 步骤 3: 查看输出
```
# 输出文件将保存在 prefix 指定的目录：
- supplementary_all_phenotypes.xlsx (10.supplementary_tables.R)
- supplementary_blood_pressure.xlsx (11.supplementary_tables.R)
- Supplementary Table S1.xlsx 到 S15.xlsx (两个脚本都会生成)
```

---

## ⚠️ 常见问题

### Q: 脚本报错说找不到文件？
**A:** 检查以下几点：
1. `prefix` 路径是否正确
2. 所有数据文件是否存在
3. 特别检查hyprcoloc文件路径：
   - 10: `outcome_hyprcoloc_processed/combined_hyprcoloc_results_3.txt`
   - 11: `hyprcoloc/hyprcoloc_all_significant_results.txt`

### Q: 可以同时运行两个脚本吗？
**A:** 可以！它们输出到不同的文件，不会冲突。

### Q: 如何修改表格样式？
**A:** 修改脚本中的 `createStyle()` 函数参数，例如：
```r
titleStyle <- createStyle(
  fontName = "Arial",      # 修改字体
  fontSize = 14,           # 修改字号
  textDecoration = "bold"  # 修改样式
)
```

### Q: 脚本运行很慢？
**A:** 这是正常的。处理大量数据和生成Excel文件需要时间。可以：
- 确保有足够内存
- 关闭其他占用内存的程序
- 耐心等待

---

## 📁 所需数据文件清单

### 两个脚本都需要：
- `model_1/obs_all_results.txt`
- `model_2/obs_all_results.txt`
- `model_3/obs_all_results.txt`
- `mr/mr_all_results.txt`
- `gsmr/gsmr_all_results.txt`
- `model_1/merged_obs_allmr_hyprcoloc.txt`
- `model_2/merged_obs_allmr_hyprcoloc.txt`
- `model_3/merged_obs_allmr_hyprcoloc.txt`
- `sg_model_1/obs_all_results.txt`
- `sg_model_2/obs_all_results.txt`
- `sg_model_3/obs_all_results.txt`

### 10.supplementary_tables.R 特有：
- `outcome_hyprcoloc_processed/combined_hyprcoloc_results_3.txt`

### 11.supplementary_tables.R 特有：
- `hyprcoloc/hyprcoloc_all_significant_results.txt`

---

## 📚 更多信息

详细文档请参考：
- **使用说明**: `README_SUPPLEMENTARY_TABLES.md`
- **代码对比**: `COMPARISON_SUPPLEMENTARY_TABLES.md`

---

## 🆘 需要帮助？

如果遇到问题：
1. 查看错误消息（脚本会显示详细的错误信息）
2. 检查数据文件路径和文件名
3. 确认所有依赖包已安装
4. 查阅详细文档

---

**最后更新**: 2026-05-09  
**版本**: 1.0
