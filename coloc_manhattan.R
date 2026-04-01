setwd("C:/Users/shyne/OneDrive/intensity/protein_screening_New/resultNew3/Result")

# =================== 1. 加载必备包 ===================
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(readxl)
library(readr)
if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)
library(VennDiagram)
library(grid)

# =================== 2. 定义原始颜色方案 ===================
original_colors <- c(
  "#918BD9",  # 颜色1
  "#B3B4DF",  # 颜色2
  "#CFCCE3",  # 颜色3
  "#F6DFD6",  # 颜色4
  "#F8B2A2",  # 颜色5
  "#F1B6BB",  # 颜色6
  "#E9687A"   # 颜色7
)

# 生成22条染色体颜色向量（按7色循环）
chr_colors <- rep(original_colors, length.out = 22)

# 打印颜色分配验证
cat("染色体颜色分配:\n")
for (i in 1:22) {
  cat(sprintf("染色体 %d: %s\n", i, chr_colors[i]))
}

# =================== 3. 数据准备与处理 ===================
# 读取合并数据
fbg_merged <- read_csv("C:/Users/shyne/OneDrive/intensity/protein_screening_New/resultNew3/Result/fbg_merged.csv", 
                       na = "", show_col_types = FALSE)
dia_merged <- read_csv("C:/Users/shyne/OneDrive/intensity/protein_screening_New/resultNew3/Result/dia_merged.csv", 
                       na = "", show_col_types = FALSE)

# 基因注释数据
annotation_data <- read_excel("D:/Standard_CMCS_Database/anno.info_updat.xlsx") %>%
  filter(Organism.y == "Human") %>%
  select(AptName, `Entrez Gene Name`) %>%
  #filter(!is.na(AptName), !is.na(`Entrez Gene Name`)) %>%
  distinct(AptName, .keep_all = TRUE)

id.exposure_human <- annotation_data$AptName

# 读取共定位数据
fbg_summary_file <- "D:/colocalization_analysis/coloc_FBG_magic/clump_200kb_p_5e-06/all_significant_results_merged.csv"
dia_summary_file <- "D:/colocalization_analysis/bbj_diabetes_results/bbj_diabetes_results/all_significant_results_merged.csv"

# # 处理FBG数据
# coloc_FBG <- read_csv(
#   fbg_summary_file, 
#   show_col_types = FALSE,
#   col_select = 1:22
# ) %>%
#   filter(protein_id %in% id.exposure_human) %>%
#   distinct(protein_id, .keep_all = TRUE) %>%
#   left_join(annotation_data, by = c("protein_id" = "AptName")) %>%
#   rename(entrez_gene = `Entrez Gene Name`)
# 
# # 处理糖尿病数据
# coloc_diabetes <- read_csv(
#   dia_summary_file, 
#   show_col_types = FALSE,
#   col_select = 1:22
# ) %>%
#   mutate(protein_id = sub(".*?_", "", protein_id)) %>%
#   filter(protein_id %in% id.exposure_human) %>%
#   distinct(protein_id, .keep_all = TRUE) %>%
#   left_join(annotation_data, by = c("protein_id" = "AptName")) %>%
#   rename(entrez_gene = `Entrez Gene Name`)

# 处理FBG数据
coloc_FBG <- read_csv(
  fbg_summary_file, 
  show_col_types = FALSE,
  col_select = 1:22
) %>%
  filter(protein_id %in% id.exposure_human) %>%
  group_by(protein_id) %>%
  slice_max(order_by = posterior_prob, n = 1, with_ties = FALSE) %>%  # 保留posterior_prob最高的那行
  ungroup() %>%
  left_join(annotation_data, by = c("protein_id" = "AptName")) %>%
  rename(entrez_gene = `Entrez Gene Name`)

# 处理糖尿病数据
coloc_diabetes <- read_csv(
  dia_summary_file, 
  show_col_types = FALSE,
  col_select = 1:22
) %>%
  mutate(protein_id = sub(".*?_", "", protein_id)) %>%
  filter(protein_id %in% id.exposure_human) %>%
  group_by(protein_id) %>%
  slice_max(order_by = posterior_prob, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(annotation_data, by = c("protein_id" = "AptName")) %>%
  rename(entrez_gene = `Entrez Gene Name`)

# 计算公共集合
fbg_exp <- na.omit(unique(fbg_merged$id.exposure))
coloc_fbg <- na.omit(unique(coloc_FBG$protein_id))
fbg_common <- intersect(fbg_exp, coloc_fbg)

dia_exp <- na.omit(unique(dia_merged$id.exposure))
coloc_diab <- na.omit(unique(coloc_diabetes$protein_id))
dia_common <- intersect(dia_exp, coloc_diab)

# 验证公共集合数量
cat("FBG公共集合数量:", length(fbg_common), "\n")
cat("Diabetes公共集合数量:", length(dia_common), "\n")

# =================== 4. 数据预处理（移除is_target标记） ===================
clean_coloc_data <- function(data, trait_name, common_set) {
  data %>%
    mutate(
      chr = as.integer(chr),
      BP = as.numeric(BP),
      posterior_prob = as.numeric(posterior_prob),
      trait = trait_name,
      is_common = protein_id %in% common_set  # 仅保留公共集合标记
    ) %>%
    filter(
      !is.na(chr), 
      !is.na(BP), 
      !is.na(posterior_prob),
      !is.na(entrez_gene),
      chr %in% 1:22
    ) %>%
    group_by(chr) %>%
    mutate(
      rel_pos = scales::rescale(BP, to = c(-0.4, 0.4)),
      plot_x = chr + rel_pos,
      chr_color = chr_colors[chr]  # 用于验证颜色
    ) %>%
    ungroup()
}

data_fbg <- clean_coloc_data(coloc_FBG, trait_name = "FBG", fbg_common)
data_diab <- clean_coloc_data(coloc_diabetes, trait_name = "Diabetes", dia_common)

# =================== 5. 绘图函数（不单独标记is_target） ===================
plot_coloc_manhattan <- function(data, main_title) {
  ggplot(data, aes(x = plot_x, y = posterior_prob)) +
    # 1. 普通点（非公共集合）
    geom_point(
      data = ~filter(., !is_common),
      aes(color = as.factor(chr)),
      alpha = 0.3,
      size = 2.5,
      stroke = 0.1
    ) +
    
    # 2. 公共集合点（保持原有样式）
    geom_point(
      data = ~filter(., is_common),
      aes(fill = as.factor(chr)),
      color = "black",             # 黑色轮廓
      size = 2.5,
      alpha = 1.0,
      stroke = 0.6,
      shape = 21                   # 实心圆
    ) +
    
    # 3. 阈值线
    geom_hline(
      yintercept = 0.5,
      color = "black",
      linetype = "dashed",
      size = 0.5,
      alpha = 0.8
    ) +
    annotate(
      "text",
      x = 22.2, y = 0.53,
      label = "Posterior Probability = 0.5",
      color = "black",
      size = 3.5,
      hjust = 1
    ) +
    
    # 4. 公共集合标签（保持样式）
    ggrepel::geom_text_repel(
      data = ~filter(., is_common),
      aes(label = entrez_gene),
      color = "black",
      size = 2.5,
      fontface = "plain",
      box.padding = 0.4,
      point.padding = 0.3,
      segment.alpha = 0.6,
      segment.size = 0.3,
      max.overlaps = 200,
      min.segment.length = 0.1,
      force = 3,
      force_pull = 1,
      direction = "both"
    ) +
    
    # 5. 颜色设置（保持同步循环）
    scale_color_manual(
      values = setNames(chr_colors, 1:22),
      guide = "none"
    ) +
    scale_fill_manual(
      values = setNames(chr_colors, 1:22),
      guide = "none"
    ) +
    
    # 6. 坐标轴设置
    scale_x_continuous(
      breaks = 1:22,
      labels = 1:22,
      limits = c(0.2, 22.8),
      name = "Chromosome",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1.1),
      breaks = seq(0.0, 1, 0.2),
      name = "Posterior Probability",
      expand = c(0, 0)
    ) +
    
    # 7. 主题设置
    theme_bw() +
    theme(
      axis.line.x = element_line(color = "black", size = 0.9),
      axis.line.y = element_line(color = "black", size = 0.9),
      axis.title.x = element_text(
        color = "black", size = 13, face = "bold", margin = margin(t = 12)
      ),
      axis.title.y = element_text(
        color = "black", size = 13, face = "bold", margin = margin(r = 12)
      ),
      axis.text.x = element_text(color = "black", size = 10.5, hjust = 0.5),
      axis.text.y = element_text(color = "black", size = 10.5),
      plot.title = element_text(
        color = "black", size = 15, face = "bold", hjust = 0.5, margin = margin(b = 15)
      ),
      panel.grid.major = element_line(color = "#ECF0F1", size = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = 10, r = 50, b = 10, l = 50, unit = "mm")
    ) +
    labs(title = main_title)
}

# =================== 6. 绘制并拼接图形 ===================
p_fbg <- plot_coloc_manhattan(
  data = data_fbg,
  main_title = "FBG Colocalization (Common Sets Highlighted)"
)

p_diab <- plot_coloc_manhattan(
  data = data_diab,
  main_title = "Diabetes Colocalization (Common Sets Highlighted)"
)

p_combined <- p_fbg + p_diab +
  plot_layout(guides = "collect", ncol = 1, heights = c(1, 1))

# =================== 7. 输出PDF ===================
pdf(
  file = "Manhattan_No_Separate_Target1.pdf",
  width = 16,
  height = 16,
  paper = "special"
)

print(p_combined)
dev.off()

message("=== 已移除is_target单独标记的曼哈顿图已保存！===")
message("文件路径：", file.path(getwd(), "Manhattan_No_Separate_Target.pdf"))
message("主要修改：移除了对is_target的特殊标记，所有点仅区分普通点和公共集合点")
message("公共集合点样式保持不变：实心圆+黑色轮廓+同染色体颜色+标签")



setwd("C:/Users/shyne/OneDrive/intensity/protein_screening_New/resultNew3/Result")

# =================== 1. 加载必备包 ===================
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(readxl)
library(readr)
if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)
library(VennDiagram)
library(grid)

# =================== 2. 定义原始颜色方案（精确还原循环逻辑） ===================
# 原始7种颜色，按此顺序循环应用于22条染色体
original_colors <- c(
  "#918BD9",  # 颜色1
  "#B3B4DF",  # 颜色2
  "#CFCCE3",  # 颜色3
  "#F6DFD6",  # 颜色4
  "#F8B2A2",  # 颜色5
  "#F1B6BB",  # 颜色6
  "#E9687A"   # 颜色7
)

# 生成与原始代码完全一致的22条染色体颜色向量（关键修复）
# 按1-7颜色循环分配给1-22号染色体
chr_colors <- rep(original_colors, length.out = 22)

# 打印颜色分配验证（方便检查）
cat("染色体颜色分配:\n")
for (i in 1:22) {
  cat(sprintf("染色体 %d: %s\n", i, chr_colors[i]))
}

# =================== 3. 数据准备与处理 ===================
# 读取合并数据
fbg_merged <- read_csv("C:/Users/shyne/OneDrive/intensity/protein_screening_New/resultNew3/Result/fbg_merged.csv", 
                       na = "", show_col_types = FALSE)
dia_merged <- read_csv("C:/Users/shyne/OneDrive/intensity/protein_screening_New/resultNew3/Result/dia_merged.csv", 
                       na = "", show_col_types = FALSE)

# 基因注释数据
annotation_data <- read_excel("D:/Standard_CMCS_Database/anno.info_updat.xlsx") %>%
  filter(Organism.y == "Human") %>%
  select(AptName, `Entrez Gene Name`) %>%
  filter(!is.na(AptName), !is.na(`Entrez Gene Name`)) %>%
  distinct(AptName, .keep_all = TRUE)

id.exposure_human <- annotation_data$AptName

# 读取共定位数据
fbg_summary_file <- "D:/colocalization_analysis/coloc_FBG_magic/clump_200kb_p_5e-06/all_significant_results_merged.csv"
dia_summary_file <- "D:/colocalization_analysis/bbj_diabetes_results/bbj_diabetes_results/all_significant_results_merged.csv"

# 处理FBG数据
coloc_FBG <- read_csv(
  fbg_summary_file, 
  show_col_types = FALSE,
  col_select = 1:22
) %>%
  filter(protein_id %in% id.exposure_human) %>%
  distinct(protein_id, .keep_all = TRUE) %>%
  left_join(annotation_data, by = c("protein_id" = "AptName")) %>%
  rename(entrez_gene = `Entrez Gene Name`)

# 处理糖尿病数据
coloc_diabetes <- read_csv(
  dia_summary_file, 
  show_col_types = FALSE,
  col_select = 1:22
) %>%
  mutate(protein_id = sub(".*?_", "", protein_id)) %>%
  filter(protein_id %in% id.exposure_human) %>%
  distinct(protein_id, .keep_all = TRUE) %>%
  left_join(annotation_data, by = c("protein_id" = "AptName")) %>%
  rename(entrez_gene = `Entrez Gene Name`)

# 计算公共集合
fbg_exp <- na.omit(unique(fbg_merged$id.exposure))
coloc_fbg <- na.omit(unique(coloc_FBG$protein_id))
fbg_common <- intersect(fbg_exp, coloc_fbg)

dia_exp <- na.omit(unique(dia_merged$id.exposure))
coloc_diab <- na.omit(unique(coloc_diabetes$protein_id))
dia_common <- intersect(dia_exp, coloc_diab)

# 验证公共集合数量
cat("FBG公共集合数量:", length(fbg_common), "\n")
cat("Diabetes公共集合数量:", length(dia_common), "\n")

# =================== 4. 数据预处理（添加染色体颜色列用于验证） ===================
clean_coloc_data <- function(data, trait_name, common_set) {
  data %>%
    mutate(
      chr = as.integer(chr),
      BP = as.numeric(BP),
      posterior_prob = as.numeric(posterior_prob),
      trait = trait_name,
      is_common = protein_id %in% common_set,
      is_target = protein_id == "seq.13109.82",
      # 添加染色体颜色列用于验证
      chr_color = chr_colors[chr]
    ) %>%
    filter(
      !is.na(chr), 
      !is.na(BP), 
      !is.na(posterior_prob),
      !is.na(entrez_gene),
      chr %in% 1:22
    ) %>%
    group_by(chr) %>%
    mutate(
      rel_pos = scales::rescale(BP, to = c(-0.4, 0.4)),
      plot_x = chr + rel_pos
    ) %>%
    ungroup()
}

data_fbg <- clean_coloc_data(coloc_FBG, trait_name = "FBG", fbg_common)
data_diab <- clean_coloc_data(coloc_diabetes, trait_name = "Diabetes", dia_common)

# 验证数据中的颜色分配
cat("\nFBG数据中染色体颜色样本:\n")
print(unique(data_fbg[, c("chr", "chr_color")]))

# =================== 5. 绘图函数（精确同步颜色循环） ===================
plot_coloc_manhattan <- function(data, main_title) {
  ggplot(data, aes(x = plot_x, y = posterior_prob)) +
    # 1. 普通点（基础层）
    geom_point(
      data = ~filter(., !is_common & !is_target),
      aes(color = as.factor(chr)),  # 使用染色体颜色
      alpha = 0.9,
      size = 1.3,
      stroke = 0.1
    ) +
    
    # 2. 公共集合点（与同染色体普通点同色）
    geom_point(
      data = ~filter(., is_common & !is_target),
      aes(fill = as.factor(chr)),  # 填充色=同染色体普通点颜色
      color = "black",             # 黑色轮廓
      size = 2.0,
      alpha = 1.0,
      stroke = 0.6,
      shape = 21                   # 可填充圆形
    ) +
    
    # 3. 目标ID点（与同染色体普通点同色）
    geom_point(
      data = ~filter(., is_target),
      aes(fill = as.factor(chr)),  # 填充色=同染色体普通点颜色
      color = "black",             # 黑色轮廓
      size = 2.5,
      alpha = 1.0,
      stroke = 0.8,
      shape = 21
    ) +
    
    # 4. 阈值线
    geom_hline(
      yintercept = 0.5,
      color = "black",
      linetype = "dashed",
      size = 0.5,
      alpha = 0.8
    ) +
    annotate(
      "text",
      x = 22.2, y = 0.53,
      label = "Posterior Probability = 0.5",
      color = "black",
      size = 3.5,
      hjust = 1
    ) +
    
    # 5. 公共集合标签
    ggrepel::geom_text_repel(
      data = ~filter(., is_common),
      aes(label = entrez_gene),
      color = "black",
      size = 2.5,
      fontface = "plain",
      box.padding = 0.4,
      point.padding = 0.3,
      segment.alpha = 0.6,
      segment.size = 0.3,
      max.overlaps = 200,
      min.segment.length = 0.1,
      force = 3,
      force_pull = 1,
      direction = "both"
    ) +
    
    # 6. 目标ID标签
    ggrepel::geom_text_repel(
      data = ~filter(., is_target),
      aes(label = entrez_gene),
      color = "black",
      size = 2.5,
      fontface = "plain",
      box.padding = 0.4,
      point.padding = 0.3,
      segment.alpha = 0.6,
      segment.size = 0.3,
      max.overlaps = 200,
      min.segment.length = 0.1
    ) +
    
    # 7. 颜色设置（核心：普通点和填充色使用完全相同的循环方案）
    scale_color_manual(
      values = setNames(chr_colors, 1:22),  # 按染色体号精确映射
      guide = "none"
    ) +
    scale_fill_manual(
      values = setNames(chr_colors, 1:22),  # 填充色与普通点颜色完全同步
      guide = "none"
    ) +
    
    # 8. 坐标轴设置
    scale_x_continuous(
      breaks = 1:22,
      labels = 1:22,
      limits = c(0.2, 22.8),
      name = "Chromosome",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1.1),
      breaks = seq(0, 1, 0.2),
      name = "Posterior Probability",
      expand = c(0, 0)
    ) +
    
    # 9. 主题设置
    theme_bw() +
    theme(
      axis.line.x = element_line(color = "black", size = 0.9),
      axis.line.y = element_line(color = "black", size = 0.9),
      axis.title.x = element_text(
        color = "black", size = 13, face = "bold", margin = margin(t = 12)
      ),
      axis.title.y = element_text(
        color = "black", size = 13, face = "bold", margin = margin(r = 12)
      ),
      axis.text.x = element_text(color = "black", size = 10.5, hjust = 0.5),
      axis.text.y = element_text(color = "black", size = 10.5),
      plot.title = element_text(
        color = "black", size = 15, face = "bold", hjust = 0.5, margin = margin(b = 15)
      ),
      panel.grid.major = element_line(color = "#ECF0F1", size = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = 10, r = 50, b = 10, l = 50, unit = "mm")
    ) +
    labs(title = main_title)
}

# =================== 6. 绘制并拼接图形 ===================
p_fbg <- plot_coloc_manhattan(
  data = data_fbg,
  main_title = "FBG Colocalization (Common Sets Highlighted)"
)

p_diab <- plot_coloc_manhattan(
  data = data_diab,
  main_title = "Diabetes Colocalization (Common Sets Highlighted)"
)

p_combined <- p_fbg + p_diab +
  plot_layout(guides = "collect", ncol = 1, heights = c(1.2, 1))

# =================== 7. 输出PDF ===================
pdf(
  file = "Manhattan_Color_Cycle_Matched.pdf",
  width = 16,
  height = 18,
  paper = "special"
)

print(p_combined)
dev.off()

message("=== 颜色循环匹配的曼哈顿图已保存！===")
message("文件路径：", file.path(getwd(), "Manhattan_Color_Cycle_Matched.pdf"))
message("颜色循环规则：每7条染色体重复一次原始颜色序列")
message("已通过名称映射确保染色体与颜色严格对应")




#====突出交集版=====

setwd("C:/Users/shyne/OneDrive/intensity/protein_screening_New/resultNew3/Result")

# =================== 1. 加载必备包 ===================
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(readxl)
library(readr)
if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)

# =================== 2. 读取数据与基因注释 ===================
annotation_data <- read_excel("D:/Standard_CMCS_Database/anno.info_updat.xlsx") %>%
  filter(Organism.y == "Human") %>%
  select(AptName, `Entrez Gene Name`) %>%
  filter(!is.na(AptName), !is.na(`Entrez Gene Name`)) %>%
  distinct(AptName, .keep_all = TRUE)

# 处理FBG数据（保留id_exposure列用于筛选）
coloc_FBG <- read_csv(
  "D:/colocalization_analysis/coloc_FBG_magic/clump_200kb_p_5e-06/all_significant_results_merged.csv",
  show_col_types = FALSE,
  col_select = 1:22
) %>%
  filter(protein_id %in% annotation_data$AptName) %>%
  distinct(protein_id, .keep_all = TRUE) %>%
  left_join(annotation_data, by = c("protein_id" = "AptName")) %>%
  rename(entrez_gene = `Entrez Gene Name`)

# 处理Diabetes数据（保留id_exposure列用于筛选）
coloc_diabetes <- read_csv(
  "D:/colocalization_analysis/bbj_diabetes_results/bbj_diabetes_results/all_significant_results_merged.csv",
  show_col_types = FALSE,
  col_select = 1:22
) %>%
  mutate(protein_id = sub(".*?_", "", protein_id)) %>%
  filter(protein_id %in% annotation_data$AptName) %>%
  distinct(protein_id, .keep_all = TRUE) %>%
  left_join(annotation_data, by = c("protein_id" = "AptName")) %>%
  rename(entrez_gene = `Entrez Gene Name`)

# =================== 3. 数据预处理（新增特定ID标记） ===================
clean_coloc_data <- function(data, trait_name) {
  data %>%
    mutate(
      chr = as.integer(chr),
      BP = as.numeric(BP),
      posterior_prob = as.numeric(posterior_prob),
      trait = trait_name,
      is_top30 = dense_rank(desc(posterior_prob)) <= 30,
      # 新增：标记id_exposure为seq.13109.82的点
      is_target = protein_id == "seq.13109.82"
    ) %>%
    filter(
      !is.na(chr), 
      !is.na(BP), 
      !is.na(posterior_prob),
      !is.na(entrez_gene),
      chr %in% 1:22
    ) %>%
    group_by(chr) %>%
    mutate(
      rel_pos = scales::rescale(BP, to = c(-0.4, 0.4)),
      plot_x = chr + rel_pos
    ) %>%
    ungroup()
}

data_fbg <- clean_coloc_data(coloc_FBG, trait_name = "FBG")
data_diab <- clean_coloc_data(coloc_diabetes, trait_name = "Diabetes")

# 验证目标ID是否存在
cat("FBG中seq.13109.82的数量:", sum(data_fbg$is_target), "\n")
cat("Diabetes中seq.13109.82的数量:", sum(data_diab$is_target), "\n")

# =================== 4. 绘图函数（新增特定ID突出显示） ===================
plot_coloc_manhattan <- function(data, main_title, y_lab = "Posterior Probability") {
  ggplot(data, aes(x = plot_x, y = posterior_prob)) +
    # 1. 普通点（非Top30且非目标ID）
    geom_point(
      data = ~filter(., !is_top30 & !is_target),
      aes(color = as.factor(chr)),
      alpha = 0.9,
      size = 1.3,
      stroke = 0.1
    ) +
    
    # 2. Top30点（非目标ID）
    geom_point(
      data = ~filter(., is_top30 & !is_target),
      aes(color = as.factor(chr)),
      size = 2.0,
      alpha = 1.0,
      stroke = 0.6,
      shape = 21,
      fill = NA
    ) +
    
    # 3. 目标ID点（seq.13109.82，用更突出的样式）
    geom_point(
      data = ~filter(., is_target),
      aes(color = as.factor(chr)),  # 保留原有染色体颜色
      size = 2.5,                   # 比Top30更大
      alpha = 1.0,
      stroke = 1.0,                 # 更粗的黑色轮廓
      shape = 21,
      fill = NA
    ) +
    
    # 4. 阈值线
    geom_hline(
      yintercept = 0.5,
      color = "black",
      linetype = "dashed",
      size = 0.5,
      alpha = 0.8
    ) +
    annotate(
      "text",
      x = 22.2, y = 0.53,
      label = "Posterior Probability = 0.5",
      color = "black",
      size = 3.5,
      hjust = 1
    ) +
    
    # 5. Top30标签
    ggrepel::geom_text_repel(
      data = ~filter(., is_top30 & !is_target),
      aes(label = entrez_gene),
      color = "black",
      size = 2.5,
      fontface = "plain",
      box.padding = 0.3,
      point.padding = 0.2,
      segment.alpha = 0.6,
      segment.size = 0.3,
      max.overlaps = 50,
      min.segment.length = 0.1
    ) +
    
    # 6. 目标ID标签（更突出）
    ggrepel::geom_text_repel(
      data = ~filter(., is_target),
      aes(label = entrez_gene),
      color = "darkred",            # 用红色区分
      size = 2.5,                   # 比Top30标签稍大
      fontface = "bold",            # 加粗
      box.padding = 0.4,
      point.padding = 0.2,
      segment.alpha = 0.6,
      segment.size = 0.3,
      max.overlaps = 50,
      min.segment.length = 0.1
    ) +
    
    # 7. 坐标轴与颜色设置
    scale_x_continuous(
      breaks = 1:22,
      labels = 1:22,
      limits = c(0.2, 22.8),
      name = "Chromosome",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1.1),
      breaks = seq(0, 1, 0.2),
      name = y_lab,
      expand = c(0, 0)
    ) +
    scale_color_manual(
      values = rep(
        c("#918BD9","#B3B4DF", "#CFCCE3", "#F6DFD6",
          "#F8B2A2", "#F1B6BB", "#E9687A"),
        4
      ),
      guide = "none"
    ) +
    
    # 8. 主题设置
    theme_bw() +
    theme(
      axis.line.x = element_line(color = "black", size = 0.9),
      axis.line.y = element_line(color = "black", size = 0.9),
      axis.title.x = element_text(
        color = "black", size = 13, face = "bold", margin = margin(t = 12)
      ),
      axis.title.y = element_text(
        color = "black", size = 13, face = "bold", margin = margin(r = 12)
      ),
      axis.ticks.x = element_line(color = "black", size = 0.8),
      axis.ticks.y = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.12, "cm"),
      axis.text.x = element_text(color = "black", size = 10.5, angle = 0, hjust = 0.5),
      axis.text.y = element_text(color = "black", size = 10.5),
      plot.title = element_text(
        color = "black", size = 15, face = "bold", hjust = 0.5, margin = margin(b = 15)
      ),
      panel.grid.major = element_line(color = "#ECF0F1", size = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 15, unit = "mm")
    ) +
    labs(
      title = main_title,
      x = "Chromosome",
      y = y_lab
    )
}

# =================== 5. 绘制并拼接子图 ===================
p_fbg <- plot_coloc_manhattan(
  data = data_fbg,
  main_title = "Colocalization with FBG (Top30 + seq.13109.82 Highlighted)",
  y_lab = "Posterior Probability"
)

p_diab <- plot_coloc_manhattan(
  data = data_diab,
  main_title = "Colocalization with Diabetes (Top30 + seq.13109.82 Highlighted)",
  y_lab = "Posterior Probability"
)

p_combined <- p_fbg + p_diab +
  plot_layout(guides = "collect", ncol = 1)

# =================== 6. 输出PDF文件 ===================
pdf(
  file = "Protein_FBG_Diabetes_Coloc_With_TargetID.pdf",
  width = 12,
  height = 14,  # 增加高度以容纳可能的额外标签
  paper = "special"
)

print(p_combined)
dev.off()

message("=== 包含特定ID突出显示的曼哈顿图已保存！===")
message("文件路径：", file.path(getwd(), "Protein_FBG_Diabetes_Coloc_With_TargetID.pdf"))
message("突出显示: id_exposure = seq.13109.82 (更大尺寸+更粗轮廓+红色标签)")



#====不突出交集版====
setwd("C:/Users/shyne/OneDrive/intensity/protein_screening_New/resultNew3/Result")

# =================== 1. 加载必备包 ===================
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(readxl)
library(readr)
if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)

# =================== 2. 读取数据与基因注释 ===================
annotation_data <- read_excel("D:/Standard_CMCS_Database/anno.info_updat.xlsx") %>%
  filter(Organism.y == "Human") %>%
  select(AptName, `Entrez Gene Name`) %>%
  filter(!is.na(AptName), !is.na(`Entrez Gene Name`)) %>%
  distinct(AptName, .keep_all = TRUE)

# 处理FBG数据
coloc_FBG <- read_csv(
  "D:/colocalization_analysis/coloc_FBG_magic/clump_200kb_p_5e-06/all_significant_results_merged.csv",
  show_col_types = FALSE,
  col_select = 1:22
) %>%
  filter(protein_id %in% annotation_data$AptName) %>%
  distinct(protein_id, .keep_all = TRUE) %>%
  left_join(annotation_data, by = c("protein_id" = "AptName")) %>%
  rename(entrez_gene = `Entrez Gene Name`)

# 处理Diabetes数据
coloc_diabetes <- read_csv(
  "D:/colocalization_analysis/bbj_diabetes_results/bbj_diabetes_results/all_significant_results_merged.csv",
  show_col_types = FALSE,
  col_select = 1:22
) %>%
  mutate(protein_id = sub(".*?_", "", protein_id)) %>%
  filter(protein_id %in% annotation_data$AptName) %>%
  distinct(protein_id, .keep_all = TRUE) %>%
  left_join(annotation_data, by = c("protein_id" = "AptName")) %>%
  rename(entrez_gene = `Entrez Gene Name`)

# =================== 3. 数据预处理 ===================
clean_coloc_data <- function(data, trait_name) {
  data %>%
    mutate(
      chr = as.integer(chr),
      BP = as.numeric(BP),
      posterior_prob = as.numeric(posterior_prob),
      trait = trait_name,
      is_top30 = dense_rank(desc(posterior_prob)) <= 30
    ) %>%
    filter(
      !is.na(chr), 
      !is.na(BP), 
      !is.na(posterior_prob),
      !is.na(entrez_gene),
      chr %in% 1:22
    ) %>%
    group_by(chr) %>%
    mutate(
      rel_pos = scales::rescale(BP, to = c(-0.4, 0.4)),
      plot_x = chr + rel_pos
    ) %>%
    ungroup()
}

data_fbg <- clean_coloc_data(coloc_FBG, trait_name = "FBG")
data_diab <- clean_coloc_data(coloc_diabetes, trait_name = "Diabetes")

# 验证数据量
cat("FBG数据总量:", nrow(data_fbg), "，Top30数量:", sum(data_fbg$is_top30), "\n")
cat("Diabetes数据总量:", nrow(data_diab), "，Top30数量:", sum(data_diab$is_top30), "\n")

# =================== 4. 绘图函数（核心调整） ===================
plot_coloc_manhattan <- function(data, main_title, y_lab = "Posterior Probability") {
  ggplot(data, aes(x = plot_x, y = posterior_prob)) +
    # 1. 普通点（非Top30）
    geom_point(
      data = ~filter(., !is_top30),
      aes(color = as.factor(chr)),
      alpha = 0.9,
      size = 1.3,
      stroke = 0.1
    ) +
    
    # 2. Top30点（核心修改：原有颜色+黑色轮廓）
    geom_point(
      data = ~filter(., is_top30),
      aes(color = as.factor(chr)),  # 保留原有染色体颜色
      size = 2.0,                   # 适当放大
      alpha = 1.0,
      stroke = 0.6,                 # 黑色轮廓
      shape = 21,                   # 使用带轮廓的形状
      fill = NA                     # 透明填充，显示原有颜色
    ) +
    
    # 3. 阈值线
    geom_hline(
      yintercept = 0.5,
      color = "black",
      linetype = "dashed",
      size = 0.5,
      alpha = 0.8
    ) +
    annotate(
      "text",
      x = 22.2, y = 0.53,
      label = "Posterior Probability = 0.5",
      color = "black",
      size = 3.5,
      hjust = 1
    ) +
    
    # 4. Top30标签（修改：非斜体+更小字体）
    ggrepel::geom_text_repel(
      data = ~filter(., is_top30),
      aes(label = entrez_gene),
      color = "black",
      size = 2.5,                   # 字体更小
      fontface = "plain",           # 非斜体
      box.padding = 0.3,
      point.padding = 0.2,
      segment.alpha = 0.6,
      segment.size = 0.3,
      max.overlaps = 50,
      min.segment.length = 0.1
    ) +
    
    # 5. 坐标轴与颜色设置
    scale_x_continuous(
      breaks = 1:22,
      labels = 1:22,
      limits = c(0.2, 22.8),
      name = "Chromosome",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1.1),
      breaks = seq(0, 1, 0.2),
      name = y_lab,
      expand = c(0, 0)
    ) +
    scale_color_manual(
      values = rep(
        c("#918BD9","#B3B4DF", "#CFCCE3", "#F6DFD6",
          "#F8B2A2", "#F1B6BB", "#E9687A"),
        4
      ),
      guide = "none"
    ) +
    
    # 6. 主题设置
    theme_bw() +
    theme(
      axis.line.x = element_line(color = "black", size = 0.9),
      axis.line.y = element_line(color = "black", size = 0.9),
      axis.title.x = element_text(
        color = "black", size = 13, face = "bold", margin = margin(t = 12)
      ),
      axis.title.y = element_text(
        color = "black", size = 13, face = "bold", margin = margin(r = 12)
      ),
      axis.ticks.x = element_line(color = "black", size = 0.8),
      axis.ticks.y = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.12, "cm"),
      axis.text.x = element_text(color = "black", size = 10.5, angle = 0, hjust = 0.5),
      axis.text.y = element_text(color = "black", size = 10.5),
      plot.title = element_text(
        color = "black", size = 15, face = "bold", hjust = 0.5, margin = margin(b = 15)
      ),
      panel.grid.major = element_line(color = "#ECF0F1", size = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 15, unit = "mm")
    ) +
    labs(
      title = main_title,
      x = "Chromosome",
      y = y_lab
    )
}

# =================== 5. 绘制并拼接子图 ===================
p_fbg <- plot_coloc_manhattan(
  data = data_fbg,
  main_title = "Colocalization with FBG (Top30 Genes Labeled)",
  y_lab = "Posterior Probability"
)

p_diab <- plot_coloc_manhattan(
  data = data_diab,
  main_title = "Colocalization with Diabetes (Top30 Genes Labeled)",
  y_lab = "Posterior Probability"
)

p_combined <- p_fbg + p_diab +
  plot_layout(guides = "collect", ncol = 1)

# =================== 6. 输出PDF文件 ===================
pdf(
  file = "Protein_FBG_Diabetes_Coloc_Top30_Adjusted.pdf",
  width = 12,
  height = 12,
  paper = "special"
)

print(p_combined)
dev.off()

message("=== 调整后的曼哈顿图已保存！===")
message("文件路径：", file.path(getwd(), "Protein_FBG_Diabetes_Coloc_Top30_Adjusted.pdf"))
message("Top30点样式：原有染色体颜色 + 黑色轮廓")
message("标签样式：非斜体，更小字体(2.5)")
