"""
HyPrColoc 血压 Circos 环形图
两层结构：外圈=蛋白，内圈=区域
同一个region对应的蛋白聚在一起
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import matplotlib.patches as patches

# =========================================================
# 0. 配置
# =========================================================
BASE_DIR = "D:/OneDrive/工作/1.工作/blood pressure"
INPUT_DIR = os.path.join(BASE_DIR, "hyprcoloc")
MODELS = ["model_1", "model_2", "model_3"]

# 配色方案
TRAIT_COLORS = {
    "ckb_sbp": "#4C78A8",
    "ckb_dbp": "#F58518",
    "ckb_pp": "#54A24B",
    "ckb_map": "#E45756",
    "ukb_hypertension": "#72B7B2"
}

BP_TRAITS = ["ckb_sbp", "ckb_dbp", "ckb_pp", "ckb_map", "ukb_hypertension"]

TRAIT_MR_MAPPING = {
    "ckb_sbp": ("sbp.obs", "ckb_sbp.union"),
    "ckb_dbp": ("dbp.obs", "ckb_dbp.union"),
    "ckb_pp": ("pp.obs", "ckb_pp.union"),
    "ckb_map": ("map.obs", "ckb_map.union"),
    "ukb_hypertension": ("hypertension.obs", "ukb_hypertension.union")
}

# =========================================================
# 1. 读取共定位数据
# =========================================================
print("读取共定位数据...")
summary_file = os.path.join(INPUT_DIR, "hyprcoloc_max_posterior_prob_by_region.txt")
if not os.path.exists(summary_file):
    raise FileNotFoundError(f"未找到文件: {summary_file}")

coloc_data_raw = pd.read_csv(summary_file, sep="\t")

# =========================================================
# 2. Circos绘图函数
# =========================================================
def draw_arc(ax, center, radius, angle_start, angle_end, width, color, alpha=1.0):
    """绘制圆弧"""
    theta = np.linspace(angle_start, angle_end, 100)
    x_outer = center[0] + (radius + width/2) * np.cos(theta)
    y_outer = center[1] + (radius + width/2) * np.sin(theta)
    x_inner = center[0] + (radius - width/2) * np.cos(theta)
    y_inner = center[1] + (radius - width/2) * np.sin(theta)
    
    verts = list(zip(x_outer, y_outer)) + list(zip(x_inner[::-1], y_inner[::-1]))
    codes = [Path.MOVETO] + [Path.LINETO] * (len(verts) - 2) + [Path.CLOSEPOLY]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor=color, alpha=alpha, edgecolor='black', linewidth=0.5)
    ax.add_patch(patch)

def draw_bezier_curve(ax, center, r1, angle1, r2, angle2, color, alpha=0.3):
    """绘制贝塞尔曲线连接两个点"""
    x1 = center[0] + r1 * np.cos(angle1)
    y1 = center[1] + r1 * np.sin(angle1)
    x2 = center[0] + r2 * np.cos(angle2)
    y2 = center[1] + r2 * np.sin(angle2)
    
    # 控制点设在中心附近
    ctrl_scale = 0.3
    ctrl_x = center[0] * (1 - ctrl_scale) + (x1 + x2) / 2 * ctrl_scale
    ctrl_y = center[1] * (1 - ctrl_scale) + (y1 + y2) / 2 * ctrl_scale
    
    verts = [(x1, y1), (ctrl_x, ctrl_y), (x2, y2)]
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', edgecolor=color, 
                             linewidth=0.8, alpha=alpha)
    ax.add_patch(patch)

def draw_straight_line(ax, center, r1, angle1, r2, angle2, color, alpha=0.3):
    """绘制直线连接两个点"""
    x1 = center[0] + r1 * np.cos(angle1)
    y1 = center[1] + r1 * np.sin(angle1)
    x2 = center[0] + r2 * np.cos(angle2)
    y2 = center[1] + r2 * np.sin(angle2)
    
    ax.plot([x1, x2], [y1, y2], color=color, alpha=alpha, linewidth=0.5)

def create_circos_plot(data, model_name, output_dir):
    """创建Circos图 - 两层结构：外圈蛋白，内圈区域"""
    print(f"\n绘制 {model_name} Circos图...")
    
    # 数据预处理
    data = data.copy()
    data['region'] = data['region'].fillna('intergenic')
    data['protein'] = data['protein'].fillna('Unknown').astype(str)
    
    # 获取唯一元素
    regions = sorted(data['region'].unique())
    
    # 重新组织蛋白：按region分组，同一个region的蛋白放在一起
    print("  组织蛋白分组...")
    region_protein_map = {}
    for region in regions:
        proteins_in_region = data[data['region'] == region]['protein'].unique()
        region_protein_map[region] = sorted(proteins_in_region)
    
    # 创建蛋白的有序列表（按region分组）
    proteins_ordered = []
    for region in regions:
        proteins_ordered.extend(region_protein_map[region])
    
    print(f"  区域数: {len(regions)}")
    print(f"  蛋白数: {len(proteins_ordered)}")
    
    # 计算角度分配
    n_regions = len(regions)
    n_proteins = len(proteins_ordered)
    
    # 为region分配角度
    region_angles = {}
    angle = 0
    gap = 0.02  # 小间隙
    angle_per_region = (2 * np.pi - gap * n_regions) / n_regions
    for region in regions:
        region_angles[region] = (angle, angle + angle_per_region)
        angle += angle_per_region + gap
    
    # 为蛋白分配角度（均匀分布在外圈）
    protein_angles = {}
    angle = 0
    angle_per_protein = (2 * np.pi - gap * len(proteins_ordered)) / len(proteins_ordered)
    for protein in proteins_ordered:
        protein_angles[protein] = (angle, angle + angle_per_protein)
        angle += angle_per_protein + gap
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(16, 16))
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    center = (0, 0)
    
    # 定义半径 - 只有两层
    r_region = 0.5     # 内圈：区域
    r_protein = 0.9    # 外圈：蛋白
    
    arc_width = 0.1
    
    # 1. 绘制内圈：区域
    print("  绘制区域圈...")
    region_colors = plt.cm.Set3(np.linspace(0, 1, len(regions)))
    region_color_map = {region: region_colors[i] for i, region in enumerate(regions)}
    
    for region in regions:
        angle_start, angle_end = region_angles[region]
        color = region_color_map[region]
        draw_arc(ax, center, r_region, angle_start, angle_end, arc_width, color, alpha=0.7)
        
        # 添加区域标签
        mid_angle = (angle_start + angle_end) / 2
        label_r = r_region + arc_width/2 + 0.08
        x = label_r * np.cos(mid_angle)
        y = label_r * np.sin(mid_angle)
        rotation = np.degrees(mid_angle)
        if rotation > 90 and rotation < 270:
            rotation += 180
        ax.text(x, y, region, ha='center', va='center', 
               fontsize=8, fontweight='bold', rotation=rotation)
    
    # 2. 绘制外圈：蛋白
    print("  绘制蛋白圈...")
    for protein in proteins_ordered:
        angle_start, angle_end = protein_angles[protein]
        # 根据所属region使用相同颜色
        protein_region = data[data['protein'] == protein]['region'].iloc[0]
        color = region_color_map[protein_region]
        draw_arc(ax, center, r_protein, angle_start, angle_end, arc_width, color, alpha=0.5)
        
        # 添加蛋白标签（显示所有蛋白）
        mid_angle = (angle_start + angle_end) / 2
        label_r = r_protein + arc_width/2 + 0.05
        x = label_r * np.cos(mid_angle)
        y = label_r * np.sin(mid_angle)
        rotation = np.degrees(mid_angle)
        if rotation > 90 and rotation < 270:
            rotation += 180
        ax.text(x, y, protein, ha='center', va='center', 
               fontsize=5, rotation=rotation)
    
    # 3. 绘制连接线：Region -> Protein（从内圈外侧直线连接到外圈内侧）
    print("  绘制连接线...")
    for _, row in data.iterrows():
        region = row['region']
        protein = row['protein']
        
        if region not in region_angles or protein not in protein_angles:
            continue
        
        angle_region = (region_angles[region][0] + region_angles[region][1]) / 2
        angle_protein = (protein_angles[protein][0] + protein_angles[protein][1]) / 2
        color = region_color_map[region]
        
        # 从内圈外侧直线连接到外圈内侧
        draw_straight_line(ax, center, r_region + arc_width/2, angle_region,
                          r_protein - arc_width/2, angle_protein, color, alpha=0.2)
    
    # 添加图例
    legend_elements = [
        mpatches.Patch(color='gray', label='Proteins (Outer)', alpha=0.5),
        mpatches.Patch(color='gray', label='Regions (Inner)', alpha=0.7),
    ]
    
    ax.legend(handles=legend_elements, loc='upper left', 
             bbox_to_anchor=(0.02, 0.98), fontsize=12)
    
    # 添加标题
    plt.title(f'Blood Pressure HyPrColoc Circos — {model_name}\nProteins (Outer) connected to Regions (Inner)', 
             fontsize=18, fontweight='bold', pad=20)
    
    # 保存图表
    output_file = os.path.join(output_dir, f'HyPrColoc_BP_Circos_{model_name}_2layers.pdf')
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  已保存: {output_file}")
    
    output_file_png = os.path.join(output_dir, f'HyPrColoc_BP_Circos_{model_name}_2layers.png')
    plt.savefig(output_file_png, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  已保存: {output_file_png}")
    
    plt.close()

# =========================================================
# 3. 逐模型处理
# =========================================================
for MODEL_NAME in MODELS:
    print(f"\n{'='*70}")
    print(f"处理 {MODEL_NAME}")
    print(f"{'='*70}")
    
    MODEL_DIR = os.path.join(BASE_DIR, MODEL_NAME)
    OUTPUT_DIR = os.path.join(MODEL_DIR, "hyprcoloc")
    
    # 读取MR数据
    mr_file = os.path.join(MODEL_DIR, "merged_obs_allmr_hyprcoloc.txt")
    if not os.path.exists(mr_file):
        print(f"跳过 {MODEL_NAME}: 未找到MR文件")
        continue
    
    mr_data = pd.read_csv(mr_file, sep="\t")
    
    # 筛选显著配对
    significant_pairs = []
    for trait, (obs_col, union_col) in TRAIT_MR_MAPPING.items():
        trait_sig = mr_data[(mr_data[obs_col] == 1) & (mr_data[union_col] == 1)]
        if not trait_sig.empty:
            for protein_id in trait_sig['protein_id']:
                significant_pairs.append((protein_id, trait))
    
    if not significant_pairs:
        print(f"{MODEL_NAME}: 没有显著配对")
        continue
    
    print(f"显著配对数: {len(significant_pairs)}")
    
    # 筛选共定位数据
    sig_pairs_df = pd.DataFrame(significant_pairs, columns=['protein_id', 'trait'])
    coloc_data = coloc_data_raw.merge(
        sig_pairs_df,
        left_on=['protein_seq', 'trait'],
        right_on=['protein_id', 'trait'],
        how='inner'
    ).drop(columns=['protein_id'])
    
    if coloc_data.empty:
        print(f"{MODEL_NAME}: 筛选后无数据")
        continue
    
    # 区域过滤（可选）
    REGION_MIN_COLOC = 2
    region_freq = coloc_data.groupby('region').size().reset_index(name='n_coloc')
    high_freq_regions = region_freq[region_freq['n_coloc'] >= REGION_MIN_COLOC]['region'].tolist()
    coloc_data_filtered = coloc_data[coloc_data['region'].isin(high_freq_regions)].copy()
    
    print(f"过滤后记录数: {len(coloc_data_filtered)}")
    
    # 绘制Circos图
    create_circos_plot(coloc_data_filtered, MODEL_NAME, OUTPUT_DIR)

print("\n" + "="*70)
print("✓ 所有模型的Circos图创建完成！")
print("="*70 + "\n")