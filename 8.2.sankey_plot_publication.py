"""
HyPrColoc 血压表型桑吉图（发表版）
为每个模型分别生成 Trait->Region->Gene 桑吉图
只包含 obs+union 都显著的蛋白
"""

import os
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

# 设置为静态图像输出
pio.kaleido.scope.mathjax = None

# =========================================================
# 0. 路径与基础配置
# =========================================================
BASE_DIR = "D:/OneDrive/工作/1.工作/blood pressure"
INPUT_DIR = os.path.join(BASE_DIR, "hyprcoloc")
MODELS = ["model_1", "model_2", "model_3"]

# =========================================================
# 1. 配色方案定义
# =========================================================
TRAIT_GROUPS = {
    "SBP": {
        "traits": ["ckb_sbp"],
        "color": "#80B1D3"
    },
    "DBP": {
        "traits": ["ckb_dbp"],
        "color": "#F58518"
    },
    "PP": {
        "traits": ["ckb_pp"],
        "color": "#66C2A5"
    },
    "MAP": {
        "traits": ["ckb_map"],
        "color": "#E45756"
    },
    "Hypertension": {
        "traits": ["ukb_hypertension"],
        "color": "#BEBADA"
    }
}

group_colors = {k: v["color"] for k, v in TRAIT_GROUPS.items()}
group_colors["other"] = "#999999"

trait_to_group = {}
for group_name, group_info in TRAIT_GROUPS.items():
    for trait in group_info["traits"]:
        trait_to_group[trait] = group_name

BP_TRAITS = ["ckb_sbp", "ckb_dbp", "ckb_pp", "ckb_map", "ukb_hypertension"]

# MR 显著性映射
TRAIT_MR_MAPPING = {
    "ckb_sbp": ("sbp.obs", "ckb_sbp.union"),
    "ckb_dbp": ("dbp.obs", "ckb_dbp.union"),
    "ckb_pp": ("pp.obs", "ckb_pp.union"),
    "ckb_map": ("map.obs", "ckb_map.union"),
    "ukb_hypertension": ("hypertension.obs", "ukb_hypertension.union")
}

NODE_COLORS = {
    "region": "#B8E6D5",
    "protein": "#98D8C8",
    "default": "#000000"
}

# =========================================================
# 2. 读取共定位数据
# =========================================================
print("读取血压共定位结果...")
summary_file = os.path.join(INPUT_DIR, "hyprcoloc_max_posterior_prob_by_region.txt")
if not os.path.exists(summary_file):
    raise FileNotFoundError(f"未找到共定位汇总文件: {summary_file}")

coloc_data_raw = pd.read_csv(summary_file, sep="\t")

# =========================================================
# 3. 辅助函数
# =========================================================
def hex_to_rgba(hex_color, alpha=0.3):
    """将hex颜色转换为rgba格式"""
    r, g, b = int(hex_color[1:3], 16), int(hex_color[3:5], 16), int(hex_color[5:7], 16)
    return f'rgba({r}, {g}, {b}, {alpha})'

def get_primary_group(series):
    """获取系列中最常见的分组"""
    counts = series.value_counts()
    if counts.empty:
        return "other"
    return counts.index[0]

# =========================================================
# 4. 逐模型处理和绘制桑吉图
# =========================================================
for MODEL_NAME in MODELS:
    print(f"\n{'='*70}")
    print(f"处理 {MODEL_NAME}")
    print(f"{'='*70}\n")
    
    MODEL_DIR = os.path.join(BASE_DIR, MODEL_NAME)
    OUTPUT_DIR = os.path.join(MODEL_DIR, "hyprcoloc_sankey")
    
    try:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        print(f"输出目录: {OUTPUT_DIR}")
    except Exception as e:
        print(f"警告: 创建输出目录失败: {e}")
        continue
    
    # 读取 MR 验证数据
    mr_file = os.path.join(MODEL_DIR, "merged_obs_allmr_hyprcoloc.txt")
    if not os.path.exists(mr_file):
        print(f"跳过 {MODEL_NAME}: 未找到 MR 文件 {mr_file}")
        continue
    
    mr_data = pd.read_csv(mr_file, sep="\t")
    
    # 为每个表型筛选 obs+union 都显著的蛋白，创建 (protein, trait) 配对
    significant_pairs = []
    for trait, (obs_col, union_col) in TRAIT_MR_MAPPING.items():
        trait_sig = mr_data[(mr_data[obs_col] == 1) & (mr_data[union_col] == 1)]
        if not trait_sig.empty:
            for protein_id in trait_sig['protein_id']:
                significant_pairs.append((protein_id, trait))
    
    if not significant_pairs:
        print(f"{MODEL_NAME}: 没有 obs+union 都显著的蛋白-表型配对")
        continue
    
    print(f"MR 显著的蛋白-表型配对数: {len(significant_pairs)}")
    
    # 创建显著配对的 DataFrame 用于筛选
    sig_pairs_df = pd.DataFrame(significant_pairs, columns=['protein_id', 'trait'])
    
    # 筛选共定位数据：只保留蛋白和表型都匹配的记录
    coloc_data = coloc_data_raw.merge(
        sig_pairs_df,
        left_on=['protein_seq', 'trait'],
        right_on=['protein_id', 'trait'],
        how='inner'
    ).drop(columns=['protein_id'])
    
    if coloc_data.empty:
        print(f"{MODEL_NAME}: 筛选后无共定位记录")
        continue
    
    # 数据预处理
    coloc_data['region'] = coloc_data['region'].fillna('intergenic')
    coloc_data['protein'] = coloc_data['protein'].fillna('Unknown').astype(str)
    coloc_data['trait_group'] = coloc_data['trait'].map(trait_to_group).fillna("other")
    
    print(f"\n数据概览:")
    print(f"  共定位记录数: {len(coloc_data)}")
    print(f"  唯一蛋白质数: {coloc_data['protein'].nunique()}")
    print(f"  唯一表型数: {coloc_data['trait'].nunique()}")
    print(f"  唯一区域数: {coloc_data['region'].nunique()}")
    
    # 保存筛选后的数据
    filtered_file = os.path.join(OUTPUT_DIR, "hyprcoloc_bp_regions_obs_union_sig.txt")
    try:
        coloc_data.to_csv(filtered_file, sep="\t", index=False)
        print(f"\n已保存筛选数据: {filtered_file}")
    except Exception as e:
        print(f"\n警告: 保存筛选数据失败: {e}")
    
    # 区域过滤
    REGION_MIN_COLOC = 2
    region_freq = coloc_data.groupby('region').size().reset_index(name='n_coloc')
    high_freq_regions = region_freq[region_freq['n_coloc'] >= REGION_MIN_COLOC]['region'].tolist()
    
    coloc_data_filtered = coloc_data[coloc_data['region'].isin(high_freq_regions)].copy()
    print(f"\n过滤后保留的区域数: {len(high_freq_regions)} (共定位次数 >= {REGION_MIN_COLOC})")
    print(f"过滤后的共定位记录数: {len(coloc_data_filtered)}")
    
    # 创建两个版本的数据：包含intergenic和不包含intergenic
    data_versions = {
        "with_intergenic": coloc_data_filtered,
        "without_intergenic": coloc_data_filtered[coloc_data_filtered['region'] != 'intergenic'].copy()
    }
    
    print(f"不包含intergenic的记录数: {len(data_versions['without_intergenic'])}")
    
    # 为每个版本生成桑吉图
    for version_name, data_version in data_versions.items():
        if data_version.empty:
            print(f"\n{version_name}: 数据为空，跳过")
            continue
        
        print(f"\n{'='*50}")
        print(f"准备桑吉图数据 - {version_name}...")
        print(f"{'='*50}")
        
        # Trait -> Region 连接
        links_trait_region = (data_version.groupby(['trait', 'region', 'trait_group'])
                              .size()
                              .reset_index(name='value'))
        
        # Region -> Gene 连接 (使用protein作为基因名)
        links_region_gene = (data_version.groupby(['region', 'protein'])
                             .agg(value=('trait', 'size'),
                                  primary_group=('trait_group', get_primary_group))
                             .reset_index())
        
        # 创建节点列表
        unique_traits = set(data_version['trait'].unique())
        traits_by_group = {}
        for group_name, group_info in TRAIT_GROUPS.items():
            group_traits = [t for t in group_info["traits"] if t in unique_traits]
            traits_by_group[group_name] = group_traits
        
        all_traits = sum(traits_by_group.values(), [])
        all_regions = data_version.groupby('region').size().sort_values(ascending=False).index.tolist()
        all_genes = sorted(data_version['protein'].unique())
        
        if not all_traits or not all_regions or not all_genes:
            print(f"{version_name}: 节点数量不足以构建桑吉图")
            continue
        
        print(f"\n节点统计:")
        print(f"  Traits: {len(all_traits)}")
        print(f"  Regions: {len(all_regions)}")
        print(f"  Genes: {len(all_genes)}")
        print(f"  Total: {len(all_traits) + len(all_regions) + len(all_genes)}")
        
        # 创建节点标签
        node_labels_display = all_traits + all_regions + all_genes
        
        # 创建索引映射
        node_to_idx = {}
        for i, trait in enumerate(all_traits):
            node_to_idx[('trait', trait)] = i
        for i, region in enumerate(all_regions):
            node_to_idx[('region', region)] = len(all_traits) + i
        for i, gene in enumerate(all_genes):
            node_to_idx[('gene', gene)] = len(all_traits) + len(all_regions) + i
        
        # 创建节点颜色
        def get_node_color(i):
            if i < len(all_traits):
                trait = all_traits[i]
                return group_colors.get(trait_to_group.get(trait, "other"), NODE_COLORS["default"])
            elif i < len(all_traits) + len(all_regions):
                return NODE_COLORS["region"]
            else:
                return NODE_COLORS["protein"]
        
        node_colors = [get_node_color(i) for i in range(len(node_labels_display))]
        
        # 准备连接数据
        source_indices = []
        target_indices = []
        values = []
        link_colors = []
        
        # Trait -> Region 连接
        for _, row in links_trait_region.iterrows():
            source_indices.append(node_to_idx[('trait', row['trait'])])
            target_indices.append(node_to_idx[('region', row['region'])])
            values.append(row['value'])
            color = group_colors.get(row['trait_group'], NODE_COLORS["default"])
            link_colors.append(hex_to_rgba(color))
        
        # Region -> Gene 连接
        for _, row in links_region_gene.iterrows():
            source_indices.append(node_to_idx[('region', row['region'])])
            target_indices.append(node_to_idx[('gene', row['protein'])])
            values.append(row['value'])
            color = group_colors.get(row['primary_group'], NODE_COLORS["default"])
            link_colors.append(hex_to_rgba(color))
        
        # 创建桑吉图
        print("\n创建桑吉图...")
        
        max_nodes = max(len(all_traits), len(all_regions), len(all_genes))
        FIGURE_CONFIG = {
            "width": 1200,
            "base_height": 800,
            "height_per_node": 20
        }
        fig_height = FIGURE_CONFIG["base_height"] + max_nodes * FIGURE_CONFIG["height_per_node"]
        
        fig = go.Figure(data=[go.Sankey(
            arrangement='snap',
            node=dict(
                pad=5,
                thickness=8,
                line=dict(color="black", width=1),
                label=node_labels_display,
                color=node_colors,
                x=[0.0] * len(all_traits) + [0.5] * len(all_regions) + [0.98] * len(all_genes),
                hovertemplate='<b>%{label}</b><br>Total: %{value}<extra></extra>',
                hoverlabel=dict(font=dict(size=14))
            ),
            link=dict(
                source=source_indices,
                target=target_indices,
                value=values,
                color=link_colors,
                hovertemplate='<b>%{source.label}</b> → <b>%{target.label}</b><br>Count: %{value}<extra></extra>',
                hoverlabel=dict(font=dict(size=13))
            )
        )])
        
        # 设置标题（根据版本不同）
        title_suffix = " (excluding intergenic)" if version_name == "without_intergenic" else ""
        fig.update_layout(
            title=dict(
                text=f"<b>Blood Pressure HyPrColoc Sankey — {MODEL_NAME}{title_suffix}</b>",
                font=dict(size=24, color='#000000'),
                x=0.5,
                xanchor='center'
            ),
            font=dict(size=20, family='Arial, sans-serif', color='#000000'),
            plot_bgcolor='white',
            paper_bgcolor='white',
            width=FIGURE_CONFIG["width"],
            height=fig_height,
            margin=dict(t=60, b=30, l=40, r=250)
        )
        
        # 保存图表
        print("保存桑吉图...")
        OUTPUT_FORMATS = {
            "pdf": {"scale": 2, "description": "For publication"},
            "png": {"scale": 2, "description": "For presentations"},
            "svg": {"scale": 1, "description": "Vector format"}
        }
        
        # 文件名后缀
        file_suffix = "_no_intergenic" if version_name == "without_intergenic" else ""
        
        saved_files = []
        for fmt, config in OUTPUT_FORMATS.items():
            filename = os.path.join(OUTPUT_DIR, f"HyPrColoc_BP_Sankey_{MODEL_NAME}{file_suffix}.{fmt}")
            fig.write_image(filename, width=FIGURE_CONFIG["width"], height=fig_height,
                           scale=config.get("scale", 1))
            saved_files.append((fmt.upper(), os.path.basename(filename)))
            print(f"  已保存 {fmt.upper()}: {os.path.basename(filename)}")
        
        print(f"\n{version_name} 完成！")
    
    print(f"\n{MODEL_NAME} 所有版本完成！保存位置: {OUTPUT_DIR}")

print("\n" + "="*70)
print("✓ 所有模型的桑吉图创建完成！")
print("="*70 + "\n")
