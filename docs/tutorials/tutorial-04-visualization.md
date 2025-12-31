# 教程 04：序列可视化

## 📋 本章导览

- **你将掌握的技能**：
  - 使用 matplotlib 进行科学可视化
  - 绘制 GC 含量分布图
  - 绘制碱基/氨基酸组成饼图和条形图
  - 绘制疏水性分布图
  - 生成序列 Logo 图
  - 创建交互式可视化

- **前置知识**：
  - 完成教程 01-03
  - 基本的 matplotlib 使用经验（可选）

- **核心挑战**：如何让抽象的序列数据变得直观可理解？

---

## 📚 理解序列可视化

### 为什么需要可视化？

你可能会想："我已经能计算 GC 含量了，为什么还要画图？"

| 场景 | 纯数字 | 可视化 |
|------|--------|--------|
| 单条序列 | 足够 | 锦上添花 |
| 比较多条序列 | 表格难读 | 一目了然 |
| 发现模式 | 几乎不可能 | 直观 |
| 报告/论文 | 不够专业 | 标准做法 |

### 常见的生物序列可视化

1. **GC 含量滑动窗口图**：显示 GC 含量沿序列的变化
2. **碱基/氨基酸组成图**：饼图或条形图
3. **疏水性分布图**：预测跨膜区域
4. **序列 Logo**：显示保守位点
5. **点阵图**：序列相似性比较

---

## 🔧 设计我们的实现

### 技术选型

| 库 | 优点 | 缺点 | 适用场景 |
|---|------|------|---------|
| matplotlib | 功能全面，标准库 | API 较繁琐 | 静态图，论文 |
| seaborn | 美观，简洁 API | 依赖 matplotlib | 统计图 |
| plotly | 交互式，Web 友好 | 文件较大 | 仪表板 |
| Rich | 终端可视化 | 功能有限 | CLI 工具 |

**我们的选择**：
- 主要使用 matplotlib（通用性强）
- 可选使用 plotly（交互式）
- 终端预览使用 Rich

### 模块设计

```
genomeflow/
├── visualization/
│   ├── __init__.py
│   ├── base.py         # 基础绑图工具
│   ├── composition.py  # 组成分析图
│   ├── profile.py      # 分布图（GC, 疏水性）
│   └── logo.py         # 序列 Logo
```

### 设计原则

1. **分离数据和展示**：计算逻辑在 analyzer 中，可视化只负责展示
2. **合理的默认值**：开箱即用，但可定制
3. **输出灵活**：支持显示、保存、返回 Figure 对象

---

## 💻 代码实现

### 步骤 1：添加依赖

更新 `pyproject.toml`：

```toml
[project]
name = "genomeflow"
version = "0.1.0"
description = "A DNA/RNA/Protein sequence analyzer"
requires-python = ">=3.12"
dependencies = [
    "matplotlib>=3.8",
]

[project.optional-dependencies]
dev = [
    "pytest>=8.0",
    "pytest-cov>=4.0",
]
interactive = [
    "plotly>=5.18",
]
all = [
    "matplotlib>=3.8",
    "plotly>=5.18",
]
```

安装依赖：

```bash
uv add matplotlib
# 可选：交互式可视化
uv add --optional interactive plotly
```

### 步骤 2：创建基础可视化模块

创建 `src/genomeflow/visualization/__init__.py`：

```python
"""
序列可视化模块。

提供各种序列分析结果的可视化功能。
"""

from genomeflow.visualization.composition import (
    plot_base_composition,
    plot_amino_acid_composition,
)
from genomeflow.visualization.profile import (
    plot_gc_content,
    plot_hydrophobicity,
)

__all__ = [
    "plot_base_composition",
    "plot_amino_acid_composition",
    "plot_gc_content",
    "plot_hydrophobicity",
]
```

创建 `src/genomeflow/visualization/base.py`：

```python
"""
可视化的基础工具和样式配置。

这个模块提供：
- 统一的颜色方案
- 图表样式配置
- 通用的辅助函数
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes


# 碱基颜色方案（常见的生物信息学配色）
BASE_COLORS = {
    "A": "#2ecc71",  # 绿色
    "T": "#e74c3c",  # 红色
    "U": "#e74c3c",  # 红色（同 T）
    "G": "#f39c12",  # 橙色
    "C": "#3498db",  # 蓝色
}

# 氨基酸颜色方案（按性质分组）
AMINO_ACID_COLORS = {
    # 疏水（棕色系）
    "A": "#8B4513", "V": "#A0522D", "I": "#D2691E",
    "L": "#CD853F", "M": "#DEB887", "F": "#F4A460",
    "W": "#BC8F8F", "P": "#D2B48C",
    # 亲水（蓝色系）
    "S": "#4169E1", "T": "#6495ED", "N": "#00CED1",
    "Q": "#20B2AA",
    # 酸性（红色）
    "D": "#DC143C", "E": "#FF6347",
    # 碱性（蓝紫色）
    "K": "#8A2BE2", "R": "#9370DB", "H": "#BA55D3",
    # 特殊（灰色/黄色）
    "G": "#808080", "C": "#FFD700", "Y": "#DAA520",
}

# 图表样式配置
STYLE_CONFIG = {
    "figure.figsize": (10, 6),
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "font.family": "sans-serif",
}


def apply_style() -> None:
    """应用统一的图表样式。"""
    import matplotlib.pyplot as plt
    plt.rcParams.update(STYLE_CONFIG)


def save_or_show(
    fig: "Figure",
    output_path: str | None = None,
    dpi: int = 150,
    show: bool = True,
) -> "Figure":
    """
    保存或显示图表。

    Args:
        fig: matplotlib Figure 对象
        output_path: 保存路径（如果为 None 则显示）
        dpi: 图像分辨率
        show: 是否显示图表

    Returns:
        Figure 对象（便于链式调用）
    """
    import matplotlib.pyplot as plt

    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        print(f"图表已保存到: {output_path}")

    if show and not output_path:
        plt.show()

    return fig
```

### 步骤 3：实现组成分析可视化

创建 `src/genomeflow/visualization/composition.py`：

```python
"""
序列组成的可视化。

包括：
- 碱基组成（DNA/RNA）
- 氨基酸组成（蛋白质）
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import matplotlib.pyplot as plt

from genomeflow.visualization.base import (
    BASE_COLORS,
    AMINO_ACID_COLORS,
    apply_style,
    save_or_show,
)

if TYPE_CHECKING:
    from matplotlib.figure import Figure
    from genomeflow.sequence import DNASequence
    from genomeflow.rna import RNASequence
    from genomeflow.protein import ProteinSequence


def plot_base_composition(
    sequence: "DNASequence | RNASequence",
    plot_type: Literal["bar", "pie"] = "bar",
    title: str | None = None,
    output_path: str | None = None,
    show: bool = True,
) -> "Figure":
    """
    绘制碱基组成图。

    Args:
        sequence: DNA 或 RNA 序列
        plot_type: 图表类型，"bar"（条形图）或 "pie"（饼图）
        title: 图表标题
        output_path: 保存路径（可选）
        show: 是否显示图表

    Returns:
        matplotlib Figure 对象

    示例：
        >>> from genomeflow import DNASequence
        >>> from genomeflow.visualization import plot_base_composition
        >>> seq = DNASequence("ATGCGATCGATCGATCG")
        >>> plot_base_composition(seq, plot_type="pie")
    """
    apply_style()

    # 统计碱基
    from collections import Counter
    counts = Counter(sequence)

    # 确定碱基顺序
    if "U" in counts:
        bases = ["A", "U", "G", "C"]
        seq_type = "RNA"
    else:
        bases = ["A", "T", "G", "C"]
        seq_type = "DNA"

    values = [counts.get(base, 0) for base in bases]
    colors = [BASE_COLORS[base] for base in bases]

    # 创建图表
    fig, ax = plt.subplots(figsize=(8, 6))

    if plot_type == "bar":
        bars = ax.bar(bases, values, color=colors, edgecolor="white", linewidth=1.5)

        # 在条形上添加数值
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax.annotate(
                f"{value}\n({value/sum(values):.1%})",
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=10,
            )

        ax.set_ylabel("计数")
        ax.set_xlabel("碱基")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    elif plot_type == "pie":
        # 饼图
        wedges, texts, autotexts = ax.pie(
            values,
            labels=bases,
            colors=colors,
            autopct="%1.1f%%",
            startangle=90,
            explode=[0.02] * len(bases),
        )

        # 美化文本
        for autotext in autotexts:
            autotext.set_fontsize(11)
            autotext.set_fontweight("bold")

    # 设置标题
    if title is None:
        title = f"{seq_type} 碱基组成 (n={len(sequence)})"
    ax.set_title(title, fontsize=14, fontweight="bold")

    plt.tight_layout()
    return save_or_show(fig, output_path, show=show)


def plot_amino_acid_composition(
    sequence: "ProteinSequence",
    top_n: int | None = None,
    plot_type: Literal["bar", "pie"] = "bar",
    group_by_property: bool = False,
    title: str | None = None,
    output_path: str | None = None,
    show: bool = True,
) -> "Figure":
    """
    绘制氨基酸组成图。

    Args:
        sequence: 蛋白质序列
        top_n: 只显示前 N 个最多的氨基酸（None 表示全部）
        plot_type: 图表类型
        group_by_property: 是否按性质分组显示
        title: 图表标题
        output_path: 保存路径
        show: 是否显示

    Returns:
        matplotlib Figure 对象

    示例：
        >>> from genomeflow import ProteinSequence
        >>> from genomeflow.visualization import plot_amino_acid_composition
        >>> protein = ProteinSequence("MKFLILLFNILCLFPVLAADNH")
        >>> plot_amino_acid_composition(protein, top_n=10)
    """
    apply_style()

    # 获取组成
    comp = sequence.amino_acid_composition()

    if group_by_property:
        return _plot_grouped_composition(comp, title, output_path, show)

    # 排序
    sorted_items = sorted(comp.counts.items(), key=lambda x: x[1], reverse=True)

    if top_n:
        sorted_items = sorted_items[:top_n]

    amino_acids = [item[0] for item in sorted_items]
    values = [item[1] for item in sorted_items]
    colors = [AMINO_ACID_COLORS.get(aa, "#808080") for aa in amino_acids]

    # 创建图表
    fig, ax = plt.subplots(figsize=(12, 6))

    if plot_type == "bar":
        bars = ax.bar(amino_acids, values, color=colors, edgecolor="white")

        # 添加数值标签
        for bar, value in zip(bars, values):
            ax.annotate(
                str(value),
                xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                fontsize=9,
            )

        ax.set_ylabel("计数")
        ax.set_xlabel("氨基酸")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    elif plot_type == "pie":
        ax.pie(
            values,
            labels=amino_acids,
            colors=colors,
            autopct="%1.1f%%",
            startangle=90,
        )

    # 标题
    if title is None:
        suffix = f"（前 {top_n} 位）" if top_n else ""
        title = f"氨基酸组成{suffix} (n={comp.total})"
    ax.set_title(title, fontsize=14, fontweight="bold")

    plt.tight_layout()
    return save_or_show(fig, output_path, show=show)


def _plot_grouped_composition(comp, title, output_path, show):
    """按性质分组绘制氨基酸组成。"""
    # 氨基酸分组
    groups = {
        "疏水": ["A", "V", "I", "L", "M", "F", "W", "P"],
        "亲水": ["S", "T", "N", "Q"],
        "酸性": ["D", "E"],
        "碱性": ["K", "R", "H"],
        "特殊": ["G", "C", "Y"],
    }

    group_colors = {
        "疏水": "#D2691E",
        "亲水": "#4169E1",
        "酸性": "#DC143C",
        "碱性": "#8A2BE2",
        "特殊": "#808080",
    }

    # 计算每组的总数
    group_counts = {}
    for group_name, amino_acids in groups.items():
        total = sum(comp.counts.get(aa, 0) for aa in amino_acids)
        group_counts[group_name] = total

    fig, ax = plt.subplots(figsize=(10, 6))

    names = list(group_counts.keys())
    values = list(group_counts.values())
    colors = [group_colors[name] for name in names]

    bars = ax.bar(names, values, color=colors, edgecolor="white", linewidth=2)

    # 添加百分比
    total = sum(values)
    for bar, value in zip(bars, values):
        pct = value / total * 100 if total > 0 else 0
        ax.annotate(
            f"{value}\n({pct:.1f}%)",
            xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            fontsize=11,
        )

    ax.set_ylabel("氨基酸数量")
    ax.set_xlabel("氨基酸性质")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if title is None:
        title = f"氨基酸性质分布 (n={comp.total})"
    ax.set_title(title, fontsize=14, fontweight="bold")

    plt.tight_layout()
    return save_or_show(fig, output_path, show=show)
```

### 步骤 4：实现分布图可视化

创建 `src/genomeflow/visualization/profile.py`：

```python
"""
序列特征分布的可视化。

包括：
- GC 含量滑动窗口图
- 疏水性分布图
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from genomeflow.visualization.base import apply_style, save_or_show

if TYPE_CHECKING:
    from matplotlib.figure import Figure
    from genomeflow.sequence import DNASequence
    from genomeflow.rna import RNASequence
    from genomeflow.protein import ProteinSequence


def plot_gc_content(
    sequence: "DNASequence | RNASequence",
    window_size: int = 100,
    step: int = 10,
    title: str | None = None,
    output_path: str | None = None,
    show: bool = True,
) -> "Figure":
    """
    绘制 GC 含量滑动窗口分布图。

    这种图可以帮助识别：
    - 基因组中 GC 含量的变化区域
    - 可能的水平基因转移区域（GC 含量异常）
    - 编码区和非编码区的边界

    Args:
        sequence: DNA 或 RNA 序列
        window_size: 滑动窗口大小
        step: 滑动步长
        title: 图表标题
        output_path: 保存路径
        show: 是否显示

    Returns:
        matplotlib Figure 对象

    示例：
        >>> from genomeflow import DNASequence
        >>> from genomeflow.visualization import plot_gc_content
        >>> seq = DNASequence("ATGC" * 100 + "GGCC" * 50 + "ATAT" * 100)
        >>> plot_gc_content(seq, window_size=20, step=5)
    """
    apply_style()

    seq_str = sequence.sequence

    # 计算滑动窗口 GC 含量
    positions: list[int] = []
    gc_values: list[float] = []

    for i in range(0, len(seq_str) - window_size + 1, step):
        window = seq_str[i:i + window_size]
        gc = (window.count("G") + window.count("C")) / window_size
        positions.append(i + window_size // 2)  # 窗口中心位置
        gc_values.append(gc * 100)  # 转换为百分比

    # 创建图表
    fig, ax = plt.subplots(figsize=(12, 5))

    # 绘制 GC 曲线
    ax.plot(positions, gc_values, color="#3498db", linewidth=1.5, label="GC 含量")

    # 添加平均线
    avg_gc = sum(gc_values) / len(gc_values) if gc_values else 0
    ax.axhline(
        y=avg_gc,
        color="#e74c3c",
        linestyle="--",
        linewidth=1.5,
        label=f"平均值 ({avg_gc:.1f}%)",
    )

    # 添加 50% 参考线
    ax.axhline(y=50, color="#95a5a6", linestyle=":", linewidth=1, alpha=0.7)

    # 填充区域
    ax.fill_between(positions, gc_values, alpha=0.3, color="#3498db")

    # 设置坐标轴
    ax.set_xlabel("序列位置 (bp)")
    ax.set_ylabel("GC 含量 (%)")
    ax.set_ylim(0, 100)
    ax.set_xlim(0, len(seq_str))

    # 添加网格
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    # 标题
    if title is None:
        title = f"GC 含量分布 (窗口={window_size}bp, 序列长度={len(seq_str)}bp)"
    ax.set_title(title, fontsize=14, fontweight="bold")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    return save_or_show(fig, output_path, show=show)


def plot_hydrophobicity(
    sequence: "ProteinSequence",
    window_size: int = 9,
    title: str | None = None,
    output_path: str | None = None,
    show: bool = True,
) -> "Figure":
    """
    绘制疏水性分布图（Kyte-Doolittle 图）。

    这种图用于：
    - 预测跨膜区域（疏水性高的区域）
    - 识别信号肽
    - 分析蛋白质结构

    Args:
        sequence: 蛋白质序列
        window_size: 滑动窗口大小（推荐 7-11）
        title: 图表标题
        output_path: 保存路径
        show: 是否显示

    Returns:
        matplotlib Figure 对象

    示例：
        >>> from genomeflow import ProteinSequence
        >>> from genomeflow.visualization import plot_hydrophobicity
        >>> # 一个包含跨膜区的假设序列
        >>> protein = ProteinSequence("MKFLILLFNILCLFPVLAADNHEKK")
        >>> plot_hydrophobicity(protein)
    """
    apply_style()

    # 获取疏水性分布
    profile = sequence.hydrophobicity_profile(window=window_size)

    if not profile:
        print("警告：序列太短，无法计算疏水性分布")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "序列太短", ha="center", va="center")
        return fig

    # 计算位置（窗口中心）
    half_window = window_size // 2
    positions = list(range(half_window, half_window + len(profile)))

    # 创建图表
    fig, ax = plt.subplots(figsize=(12, 5))

    # 绘制疏水性曲线
    ax.plot(positions, profile, color="#2c3e50", linewidth=1.5)

    # 填充正负区域
    profile_array = np.array(profile)
    positions_array = np.array(positions)

    # 疏水区域（正值）- 红色
    ax.fill_between(
        positions_array,
        profile_array,
        0,
        where=(profile_array > 0),
        color="#e74c3c",
        alpha=0.4,
        label="疏水",
    )

    # 亲水区域（负值）- 蓝色
    ax.fill_between(
        positions_array,
        profile_array,
        0,
        where=(profile_array < 0),
        color="#3498db",
        alpha=0.4,
        label="亲水",
    )

    # 零线
    ax.axhline(y=0, color="#2c3e50", linestyle="-", linewidth=0.5)

    # 跨膜阈值线（通常 1.6）
    ax.axhline(
        y=1.6,
        color="#e74c3c",
        linestyle="--",
        linewidth=1,
        alpha=0.7,
        label="跨膜阈值 (1.6)",
    )

    # 设置坐标轴
    ax.set_xlabel("氨基酸位置")
    ax.set_ylabel("疏水性 (Kyte-Doolittle)")
    ax.set_xlim(0, len(sequence))

    # 网格和图例
    ax.grid(True, alpha=0.3, axis="y")
    ax.legend(loc="upper right")

    # 标题
    if title is None:
        title = f"疏水性分布 (窗口={window_size}, GRAVY={sequence.gravy():.2f})"
    ax.set_title(title, fontsize=14, fontweight="bold")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    return save_or_show(fig, output_path, show=show)


def plot_sequence_overview(
    dna_sequence: "DNASequence",
    output_path: str | None = None,
    show: bool = True,
) -> "Figure":
    """
    绘制序列综合分析图（多子图）。

    包含：
    - GC 含量分布
    - 碱基组成
    - 基本统计信息

    Args:
        dna_sequence: DNA 序列
        output_path: 保存路径
        show: 是否显示

    Returns:
        matplotlib Figure 对象
    """
    apply_style()

    from collections import Counter

    fig = plt.figure(figsize=(14, 8))

    # 创建网格布局
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

    # 1. 碱基组成饼图
    ax1 = fig.add_subplot(gs[0, 0])
    counts = Counter(dna_sequence)
    bases = ["A", "T", "G", "C"]
    values = [counts.get(b, 0) for b in bases]
    colors = ["#2ecc71", "#e74c3c", "#f39c12", "#3498db"]
    ax1.pie(values, labels=bases, colors=colors, autopct="%1.1f%%", startangle=90)
    ax1.set_title("碱基组成")

    # 2. 碱基条形图
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.bar(bases, values, color=colors)
    ax2.set_ylabel("计数")
    ax2.set_title("碱基分布")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # 3. 统计信息文本
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.axis("off")

    gc = dna_sequence.gc_content()
    stats_text = f"""
序列统计信息
━━━━━━━━━━━━━━━━━━━━━━
长度:     {len(dna_sequence):,} bp
GC 含量:  {gc:.2%}

碱基计数:
  A: {counts.get('A', 0):,}
  T: {counts.get('T', 0):,}
  G: {counts.get('G', 0):,}
  C: {counts.get('C', 0):,}

A+T: {counts.get('A', 0) + counts.get('T', 0):,} ({(counts.get('A', 0) + counts.get('T', 0))/len(dna_sequence):.1%})
G+C: {counts.get('G', 0) + counts.get('C', 0):,} ({gc:.1%})
"""
    ax3.text(
        0.1, 0.9, stats_text,
        transform=ax3.transAxes,
        fontsize=11,
        verticalalignment="top",
        fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    # 4. GC 含量分布（下方横跨）
    ax4 = fig.add_subplot(gs[1, :])

    window_size = min(100, len(dna_sequence) // 5)
    step = max(1, window_size // 10)

    if window_size >= 10:
        seq_str = dna_sequence.sequence
        positions = []
        gc_values = []

        for i in range(0, len(seq_str) - window_size + 1, step):
            window = seq_str[i:i + window_size]
            gc = (window.count("G") + window.count("C")) / window_size * 100
            positions.append(i + window_size // 2)
            gc_values.append(gc)

        ax4.plot(positions, gc_values, color="#3498db", linewidth=1)
        ax4.fill_between(positions, gc_values, alpha=0.3, color="#3498db")
        ax4.axhline(y=gc * 100, color="#e74c3c", linestyle="--", label=f"平均 {gc*100:.1f}%")
        ax4.set_ylim(0, 100)
        ax4.set_xlim(0, len(seq_str))
        ax4.legend()

    ax4.set_xlabel("序列位置 (bp)")
    ax4.set_ylabel("GC 含量 (%)")
    ax4.set_title(f"GC 含量分布 (窗口大小={window_size})")
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)
    ax4.grid(True, alpha=0.3)

    plt.suptitle(
        f"DNA 序列分析报告",
        fontsize=16,
        fontweight="bold",
        y=1.02,
    )

    plt.tight_layout()
    return save_or_show(fig, output_path, show=show)
```

### 步骤 5：更新包导出

更新 `src/genomeflow/__init__.py`，添加可视化模块：

```python
# ... 现有导入 ...

# 可视化（可选导入，避免强制依赖 matplotlib）
try:
    from genomeflow.visualization import (
        plot_base_composition,
        plot_amino_acid_composition,
        plot_gc_content,
        plot_hydrophobicity,
    )
    _HAS_VISUALIZATION = True
except ImportError:
    _HAS_VISUALIZATION = False

# 更新 __all__
if _HAS_VISUALIZATION:
    __all__ += [
        "plot_base_composition",
        "plot_amino_acid_composition",
        "plot_gc_content",
        "plot_hydrophobicity",
    ]
```

### 步骤 6：编写测试

创建 `tests/test_visualization.py`：

```python
"""可视化模块的测试。"""

import pytest

# 跳过如果没有安装 matplotlib
pytest.importorskip("matplotlib")

from genomeflow import DNASequence, ProteinSequence
from genomeflow.visualization import (
    plot_base_composition,
    plot_gc_content,
    plot_amino_acid_composition,
    plot_hydrophobicity,
)


class TestBaseComposition:
    """测试碱基组成图。"""

    def test_bar_chart(self):
        seq = DNASequence("ATGCGATCGATCGATCG")
        fig = plot_base_composition(seq, plot_type="bar", show=False)
        assert fig is not None

    def test_pie_chart(self):
        seq = DNASequence("ATGCGATCGATCGATCG")
        fig = plot_base_composition(seq, plot_type="pie", show=False)
        assert fig is not None


class TestGCContent:
    """测试 GC 含量分布图。"""

    def test_gc_plot(self):
        seq = DNASequence("ATGC" * 50)
        fig = plot_gc_content(seq, window_size=20, step=5, show=False)
        assert fig is not None


class TestAminoAcidComposition:
    """测试氨基酸组成图。"""

    def test_bar_chart(self):
        protein = ProteinSequence("MKFLILLFNILCLFPVLAADNH")
        fig = plot_amino_acid_composition(protein, show=False)
        assert fig is not None

    def test_top_n(self):
        protein = ProteinSequence("MKFLILLFNILCLFPVLAADNH")
        fig = plot_amino_acid_composition(protein, top_n=5, show=False)
        assert fig is not None


class TestHydrophobicity:
    """测试疏水性分布图。"""

    def test_hydrophobicity_plot(self):
        protein = ProteinSequence("MKFLILLFNILCLFPVLAADNHEKKERR")
        fig = plot_hydrophobicity(protein, window_size=5, show=False)
        assert fig is not None
```

---

## ✅ 使用示例

### 基本使用

```python
from genomeflow import DNASequence, ProteinSequence
from genomeflow.visualization import (
    plot_base_composition,
    plot_gc_content,
    plot_amino_acid_composition,
    plot_hydrophobicity,
)

# 1. DNA 碱基组成
dna = DNASequence("ATGCGATCGATCGATCGATCGATCGATCG")
plot_base_composition(dna, plot_type="pie")

# 2. GC 含量分布
long_dna = DNASequence("ATGC" * 100 + "GGCC" * 50 + "ATAT" * 100)
plot_gc_content(long_dna, window_size=50)

# 3. 蛋白质氨基酸组成
protein = ProteinSequence("MKFLILLFNILCLFPVLAADNHEKKERR")
plot_amino_acid_composition(protein, group_by_property=True)

# 4. 疏水性分布
plot_hydrophobicity(protein, window_size=7)
```

### 保存图表

```python
# 保存为 PNG
plot_base_composition(dna, output_path="composition.png")

# 保存为 PDF（适合论文）
plot_gc_content(long_dna, output_path="gc_content.pdf")

# 保存为 SVG（可编辑矢量图）
plot_hydrophobicity(protein, output_path="hydrophobicity.svg")
```

### 自定义样式

```python
import matplotlib.pyplot as plt

# 使用不同的样式
plt.style.use("seaborn-v0_8-whitegrid")
plot_gc_content(long_dna)

# 恢复默认
plt.style.use("default")
```

---

## 🤔 深入思考

<details>
<summary>为什么 GC 含量要用滑动窗口？</summary>

**问题**：单一的 GC 含量值丢失了位置信息。

**示例**：
```
序列 A: GGGGCCCCATAT  GC = 66.7%
序列 B: GCGCATGCATAT  GC = 50%
```

滑动窗口可以揭示：
- 哪些区域 GC 富集？
- 是否有外源 DNA 插入（GC 跳变）？
- 编码区和非编码区的边界在哪？

**窗口大小选择**：
- 太小：噪声大
- 太大：细节丢失
- 推荐：50-200 bp（基因组），10-50 bp（短序列）

</details>

<details>
<summary>如何解读疏水性分布图？</summary>

**Kyte-Doolittle 量表**：
- 正值 → 疏水（蛋白质内部/膜区）
- 负值 → 亲水（蛋白质表面）

**跨膜区特征**：
- 连续 15-25 个氨基酸
- 平均疏水性 > 1.6
- 两侧有带电氨基酸（锚定）

**注意**：
- 这只是预测，需要实验验证
- 信号肽也会显示疏水峰
- 专业工具如 TMHMM 更准确

</details>

---

## 📝 总结

通过这个教程，你学会了：

1. **matplotlib 基础**：创建科学图表
2. **碱基/氨基酸组成图**：饼图和条形图
3. **GC 含量分布图**：滑动窗口分析
4. **疏水性分布图**：蛋白质结构预测

### 可视化最佳实践

| 原则 | 说明 |
|------|------|
| 简洁 | 移除不必要的元素 |
| 清晰 | 标签、图例完整 |
| 一致 | 统一的颜色方案 |
| 可复现 | 代码保存，参数记录 |

---

## 🎉 教程系列完结

恭喜你完成了 GenomeFlow 教程系列！

你已经学会了：
- 教程 01：DNA 序列分析基础
- 教程 02：RNA 序列支持
- 教程 03：蛋白质序列支持
- 教程 04：序列可视化

接下来你可以：
- 扩展更多分析功能
- 添加序列比对算法
- 构建 Web 界面
- 贡献到开源社区

**祝你在生物信息学的道路上越走越远！**
