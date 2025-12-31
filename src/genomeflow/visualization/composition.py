"""
序列组成可视化。

提供碱基/氨基酸组成的图表绑制功能。
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt

if TYPE_CHECKING:
    from genomeflow.sequence import DNASequence
    from genomeflow.io import FastaRecord


def plot_base_composition(
    seq: "DNASequence",
    title: str = "碱基组成",
    save_path: str | Path | None = None,
    show: bool = True,
) -> plt.Figure:
    """
    绑制 DNA 序列的碱基组成饼图。

    Args:
        seq: DNA 序列对象
        title: 图表标题
        save_path: 保存路径（可选）
        show: 是否显示图表

    Returns:
        matplotlib Figure 对象

    示例：
        >>> from genomeflow.sequence import DNASequence
        >>> seq = DNASequence("ATGCGATCGATCG")
        >>> fig = plot_base_composition(seq, show=False)
    """
    # 统计碱基
    counts = Counter(seq)

    # 准备数据
    bases = ["A", "T", "G", "C"]
    values = [counts.get(base, 0) for base in bases]
    colors = ["#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4"]

    # 创建图表
    fig, ax = plt.subplots(figsize=(8, 8))

    # 绑制饼图
    wedges, texts, autotexts = ax.pie(
        values,
        labels=bases,
        autopct="%1.1f%%",
        colors=colors,
        explode=(0.02, 0.02, 0.02, 0.02),
        shadow=True,
        startangle=90,
    )

    # 美化文字
    for text in texts:
        text.set_fontsize(14)
        text.set_fontweight("bold")
    for autotext in autotexts:
        autotext.set_fontsize(12)
        autotext.set_color("white")

    ax.set_title(title, fontsize=16, fontweight="bold", pad=20)

    # 添加图例
    legend_labels = [f"{base}: {count}" for base, count in zip(bases, values)]
    ax.legend(
        wedges,
        legend_labels,
        title="碱基计数",
        loc="center left",
        bbox_to_anchor=(1, 0, 0.5, 1),
    )

    plt.tight_layout()

    # 保存图表
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_gc_distribution(
    records: list["FastaRecord"],
    title: str = "GC 含量分布",
    save_path: str | Path | None = None,
    show: bool = True,
) -> plt.Figure:
    """
    绑制多条序列的 GC 含量分布直方图。

    Args:
        records: FastaRecord 列表
        title: 图表标题
        save_path: 保存路径（可选）
        show: 是否显示图表

    Returns:
        matplotlib Figure 对象

    示例：
        >>> from genomeflow.io import read_fasta
        >>> records = list(read_fasta("sequences.fasta"))
        >>> fig = plot_gc_distribution(records, show=False)
    """
    from genomeflow.analyzer import gc_content

    # 计算每条序列的 GC 含量
    gc_values = [gc_content(record.sequence) * 100 for record in records]

    # 创建图表
    fig, ax = plt.subplots(figsize=(10, 6))

    # 绑制直方图
    n, bins, patches = ax.hist(
        gc_values,
        bins=20,
        color="#45B7D1",
        edgecolor="white",
        alpha=0.8,
    )

    # 添加平均线
    mean_gc = sum(gc_values) / len(gc_values) if gc_values else 0
    ax.axvline(
        mean_gc,
        color="#FF6B6B",
        linestyle="--",
        linewidth=2,
        label=f"平均值: {mean_gc:.1f}%",
    )

    # 设置标签
    ax.set_xlabel("GC 含量 (%)", fontsize=12)
    ax.set_ylabel("序列数量", fontsize=12)
    ax.set_title(title, fontsize=16, fontweight="bold")
    ax.legend(fontsize=11)

    # 设置 x 轴范围
    ax.set_xlim(0, 100)

    # 添加网格
    ax.grid(True, alpha=0.3, linestyle="--")

    plt.tight_layout()

    # 保存图表
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_amino_acid_composition(
    protein_seq: str,
    title: str = "氨基酸组成",
    save_path: str | Path | None = None,
    show: bool = True,
) -> plt.Figure:
    """
    绑制蛋白质序列的氨基酸组成条形图。

    Args:
        protein_seq: 蛋白质序列字符串
        title: 图表标题
        save_path: 保存路径（可选）
        show: 是否显示图表

    Returns:
        matplotlib Figure 对象

    示例：
        >>> fig = plot_amino_acid_composition("MKFLILLFNILCLFPVLAADNH", show=False)
    """
    # 统计氨基酸
    counts = Counter(protein_seq.upper())

    # 按字母顺序排列
    amino_acids = sorted(counts.keys())
    values = [counts[aa] for aa in amino_acids]

    # 创建图表
    fig, ax = plt.subplots(figsize=(12, 6))

    # 根据氨基酸类型着色
    # 疏水: 蓝色, 极性: 绿色, 酸性: 红色, 碱性: 紫色
    colors = []
    hydrophobic = set("AILMFWVY")
    polar = set("STNQCGP")
    acidic = set("DE")
    basic = set("RKH")

    for aa in amino_acids:
        if aa in hydrophobic:
            colors.append("#4ECDC4")  # 青色
        elif aa in polar:
            colors.append("#96CEB4")  # 绿色
        elif aa in acidic:
            colors.append("#FF6B6B")  # 红色
        elif aa in basic:
            colors.append("#9B59B6")  # 紫色
        else:
            colors.append("#95A5A6")  # 灰色

    # 绑制条形图
    bars = ax.bar(amino_acids, values, color=colors, edgecolor="white")

    # 在条形上方添加数值
    for bar, value in zip(bars, values):
        if value > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.1,
                str(value),
                ha="center",
                va="bottom",
                fontsize=9,
            )

    # 设置标签
    ax.set_xlabel("氨基酸", fontsize=12)
    ax.set_ylabel("数量", fontsize=12)
    ax.set_title(title, fontsize=16, fontweight="bold")

    # 添加图例
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor="#4ECDC4", label="疏水性"),
        Patch(facecolor="#96CEB4", label="极性"),
        Patch(facecolor="#FF6B6B", label="酸性"),
        Patch(facecolor="#9B59B6", label="碱性"),
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    plt.tight_layout()

    # 保存图表
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()

    return fig
