"""
序列特征分布可视化。

提供沿序列位置变化的特征图绑制功能。
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt

if TYPE_CHECKING:
    from genomeflow.protein import ProteinSequence
    from genomeflow.sequence import DNASequence


def plot_hydrophobicity_profile(
    protein: ProteinSequence,
    window: int = 9,
    title: str = "疏水性分布",
    save_path: str | Path | None = None,
    show: bool = True,
) -> plt.Figure:
    """
    绑制蛋白质疏水性分布图。

    使用 Kyte-Doolittle 疏水性量表，通过滑动窗口计算
    每个位置的平均疏水性。

    正值区域可能是：
    - 跨膜区域
    - 蛋白质疏水核心

    负值区域可能是：
    - 蛋白质表面
    - 与水接触的区域

    Args:
        protein: 蛋白质序列对象
        window: 滑动窗口大小，推荐 7-11
        title: 图表标题
        save_path: 保存路径（可选）
        show: 是否显示图表

    Returns:
        matplotlib Figure 对象

    示例：
        >>> from genomeflow.protein import ProteinSequence
        >>> protein = ProteinSequence("MKFLILLFNILCLFPVLAADNH")
        >>> fig = plot_hydrophobicity_profile(protein, show=False)
    """
    # 计算疏水性分布
    profile = protein.hydrophobicity_profile(window=window)

    if not profile:
        raise ValueError(f"序列太短，无法计算 window={window} 的疏水性分布")

    # 创建图表
    fig, ax = plt.subplots(figsize=(12, 5))

    # 位置（从窗口中心开始）
    half_window = window // 2
    positions = list(range(half_window + 1, half_window + len(profile) + 1))

    # 绑制折线
    ax.plot(positions, profile, color="#45B7D1", linewidth=1.5)

    # 填充正负区域
    ax.fill_between(
        positions,
        profile,
        0,
        where=[v > 0 for v in profile],
        color="#FF6B6B",
        alpha=0.3,
        label="疏水区域",
    )
    ax.fill_between(
        positions,
        profile,
        0,
        where=[v <= 0 for v in profile],
        color="#4ECDC4",
        alpha=0.3,
        label="亲水区域",
    )

    # 添加零线
    ax.axhline(y=0, color="black", linestyle="-", linewidth=0.5)

    # 添加跨膜域阈值线（经验值）
    ax.axhline(
        y=1.6,
        color="#FF6B6B",
        linestyle="--",
        linewidth=1,
        alpha=0.7,
        label="跨膜阈值 (1.6)",
    )

    # 设置标签
    ax.set_xlabel("氨基酸位置", fontsize=12)
    ax.set_ylabel(f"疏水性 (窗口={window})", fontsize=12)
    ax.set_title(title, fontsize=16, fontweight="bold")
    ax.legend(loc="upper right")

    # 添加网格
    ax.grid(True, alpha=0.3, linestyle="--")

    plt.tight_layout()

    # 保存图表
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_gc_content_window(
    seq: DNASequence,
    window: int = 100,
    step: int = 10,
    title: str = "GC 含量分布",
    save_path: str | Path | None = None,
    show: bool = True,
) -> plt.Figure:
    """
    绑制沿序列位置的 GC 含量变化图。

    使用滑动窗口计算每个位置的局部 GC 含量，
    可以揭示基因组的区域特征。

    Args:
        seq: DNA 序列对象
        window: 滑动窗口大小
        step: 滑动步长
        title: 图表标题
        save_path: 保存路径（可选）
        show: 是否显示图表

    Returns:
        matplotlib Figure 对象

    示例：
        >>> from genomeflow.sequence import DNASequence
        >>> seq = DNASequence("ATGC" * 100)
        >>> fig = plot_gc_content_window(seq, window=20, show=False)
    """
    seq_str = seq.sequence

    if len(seq_str) < window:
        raise ValueError(f"序列长度 ({len(seq_str)}) 小于窗口大小 ({window})")

    # 计算滑动窗口 GC 含量
    positions: list[int] = []
    gc_values: list[float] = []

    for i in range(0, len(seq_str) - window + 1, step):
        window_seq = seq_str[i : i + window]
        gc_count = window_seq.count("G") + window_seq.count("C")
        gc_ratio = gc_count / window * 100

        positions.append(i + window // 2)  # 使用窗口中心位置
        gc_values.append(gc_ratio)

    # 创建图表
    fig, ax = plt.subplots(figsize=(12, 5))

    # 绑制折线
    ax.plot(positions, gc_values, color="#45B7D1", linewidth=1.5)

    # 填充区域
    ax.fill_between(positions, gc_values, alpha=0.3, color="#45B7D1")

    # 添加平均线
    mean_gc = sum(gc_values) / len(gc_values)
    ax.axhline(
        y=mean_gc,
        color="#FF6B6B",
        linestyle="--",
        linewidth=2,
        label=f"平均值: {mean_gc:.1f}%",
    )

    # 添加参考线
    ax.axhline(y=50, color="gray", linestyle=":", linewidth=1, alpha=0.5)

    # 设置标签
    ax.set_xlabel("序列位置 (bp)", fontsize=12)
    ax.set_ylabel(f"GC 含量 % (窗口={window}bp)", fontsize=12)
    ax.set_title(title, fontsize=16, fontweight="bold")
    ax.legend(loc="upper right")

    # 设置 y 轴范围
    ax.set_ylim(0, 100)

    # 添加网格
    ax.grid(True, alpha=0.3, linestyle="--")

    plt.tight_layout()

    # 保存图表
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_codon_usage(
    codons: dict[str, int],
    title: str = "密码子使用频率",
    save_path: str | Path | None = None,
    show: bool = True,
) -> plt.Figure:
    """
    绑制密码子使用频率热力图。

    Args:
        codons: 密码子计数字典 {密码子: 计数}
        title: 图表标题
        save_path: 保存路径（可选）
        show: 是否显示图表

    Returns:
        matplotlib Figure 对象

    示例：
        >>> codons = {"AUG": 10, "UUU": 5, "UUC": 8}
        >>> fig = plot_codon_usage(codons, show=False)
    """
    import numpy as np

    # 密码子按第一、二位碱基分组
    bases = ["U", "C", "A", "G"]

    # 创建 16x4 的矩阵（16 个二碱基组合 x 4 个第三位碱基）
    # 行: 第一位 + 第二位 (UU, UC, UA, UG, CU, ...)
    # 列: 第三位 (U, C, A, G)

    data = np.zeros((16, 4))
    row_labels = []
    col_labels = bases

    idx = 0
    for b1 in bases:
        for b2 in bases:
            row_labels.append(f"{b1}{b2}")
            for j, b3 in enumerate(bases):
                codon = f"{b1}{b2}{b3}"
                data[idx, j] = codons.get(codon, 0)
            idx += 1

    # 创建图表
    fig, ax = plt.subplots(figsize=(8, 12))

    # 绑制热力图
    im = ax.imshow(data, cmap="YlOrRd", aspect="auto")

    # 设置刻度
    ax.set_xticks(range(4))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(range(16))
    ax.set_yticklabels(row_labels)

    # 添加数值标注
    for i in range(16):
        for j in range(4):
            value = int(data[i, j])
            if value > 0:
                text_color = "white" if value > data.max() / 2 else "black"
                ax.text(j, i, str(value), ha="center", va="center", color=text_color)

    # 设置标签
    ax.set_xlabel("第三位碱基", fontsize=12)
    ax.set_ylabel("第一、二位碱基", fontsize=12)
    ax.set_title(title, fontsize=16, fontweight="bold")

    # 添加颜色条
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("计数", fontsize=11)

    plt.tight_layout()

    # 保存图表
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()

    return fig
