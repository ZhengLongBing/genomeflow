"""
DNA 序列分析功能。

设计决策：为什么用独立函数而不是 DNASequence 的方法？
1. 单一职责：DNASequence 负责数据，analyzer 负责分析
2. 可测试性：函数更容易单独测试
3. 可扩展性：添加新分析功能不需要修改 DNASequence
4. 函数式风格：更容易组合和复用
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from genomeflow.sequence import DNASequence


@dataclass(frozen=True)
class BaseFrequency:
    """
    碱基频率统计结果。

    使用 dataclass 的好处：
    1. 自动生成 __init__, __repr__, __eq__
    2. frozen=True 使其不可变，更安全
    3. 作为返回类型比字典更清晰
    """

    a_count: int
    t_count: int
    g_count: int
    c_count: int
    total: int

    @property
    def a_ratio(self) -> float:
        """A 的比例。"""
        return self.a_count / self.total if self.total > 0 else 0.0

    @property
    def t_ratio(self) -> float:
        """T 的比例。"""
        return self.t_count / self.total if self.total > 0 else 0.0

    @property
    def g_ratio(self) -> float:
        """G 的比例。"""
        return self.g_count / self.total if self.total > 0 else 0.0

    @property
    def c_ratio(self) -> float:
        """C 的比例。"""
        return self.c_count / self.total if self.total > 0 else 0.0


def gc_content(seq: DNASequence) -> float:
    """
    计算 GC 含量（G 和 C 碱基的百分比）。

    GC 含量是生物信息学中最常用的指标之一：
    - 物种鉴定：不同物种有特征性的 GC 含量
    - 基因预测：编码区和非编码区 GC 含量不同
    - 实验设计：影响 DNA 熔解温度和 PCR 条件

    Args:
        seq: DNA 序列

    Returns:
        GC 含量，范围 0.0-1.0（0%-100%）
        空序列返回 0.0

    示例：
        >>> gc_content(DNASequence("ATGC"))
        0.5
        >>> gc_content(DNASequence("GGCC"))
        1.0
    """
    if len(seq) == 0:
        return 0.0

    # 使用 Counter 统计，比手动循环更 Pythonic
    # Counter 返回一个字典，如 {'A': 2, 'T': 2, 'G': 1, 'C': 1}
    counts = Counter(seq)
    gc_count = counts.get("G", 0) + counts.get("C", 0)

    return gc_count / len(seq)


def base_frequency(seq: DNASequence) -> BaseFrequency:
    """
    统计各碱基的出现次数和比例。

    Args:
        seq: DNA 序列

    Returns:
        BaseFrequency 对象，包含各碱基的计数和比例

    示例：
        >>> freq = base_frequency(DNASequence("AATGC"))
        >>> freq.a_count
        2
        >>> freq.a_ratio
        0.4
    """
    counts = Counter(seq)

    return BaseFrequency(
        a_count=counts.get("A", 0),
        t_count=counts.get("T", 0),
        g_count=counts.get("G", 0),
        c_count=counts.get("C", 0),
        total=len(seq),
    )


def find_motif(seq: DNASequence, motif: str) -> list[int]:
    """
    在序列中查找 motif（短序列模式）的所有出现位置。

    Motif 是具有生物学意义的短序列模式，例如：
    - 转录因子结合位点
    - 限制性内切酶识别位点
    - 启动子序列

    Args:
        seq: 要搜索的 DNA 序列
        motif: 要查找的短序列模式

    Returns:
        所有匹配位置的列表（0-indexed）

    示例：
        >>> find_motif(DNASequence("ATGATGATG"), "ATG")
        [0, 3, 6]
    """
    # 验证 motif
    motif = motif.upper()
    invalid = set(motif) - DNASequence("A").VALID_BASES
    if invalid:
        raise ValueError(f"Motif 包含无效字符: {invalid}")

    positions: list[int] = []
    seq_str = seq.sequence
    motif_len = len(motif)

    # 滑动窗口查找
    # 为什么不用 str.find()？因为我们要找所有位置，包括重叠的
    for i in range(len(seq_str) - motif_len + 1):
        if seq_str[i : i + motif_len] == motif:
            positions.append(i)

    return positions


def transcribe(seq: DNASequence) -> str:
    """
    将 DNA 序列转录为 RNA 序列。

    转录是基因表达的第一步：DNA -> RNA -> 蛋白质
    规则很简单：将 T（胸腺嘧啶）替换为 U（尿嘧啶）

    Args:
        seq: DNA 序列

    Returns:
        RNA 序列字符串（包含 A、U、G、C）

    示例：
        >>> transcribe(DNASequence("ATGC"))
        'AUGC'
    """
    return seq.sequence.replace("T", "U")


def calculate_molecular_weight(seq: DNASequence) -> float:
    """
    计算单链 DNA 的近似分子量（单位：道尔顿 Da）。

    分子量在实验中很重要：
    - 凝胶电泳迁移率与分子量相关
    - 摩尔浓度计算需要分子量

    使用各碱基的平均分子量（已减去水分子，因为形成磷酸二酯键会脱水）

    Args:
        seq: DNA 序列

    Returns:
        分子量（道尔顿）
    """
    # 各碱基核苷酸的分子量（Da）
    # 这些是形成多核苷酸链时的有效分子量
    weights = {
        "A": 331.2,  # dAMP
        "T": 322.2,  # dTMP
        "G": 347.2,  # dGMP
        "C": 307.2,  # dCMP
    }

    return sum(weights[base] for base in seq)
