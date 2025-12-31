# 教程 03：蛋白质序列支持

## 📋 本章导览

- **你将掌握的技能**：
  - 理解蛋白质序列的特点
  - 实现 ProteinSequence 类
  - 计算蛋白质理化性质（分子量、等电点等）
  - 分析氨基酸组成和疏水性
  - 构建完整的"中心法则"工作流

- **前置知识**：
  - 完成教程 01 和 02
  - 理解氨基酸的基本概念

- **核心挑战**：蛋白质有 20 种氨基酸，如何高效地处理和分析？

---

## 📚 理解蛋白质序列

### 什么是蛋白质？

蛋白质是生命活动的执行者，由氨基酸链组成：

```
DNA → RNA → 蛋白质
         ↑
      我们在这里
```

### 20 种标准氨基酸

| 符号 | 三字母 | 名称 | 性质 |
|------|--------|------|------|
| A | Ala | 丙氨酸 | 疏水，小 |
| C | Cys | 半胱氨酸 | 可形成二硫键 |
| D | Asp | 天冬氨酸 | 酸性，带负电 |
| E | Glu | 谷氨酸 | 酸性，带负电 |
| F | Phe | 苯丙氨酸 | 疏水，芳香族 |
| G | Gly | 甘氨酸 | 最小，灵活 |
| H | His | 组氨酸 | 碱性，可带正电 |
| I | Ile | 异亮氨酸 | 疏水 |
| K | Lys | 赖氨酸 | 碱性，带正电 |
| L | Leu | 亮氨酸 | 疏水 |
| M | Met | 甲硫氨酸 | 起始氨基酸 |
| N | Asn | 天冬酰胺 | 亲水 |
| P | Pro | 脯氨酸 | 结构限制 |
| Q | Gln | 谷氨酰胺 | 亲水 |
| R | Arg | 精氨酸 | 碱性，带正电 |
| S | Ser | 丝氨酸 | 亲水，可磷酸化 |
| T | Thr | 苏氨酸 | 亲水，可磷酸化 |
| V | Val | 缬氨酸 | 疏水 |
| W | Trp | 色氨酸 | 疏水，芳香族 |
| Y | Tyr | 酪氨酸 | 芳香族，可磷酸化 |

### 蛋白质分析的意义

**为什么要分析蛋白质序列？**

1. **功能预测**：序列 → 结构 → 功能
2. **蛋白质工程**：改造蛋白质性质
3. **药物设计**：靶点分析
4. **进化研究**：序列相似性分析

**常见分析指标**：

| 指标 | 说明 | 用途 |
|------|------|------|
| 分子量 | 蛋白质质量 | 电泳、质谱 |
| 等电点 (pI) | 净电荷为 0 的 pH | 纯化、电泳 |
| 疏水性 | 亲/疏水程度 | 膜蛋白预测 |
| 消光系数 | 吸光能力 | 浓度测定 |

---

## 🔧 设计我们的实现

### 与核酸序列的区别

| 特性 | DNA/RNA | 蛋白质 |
|------|---------|--------|
| 单元数量 | 4 种碱基 | 20 种氨基酸 |
| 互补概念 | 有（配对） | 无 |
| 方向 | 5' → 3' | N端 → C端 |
| 化学性质 | 较统一 | 差异大 |

**设计决策**：蛋白质序列不继承 `BaseSequence`，因为没有"互补"概念。

### 架构图

```
┌─────────────────────────────────────────────────────────┐
│              ProteinSequence                            │
│  ├── _sequence: str (氨基酸序列)                         │
│  ├── molecular_weight() → float                        │
│  ├── isoelectric_point() → float                       │
│  ├── amino_acid_composition() → dict                   │
│  ├── hydrophobicity_profile() → list[float]            │
│  └── extinction_coefficient() → float                  │
└─────────────────────────────────────────────────────────┘
                          ↑
                          │ translate()
┌─────────────────────────────────────────────────────────┐
│                    RNASequence                          │
└─────────────────────────────────────────────────────────┘
```

---

## 💻 代码实现

### 步骤 1：创建 ProteinSequence 类

创建 `src/genomeflow/protein.py`：

```python
"""
蛋白质序列类。

蛋白质是由氨基酸组成的生物大分子，
是基因表达的最终产物，执行各种生物学功能。
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterator


class InvalidProteinError(ValueError):
    """蛋白质序列包含无效字符时抛出。"""
    pass


# 20 种标准氨基酸
AMINO_ACIDS: frozenset[str] = frozenset("ACDEFGHIKLMNPQRSTVWY")

# 特殊符号
# * = 终止密码子产生的终止符
# X = 未知氨基酸
# - = 比对中的间隙
EXTENDED_AMINO_ACIDS: frozenset[str] = frozenset("ACDEFGHIKLMNPQRSTVWY*X-")


# 氨基酸分子量（Da）
# 这是残基分子量（已减去水，因为肽键形成会脱水）
AMINO_ACID_WEIGHTS: dict[str, float] = {
    "A": 89.09,   # Alanine
    "C": 121.15,  # Cysteine
    "D": 133.10,  # Aspartic acid
    "E": 147.13,  # Glutamic acid
    "F": 165.19,  # Phenylalanine
    "G": 75.07,   # Glycine
    "H": 155.16,  # Histidine
    "I": 131.17,  # Isoleucine
    "K": 146.19,  # Lysine
    "L": 131.17,  # Leucine
    "M": 149.21,  # Methionine
    "N": 132.12,  # Asparagine
    "P": 115.13,  # Proline
    "Q": 146.15,  # Glutamine
    "R": 174.20,  # Arginine
    "S": 105.09,  # Serine
    "T": 119.12,  # Threonine
    "V": 117.15,  # Valine
    "W": 204.23,  # Tryptophan
    "Y": 181.19,  # Tyrosine
}

# 疏水性指数（Kyte-Doolittle 量表）
# 正值表示疏水，负值表示亲水
HYDROPHOBICITY: dict[str, float] = {
    "A": 1.8,    # 疏水
    "C": 2.5,
    "D": -3.5,   # 亲水（酸性）
    "E": -3.5,
    "F": 2.8,    # 强疏水（芳香族）
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,    # 强疏水
    "K": -3.9,   # 亲水（碱性）
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,   # 最亲水
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,    # 强疏水
    "W": -0.9,
    "Y": -1.3,
}

# 氨基酸 pKa 值（用于计算等电点）
# pKa_COOH: C端羧基
# pKa_NH2: N端氨基
# pKa_side: 侧链（如果可解离）
PKA_VALUES: dict[str, dict[str, float]] = {
    "D": {"side": 3.9},   # 天冬氨酸侧链羧基
    "E": {"side": 4.1},   # 谷氨酸侧链羧基
    "H": {"side": 6.0},   # 组氨酸侧链咪唑
    "C": {"side": 8.3},   # 半胱氨酸侧链巯基
    "Y": {"side": 10.1},  # 酪氨酸侧链酚羟基
    "K": {"side": 10.5},  # 赖氨酸侧链氨基
    "R": {"side": 12.5},  # 精氨酸侧链胍基
}

# N端和C端的 pKa
PKA_N_TERMINUS = 9.69
PKA_C_TERMINUS = 2.34


@dataclass(frozen=True)
class AminoAcidComposition:
    """氨基酸组成分析结果。"""
    counts: dict[str, int]
    frequencies: dict[str, float]
    total: int

    def __str__(self) -> str:
        lines = [f"总氨基酸数: {self.total}"]
        for aa in sorted(self.counts.keys()):
            count = self.counts[aa]
            freq = self.frequencies[aa]
            lines.append(f"  {aa}: {count:4d} ({freq:5.1%})")
        return "\n".join(lines)


@dataclass(frozen=True)
class ProteinProperties:
    """蛋白质理化性质。"""
    molecular_weight: float
    isoelectric_point: float
    charge_at_ph7: float
    gravy: float  # Grand Average of Hydropathy
    extinction_coefficient: tuple[float, float]  # (还原态, 氧化态)


class ProteinSequence:
    """
    表示一条蛋白质（氨基酸）序列。

    蛋白质序列使用单字母氨基酸代码表示。
    支持 20 种标准氨基酸，以及一些特殊符号。

    示例：
        >>> protein = ProteinSequence("MKFLILLFNILCLFPVLAADNH")
        >>> protein.molecular_weight()
        2567.89
    """

    # 允许的字符：标准氨基酸 + 特殊符号
    VALID_CHARS: frozenset[str] = EXTENDED_AMINO_ACIDS

    def __init__(self, sequence: str, strict: bool = True) -> None:
        """
        创建蛋白质序列对象。

        Args:
            sequence: 氨基酸序列（单字母代码）
            strict: 严格模式。如果为 True，只接受 20 种标准氨基酸；
                   如果为 False，接受 X（未知）和 *（终止）等符号。

        Raises:
            InvalidProteinError: 如果序列包含无效字符
        """
        normalized = sequence.upper()

        # 验证
        valid_set = AMINO_ACIDS if strict else self.VALID_CHARS
        invalid_chars = set(normalized) - valid_set
        if invalid_chars:
            raise InvalidProteinError(
                f"序列包含无效字符: {invalid_chars}。"
                f"有效氨基酸: {sorted(valid_set)}"
            )

        self._sequence = normalized
        self._strict = strict

    @property
    def sequence(self) -> str:
        """返回序列字符串。"""
        return self._sequence

    def __len__(self) -> int:
        return len(self._sequence)

    def __getitem__(self, index: int | slice) -> str:
        return self._sequence[index]

    def __iter__(self) -> Iterator[str]:
        return iter(self._sequence)

    def __str__(self) -> str:
        return self._sequence

    def __repr__(self) -> str:
        if len(self._sequence) > 50:
            display = f"{self._sequence[:25]}...{self._sequence[-25:]}"
        else:
            display = self._sequence
        return f"ProteinSequence('{display}')"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, ProteinSequence):
            return self._sequence == other._sequence
        if isinstance(other, str):
            return self._sequence == other.upper()
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self._sequence)

    def molecular_weight(self) -> float:
        """
        计算蛋白质分子量（Da）。

        分子量 = Σ(各氨基酸残基分子量) + 水分子分子量
        （N端和C端各保留一个H/OH，形成完整多肽）

        Returns:
            分子量（道尔顿）

        示例：
            >>> ProteinSequence("MK").molecular_weight()
            277.40  # 149.21 + 146.19 - 18.0 (肽键脱水)
        """
        # 计算残基分子量之和
        weight = sum(
            AMINO_ACID_WEIGHTS.get(aa, 0.0)
            for aa in self._sequence
            if aa in AMINO_ACID_WEIGHTS
        )

        # 减去形成肽键脱去的水（n-1 个肽键）
        num_peptide_bonds = max(0, len(self._sequence) - 1)
        water_weight = 18.015  # H2O 分子量
        weight -= num_peptide_bonds * water_weight

        # 加上一个水（N端的H和C端的OH）
        weight += water_weight

        return round(weight, 2)

    def amino_acid_composition(self) -> AminoAcidComposition:
        """
        分析氨基酸组成。

        Returns:
            AminoAcidComposition 对象，包含计数和频率

        示例：
            >>> protein = ProteinSequence("AAAMKKK")
            >>> comp = protein.amino_acid_composition()
            >>> comp.counts["A"]
            3
        """
        # 只统计标准氨基酸
        standard_seq = [aa for aa in self._sequence if aa in AMINO_ACIDS]
        counts = Counter(standard_seq)
        total = len(standard_seq)

        frequencies = {
            aa: count / total if total > 0 else 0.0
            for aa, count in counts.items()
        }

        return AminoAcidComposition(
            counts=dict(counts),
            frequencies=frequencies,
            total=total,
        )

    def hydrophobicity_profile(self, window: int = 9) -> list[float]:
        """
        计算疏水性分布（滑动窗口平均）。

        使用 Kyte-Doolittle 量表计算每个位置的疏水性。
        正值表示疏水区域，可能是跨膜区或蛋白质内部。

        Args:
            window: 滑动窗口大小，通常用 7-11

        Returns:
            每个位置的疏水性值列表

        示例：
            >>> protein = ProteinSequence("MILVAGFYW")
            >>> profile = protein.hydrophobicity_profile(window=3)
        """
        if len(self._sequence) < window:
            return []

        profile: list[float] = []
        half_window = window // 2

        for i in range(half_window, len(self._sequence) - half_window):
            window_seq = self._sequence[i - half_window:i + half_window + 1]
            values = [
                HYDROPHOBICITY.get(aa, 0.0)
                for aa in window_seq
                if aa in HYDROPHOBICITY
            ]
            if values:
                avg = sum(values) / len(values)
                profile.append(round(avg, 2))

        return profile

    def gravy(self) -> float:
        """
        计算 GRAVY（Grand Average of Hydropathy）。

        GRAVY 是所有氨基酸疏水性的平均值。
        正值表示蛋白质整体疏水（可能是膜蛋白）。
        负值表示蛋白质整体亲水（可能是可溶蛋白）。

        Returns:
            GRAVY 值

        示例：
            >>> ProteinSequence("ILVAG").gravy()  # 疏水残基
            3.16
            >>> ProteinSequence("DEKER").gravy()  # 亲水残基
            -3.5
        """
        values = [
            HYDROPHOBICITY.get(aa, 0.0)
            for aa in self._sequence
            if aa in HYDROPHOBICITY
        ]
        if not values:
            return 0.0
        return round(sum(values) / len(values), 2)

    def _calculate_charge(self, ph: float) -> float:
        """
        计算在指定 pH 下的净电荷。

        使用 Henderson-Hasselbalch 方程计算每个可解离基团的电荷。
        """
        charge = 0.0

        # N端氨基（正电荷）
        charge += 1.0 / (1.0 + 10 ** (ph - PKA_N_TERMINUS))

        # C端羧基（负电荷）
        charge -= 1.0 / (1.0 + 10 ** (PKA_C_TERMINUS - ph))

        # 侧链
        for aa in self._sequence:
            if aa not in PKA_VALUES:
                continue

            pka = PKA_VALUES[aa]["side"]

            if aa in ("D", "E", "C", "Y"):
                # 酸性侧链（失去H带负电）
                charge -= 1.0 / (1.0 + 10 ** (pka - ph))
            else:
                # 碱性侧链（获得H带正电）
                charge += 1.0 / (1.0 + 10 ** (ph - pka))

        return charge

    def isoelectric_point(self) -> float:
        """
        计算等电点（pI）。

        等电点是蛋白质净电荷为零时的 pH 值。
        在等电点时，蛋白质在电场中不移动。

        使用二分法搜索净电荷为 0 的 pH 值。

        Returns:
            等电点 pH 值

        示例：
            >>> ProteinSequence("MKFLK").isoelectric_point()
            10.0  # 富含碱性氨基酸
        """
        # 二分查找
        ph_low, ph_high = 0.0, 14.0

        while ph_high - ph_low > 0.01:
            ph_mid = (ph_low + ph_high) / 2
            charge = self._calculate_charge(ph_mid)

            if charge > 0:
                ph_low = ph_mid
            else:
                ph_high = ph_mid

        return round((ph_low + ph_high) / 2, 2)

    def charge_at_ph(self, ph: float = 7.0) -> float:
        """
        计算在指定 pH 下的净电荷。

        Args:
            ph: pH 值，默认 7.0（生理 pH）

        Returns:
            净电荷

        示例:
            >>> ProteinSequence("DDDKKK").charge_at_ph(7.0)
            约等于 0（3个负电荷 + 3个正电荷）
        """
        return round(self._calculate_charge(ph), 2)

    def extinction_coefficient(self) -> tuple[float, float]:
        """
        计算消光系数（280nm）。

        消光系数用于通过紫外吸收测定蛋白质浓度。
        返回两个值：
        - 还原态（二硫键被还原）
        - 氧化态（二硫键形成）

        计算公式（Pace et al., 1995）：
        ε = nTyr × 1490 + nTrp × 5500 + nCys × 125（氧化态）

        Returns:
            (还原态消光系数, 氧化态消光系数)，单位 M⁻¹ cm⁻¹

        示例：
            >>> ProteinSequence("WWYY").extinction_coefficient()
            (13960, 13960)  # 无 Cys，两个值相同
        """
        n_trp = self._sequence.count("W")
        n_tyr = self._sequence.count("Y")
        n_cys = self._sequence.count("C")

        # 还原态：Cys 不贡献
        reduced = n_trp * 5500 + n_tyr * 1490

        # 氧化态：假设所有 Cys 形成二硫键
        # 每对 Cys 贡献 125
        oxidized = reduced + (n_cys // 2) * 125

        return (reduced, oxidized)

    def get_properties(self) -> ProteinProperties:
        """
        获取蛋白质的主要理化性质。

        Returns:
            ProteinProperties 对象

        示例：
            >>> props = ProteinSequence("MKFLILLFNILCLFPVLAADNH").get_properties()
            >>> print(f"MW: {props.molecular_weight:.1f} Da")
            >>> print(f"pI: {props.isoelectric_point:.2f}")
        """
        return ProteinProperties(
            molecular_weight=self.molecular_weight(),
            isoelectric_point=self.isoelectric_point(),
            charge_at_ph7=self.charge_at_ph(7.0),
            gravy=self.gravy(),
            extinction_coefficient=self.extinction_coefficient(),
        )
```

### 步骤 2：完善 RNA 到蛋白质的工作流

更新 `src/genomeflow/rna.py`，修改 `translate` 方法返回 `ProteinSequence`：

```python
# 在 rna.py 中添加新的翻译方法

def translate_to_protein(
    self,
    start_codon: bool = True,
    stop_at_stop: bool = True,
) -> "ProteinSequence":
    """
    将 RNA 翻译为蛋白质序列对象。

    Args:
        start_codon: 是否从起始密码子开始
        stop_at_stop: 是否在终止密码子处停止

    Returns:
        ProteinSequence 对象

    示例：
        >>> rna = RNASequence("AUGUUUUAA")
        >>> protein = rna.translate_to_protein()
        >>> protein.molecular_weight()
    """
    from genomeflow.protein import ProteinSequence

    result = self.translate(start_codon, stop_at_stop)
    # 移除终止符号（如果有）
    protein_seq = result.protein.replace("*", "")
    return ProteinSequence(protein_seq, strict=True)
```

### 步骤 3：更新包导出

更新 `src/genomeflow/__init__.py`：

```python
"""
GenomeFlow - 序列分析工具

支持 DNA、RNA 和蛋白质序列的分析。
"""

from genomeflow.base import BaseSequence, InvalidSequenceError
from genomeflow.sequence import DNASequence
from genomeflow.rna import RNASequence, TranslationResult, CODON_TABLE
from genomeflow.protein import (
    ProteinSequence,
    InvalidProteinError,
    AminoAcidComposition,
    ProteinProperties,
    AMINO_ACIDS,
    AMINO_ACID_WEIGHTS,
    HYDROPHOBICITY,
)
from genomeflow.analyzer import (
    gc_content,
    base_frequency,
    find_motif,
    transcribe,
    calculate_molecular_weight,
    BaseFrequency,
)
from genomeflow.io import read_fasta, write_fasta, FastaRecord

__version__ = "0.1.0"

__all__ = [
    # 基类
    "BaseSequence",
    "InvalidSequenceError",
    # DNA
    "DNASequence",
    # RNA
    "RNASequence",
    "TranslationResult",
    "CODON_TABLE",
    # 蛋白质
    "ProteinSequence",
    "InvalidProteinError",
    "AminoAcidComposition",
    "ProteinProperties",
    "AMINO_ACIDS",
    "AMINO_ACID_WEIGHTS",
    "HYDROPHOBICITY",
    # 分析
    "gc_content",
    "base_frequency",
    "find_motif",
    "transcribe",
    "calculate_molecular_weight",
    "BaseFrequency",
    # IO
    "read_fasta",
    "write_fasta",
    "FastaRecord",
]
```

### 步骤 4：编写测试

创建 `tests/test_protein.py`：

```python
"""蛋白质序列的单元测试。"""

import pytest

from genomeflow.protein import (
    ProteinSequence,
    InvalidProteinError,
)


class TestProteinSequenceCreation:
    """测试 ProteinSequence 的创建。"""

    def test_create_valid_sequence(self):
        protein = ProteinSequence("MKFLILLFNILCLFPVLAADNH")
        assert len(protein) == 22

    def test_lowercase_converted(self):
        protein = ProteinSequence("mkfl")
        assert str(protein) == "MKFL"

    def test_invalid_character_rejected(self):
        with pytest.raises(InvalidProteinError):
            ProteinSequence("MKXFL")  # X 在严格模式下无效

    def test_non_strict_mode(self):
        protein = ProteinSequence("MKX*", strict=False)
        assert "X" in str(protein)
        assert "*" in str(protein)

    def test_invalid_base_rejected(self):
        with pytest.raises(InvalidProteinError):
            ProteinSequence("ATGC")  # 这是 DNA，不是蛋白质


class TestMolecularWeight:
    """测试分子量计算。"""

    def test_single_amino_acid(self):
        # Glycine (G) 分子量 75.07
        protein = ProteinSequence("G")
        mw = protein.molecular_weight()
        assert mw == pytest.approx(75.07, rel=0.01)

    def test_dipeptide(self):
        # M + K - H2O（肽键脱水）+ H2O（末端）
        protein = ProteinSequence("MK")
        mw = protein.molecular_weight()
        # 149.21 + 146.19 - 18.015 + 18.015 = 295.40
        expected = 149.21 + 146.19
        assert mw == pytest.approx(expected, rel=0.01)


class TestAminoAcidComposition:
    """测试氨基酸组成分析。"""

    def test_composition(self):
        protein = ProteinSequence("AAAKKK")
        comp = protein.amino_acid_composition()
        assert comp.counts["A"] == 3
        assert comp.counts["K"] == 3
        assert comp.frequencies["A"] == pytest.approx(0.5)

    def test_total_count(self):
        protein = ProteinSequence("MKFLILLFNILCLFPVLAADNH")
        comp = protein.amino_acid_composition()
        assert comp.total == 22


class TestHydrophobicity:
    """测试疏水性分析。"""

    def test_hydrophobic_protein(self):
        # I, L, V 是强疏水氨基酸
        protein = ProteinSequence("ILVILVILVILVILVIL")
        gravy = protein.gravy()
        assert gravy > 3.0  # 应该是正值（疏水）

    def test_hydrophilic_protein(self):
        # D, E, K, R 是亲水氨基酸
        protein = ProteinSequence("DEKRDEKRDEKR")
        gravy = protein.gravy()
        assert gravy < -3.0  # 应该是负值（亲水）

    def test_hydrophobicity_profile(self):
        protein = ProteinSequence("MILVAGFYWDEKR")
        profile = protein.hydrophobicity_profile(window=5)
        assert len(profile) > 0
        # 开头疏水，结尾亲水
        assert profile[0] > profile[-1]


class TestIsoelectricPoint:
    """测试等电点计算。"""

    def test_acidic_protein(self):
        # 富含酸性氨基酸
        protein = ProteinSequence("DDDEEEMMM")
        pi = protein.isoelectric_point()
        assert pi < 5.0  # 酸性蛋白质 pI 低

    def test_basic_protein(self):
        # 富含碱性氨基酸
        protein = ProteinSequence("KKKRRRMMMM")
        pi = protein.isoelectric_point()
        assert pi > 9.0  # 碱性蛋白质 pI 高


class TestExtinctionCoefficient:
    """测试消光系数计算。"""

    def test_no_aromatic(self):
        # 没有 W, Y, C
        protein = ProteinSequence("AAAAKKKKK")
        reduced, oxidized = protein.extinction_coefficient()
        assert reduced == 0
        assert oxidized == 0

    def test_with_tryptophan(self):
        protein = ProteinSequence("WWAAA")
        reduced, oxidized = protein.extinction_coefficient()
        assert reduced == 2 * 5500  # 2 个 Trp

    def test_with_cysteine(self):
        protein = ProteinSequence("CCAAA")
        reduced, oxidized = protein.extinction_coefficient()
        assert oxidized > reduced  # 氧化态考虑二硫键


class TestGetProperties:
    """测试综合属性获取。"""

    def test_get_properties(self):
        protein = ProteinSequence("MKFLILLFNILCLFPVLAADNH")
        props = protein.get_properties()

        assert props.molecular_weight > 0
        assert 0 < props.isoelectric_point < 14
        assert isinstance(props.gravy, float)
        assert len(props.extinction_coefficient) == 2
```

---

## ✅ 完整工作流示例

现在我们可以实现完整的"中心法则"工作流：

```python
from genomeflow import DNASequence

# 1. DNA 序列
dna = DNASequence("ATGAAATTTGGGCCCTAGATG")
print(f"DNA: {dna}")
print(f"长度: {len(dna)} bp")
print(f"GC 含量: {dna.gc_content():.1%}")

# 2. 转录为 RNA
rna = dna.transcribe()
print(f"\nRNA: {rna}")

# 3. 翻译为蛋白质
result = rna.translate()
print(f"\n翻译结果:")
print(f"  蛋白质: {result.protein}")
print(f"  遇到终止密码子: {result.stop_codon}")

# 4. 分析蛋白质性质
if result.protein:
    from genomeflow import ProteinSequence
    protein = ProteinSequence(result.protein.replace("*", ""))

    props = protein.get_properties()
    print(f"\n蛋白质性质:")
    print(f"  分子量: {props.molecular_weight:.1f} Da")
    print(f"  等电点: {props.isoelectric_point:.2f}")
    print(f"  pH 7.0 净电荷: {props.charge_at_ph7:+.1f}")
    print(f"  GRAVY: {props.gravy:.2f}")

    comp = protein.amino_acid_composition()
    print(f"\n氨基酸组成:")
    print(comp)
```

---

## 🤔 深入思考

<details>
<summary>为什么 ProteinSequence 不继承 BaseSequence？</summary>

虽然蛋白质序列和核酸序列有相似之处（都是线性序列），但：

1. **核心概念不同**：
   - 核酸有互补配对，蛋白质没有
   - 核酸有方向（5'→3'），蛋白质也有（N→C），但意义不同

2. **避免不恰当的继承**：
   - 如果继承 BaseSequence，会继承 `complement()` 方法
   - 蛋白质没有"互补链"概念，这会造成混淆

3. **设计原则**：
   - 继承应该表示"is-a"关系
   - 蛋白质不是"一种核酸"

如果需要统一接口，可以使用 Protocol：

```python
from typing import Protocol

class SequenceProtocol(Protocol):
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> str: ...
```

</details>

<details>
<summary>等电点计算的原理是什么？</summary>

等电点（pI）是蛋白质净电荷为零时的 pH。

**可解离基团**：
- N端氨基（pKa ≈ 9.7）
- C端羧基（pKa ≈ 2.3）
- 侧链：D, E（酸性），H, K, R（碱性），C, Y（特殊）

**Henderson-Hasselbalch 方程**：
```
pH = pKa + log([A⁻]/[HA])
```

对于酸性基团（失去 H 带负电）：
```
电荷 = -1 / (1 + 10^(pKa - pH))
```

对于碱性基团（获得 H 带正电）：
```
电荷 = +1 / (1 + 10^(pH - pKa))
```

**二分查找**：在 pH 0-14 范围内搜索净电荷为 0 的点。

</details>

---

## 📝 总结

通过这个教程，你学会了：

1. **蛋白质基础**：20 种氨基酸的性质和分类
2. **ProteinSequence 类**：独立于核酸的序列类
3. **理化性质计算**：分子量、等电点、疏水性、消光系数
4. **完整工作流**：DNA → RNA → 蛋白质 → 性质分析

下一步：[教程 04 - 序列可视化](tutorial-04-visualization.md)
