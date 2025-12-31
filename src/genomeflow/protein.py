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
    "A": 89.09,  # Alanine
    "C": 121.15,  # Cysteine
    "D": 133.10,  # Aspartic acid
    "E": 147.13,  # Glutamic acid
    "F": 165.19,  # Phenylalanine
    "G": 75.07,  # Glycine
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
    "A": 1.8,  # 疏水
    "C": 2.5,
    "D": -3.5,  # 亲水（酸性）
    "E": -3.5,
    "F": 2.8,  # 强疏水（芳香族）
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,  # 强疏水
    "K": -3.9,  # 亲水（碱性）
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,  # 最亲水
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,  # 强疏水
    "W": -0.9,
    "Y": -1.3,
}

# 氨基酸 pKa 值（用于计算等电点）
# pKa_COOH: C端羧基
# pKa_NH2: N端氨基
# pKa_side: 侧链（如果可解离）
PKA_VALUES: dict[str, dict[str, float]] = {
    "D": {"side": 3.9},  # 天冬氨酸侧链羧基
    "E": {"side": 4.1},  # 谷氨酸侧链羧基
    "H": {"side": 6.0},  # 组氨酸侧链咪唑
    "C": {"side": 8.3},  # 半胱氨酸侧链巯基
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
                f"序列包含无效字符: {invalid_chars}。" f"有效氨基酸: {sorted(valid_set)}"
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
            aa: count / total if total > 0 else 0.0 for aa, count in counts.items()
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
            window_seq = self._sequence[i - half_window : i + half_window + 1]
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
