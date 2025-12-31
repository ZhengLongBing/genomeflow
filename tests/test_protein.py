"""
ProteinSequence 类的单元测试。

测试蛋白质序列的创建、理化性质计算等功能。
"""

import pytest

from genomeflow.protein import (
    AMINO_ACIDS,
    InvalidProteinError,
    ProteinSequence,
)


class TestProteinSequenceCreation:
    """测试 ProteinSequence 的创建和验证。"""

    def test_valid_sequence(self):
        """有效序列应该成功创建。"""
        seq = ProteinSequence("MKFLK")
        assert str(seq) == "MKFLK"

    def test_lowercase_converted(self):
        """小写字母应该转换为大写。"""
        seq = ProteinSequence("mkflk")
        assert str(seq) == "MKFLK"

    def test_all_amino_acids(self):
        """应该支持所有 20 种标准氨基酸。"""
        all_aa = "".join(sorted(AMINO_ACIDS))
        seq = ProteinSequence(all_aa)
        assert len(seq) == 20

    def test_invalid_character_strict(self):
        """严格模式下无效字符应该抛出异常。"""
        with pytest.raises(InvalidProteinError):
            ProteinSequence("MK*", strict=True)

    def test_extended_characters_non_strict(self):
        """非严格模式应该接受 * 和 X。"""
        seq = ProteinSequence("MK*X", strict=False)
        assert str(seq) == "MK*X"

    def test_invalid_character_non_strict(self):
        """非严格模式下无效字符仍应抛出异常。"""
        with pytest.raises(InvalidProteinError):
            ProteinSequence("MKZ", strict=False)  # Z 不是有效字符


class TestProteinSequenceLength:
    """测试序列长度。"""

    def test_length(self):
        """len() 应该返回正确的长度。"""
        assert len(ProteinSequence("MKFLK")) == 5
        assert len(ProteinSequence("M")) == 1


class TestProteinMolecularWeight:
    """测试分子量计算。"""

    def test_single_amino_acid(self):
        """单个氨基酸的分子量。"""
        # 单个氨基酸 = 残基分子量 + 水
        met = ProteinSequence("M")
        mw = met.molecular_weight()
        # M 的残基分子量是 149.21，加上水 18.015
        expected = 149.21 + 18.015 - 0 * 18.015  # 0 个肽键
        assert abs(mw - expected) < 0.1

    def test_dipeptide(self):
        """二肽的分子量。"""
        # 两个氨基酸，一个肽键
        dipeptide = ProteinSequence("MK")
        mw = dipeptide.molecular_weight()
        # 149.21 + 146.19 - 18.015 (一个肽键脱水) + 18.015 (端基水)
        expected = 149.21 + 146.19
        assert abs(mw - expected) < 0.1

    def test_molecular_weight_positive(self):
        """分子量应该是正数。"""
        seq = ProteinSequence("MKFLILLFNILCLFPVLAADNH")
        mw = seq.molecular_weight()
        assert mw > 0


class TestProteinAminoAcidComposition:
    """测试氨基酸组成分析。"""

    def test_composition_counts(self):
        """应该正确统计各氨基酸数量。"""
        seq = ProteinSequence("AAAMKK")
        comp = seq.amino_acid_composition()

        assert comp.counts["A"] == 3
        assert comp.counts["M"] == 1
        assert comp.counts["K"] == 2
        assert comp.total == 6

    def test_composition_frequencies(self):
        """应该正确计算频率。"""
        seq = ProteinSequence("AAMM")
        comp = seq.amino_acid_composition()

        assert abs(comp.frequencies["A"] - 0.5) < 0.01
        assert abs(comp.frequencies["M"] - 0.5) < 0.01


class TestProteinHydrophobicity:
    """测试疏水性相关功能。"""

    def test_gravy_hydrophobic(self):
        """疏水残基应该有正的 GRAVY。"""
        # I, L, V 都是强疏水
        seq = ProteinSequence("ILVAG")
        gravy = seq.gravy()
        assert gravy > 0

    def test_gravy_hydrophilic(self):
        """亲水残基应该有负的 GRAVY。"""
        # D, E, K, R 都是亲水
        seq = ProteinSequence("DEKRN")
        gravy = seq.gravy()
        assert gravy < 0

    def test_hydrophobicity_profile_length(self):
        """疏水性分布应该返回正确长度。"""
        seq = ProteinSequence("MKFLILLFNILCLFPVLAADNH")
        profile = seq.hydrophobicity_profile(window=9)
        expected_length = len(seq) - 9 + 1  # 窗口滑动次数
        assert len(profile) == expected_length

    def test_hydrophobicity_profile_short_sequence(self):
        """序列太短应该返回空列表。"""
        seq = ProteinSequence("MKF")
        profile = seq.hydrophobicity_profile(window=9)
        assert profile == []


class TestProteinIsoelectricPoint:
    """测试等电点计算。"""

    def test_pi_basic_protein(self):
        """富含碱性氨基酸的蛋白质应该有高 pI。"""
        # 富含 K, R 的蛋白质
        seq = ProteinSequence("MKKKRRR")
        pi = seq.isoelectric_point()
        assert pi > 9.0

    def test_pi_acidic_protein(self):
        """富含酸性氨基酸的蛋白质应该有低 pI。"""
        # 富含 D, E 的蛋白质
        seq = ProteinSequence("MDDDEEE")
        pi = seq.isoelectric_point()
        assert pi < 5.0

    def test_pi_neutral_protein(self):
        """中性蛋白质应该有中等 pI。"""
        seq = ProteinSequence("MKDE")  # 1 K, 1 D, 1 E
        pi = seq.isoelectric_point()
        assert 4.0 < pi < 10.0


class TestProteinCharge:
    """测试电荷计算。"""

    def test_charge_at_ph7(self):
        """在 pH 7 时计算电荷。"""
        seq = ProteinSequence("MKDE")
        charge = seq.charge_at_ph(7.0)
        # K 带正电，D 和 E 带负电
        # 结果应该接近 -1
        assert -2 < charge < 0


class TestProteinExtinctionCoefficient:
    """测试消光系数计算。"""

    def test_extinction_with_trp_tyr(self):
        """含有 W 和 Y 应该有非零消光系数。"""
        seq = ProteinSequence("MWWY")
        reduced, oxidized = seq.extinction_coefficient()

        # 2 W * 5500 + 1 Y * 1490 = 12490
        assert reduced == 2 * 5500 + 1 * 1490
        assert oxidized == reduced  # 没有 Cys

    def test_extinction_with_cys(self):
        """含有 Cys 时氧化态和还原态应该不同。"""
        seq = ProteinSequence("MWCC")
        reduced, oxidized = seq.extinction_coefficient()

        assert reduced == 1 * 5500  # 只有 W
        assert oxidized == reduced + 1 * 125  # 一对 Cys 形成二硫键


class TestProteinProperties:
    """测试获取综合理化性质。"""

    def test_get_properties(self):
        """应该返回包含所有性质的对象。"""
        seq = ProteinSequence("MKFLILLFNILCLFPVLAADNH")
        props = seq.get_properties()

        assert props.molecular_weight > 0
        assert 0 < props.isoelectric_point < 14
        assert isinstance(props.gravy, float)
        assert isinstance(props.extinction_coefficient, tuple)
        assert len(props.extinction_coefficient) == 2
