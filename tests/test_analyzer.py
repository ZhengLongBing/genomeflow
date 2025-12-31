"""
analyzer 模块的单元测试。

测试 GC 含量、碱基频率、motif 查找等分析功能。
"""

import pytest

from genomeflow.analyzer import (
    gc_content,
    base_frequency,
    find_motif,
    transcribe,
    calculate_molecular_weight,
)
from genomeflow.sequence import DNASequence


class TestGCContent:
    """测试 GC 含量计算。"""

    def test_all_gc(self):
        """全 GC 序列应该返回 1.0。"""
        seq = DNASequence("GGCC")
        assert gc_content(seq) == 1.0

    def test_no_gc(self):
        """没有 GC 应该返回 0.0。"""
        seq = DNASequence("AATT")
        assert gc_content(seq) == 0.0

    def test_half_gc(self):
        """一半 GC 应该返回 0.5。"""
        seq = DNASequence("ATGC")
        assert gc_content(seq) == 0.5

    def test_empty_sequence(self):
        """空序列应该返回 0.0。"""
        seq = DNASequence("")
        assert gc_content(seq) == 0.0

    def test_typical_sequence(self):
        """典型序列的 GC 含量计算。"""
        seq = DNASequence("ATGCGATCGATCG")
        # G: 4, C: 3, A: 3, T: 3 => GC = 7/13 ≈ 0.538
        gc = gc_content(seq)
        assert 0.53 < gc < 0.55


class TestBaseFrequency:
    """测试碱基频率统计。"""

    def test_equal_bases(self):
        """等量碱基应该各占 25%。"""
        seq = DNASequence("ATGC")
        freq = base_frequency(seq)

        assert freq.a_count == 1
        assert freq.t_count == 1
        assert freq.g_count == 1
        assert freq.c_count == 1
        assert freq.total == 4

        assert abs(freq.a_ratio - 0.25) < 0.001
        assert abs(freq.t_ratio - 0.25) < 0.001
        assert abs(freq.g_ratio - 0.25) < 0.001
        assert abs(freq.c_ratio - 0.25) < 0.001

    def test_unequal_bases(self):
        """不等量碱基应该正确计算比例。"""
        seq = DNASequence("AAATG")  # A: 3, T: 1, G: 1, C: 0
        freq = base_frequency(seq)

        assert freq.a_count == 3
        assert freq.t_count == 1
        assert freq.g_count == 1
        assert freq.c_count == 0
        assert freq.total == 5

        assert abs(freq.a_ratio - 0.6) < 0.001

    def test_empty_sequence(self):
        """空序列所有比例应该是 0。"""
        seq = DNASequence("")
        freq = base_frequency(seq)

        assert freq.total == 0
        assert freq.a_ratio == 0.0
        assert freq.t_ratio == 0.0


class TestFindMotif:
    """测试 motif 查找。"""

    def test_single_occurrence(self):
        """单个 motif 出现。"""
        seq = DNASequence("ATGATGATG")
        positions = find_motif(seq, "ATG")
        assert positions == [0, 3, 6]

    def test_no_occurrence(self):
        """没有 motif 应该返回空列表。"""
        seq = DNASequence("ATGATGATG")
        positions = find_motif(seq, "CCC")
        assert positions == []

    def test_overlapping_motif(self):
        """重叠的 motif 应该都被找到。"""
        seq = DNASequence("AAAA")
        positions = find_motif(seq, "AA")
        assert positions == [0, 1, 2]

    def test_case_insensitive(self):
        """motif 搜索应该不区分大小写。"""
        seq = DNASequence("ATGATG")
        positions = find_motif(seq, "atg")
        assert positions == [0, 3]

    def test_invalid_motif(self):
        """包含无效字符的 motif 应该抛出异常。"""
        seq = DNASequence("ATGATG")
        with pytest.raises(ValueError) as excinfo:
            find_motif(seq, "ATX")
        assert "无效字符" in str(excinfo.value)


class TestTranscribe:
    """测试转录功能。"""

    def test_basic_transcription(self):
        """基本转录：T -> U。"""
        seq = DNASequence("ATGC")
        rna = transcribe(seq)
        assert rna == "AUGC"

    def test_all_t(self):
        """全 T 序列。"""
        seq = DNASequence("TTTT")
        rna = transcribe(seq)
        assert rna == "UUUU"

    def test_no_t(self):
        """没有 T 的序列。"""
        seq = DNASequence("AGCC")
        rna = transcribe(seq)
        assert rna == "AGCC"

    def test_empty(self):
        """空序列。"""
        seq = DNASequence("")
        rna = transcribe(seq)
        assert rna == ""


class TestMolecularWeight:
    """测试分子量计算。"""

    def test_single_base(self):
        """单碱基分子量。"""
        seq = DNASequence("A")
        mw = calculate_molecular_weight(seq)
        assert mw == 331.2  # dAMP

    def test_atgc(self):
        """ATGC 的分子量。"""
        seq = DNASequence("ATGC")
        mw = calculate_molecular_weight(seq)
        # 331.2 + 322.2 + 347.2 + 307.2 = 1307.8
        expected = 331.2 + 322.2 + 347.2 + 307.2
        assert abs(mw - expected) < 0.1

    def test_longer_sequence(self):
        """较长序列的分子量应该是正数。"""
        seq = DNASequence("ATGCGATCGATCG")
        mw = calculate_molecular_weight(seq)
        assert mw > 0
