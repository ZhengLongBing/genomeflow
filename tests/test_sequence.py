"""
DNASequence 类的单元测试。

测试 DNA 序列的创建、验证、操作等功能。
"""

import pytest

from genomeflow.sequence import DNASequence, InvalidSequenceError


class TestDNASequenceCreation:
    """测试 DNASequence 的创建和验证。"""

    def test_valid_sequence(self):
        """有效序列应该成功创建。"""
        seq = DNASequence("ATGC")
        assert str(seq) == "ATGC"

    def test_lowercase_converted(self):
        """小写字母应该转换为大写。"""
        seq = DNASequence("atgc")
        assert str(seq) == "ATGC"

    def test_mixed_case(self):
        """混合大小写应该正常处理。"""
        seq = DNASequence("AtGc")
        assert str(seq) == "ATGC"

    def test_empty_sequence(self):
        """空序列应该可以创建。"""
        seq = DNASequence("")
        assert len(seq) == 0

    def test_invalid_character_raises_error(self):
        """包含无效字符应该抛出异常。"""
        with pytest.raises(InvalidSequenceError) as excinfo:
            DNASequence("ATGX")
        assert "X" in str(excinfo.value)

    def test_invalid_multiple_characters(self):
        """多个无效字符都应该被报告。"""
        with pytest.raises(InvalidSequenceError) as excinfo:
            DNASequence("ATGXYZ")
        assert "X" in str(excinfo.value) or "Y" in str(excinfo.value)


class TestDNASequenceLength:
    """测试序列长度相关功能。"""

    def test_length(self):
        """len() 应该返回正确的长度。"""
        assert len(DNASequence("ATGC")) == 4
        assert len(DNASequence("")) == 0
        assert len(DNASequence("A" * 1000)) == 1000


class TestDNASequenceIndexing:
    """测试序列索引和切片。"""

    def test_single_index(self):
        """单个索引应该返回字符。"""
        seq = DNASequence("ATGC")
        assert seq[0] == "A"
        assert seq[1] == "T"
        assert seq[-1] == "C"

    def test_slice(self):
        """切片应该返回字符串。"""
        seq = DNASequence("ATGCGATC")
        assert seq[0:4] == "ATGC"
        # ATGCGATC[::2] = A, G, G, T (位置 0, 2, 4, 6)
        assert seq[::2] == "AGGT"


class TestDNASequenceIteration:
    """测试序列迭代。"""

    def test_iteration(self):
        """应该能够迭代序列。"""
        seq = DNASequence("ATGC")
        bases = list(seq)
        assert bases == ["A", "T", "G", "C"]

    def test_in_operator(self):
        """in 操作符应该正常工作。"""
        seq = DNASequence("ATGC")
        assert "A" in seq
        assert "X" not in seq


class TestDNASequenceComplement:
    """测试互补链功能。"""

    def test_complement(self):
        """互补链应该正确计算。"""
        seq = DNASequence("ATGC")
        comp = seq.complement()
        assert str(comp) == "TACG"

    def test_reverse_complement(self):
        """反向互补链应该正确计算。"""
        seq = DNASequence("ATGC")
        rev_comp = seq.reverse_complement()
        assert str(rev_comp) == "GCAT"

    def test_complement_empty(self):
        """空序列的互补链也是空的。"""
        seq = DNASequence("")
        assert str(seq.complement()) == ""

    def test_double_complement(self):
        """两次互补应该回到原序列。"""
        seq = DNASequence("ATGCGATCGATCG")
        assert str(seq.complement().complement()) == str(seq)


class TestDNASequenceTranscription:
    """测试转录功能。"""

    def test_transcribe(self):
        """转录应该将 T 替换为 U。"""
        seq = DNASequence("ATGC")
        rna = seq.transcribe()
        assert str(rna) == "AUGC"

    def test_transcribe_no_t(self):
        """没有 T 的序列转录后不变。"""
        seq = DNASequence("AGCC")
        rna = seq.transcribe()
        assert str(rna) == "AGCC"

    def test_transcribe_all_t(self):
        """全 T 序列应该变成全 U。"""
        seq = DNASequence("TTTT")
        rna = seq.transcribe()
        assert str(rna) == "UUUU"


class TestDNASequenceEquality:
    """测试序列相等性比较。"""

    def test_equal_sequences(self):
        """相同序列应该相等。"""
        seq1 = DNASequence("ATGC")
        seq2 = DNASequence("ATGC")
        assert seq1 == seq2

    def test_equal_different_case(self):
        """大小写不同但内容相同应该相等。"""
        seq1 = DNASequence("ATGC")
        seq2 = DNASequence("atgc")
        assert seq1 == seq2

    def test_not_equal(self):
        """不同序列应该不相等。"""
        seq1 = DNASequence("ATGC")
        seq2 = DNASequence("GCTA")
        assert seq1 != seq2

    def test_equal_to_string(self):
        """应该能与字符串比较。"""
        seq = DNASequence("ATGC")
        assert seq == "ATGC"
        assert seq == "atgc"


class TestDNASequenceHash:
    """测试序列哈希功能。"""

    def test_hashable(self):
        """序列应该可哈希。"""
        seq = DNASequence("ATGC")
        hash(seq)  # 不应该抛出异常

    def test_same_hash(self):
        """相同序列应该有相同的哈希值。"""
        seq1 = DNASequence("ATGC")
        seq2 = DNASequence("ATGC")
        assert hash(seq1) == hash(seq2)

    def test_usable_in_set(self):
        """序列应该能放入集合。"""
        seq1 = DNASequence("ATGC")
        seq2 = DNASequence("ATGC")
        seq3 = DNASequence("GCTA")
        s = {seq1, seq2, seq3}
        assert len(s) == 2  # seq1 和 seq2 相同


class TestDNASequenceRepr:
    """测试序列的字符串表示。"""

    def test_short_repr(self):
        """短序列应该完整显示。"""
        seq = DNASequence("ATGC")
        repr_str = repr(seq)
        assert "ATGC" in repr_str
        assert "DNASequence" in repr_str

    def test_long_repr_truncated(self):
        """长序列应该被截断。"""
        seq = DNASequence("A" * 100)
        repr_str = repr(seq)
        assert "..." in repr_str
