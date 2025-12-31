"""
RNASequence 类的单元测试。

测试 RNA 序列的创建、翻译、ORF 查找等功能。
"""

import pytest

from genomeflow.rna import CODON_TABLE, RNASequence
from genomeflow.sequence import InvalidSequenceError


class TestRNASequenceCreation:
    """测试 RNASequence 的创建和验证。"""

    def test_valid_sequence(self):
        """有效序列应该成功创建。"""
        seq = RNASequence("AUGC")
        assert str(seq) == "AUGC"

    def test_lowercase_converted(self):
        """小写字母应该转换为大写。"""
        seq = RNASequence("augc")
        assert str(seq) == "AUGC"

    def test_invalid_t_raises_error(self):
        """T 应该被拒绝（RNA 用 U）。"""
        with pytest.raises(InvalidSequenceError) as excinfo:
            RNASequence("ATGC")
        assert "T" in str(excinfo.value)

    def test_invalid_character(self):
        """其他无效字符应该抛出异常。"""
        with pytest.raises(InvalidSequenceError):
            RNASequence("AUGX")


class TestRNASequenceComplement:
    """测试 RNA 互补链功能。"""

    def test_complement(self):
        """RNA 互补链：A-U, G-C。"""
        seq = RNASequence("AUGC")
        comp = seq.complement()
        assert str(comp) == "UACG"

    def test_reverse_complement(self):
        """反向互补链。"""
        seq = RNASequence("AUGC")
        rev_comp = seq.reverse_complement()
        assert str(rev_comp) == "GCAU"


class TestRNASequenceReverseTranscription:
    """测试反向转录功能。"""

    def test_reverse_transcribe(self):
        """反向转录应该将 U 替换为 T。"""
        rna = RNASequence("AUGC")
        dna = rna.reverse_transcribe()
        assert str(dna) == "ATGC"

    def test_reverse_transcribe_all_u(self):
        """全 U 序列应该变成全 T。"""
        rna = RNASequence("UUUU")
        dna = rna.reverse_transcribe()
        assert str(dna) == "TTTT"


class TestRNATranslation:
    """测试 RNA 翻译功能。"""

    def test_simple_translation(self):
        """简单翻译：AUG (M) + UUU (F) + UAA (stop)。"""
        rna = RNASequence("AUGUUUUAA")
        result = rna.translate()

        assert result.protein == "MF"
        assert result.stop_codon is True
        assert result.remaining_bases == 0

    def test_translation_no_start_codon(self):
        """没有起始密码子应该返回空。"""
        rna = RNASequence("UUUUUUUAA")
        result = rna.translate(start_codon=True)

        assert result.protein == ""
        assert result.stop_codon is False

    def test_translation_ignore_start_codon(self):
        """不要求起始密码子时应该从头开始翻译。"""
        rna = RNASequence("UUUUUAUAA")
        result = rna.translate(start_codon=False)

        # UUU=F, UUA=L, UAA=stop
        assert result.protein == "FL"
        assert result.stop_codon is True

    def test_translation_no_stop_codon(self):
        """没有终止密码子应该翻译到底。"""
        rna = RNASequence("AUGUUUGCC")
        result = rna.translate()

        assert result.protein == "MFA"
        assert result.stop_codon is False
        assert result.remaining_bases == 0

    def test_translation_remaining_bases(self):
        """不完整密码子应该记录为剩余碱基。"""
        rna = RNASequence("AUGUUUAA")  # 8 碱基，最后 2 个不足一个密码子
        result = rna.translate()

        assert result.protein == "MF"
        assert result.remaining_bases == 2

    def test_translation_continue_past_stop(self):
        """stop_at_stop=False 时应该继续翻译。"""
        rna = RNASequence("AUGUAAUUU")  # M, *, F
        result = rna.translate(stop_at_stop=False)

        assert "*" in result.protein
        assert result.stop_codon is True

    def test_all_amino_acids(self):
        """验证密码子表包含所有 20 种氨基酸。"""
        amino_acids = set(CODON_TABLE.values()) - {"*"}
        assert len(amino_acids) == 20


class TestRNATranslateToProtein:
    """测试翻译为蛋白质对象。"""

    def test_translate_to_protein(self):
        """应该返回 ProteinSequence 对象。"""
        rna = RNASequence("AUGUUUUAA")
        protein = rna.translate_to_protein()

        assert str(protein) == "MF"
        # 验证是 ProteinSequence 类型
        from genomeflow.protein import ProteinSequence

        assert isinstance(protein, ProteinSequence)


class TestRNAFindORFs:
    """测试 ORF 查找功能。"""

    def test_find_single_orf(self):
        """应该找到单个 ORF。"""
        rna = RNASequence("AUGAAAUAA")  # MK*
        orfs = rna.find_orfs(min_length=1)

        assert len(orfs) == 1
        start, end, protein = orfs[0]
        assert start == 0
        assert end == 9
        assert protein == "MK"

    def test_find_multiple_orfs(self):
        """应该找到多个 ORF。"""
        # 两个不重叠的 ORF
        rna = RNASequence("AUGAAAUAAAUGCCUUGA")
        orfs = rna.find_orfs(min_length=1)

        assert len(orfs) >= 1

    def test_min_length_filter(self):
        """应该过滤掉太短的 ORF。"""
        rna = RNASequence("AUGAAAUAA")  # MK* (2 氨基酸)
        orfs = rna.find_orfs(min_length=5)

        assert len(orfs) == 0

    def test_no_orf(self):
        """没有起始密码子时应该返回空。"""
        rna = RNASequence("UUUUUUUUUUAA")
        orfs = rna.find_orfs(min_length=1)

        assert len(orfs) == 0

    def test_orf_no_stop(self):
        """没有终止密码子的不算完整 ORF。"""
        rna = RNASequence("AUGAAAGCC")
        orfs = rna.find_orfs(min_length=1)

        assert len(orfs) == 0
