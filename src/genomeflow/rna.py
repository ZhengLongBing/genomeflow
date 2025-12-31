"""
RNA 序列类。

RNA（核糖核酸）是遗传信息传递的中间分子。
与 DNA 的主要区别是用 U（尿嘧啶）代替 T（胸腺嘧啶）。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from genomeflow.base import BaseSequence

if TYPE_CHECKING:
    from genomeflow.protein import ProteinSequence
    from genomeflow.sequence import DNASequence


# 标准遗传密码表
# 密码子 -> 氨基酸单字母代码
CODON_TABLE: dict[str, str] = {
    # 苯丙氨酸 (Phe, F)
    "UUU": "F",
    "UUC": "F",
    # 亮氨酸 (Leu, L)
    "UUA": "L",
    "UUG": "L",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    # 异亮氨酸 (Ile, I)
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    # 甲硫氨酸 (Met, M) - 起始密码子
    "AUG": "M",
    # 缬氨酸 (Val, V)
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    # 丝氨酸 (Ser, S)
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "AGU": "S",
    "AGC": "S",
    # 脯氨酸 (Pro, P)
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    # 苏氨酸 (Thr, T)
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    # 丙氨酸 (Ala, A)
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    # 酪氨酸 (Tyr, Y)
    "UAU": "Y",
    "UAC": "Y",
    # 终止密码子
    "UAA": "*",
    "UAG": "*",
    "UGA": "*",
    # 组氨酸 (His, H)
    "CAU": "H",
    "CAC": "H",
    # 谷氨酰胺 (Gln, Q)
    "CAA": "Q",
    "CAG": "Q",
    # 天冬酰胺 (Asn, N)
    "AAU": "N",
    "AAC": "N",
    # 赖氨酸 (Lys, K)
    "AAA": "K",
    "AAG": "K",
    # 天冬氨酸 (Asp, D)
    "GAU": "D",
    "GAC": "D",
    # 谷氨酸 (Glu, E)
    "GAA": "E",
    "GAG": "E",
    # 半胱氨酸 (Cys, C)
    "UGU": "C",
    "UGC": "C",
    # 色氨酸 (Trp, W)
    "UGG": "W",
    # 精氨酸 (Arg, R)
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    # 甘氨酸 (Gly, G)
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


@dataclass(frozen=True)
class TranslationResult:
    """
    翻译结果。

    Attributes:
        protein: 氨基酸序列（单字母代码）
        stop_codon: 是否遇到终止密码子
        remaining_bases: 未翻译的碱基数（不足3个）
    """

    protein: str
    stop_codon: bool
    remaining_bases: int


class RNASequence(BaseSequence):
    """
    表示一条 RNA 序列。

    RNA 序列只能包含四种碱基：A、U、G、C（不区分大小写）。

    示例：
        >>> seq = RNASequence("AUGC")
        >>> len(seq)
        4
        >>> seq.complement()
        RNASequence('UACG')
    """

    @property
    def VALID_BASES(self) -> frozenset[str]:
        return frozenset("AUGC")

    @property
    def COMPLEMENT_MAP(self) -> dict[str, str]:
        return {
            "A": "U",
            "U": "A",
            "G": "C",
            "C": "G",
        }

    def reverse_transcribe(self) -> DNASequence:
        """
        将 RNA 反向转录为 DNA。

        反向转录是转录的逆过程：U → T

        这在某些病毒（如 HIV）中自然发生，
        也是分子生物学实验中的常见操作。

        Returns:
            对应的 DNA 序列

        示例：
            >>> RNASequence("AUGC").reverse_transcribe()
            DNASequence('ATGC')
        """
        from genomeflow.sequence import DNASequence

        dna_seq = self._sequence.replace("U", "T")
        return DNASequence(dna_seq)

    def translate(
        self,
        start_codon: bool = True,
        stop_at_stop: bool = True,
    ) -> TranslationResult:
        """
        将 RNA 翻译为蛋白质（氨基酸序列）。

        翻译是基因表达的第二步：RNA → 蛋白质。
        每三个碱基（密码子）编码一个氨基酸。

        Args:
            start_codon: 是否从起始密码子（AUG）开始。
                        如果为 True，会跳过 AUG 之前的序列。
            stop_at_stop: 是否在遇到终止密码子时停止。
                         如果为 False，终止密码子翻译为 '*'。

        Returns:
            TranslationResult 包含翻译结果

        示例：
            >>> RNASequence("AUGUUUUAA").translate()
            TranslationResult(protein='MF', stop_codon=True, remaining_bases=0)
        """
        seq = self._sequence

        # 如果需要从起始密码子开始，找到 AUG
        if start_codon:
            start_pos = seq.find("AUG")
            if start_pos == -1:
                return TranslationResult(
                    protein="",
                    stop_codon=False,
                    remaining_bases=len(seq),
                )
            seq = seq[start_pos:]

        # 翻译
        protein_parts: list[str] = []
        stop_found = False

        for i in range(0, len(seq) - 2, 3):
            codon = seq[i : i + 3]
            amino_acid = CODON_TABLE.get(codon, "X")  # X 表示未知

            if amino_acid == "*":
                stop_found = True
                if stop_at_stop:
                    break
                protein_parts.append(amino_acid)
            else:
                protein_parts.append(amino_acid)

        # 计算剩余碱基
        translated_length = len(protein_parts) * 3
        if stop_found and stop_at_stop:
            translated_length += 3  # 终止密码子也消耗3个碱基
        remaining = len(seq) - translated_length

        return TranslationResult(
            protein="".join(protein_parts),
            stop_codon=stop_found,
            remaining_bases=max(0, remaining),
        )

    def translate_to_protein(
        self,
        start_codon: bool = True,
        stop_at_stop: bool = True,
    ) -> ProteinSequence:
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

    def find_orfs(self, min_length: int = 30) -> list[tuple[int, int, str]]:
        """
        查找所有开放阅读框（Open Reading Frame, ORF）。

        ORF 是从起始密码子（AUG）到终止密码子的序列，
        可能编码蛋白质。

        Args:
            min_length: 最小 ORF 长度（氨基酸数），默认 30

        Returns:
            列表，每个元素是 (起始位置, 结束位置, 蛋白质序列)

        示例：
            >>> rna = RNASequence("AUGAAAUAA")
            >>> rna.find_orfs(min_length=1)
            [(0, 9, 'MK')]
        """
        orfs: list[tuple[int, int, str]] = []
        seq = self._sequence

        # 在三个阅读框中搜索
        for frame in range(3):
            i = frame
            while i < len(seq) - 2:
                codon = seq[i : i + 3]

                # 找到起始密码子
                if codon == "AUG":
                    start = i
                    protein_parts: list[str] = []

                    # 继续翻译直到终止密码子
                    j = i
                    while j < len(seq) - 2:
                        codon = seq[j : j + 3]
                        amino_acid = CODON_TABLE.get(codon, "X")

                        if amino_acid == "*":
                            # 找到完整的 ORF
                            if len(protein_parts) >= min_length:
                                orfs.append(
                                    (
                                        start,
                                        j + 3,
                                        "".join(protein_parts),
                                    )
                                )
                            break
                        else:
                            protein_parts.append(amino_acid)
                        j += 3

                    i = j + 3  # 跳过已处理的部分
                else:
                    i += 3

        return orfs
