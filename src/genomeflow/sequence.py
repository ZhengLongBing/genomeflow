"""
DNA 序列类。

继承自 BaseSequence，实现 DNA 特定的碱基规则和操作。
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from genomeflow.base import BaseSequence, InvalidSequenceError

if TYPE_CHECKING:
    from genomeflow.rna import RNASequence

# 重新导出，保持向后兼容
__all__ = ["DNASequence", "InvalidSequenceError"]


class DNASequence(BaseSequence):
    """
    表示一条 DNA 序列。

    DNA 序列只能包含四种碱基：A、T、G、C（不区分大小写）。

    示例：
        >>> seq = DNASequence("ATGC")
        >>> len(seq)
        4
        >>> seq.complement()
        DNASequence('TACG')
    """

    @property
    def VALID_BASES(self) -> frozenset[str]:
        return frozenset("ATGC")

    @property
    def COMPLEMENT_MAP(self) -> dict[str, str]:
        return {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
        }

    def transcribe(self) -> "RNASequence":
        """
        将 DNA 转录为 RNA。

        转录规则：T → U，其他碱基不变。

        这模拟了细胞中 DNA 转录为 mRNA 的过程。
        实际上，转录是基于模板链（反义链）进行的，
        这里简化为直接替换。

        Returns:
            对应的 RNA 序列

        示例：
            >>> DNASequence("ATGC").transcribe()
            RNASequence('AUGC')
        """
        # 延迟导入避免循环依赖
        from genomeflow.rna import RNASequence

        rna_seq = self._sequence.replace("T", "U")
        return RNASequence(rna_seq)
