"""
GenomeFlow - DNA 序列分析工具包。

一个面向生物信息学初学者的 Python 库，提供 DNA、RNA 和蛋白质序列
的分析功能。

快速开始：
    >>> from genomeflow import DNASequence, gc_content
    >>> seq = DNASequence("ATGCGATCGATCG")
    >>> print(f"GC 含量: {gc_content(seq):.2%}")
    GC 含量: 53.85%

主要功能：
    - DNA/RNA/蛋白质序列表示
    - GC 含量、碱基频率分析
    - 序列转录和翻译
    - FASTA 文件读写
    - 序列可视化
"""

from genomeflow.analyzer import (
    BaseFrequency,
    base_frequency,
    calculate_molecular_weight,
    find_motif,
    gc_content,
    transcribe,
)
from genomeflow.io import FastaRecord, read_fasta, write_fasta
from genomeflow.protein import (
    AminoAcidComposition,
    InvalidProteinError,
    ProteinProperties,
    ProteinSequence,
)
from genomeflow.rna import CODON_TABLE, RNASequence, TranslationResult
from genomeflow.sequence import DNASequence, InvalidSequenceError

__version__ = "0.1.0"

__all__ = [
    # 版本
    "__version__",
    # 序列类
    "DNASequence",
    "RNASequence",
    "ProteinSequence",
    # 异常
    "InvalidSequenceError",
    "InvalidProteinError",
    # 数据类
    "TranslationResult",
    "ProteinProperties",
    "AminoAcidComposition",
    "BaseFrequency",
    "FastaRecord",
    # 分析函数
    "gc_content",
    "base_frequency",
    "find_motif",
    "transcribe",
    "calculate_molecular_weight",
    # IO 函数
    "read_fasta",
    "write_fasta",
    # 常量
    "CODON_TABLE",
]
