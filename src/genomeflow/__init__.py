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

from genomeflow.sequence import DNASequence, InvalidSequenceError
from genomeflow.rna import RNASequence, TranslationResult, CODON_TABLE
from genomeflow.protein import (
    ProteinSequence,
    InvalidProteinError,
    ProteinProperties,
    AminoAcidComposition,
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
