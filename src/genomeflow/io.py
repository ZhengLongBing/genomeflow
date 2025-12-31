"""
序列文件的读写功能。

FASTA 格式是最常见的序列存储格式：
- 以 > 开头的行是序列标题（header）
- 后续行是序列内容
- 可以包含多条序列

示例 FASTA 文件：
    >seq1 Homo sapiens gene
    ATGCGATCGATCG
    ATCGATCGATCGA
    >seq2 Another sequence
    GGCCAATTGGCC
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from genomeflow.sequence import DNASequence


@dataclass
class FastaRecord:
    """
    表示 FASTA 文件中的一条记录。

    Attributes:
        id: 序列标识符（> 后的第一个词）
        description: 完整描述行（不含 >）
        sequence: DNA 序列对象
    """

    id: str
    description: str
    sequence: DNASequence


def read_fasta(file_path: str | Path) -> Iterator[FastaRecord]:
    """
    读取 FASTA 文件，逐条返回序列记录。

    为什么用 Iterator（生成器）而不是返回列表？
    - 内存效率：大文件可能包含数百万条序列，一次性加载会耗尽内存
    - 惰性求值：只在需要时才处理每条序列
    - 管道处理：可以边读边处理，不需要等待全部读完

    Args:
        file_path: FASTA 文件路径

    Yields:
        FastaRecord 对象

    Raises:
        FileNotFoundError: 文件不存在
        InvalidSequenceError: 序列包含无效字符

    示例：
        >>> for record in read_fasta("sequences.fasta"):
        ...     print(f"{record.id}: {len(record.sequence)} bp")
    """
    path = Path(file_path)

    current_header: str | None = None
    current_sequence_parts: list[str] = []

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            # 跳过空行
            if not line:
                continue

            # 遇到新的序列标题
            if line.startswith(">"):
                # 先输出之前积累的序列（如果有）
                if current_header is not None:
                    yield _create_record(current_header, current_sequence_parts)

                # 开始新序列
                current_header = line[1:]  # 去掉 >
                current_sequence_parts = []
            else:
                # 序列内容行
                # FASTA 允许序列分多行，所以要累积
                current_sequence_parts.append(line)

    # 别忘了最后一条序列
    if current_header is not None:
        yield _create_record(current_header, current_sequence_parts)


def _create_record(header: str, sequence_parts: list[str]) -> FastaRecord:
    """从解析的数据创建 FastaRecord 对象。"""
    # 序列 ID 是 header 的第一个词
    parts = header.split(maxsplit=1)
    seq_id = parts[0]
    description = header

    # 合并所有序列行
    sequence_str = "".join(sequence_parts)

    return FastaRecord(
        id=seq_id,
        description=description,
        sequence=DNASequence(sequence_str),
    )


def write_fasta(
    records: Iterator[FastaRecord] | list[FastaRecord],
    file_path: str | Path,
    line_width: int = 60,
) -> None:
    """
    将序列记录写入 FASTA 文件。

    Args:
        records: 要写入的序列记录
        file_path: 输出文件路径
        line_width: 每行序列的最大宽度（FASTA 惯例是 60 或 80）

    示例：
        >>> records = [FastaRecord("seq1", "My sequence", DNASequence("ATGC"))]
        >>> write_fasta(records, "output.fasta")
    """
    path = Path(file_path)

    with path.open("w", encoding="utf-8") as f:
        for record in records:
            # 写入标题行
            f.write(f">{record.description}\n")

            # 写入序列，按 line_width 换行
            seq_str = record.sequence.sequence
            for i in range(0, len(seq_str), line_width):
                f.write(seq_str[i : i + line_width] + "\n")
