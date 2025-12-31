"""
io 模块的单元测试。

测试 FASTA 文件的读写功能。
"""

from pathlib import Path

import pytest

from genomeflow.io import FastaRecord, read_fasta, write_fasta
from genomeflow.sequence import DNASequence


class TestReadFasta:
    """测试 FASTA 文件读取。"""

    def test_read_single_sequence(self, tmp_path: Path):
        """读取单条序列。"""
        fasta_file = tmp_path / "single.fasta"
        fasta_file.write_text(">seq1 Test sequence\nATGCGATCGATCG\n")

        records = list(read_fasta(fasta_file))

        assert len(records) == 1
        assert records[0].id == "seq1"
        assert records[0].description == "seq1 Test sequence"
        assert str(records[0].sequence) == "ATGCGATCGATCG"

    def test_read_multiple_sequences(self, tmp_path: Path):
        """读取多条序列。"""
        fasta_file = tmp_path / "multiple.fasta"
        fasta_file.write_text(
            ">seq1 First\nATGC\n>seq2 Second\nGCTA\n>seq3 Third\nAAAA\n"
        )

        records = list(read_fasta(fasta_file))

        assert len(records) == 3
        assert records[0].id == "seq1"
        assert records[1].id == "seq2"
        assert records[2].id == "seq3"

    def test_read_multiline_sequence(self, tmp_path: Path):
        """读取多行序列。"""
        fasta_file = tmp_path / "multiline.fasta"
        fasta_file.write_text(">seq1 Multi-line\nATGC\nGCTA\nAAAA\n")

        records = list(read_fasta(fasta_file))

        assert len(records) == 1
        assert str(records[0].sequence) == "ATGCGCTAAAAA"

    def test_skip_empty_lines(self, tmp_path: Path):
        """跳过空行。"""
        fasta_file = tmp_path / "empty_lines.fasta"
        fasta_file.write_text(">seq1 Test\n\nATGC\n\nGCTA\n\n")

        records = list(read_fasta(fasta_file))

        assert len(records) == 1
        assert str(records[0].sequence) == "ATGCGCTA"

    def test_read_empty_file(self, tmp_path: Path):
        """读取空文件。"""
        fasta_file = tmp_path / "empty.fasta"
        fasta_file.write_text("")

        records = list(read_fasta(fasta_file))

        assert len(records) == 0

    def test_file_not_found(self, tmp_path: Path):
        """文件不存在应该抛出异常。"""
        with pytest.raises(FileNotFoundError):
            list(read_fasta(tmp_path / "nonexistent.fasta"))

    def test_id_extraction(self, tmp_path: Path):
        """ID 应该是 header 的第一个词。"""
        fasta_file = tmp_path / "header.fasta"
        fasta_file.write_text(">seq1|gene|organism Description here\nATGC\n")

        records = list(read_fasta(fasta_file))

        assert records[0].id == "seq1|gene|organism"
        assert records[0].description == "seq1|gene|organism Description here"


class TestWriteFasta:
    """测试 FASTA 文件写入。"""

    def test_write_single_record(self, tmp_path: Path):
        """写入单条记录。"""
        output_file = tmp_path / "output.fasta"
        records = [
            FastaRecord(
                id="seq1",
                description="seq1 Test sequence",
                sequence=DNASequence("ATGC"),
            )
        ]

        write_fasta(records, output_file)

        content = output_file.read_text()
        assert ">seq1 Test sequence" in content
        assert "ATGC" in content

    def test_write_multiple_records(self, tmp_path: Path):
        """写入多条记录。"""
        output_file = tmp_path / "output.fasta"
        records = [
            FastaRecord("seq1", "seq1 First", DNASequence("ATGC")),
            FastaRecord("seq2", "seq2 Second", DNASequence("GCTA")),
        ]

        write_fasta(records, output_file)

        content = output_file.read_text()
        assert ">seq1 First" in content
        assert ">seq2 Second" in content
        assert "ATGC" in content
        assert "GCTA" in content

    def test_line_wrapping(self, tmp_path: Path):
        """长序列应该按指定宽度换行。"""
        output_file = tmp_path / "output.fasta"
        long_seq = "A" * 100
        records = [
            FastaRecord("seq1", "seq1 Long sequence", DNASequence(long_seq))
        ]

        write_fasta(records, output_file, line_width=60)

        lines = output_file.read_text().strip().split("\n")
        # 第一行是 header，后面是序列
        assert lines[0].startswith(">")
        assert len(lines[1]) == 60  # 第一行序列 60 字符
        assert len(lines[2]) == 40  # 第二行序列 40 字符

    def test_custom_line_width(self, tmp_path: Path):
        """自定义行宽。"""
        output_file = tmp_path / "output.fasta"
        records = [
            FastaRecord("seq1", "seq1 Test", DNASequence("A" * 100))
        ]

        write_fasta(records, output_file, line_width=20)

        lines = output_file.read_text().strip().split("\n")
        # 应该有 5 行序列（100 / 20 = 5）
        seq_lines = [line for line in lines if not line.startswith(">")]
        assert len(seq_lines) == 5
        assert all(len(line) == 20 for line in seq_lines)


class TestFastaRoundTrip:
    """测试 FASTA 文件的读写往返。"""

    def test_roundtrip(self, tmp_path: Path):
        """写入后读取应该得到相同内容。"""
        output_file = tmp_path / "roundtrip.fasta"

        original_records = [
            FastaRecord("seq1", "seq1 First sequence", DNASequence("ATGCGATCGATCG")),
            FastaRecord("seq2", "seq2 Second sequence", DNASequence("GCTAGCTAGCTA")),
        ]

        # 写入
        write_fasta(original_records, output_file)

        # 读取
        read_records = list(read_fasta(output_file))

        # 验证
        assert len(read_records) == 2
        for orig, read in zip(original_records, read_records, strict=False):
            assert orig.id == read.id
            assert orig.description == read.description
            assert str(orig.sequence) == str(read.sequence)
