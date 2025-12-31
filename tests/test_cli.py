"""
CLI 模块的单元测试。

使用 Click 的测试工具 CliRunner 测试命令行接口。
"""

from pathlib import Path

import pytest
from click.testing import CliRunner

from genomeflow.cli import cli


@pytest.fixture
def runner():
    """创建 CLI 测试运行器。"""
    return CliRunner()


class TestQuickCommand:
    """测试 quick 命令。"""

    def test_quick_valid_sequence(self, runner: CliRunner):
        """有效序列应该正常分析。"""
        result = runner.invoke(cli, ["quick", "ATGCGATCG"])

        assert result.exit_code == 0
        assert "序列:" in result.output
        assert "长度:" in result.output
        assert "GC 含量:" in result.output
        assert "碱基组成:" in result.output

    def test_quick_shows_transcription(self, runner: CliRunner):
        """应该显示 RNA 转录结果。"""
        result = runner.invoke(cli, ["quick", "ATGC"])

        assert result.exit_code == 0
        assert "RNA 转录:" in result.output
        assert "AUGC" in result.output

    def test_quick_shows_reverse_complement(self, runner: CliRunner):
        """应该显示反向互补链。"""
        result = runner.invoke(cli, ["quick", "ATGC"])

        assert result.exit_code == 0
        assert "反向互补:" in result.output
        assert "GCAT" in result.output

    def test_quick_invalid_sequence(self, runner: CliRunner):
        """无效序列应该报错。"""
        result = runner.invoke(cli, ["quick", "ATGX"])

        assert result.exit_code != 0
        assert "错误" in result.output


class TestComplementCommand:
    """测试 complement 命令。"""

    def test_complement(self, runner: CliRunner):
        """应该返回互补链。"""
        result = runner.invoke(cli, ["complement", "ATGC"])

        assert result.exit_code == 0
        assert "互补链:" in result.output
        assert "TACG" in result.output

    def test_reverse_complement(self, runner: CliRunner):
        """使用 -r 应该返回反向互补链。"""
        result = runner.invoke(cli, ["complement", "ATGC", "-r"])

        assert result.exit_code == 0
        assert "反向互补链:" in result.output
        assert "GCAT" in result.output

    def test_reverse_complement_long_option(self, runner: CliRunner):
        """使用 --reverse 应该返回反向互补链。"""
        result = runner.invoke(cli, ["complement", "ATGC", "--reverse"])

        assert result.exit_code == 0
        assert "GCAT" in result.output

    def test_complement_invalid_sequence(self, runner: CliRunner):
        """无效序列应该报错。"""
        result = runner.invoke(cli, ["complement", "ATGX"])

        assert result.exit_code != 0


class TestTranslateCommand:
    """测试 translate 命令。"""

    def test_translate_with_stop(self, runner: CliRunner):
        """翻译到终止密码子。"""
        # ATGTTTTAA = M, F, Stop
        result = runner.invoke(cli, ["translate", "ATGTTTTAA"])

        assert result.exit_code == 0
        assert "DNA:" in result.output
        assert "RNA:" in result.output
        assert "蛋白质:" in result.output
        assert "MF" in result.output

    def test_translate_shows_stop_codon(self, runner: CliRunner):
        """应该显示是否遇到终止密码子。"""
        result = runner.invoke(cli, ["translate", "ATGTTTTAA"])

        assert result.exit_code == 0
        assert "遇到终止密码子:" in result.output

    def test_translate_shows_properties(self, runner: CliRunner):
        """应该显示蛋白质性质。"""
        result = runner.invoke(cli, ["translate", "ATGAAATTT"])

        assert result.exit_code == 0
        assert "蛋白质性质:" in result.output
        assert "分子量:" in result.output

    def test_translate_invalid_sequence(self, runner: CliRunner):
        """无效序列应该报错。"""
        result = runner.invoke(cli, ["translate", "ATGX"])

        assert result.exit_code != 0


class TestAnalyzeCommand:
    """测试 analyze 命令。"""

    def test_analyze_fasta(self, runner: CliRunner, tmp_path: Path):
        """分析 FASTA 文件。"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1 Test\nATGCGATCGATCG\n")

        result = runner.invoke(cli, ["analyze", "-f", str(fasta_file)])

        assert result.exit_code == 0
        assert "序列: seq1" in result.output
        assert "长度:" in result.output
        assert "GC 含量:" in result.output

    def test_analyze_multiple_sequences(self, runner: CliRunner, tmp_path: Path):
        """分析多条序列。"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1 First\nATGC\n>seq2 Second\nGCTA\n")

        result = runner.invoke(cli, ["analyze", "-f", str(fasta_file)])

        assert result.exit_code == 0
        assert "seq1" in result.output
        assert "seq2" in result.output

    def test_analyze_with_motif(self, runner: CliRunner, tmp_path: Path):
        """搜索 motif。"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1 Test\nATGATGATG\n")

        result = runner.invoke(cli, ["analyze", "-f", str(fasta_file), "-m", "ATG"])

        assert result.exit_code == 0
        assert "Motif 'ATG'" in result.output

    def test_analyze_motif_not_found(self, runner: CliRunner, tmp_path: Path):
        """未找到 motif。"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1 Test\nAAAAAAAA\n")

        result = runner.invoke(cli, ["analyze", "-f", str(fasta_file), "-m", "GGG"])

        assert result.exit_code == 0
        assert "未找到 motif" in result.output

    def test_analyze_file_not_found(self, runner: CliRunner, tmp_path: Path):
        """文件不存在应该报错。"""
        result = runner.invoke(
            cli, ["analyze", "-f", str(tmp_path / "nonexistent.fasta")]
        )

        assert result.exit_code != 0


class TestVersionOption:
    """测试版本选项。"""

    def test_version(self, runner: CliRunner):
        """--version 应该显示版本。"""
        result = runner.invoke(cli, ["--version"])

        assert result.exit_code == 0
        assert "genomeflow" in result.output
        assert "0.1.0" in result.output


class TestHelpOption:
    """测试帮助选项。"""

    def test_main_help(self, runner: CliRunner):
        """主命令帮助。"""
        result = runner.invoke(cli, ["--help"])

        assert result.exit_code == 0
        assert "GenomeFlow" in result.output
        assert "quick" in result.output
        assert "analyze" in result.output

    def test_quick_help(self, runner: CliRunner):
        """quick 命令帮助。"""
        result = runner.invoke(cli, ["quick", "--help"])

        assert result.exit_code == 0
        assert "快速分析" in result.output

    def test_analyze_help(self, runner: CliRunner):
        """analyze 命令帮助。"""
        result = runner.invoke(cli, ["analyze", "--help"])

        assert result.exit_code == 0
        assert "FASTA" in result.output
