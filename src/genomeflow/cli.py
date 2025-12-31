"""
命令行接口。

使用 Click 框架构建 CLI，它是 Python 中最流行的 CLI 库之一。

为什么选择 Click？
1. 装饰器语法简洁直观
2. 自动生成帮助信息
3. 参数验证和类型转换
4. 支持命令组和子命令
5. 彩色终端输出
6. 被 Flask 等知名项目采用
"""

from __future__ import annotations

from pathlib import Path

import click

from genomeflow.analyzer import (
    base_frequency,
    calculate_molecular_weight,
    find_motif,
    gc_content,
    transcribe,
)
from genomeflow.io import read_fasta
from genomeflow.sequence import DNASequence, InvalidSequenceError


# 创建命令组
# @click.group() 将函数变成一个可以包含多个子命令的组
@click.group()
@click.version_option(version="0.1.0", prog_name="genomeflow")
def cli() -> None:
    """GenomeFlow - DNA 序列分析工具

    一个面向生物信息学初学者的命令行工具。

    示例：

        genomeflow quick ATGCGATCG

        genomeflow analyze -f sequences.fasta

        genomeflow complement ATGC --reverse
    """
    pass


@cli.command()
@click.argument("sequence")
def quick(sequence: str) -> None:
    """快速分析单条 DNA 序列。

    SEQUENCE: 要分析的 DNA 序列字符串
    """
    try:
        seq = DNASequence(sequence)

        click.echo(f"序列: {seq}")
        click.echo(f"长度: {len(seq)} bp")
        click.echo(f"GC 含量: {gc_content(seq):.2%}")

        freq = base_frequency(seq)
        click.echo(
            f"碱基组成: A={freq.a_count} ({freq.a_ratio:.1%}) "
            f"T={freq.t_count} ({freq.t_ratio:.1%}) "
            f"G={freq.g_count} ({freq.g_ratio:.1%}) "
            f"C={freq.c_count} ({freq.c_ratio:.1%})"
        )

        click.echo(f"RNA 转录: {transcribe(seq)}")
        click.echo(f"反向互补: {seq.reverse_complement()}")

    except InvalidSequenceError as e:
        # click.style 可以给输出添加颜色
        click.echo(click.style(f"错误: {e}", fg="red"), err=True)
        raise SystemExit(1)


@cli.command()
@click.option(
    "-f",
    "--file",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="FASTA 文件路径",
)
@click.option(
    "-m",
    "--motif",
    type=str,
    default=None,
    help="要搜索的 motif（可选）",
)
def analyze(file: Path, motif: str | None) -> None:
    """分析 FASTA 文件中的序列。

    读取 FASTA 文件，对每条序列进行分析，
    输出 GC 含量、碱基组成、分子量等信息。
    """
    try:
        for record in read_fasta(file):
            click.echo(f"\n{'=' * 50}")
            click.echo(click.style(f"序列: {record.id}", fg="cyan", bold=True))
            click.echo(f"描述: {record.description}")
            click.echo(f"长度: {len(record.sequence)} bp")

            # GC 含量
            gc = gc_content(record.sequence)
            # 根据 GC 含量显示不同颜色
            gc_color = "green" if 0.4 <= gc <= 0.6 else "yellow"
            click.echo(f"GC 含量: " + click.style(f"{gc:.2%}", fg=gc_color))

            # 碱基频率
            freq = base_frequency(record.sequence)
            click.echo(
                f"碱基组成: A={freq.a_count} T={freq.t_count} "
                f"G={freq.g_count} C={freq.c_count}"
            )

            # 分子量
            mw = calculate_molecular_weight(record.sequence)
            click.echo(f"分子量: {mw:.1f} Da")

            # 如果指定了 motif，搜索它
            if motif:
                positions = find_motif(record.sequence, motif)
                if positions:
                    click.echo(
                        f"Motif '{motif}' 位置: "
                        + click.style(str(positions), fg="green")
                    )
                else:
                    click.echo(f"未找到 motif '{motif}'")

        click.echo(f"\n{'=' * 50}")

    except InvalidSequenceError as e:
        click.echo(click.style(f"错误: {e}", fg="red"), err=True)
        raise SystemExit(1)


@cli.command()
@click.argument("sequence")
@click.option(
    "-r",
    "--reverse",
    is_flag=True,
    help="返回反向互补链（而非互补链）",
)
def complement(sequence: str, reverse: bool) -> None:
    """获取 DNA 序列的互补链或反向互补链。

    SEQUENCE: 要处理的 DNA 序列字符串
    """
    try:
        seq = DNASequence(sequence)

        if reverse:
            result = seq.reverse_complement()
            click.echo(f"反向互补链: {result}")
        else:
            result = seq.complement()
            click.echo(f"互补链: {result}")

    except InvalidSequenceError as e:
        click.echo(click.style(f"错误: {e}", fg="red"), err=True)
        raise SystemExit(1)


@cli.command()
@click.argument("sequence")
def translate(sequence: str) -> None:
    """将 DNA 序列翻译为蛋白质。

    SEQUENCE: 要翻译的 DNA 序列字符串
    """
    try:
        dna = DNASequence(sequence)
        rna = dna.transcribe()
        result = rna.translate()

        click.echo(f"DNA: {dna}")
        click.echo(f"RNA: {rna}")
        click.echo(f"蛋白质: {click.style(result.protein, fg='green', bold=True)}")
        click.echo(f"遇到终止密码子: {'是' if result.stop_codon else '否'}")

        if result.protein:
            from genomeflow.protein import ProteinSequence

            protein = ProteinSequence(result.protein.replace("*", ""))
            props = protein.get_properties()
            click.echo(f"\n蛋白质性质:")
            click.echo(f"  分子量: {props.molecular_weight:.1f} Da")
            click.echo(f"  等电点: {props.isoelectric_point:.2f}")
            click.echo(f"  GRAVY: {props.gravy:.2f}")

    except InvalidSequenceError as e:
        click.echo(click.style(f"错误: {e}", fg="red"), err=True)
        raise SystemExit(1)


# 入口点函数
def main() -> None:
    """CLI 入口点。"""
    cli()


if __name__ == "__main__":
    main()
