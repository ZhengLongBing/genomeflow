# 教程 02：RNA 序列支持

## 📋 本章导览

- **你将掌握的技能**：
  - 理解 RNA 与 DNA 的区别
  - 设计可复用的序列基类
  - 实现 RNASequence 类
  - 实现 DNA 到 RNA 的转录
  - 理解密码子和翻译过程

- **前置知识**：
  - 完成教程 01
  - 理解 Python 类继承

- **核心挑战**：如何设计统一的序列抽象，同时保持各序列类型的特性？

---

## 📚 理解 RNA

### RNA 与 DNA 的区别

你可能会想："RNA 和 DNA 有什么不同？为什么要单独实现？"

| 特性 | DNA | RNA |
|------|-----|-----|
| 糖 | 脱氧核糖 | 核糖 |
| 碱基 | A, T, G, C | A, **U**, G, C |
| 结构 | 通常双链 | 通常单链 |
| 稳定性 | 稳定 | 较不稳定 |
| 功能 | 遗传信息存储 | 信息传递、催化 |

关键区别：**RNA 用 U（尿嘧啶）代替 T（胸腺嘧啶）**

```
DNA: ATGCGATCG
     ↓ 转录
RNA: AUGCGAUCG
```

### RNA 的类型

| 类型 | 全名 | 功能 |
|------|------|------|
| mRNA | 信使 RNA | 携带遗传信息到核糖体 |
| tRNA | 转运 RNA | 运输氨基酸 |
| rRNA | 核糖体 RNA | 组成核糖体 |
| miRNA | 微小 RNA | 基因调控 |

### 密码子表

RNA 以三个碱基（密码子）为单位编码氨基酸：

```
        第二位
      U     C     A     G
   ┌─────┬─────┬─────┬─────┐
 U │ Phe │ Ser │ Tyr │ Cys │ U
   │ Phe │ Ser │ Tyr │ Cys │ C
第 │ Leu │ Ser │ STOP│ STOP│ A  第
一 │ Leu │ Ser │ STOP│ Trp │ G  三
位 ├─────┼─────┼─────┼─────┤    位
 C │ Leu │ Pro │ His │ Arg │ U
   │ Leu │ Pro │ His │ Arg │ C
   │ Leu │ Pro │ Gln │ Arg │ A
   │ Leu │ Pro │ Gln │ Arg │ G
   ├─────┼─────┼─────┼─────┤
 A │ Ile │ Thr │ Asn │ Ser │ U
   │ Ile │ Thr │ Asn │ Ser │ C
   │ Ile │ Thr │ Lys │ Arg │ A
   │ Met*│ Thr │ Lys │ Arg │ G
   ├─────┼─────┼─────┼─────┤
 G │ Val │ Ala │ Asp │ Gly │ U
   │ Val │ Ala │ Asp │ Gly │ C
   │ Val │ Ala │ Glu │ Gly │ A
   │ Val │ Ala │ Glu │ Gly │ G
   └─────┴─────┴─────┴─────┘

* AUG 是起始密码子，编码 Met（甲硫氨酸）
STOP: UAA, UAG, UGA 是终止密码子
```

---

## 🔧 设计我们的实现

### 架构思考

**问题**：DNA 和 RNA 有很多共同点（都是核酸序列），如何避免代码重复？

**方案对比**：

| 方案 | 优点 | 缺点 |
|------|------|------|
| 复制 DNASequence 代码 | 简单直接 | 重复代码，维护困难 |
| 继承自共同基类 | 复用代码，统一接口 | 需要设计抽象 |
| 组合模式 | 灵活 | 可能过度设计 |
| 泛型/协议 | 类型安全 | Python 中较复杂 |

**我们的选择**：使用抽象基类 + 继承

```
           BaseSequence (抽象基类)
           /           \
    DNASequence     RNASequence
```

### 设计决策

**1. 哪些应该放在基类？**
- 序列存储和验证框架
- 长度、索引、迭代等通用操作
- 相等性比较

**2. 哪些应该在子类实现？**
- 有效碱基定义（ATGC vs AUGC）
- 互补规则
- 特定转换（DNA→RNA 转录）

### 架构图

```
┌─────────────────────────────────────────────────────────┐
│              BaseSequence (抽象基类)                     │
│  ├── _sequence: str                                     │
│  ├── __len__(), __getitem__(), __iter__()              │
│  ├── __eq__(), __hash__()                              │
│  └── @abstractmethod: VALID_BASES, complement()        │
└─────────────────────────┬───────────────────────────────┘
                          │
          ┌───────────────┴───────────────┐
          ▼                               ▼
┌─────────────────────┐       ┌─────────────────────┐
│    DNASequence      │       │    RNASequence      │
│  VALID_BASES: ATGC  │       │  VALID_BASES: AUGC  │
│  complement()       │       │  complement()       │
│  transcribe()→RNA   │       │  reverse_transcribe()│
└─────────────────────┘       └─────────────────────┘
```

---

## 💻 代码实现

### 步骤 1：重构为抽象基类

创建 `src/genomeflow/base.py`：

```python
"""
序列的抽象基类。

设计理念：
- 提取 DNA 和 RNA 的共同行为到基类
- 使用 Python 的 ABC（Abstract Base Class）机制
- 子类只需实现特定的碱基规则和操作
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Iterator


class InvalidSequenceError(ValueError):
    """当序列包含无效字符时抛出此异常。"""
    pass


class BaseSequence(ABC):
    """
    核酸序列的抽象基类。

    这是一个抽象类，不能直接实例化。
    子类必须实现 VALID_BASES 和 complement() 方法。

    为什么使用抽象基类？
    1. 强制子类实现必要的方法
    2. 共享通用实现，减少代码重复
    3. 提供统一的接口，便于多态使用
    """

    @property
    @abstractmethod
    def VALID_BASES(self) -> frozenset[str]:
        """
        有效的碱基集合。

        子类必须实现此属性。
        DNA 返回 frozenset("ATGC")
        RNA 返回 frozenset("AUGC")
        """
        ...

    @property
    @abstractmethod
    def COMPLEMENT_MAP(self) -> dict[str, str]:
        """
        碱基配对映射。

        子类必须实现此属性。
        DNA: A↔T, G↔C
        RNA: A↔U, G↔C
        """
        ...

    def __init__(self, sequence: str) -> None:
        """
        创建序列对象。

        Args:
            sequence: 序列字符串（不区分大小写）

        Raises:
            InvalidSequenceError: 如果序列包含无效字符
        """
        normalized = sequence.upper()
        self._validate(normalized)
        self._sequence = normalized

    def _validate(self, sequence: str) -> None:
        """验证序列是否只包含有效碱基。"""
        invalid_chars = set(sequence) - self.VALID_BASES
        if invalid_chars:
            raise InvalidSequenceError(
                f"序列包含无效字符: {invalid_chars}。"
                f"有效字符: {self.VALID_BASES}"
            )

    @property
    def sequence(self) -> str:
        """返回原始序列字符串。"""
        return self._sequence

    def __len__(self) -> int:
        return len(self._sequence)

    def __getitem__(self, index: int | slice) -> str:
        return self._sequence[index]

    def __iter__(self) -> Iterator[str]:
        return iter(self._sequence)

    def __str__(self) -> str:
        return self._sequence

    def __repr__(self) -> str:
        class_name = self.__class__.__name__
        if len(self._sequence) > 50:
            display = f"{self._sequence[:25]}...{self._sequence[-25:]}"
        else:
            display = self._sequence
        return f"{class_name}('{display}')"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, BaseSequence):
            return (
                type(self) == type(other) and
                self._sequence == other._sequence
            )
        if isinstance(other, str):
            return self._sequence == other.upper()
        return NotImplemented

    def __hash__(self) -> int:
        return hash((type(self).__name__, self._sequence))

    def complement(self) -> BaseSequence:
        """
        返回互补链。

        Returns:
            同类型的互补序列
        """
        comp_seq = "".join(self.COMPLEMENT_MAP[base] for base in self._sequence)
        return type(self)(comp_seq)

    def reverse_complement(self) -> BaseSequence:
        """
        返回反向互补链。

        Returns:
            同类型的反向互补序列
        """
        comp_seq = "".join(self.COMPLEMENT_MAP[base] for base in self._sequence)
        return type(self)(comp_seq[::-1])

    def gc_content(self) -> float:
        """
        计算 GC 含量。

        GC 含量是 G 和 C 碱基的百分比，对 DNA 和 RNA 都有意义。

        Returns:
            GC 含量，范围 0.0-1.0
        """
        if len(self._sequence) == 0:
            return 0.0
        gc_count = self._sequence.count("G") + self._sequence.count("C")
        return gc_count / len(self._sequence)
```

### 步骤 2：重构 DNASequence

更新 `src/genomeflow/sequence.py`：

```python
"""
DNA 序列类。

继承自 BaseSequence，实现 DNA 特定的碱基规则和操作。
"""

from __future__ import annotations

from genomeflow.base import BaseSequence, InvalidSequenceError

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
```

### 步骤 3：实现 RNASequence

创建 `src/genomeflow/rna.py`：

```python
"""
RNA 序列类。

RNA（核糖核酸）是遗传信息传递的中间分子。
与 DNA 的主要区别是用 U（尿嘧啶）代替 T（胸腺嘧啶）。
"""

from __future__ import annotations

from dataclasses import dataclass

from genomeflow.base import BaseSequence


# 标准遗传密码表
# 密码子 -> 氨基酸单字母代码
CODON_TABLE: dict[str, str] = {
    # 苯丙氨酸 (Phe, F)
    "UUU": "F", "UUC": "F",
    # 亮氨酸 (Leu, L)
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    # 异亮氨酸 (Ile, I)
    "AUU": "I", "AUC": "I", "AUA": "I",
    # 甲硫氨酸 (Met, M) - 起始密码子
    "AUG": "M",
    # 缬氨酸 (Val, V)
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    # 丝氨酸 (Ser, S)
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    # 脯氨酸 (Pro, P)
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # 苏氨酸 (Thr, T)
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # 丙氨酸 (Ala, A)
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # 酪氨酸 (Tyr, Y)
    "UAU": "Y", "UAC": "Y",
    # 终止密码子
    "UAA": "*", "UAG": "*", "UGA": "*",
    # 组氨酸 (His, H)
    "CAU": "H", "CAC": "H",
    # 谷氨酰胺 (Gln, Q)
    "CAA": "Q", "CAG": "Q",
    # 天冬酰胺 (Asn, N)
    "AAU": "N", "AAC": "N",
    # 赖氨酸 (Lys, K)
    "AAA": "K", "AAG": "K",
    # 天冬氨酸 (Asp, D)
    "GAU": "D", "GAC": "D",
    # 谷氨酸 (Glu, E)
    "GAA": "E", "GAG": "E",
    # 半胱氨酸 (Cys, C)
    "UGU": "C", "UGC": "C",
    # 色氨酸 (Trp, W)
    "UGG": "W",
    # 精氨酸 (Arg, R)
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # 甘氨酸 (Gly, G)
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
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

    def reverse_transcribe(self) -> "DNASequence":
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
            codon = seq[i:i + 3]
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
                codon = seq[i:i + 3]

                # 找到起始密码子
                if codon == "AUG":
                    start = i
                    protein_parts: list[str] = []

                    # 继续翻译直到终止密码子
                    j = i
                    while j < len(seq) - 2:
                        codon = seq[j:j + 3]
                        amino_acid = CODON_TABLE.get(codon, "X")

                        if amino_acid == "*":
                            # 找到完整的 ORF
                            if len(protein_parts) >= min_length:
                                orfs.append((
                                    start,
                                    j + 3,
                                    "".join(protein_parts),
                                ))
                            break
                        else:
                            protein_parts.append(amino_acid)
                        j += 3

                    i = j + 3  # 跳过已处理的部分
                else:
                    i += 3

        return orfs
```

### 步骤 4：更新包导出

更新 `src/genomeflow/__init__.py`：

```python
"""
GenomeFlow - 序列分析工具

支持 DNA、RNA 和蛋白质序列的分析。
"""

from genomeflow.base import BaseSequence, InvalidSequenceError
from genomeflow.sequence import DNASequence
from genomeflow.rna import RNASequence, TranslationResult, CODON_TABLE
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
    # 基类
    "BaseSequence",
    "InvalidSequenceError",
    # DNA
    "DNASequence",
    # RNA
    "RNASequence",
    "TranslationResult",
    "CODON_TABLE",
    # 分析
    "gc_content",
    "base_frequency",
    "find_motif",
    "transcribe",
    "calculate_molecular_weight",
    "BaseFrequency",
    # IO
    "read_fasta",
    "write_fasta",
    "FastaRecord",
]
```

### 步骤 5：编写测试

创建 `tests/test_rna.py`：

```python
"""RNA 序列的单元测试。"""

import pytest

from genomeflow.rna import RNASequence, TranslationResult
from genomeflow.sequence import DNASequence
from genomeflow.base import InvalidSequenceError


class TestRNASequenceCreation:
    """测试 RNASequence 的创建。"""

    def test_create_valid_sequence(self):
        seq = RNASequence("AUGC")
        assert str(seq) == "AUGC"

    def test_lowercase_converted(self):
        seq = RNASequence("augc")
        assert str(seq) == "AUGC"

    def test_invalid_thymine_rejected(self):
        """RNA 不应该包含 T。"""
        with pytest.raises(InvalidSequenceError) as exc_info:
            RNASequence("ATGC")  # T 是无效的
        assert "T" in str(exc_info.value)

    def test_valid_uracil(self):
        """RNA 应该接受 U。"""
        seq = RNASequence("AUUU")
        assert "U" in str(seq)


class TestRNAComplement:
    """测试 RNA 互补链。"""

    def test_complement(self):
        seq = RNASequence("AUGC")
        comp = seq.complement()
        assert str(comp) == "UACG"
        assert isinstance(comp, RNASequence)

    def test_reverse_complement(self):
        seq = RNASequence("AUGC")
        rev_comp = seq.reverse_complement()
        assert str(rev_comp) == "GCAU"


class TestTranscription:
    """测试 DNA 和 RNA 之间的转换。"""

    def test_dna_to_rna(self):
        dna = DNASequence("ATGC")
        rna = dna.transcribe()
        assert str(rna) == "AUGC"
        assert isinstance(rna, RNASequence)

    def test_rna_to_dna(self):
        rna = RNASequence("AUGC")
        dna = rna.reverse_transcribe()
        assert str(dna) == "ATGC"
        assert isinstance(dna, DNASequence)

    def test_round_trip(self):
        """DNA → RNA → DNA 应该得到原序列。"""
        original = DNASequence("ATGCGATCGATCG")
        result = original.transcribe().reverse_transcribe()
        assert result == original


class TestTranslation:
    """测试 RNA 翻译。"""

    def test_simple_translation(self):
        # AUG = Met, UUU = Phe, UAA = Stop
        rna = RNASequence("AUGUUUUAA")
        result = rna.translate()
        assert result.protein == "MF"
        assert result.stop_codon is True

    def test_no_start_codon(self):
        rna = RNASequence("UUUUUU")
        result = rna.translate(start_codon=True)
        assert result.protein == ""

    def test_translate_without_start_codon_check(self):
        # UUU = Phe
        rna = RNASequence("UUUUUU")
        result = rna.translate(start_codon=False)
        assert result.protein == "FF"

    def test_translate_with_stop_in_middle(self):
        # AUG=M, UAA=Stop, UUU=F
        rna = RNASequence("AUGUAAUUU")
        result = rna.translate(stop_at_stop=True)
        assert result.protein == "M"
        assert result.stop_codon is True

    def test_translate_through_stop(self):
        rna = RNASequence("AUGUAAUUU")
        result = rna.translate(stop_at_stop=False)
        assert "*" in result.protein


class TestORFFinding:
    """测试 ORF 查找。"""

    def test_find_single_orf(self):
        # AUG AAA UAA = Met Lys Stop
        rna = RNASequence("AUGAAAUAA")
        orfs = rna.find_orfs(min_length=1)
        assert len(orfs) == 1
        assert orfs[0][2] == "MK"

    def test_no_orf_without_stop(self):
        rna = RNASequence("AUGAAAAAA")
        orfs = rna.find_orfs(min_length=1)
        assert len(orfs) == 0  # 没有终止密码子

    def test_min_length_filter(self):
        rna = RNASequence("AUGAAAUAA")  # 只有 2 个氨基酸
        orfs = rna.find_orfs(min_length=5)
        assert len(orfs) == 0
```

---

## ✅ 测试和验证

### 运行测试

```bash
# 运行所有测试
uv run pytest

# 只运行 RNA 测试
uv run pytest tests/test_rna.py -v

# 查看测试覆盖率
uv run pytest --cov=genomeflow
```

### 手动验证

```python
from genomeflow import DNASequence, RNASequence

# DNA 转录为 RNA
dna = DNASequence("ATGAAATTTGGG")
rna = dna.transcribe()
print(f"DNA: {dna}")
print(f"RNA: {rna}")

# RNA 翻译为蛋白质
result = rna.translate()
print(f"蛋白质: {result.protein}")
print(f"遇到终止密码子: {result.stop_codon}")

# 查找 ORF
long_rna = RNASequence("UUUAUGAAAGGGCCCUAAUUUAUGCCCUAG")
orfs = long_rna.find_orfs(min_length=2)
for start, end, protein in orfs:
    print(f"ORF [{start}:{end}]: {protein}")
```

---

## 🤔 深入思考

<details>
<summary>为什么使用抽象基类而不是 Protocol？</summary>

Python 3.8+ 引入了 `typing.Protocol`，支持结构化子类型（鸭子类型的形式化）。

**Protocol 方式**：
```python
class SequenceProtocol(Protocol):
    def __len__(self) -> int: ...
    def complement(self) -> "SequenceProtocol": ...
```

**ABC 方式**（我们的选择）：
```python
class BaseSequence(ABC):
    @abstractmethod
    def complement(self) -> "BaseSequence": ...
```

**选择 ABC 的原因**：
1. 可以共享实现代码（`__len__`, `__getitem__` 等）
2. 强制继承关系，`isinstance()` 检查有效
3. 对于新手更容易理解

**何时用 Protocol**：
- 需要支持不能修改的第三方类
- 更看重灵活性而非代码共享

</details>

<details>
<summary>为什么 translate() 返回 dataclass 而不是字符串？</summary>

翻译结果不仅仅是氨基酸序列，还包含：
- 是否遇到终止密码子
- 剩余未翻译的碱基数

这些信息对于分析很重要。使用 dataclass：
- 类型安全，IDE 提示友好
- 不可变（frozen=True）
- 容易扩展（将来可以加更多字段）

</details>

---

## 📝 总结

通过这个教程，你学会了：

1. **重构技巧**：将共同行为提取到抽象基类
2. **RNA 特性**：U 代替 T，以及与 DNA 的转换
3. **翻译机制**：密码子表和蛋白质合成
4. **ORF 查找**：识别可能的编码区域

下一步：[教程 03 - 蛋白质序列支持](tutorial-03-protein-sequence.md)
