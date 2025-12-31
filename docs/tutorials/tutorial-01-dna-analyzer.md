# 教程：使用 Python 编写 DNA 序列分析器

## 📋 本章导览

- **你将掌握的技能**：
  - 理解 DNA 序列的基本概念和分析意义
  - 使用 uv 管理 Python 项目
  - 编写模块化、类型安全的 Python 代码
  - 实现 GC 含量计算、碱基统计、互补链生成等核心功能
  - 使用现代 Python 工具链进行测试

- **预计学习时间**：2-3 小时

- **前置知识**：
  - 基础 Python 语法（变量、函数、循环）
  - 会使用命令行终端
  - 不需要生物学背景（本教程会讲解必要概念）

- **核心挑战**：如何将生物学概念转化为清晰的代码结构？

---

## 📚 理解 DNA 序列分析

### 什么是 DNA？为什么要分析它？

你可能会想："DNA 是什么？为什么程序员要关心它？"

想象一下，DNA 就像是生命的"源代码"。就像软件由 0 和 1 组成，生命由 4 个"字符"组成：

| 碱基 | 符号 | 全名 |
|------|------|------|
| 腺嘌呤 | A | Adenine |
| 胸腺嘧啶 | T | Thymine |
| 鸟嘌呤 | G | Guanine |
| 胞嘧啶 | C | Cytosine |

一段 DNA 序列看起来就像这样：
```
ATGCGATCGATCGATCGATCG
```

**为什么要分析 DNA 序列？**

1. **疾病诊断**：某些基因突变与疾病相关
2. **物种鉴定**：通过 DNA 判断物种
3. **进化研究**：比较不同物种的 DNA 相似度
4. **药物开发**：理解基因如何影响药物反应

### 核心概念：我们要分析什么？

作为生物信息学入门，我们将实现以下分析功能：

#### 1. GC 含量（GC Content）

GC 含量是序列中 G 和 C 碱基的百分比。

```
序列：ATGCGC
G 的数量：2
C 的数量：2
总长度：6
GC 含量：(2+2)/6 = 66.7%
```

**为什么重要？**
- GC 含量高的 DNA 更稳定（G-C 配对有 3 个氢键，A-T 只有 2 个）
- 不同物种有特征性的 GC 含量
- PCR 实验设计需要考虑 GC 含量

#### 2. 碱基频率统计

统计每种碱基出现的次数和比例：

```
序列：ATGCGATCG
A: 2 (22.2%)
T: 2 (22.2%)
G: 3 (33.3%)
C: 2 (22.2%)
```

#### 3. 互补链（Complement）和反向互补链（Reverse Complement）

DNA 是双链结构，两条链通过碱基配对连接：
- A 配对 T
- G 配对 C

```
原始链：    5'-ATGC-3'
互补链：    3'-TACG-5'
反向互补链：5'-GCAT-3'（反向互补链方向翻转）
```

**为什么重要？**
- DNA 复制和转录依赖互补配对
- 引物设计必须考虑方向性
- 数据库搜索常用反向互补序列

<details>
<summary>💡 思考：为什么反向互补链要"反向"？</summary>

DNA 有方向性，称为 5' 端和 3' 端。两条链是反向平行的，就像双向车道。当我们说"互补链"时，如果只是简单替换碱基，方向是相反的。"反向互补"则是把方向调整成一致，这在序列比对中非常重要。

</details>

---

## 🔧 设计我们的实现

### 架构思考

在写代码之前，让我们思考几个问题：

**1. 如何表示 DNA 序列？**

| 方案 | 优点 | 缺点 |
|------|------|------|
| 字符串 `str` | 简单直观 | 无法保证只含 ATGC |
| 自定义类 | 可以验证、添加方法 | 稍复杂 |
| 枚举 + 列表 | 类型安全 | 使用不便 |

**我们的选择**：使用自定义类包装字符串，在构造时验证，兼顾简洁和安全。

**2. 如何组织代码模块？**

```
genomeflow/
├── __init__.py
├── sequence.py      # DNA 序列类
├── analyzer.py      # 分析功能
├── io.py           # 文件读写（FASTA 格式）
└── cli.py          # 命令行接口
```

**3. 错误处理策略**

- 无效碱基：抛出明确的异常
- 空序列：根据场景返回 0 或抛异常
- 文件不存在：使用标准 FileNotFoundError

### 架构图

```
┌─────────────────────────────────────────────────────────┐
│                      CLI (cli.py)                       │
│                   命令行参数解析                          │
└─────────────────────────┬───────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│                   Analyzer (analyzer.py)                │
│        gc_content() | base_frequency() | ...           │
└─────────────────────────┬───────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│                   DNASequence (sequence.py)             │
│              序列存储 | 验证 | 基本操作                    │
└─────────────────────────┬───────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│                      IO (io.py)                         │
│                  FASTA 文件读写                          │
└─────────────────────────────────────────────────────────┘
```

### 设计权衡

**面向对象 vs 函数式？**

我们采用混合方式：
- `DNASequence` 类：封装数据和基本操作
- 分析函数：独立函数，接收序列作为参数

这样做的好处是：
- 类提供清晰的数据边界和验证
- 函数便于组合和测试
- 符合 Python 的实用主义风格

---

## 💻 代码实现

### 步骤 1：使用 uv 初始化项目

首先，确保你已安装 [uv](https://github.com/astral-sh/uv)——它是现代 Python 的项目管理工具，比 pip 快得多。

```bash
# 安装 uv（如果还没安装）
curl -LsSf https://astral.sh/uv/install.sh | sh

# 创建新项目（如果从零开始）
uv init genomeflow
cd genomeflow

# 或者在现有目录初始化
uv init
```

<details>
<summary>💡 为什么用 uv 而不是 pip？</summary>

1. **速度快**：用 Rust 编写，比 pip 快 10-100 倍
2. **依赖解析更好**：使用现代算法，减少冲突
3. **虚拟环境管理**：自动创建和管理 .venv
4. **锁文件支持**：确保环境可复现
5. **一体化**：替代 pip、pip-tools、virtualenv、pyenv 等多个工具

</details>

### 步骤 2：配置项目

编辑 `pyproject.toml`：

```toml
[project]
name = "genomeflow"
version = "0.1.0"
description = "A DNA sequence analyzer for bioinformatics beginners"
requires-python = ">=3.12"
dependencies = [
    "click>=8.1",  # CLI 框架
]

[project.optional-dependencies]
dev = [
    "pytest>=8.0",
    "pytest-cov>=4.0",
]

[project.scripts]
# 命令行入口点：输入 genomeflow 即可运行
genomeflow = "genomeflow.cli:main"

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"
```

安装依赖：

```bash
uv sync --extra dev
```

### 步骤 3：创建项目结构

```bash
mkdir -p src/genomeflow tests
touch src/genomeflow/__init__.py
touch src/genomeflow/sequence.py
touch src/genomeflow/analyzer.py
touch src/genomeflow/io.py
touch src/genomeflow/cli.py
touch tests/__init__.py
touch tests/test_sequence.py
touch tests/test_analyzer.py
```

更新 `pyproject.toml` 添加 src 布局支持：

```toml
[tool.hatch.build.targets.wheel]
packages = ["src/genomeflow"]
```

### 步骤 4：实现 DNASequence 类

创建 `src/genomeflow/sequence.py`：

```python
"""
DNA 序列的核心数据结构。

为什么需要自定义类而不是直接用字符串？
1. 验证：确保序列只包含有效碱基
2. 封装：提供统一的操作接口
3. 类型安全：IDE 和类型检查器能提供更好的支持
"""

from __future__ import annotations

# 使用 typing 模块提供类型提示
# 类型提示的好处：IDE 自动补全、提前发现错误、代码即文档
from typing import Iterator


class InvalidSequenceError(ValueError):
    """当 DNA 序列包含无效字符时抛出此异常。"""
    pass


class DNASequence:
    """
    表示一条 DNA 序列。

    DNA 序列只能包含四种碱基：A、T、G、C（不区分大小写）。
    内部统一存储为大写形式。

    示例：
        >>> seq = DNASequence("ATGC")
        >>> len(seq)
        4
        >>> seq[0]
        'A'
    """

    # 有效的 DNA 碱基
    # 使用 frozenset 因为：1) 查找是 O(1)  2) 不可变，更安全
    VALID_BASES: frozenset[str] = frozenset("ATGC")

    # 碱基配对规则：A-T, G-C
    COMPLEMENT_MAP: dict[str, str] = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
    }

    def __init__(self, sequence: str) -> None:
        """
        创建 DNA 序列对象。

        Args:
            sequence: DNA 序列字符串，只能包含 A、T、G、C（不区分大小写）

        Raises:
            InvalidSequenceError: 如果序列包含无效字符

        示例：
            >>> DNASequence("atgc")  # 小写会自动转大写
            DNASequence('ATGC')
        """
        # 统一转为大写，方便后续处理
        normalized = sequence.upper()

        # 验证序列
        # 这里我们选择在构造时验证，而非每次操作时验证
        # 权衡：构造时多花一点时间，换取后续操作的安全和高效
        self._validate(normalized)

        self._sequence = normalized

    def _validate(self, sequence: str) -> None:
        """验证序列是否只包含有效碱基。"""
        invalid_chars = set(sequence) - self.VALID_BASES
        if invalid_chars:
            raise InvalidSequenceError(
                f"序列包含无效字符: {invalid_chars}。"
                f"DNA 序列只能包含 A、T、G、C。"
            )

    @property
    def sequence(self) -> str:
        """返回原始序列字符串。"""
        return self._sequence

    def __len__(self) -> int:
        """返回序列长度。"""
        return len(self._sequence)

    def __getitem__(self, index: int | slice) -> str:
        """支持索引和切片访问。"""
        return self._sequence[index]

    def __iter__(self) -> Iterator[str]:
        """支持迭代访问每个碱基。"""
        return iter(self._sequence)

    def __str__(self) -> str:
        """返回序列字符串（用于 print）。"""
        return self._sequence

    def __repr__(self) -> str:
        """返回可用于重建对象的表示（用于调试）。"""
        # 如果序列太长，截断显示
        if len(self._sequence) > 50:
            display = f"{self._sequence[:25]}...{self._sequence[-25:]}"
        else:
            display = self._sequence
        return f"DNASequence('{display}')"

    def __eq__(self, other: object) -> bool:
        """判断两个序列是否相等。"""
        if isinstance(other, DNASequence):
            return self._sequence == other._sequence
        if isinstance(other, str):
            return self._sequence == other.upper()
        return NotImplemented

    def __hash__(self) -> int:
        """使序列可哈希，可用作字典键。"""
        return hash(self._sequence)

    def complement(self) -> DNASequence:
        """
        返回互补链。

        DNA 双链中，A 与 T 配对，G 与 C 配对。

        示例：
            >>> DNASequence("ATGC").complement()
            DNASequence('TACG')
        """
        comp_seq = "".join(self.COMPLEMENT_MAP[base] for base in self._sequence)
        return DNASequence(comp_seq)

    def reverse_complement(self) -> DNASequence:
        """
        返回反向互补链。

        先取互补，再反转方向。这在序列比对中非常常用，
        因为 DNA 双链方向相反（5'->3' 和 3'->5'）。

        示例：
            >>> DNASequence("ATGC").reverse_complement()
            DNASequence('GCAT')
        """
        comp_seq = "".join(self.COMPLEMENT_MAP[base] for base in self._sequence)
        return DNASequence(comp_seq[::-1])
```

<details>
<summary>💡 思考：为什么要实现这么多魔术方法（__len__, __getitem__ 等）？</summary>

这些魔术方法让 `DNASequence` 的行为像原生 Python 类型一样自然：

```python
seq = DNASequence("ATGC")

# 因为实现了 __len__
print(len(seq))  # 4

# 因为实现了 __getitem__
print(seq[0])    # 'A'
print(seq[1:3])  # 'TG'

# 因为实现了 __iter__
for base in seq:
    print(base)

# 因为实现了 __eq__ 和 __hash__
sequences = {seq: "my sequence"}
```

这是 Python 的"鸭子类型"哲学：如果它走起来像鸭子，叫起来像鸭子，那它就是鸭子。

</details>

### 步骤 5：实现分析功能

创建 `src/genomeflow/analyzer.py`：

```python
"""
DNA 序列分析功能。

设计决策：为什么用独立函数而不是 DNASequence 的方法？
1. 单一职责：DNASequence 负责数据，analyzer 负责分析
2. 可测试性：函数更容易单独测试
3. 可扩展性：添加新分析功能不需要修改 DNASequence
4. 函数式风格：更容易组合和复用
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from genomeflow.sequence import DNASequence


@dataclass(frozen=True)
class BaseFrequency:
    """
    碱基频率统计结果。

    使用 dataclass 的好处：
    1. 自动生成 __init__, __repr__, __eq__
    2. frozen=True 使其不可变，更安全
    3. 作为返回类型比字典更清晰
    """
    a_count: int
    t_count: int
    g_count: int
    c_count: int
    total: int

    @property
    def a_ratio(self) -> float:
        """A 的比例。"""
        return self.a_count / self.total if self.total > 0 else 0.0

    @property
    def t_ratio(self) -> float:
        """T 的比例。"""
        return self.t_count / self.total if self.total > 0 else 0.0

    @property
    def g_ratio(self) -> float:
        """G 的比例。"""
        return self.g_count / self.total if self.total > 0 else 0.0

    @property
    def c_ratio(self) -> float:
        """C 的比例。"""
        return self.c_count / self.total if self.total > 0 else 0.0


def gc_content(seq: DNASequence) -> float:
    """
    计算 GC 含量（G 和 C 碱基的百分比）。

    GC 含量是生物信息学中最常用的指标之一：
    - 物种鉴定：不同物种有特征性的 GC 含量
    - 基因预测：编码区和非编码区 GC 含量不同
    - 实验设计：影响 DNA 熔解温度和 PCR 条件

    Args:
        seq: DNA 序列

    Returns:
        GC 含量，范围 0.0-1.0（0%-100%）
        空序列返回 0.0

    示例：
        >>> gc_content(DNASequence("ATGC"))
        0.5
        >>> gc_content(DNASequence("GGCC"))
        1.0
    """
    if len(seq) == 0:
        return 0.0

    # 使用 Counter 统计，比手动循环更 Pythonic
    # Counter 返回一个字典，如 {'A': 2, 'T': 2, 'G': 1, 'C': 1}
    counts = Counter(seq)
    gc_count = counts.get("G", 0) + counts.get("C", 0)

    return gc_count / len(seq)


def base_frequency(seq: DNASequence) -> BaseFrequency:
    """
    统计各碱基的出现次数和比例。

    Args:
        seq: DNA 序列

    Returns:
        BaseFrequency 对象，包含各碱基的计数和比例

    示例：
        >>> freq = base_frequency(DNASequence("AATGC"))
        >>> freq.a_count
        2
        >>> freq.a_ratio
        0.4
    """
    counts = Counter(seq)

    return BaseFrequency(
        a_count=counts.get("A", 0),
        t_count=counts.get("T", 0),
        g_count=counts.get("G", 0),
        c_count=counts.get("C", 0),
        total=len(seq),
    )


def find_motif(seq: DNASequence, motif: str) -> list[int]:
    """
    在序列中查找 motif（短序列模式）的所有出现位置。

    Motif 是具有生物学意义的短序列模式，例如：
    - 转录因子结合位点
    - 限制性内切酶识别位点
    - 启动子序列

    Args:
        seq: 要搜索的 DNA 序列
        motif: 要查找的短序列模式

    Returns:
        所有匹配位置的列表（0-indexed）

    示例：
        >>> find_motif(DNASequence("ATGATGATG"), "ATG")
        [0, 3, 6]
    """
    # 验证 motif
    motif = motif.upper()
    invalid = set(motif) - DNASequence.VALID_BASES
    if invalid:
        raise ValueError(f"Motif 包含无效字符: {invalid}")

    positions: list[int] = []
    seq_str = seq.sequence
    motif_len = len(motif)

    # 滑动窗口查找
    # 为什么不用 str.find()？因为我们要找所有位置，包括重叠的
    for i in range(len(seq_str) - motif_len + 1):
        if seq_str[i:i + motif_len] == motif:
            positions.append(i)

    return positions


def transcribe(seq: DNASequence) -> str:
    """
    将 DNA 序列转录为 RNA 序列。

    转录是基因表达的第一步：DNA -> RNA -> 蛋白质
    规则很简单：将 T（胸腺嘧啶）替换为 U（尿嘧啶）

    Args:
        seq: DNA 序列

    Returns:
        RNA 序列字符串（包含 A、U、G、C）

    示例：
        >>> transcribe(DNASequence("ATGC"))
        'AUGC'
    """
    return seq.sequence.replace("T", "U")


def calculate_molecular_weight(seq: DNASequence) -> float:
    """
    计算单链 DNA 的近似分子量（单位：道尔顿 Da）。

    分子量在实验中很重要：
    - 凝胶电泳迁移率与分子量相关
    - 摩尔浓度计算需要分子量

    使用各碱基的平均分子量（已减去水分子，因为形成磷酸二酯键会脱水）

    Args:
        seq: DNA 序列

    Returns:
        分子量（道尔顿）
    """
    # 各碱基核苷酸的分子量（Da）
    # 这些是形成多核苷酸链时的有效分子量
    weights = {
        "A": 331.2,  # dAMP
        "T": 322.2,  # dTMP
        "G": 347.2,  # dGMP
        "C": 307.2,  # dCMP
    }

    return sum(weights[base] for base in seq)
```

<details>
<summary>💡 思考：为什么 BaseFrequency 用 @dataclass 而不是普通字典？</summary>

比较两种方式：

```python
# 方式 1：返回字典
def base_frequency_dict(seq):
    return {"a": 2, "t": 1, "g": 1, "c": 1, "total": 5}

result = base_frequency_dict(seq)
print(result["a"])  # 没有自动补全，可能拼错键名

# 方式 2：返回 dataclass
@dataclass
class BaseFrequency:
    a_count: int
    ...

result = base_frequency(seq)
print(result.a_count)  # IDE 自动补全，类型检查
```

dataclass 的优势：
1. **类型安全**：IDE 和 mypy 能检查属性名
2. **自动补全**：写 `result.` 时 IDE 会提示所有属性
3. **不可变性**：`frozen=True` 防止意外修改
4. **自文档化**：清晰地表达返回值结构

</details>

### 步骤 6：实现 FASTA 文件读写

FASTA 是生物信息学最常用的序列文件格式。创建 `src/genomeflow/io.py`：

```python
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
                f.write(seq_str[i:i + line_width] + "\n")
```

### 步骤 7：实现命令行接口

首先添加 click 依赖：

```bash
uv add click
```

创建 `src/genomeflow/cli.py`：

```python
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
    gc_content,
    find_motif,
    transcribe,
    calculate_molecular_weight,
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
    "-f", "--file",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="FASTA 文件路径",
)
@click.option(
    "-m", "--motif",
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
    "-r", "--reverse",
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


# 入口点函数
def main() -> None:
    """CLI 入口点。"""
    cli()


if __name__ == "__main__":
    main()
```

<details>
<summary>Click vs argparse：有什么区别？</summary>

| 特性 | argparse | Click |
|------|----------|-------|
| 语法 | 命令式 | 装饰器式 |
| 子命令 | 需要手动设置 | `@group` + `@command` |
| 参数验证 | 基础 | 丰富的内置类型 |
| 彩色输出 | 需要额外库 | 内置 `click.style()` |
| 帮助信息 | 自动生成 | 自动生成，更美观 |
| 测试支持 | 较难 | 内置 `CliRunner` |

**argparse 示例**：
```python
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", required=True)
args = parser.parse_args()
print(args.file)
```

**Click 示例**：
```python
@click.command()
@click.option("-f", "--file", required=True)
def main(file):
    print(file)
```

Click 的装饰器语法更接近函数签名，代码更易读。

</details>

### 步骤 8：完善包初始化

编辑 `src/genomeflow/__init__.py`：

```python
"""
GenomeFlow - DNA 序列分析工具

一个面向生物信息学初学者的 Python 库，用于：
- DNA 序列的基本操作（互补链、反向互补等）
- 序列分析（GC 含量、碱基频率、motif 搜索等）
- FASTA 文件读写
"""

from genomeflow.sequence import DNASequence, InvalidSequenceError
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

# 定义公开 API
# 这告诉 IDE 和工具哪些是公开接口
__all__ = [
    # 核心类
    "DNASequence",
    "InvalidSequenceError",
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
]
```

---

## ✅ 测试和验证

### 步骤 9：编写单元测试

创建 `tests/test_sequence.py`：

```python
"""DNASequence 类的单元测试。"""

import pytest

from genomeflow.sequence import DNASequence, InvalidSequenceError


class TestDNASequenceCreation:
    """测试 DNASequence 的创建和验证。"""

    def test_create_valid_sequence(self):
        """正常创建应该成功。"""
        seq = DNASequence("ATGC")
        assert str(seq) == "ATGC"

    def test_lowercase_converted_to_uppercase(self):
        """小写字母应该自动转为大写。"""
        seq = DNASequence("atgc")
        assert str(seq) == "ATGC"

    def test_mixed_case(self):
        """混合大小写应该正常处理。"""
        seq = DNASequence("AtGc")
        assert str(seq) == "ATGC"

    def test_invalid_character_raises_error(self):
        """包含无效字符应该抛出异常。"""
        with pytest.raises(InvalidSequenceError) as exc_info:
            DNASequence("ATGX")

        assert "无效字符" in str(exc_info.value)
        assert "X" in str(exc_info.value)

    def test_empty_sequence(self):
        """空序列应该是允许的。"""
        seq = DNASequence("")
        assert len(seq) == 0


class TestDNASequenceOperations:
    """测试 DNASequence 的各种操作。"""

    def test_len(self):
        """测试长度计算。"""
        assert len(DNASequence("ATGC")) == 4
        assert len(DNASequence("")) == 0

    def test_getitem(self):
        """测试索引访问。"""
        seq = DNASequence("ATGC")
        assert seq[0] == "A"
        assert seq[-1] == "C"
        assert seq[1:3] == "TG"

    def test_iteration(self):
        """测试迭代。"""
        seq = DNASequence("ATGC")
        bases = list(seq)
        assert bases == ["A", "T", "G", "C"]

    def test_equality(self):
        """测试相等性比较。"""
        seq1 = DNASequence("ATGC")
        seq2 = DNASequence("ATGC")
        seq3 = DNASequence("GGCC")

        assert seq1 == seq2
        assert seq1 != seq3
        assert seq1 == "ATGC"  # 与字符串比较

    def test_hash(self):
        """测试可哈希性。"""
        seq1 = DNASequence("ATGC")
        seq2 = DNASequence("ATGC")

        # 相同序列应该有相同哈希
        assert hash(seq1) == hash(seq2)

        # 可以用作字典键
        d = {seq1: "value"}
        assert d[seq2] == "value"


class TestComplement:
    """测试互补链功能。"""

    def test_complement(self):
        """测试互补链。"""
        seq = DNASequence("ATGC")
        comp = seq.complement()
        assert str(comp) == "TACG"

    def test_reverse_complement(self):
        """测试反向互补链。"""
        seq = DNASequence("ATGC")
        rev_comp = seq.reverse_complement()
        assert str(rev_comp) == "GCAT"

    def test_complement_is_idempotent(self):
        """互补的互补应该等于原序列。"""
        seq = DNASequence("ATGCGATCGA")
        assert seq.complement().complement() == seq

    def test_reverse_complement_twice(self):
        """反向互补两次应该等于原序列。"""
        seq = DNASequence("ATGCGATCGA")
        assert seq.reverse_complement().reverse_complement() == seq
```

创建 `tests/test_analyzer.py`：

```python
"""分析功能的单元测试。"""

import pytest

from genomeflow.sequence import DNASequence
from genomeflow.analyzer import (
    gc_content,
    base_frequency,
    find_motif,
    transcribe,
    calculate_molecular_weight,
)


class TestGCContent:
    """测试 GC 含量计算。"""

    def test_fifty_percent(self):
        """50% GC 含量。"""
        seq = DNASequence("ATGC")
        assert gc_content(seq) == 0.5

    def test_all_gc(self):
        """100% GC 含量。"""
        seq = DNASequence("GGCC")
        assert gc_content(seq) == 1.0

    def test_no_gc(self):
        """0% GC 含量。"""
        seq = DNASequence("AATT")
        assert gc_content(seq) == 0.0

    def test_empty_sequence(self):
        """空序列应返回 0。"""
        seq = DNASequence("")
        assert gc_content(seq) == 0.0

    def test_realistic_sequence(self):
        """测试更真实的序列。"""
        # 这是一段典型的人类基因序列片段
        seq = DNASequence("ATGCGATCGATCGATCGATCG")
        gc = gc_content(seq)
        # 应该在 40-60% 之间
        assert 0.4 <= gc <= 0.6


class TestBaseFrequency:
    """测试碱基频率统计。"""

    def test_basic_frequency(self):
        """基本频率统计。"""
        seq = DNASequence("AATGC")
        freq = base_frequency(seq)

        assert freq.a_count == 2
        assert freq.t_count == 1
        assert freq.g_count == 1
        assert freq.c_count == 1
        assert freq.total == 5

    def test_ratio_calculation(self):
        """比例计算。"""
        seq = DNASequence("AAAA")
        freq = base_frequency(seq)

        assert freq.a_ratio == 1.0
        assert freq.t_ratio == 0.0

    def test_empty_sequence(self):
        """空序列的频率。"""
        seq = DNASequence("")
        freq = base_frequency(seq)

        assert freq.total == 0
        assert freq.a_ratio == 0.0  # 不应该除零错误


class TestFindMotif:
    """测试 motif 搜索。"""

    def test_find_single_occurrence(self):
        """找到单个匹配。"""
        seq = DNASequence("ATGCATGC")
        positions = find_motif(seq, "ATG")
        assert positions == [0, 4]

    def test_overlapping_matches(self):
        """重叠匹配。"""
        seq = DNASequence("AAAA")
        positions = find_motif(seq, "AA")
        assert positions == [0, 1, 2]

    def test_no_match(self):
        """无匹配。"""
        seq = DNASequence("AAAA")
        positions = find_motif(seq, "GGG")
        assert positions == []

    def test_case_insensitive(self):
        """不区分大小写。"""
        seq = DNASequence("ATGC")
        positions = find_motif(seq, "atg")
        assert positions == [0]

    def test_invalid_motif(self):
        """无效 motif 应该报错。"""
        seq = DNASequence("ATGC")
        with pytest.raises(ValueError):
            find_motif(seq, "ATX")


class TestTranscribe:
    """测试转录功能。"""

    def test_transcribe(self):
        """基本转录。"""
        seq = DNASequence("ATGC")
        rna = transcribe(seq)
        assert rna == "AUGC"

    def test_no_thymine(self):
        """没有 T 的序列转录后不变。"""
        seq = DNASequence("GGCC")
        rna = transcribe(seq)
        assert rna == "GGCC"


class TestMolecularWeight:
    """测试分子量计算。"""

    def test_molecular_weight(self):
        """基本分子量计算。"""
        seq = DNASequence("A")
        mw = calculate_molecular_weight(seq)
        assert mw == pytest.approx(331.2, rel=0.01)

    def test_longer_sequence(self):
        """较长序列的分子量。"""
        seq = DNASequence("ATGC")
        mw = calculate_molecular_weight(seq)
        # 331.2 + 322.2 + 347.2 + 307.2 = 1307.8
        assert mw == pytest.approx(1307.8, rel=0.01)
```

### 运行测试

```bash
# 进入项目目录
cd /path/to/genomeflow

# 运行所有测试
uv run pytest

# 运行测试并显示覆盖率
uv run pytest --cov=genomeflow --cov-report=term-missing

# 只运行特定测试
uv run pytest tests/test_sequence.py -v
```

### 手动测试

创建一个测试用的 FASTA 文件 `test.fasta`：

```
>seq1 Human BRCA1 gene fragment
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
>seq2 E. coli gene fragment
GGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCC
```

然后运行：

```bash
# 安装项目（可编辑模式）
uv pip install -e .

# 快速分析单条序列
genomeflow quick ATGCGATCGATCG

# 分析 FASTA 文件
genomeflow analyze -f test.fasta

# 搜索特定 motif
genomeflow analyze -f test.fasta -m GATC

# 获取互补链
genomeflow complement ATGC

# 获取反向互补链
genomeflow complement ATGC -r
```

预期输出示例：

```
$ genomeflow quick ATGCGATCGATCG

序列: ATGCGATCGATCG
长度: 13 bp
GC 含量: 53.85%
碱基组成: A=3 (23.1%) T=3 (23.1%) G=4 (30.8%) C=3 (23.1%)
RNA 转录: AUGCGAUCGAUCG
反向互补: CGATCGATCGCAT
```

---

## 🤔 深入思考

### 常见问题解答

**Q: 为什么选择 uv 而不是 pip 或 poetry？**

A: uv 是 2024 年兴起的现代工具，有几个优势：
- 速度极快（Rust 编写）
- 自动管理虚拟环境
- 依赖解析更智能
- 一个工具替代多个（pip、virtualenv、pip-tools）

**Q: 为什么用 `src/` 布局？**

A: src 布局是 Python 社区推荐的最佳实践：
- 避免意外导入未安装的本地代码
- 测试运行的是安装后的包，更接近真实环境
- 更清晰的项目结构

**Q: 什么时候应该加类型注解？**

A: 本教程在以下场景使用类型注解：
- 函数签名（参数和返回值）
- 类属性
- 复杂的数据结构

不需要到处都加：局部变量如果类型显而易见就不用加。

**Q: 为什么 `read_fasta` 返回生成器而不是列表？**

A: 生物信息学文件可能非常大（GB 级别）。生成器允许：
- 边读边处理，内存占用恒定
- 可以提前终止（只要前 10 条）
- 支持管道式处理

### 设计模式分析

本项目用到的几个关键模式：

1. **数据类（Data Class）**：`BaseFrequency`, `FastaRecord`
   - 用于纯数据容器
   - 自动生成常用方法

2. **构建时验证**：`DNASequence.__init__`
   - 确保对象始终处于有效状态
   - 后续操作不需要重复验证

3. **生成器模式**：`read_fasta`
   - 惰性求值，内存高效
   - 适合处理大数据

4. **命令模式**：CLI 子命令
   - 每个命令一个处理函数
   - 易于扩展新命令

---

## 📝 总结与展望

### 你学到了什么

通过这个教程，你掌握了：

1. **生物信息学基础**
   - DNA 的结构和碱基配对规则
   - GC 含量的意义
   - FASTA 文件格式

2. **现代 Python 工程实践**
   - 使用 uv 管理项目和依赖
   - src 布局的项目结构
   - 类型注解的合理使用
   - dataclass 简化数据类定义

3. **软件设计原则**
   - 构建时验证确保数据有效性
   - 分离数据和行为（DNASequence vs analyzer）
   - 生成器处理大文件
   - 模块化组织代码

4. **测试驱动开发**
   - 使用 pytest 编写单元测试
   - 测试边界情况（空序列、无效输入）

### 下一步方向

如果你想继续深入，可以尝试：

1. **添加更多分析功能**
   - 序列比对（简单的点阵图）
   - 开放阅读框（ORF）查找
   - 密码子使用统计

2. **支持更多格式**
   - FASTQ（带质量分数的序列）
   - GenBank 格式

3. **性能优化**
   - 使用 NumPy 加速统计计算
   - 多线程处理大文件

4. **学习真正的生物信息学工具**
   - Biopython：Python 生物信息学库
   - Bioconda：生物信息学软件包管理

---

## 📚 参考资源

- [uv 官方文档](https://github.com/astral-sh/uv)
- [Python 类型注解指南](https://docs.python.org/zh-cn/3/library/typing.html)
- [FASTA 格式规范](https://en.wikipedia.org/wiki/FASTA_format)
- [Biopython 教程](https://biopython.org/wiki/Documentation)

---

**恭喜你完成了这个教程！** 你现在已经具备了编写生物信息学工具的基础能力。记住，最好的学习方式是动手实践——试着扩展这个项目，添加你感兴趣的功能吧！
