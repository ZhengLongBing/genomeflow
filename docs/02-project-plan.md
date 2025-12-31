# 项目规划和架构设计

本文档描述 GenomeFlow 项目的整体架构、模块划分和开发路线图。

---

## 目录

- [项目愿景](#项目愿景)
- [架构设计](#架构设计)
- [模块划分](#模块划分)
- [API 设计](#api-设计)
- [开发路线图](#开发路线图)

---

## 项目愿景

### 目标用户

- 生物信息学入门学习者
- 需要快速分析 DNA 序列的研究人员
- 希望学习 Python 工程实践的开发者

### 设计原则

1. **简单易用**：API 直观，无需深厚生物学背景
2. **渐进学习**：代码结构清晰，适合学习和扩展
3. **工程规范**：采用现代 Python 最佳实践
4. **实用优先**：覆盖最常用的分析需求

### 核心功能

| 功能 | 描述 | 优先级 | 状态 |
|------|------|--------|------|
| DNA 序列操作 | 互补链、反向互补、转录 | P0 | ✅ |
| 基础分析 | GC 含量、碱基频率 | P0 | ✅ |
| 文件读写 | FASTA 格式支持 | P0 | ✅ |
| 命令行工具 | 快速分析序列 | P0 | ✅ |
| RNA 序列支持 | 转录、翻译、ORF 查找 | P0 | ✅ |
| 蛋白质序列支持 | 分子量、等电点、疏水性 | P0 | ✅ |
| 可视化 | GC 分布、组成图、疏水性图 | P1 | ✅ |
| Motif 搜索 | 查找序列模式 | P1 | ✅ |
| 分子量计算 | DNA/蛋白质质量估算 | P1 | ✅ |
| ORF 查找 | 开放阅读框识别 | P2 | ✅ |
| 序列比对 | 简单的序列相似度比较 | P2 | 计划中 |

---

## 架构设计

### 分层架构

```
┌─────────────────────────────────────────────────────────┐
│                    CLI 层 (cli.py)                      │
│            命令行参数解析、用户交互、结果展示              │
└─────────────────────────┬───────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│                  业务逻辑层 (analyzer.py)                │
│         GC含量、碱基频率、Motif搜索、分子量计算            │
└─────────────────────────┬───────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│                  核心数据层 (sequence.py)                │
│            DNASequence 类、验证、基本操作                 │
└─────────────────────────┬───────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│                    IO 层 (io.py)                        │
│                FASTA/FASTQ 文件读写                      │
└─────────────────────────────────────────────────────────┘
```

### 设计决策

**为什么分层？**

| 层 | 职责 | 好处 |
|---|------|------|
| CLI | 用户交互 | 可替换为 GUI 或 Web |
| 业务逻辑 | 分析算法 | 独立测试，易于扩展 |
| 核心数据 | 数据表示 | 统一验证，类型安全 |
| IO | 文件处理 | 支持多种格式 |

**数据流向**

```
用户输入 → CLI 解析 → 读取文件/创建序列 → 分析处理 → 格式化输出
```

### 错误处理策略

| 错误类型 | 处理方式 |
|----------|----------|
| 无效碱基 | 抛出 `InvalidSequenceError` |
| 文件不存在 | 抛出 `FileNotFoundError` |
| 空序列 | 返回合理默认值（如 GC=0） |
| 参数错误 | CLI 层显示帮助信息 |

---

## 模块划分

### 目录结构

```
genomeflow/
├── src/
│   └── genomeflow/
│       ├── __init__.py      # 包入口，导出公共 API
│       ├── base.py          # 抽象基类：BaseSequence
│       ├── sequence.py      # DNA：DNASequence 类
│       ├── rna.py           # RNA：RNASequence 类、翻译功能
│       ├── protein.py       # 蛋白质：ProteinSequence 类
│       ├── analyzer.py      # 分析：各种分析函数
│       ├── io.py            # 读写：FASTA 文件处理
│       ├── cli.py           # 命令行：用户交互
│       └── visualization/   # 可视化模块
│           ├── __init__.py
│           ├── base.py      # 基础样式配置
│           ├── composition.py  # 组成分析图
│           └── profile.py   # 分布图
├── tests/
│   ├── __init__.py
│   ├── test_sequence.py    # DNA 序列测试
│   ├── test_rna.py         # RNA 序列测试
│   ├── test_protein.py     # 蛋白质序列测试
│   ├── test_analyzer.py    # 分析函数测试
│   ├── test_io.py          # IO 函数测试
│   └── test_visualization.py  # 可视化测试
├── docs/
│   ├── 01-concepts.md      # 核心概念
│   ├── 02-project-plan.md  # 项目规划（本文档）
│   ├── 03-tech-stack.md    # 技术选型
│   └── tutorials/          # 分步教程
│       ├── README.md
│       ├── tutorial-01-dna-analyzer.md
│       ├── tutorial-02-rna-sequence.md
│       ├── tutorial-03-protein-sequence.md
│       └── tutorial-04-visualization.md
├── pyproject.toml          # 项目配置
└── README.md               # 项目说明
```

### 模块职责

#### base.py - 抽象基类

```python
class BaseSequence(ABC):
    """核酸序列的抽象基类"""
    - 共享的序列操作（长度、索引、迭代）
    - 抽象方法：VALID_BASES, complement()
    - GC 含量计算

class InvalidSequenceError(ValueError):
    """序列验证失败时抛出"""
```

#### sequence.py - DNA 序列

```python
class DNASequence(BaseSequence):
    """DNA 序列"""
    - VALID_BASES = ATGC
    - 互补链、反向互补
    - transcribe() → RNASequence
```

#### rna.py - RNA 序列

```python
class RNASequence(BaseSequence):
    """RNA 序列"""
    - VALID_BASES = AUGC
    - translate() → 蛋白质序列
    - find_orfs() → ORF 列表
    - reverse_transcribe() → DNASequence
```

#### protein.py - 蛋白质序列

```python
class ProteinSequence:
    """蛋白质序列"""
    - 20 种标准氨基酸
    - molecular_weight() → 分子量
    - isoelectric_point() → 等电点
    - hydrophobicity_profile() → 疏水性分布
    - amino_acid_composition() → 氨基酸组成
```

#### analyzer.py - 分析功能

```python
# 统计类
@dataclass
class BaseFrequency:
    """碱基频率统计结果"""

# 分析函数
def gc_content(seq: DNASequence) -> float: ...
def base_frequency(seq: DNASequence) -> BaseFrequency: ...
def find_motif(seq: DNASequence, motif: str) -> list[int]: ...
def transcribe(seq: DNASequence) -> str: ...
def calculate_molecular_weight(seq: DNASequence) -> float: ...
```

#### io.py - 文件读写

```python
@dataclass
class FastaRecord:
    """FASTA 记录"""
    id: str
    description: str
    sequence: DNASequence

def read_fasta(file_path) -> Iterator[FastaRecord]: ...
def write_fasta(records, file_path) -> None: ...
```

#### visualization/ - 可视化模块

```python
# 碱基/氨基酸组成图
def plot_base_composition(seq, plot_type="bar"): ...
def plot_amino_acid_composition(protein): ...

# 分布图
def plot_gc_content(seq, window_size=100): ...
def plot_hydrophobicity(protein, window_size=9): ...

# 综合报告
def plot_sequence_overview(dna_sequence): ...
```

#### cli.py - 命令行接口

```bash
genomeflow analyze -f sequences.fasta [-m MOTIF]
genomeflow quick SEQUENCE
genomeflow complement SEQUENCE [-r]
```

---

## API 设计

### 公共 API

通过 `__init__.py` 导出的接口：

```python
from genomeflow import (
    # 核心类
    DNASequence,
    InvalidSequenceError,
    BaseFrequency,
    FastaRecord,

    # 分析函数
    gc_content,
    base_frequency,
    find_motif,
    transcribe,
    calculate_molecular_weight,

    # IO 函数
    read_fasta,
    write_fasta,
)
```

### 使用示例

```python
from genomeflow import DNASequence, gc_content, read_fasta

# 方式 1：直接创建序列
seq = DNASequence("ATGCGATCGATCG")
print(f"GC 含量: {gc_content(seq):.2%}")
print(f"反向互补: {seq.reverse_complement()}")

# 方式 2：从文件读取
for record in read_fasta("sequences.fasta"):
    print(f"{record.id}: {len(record.sequence)} bp")
```

### 设计原则

1. **类型安全**：函数签名使用类型注解
2. **不可变性**：`BaseFrequency` 等数据类是 frozen
3. **惰性求值**：`read_fasta` 返回生成器
4. **显式优于隐式**：不使用魔法，行为可预测

---

## 开发路线图

### 里程碑 1：基础功能（MVP）

**目标**：实现最小可用版本

- [x] 项目初始化（uv + pyproject.toml）
- [ ] DNASequence 类实现
- [ ] 基础分析函数（GC、碱基频率）
- [ ] FASTA 文件读写
- [ ] 基础 CLI
- [ ] 单元测试

### 里程碑 2：功能增强

**目标**：覆盖常用分析需求

- [ ] Motif 搜索功能
- [ ] 分子量计算
- [ ] 转录功能
- [ ] CLI 增强（更多选项）
- [ ] 测试覆盖率 > 90%

### 里程碑 3：高级功能

**目标**：支持复杂分析场景

- [ ] ORF（开放阅读框）查找
- [ ] 密码子使用统计
- [ ] 简单序列比对
- [ ] FASTQ 格式支持
- [ ] 性能优化（大文件处理）

### 里程碑 4：生态完善

**目标**：便于分享和使用

- [ ] 完整文档
- [ ] PyPI 发布
- [ ] GitHub Actions CI
- [ ] 示例数据和教程

---

## 扩展考虑

### 已完成的扩展

1. **RNA 序列支持** ✅
   - RNASequence 类，继承自 BaseSequence
   - 支持 U（尿嘧啶）
   - 翻译为蛋白质、ORF 查找

2. **蛋白质序列支持** ✅
   - ProteinSequence 类
   - 20 种氨基酸
   - 理化性质计算（分子量、等电点、疏水性）

3. **可视化** ✅
   - 碱基/氨基酸组成图
   - GC 含量分布图
   - 疏水性分布图

### 未来可能的方向

1. **序列比对**
   - 点阵图
   - 简单的全局/局部比对

2. **更多文件格式**
   - FASTQ（带质量分数）
   - GenBank 格式

3. **与 Biopython 集成**
   - 互操作接口
   - 格式转换

4. **Web 界面**
   - 使用 Streamlit 或 Gradio
   - 交互式分析

### 不会做的事情

- 复杂的序列比对算法（使用 BLAST）
- 基因组组装（使用专业工具）
- 结构预测（使用 AlphaFold）

---

下一步：阅读 [03-tech-stack.md](03-tech-stack.md) 了解技术选型。
