# GenomeFlow

DNA 序列分析工具包 - 面向生物信息学初学者的 Python 库。

## 功能特性

- **DNA 序列分析**：GC 含量、碱基频率、motif 查找
- **RNA 序列支持**：转录、翻译、ORF 查找
- **蛋白质分析**：分子量、等电点、疏水性分析
- **文件 I/O**：FASTA 格式读写
- **可视化**：碱基组成图、疏水性分布图
- **命令行工具**：快速分析序列

## 安装

```bash
# 使用 uv
uv add genomeflow

# 或使用 pip
pip install genomeflow
```

## 快速开始

### Python API

```python
from genomeflow import DNASequence, gc_content, base_frequency

# 创建 DNA 序列
seq = DNASequence("ATGCGATCGATCG")

# 分析 GC 含量
print(f"GC 含量: {gc_content(seq):.2%}")

# 碱基频率
freq = base_frequency(seq)
print(f"A: {freq.a_count}, T: {freq.t_count}, G: {freq.g_count}, C: {freq.c_count}")

# 转录为 RNA
rna = seq.transcribe()
print(f"RNA: {rna}")

# 翻译为蛋白质
protein = rna.translate_to_protein()
print(f"蛋白质: {protein}")
```

### 命令行

```bash
# 快速分析序列
genomeflow quick ATGCGATCGATCG

# 分析 FASTA 文件
genomeflow analyze -f sequences.fasta

# 获取反向互补链
genomeflow complement ATGC --reverse

# 翻译为蛋白质
genomeflow translate ATGTTTTAA
```

## 开发

```bash
# 克隆仓库
git clone https://github.com/yourusername/genomeflow.git
cd genomeflow

# 安装开发依赖
uv sync --dev

# 运行测试
uv run pytest

# 运行类型检查
uv run mypy src

# 运行代码检查
uv run ruff check src
```

## 许可证

MIT License
