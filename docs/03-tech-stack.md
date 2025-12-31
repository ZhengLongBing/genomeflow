# 技术选型说明

本文档解释 GenomeFlow 项目的技术选择，以及为什么做出这些决策。

---

## 目录

- [Python 版本](#python-版本)
- [项目管理工具](#项目管理工具)
- [项目结构](#项目结构)
- [代码风格和工具](#代码风格和工具)
- [测试框架](#测试框架)
- [第三方库评估](#第三方库评估)

---

## Python 版本

### 选择：Python 3.12+

**原因**：

1. **类型系统增强**
   - 内置泛型语法：`list[int]` 而非 `List[int]`
   - 更好的类型推断

2. **性能提升**
   - 解释器优化，运行更快
   - 更低的内存占用

3. **新语法特性**
   - 改进的错误消息
   - 更好的 f-string 支持

**权衡**：

| 方面 | Python 3.12 | Python 3.9 |
|------|-------------|------------|
| 新特性 | 更多 | 较少 |
| 兼容性 | 部分库可能不支持 | 广泛支持 |
| 维护周期 | 到 2028 年 | 到 2025 年 |

对于新项目，选择较新版本是更好的长期投资。

---

## 项目管理工具

### 选择：uv

**什么是 uv？**

[uv](https://github.com/astral-sh/uv) 是 Astral 公司开发的 Python 包管理工具，用 Rust 编写。

**对比其他工具**：

| 工具 | 优点 | 缺点 |
|------|------|------|
| **pip** | 标准工具，无需安装 | 慢，依赖解析有限 |
| **poetry** | 功能完善，锁文件支持 | 较慢，配置复杂 |
| **pdm** | PEP 标准兼容 | 社区较小 |
| **uv** | 极快，一体化 | 较新，生态还在发展 |

**为什么选 uv？**

1. **速度**
   ```bash
   # 安装 100 个包
   pip:     45 秒
   poetry:  30 秒
   uv:      3 秒
   ```

2. **一体化**
   ```bash
   # uv 替代了多个工具
   pip        → uv pip
   pip-tools  → uv pip compile
   virtualenv → uv venv
   pyenv      → uv python
   ```

3. **现代体验**
   - 自动创建虚拟环境
   - 智能依赖解析
   - 锁文件支持

**基本用法**：

```bash
# 安装 uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# 初始化项目
uv init myproject

# 添加依赖
uv add requests

# 添加开发依赖
uv add --dev pytest

# 同步环境
uv sync

# 运行命令
uv run python script.py
uv run pytest
```

---

## 项目结构

### 选择：src 布局

```
项目根目录/
├── src/
│   └── genomeflow/    # 包代码
├── tests/             # 测试代码
└── pyproject.toml
```

**对比平铺布局**：

```
# 平铺布局（不推荐）
项目根目录/
├── genomeflow/        # 包代码直接在根目录
├── tests/
└── pyproject.toml
```

**为什么用 src 布局？**

| 问题 | 平铺布局 | src 布局 |
|------|----------|----------|
| 测试未安装的代码 | 可能发生 | 不可能 |
| 导入混淆 | 容易发生 | 清晰分离 |
| 打包问题 | 容易遗漏文件 | 边界清晰 |

**实际例子**：

```bash
# 平铺布局的问题
cd myproject
python -c "import genomeflow"  # 导入的是本地目录，不是安装的包！

# src 布局
cd myproject
python -c "import genomeflow"  # 必须先安装才能导入
```

---

## 代码风格和工具

### 类型注解策略

**原则**：在有价值的地方使用类型注解

| 场景 | 是否使用 | 原因 |
|------|----------|------|
| 函数签名 | 是 | API 文档化，IDE 支持 |
| 类属性 | 是 | 数据结构清晰 |
| 局部变量 | 视情况 | 显而易见时不加 |
| 返回复杂结构 | 是 | 帮助理解返回值 |

**示例**：

```python
# 好：函数签名清晰
def gc_content(seq: DNASequence) -> float:
    counts = Counter(seq)  # 局部变量，类型显而易见
    return (counts["G"] + counts["C"]) / len(seq)

# 好：复杂返回值
def analyze(seq: DNASequence) -> tuple[float, BaseFrequency]:
    ...

# 不必要：显而易见的局部变量
name: str = "ATGC"  # 多余，直接 name = "ATGC"
```

### 代码格式化

**工具**：[Ruff](https://github.com/astral-sh/ruff)

Ruff 是 Python 的超快代码检查和格式化工具。

```toml
# pyproject.toml
[tool.ruff]
line-length = 88
target-version = "py312"

[tool.ruff.lint]
select = ["E", "F", "I", "UP"]
```

**常用规则**：
- `E`: pycodestyle 错误
- `F`: pyflakes
- `I`: isort（导入排序）
- `UP`: pyupgrade（现代语法）

### 文档字符串

**风格**：Google 风格

```python
def gc_content(seq: DNASequence) -> float:
    """
    计算 GC 含量。

    Args:
        seq: DNA 序列

    Returns:
        GC 含量，范围 0.0-1.0

    Raises:
        ValueError: 如果序列为空

    Example:
        >>> gc_content(DNASequence("ATGC"))
        0.5
    """
```

---

## 测试框架

### 选择：pytest

**为什么不用 unittest？**

| 特性 | unittest | pytest |
|------|----------|--------|
| 语法 | 需要类和方法 | 简单函数 |
| 断言 | `self.assertEqual(a, b)` | `assert a == b` |
| 夹具 (fixture) | setUp/tearDown | 灵活的 fixture |
| 插件生态 | 有限 | 丰富 |
| 输出 | 简单 | 详细，美观 |

**pytest 示例**：

```python
# 简洁的测试函数
def test_gc_content():
    seq = DNASequence("ATGC")
    assert gc_content(seq) == 0.5

# 参数化测试
@pytest.mark.parametrize("sequence,expected", [
    ("ATGC", 0.5),
    ("GGCC", 1.0),
    ("AATT", 0.0),
])
def test_gc_content_various(sequence, expected):
    assert gc_content(DNASequence(sequence)) == expected

# 夹具
@pytest.fixture
def sample_sequence():
    return DNASequence("ATGCGATCGATCG")

def test_with_fixture(sample_sequence):
    assert len(sample_sequence) == 13
```

### 测试配置

```toml
# pyproject.toml
[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-v --tb=short"
```

---

## 第三方库评估

### 不使用 Biopython

**Biopython** 是 Python 生物信息学的标准库，功能非常强大。

**为什么不用？**

| 考虑因素 | 使用 Biopython | 自己实现 |
|----------|---------------|----------|
| 学习价值 | 低（用现成的） | 高（理解原理） |
| 依赖大小 | 大（~50MB） | 无依赖 |
| 灵活性 | 按库设计 | 完全控制 |
| 适用场景 | 生产项目 | 教学项目 |

**本项目定位**：教学和学习，因此选择从零实现核心功能。

**如果是生产项目**，建议直接使用 Biopython：

```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

for record in SeqIO.parse("sequences.fasta", "fasta"):
    print(f"{record.id}: GC={gc_fraction(record.seq):.2%}")
```

### CLI 框架评估

| 框架 | 优点 | 缺点 |
|------|------|------|
| argparse | 标准库，无依赖 | 语法繁琐，子命令麻烦 |
| **click** | 装饰器语法，功能丰富 | 外部依赖 |
| typer | 类型注解驱动，最现代 | 依赖 click，较新 |

**选择**：Click

理由：
- **装饰器语法直观**：代码即文档，参数定义清晰
- **功能丰富**：彩色输出、进度条、自动补全
- **测试友好**：内置 `CliRunner` 便于测试
- **生态成熟**：Flask 等知名项目采用，稳定可靠
- **子命令支持好**：`@group` + `@command` 简洁明了

**Click vs argparse 对比**：

```python
# argparse 写法
parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(dest="command")
quick_parser = subparsers.add_parser("quick")
quick_parser.add_argument("sequence")
args = parser.parse_args()
if args.command == "quick":
    do_quick(args.sequence)

# Click 写法
@click.group()
def cli(): pass

@cli.command()
@click.argument("sequence")
def quick(sequence):
    do_quick(sequence)
```

Click 的装饰器风格更 Pythonic，代码量更少，可读性更好。

**为什么不选 typer？**

Typer 基于 Click，使用类型注解定义参数，更加现代：

```python
import typer

def quick(sequence: str):
    """快速分析序列"""
    ...
```

虽然 Typer 更简洁，但：
- Click 更成熟稳定
- Click 文档和社区资源更丰富
- 对于入门学习，Click 的显式装饰器更容易理解

### 数据处理库

| 需求 | 选择 | 备注 |
|------|------|------|
| 碱基统计 | `collections.Counter` | 标准库足够 |
| 数据类 | `dataclasses` | 标准库 |
| 大规模数据 | 未来考虑 NumPy | 当前不需要 |

---

## 开发环境建议

### 推荐 IDE

1. **VS Code** + Python 扩展
   - 免费，轻量
   - 优秀的 Python 支持

2. **PyCharm**
   - 功能最全面
   - 社区版免费

### 推荐配置

```json
// VS Code settings.json
{
  "python.defaultInterpreterPath": ".venv/bin/python",
  "python.analysis.typeCheckingMode": "basic",
  "[python]": {
    "editor.formatOnSave": true,
    "editor.defaultFormatter": "charliermarsh.ruff"
  }
}
```

### Git 配置

```gitignore
# .gitignore
.venv/
__pycache__/
*.pyc
.pytest_cache/
*.egg-info/
dist/
build/
.coverage
```

---

## 总结

| 类别 | 选择 | 关键理由 |
|------|------|----------|
| Python 版本 | 3.12+ | 新特性，长期支持 |
| 包管理 | uv | 速度快，一体化 |
| 项目结构 | src 布局 | 避免导入问题 |
| 类型注解 | 适度使用 | 平衡清晰度和简洁性 |
| 格式化 | Ruff | 快速，功能全面 |
| 测试 | pytest | 简洁，生态丰富 |
| CLI | Click | 装饰器语法，功能丰富 |
| 可视化 | matplotlib | 科学绑图标准 |

这些选择代表了 2024-2025 年 Python 社区的最佳实践，适合新项目采用。

---

下一步：进入 [tutorials/](tutorials/README.md) 开始动手实践。
