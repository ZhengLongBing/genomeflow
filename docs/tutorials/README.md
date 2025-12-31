# GenomeFlow 教程索引

欢迎来到 GenomeFlow 教程系列！本系列教程将带你从零开始构建一个完整的序列分析工具，支持 DNA、RNA 和蛋白质。

---

## 学习路径

```
开始
  │
  ▼
┌─────────────────────────────────────┐
│  阅读背景知识（可选但推荐）            │
│  ├── 01-concepts.md   核心概念       │
│  ├── 02-project-plan.md 项目规划     │
│  └── 03-tech-stack.md 技术选型       │
└─────────────────────────────────────┘
  │
  ▼
┌─────────────────────────────────────┐
│  教程 01：DNA 序列分析器              │
│  DNASequence 类 + 基础分析 + CLI     │
│  预计时间：2-3 小时                   │
└─────────────────────────────────────┘
  │
  ▼
┌─────────────────────────────────────┐
│  教程 02：RNA 序列支持                │
│  RNASequence + 转录 + 翻译 + ORF     │
│  预计时间：1.5-2 小时                 │
└─────────────────────────────────────┘
  │
  ▼
┌─────────────────────────────────────┐
│  教程 03：蛋白质序列支持              │
│  ProteinSequence + 理化性质分析       │
│  预计时间：1.5-2 小时                 │
└─────────────────────────────────────┘
  │
  ▼
┌─────────────────────────────────────┐
│  教程 04：序列可视化                  │
│  matplotlib + 组成图 + 分布图        │
│  预计时间：1-1.5 小时                 │
└─────────────────────────────────────┘
  │
  ▼
┌─────────────────────────────────────┐
│  完成！你已掌握：                     │
│  ✅ 中心法则完整工作流                │
│  ✅ 现代 Python 工程实践              │
│  ✅ 序列分析与可视化                  │
└─────────────────────────────────────┘
```

---

## 教程列表

### 教程 01：DNA 序列分析器基础

**文件**：[tutorial-01-dna-analyzer.md](tutorial-01-dna-analyzer.md)

**你将学到**：
- 使用 uv 管理 Python 项目
- 实现 DNASequence 类
- 计算 GC 含量和碱基频率
- 生成互补链和反向互补链
- 实现 FASTA 文件读写
- 编写单元测试
- 构建命令行接口

**前置知识**：
- 基础 Python 语法
- 命令行基本操作
- 不需要生物学背景

**产出**：一个可用的 DNA 分析 CLI 工具

---

### 教程 02：RNA 序列支持

**文件**：[tutorial-02-rna-sequence.md](tutorial-02-rna-sequence.md)

**你将学到**：
- 理解 RNA 与 DNA 的区别
- 设计可复用的抽象基类
- 实现 RNASequence 类
- DNA 转录为 RNA
- RNA 翻译为蛋白质（密码子表）
- 查找开放阅读框（ORF）

**前置知识**：
- 完成教程 01
- 理解 Python 类继承

**产出**：支持 DNA→RNA→蛋白质的完整转换

---

### 教程 03：蛋白质序列支持

**文件**：[tutorial-03-protein-sequence.md](tutorial-03-protein-sequence.md)

**你将学到**：
- 20 种标准氨基酸的性质
- 实现 ProteinSequence 类
- 计算蛋白质分子量
- 计算等电点（pI）
- 分析疏水性分布
- 计算消光系数

**前置知识**：
- 完成教程 01 和 02
- 理解氨基酸的基本概念

**产出**：完整的"中心法则"工作流

---

### 教程 04：序列可视化

**文件**：[tutorial-04-visualization.md](tutorial-04-visualization.md)

**你将学到**：
- 使用 matplotlib 进行科学可视化
- 绘制碱基/氨基酸组成图（饼图、条形图）
- 绘制 GC 含量滑动窗口分布图
- 绘制疏水性分布图（跨膜区预测）
- 创建综合分析报告

**前置知识**：
- 完成教程 01-03
- 基本的 matplotlib 使用经验（可选）

**产出**：专业的序列分析图表

---

## 如何使用这些教程

### 推荐方式：边学边做

1. **阅读一个小节**
2. **动手敲代码**（不要复制粘贴）
3. **运行测试验证**
4. **遇到问题先思考**
5. **继续下一节**

### 代码仓库

每个教程结束后，你的代码应该与教程中展示的一致。如果遇到困难，可以参考完整代码。

### 提问和反馈

如果教程中有不清楚的地方，欢迎：
- 在代码中添加注释记录疑问
- 尝试修改代码验证理解
- 查阅参考资源

---

## 前置准备

### 环境要求

- Python 3.12 或更高版本
- 终端/命令行
- 代码编辑器（推荐 VS Code 或 PyCharm）

### 安装 uv

```bash
# macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows (PowerShell)
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"

# 验证安装
uv --version
```

### 创建项目

```bash
# 创建新项目
uv init genomeflow
cd genomeflow

# 安装开发依赖
uv add --dev pytest pytest-cov

# 安装可视化依赖（教程04需要）
uv add matplotlib

# 验证环境
uv run python --version
```

---

## 教程概览表

| 教程 | 主题 | 新增模块 | 关键技能 |
|------|------|----------|----------|
| 01 | DNA 分析 | sequence.py, analyzer.py, io.py, cli.py | 类设计、测试、CLI |
| 02 | RNA 支持 | base.py, rna.py | 抽象基类、翻译 |
| 03 | 蛋白质支持 | protein.py | 理化性质计算 |
| 04 | 可视化 | visualization/ | matplotlib 绑图 |

---

## 参考资源

### 生物学背景

- [01-concepts.md](../01-concepts.md) - 本项目的核心概念
- [NCBI 教程](https://www.ncbi.nlm.nih.gov/guide/all/) - 官方生物信息学资源

### Python 技术

- [03-tech-stack.md](../03-tech-stack.md) - 本项目的技术选型
- [Real Python](https://realpython.com/) - 高质量 Python 教程
- [Python 官方文档](https://docs.python.org/zh-cn/3/)

### 工具文档

- [uv 文档](https://github.com/astral-sh/uv)
- [pytest 文档](https://docs.pytest.org/)
- [matplotlib 文档](https://matplotlib.org/stable/)
- [Ruff 文档](https://docs.astral.sh/ruff/)

---

## 学习建议

### 对于编程初学者

1. 先完成 Python 基础教程
2. 不要跳过测试部分
3. 遇到错误是好事，认真阅读错误信息

### 对于有经验的开发者

1. 重点关注生物学概念
2. 可以快速浏览代码，重点看设计决策
3. 尝试扩展功能作为练习

### 对于生物学背景的学习者

1. 代码部分可能需要多花时间
2. 利用你的领域知识理解"为什么"
3. 思考如何应用到自己的研究中

---

**准备好了吗？让我们从 [教程 01](tutorial-01-dna-analyzer.md) 开始吧！**
