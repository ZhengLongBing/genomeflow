# 核心概念和原理

本文档介绍 DNA 序列分析的基础知识，帮助你建立全局认知。即使没有生物学背景，也能理解我们要解决的问题。

---

## 目录

- [DNA 是什么](#dna-是什么)
- [为什么要用程序分析 DNA](#为什么要用程序分析-dna)
- [核心分析指标](#核心分析指标)
- [常用文件格式](#常用文件格式)
- [关键术语表](#关键术语表)

---

## DNA 是什么

### 生命的"源代码"

如果把生命比作一台计算机：
- **DNA** 是存储在硬盘里的程序源代码
- **RNA** 是从源代码编译出的中间代码
- **蛋白质** 是最终运行的程序

DNA（脱氧核糖核酸）是一种长链分子，由四种碱基组成：

| 碱基 | 符号 | 英文全名 | 配对碱基 |
|------|------|----------|----------|
| 腺嘌呤 | A | Adenine | T |
| 胸腺嘧啶 | T | Thymine | A |
| 鸟嘌呤 | G | Guanine | C |
| 胞嘧啶 | C | Cytosine | G |

### 双螺旋结构

DNA 由两条链组成，像拉链一样配对：

```
5' ----A---T---G---C---A---T---- 3'
       |   |   |   |   |   |
3' ----T---A---C---G---T---A---- 5'
```

关键规则：
- **A 永远与 T 配对**（2 个氢键）
- **G 永远与 C 配对**（3 个氢键）
- 两条链**方向相反**（5'→3' 和 3'→5'）

<details>
<summary>什么是 5' 和 3' 端？</summary>

这是化学结构的命名。DNA 链有方向性，就像单行道：
- 5' 端是起点（带磷酸基团）
- 3' 端是终点（带羟基）

细胞里的酶（如 DNA 聚合酶）只能沿 5'→3' 方向工作，这就是为什么方向很重要。

</details>

### 从 DNA 到蛋白质

```
DNA  →  RNA  →  蛋白质
    转录     翻译
```

1. **转录（Transcription）**：DNA → RNA
   - 将 T 替换为 U（尿嘧啶）
   - 产生信使 RNA（mRNA）

2. **翻译（Translation）**：RNA → 蛋白质
   - 每 3 个碱基（密码子）对应 1 个氨基酸
   - 氨基酸链折叠成蛋白质

---

## 为什么要用程序分析 DNA

### 数据规模

| 生物 | 基因组大小 |
|------|-----------|
| 大肠杆菌 | 460 万碱基 |
| 人类 | 32 亿碱基 |
| 小麦 | 170 亿碱基 |

手工分析这么多数据？不可能。

### 实际应用场景

**1. 医学诊断**
- 检测致病基因突变
- 癌症基因分型
- 药物敏感性预测

**2. 科学研究**
- 新物种基因组测序
- 进化关系分析
- 基因功能研究

**3. 农业育种**
- 作物性状关联分析
- 分子标记辅助选择

**4. 法医鉴定**
- DNA 指纹比对
- 亲子鉴定

---

## 核心分析指标

### 1. GC 含量（GC Content）

**定义**：序列中 G 和 C 碱基的百分比

```
GC% = (G + C) / 总长度 × 100%
```

**为什么重要？**

| 特性 | 说明 |
|------|------|
| 稳定性 | GC 含量高的 DNA 更稳定（G-C 有 3 个氢键） |
| 物种特征 | 不同物种有特征性 GC 含量 |
| 基因预测 | 编码区通常 GC 含量较高 |
| 实验设计 | PCR 引物设计需考虑 GC 含量 |

**典型值参考**：
- 大肠杆菌：~51%
- 人类：~41%
- 疟原虫：~19%（极端低 GC）

### 2. 碱基频率（Base Frequency）

统计每种碱基的数量和比例。

**查加夫法则**（Chargaff's Rules）：
- 在双链 DNA 中：A% ≈ T%，G% ≈ C%
- 这是碱基配对的必然结果

如果统计结果严重偏离这个规律，可能：
- 序列是单链的
- 数据有质量问题
- 序列来自特殊生物

### 3. 互补链与反向互补链

| 概念 | 操作 | 示例 |
|------|------|------|
| 互补链 | A↔T, G↔C | ATGC → TACG |
| 反向互补 | 互补 + 反转 | ATGC → GCAT |

**为什么需要反向互补？**

因为 DNA 双链方向相反：
```
原始链：   5'-ATGC-3'
互补链：   3'-TACG-5'  (方向相反)
反向互补： 5'-GCAT-3'  (方向一致，便于比对)
```

在序列比对、引物设计中，经常需要用到反向互补。

### 4. 序列模式（Motif）

Motif 是具有生物学意义的短序列模式。

**常见类型**：
- **转录因子结合位点**：基因调控区域
- **限制性酶切位点**：EcoRI 识别 GAATTC
- **起始密码子**：ATG（翻译开始）
- **终止密码子**：TAA、TAG、TGA

### 5. 分子量（Molecular Weight）

DNA 片段的质量，单位是道尔顿（Da）或千道尔顿（kDa）。

**用途**：
- 凝胶电泳：迁移率与分子量相关
- 摩尔浓度计算
- 测序文库质控

---

## 常用文件格式

### FASTA 格式

最常用的序列存储格式，简单直观。

```
>seq_id description text
ATGCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGA
>another_seq Another description
GGCCAATTGGCCAATTGGCCAATTGGCCAATT
```

**格式规则**：
- `>` 开头的行是标题（header）
- 第一个空格前是序列 ID
- 后续行是序列内容
- 序列可以换行（通常每行 60-80 字符）
- 一个文件可包含多条序列

### FASTQ 格式

带质量分数的序列格式，用于测序原始数据。

```
@read_id
ATGCGATCGATCG
+
IIIIIIIIIIIII
```

**四行一组**：
1. `@` 开头的 ID 行
2. 序列
3. `+`（可选重复 ID）
4. 质量分数（ASCII 编码）

### GenBank 格式

注释丰富的格式，包含基因位置、功能等信息。

```
LOCUS       SCU49845     5028 bp    DNA     PLN       21-JUN-1999
DEFINITION  Saccharomyces cerevisiae TCP1-beta gene, partial cds
...
ORIGIN
        1 gatcctccat atacaacggt atctccacct caggtttaga tctcaacaac ggaaccattg
...
```

---

## 关键术语表

| 术语 | 英文 | 解释 |
|------|------|------|
| 碱基 | Base | A、T、G、C 四种核苷酸的简称 |
| 碱基对 | Base Pair (bp) | 两条链配对的碱基，也是长度单位 |
| 核苷酸 | Nucleotide | 碱基 + 糖 + 磷酸的完整单元 |
| 互补链 | Complement | 按配对规则对应的另一条链 |
| 正链/负链 | Sense/Antisense | 编码蛋白的链/模板链 |
| 密码子 | Codon | 3 个碱基为一组，编码氨基酸 |
| 启动子 | Promoter | 基因转录的起始信号区域 |
| ORF | Open Reading Frame | 开放阅读框，可能编码蛋白的区域 |
| PCR | Polymerase Chain Reaction | 聚合酶链式反应，扩增 DNA 片段 |
| 引物 | Primer | PCR 中指定扩增区域的短序列 |

---

## 延伸阅读

如果你想深入了解：

1. **入门书籍**
   - 《分子生物学》（Molecular Biology of the Cell）
   - 《生物信息学》（Bioinformatics）

2. **在线课程**
   - Coursera: Bioinformatics Specialization
   - edX: Data Analysis for Life Sciences

3. **工具和数据库**
   - NCBI（美国国家生物技术信息中心）
   - Ensembl（基因组数据库）
   - UniProt（蛋白质数据库）

---

下一步：阅读 [02-project-plan.md](02-project-plan.md) 了解项目架构设计。
