# Contributing to GenomeFlow

感谢你对 GenomeFlow 项目的兴趣！我们欢迎各种形式的贡献。

## 如何贡献

### 报告 Bug

如果你发现了 bug，请在 GitHub Issues 中提交，并包含：

- 问题的清晰描述
- 复现步骤
- 期望行为与实际行为
- 你的 Python 版本和操作系统
- 相关的代码片段或错误信息

### 提交功能建议

我们欢迎新功能建议！请在 Issues 中描述：

- 功能的目的和用例
- 期望的 API 设计
- 可能的实现方式（可选）

### 提交代码

1. **Fork 仓库**

2. **克隆你的 fork**
   ```bash
   git clone https://github.com/your-username/GenomeFlow.git
   cd GenomeFlow
   ```

3. **创建新分支**
   ```bash
   git checkout -b feature/your-feature-name
   ```

4. **设置开发环境**
   ```bash
   # 使用 uv 安装依赖
   uv sync

   # 或使用 pip
   pip install -e ".[dev]"
   ```

5. **进行修改**

6. **运行测试**
   ```bash
   uv run pytest
   ```

7. **运行代码检查**
   ```bash
   uv run ruff check src tests
   uv run mypy src
   ```

8. **提交更改**
   ```bash
   git add .
   git commit -m "feat: 描述你的更改"
   ```

9. **推送并创建 Pull Request**
   ```bash
   git push origin feature/your-feature-name
   ```

## 代码规范

### Python 风格

- 遵循 PEP 8 规范
- 使用 ruff 进行代码格式化和检查
- 行长度限制为 88 字符
- 使用类型注解

### 提交信息规范

使用 [Conventional Commits](https://www.conventionalcommits.org/) 格式：

- `feat:` 新功能
- `fix:` Bug 修复
- `docs:` 文档更新
- `test:` 测试相关
- `refactor:` 代码重构
- `chore:` 构建/工具相关

示例：
```
feat: 添加蛋白质二级结构预测功能
fix: 修复 GC 含量计算的边界情况
docs: 更新 RNA 分析教程
```

### 测试要求

- 所有新功能必须包含测试
- 修复 bug 时请添加回归测试
- 保持测试覆盖率不低于现有水平

## 开发指南

### 项目结构

```
GenomeFlow/
├── src/genomeflow/    # 源代码
├── tests/             # 测试文件
├── docs/              # 文档
└── pyproject.toml     # 项目配置
```

### 添加新功能

1. 在 `src/genomeflow/` 中添加实现
2. 在 `__init__.py` 中导出公共 API
3. 在 `tests/` 中添加测试
4. 更新相关文档

## 行为准则

请保持友善和尊重。我们致力于为所有人提供一个友好、安全和欢迎的环境。

## 问题？

如有任何问题，欢迎在 Issues 中提问或讨论。
