# 运筹学学习笔记 · Lecture Notes on Operations Research

[English](README.md) | **简体中文** | [Read online](https://betterweizh.github.io/lecture-notes-operation-research/en/) | [在线阅读](https://betterweizh.github.io/lecture-notes-operation-research/zh)

本项目汇总了个人运筹学学习笔记，内容包括线性规划、对偶理论与灵敏度分析、网络优化、动态规划、非线性规划和约束线性化，并以线性代数、矩阵微积分和凸性作为预备知识。每种方法都从基本原理推导，并配有约一百道完整的例题；书中另收录 160 余道例题和习题，附可折叠的完整解答；所有数值均已用代码重新计算；最后一部分用 Python（NumPy、SciPy/HiGHS、PuLP、NetworkX、CVXPY）实现书中的算法。

全书使用 MyST Markdown 编写，由 Sphinx 构建，通过 GitHub Pages 发布，提供英文版和简体中文版。

## 语言切换


| 版本     | 源文件     | 构建结果          |
| -------- | ---------- | ----------------- |
| English  | `docs/en/` | `_build/html/en/` |
| 简体中文 | `docs/zh/` | `_build/html/zh/` |

点击任意页面右上角的翻译图标，即可打开另一种语言中的同一章节；根页面 `_build/html/index.html` 会记住读者的选择。两个版本的文件名、标签、公式、代码和图片名完全相同；图片按版本分别绘制，标注使用各自的语言。`scripts/glossary.py` 是中英文术语的唯一来源，并为两个版本生成 `glossary.md`。

## 目录


| 部分                     | 章节                                                                                                                                                                                                   |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **第一部分　预备知识**   | 1.1 线性代数基础 · 1.2 向量与矩阵微积分 · 1.3 运筹学导论                                                                                                                                             |
| **第二部分　凸优化**     | 2.1 凸集 · 2.2 凸函数 · 2.3 凸优化与拉格朗日对偶                                                                                                                                                     |
| **第三部分　线性规划**   | 3.1 线性规划建模 · 3.2 线性规划的几何 · 3.3 单纯形法 · 3.4 人工变量：两阶段法与大 M 法 · 3.5 对偶理论 · 3.6 修正单纯形法与对偶单纯形法 · 3.7 灵敏度分析                                          |
| **第四部分　网络优化**   | 4.1 网络与最小生成树 · 4.2 最短路问题 · 4.3 最大流 · 4.4 最小费用流                                                                                                                                 |
| **第五部分　动态规划**   | 5.1 动态规划基础 · 5.2 设备更新 · 5.3 资源分配与资本预算 · 5.4 最短路的动态规划方法                                                                                                                 |
| **第六部分　非线性规划** | 6.1 非线性规划模型 · 6.2 最优性条件 · 6.3 KKT 条件 · 6.4 线搜索方法 · 6.5 无约束优化算法 · 6.6 约束优化算法                                                                                       |
| **第七部分　线性化**     | 7.1 逻辑约束与大 M 法 · 7.2 线性化技术 · 7.3 分段线性近似                                                                                                                                            |
| **第八部分　例题与习题** | 8.1 线性规划建模与几何 · 8.2 单纯形法 · 8.3 对偶 · 8.4 修正与对偶单纯形法、灵敏度分析 · 8.5 网络优化 · 8.6 动态规划基础与最短路 · 8.7 动态规划应用 · 8.8 最优性条件与 KKT · 8.9 非线性规划算法 |
| **第九部分　实现**       | 9.1 软件环境与求解器 · 9.2 用 Python 建立 LP 与 MILP 模型 · 9.3 线性规划算法的 Python 实现 · 9.4 网络算法的 Python 实现 · 9.5 动态规划的 Python 实现 · 9.6 非线性规划的 Python 实现               |
| **附录**                 | 符号表 · 参考文献                                                                                                                                                                                     |

## 仓库结构

```text
docs/
  _shared/            各语言版本共用的配置、扩展、模板与样式
    conf_common.py    Sphinx 设置、MathJax 宏、主题
    _ext/             按部分编号（带部分号前缀）、语言切换
    figures/          现成（非生成）图片及其来源说明 SOURCES.md
    locales/          证明类环境的中文名称（定义、定理、例……）
  en/                 英文版：conf.py、index.md、PP-CC-*.md 各章（部分号-章号）、notation.md、glossary.md、references.md
  zh/                 简体中文版：文件与英文版相同，内容为中文
scripts/
  build_docs.py       构建网站（先生成图片）
  make_figures.py     绘制生成图片
  figures/            每章一个绘图模块（pPPcCC.py）及 common.py
  run_code_cells.py   运行代码章节中的 Python 代码块
  glossary.py         中英文术语表，为两个版本生成 glossary.md
.github/workflows/    构建并发布到 GitHub Pages
requirements.txt      Sphinx 工具链及代码章节所需的软件包
```

`sources/`（笔记所依据的课程资料）只保存在本地，不纳入仓库（已列入 `.gitignore`）。

## 图片

所有图片均由代码生成：每张图由 `scripts/figures/pPPcCC.py`（`PP`、`CC` 为首次出现该图的章节所在的部分号和章号）中一个带注释的 Python 函数绘制，名称以部分号和章号开头（如 `03-02-product-mix-region`），构建网站时写入 `docs/en/figures/`。生成的文件不纳入 git，请勿手工修改。现成（非生成）图片应放在 `docs/_shared/figures/`，并在其 `SOURCES.md` 中注明来源。绘图辅助函数见 [`scripts/figures/README.md`](scripts/figures/README.md)。

## 本地构建

```bash
python -m venv .venv
source .venv/bin/activate            # Windows：.venv\Scripts\activate
pip install -r requirements.txt

python scripts/build_docs.py --serve        # 构建并打开 http://127.0.0.1:8000/
python scripts/build_docs.py --strict       # 警告视为错误（与 GitHub 上一致）
python scripts/build_docs.py --clean        # 从头重新构建
python scripts/build_docs.py --figures-only # 只绘制图片
python scripts/build_docs.py --lang zh       # 只构建中文版
python scripts/run_code_cells.py docs/en/09-*.md   # 运行全部代码块
```

网站输出到 `_build/html/`。缺失或比绘图模块旧的图片会自动重新绘制。

## 编辑流程

1. 同时修改 `docs/en/` 中的章节和 `docs/zh/` 中的同名文件（若图片有变化，同时修改 `scripts/figures/` 中对应的绘图模块；含文字的标注有英文和中文两个版本）。保持标签（`eq-03-02-…`、`thm-03-05-…`、`fig-04-03-…`）不变，且两个版本一致，其他章节会引用它们。
2. 新增章节需要同时加入 `docs/en/index.md` 和 `docs/zh/index.md` 中相应的 toctree；新术语加入 `scripts/glossary.py` 后运行该脚本。
3. 运行 `python scripts/build_docs.py --strict`；修改代码章节时，再运行 `scripts/run_code_cells.py`。
4. 推送到 `main` 分支后，GitHub Actions 会构建网站并发布到 GitHub Pages（首次需在仓库设置中把 Pages 的来源设为 "GitHub Actions"）。

## 写作约定

- 向量用粗体小写（$\mathbf{x}$），矩阵用粗体大写（$\mathbf{A}$），统一使用 `docs/_shared/conf_common.py` 中定义的 MathJax 宏（`\vx`、`\mA`、`\T` 等）；符号表列出了全部符号。
- 单纯形表的目标行放在最后，记为 $-z$，其中为检验数 $\bar c_j$；对偶变量记为 $\boldsymbol{\pi}$。
- 各章采用统一结构：概述、学习目标、正文各节、小结与常见错误、参考文献。例题使用 `{prf:example}`，其后为默认折叠的 `{dropdown} Solution`；证明同样默认折叠。
- 章节在各部分内编号（3.2 即第三部分第 2 章），图、表、例题和公式带同样的前缀（3.2.1）。
- 标签带章节号；跨章引用使用 `{doc}` 和标签，不使用标题锚点。
- 参考文献采用 APA 第 7 版格式，每章列出本章引用的文献。

## 作者与许可

Better W.Z.，以 [MIT 许可证](LICENSE) 发布。笔记所依据的课程资料及所有被引用的作品，其权利归原作者所有。

## 引用方式

Better W.Z. (2026). *Lecture Notes on Operations Research*. GitHub.
https://github.com/betterweizh/lecture-notes-operation-research

```bibtex
@misc{betterwz2026lecture,
  author       = {Better W.Z.},
  title        = {Lecture Notes on Operations Research},
  year         = {2026},
  howpublished = {\url{https://github.com/betterweizh/lecture-notes-operation-research}},
  note         = {Lecture notes}
}
```
