# 运筹学学习笔记

本学习笔记系统阐述运筹学的核心模型与方法：线性规划及其对偶理论与灵敏度分析、网络优化、动态规划、非线性规划，以及将逻辑关系与非线性关系转化为混合整数线性模型的线性化技术。每种方法均先陈述其假设，再从基本原理出发加以推导。所需的线性代数与矩阵微积分知识在预备知识部分复习；凸集、凸函数与拉格朗日对偶则单独成章。

部分方法配有附完整解答的例题与习题，全部数值均已通过代码重新计算核实。学习笔记的最后一部分以 Python 实现各算法，并演示如何借助 SciPy（HiGHS）、PuLP、NetworkX 与 CVXPY 求解同一模型。

本学习笔记以新加坡国立大学 CE6001 Operations and Management of Infrastructure Systems（基础设施系统的运营与管理）课程讲义（主讲教师：孟强教授，作者于 2021/22 学年修读）为基础，并结合个人学习笔记与经典参考书，尤其是 Taha 的 Operations Research: An Introduction、Boyd 与 Vandenberghe 的 Convex Optimization，以及 Bazaraa、Sherali 与 Shetty 的 Nonlinear Programming，将全部材料重新组织为符号统一的学习笔记。

本笔记于 2026 年 11 月在 Claude Opus 5.5（Anthropic）的协助下完成内容审阅、修订与整理，译为简体中文，并构建了网站。

作者水平有限，书中难免存在疏漏与错误，恳请读者不吝指正。如发现任何问题，欢迎通过 [GitHub Issue](https://github.com/betterweizh/lecture-notes-operation-research/issues) 提交反馈，共同完善这份笔记。

:::{admonition} 英文版
:class: tip

本书提供完整的英文版。点击任意页面右上角的翻译图标 {fas}`language`，即可打开另一种语言中的同一章节。
:::

## 全书结构

- **第一部分 — 预备知识** 复习全书所用的线性代数（线性方程组、主元变换、特征值、有定性）和向量与矩阵微积分（梯度、黑塞矩阵、泰勒展开），并介绍运筹学：它的历史、建模过程以及本书涉及的问题类别。
- **第二部分 — 凸优化** 讨论凸集、凸函数、凸优化问题与拉格朗日对偶。
- **第三部分 — 线性规划** 依次讲述模型建立、线性规划的几何、含人工变量的单纯形法、对偶理论、修正单纯形法与对偶单纯形法，以及灵敏度分析。
- **第四部分 — 网络优化** 讨论最小生成树、最短路、最大流与最小费用流。
- **第五部分 — 动态规划** 介绍最优性原理及其在设备更新、资源分配和最短路问题中的应用。
- **第六部分 — 非线性规划** 讨论模型、最优性条件与 Karush–Kuhn–Tucker（KKT）条件、线搜索，以及无约束问题和约束问题的算法。
- **第七部分 — 线性化** 说明如何把逻辑条件、乘积项和分段线性函数转化为混合整数线性约束。
- **第八部分 — 例题与习题** 按主题逐一讲解第三至第六部分的例题与习题。
- **第九部分 — 实现** 介绍软件环境的搭建，比较各种求解器与建模工具，并用 Python 实现各算法，每个实现都与库中的求解器对照检验。

各章在其所属部分内编号：第 3.2 章是第三部分的第二章，该章中的图、表、例和公式依次编号为 3.2.1、3.2.2、…。证明以及例题与习题的解答默认折叠：点击 **证明** 或 **解答** 栏即可展开。建议先独立完成每道习题，再查看解答。

附录收录了各章通用的{doc}`notation`、英汉对照的{doc}`glossary`以及全部{doc}`references`。

## 许可

本笔记的文字、图片与代码以 MIT 许可证发布：在所有副本或其主要部分中保留版权声明“Copyright (c) 2026 Better W.Z.”及许可证文本的前提下，可以使用、复制、修改、合并、出版、发行、再许可和销售本笔记的副本。本笔记按“原样”提供，不附带任何形式的担保。本笔记所依据的课程材料以及所有引用的著作，其权利仍归各自的所有者所有。

## 引用方式

Better W.Z. (2026). *Lecture Notes on Operations Research*.
GitHub. https://github.com/betterweizh/lecture-notes-operation-research

```bibtex
@misc{betterwz2026lecture,
  author       = {Better W.Z.},
  title        = {Lecture Notes on Operations Research},
  year         = {2026},
  howpublished = {\url{https://github.com/betterweizh/lecture-notes-operation-research}},
  note         = {Lecture notes}
}
```

使用本笔记中的结论、例题或图片时，也请引用该章参考文献中给出的原始出处。

## 目录

```{toctree}
:caption: 第一部分 — 预备知识
:numbered:
:maxdepth: 1

01-01-linear-algebra
01-02-matrix-calculus
01-03-introduction
```

```{toctree}
:caption: 第二部分 — 凸优化
:numbered:
:maxdepth: 1

02-01-convex-sets
02-02-convex-functions
02-03-convex-optimization-duality
```

```{toctree}
:caption: 第三部分 — 线性规划
:numbered:
:maxdepth: 1

03-01-lp-model-building
03-02-lp-geometry
03-03-simplex-method
03-04-two-phase-big-m
03-05-lp-duality
03-06-revised-dual-simplex
03-07-sensitivity-analysis
```

```{toctree}
:caption: 第四部分 — 网络优化
:numbered:
:maxdepth: 1

04-01-networks-mst
04-02-shortest-paths
04-03-maximum-flow
04-04-minimum-cost-flow
```

```{toctree}
:caption: 第五部分 — 动态规划
:numbered:
:maxdepth: 1

05-01-dp-fundamentals
05-02-dp-equipment-replacement
05-03-dp-resource-allocation
05-04-dp-shortest-paths
```

```{toctree}
:caption: 第六部分 — 非线性规划
:numbered:
:maxdepth: 1

06-01-nlp-models
06-02-optimality-conditions
06-03-kkt-conditions
06-04-line-search
06-05-unconstrained-algorithms
06-06-constrained-algorithms
```

```{toctree}
:caption: 第七部分 — 线性化
:numbered:
:maxdepth: 1

07-01-logical-constraints
07-02-linearization-techniques
07-03-piecewise-linear-approximation
```

```{toctree}
:caption: 第八部分 — 例题与习题
:numbered:
:maxdepth: 1

08-01-examples-lp-modeling
08-02-examples-simplex
08-03-examples-duality
08-04-examples-revised-sensitivity
08-05-examples-networks
08-06-examples-dp-fundamentals
08-07-examples-dp-applications
08-08-examples-nlp-optimality
08-09-examples-nlp-algorithms
```

```{toctree}
:caption: 第九部分 — 实现
:numbered:
:maxdepth: 1

09-01-code-setup-solvers
09-02-code-lp-modeling
09-03-code-lp-algorithms
09-04-code-networks
09-05-code-dynamic-programming
09-06-code-nlp
```

```{toctree}
:caption: 附录
:maxdepth: 1

notation
glossary
references
```
