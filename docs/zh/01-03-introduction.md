# 运筹学导论

## 概述

运筹学（operations research, OR）是用数学模型进行决策的学科：把一个实际的决策问题写成由决策变量、目标和约束组成的**最优化模型**（optimization model），用算法求解该模型，并在使用之前对照实际情况检验所得的答案。本章是第一部分（预备知识）的最后一章：在{doc}`01-01-linear-algebra`和{doc}`01-02-matrix-calculus`介绍了数学工具之后，本章说明运筹学是什么、从何而来，以交通、供水和能源基础设施的规划与运营为例说明它的应用领域，并描述从问题定义到实施的建模过程。随后本章确定全书其余部分所用的术语：最优化问题的一般形式 {eq}`eq-01-03-general-form`、可行解与最优解、最优值、不可行问题与无界问题，以及局部最优与全局最优。本章最后介绍本书讨论的问题类别，并给出全书九个部分的结构图（{numref}`fig-01-03-roadmap`）和阅读建议。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 解释什么是运筹学，并概述其历史；
- 描述一项运筹学研究从问题定义到实施的各个阶段；
- 从文字描述的决策问题中识别出决策变量、参数、目标和约束，并把它写成一般形式 {eq}`eq-01-03-general-form`；
- 正确使用可行、不可行、无界、最优解、最优值、局部最优和全局最优等术语；
- 说出本书涉及的问题类别，并找到讨论它们的部分和章节。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\vx\in\R^n$ | 决策变量向量 |
| $f$ | 目标函数 |
| $g_i$, $i=1,\dots,m$ | 不等式约束函数，$g_i(\vx)\le0$ |
| $h_j$, $j=1,\dots,\ell$ | 等式约束函数，$h_j(\vx)=0$ |
| $\Omega$ | 可行集 |
| $\vx^*$ | 最优解 |
| $p^*$ | 最小化问题的最优值，$p^*=\inf\{f(\vx):\vx\in\Omega\}$ |

(sec-01-03-what-is-or)=
## 什么是运筹学

每个运营某种系统的组织都必须决定如何使用有限的资源：资金、人员、车辆、机器、能源、水和时间。公交公司要决定每条线路开行多少辆公交车以及何时开行；供水公司要决定从每个水源取多少水，以及如何安排水泵的运行；电力系统调度机构要决定开启哪些发电厂，以及每个电厂每小时发多少电；道路管理部门要决定今年翻修哪些路面，哪些可以暂缓。这类决策有三个共同特征：可能的选择很多，往往有无穷多个；选择受到物理、技术、财务或法律上的限制；就成本、出行时间、可靠性或排放等可度量的目标而言，某些选择明显优于其他选择。

**运筹学**就是把数学模型和算法应用于这类决策。这一名称源于它最初对军事*作战*（operations）的研究（见{ref}`sec-01-03-history`）；*管理科学*（management science）和*决策科学*（decision science）指的是同一领域，而 *operational research* 是英国的叫法。运筹学方法有三个特点：

1. **建模。** 用数学模型描述决策问题，保留本质性的权衡，舍去不影响决策的细节。模型把假设明确写出来，使之可以讨论和检验。
2. **优化。** 在模型允许的所有选择中，运筹学寻找对给定目标而言最好的一个，或者至少可以证明接近最好的一个。这使运筹学有别于仿真和经验法则：后两者评价一个给定的方案，但不去寻找更好的方案。
3. **算法。** 实际中的模型有成千上万乃至数百万个变量。求解它们需要正确性可以证明、计算量随问题规模合理增长的算法，并在软件（求解器）中实现。

因此，运筹学对一个模型提出三类问题：最优解在*哪里*（问题的几何结构与最优性条件），*如何*找到它们（算法），以及它们*意味着什么*（答案对数据的灵敏度，以及结果的经济解释，例如稀缺资源的价格）。本书的各个部分依次对每一类模型回答这三个问题。

运筹学关注的是系统**运营与管理**方面的决策，因此基础设施是它的一个天然应用领域。交通网络、给排水系统、电网、港口和通信网络规模大、造价高，并由众多用户共享；在规划和运营方式上的微小改进，就能节省大量的资金、能源和时间。

(sec-01-03-history)=
## 简史

今天归入运筹学的一些数学思想比这个名称出现得更早。Agner Krarup Erlang 约在 1909 年分析了哥本哈根的电话话务，为确定电话交换机的规模创立了排队论，这是一个早期的基础设施规划问题。1939 年，Leonid Kantorovich 把生产计划问题表述成今天所说的线性规划；1941 年，Frank Hitchcock 提出了运输问题（Hillier & Lieberman, 2021, Chapter 1）。

这一名称本身出现于 20 世纪 30 年代末，当时英国科学家受命研究如何在防空体系中最好地使用新建的雷达站。第二次世界大战期间，由物理学家、数学家、生物学家和工程师组成的混合团队分析了各种军事行动：雷达和高射炮的部署、商船护航队的规模、反潜巡逻的搜索方式。他们的工作在英国称为 *operational research*，在美国称为 *operations research*——美国海军和陆军航空兵也成立了类似的小组。这些团队之所以受到重视，是因为他们把作战行动当作系统来看待，收集数据，建立简单的定量模型，并在实战中检验自己的建议。

战后，这些方法转向工业和政府部门，这一领域也形成了它的数学核心。

- 1947 年，George B. Dantzig 在研究美国空军的计划问题时，提出了一般的**线性规划**（linear programming, LP）问题，并发明了求解它的**单纯形法**（simplex method）；同一年，John von Neumann 勾勒出线性规划的对偶理论。Dantzig (1963) 记述了这段起源。线性规划是{doc}`03-01-lp-model-building`至{doc}`03-07-sensitivity-analysis`的主题。
- 1951 年，Harold W. Kuhn 和 Albert W. Tucker 发表了**非线性规划**（nonlinear programming）的最优性条件（Kuhn & Tucker, 1951）。William Karush 在他 1939 年的硕士论文中已推导出同样的条件（Karush, 1939），如今这些条件以三人共同命名。这些 **Karush–Kuhn–Tucker（KKT）条件**出现在{doc}`02-03-convex-optimization-duality`和{doc}`06-03-kkt-conditions`中。
- 20 世纪 50 年代初，Richard Bellman 基于他的最优性原理，为多阶段决策问题发展了**动态规划**（dynamic programming）（Bellman, 1957）；见{doc}`05-01-dp-fundamentals`。
- 网络优化也在同一个十年中成形：Ford 和 Fulkerson (1956) 的最大流算法，Kruskal (1956) 和 Prim (1957) 的最小生成树算法，以及 Dijkstra 的最短路算法（Dijkstra, 1959）（{doc}`04-01-networks-mst`至{doc}`04-04-minimum-cost-flow`）。

随后出现了专业学会和期刊（美国运筹学会成立于 1952 年）；随着计算机的普及，运筹学成为工程、商业和公共管理的组成部分。后来的里程碑包括整数规划方法（割平面与分支定界，约 1960 年）、20 世纪 80 年代的线性规划内点法，以及自 20 世纪 90 年代以来能在台式计算机上求解含数百万变量模型的商业和开源求解器（{doc}`09-01-code-setup-solvers`）。今天，运筹学方法也是机器学习（训练模型就是一个最优化问题）以及用于实时运营基础设施的数字孪生的组成部分。

(sec-01-03-applications)=
## 运筹学的应用领域

凡是决策反复进行、涉及大笔资金或必须协调众多相互作用部分的地方，都会用到运筹学。{numref}`tab-01-03-applications`列出了基础设施系统中的典型决策、描述它们的模型类别，以及本书中讨论这些类别的章节。

```{list-table} 基础设施系统中典型的运筹学决策
:name: tab-01-03-applications
:header-rows: 1
:widths: 18 42 22 18

* - 领域
  - 决策
  - 模型类别
  - 章节
* - 交通
  - 公交和列车的线路与时刻表；乘务人员排班
  - 线性规划，整数规划
  - {doc}`3.1 <03-01-lp-model-building>`, {doc}`7.1 <07-01-logical-constraints>`–{doc}`7.2 <07-02-linearization-techniques>`
* - 交通
  - 车辆、应急服务和货运的最快或最省路线
  - 最短路
  - {doc}`4.2 <04-02-shortest-paths>`, {doc}`5.4 <05-04-dp-shortest-paths>`
* - 交通
  - 交通分配：交通量如何在拥堵的道路网络上分布，以及应新建或拓宽哪些路段
  - 凸规划与非线性规划
  - {doc}`2.2 <02-02-convex-functions>`–{doc}`2.3 <02-03-convex-optimization-duality>`, {doc}`6.1 <06-01-nlp-models>`–{doc}`6.6 <06-06-constrained-algorithms>`
* - 供水
  - 每日从水库、地下水和海水淡化厂的取水量；水泵调度
  - 线性规划，非线性规划
  - {doc}`3.1 <03-01-lp-model-building>`, {doc}`6.1 <06-01-nlp-models>`
* - 供水
  - 来水不确定时水库在整个季节中的放水策略
  - 动态规划
  - {doc}`5.1 <05-01-dp-fundamentals>`–{doc}`5.3 <05-03-dp-resource-allocation>`
* - 能源
  - 以最低成本满足需求的各发电厂出力（经济调度）；开启哪些机组（机组组合）
  - 凸规划，混合整数规划
  - {doc}`2.3 <02-03-convex-optimization-duality>`, {doc}`7.1 <07-01-logical-constraints>`–{doc}`7.3 <07-03-piecewise-linear-approximation>`
* - 能源与供水
  - 网络中管道、电缆和输电线路的容量
  - 网络流，混合整数规划
  - {doc}`4.3 <04-03-maximum-flow>`–{doc}`4.4 <04-04-minimum-cost-flow>`, {doc}`7.1 <07-01-logical-constraints>`
* - 所有领域
  - 何时维修或更换路面、桥梁、水泵和车辆
  - 动态规划
  - {doc}`5.2 <05-02-dp-equipment-replacement>`
* - 所有领域
  - 把有限的预算分配给各个项目
  - 整数规划与动态规划
  - {doc}`5.3 <05-03-dp-resource-allocation>`, {doc}`7.1 <07-01-logical-constraints>`
```

在基础设施之外，同样的模型还用于制造业的生产与库存计划、物流中的包裹路径规划、医院手术室的排程、金融中的投资组合配置以及电信中的频率分配。Hillier and Lieberman (2021, Chapter 1) 和 Taha (2017, Chapter 1) 介绍了许多这样的应用。表中提到的交通分配模型在 Sheffi (1985) 中有深入讨论。

(sec-01-03-modelling-process)=
## 建模过程

一项运筹学研究远不止求解一个方程。它通常分为五个阶段（{numref}`fig-01-03-modelling-process`）。

```{figure} figures/01-03-modelling-process.svg
:name: fig-01-03-modelling-process
:width: 680px
:class: dark-light
:alt: 从问题定义到实施的五个方框，由箭头相连，并有从验证返回模型建立和问题定义的反馈箭头

运筹学研究的各个阶段。验证常常使分析人员回到模型建立甚至问题定义阶段；整个过程是一个循环，而不是一条直线。
```

1. **问题定义。** 弄清楚要决定什么、由谁决定、目标是什么、限制条件是什么。决策者最初的描述很少是精确的：“降低供水管网的运行费用”可能指水泵的能耗费用、水处理的药剂费用或两者兼有，也可能允许或不允许改变水库的调度规则。这一阶段确定范围、目标和备选方案。
2. **模型建立。** 把问题翻译成数学：选择决策变量，把目标和约束表示为决策变量的函数，并收集数据（参数）。其中的技巧在于选择合适的详细程度：既要足以刻画各种权衡，又不能详细到模型无法求解或数据无法获得。每个模型都建立在假设之上，例如线性规划中的比例性和可加性（{ref}`sec-03-01-assumptions`），这些假设应当明确写下来。
3. **求解。** 用合适的算法求解模型，通常借助求解器。本书第二至第七部分讲解这些算法，因为了解算法的工作原理，才能知道哪些模型可以求解、模型可以有多大，以及输出结果的含义（例如对偶值和最优性范围，见{doc}`03-05-lp-duality`和{doc}`03-07-sensitivity-analysis`）。
4. **验证。** 检查模型给出的答案是否合理。输入过去的数据时，模型能否复现过去的决策？数据变化时，解的变化是否合理？实际中重要的约束是否都已包括？一个看起来好得令人难以置信的解，或者一个无界或不可行的模型（见{ref}`sec-01-03-basic-terms`），通常意味着模型的建立有错误。
5. **实施。** 把决策付诸实践，这需要取得使用者的信任；并随着系统的变化维护模型。许多模型是例行运行的，例如每隔几分钟求解一次的调度模型。

本书通篇使用这一循环的简化版本，共四步：陈述决策问题；建立模型；设计或选择求解算法；在例题和案例研究上检验模型和算法。本书着重于第二步和第三步，即模型建立和算法设计，因为这两种技能可以从一个应用迁移到另一个应用。

(sec-01-03-model-parts)=
## 最优化模型的组成部分

本书中的每个最优化模型都由四部分组成。

- **决策变量**（decision variable）是决策者能够控制的量，汇集成向量 $\vx=(x_1,\dots,x_n)^\T\in\R^n$：从每个水源取用的水量、每台发电机的出力、每条线路上的公交车数量，以及诸如“是否修建这条路段”之类的 0–1 选择。
- **参数**（parameter，即数据）是决策者在模型中无法改变的量：需求、容量、单位费用、行程时间。它们被视为已知常数；最优决策如何依赖于它们，是灵敏度分析的主题。
- **目标函数**（objective function）$f(\vx)$ 衡量决策的优劣：要最小化的费用或出行时间，或者要最大化的利润或通过量。
- **约束**（constraint）描述哪些决策是允许的：物理限制（容量、流量守恒）、技术规则（压力限值、安全裕度）、政策和预算。

:::{prf:example} 为一个片区供水
:label: ex-01-03-water-supply

一家供水公司每天必须向某个片区输送至少 120 ML（百万升）的净化水。它可以从水库取水，抽水和处理费用为每 ML 250 美元，但水库的安全出水量每天至多为 80 ML。它也可以从海水淡化厂购水，价格为每 ML 1 100 美元，以该厂每天 100 ML 的产能为上限。它应从每个水源各取多少水？把这一决策表述为最优化模型，指出其决策变量、参数、目标和约束，并求解。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*决策变量。* $x_1$ = 从水库取的水量，$x_2$ = 从海水淡化厂购买的水量，单位均为 ML/日。

*参数。* 需求 120 ML/日，产能 80 ML/日和 100 ML/日，单位费用每 ML 250 美元和 1 100 美元。

*目标。* 每日费用 $f(\vx)=250x_1+1100x_2$（美元/日），求其最小值。

*约束。* 必须满足需求，$x_1+x_2\ge120$；产能限制 $x_1\le80$ 和 $x_2\le100$；取水量不能为负，$x_1,x_2\ge0$。模型为

$$
\begin{aligned}
\min_{\vx}\ & 250x_1+1100x_2\\
\st\ & x_1+x_2\ge120,\\
& x_1\le80,\\
& x_2\le100,\\
& x_1,x_2\ge0.
\end{aligned}
$$ (eq-01-03-water-model)

*最优方案。* 水库更便宜，因此应把它用到上限，其余部分由海水淡化厂补足：$\vx^*=(80,40)$，费用为每天 $250\cdot80+1100\cdot40=\$64{,}000$。这一直觉可以变成证明。对每个可行的 $\vx$，由需求约束得 $x_2\ge120-x_1$，因此

$$
\begin{aligned}
250x_1+1100x_2
&\ge 250x_1+1100\,(120-x_1)\\
&= 132{,}000-850\,x_1\\
&\ge 132{,}000-850\cdot80\\
&= 64{,}000,
\end{aligned}
$$

其中最后一个不等式用到了 $x_1\le80$。没有哪个可行方案的费用低于每天 64 000 美元，而 $\vx^*$ 达到了这个下界。
:::

模型 {eq}`eq-01-03-water-model` 是一个线性规划问题：目标和所有约束关于 $\vx$ 都是线性的。它虽小，却显示了典型的结构：一种便宜但有限的资源，一种昂贵但充足的资源，以及一个迫使人们进行权衡的需求条件。{prf:ref}`ex-01-03-water-supply`解答末尾的论证把约束组合起来，得到目标函数的一个下界，这正是对偶理论的萌芽（{doc}`02-03-convex-optimization-duality`、{doc}`03-05-lp-duality`）。

(sec-01-03-general-form)=
## 最优化问题的一般形式

把模型的约束写成 $g_i(\vx)\le0$，$i=1,\dots,m$（不等式）和 $h_j(\vx)=0$，$j=1,\dots,\ell$（等式）。满足所有约束的决策构成的集合

$$
\begin{aligned}
\Omega=\bigl\{\vx\in\R^n:\ & g_i(\vx)\le0,\ i=1,\dots,m,\\
& h_j(\vx)=0,\ j=1,\dots,\ell\bigr\},
\end{aligned}
$$ (eq-01-03-feasible-set)

称为**可行集**（feasible set，或可行域，feasible region）。一般的最优化问题为

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)\\
\st\ & \vx\in\Omega .
\end{aligned}
$$ (eq-01-03-general-form)

每个约束都可以化为这两种形式之一：$a(\vx)\ge b(\vx)$ 变为 $b(\vx)-a(\vx)\le0$。在供水模型 {eq}`eq-01-03-water-model` 中，$g_1(\vx)=120-x_1-x_2$，$g_2(\vx)=x_1-80$，$g_3(\vx)=x_2-100$，$g_4(\vx)=-x_1$，$g_5(\vx)=-x_2$，其中 $m=5$，$\ell=0$。

最大化问题不需要单独的理论，因为

$$
\max_{\vx\in\Omega} f(\vx) = -\min_{\vx\in\Omega}\bigl(-f(\vx)\bigr),
$$ (eq-01-03-max-min)

并且两个问题有相同的最优解。本书对最小化问题叙述一般理论；线性规划各章在有帮助时同时给出两种优化方向下的规则。

(sec-01-03-basic-terms)=
## 基本术语

以下术语在全书中通用。

:::{prf:definition} 可行解与最优解
:label: def-01-03-feasible-optimal

考虑问题 {eq}`eq-01-03-general-form`。

1. 点 $\vx\in\Omega$ 称为**可行解**（feasible solution）；点 $\vx\notin\Omega$ 称为**不可行**（infeasible）的。
2. 若 $\Omega=\emptyset$，则称该问题**不可行**。
3. **最优值**（optimal value）为 $p^*=\inf\{f(\vx):\vx\in\Omega\}$，并约定：若问题不可行，则 $p^*=+\infty$；若 $f$ 在 $\Omega$ 上可以取任意小的值，则 $p^*=-\infty$，此时称该问题**无界**（unbounded）。
4. 满足 $f(\vx^*)=p^*$ 的可行点 $\vx^*$ 称为**最优解**（optimal solution，即极小点）。全部最优解构成的集合可能是空集、单点集，也可能含有无穷多个点。
:::

区分最优*解* $\vx^*$（一个决策）与最优*值* $p^*$（一个数）是重要的。在{prf:ref}`ex-01-03-water-supply`中，最优解是 $\vx^*=(80,40)$ ML/日，最优值是每天 $p^*=\$64{,}000$。如果海水淡化的费用也是每 ML 250 美元，则最优值为每天 30 000 美元，并且每个满足 $x_1+x_2=120$、$20\le x_1\le80$ 的方案都是最优的：最优值只有一个，最优解却有无穷多个（多重最优解，{ref}`sec-03-03-special-cases`）。

可能出现三种问题，每一种在实践中都有其典型原因。

- **不可行。** 如果该片区每天需要 200 ML，则约束 $x_1+x_2\ge200$、$x_1\le80$、$x_2\le100$ 不可能同时成立，因为两个水源每天至多供水 180 ML。这时模型没有可行解。在实践中，模型不可行意味着要求过于严格，或者数据有错误。
- **无界。** 假设模型还允许供水公司以每 ML 1 300 美元的价格把多余的水卖给相邻片区，而建模者忘记了海水淡化厂的产能。那么每多以 1 100 美元买进、以 1 300 美元卖出 1 ML，费用就降低 200 美元，没有止境，于是 $p^*=-\infty$。实际系统绝不会无界；无界的模型一定遗漏了某个约束。
- **最优值取不到。** 下确界可以有限却取不到。在 $\Omega=\{x\in\R:x\ge1\}$ 上极小化 $f(x)=1/x$，得 $p^*=0$，但对每个可行的 $x$ 都有 $f(x)>0$，因此不存在最优解。魏尔斯特拉斯定理（{prf:ref}`thm-02-01-weierstrass`）给出了最小值必定存在的条件（在非空、闭且有界的可行集上的连续目标函数）。

最后，一个解可能在其邻近的解中最好，却不是整体上最好的。

:::{prf:definition} 局部极小点与全局极小点
:label: def-01-03-local-global

若对所有 $\vx\in\Omega$ 有 $f(\vx^*)\le f(\vx)$，则称可行点 $\vx^*$ 为 {eq}`eq-01-03-general-form` 的**全局极小点**（global minimum）。若存在 $\varepsilon>0$，使得对所有满足 $\norm{\vx-\vx^*}<\varepsilon$ 的 $\vx\in\Omega$ 都有 $f(\vx^*)\le f(\vx)$，则称它为**局部极小点**（local minimum）。
:::

每个全局极小点都是局部极小点，反之则不然。如果海水淡化厂有很强的规模经济，其每 ML 的费用随购买量增加而急剧下降，那么费用函数就不再是线性的，模型可能有一个主要使用水库水的局部极小点，以及另一个主要使用淡化水的全局极小点。只会一步步改进解的算法可能停在错误的那一个上。对于**凸**问题，这种情况不会发生：每个局部极小点都是全局极小点（{prf:ref}`thm-02-03-local-global`）。这就是本书在预备知识之后紧接着就在第二部分研究凸性（{doc}`02-01-convex-sets`、{doc}`02-02-convex-functions`）的主要原因。一般非线性规划的局部最优与全局最优是{doc}`06-01-nlp-models`和{doc}`06-02-optimality-conditions`的主题。

(sec-01-03-problem-classes)=
## 本书涉及的问题类别

一般形式 {eq}`eq-01-03-general-form` 涵盖了本书中的所有模型。各个类别的区别在于 $f$、$g_i$、$h_j$ 以及变量的类型，而这种类型决定了可以使用哪些算法。

```{list-table} 本书讨论的问题类别
:name: tab-01-03-problem-classes
:header-rows: 1
:widths: 24 46 30

* - 类别
  - 结构
  - 主要方法（部分，章节）
* - 凸优化
  - $f$、$g_i$ 为凸函数，$h_j$ 为仿射函数
  - 凸性、拉格朗日对偶、KKT 条件（第二部分，{doc}`2.1 <02-01-convex-sets>`–{doc}`2.3 <02-03-convex-optimization-duality>`）
* - 线性规划（LP）
  - $f$、$g_i$、$h_j$ 为线性（仿射）函数；连续变量
  - 单纯形法、对偶、灵敏度分析（第三部分，{doc}`3.1 <03-01-lp-model-building>`–{doc}`3.7 <03-07-sensitivity-analysis>`）
* - 网络优化
  - 图上的线性规划：节点处的流量守恒，弧的容量和费用
  - 生成树、最短路、最大流、最小费用流（第四部分，{doc}`4.1 <04-01-networks-mst>`–{doc}`4.4 <04-04-minimum-cost-flow>`）
* - 动态规划（DP）
  - 分阶段决策；目标函数在各阶段上累加
  - 值函数的递推（第五部分，{doc}`5.1 <05-01-dp-fundamentals>`–{doc}`5.4 <05-04-dp-shortest-paths>`）
* - 非线性规划（NLP）
  - $f$、$g_i$、$h_j$ 光滑，可能非凸
  - 最优性条件、线搜索、下降法与罚函数法（第六部分，{doc}`6.1 <06-01-nlp-models>`–{doc}`6.6 <06-06-constrained-algorithms>`）
* - 混合整数规划（MIP）
  - 部分变量限定为整数或取值于 $\{0,1\}$
  - 用 0–1 变量建模、线性化（第七部分，{doc}`7.1 <07-01-logical-constraints>`–{doc}`7.3 <07-03-piecewise-linear-approximation>`）
```

这些类别相互重叠。线性规划是凸问题；网络流问题是具有特殊结构的线性规划，更快的算法利用了这种结构；最短路问题可以作为线性规划求解，也可以用 Dijkstra 算法或动态规划求解；非线性关系或逻辑关系常常可以转化为混合整数*线性*模型（{doc}`07-02-linearization-techniques`、{doc}`07-03-piecewise-linear-approximation`），交给 MIP 求解器处理。

(sec-01-03-roadmap)=
## 全书结构图

{numref}`fig-01-03-roadmap`展示了全书的九个部分及其相互依赖关系。

```{figure} figures/01-03-roadmap.svg
:name: fig-01-03-roadmap
:width: 760px
:class: dark-light
:alt: 全书的框图：顶部是第一部分预备知识的横条；中间是第二至第七部分（凸优化、线性规划、网络优化、动态规划、非线性规划、线性化）的六个方框，由依赖关系箭头相连；底部是第八部分例题与习题和第九部分实现的横条

全书的九个部分及其依赖关系。实线箭头表示后面的部分用到前面的部分；虚线连接从两个侧面讨论同一问题的部分（由线性规划对偶推导 KKT 条件；用标号法和动态规划求最短路）。第八部分汇集了第三至第六部分的例题与习题，第九部分用 Python 实现各方法。
```

**第一部分**（预备知识）汇集了后面要用的工具：矩阵与线性方程组、特征值与二次型的有定性（{doc}`01-01-linear-algebra`），梯度、黑塞矩阵与泰勒展开（{doc}`01-02-matrix-calculus`），以及本章中最优化模型的术语。前两章可以先读，也可以在需要时查阅；本章则值得在阅读任何后续部分之前先读。

**第二部分**（凸优化），即{doc}`02-01-convex-sets`至{doc}`02-03-convex-optimization-duality`，讨论凸集、凸函数和拉格朗日对偶。它为线性规划的几何（{doc}`03-02-lp-geometry`）、线性规划对偶（{doc}`03-05-lp-duality`）以及非线性规划的最优性条件（{doc}`06-02-optimality-conditions`、{doc}`06-03-kkt-conditions`）提供了语言。

**第三部分**（线性规划），即{doc}`03-01-lp-model-building`至{doc}`03-07-sensitivity-analysis`，建立线性规划模型，用单纯形法及其变形求解，并分析对偶性与灵敏度。**第四部分**（网络优化），即{doc}`04-01-networks-mst`至{doc}`04-04-minimum-cost-flow`，用到线性规划理论，但有其自身更快的算法。**第五部分**（动态规划），即{doc}`05-01-dp-fundamentals`至{doc}`05-04-dp-shortest-paths`，除第一部分外几乎不需要其他预备知识，可以随时阅读；其最后一章重新讨论{doc}`04-02-shortest-paths`。**第六部分**（非线性规划），即{doc}`06-01-nlp-models`至{doc}`06-06-constrained-algorithms`，建立在凸性与对偶的基础上。**第七部分**（线性化），即{doc}`07-01-logical-constraints`至{doc}`07-03-piecewise-linear-approximation`，说明如何用 0–1 变量对逻辑条件和非线性函数建模，使混合整数线性求解器能够处理它们。

**第八部分**（例题与习题）为第三至第六部分提供了更多附折叠解答的例题与习题。理论各章只保留发展方法所必需的例题，其解答同样是折叠的；其余的都在第八部分。无论在哪里，都请先自己尝试解题，再展开解答。**第九部分**（实现）说明如何用 Python 求解这些模型，以及如何编程实现第三至第六部分的算法；代码复现了例题中的数值。

主要关心线性模型的读者可以依次阅读第三和第四部分（{doc}`03-01-lp-model-building`至{doc}`04-04-minimum-cost-flow`），然后阅读第七部分（{doc}`07-01-logical-constraints`至{doc}`07-03-piecewise-linear-approximation`）。主要关心非线性模型的读者应先读第二部分（{doc}`02-01-convex-sets`至{doc}`02-03-convex-optimization-duality`），再读第六部分（{doc}`06-01-nlp-models`至{doc}`06-06-constrained-algorithms`），并以{doc}`03-05-lp-duality`作为 KKT 条件的背景知识。无论哪种方式，都可以在阅读每个理论章节的同时阅读第八和第九部分中相应的章节。

## 小结

- 运筹学用数学模型和算法为资源有限的系统做决策；基础设施系统（交通、供水、能源）是其主要应用领域（{numref}`tab-01-03-applications`）。
- 它的数学核心形成于 1939 年至 1960 年间：线性规划与单纯形法（1947 年）、KKT 条件（1939 年、1951 年）、动态规划与网络算法（20 世纪 50 年代）。
- 一项运筹学研究要经历问题定义、模型建立、求解、验证和实施等阶段，通常不止一轮（{numref}`fig-01-03-modelling-process`）。
- 模型由决策变量、参数、目标和约束组成；本书中的每个模型都具有一般形式 {eq}`eq-01-03-general-form`，其可行集为 {eq}`eq-01-03-feasible-set`，最大化问题可以通过 {eq}`eq-01-03-max-min` 化为最小化问题。
- 一个问题可能不可行、无界、有有限但取不到的最优值，或者有一个或多个最优解（{prf:ref}`def-01-03-feasible-optimal`）；局部极小点不一定是全局极小点（{prf:ref}`def-01-03-local-global`），凸问题除外。
- 全书共九个部分：预备知识（第一部分）；凸优化、线性规划、网络优化、动态规划、非线性规划以及面向混合整数模型的线性化（第二至第七部分）；例题与习题（第八部分）；以及 Python 实现（第九部分）（{numref}`tab-01-03-problem-classes`、{numref}`fig-01-03-roadmap`）。

## 参考文献

- Bellman, R. (1957). *Dynamic programming*. Princeton University Press.
  (Preface and Chapter 3, the principle of optimality.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Chapter 1,
  mathematical optimization and its problem classes.)
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Chapters 1–2, the origins of linear programming.)
- Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
  *Numerische Mathematik, 1*(1), 269–271. <https://doi.org/10.1007/BF01386390>
- Ford, L. R., & Fulkerson, D. R. (1956). Maximal flow through a network.
  *Canadian Journal of Mathematics, 8*, 399–404.
  <https://doi.org/10.4153/CJM-1956-045-5>
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 1, the origins and nature of
  OR; Chapter 2, overview of the modeling approach.)
- Karush, W. (1939). *Minima of functions of several variables with
  inequalities as side constraints* [Master's thesis, University of
  Chicago].
- Kruskal, J. B. (1956). On the shortest spanning subtree of a graph and
  the traveling salesman problem. *Proceedings of the American Mathematical
  Society, 7*(1), 48–50. <https://doi.org/10.1090/S0002-9939-1956-0078686-7>
- Kuhn, H. W., & Tucker, A. W. (1951). Nonlinear programming. In J. Neyman
  (Ed.), *Proceedings of the Second Berkeley Symposium on Mathematical
  Statistics and Probability* (pp. 481–492). University of California
  Press.
- Prim, R. C. (1957). Shortest connection networks and some
  generalizations. *Bell System Technical Journal, 36*(6), 1389–1401.
  <https://doi.org/10.1002/j.1538-7305.1957.tb01515.x>
- Sheffi, Y. (1985). *Urban transportation networks: Equilibrium analysis
  with mathematical programming methods*. Prentice-Hall. (Chapter 1,
  network representation and equilibrium.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 1, what is operations research.)
- Williams, H. P. (2013). *Model building in mathematical programming*
  (5th ed.). Wiley. (Chapter 1, the concept of a model and mathematical
  programming models.)
