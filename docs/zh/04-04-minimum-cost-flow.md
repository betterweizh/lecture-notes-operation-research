# 最小费用流

## 概述

{doc}`04-03-maximum-flow`中的最大流问题只问能*送多少*，而不考虑运送的代价。大多数物流问题恰好相反：已知每个工厂的供应量和每个市场的需求量，要通过由道路、仓库和管道组成的网络运送货物，在满足弧容量以及（有时）合同规定的最小流量的前提下，使**总费用最小**。这就是**最小费用流问题**（minimum-cost flow problem）。它把最短路问题、最大流问题、运输问题、指派问题和转运问题都作为特例包含在内，而它本身仍是一个线性规划，只是约束矩阵非常特殊：节点–弧关联矩阵（node–arc incidence matrix）。这个矩阵是**全单模**（totally unimodular）的，因此整数的供应量和容量给出整数的最优流。本章先给出模型，说明如何处理不平衡的数据和流量下界，证明整数性，再从两种等价的形式推导最优性条件（带节点势的检验数，以及不存在负费用圈），并由此得到算法：连续最短路算法、消圈算法，以及网络单纯形法的概要。本章用到{doc}`04-03-maximum-flow`中的残量网络和增广路、{doc}`04-02-shortest-paths`中带标号的最短路、{doc}`03-02-lp-geometry`中的基本可行解，以及{doc}`03-05-lp-duality`中的对偶线性规划和互补松弛性。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 把带有供应量、需求量、单位费用以及流量上下界的最小费用流问题写成线性规划；
- 识别不可行或不平衡的数据，添加虚拟节点，并通过变量替换去掉流量下界；
- 证明节点–弧关联矩阵是全单模的，并解释为什么整数数据给出整数最优流；
- 把最短路问题、最大流问题、运输问题、指派问题和转运问题写成最小费用流问题；
- 用节点势和检验数，或者通过在残量网络中搜索负圈，检验一个流的最优性；
- 用连续最短路算法或消圈算法求解小规模实例，并描述网络单纯形法。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $b_i$ | 节点 $i$ 的净供应量：$b_i>0$ 为供应，$b_i<0$ 为需求，$b_i=0$ 为转运节点（单位） |
| $c_{ij}$ | 弧 $(i,j)$ 上单位流量的费用（美元/单位） |
| $l_{ij}$，$u_{ij}$ | 流量 $x_{ij}$ 的下界和上界（单位） |
| $z$ | 总费用 $\sum c_{ij}x_{ij}$（美元） |
| $\pi_i$ | 节点 $i$ 的势（对偶变量）（美元/单位） |
| $y_{ij}$ | 上界约束 $x_{ij}\le u_{ij}$ 的对偶变量 |
| $\bar c_{ij}=c_{ij}-\pi_i+\pi_j$ | 弧 $(i,j)$ 的检验数 |
| $G(\vx)$ | $\vx$ 的残量网络，前向残量弧的费用为 $c_{ij}$，后向残量弧的费用为 $-c_{ij}$ |
| $\mathcal W$，$c(\mathcal W)$ | $G(\vx)$ 中的一个有向圈及其费用 |

(sec-04-04-model)=
## 最小费用流问题

设 $G=(\mathcal N,\mathcal A)$ 是一个有向网络。每个节点 $i$ 有**净供应量**（net supply）$b_i$：若 $b_i>0$，它生产 $b_i$ 个单位；若 $b_i<0$，它消耗 $-b_i$ 个单位；若 $b_i=0$，它只让流量通过。每条弧有单位费用 $c_{ij}$ 和界 $l_{ij}\le x_{ij}\le u_{ij}$（未给出界时取 $l_{ij}=0$，$u_{ij}=\infty$）。

:::{prf:definition} 最小费用流问题
:label: def-04-04-mcf

**最小费用流问题**是如下线性规划：

$$
\begin{aligned}
\min_{\vx}\ & z=\sum_{(i,j)\in\mathcal A}c_{ij}x_{ij} \\
\st\ & \sum_{j:(i,j)\in\mathcal A}x_{ij}-\sum_{j:(j,i)\in\mathcal A}x_{ji}=b_i, \quad i\in\mathcal N,\\
& l_{ij}\le x_{ij}\le u_{ij}, \quad (i,j)\in\mathcal A.
\end{aligned}
$$ (eq-04-04-mcf)

其中的等式约束是**流量平衡**（flow-balance，即流量守恒）方程：从节点 $i$ 流出的量减去流入的量等于它的净供应量。写成矩阵形式即 $\mA\vx=\vb$，其中 $\mA$ 是{ref}`sec-04-01-representations`中的节点–弧关联矩阵，弧 $(i,j)$ 对应的列为 $\ve_i-\ve_j$。
:::

:::{prf:example} 弧上有界的供应链
:label: ex-04-04-supply-chain

三个工厂（节点 1、2、3）分别供应某产品 100、200 和 50 个单位；三个市场（节点 4、5、6）分别需要 150、80 和 120 个单位。节点 3 既是工厂，也是其他工厂可以经由它转运的仓库。各条弧、以美元计的单位费用以及由合同和线路能力规定的界如下：

| 弧 | $(1,2)$ | $(1,3)$ | $(1,4)$ | $(2,3)$ | $(2,5)$ |
|---|---|---|---|---|---|
| $c_{ij}$（美元） | 3 | 4 | 1 | 5 | 6 |
| $[l_{ij},u_{ij}]$ | – | – | $[50,80]$ | – | – |

| 弧 | $(3,4)$ | $(3,5)$ | $(4,6)$ | $(5,6)$ |
|---|---|---|---|---|
| $c_{ij}$（美元） | 1 | 2 | 2 | 4 |
| $[l_{ij},u_{ij}]$ | $[70,120]$ | – | $[100,120]$ | – |

（短横线表示 $[0,\infty)$。）求费用最小的运输方案。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

取 $\vb=(100,200,50,-150,-80,-120)$，模型 {eq}`eq-04-04-mcf` 为

$$
\begin{aligned}
\min\ & x_{14}+4x_{13}+3x_{12}+5x_{23}+6x_{25}+x_{34}+2x_{35}+2x_{46}+4x_{56}\\
\st\ & x_{12}+x_{13}+x_{14}=100\\
& x_{23}+x_{25}-x_{12}=200\\
& x_{34}+x_{35}-x_{13}-x_{23}=50\\
& x_{46}-x_{14}-x_{34}=-150\\
& x_{56}-x_{25}-x_{35}=-80\\
& -x_{46}-x_{56}=-120\\
& 50\le x_{14}\le80\\
& 70\le x_{34}\le120\\
& 100\le x_{46}\le120\\
& \text{其余 } x_{ij}\ge0.
\end{aligned}
$$

*（节点 2 的平衡方程中减去的是流入量 $x_{12}$；若该行写成 $x_{13}$，模型就不再描述这个网络。）*

这组数据实际上是不可行的（{ref}`sec-04-04-feasibility`）。去掉弧 $(4,6)$ 上的下界后，{numref}`tab-04-04-supply-chain-solution`中的最优方案为 $x_{13}=20$，$x_{14}=80$，$x_{23}=50$，$x_{25}=150$，$x_{34}=120$，$x_{46}=50$，$x_{56}=70$，弧 $(1,2)$ 和 $(3,5)$ 上不运货，费用为 1 810 美元（{ref}`sec-04-04-supply-chain-solved`）。
:::

在求解{prf:ref}`ex-04-04-supply-chain`的模型之前，下一节先检查解是否存在。

(sec-04-04-feasibility)=
## 可行性、平衡与虚拟节点

**平衡。** 把所有平衡方程相加得到 $0=\sum_i b_i$，因为每个 $x_{ij}$ 恰好以 $+$ 号和 $-$ 号各出现一次。因此 $\sum_ib_i\ne0$ 的问题按原样是不可行的，必须先确定这种不平衡意味着什么：

- **供过于求**（$\sum_i b_i>0$）：生产的产品不能全部售出。添加一个**虚拟需求节点**（dummy demand node）$d$，令 $b_d=-\sum_ib_i$，并从每个供应节点向 $d$ 引一条费用为 0（或为储存、处置费用）的弧。流入 $d$ 的流量就是各工厂未用完的产能。
- **供不应求**（$\sum_ib_i<0$）：有些需求无法满足。添加一个**虚拟供应节点**，令 $b=-\sum_ib_i$，并向每个需求节点引一条弧，其费用等于每短缺一个单位的惩罚（损失的利润）。从它流出的流量就是各市场的短缺量。

两种情形下方程组都变得相容；在连通网络上 $\mA$ 的秩为 $\abs{\mathcal N}-1$，所以有一个平衡方程是冗余的。

**容量。** 即使数据平衡，也可能因为界而不可行。对节点集 $S$，流出 $S$ 的流量减去流入 $S$ 的流量必须等于净供应量 $\sum_{i\in S}b_i$；这只有在下式成立时才可能：

$$
\sum_{i\in S}b_i\le\sum_{(i,j):\,i\in S,\,j\notin S}u_{ij}-\sum_{(j,i):\,j\notin S,\,i\in S}l_{ji}
\quad\text{对每个 } S\subseteq\mathcal N,
$$ (eq-04-04-gale)

也就是说，$S$ 的净供应量在扣除流入 $S$ 的弧上的强制流量之后，必须能通过从 $S$ 流出的弧送出。Gale（1957）和 Hoffman（1960）证明了这些条件（连同 $\vl\le\vu$ 和平衡条件）也是充分的。实践中只需一次最大流计算即可检验可行性：去掉流量下界（{ref}`sec-04-04-lower-bounds`）后，用容量为 $b_i$ 的弧把超级源点（super-source）连到每个供应节点 $i$，用容量为 $-b_j$ 的弧把每个需求节点 $j$ 连到超级汇点（super-sink）；问题可行当且仅当最大流使所有源点弧饱和。

**供应链例子是不可行的。** 在 {eq}`eq-04-04-gale` 中取 $S=\mathcal N\setminus\{4\}$：它的净供应量为 $-b_4=150$，从它流出的弧 $(1,4)$ 和 $(3,4)$ 最多能运 $80+120=200$，而流入它的弧 $(4,6)$ 至少要运 100，于是条件 $150\le200-100$ 不成立。换句话说：节点 4 必须向它的市场交付 150 个单位，并至少向节点 6 送出 $l_{46}=100$ 个单位，因此它至少需要 250 个单位的流入量，但它的两条入弧最多只能运 200。线性规划求解器也确认，在这组数据下 {eq}`eq-04-04-mcf` 没有可行解。必须放松某个界。改动最小的是弧 $(4,6)$ 的下界：任何 $l_{46}\le50$ 都使问题可行。本章其余部分采用修正后的界 $0\le x_{46}\le120$（对每个 $l_{46}\le50$，最优解都相同）。

(sec-04-04-lower-bounds)=
## 去掉流量下界

去掉流量下界的方法是从下界起算每条弧的流量。作替换

$$
\begin{aligned}
x_{ij}&=l_{ij}+x'_{ij},\\
0\le x'_{ij}&\le u_{ij}-l_{ij}.
\end{aligned}
$$ (eq-04-04-shift)

强制的 $l_{ij}$ 个单位无论如何都要离开节点 $i$ 到达节点 $j$，所以净供应量变为

$$
b'_i=b_i-\sum_{j:(i,j)\in\mathcal A}l_{ij}+\sum_{j:(j,i)\in\mathcal A}l_{ji},
$$ (eq-04-04-shifted-supply)

而费用变为 $z=\sum c_{ij}x'_{ij}+\sum c_{ij}l_{ij}$，即一个常数加上一个下界为零的问题。因此下面的每个算法都可以假设 $l_{ij}=0$。

对修正后的供应链，有界的弧是 $(1,4)$（界为 $[50,80]$）和 $(3,4)$（界为 $[70,120]$）。平移后得到

- $b'_1=100-50=50$，$b'_3=50-70=-20$，$b'_4=-150+50+70=-30$，其余 $b'_i$ 不变；
- $0\le x'_{14}\le30$ 和 $0\le x'_{34}\le50$；
- 常数 $1\cdot50+1\cdot70=\$120$。

节点 3 从供应节点变成了需求节点：它自己的 50 个单位再加 20 个单位已经被弧 $(3,4)$ 占用。

(sec-04-04-integrality)=
## 全单模性与整数流

:::{prf:definition} 全单模矩阵
:label: def-04-04-tu

若一个矩阵的每个方子矩阵的行列式都等于 $0$、$+1$ 或 $-1$，则称它是**全单模**（TU）的。特别地，它的每个元素都是 $0$ 或 $\pm1$。
:::

:::{prf:theorem} 关联矩阵是全单模的
:label: thm-04-04-tu

有向网络的节点–弧关联矩阵 $\mA$ 是全单模的。约束 $\mA\vx=\vb$，$\vx+\vw=\vu$（带松弛变量 $\vw$ 的上界约束）的矩阵 $\begin{pmatrix}\mA&\vzero\\ \mI&\mI\end{pmatrix}$ 也是全单模的。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

设 $\mQ$ 是 $\mA$ 的一个 $k\times k$ 子矩阵，对 $k$ 用归纳法。$k=1$ 时，该元素为 $0$ 或 $\pm1$。$k>1$ 时，$\mQ$ 的每一列至多有一个 $+1$ 和至多一个 $-1$，因为 $\mA$ 的每一列都是如此。于是只剩三种情形。若 $\mQ$ 有零列，则 $\det\mQ=0$。若某一列只有一个非零元素，就按这一列展开 $\det\mQ$：$\det\mQ=\pm\det\mQ'$，其中 $\mQ'$ 是 $(k-1)\times(k-1)$ 子矩阵，由归纳假设其行列式为 $0$ 或 $\pm1$。否则每一列恰有一个 $+1$ 和一个 $-1$，$\mQ$ 的各行相加得零向量，故 $\det\mQ=0$。

对第二个矩阵，任取一个方子矩阵。属于 $\vw$ 的列在整个矩阵中只有一个 1，所以在子矩阵中它要么是零列（行列式为 0），要么是单位列，可以沿它展开。反复这样做，就去掉了所有 $\vw$ 列，剩下 $\begin{pmatrix}\mA\\\mI\end{pmatrix}$ 的一个方子矩阵；其中来自 $\mI$ 的行是零行或单位行，沿它们展开后剩下 $\mA$ 的一个方子矩阵。
:::

:::{prf:theorem} 最小费用流的整数性
:label: thm-04-04-integrality

若所有 $b_i$、$l_{ij}$ 和 $u_{ij}$ 都是整数，则 {eq}`eq-04-04-mcf` 的每个基本可行解都是整数解。因此，若问题有最优解，则它有整数最优解，而且单纯形法能找到它。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

经过保持数据为整数的平移 {eq}`eq-04-04-shift` 之后，像{prf:ref}`thm-04-04-tu`那样用松弛变量写出上界约束，并删去一个冗余的平衡方程行。基本可行解（{prf:ref}`def-03-02-basic-solution`）满足 $\mB\vx_B=\vb'$，其中 $\mB$ 是某个 TU 矩阵的非奇异方子矩阵，因此 $\det\mB=\pm1$。由克拉默法则（Cramer's rule），$\vx_B$ 的每个分量都是一个整数行列式与 $\pm1$ 之比，即整数。非基变量为 0。有最优解的标准形线性规划必有最优基本可行解（{prf:ref}`thm-03-02-fundamental`）。
:::

这就是网络问题“虽然是线性规划，却具有组合性质”的原因：下面的指派问题需要 0–1 决策，{ref}`sec-04-03-matching`中的匹配问题需要完整的配对，而两者都作为普通的线性规划来求解。费用 $c_{ij}$ 不必是整数。

(sec-04-04-special-cases)=
## 特例

**最短路。** 从 $s$ 向 $t$ 运送一个单位：$b_s=1$，$b_t=-1$，其余 $b_i=0$，$u_{ij}=\infty$，$c_{ij}$ 为弧长。整数最优流是一个 0–1 向量，其中取值为 1 的弧构成一条最短的 $s$–$t$ 路（前提是不存在负圈）。若取 $b_s=\abs{\mathcal N}-1$，其余每个节点 $b_i=-1$，则解是一棵最短路树（{prf:ref}`prop-04-02-tree`）。

**最大流。** 令所有费用为 0、所有 $b_i=0$，并添加一条容量无穷、费用为 $-1$ 的弧 $(t,s)$。这时流是一个循环流（circulation），其费用等于从 $t$ 返回 $s$ 的流量的相反数，即 $s$–$t$ 流值的相反数；最小费用为 $-v^*$。

**运输问题。** 节点分为供应点（工厂，$b_i>0$）和需求点（市场，$b_j<0$），每条弧都从工厂指向市场。平衡方程表示每个工厂运出其供应量、每个市场收到其需求量：

$$
\begin{aligned}
\sum_j x_{ij}&=b_i &&\text{（工厂）},\\
\sum_i x_{ij}&=-b_j &&\text{（市场）}.
\end{aligned}
$$ (eq-04-04-transportation)

**指派问题。** 这是有 $n$ 个工厂（工人）和 $n$ 个市场（任务）、每个 $b_i=1$ 且每个 $b_j=-1$ 的运输问题。由{prf:ref}`thm-04-04-integrality`，最优基本解满足 $x_{ij}\in\{0,1\}$，每行每列恰有一个 1，即一个置换。虽然有专门的算法（Kuhn 于 1955 年提出的匈牙利算法），但线性规划本身已经给出精确解。

**转运问题。** 存在 $b_i=0$ 的中间节点（仓库、港口、越库配送中心），货物可以经由它们中转。任何无容量限制的最小费用流问题都是转运问题；{prf:ref}`ex-08-05-transshipment`在此基础上增加了仓库容量和一个虚拟市场。

:::{prf:example} 一个运输问题
:label: ex-04-04-transportation

两个工厂 A 和 B 每周分别生产某产品 40 吨和 50 吨；三个市场分别需要 20、25 和 45 吨。每吨的运费（美元）为

| | 市场 1 | 市场 2 | 市场 3 | 供应量 |
|---|---|---|---|---|
| 工厂 A | 4 | 6 | 9 | 40 |
| 工厂 B | 2 | 7 | 5 | 50 |
| 需求量 | 20 | 25 | 45 | 90 |

求费用最小的运输方案。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

数据是平衡的（$40+50=20+25+45=90$）；网络就是{numref}`fig-04-04-transportation`中的二部网络，其中 $b_A=40$，$b_B=50$，$b_1=-20$，$b_2=-25$，$b_3=-45$，共有六条无容量限制的弧。连续最短路算法（{ref}`sec-04-04-ssp`，{numref}`tab-04-04-ssp`）求得最优方案 $x_{A1}=15$，$x_{A2}=25$，$x_{B1}=5$，$x_{B3}=45$，费用为每周 445 美元。
:::

:::{prf:example} 一个指派问题
:label: ex-04-04-assignment

要给三名工人 W1、W2、W3 在任务 J1、J2、J3 中各分配一项任务。每名工人完成每项任务所需的时间（小时）为

| | J1 | J2 | J3 |
|---|---|---|---|
| W1 | 9 | 4 | 7 |
| W2 | 3 | 6 | 5 |
| W3 | 8 | 7 | 4 |

使总时间最短。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

作为最小费用流问题，每名工人是 $b=1$ 的节点，每项任务是 $b=-1$ 的节点，九条弧（工人，任务）的费用就是上表中的时间。不加任何整数约束求解线性规划松弛，得到 $x_{W1,J2}=x_{W2,J1}=x_{W3,J3}=1$，其余 $x=0$：W1 做 J2，W2 做 J1，W3 做 J3，共用 $4+3+4=11$ 小时。列出全部 $3!=6$ 个置换即可验证（它们的总时间分别为 11、17、17、19、21 和 21 小时）。线性规划不需要被告知答案必须是一个置换；全单模性保证了这一点。下一节的节点势可以证明这个结果：取 $\pi_{W1}=1$，$\pi_{W2}=3$，$\pi_{W3}=4$，$\pi_{J1}=0$，$\pi_{J2}=-3$，$\pi_{J3}=0$，则每个检验数 $c_{ij}-\pi_i+\pi_j$ 都非负（例如 W1–J1 的检验数为 $9-1+0=8$），而选中的三条弧的检验数为 0。
:::

(sec-04-04-optimality)=
## 最优性条件

### 检验数与节点势

假设 $\vl=\vzero$（经过 {eq}`eq-04-04-shift` 之后），给节点 $i$ 的平衡方程配上对偶变量 $\pi_i$，给上界约束 $x_{ij}\le u_{ij}$ 配上 $y_{ij}\ge0$。{eq}`eq-04-04-mcf` 的对偶问题（{doc}`03-05-lp-duality`）为

$$
\begin{aligned}
\max_{\vpi,\,\vy}\ & \sum_i b_i\pi_i-\sum_{(i,j)}u_{ij}y_{ij}\\
\st\ & \pi_i-\pi_j-y_{ij}\le c_{ij}, \quad (i,j)\in\mathcal A,\\
& y_{ij}\ge0,
\end{aligned}
$$ (eq-04-04-dual)

其中 $\vpi$ 无符号限制。$\pi_i$ 称为**节点势**（node potential），而

$$
\bar c_{ij}=c_{ij}-\pi_i+\pi_j
$$ (eq-04-04-reduced-cost)

称为弧 $(i,j)$ 的**检验数**（reduced cost），它正是单纯形法中列 $\ve_i-\ve_j$ 的检验数 $c_j-\vpi^\T\mA_{\cdot j}$。给所有节点势加上同一个常数不改变任何检验数（因为有一个平衡方程是冗余的），所以可以固定一个节点势，例如取为 0。

:::{prf:theorem} 检验数最优性条件
:label: thm-04-04-reduced-cost

{eq}`eq-04-04-mcf` 的可行流 $\vx$ 是最优的，当且仅当存在节点势 $\vpi$，使得对每条弧都有

$$
\begin{aligned}
\bar c_{ij}>0 &\ \Rightarrow\ x_{ij}=l_{ij},\\
\bar c_{ij}<0 &\ \Rightarrow\ x_{ij}=u_{ij},\\
l_{ij}<x_{ij}<u_{ij} &\ \Rightarrow\ \bar c_{ij}=0.
\end{aligned}
$$ (eq-04-04-cs)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

在 $\vl=\vzero$ 时取 $y_{ij}=\max\{0,-\bar c_{ij}\}$，这是使 $(\vpi,\vy)$ 对偶可行的最小 $y_{ij}$。对问题对 {eq}`eq-04-04-mcf`–{eq}`eq-04-04-dual`，互补松弛性要求 $x_{ij}>0\Rightarrow \pi_i-\pi_j-y_{ij}=c_{ij}$ 以及 $y_{ij}>0\Rightarrow x_{ij}=u_{ij}$，在这样取的 $\vy$ 下，它们恰好就是 {eq}`eq-04-04-cs` 中的三个蕴含关系。由线性规划对偶理论（{prf:ref}`thm-03-05-complementary-slackness`），原问题的可行解是最优的，当且仅当存在与之互补的对偶解。当 $\vl\ne\vzero$ 时，把这一论证用于平移后的问题即可；检验数在平移下不变。
:::

这些条件有经济上的解释。把 $\pi_i$ 看作一单位产品在节点 $i$ 处的价值（在拥有产品的工厂处高，在需要产品的市场处低）。从 $i$ 运到 $j$，就是以费用 $c_{ij}$ 把一个价值 $\pi_i$ 的单位运到价值为 $\pi_j$ 的地方。若 $\bar c_{ij}>0$，这次运输得不偿失，该弧的使用量应尽可能小；若 $\bar c_{ij}<0$，运输有利可图，该弧应运满；介于两者之间的弧恰好盈亏平衡。

### 负圈

{doc}`04-03-maximum-flow`中的残量网络可以推广到带费用的情形。对满足 $\vl\le\vx\le\vu$ 的流 $\vx$，$G(\vx)$ 包含

- 当 $x_{ij}<u_{ij}$ 时，一条残余容量为 $u_{ij}-x_{ij}$、费用为 $c_{ij}$ 的前向弧 $(i,j)$；
- 当 $x_{ij}>l_{ij}$ 时，一条残余容量为 $x_{ij}-l_{ij}$、费用为 $-c_{ij}$ 的后向弧 $(j,i)$（把流量推回去，每单位节省 $c_{ij}$）。

$G(\vx)$ 中的一个有向圈 $\mathcal W$ 可以在不改变任何净供应量的情况下运送 $\Delta=\min$（$\mathcal W$ 上的残余容量）$>0$ 个单位，总费用的改变量为 $\Delta\,c(\mathcal W)$，其中 $c(\mathcal W)$ 是沿圈各弧费用之和。

:::{prf:theorem} 负圈最优性条件
:label: thm-04-04-negative-cycle

可行流 $\vx$ 是最优的，当且仅当它的残量网络 $G(\vx)$ 中不含费用为负的有向圈。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

若 $G(\vx)$ 中有负圈，沿它推送流量就得到一个费用更低的可行流，所以 $\vx$ 不是最优的。反之，设不存在负圈。添加一个节点 $r$，从它向每个节点引一条费用为零的弧；于是 $G(\vx)$ 中从 $r$ 出发的最短路长度 $d_i$ 是良定义的，并且对每条残量弧满足 $d_j\le d_i+(\text{弧 }(i,j)\text{ 的费用})$（{prf:ref}`thm-04-02-optimality`）。令 $\pi_i=-d_i$。对前向残量弧（$x_{ij}<u_{ij}$），得 $c_{ij}-\pi_i+\pi_j=c_{ij}+d_i-d_j\ge0$，故 $\bar c_{ij}\ge0$；对后向残量弧（$x_{ij}>l_{ij}$），得 $-c_{ij}+d_j-d_i\ge0$，故 $\bar c_{ij}\le0$。因此 $\bar c_{ij}>0$ 迫使 $x_{ij}=l_{ij}$，$\bar c_{ij}<0$ 迫使 $x_{ij}=u_{ij}$，而严格介于上下界之间的弧满足 $\bar c_{ij}=0$：这正是条件 {eq}`eq-04-04-cs`，由{prf:ref}`thm-04-04-reduced-cost`知 $\vx$ 是最优的。
:::

(sec-04-04-supply-chain-solved)=
### 供应链问题的求解

修正后的供应链问题的最优流见{numref}`tab-04-04-supply-chain-solution`，其费用为

$$
\begin{aligned}
z^*&=1\cdot80+4\cdot20+3\cdot0+5\cdot50+6\cdot150+1\cdot120+2\cdot0+2\cdot50+4\cdot70\\
&=\$1\,810.
\end{aligned}
$$

表中的节点势（归一化为 $\pi_6=0$）对节点 1 到 6 依次为 $\pi=(9,10,5,2,4,0)$。每条弧都满足 {eq}`eq-04-04-cs`：严格位于上下界之间的五条弧有 $\bar c_{ij}=0$；$(1,2)$ 和 $(3,5)$ 处于下界 0，且 $\bar c_{ij}>0$；$(1,4)$ 和 $(3,4)$ 处于上界，且 $\bar c_{ij}<0$。所以这个流是最优的；又因为没有一条处于界上的弧的检验数为 0，它是唯一的最优流。正如{prf:ref}`thm-04-04-integrality`所预言，各流量都是整数。平衡条件的验证是直接的；例如，节点 3 收到 $20+50$，送出 $120+0$，净供应量为 $50$。

```{list-table} 修正后供应链的最优流，节点势 $\pi=(9,10,5,2,4,0)$
:name: tab-04-04-supply-chain-solution
:header-rows: 1
:stub-columns: 1

* - 弧
  - $c_{ij}$（美元）
  - $[l_{ij},u_{ij}]$
  - $x^*_{ij}$
  - $\bar c_{ij}=c_{ij}-\pi_i+\pi_j$
  - 状态
* - $(1,2)$
  - 3
  - $[0,\infty)$
  - 0
  - 4
  - 处于下界
* - $(1,3)$
  - 4
  - $[0,\infty)$
  - 20
  - 0
  - 基变量
* - $(1,4)$
  - 1
  - $[50,80]$
  - 80
  - $-6$
  - 处于上界
* - $(2,3)$
  - 5
  - $[0,\infty)$
  - 50
  - 0
  - 基变量
* - $(2,5)$
  - 6
  - $[0,\infty)$
  - 150
  - 0
  - 基变量
* - $(3,4)$
  - 1
  - $[70,120]$
  - 120
  - $-2$
  - 处于上界
* - $(3,5)$
  - 2
  - $[0,\infty)$
  - 0
  - 1
  - 处于下界
* - $(4,6)$
  - 2
  - $[0,120]$
  - 50
  - 0
  - 基变量
* - $(5,6)$
  - 4
  - $[0,\infty)$
  - 70
  - 0
  - 基变量
```

```{figure} figures/04-04-supply-chain.svg
:name: fig-04-04-supply-chain
:width: 600px
:class: dark-light
:alt: 供应链网络，标出供应量、需求量、弧费用和界；右图为最优流，五条基弧突出显示

(a) {prf:ref}`ex-04-04-supply-chain` 中的供应链：方括号内为净供应量，有界弧的标注为“费用（下界，上界）”。弧 $(4,6)$ 上的下界 100 使数据不可行，因而被去掉。(b) 最优流 $x^*_{ij}$；突出显示的弧严格位于上下界之间，构成一棵生成树，$(1,4)$ 和 $(3,4)$ 处于上界，$(1,2)$ 和 $(3,5)$ 上没有流量。
```

检验数也是价格。两个起作用的上界对应 $\bar c_{14}=-6$ 和 $\bar c_{34}=-2$：弧 $(1,4)$ 的容量每增加一个单位，费用降低 6 美元；弧 $(3,4)$ 的容量每增加一个单位，费用降低 2 美元。节点势之差 $\pi_1-\pi_6=9$ 是在工厂 1 多生产一个单位并在市场 6 多售出一个单位的边际成本：这个单位沿 $1\to3$ 运输（4 美元），顶替 $2\to3$ 上的一个单位（节省 5 美元），被顶替的单位沿 $2\to5\to6$ 运输（6 美元 + 4 美元），合计 $4-5+6+4=\$9$。按这些改动重新求解线性规划，恰好得到 1 804 美元、1 808 美元和 1 819 美元。

(sec-04-04-ssp)=
## 连续最短路算法

第一个算法的思路是始终保持费用最优，逐步趋向可行。添加一个超级源点 $s$，从它向每个供应节点引一条容量为 $b_i$ 的弧 $(s,i)$，再添加一个超级汇点 $t$，从每个需求节点向它引一条容量为 $-b_j$ 的弧 $(j,t)$，这些弧的费用都为 0，并假设 $\vl=\vzero$。从零流出发（当所有 $c_{ij}\ge0$ 时零流没有负圈），反复沿一条*最便宜的*增广路推送流量。

:::{prf:algorithm} 连续最短路算法
:label: alg-04-04-ssp

**输入：**一个平衡的问题 {eq}`eq-04-04-mcf`，其中 $\vl=\vzero$，$\vc\ge\vzero$，并按上述方式添加了 $s$ 和 $t$。
**输出：**一个最小费用流，或者问题不可行的结论。

1. **初始化。** $\vx=\vzero$，$z=0$。
2. **最短路。** 在残量网络 $G(\vx)$（前向弧费用为 $c_{ij}$，后向弧费用为 $-c_{ij}$）中找一条从 $s$ 到 $t$ 的最短有向路 $\mathcal P$。若不存在，转第 4 步。
3. **增广。** 令 $\Delta$ 为 $\mathcal P$ 上最小的残余容量。沿 $\mathcal P$ 推送 $\Delta$ 个单位（增加前向弧上的流量，抵消后向弧上的流量），把 $\Delta\cdot c(\mathcal P)$ 加到 $z$ 上，转第 2 步。
4. **停止。** 若所有弧 $(s,i)$ 都已饱和，则 $\vx$ 限制在原有弧上即为最优解，费用为 $z$；否则问题不可行。
:::

**为什么算法有效。** 设 $G(\vx)$ 没有负圈，$\mathcal P$ 是一条最短 $s$–$t$ 路。增广之后新出现的残量弧是 $\mathcal P$ 上各弧的反向弧。记增广前的最短路长度为 $d$，则 $\mathcal P$ 上的每条弧满足 $d_j=d_i+c_{ij}$，所以它的反向弧的费用 $-c_{ij}=d_i-d_j$，于是增广后的每条残量弧仍满足 $c_{ij}+d_i-d_j\ge0$。沿任意圈求和可知没有费用为负的圈：最优性（{prf:ref}`thm-04-04-negative-cycle`）得以保持。在整数数据下，每次增广至少运送一个单位，所以最多需要 $\sum_{b_i>0}b_i$ 次增广。由于反向弧的费用为负，第 2 步需要允许负费用的方法（Bellman–Ford 算法，{ref}`sec-05-04-successive-approximation`）；另一种做法是用距离 $d$ 定义节点势 $\pi=-d$，使所有检验数非负，再在检验数上运行 Dijkstra 算法。

**运输问题的例子。** {numref}`tab-04-04-ssp`列出了{prf:ref}`ex-04-04-transportation`的迭代过程。

```{list-table} 运输问题例子的连续最短路算法
:name: tab-04-04-ssp
:header-rows: 1

* - 迭代
  - 最短路
  - 路的费用（美元/吨）
  - $\Delta$（吨）
  - 累计费用（美元）
* - 1
  - $s\to B\to1\to t$
  - 2
  - 20
  - 40
* - 2
  - $s\to B\to3\to t$
  - 5
  - 30
  - 190
* - 3
  - $s\to A\to2\to t$
  - 6
  - 25
  - 340
* - 4
  - $s\to A\to1\to B\to3\to t$
  - $4-2+5=7$
  - 15
  - 445
```

1. 所有路线中最便宜的是 B → 1，每吨 2 美元；市场 1 需要 20 吨，所以 $\Delta=20$。
2. 市场 1 已满足。剩下最便宜的路线是 B → 3，每吨 5 美元；B 还剩 30 吨，所以 $\Delta=30$，工厂 B 的产品运完。
3. 只有工厂 A 还有产品。它运往市场 2 每吨 6 美元，运往市场 3 每吨 9 美元，经残量路 $A\to1\to B\to3$ 运往市场 3 每吨 7 美元。最便宜的是 $A\to2$；$\Delta=25$ 满足了市场 2。
4. 市场 3 还需要 15 吨。直达弧 $A\to3$ 每吨 9 美元，而路 $A\to1\to B\to3$ 每吨 $4-2+5=\$7$：工厂 A 从工厂 B 手中接过市场 1 的 15 吨（后向弧 $1\to B$ 抵消了 $x_{B1}$ 中的 15 吨），B 把这 15 吨运往市场 3。$\Delta=\min\{15,20,15\}=15$。

运送 90 吨后所有供应弧都已饱和。最优方案为 $x_{A1}=15$，$x_{A2}=25$，$x_{B1}=5$，$x_{B3}=45$，$x_{A3}=x_{B2}=0$，每周费用为 $4\cdot15+6\cdot25+2\cdot5+5\cdot45=\$445$，即表中最后的累计费用。节点势 $\pi_A=7$，$\pi_B=5$，$\pi_1=3$，$\pi_2=1$，$\pi_3=0$ 使四条使用的弧的检验数为 0，未使用的弧的检验数为 $\bar c_{A3}=9-7+0=2$，$\bar c_{B2}=7-5+1=3$，这就证明了最优性（{prf:ref}`thm-04-04-reduced-cost`）；两者都为正，所以方案唯一。第 4 次迭代说明了为什么需要残量网络：从不修改已有运量的贪心规则会选用 $A\to3$，多付 30 美元。

```{figure} figures/04-04-transportation.svg
:name: fig-04-04-transportation
:width: 420px
:class: dark-light
:alt: 二部网络，左侧为工厂 A 和 B，右侧为市场 1 到 3，弧上标有单位费用，最优方案的四条弧突出显示

{prf:ref}`ex-04-04-transportation` 的运输网络，标出单位费用（美元/吨），方括号内为净供应量。突出显示的弧在最优方案中有流量：$A\to1$ 15 吨，$A\to2$ 25 吨，$B\to1$ 5 吨，$B\to3$ 45 吨。
```

(sec-04-04-cycle-cancelling)=
## 消圈算法

第二个算法的思路相反：从任意一个可行流出发，不断降低费用，直到{prf:ref}`thm-04-04-negative-cycle`证明其最优性。

:::{prf:algorithm} 消圈算法
:label: alg-04-04-cycle-cancelling

**输入：**一个最小费用流问题 {eq}`eq-04-04-mcf`。
**输出：**一个最小费用流。

1. **可行流。** 求一个可行流 $\vx$（例如用一次最大流计算，见{ref}`sec-04-04-feasibility`）。若不存在，停止。
2. **检测。** 在 $G(\vx)$ 中搜索满足 $c(\mathcal W)<0$ 的有向圈 $\mathcal W$（Bellman–Ford 算法或 Floyd–Warshall 算法都能检测负圈）。若不存在，停止：$\vx$ 是最优的。
3. **消圈。** 令 $\Delta$ 为 $\mathcal W$ 上最小的残余容量；沿 $\mathcal W$ 推送 $\Delta$ 个单位，费用降低 $\Delta\abs{c(\mathcal W)}$。转第 2 步。
:::

在整数数据下，每次消圈至少使费用降低 1，所以算法必然终止。若每次选取每条弧的*平均*费用最负的圈，算法就成为多项式时间算法（Goldberg & Tarjan, 1989）。消圈算法也解释了上面第 4 次迭代中的“推回”：如果在该次迭代之前的方案上再在 $A\to3$ 上运送 15 吨，就会留下残量圈 $A\to1\to B\to3\to A$，其费用为 $4-2+5-9=-2$，消去这个圈就得到最优解。完整的求解过程见{prf:ref}`ex-08-05-cycle-cancelling`。

(sec-04-04-network-simplex)=
## 网络单纯形法（概要）

把单纯形法用于 {eq}`eq-04-04-mcf` 时，可以不用单纯形表而直接在网络上进行，因为基有图论上的含义。

:::{prf:theorem} 基就是生成树
:label: thm-04-04-spanning-tree-basis

对连通网络，$\mA$（删去一个冗余行后）的 $\abs{\mathcal N}-1$ 列构成一个基，当且仅当对应的弧（不计方向）构成一棵生成树。
:::

理由是：一个圈上各弧对应的列线性相关（按行进方向取符号 $\pm1$ 相加，它们恰好抵消），所以基中不含圈；而 $\abs{\mathcal N}-1$ 条不含圈的弧构成一棵生成树（{prf:ref}`thm-04-01-tree`）。因此在基本解中，每条非树弧都处于其下界或上界，树弧上的流量由平衡方程确定，可以从叶节点向内逐步算出。**网络单纯形法**（network simplex method）的一次迭代如下：

1. **节点势。** 在根节点固定 $\pi_r=0$，从根向外沿树弧求解 $\bar c_{ij}=c_{ij}-\pi_i+\pi_j=0$。
2. **定价。** 计算非树弧的 $\bar c_{ij}$。若每条处于下界的弧都满足 $\bar c_{ij}\ge0$，每条处于上界的弧都满足 $\bar c_{ij}\le0$，则条件 {eq}`eq-04-04-cs` 成立：停止。
3. **进基弧。** 否则选一条违反条件的弧。它与树一起恰好构成一个圈。
4. **比值检验与主元变换。** 沿使费用下降的方向在这个圈上推送流量，直到圈上某条弧达到其界；这条弧离开树，进基弧加入树，并更新节点势。

所有运算都只是沿树上路径的加法和减法，这就是该方法在实践中速度快、并且保持整数性的原因。退化的主元变换（步长为零）可能导致循环；Cunningham（1976）提出的“强可行树”规则可以防止循环。对供应链问题，最优树是{numref}`fig-04-04-supply-chain`(b) 中突出显示的弧集 $\{(1,3),(2,3),(2,5), (5,6),(4,6)\}$；取 $\pi_6=0$，第 1 步得到 $\pi_5=4$，$\pi_2=10$，$\pi_3=5$，$\pi_1=9$，$\pi_4=2$，第 2 步立即停止，正如{numref}`tab-04-04-supply-chain-solution`所示。

## 本章回顾

### 小结

- 最小费用流问题 {eq}`eq-04-04-mcf` 在流量平衡 $\mA\vx=\vb$ 和界 $\vl\le\vx\le\vu$ 的约束下最小化 $\sum c_{ij}x_{ij}$（{prf:ref}`def-04-04-mcf`）。
- 它要求 $\sum_ib_i=0$（否则添加虚拟节点），并要求每个节点集上有足够的容量 {eq}`eq-04-04-gale`；{prf:ref}`ex-04-04-supply-chain`中的供应链在节点 4 处不满足这一条件，去掉一个下界后得到修正。流量下界可以通过 {eq}`eq-04-04-shift`–{eq}`eq-04-04-shifted-supply` 去掉。
- 关联矩阵是全单模的（{prf:ref}`thm-04-04-tu`），所以整数数据给出整数的最优基本流（{prf:ref}`thm-04-04-integrality`）。
- 最短路、最大流、运输、指派和转运问题都是它的特例（{ref}`sec-04-04-special-cases`）。
- 一个流是最优的，当且仅当存在满足 {eq}`eq-04-04-cs` 的节点势（{prf:ref}`thm-04-04-reduced-cost`），等价地，当且仅当它的残量网络中没有负圈（{prf:ref}`thm-04-04-negative-cycle`）。供应链问题的最优费用为 1 810 美元（{numref}`tab-04-04-supply-chain-solution`）。
- 连续最短路算法（{prf:ref}`alg-04-04-ssp`）经四次增广求得运输问题例子的最优费用 445 美元（{numref}`tab-04-04-ssp`）；消圈算法（{prf:ref}`alg-04-04-cycle-cancelling`）和网络单纯形法（基即生成树，{prf:ref}`thm-04-04-spanning-tree-basis`）是另外两种标准方法。

### 常见错误

- **数据不平衡。** 若 $\sum b_i\ne0$，无论网络如何，线性规划都不可行；应添加虚拟节点，并确定其弧的费用。
- **在供应量中忘记流量下界。** 作替换 $x_{ij}=l_{ij}+x'_{ij}$ 后，弧的两个端点的净供应量都要改变 {eq}`eq-04-04-shifted-supply`，而且常数 $\sum c_{ij}l_{ij}$ 必须加回到费用中。
- **检验数的符号。** 在 $\bar c_{ij}=c_{ij}-\pi_i+\pi_j$ 下，处于*上*界的弧在 $\bar c_{ij}\le0$（而不是 $\ge0$）时才满足最优性。
- **忽略残量网络中的后向弧。** 后向弧的费用为负，所以不借助节点势就不能对它们使用 Dijkstra 算法；而没有后向弧，连续最短路算法可能停在非最优的方案上。
- **以为整数性是自动成立的。** 它要求 $\vb$、$\vl$ 和 $\vu$ 为整数，并且最优解是*基本*解；内点法求解器可能返回最优面上的一个分数点。

## 参考文献

- Ahuja, R. K., Magnanti, T. L., & Orlin, J. B. (1993). *Network flows:
  Theory, algorithms, and applications*. Prentice Hall. (Chapter 9:
  optimality conditions, duality, cycle-canceling and successive
  shortest path algorithms; Chapter 11: network simplex algorithm; Chapter
  12: assignment problems.)
- Bazaraa, M. S., Jarvis, J. J., & Sherali, H. D. (2010). *Linear programming
  and network flows* (4th ed.). Wiley. (Chapter 9: minimal-cost network flow
  problems, total unimodularity and the network simplex method; Chapter 10:
  the transportation and assignment problems.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Chapter 7: network flow problems, the
  network simplex algorithm and integrality.)
- Cunningham, W. H. (1976). A network simplex method. *Mathematical
  Programming, 11*(1), 105–116.
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Chapters 14–17 and 19–20: the transportation, assignment
  and transshipment problems, networks, maximal flows and the primal–dual
  method for transportation problems.)
- Gale, D. (1957). A theorem on flows in networks. *Pacific Journal of
  Mathematics, 7*(2), 1073–1082.
- Goldberg, A. V., & Tarjan, R. E. (1989). Finding minimum-cost circulations
  by canceling negative cycles. *Journal of the ACM, 36*(4), 873–886.
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 9: transportation and
  assignment problems; Sections 10.6–10.7: minimum cost flow problem and the
  network simplex method.)
- Hoffman, A. J. (1960). Some recent applications of the theory of linear
  inequalities to extremal combinatorial analysis. In R. Bellman & M. Hall
  (Eds.), *Combinatorial analysis* (pp. 113–127). American Mathematical
  Society.
- Kuhn, H. W. (1955). The Hungarian method for the assignment problem.
  *Naval Research Logistics Quarterly, 2*(1–2), 83–97.
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 5: transportation model and its variants; Chapter 6:
  network models.)
