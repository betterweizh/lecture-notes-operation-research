# 网络与最小生成树

## 概述

许多基础设施系统都由地点和地点之间的连接构成：道路交叉口与路段、泵站与管道、交换机与电缆。**网络**（network）记录的正是这种结构，一大类最优化模型（最短路、最大流、最小费用流）都建立在网络之上。本章确定网络各章所使用的术语（节点、弧、路、圈、连通性、树），证明关于树的基本事实（$\abs{\mathcal N}$ 个节点上的树有 $\abs{\mathcal N}-1$ 条弧，增加或删去一条弧的效果是可以预见的），并说明网络的存储方式。然后求解第一个网络最优化问题：以最小的总费用把所有节点连接起来，即**最小生成树**（minimum spanning tree, MST）问题。割性质与圈性质这两个最优性结果解释了为什么在这里贪婪的选择是最优的；由它们可以推出 **Prim 算法**（从一个起始节点出发生长出一棵树）和 **Kruskal 算法**（每次加入最便宜的连接）的正确性。两种算法都在一个有七个站点的污水网络上完整地执行一遍。这些算法是组合性的，不需要最优化理论；只有{ref}`sec-04-01-representations`用到{doc}`01-01-linear-algebra`中的矩阵与秩。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 定义有向网络与无向网络、路、链、圈、连通性、树和生成树，并能在图中辨认它们；
- 证明树的各种等价刻画，并运用生成树的交换性质；
- 将网络存储为弧表、节点–弧关联矩阵、邻接矩阵或邻接表，并在它们之间作出选择；
- 叙述并证明割性质、圈性质以及最小生成树的最优性条件；
- 手工执行 Prim 算法和 Kruskal 算法，论证其正确性，并给出其运行时间；
- 判断最小生成树是否唯一，并在不唯一时列出其他最优树。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $G=(\mathcal N,\mathcal A)$ | 节点集为 $\mathcal N$、弧（边）集为 $\mathcal A$ 的网络 |
| $(i,j)$ | 弧尾为 $i$、弧头为 $j$ 的弧；在无向网络中指连接 $i$ 与 $j$ 的边，也记作 $i$–$j$ |
| $\abs{\mathcal N}$, $\abs{\mathcal A}$ | 节点数，弧数 |
| $c_{ij}$ | 弧 $(i,j)$ 的长度或费用 |
| $\deg(i)$ | 节点 $i$ 的度（与 $i$ 关联的边数） |
| $\mathcal T$, $c(\mathcal T)$ | 一棵（生成）树（以其边集表示）及其总费用 $\sum_{(i,j)\in\mathcal T}c_{ij}$ |
| $[S,\bar S]$ | 割：一端在 $S\subset\mathcal N$ 中、另一端在 $\bar S=\mathcal N\setminus S$ 中的边 |
| $\mathcal W$ | 一个圈 |
| $C_k$, $\bar C_k$ | Prim 算法第 $k$ 次迭代中已连通与未连通的节点集 |
| $[\gamma(j),P(j)]$ | Prim 算法中未连通节点的标号：它到 $C_k$ 的最便宜边的费用，以及该边在 $C_k$ 中的另一端节点 |
| $\mA$ | 节点–弧关联矩阵（网络流模型的约束矩阵） |

(sec-04-01-graphs)=
## 无向图与有向图

:::{prf:definition} 网络
:label: def-04-01-network

**网络**（network，或称图，graph）$G=(\mathcal N,\mathcal A)$ 由有限的**节点**（node，或称顶点，vertex）集 $\mathcal N$ 和**弧**（arc，或称连接，link）集 $\mathcal A$ 组成，每条弧连接两个不同的节点。

- 在**有向网络**（directed network，又称有向图，digraph）中，每条弧都是一个有序对 $(i,j)$，其**弧尾**（tail）为 $i$，**弧头**（head）为 $j$；它只能从 $i$ 走向 $j$。
- 在**无向网络**（undirected network）中，每条弧都是一个无序对，称为**边**（edge）；边 $(i,j)$ 与 $(j,i)$ 相同，两个方向都可以通行。

附加在弧上的数值，如长度、费用或容量，称为弧的**属性**（attribute）。弧 $(i,j)$ 的长度或费用记作 $c_{ij}$。
:::

在网络各章中，从 $i$ 到 $j$ 至多有一条弧，且没有从某节点到其自身的弧，因此一条弧由它的两个端点唯一确定。由一条弧连接的两个节点称为**相邻**（adjacent）的，这条弧与这两个节点都**关联**（incident）。无向网络中节点的**度**（degree）$\deg(i)$ 是与它关联的边数。每条边都有两个端点，因此

$$
\sum_{i\in\mathcal N}\deg(i)=2\abs{\mathcal A}.
$$ (eq-04-01-handshake)

在有向网络中，$i$ 的**出度**（out-degree）是以 $i$ 为弧尾的弧数，**入度**（in-degree）是以 $i$ 为弧头的弧数；对所有节点求和，两者都等于 $\abs{\mathcal A}$。把每条边 $(i,j)$ 换成两条弧 $(i,j)$ 和 $(j,i)$，两者都带有该边的属性，就总能把无向网络变成有向网络；{doc}`04-02-shortest-paths`和{doc}`04-03-maximum-flow`中的最短路算法和流算法都在这种有向形式上运行。

:::{prf:example} 一个小型道路网络
:label: ex-04-01-road-network

五个道路交叉口由八条路段相连（{numref}`fig-04-01-road-network`）。将交叉口编号为 1 到 5，

$$
\begin{aligned}
\mathcal N&=\{1,2,3,4,5\},\\
\mathcal A&=\{(3,1),(3,4),(1,2),(1,4),(1,5),(2,5),(2,4),(4,5)\},
\end{aligned}
$$

各路段的长度 $c_{ij}$（单位：英里）为

| 弧 $(i,j)$ | $(3,1)$ | $(3,4)$ | $(1,2)$ | $(1,4)$ | $(1,5)$ | $(2,5)$ | $(2,4)$ | $(4,5)$ |
|---|---|---|---|---|---|---|---|---|
| $c_{ij}$ | 1 | 2 | 0.5 | 2.5 | 5 | 1.2 | 1.7 | 3 |

求所有道路都是双向通行时（图 (a)）各节点的度，以及每条道路都只能沿其弧的列出方向单向通行时（图 (b)）各节点的出度和入度。用弧数对两种结果加以验证。

*（这一网络的某些画法把节点 2 与 4 之间的弧画成从 4 指向 2；这里始终取为 $(2,4)$，与上面的弧集一致。）*

```{figure} figures/04-01-road-network.svg
:name: fig-04-01-road-network
:width: 700px
:class: dark-light
:alt: 由八条带长度的连接相连的五个节点，分别画成不带箭头和带箭头的形式

{prf:ref}`ex-04-01-road-network` 的道路网络：(a) 无向网络，标出各路段长度（英里）；(b) 有向网络，每条弧都按 $\mathcal A$ 中列出的方向定向。
```
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

作为无向网络（所有道路双向通行，图 (a)），各节点的度为 $\deg(1)=4$、$\deg(2)=3$、$\deg(3)=2$、$\deg(4)=4$ 和 $\deg(5)=3$；其和为 $16=2\times 8$，与 {eq}`eq-04-01-handshake` 的要求一致。

如果每条道路都只能沿其弧的列出方向单向通行（图 (b)），网络就是有向的。节点 1 到 5 的出度依次为 3、2、2、1、0，入度依次为 1、1、0、3、3；两组数之和都等于 $\abs{\mathcal A}=8$。节点 3 没有进入的弧，节点 5 没有离开的弧。
:::

(sec-04-01-paths)=
## 路、圈与连通性

:::{prf:definition} 路、链与圈
:label: def-04-01-path

设 $G=(\mathcal N,\mathcal A)$ 是一个网络。

- 在无向网络中，从 $i$ 到 $j$ 的**路**（path）是一列互不相同的节点 $i=i_1,i_2,\dots,i_p=j$，其中每一对相继节点 $i_k,i_{k+1}$ 都由一条边相连。路上的边互不相同，共有 $p-1$ 条。
- 在有向网络中，具有同样性质、但弧可以沿任一方向使用的节点序列称为**链**（chain，即无向的路）。如果每条弧都沿其自身方向使用，即对所有 $k$ 有 $(i_k,i_{k+1})\in\mathcal A$，就称为**有向路**（directed path）。
- **圈**（cycle）是回到起始节点的路：互不相同的节点 $i_1,\dots,i_p$，$p\ge3$，相继节点之间有边，并且有边 $(i_p,i_1)$。**有向圈**（directed cycle）沿各弧自身的方向使用弧 $(i_1,i_2),\dots, (i_{p-1},i_p),(i_p,i_1)$（$p\ge2$）。没有有向圈的有向网络称为**无圈的**（acyclic）。

路、链或圈的**长度**是其各弧属性之和。
:::

:::{prf:definition} 连通性与子网络
:label: def-04-01-connected

如果任意两个不同的节点之间都至少有一条路相连，则称无向网络是**连通的**（connected）。如果有向网络的无向形式是连通的，则称它是（弱）连通的；如果从每个节点到其他每个节点都有有向路，则称它是**强连通的**（strongly connected）。$G$ 的**子网络**（sub-network）是满足 $\mathcal N'\subseteq\mathcal N$、$\mathcal A'\subseteq\mathcal A$ 且 $\mathcal A'$ 中每条弧的两个端点都在 $\mathcal N'$ 中的网络 $G'=(\mathcal N',\mathcal A')$；若 $\mathcal N'=\mathcal N$，则称它是**生成的**（spanning）。$G$ 的极大连通子网络称为它的**连通分支**（component）。
:::

把{prf:ref}`ex-04-01-road-network`的道路网络看作无向网络（{numref}`fig-04-01-terminology`）：

- 3–1–5 是从节点 3 到节点 5 的一条路，长度为 $1+5=6$ 英里（图 a）；
- 1–2–4–1 是一个圈，长度为 $0.5+1.7+2.5=4.7$ 英里（图 b）。路线 3–1–2–4–1–5 连接节点 3 和 5，但两次经过节点 1，因此不是路；它包含圈 1–2–4–1，去掉这个圈后剩下路 3–1–5；
- 该网络是连通的，但删去节点 5 处三条边后得到的生成子网络不连通：节点 5 自成一个连通分支（图 c）。

现在把该网络看作有向网络（{numref}`fig-04-01-road-network`的图 (b)）。序列 3、1、2、5 是一条有向路，长度为 $1+0.5+1.2=2.7$ 英里（图 f）；从 3 到 5 还有另外四条有向路：3→1→2→4→5、3→1→4→5、3→1→5 和 3→4→5。圈 1–2–4–1 不是有向圈，因为封闭它需要弧 $(4,1)$，而 $\mathcal A$ 中有的是 $(1,4)$；节点、边和弧都相同，只是方向不同。事实上这个有向网络是无圈的：按 3、1、2、4、5 的顺序，每条弧都从较前的节点指向较后的节点，因此任何有向路都不可能回到起点。它是连通的，但不是强连通的，因为没有弧进入节点 3，也没有弧离开节点 5。

```{figure} figures/04-01-terminology.svg
:name: fig-04-01-terminology
:width: 700px
:class: dark-light
:alt: 五节点网络的六个副本，每个副本突出显示一个不同的子网络

道路网络上的术语；突出显示的边和节点构成所讨论的对象，灰色虚线边是网络的其余部分。(a) 路 3–1–5。(b) 圈 1–2–4–1。(c) 一个不连通的生成子网络（节点 5 是孤立的）。(d) 一棵不能生成整个网络的树（缺少节点 5）。(e) 一棵生成树。(f) 有向形式中的有向路 3→1→2→5。
```

(sec-04-01-trees)=
## 树

:::{prf:definition} 树、森林与生成树
:label: def-04-01-tree

没有圈的无向网络称为**森林**（forest），连通的森林称为**树**（tree）。网络 $G$ 的**生成树**（spanning tree）是 $G$ 的一个本身是树的生成子网络。树中度为 1 的节点称为**叶节点**（leaf）。树记作 $\mathcal T$，并与其边集等同看待。
:::

在{numref}`fig-04-01-terminology`中，图 (d) 是节点 1 到 4 上的一棵树，但不是道路网络的生成树，而图 (e) 是一棵生成树。树是使一组节点保持连通的最便宜的方式，下面的结果完整地刻画了树。

:::{prf:lemma} 叶节点
:label: lem-04-01-leaves

至少有两个节点的树至少有两个叶节点。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

有两个或更多节点的树是连通的，因此有边。取一条边数最多的路 $i_1,\dots,i_p$（$p\ge2$；由于网络有限，这样的路存在）。节点 $i_1$ 与 $i_2$ 相邻。它在路外没有相邻节点，因为这样的相邻节点可以延长这条路；它也没有满足 $k\ge3$ 的相邻节点 $i_k$，因为边 $(i_1,i_k)$ 与路 $i_1,\dots,i_k$ 会构成一个圈。因此 $\deg(i_1)=1$，同理 $\deg(i_p)=1$。
:::

:::{prf:theorem} 树的刻画
:label: thm-04-01-tree

对无向网络 $G=(\mathcal N,\mathcal A)$，以下各条等价：

1. $G$ 是树（连通且无圈）；
2. $G$ 连通且 $\abs{\mathcal A}=\abs{\mathcal N}-1$；
3. $G$ 无圈且 $\abs{\mathcal A}=\abs{\mathcal N}-1$；
4. $G$ 的任意两个节点之间恰有一条路相连。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

*树有 $\abs{\mathcal N}-1$ 条边。* 对 $\abs{\mathcal N}$ 用归纳法。只有一个节点的树没有边。若 $\abs{\mathcal N}\ge2$，删去一个叶节点及其边（{prf:ref}`lem-04-01-leaves`）。剩下的部分没有圈，而且是连通的，因为叶节点不可能是一条路的内部节点（内部节点的度至少为 2）。由归纳假设，它有 $\abs{\mathcal N}-2$ 条边，所以原来的树有 $\abs{\mathcal N}-1$ 条边。这就证明了 1 ⇒ 2 和 1 ⇒ 3。将此结论应用于每个连通分支，还可知有 $q$ 个连通分支的森林有 $\abs{\mathcal N}-q$ 条边。

*3 ⇒ 1。* 无圈的网络是森林；若它有 $q$ 个连通分支，则有 $\abs{\mathcal N}-q$ 条边。于是由 $\abs{\mathcal N}-q=\abs{\mathcal N}-1$ 得 $q=1$：网络是连通的。

*2 ⇒ 1。* 若 $G$ 有圈，删去圈上的一条边后网络仍连通（圈的其余部分仍连接被删边的两个端点）。重复这一步骤直到不再有圈，就得到一棵边数少于 $\abs{\mathcal N}-1$ 的树，这与 1 ⇒ 2 矛盾。

*1 ⇒ 4。* 树是连通的，因此至少有一条路。假设有两条不同的路 $P$ 和 $Q$ 连接 $i$ 与 $j$。从 $i$ 出发，设 $u$ 是两条路分开之前的最后一个节点，$w$ 是 $P$ 上在 $u$ 之后第一个也位于 $Q$ 上的节点（由于 $j$ 同时在两条路上，$w$ 存在）。$P$ 从 $u$ 到 $w$ 的部分与 $Q$ 从 $w$ 回到 $u$ 的部分只有 $u$ 和 $w$ 两个公共节点，而且它们不是同一条边，因此合起来构成一个圈，矛盾。

*4 ⇒ 1。* 网络是连通的。如果有圈，则圈上两个相继节点之间会有两条路（它们之间的边以及圈的其余部分），这与唯一性矛盾。
:::

这个定理有两个推论，它们是所有生成树算法的基础。

:::{prf:corollary} 生成树的交换性质
:label: cor-04-01-exchange

设 $\mathcal T$ 是连通网络 $G$ 的一棵生成树。

1. 把边 $(k,l)\notin\mathcal T$ 加入 $\mathcal T$，恰好产生一个圈：边 $(k,l)$ 连同树中从 $k$ 到 $l$ 的唯一一条路。删去这个圈上的任意一条边，得到另一棵生成树。
2. 删去边 $(k,l)\in\mathcal T$，把 $\mathcal T$ 恰好分成两棵树，其节点集分别为 $S\ni k$ 和 $\bar S\ni l$。加入 $G$ 中任意一条一端在 $S$、另一端在 $\bar S$ 中的边，得到另一棵生成树。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

1. 经过新边的圈由 $(k,l)$ 和 $\mathcal T$ 中一条从 $l$ 回到 $k$ 的路组成，而这样的路只有一条（{prf:ref}`thm-04-01-tree`第 4 条）。删去这个圈上的一条边，剩下一个有 $\abs{\mathcal N}-1$ 条边的连通网络，由第 2 条，它是一棵树。
2. 删边之后有 $\abs{\mathcal N}-2$ 条边且没有圈，因此这个森林有 $\abs{\mathcal N}-(\abs{\mathcal N}-2)=2$ 个连通分支。它们之间的一条边重新把网络连通起来，边数为 $\abs{\mathcal N}-1$。
:::

每个连通网络都有生成树：逐条删去圈上的边，直到不再有圈（即 2 ⇒ 1 的论证）。生成树的数目通常非常多。由**基尔霍夫矩阵树定理**（Kirchhoff's matrix-tree theorem），生成树的数目等于拉普拉斯矩阵（对角线上为各节点的度，每对相邻节点对应的位置为 $-1$）的任一代数余子式；对道路网络，由此得到 40 棵生成树。对 $\abs{\mathcal N}$ 个节点上的完全网络（每对节点都相连），这个数目为 $\abs{\mathcal N}^{\abs{\mathcal N}-2}$（Cayley 公式），20 个节点时已约为 $2.6\times10^{23}$，因此生成树问题不可能靠逐一列举生成树来求解。

(sec-04-01-representations)=
## 网络的表示

图形便于人阅读；算法则需要以下某种数据结构。

**弧表**（arc list）。各条弧以三元组 $(i,j,c_{ij})$ 的形式列出。所需存储量与 $\abs{\mathcal A}$ 成正比，这是常用的输入格式，也是 Kruskal 算法进行排序的对象。

**节点–弧关联矩阵**（node–arc incidence matrix）。矩阵 $\mA\in\R^{\abs{\mathcal N}\times \abs{\mathcal A}}$ 每个节点对应一行，每条弧对应一列。弧 $(i,j)$ 对应的列为 $\ve_i-\ve_j$：在第 $i$ 行（弧尾）为 $+1$，在第 $j$ 行（弧头）为 $-1$，其余为零。{numref}`tab-04-01-incidence`给出了有向道路网络的关联矩阵。第 $i$ 行中 $+1$ 的个数等于 $i$ 的出度，$-1$ 的个数等于其入度。

```{list-table} 有向道路网络的节点–弧关联矩阵
:name: tab-04-01-incidence
:header-rows: 1
:stub-columns: 1

* - 节点
  - $(1,2)$
  - $(1,4)$
  - $(1,5)$
  - $(2,4)$
  - $(2,5)$
  - $(3,1)$
  - $(3,4)$
  - $(4,5)$
* - 1
  - 1
  - 1
  - 1
  - 0
  - 0
  - −1
  - 0
  - 0
* - 2
  - −1
  - 0
  - 0
  - 1
  - 1
  - 0
  - 0
  - 0
* - 3
  - 0
  - 0
  - 0
  - 0
  - 0
  - 1
  - 1
  - 0
* - 4
  - 0
  - −1
  - 0
  - −1
  - 0
  - 0
  - −1
  - 1
* - 5
  - 0
  - 0
  - −1
  - 0
  - −1
  - 0
  - 0
  - −1
```

每一列恰有一个 $+1$ 和一个 $-1$，因此各行之和为零，$\vone^\T\mA=\vzero^\T$，从而 $\rank\mA\le\abs{\mathcal N}-1$。对连通网络，秩恰为 $\abs{\mathcal N}-1$（这里为 4）。该矩阵是每个网络流模型中流量守恒方程 $\mA\vx=\vb$ 的约束矩阵（{doc}`04-02-shortest-paths`，{doc}`04-03-maximum-flow`，{doc}`04-04-minimum-cost-flow`），它的特殊结构正是这些线性规划问题具有整数最优解的原因（{prf:ref}`thm-04-04-tu`，{prf:ref}`thm-04-04-integrality`）。但作为存储方式它很浪费：其 $\abs{\mathcal N}\abs{\mathcal A}$ 个元素中只有 $2\abs{\mathcal A}$ 个非零。对无向网络，每列的两个非零元素都是 $+1$。

**邻接矩阵**（adjacency matrix）。这是一个 $\abs{\mathcal N}\times\abs{\mathcal N}$ 矩阵，当 $(i,j)\in\mathcal A$ 时位置 $(i,j)$ 上的元素为 1，否则为 0；对无向网络它是对称的。它的赋权形式在弧对应的位置上存放 $c_{ij}$，对角线上为 0，不存在的弧对应 $+\infty$；这种**长度矩阵**（length matrix）是 Floyd–Warshall 算法（{prf:ref}`alg-04-02-floyd`）的出发点。所需存储量为 $\abs{\mathcal N}^2$，检验一条弧是否存在只需常数时间，适合稠密网络。

**邻接表**（adjacency list）。对每个节点 $i$，列出由 $i$ 出发的弧所到达的节点 $\mathcal A(i)=\{j:(i,j)\in \mathcal A\}$，并附上弧的属性。对有向道路网络，

$$
\begin{aligned}
\mathcal A(1)&=\{2,4,5\},\\
\mathcal A(2)&=\{4,5\},\\
\mathcal A(3)&=\{1,4\},\\
\mathcal A(4)&=\{5\},\\
\mathcal A(5)&=\varnothing .
\end{aligned}
$$

所需存储量与 $\abs{\mathcal N}+\abs{\mathcal A}$ 成正比，扫描由 $i$ 出发的弧所需时间与 $\abs{\mathcal A(i)}$ 成正比。若把邻接表按弧尾节点排序存放在一个数组中，并为每个节点设一个指向其第一条弧的指针，就得到**前向星**（forward star）表示。Prim 算法和 Dijkstra 算法要反复扫描某一个节点的弧，因此都用邻接表实现；相应的 Python 数据结构见{doc}`09-04-code-networks`。

(sec-04-01-mst-problem)=
## 最小生成树问题

:::{prf:definition} 最小生成树
:label: def-04-01-mst

设 $G=(\mathcal N,\mathcal A)$ 是边费用为 $c_{ij}$ 的连通无向网络。**最小生成树问题**是求 $G$ 的一棵总费用最小的生成树 $\mathcal T$：

$$
c(\mathcal T)=\sum_{(i,j)\in\mathcal T}c_{ij}.
$$ (eq-04-01-tree-cost)
:::

“以尽可能低的费用连接所有节点”这一任务要求的是费用最低的连通生成子网络；当费用为正时，这样的子网络没有圈（圈上的一条边可以去掉而不破坏连通性），因此它是一棵生成树。另有两点说明以后会用到。第一，每棵生成树都有同样多的边，即 $\abs{\mathcal N}-1$ 条，因此给每个费用加上同一个常数，所有树的费用都改变同样的量，最优树不变；所以负费用并无妨碍。第二，下面的算法只比较费用的大小，因此费用可以是任意实数。

典型的应用是只要求连通性的网络设计：

- **基础设施布局**：从水库出发或通往处理厂的管道、电缆和光纤网络、输电线路、乡村道路，其中一条连接的费用就是它的建设费用；
- **聚类**（clustering）：当费用度量的是差异程度时，删去最小生成树中费用最高的 $q-1$ 条边，就得到 $q$ 个类（单链接聚类，single-linkage clustering）；{prf:ref}`ex-08-05-machine-cells`用这种方法把机器划分为制造单元；
- **近似**：最小生成树是求解更难问题（如旅行商问题）的启发式方法的出发点。

本问题是对无向网络提出的。它的有向版本，即从一个根出发能到达每个节点的费用最低的弧集（最小费用树形图，minimum-cost arborescence），需要不同的算法，本书不作讨论。

:::{prf:example} 一个污水收集网络
:label: ex-04-01-sewage

一个深层隧道污水系统的七个站点必须用隧道连接起来，使每个站点的污水都能到达主处理站，即节点 7。技术上可行的隧道共有十四条（{numref}`fig-04-01-sewage-network`）；没有连接的站点对若直接相连则费用过高。各条隧道的建设费用（单位：百万美元）为

| 连接 | 1–2 | 1–3 | 1–4 | 1–5 | 2–3 | 2–4 | 2–5 | 3–4 | 3–6 | 3–7 | 4–5 | 4–6 | 5–7 | 6–7 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 费用 | 1 | 5 | 7 | 9 | 6 | 4 | 3 | 5 | 10 | 8 | 8 | 3 | 3 | 4 |

应修建哪些隧道，才能以最小的总费用把所有站点连接起来？

```{figure} figures/04-01-sewage-network.svg
:name: fig-04-01-sewage-network
:width: 560px
:class: dark-light
:alt: 由十四条可能的隧道连接相连的七个站点及其费用；节点 7 被特别标出

{prf:ref}`ex-04-01-sewage` 的污水网络。标注为建设费用（单位：百万美元）；节点 7 是主处理站。
```
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

任何连通的布局都能使每个站点到达节点 7，因此费用最低的布局就是一棵最小生成树。它有 $7-1=6$ 条连接，从 14 条中选出。在 $\binom{14}{6}=3003$ 个六条连接的集合中，由矩阵树定理可知有 1 020 个是生成树。其中之一为

$$
\begin{gathered}
\mathcal T=\{1\text{–}2,\ 2\text{–}5,\ 5\text{–}7,\ 1\text{–}3,\ 3\text{–}6,\ 4\text{–}5\},\\
c(\mathcal T)=1+3+3+5+10+8=30,
\end{gathered}
$$

这是一个费用为 3 000 万美元的可行布局。Prim 算法（{prf:ref}`ex-04-01-sewage-prim`）和 Kruskal 算法（{prf:ref}`ex-04-01-sewage-kruskal`）无需列举所有树就能求出最优解。
:::

(sec-04-01-optimality)=
## 最优性条件：割性质与圈性质

对满足 $\varnothing\ne S\ne\mathcal N$ 的节点集 $S$，**割**（cut）$[S,\bar S]$ 是一端在 $S$ 中、另一端在 $\bar S=\mathcal N\setminus S$ 中的边的集合。每棵生成树都必须把 $S$ 与 $\bar S$ 连接起来，因此都至少包含每个割中的一条边。两个简单的交换论证说明了哪些边属于最小生成树。

:::{prf:theorem} 割性质
:label: thm-04-01-cut-property

设 $F$ 是包含于某棵最小生成树中的边集，$[S,\bar S]$ 是不含 $F$ 中任何边的一个割。若 $(k,l)$ 是 $[S,\bar S]$ 中费用最小的边，则 $F\cup\{(k,l)\}$ 包含于某棵最小生成树中。若 $(k,l)$ 是该割中唯一的费用最小的边，则每棵最小生成树都包含它。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

设 $\mathcal T^*$ 是满足 $F\subseteq\mathcal T^*$ 的最小生成树。若 $(k,l)\in\mathcal T^*$，则无需证明。否则，把 $(k,l)$ 加入 $\mathcal T^*$ 就形成一个圈（{prf:ref}`cor-04-01-exchange`第 1 条），它经过 $(k,l)$ 离开 $S$，必须经过该割的另一条边 $(p,q)$ 返回；$(p,q)$ 是树边，并且由于 $F$ 中没有边位于该割中，$(p,q)\notin F$。树 $\mathcal T'=\mathcal T^*\cup\{(k,l)\} \setminus\{(p,q)\}$ 包含 $F\cup\{(k,l)\}$，其费用为

$$
c(\mathcal T')=c(\mathcal T^*)+c_{kl}-c_{pq}\le c(\mathcal T^*),
$$

因为 $(k,l)$ 是该割中最便宜的边。因此 $\mathcal T'$ 也是最小生成树。若 $(k,l)$ 严格比该割中所有其他边都便宜，则对每棵不含 $(k,l)$ 的最小生成树 $\mathcal T^*$（取 $F=\varnothing$），上述不等式都是严格的，这是不可能的；因此每棵最小生成树都包含 $(k,l)$。
:::

:::{prf:theorem} 圈性质
:label: thm-04-01-cycle-property

设 $\mathcal W$ 是 $G$ 中的一个圈，$(k,l)$ 是 $\mathcal W$ 上一条严格比 $\mathcal W$ 上所有其他边都贵的边。则没有任何最小生成树包含 $(k,l)$。若 $(k,l)$ 只是 $\mathcal W$ 上最贵的边之一（允许并列），则存在某棵最小生成树不包含它。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

设 $\mathcal T^*$ 是包含 $(k,l)$ 的最小生成树。删去 $(k,l)$ 把它分成两棵树，其节点集分别为 $S\ni k$ 和 $\bar S\ni l$（{prf:ref}`cor-04-01-exchange`第 2 条）。圈的其余部分，即去掉 $(k,l)$ 的 $\mathcal W$，是一条从 $k$ 到 $l$ 的路，因此它包含一条一端在 $S$、另一端在 $\bar S$ 中的边 $(p,q)\ne(k,l)$。于是 $\mathcal T'=\mathcal T^*\setminus\{(k,l)\}\cup\{(p,q)\}$ 是一棵生成树，且 $c(\mathcal T')=c(\mathcal T^*)-c_{kl}+c_{pq}$。若 $c_{pq}<c_{kl}$，这与 $\mathcal T^*$ 的最优性矛盾；若只有 $c_{pq}\le c_{kl}$，则 $\mathcal T'$ 是一棵不含 $(k,l)$ 的最小生成树。
:::

这两个性质合在一起，给出了检验一棵给定的树是否最优的方法，对每条非树边只需考察一个圈。

:::{prf:theorem} 最优性条件
:label: thm-04-01-mst-optimality

生成树 $\mathcal T^*$ 是最小生成树，当且仅当对每条边 $(k,l)\notin\mathcal T^*$，$k$ 与 $l$ 之间的树路上的每条边 $(i,j)$ 都满足 $c_{ij}\le c_{kl}$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

*必要性。* 若非树边 $(k,l)$ 的树路上有一条边 $(i,j)$ 满足 $c_{ij}>c_{kl}$，则 $\mathcal T^*\cup\{(k,l)\}\setminus\{(i,j)\}$ 是一棵比 $\mathcal T^*$ 更便宜的生成树（{prf:ref}`cor-04-01-exchange`第 1 条）。

*充分性。* 设 $\mathcal T^*$ 满足上述条件，并在所有最小生成树中取与 $\mathcal T^*$ 公共边最多的一棵 $\mathcal T^o$。假设 $\mathcal T^o\ne\mathcal T^*$，取 $(i,j)\in\mathcal T^*\setminus\mathcal T^o$。从 $\mathcal T^*$ 中删去 $(i,j)$ 确定了一个割 $[S,\bar S]$，$(i,j)$ 是 $\mathcal T^*$ 在其中的唯一一条边。把 $(i,j)$ 加入 $\mathcal T^o$ 形成一个圈，它在第二条边 $(k,l)\in\mathcal T^o$ 处穿过这个割，并且 $(k,l)\notin\mathcal T^*$。$\mathcal T^*$ 中 $k$ 与 $l$ 之间的树路必须穿过这个割，因此要经过 $(i,j)$，由条件得 $c_{ij}\le c_{kl}$。对 $\mathcal T^o$ 应用必要性部分（非树边为 $(i,j)$，$(k,l)$ 是其圈上的边），得 $c_{kl}\le c_{ij}$。因此 $c_{ij}=c_{kl}$，而 $\mathcal T^o\cup\{(i,j)\} \setminus\{(k,l)\}$ 是一棵与 $\mathcal T^*$ 多一条公共边的最小生成树，矛盾。所以 $\mathcal T^*=\mathcal T^o$ 是最优的。
:::

:::{prf:corollary} 唯一性
:label: cor-04-01-unique

若每条非树边 $(k,l)$ 都严格比其树路上的每条边贵，则 $\mathcal T^*$ 是唯一的最小生成树；特别地，当所有边的费用互不相同时，最小生成树是唯一的。若某条非树边 $(k,l)$ 与其树路上的某条边 $(i,j)$ 费用相同，则 $\mathcal T^*\cup\{(k,l)\}\setminus\{(i,j)\}$ 是另一棵最小生成树。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

{prf:ref}`thm-04-01-mst-optimality`证明中的充分性部分表明，一棵最小生成树 $\mathcal T^o\ne\mathcal T^*$ 会给出一条非树边 $(k,l)$，它与其树路上的某条边 $(i,j)$ 满足 $c_{kl}=c_{ij}$，而严格不等式排除了这种情况。最后一个结论就是{prf:ref}`cor-04-01-exchange`中的交换，它不改变费用。
:::

(sec-04-01-prim)=
## Prim 算法

割性质提示我们只生长一棵树。从任一节点出发；在每次迭代中，考察已连通节点与其余节点之间的割，并加入该割中一条最便宜的边。这种“已连通集合”方法由 Prim (1957) 发表；Jarník (1930) 更早已经发现了它。

:::{prf:algorithm} Prim 算法
:label: alg-04-01-prim

**输入：** 费用为 $c_{ij}$ 的连通无向网络 $G=(\mathcal N,\mathcal A)$，以及起始节点 $r$。
**输出：** 一棵最小生成树 $\mathcal T$。

1. **初始化。** 令 $C_1=\{r\}$，$\bar C_1=\mathcal N\setminus\{r\}$，$\mathcal T=\varnothing$，$k=1$。
2. **终止检验。** 若 $\bar C_k=\varnothing$，则停止：$\mathcal T$ 是最小生成树。
3. **选择。** 选取一条满足 $i^*\in C_k$、$j^*\in\bar C_k$ 且费用最小的边 $(i^*,j^*)$，$c_{i^*j^*}=\min\{c_{ij}:(i,j)\in[C_k,\bar C_k]\}$（并列时任选其一）。
4. **连接。** 令 $\mathcal T\leftarrow\mathcal T\cup\{(i^*,j^*)\}$，$C_{k+1}=C_k\cup\{j^*\}$，$\bar C_{k+1}=\bar C_k\setminus\{j^*\}$，$k\leftarrow k+1$，转第 2 步。
:::

无需检验是否成圈：新边恰有一个端点 $j^*$ 不在已连通集合中，而经过 $j^*$ 的圈需要 $j^*$ 处有两条树边。*（有些对该算法的叙述还要求“加入 $j^*$ 不得产生圈”，这一条件总是成立的。）* 每次迭代把一个节点从 $\bar C_k$ 移到 $C_k$，因此算法恰好在 $\abs{\mathcal N}-1$ 次迭代后停止，得到 $\abs{\mathcal N}-1$ 条边，与{prf:ref}`thm-04-01-tree`的要求一致。

:::{prf:theorem} Prim 算法的正确性
:label: thm-04-01-prim

在连通网络上，无论起始节点和并列时的选择如何，Prim 算法都给出一棵最小生成树。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

我们用归纳法证明：在每次迭代开始时，当前的边集 $\mathcal T$ 包含于某棵最小生成树中。对 $\mathcal T=\varnothing$ 这显然成立。在第 $k$ 次迭代中，$\mathcal T$ 的每条边的两个端点都在 $C_k$ 中，因此 $\mathcal T$ 中没有边位于割 $[C_k,\bar C_k]$ 中，而选出的边是该割中最便宜的边。由割性质（{prf:ref}`thm-04-01-cut-property`，取 $F=\mathcal T$，$S=C_k$）知，$\mathcal T\cup\{(i^*,j^*)\}$ 仍包含于某棵最小生成树中。由于 $G$ 连通，这个割永远不会是空的。最后，$\mathcal T$ 是一棵生成树（连通，有 $\abs{\mathcal N}-1$ 条边），并包含于一棵同样有 $\abs{\mathcal N}-1$ 条边的最小生成树中；两者相等。
:::

**运行时间。** 若每次迭代都搜索整个割，每次迭代的代价为 $O(\abs{\mathcal A})$，总代价为 $O(\abs{\mathcal N}\abs{\mathcal A})$。更好的做法是为每个未连通节点 $j$ 保存一个标号 $[\gamma(j),P(j)]$：从 $j$ 到已连通集合的最便宜边的费用 $\gamma(j)$，以及该边另一端的节点 $P(j)\in C_k$（若没有这样的边，则 $\gamma(j)=\infty$）。于是第 3 步选取 $\gamma$ 最小的未连通节点；$j^*$ 加入之后，只需检查 $j^*$ 处的边：对每个满足 $c_{j^*j}<\gamma(j)$ 的相邻节点 $j\in\bar C_{k+1}$，令 $\gamma(j)=c_{j^*j}$，$P(j)=j^*$。若标号存放在数组中，每次选择需要 $O(\abs{\mathcal N})$ 时间，整个算法需要 $O(\abs{\mathcal N}^2)$ 时间，这对稠密网络已是最优的。若标号存放在二叉堆中，时间为 $O(\abs{\mathcal A}\log\abs{\mathcal N})$；若用斐波那契堆，则为 $O(\abs{\mathcal A}+\abs{\mathcal N}\log\abs{\mathcal N})$（Fredman & Tarjan, 1987）。这种标号方案与 Dijkstra 最短路算法（{prf:ref}`alg-04-02-dijkstra`）中的相同；只是更新不同，这里用 $c_{j^*j}$，那里用 $V(j^*)+c_{j^*j}$。

:::{prf:example} 在污水网络上执行 Prim 算法
:label: ex-04-01-sewage-prim

从主处理站即节点 7 出发，对{prf:ref}`ex-04-01-sewage`的污水网络执行 Prim 算法。当有多条连接并列费用最小时，选取引入编号较小站点的那一条；若仍并列，再选取在 $C_k$ 中一端编号较小的那一条。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

{numref}`tab-04-01-prim`列出了各次迭代。在第 1 次迭代中，割 $[\{7\},\{1,\dots,6\}]$ 包含 5–7 (3)、6–7 (4) 和 3–7 (8)；最便宜的是 5–7。在第 2 次迭代中，$C_2=\{5,7\}$ 的割包含 1–5 (9)、2–5 (3)、4–5 (8)、3–7 (8) 和 6–7 (4)，2–5 被加入。在第 3 次迭代中，$C_3=\{2,5,7\}$（{numref}`fig-04-01-prim-steps`a），割包含 1–2 (1)、2–3 (6)、2–4 (4)、1–5 (9)、4–5 (8)、3–7 (8) 和 6–7 (4)，1–2 被加入。在第 4 次迭代中，有两条费用为 4 的连接 2–4 和 6–7 并列；按规则选 2–4（它引入站点 4）。在第 5 次迭代中（{numref}`fig-04-01-prim-steps`b），割中最便宜的连接是 4–6 (3)。在第 6 次迭代中只剩下站点 3；它到已连通集合的连接为 1–3 (5)、3–4 (5)、2–3 (6)、3–7 (8) 和 3–6 (10)，1–3 与 3–4 之间的并列按规则选择 1–3。

```{list-table} 从节点 7 出发在污水网络上执行 Prim 算法（费用单位：百万美元）
:name: tab-04-01-prim
:header-rows: 1

* - 迭代 $k$
  - 已连通集合 $C_k$
  - $[C_k,\bar C_k]$ 中最便宜的连接
  - 加入的连接
  - 累计费用
* - 1
  - $\{7\}$
  - 5–7 (3)
  - 5–7
  - 3
* - 2
  - $\{5,7\}$
  - 2–5 (3)
  - 2–5
  - 6
* - 3
  - $\{2,5,7\}$
  - 1–2 (1)
  - 1–2
  - 7
* - 4
  - $\{1,2,5,7\}$
  - 2–4 (4), 6–7 (4)
  - 2–4
  - 11
* - 5
  - $\{1,2,4,5,7\}$
  - 4–6 (3)
  - 4–6
  - 14
* - 6
  - $\{1,2,4,5,6,7\}$
  - 1–3 (5), 3–4 (5)
  - 1–3
  - 19
```

六次迭代后 $\bar C_7=\varnothing$。应修建的隧道为

$$
\begin{gathered}
\mathcal T^*=\{5\text{–}7,\ 2\text{–}5,\ 1\text{–}2,\ 2\text{–}4,\ 4\text{–}6,\ 1\text{–}3\},\\
c(\mathcal T^*)=3+3+1+4+3+5=19,
\end{gathered}
$$

总费用为 1 900 万美元（{numref}`fig-04-01-sewage-mst`），比{prf:ref}`ex-04-01-sewage`中的可行布局少 1 100 万美元。两处并列已经表明最优布局不一定唯一：在第 4 次迭代中选择 6–7，或在第 6 次迭代中选择 3–4，会得到费用相同的其他树（见{prf:ref}`ex-04-01-sewage-kruskal`中的讨论）。

```{figure} figures/04-01-prim-steps.svg
:name: fig-04-01-prim-steps
:width: 560px
:class: dark-light
:alt: 污水网络的两个副本，已连通集合的节点被填充，树的连接以粗线表示，割中的连接以虚线表示

在污水网络上执行 Prim 算法。填充的节点构成已连通集合，粗线为当前的树，虚线为割 $[C_k,\bar C_k]$ 中的连接及其费用；灰色连接不在割中。(a) $C_3=\{2,5,7\}$：接下来加入连接 1–2（费用 1）。(b) $C_5=\{1,2,4,5,7\}$：接下来加入连接 4–6（费用 3）。
```

```{figure} figures/04-01-sewage-mst.svg
:name: fig-04-01-sewage-mst
:width: 560px
:class: dark-light
:alt: 污水网络，六条树连接以粗线表示，两条可替换的连接以虚线表示

污水网络的最小生成树（粗线，总费用 1 900 万美元）。虚线连接 6–7 和 3–4 分别与其圈上的树连接 2–4 和 1–3 费用并列；交换它们可以得到另外三棵最优树。
```
:::

(sec-04-01-kruskal)=
## Kruskal 算法

Kruskal 算法（Kruskal, 1956）考察的是边而不是节点：它按费用从低到高依次考察各条边，只要某条边不会与已保留的边构成圈，就保留它。在算法进行过程中，保留的边构成一个森林，其连通分支不断生长、合并。

:::{prf:algorithm} Kruskal 算法
:label: alg-04-01-kruskal

**输入：** 费用为 $c_{ij}$ 的无向网络 $G=(\mathcal N,\mathcal A)$。
**输出：** 一棵最小生成树 $\mathcal T$（若 $G$ 不连通，则为最小生成森林）。

1. **排序。** 将各边按费用非减的顺序排列为 $e_1,e_2,\dots,e_{\abs{\mathcal A}}$。令 $\mathcal T=\varnothing$（每个节点自成一个连通分支），$k=1$。
2. **终止检验。** 若 $\abs{\mathcal T}=\abs{\mathcal N}-1$ 或 $k>\abs{\mathcal A}$，则停止。
3. **检验该边。** 设 $e_k=(i,j)$。若 $i$ 与 $j$ 位于 $(\mathcal N,\mathcal T)$ 的不同连通分支中，则令 $\mathcal T\leftarrow\mathcal T\cup\{e_k\}$，并合并这两个连通分支；否则拒绝 $e_k$，因为它会构成圈。
4. **下一条边。** 令 $k\leftarrow k+1$，转第 2 步。
:::

:::{prf:theorem} Kruskal 算法的正确性
:label: thm-04-01-kruskal

在连通网络上，Kruskal 算法给出一棵最小生成树。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

被接受的边永远不会构成圈，因此 $\mathcal T$ 始终是森林。设 $e_k=(i,j)$ 被接受，$S$ 是此前 $i$ 所在的连通分支。$\mathcal T$ 中没有边位于割 $[S,\bar S]$ 中。该割中任何比 $e_k$ 便宜的边都在之前被考察过；当时它的两个端点也位于不同的连通分支中（连通分支只会合并），因此它会被接受，现在就是 $\mathcal T$ 中位于该割中的一条边，这是不可能的。因此 $e_k$ 是该割中最便宜的边，由割性质（{prf:ref}`thm-04-01-cut-property`）并用归纳法可知，$\mathcal T$ 始终包含于某棵最小生成树中。若算法停止时边数少于 $\abs{\mathcal N}-1$，则 $G$ 中有某条边连接最终的两个连通分支，而这条边在被考察时就会被接受。所以算法结束时得到一棵生成树。
:::

**运行时间与并查集。** 排序需要 $O(\abs{\mathcal A}\log\abs{\mathcal A})=O(\abs{\mathcal A}\log \abs{\mathcal N})$ 时间。第 3 步需要判断两个节点是否在同一连通分支中，并合并连通分支。一种简单的**并查集**（union–find）方案为每个节点存储其所在连通分支的名称：检验时比较两个名称，合并时把较小连通分支中的节点重新命名。一个节点只有在其所在连通分支的规模至少翻倍时才会被重新命名，因此至多被重新命名 $\log_2\abs{\mathcal N}$ 次，所有重新命名的总代价为 $O(\abs{\mathcal N}\log\abs{\mathcal N})$。排序占主导地位，Kruskal 算法的运行时间为 $O(\abs{\mathcal A}\log \abs{\mathcal N})$。对以弧表给出的稀疏网络，它是自然的选择，而 Prim 算法的数组实现适合稠密网络。（若使用并查集树和路径压缩，合并的代价几乎是线性的；见 Ahuja et al., 1993, Section 13.4。）

:::{prf:example} 污水网络的 Kruskal 算法与多重最优解
:label: ex-04-01-sewage-kruskal

对{prf:ref}`ex-04-01-sewage`的污水网络执行 Kruskal 算法，费用并列的连接按其端点编号递增的顺序排列，并检验最优布局是否唯一。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

排序后的列表开头为 1–2 (1)、2–5 (3)、4–6 (3)、5–7 (3)、2–4 (4)、6–7 (4)、1–3 (5)、3–4 (5)。{numref}`tab-04-01-kruskal`给出了各步决策。接受 2–4 之后，节点 6 和 7 已经在同一连通分支中，因此 6–7 会构成圈 6–7–5–2–4–6，被拒绝。接受 1–3 后已有六条连接，算法停止，其余七条连接无需考察。

```{list-table} 在污水网络上执行 Kruskal 算法（费用单位：百万美元）
:name: tab-04-01-kruskal
:header-rows: 1

* - 连接（费用）
  - 决策
  - 该步之后的连通分支
* - 1–2 (1)
  - 接受
  - $\{1,2\},\{3\},\{4\},\{5\},\{6\},\{7\}$
* - 2–5 (3)
  - 接受
  - $\{1,2,5\},\{3\},\{4\},\{6\},\{7\}$
* - 4–6 (3)
  - 接受
  - $\{1,2,5\},\{3\},\{4,6\},\{7\}$
* - 5–7 (3)
  - 接受
  - $\{1,2,5,7\},\{3\},\{4,6\}$
* - 2–4 (4)
  - 接受
  - $\{1,2,4,5,6,7\},\{3\}$
* - 6–7 (4)
  - 拒绝（成圈）
  - 不变
* - 1–3 (5)
  - 接受
  - $\{1,2,3,4,5,6,7\}$
```

被接受的连接恰好就是 Prim 算法求得的树 $\mathcal T^*$，费用为 $1+3+3+3+4+5=19$，即 1 900 万美元。

*最优性与唯一性。* {numref}`tab-04-01-sewage-check`应用了{prf:ref}`thm-04-01-mst-optimality`：对 $\mathcal T^*$ 之外的每条连接，给出其两端点之间的树路以及该路上最贵的连接。

```{list-table} 污水网络的非树连接及其树路
:name: tab-04-01-sewage-check
:header-rows: 1

* - 非树连接（费用）
  - 树路
  - 路上的最大费用
  - 比较
* - 1–4 (7)
  - 1–2–4
  - 4
  - $7>4$
* - 1–5 (9)
  - 1–2–5
  - 3
  - $9>3$
* - 2–3 (6)
  - 2–1–3
  - 5
  - $6>5$
* - 3–4 (5)
  - 3–1–2–4
  - 5（连接 1–3）
  - 并列
* - 3–6 (10)
  - 3–1–2–4–6
  - 5
  - $10>5$
* - 3–7 (8)
  - 3–1–2–5–7
  - 5
  - $8>5$
* - 4–5 (8)
  - 4–2–5
  - 4
  - $8>4$
* - 6–7 (4)
  - 6–4–2–5–7
  - 4（连接 2–4）
  - 并列
```

没有任何非树连接比其树路上的某条连接便宜，因此 $\mathcal T^*$ 是最优的。有两处并列，因此由{prf:ref}`cor-04-01-unique`，最优解不唯一：6–7 可以替换 2–4，3–4 可以替换 1–3。这两种交换是相互独立的（做了第一种交换之后，3–4 的树路为 3–1–2–5–7–6–4，其上最贵的连接仍是 1–3），因此共有四种最优布局，费用都是 1 900 万美元：

$$
\begin{gathered}
\{1\text{–}2,\ 1\text{–}3,\ 2\text{–}4,\ 2\text{–}5,\ 4\text{–}6,\ 5\text{–}7\},\\
\{1\text{–}2,\ 1\text{–}3,\ 2\text{–}5,\ 4\text{–}6,\ 5\text{–}7,\ 6\text{–}7\},\\
\{1\text{–}2,\ 2\text{–}4,\ 2\text{–}5,\ 3\text{–}4,\ 4\text{–}6,\ 5\text{–}7\},\\
\{1\text{–}2,\ 2\text{–}5,\ 3\text{–}4,\ 4\text{–}6,\ 5\text{–}7,\ 6\text{–}7\}.
\end{gathered}
$$

对全部 1 020 棵生成树的检查证实，费用为 1 900 万美元的恰好是这四棵，没有其他。只要并列时的选择适当，Prim 算法可以得到其中任何一棵；Kruskal 算法在并列连接按其他顺序排列时也是如此。
:::

(sec-04-01-boruvka)=
## Borůvka 算法

最古老的最小生成树方法由 Borůvka (1926) 为设计电力网络而发表，它同时作用于所有连通分支。在每一轮中，当前森林的每个连通分支（开始时即每个节点）选取一条离开它的最便宜的边，所有选出的边一起加入；每条边的合理性都由割性质对其自身连通分支的应用来保证。并列必须按一个固定的规则打破，例如按边的编号，使选出的边不会构成圈。每一轮至少使连通分支数减半，因此至多有 $\lceil\log_2\abs{\mathcal N}\rceil$ 轮，每轮的工作量为 $O(\abs{\mathcal A})$，总计 $O(\abs{\mathcal A}\log\abs{\mathcal N})$。在污水网络上，若并列时选取另一端编号较小的边，第一轮选出 1–2、1–3、2–5、4–6 和 5–7，留下连通分支 $\{1,2,3,5,7\}$ 和 $\{4,6\}$；第二轮用 2–4 把它们连起来，再次得到费用为 1 900 万美元的树。由于每个连通分支都独立行动，该方法是并行和分布式 MST 算法的基础。

## 本章回顾

### 小结

- 网络 $G=(\mathcal N,\mathcal A)$ 由节点以及弧（有向）或边（无向）组成；路、链、有向路、圈和连通性的定义见{prf:ref}`def-04-01-path`和{prf:ref}`def-04-01-connected`，各节点的度满足 {eq}`eq-04-01-handshake`。
- 树是没有圈的连通网络；等价地，它是有 $\abs{\mathcal N}-1$ 条边的连通网络，或有 $\abs{\mathcal N}-1$ 条边的无圈网络，或任意两个节点之间有唯一一条路的网络（{prf:ref}`thm-04-01-tree`）。加入一条非树边产生一个圈，删去一条树边产生一个割（{prf:ref}`cor-04-01-exchange`）。
- 网络可以存储为弧表、节点–弧关联矩阵（连通时秩为 $\abs{\mathcal N}-1$，{numref}`tab-04-01-incidence`）、邻接矩阵或长度矩阵，或者邻接表。
- 最小生成树问题 {eq}`eq-04-01-tree-cost` 要求费用最低的生成树。任何割中最便宜的边属于某棵最优树（{prf:ref}`thm-04-01-cut-property`）；圈上严格最贵的边不属于任何最优树（{prf:ref}`thm-04-01-cycle-property`）；一棵树是最优的，当且仅当没有任何非树边比其树路上的某条边便宜（{prf:ref}`thm-04-01-mst-optimality`），而并列情况决定唯一性（{prf:ref}`cor-04-01-unique`）。
- Prim 算法（{prf:ref}`alg-04-01-prim`）经 $\abs{\mathcal N}-1$ 次迭代生长出一棵树，使用标号时需 $O(\abs{\mathcal N}^2)$ 时间，使用堆时需 $O(\abs{\mathcal A}\log\abs{\mathcal N})$ 时间；Kruskal 算法（{prf:ref}`alg-04-01-kruskal`）借助并查集扫描排好序的边，需 $O(\abs{\mathcal A}\log\abs{\mathcal N})$ 时间；Borůvka 算法在并行的各轮中合并所有连通分支。
- 对污水网络，两种算法都给出费用为 1 900 万美元的布局（{numref}`fig-04-01-sewage-mst`），它是 1 020 棵生成树中四棵最优树之一。

### 常见错误

- **把最小生成树当作一组最短路线。** 最小生成树最小化的是各条连接的总费用，而不是从一个节点到其他节点的距离。在污水网络中，站点 3 沿 3–1–2–5–7 到达节点 7，费用为 12，而直接连接 3–7 的费用只有 8。最短路线是{doc}`04-02-shortest-paths`的主题（对比{prf:ref}`prop-04-02-tree`）。
- **只考察最后加入的节点。** Prim 算法从整个已连通集合 $C_k$ 出发选取最便宜的边，而不是只从最近加入的节点出发；后一种“最近邻”规则可能错过较早加入的节点处的便宜边。
- **在 Kruskal 算法中因为某条边的一个端点已在树中就拒绝它。** 只有当一条边的两个端点在同一连通分支中时才拒绝它。算法早期森林有许多连通分支，一条与已连通节点相接的边仍可能连接其中两个连通分支。
- **执行一次算法后就断言唯一性。** 执行一次只得到一棵最优树。它是否唯一要由{prf:ref}`cor-04-01-unique`判断：寻找与其树路上某条边费用并列的非树边，如{prf:ref}`ex-04-01-sewage-kruskal`所示。
- **把迭代次数或边数算错。** $\abs{\mathcal N}$ 个节点上的生成树有 $\abs{\mathcal N}-1$ 条边，因此 Prim 算法执行 $\abs{\mathcal N}-1$ 次迭代，Kruskal 算法在接受 $\abs{\mathcal N}-1$ 条边后即可停止。
- **把这些算法用于有向网络。** Prim 算法和 Kruskal 算法都忽略弧的方向；沿有向弧从一个根出发到达每个节点的最便宜结构是另一个问题。

## 参考文献

- Ahuja, R. K., Magnanti, T. L., & Orlin, J. B. (1993). *Network flows: Theory,
  algorithms, and applications*. Prentice Hall. (Chapter 2: notation,
  definitions and network representations; Chapter 13: minimum spanning trees,
  optimality conditions, Kruskal's, Prim's and Sollin's algorithms.)
- Borůvka, O. (1926). O jistém problému minimálním [On a certain minimal
  problem]. *Práce Moravské přírodovědecké společnosti, 3*, 37–58.
- Fredman, M. L., & Tarjan, R. E. (1987). Fibonacci heaps and their uses in
  improved network optimization algorithms. *Journal of the ACM, 34*(3),
  596–615.
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Sections 10.2 and 10.4: terminology of
  networks and the minimum spanning tree problem.)
- Jarník, V. (1930). O jistém problému minimálním [About a certain minimal
  problem]. *Práce Moravské přírodovědecké společnosti, 6*(4), 57–63.
- Kruskal, J. B. (1956). On the shortest spanning subtree of a graph and the
  traveling salesman problem. *Proceedings of the American Mathematical
  Society, 7*(1), 48–50. <https://doi.org/10.1090/S0002-9939-1956-0078686-7>
- Prim, R. C. (1957). Shortest connection networks and some generalizations.
  *Bell System Technical Journal, 36*(6), 1389–1401.
  <https://doi.org/10.1002/j.1538-7305.1957.tb01515.x>
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Sections 6.1–6.2: scope and definition of network models, minimal
  spanning tree algorithm.)
