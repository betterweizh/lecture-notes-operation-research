# 最短路问题

## 概述

在网络中寻找最便宜、最快或最可靠的路线，是求解次数最多的网络最优化问题：路线规划软件、通信网络中的分组路由以及交通分配模型每天要求解数以百万计的最短路问题；而许多看起来与路线无关的规划问题（设备更新、采购计划、项目进度安排），只要画出合适的网络，也都可以化为最短路问题。本章叙述最短路问题（shortest-path problem）及其变形，并从一个观察出发推导出全部内容：最短路的一部分本身也是最短路。由这一观察可以得到最短路长度所满足的 **Bellman 方程**（Bellman's equations）。当所有弧长都非负时，**Dijkstra 算法**按距离递增的顺序把节点标号逐个变为永久标号，从而对一个起点求解这些方程；我们证明它的正确性，计算它的工作量，由前驱标号恢复出路，并用一个反例说明为什么负弧会使它失效。**Floyd–Warshall 算法**求出所有节点对之间的最短路，允许有负弧，并能检测负圈。最后，把最短路问题写成线性规划问题，其对偶变量就是节点标号，并把它用作建模工具。本章以{doc}`04-01-networks-mst`（网络、路、圈、树以及节点–弧关联矩阵）为基础；{ref}`sec-04-02-lp`还用到{doc}`03-05-lp-duality`中线性规划的对偶问题和互补松弛性。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 叙述单点对、单源和所有点对最短路问题，以及每种方法适用的假设（弧长非负、不存在负圈）；
- 由子路的最优性推导 Bellman 方程；
- 用标号 $[V(j),P(j)]$ 手工执行 Dijkstra 算法，恢复最短路和最短路树，并证明该算法在弧长非负时的正确性；
- 执行 Floyd–Warshall 算法，由前驱矩阵恢复路，并检测负圈；
- 把最短路问题表述为线性规划问题，把其对偶问题解释为节点势，并用约化长度验证最优性；
- 把设备更新等决策问题建模为最短路问题。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $s$, $t$ | 最短路问题的起点和终点 |
| $c_{ij}$ | 弧 $(i,j)$ 的长度（距离、时间或费用） |
| $c(Q)$ | 路 $Q$ 的长度，即其各弧长度之和 |
| $V(j)$ | 节点 $j$ 的标号：到目前为止找到的从 $s$ 到 $j$ 的最好路的长度；算法终止时即最短路长度 |
| $P(j)$ | 节点 $j$ 在该路上的前驱节点 |
| $[V(j),P(j)]$ | Dijkstra 算法中节点 $j$ 的标号 |
| $S$, $T$ | 具有永久标号与具有临时标号的节点集 |
| $\mD_k=(d_{ij}^{(k)})$ | 第 $k$ 次迭代后的 Floyd–Warshall 距离矩阵 |
| $\mP_k=(p_{ij}^{(k)})$ | 第 $k$ 次迭代后的 Floyd–Warshall 前驱矩阵 |
| $\bar c_{ij}=c_{ij}+V(i)-V(j)$ | 弧 $(i,j)$ 关于标号 $V$ 的约化长度 |
| $\mathcal W$, $c(\mathcal W)$ | 一个有向圈及其长度 |

(sec-04-02-problem)=
## 最短路问题

设 $G=(\mathcal N,\mathcal A)$ 是一个有向网络，其中每条弧 $(i,j)$ 都有长度 $c_{ij}$。有向路 $Q=(i_1,i_2,\dots,i_p)$ 的**长度**为

$$
c(Q)=\sum_{k=1}^{p-1}c_{i_ki_{k+1}}.
$$ (eq-04-02-path-length)

:::{prf:definition} 最短路问题
:label: def-04-02-shortest-path

给定起点 $s$ 和终点 $t$，**最短路问题**要求一条从 $s$ 到 $t$ 的长度最小的有向路。它的常见变形有：

- **单点对**：一个起点 $s$ 和一个终点 $t$；
- **单源**（single-source shortest paths）：从 $s$ 到其他每个节点的最短路；
- **单终点**：从每个节点到 $t$ 的最短路（把所有弧反向，再从 $t$ 出发求解单源问题）；
- **所有点对**（all-pairs shortest paths）：每个有序节点对之间的最短路。

从 $s$ 到 $j$ 的最短路的长度称为从 $s$ 到 $j$ 的**距离**（distance）；若不存在有向路，则距离为 $+\infty$。
:::

单点对问题通常作为提前停止的单源问题来求解，因为这些方法在到达 $t$ 的过程中无论如何都会求出到 $s$ 附近各节点的距离。以下两点说明确定了本章的范围。

- **无向网络。** 一条两个方向都可以使用的边被替换为两条长度相同的弧 $(i,j)$ 和 $(j,i)$（{ref}`sec-04-01-graphs`）。于是长度为负的边会构成一个长度为负的圈 $i\to j\to i$，因此对无向网络，下面的方法都假设长度非负。
- **负圈。** 如果一个总长度为负的有向圈可以从 $s$ 到达，并且可以从它到达 $t$，那么沿着它一圈又一圈地走，就能使从 $s$ 到 $t$ 的路线任意短；这时在路线上提出的问题是无界的，而在路（不重复经过节点）上提出的问题则变得非常难解。因此本章所有方法都假设**不存在负圈**（negative cycle）。Dijkstra 算法需要更强的假设，即对所有弧 $c_{ij}\ge0$；Floyd–Warshall 算法只要求不存在负圈，并且在负圈存在时能把它检测出来。

:::{prf:example} 从节点 1 到节点 7 的最短路线
:label: ex-04-02-route

{numref}`fig-04-02-dijkstra-network`给出一个有七个节点、十二条弧的有向网络；弧上的标注为长度，例如以分钟计的行驶时间。求从节点 1 到节点 7 的一条最短路。

```{figure} figures/04-02-dijkstra-network.svg
:name: fig-04-02-dijkstra-network
:width: 540px
:class: dark-light
:alt: 由十二条有向弧相连的七个节点，弧长在 1 到 11 之间

{prf:ref}`ex-04-02-route` 的网络。弧上的标注为长度 $c_{ij}$（单位：分钟）；弧 $(3,1)$、$(3,2)$ 和 $(5,3)$ 指回起点方向。
```
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

最直接的方法是**枚举**（enumeration）：列出从节点 1 到节点 7 的每一条有向路，保留其中最短的一条。这个网络中共有七条这样的路（{numref}`tab-04-02-paths`），最短的是 $1\to4\to3\to7$，长度为 $4+3+5=12$ 分钟。

```{list-table} 从节点 1 到节点 7 的全部有向路，按长度排列（单位：分钟）
:name: tab-04-02-paths
:header-rows: 1

* - 路
  - 长度
* - $1\to4\to3\to7$
  - $4+3+5=12$
* - $1\to2\to5\to7$
  - $2+5+6=13$
* - $1\to4\to6\to7$
  - $4+1+9=14$
* - $1\to2\to5\to3\to7$
  - $2+5+8+5=20$
* - $1\to4\to3\to6\to7$
  - $4+3+7+9=23$
* - $1\to4\to3\to2\to5\to7$
  - $4+3+11+5+6=29$
* - $1\to2\to5\to3\to6\to7$
  - $2+5+8+7+9=31$
```
:::

像{prf:ref}`ex-04-02-route`的解答那样的枚举无法推广到大规模问题。在一个有 $q$ 层、每层 $r$ 个节点、每个节点都与下一层每个节点相连的网络中，从第一层到最后一层共有 $r^q$ 条路，20 层、每层 10 个节点时就已达 $10^{20}$ 条。下面的方法转而利用问题的结构。

(sec-04-02-bellman)=
## 子路的最优性与 Bellman 方程

{numref}`tab-04-02-paths`中的枚举做了许多无用功：路 $1\to4\to3\to7$ 和 $1\to4\to3\to6\to7$ 有公共的子路 $1\to4\to3$，其长度被计算了两次；而一旦知道 $1\to4\to3$ 是到达节点 3 的最好方式，任何以其他方式到达节点 3 的路都无需再考虑。下面的引理把这一点表述得更加精确。

:::{prf:lemma} 子路的最优性
:label: lem-04-02-subpath

设 $G$ 中没有负圈。若 $Q=(s,\dots,i,\dots,j)$ 是从 $s$ 到 $j$ 的最短路，则它从 $s$ 到 $i$ 的部分是从 $s$ 到 $i$ 的最短路。特别地，若 $(i,j)$ 是 $Q$ 的最后一条弧，则到 $j$ 的距离等于到 $i$ 的距离加上 $c_{ij}$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

把 $Q$ 写成 $Q_1$（从 $s$ 到 $i$）后接 $Q_2$（从 $i$ 到 $j$）。假设某条从 $s$ 到 $i$ 的路 $R$ 比 $Q_1$ 短。则 $R$ 后接 $Q_2$ 是一条从 $s$ 到 $j$ 的路线，长度为 $c(R)+c(Q_2)<c(Q)$。如果它重复经过某个节点，它就由一条从 $s$ 到 $j$ 的路加上一个或多个圈组成；这些圈的长度非负，因此去掉它们后得到一条从 $s$ 到 $j$ 的路，其长度不超过该路线的长度，从而比 $Q$ 短。这与 $Q$ 的最优性矛盾。
:::

这就是动态规划的**最优性原理**（principle of optimality）的最简单形式：最优决策序列的剩余部分对于它所出发的状态必须是最优的（{prf:ref}`thm-05-01-principle`）。把这个引理应用于每条最短路的最后一条弧，就得到关于距离的一个方程组。

:::{prf:theorem} Bellman 方程
:label: thm-04-02-bellman

设 $G$ 中没有负圈，并且每个节点都可以从 $s$ 到达。则从 $s$ 出发的距离 $V(j)$ 满足

$$
\begin{aligned}
V(s)&=0,\\
V(j)&=\min_{i:\,(i,j)\in\mathcal A}\bigl\{V(i)+c_{ij}\bigr\}
\qquad\text{对所有 } j\ne s.
\end{aligned}
$$ (eq-04-02-bellman)

此外，若每个圈的长度都为正，则这些距离是 {eq}`eq-04-02-bellman` 的唯一解。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

*距离满足这些方程。* $V(s)=0$，因为只由 $s$ 构成的路长度为 0，而且任何回到 $s$ 的路线都不是负的。对 $j\ne s$，每条从 $s$ 到 $j$ 的路都以某条弧 $(i,j)$ 结尾，且长度至少为 $V(i)+c_{ij}$，因此 $V(j)\ge\min_i\{V(i)+c_{ij}\}$。反之，对每条弧 $(i,j)$，到 $i$ 的最短路后接 $(i,j)$ 是一条到 $j$ 的长度为 $V(i)+c_{ij}$ 的路线，从中去掉圈不会使它变长，因此 $V(j)\le V(i)+c_{ij}$。

*唯一性。* 设 $U$ 是任意一个解。对每个 $j\ne s$，取一条达到最小值的弧 $(q_j,j)$，于是 $U(j)=U(q_j)+c_{q_jj}$。从 $j$ 出发沿这些弧向后追溯，永远不会回到已经访问过的节点：若沿选定的弧构成一个圈 $j_1\to\dots\to j_p\to j_1$，把相应的方程相加就会得到圈长为 0，而这已被排除。因此向后追溯的途径终止于 $s$，并描出一条从 $s$ 到 $j$、长度为 $U(j)$ 的路；所以 $U(j)\ge V(j)$。另一方面，沿一条最短路 $s=i_1,\dots,i_p=j$，把不等式 $U(i_{k+1})\le U(i_k)+c_{i_ki_{k+1}}$ 相加得到 $U(j)\le V(j)$。因此 $U=V$。
:::

对{prf:ref}`ex-04-02-route`的网络，利用{prf:ref}`ex-04-02-dijkstra`中求得的距离 $V(3)=7$、$V(5)=7$ 和 $V(6)=5$，节点 7 的方程为

$$
\begin{aligned}
V(7)&=\min\{V(3)+5,\ V(5)+6,\ V(6)+9\}\\
&=\min\{7+5,\ 7+6,\ 5+9\}=12,
\end{aligned}
$$

这里可以看出困难所在：这些方程一般是相互循环的（节点 3 的方程中含有 $V(5)$，节点 5 的方程中含有 $V(2)$，而节点 3 又有一条弧指回节点 2）。只有按照合适的顺序，才能逐个节点地求解它们。在无圈网络中，拓扑序（topological order）就能做到这一点（{ref}`sec-05-04-acyclic`）；在弧长非负时，Dijkstra 算法在运行过程中找出这一顺序，即距离递增的顺序。

(sec-04-02-dijkstra)=
## Dijkstra 算法

Dijkstra (1959) 发表的方法是一种**标号设定算法**（label-setting algorithm）。每个节点 $j$ 都带有一个标号（label）$[V(j),P(j)]$：$V(j)$ 是到目前为止找到的从 $s$ 到 $j$ 的最好路的长度，$P(j)$ 是该路上 $j$ 之前的节点。标号在还可能改进时是**临时的**（temporary），在已知其最优时变为**永久的**（permanent）。每次迭代**扫描**最近获得永久标号的节点 $i$：利用每条通向尚无永久标号节点的弧 $(i,j)$，尝试改进 $j$ 的标号。然后把 $V$ 最小的临时标号变为永久标号。

:::{prf:algorithm} Dijkstra 算法
:label: alg-04-02-dijkstra

**输入：** 长度为 $c_{ij}\ge0$ 的有向网络 $G=(\mathcal N,\mathcal A)$、起点 $s$，以及（可选的）终点 $t$。
**输出：** 对从 $s$ 可以到达的每个节点 $j$（若给定了终点，则对直到 $t$ 为止的节点），标号 $[V(j),P(j)]$：从 $s$ 到 $j$ 的距离以及 $j$ 在一条最短路上的前驱节点。

1. **初始化。** 给 $s$ 永久标号 $[V(s),P(s)]=[0,-]$。令 $S=\{s\}$，$T=\varnothing$，$i=s$。
2. **终止检验。** 若 $i=t$，或每个节点都已有永久标号，则停止。
3. **扫描节点 $i$。** 对每条满足 $j\notin S$ 的弧 $(i,j)\in\mathcal A$：
   - 若 $j$ 尚无标号，则给它临时标号 $[V(i)+c_{ij},\,i]$，并令 $T\leftarrow T\cup\{j\}$；
   - 若 $j\in T$ 且 $V(i)+c_{ij}<V(j)$，则令 $V(j)\leftarrow V(i)+c_{ij}$，$P(j)\leftarrow i$；
   - 否则 $j$ 的标号不变。
4. **设定永久标号。** 若 $T=\varnothing$，则停止：其余节点无法从 $s$ 到达。否则选取 $V(r)$ 最小的 $r\in T$（并列时任选其一），令 $S\leftarrow S\cup\{r\}$，$T\leftarrow T\setminus\{r\}$，$i\leftarrow r$，转第 2 步。
:::

到任何永久节点 $j$ 的路都由前驱标号**向后**读出：$j\leftarrow P(j)\leftarrow P(P(j))\leftarrow\cdots\leftarrow s$。该算法也适用于无向网络，此时每条边都要沿两个方向扫描。

:::{prf:example} 在七节点网络上执行 Dijkstra 算法
:label: ex-04-02-dijkstra

对{prf:ref}`ex-04-02-route`的网络，从 $s=1$ 到 $t=7$ 应用{prf:ref}`alg-04-02-dijkstra`。当两个临时标号并列时，把编号较小的节点设为永久节点。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**第 0 次迭代。** 节点 1 得到永久标号 $[0,-]$；$S=\{1\}$，$T=\varnothing$，$i=1$。

**第 1 次迭代。** 节点 1 不是终点。扫描它的弧，节点 2 得到标号 $[0+2,1]=[2,1]$，节点 4 得到标号 $[0+4,1]=[4,1]$，于是 $T=\{2,4\}$。最小的临时标号是 $V(2)=2$，节点 2 变为永久节点。

**第 2 次迭代。** 扫描节点 2：弧 $(2,5)$ 给节点 5 标号 $[2+5,2]=[7,2]$。现在 $T=\{4,5\}$，$V(4)=4<V(5)=7$，节点 4 变为永久节点。

**第 3 次迭代。** 扫描节点 4：节点 3 得到 $[4+3,4]=[7,4]$，节点 6 得到 $[4+1,4]=[5,4]$。在 $T=\{3,5,6\}$ 中最小的标号是 $V(6)=5$，节点 6 变为永久节点。

**第 4 次迭代。** 扫描节点 6：节点 7 得到 $[5+9,6]=[14,6]$。现在 $T=\{3,5,7\}$，$V(3)=V(5)=7$；按并列规则，节点 3 变为永久节点。

**第 5 次迭代。** 扫描节点 3：弧 $(3,1)$、$(3,2)$ 和 $(3,6)$ 通向永久节点，跳过。弧 $(3,7)$ 给出 $V(3)+5=12<14$，于是节点 7 的标号改进为 $[12,3]$。节点 5（$V(5)=7$）变为永久节点。

**第 6 次迭代。** 扫描节点 5：弧 $(5,3)$ 通向永久节点。弧 $(5,7)$ 给出 $7+6=13\ge12$，因此节点 7 的标号仍为 $[12,3]$。节点 7 是唯一的临时节点，变为永久节点。

**第 7 次迭代。** 刚刚变为永久节点的是终点，因此算法停止，$V(7)=12$。

{numref}`tab-04-02-dijkstra`汇总了这次运行。从节点 7 回溯得到 $7\leftarrow P(7)=3\leftarrow P(3)=4\leftarrow P(4)=1$，因此最短路为 $1\to4\to3\to7$，长度为 $4+3+5=12$ 分钟，与{numref}`tab-04-02-paths`一致。第 4 次迭代中的两个并列标号不影响结果：先把节点 5 设为永久节点，得到的最终标号相同。

```{list-table} 从节点 1 出发执行 Dijkstra 算法（长度单位：分钟）；临时标号为扫描之后的值
:name: tab-04-02-dijkstra
:header-rows: 1

* - 迭代
  - 扫描的节点 $i$
  - 标号的变化
  - 临时标号（$T$）
  - 设为永久的节点
* - 0
  - –
  - $1:[0,-]$
  - –
  - 1
* - 1
  - 1
  - $2:[2,1]$, $4:[4,1]$
  - $2:[2,1]$, $4:[4,1]$
  - 2
* - 2
  - 2
  - $5:[7,2]$
  - $4:[4,1]$, $5:[7,2]$
  - 4
* - 3
  - 4
  - $3:[7,4]$, $6:[5,4]$
  - $3:[7,4]$, $5:[7,2]$, $6:[5,4]$
  - 6
* - 4
  - 6
  - $7:[14,6]$
  - $3:[7,4]$, $5:[7,2]$, $7:[14,6]$
  - 3（与 5 并列）
* - 5
  - 3
  - $7:[14,6]\to[12,3]$
  - $5:[7,2]$, $7:[12,3]$
  - 5
* - 6
  - 5
  - 无（$13\ge12$）
  - $7:[12,3]$
  - 7
* - 7
  - –
  - 停止：$i=t=7$
  - –
  - –
```

*（起点的初始标号有时写作 $[0,0]$ 或 $(0,0)$；其第二个分量并不是节点，这里 $[0,-]$ 中的短横线表示起点没有前驱节点。）*

```{figure} figures/04-02-dijkstra-tree.svg
:name: fig-04-02-dijkstra-tree
:width: 560px
:class: dark-light
:alt: 七节点网络，每个节点旁标有标号 [V, P]，六条粗线弧构成以节点 1 为根的树

{prf:ref}`ex-04-02-dijkstra` 的最终标号 $[V(j),P(j)]$。粗线弧 $(P(j),j)$ 构成以节点 1 为根的最短路树；填充的节点位于长度为 12 的最短路 $1\to4\to3\to7$ 上。
```
:::

如{prf:ref}`ex-04-02-dijkstra`的解答所示，这次运行所做的比要求的更多：由于节点 7 恰好是最远的节点，每个节点都得到了它到节点 1 的距离。所有永久节点的弧 $(P(j),j)$ 在{numref}`fig-04-02-dijkstra-tree`中以粗线画出；它们构成一棵树，下一节将解释其原因。

(sec-04-02-dijkstra-correctness)=
## 正确性、最短路树与运行时间

该算法是贪婪的：一个标号一旦成为最小的临时标号，就被宣布为最终标号，以后再也不会重新考虑。下面的定理说明，当没有负弧时这样做为什么是安全的。

:::{prf:theorem} Dijkstra 算法的正确性
:label: thm-04-02-dijkstra

设对所有弧 $c_{ij}\ge0$。每当节点 $r$ 被设为永久节点时，$V(r)$ 就是从 $s$ 到 $r$ 的距离，并且从 $r$ 沿前驱标号回溯到 $s$ 所描出的是一条最短路。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

每个标号 $[V(j),P(j)]$ 都描述了一条实际存在的路：在设定该标号时有 $V(j)=V(P(j))+c_{P(j)j}$，且 $P(j)$ 当时已是永久节点。因此 $V(j)$ 永远不小于从 $s$ 到 $j$ 的距离 $d(j)$，我们只需证明 $V(r)\le d(r)$。

对节点变为永久节点的顺序用归纳法；结论对 $s$ 成立。设 $r$ 是下一个变为永久节点的节点，并假设 $S$ 中每个节点都已有正确的距离。任取一条从 $s$ 到 $r$ 的路 $Q$。它从 $S$ 内出发，终止于 $S$ 外，因此有第一条满足 $u\in S$、$w\notin S$ 的弧 $(u,w)$。扫描 $u$ 时，$w$ 得到了满足 $V(w)\le V(u)+c_{uw}=d(u)+c_{uw}$ 的标号，此后标号只会减小。由于 $Q$ 在 $w$ 之后的弧长度非负，

$$
c(Q)\ \ge\ d(u)+c_{uw}\ \ge\ V(w)\ \ge\ V(r),
$$

其中最后一个不等式成立，是因为 $w\in T$，而 $r$ 具有最小的临时标号。这对每条路 $Q$ 都成立，因此 $V(r)\le d(r)$。
:::

非负性假设恰好在一处起作用：$Q$ 在 $S$ 外第一个节点之后的部分不能使 $Q$ 变短。这一论证还表明，节点按距离非减的顺序变为永久节点，因此给定终点 $t$ 时，一旦 $t$ 成为永久节点，算法即可停止。

:::{prf:proposition} 最短路树
:label: prop-04-02-tree

当 Dijkstra 算法在所有可到达的节点都成为永久节点后停止时，弧 $(P(j),j)$，$j\ne s$，构成一棵以 $s$ 为根的树，其中从 $s$ 到每个节点 $j$ 的唯一一条路都是最短路。这棵树称为**最短路树**（shortest-path tree）。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

每个节点 $j\ne s$ 恰有一条进入它的树弧，来自 $P(j)$，并且 $P(j)$ 比 $j$ 更早成为永久节点。因此沿前驱节点追溯时，所访问的节点成为永久节点的时间严格地越来越早，最终必然到达 $s$，所以这些弧不含圈，并把每个节点都与 $s$ 连通：在可到达的节点上有 $\abs{\mathcal N}-1$ 条弧，它们构成一棵树（{prf:ref}`thm-04-01-tree`）。由{prf:ref}`thm-04-02-dijkstra`，到 $j$ 的树路长度为 $V(j)=d(j)$。
:::

因此单源问题有一个紧凑的答案：每个节点一个前驱节点。在{numref}`fig-04-02-dijkstra-tree`中，这棵树由弧 $(1,2)$、$(1,4)$、$(2,5)$、$(4,3)$、$(4,6)$ 和 $(3,7)$ 组成。注意它与最小生成树（{doc}`04-01-networks-mst`）的区别：最短路树最小化的是从根到每个节点的路的长度，而不是树弧的总长度。

**运行时间。** 每个节点被扫描一次，每条弧在扫描过程中被考察一次，总代价为 $O(\abs{\mathcal A})$。其余的工作是在第 4 步中寻找最小的临时标号，共 $\abs{\mathcal N}-1$ 次。若标号存放在数组中，每次搜索的代价为 $O(\abs{\mathcal N})$，因此算法的运行时间为 $O(\abs{\mathcal N}^2)$，这对稠密网络已是最优的。若临时标号存放在二叉堆中，每次选择和每次减小标号的代价为 $O(\log\abs{\mathcal N})$，总时间为 $O(\abs{\mathcal A}\log\abs{\mathcal N})$，这对道路网络这样的稀疏网络要快得多。斐波那契堆可以把它降低到 $O(\abs{\mathcal A}+\abs{\mathcal N}\log\abs{\mathcal N})$（Ahuja et al., 1993, Section 4.7）。这一方案与使用标号的 Prim 算法（{ref}`sec-04-01-prim`）相同，只是更新方式不同。

:::{prf:example} 负弧使 Dijkstra 算法失效
:label: ex-04-02-negative-arc

某网络有弧 $(1,2)$（长度为 2）、$(1,3)$（长度为 3）和 $(3,2)$（长度为 $-2$），完全没有圈。从节点 1 出发执行 Dijkstra 算法，并把它给节点 2 的标号与真实距离相比较。给所有长度加上一个常数使其非负，能否修复这一方法？
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

从节点 1 出发的 Dijkstra 算法（{numref}`fig-04-02-negative-arc`）给节点 2 标号 $[2,1]$，给节点 3 标号 $[3,1]$，并由于 $2<3$ 而把节点 2 设为永久节点。扫描节点 2 不改变任何标号。节点 3 以 $[3,1]$ 成为永久节点，扫描它时弧 $(3,2)$ 通向一个永久节点，因而被忽略。算法报告 $V(2)=2$，但路 $1\to3\to2$ 的长度为 $3-2=1$。

{prf:ref}`thm-04-02-dijkstra`的证明在 $S$ 外第一个节点之后的弧处失效：路 $1\to3\to2$ 经节点 3 离开 $S=\{1\}$ 时长度为 3，此后反而变短。给所有长度加上一个常数使其非负也无济于事：加上 2 后长度变为 4、5 和 0，最短路变成 $1\to2$（长度 4 对 5），因为加到每条弧上的常数对弧数多的路惩罚得比弧数少的路更重。

```{figure} figures/04-02-negative-arc.svg
:name: fig-04-02-negative-arc
:width: 360px
:class: dark-light
:alt: 三个节点；弧 1 到 2 长度为 2，弧 1 到 3 长度为 3，弧 3 到 2 长度为负 2，并标有 Dijkstra 算法给出的标号

{prf:ref}`ex-04-02-negative-arc` 的反例。Dijkstra 算法过早地把节点 2 的标号 $[2,1]$ 设为永久标号；粗线所示的路 $1\to3\to2$ 长度为 1。
```
:::

对于有负弧但没有负圈的情形，单源问题用**标号修正算法**（label-correcting algorithm）求解，这类方法允许标号在被使用之后再改变。其中最著名的是 Bellman–Ford 算法（Bellman, 1958; Ford, 1956），它分轮对所有弧施行更新 $V(j)\leftarrow\min\{V(j),V(i)+c_{ij}\}$；至多需要 $\abs{\mathcal N}-1$ 轮，总工作量为 $O(\abs{\mathcal N}\abs{\mathcal A})$，若第 $\abs{\mathcal N}$ 轮仍有改进，就说明存在负圈。在{ref}`sec-05-04-successive-approximation`中，它作为动态规划的逐次逼近法加以展开。

(sec-04-02-floyd)=
## 所有点对：Floyd–Warshall 算法

当所有长度都非负时，从每个节点出发执行一次 Dijkstra 算法，就能求出所有节点对之间的距离，用数组实现时需要 $O(\abs{\mathcal N}^3)$ 时间。Floyd–Warshall 算法（Floyd, 1962; Warshall, 1962）所需时间相同，但更为简单，而且允许有负弧。它作用于{ref}`sec-04-01-representations`中的**长度矩阵** $\mD_0=(d_{ij}^{(0)})$：$d_{ii}^{(0)}=0$，对于弧有 $d_{ij}^{(0)}=c_{ij}$，对于不存在的弧有 $d_{ij}^{(0)}=+\infty$。把节点编号为 $1,\dots, \abs{\mathcal N}$，并令

$$
d_{ij}^{(k)}=\text{从 } i \text{ 到 } j \text{、所有中间节点都在 } \{1,\dots,k\}
\text{ 中的最短路的长度}.
$$ (eq-04-02-floyd-def)

$k=0$ 时不允许有任何中间节点，得到的是长度矩阵；$k=\abs{\mathcal N}$ 时允许任何路，得到的是距离。

:::{prf:theorem} Floyd–Warshall 递推式
:label: thm-04-02-floyd

若 $G$ 中没有负圈，则对 $k=1,\dots,\abs{\mathcal N}$ 以及所有 $i$、$j$，

$$
d_{ij}^{(k)}=\min\bigl\{d_{ij}^{(k-1)},\ d_{ik}^{(k-1)}+d_{kj}^{(k-1)}\bigr\}.
$$ (eq-04-02-floyd)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

如果一条从 $i$ 到 $j$ 的路的中间节点都在 $\{1,\dots,k\}$ 中，就称它是*容许的*，并考虑一条最短的容许路。若它不经过 $k$，则其中间节点都在 $\{1,\dots,k-1\}$ 中，其长度至少为 $d_{ij}^{(k-1)}$。若它经过 $k$，则只经过一次，它从 $i$ 到 $k$ 和从 $k$ 到 $j$ 的两部分的中间节点都在 $\{1,\dots,k-1\}$ 中，因此其长度至少为 $d_{ik}^{(k-1)}+d_{kj}^{(k-1)}$。所以 $d_{ij}^{(k)}$ 至少等于 {eq}`eq-04-02-floyd` 的右端。反之，右端的两个数都是容许路线的长度：前者是一条容许路的长度，后者是一条从 $i$ 经 $k$ 到 $j$ 的路线的长度。若该路线的两部分有公共节点，则它含有一个圈，其长度非负，可以去掉（即{prf:ref}`lem-04-02-subpath`中的论证），剩下一条不更长的容许路。因此两个数都不小于 $d_{ij}^{(k)}$，等式成立。
:::

更新式 {eq}`eq-04-02-floyd` 称为关于 $i$、$k$ 和 $j$ 的**三元运算**（triple operation）。在第 $k$ 次迭代中，矩阵的第 $k$ 行和第 $k$ 列不变（因为 $d_{kk}^{(k-1)}=0$），因此可以就地覆盖矩阵。为了恢复出路，用**前驱矩阵**（predecessor matrix）$\mP_k=(p_{ij}^{(k)})$ 存储当前从 $i$ 到 $j$ 的最好路上 $j$ 之前的节点，与 Dijkstra 算法中的标号 $P(j)$ 完全相同。

:::{prf:algorithm} Floyd–Warshall 算法
:label: alg-04-02-floyd

**输入：** 节点为 $1,\dots,\abs{\mathcal N}$、长度为 $c_{ij}$（符号任意）的有向网络。
**输出：** 距离矩阵 $\mD$ 和前驱矩阵 $\mP$，或者一个负圈。

1. **初始化。** 令 $d_{ii}=0$；对每条弧 $(i,j)$ 令 $d_{ij}=c_{ij}$，$p_{ij}=i$；对其他每个节点对令 $d_{ij}=+\infty$，$p_{ij}=-$。
2. **对中间节点迭代。** 对 $k=1,\dots,\abs{\mathcal N}$（外层循环），以及每个节点对 $i$、$j$（内层循环）：若 $d_{ik}+d_{kj}<d_{ij}$，则令 $d_{ij}\leftarrow d_{ik}+d_{kj}$，$p_{ij}\leftarrow p_{kj}$。
3. **检验负圈。** 若对某个节点 $i$ 有 $d_{ii}<0$（这可以在每次迭代 $k$ 之后检验），则停止：有一个负圈经过 $i$。否则 $\mD$ 中存放的就是全部距离。
4. **恢复路。** 从 $i$ 到 $j$ 的最短路向后读出：$j\leftarrow p_{ij}\leftarrow p_{i,p_{ij}}\leftarrow\cdots\leftarrow i$。
:::

关于 $k$ 的循环必须是**外层**循环：定理是由完整的矩阵 $\mD_{k-1}$ 构造 $\mD_k$ 的。算法对 $\abs{\mathcal N}^2$ 个节点对中的每一个执行 $\abs{\mathcal N}$ 次三元运算，工作量为 $O(\abs{\mathcal N}^3)$，存储量为 $O(\abs{\mathcal N}^2)$，与弧数无关。对角线上出现负元素为什么意味着存在负圈：每个元素 $d_{ij}$ 始终是某条从 $i$ 到 $j$ 的路线的长度，因此 $d_{ii}<0$ 意味着有一条经过 $i$ 的长度为负的闭路线，而这样的路线含有一个负圈。反之，对 $k$ 用归纳法可知，第 $k$ 次迭代之后，每个元素 $d_{uv}$ 都不超过任何一条中间节点在 $\{1,\dots,k\}$ 中的从 $u$ 到 $v$ 的路的长度。若存在负圈 $\mathcal W$，设 $h$ 是其中编号最大的节点，$g$ 是其中另一个节点。$\mathcal W$ 的两部分，即从 $g$ 到 $h$ 和从 $h$ 到 $g$ 的部分，在 $k=h-1$ 时都是这样的路，因此最迟在第 $h$ 次迭代就会使 $d_{gg}\le d_{gh}+d_{hg}\le c(\mathcal W)<0$。

:::{prf:example} 在四节点网络上执行 Floyd–Warshall 算法
:label: ex-04-02-floyd

一家快递公司在四个仓库之间经营{numref}`fig-04-02-floyd-network`中的六条线路。弧上的标注为净费用（单位：百美元）；从仓库 1 到仓库 4 的线路可以带一批有偿回程货物，因此其净费用为 $-3$（即盈利 300 美元）。求每对仓库之间最便宜的路线。

```{figure} figures/04-02-floyd-network.svg
:name: fig-04-02-floyd-network
:width: 330px
:class: dark-light
:alt: 排成正方形的四个节点和六条有向弧；从 1 到 4 的弧长度为负 3

{prf:ref}`ex-04-02-floyd` 的快递网络。弧上的标注为净费用（单位：百美元）；十二条最便宜路线中有四条使用负弧 $(1,4)$。
```
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

该网络有一条负弧以及圈 $1\to4\to3\to1$（费用为 $-3+2+3=2$）、$1\to2\to3\to1$（11）和 $3\to4\to3$（9），都不是负圈，因此可以应用{prf:ref}`alg-04-02-floyd`。在下面的矩阵中，第 $i$ 行第 $j$ 列对应从 $i$ 到 $j$ 的路线，加框的元素是在该次迭代中发生改变的元素。

$$
\begin{gathered}
\mD_0=\begin{pmatrix}0&5&\infty&-3\\ \infty&0&3&\infty\\ 3&\infty&0&7\\ \infty&\infty&2&0\end{pmatrix},\\
\mP_0=\begin{pmatrix}-&1&-&1\\ -&-&2&-\\ 3&-&-&3\\ -&-&4&-\end{pmatrix}.
\end{gathered}
$$

**第 $k=1$ 次迭代**（经过仓库 1 的路线）。只有满足 $d_{i1}<\infty$ 的行（第 3 行）和满足 $d_{1j}<\infty$ 的列（第 2 列和第 4 列）可能改进：$d_{32}=d_{31}+d_{12}=3+5=8<\infty$，$d_{34}=3+(-3)=0<7$；两条新路都以从仓库 1 出发的弧结尾，因此 $p_{32}=p_{34}=p_{12}=p_{14}=1$。

$$
\begin{gathered}
\mD_1=\begin{pmatrix}0&5&\infty&-3\\ \infty&0&3&\infty\\ 3&\boxed{8}&0&\boxed{0}\\ \infty&\infty&2&0\end{pmatrix},\\
\mP_1=\begin{pmatrix}-&1&-&1\\ -&-&2&-\\ 3&\boxed{1}&-&\boxed{1}\\ -&-&4&-\end{pmatrix}.
\end{gathered}
$$

**第 $k=2$ 次迭代。** 第 2 列在第 1 行和第 3 行有限，第 2 行只在第 3 列有限：$d_{13}=5+3=8<\infty$，$p_{13}=p_{23}=2$；$d_{33}$ 将变为 $8+3=11>0$，因此不改变。

$$
\begin{gathered}
\mD_2=\begin{pmatrix}0&5&\boxed{8}&-3\\ \infty&0&3&\infty\\ 3&8&0&0\\ \infty&\infty&2&0\end{pmatrix},\\
\mP_2=\begin{pmatrix}-&1&\boxed{2}&1\\ -&-&2&-\\ 3&1&-&1\\ -&-&4&-\end{pmatrix}.
\end{gathered}
$$

**第 $k=3$ 次迭代。** 从每个仓库都能到达仓库 3，从仓库 3 也能到达每个仓库，有四个元素得到改进：$d_{21}=3+3=6$、$d_{24}=3+0=3$、$d_{41}=2+3=5$ 和 $d_{42}=2+8=10$，相应地 $p_{21}=p_{41}=p_{31}=3$，$p_{24}=p_{34}=1$，$p_{42}=p_{32}=1$。从仓库 1 出发经过仓库 3 的路线没有改进；例如 $d_{13}+d_{32}=8+8=16>5=d_{12}$。

$$
\begin{gathered}
\mD_3=\begin{pmatrix}0&5&8&-3\\ \boxed{6}&0&3&\boxed{3}\\ 3&8&0&0\\ \boxed{5}&\boxed{10}&2&0\end{pmatrix},\\
\mP_3=\begin{pmatrix}-&1&2&1\\ \boxed{3}&-&2&\boxed{1}\\ 3&1&-&1\\ \boxed{3}&\boxed{1}&4&-\end{pmatrix}.
\end{gathered}
$$

**第 $k=4$ 次迭代。** 经过仓库 4 的路线中只有一条得到改进：$d_{13}=d_{14}+d_{43}=-3+2=-1<8$，$p_{13}=p_{43}=4$。

$$
\begin{gathered}
\mD_4=\begin{pmatrix}0&5&\boxed{-1}&-3\\ 6&0&3&3\\ 3&8&0&0\\ 5&10&2&0\end{pmatrix},\\
\mP_4=\begin{pmatrix}-&1&\boxed{4}&1\\ 3&-&2&1\\ 3&1&-&1\\ 3&1&4&-\end{pmatrix}.
\end{gathered}
$$

$\mD_4$ 的对角线为零，这证实了不存在负圈，$\mD_4$ 中存放的就是全部最便宜路线的费用。例如，从仓库 4 到仓库 2 最便宜的路线费用为 $d_{42}=10$，即 1 000 美元；在 $\mP_4$ 的第 4 行中回溯得到 $2\leftarrow p_{42}=1\leftarrow p_{41}=3\leftarrow p_{43}=4$，即路线 $4\to3\to1\to2$，费用为 $2+3+5=10$。从仓库 2 到仓库 4 的路线为 $2\to3\to1\to4$，费用为 $3+3-3=3$（300 美元），它使用了盈利的线路；从仓库 1 到仓库 3 的路线为 $1\to4\to3$，净盈利 100 美元。用一个独立的程序实现加以验证，得到相同的矩阵和路。

*负圈。* 如果线路 $(1,4)$ 上的回程货物报酬更高，比如 $c_{14}=-6$，则圈 $1\to4\to3\to1$ 的费用为 $-6+2+3=-1$。重新执行算法，第 1 次迭代令 $d_{34}=3-6=-3$，第 3 次迭代给出 $d_{44}=d_{43}+d_{34}=2-3=-1<0$：对角线上的负元素揭示了圈 $4\to3\to1\to4$，在能够使用这个圈的仓库之间不存在最便宜的路线。
:::

:::{prf:remark} 序列矩阵
:label: rem-04-02-sequence-matrix

有些教科书（Taha, 2017, Section 6.3）记录的不是前驱节点而是中间节点：**序列矩阵**（sequence matrix）开始时对每条弧取 $s_{ij}=j$，每当第 $k$ 次迭代改进了 $d_{ij}$，就令 $s_{ij}\leftarrow k$。于是 $s_{ij}=j$ 表示直接使用弧 $(i,j)$，否则 $s_{ij}=k$ 表示该路线经过 $k$，再用同样的方法恢复从 $i$ 到 $k$ 和从 $k$ 到 $j$ 的路线。在{prf:ref}`ex-04-02-floyd`中，最终的序列矩阵各行为 $(-,2,4,4)$、$(3,-,3,3)$、$(1,1,-,1)$ 和 $(3,3,3,-)$；例如由 $s_{42}=3$、$s_{43}=3$（直接）以及 $s_{32}=1$、$s_{31}=1$ 和 $s_{12}=2$（直接）可以重建 $4\to3\to1\to2$。两种矩阵描述的是相同的路；前驱矩阵不需要递归，这就是本书采用它的原因。
:::

(sec-04-02-lp)=
## 线性规划形式与节点势

最短路问题是一种特殊的最小费用流问题（{ref}`sec-04-04-special-cases`）：以最小的费用把一个单位的流从 $s$ 送到 $t$，其中 $x_{ij}$ 是弧 $(i,j)$ 上的流量：

$$
\begin{aligned}
\min_{\vx}\ & \sum_{(i,j)\in\mathcal A}c_{ij}x_{ij}\\
\st\ & \sum_{j:(i,j)\in\mathcal A}x_{ij}-\sum_{j:(j,i)\in\mathcal A}x_{ji}
=\begin{cases}1,& i=s,\\-1,& i=t,\\0,&\text{其他情形，}\end{cases}\\
& x_{ij}\ge0\qquad\text{对所有 }(i,j)\in\mathcal A.
\end{aligned}
$$ (eq-04-02-lp)

左端各式就是 $\mA\vx$ 的各行，其中 $\mA$ 是节点–弧关联矩阵（{numref}`tab-04-01-incidence`）。每条从 $s$ 到 $t$ 的有向路都给出一个可行解：在其各弧上 $x_{ij}=1$，费用等于路的长度。反之，由于关联矩阵是全单模的（{prf:ref}`thm-04-04-tu`），{eq}`eq-04-02-lp` 的基本最优解是整数解，并且在不存在负圈时可以取为这样一条路。（若存在负圈，则该线性规划问题无界：流可以以负费用沿圈循环。）若要用一个线性规划问题求解单源问题，可令 $s$ 的供应量为 $\abs{\mathcal N}-1$，其他每个节点的需求量为 1；最优解沿一棵最短路树输送流量。

为节点 $i$ 的守恒方程配置一个自由对偶变量 $\pi_i$。$x_{ij}$ 的列在第 $i$ 行为 $+1$，在第 $j$ 行为 $-1$，因此 {eq}`eq-04-02-lp` 的对偶问题（{doc}`03-05-lp-duality`）为 $\max\{\pi_s-\pi_t:\ \pi_i-\pi_j\le c_{ij}\ \text{对所有 }(i,j)\}$。由于 $\mA$ 的各行之和为零，$\pi_i$ 只能确定到相差一个公共常数；记 $V(j)=\pi_s-\pi_j$（从而 $V(s)=0$），对偶问题就变为

$$
\begin{aligned}
\max_{V}\ & V(t)-V(s)\\
\st\ & V(j)-V(i)\le c_{ij}\qquad\text{对所有 }(i,j)\in\mathcal A.
\end{aligned}
$$ (eq-04-02-dual)

对偶变量是**节点势**（node potential）：每条弧限制了势沿它可以升高多少。一个物理图景是：一组珠子（节点）由长度为 $c_{ij}$ 的绳子（弧）连在一起；把 $t$ 尽可能远地拉离 $s$，绷紧的绳子就构成一条最短路。定义弧关于势的**约化长度**（reduced length）

$$
\bar c_{ij}=c_{ij}+V(i)-V(j).
$$ (eq-04-02-reduced-length)

对偶可行性意味着每条弧都有 $\bar c_{ij}\ge0$，而互补松弛性（{prf:ref}`thm-03-05-complementary-slackness`）意味着每条有流量的弧上 $\bar c_{ij}=0$。这就给出了一个与标号如何算出无关的最优性证书。

:::{prf:theorem} 最短路最优性条件
:label: thm-04-02-optimality

设 $V(s)=0$，并且对每个节点 $j$，$V(j)$ 是某条从 $s$ 到 $j$ 的有向路的长度。则各 $V(j)$ 是从 $s$ 出发的距离，当且仅当

$$
V(j)\le V(i)+c_{ij}\qquad\text{对所有 }(i,j)\in\mathcal A,
$$ (eq-04-02-optimality)

即当且仅当所有约化长度 {eq}`eq-04-02-reduced-length` 都非负。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

若某条弧有 $V(j)>V(i)+c_{ij}$，则到 $i$ 的路后接 $(i,j)$ 比 $V(j)$ 短（像前面一样去掉圈之后），因此 $V(j)$ 不是距离。反之，对任意一条路 $s=i_1,\dots,i_p=j$，把不等式 {eq}`eq-04-02-optimality` 沿其各弧相加，得到 $V(j)-V(s)\le c(Q)$，因此没有任何路比 $V(j)$ 短。
:::

在{prf:ref}`ex-04-02-route`的网络中，{prf:ref}`ex-04-02-dijkstra`的标号给出{numref}`tab-04-02-reduced-costs`中的约化长度。它们全都非负，因此这些标号就是距离；并且路 $1\to4\to3\to7$ 上的弧都有 $\bar c_{ij}=0$，这正是互补松弛性的要求。事实上，线性规划问题 {eq}`eq-04-02-lp` 与其对偶问题 {eq}`eq-04-02-dual` 的最优值都是 12。约化长度为零的弧恰好就是最短路树的弧。对单点对问题，最优的势不唯一：固定其他标号，3 与 5 之间的任何 $V(6)$ 都满足 {eq}`eq-04-02-optimality`，因为节点 6 不在最优路上，所以对于这样的节点，求解器可能返回不同的势。

```{list-table} 节点 1 到 7 的标号为 $V=(0,2,7,4,7,5,12)$ 时的约化长度 $\bar c_{ij}=c_{ij}+V(i)-V(j)$
:name: tab-04-02-reduced-costs
:header-rows: 1

* - 弧
  - $c_{ij}$
  - $\bar c_{ij}$
  - 弧
  - $c_{ij}$
  - $\bar c_{ij}$
* - $(1,2)$
  - 2
  - 0
  - $(3,7)$
  - 5
  - 0
* - $(1,4)$
  - 4
  - 0
  - $(4,3)$
  - 3
  - 0
* - $(2,5)$
  - 5
  - 0
  - $(4,6)$
  - 1
  - 0
* - $(3,1)$
  - 10
  - 17
  - $(5,3)$
  - 8
  - 8
* - $(3,2)$
  - 11
  - 16
  - $(5,7)$
  - 6
  - 1
* - $(3,6)$
  - 7
  - 9
  - $(6,7)$
  - 9
  - 2
```

约化长度也可以用来解释标号修正算法：Bellman–Ford 算法不断修复满足 $\bar c_{ij}<0$ 的弧，直到不再有这样的弧为止；而对某组可行的势，把 $c_{ij}$ 换成 $\bar c_{ij}$，就把一个有负弧（但没有负圈）的网络变成一个长度非负、最短路不变的网络，因为每条从 $i$ 到 $j$ 的路的长度都改变同一个量 $V(i)-V(j)$。在有负弧的稀疏网络上，所有点对问题正是这样求解的：先执行一次 Bellman–Ford 算法，再从每个节点出发执行 Dijkstra 算法。

(sec-04-02-modelling)=
## 用最短路建模

许多序贯决策问题都是某个网络上的最短路问题，该网络的节点是**状态**（state），弧是**决策**（decision）。由于状态按时间排序，这样的网络往往是无圈的。

**设备更新**（equipment replacement）。未来 $n$ 年需要一台机器。令节点 $i$（$i=1,\dots,n+1$）表示第 $i$ 年年初刚刚买入一台新机器（节点 $n+1$ 表示计划期末）。满足 $i<j$ 的弧 $(i,j)$ 表示“在第 $i$ 年年初买入一台机器，使用 $j-i$ 年，在第 $j$ 年年初折价换新”，其长度为

$$
c_{ij}=\text{购置价格}+\sum_{a=0}^{j-i-1}\text{役龄为 }a\text{ 时的运行费用}
-\text{役龄为 }j-i\text{ 时的旧设备折价} .
$$ (eq-04-02-replacement)

每个更新方案都是一条从节点 1 到节点 $n+1$ 的路，其长度就是该方案的总费用，因此最便宜的方案就是一条最短路。该网络是无圈的，因此其弧长可以取任意符号。{doc}`05-02-dp-equipment-replacement`用动态规划对带数值数据的这一模型进行了求解，并画出了同样的网络（{ref}`sec-05-02-network`）。

**最可靠路线。** 若弧 $(i,j)$ 以概率 $q_{ij}\in(0,1]$ 正常工作，且与其他弧相互独立，则一条路线正常工作的概率为 $\prod q_{ij}$。最大化这个乘积等价于最小化 $\sum(-\ln q_{ij})$，这是一个长度 $-\ln q_{ij}$ 非负的最短路问题，可以用 Dijkstra 算法求解。

更多此类模型见{doc}`08-05-examples-networks`：采购与库存计划（{prf:ref}`ex-08-05-purchasing`）、写成路问题的背包问题（{prf:ref}`ex-08-05-knapsack-path`），以及社交网络中的联系链（{prf:ref}`ex-08-05-hide-and-seek`）。

(sec-04-02-longest)=
## 无圈网络中的最长路

**最长路**（longest path）问题要求总长度最大的路。它是长度为 $-c_{ij}$ 的最短路问题，因此在含有正长度圈的网络中，它是一个含有负圈的最短路问题，一般是困难的。在**无圈**网络中根本没有圈，问题就变得容易了：按拓扑序给节点编号（每条弧都从编号较小的节点指向编号较大的节点），然后用 $\max$ 代替 $\min$ 求解 Bellman 方程 $V(j)=\max_{i:(i,j)\in\mathcal A}\{V(i)+c_{ij}\}$，按这一顺序逐个节点求解，需要 $O(\abs{\mathcal A})$ 时间。在关键路线法（critical path method，以活动为弧、事件为节点）中，无圈网络中的最长路决定了项目的工期；它还可以求解画成阶段–状态网络的背包型分配问题；具体过程见{ref}`sec-05-04-acyclic`。

## 本章回顾

### 小结

- 最短路问题要求给定节点之间长度最小的有向路（{prf:ref}`def-04-02-shortest-path`）；在不存在负圈时它是适定的。无向边变成一对方向相反的弧。
- 最短路的子路也是最短路（{prf:ref}`lem-04-02-subpath`），因此距离满足 Bellman 方程 {eq}`eq-04-02-bellman`（{prf:ref}`thm-04-02-bellman`）。
- Dijkstra 算法（{prf:ref}`alg-04-02-dijkstra`）按距离递增的顺序把标号 $[V(j),P(j)]$ 设为永久标号；它在 $c_{ij}\ge0$ 时是正确的（{prf:ref}`thm-04-02-dijkstra`），运行时间为 $O(\abs{\mathcal N}^2)$ 或 $O(\abs{\mathcal A}\log\abs{\mathcal N})$，给出一棵最短路树（{prf:ref}`prop-04-02-tree`），在有负弧时会失效（{prf:ref}`ex-04-02-negative-arc`）。在{prf:ref}`ex-04-02-dijkstra`中，从节点 1 到节点 7 的最短路为 $1\to4\to3\to7$，长度为 12（{numref}`fig-04-02-dijkstra-tree`）。
- Floyd–Warshall 算法（{prf:ref}`alg-04-02-floyd`）在外层循环中对 $k=1,\dots,\abs{\mathcal N}$ 施行三元运算 {eq}`eq-04-02-floyd`，运行时间为 $O(\abs{\mathcal N}^3)$；它允许有负弧，并通过对角线上的负元素揭示负圈（{prf:ref}`ex-04-02-floyd`）。
- 最短路问题就是线性规划问题 {eq}`eq-04-02-lp`；其对偶问题 {eq}`eq-04-02-dual` 以节点势 $V(j)$ 为变量，标号是最优的当且仅当所有约化长度 {eq}`eq-04-02-reduced-length` 都非负（{prf:ref}`thm-04-02-optimality`，{numref}`tab-04-02-reduced-costs`）。
- 设备更新问题和可靠性问题都是最短路模型（{ref}`sec-04-02-modelling`）；在无圈网络中求最长路是容易的（{ref}`sec-04-02-longest`）。

### 常见错误

- **在有负弧时使用 Dijkstra 算法，或通过加常数来“修正”负弧。** {prf:ref}`thm-04-02-dijkstra`的证明需要 $c_{ij}\ge0$，{prf:ref}`ex-04-02-negative-arc`展示了失效的情形。由于各条路的弧数不同，给每条弧加上一个常数会改变哪条路最短。应改用 Bellman–Ford 算法或 Floyd–Warshall 算法。
- **过早地设定永久标号，或在终点第一次得到标号时就停止。** 只有最小的临时标号才可以变为永久标号。在{prf:ref}`ex-04-02-dijkstra`中，节点 7 第一次得到的标号是 $[14,6]$；若在那里停止，就会报告长度为 14 而不是 12 的路 $1\to4\to6\to7$。
- **把 Floyd–Warshall 算法的循环顺序写错。** 中间节点 $k$ 必须是最外层循环（{prf:ref}`alg-04-02-floyd`）。若 $k$ 是最内层循环，元素 $d_{ij}$ 在对 $k$ 的一遍扫描之后就被确定下来，而所用的元素 $d_{ik}$、$d_{kj}$ 尚未达到最优，于是有些距离会偏大。
- **把三元运算的下标弄混。** 检验条件是 $d_{ik}+d_{kj}<d_{ij}$：路线先从 $i$ 到 $k$，再从 $k$ 到 $j$。写成 $d_{ik}+d_{jk}$ 就用了从 $j$ 到 $k$ 的路线，在有向网络中这是另一个量。
- **用 $k$ 更新前驱矩阵。** 在前驱矩阵中，更新式是 $p_{ij}\leftarrow p_{kj}$（从 $k$ 到 $j$ 的最好路线上 $j$ 之前的节点）。若存储 $k$，得到的是{prf:ref}`rem-04-02-sequence-matrix`中的序列矩阵，它必须递归地解码。
- **从前驱标号向前读出路。** $P(j)$ 指向 $j$ *之前*的节点；应从终点出发回溯到起点。
- **把最短路树与最小生成树混为一谈。** 最小生成树最小化其边的总长度，与任何起点无关；最短路树最小化从根出发的每一条路（{prf:ref}`prop-04-02-tree`）。

## 参考文献

- Ahuja, R. K., Magnanti, T. L., & Orlin, J. B. (1993). *Network flows: Theory,
  algorithms, and applications*. Prentice Hall. (Chapter 4: shortest paths with
  label-setting algorithms, acyclic networks, Dijkstra's algorithm and its heap
  implementations; Chapter 5: optimality conditions, label-correcting
  algorithms, negative cycles and the all-pairs problem.)
- Bellman, R. (1957). *Dynamic programming*. Princeton University Press.
  (Chapter 3: the principle of optimality.)
- Bellman, R. (1958). On a routing problem. *Quarterly of Applied
  Mathematics, 16*(1), 87–90. <https://doi.org/10.1090/qam/102435>
- Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
  *Numerische Mathematik, 1*(1), 269–271. <https://doi.org/10.1007/BF01386390>
- Floyd, R. W. (1962). Algorithm 97: Shortest path. *Communications of the
  ACM, 5*(6), 345. <https://doi.org/10.1145/367766.368168>
- Ford, L. R., Jr. (1956). *Network flow theory* (Report P-923). RAND
  Corporation.
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Section 10.3: the shortest-path problem.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Section 6.3: shortest-route problem, examples of applications,
  Dijkstra's and Floyd's algorithms, linear programming formulation.)
- Warshall, S. (1962). A theorem on Boolean matrices. *Journal of the ACM,
  9*(1), 11–12.
