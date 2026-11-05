# 单纯形法

## 概述

线性规划基本定理（{doc}`03-02-lp-geometry`）指出：有最优解的线性规划问题（linear program, LP）必有一个最优的**基本可行解**（basic feasible solution, BFS），也就是一个最优极点。单纯形法（simplex method）由 G. B. Dantzig 于 1947 年提出，它把这一结论变成了算法：从某个 BFS 出发，只要存在更好的相邻 BFS，就移动过去。本章从基的代数出发建立这一方法。首先把 LP 写成关于某个基的**典式**（canonical form），从中可以一眼读出当前的 BFS、它的目标函数值以及各非基变量的**检验数**（reduced cost）。检验数给出最优性检验；检验不通过时，令一个非基变量增大，由**最小比值检验**（minimum-ratio test）决定增大多少，并有一个基变量离开基。随后，同样的一步在**单纯形表**（simplex tableau）上通过一次高斯–若尔当**主元变换**（pivot）完成。整个方法在产品组合问题上完整地运行一遍，并画在其可行域上。本章最后讨论迭代结束时可能出现的异常或特殊情形：无界、多重最优解、退化与**循环**（cycling），以及防止循环的 **Bland 规则**（Bland's rule）。对于没有现成初始基的 LP，如何求得初始基是{doc}`03-04-two-phase-big-m`的内容。全章都要用到{doc}`01-01-linear-algebra`中的逆矩阵和高斯–若尔当消元法。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 把 LP 写成关于某个基的典式（方程形式、矩阵形式和表格形式），并从中读出 BFS 及其目标函数值；
- 计算检验数和单纯形乘子，并对极小化和极大化问题进行最优性检验；
- 由“下一个解必须保持可行”这一要求推导最小比值检验，并解释为什么所得的解仍是 BFS；
- 用精确分数在单纯形表上进行主元变换，并在二维可行域上追踪该方法的迭代路径；
- 证明当所有 BFS 均非退化时该方法有限终止，并在单纯形表中识别无界、多重最优解和退化；
- 解释循环是如何发生的，并应用 Bland 规则加以防止。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\mathcal B$ | 基变量的下标集，$\abs{\mathcal B}=m$ |
| $B(i)$ | 单纯形表第 $i$ 行基变量的下标，于是 $\vx_B=(x_{B(1)},\dots,x_{B(m)})^\T$ |
| $\bar{\va}_j=\mB^{-1}\mA_{\cdot j}$ | $x_j$ 的更新列，其元素为 $\bar a_{ij}$ |
| $s$, $r$ | 进基变量的下标、出基变量所在的行（主元行） |
| $\bar a_{rs}$ | 主元 |
| $\vd$ | 单纯形迭代的棱方向，$\vd_B=-\bar{\va}_s$，$d_s=1$ |
| $\theta$ | 一次主元变换的步长，即最小比值 |

其余符号（$\mA$、$\vb$、$\vc$、$\mB$、$\mN$、$\vx_B$、$\vx_N$、$\bar{\vb}$、$\bar c_j$、$\vpi$、$z$）与本书符号页中的相同。

(sec-03-03-canonical-form)=
## 基、基本可行解与典式

### 标准形与基

单纯形法作用于**标准形**（standard form）的 LP

$$
\begin{aligned}
\min_{\vx}\ & z=\vc^\T\vx\\
\st\ & \mA\vx=\vb,\\
& \vx\ge\vzero,
\end{aligned}
$$ (eq-03-03-standard-form)

其中 $\mA\in\R^{m\times n}$ 行满秩（秩为 $m$），且 $\vb\ge\vzero$。{doc}`03-02-lp-geometry`说明了如何把任意 LP 化为这种形式：不等式约束引入松弛变量或剩余变量，自由变量拆成两个非负部分之差，下界平移到零，右端项为负的约束两边乘以 $-1$，而 $\max\vc^\T\vx$ 等同于 $-\min(-\vc)^\T\vx$。具体的变换见{prf:ref}`ex-08-02-slack-surplus`和{prf:ref}`ex-08-02-bounds-free`。本书对极大化问题保持其原样，其规则与极小化的规则并列给出。

**基**（basis）是由 $m$ 个列下标组成的集合 $\mathcal B$，相应的列构成非奇异的 $m\times m$ 矩阵 $\mB=[\mA_{\cdot B(1)},\dots,\mA_{\cdot B(m)}]$。其余的列构成非基矩阵 $\mN$，变量相应地分成 $\vx_B$ 和 $\vx_N$，费用分成 $\vc_B$ 和 $\vc_N$。令 $\vx_N=\vzero$ 并解 $\mB\vx_B=\vb$，得到**基本解**（basic solution）$\vx_B=\bar{\vb}=\mB^{-1}\vb$，$\vx_N=\vzero$。当 $\bar{\vb}\ge\vzero$ 时它是**基本可行解**；当某个基变量取零值时它是**退化**（degenerate）的。BFS 恰好是可行多面体的极点（{prf:ref}`thm-03-02-extreme-bfs`）；只相差一个下标的两个基称为**相邻**（adjacent）的，它们的 BFS 或者相同，或者由多面体的一条棱相连（{doc}`03-02-lp-geometry`）。

### 典式

把 $\mA\vx=\vb$ 分块写成 $\mB\vx_B+\mN\vx_N=\vb$，再左乘 $\mB^{-1}$：

$$
\vx_B+\mB^{-1}\mN\vx_N=\mB^{-1}\vb=\bar{\vb}.
$$ (eq-03-03-canonical-constraints)

现在每个基变量恰好出现在一行中，且系数为 $+1$。把 $\vx_B=\bar{\vb}-\mB^{-1}\mN\vx_N$ 代入目标函数，得到

$$
\begin{aligned}
z&=\vc_B^\T\vx_B+\vc_N^\T\vx_N\\
&=\vc_B^\T\mB^{-1}\vb+\left(\vc_N^\T-\vc_B^\T\mB^{-1}\mN\right)\vx_N\\
&=\vc_B^\T\bar{\vb}+\sum_{j\notin\mathcal B}\bar c_j x_j ,
\end{aligned}
$$ (eq-03-03-objective-reduced)

其中**单纯形乘子**（simplex multipliers）和**检验数**（相对费用系数）为

$$
\begin{aligned}
\vpi^\T&=\vc_B^\T\mB^{-1},\\
\bar c_j&=c_j-\vpi^\T\mA_{\cdot j}, \qquad j=1,\dots,n.
\end{aligned}
$$ (eq-03-03-reduced-cost)

对第 $i$ 行的基变量 $j=B(i)$，其列为 $\mA_{\cdot j}=\mB\ve_i$，所以 $\bar c_j=c_j-\vc_B^\T\mB^{-1}\mB\ve_i=c_j-c_{B(i)}=0$。乘子 $\vpi$ 就是{doc}`03-05-lp-duality`中的对偶变量（{ref}`sec-03-05-reading-duals`）；在这里它只是一次算出全部 $\bar c_j$ 的便捷工具。下面定理的第 1 部分就是针对典式重新表述的{prf:ref}`thm-03-02-optimality-test`。

:::{prf:definition} 典式
:label: def-03-03-canonical-form

称标准形 LP 关于基 $\mathcal B$ 处于**典式**，如果它写成

$$
\begin{aligned}
\vx_B+\mB^{-1}\mN\vx_N&=\bar{\vb},\\
-z+\bar{\vc}_N^\T\vx_N&=-\vc_B^\T\bar{\vb},
\end{aligned}
$$

且 $\bar{\vb}=\mB^{-1}\vb\ge\vzero$：每个约束含有一个系数为 $+1$ 的基变量，该变量不出现在其他约束中，也不出现在目标行中，并且各右端项非负。
:::

这里把目标函数写成了多出来的一个方程 $-z+\vc^\T\vx=0$，其中 $-z$ 扮演一个永不出基的基变量。在典式中，当前 BFS 无须任何计算就能读出：非基变量为零，第 $i$ 行的基变量等于 $\bar b_i$，目标行的右端项等于 $-z$。

:::{prf:example} 产品组合问题的典式
:label: ex-03-03-product-mix

{prf:ref}`ex-03-01-product-mix`中的产品组合问题为一个成型车间安排每周产量：$x_1$ 和 $x_2$ 分别是每周生产的果汁杯和鸡尾酒杯的数量（单位：百箱），每百箱分别贡献 500 美元和 450 美元的利润。每百箱分别需要 6 小时和 5 小时的机时（可用 60 小时），以及 10 立方米和 20 立方米的仓储空间（可用 150 立方米），而果汁杯至多能售出 8 百箱：

$$
\begin{aligned}
\max_{\vx}\ & z=500x_1+450x_2\\
\st\ & 6x_1+5x_2\le 60,\\
& 10x_1+20x_2\le 150,\\
& x_1\le 8,\\
& x_1,x_2\ge 0.
\end{aligned}
$$

把这个 LP 写成标准形，并分别写出关于松弛基 $\mathcal B=\{3,4,5\}$ 和关于基 $\mathcal B=\{1,3,4\}$ 的典式。对每个基给出 BFS、其目标函数值以及非基变量的检验数。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

引入松弛变量 $x_3$（闲置机时）、$x_4$（未用仓储空间，立方米）和 $x_5$（未满足的果汁杯需求，百箱），标准形 {eq}`eq-03-02-product-mix-standard` 为

$$
\begin{gathered}
\mA=\begin{pmatrix}6&5&1&0&0\\10&20&0&1&0\\1&0&0&0&1\end{pmatrix},\\
\vb=\begin{pmatrix}60\\150\\8\end{pmatrix},\\
\vc=(500,450,0,0,0)^\T .
\end{gathered}
$$

松弛变量的列构成单位矩阵，所以该方程组已经是关于 $\mathcal B=\{3,4,5\}$ 的典式：

$$
\begin{aligned}
6x_1+5x_2+x_3&=60,\\
10x_1+20x_2+x_4&=150,\\
x_1+x_5&=8,\\
-z+500x_1+450x_2&=0 .
\end{aligned}
$$

其 BFS 为 $\vx=(0,0,60,150,8)^\T$，$z=0$：什么都不生产，所有资源都闲置。这里 $\vc_B=\vzero$，所以 $\vpi=\vzero$，检验数就是 $\bar c_1=500$ 和 $\bar c_2=450$。

第二个基（在{ref}`sec-03-03-matrix-step`中使用）是 $\mathcal B=\{1,3,4\}$，$\vx_B=(x_1,x_3,x_4)^\T$：

$$
\begin{gathered}
\mB=\begin{pmatrix}6&1&0\\10&0&1\\1&0&0\end{pmatrix},\\
\mB^{-1}=\begin{pmatrix}0&0&1\\1&0&-6\\0&1&-10\end{pmatrix},\\
\bar{\vb}=\mB^{-1}\vb=\begin{pmatrix}8\\12\\70\end{pmatrix}\ge\vzero ,
\end{gathered}
$$

所以这个基是可行的。由 $\vc_B=(500,0,0)^\T$ 得乘子 $\vpi^\T=\vc_B^\T\mB^{-1}=(0,0,500)$，对非基列 $\mN=[\mA_{\cdot2},\mA_{\cdot5}]$ 有

$$
\begin{gathered}
\mB^{-1}\mN=\begin{pmatrix}0&1\\5&-6\\20&-10\end{pmatrix},\\
\bar c_2=450-\vpi^\T\mA_{\cdot2}=450,\\
\bar c_5=0-\vpi^\T\mA_{\cdot5}=-500,
\end{gathered}
$$

且 $z=\vc_B^\T\bar{\vb}=500\cdot8=4000$。关于这个基的典式为

$$
\begin{aligned}
x_1+x_5&=8,\\
5x_2+x_3-6x_5&=12,\\
20x_2+x_4-10x_5&=70,\\
-z+450x_2-500x_5&=-4000 :
\end{aligned}
$$

BFS $\vx=(8,0,12,70,0)^\T$ 生产 8 百箱果汁杯，剩余 12 个机时和 70 立方米仓储空间未用，每周获利 4000 美元。{prf:ref}`ex-08-02-product-mix-bfs`用高斯–若尔当消元法得到了同样的方程。
:::

(sec-03-03-matrix-step)=
## 矩阵形式的单纯形迭代

### 最优性检验

式 {eq}`eq-03-03-objective-reduced` 对 $\mA\vx=\vb$ 的**每一个**解都成立，而不仅仅对 BFS 成立。它说明了当非基变量从零开始增大时目标函数如何变化：$\bar c_j$ 是当基变量随之调整以保持 $\mA\vx=\vb$ 时，$x_j$ 每增加一个单位 $z$ 的变化率。

:::{prf:theorem} 最优性检验
:label: thm-03-03-optimality

设 $\mathcal B$ 是极小化问题 {eq}`eq-03-03-standard-form` 的一个可行基，其 BFS 为 $\bar\vx$。

1. 若对所有 $j$ 有 $\bar c_j\ge 0$，则 $\bar\vx$ 是最优解。
2. 若 $\bar\vx$ 非退化且对某个 $s$ 有 $\bar c_s<0$，则 $\bar\vx$ 不是最优解。

对极大化问题，第 1 部分改为 $\bar c_j\le 0$，第 2 部分改为 $\bar c_s>0$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

1. 设 $\vx$ 是任一可行解。由 {eq}`eq-03-03-objective-reduced`，因为 $\bar c_j\ge0$ 且 $x_j\ge 0$，有 $\vc^\T\vx=\vc_B^\T\bar{\vb}+\sum_{j\notin\mathcal B}\bar c_jx_j\ge\vc_B^\T\bar{\vb}=\vc^\T\bar\vx$。
2. 把 $x_s$ 增大到一个小的 $\theta>0$，其他非基变量保持为零，并令 $\vx_B=\bar{\vb}-\theta\bar{\va}_s$。它满足 $\mA\vx=\vb$；又因为 $\bar{\vb}>\vzero$，当 $\theta$ 足够小时 $\vx_B\ge\vzero$。它的费用为 $\vc^\T\bar\vx+\theta\bar c_s<\vc^\T\bar\vx$。

对极大化问题，把上述结论应用于 $\min(-\vc)^\T\vx$，其检验数为 $-\bar c_j$。
:::

第 1 部分是*充分*条件。在退化的 BFS 处它并不必要：这时第 2 部分的步长可能为零，一个最优 BFS 可以由一个检验数不满足检验的基来表示。{prf:ref}`ex-08-02-a3-q16`给出了这样的基。最优点不变，只需通过步长为零的主元变换更换它的基，直到检验成立为止。（有些教材把该检验表述为“当且仅当”，这只在非退化的假设下成立。）

### 沿棱移动

设检验不通过：某个非基变量 $x_s$ 有 $\bar c_s<0$（对极大化问题为 $\bar c_s>0$）。单纯形法以最简单的方式改变 BFS：把 $x_s$ 从零增大到某个 $\theta\ge0$，其他非基变量都保持为零。于是基变量由 {eq}`eq-03-03-canonical-constraints` 确定：

$$
\begin{gathered}
\vx_B(\theta)=\bar{\vb}-\theta\,\bar{\va}_s,\\
x_s(\theta)=\theta,\\
z(\theta)=z+\theta\,\bar c_s .
\end{gathered}
$$ (eq-03-03-edge)

写成向量形式为 $\vx(\theta)=\bar\vx+\theta\vd$，其中 $d_s=1$，$\vd_B=-\bar{\va}_s$，其余分量为零。从几何上看，$\vd$ 指向可行多面体上从 $\bar\vx$ 出发的一条棱；$\theta$ 每增加一个单位，目标函数改善 $\abs{\bar c_s}$。

能走多远？目标函数一直在改善，所以唯一的限制是可行性 $\vx_B(\theta)\ge\vzero$，即对每一行 $i$ 有 $\bar b_i-\theta\bar a_{is}\ge 0$。$\bar a_{is}\le 0$ 的行不构成限制：其基变量保持不变或增大。$\bar a_{is}>0$ 的行要求 $\theta\le\bar b_i/\bar a_{is}$。因此最大可行步长由**最小比值检验**给出：

$$
\theta=\min\left\{\frac{\bar b_i}{\bar a_{is}}:\ \bar a_{is}>0\right\}
=\frac{\bar b_r}{\bar a_{rs}} .
$$ (eq-03-03-ratio-test)

在 $\theta$ 处，达到最小值的行 $r$ 的基变量 $x_{B(r)}$ 下降到零。它**出基**，而 $x_s$ **进基**：新基为 $\mathcal B'=\mathcal B\setminus\{B(r)\}\cup\{s\}$。还剩两个问题。

*新的点仍是 BFS 吗？* 点 $\vx(\theta)$ 按构造是可行的，且 $x_{B(r)}=0$，除 $x_s$ 外的非基变量均为零，所以只要 $\mathcal B'$ 是一个基，它就是 $\mathcal B'$ 的基本解。而 $\mathcal B'$ 确实是基：$\mB^{-1}\mB'$ 等于把单位矩阵的第 $r$ 列换成 $\bar{\va}_s$ 所得的矩阵，其行列式为 $\bar a_{rs}\neq 0$（{prf:ref}`prop-01-01-basis-exchange`）。所以 $\mB'$ 非奇异，$\mathcal B'$ 是与 $\mathcal B$ 相邻的可行基。

*如果 $\bar{\va}_s$ 没有正元素呢？* 这时对每个 $\theta\ge0$，$\vx(\theta)$ 都可行，而 $z(\theta)=z+\theta\bar c_s$ 无限改善：该 LP **无界**（unbounded），$\{\bar\vx+\theta\vd:\theta\ge0\}$ 就是使目标函数无界的一条可行解射线（{ref}`sec-03-03-special-cases`）。

因此，“哪个变量出基？”的答案是唯一确定的。正如{prf:ref}`ex-03-03-matrix-step`的解答所示，任何其他选择要么得不到基，要么得到不可行的基本解。

:::{prf:example} 用矩阵运算完成一次单纯形迭代
:label: ex-03-03-matrix-step

接着{prf:ref}`ex-03-03-product-mix`，从基 $\mathcal B=\{1,3,4\}$ 出发，BFS 为 $(x_1,x_2)=(8,0)$，$z=\$4000$。检验数为 $\bar c_2=450>0$ 和 $\bar c_5=-500$。用矩阵运算完成一次单纯形迭代：选取进基变量，逐一考察每个可能出基的基变量并与比值检验 {eq}`eq-03-03-ratio-test` 比较，再计算新的基逆矩阵、乘子、检验数和目标函数值。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

这是极大化问题，所以 $x_2$ 是进基变量：在基变量随之调整的情况下，每多生产一百箱鸡尾酒杯，利润增加 450 美元。它的更新列和当前值为

$$
\begin{gathered}
\bar{\va}_2=\mB^{-1}\mA_{\cdot2}=\begin{pmatrix}0\\5\\20\end{pmatrix},\\
\bar{\vb}=\begin{pmatrix}8\\12\\70\end{pmatrix} \quad\text{（依次为 }x_1,\ x_3,\ x_4\text{ 所在行）}.
\end{gathered}
$$

原则上有三个变量可能出基。

- **$x_1$ 出基。** $x_2,x_3,x_4$ 的列为 $(5,20,0)^\T$、$(1,0,0)^\T$ 和 $(0,1,0)^\T$，它们的第三个分量都为零，所以 $\mB$ 奇异，$\{2,3,4\}$ 不是基。这对应于 $\bar a_{12}=0$ 的情形：增大 $x_2$ 根本不改变 $x_1$。
- **$x_4$ 出基。** 基 $\{1,3,2\}$ 给出 $(x_1,x_3,x_2)=(8,-\tfrac{11}{2},\tfrac72)$。这是基本解，但不可行：点 $(x_1,x_2)=(8,\tfrac72)$ 位于仓储约束线和需求约束线上，需要 $6\cdot8+5\cdot\tfrac72=65.5$ 小时机时，超过了可用的 60 小时。步长 $\theta=\tfrac72$ 太长了。
- **$x_3$ 出基。** 基 $\{1,2,4\}$ 给出 $(x_1,x_2,x_4)=(8,\tfrac{12}{5},22)\ge\vzero$，是一个 BFS。

比值检验 {eq}`eq-03-03-ratio-test` 直接找到第三种情形：比值为 $12/5$（$x_3$ 所在行）和 $70/20=7/2$（$x_4$ 所在行），而 $x_1$ 所在行因为 $\bar a_{12}=0$ 没有比值。最小值 $\theta=\tfrac{12}{5}$ 属于 $x_3$。沿着棱 {eq}`eq-03-03-edge`，$x_3=12-5\theta$，$x_4=70-20\theta$；机时最先用完。

对新基，按 $(x_1,x_2,x_4)$ 的顺序，

$$
\begin{gathered}
\mB'^{-1}=\begin{pmatrix}0&0&1\\\tfrac15&0&-\tfrac65\\-4&1&14\end{pmatrix},\\
\vpi^\T=(500,450,0)\,\mB'^{-1}=(90,0,-40),
\end{gathered}
$$

所以 $\bar c_3=0-90=-90$，$\bar c_5=0-(-40)=40$，且 $z=500\cdot8+450\cdot\tfrac{12}{5}=5080$。正如 {eq}`eq-03-03-edge` 所预言的，$z$ 增加了 $\theta\bar c_2=\tfrac{12}{5}\cdot450=1080$。新的 BFS $(8,\tfrac{12}{5})$ 还不是最优的：$\bar c_5=40>0$ 说明少生产果汁杯、多生产鸡尾酒杯仍然有利可图。
:::

每一步都重新计算 $\mB^{-1}$ 是一种浪费。下面的单纯形表只用一次消元就更新典式中的全部量；{doc}`03-06-revised-dual-simplex`中的修正单纯形法又回到矩阵形式，并以乘积形式更新 $\mB^{-1}$（{ref}`sec-03-06-product-form`，{eq}`eq-01-01-product-form`）。

(sec-03-03-tableau)=
## 单纯形表

### 格式与符号约定

**单纯形表**以分离系数表的形式存储典式（{prf:ref}`def-03-03-canonical-form`）。在本书中，每个约束占一行，以其基变量为行标；目标行放在最后，行标为 $-z$：

| 基变量 | $\vx_B$ | $\vx_N$ | 右端项 |
|:--|:-:|:-:|:-:|
| $\vx_B$ | $\mI$ | $\mB^{-1}\mN$ | $\bar{\vb}=\mB^{-1}\vb$ |
| $-z$ | $\vzero^\T$ | $\bar{\vc}_N^\T=\vc_N^\T-\vpi^\T\mN$ | $-\vc_B^\T\bar{\vb}$ |

约束行的元素是 $\bar a_{ij}$，“右端项”列中是 $\bar b_i$，目标行存放检验数 $\bar c_j$（在基变量列中为零），其右端项为 $-z$。行标 $-z$ 提醒我们这一行是方程 $-z+\sum_j\bar c_jx_j=-z_0$，其中 $z_0$ 是当前的目标函数值。

各条规则只取决于 $\bar c_j$ 的符号：

```{list-table} 表格单纯形法的规则
:name: tab-03-03-rules
:header-rows: 1
:stub-columns: 1

* -
  - 极小化
  - 极大化
* - 最优的条件
  - 所有 $\bar c_j\ge0$
  - 所有 $\bar c_j\le0$
* - 进基变量（Dantzig 规则）
  - 最负的 $\bar c_j$
  - 最正的 $\bar c_j$
* - 无界的条件
  - $\bar c_s<0$ 且对所有 $i$ 有 $\bar a_{is}\le0$
  - $\bar c_s>0$ 且对所有 $i$ 有 $\bar a_{is}\le0$
* - 出基变量
  - 最小比值 {eq}`eq-03-03-ratio-test`
  - 最小比值 {eq}`eq-03-03-ratio-test`
* - 一次主元变换中 $z$ 的变化
  - $\theta\bar c_s\le0$
  - $\theta\bar c_s\ge0$
```

有些教材（例如 Taha (2017)）改用 **$z$ 行**：目标函数写成 $z-\vc^\T\vx=0$，于是该行存放 $-\bar c_j$，右端项为 $+z$，上面所有符号检验都要反过来（对极大化问题，当 $z$ 行所有元素都 $\ge0$ 时最优，进基变量取最负者）。两种约定给出相同的主元变换；本书始终使用 $-z$ 行。

### 主元变换

从基 $\mathcal B$ 转到 $\mathcal B'$，就是把 $x_s$ 的列变成第 $r$ 行为 $1$ 的单位向量，同时保持其他基变量的列仍为单位向量。这正是以**主元**（pivot element）$\bar a_{rs}$ 进行的一步高斯–若尔当消元（{prf:ref}`alg-01-01-pivot`）：

$$
\begin{aligned}
\text{第 } r \text{ 行：}&\quad \bar a'_{rj}=\frac{\bar a_{rj}}{\bar a_{rs}},\\
&\quad \bar b'_r=\frac{\bar b_r}{\bar a_{rs}}=\theta,\\
\text{第 } i \text{ 行（}i\ne r\text{）：}&\quad \bar a'_{ij}=\bar a_{ij}-\bar a_{is}\,\bar a'_{rj},\\
&\quad \bar b'_i=\bar b_i-\bar a_{is}\,\theta,\\
\text{目标行：}&\quad \bar c'_j=\bar c_j-\bar c_s\,\bar a'_{rj},\\
&\quad -z'=-z-\bar c_s\,\theta .
\end{aligned}
$$ (eq-03-03-pivot)

目标行与其他行一样，用同样的行变换更新。结果仍是典式：行变换把方程组换成等价的方程组，基变量的列是单位向量，目标行在这些列上为零。由于关于一个基的典式是唯一的（它就是 $\mB'^{-1}\mA\vx=\mB'^{-1}\vb$ 和 {eq}`eq-03-03-objective-reduced`），更新后的目标行恰好存放新的检验数，其右端项为 $-z'$，其中

$$
z'=z+\theta\,\bar c_s .
$$ (eq-03-03-z-change)

由于 $\theta\ge0$，主元变换永远不会使目标函数变差，而当 $\theta>0$ 时严格改善目标函数。$\theta=0$ 恰好出现在出基变量已经为零的时候，即在退化的 BFS 处（{ref}`sec-03-03-degeneracy`）。

主元变换还顺带给出了 $\mB'^{-1}$：在初始表中构成单位矩阵的那些列（这里是松弛变量的列）在以后每张表中都存放着 $\mB^{-1}$，因为它们记录了到目前为止所做的行变换。当初始基变量的费用为零时，目标行的这些列存放 $-\vpi^\T$。

(sec-03-03-algorithm)=
## 单纯形算法

:::{prf:algorithm} 单纯形法（表格形式）
:label: alg-03-03-simplex

**输入：** 一个标准形 LP，以及关于可行基 $\mathcal B$ 处于典式的单纯形表。
**输出：** 一个最优 BFS，或一条使目标函数无界的射线。

1. **最优性检验。** 若对所有 $j$ 有 $\bar c_j\ge0$（极大化：$\bar c_j\le0$），停止：当前 BFS 是最优解（{prf:ref}`thm-03-03-optimality`）。
2. **进基变量。** 选取一个 $\bar c_s<0$（极大化：$\bar c_s>0$）的非基变量 $x_s$。Dantzig 规则取最负（最正）的 $\bar c_s$；Bland 规则（{prf:ref}`alg-03-03-bland`）取下标最小者。
3. **无界性检验。** 若对每一行 $i$ 有 $\bar a_{is}\le0$，停止：LP 沿 {eq}`eq-03-03-edge` 的射线 $\bar\vx+\theta\vd$，$\theta\ge0$ 无界。
4. **比值检验。** 计算 $\theta=\min\{\bar b_i/\bar a_{is}:\bar a_{is}>0\}$，并选取一个达到该最小值的行 $r$（并列时取基变量下标最小的行）。变量 $x_{B(r)}$ 出基。
5. **主元变换。** 按 {eq}`eq-03-03-pivot` 以 $\bar a_{rs}$ 为主元进行主元变换，把第 $r$ 行的行标改为 $x_s$，转第 1 步。
:::

第 2 步留有选择的余地。Dantzig 规则是经典规则，通常只需较少的迭代次数，但并不总能使主元变换次数最少（{prf:ref}`ex-08-02-a3-q8`把它与“最小正检验数”规则作了比较）。任何选取改善变量的规则都是正确的；选择方式影响速度，在退化顶点处还影响能否终止。

:::{prf:example} 用单纯形法求解产品组合问题
:label: ex-03-03-simplex-run

从松弛基出发，用 Dantzig 规则求解{prf:ref}`ex-03-03-product-mix`。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

每张表的主元用粗体表示。

```{table} 产品组合问题的初始表（基 $\{x_3,x_4,x_5\}$，顶点 $(0,0)$）。
:name: tab-03-03-pm-0

| 基变量 | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | 右端项 |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $6$ | $5$ | $1$ | $0$ | $0$ | $60$ |
| $x_4$ | $10$ | $20$ | $0$ | $1$ | $0$ | $150$ |
| $x_5$ | $\boldsymbol{1}$ | $0$ | $0$ | $0$ | $1$ | $8$ |
| $-z$ | $500$ | $450$ | $0$ | $0$ | $0$ | $0$ |
```

**第 1 次迭代。** 最大的检验数是 $\bar c_1=500$，所以 $x_1$ 进基。比值为 $60/6=10$、$150/10=15$ 和 $8/1=8$；最小值 $\theta=8$ 在 $x_5$ 所在行，所以 $x_5$ 出基：需求上限先于机时和仓储空间达到。以 $1$ 为主元进行主元变换，即从其他各行分别减去 $x_5$ 行的 6 倍、10 倍和 500 倍，$z$ 增加 $8\cdot500=4000$。

```{table} 表 1（基 $\{x_3,x_4,x_1\}$，顶点 $(8,0)$）。
:name: tab-03-03-pm-1

| 基变量 | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | 右端项 |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $0$ | $\boldsymbol{5}$ | $1$ | $0$ | $-6$ | $12$ |
| $x_4$ | $0$ | $20$ | $0$ | $1$ | $-10$ | $70$ |
| $x_1$ | $1$ | $0$ | $0$ | $0$ | $1$ | $8$ |
| $-z$ | $0$ | $450$ | $0$ | $0$ | $-500$ | $-4000$ |
```

这正是{prf:ref}`ex-03-03-product-mix`中用矩阵运算求得的典式，只是行的顺序不同。

**第 2 次迭代。** 现在 $\bar c_2=450$ 是唯一的正元素，所以 $x_2$ 进基。比值为 $12/5$ 和 $70/20=7/2$（$x_1$ 所在行 $\bar a_{12}=0$，没有比值）；与{prf:ref}`ex-03-03-matrix-step`一样，$x_3$ 以 $\theta=\tfrac{12}{5}$ 出基，$z$ 增加 $\tfrac{12}{5}\cdot450=1080$，达到 5080。

```{table} 表 2（基 $\{x_2,x_4,x_1\}$，顶点 $(8,\tfrac{12}{5})$）。
:name: tab-03-03-pm-2

| 基变量 | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | 右端项 |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_2$ | $0$ | $1$ | $\tfrac{1}{5}$ | $0$ | $-\tfrac{6}{5}$ | $\tfrac{12}{5}$ |
| $x_4$ | $0$ | $0$ | $-4$ | $1$ | $\boldsymbol{14}$ | $22$ |
| $x_1$ | $1$ | $0$ | $0$ | $0$ | $1$ | $8$ |
| $-z$ | $0$ | $0$ | $-90$ | $0$ | $40$ | $-5080$ |
```

**第 3 次迭代。** 需求上限的松弛变量 $x_5$ 有 $\bar c_5=40>0$，于是进基：车间应当*少*生产果汁杯。$x_5$ 的列为 $(-\tfrac65,14,1)^\T$；第一个元素为负（$x_5$ 增大时 $x_2$ 也增大），比值为 $22/14=\tfrac{11}{7}$ 和 $8/1=8$。所以 $x_4$ 以 $\theta=\tfrac{11}{7}$ 出基，$z$ 增加 $\tfrac{11}{7}\cdot40=\tfrac{440}{7}$。

```{table} 表 3，最优表（基 $\{x_2,x_5,x_1\}$，顶点 $(\tfrac{45}{7},\tfrac{30}{7})$）。
:name: tab-03-03-pm-3

| 基变量 | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | 右端项 |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_2$ | $0$ | $1$ | $-\tfrac{1}{7}$ | $\tfrac{3}{35}$ | $0$ | $\tfrac{30}{7}$ |
| $x_5$ | $0$ | $0$ | $-\tfrac{2}{7}$ | $\tfrac{1}{14}$ | $1$ | $\tfrac{11}{7}$ |
| $x_1$ | $1$ | $0$ | $\tfrac{2}{7}$ | $-\tfrac{1}{14}$ | $0$ | $\tfrac{45}{7}$ |
| $-z$ | $0$ | $0$ | $-\tfrac{550}{7}$ | $-\tfrac{20}{7}$ | $0$ | $-\tfrac{36000}{7}$ |
```

**最优性。** 所有检验数都 $\le0$，所以表 3 是最优表：每周生产 $x_1^*=\tfrac{45}{7}\approx6.429$ 百箱和 $x_2^*=\tfrac{30}{7}\approx4.286$ 百箱，每周利润为 $z^*=\tfrac{36000}{7}\approx\$5142.86$。机时和仓储空间全部用完（$x_3=x_4=0$），需求上限还剩 $x_5=\tfrac{11}{7}\approx1.571$ 百箱的余量。非基变量的检验数都严格为负，所以这个最优解是唯一的（{ref}`sec-03-03-special-cases`）。

*（此表的较早版本把 $x_2$ 行、$x_3$ 列的元素写成了 $\tfrac17$；正确的元素是 $-\tfrac17$，这也是唯一与目标行一致的值：$0-\left(450\cdot(-\tfrac17)+500\cdot\tfrac27\right)=-\tfrac{550}{7}$。）*

**验算。** 该解满足 $6\cdot\tfrac{45}{7}+5\cdot\tfrac{30}{7}=60$、$10\cdot\tfrac{45}{7}+20\cdot\tfrac{30}{7}=150$ 和 $\tfrac{45}{7}\le8$，且 $500\cdot\tfrac{45}{7}+450\cdot\tfrac{30}{7}=\tfrac{36000}{7}$；各步增量之和 $4000+1080+\tfrac{440}{7}=\tfrac{36000}{7}$ 也吻合。目标行的松弛变量列中是 $-\pi_1=-\tfrac{550}{7}$ 和 $-\pi_2=-\tfrac{20}{7}$：每周多一个机时价值 $\tfrac{550}{7}\approx\$78.57$，多一立方米仓储空间价值 $\tfrac{20}{7}\approx\$2.86$。这些影子价格是{doc}`03-05-lp-duality`（{prf:ref}`ex-03-05-product-mix-dual`，{prf:ref}`thm-03-05-shadow-price`）和{doc}`03-07-sensitivity-analysis`的主题。

```{figure} figures/03-03-simplex-path.svg
:name: fig-03-03-simplex-path
:width: 520px
:class: dark-light
:alt: 产品组合问题的可行域及单纯形法路径：从原点沿 x1 轴到 (8,0)，向上到 (8,12/5)，再到最优点 (45/7,30/7)

单纯形法在产品组合问题上的路径。每张表对应一个顶点：$\vx^{(0)}=(0,0)$、$\vx^{(1)}=(8,0)$、$\vx^{(2)}=(8,\tfrac{12}{5})$，以及最优点 $\vx^{(3)}=(\tfrac{45}{7},\tfrac{30}{7})$，虚线等值线 $500x_1+450x_2=\tfrac{36000}{7}$ 在该点与可行域相切。每次主元变换都沿一条棱移动到相邻顶点。空心点 $(8,\tfrac72)$ 是第二步中若改由 $x_4$ 而非 $x_3$ 出基所得到的不可行基本解（{prf:ref}`ex-03-03-matrix-step`）。
```
:::

在{numref}`fig-03-03-simplex-path`的路径上，Dantzig 规则因为 $x_1$ 的检验数较大而先选了 $x_1$，该方法用了三次主元变换。如果先让 $x_2$ 进基，路径 $(0,0)\to(0,\tfrac{15}{2})\to(\tfrac{45}{7},\tfrac{30}{7})$ 只需两次就能到达最优点：最大的检验数衡量的是进基变量每增加一个单位带来的改善，而不是每次主元变换带来的改善。

(sec-03-03-special-cases)=
## 终止、无界与多重最优解

### 有限终止

:::{prf:theorem} 无退化时的有限终止
:label: thm-03-03-finite

若 LP 的每个 BFS 都非退化，则{prf:ref}`alg-03-03-simplex`无论采用何种选取改善进基变量的规则，都在有限次迭代后停止：要么得到最优基，要么得到一条使目标函数无界的射线。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

考虑极小化问题。在每次不停止的迭代中，$\bar c_s<0$，且由非退化性 $\theta=\bar b_r/\bar a_{rs}>0$，所以由 {eq}`eq-03-03-z-change`，目标函数严格下降。BFS 的目标函数值由其基决定，因此没有哪个基会出现两次。基至多有 $\binom{n}{m}$ 个，因而迭代次数至多为这么多，并且算法只能在第 1 步（最优）或第 3 步（无界）停止。若在第 1 步停止，则由{prf:ref}`thm-03-03-optimality`，该 BFS 是最优的。
:::

有了一个初始可行基（{doc}`03-04-two-phase-big-m`说明如何求得），这就对非退化问题构造性地证明了线性规划基本定理（{prf:ref}`thm-03-02-fundamental`）：可行且有界的 LP 有最优 BFS。界 $\binom{n}{m}$ 非常大，而且已知存在一些例子（Klee & Minty, 1972），在其上 Dantzig 规则要访问指数多个顶点；但在实际中，主元变换的次数通常只是 $m$ 的一个不大的倍数。

### 无界

当进基列没有正元素时，算法在第 3 步停止。由 {eq}`eq-03-03-edge`，点 $\bar\vx+\theta\vd$ 对所有 $\theta\ge0$ 都可行，而目标函数改变 $\theta\bar c_s$，没有界限。可行域包含整条射线，该 LP 没有最优解。有两点说明：

- 可行域*无界*并不意味着*目标函数*无界。LP 的无界性是可行域与费用向量共同决定的性质，而单纯形表只对它所考察的那个进基列检测无界性；另一个改善列可能仍有正元素。
- 在建模中，无界的 LP 几乎总意味着漏掉了某个约束（资源没有上限，忘记了需求上限等）。

{prf:ref}`ex-08-02-unbounded-ray`和{prf:ref}`ex-08-02-a3-q15`详细讨论了无界的 LP。

### 多重最优解

在极小化问题的最优表中，设某个非基变量 $x_j$ 有 $\bar c_j=0$。由 {eq}`eq-03-03-edge`，让 $x_j$ 进基不改变 $z$，所以到达的每个点也都是最优的：

- 若比值检验给出 $\theta>0$，主元变换得到另一个最优 BFS，两个顶点之间的整条棱都是最优的；
- 若 $\theta=0$（退化），主元变换改变了基，但没有改变点；
- 若 $\bar{\va}_j\le\vzero$，则射线 $\bar\vx+\theta\vd$，$\theta\ge0$ 由最优解组成，其中除 $\bar\vx$ 本身外都不是顶点（{prf:ref}`ex-08-02-a3-q14`）。

反之，若每个非基变量的 $\bar c_j$ 都严格为正（对极大化问题为严格为负），则最优解**唯一**：由 {eq}`eq-03-03-objective-reduced`，任何满足 $\vx_N\neq\vzero$ 的可行解 $\vx$ 的费用都比 $\bar\vx$ 大，而 $\vx_N=\vzero$ 迫使 $\vx_B=\bar{\vb}$。在退化的最优解处，其逆命题不成立：即使某个 $\bar c_j=0$，最优解也可能是唯一的。

(sec-03-03-degeneracy)=
## 退化与循环

### 退化的主元变换

当某个基变量为零时，BFS 是退化的。从几何上看，该顶点处的有效约束多于确定该顶点所需的约束个数：在平面上，有三条或更多条约束直线（包括坐标轴）通过该点，并且有多个基描述同一个点。在算法中，退化源于比值检验中的**并列**（tie）：当两行同时达到最小比值时，只有其中一行的基变量出基，另一行的基变量仍留在基中，取值为零。

在退化的 BFS 处，下一次比值检验可能给出 $\theta=0$：出基变量是一个已经为零的基变量，主元变换改变了基，但没有移动点，也没有改变 $z$。这样的**退化的主元变换**（degenerate pivot）并非白费。它可能是到达一个通过最优性检验的基所必需的，或者是找到一条可以继续移动的棱所必需的。{prf:ref}`ex-08-02-a3-q13`展示了这种*暂时退化*，{prf:ref}`ex-08-02-a3-q12`统计了可以描述同一个退化顶点的基的个数。

### 循环

{prf:ref}`thm-03-03-finite`的证明在退化顶点处失效：若 $\theta=0$，目标函数不下降，没有什么能阻止该方法回到已经访问过的基。一旦在确定性规则下发生这种情况，同样的序列就会永远重复下去。这就是**循环**。它在实际中很少见，但确实可能发生，下面的经典例子就说明了这一点；它是 Bertsimas 和 Tsitsiklis（1997, Example 3.6）给出的 Beale（1955）例子的版本。

:::{prf:example} Dantzig 规则下的循环
:label: ex-03-03-cycling

考虑

$$
\begin{aligned}
\min_{\vx}\ & z=-\tfrac34x_4+20x_5-\tfrac12x_6+6x_7\\
\st\ & x_1+\tfrac14x_4-8x_5-x_6+9x_7=0,\\
& x_2+\tfrac12x_4-12x_5-\tfrac12x_6+3x_7=0,\\
& x_3+x_6=1,\\
& x_1,\dots,x_7\ge0,
\end{aligned}
$$

它是关于 $\{x_1,x_2,x_3\}$ 的典式：

| 基变量 | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | 右端项 |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $0$ | $\boldsymbol{\tfrac{1}{4}}$ | $-8$ | $-1$ | $9$ | $0$ |
| $x_2$ | $0$ | $1$ | $0$ | $\tfrac{1}{2}$ | $-12$ | $-\tfrac{1}{2}$ | $3$ | $0$ |
| $x_3$ | $0$ | $0$ | $1$ | $0$ | $0$ | $1$ | $0$ | $1$ |
| $-z$ | $0$ | $0$ | $0$ | $-\tfrac{3}{4}$ | $20$ | $-\tfrac{1}{2}$ | $6$ | $0$ |

BFS $(0,0,1,0,0,0,0)$ 是退化的。使用 Dantzig 规则（最负的 $\bar c_j$ 进基），比值检验中出现并列时取出基变量下标最小者。证明单纯形法此时会发生循环。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

精确的主元变换给出如下序列：

```{list-table} 回到初始基的六次退化主元变换
:name: tab-03-03-cycling
:header-rows: 1

* - 迭代
  - 进基变量（$\bar c_s$）
  - 并列的行
  - 出基变量
  - 主元
  - 主元变换后的基
* - 1
  - $x_4$ ($-\tfrac34$)
  - $x_1,x_2$
  - $x_1$
  - $\tfrac14$
  - $\{x_4,x_2,x_3\}$
* - 2
  - $x_5$ ($-4$)
  - $x_2$
  - $x_2$
  - $4$
  - $\{x_4,x_5,x_3\}$
* - 3
  - $x_6$ ($-2$)
  - $x_4,x_5$
  - $x_4$
  - $8$
  - $\{x_6,x_5,x_3\}$
* - 4
  - $x_7$ ($-3$)
  - $x_5$
  - $x_5$
  - $\tfrac{3}{16}$
  - $\{x_6,x_7,x_3\}$
* - 5
  - $x_1$ ($-1$)
  - $x_6,x_7$
  - $x_6$
  - $2$
  - $\{x_1,x_7,x_3\}$
* - 6
  - $x_2$ ($-2$)
  - $x_7$
  - $x_7$
  - $\tfrac13$
  - $\{x_1,x_2,x_3\}$
```

每次主元变换都有 $\theta=0$，点和目标函数值 $z=0$ 始终不变，六次主元变换之后单纯形表又回到初始表。尽管该 LP 有最优解 $\vx^*=(\tfrac34,0,0,1,0,1,0)$，$z^*=-\tfrac54$，该方法却会永远重复这六次主元变换。
:::

### Bland 规则

通过仔细选取进基变量和出基变量可以防止循环。最简单的这种选法由 Bland（1977）提出。

:::{prf:algorithm} Bland 最小下标规则
:label: alg-03-03-bland

在{prf:ref}`alg-03-03-simplex`的第 2 步和第 4 步中：

1. **进基变量。** 在 $\bar c_j<0$（极大化：$\bar c_j>0$）的非基变量中，选取**下标最小**的 $j$。
2. **出基变量。** 在达到最小比值的各行中，选取其基变量**下标最小**的行。
:::

:::{prf:theorem} Bland（1977）
:label: thm-03-03-bland

采用 Bland 规则时，单纯形法永不循环；无论 LP 是否退化，它都在有限次迭代后停止，得到最优基或一条无界射线。
:::

其证明（Bland, 1977；Bertsimas & Tsitsiklis, 1997, §3.4）先假设出现一个循环，并考察在该循环中进基和出基的变量里下标最大的那个变量。比较它出基时的单纯形表与它再次进基时的单纯形表，可以在检验数的符号上导出矛盾。把 Bland 规则应用于{prf:ref}`ex-03-03-cycling`，前三次主元变换与前面相同。但在第四张表中，$x_1$（$\bar c_1=-2$）和 $x_7$（$\bar c_7=-3$）都是改善变量，Bland 规则取 $x_1$，而 Dantzig 规则取的是 $x_7$：

| 迭代 | 进基变量 | 出基变量 | $\theta$ | 主元变换后的 $z$ |
|:-:|:-:|:-:|:-:|:-:|
| 4 | $x_1$ | $x_5$ | $0$ | $0$ |
| 5 | $x_2$ | $x_3$ | $\tfrac12$ | $-\tfrac12$ |
| 6 | $x_4$ | $x_2$ | $1$ | $-\tfrac54$ |

第六张表

| 基变量 | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | 右端项 |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_6$ | $0$ | $0$ | $1$ | $0$ | $0$ | $1$ | $0$ | $1$ |
| $x_1$ | $1$ | $-\tfrac{1}{2}$ | $\tfrac{3}{4}$ | $0$ | $-2$ | $0$ | $\tfrac{15}{2}$ | $\tfrac{3}{4}$ |
| $x_4$ | $0$ | $2$ | $1$ | $1$ | $-24$ | $0$ | $6$ | $1$ |
| $-z$ | $0$ | $\tfrac{3}{2}$ | $\tfrac{5}{4}$ | $0$ | $2$ | $0$ | $\tfrac{21}{2}$ | $\tfrac{5}{4}$ |

的所有检验数都 $\ge0$，给出最优解 $\vx^*=(\tfrac34,0,0,1,0,1,0)$，$z^*=-\tfrac54$（已由 LP 求解器验证）。

Bland 规则很简单，但往往很慢，因为它不考虑检验数的大小。另有两种经典的补救办法。**字典序规则**（lexicographic rule）在比值检验出现并列时，把 $[\bar{\vb}\ \ \mB^{-1}]$ 的各整行除以主元列中相应的元素后按字典序比较；它等价于 Charnes（1952）的思想，即把右端项**摄动**（perturbation）为 $b_i+\varepsilon^i$（$\varepsilon>0$ 极小），使每个 BFS 都非退化（Dantzig, 1963, ch. 10；Bertsimas & Tsitsiklis, 1997, §3.4）。商业求解器对 $\vb$ 或变量界作微小的随机摄动，在浮点运算中具有同样的效果。

## 本章回顾

### 小结

- 相对于一个可行基，LP 可写成典式（{prf:ref}`def-03-03-canonical-form`）：$\vx_B=\bar{\vb}-\mB^{-1}\mN\vx_N$，$z=\vc_B^\T\bar{\vb}+\sum_{j\notin\mathcal B}\bar c_jx_j$，其中检验数 $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$，$\vpi^\T=\vc_B^\T\mB^{-1}$（{eq}`eq-03-03-objective-reduced`，{eq}`eq-03-03-reduced-cost`）。
- 对极小化问题，若所有 $\bar c_j\ge0$，则 BFS 是最优解（对极大化问题为所有 $\bar c_j\le0$）；在非退化的 BFS 处该检验也是必要的（{prf:ref}`thm-03-03-optimality`）。
- 改善变量 $x_s$ 沿棱 {eq}`eq-03-03-edge` 增大；最小比值检验 {eq}`eq-03-03-ratio-test` 给出最长可行步长和出基变量，新基与原基相邻，且 $z'=z+\theta\bar c_s$（{eq}`eq-03-03-z-change`）。
- 在单纯形表中（目标行放在最后，行标为 $-z$，存放 $\bar c_j$ 和 $-z$），一次高斯–若尔当主元变换 {eq}`eq-03-03-pivot` 完成整个迭代步（{prf:ref}`alg-03-03-simplex`，{numref}`tab-03-03-rules`）。
- 产品组合问题沿{numref}`fig-03-03-simplex-path`的路径经三次主元变换求得最优解 $(\tfrac{45}{7},\tfrac{30}{7})$，每周 $z^*=\tfrac{36000}{7}\approx\$5142.86$。
- 无退化时该方法有限终止（{prf:ref}`thm-03-03-finite`）。没有正元素的改善列证明了无界性；最优时非基变量检验数为零预示着多重最优解。
- 退化的主元变换（$\theta=0$）可能使该方法循环（{prf:ref}`ex-03-03-cycling`）；Bland 最小下标规则（{prf:ref}`alg-03-03-bland`，{prf:ref}`thm-03-03-bland`）以及字典序规则或摄动规则可以防止循环。

### 常见错误

- **符号检验用反。** 在本书的 $-z$ 行中，*极大化*问题当所有 $\bar c_j\le0$ 时最优，进基变量取 $\bar c_j$ 最*正*者；用 Taha 的 $z$ 行时符号相反。在同一个问题中混用两种约定，会导致使 $z$ 变差的主元变换。
- **对非正元素计算比值。** $\bar a_{is}\le0$ 的行不参加比值检验。负的比值或分母为零的比值都不限制步长。
- **按最大比值或最大主元选取出基变量。** 除达到最小比值的行之外，任何其他行都会给出不可行的基本解（{prf:ref}`ex-03-03-matrix-step`）。
- **主元变换时忘记目标行**，或用不同的倍数更新目标行。目标行只是消元步骤 {eq}`eq-03-03-pivot` 中的又一行而已。
- **从目标行读取目标函数值时弄错符号。** $-z$ 行的右端项是 $-z$。
- **因为某个 $\bar c_j$ 不满足检验就断定退化的单纯形表不是最优的**，或者因为检验严格成立就断定最优解唯一。这两个结论都需要非退化性。
- **舍入。** 单纯形表中的小数舍入误差积累很快；手算时应保持精确分数。

## 参考文献

- Beale, E. M. L. (1955). Cycling in the dual simplex algorithm. *Naval
  Research Logistics Quarterly, 2*(4), 269–275. (The cycling example.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Sections 3.1–3.3: optimality
  conditions, development of the simplex method, the full tableau and
  Example 3.6 on cycling; Section 3.4: lexicographic and Bland's rules;
  Section 3.7: computational efficiency.)
- Bland, R. G. (1977). New finite pivoting rules for the simplex method.
  *Mathematics of Operations Research, 2*(2), 103–107.
  <https://doi.org/10.1287/moor.2.2.103>
- Charnes, A. (1952). Optimality and degeneracy in linear programming.
  *Econometrica, 20*(2), 160–170. <https://doi.org/10.2307/1907845>
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Chapter 5: the simplex method; Chapter 10: resolving
  degeneracy.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 4: the simplex method, its
  algebra and tableau form, tie breaking.)
- Klee, V., & Minty, G. J. (1972). How good is the simplex algorithm? In
  O. Shisha (Ed.), *Inequalities III* (pp. 159–175). Academic Press.
  (Exponentially many pivots under Dantzig's rule.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Sections 3.1–3.3: equation form, from the graphical to the
  algebraic solution, the simplex method with the $z$-row convention;
  Section 3.5: degeneracy, alternative optima and unboundedness.)
