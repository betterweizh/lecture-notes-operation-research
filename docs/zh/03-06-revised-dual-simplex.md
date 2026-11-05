# 修正单纯形法与对偶单纯形法

## 概述

{doc}`03-03-simplex-method`中的表格方法在每次主元变换时都要重新计算单纯形表中的每一个元素，尽管一次迭代只需要检验数、一列和右端项。**修正单纯形法**（revised simplex method）保持原始数据 $\mA,\vb,\vc$ 不变，只携带当前基矩阵的逆 $\mB^{-1}$；其他一切都在需要时由它计算出来。本章把单纯形迭代写成矩阵形式，说明一次主元变换相当于用一个**初等矩阵（eta 矩阵）**（elementary (eta) matrix）左乘 $\mB^{-1}$，从而 $\mB^{-1}$ 是这类矩阵的乘积（**逆矩阵的乘积形式**，product form of the inverse），并给出修正算法、带加边逆矩阵 $\mT$ 的两阶段版本及其运算量。本章后半部分把单纯形法反过来：**对偶单纯形法**（dual simplex method）从检验数已达最优但基本解不可行的基出发，在保持检验数最优的同时通过主元变换走向可行。它是在增加约束或改变右端项之后进行重新优化的工具，这也是{doc}`03-07-sensitivity-analysis`依赖它的原因。本章完整地求解了三个例子：一次乘积形式的修正单纯形法求解、一次两阶段修正单纯形法求解和一次对偶单纯形主元变换。本章要用到{doc}`01-01-linear-algebra`中的分块矩阵和初等矩阵、{doc}`03-04-two-phase-big-m`中的两阶段法，以及{doc}`03-05-lp-duality`中的对偶线性规划问题（linear program, LP）与对偶可行性。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 用 $\mB^{-1}$ 和原始数据表示单纯形迭代中的各个量（$\bar\vb$、$\vpi$、$\bar c_j$、$\bar\mA_{\cdot s}$）；
- 构造一次主元变换的 eta 矩阵，并以乘积形式更新 $\mB^{-1}$；
- 手工执行修正单纯形法，包括带加边逆矩阵 $\mT$ 的两阶段版本；
- 比较表格方法与修正方法每次迭代的工作量；
- 推导对偶单纯形法的出基规则和进基规则，并以表格形式和修正形式执行该方法；
- 在增加约束后进行重新优化，并在没有现成的对偶可行基时借助人工界约束启动对偶单纯形法。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\mB_k$ | 第 $k$ 次迭代的基矩阵 |
| $\mE_k$ | 第 $k$ 次主元变换的 eta 矩阵（初等矩阵、主元变换矩阵）：$\mB_{k+1}^{-1}=\mE_k\mB_k^{-1}$ |
| $\veta$ | eta 列，即 $\mE_k$ 中与单位矩阵不同的那一列 |
| $\bar\mA_{\cdot s}=\mB^{-1}\mA_{\cdot s}$ | 进基变量 $x_s$ 的更新列 |
| $\bar\mA_{r\cdot}$, $\bar a_{rj}$ | $\mB^{-1}\mA$ 的第 $r$ 行（对偶单纯形法的主元行）及其元素 |
| $\ve_r$ | 第 $r$ 个单位向量 |
| $\hat\mB$ | 加边基矩阵 $\begin{pmatrix}\mB&\vzero\\\vc_B^\T&1\end{pmatrix}$；$\mT=\hat\mB^{-1}$ |
| $\theta_D$ | 对偶单纯形法的步长（对偶比值） |

(sec-03-06-why-revised)=
## 为什么需要修正方法

考虑标准形 LP

$$
\begin{aligned}
\min_{\vx}\ & \vc^\T\vx \\
\st\ & \mA\vx=\vb,\\
& \vx\ge\vzero,
\end{aligned}
$$ (eq-03-06-lp)

其中 $\mA\in\R^{m\times n}$ 的秩为 $m$。对基 $\mB$（基变量 $\vx_B$ 的列）和非基矩阵 $\mN$，{doc}`03-03-simplex-method`中的单纯形表恰好存放以下各量：

$$
\begin{aligned}
\bar\vb &= \mB^{-1}\vb && \text{基变量的值，}\\
\vpi^\T &= \vc_B^\T\mB^{-1} && \text{单纯形乘子（对偶解），}\\
\bar c_j &= c_j-\vpi^\T\mA_{\cdot j} && x_j\ \text{的检验数，}\\
\bar\mA_{\cdot j} &= \mB^{-1}\mA_{\cdot j} && x_j\ \text{的更新列，}\\
-z &= -\vc_B^\T\bar\vb = -\vpi^\T\vb && \text{目标元素。}
\end{aligned}
$$ (eq-03-06-matrix-quantities)

该方法的一次迭代（对极小化问题）只用到其中一部分：

1. 非基变量的检验数 $\bar c_j$，用于检验最优性并选取 $\bar c_s<0$ 的进基变量 $x_s$；
2. 单独一列 $\bar\mA_{\cdot s}$，用于比值检验 $\theta=\min\{\bar b_i/\bar a_{is}:\bar a_{is}>0\}$ 和无界性检验 $\bar\mA_{\cdot s}\le\vzero$；
3. 右端项 $\bar\vb$。

其余 $n-m-1$ 个更新列在每次主元变换中都被计算出来然后丢弃。有两点观察使得这种浪费在实际中很大。第一，计算经验表明，单纯形法通常在 $m$ 量级的主元变换次数内达到最优（常见的报告数字大约在 $m$ 到 $3m$ 之间），所以当 $n\gg m$ 时，$n$ 列中的大多数从不进基。第二，实际的约束矩阵是**稀疏的**（sparse，大多数 $a_{ij}$ 为零），而单纯形表会破坏稀疏性：即使 $\mA$ 是稀疏的，经过几次主元变换后 $\mB^{-1}\mA$ 也会变得稠密。

因此，**修正单纯形法**存储永不改变的原始数据和 $\mB^{-1}$ 的某种表示，只在迭代需要的地方计算 {eq}`eq-03-06-matrix-quantities`。剩下的问题只是在主元变换后如何低成本地更新 $\mB^{-1}$。

(sec-03-06-product-form)=
## 初等矩阵与逆矩阵的乘积形式

一次主元变换替换一个基列：进基列 $\mA_{\cdot s}$ 取代出基变量的列，后者位于基的第 $r$ 个位置。新的逆矩阵等于用一个与单位矩阵只有一列不同的矩阵左乘旧的逆矩阵。

:::{prf:definition} 主元变换的 eta 矩阵
:label: def-03-06-eta

设 $\bar\mA_{\cdot s}=\mB^{-1}\mA_{\cdot s}$，主元为 $\bar a_{rs}\ne0$。以 $(r,s)$ 为主元的主元变换的 **eta 矩阵**（eta matrix，也称主元变换矩阵、初等矩阵）$\mE$ 是把 $m\times m$ 单位矩阵的第 $r$ 列换成 **eta 列**（eta column）

$$
\veta=\left(-\frac{\bar a_{1s}}{\bar a_{rs}},\ \dots,\
-\frac{\bar a_{r-1,s}}{\bar a_{rs}},\ \frac{1}{\bar a_{rs}},\
-\frac{\bar a_{r+1,s}}{\bar a_{rs}},\ \dots,\
-\frac{\bar a_{ms}}{\bar a_{rs}}\right)^\T .
$$ (eq-03-06-eta)

所得的矩阵。
:::

用 $\mE$ 左乘向量 $\vv$，恰好完成单纯形表主元变换中的高斯–若尔当行变换：第 $r$ 行除以 $\bar a_{rs}$，其他每一行 $i$ 减去旧的第 $r$ 行的 $\bar a_{is}/\bar a_{rs}$ 倍，

$$
\begin{gathered}
(\mE\vv)_r=\frac{v_r}{\bar a_{rs}},\\
(\mE\vv)_i=v_i-\frac{\bar a_{is}}{\bar a_{rs}}\,v_r\quad(i\ne r).
\end{gathered}
$$ (eq-03-06-eta-action)

特别地，$\mE\bar\mA_{\cdot s}=\ve_r$：主元变换把进基列变成单位列。（这就是{prf:ref}`prop-01-01-pivot-matrix`中的主元变换矩阵 {eq}`eq-01-01-eta-matrix`，其中 $\vw=\bar\mA_{\cdot s}$；写成初等矩阵之积，它是 $m-1$ 次行倍加和一次行倍乘。）

:::{prf:theorem} 基逆矩阵的更新
:label: thm-03-06-update

设 $\mB'$ 是把 $\mB$ 的第 $r$ 列换成 $\mA_{\cdot s}$ 所得的矩阵，其中 $\bar a_{rs}\ne0$。则 $\mB'$ 非奇异，并且

$$
\mB'^{-1}=\mE\,\mB^{-1},
$$ (eq-03-06-update)

其中 $\mE$ 是以 $(r,s)$ 为主元的主元变换的 eta 矩阵 {eq}`eq-03-06-eta`。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

记 $\mF$ 为把单位矩阵的第 $r$ 列换成 $\bar\mA_{\cdot s}$ 所得的矩阵。$\mB\mF$ 的第 $i\ne r$ 列是 $\mB\ve_i$，即 $\mB$ 的第 $i$ 列，而第 $r$ 列是 $\mB\bar\mA_{\cdot s}=\mA_{\cdot s}$。因此 $\mB'=\mB\mF$。矩阵 $\mF$ 除第 $r$ 列外与单位矩阵相同，所以 $\det\mF=\bar a_{rs}\ne0$，$\mB'$ 非奇异。由 {eq}`eq-03-06-eta-action`，对 $i\ne r$ 有 $\mE\ve_i=\ve_i$，且 $\mE\bar\mA_{\cdot s}=\ve_r$，所以 $\mE\mF=\mI$，即 $\mE=\mF^{-1}$。因此 $\mB'^{-1}=\mF^{-1}\mB^{-1}=\mE\mB^{-1}$。
:::

若该方法从松弛基 $\mB_1=\mI$ 出发，则经过 $k-1$ 次主元变换后

$$
\mB_k^{-1}=\mE_{k-1}\mE_{k-2}\cdots\mE_1 .
$$ (eq-03-06-product-form)

这就是**逆矩阵的乘积形式**（Dantzig & Orchard-Hays, 1954），在 {eq}`eq-01-01-product-form` 中已经出现过（那里从第一次换基时的 $\mE_1$ 开始计数）。一个 eta 矩阵只需存储其位置 $r$ 和列 $\veta$，即 $m+1$ 个数；由于 $\bar\mA_{\cdot s}$ 是稀疏的，往往还更少。{eq}`eq-03-06-matrix-quantities` 中的乘积无须形成 $\mB^{-1}$ 就能计算：

* $\bar\mA_{\cdot s}=\mE_{k-1}(\cdots(\mE_1\mA_{\cdot s}))$，从第一个到最后一个依次作用各 eta 矩阵（由 {eq}`eq-03-06-eta-action`，每次作用约需 $m$ 次乘法）；
* $\vpi^\T=(\cdots(\vc_B^\T\mE_{k-1})\cdots)\mE_1$，从最后一个到第一个依次作用；行向量 $\vu^\T$ 右乘 $\mE$ 只改变其第 $r$ 个元素，该元素变为 $\vu^\T\veta$。

eta 矩阵的列表每次主元变换增加一个，所以这些乘积的计算量也随之增长，舍入误差也会积累。因此程序会定期**重新求逆**（reinversion）：重新分解当前的 $\mB$（如今采用 LU 分解，而不是 eta 矩阵之积），并重新开始列表。像下面这样的小型手算例子，则直接显式地携带 $\mB^{-1}$，并用 {eq}`eq-03-06-update` 更新它。

(sec-03-06-revised-algorithm)=
## 修正单纯形算法

:::{prf:algorithm} 修正单纯形法（极小化）
:label: alg-03-06-revised-simplex

**输入：** {eq}`eq-03-06-lp` 的数据 $\mA,\vb,\vc$，一个可行基 $\mB$ 及其 $\mB^{-1}$（或其 eta 因子），且 $\vx_B=\mB^{-1}\vb\ge\vzero$。
**输出：** 一个最优基，或无界性的证明。

1. **定价。** 计算乘子 $\vpi^\T=\vc_B^\T\mB^{-1}$ 和非基变量的检验数 $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$。若它们都满足 $\bar c_j\ge0$，停止：$\vx_B=\mB^{-1}\vb$，$\vx_N=\vzero$ 是最优解，且 $z^*=\vpi^\T\vb$。
2. **选取进基变量。** 选取 $\bar c_s<0$ 的 $x_s$（例如取最负者），并计算 $\bar\mA_{\cdot s}=\mB^{-1}\mA_{\cdot s}$。
3. **比值检验。** 若 $\bar\mA_{\cdot s}\le\vzero$，停止：LP 沿增大 $x_s$ 的射线无界。否则，以 $\bar\vb=\vx_B$ 计算 $\theta=\min\{\bar b_i/\bar a_{is}:\bar a_{is}>0\}$，并令 $r$ 为达到最小值的一行。第 $r$ 行的基变量出基。
4. **更新。** 按 {eq}`eq-03-06-eta` 构造以 $(r,s)$ 为主元的 eta 矩阵 $\mE$，令 $\mB^{-1}\leftarrow\mE\mB^{-1}$，$\vx_B\leftarrow\mE\vx_B$（新值为：第 $r$ 行的 $x_s$ 取 $\theta$，其他各行取 $\bar b_i-\theta\bar a_{is}$），把第 $r$ 行的基变量换成 $x_s$，转第 1 步。
:::

```{figure} figures/03-06-revised-iteration.svg
:name: fig-03-06-revised-iteration
:width: 560px
:class: dark-light
:alt: 由四个方框组成的循环：定价、选取进基变量、比值检验、更新逆矩阵，并有通往最优和无界两种情形的出口

修正单纯形法的一次迭代（{prf:ref}`alg-03-06-revised-simplex`）。从一次迭代传到下一次迭代的只有 $\mB^{-1}$（或其 eta 因子）；$\mA$、$\vb$ 和 $\vc$ 从不改变。
```

{numref}`fig-03-06-revised-iteration`把四个步骤画成一个循环。对极大化问题，最优性检验为 $\bar c_j\le0$，进基变量满足 $\bar c_s>0$；其他步骤不变。当采用相同的进基和出基规则时，该方法访问的基与表格方法相同，所以有限性、退化和 Bland 规则（{prf:ref}`alg-03-03-bland`）的结论都原样适用。手算时，用一张列为“基变量 $\mid\mB^{-1}\mid\vx_B\mid\bar\mA_{\cdot s}\mid$ 比值”的工作表代替完整的单纯形表。

:::{prf:example} 乘积形式的修正单纯形法求解
:label: ex-03-06-revised-run

从松弛基 $(x_4,x_5,x_6)$ 出发，用修正单纯形法求解

$$
\begin{aligned}
\min_{\vx}\ & -3x_1-x_2-3x_3 \\
\st\ & 2x_1+x_2+x_3+x_4=2,\\
& x_1+2x_2+3x_3+x_5=5,\\
& 2x_1+2x_2+x_3+x_6=6,\\
& x_1,\dots,x_6\ge0
\end{aligned}
$$

并把每个 $\mB_k^{-1}$ 写成 eta 矩阵之积。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

数据为 $\mA_{\cdot1}=(2,1,2)^\T$、$\mA_{\cdot2}=(1,2,2)^\T$、$\mA_{\cdot3}=(1,3,1)^\T$，松弛变量的列 $\ve_1,\ve_2,\ve_3$，$\vb=(2,5,6)^\T$ 以及 $\vc=(-3,-1,-3,0,0,0)^\T$。下面的运算都是精确的。

**第 1 次迭代。** $\mB_1=\mI$，所以 $\mB_1^{-1}=\mI$，$\vx_B=(2,5,6)^\T$，$z=0$，$\vpi^\T=(0,0,0)$。检验数为 $(\bar c_1,\bar c_2,\bar c_3)=(-3,-1,-3)$。任何检验数为负的变量都可以进基；我们让 $x_2$ 进基，因为它的主元为 1，可以使手算保持为整数（按最负规则会选 $x_1$；见本解答后面的说明）。

| 基变量 | $\mB_1^{-1}$ | | | $\vx_B$ | $\bar\mA_{\cdot2}$ | 比值 |
|---|---|---|---|---|---|---|
| $x_4$ | 1 | 0 | 0 | 2 | **1** | 2 |
| $x_5$ | 0 | 1 | 0 | 5 | 2 | 5/2 |
| $x_6$ | 0 | 0 | 1 | 6 | 2 | 3 |

最小比值为第 1 行的 2，所以 $x_4$ 出基。由 {eq}`eq-03-06-eta`，eta 列为 $(1/1,-2/1,-2/1)^\T$，位于第 1 个位置：

$$
\begin{gathered}
\mE_1=\begin{pmatrix}1&0&0\\-2&1&0\\-2&0&1\end{pmatrix},\\
\mB_2^{-1}=\mE_1\mB_1^{-1}=\mE_1 .
\end{gathered}
$$

**第 2 次迭代。** 基为 $(x_2,x_5,x_6)$，$\vx_B=\mB_2^{-1}\vb=(2,1,2)^\T$，$z=-2$。定价：$\vpi^\T=\vc_B^\T\mB_2^{-1}=(-1,0,0)\mB_2^{-1}=(-1,0,0)$，所以 $\bar c_1=-3-(-1)(2)=-1$，$\bar c_3=-3-(-1)(1)=-2$，$\bar c_4=0-(-1)(1)=1$。最负的是 $\bar c_3$；$x_3$ 进基，$\bar\mA_{\cdot3}=\mB_2^{-1}(1,3,1)^\T=(1,1,-1)^\T$。

| 基变量 | $\mB_2^{-1}$ | | | $\vx_B$ | $\bar\mA_{\cdot3}$ | 比值 |
|---|---|---|---|---|---|---|
| $x_2$ | 1 | 0 | 0 | 2 | 1 | 2 |
| $x_5$ | −2 | 1 | 0 | 1 | **1** | 1 |
| $x_6$ | −2 | 0 | 1 | 2 | −1 | — |

第 2 行胜出，$x_5$ 出基，eta 列 $(-1,1,1)^\T$ 位于第 2 个位置：

$$
\begin{gathered}
\mE_2=\begin{pmatrix}1&-1&0\\0&1&0\\0&1&1\end{pmatrix},\\
\mB_3^{-1}=\mE_2\mE_1=\begin{pmatrix}3&-1&0\\-2&1&0\\-4&1&1\end{pmatrix}.
\end{gathered}
$$

**第 3 次迭代。** 基为 $(x_2,x_3,x_6)$，$\vx_B=(1,1,3)^\T$，$z=-4$。定价：$\vpi^\T=(-1,-3,0)\mB_3^{-1}=(3,-2,0)$，所以 $\bar c_1=-3-(6-2)=-7$，$\bar c_4=0-3=-3$，$\bar c_5=0-(-2)=2$。现在 $x_1$ 进基，$\bar\mA_{\cdot1}=\mB_3^{-1}(2,1,2)^\T=(5,-3,-5)^\T$。

| 基变量 | $\mB_3^{-1}$ | | | $\vx_B$ | $\bar\mA_{\cdot1}$ | 比值 |
|---|---|---|---|---|---|---|
| $x_2$ | 3 | −1 | 0 | 1 | **5** | 1/5 |
| $x_3$ | −2 | 1 | 0 | 1 | −3 | — |
| $x_6$ | −4 | 1 | 1 | 3 | −5 | — |

只有第 1 行的元素为正，所以 $x_2$ 以 $\theta=1/5$ 出基。eta 列为 $(1/5,\,3/5,\,5/5)^\T$，位于第 1 个位置：

$$
\begin{gathered}
\mE_3=\begin{pmatrix}1/5&0&0\\3/5&1&0\\1&0&1\end{pmatrix},\\
\mB_4^{-1}=\mE_3\mE_2\mE_1=\begin{pmatrix}3/5&-1/5&0\\-1/5&2/5&0\\-1&0&1\end{pmatrix}.
\end{gathered}
$$

**第 4 次迭代。** 基为 $(x_1,x_3,x_6)$，$\vx_B=\mB_4^{-1}\vb=(1/5,8/5,4)^\T$，$z=-3\cdot\tfrac15-3\cdot\tfrac85=-\tfrac{27}{5}$。定价：$\vpi^\T=(-3,-3,0)\mB_4^{-1}=(-6/5,-3/5,0)$，且

$$
\begin{gathered}
\bar c_2=-1-\left(-\tfrac65-\tfrac65\right)=\tfrac75,\\
\bar c_4=\tfrac65,\\
\bar c_5=\tfrac35 .
\end{gathered}
$$

所有检验数都非负，所以该基是最优的：

$$
\begin{gathered}
\vx^*=\left(\tfrac15,\,0,\,\tfrac85,\,0,\,0,\,4\right)^\T,\\
z^*=-\tfrac{27}{5}=-5.4 .
\end{gathered}
$$

**验算。** 乘积 $\mB_4^{-1}\mB_4$（其中 $\mB_4=(\mA_{\cdot1},\mA_{\cdot3},\mA_{\cdot6})$）是单位矩阵，对偶目标函数值 $\vpi^\T\vb=-\tfrac65\cdot2-\tfrac35\cdot5=-\tfrac{27}{5}$ 等于 $z^*$，正如强对偶性（{prf:ref}`thm-03-05-strong-duality`）所要求的。求解器（HiGHS）给出同样的解。作为比较，{numref}`tab-03-06-final-tableau`是表格方法会得到的完整最终表 $\mB_4^{-1}\mA$；修正方法从未形成这张表，但它的检验数行和右端项恰好就是第 1 步所算出的。

```{table} 修正单纯形法求解的最终表，为便于比较按 $\mB_4^{-1}\mA$ 和 $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ 计算。
:name: tab-03-06-final-tableau

| 基变量 | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | 右端项 |
|---|---|---|---|---|---|---|---|
| $x_1$ | 1 | 1/5 | 0 | 3/5 | −1/5 | 0 | 1/5 |
| $x_3$ | 0 | 3/5 | 1 | −1/5 | 2/5 | 0 | 8/5 |
| $x_6$ | 0 | 1 | 0 | −1 | 0 | 1 | 4 |
| $-z$ | 0 | 7/5 | 0 | 6/5 | 3/5 | 0 | 27/5 |
```

*（本例的较早版本把 $\mE_3$ 的第三行印成 $(-1,0,1)$，把 $\mB_4^{-1}$ 的第三行印成 $(-1,1,1)$，把最终表的两处元素印成 $-1/5$ 和 $(0,3,0,-1,1,1)$，并且给出的解有五个分量。上面的值才是正确的：$\mE_3$ 第三行的元素为 $-\bar a_{31}/\bar a_{11}=-(-5)/5=1$；若使用印错的 $\mB_4^{-1}$，$x_6$ 的值将是 9 而不是 4。）*
:::

{prf:ref}`ex-03-06-revised-run`的解答在第 1 次迭代中有意违背了最负规则。若采用 Dantzig 规则（先选 $x_1$，并列时取较小下标），同一个 LP 只需两次主元变换即可求解：先 $x_1$ 取代 $x_4$，再 $x_3$ 取代 $x_5$，最终得到同一个基 $(x_1,x_3,x_6)$。没有哪一种规则在每个问题上都是最好的。

(sec-03-06-operation-counts)=
## 运算量

对稠密数据，统计每次迭代的乘法次数，这是主要的计算开销。

* **表格方法。** 一次主元变换要更新单纯形表的全部 $(m+1)(n+1)$ 个元素，每个元素约一次乘法：大约 $mn$ 次。
* **显式存储 $\mB^{-1}$ 的修正方法。** 定价需要计算 $\vpi^\T=\vc_B^\T\mB^{-1}$（$m^2$ 次）以及 $n-m$ 个非基列的 $\vpi^\T\mA_{\cdot j}$（$m(n-m)$ 次）；进基列 $\mB^{-1}\mA_{\cdot s}$ 需要 $m^2$ 次；比值检验需要 $m$ 次；按 {eq}`eq-03-06-eta-action` 更新 $\mE\mB^{-1}$ 约需 $m^2$ 次。总计约 $3m^2+m(n-m)$ 次。

对稠密矩阵，两者是同一量级的：$m=3$，$n=6$ 时（{prf:ref}`ex-03-06-revised-run`）约为 28 次对 36 次乘法；$m=100$，$n=1000$ 时约为 $1.0\times10^5$ 次对 $1.2\times10^5$ 次。修正方法胜出有三个原因，而稠密情形下的计数掩盖了这些原因。

1. **稀疏性。** 定价作用于原始的 $\mA$。若其元素中只有一小部分非零，定价的开销也只有 $m(n-m)$ 的那一小部分。当 $m=1000$，$n=10\,000$，非零元素占 1% 时，表格方法的一次主元变换约需 $10^7$ 次乘法，修正方法的一次迭代约需 $3.1\times10^6$ 次；若再对 $\mB$ 采用稀疏分解，$m^2$ 项也会缩小。
2. **部分定价**（partial pricing）。第 1 步可以在找到最先几个有吸引力的列后就停止，而不必对全部 $n-m$ 列定价；表格方法无法利用这一点，因为它无论如何都要更新每一列。
3. **存储与精度。** 修正方法存储（稀疏的）$\mA$ 加上 $\mB^{-1}$ 或其因子，而不是稠密的 $m\times n$ 单纯形表；而且每个量都由原始数据重新计算，所以舍入误差不会在各列中积累。

由于这些原因，所有实用的单纯形程序都是修正单纯形程序。单纯形表仍然是学习和手算小型例子最清楚的方式。

(sec-03-06-two-phase-revised)=
## 两阶段修正单纯形法与加边逆矩阵

当没有现成的松弛基时，两阶段法（{prf:ref}`alg-03-04-two-phase`）的第一阶段极小化人工变量 $R_i$ 之和 $\phi$。在修正形式中，把原目标函数作为一个额外的约束行保留下来是很方便的，

$$
\vc^\T\vx+(-z)=0,
$$ (eq-03-06-cost-row)

其中 $-z$ 被看作一个始终在基中的自由变量。这个扩大后的方程组的基是**加边基矩阵**（bordered basis matrix）

$$
\begin{aligned}
\hat\mB&=\begin{pmatrix}\mB&\vzero\\\vc_B^\T&1\end{pmatrix},\\
\mT=\hat\mB^{-1}&=\begin{pmatrix}\mB^{-1}&\vzero\\-\vc_B^\T\mB^{-1}&1\end{pmatrix}
=\begin{pmatrix}\mB^{-1}&\vzero\\-\vpi^\T&1\end{pmatrix}.
\end{aligned}
$$ (eq-03-06-bordered-inverse)

$\mT$ 的公式可以通过乘法验证：$\hat\mB\mT=\begin{pmatrix}\mB\mB^{-1}&\vzero\\\vc_B^\T\mB^{-1}-\vpi^\T&1\end{pmatrix}=\mI$。把 $\mT$ 作用于一个变量的扩展列和扩展右端项，就得到单纯形表的一整列及其检验数：

$$
\begin{aligned}
\mT\begin{pmatrix}\mA_{\cdot j}\\c_j\end{pmatrix}
&=\begin{pmatrix}\bar\mA_{\cdot j}\\\bar c_j\end{pmatrix},\\
\mT\begin{pmatrix}\vb\\0\end{pmatrix}
&=\begin{pmatrix}\bar\vb\\-z\end{pmatrix}.
\end{aligned}
$$ (eq-03-06-T-column)

所以 $\mT$ 的最后一行存放 $-\vpi^\T$，而 $\mT(\vb;0)$ 的最后一个元素就是单纯形表的目标元素 $-z$。矩阵 $\mT$（**T 矩阵**，T-matrix）与 $\mB^{-1}$ 完全一样，用 $(m+1)\times(m+1)$ 的 eta 矩阵更新；多出的一行只是主元变换中的又一行而已。

在第一阶段中，基由人工变量和 $-z$ 组成。定价使用第一阶段的费用：每个人工变量为 1，其他一切都为 0，*包括 $-z$ 的列*。当 $\phi=0$ 时第一阶段结束。若此时（通常如此）没有人工变量在基中，则基列都是原 $\mA$ 的列，其取值满足 $\mA\vx=\vb$，$\vx\ge\vzero$：它们构成原 LP 的一个基本可行解。（取值为 0 但仍在基中的人工变量，按{prf:ref}`thm-03-04-phase-one`那样被移出基，或者把其所在行作为冗余行删去。）由于 $-z$ 始终留在基中，当前的 $\mT$ *就是*这个基的加边逆矩阵：它的左上块是 $\mB^{-1}$，最后一行是对应于原始费用的 $(-\vpi^\T,1)$。第二阶段可以立即开始，无须重新计算任何东西。

:::{prf:example} 两阶段修正单纯形法
:label: ex-03-06-two-phase-revised

用带加边逆矩阵 $\mT$ 的两阶段修正单纯形法求解

$$
\begin{aligned}
\min_{\vx}\ & 3x_1-x_2-7x_3+3x_4+x_5 \\
\st\ & 5x_1-4x_2+13x_3-2x_4+x_5=20,\\
& x_1-x_2+5x_3-x_4+x_5=8,\\
& x_1,\dots,x_5\ge0
\end{aligned}
$$
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**准备。** 对两行分别加入人工变量 $x_6,x_7$，并加入费用行 {eq}`eq-03-06-cost-row`。第一阶段极小化 $\phi=x_6+x_7$。扩展列 $(\mA_{\cdot j};c_j)$ 为

| | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | $-z$ | 右端项 |
|---|---|---|---|---|---|---|---|---|---|
| 第 1 行 | 5 | −4 | 13 | −2 | 1 | 1 | 0 | 0 | 20 |
| 第 2 行 | 1 | −1 | 5 | −1 | 1 | 0 | 1 | 0 | 8 |
| 费用行 | 3 | −1 | −7 | 3 | 1 | 0 | 0 | 1 | 0 |
| 第一阶段费用 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | |

初始基为 $(x_6,x_7,-z)$，$\mT_1=\mI$。

**第一阶段，第 1 次迭代。** 第一阶段乘子为 $(1,1,0)\mT_1=(1,1,0)$（第三个元素是 $-z$ 的第一阶段费用 0），$x_1,\dots,x_5$ 的第一阶段检验数 $0-(1,1,0)(\mA_{\cdot j};c_j)$ 为 $(-6,\,5,\,-18,\,3,\,-2)$。最负的是 $-18$：$x_3$ 进基，扩展列为 $\mT_1(13,5,-7)^\T=(13,5,-7)^\T$。

| 基变量 | $\mT_1$ | | | 取值 | $\bar\mA_{\cdot3}$ | 比值 |
|---|---|---|---|---|---|---|
| $x_6$ | 1 | 0 | 0 | 20 | **13** | 20/13 |
| $x_7$ | 0 | 1 | 0 | 8 | 5 | 8/5 |
| $-z$ | 0 | 0 | 1 | 0 | −7 | — |

$-z$ 所在行从不参加比值检验（$-z$ 是自由变量）。由于 $20/13<8/5$，$x_6$ 出基，且

$$
\begin{gathered}
\mE_1=\begin{pmatrix}1/13&0&0\\-5/13&1&0\\7/13&0&1\end{pmatrix},\\
\mT_2=\mE_1 .
\end{gathered}
$$

**第一阶段，第 2 次迭代。** 取值为 $\mT_2(20,8,0)^\T=(20/13,\,4/13,\,140/13)^\T$，所以 $\phi=4/13$。乘子为 $(0,1,0)\mT_2=(-5/13,\,1,\,0)$；$x_1,x_2,x_4,x_5,x_6$ 的检验数为 $(12/13,\,-7/13,\,3/13,\,-8/13,\,18/13)$。最负的是 $-8/13$：$x_5$ 进基，其列为 $\mT_2(1,1,1)^\T=(1/13,\,8/13,\,20/13)^\T$。

| 基变量 | $\mT_2$ | | | 取值 | $\bar\mA_{\cdot5}$ | 比值 |
|---|---|---|---|---|---|---|
| $x_3$ | 1/13 | 0 | 0 | 20/13 | 1/13 | 20 |
| $x_7$ | −5/13 | 1 | 0 | 4/13 | **8/13** | 1/2 |
| $-z$ | 7/13 | 0 | 1 | 140/13 | 20/13 | — |

$x_7$ 出基，且

$$
\begin{gathered}
\mE_2=\begin{pmatrix}1&-1/8&0\\0&13/8&0\\0&-5/2&1\end{pmatrix},\\
\mT_3=\mE_2\mT_2=\begin{pmatrix}1/8&-1/8&0\\-5/8&13/8&0\\3/2&-5/2&1\end{pmatrix}.
\end{gathered}
$$

**第一阶段结束。** $(x_3,x_5,-z)$ 的取值为 $\mT_3(20,8,0)^\T=(3/2,\,1/2,\,10)^\T$，所以 $\phi=0$，两个人工变量都是非基变量：现在第一阶段检验数为 $(0,0,0,0,0,1,1)$，而 $(x_3,x_5)=(3/2,1/2)$ 是原 LP 的一个基本可行解，$z=-10$。“第二阶段的 $\mB^{-1}$ 是什么？”这个问题的答案可以从 $\mT_3$ 中读出：它的左上块是基 $(x_3,x_5)$ 的 $\mB^{-1}$，最后一行是对应于原始费用的 $-\vpi^\T$，

$$
\begin{gathered}
\mB^{-1}=\begin{pmatrix}1/8&-1/8\\-5/8&13/8\end{pmatrix},\\
\vpi^\T=\left(-\tfrac32,\ \tfrac52\right).
\end{gathered}
$$

从此以后人工变量的列都被删去。

**第二阶段，第 3 次迭代。** 由 $\mA_{\cdot1}=(5,1)^\T$、$\mA_{\cdot2}=(-4,-1)^\T$、$\mA_{\cdot4}=(-2,-1)^\T$，$x_1,x_2,x_4$ 的检验数为

$$
\begin{aligned}
(\bar c_1,\bar c_2,\bar c_4)&=(3,-1,3)-\left(-\tfrac32,\tfrac52\right)
\begin{pmatrix}5&-4&-2\\1&-1&-1\end{pmatrix}\\
&=\left(8,\,-\tfrac92,\,\tfrac52\right).
\end{aligned}
$$

$x_2$ 进基，$\mT_3(-4,-1,-1)^\T=(-3/8,\,7/8,\,-9/2)^\T$；正如 {eq}`eq-03-06-T-column` 所保证的，最后一个元素就是 $\bar c_2$。

| 基变量 | $\mB^{-1}$ | | 取值 | $\bar\mA_{\cdot2}$ | 比值 |
|---|---|---|---|---|---|
| $x_3$ | 1/8 | −1/8 | 3/2 | −3/8 | — |
| $x_5$ | −5/8 | 13/8 | 1/2 | **7/8** | 4/7 |

$x_5$ 以 $\theta=4/7$ 出基。在 $2\times2$ 块上，eta 矩阵为 $\begin{pmatrix}1&3/7\\0&8/7\end{pmatrix}$；在加边的 $3\times3$ 逆矩阵上，它是同一个矩阵再加上额外的一行 $(0,\,36/7,\,1)$，并且

$$
\begin{gathered}
\mT_4=\begin{pmatrix}-1/7&4/7&0\\-5/7&13/7&0\\-12/7&41/7&1\end{pmatrix},\\
\mT_4\begin{pmatrix}20\\8\\0\end{pmatrix}=\begin{pmatrix}12/7\\4/7\\88/7\end{pmatrix}.
\end{gathered}
$$

**第二阶段，第 4 次迭代。** 基为 $(x_3,x_2)$，取值为 $(12/7,\,4/7)$，$-z=88/7$。由 $\mT_4$ 的最后一行得 $\vpi^\T=(12/7,\,-41/7)$，且

$$
(\bar c_1,\bar c_4,\bar c_5)=\left(\tfrac27,\ \tfrac47,\ \tfrac{36}{7}\right)\ge0 .
$$

该基是最优的：

$$
\begin{gathered}
\vx^*=\left(0,\,\tfrac47,\,\tfrac{12}{7},\,0,\,0\right)^\T,\\
z^*=-\tfrac{88}{7}\approx-12.57 .
\end{gathered}
$$

验算：$5\cdot0-4\cdot\tfrac47+13\cdot\tfrac{12}{7}=\tfrac{140}{7}=20$，$-\tfrac47+5\cdot\tfrac{12}{7}=8$；目标函数值为 $-\tfrac47-7\cdot\tfrac{12}{7}=-\tfrac{88}{7}$，且 $\vpi^\T\vb=\tfrac{12}{7}\cdot20-\tfrac{41}{7}\cdot8=-\tfrac{88}{7}$。HiGHS 给出同样的最优解。

*（本例的较早版本把初始基变量列为 $x_5,x_6,-z$ 而不是 $x_6,x_7,-z$，在第二阶段定价中对 $x_4$ 使用了费用 $-7$ 和列 $(2,-1)^\T$，并且在没有给出最优解的情况下就结束了。$x_4$ 的费用是 3，其列为 $(-2,-1)^\T$；检验数 $(8,-9/2,5/2)$ 对这些数据是正确的。）*
:::

(sec-03-06-dual-simplex)=
## 对偶单纯形法

当 $\bar\vb=\mB^{-1}\vb\ge\vzero$ 时，基是**原始可行**（primal feasible）的；当所有检验数都具有最优符号（对极小化问题为 $\bar c_j\ge0$）时，基是**对偶可行**（dual feasible）的。对偶可行意味着 $\vpi^\T=\vc_B^\T\mB^{-1}$ 对每个 $j$ 满足 $\vpi^\T\mA_{\cdot j}\le c_j$，即 $\vpi$ 是对偶 LP $\max\{\vb^\T\vpi:\mA^\T\vpi\le\vc\}$ 的可行解（{prf:ref}`def-03-05-dual-feasible-basis`）。两者兼备的基是最优的（{prf:ref}`prop-03-05-basis-optimal`）。原始单纯形法保持原始可行性并朝着对偶可行性推进；**对偶单纯形法**（Lemke, 1954）则相反：

* 它从一个对偶可行但原始不可行（某个 $\bar b_r<0$）的基出发；
* 它通过主元变换使基保持对偶可行；
* 当基变为原始可行（因而最优）时，或者当某一行证明该 LP 没有可行解时，它停止。

这样的初始基是自然出现的：费用非负的 $\ge$ 约束给出右端项为负的松弛基（{prf:ref}`ex-08-04-dual-infeasible`）；而在改变右端项或增加约束之后，原来的最优基保持其检验数不变，但可能失去可行性（{doc}`03-07-sensitivity-analysis`）。

### 主元变换规则的推导

取一张对偶可行的单纯形表，其第 $r$ 行有 $\bar b_r<0$，该行的基变量 $x_{B(r)}$ 必须出基（它是负的）。记 $\mathcal B$ 为基变量的下标集，第 $r$ 行为

$$
x_{B(r)}+\sum_{j\notin\mathcal B}\bar a_{rj}x_j=\bar b_r<0 .
$$ (eq-03-06-pivot-row)

**不可行性检验。** 若对每个非基变量 $j$ 都有 $\bar a_{rj}\ge0$，则对每个 $\vx\ge\vzero$，{eq}`eq-03-06-pivot-row` 的左端都非负，而右端为负。因此不存在可行的 $\vx$。

**进基变量。** 否则，以某个 $\bar a_{rs}<0$ 为主元进行主元变换；需要负的主元，才能使 $x_s$ 的新值 $\bar b_r/\bar a_{rs}$ 为正。主元变换从 $-z$ 行中减去第 $r$ 行的 $\bar c_s/\bar a_{rs}$ 倍。记 $\theta_D=-\bar c_s/\bar a_{rs}\ge0$，新的检验数为

$$
\bar c_j'=\bar c_j+\theta_D\,\bar a_{rj}\qquad(j\notin\mathcal B),
$$ (eq-03-06-dual-update)

而出基变量得到 $\bar c'=\theta_D$。对 $\bar a_{rj}\ge0$，检验数不会减小。对 $\bar a_{rj}<0$，检验数保持非负当且仅当 $\theta_D\le-\bar c_j/\bar a_{rj}$。因此保持对偶可行性当且仅当 $\theta_D$ 是**对偶比值**（dual ratio）

$$
\theta_D=\min\left\{\frac{\bar c_j}{-\bar a_{rj}}:\ j\notin\mathcal B,\ \bar a_{rj}<0\right\},
$$ (eq-03-06-dual-ratio)

且进基变量 $x_s$ 是达到该最小值的变量。（步长再大，就会使达到最小值的那个变量的检验数变为负。）

**目标函数。** 目标元素从 $-z$ 变为 $-z+\theta_D\bar b_r$，所以 $z$ 增加 $\theta_D\abs{\bar b_r}\ge0$。由于 $z=\vc_B^\T\bar\vb=\vpi^\T\vb$ 是对偶可行解 $\vpi$ 的对偶目标函数值，该方法使对偶目标函数值单调增加：它就是在原始单纯形表上执行的、应用于对偶 LP 的原始单纯形法。若每次主元变换都有 $\theta_D>0$（对偶非退化），则没有基会重复，该方法在有限次主元变换后停止。

对极大化问题，对偶可行性为 $\bar c_j\le0$，同样的论证给出比值 $\min\{\abs{\bar c_j/\bar a_{rj}}:\bar a_{rj}<0\}$；取绝对值的形式涵盖了两种情形。

:::{prf:algorithm} 对偶单纯形法（极小化）
:label: alg-03-06-dual-simplex

**输入：** 一个对偶可行的基 $\mB$，即对所有非基变量 $j$ 有 $\bar c_j\ge0$。
**输出：** 一个最优基，或该 LP 不可行的证明。

1. **最优性检验。** 若 $\bar\vb=\mB^{-1}\vb\ge\vzero$，停止：该基是最优的。
2. **出基变量。** 选取一个 $\bar b_r<0$ 的行 $r$（通常取最负者）。其基变量出基。
3. **不可行性检验。** 计算 $\mB^{-1}\mA$ 的第 $r$ 行，即对非基变量 $j$ 计算 $\bar a_{rj}=(\ve_r^\T\mB^{-1})\mA_{\cdot j}$。若所有 $\bar a_{rj}\ge0$，停止：该 LP 没有可行解。
4. **进基变量。** 选取达到对偶比值 {eq}`eq-03-06-dual-ratio` 的 $x_s$。
5. **主元变换。** 以 $\bar a_{rs}$ 为主元进行主元变换（表格形式）；或者在修正形式中计算 $\bar\mA_{\cdot s}=\mB^{-1}\mA_{\cdot s}$，并用以 $(r,s)$ 为主元的 eta 矩阵令 $\mB^{-1}\leftarrow\mE\mB^{-1}$。转第 1 步。
:::

在修正形式中，第 3 步只需要行向量 $\ve_r^\T\mB^{-1}$（逆矩阵的第 $r$ 行）以及对每个非基列计算一次内积；第 5 步只需要一列。乘子 $\vpi$ 按{prf:ref}`alg-03-06-revised-simplex`的方式更新或重新计算，以得到检验数。完整的修正形式求解见{prf:ref}`ex-08-04-a5-q4-revised-dual`。

:::{prf:example} 一次对偶单纯形主元变换
:label: ex-03-06-dual-simplex

某极小化问题有如下典式单纯形表，基变量为 $x_1,x_2,x_3$，目标函数值为 0。对它进行重新优化。

| 基变量 | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | 右端项 |
|---|---|---|---|---|---|---|---|
| $x_1$ | 1 | 0 | 0 | 4 | −5 | 7 | 8 |
| $x_2$ | 0 | 1 | 0 | −2 | 4 | −2 | −2 |
| $x_3$ | 0 | 0 | 1 | 1 | −3 | 2 | 2 |
| $-z$ | 0 | 0 | 0 | 1 | 3 | 2 | 0 |
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

非基变量的检验数 $(1,3,2)$ 都非负，所以该基是对偶可行的，但 $x_2=-2<0$：该基原始不可行，适用对偶单纯形法。

*出基变量。* 第 2 行是唯一右端项为负的行，所以 $x_2$ 出基。

*进基变量。* 第 2 行的负元素为 $\bar a_{24}=-2$ 和 $\bar a_{26}=-2$。把第 2 行的 $\theta_D$ 倍加到 $-z$ 行上，$x_4$ 和 $x_6$ 的检验数分别变为 $1-2\theta_D$ 和 $2-2\theta_D$；两者都保持非负，仅当 $\theta_D\le\min\{1/2,\,1\}=1/2$。这就是对偶比值 {eq}`eq-03-06-dual-ratio`：$x_4$ 为 $1/2$，$x_6$ 为 $2/2=1$。所以 $x_4$ 进基，主元为 $\bar a_{24}=-2$。

*主元变换。*

| 基变量 | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | 右端项 |
|---|---|---|---|---|---|---|---|
| $x_1$ | 1 | 2 | 0 | 0 | 3 | 3 | 4 |
| $x_4$ | 0 | −1/2 | 0 | 1 | −2 | 1 | 1 |
| $x_3$ | 0 | 1/2 | 1 | 0 | −1 | 1 | 1 |
| $-z$ | 0 | 1/2 | 0 | 0 | 5 | 1 | −1 |

现在右端项非负，而检验数 $(1/2,5,1)$ 仍然非负：该表是最优表，$x_1=4$，$x_3=1$，$x_4=1$，其他变量为 0，$z^*=1$。目标函数值从 0 增加到 1，即 $\theta_D\abs{\bar b_2}=\tfrac12\cdot2$。如果改让 $x_6$ 进基（比值为 1），这次主元变换会使 $x_4$ 的检验数变为 $1-2\cdot1=-1<0$，破坏对偶可行性；这就是规则只取负元素并取最小比值的原因。用 HiGHS 求解该表对应的方程组，得到同样的最优解。
:::

(sec-03-06-adding-constraint)=
## 增加约束与重新优化

设某个 LP 已经求解，现在必须增加一个新约束 $\va^\T\vx\le\beta$（一个新的限制、一个割，或者一个为节省工作量而先前略去的约束，见{prf:ref}`ex-08-04-a6-q4-secondary`）。若最优解 $\vx^*$ 满足它，则 $\vx^*$ 仍然最优，因为增加约束只会缩小可行集。否则：

1. **添加该行。** 加入松弛变量 $x_{n+1}\ge0$，使 $\va^\T\vx+x_{n+1}=\beta$，并让 $x_{n+1}$ 成为新行的基变量。
2. **恢复典式。** 新行中含有原单纯形表的基变量。从新行中减去每个基变量 $x_j$ 所在行的 $a_j$ 倍。结果就是 $x_{n+1}$ 所在的行，其右端项为 $\beta-\va^\T\vx^*<0$。
3. **重新优化。** 检验数不变（松弛变量的费用为 0，原有各行不变），所以扩大后的基是对偶可行且原始不可行的：以新行为主元行开始，执行{prf:ref}`alg-03-06-dual-simplex`。

用矩阵语言表述，新的基矩阵及其逆为

$$
\begin{gathered}
\mB'=\begin{pmatrix}\mB&\vzero\\\va_B^\T&1\end{pmatrix},\\
\mB'^{-1}=\begin{pmatrix}\mB^{-1}&\vzero\\-\va_B^\T\mB^{-1}&1\end{pmatrix},
\end{gathered}
$$ (eq-03-06-new-row)

它与 {eq}`eq-03-06-bordered-inverse` 具有相同的加边结构；$\mB'^{-1}\mA$ 的新行为 $\va^\T-\va_B^\T\mB^{-1}\mA$，这就是第 2 步。等式约束可以当作两个不等式处理，或者在其所在行中加入一个人工变量。通常只需几次对偶单纯形主元变换，远少于重新求解。{prf:ref}`ex-03-07-new-constraint`把这一过程应用于{doc}`03-07-sensitivity-analysis`的贯穿示例。

(sec-03-06-artificial-constraint)=
## 人工约束法

对偶单纯形法需要一个对偶可行的初始基。当显然的松弛基既不原始可行也不对偶可行时（例如费用为正且带 $\ge$ 约束的极大化问题），可以用一个额外的约束构造出一个对偶可行基。

设 $Q$ 是检验数符号错误（极小化问题中 $\bar c_j<0$，极大化问题中 $\bar c_j>0$）的非基变量的集合，或者简单地取全部非基变量。加入**人工约束**（artificial constraint）

$$
\sum_{j\in Q}x_j+x_{n+1}=M,
$$ (eq-03-06-artificial-constraint)

其中 $M$ 是一个足够大的常数，使该约束不切掉任何有意义的可行点，$x_{n+1}$ 是它的松弛变量，作为基变量取值为 $M$。在这一行中，以检验数符号错误程度最大的列 $x_s$（对极小化问题，$\bar c_s=\min_{j\in Q}\bar c_j<0$）为主元列进行主元变换。主元为 1，新的检验数为

$$
\begin{gathered}
\bar c_j'=\bar c_j-\bar c_s\ \ (j\in Q),\\
\bar c_{n+1}'=-\bar c_s>0,
\end{gathered}
$$

而 $Q$ 以外的检验数不变。由于 $\bar c_s$ 最小，每个 $\bar c_j'\ge0$：该基是对偶可行的。接下来执行对偶单纯形法，与{ref}`sec-03-04-big-m`中的大 M 法一样，把 $M$ 作为符号携带，比较表达式 $aM+b$ 时先比较 $a$，再比较 $b$。

最终，要么松弛变量 $x_{n+1}$ 在基中且为正，人工约束不是紧约束，所得最优解就是原 LP 的最优解；要么人工约束是紧约束。在第二种情形下，若最优值依赖于 $M$（随 $M$ 无限增长），则原 LP 无界；若不依赖于 $M$，则 $x_{n+1}$ 的检验数为零，该值就是原 LP 的最优值，再用一次原始主元变换让 $x_{n+1}$ 进基，就得到一个不含 $M$ 的最优解。若对偶单纯形法因某行不可行而停止，则原 LP 不可行。该方法的例题见{prf:ref}`ex-08-04-a5-q5-artificial-constraint`。

(sec-03-06-primal-vs-dual)=
## 原始单纯形法与对偶单纯形法的比较

```{table} 原始单纯形法与对偶单纯形法的对照（极小化）。
:name: tab-03-06-primal-vs-dual

| | 原始单纯形法 | 对偶单纯形法 |
|---|---|---|
| 初始基 | 原始可行，$\bar\vb\ge\vzero$ | 对偶可行，$\bar c_j\ge0$ |
| 每次主元变换中保持 | 原始可行性 | 对偶可行性 |
| 目标 | 对偶可行性 | 原始可行性 |
| 先选取 | 进基变量（$\bar c_s<0$） | 出基变量（$\bar b_r<0$） |
| 比值检验 | 在进基列中对 $\bar a_{is}>0$ 进行 | 在出基行中对 $\bar a_{rj}<0$ 进行 |
| 目标函数值 $z$ | 不增加 | 不减少 |
| 失败检验 | 列 $\bar\mA_{\cdot s}\le\vzero$：无界 | 行 $\bar\mA_{r\cdot}\ge\vzero$ 且 $\bar b_r<0$：不可行 |
| 典型用途 | 从松弛基或第一阶段的基出发求解 | 在 $\vb$ 改变或增加新约束后重新优化 |
```

两种失败检验互为对偶：原问题无界意味着对偶问题不可行，反之亦然（{numref}`tab-03-05-status`）。在灵敏度分析中，选择是机械的：破坏对偶可行性的改变（费用、新变量）用原始单纯形主元变换修复；破坏原始可行性的改变（右端项、新约束）用对偶单纯形主元变换修复。

## 本章回顾

### 小结

- 一次单纯形迭代只需要 $\vpi^\T=\vc_B^\T\mB^{-1}$、检验数、一列 $\mB^{-1}\mA_{\cdot s}$ 和 $\bar\vb=\mB^{-1}\vb$（{eq}`eq-03-06-matrix-quantities`）；修正方法由原始数据计算这些量（{prf:ref}`alg-03-06-revised-simplex`）。
- 一次主元变换相当于用一个 eta 矩阵左乘逆矩阵（{prf:ref}`def-03-06-eta`，{prf:ref}`thm-03-06-update`），所以 $\mB_k^{-1}=\mE_{k-1}\cdots\mE_1$（{eq}`eq-03-06-product-form`）。
- 对稠密数据，修正方法每次迭代约需 $3m^2+m(n-m)$ 次乘法，而表格方法约需 $mn$ 次；稀疏性、部分定价和精度使它成为实际中的首选方法。
- 加边逆矩阵 $\mT$（{eq}`eq-03-06-bordered-inverse`）同时携带 $\mB^{-1}$ 和 $-\vpi^\T$，并给出单纯形表的整列及其检验数（{eq}`eq-03-06-T-column`）；在第一阶段中一直保留它，第二阶段就可以无须额外工作而直接开始（{prf:ref}`ex-03-06-two-phase-revised`）。
- 对偶单纯形法保持 $\bar c\ge\vzero$，移除一个取负值的基变量，并让对偶比值最小的变量进基（{eq}`eq-03-06-dual-ratio`，{prf:ref}`alg-03-06-dual-simplex`）；没有负元素的主元行证明了不可行性。
- 新约束通过对基加边（{eq}`eq-03-06-new-row`）加入，再用对偶单纯形法重新优化；当没有现成的对偶可行基时，人工界约束（{eq}`eq-03-06-artificial-constraint`）可以构造出一个对偶可行的起点。

### 常见错误

- 用右乘 $\mB^{-1}\mE$ 来更新 $\mB^{-1}$；eta 矩阵应当左乘，即 $\mE\mB^{-1}$。
- 忘记 eta 列的符号：非主元位置的元素为 $-\bar a_{is}/\bar a_{rs}$，主元位置的元素为 $1/\bar a_{rs}$。
- 在加边方法的比值检验中包括 $-z$（或 $\phi$）所在的行；$-z$ 是自由变量。
- 在第一阶段中以费用 1 对 $-z$ 定价；它的第一阶段费用为 0。
- 在对偶单纯形法中，在比值检验里取正的 $\bar a_{rj}$，或者取最大比值：两者都会破坏对偶可行性。
- 从一个非对偶可行的基开始对偶单纯形法；应先使用人工约束法。
- 由一个 $\bar b_r<0$ 但仍有负元素的行断定不可行；只有所有元素都非负的行才能证明不可行。

## 参考文献

- Bazaraa, M. S., Jarvis, J. J., & Sherali, H. D. (2010). *Linear programming
  and network flows* (4th ed.). Wiley. (Section 5.1: the revised simplex
  method; Sections 6.4 and 6.6: the dual simplex method and the
  artificial-constraint technique.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Section 3.3: implementations of the
  simplex method, revised simplex and operation counts; Section 4.5: the
  dual simplex method.)
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Chapters 8–9: pivoting and inverses in product form, and
  the simplex method using multipliers, that is, the revised simplex method.)
- Dantzig, G. B., & Orchard-Hays, W. (1954). The product form for the
  inverse in the simplex method. *Mathematical Tables and Other Aids to
  Computation, 8*(46), 64–67.
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 5: the matrix form of the
  simplex method and the revised simplex method; Section 8.1: the dual
  simplex method.)
- Lemke, C. E. (1954). The dual method of solving the linear programming
  problem. *Naval Research Logistics Quarterly, 1*(1), 36–47.
  <https://doi.org/10.1002/nav.3800010107>
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Section 4.4: the dual simplex method; Section 7.2: the revised
  simplex method.)
