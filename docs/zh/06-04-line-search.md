# 线搜索方法

## 概述

{doc}`06-02-optimality-conditions`和{doc}`06-03-kkt-conditions`中的最优性条件告诉我们如何*识别*极小点；但它们很少告诉我们如何*找到*极小点，因为方程 $\nabla f(\vx)=\vzero$ 或 Karush–Kuhn–Tucker（KKT）方程组很少能求出闭式解。数值方法把一个困难问题替换为一系列容易的问题：从当前点出发，选择一个方向，沿它移动一段距离，然后重复。本章首先把**迭代算法**（iterative algorithm）的思想精确化（算法映射、下降函数、终止准则和收敛速度），然后讨论每种这类方法都需要的一步，即**线搜索**（line search）：沿一个方向最小化一元函数。本章建立三种经典方法，并在同一个函数——$[3,6]$ 上的 $f(x)=\sin x$——上运行它们：**二分法**（bisection search），利用 $f'$ 的符号把区间减半；**黄金分割法**（golden-section search），只需要函数值，并且借助比值 $\tau=(\sqrt5-1)/2$ 在每一步重复利用一个函数值；**牛顿法**（Newton's method），用二次函数拟合，具有二次收敛性，但只有初始点足够接近极小点时才收敛。本章最后介绍二次插值法以及非精确的 **Armijo** 条件和 **Wolfe** 条件，实用程序在{doc}`06-05-unconstrained-algorithms`和{doc}`06-06-constrained-algorithms`的多维方法内部使用它们。全章用到{doc}`01-02-matrix-calculus`中的泰勒展开和方向导数，以及{doc}`02-02-convex-functions`中函数的凸性。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 用算法映射、下降函数和终止准则描述迭代算法，并把收敛序列分类为次线性、线性、超线性或二次收敛；
- 把沿某方向的线搜索化为一维最小化问题，并说明一元函数何时是单峰的；
- 手工执行二分法和黄金分割法，推导黄金比 $\tau$，并计算达到给定精度所需的函数值计算次数；
- 由二次模型推导牛顿法，证明其局部二次收敛性，并识别其失效方式；
- 过三点拟合抛物线，并以其顶点作为下一个试探点；
- 用 Armijo 条件和 Wolfe 条件检验步长，并执行回溯线搜索。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\vx^{(k)}$, $x_k$ | 算法在 $\R^n$ 中的第 $k$ 个迭代点；一维情形的标量迭代点 |
| $\vd^{(k)}$, $\alpha_k$ | 第 $k$ 次迭代的搜索方向和步长 |
| $\varphi(\alpha)$ | 沿直线的函数，$\varphi(\alpha)=f(\vx+\alpha\vd)$ |
| $[a_k,b_k]$ | 第 $k$ 次迭代的不确定区间（一维搜索） |
| $p_k<q_k$ | 黄金分割法的两个试探点 |
| $\tau$ | 黄金分割比 $(\sqrt5-1)/2\approx0.6180$ |
| $e_k$ | 第 $k$ 个迭代点的误差，$e_k=\norm{\vx^{(k)}-\vx^*}$（或 $\abs{x_k-x^*}$） |
| $\varepsilon$ | 终止准则的容差 |
| $c_1$, $c_2$ | Armijo 条件和曲率（Wolfe）条件中的常数，$0<c_1<c_2<1$ |
| $\gamma$ | 回溯因子，$0<\gamma<1$ |
| $F_N$ | 斐波那契数，$F_0=F_1=1$，$F_N=F_{N-1}+F_{N-2}$ |

(sec-06-04-iterative-algorithms)=
## 迭代算法

考虑非线性规划问题

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)\\
\st\ & g_i(\vx)\le0, && i=1,\dots,m,\\
& h_j(\vx)=0, && j=1,\dots,\ell,
\end{aligned}
$$ (eq-06-04-nlp)

其可行集为 $\Omega$（{eq}`eq-06-01-feasible-set`）。怎样设计一种方法来找到 {eq}`eq-06-04-nlp` 的局部极小点，或者至少找到一个 KKT 点？下面三个事实决定了所有答案的形态。

1. 可行集 $\Omega$ 包含无穷多个点，因此不能穷举搜索；方法必须从一个点*移动*到另一个点。
2. 在给定点处，方法只能计算局部信息：$f$、$g_i$、$h_j$ 的值及其一阶（有时二阶）导数。下一步移动的方向必须由这些信息构造。
3. 最优性条件（{doc}`06-02-optimality-conditions`，{doc}`06-03-kkt-conditions`）通常不能直接求解，但可以在给定点处*检验*。因此它们用作终止准则。

基于这些事实的方法产生一个序列 $\vx^{(0)},\vx^{(1)},\vx^{(2)},\dots$，其中每个点都按固定规则由前一个点算出。

:::{prf:definition} 算法映射
:label: def-06-04-algorithm-map

**迭代算法**由映射 $A$ 定义，它给其定义域中的每个点 $\vx$ 指定一个新点或一个候选点集合。

- **点到点映射**（point-to-point map）给出唯一的后继点，$\vx^{(k+1)}=A(\vx^{(k)})$，$k=0,1,2,\dots$
- **点到集映射**（point-to-set map）给出可接受后继点的集合 $A(\vx^{(k)})$，算法可以选取任意 $\vx^{(k+1)}\in A(\vx^{(k)})$。
:::

点到集映射是描述大多数方法的自然方式，因为一次迭代通常是若干较简单映射的*复合*：选择一个方向，再沿该方向选择步长，然后（可能）进行取整或投影。当线搜索不精确，或沿直线的极小点不唯一时，第二个映射返回的是可接受步长的集合，而不是单个步长。把整个迭代看作点到集映射，才使一个收敛定理能够涵盖许多变体（Luenberger & Ye, 2016, Ch. 7；Bazaraa et al., 2006, Ch. 7）。

当某种质量度量在相邻迭代点之间得到改进时，方法就取得了进展。

:::{prf:definition} 下降函数
:label: def-06-04-descent-function

设 $S$ 是被视为解的点的集合（例如 $f$ 的驻点）。若对每个 $\vx\notin S$ 和每个 $\vy\in A(\vx)$ 有 $M(\vy)<M(\vx)$，且对每个 $\vx\in S$ 和每个 $\vy\in A(\vx)$ 有 $M(\vy)\le M(\vx)$，则称连续函数 $M$ 为算法 $A$ 和解集 $S$ 的**下降函数**（descent function，又称价值函数）。具有下降函数的算法称为**下降算法**（descent algorithm）。
:::

无约束最小化中 $M$ 的常见选择有

$$
\begin{gathered}
M(\vx)=f(\vx),\\
M(\vx)=\norm{\vx-\vx^*}^2,\\
M(\vx)=\norm{\nabla f(\vx)},\\
M(\vx)=-\min_{\norm{\vd}\le1}\nabla f(\vx)^\T\vd .
\end{gathered}
$$

第二个在证明中有用，但不能用于计算，因为 $\vx^*$ 未知。最后一个是沿单位方向所能得到的最陡下降率；它非负，且恰在驻点处为零。方向必须归一化：没有界 $\norm{\vd}\le1$ 时，只要 $\nabla f(\vx)\ne\vzero$，线性函数 $\nabla f(\vx)^\T\vd$ 的最小值就是 $-\infty$。没有约束时，最小值在 $\vd=-\nabla f(\vx)/\norm{\nabla f(\vx)}$ 处取得，最后一个度量等于 $\norm{\nabla f(\vx)}$；当 $\vd$ 被限制为约束问题的可行方向时，它就成为一个真正不同的平稳性度量。对于约束问题，价值函数把目标函数与不可行性的度量结合起来，{ref}`sec-06-06-penalty`中的罚函数就是这样做的。

**收敛性。** 若从所考虑的初始点出发，算法产生的序列收敛到一个解——局部极小点或 KKT 点——则称算法是**收敛的**。大多数定理给出的保证要弱一些，即有界序列的每个*极限点*都是驻点。当这样的保证对每个初始点都成立时，称为**全局收敛性**（global convergence）；它不涉及收敛速度。需要初始点“足够接近”$\vx^*$ 的定理是**局部的**。

(sec-06-04-descent-step)=
### 方向与步长

本章及后两章的几乎所有方法都按如下方式更新迭代点：

$$
\vx^{(k+1)}=\vx^{(k)}+\alpha_k\vd^{(k)},\qquad \alpha_k>0,
$$ (eq-06-04-update)

其中 $\vd^{(k)}$ 是**搜索方向**（search direction），$\alpha_k$ 是**步长**（step length）。若

$$
\nabla f(\vx)^\T\vd<0.
$$ (eq-06-04-descent-direction)

则称方向是 $\vx$ 处的**下降方向**（{prf:ref}`def-06-02-descent-direction`）。这个条件之所以恰当，是因为一阶泰勒展开 $f(\vx+\alpha\vd)=f(\vx)+\alpha\nabla f(\vx)^\T\vd+o(\alpha)$：若 $\nabla f(\vx)^\T\vd<0$，则对小的 $\alpha>0$ 线性项压倒余项，因此对所有充分小的 $\alpha>0$ 有 $f(\vx+\alpha\vd)<f(\vx)$（完整论证见{prf:ref}`lem-06-02-descent`）。$\vd^{(k)}$ 的选择区分了{doc}`06-05-unconstrained-algorithms`中的各种方法；$\alpha_k$ 的选择就是**线搜索**，它是本章其余部分的主题。

**终止准则。** 对于用户选定的容差 $\varepsilon$，实用算法在以下任一情况成立时停止：

- 最优性条件近似成立，$\norm{\nabla f(\vx^{(k)})}\le\varepsilon$（或 KKT 条件在 $\varepsilon$ 精度内成立）；
- 迭代点或目标函数值停滞，$\norm{\vx^{(k+1)}-\vx^{(k)}}\le\varepsilon$ 或 $\abs{f(\vx^{(k+1)})-f(\vx^{(k)})}\le\varepsilon$，通常相对于 $\norm{\vx^{(k)}}$ 或 $\abs{f(\vx^{(k)})}$ 来度量；
- 包含解的区间足够短，$b_k-a_k\le\varepsilon$（一维搜索）；
- 达到最大迭代次数或最大函数值计算次数（这是一道安全网，而不是成功的标志）。

只用停滞检验是有风险的：收敛缓慢的方法可能在离解很远的地方只迈出极小的步子。

(sec-06-04-rates)=
### 收敛速度

为了比较方法，我们度量收敛序列的误差 $e_k=\norm{\vx^{(k)}-\vx^*}$ 趋于零的快慢。

:::{prf:definition} 收敛速度
:label: def-06-04-convergence-rates

设 $\vx^{(k)}\to\vx^*$，$e_k=\norm{\vx^{(k)}-\vx^*}>0$。称收敛是

- 速率为 $\beta\in(0,1)$ 的**线性收敛**，若 $\limsup_{k\to\infty}e_{k+1}/e_k=\beta$；
- **超线性收敛**，若 $\lim_{k\to\infty}e_{k+1}/e_k=0$；
- $r>1$ **阶收敛**，若对某个常数 $C$ 和所有充分大的 $k$ 有 $e_{k+1}\le Ce_k^{\,r}$；$2$ 阶收敛即**二次收敛**；
- **次线性收敛**，若 $e_{k+1}/e_k\to1$。
:::

序列 $e_k=2^{-k}$、$e_k=k^{-k}$、$e_k=2^{-2^k}$ 和 $e_k=1/k$ 分别是线性（速率 $\tfrac12$）、超线性、二次和次线性收敛的。线性收敛的方法每次迭代获得固定位数的有效数字（约 $-\log_{10}\beta$ 位）；二次收敛的方法一旦接近解，每次迭代有效数字位数大致翻倍；当 $e_k\approx C/k$ 时，次线性收敛的方法每多得一位有效数字，大约需要十倍的迭代次数。下面的二分法（线性，$\beta=\tfrac12$）、黄金分割法（线性，$\beta=\tau$）和牛顿法（二次）是最初的例子；最速下降法（线性）和 $\R^n$ 中的牛顿法（二次）见{doc}`06-05-unconstrained-algorithms`，Frank–Wolfe 方法（次线性）见{doc}`06-06-constrained-algorithms`。

(sec-06-04-one-dimensional)=
## 一维问题

一旦固定了 $\vx=\vx^{(k)}$ 和 $\vd=\vd^{(k)}$，步长就通过在步长区间 $0\le\alpha\le\bar\alpha$（无约束问题中 $\bar\alpha=\infty$）上最小化一元函数

$$
\begin{aligned}
\varphi(\alpha)&=f(\vx+\alpha\vd),\\
\varphi'(\alpha)&=\nabla f(\vx+\alpha\vd)^\T\vd,
\end{aligned}
$$ (eq-06-04-phi)

来求得。导数由链式法则（{prf:ref}`thm-01-02-chain-rule`）给出；特别地，对下降方向有 $\varphi'(0)=\nabla f(\vx)^\T\vd<0$。精确地最小化 $\varphi$ 称为**精确线搜索**（exact line search）；接受任何使 $\varphi$ 减小得“足够多”的步长称为**非精确线搜索**（inexact line search，{ref}`sec-06-04-inexact`）。

当移动朝向目标点 $\vy$，即 $\vd=\vy-\vx$ 时，直线上的点是凸组合，

$$
\vx+\alpha(\vy-\vx)=(1-\alpha)\vx+\alpha\vy,
$$

而 $\alpha\in[0,1]$ 恰好覆盖从 $\vx$（$\alpha=0$）到 $\vy$（$\alpha=1$）的线段。Frank–Wolfe 方法（{ref}`sec-06-06-frank-wolfe`）就是这样工作的。

本章其余部分把线搜索写成一般问题

$$
\min_{a\le x\le b} f(x)
$$ (eq-06-04-1d)

其中 $f$ 是区间 $[a,b]$ 上的一元实函数，$[a,b]$ 称为**不确定区间**（interval of uncertainty）。所有方法要么缩小这个区间，要么把单个点向极小点 $x^*$ 移动。它们需要某种结构：任意函数在 $[a,b]$ 中可能有许多局部极小点，有限次函数值计算无法确定全局极小点。所用的结构是单峰性。

:::{prf:definition} 单峰（严格拟凸）函数
:label: def-06-04-unimodal

若对 $[a,b]$ 中任意两点 $x_1\ne x_2$ 和任意 $\theta\in(0,1)$ 有

$$
f(\theta x_1+(1-\theta)x_2)<\max\{f(x_1),f(x_2)\}.
$$

则称 $f$ 在 $[a,b]$ 上**严格拟凸**（strictly quasiconvex）。对于在 $[a,b]$ 上取得最小值的函数（例如连续函数），这与**单峰性**（unimodality）相同：$f$ 在 $[a,b]$ 上有唯一的极小点 $x^*$，在 $[a,x^*]$ 上严格递减，在 $[x^*,b]$ 上严格递增。
:::

严格拟凸性加强了{prf:ref}`def-02-02-quasiconvex`中的拟凸性。严格凸函数是单峰的，但单峰函数不一定是凸的或可微的：$[-1,1]$ 上的 $f(x)=\sqrt{\abs{x}}$ 是单峰的，但既不凸也不可微；$\max\{x^2,(x-1)^2\}$ 是单峰的，但在其极小点 $x=\tfrac12$ 处不可微。下面的结果是所有区间消去方法的基础：两个函数值就足以舍去区间的一部分。

:::{prf:theorem} 区间消去
:label: thm-06-04-interval-elimination

设 $f$ 在 $[a,b]$ 上严格拟凸，$a\le p<q\le b$。

1. 若 $f(p)>f(q)$，则对每个 $x\in[a,p]$ 有 $f(x)>f(q)$，因此极小点位于 $(p,b]$ 中。
2. 若 $f(p)<f(q)$，则对每个 $x\in[q,b]$ 有 $f(x)>f(p)$，因此极小点位于 $[a,q)$ 中。
3. 若 $f(p)=f(q)$，则极小点位于 $[p,q]$ 中。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

(1) 设 $x\in[a,p)$。则 $p$ 严格位于 $x$ 与 $q$ 之间，记 $p=\theta x+(1-\theta)q$，$\theta\in(0,1)$。由严格拟凸性得 $f(p)<\max\{f(x),f(q)\}$。由于 $f(p)>f(q)$，最大者必为 $f(x)$，因此 $f(x)>f(p)>f(q)$；在 $x=p$ 处结论就是假设本身。因此 $[a,p]$ 中没有极小点。(2) 是镜像情形：对 $x\in(q,b]$，点 $q$ 位于 $p$ 与 $x$ 之间，因此 $f(q)<\max\{f(p),f(x)\}$，而 $f(q)>f(p)$ 迫使 $f(x)>f(q)>f(p)$。(3) 若 $f(p)=f(q)$，对 $x\in[a,p)$ 应用 (1) 的论证得 $f(p)<\max\{f(x),f(q)\}=\max\{f(x),f(p)\}$，因此 $f(x)>f(p)$；类似地，对 $x\in(q,b]$ 有 $f(x)>f(q)$。
:::

用文字表述：舍去函数值**较大**的试探点之外的那部分区间。在情形 2 中，证明给出了更强的界：对 $x\in(q,b]$ 有 $f(x)>f(q)$，这是情形 1 的镜像。

(sec-06-04-bisection)=
## 二分法

二分法使用导数。若 $f$ 在 $[a,b]$ 上可微且极小点位于内部，则 $f'(x^*)=0$，寻找 $x^*$ 就变成寻找 $f'$ 的根。假设

$$
f'(a)<0<f'(b),
$$ (eq-06-04-bisection-bracket)

且 $f'$ 在 $[a,b]$ 中只变号一次（由负变正）；例如当 $f$ 是连续可微的凸函数时，这一点成立。于是 $f'$ 在中点处的符号指明哪一半包含 $x^*$：若 $f'(x_k)>0$，则函数在 $x_k$ 处递增，因此 $x^*<x_k$；若 $f'(x_k)<0$，则 $x^*>x_k$。

:::{prf:algorithm} 二分法
:label: alg-06-04-bisection

**输入：** $[a,b]$ 上满足 $f'(a)<0<f'(b)$ 的可微函数 $f$；容差 $\varepsilon>0$。
**输出：** 满足 $\abs{\hat x-x^*}\le\varepsilon/2$ 的估计值 $\hat x$。

1. **初始化。** 令 $a_0=a$，$b_0=b$，$k=0$。
2. **终止检验。** 若 $b_k-a_k\le\varepsilon$，则停止并返回 $\hat x=(a_k+b_k)/2$。
3. **中点。** 计算 $x_k=(a_k+b_k)/2$ 和 $f'(x_k)$。
4. **减半。** 若 $f'(x_k)>0$，令 $a_{k+1}=a_k$，$b_{k+1}=x_k$。若 $f'(x_k)<0$，令 $a_{k+1}=x_k$，$b_{k+1}=b_k$。若 $f'(x_k)=0$，则停止，$\hat x=x_k$。
5. **下一步。** 令 $k\leftarrow k+1$，转第 2 步。
:::

每次迭代把区间减半，$b_{k+1}-a_{k+1}=\tfrac12(b_k-a_k)$，因此 $b_k-a_k=(b-a)/2^k$。一旦 $(b-a)/2^k\le\varepsilon$，终止检验即成立，即经过

$$
k=\left\lceil\log_2\frac{b-a}{\varepsilon}\right\rceil
$$ (eq-06-04-bisection-count)

次迭代，每次迭代计算一次 $f'$。由于始终有 $x^*\in[a_k,b_k]$，返回的中点与 $x^*$ 的距离不超过 $\varepsilon/2$。误差界每步减半：二分法以速率 $\tfrac12$ 线性收敛，并且对满足假设的每个函数都如此，对初始点没有要求。

:::{prf:example} 用二分法求 $\sin x$ 的极小点
:label: ex-06-04-bisection-sine

取 $\varepsilon=0.05$，用二分法在 $[3,6]$ 上最小化 $f(x)=\sin x$。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

此处 $f'(x)=\cos x$，$f'(3)=-0.9900<0$，$f'(6)=0.9602>0$。在 $[3,6]$ 上，余弦函数在 $3\pi/2$ 左侧为负、右侧为正，因此 $f'$ 在该区间内只有一个根 $x^*=3\pi/2=4.712389$，$f$ 在该区间上是单峰的。公式 {eq}`eq-06-04-bisection-count` 给出 $k=\lceil\log_2(3/0.05)\rceil=\lceil5.91\rceil=6$ 次减半。{numref}`tab-06-04-bisection`列出了各次迭代。

```{list-table} 对 $f'(x)=\cos x$ 用二分法，$[a_0,b_0]=[3,6]$，$\varepsilon=0.05$
:name: tab-06-04-bisection
:header-rows: 1

* - $k$
  - $a_k$
  - $b_k$
  - $b_k-a_k$
  - $x_k$
  - $f'(x_k)$
  - 保留的一半
* - 0
  - 3
  - 6
  - 3
  - 4.5
  - $-0.2108$
  - 右
* - 1
  - 4.5
  - 6
  - 1.5
  - 5.25
  - $0.5121$
  - 左
* - 2
  - 4.5
  - 5.25
  - 0.75
  - 4.875
  - $0.1619$
  - 左
* - 3
  - 4.5
  - 4.875
  - 0.375
  - 4.6875
  - $-0.02489$
  - 右
* - 4
  - 4.6875
  - 4.875
  - 0.1875
  - 4.78125
  - $0.06881$
  - 左
* - 5
  - 4.6875
  - 4.78125
  - 0.09375
  - 4.734375
  - $0.02198$
  - 左
* - 6
  - 4.6875
  - 4.734375
  - 0.046875
  - 停止
  -
  -
```

在 $k=6$ 时宽度 $0.046875\le0.05$ 触发终止检验，估计值为中点 $\hat x=(4.6875+4.734375)/2=4.7109375$，$f'(\hat x)=-0.001451$。其误差为 $\abs{\hat x-x^*}=0.001451$，远在保证的 $\varepsilon/2=0.025$ 之内；共计算了六次导数。*（本表的某些印刷版本把最后一个区间的右端点写成 $4.73475$，这是 $4.734375$ 的笔误，并由此算出最后的中点和导数为 $4.71113$ 和 $-0.0013$；正确的值是 $4.710938$ 和 $-0.001451$。）*
:::

```{figure} figures/06-04-interval-reduction.svg
:name: fig-06-04-interval-reduction
:width: 640px
:class: dark-light
:alt: 两幅子图，显示逐层嵌套的水平区间向 sin x 的极小点 3 pi / 2 收缩；左图为二分法及其中点，右图为黄金分割法，每行两个试探点

在 $[3,6]$ 上求 $\min\sin x$ 时的不确定区间。(a) 二分法（{prf:ref}`ex-06-04-bisection-sine`）：每一行在中点（圆点）处把区间减半。(b) 黄金分割法（{prf:ref}`ex-06-04-golden-sine`）：每一行保留区间的 $\tau\approx0.618$ 倍；实心点是从上一行重复利用的试探点，空心点是唯一的一次新计算。虚线为 $x^*=3\pi/2$。
```

(sec-06-04-golden-section)=
## 黄金分割法

二分法需要导数。**黄金分割法**只需要函数值，并且只假设单峰性（{prf:ref}`def-06-04-unimodal`），因此它也适用于不可微函数。它在 $[a_k,b_k]$ 中放置两个试探点 $p_k<q_k$，并按{prf:ref}`thm-06-04-interval-elimination`舍去一端。试探点的位置按两个目标来选取。

- **对称性。** 两个可能的新区间 $[a_k,q_k]$ 和 $[p_k,b_k]$ 应当等长，使最坏情形与最好情形一样好。这要求 $p_k-a_k=b_k-q_k$。把保留的区间所占的公共比例记为 $\tau$，

  $$
  \begin{aligned}
  p_k&=a_k+(1-\tau)(b_k-a_k),\\
  q_k&=a_k+\tau(b_k-a_k),
  \end{aligned}
  $$ (eq-06-04-golden-points)

  新区间的长度为 $b_{k+1}-a_{k+1}=\tau(b_k-a_k)$。
- **重复利用。** 保留下来的试探点应当是下一个区间的两个试探点之一，使每次迭代只需要计算**一个**新的函数值。

假设保留 $[a_k,q_k]$（定理的情形 2）。它的两个试探点位于从 $a_k$ 量起、其长度 $\tau(b_k-a_k)$ 的 $1-\tau$ 和 $\tau$ 比例处。旧点 $p_k$ 与 $a_k$ 的距离为 $(1-\tau)(b_k-a_k)$，它与新的*右*试探点 $q_{k+1}$ 重合的条件是

$$
\begin{aligned}
(1-\tau)(b_k-a_k)&=\tau\cdot\tau(b_k-a_k),\\
\text{即}\quad\tau^2&=1-\tau .
\end{aligned}
$$ (eq-06-04-golden-equation)

由对称性，当保留 $[p_k,b_k]$ 时，同一条件使 $q_k$ 成为新的左试探点 $p_{k+1}$。$\tau^2+\tau-1=0$ 的正根为

$$
\begin{aligned}
\tau&=\frac{\sqrt5-1}{2}\approx0.618034,\\
1-\tau&=\tau^2\approx0.381966,
\end{aligned}
$$ (eq-06-04-tau)

它是黄金比 $(1+\sqrt5)/2$ 的倒数；方法因此得名。{numref}`fig-06-04-golden-section-geometry`显示了相邻的两层。

```{figure} figures/06-04-golden-section-geometry.svg
:name: fig-06-04-golden-section-geometry
:width: 560px
:class: dark-light
:alt: 上下两条数轴；上面是带试探点 p 和 q 的区间 a 到 b，下面是区间 a 到 q，其中旧点 p 成为新的右试探点

黄金分割试探点的位置。上方数轴：$[a_k,b_k]$，$p_k$ 和 $q_k$ 位于其长度的 $1-\tau=\tau^2$ 和 $\tau$ 比例处。若 $f(p_k)<f(q_k)$，则舍去 $(q_k,b_k]$ 部分。下方数轴：在长度为 $\tau(b_k-a_k)$ 的新区间 $[a_k,q_k]$ 中，旧点 $p_k$ 恰好位于 $\tau$ 比例处，因此它成为 $q_{k+1}$，只需计算 $p_{k+1}$ 处的函数值。
```

:::{prf:algorithm} 黄金分割法
:label: alg-06-04-golden-section

**输入：** $[a,b]$ 上的严格拟凸函数 $f$；容差 $\varepsilon>0$。
**输出：** 满足 $\abs{\hat x-x^*}\le\varepsilon/2$ 的估计值 $\hat x$。

1. **初始化。** 令 $a_0=a$，$b_0=b$，$p_0=a+(1-\tau)(b-a)$，$q_0=a+\tau(b-a)$；计算 $f(p_0)$ 和 $f(q_0)$；令 $k=0$。
2. **终止检验。** 若 $b_k-a_k\le\varepsilon$，则停止并返回 $\hat x=(a_k+b_k)/2$。
3. **消去并重复利用。**
   - 若 $f(p_k)>f(q_k)$：令 $a_{k+1}=p_k$，$b_{k+1}=b_k$，$p_{k+1}=q_k$（函数值已知），并计算新点 $q_{k+1}=a_{k+1}+\tau(b_{k+1}-a_{k+1})$ 处的 $f$ 值。
   - 若 $f(p_k)\le f(q_k)$：令 $a_{k+1}=a_k$，$b_{k+1}=q_k$，$q_{k+1}=p_k$（函数值已知），并计算新点 $p_{k+1}=a_{k+1}+(1-\tau)(b_{k+1}-a_{k+1})$ 处的 $f$ 值。
4. **下一步。** 令 $k\leftarrow k+1$，转第 2 步。
:::

相等的情形 $f(p_k)=f(q_k)$ 归入第二种情况；由{prf:ref}`thm-06-04-interval-elimination`的情形 3，极小点位于 $[p_k,q_k]\subset[a_k,q_k]$ 中，因此不会丢失任何东西。

**需要计算多少次函数值？** 经过 $k$ 次缩减后，区间长度为 $\tau^k(b-a)$，而前 $k$ 次缩减需要 $2+(k-1)=k+1$ 个函数值。因此，用 $N$ 次函数值计算，最终区间长度为 $\tau^{N-1}(b-a)$，而要达到精度 $b_k-a_k\le\varepsilon$ 需要

$$
\begin{aligned}
k&=\left\lceil\frac{\ln(\varepsilon/(b-a))}{\ln\tau}\right\rceil\ \text{次缩减},\\
N&=k+1\ \text{次函数值计算。}
\end{aligned}
$$ (eq-06-04-golden-count)

该方法以速率 $\tau\approx0.618$ 线性收敛，每步比二分法（$\tfrac12$）慢：把区间缩小 $10^{6}$ 倍需要 $30$ 个函数值，而二分法需要 $20$ 个导数值。作为交换，黄金分割法不使用导数，并且能容忍折点。

:::{prf:example} 用黄金分割法求 $\sin x$ 的极小点
:label: ex-06-04-golden-sine

取 $\varepsilon=0.05$，用黄金分割法在 $[3,6]$ 上最小化 $f(x)=\sin x$。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

公式 {eq}`eq-06-04-golden-count` 给出 $k=\lceil\ln(0.05/3)/\ln0.618\rceil=\lceil8.51\rceil=9$ 次缩减和 $N=10$ 个函数值。初始试探点为 $p_0=3+0.381966\cdot3=4.1459$ 和 $q_0=3+0.618034\cdot3=4.8541$。由于在极小点附近函数值只在小数点后第四位才有差别，{numref}`tab-06-04-golden`中的函数值保留六位小数。

```{list-table} 对 $\sin x$ 用黄金分割法，$[a_0,b_0]=[3,6]$，$\varepsilon=0.05$（每行中重复利用的点以粗体表示）
:name: tab-06-04-golden
:header-rows: 1

* - $k$
  - $a_k$
  - $b_k$
  - $b_k-a_k$
  - $p_k$
  - $q_k$
  - $f(p_k)$
  - $f(q_k)$
* - 0
  - 3
  - 6
  - 3
  - 4.1459
  - 4.8541
  - $-0.843789$
  - $-0.989976$
* - 1
  - 4.1459
  - 6
  - 1.8541
  - **4.8541**
  - 5.2918
  - $-0.989976$
  - $-0.836787$
* - 2
  - 4.1459
  - 5.2918
  - 1.1459
  - 4.5836
  - **4.8541**
  - $-0.991717$
  - $-0.989976$
* - 3
  - 4.1459
  - 4.8541
  - 0.7082
  - 4.4164
  - **4.5836**
  - $-0.956516$
  - $-0.991717$
* - 4
  - 4.4164
  - 4.8541
  - 0.4377
  - **4.5836**
  - 4.6869
  - $-0.991717$
  - $-0.999676$
* - 5
  - 4.5836
  - 4.8541
  - 0.2705
  - **4.6869**
  - 4.7508
  - $-0.999676$
  - $-0.999263$
* - 6
  - 4.5836
  - 4.7508
  - 0.1672
  - 4.6475
  - **4.6869**
  - $-0.997892$
  - $-0.999676$
* - 7
  - 4.6475
  - 4.7508
  - 0.1033
  - **4.6869**
  - 4.7113
  - $-0.999676$
  - $-0.999999$
* - 8
  - 4.6869
  - 4.7508
  - 0.06386
  - **4.7113**
  - 4.7264
  - $-0.999999$
  - $-0.999902$
* - 9
  - 4.6869
  - 4.7264
  - 0.03947
  - 停止
  -
  -
  -
```

第 0 行：$f(p_0)>f(q_0)$，因此舍去 $[3,p_0]$，新区间为 $[4.1459,6]$，旧的 $q_0$ 成为 $p_1$，只有 $q_1=4.1459+0.618034\cdot1.8541=5.2918$ 是新点。第 1 行：$f(p_1)<f(q_1)$，因此舍去 $(q_1,6]$，$p_1$ 成为 $q_2$。以后各行都遵循同样的模式，每个宽度都是前一个的 $\tau$ 倍（$3\tau^8=0.06386$）。搜索在 $k=9$ 时停止，$[a_9,b_9]=[4.6869,4.7264]$，它包含 $x^*=4.7124$；估计值为 $\hat x=4.7067$，误差 $0.0057$，共用了 $10$ 个函数值。二分法用六个导数值得到了更小的最终区间；黄金分割法则完全没有使用导数。
:::

**斐波那契法。** 若函数值的个数 $N$ 事先固定，则可以在全部 $N$ 步上整体优化带重复利用的对称放置方式，而不是每步都相同。Kiefer（1953）找到的最优策略把试探点放在相邻斐波那契数之比 $F_{N-1}/F_N$、$F_{N-2}/F_{N-1}$、……处（$F_0=F_1=1$，$F_N=F_{N-1}+F_{N-2}$），并把区间缩小为原来的 $1/F_N$，这是用 $N$ 次计算所能达到的最好结果（忽略最后一步中一个小的可分辨性常数）。由于 $F_{N-1}/F_N\to\tau$，黄金分割法是斐波那契法的极限，而且损失很小：计算 $N=10$ 次后，缩小因子分别为 $1/F_{10}=1/89=0.01124$ 和 $\tau^9=0.01316$；当 $N$ 很大时，黄金分割法的区间只长约 $17\,\%$。实践中更常用黄金分割法，因为它不需要事先知道 $N$。

(sec-06-04-newton)=
## 牛顿法

二分法和黄金分割法每计算一次就把区间缩小一个固定的倍数。若 $f$ 二阶可微，则通过**曲线拟合法**（curve fitting）可以得到快得多的方法：在当前点 $x_k$ 附近，用在函数值、斜率和曲率上都与 $f$ 一致的二次函数

$$
q(x)=f(x_k)+f'(x_k)(x-x_k)+\tfrac12f''(x_k)(x-x_k)^2,
$$ (eq-06-04-newton-model)

代替 $f$，并取其驻点作为下一个迭代点。令 $q'(x)=f'(x_k)+f''(x_k)(x-x_k)=0$，得

$$
x_{k+1}=x_k-\frac{f'(x_k)}{f''(x_k)}.
$$ (eq-06-04-newton-1d)

若 $f''(x_k)>0$，则模型是开口向上的抛物线，$x_{k+1}$ 是它的极小点。方程 {eq}`eq-06-04-newton-1d` 也是求解 $f'(x)=0$ 的牛顿法：$x_{k+1}$ 是 $f'$ 的图像在 $x_k$ 处的切线与横轴的交点（{numref}`fig-06-04-newton-tangents`）。

:::{prf:algorithm} 一维牛顿法
:label: alg-06-04-newton-1d

**输入：** 二阶可微函数 $f$；初始点 $x_0$；容差 $\varepsilon>0$。
**输出：** $f$ 的近似驻点。

1. **初始化。** 令 $k=0$。
2. **终止检验。** 若 $\abs{f'(x_k)}\le\varepsilon$，则停止并返回 $x_k$。
3. **牛顿步。** 令 $x_{k+1}=x_k-f'(x_k)/f''(x_k)$。
4. **下一步。** 令 $k\leftarrow k+1$，转第 2 步。
:::

:::{prf:theorem} 牛顿法的局部二次收敛性
:label: thm-06-04-newton-quadratic

设 $f$ 在满足 $f'(x^*)=0$ 和 $f''(x^*)\ne0$ 的点 $x^*$ 附近三阶连续可微。则存在 $\delta>0$，使得对每个满足 $\abs{x_0-x^*}<\delta$ 的 $x_0$，迭代 {eq}`eq-06-04-newton-1d` 有定义，收敛到 $x^*$，并满足

$$
\begin{aligned}
&\abs{x_{k+1}-x^*}\le C\abs{x_k-x^*}^2,\\
&\lim_{k\to\infty}\frac{x_{k+1}-x^*}{(x_k-x^*)^2}=\frac{f'''(x^*)}{2f''(x^*)}.
\end{aligned}
$$ (eq-06-04-newton-error)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

选取 $\delta_0>0$，使得在 $I=[x^*-\delta_0,x^*+\delta_0]$ 上 $\abs{f''(x)}\ge m>0$ 且 $\abs{f'''(x)}\le L$（由连续性和 $f''(x^*)\ne0$，这是可以做到的）。对 $x_k\in I$，在 $x_k$ 处应用带余项的泰勒定理，对介于 $x_k$ 与 $x^*$ 之间的某个 $\xi_k$ 有

$$
\begin{aligned}
0=f'(x^*)&=f'(x_k)+f''(x_k)(x^*-x_k)\\
&\quad+\tfrac12f'''(\xi_k)(x^*-x_k)^2.
\end{aligned}
$$

除以 $f''(x_k)$ 并利用 {eq}`eq-06-04-newton-1d`，即 $x_k-f'(x_k)/f''(x_k)=x_{k+1}$：

$$
x_{k+1}-x^*=\frac{f'''(\xi_k)}{2f''(x_k)}(x_k-x^*)^2 .
$$

因此 $\abs{x_{k+1}-x^*}\le C\abs{x_k-x^*}^2$，其中 $C=L/(2m)$。若再取 $\delta=\min\{\delta_0,1/(2C)\}$ 并且 $\abs{x_k-x^*}<\delta$，则 $\abs{x_{k+1}-x^*}\le C\delta\abs{x_k-x^*}\le\tfrac12\abs{x_k-x^*}$，因此 $x_{k+1}$ 仍在 $I$ 中，由归纳法得 $x_k\to x^*$。由于 $\xi_k\to x^*$ 且 $x_k\to x^*$，该比值收敛到 $f'''(x^*)/(2f''(x^*))$。
:::

该定理是局部的，并且没有说明找到的是*哪一个*驻点：牛顿法寻找的是 $f'$ 的零点，极大点（$f''(x^*)<0$）和极小点同样会吸引迭代点。下面是典型的失效情形。

- **负曲率。** 若 $f''(x_k)<0$，模型 {eq}`eq-06-04-newton-model` 是开口向下的抛物线，“牛顿步”走向它的极大点，往往是上坡方向。
- **曲率很小。** 若 $f''(x_k)\approx0$，步长 $-f'(x_k)/f''(x_k)$ 极大，迭代点被抛到很远的地方；若 $f''(x_k)=0$，则步长无定义。
- **离开区间。** 方程 {eq}`eq-06-04-newton-1d` 不考虑界 $a\le x\le b$。
- **循环或发散**：从不好的初始点出发，即使对凸函数也可能出现；例如，用牛顿法求 $f'(x)=\arctan x$（某个严格凸函数 $f$ 的导数）的根，从 $\abs{x_0}>1.39$ 出发会发散。

因此，实用的一维程序会对牛顿步加以**保护措施**（safeguard）：它们维持一个满足 $f'(a_k)<0<f'(b_k)$ 的区间 $[a_k,b_k]$，每当牛顿点落在该区间之外或曲率不为正时，就退回二分步。这样在极小点附近牛顿步会被接受，快速的局部收敛性得以保留。

:::{prf:example} 从两个初始点出发用牛顿法求 $\sin x$ 的极小点
:label: ex-06-04-newton-sine

取 $\varepsilon=10^{-6}$，用牛顿法在 $[3,6]$ 上最小化 $f(x)=\sin x$，(a) 从 $x_0=4$ 出发，(b) 从 $x_0=3$ 出发。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

导数为 $f'(x)=\cos x$ 和 $f''(x)=-\sin x$，因此 {eq}`eq-06-04-newton-1d` 成为 $x_{k+1}=x_k+\cos x_k/\sin x_k$。

```{list-table} 对 $\sin x$ 用牛顿法（数值保留 4 位有效数字）
:name: tab-06-04-newton
:header-rows: 1

* - 初始点
  - $k$
  - $x_k$
  - $f'(x_k)$
  - $f''(x_k)$
  - $x_{k+1}$
  - $\abs{x_k-3\pi/2}$
* - (a) $x_0=4$
  - 0
  - 4.000
  - $-0.6536$
  - $0.7568$
  - 4.864
  - $0.7124$
* -
  - 1
  - 4.864
  - $0.1507$
  - $0.9886$
  - 4.711
  - $0.1513$
* -
  - 2
  - 4.711
  - $-1.165\times10^{-3}$
  - $1.000$
  - 4.712
  - $1.165\times10^{-3}$
* -
  - 3
  - 4.712
  - $5.274\times10^{-10}$
  - $1.000$
  - 停止
  - $5.274\times10^{-10}$
* - (b) $x_0=3$
  - 0
  - 3.000
  - $-0.9900$
  - $-0.1411$
  - $-4.015$
  - $1.712$
* -
  - 1
  - $-4.015$
  - $-0.6420$
  - $-0.7667$
  - $-4.853$
  -
* -
  - 2
  - $-4.853$
  - $0.1398$
  - $-0.9902$
  - $-4.711$
  -
* -
  - 3
  - $-4.711$
  - $-9.272\times10^{-4}$
  - $-1.000$
  - $-4.712$
  -
* -
  - 4
  - $-4.712$
  - $2.657\times10^{-10}$
  - $-1.000$
  - 停止
  -
```

**(a)** 从 $x_0=4$ 出发，曲率 $f''(4)=0.7568$ 为正，模型是凸的，迭代点 $4$ → $4.8637$ → $4.7112$ → $4.712389$ 在三步内达到 $x^*=3\pi/2$，精确到十位数字。误差 $0.7124$、$0.1513$、$1.165\times10^{-3}$、$5.3\times10^{-10}$ 的缩小甚至比二次收敛还快：这里 $f'''(x^*)=-\cos(3\pi/2)=0$，因此 {eq}`eq-06-04-newton-error` 中的极限为零，事实上 $x_{k+1}-x^*\approx-(x_k-x^*)^3/3$（比值 $e_{k+1}/e_k^3$ 依次为 $0.42,\ 0.34,\ 0.333$）。四对导数值所达到的精度远远超过区间方法的六次或十次计算。

**(b)** 从 $x_0=3$ 出发，曲率 $f''(3)=-0.1411$ 为负且很小。模型是一条平坦的开口向下的抛物线，其驻点 $x_1=3-(-0.9900)/(-0.1411)=-4.0153$ 远在 $[3,6]$ 之外。从那里出发，迭代再次非常快地收敛到 $-4.712389=-3\pi/2$，在该点 $f''=-1<0$，$f=\sin(-3\pi/2)=1$：这是 $\sin x$ 的**极大点**，不是极小点，甚至不在区间内。终止检验 $\abs{f'(x_k)}\le\varepsilon$ 得到满足，因此算法报告成功。
:::

{prf:ref}`ex-06-04-newton-sine`的 (b) 同时展示了三种失效方式，并说明了为什么必须检验 $f''(x_k)$ 的符号和界；从 $x_0=3$ 出发的带保护措施的方法会拒绝这一步，转而对区间 $[3,6]$ 进行二分。

```{figure} figures/06-04-newton-tangents.svg
:name: fig-06-04-newton-tangents
:width: 640px
:class: dark-light
:alt: 两幅子图，显示 cos x 的图像及切线；左图中从 x0 = 4 出发的切线很快到达根 3 pi / 2；右图中 x0 = 3 处几乎水平的切线在 -4 附近与横轴相交

牛顿法作为 $f'(x)=\cos x$ 上的切线步（{prf:ref}`ex-06-04-newton-sine`）。(a) 从 $x_0=4$ 出发，切线（虚线）与横轴交于 $x_1=4.864$ 和 $x_2=4.711$，紧邻根 $3\pi/2$（$\sin x$ 的极小点）。(b) 从 $x_0=3$ 出发，切线几乎水平（$f''(3)=-0.141$），与横轴交于 $x_1=-4.015$；随后迭代点收敛到 $-3\pi/2$，这是 $\sin x$ 在 $[3,6]$ 之外的一个极大点。
```

(sec-06-04-curve-fitting)=
## 二次插值法

牛顿法利用一点处的导数拟合二次函数。当导数不可用或计算代价高时，可以过三个点 $x_1<x_2<x_3$ 处的*函数值*拟合二次函数，其中 $f(x_2)<\min\{f(x_1),f(x_3)\}$（即极小点的一个括区间）。过这三点的抛物线的顶点为

$$
x_{\text{new}}=x_2-\frac12\,
\frac{(x_2-x_1)^2\bigl(f(x_2)-f(x_3)\bigr)-(x_2-x_3)^2\bigl(f(x_2)-f(x_1)\bigr)}
{(x_2-x_1)\bigl(f(x_2)-f(x_3)\bigr)-(x_2-x_3)\bigl(f(x_2)-f(x_1)\bigr)} ,
$$ (eq-06-04-quadratic-vertex)

由于括区间条件，它位于 $(x_1,x_3)$ 中。新点替换旧点中的一个，使最好的点及其两个相邻点再次构成括区间，然后重复这一过程。对于 $\sin x$，取初始点 $4$、$4.5$、$5$，各次顶点为 $4.711130$、$4.712581$、$4.712393$ 和 $4.712389$，误差为 $1.3\times10^{-3}$、$1.9\times10^{-4}$、$3.7\times10^{-6}$ 和 $1.1\times10^{-11}$：每次迭代一个函数值，且超线性收敛（对光滑函数，收敛阶约为 $1.3$；Luenberger & Ye, 2016, Ch. 8）。与牛顿法一样，二次插值法在远离极小点时可能表现不佳，因此稳健的程序把它与黄金分割步结合起来（Brent 方法；Brent, 1973）。

(sec-06-04-inexact)=
## 非精确线搜索：Armijo 条件与 Wolfe 条件

在 $\R^n$ 中的方法内部，精确线搜索很少值得付出其代价：方向 $\vd^{(k)}$ 在下一次迭代中反正会改变，重要的是一个能取得合理进展的步长。非精确线搜索接受通过两项检验的任何 $\alpha$。取 {eq}`eq-06-04-phi` 中的 $\varphi$ 和常数 $0<c_1<c_2<1$：

$$
\varphi(\alpha)\le\varphi(0)+c_1\alpha\varphi'(0)
\qquad\text{（充分下降条件，即 Armijo 条件）},
$$ (eq-06-04-armijo)

$$
\varphi'(\alpha)\ge c_2\varphi'(0)
\qquad\text{（曲率条件）}.
$$ (eq-06-04-curvature)

两者合起来就是 **Wolfe 条件**（Wolfe conditions；Armijo, 1966；Wolfe, 1969）。条件 {eq}`eq-06-04-armijo` 说实际下降量至少是切线所预测下降量的 $c_1$ 倍；它排除过长的步长。由于 $\varphi'(0)<0$，右端是一条斜率为 $c_1\varphi'(0)$ 的直线，比切线平缓，因此所有充分小的 $\alpha>0$ 都满足它。条件 {eq}`eq-06-04-curvature` 排除过短的步长：它要求斜率已增大到至少 $c_2\varphi'(0)$，从而步长已越过陡峭的初始段。**强 Wolfe 条件**（strong Wolfe conditions）把 {eq}`eq-06-04-curvature` 替换为 $\abs{\varphi'(\alpha)}\le c_2\abs{\varphi'(0)}$，它还排除了 $\varphi$ 已经陡峭上升的步长。典型取值为 $c_1=10^{-4}$ 和 $c_2=0.9$（牛顿法和拟牛顿法）或 $c_2=0.1$（共轭梯度法）。

若 $\varphi$ 在 $\alpha>0$ 上连续可微且有下界，则存在满足 Wolfe 条件和强 Wolfe 条件的步长区间（Nocedal & Wright, 2006, Lemma 3.1）。仅满足 Armijo 条件的一种简单方法是回溯。

:::{prf:algorithm} 回溯（Armijo）线搜索
:label: alg-06-04-backtracking

**输入：** $\vx$，下降方向 $\vd$（$\nabla f(\vx)^\T\vd<0$），初始步长 $\bar\alpha>0$，常数 $c_1\in(0,1)$ 和 $\gamma\in(0,1)$（例如 $\gamma=\tfrac12$）。
**输出：** 满足 {eq}`eq-06-04-armijo` 的步长 $\alpha$。

1. **初始化。** 令 $\alpha=\bar\alpha$。
2. **检验。** 若 $f(\vx+\alpha\vd)\le f(\vx)+c_1\alpha\nabla f(\vx)^\T\vd$，则停止并返回 $\alpha$。
3. **缩小。** 令 $\alpha\leftarrow\gamma\alpha$，转第 2 步。
:::

由于 Armijo 条件对所有充分小的 $\alpha$ 成立，循环必定终止。对于牛顿型方向，从 $\bar\alpha=1$ 开始是自然的，因为在解附近单位步长是正确的（{doc}`06-05-unconstrained-algorithms`）。除非被接受的步长就是 $\bar\alpha$ 本身，否则它是某个因过长而被拒绝的步长的 $\gamma$ 倍，因此不会不必要地短；这就是回溯可以不用曲率条件的原因。

作为例示，再取 $f(x)=\sin x$，现在把它看作 $\R$ 上的函数，在 $x=3$ 处取下降方向 $d=-f'(3)=0.9900$，于是 $\varphi(\alpha)=\sin(3+0.99\alpha)$，$\varphi(0)=0.1411$，$\varphi'(0)=-0.9801$。取 $c_1=0.25$，$\bar\alpha=4$，$\gamma=\tfrac12$：

- $\alpha=4$：$\varphi(4)=0.6263>\varphi(0)+c_1\cdot4\cdot\varphi'(0)=-0.8390$；拒绝（这一步越过极小点，进入了 $\sin x$ 的下一段上升部分）；
- $\alpha=2$：$\varphi(2)=-0.9644\le0.1411-0.25\cdot2\cdot0.9801=-0.3489$；接受，得 $x=4.9800$。

被接受的步长在 $c_2=0.5$ 时也满足强 Wolfe 条件：$\abs{\varphi'(2)}=0.2618\le0.5\cdot0.9801=0.4900$。{numref}`fig-06-04-armijo-wolfe`显示了可接受的集合：Armijo 条件在 $0<\alpha\le2.749$ 上成立，曲率条件在 $1.207\le\alpha\le5.426$ 上成立（只有越过下一个局部极大点之后，斜率才再次低于 $c_2\varphi'(0)$），因此 Wolfe 步长构成 $[1.207,2.749]$，强 Wolfe 步长构成 $[1.207,2.253]$；精确极小点 $\alpha^*=1.730$（它到达 $x^*=3\pi/2$）位于两者之中。

```{figure} figures/06-04-armijo-wolfe.svg
:name: fig-06-04-armijo-wolfe
:width: 560px
:class: dark-light
:alt: phi(alpha) = sin(3 + 0.99 alpha) 的图像，附零点处的切线、较平缓的 Armijo 直线、表示可接受步长的阴影区间，并标出两次回溯试探

$\varphi(\alpha)=\sin(3+0.99\alpha)$ 在 $c_1=0.25$，$c_2=0.5$ 时的 Armijo 条件与 Wolfe 条件。位于 Armijo 直线 $\varphi(0)+c_1\alpha\varphi'(0)$（虚线）下方时下降量是充分的。条形：满足 Armijo 条件、Wolfe 条件和强 Wolfe 条件的步长。从 $\bar\alpha=4$ 开始的回溯拒绝 $\alpha=4$（空心标记），接受 $\alpha=2$（实心标记）。
```

## 本章回顾

### 小结

- 迭代算法是一个（点到点或点到集的）映射（{prf:ref}`def-06-04-algorithm-map`），它使下降函数减小（{prf:ref}`def-06-04-descent-function`）；大多数方法沿下降方向 {eq}`eq-06-04-descent-direction` 迈步 $\vx^{(k+1)}=\vx^{(k)}+\alpha_k\vd^{(k)}$，并在近似最优时停止。
- 收敛速度分为线性、超线性、二次和次线性（{prf:ref}`def-06-04-convergence-rates`）。
- 线搜索最小化 $\varphi(\alpha)=f(\vx+\alpha\vd)$ {eq}`eq-06-04-phi`；对于单峰函数，两个函数值就足以舍去区间的一部分（{prf:ref}`thm-06-04-interval-elimination`）。
- 二分法（{prf:ref}`alg-06-04-bisection`）每步用一个导数值把区间减半，需要 $\lceil\log_2((b-a)/\varepsilon)\rceil$ 步 {eq}`eq-06-04-bisection-count`。
- 黄金分割法（{prf:ref}`alg-06-04-golden-section`）每步保留 $\tau=(\sqrt5-1)/2$ 比例的区间，$\tau$ 由 $\tau^2=1-\tau$ {eq}`eq-06-04-golden-equation` 得到，每步需要一个函数值；斐波那契法是它的最优有限步版本。
- 牛顿法 $x_{k+1}=x_k-f'(x_k)/f''(x_k)$ {eq}`eq-06-04-newton-1d` 在满足 $f''(x^*)\ne0$ 的点附近二次收敛（{prf:ref}`thm-06-04-newton-quadratic`），但可能走向极大点或离开区间（{numref}`fig-06-04-newton-tangents`）。
- 二次插值法 {eq}`eq-06-04-quadratic-vertex` 只使用函数值，且超线性收敛。
- 非精确线搜索接受满足 Armijo 条件 {eq}`eq-06-04-armijo` 和曲率条件 {eq}`eq-06-04-curvature` 的步长；回溯（{prf:ref}`alg-06-04-backtracking`）保证前者成立。

### 常见错误

- **把直线写成 $\alpha\vx+(1-\alpha)\vy$。** 取 $\vd=\vy-\vx$ 时，点 $\vx+\alpha\vd$ 等于 $(1-\alpha)\vx+\alpha\vy$；权重恰好相反，并且 $\alpha=0$ 必须给出当前点。
- **舍错了一端。** 在区间消去中，舍去的是函数值*较大*的试探点之外的部分；两个函数值相等时，极小点位于两点之间。
- **重新计算两个黄金分割点。** 取 $\tau=(\sqrt5-1)/2$ 时，总有一个旧点被重复利用；把 $\tau$ 舍入为 $0.618$ 并重新计算两个点，会使工作量加倍，并使试探点逐渐漂移。
- **轻信牛顿法的终止检验。** $\abs{f'(x_k)}\le\varepsilon$ 在极大点处同样成立；应检验 $f''(x_k)>0$ 和界（{prf:ref}`ex-06-04-newton-sine`，初始点 $x_0=3$）。
- **只因步长很小就停止。** 缓慢的方法在离解很远处也会迈出小步；应把停滞检验与最优性检验结合起来。
- **在没有变号时使用二分法。** 二分法需要 $f'(a)<0<f'(b)$；若两端符号相同，则极小点在端点处，或者括区间有误。


## 参考文献

- Armijo, L. (1966). Minimization of functions having Lipschitz continuous
  first partial derivatives. *Pacific Journal of Mathematics, 16*(1), 1–3.
  <https://doi.org/10.2140/pjm.1966.16.1>
- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 7: the
  concept of an algorithm, algorithmic maps, descent functions and
  convergence; Sections 8.1–8.3: line search without and with derivatives,
  golden-section and Fibonacci search, bisection, Newton's method, quadratic
  fit and the Armijo rule.)
- Brent, R. P. (1973). *Algorithms for minimization without derivatives*.
  Prentice-Hall. (Chapter 5: one-dimensional minimization combining
  golden-section search and parabolic interpolation.)
- Kiefer, J. (1953). Sequential minimax search for a maximum. *Proceedings of
  the American Mathematical Society, 4*(3), 502–506.
  <https://doi.org/10.1090/S0002-9939-1953-0055639-3>
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming* (4th
  ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3> (Chapter 7:
  algorithms, global convergence and speed of convergence; Chapter 8: line
  search by Fibonacci and golden-section search, Newton's method, curve
  fitting and inexact line search.)
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Chapter 2:
  overview of line-search methods; Appendix A.2: rates of convergence;
  Sections 3.1 and 3.5: the Wolfe conditions, backtracking and step-length selection by
  interpolation.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 21: nonlinear programming algorithms, unconstrained
  one-variable search.)
- Wolfe, P. (1969). Convergence conditions for ascent methods. *SIAM Review,
  11*(2), 226–235. <https://doi.org/10.1137/1011036>
