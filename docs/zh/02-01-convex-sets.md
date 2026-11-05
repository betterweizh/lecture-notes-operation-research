# 凸集

## 概述

最优化问题的可行集 $\Omega$ 是 $\R^n$ 的一个子集，它的形状在很大程度上决定了关于该问题能得出哪些结论。本章汇集全书其余部分要用到的关于集合的事实。第一部分是拓扑方面的：开集、闭集、有界集和紧集，序列及其极限点，以及魏尔斯特拉斯定理（{prf:ref}`thm-02-01-weierstrass`），它保证连续函数在非空紧集上能取到最小值。第二部分是几何方面的：仿射集与凸集（convex set）、凸包与锥、超平面与半空间、多面体、球与椭球、半正定矩阵锥、保凸运算以及极点。本章用七个示例集合（{prf:ref}`ex-02-01-seven-sets`）逐一检验所有这些性质。本章最后给出分离超平面定理和支撑超平面定理（{prf:ref}`thm-02-01-separation-point`、{prf:ref}`thm-02-01-supporting-hyperplane`），它们是对偶理论和最优性条件背后的几何事实。本章用到{doc}`01-01-linear-algebra`中的向量、范数、内积和（半）正定矩阵，以及{doc}`01-03-introduction`中的可行集和最优值。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 判断 $\R^n$ 中的集合是否为开集、闭集、有界集或紧集，并求序列的极限点；
- 应用魏尔斯特拉斯定理证明最优化问题有最优解；
- 识别仿射集、凸集、锥、超平面、半空间、多面体、球、椭球和范数锥，并证明它们是凸的；
- 构造凸包和锥包，并识别极点；
- 利用交集和仿射映射证明一个集合是凸的；
- 叙述分离超平面定理和支撑超平面定理，并解释前者的证明。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\vx,\vy,\vz$; $\vx_1,\dots,\vx_k$ | $\R^n$ 中的点（带下标的粗体字母表示点，而不是分量） |
| $\vx^{(k)}$ | 点列的第 $k$ 项 |
| $B(\vx,r)$, $\bar B(\vx,r)$ | 开球 $\{\vy:\norm{\vy-\vx}<r\}$ 和闭球 $\{\vy:\norm{\vy-\vx}\le r\}$ |
| $\intr\Omega$, $\operatorname{cl}\Omega$, $\operatorname{bd}\Omega$ | $\Omega$ 的内部、闭包和边界 |
| $\theta$, $\theta_i$ | 仿射组合、凸组合或锥组合的权重 |
| $\aff S$, $\conv S$, $\cone S$ | 集合 $S$ 的仿射包、凸包和锥包 |
| $C$, $K$ | 凸集；锥 |
| $H$, $P$ | 超平面；多面体 |
| $\R^n_+$ | 非负卦限 $\{\vx\in\R^n:\vx\ge\vzero\}$ |
| $\mathcal E$ | 椭球 |
| $\Sym^n$, $\Sym^n_+$, $\Sym^n_{++}$ | $n\times n$ 对称矩阵、半正定矩阵和正定矩阵 |
| $\operatorname{proj}_C(\vy)$ | $C$ 中距 $\vy$ 最近的点 |

除非另有说明，$\norm{\cdot}$ 表示欧氏范数 $\norm{\vx}=\sqrt{\vx^\T\vx}$。

(sec-02-01-topology)=
## 开集、闭集、有界集与紧集

最优化问题要求集合中的最好点，而最好点是否存在，取决于集合是否包含它的“边缘”以及它是否延伸到无穷远。球使这些概念变得精确。以 $\vx$ 为中心、半径为 $r>0$ 的**开球**（open ball）为 $B(\vx,r)=\{\vy\in\R^n:\norm{\vy-\vx}<r\}$，**闭球**（closed ball）为 $\bar B(\vx,r)=\{\vy:\norm{\vy-\vx}\le r\}$。

:::{prf:definition} 内部、开集与闭集
:label: def-02-01-open-closed

设 $\Omega\subseteq\R^n$。

1. 若存在某个球 $B(\vx,\varepsilon)$，$\varepsilon>0$，完全包含在 $\Omega$ 中，则称点 $\vx\in\Omega$ 为 $\Omega$ 的**内点**（interior point）。全体内点构成的集合称为**内部**（interior）$\intr\Omega$。
2. 若 $\Omega$ 的每个点都是内点，即 $\Omega=\intr\Omega$，则称 $\Omega$ 为**开集**（open set）。
3. 若 $\Omega$ 的补集 $\R^n\setminus\Omega$ 是开集，则称 $\Omega$ 为**闭集**（closed set）。
4. 若以点 $\vx$（不论是否属于 $\Omega$）为中心的每个球都既包含 $\Omega$ 中的点，又包含不在 $\Omega$ 中的点，则称 $\vx$ 为 $\Omega$ 的**边界点**（boundary point）。全体边界点构成的集合称为**边界**（boundary）$\operatorname{bd}\Omega$，而 $\operatorname{cl}\Omega=\Omega\cup \operatorname{bd}\Omega$ 称为 $\Omega$ 的**闭包**（closure）。
:::

开集不包含它的任何边界点，闭集包含它的全部边界点。大多数集合既不是开集也不是闭集：区间 $[0,1)$ 包含边界点 $0$，却不包含边界点 $1$。因此“闭”并不是“开”的反面。空集和 $\R^n$ 既是开集又是闭集。在最优化中，由有限个不等式 $g_i(\vx)\le0$（$g_i$ 连续）定义的集合是闭集，而严格不等式 $g_i(\vx)<0$ 定义的是开集。

### 序列与极限点

$\R^n$ 中的**序列**（sequence）是一个由点组成的无穷列表 $\vx^{(0)},\vx^{(1)},\vx^{(2)},\dots$，记为 $\{\vx^{(k)}\}$。迭代算法（{doc}`06-04-line-search`、{doc}`06-05-unconstrained-algorithms`）会产生这样的序列，对算法的分析就是要回答序列是否趋于稳定以及稳定在何处。

:::{prf:definition} 收敛、子序列、极限点
:label: def-02-01-limit-point

1. 若当 $k\to\infty$ 时 $\norm{\vx^{(k)}-\vx^*}\to0$，则称 $\{\vx^{(k)}\}$ **收敛**（converge）于 $\vx^*$，记为 $\vx^{(k)}\to\vx^*$。
2. 若 $k_1<k_2<k_3<\cdots$ 是递增的指标，则序列 $\vx^{(k_1)},\vx^{(k_2)},\dots$ 称为 $\{\vx^{(k)}\}$ 的**子序列**（subsequence）；例如由偶数指标的项组成的序列。
3. 若某个子序列收敛于 $\vx^*$，则称 $\vx^*$ 为 $\{\vx^{(k)}\}$ 的**极限点**（limit point）。
:::

收敛序列恰有一个极限点，即它的极限。一个序列也可能有多个极限点，或者没有极限点。

:::{prf:example} 有两个极限点的序列
:label: ex-02-01-limit-points

求 $\R^2$ 中以下列各项构成的序列的极限点：

$$
\vx^{(k)}=\left(\frac{1}{k+1},\ (-1)^k\,2+\frac{1}{k^2-0.5}\right)^\T,
\qquad k=0,1,2,\dots
$$
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

第一个分量趋于 $0$。第二个分量取决于 $k$ 的奇偶性。当 $k=2s$ 为偶数时，它等于 $2+1/(4s^2-0.5)$，趋于 $2$；当 $k=2s+1$ 为奇数时，它等于 $-2+1/((2s+1)^2-0.5)$，趋于 $-2$。因此偶数指标项构成的子序列收敛于 $(0,2)$，奇数指标项构成的子序列收敛于 $(0,-2)$：例如 $\vx^{(10)}=(0.0909,\ 2.0101)$，$\vx^{(11)}=(0.0833,\ -1.9917)$。

没有其他极限点。任何收敛的子序列都含有无穷多个偶数指标项或无穷多个奇数指标项，而这些项趋近于 $(0,2)$ 或 $(0,-2)$，因此子序列的极限必定是这两点之一。所以极限点是 $(0,2)$ 和 $(0,-2)$。由于有两个极限点，序列本身不收敛。
:::

闭集可以用序列来刻画，而且这种刻画往往比定义更便于使用。

:::{prf:theorem} 闭集包含其中序列的极限
:label: thm-02-01-closed-sequences

集合 $\Omega\subseteq\R^n$ 是闭集，当且仅当由 $\Omega$ 中的点构成的每个收敛序列的极限都属于 $\Omega$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

设 $\Omega$ 是闭集，$\vx^{(k)}\in\Omega$，$\vx^{(k)}\to\vx^*$。若 $\vx^*\notin\Omega$，则 $\vx^*$ 位于开的补集中，因此存在某个球 $B(\vx^*,\varepsilon)$ 不含 $\Omega$ 中的点。但对所有充分大的 $k$ 有 $\norm{\vx^{(k)}-\vx^*}<\varepsilon$，矛盾。

反之，设 $\Omega$ 包含其中收敛序列的极限，并设 $\vy\notin\Omega$。若每个球 $B(\vy,1/k)$ 都含有一点 $\vx^{(k)}\in\Omega$，则 $\vx^{(k)}\to\vy$，从而 $\vy$ 将属于 $\Omega$，矛盾。因此以 $\vy$ 为中心的某个球与 $\Omega$ 不相交，即补集是开集。
:::

### 有界集与紧集

:::{prf:definition} 有界集与紧集
:label: def-02-01-bounded-compact

若存在数 $M$，使得对所有 $\vx\in\Omega$ 有 $\norm{\vx}\le M$（等价地，对所有分量和某个 $M'$ 有 $\abs{x_j}\le M'$），则称集合 $\Omega\subseteq\R^n$ **有界**（bounded）。若它既是闭集又有界，则称它为**紧集**（compact set）。
:::

紧集的关键性质是其中的序列无法“逃逸”。由**波尔查诺–魏尔斯特拉斯定理**（Bolzano–Weierstrass theorem），$\R^n$ 中每个有界序列都至少有一个极限点；若序列位于紧集 $\Omega$ 中，则由{prf:ref}`thm-02-01-closed-sequences`，该极限点也属于 $\Omega$。{prf:ref}`ex-02-01-limit-points`中的序列是有界的，它确实有极限点；序列 $\vx^{(k)}=(k,0)$ 是无界的，它没有极限点。

(sec-02-01-weierstrass)=
### 最小值的存在性

:::{prf:theorem} 魏尔斯特拉斯定理
:label: thm-02-01-weierstrass

设 $\Omega\subseteq\R^n$ 非空且紧，$f:\Omega\to\R$ 连续。则 $f$ 在 $\Omega$ 上取到最小值和最大值：存在点 $\vx^*,\vx^\sharp\in\Omega$，使得对所有 $\vx\in\Omega$ 有 $f(\vx^*)\le f(\vx)\le f(\vx^\sharp)$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

令 $p^*=\inf\{f(\vx):\vx\in\Omega\}$，它可能为 $-\infty$。选取点 $\vx^{(k)}\in\Omega$，使 $f(\vx^{(k)})\to p^*$。由于 $\Omega$ 是紧集，存在一个子序列收敛于某个 $\vx^*\in\Omega$，再由 $f$ 的连续性，沿该子序列有 $f(\vx^*)=\lim f(\vx^{(k_s)})=p^*$。因此 $p^*$ 有限且能取到。对 $-f$ 用同样的论证即得最大值。
:::

如下面三个小例子所示，每个假设都是必需的。在无界集 $[1,\infty)$ 上，函数 $1/x$ 的下确界为 $0$，但永远取不到（{doc}`01-03-introduction`）。在不是闭集的 $(0,1]$ 上，$f(x)=x$ 的下确界为 $0$，取不到。在紧集 $[0,1]$ 上，满足 $f(0)=1$、且当 $x>0$ 时 $f(x)=x$ 的不连续函数的下确界同样为 $0$，取不到。

实际中的许多可行集是闭的但无界，例如 $\{\vx:\vx\ge\vzero\}$。如果目标函数在远处无限增大，定理仍然适用。若当 $\vx\in\Omega$ 且 $\norm{\vx}\to\infty$ 时总有 $f(\vx)\to\infty$，则称 $f$ 在 $\Omega$ 上是**强制的**（coercive）。

:::{prf:corollary} 强制函数能取到最小值
:label: cor-02-01-coercive

若 $\Omega$ 非空且闭，$f$ 在 $\Omega$ 上连续且强制，则 $f$ 在 $\Omega$ 上取到最小值。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

取 $\vx^{(0)}\in\Omega$。集合 $\{\vx\in\Omega:f(\vx)\le f(\vx^{(0)})\}$ 非空、闭（闭集中的一个连续函数不等式）且有界（由强制性），因此是紧集。由{prf:ref}`thm-02-01-weierstrass`，$f$ 在这个集合上存在极小点，它也是 $f$ 在整个 $\Omega$ 上的极小点，因为被排除在外的点的函数值更大。
:::

例如，到固定点 $\vy$ 的距离的平方 $\norm{\vx-\vy}^2$ 是强制的，因此每个非空闭集都含有一个距 $\vy$ 最近的点。这一事实将在{ref}`sec-02-01-separation`中用到。

(sec-02-01-affine-sets)=
## 直线、线段与仿射集

设 $\vx\ne\vy$ 是 $\R^n$ 中的两点。点

$$
\vz=\theta\vx+(1-\theta)\vy=\vy+\theta(\vx-\vy),\qquad\theta\in\R,
$$ (eq-02-01-line)

构成过 $\vx$ 和 $\vy$ 的**直线**（line）：从 $\vy$（$\theta=0$）出发，沿方向 $\vx-\vy$ 移动，在 $\theta=1$ 时到达 $\vx$。把 $\theta$ 限制在 $[0,1]$ 中，就得到 $\vx$ 与 $\vy$ 之间的**线段**（line segment），记为 $[\vx,\vy]$。

:::{prf:definition} 仿射集
:label: def-02-01-affine-set

若集合 $A\subseteq\R^n$ 包含过其中任意两个不同点的直线，则称 $A$ 为**仿射集**（affine set）。满足 $\theta_1+\dots+\theta_k=1$ 的点 $\theta_1\vx_1+\dots+\theta_k\vx_k$ 称为 $\vx_1,\dots,\vx_k$ 的**仿射组合**（affine combination）。
:::

仿射集包含其中各点的每一个仿射组合（对 $k$ 用归纳法，如下面的{prf:ref}`thm-02-01-convex-combinations`）。最重要的例子是线性方程组的解集。

:::{prf:proposition} 线性方程组的解集是仿射集
:label: prop-02-01-linear-equations-affine

对 $\mA\in\R^{m\times n}$ 和 $\vb\in\R^m$，集合 $\{\vx:\mA\vx=\vb\}$ 是仿射集。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

若 $\mA\vx=\vb$，$\mA\vy=\vb$，且 $\theta\in\R$，则 $\mA(\theta\vx+(1-\theta)\vy)=\theta\vb+(1-\theta)\vb=\vb$。
:::

反之亦然：每个仿射集都可以对某个 $\mA$ 和 $\vb$ 写成 $\{\vx:\mA\vx=\vb\}$。若 $A$ 是仿射集且 $\vx_0\in A$，则 $V=A-\vx_0=\{\vx-\vx_0:\vx\in A\}$ 是线性子空间，而 $A=\vx_0+V$ 是平移后的子空间；它的**维数**定义为 $V$ 的维数。用{doc}`01-01-linear-algebra`的语言来说，$V$ 是 $\mA$ 的零空间，维数为 $n-\rank\mA$。集合 $S$ 的**仿射包**（affine hull）$\aff S$ 是 $S$ 中各点的全部仿射组合构成的集合，即包含 $S$ 的最小仿射集。两个不同点的仿射包是过它们的直线；不共线的三个点的仿射包是一个平面。

(sec-02-01-convex-sets)=
## 凸集

:::{prf:definition} 凸集
:label: def-02-01-convex-set

若集合 $C\subseteq\R^n$ 包含其中任意两点之间的线段：

$$
\vx,\vy\in C,\ \theta\in[0,1]\ \Longrightarrow\ \theta\vx+(1-\theta)\vy\in C .
$$ (eq-02-01-convex-set)

则称 $C$ 是**凸的**（convex）。
:::

每个仿射集都是凸的，因为它包含整条直线，而不仅仅是线段。空集、单点、射线 $\{\vx_0+t\vd:t\ge0\}$ 以及 $\R^n$ 都是凸的。{numref}`fig-02-01-convex-nonconvex`给出了一个凸集和三个非凸集：对后者中的每一个，都能在集合中找到两点，使它们之间的线段离开该集合。

```{figure} figures/02-01-convex-nonconvex.svg
:name: fig-02-01-convex-nonconvex
:width: 640px
:class: dark-light
:alt: 一个内部画有线段的阴影六边形，以及三个非凸集（肾形、圆环、在一个角上相接的两个六边形），每个都画有一条离开集合的线段

(a) 凸六边形：其中任意两点之间的线段都留在其中。(b) 肾形集合、(c) 圆环以及 (d) 在一个顶点处相接的两个六边形都不是凸的：线段的虚线部分位于集合之外。
```

### 凸组合与凸包

满足 $\theta_i\ge0$ 且 $\theta_1+\dots+\theta_k=1$ 的点 $\theta_1\vx_1+\dots+\theta_k\vx_k$ 称为 $\vx_1,\dots,\vx_k$ 的**凸组合**（convex combination），即这些点的加权平均。当 $k=2$ 时，凸组合就是线段 $[\vx_1,\vx_2]$ 上的点；对不共线的三个点，凸组合填满以它们为顶点的三角形。

:::{prf:theorem} 凸集包含其中各点的所有凸组合
:label: thm-02-01-convex-combinations

集合 $C$ 是凸的，当且仅当 $C$ 中有限个点的每个凸组合都属于 $C$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

若 $C$ 包含所有凸组合，则它包含任意两点的凸组合，因此是凸的。反之，设 $C$ 是凸的；我们对 $k$ 用归纳法证明 $C$ 中 $k$ 个点的凸组合属于 $C$。$k=1$ 时无需证明，$k=2$ 就是定义。设 $\vz=\sum_{i=1}^k\theta_i\vx_i$，$k\ge3$。若 $\theta_k=1$，则 $\vz=\vx_k\in C$。否则

$$
\vz=(1-\theta_k)\sum_{i=1}^{k-1}\frac{\theta_i}{1-\theta_k}\vx_i+\theta_k\vx_k .
$$

权重 $\theta_i/(1-\theta_k)$，$i<k$，非负且和为 $1$，因此第一个和式是 $k-1$ 个点的凸组合，由归纳假设属于 $C$。于是 $\vz$ 是 $C$ 中两点的凸组合，属于 $C$。
:::

:::{prf:definition} 凸包
:label: def-02-01-convex-hull

集合 $S\subseteq\R^n$ 的**凸包**（convex hull）$\conv S$ 是 $S$ 中有限个点的全部凸组合构成的集合。
:::

$\conv S$ 是凸的，因为两个凸组合的凸组合仍是 $S$ 中各点的凸组合；并且由{prf:ref}`thm-02-01-convex-combinations`，每个包含 $S$ 的凸集都包含 $\conv S$。因此 $\conv S$ 是包含 $S$ 的*最小*凸集。对有限集 $S$，它是一个以 $S$ 中的点为角点的多边形（在 $\R^2$ 中）或多胞形（一般情形）；见{numref}`fig-02-01-hulls`(a)。**Carathéodory 定理**（Carathéodory's theorem）进一步指出，在 $\R^n$ 中，$\conv S$ 的每个点都是 $S$ 中至多 $n+1$ 个点的凸组合；在平面上，凸包中的每个点都位于某个角点属于 $S$ 的三角形中。

```{figure} figures/02-01-hulls.svg
:name: fig-02-01-hulls
:width: 600px
:class: dark-light
:alt: 左图：九个点及作为其凸包的阴影多边形，角点已标出。右图：从原点出发的两个向量及它们之间作为其锥包的阴影楔形

(a) 九个点的凸包是阴影多边形；其角点（实心）是极点，其余各点（空心）位于内部或边上。(b) $\vx_1$ 和 $\vx_2$ 的锥包是全部组合 $\theta_1\vx_1+\theta_2\vx_2$（$\theta_1,\theta_2\ge0$）构成的楔形。
```

(sec-02-01-cones)=
## 锥

:::{prf:definition} 锥与凸锥
:label: def-02-01-cone

若 $\vx\in K$ 且 $\theta\ge0$ 蕴含 $\theta\vx\in K$，则称集合 $K\subseteq\R^n$ 为**锥**（cone）。若它还是凸的，则称为**凸锥**（convex cone）；等价地，对所有 $\vx_1,\vx_2\in K$ 和 $\theta_1,\theta_2\ge0$ 有 $\theta_1\vx_1+\theta_2\vx_2\in K$。
:::

所有 $\theta_i\ge0$（但对其和不加限制）的点 $\theta_1\vx_1+\dots+\theta_k\vx_k$ 称为各 $\vx_i$ 的**锥组合**（conic combination）。**锥包**（conic hull）$\cone S$ 是 $S$ 中各点的全部锥组合构成的集合，即包含 $S$ 的最小凸锥（{numref}`fig-02-01-hulls`(b)）。从几何上看，锥是从原点出发的射线的并；凸锥是没有“缺口”的锥。

例如：

- 非负卦限 $\R^n_+=\{\vx:\vx\ge\vzero\}$；
- 过原点的半空间 $\{\vx:\va^\T\vx\le0\}$；
- 由矩阵各列生成的锥 $\{\mA\vw:\vw\ge\vzero\}$，它是各列 $\mA_{\cdot1},\dots,\mA_{\cdot n}$ 的锥包。这类锥出现在 Farkas 引理（{prf:ref}`thm-03-05-farkas`）以及 Karush–Kuhn–Tucker（KKT）条件的几何形式 {eq}`eq-06-03-kkt-cone` 中，后者要求目标函数的负梯度位于有效约束的梯度所生成的锥中；
- 后面几节中的范数锥和半正定锥。

两条半轴的并 $\{(t,0):t\ge0\}\cup\{(0,t):t\ge0\}$ 是锥，但不是凸锥：缺少 $(1,0)+(0,1)$。

(sec-02-01-hyperplanes)=
## 超平面与半空间

:::{prf:definition} 超平面与半空间
:label: def-02-01-hyperplane

对 $\va\in\R^n$，$\va\ne\vzero$，以及 $b\in\R$，集合 $H=\{\vx:\va^\T\vx=b\}$ 称为以 $\va$ 为**法向量**（normal vector）的**超平面**（hyperplane）。集合 $\{\vx:\va^\T\vx\le b\}$ 和 $\{\vx:\va^\T\vx\ge b\}$ 是以 $H$ 为边界的两个闭**半空间**（half-space）；换成严格不等式，它们就是开半空间。
:::

在 $\R^2$ 中超平面是直线，在 $\R^3$ 中超平面是平面。若 $\vx_0$ 是 $H$ 上的任意一点，即 $\va^\T\vx_0=b$，则

$$
\begin{aligned}
H&=\{\vx:\va^\T(\vx-\vx_0)=0\},\\
\{\vx:\va^\T\vx\le b\}&=\{\vx:\va^\T(\vx-\vx_0)\le0\}.
\end{aligned}
$$ (eq-02-01-hyperplane-offset)

第一个集合由使 $\vx-\vx_0$ 与 $\va$ 正交的点 $\vx$ 组成：$H$ 是维数为 $n-1$ 的仿射集 $\vx_0+\{\va\}^\perp$。半空间 $\va^\T\vx\le b$ 由使 $\vx-\vx_0$ 与 $\va$ 成钝角或直角的点组成，即位于 $H$ *背离* $\va$ 一侧的点（{numref}`fig-02-01-hyperplane`）。从 $H$ 出发沿 $\va$ 的方向移动，$\va^\T\vx$ 增大；这就是“线性函数 $\va^\T\vx$ 的梯度 $\va$ 指向函数增大的方向”这一命题的几何内容（{prf:ref}`thm-01-02-steepest-ascent`）。

```{figure} figures/02-01-hyperplane.svg
:name: fig-02-01-hyperplane
:width: 460px
:class: dark-light
:alt: 平面上的一条直线，在其上一点 x0 处画出法向量 a；背离 a 的一侧标注为 a 的转置乘 x 至多为 b，另一侧至少为 b

过 $\vx_0$、以 $\va$ 为法向量的超平面 $\va^\T\vx=b$。对阴影半空间 $\va^\T\vx\le b$ 中的点 $\vx$，向量 $\vx-\vx_0$ 与 $\va$ 成钝角。
```

超平面是仿射集，因而是凸的。半空间是凸的：若 $\va^\T\vx\le b$，$\va^\T\vy\le b$，且 $\theta\in[0,1]$，则 $\va^\T(\theta\vx+(1-\theta)\vy)=\theta\va^\T\vx+(1-\theta)\va^\T\vy\le\theta b+(1-\theta)b=b$。半空间是线性规划（linear program, LP）可行集的基本构件：每个线性不等式约束定义一个半空间。

(sec-02-01-polyhedra)=
## 多面体

:::{prf:definition} 多面体与多胞形
:label: def-02-01-polyhedron

**多面体**（polyhedron）是有限个闭半空间的交：

$$
P=\{\vx\in\R^n:\va_i^\T\vx\le b_i,\ i=1,\dots,m\}=\{\vx:\mA\vx\le\vb\},
$$ (eq-02-01-polyhedron)

其中 $\mA$ 的各行为 $\va_i^\T$，向量之间的 $\le$ 按分量成立。有界的多面体称为**多胞形**（polytope）。
:::

等式约束是允许的，因为每个等式都是一对不等式：$\va^\T\vx=b$ 成立当且仅当 $\va^\T\vx\le b$ 且 $-\va^\T\vx\le-b$。因此，标准形线性规划的可行集 $\{\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$ 是多面体，由有限个线性等式和不等式描述的每个集合也都是多面体。其他例子有非负卦限 $\R^n_+$、长方体 $\{\vx:l_j\le x_j\le u_j\}$、半空间、仿射集，以及**单位单纯形**（unit simplex）$\{\vx:\vx\ge\vzero,\ \vone^\T\vx=1\}$（概率向量的集合）。

每个多面体都是凸的（作为凸集的交，见{ref}`sec-02-01-operations`），也是闭的（作为闭集的交）。因此多胞形是紧集，由{prf:ref}`thm-02-01-weierstrass`，连续函数（特别是线性函数）在非空多胞形上能取到最小值。

多面体还有一种“内部”描述。**Minkowski–Weyl 定理**（Minkowski–Weyl theorem）指出：一个集合是多胞形，当且仅当它是有限个点的凸包；并且每个多面体都是有限个点的凸包与有限个方向的锥包之和。线性规划基本定理（{prf:ref}`thm-03-02-fundamental`）表明，有最优解的标准形线性规划必有一个最优角点。

(sec-02-01-balls)=
## 球、椭球与范数锥

以 $\vx_c$ 为中心、半径为 $r>0$ 的欧氏闭球为

$$
\bar B(\vx_c,r)=\{\vx:\norm{\vx-\vx_c}\le r\}=\{\vx_c+r\vu:\norm{\vu}\le1\}.
$$ (eq-02-01-ball)

它是凸的：若 $\norm{\vx-\vx_c}\le r$，$\norm{\vy-\vx_c}\le r$，且 $\theta\in[0,1]$，则由三角不等式得

$$
\begin{aligned}
\norm{\theta\vx+(1-\theta)\vy-\vx_c}
&=\norm{\theta(\vx-\vx_c)+(1-\theta)(\vy-\vx_c)}\\
&\le\theta\norm{\vx-\vx_c}+(1-\theta)\norm{\vy-\vx_c}\\
&\le r .
\end{aligned}
$$

证明只用到了三角不等式和齐次性，因此*任何*范数的球 $\{\vx:\norm{\vx-\vx_c}\le r\}$ 都是凸的。在 $\R^2$ 中，$\norm{\cdot}_1$ 的单位球是菱形，$\norm{\cdot}_\infty$ 的单位球是正方形（{numref}`fig-02-01-balls`(a)）。两者都是多胞形：$\ell_1$ 球是 $2n$ 个点 $\pm\ve_j$ 的凸包，$\ell_\infty$ 球是有 $2^n$ 个角点的长方体 $\{-1\le x_j\le1\}$。欧氏球不是多面体：它的边界是弯曲的，有限个半空间无法构成它。

:::{prf:definition} 椭球
:label: def-02-01-ellipsoid

对 $\vx_c\in\R^n$ 和对称正定矩阵 $\mP\in\Sym^n_{++}$，集合

$$
\mathcal E=\{\vx:(\vx-\vx_c)^\T\mP^{-1}(\vx-\vx_c)\le1\}
$$ (eq-02-01-ellipsoid)

称为以 $\vx_c$ 为中心的**椭球**（ellipsoid）。
:::

当 $\mP=r^2\mI$ 时，椭球就是球 $\bar B(\vx_c,r)$。一般地，设 $\mP=\mV\diag(\lambda_1(\mP),\dots,\lambda_n(\mP))\mV^\T$ 是 $\mP$ 的谱分解（{doc}`01-01-linear-algebra`）。$\mathcal E$ 的半轴沿特征向量（$\mV$ 的各列）方向，长度为 $\sqrt{\lambda_i(\mP)}$（{numref}`fig-02-01-balls`(b)）。等价地，$\mathcal E=\{\vx_c+\mP^{1/2}\vu:\norm{\vu}\le1\}$ 是单位球在一个仿射映射下的像，由{prf:ref}`thm-02-01-affine-maps`即证明了它的凸性。例如，椭球可以描述需求估计值的不确定区域；它们也是凸二次函数的下水平集（{prf:ref}`thm-02-02-sublevel`）。

```{figure} figures/02-01-balls.svg
:name: fig-02-01-balls
:width: 600px
:class: dark-light
:alt: 左图：叠画在一起的 1-范数（菱形）、2-范数（圆盘）和无穷范数（正方形）的单位球。右图：一个椭圆及其中心和两条相互垂直的半轴

(a) $\R^2$ 中 $\norm{\cdot}_1$（菱形）、$\norm{\cdot}_2$（圆盘）和 $\norm{\cdot}_\infty$（正方形）的单位球。(b) 椭球 {eq}`eq-02-01-ellipsoid`，其中 $\vx_c=(2,1)$，$\mP$ 的特征值为 $4$ 和 $1$；半轴为 $2\vv_1$ 和 $\vv_2$，其中 $\vv_1$ 和 $\vv_2$ 是 $\mP$ 对应于这两个特征值的单位特征向量。
```

把半径也作为变量，就得到一个锥。范数 $\norm{\cdot}$ 的**范数锥**（norm cone）为

$$
K=\{(\vx,t)\in\R^{n+1}:\norm{\vx}\le t\}.
$$ (eq-02-01-norm-cone)

它是凸锥：用 $\theta\ge0$ 缩放 $(\vx,t)$ 保持 $\norm{\vx}\le t$ 不变，而凸性与球的情形一样由三角不等式推出。对欧氏范数，$K$ 称为**二阶锥**（second-order cone，又称“冰淇淋锥”）；在 $\R^3$ 中，它是实心圆锥 $\sqrt{x_1^2+x_2^2}\le t$。含有形如 $\norm{\mA\vx+\vb}\le\vc^\T\vx+d$ 的约束的最优化问题（二阶锥规划）几乎可以像线性规划一样高效地求解（Boyd & Vandenberghe, 2004, §4.4.2）。

(sec-02-01-psd-cone)=
## 半正定锥

矩阵的集合也可以是凸的。设 $\Sym^n$ 为 $n\times n$ 对称矩阵的集合，$\Sym^n_+$ 为其中半正定矩阵的集合（$\mA\succeq0$，即对所有 $\vz$ 有 $\vz^\T\mA\vz\ge0$），$\Sym^n_{++}$ 为其中正定矩阵的集合（{doc}`01-01-linear-algebra`）。

$\Sym^n_+$ 是凸锥：若 $\mA,\mB\succeq0$ 且 $\theta_1,\theta_2\ge0$，则对每个 $\vz$ 有

$$
\vz^\T(\theta_1\mA+\theta_2\mB)\vz=\theta_1\vz^\T\mA\vz+\theta_2\vz^\T\mB\vz\ge0 .
$$

当 $n=2$ 时，矩阵 $\begin{pmatrix}u&v\\v&w\end{pmatrix}$ 半正定当且仅当 $u\ge0$、$w\ge0$ 且 $uw-v^2\ge0$。在坐标 $(u,v,w)\in\R^3$ 中，这个集合是一个边界弯曲的实心锥，因此 $\Sym^2_+$ 是凸锥，但不是多面体。形如“一个仿射依赖于 $\vx$ 的矩阵半正定”的约束定义了半定规划（Boyd & Vandenberghe, 2004, §4.6.2）；在本书中，这个锥主要通过凸函数的二阶判别法出现（{prf:ref}`thm-02-02-second-order`）。

(sec-02-01-operations)=
## 保凸运算

直接验证定义 {eq}`eq-02-01-convex-set` 往往很麻烦。通常更容易的做法是：用保持凸性的运算，从简单的凸集构造出所要的集合。

:::{prf:theorem} 交集
:label: thm-02-01-intersection

任意一族（有限或无限个）凸集的交是凸的。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

设 $\vx,\vy$ 属于这一族中的每个集合，$\theta\in[0,1]$。每个集合都是凸的，因此每个集合都包含 $\theta\vx+(1-\theta)\vy$，交集也包含它。
:::

多面体是半空间的交，因而是凸的。该定理也适用于无限个集合的交：$\Sym^n_+=\bigcap_{\vz\in\R^n}\{\mA\in\Sym^n:\vz^\T\mA\vz\ge0\}$ 是 $\Sym^n$ 中无穷多个半空间的交，每个 $\vz$ 对应一个。只要每个集合 $\{\vx:g_i(\vx)\le0\}$ 都是凸的，问题的可行集 $\{\vx:g_i(\vx)\le0,\ i=1,\dots,m\}$ 就是凸的；当 $g_i$ 是凸函数时情况正是如此（{prf:ref}`thm-02-02-sublevel`）。

凸集的**并**一般*不是*凸的；{prf:ref}`ex-02-01-seven-sets`中的集合 $\Omega_3$ 是两个凸多边形的并，它不是凸的。

:::{prf:theorem} 仿射映射
:label: thm-02-01-affine-maps

设 $\phi(\vx)=\mA\vx+\vb$，其中 $\mA\in\R^{m\times n}$，$\vb\in\R^m$。若 $C\subseteq\R^n$ 是凸集，则它的像 $\phi(C)=\{\mA\vx+\vb:\vx\in C\}$ 是凸的。若 $D\subseteq\R^m$ 是凸集，则它的原像 $\phi^{-1}(D)=\{\vx:\mA\vx+\vb\in D\}$ 是凸的。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

仿射映射保持凸组合：由于 $\theta\vb+(1-\theta)\vb=\vb$，有 $\phi(\theta\vx+(1-\theta)\vy)=\theta\phi(\vx)+(1-\theta)\phi(\vy)$。若 $\phi(\vx),\phi(\vy)\in\phi(C)$，其中 $\vx,\vy\in C$，则组合 $\theta\phi(\vx)+(1-\theta)\phi(\vy)$ 是点 $\theta\vx+(1-\theta)\vy\in C$ 的像。若 $\phi(\vx),\phi(\vy)\in D$，则由于 $D$ 是凸的，$\phi(\theta\vx+(1-\theta)\vy)\in D$。
:::

推论：

- 凸集的平移 $C+\vx_0$、缩放 $\alpha C$ 和旋转都是凸的；
- 凸集在其部分坐标上的投影 $\{\vx_1:(\vx_1,\vx_2)\in C\}$ 是凸的（线性映射下的像）；
- 两个凸集之和 $C_1+C_2=\{\vx+\vy:\vx\in C_1,\vy\in C_2\}$ 是凸的（凸集 $C_1\times C_2$ 在 $(\vx,\vy)\mapsto\vx+\vy$ 下的像）；
- 集合 $\{\vx:\norm{\mA\vx-\vb}\le r\}$ 是凸的（球的原像）；
- 椭球 {eq}`eq-02-01-ellipsoid` 是凸的（单位球的像）。

另外两种运算，即透视映射和线性分式映射，也保持凸性（Boyd & Vandenberghe, 2004, §2.3.3）；本书不需要用到它们。

(sec-02-01-seven-sets)=
## 七个示例集合

下面的例子把本章的定义应用于一族集合，其中五个位于平面上（{numref}`fig-02-01-example-sets`）。

```{figure} figures/02-01-example-sets.svg
:name: fig-02-01-example-sets
:width: 720px
:class: dark-light
:alt: 五个小图，分别显示四边形 Omega 1、圆盘 Omega 2、带有一条位于集合外的虚线线段的非凸四边形 Omega 3、有两条虚线开边的 Omega 4，以及半平面 Omega 5

{prf:ref}`ex-02-01-seven-sets` 中的集合 $\Omega_1$ 至 $\Omega_5$。$\Omega_4$ 的虚线边不属于该集合。在 $\Omega_3$ 中，从 $(0,0)$ 到 $(3,0)$ 的线段离开了集合，因此 $\Omega_3$ 不是凸的。$\Omega_5$ 向右上方无限延伸。
```

:::{prf:example} 七个集合的性质
:label: ex-02-01-seven-sets

对每个集合，判断它是否为开集、闭集、有界集、紧集、凸集和多面体。

$$
\begin{aligned}
\Omega_1&=\{\vx\in\R^2:x_1+x_2\le2,\ x_1-x_2\le1,\ x_1\ge0,\ x_2\ge0\},\\
\Omega_2&=\{\vx\in\R^2:(x_1-1)^2+(x_2-1)^2\le1\},\\
\Omega_3&=\{\vx\in\R^2:x_1+x_2\le3,\ g(\vx)\le0,\ x_1\ge0,\ x_2\ge0\},\\
\Omega_4&=\{\vx\in\R^2:x_1+x_2<2,\ x_1-x_2<1,\ x_1\ge0,\ x_2\ge0\},\\
\Omega_5&=\{\vx\in\R^2:x_1+x_2\ge2\},\\
\Omega_6&=\R^n_+=\{\vx\in\R^n:x_j\ge0,\ j=1,\dots,n\},\\
\Omega_7&=\R^n,
\end{aligned}
$$

其中当 $0\le x_1\le1$ 时 $g(\vx)=x_1-x_2$，当 $x_1\ge1$ 时 $g(\vx)=3-x_1-2x_2$。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**$\Omega_1$** 是以 $(0,0)$、$(1,0)$、$(\tfrac32,\tfrac12)$ 和 $(0,2)$ 为角点的四边形；角点 $(\tfrac32,\tfrac12)$ 是 $x_1+x_2=2$、$x_1-x_2=1$ 的解。它是四个闭半空间的交，因此是多面体，是闭的、凸的。它有界（其上 $\norm{\vx}\le2$），因此是多胞形，也是紧集。它不是开集：边界点 $(0,0)$ 属于它。

**$\Omega_2$** 是以 $(1,1)$ 为圆心、半径为 $1$、与两条坐标轴都相切的闭圆盘。它是闭的、有界的、紧的和凸的（{ref}`sec-02-01-balls`），不是开集，也不是多面体，因为它的边界是圆。

**$\Omega_3$** 上方以直线 $x_1+x_2=3$ 为界，下方以从 $(0,0)$ 经 $(1,1)$ 到 $(3,0)$ 的折线（即 $g(\vx)=0$ 处）为界：它是以 $(0,0)$、$(1,1)$、$(3,0)$、$(0,3)$ 为角点的四边形，在 $(1,1)$ 处的内角约为 $252^\circ$。它是闭的（是两个闭多边形，即 $x_1\le1$ 部分和 $x_1\ge1$ 部分的并）且有界，因此是紧集。它**不是凸的**：$(0,0)$ 和 $(3,0)$ 属于 $\Omega_3$，但它们的中点 $(1.5,0)$ 不属于，因为在 $x_1=1.5$ 处该集合要求 $x_2\ge(3-1.5)/2=0.75$。因此它不是多面体。它不是开集。

**$\Omega_4$** 是去掉两条斜边的 $\Omega_1$。它是凸的（作为两个开半空间和两个闭半空间的交），并且有界。它既不是开集也不是闭集。它不是开集，因为点 $(0.5,0)$ 属于它，但以该点为中心的每个球都含有 $x_2<0$ 的点。它不是闭集，因为点 $(\tfrac32-\tfrac1k,\tfrac12)$，$k=2,3,\dots$，属于 $\Omega_4$，而它们的极限 $(\tfrac32,\tfrac12)$ 不属于（{prf:ref}`thm-02-01-closed-sequences`）。因此它不是紧集，也不是多面体（多面体是闭的）。它的闭包是 $\Omega_1$。

**$\Omega_5$** 是闭半平面：闭的、凸的、是多面体、不是开集。它无界（对所有 $t\ge1$ 包含 $(t,t)$），因此不是紧集。

**$\Omega_6$** 即非负卦限，是闭的、凸的、是多面体，也是凸锥；它无界，不是开集。

**$\Omega_7=\R^n$** 既是开集又是闭集，是凸的，且无界。它是多面体，因为它具有形式 {eq}`eq-02-01-polyhedron`，其中不等式个数 $m=0$（或者只有一行 $\vzero^\T\vx\le0$）；它也是仿射集和锥。

{numref}`tab-02-01-seven-sets`汇总了这些答案。

```{list-table} 七个示例集合 $\Omega_1$ 至 $\Omega_7$ 的性质
:name: tab-02-01-seven-sets
:header-rows: 1

* - 集合
  - 开集
  - 闭集
  - 有界
  - 紧集
  - 凸集
  - 多面体
* - $\Omega_1$
  - 否
  - 是
  - 是
  - 是
  - 是
  - 是（多胞形）
* - $\Omega_2$
  - 否
  - 是
  - 是
  - 是
  - 是
  - 否
* - $\Omega_3$
  - 否
  - 是
  - 是
  - 是
  - 否
  - 否
* - $\Omega_4$
  - 否
  - 否
  - 是
  - 否
  - 是
  - 否
* - $\Omega_5$
  - 否
  - 是
  - 否
  - 否
  - 是
  - 是
* - $\Omega_6$
  - 否
  - 是
  - 否
  - 否
  - 是
  - 是
* - $\Omega_7$
  - 是
  - 是
  - 否
  - 否
  - 是
  - 是
```
:::

由{prf:ref}`thm-02-01-weierstrass`和{numref}`tab-02-01-seven-sets`中的答案，每个连续的目标函数都在 $\Omega_1$、$\Omega_2$ 和 $\Omega_3$ 上取到最小值。在 $\Omega_4$ 上则未必：$f(\vx)=-x_1$ 趋近于 $-\tfrac32$，但永远取不到。在 $\Omega_5$ 至 $\Omega_7$ 上，最小值是否存在取决于目标函数，如{prf:ref}`cor-02-01-coercive`所示。

(sec-02-01-extreme-points)=
## 极点

:::{prf:definition} 极点
:label: def-02-01-extreme-point

若凸集 $C$ 中的点 $\vx$ 不是 $C$ 中另外两点的凸组合，即只要 $\vx=\theta\vy+(1-\theta)\vz$，其中 $\vy,\vz\in C$，$\theta\in(0,1)$，就有 $\vy=\vz=\vx$，则称 $\vx$ 为 $C$ 的**极点**（extreme point）。
:::

极点不严格位于 $C$ 中任何线段的内部。多边形的极点是它的角点；有限集的凸包的极点是凸包的角点（{numref}`fig-02-01-hulls`(a)）。圆盘的每个边界点都是极点，因此 $\Omega_2$ 有无穷多个极点。半空间、直线和 $\R^n$ 都没有极点。卦限 $\R^n_+$ 恰有一个极点，即原点。在{prf:ref}`ex-02-01-seven-sets`中，$\Omega_1$ 有四个极点 $(0,0)$、$(1,0)$、$(\tfrac32,\tfrac12)$、$(0,2)$，而 $\Omega_4$ 只保留了 $(0,0)$：其余三个角点位于被去掉的边上。

极点之所以重要，是因为它们足以描述一个紧凸集，也足以求出线性函数在其上的最小值。Minkowski 定理（在无穷维情形下为 Krein–Milman 定理）指出，$\R^n$ 中的紧凸集是其极点的凸包。因此，若 $\vx^*$ 在这样的集合上极小化 $\vc^\T\vx$，且 $\vx^*=\sum_i\theta_i\vv_i$ 是极点 $\vv_i$ 的凸组合，则

$$
\vc^\T\vx^*=\sum_i\theta_i\,\vc^\T\vv_i\ \ge\ \min_i\vc^\T\vv_i ,
$$

因此某个极点也是最优的。对于线性规划的可行集（它们是多面体，可能无界），{doc}`03-02-lp-geometry`证明了相应的结论（{prf:ref}`thm-03-02-fundamental`），并把极点等同于单纯形法所处理的基本可行解（{prf:ref}`thm-03-02-extreme-bfs`）。

(sec-02-01-separation)=
## 分离超平面与支撑超平面

凸集与其外部的一点总可以用一个超平面分离。这个看似简单的事实是对偶理论（{prf:ref}`thm-02-03-slater`）、Farkas 引理（{prf:ref}`thm-03-05-farkas`）和 KKT 条件（{prf:ref}`thm-06-03-kkt`）的基础。证明要用到集合中的最近点。

:::{prf:theorem} 到闭凸集上的投影
:label: thm-02-01-projection

设 $C\subseteq\R^n$ 非空、闭且凸，$\vy\in\R^n$。则 $C$ 中恰有一点 $\vx^*$ 距 $\vy$ 最近，记为 $\vx^*=\operatorname{proj}_C(\vy)$；并且点 $\vx^*\in C$ 是这个点当且仅当

$$
(\vy-\vx^*)^\T(\vx-\vx^*)\le0\quad\text{对所有 }\vx\in C .
$$ (eq-02-01-projection)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

*存在性。* 函数 $\norm{\vx-\vy}^2$ 连续且强制，因此由{prf:ref}`cor-02-01-coercive`，它在闭集 $C$ 上取到最小值。

*刻画。* 设 $\vx^*$ 是一个最近点，$\vx\in C$，$t\in(0,1]$。由凸性，点 $\vx^*+t(\vx-\vx^*)$ 属于 $C$，因此

$$
\begin{aligned}
\norm{\vy-\vx^*}^2
&\le\norm{\vy-\vx^*-t(\vx-\vx^*)}^2\\
&=\norm{\vy-\vx^*}^2-2t(\vy-\vx^*)^\T(\vx-\vx^*)+t^2\norm{\vx-\vx^*}^2 .
\end{aligned}
$$

消去相同的项，除以 $2t$，再令 $t\to0$，即得 {eq}`eq-02-01-projection`。反之，若 $\vx^*\in C$ 满足 {eq}`eq-02-01-projection`，则对每个 $\vx\in C$ 有

$$
\begin{aligned}
\norm{\vy-\vx}^2&=\norm{\vy-\vx^*}^2+\norm{\vx^*-\vx}^2-2(\vy-\vx^*)^\T(\vx-\vx^*)\\
&\ge\norm{\vy-\vx^*}^2 ,
\end{aligned}
$$

因此 $\vx^*$ 是最近点。

*唯一性。* 若 $\vx^*$ 和 $\vx^{**}$ 都满足 {eq}`eq-02-01-projection`，则对 $\vx^*$ 取 $\vx=\vx^{**}$、对 $\vx^{**}$ 取 $\vx=\vx^*$ 应用该条件，再相加得 $\norm{\vx^{**}-\vx^*}^2\le0$。
:::

条件 {eq}`eq-02-01-projection` 说的是：$\vy-\vx^*$ 与每个指向集合内部的方向 $\vx-\vx^*$ 之间的夹角至少为 $90^\circ$。其中已经蕴含了分离超平面。

:::{prf:theorem} 点与闭凸集的分离
:label: thm-02-01-separation-point

设 $C\subseteq\R^n$ 非空、闭且凸，$\vy\notin C$。则存在 $\va\ne\vzero$ 和 $b$，使得

$$
\va^\T\vx<b<\va^\T\vy\quad\text{对所有 }\vx\in C .
$$ (eq-02-01-strict-separation)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

令 $\vx^*=\operatorname{proj}_C(\vy)$，$\va=\vy-\vx^*$；由于 $\vy\notin C$，$\va$ 非零。取过 $\vx^*$ 与 $\vy$ 中点的超平面，$b=\va^\T(\vx^*+\vy)/2$。则 $\va^\T\vy=b+\tfrac12\norm{\va}^2>b$，且 $\va^\T\vx^*=b-\tfrac12\norm{\va}^2$。对 $\vx\in C$，由 {eq}`eq-02-01-projection` 得 $\va^\T(\vx-\vx^*)\le0$，因此 $\va^\T\vx\le\va^\T\vx^*=b-\tfrac12\norm{\va}^2<b$。
:::

{numref}`fig-02-01-separation`(a) 展示了这一构造。由此可知，闭凸集是包含它的所有闭半空间的交：集合外的每个点都被其中某个半空间切除。这就是凸集的“外部”描述，多面体是其中半空间个数有限的情形。

```{figure} figures/02-01-separation.svg
:name: fig-02-01-separation
:width: 640px
:class: dark-light
:alt: 左图：圆盘 C、圆盘外的一点 y、圆周上的最近点 x-star，以及二者之间一条与 y 减 x-star 垂直的直线。右图：一个凸多边形，一条在边界点 x0 处与之相切的直线，以及背离集合的法向量 a

(a) 圆盘 $C=\bar B((1,1),1)$ 与点 $\vy=(3,2)$ 的分离超平面：法向量为 $\va=\vy-\operatorname{proj}_C(\vy)$，超平面经过虚线线段的中点。(b) 边界点 $\vx_0$ 处的支撑超平面：整个集合位于半空间 $\va^\T\vx\le\va^\T\vx_0$ 中。
```

两个不相交的凸集也可以分离。

:::{prf:theorem} 分离超平面定理
:label: thm-02-01-separating-hyperplane

若 $C$ 和 $D$ 非空、凸且不相交，则存在 $\va\ne\vzero$ 和 $b$，使得对所有 $\vx\in C$ 有 $\va^\T\vx\le b$，对所有 $\vx\in D$ 有 $\va^\T\vx\ge b$。若此外 $C$ 是闭集、$D$ 是紧集，则不等式可以取为严格的。
:::

严格的情形由{prf:ref}`thm-02-01-separation-point`推出：集合 $C-D=\{\vx-\vy:\vx\in C,\vy\in D\}$ 是凸的（{prf:ref}`thm-02-01-affine-maps`）、闭的（因为 $D$ 是紧集），且不含 $\vzero$，因此 $\vzero$ 可以与它分离。一般情形需要一个取极限的论证（Boyd & Vandenberghe, 2004, §2.5.1）。对于不是紧集的闭集，严格分离可能不成立：闭凸集 $\{\vx:x_2\le0\}$ 和 $\{\vx:x_1>0,\ x_1x_2\ge1\}$ 不相交，但它们可以任意接近，唯一的分离超平面是 $x_2=0$，它与两者都不接触。

:::{prf:theorem} 支撑超平面定理
:label: thm-02-01-supporting-hyperplane

设 $C$ 是凸集，$\vx_0\in\operatorname{bd}C$。则存在 $\va\ne\vzero$，使得对所有 $\vx\in C$ 有 $\va^\T\vx\le\va^\T\vx_0$。超平面 $\{\vx:\va^\T\vx=\va^\T\vx_0\}$ 称为 $C$ 在 $\vx_0$ 处的**支撑超平面**（supporting hyperplane）。
:::

:::{dropdown} 证明思路
:class-container: proof
:animate: fade-in

由于 $\vx_0$ 是凸集 $C$ 的边界点，存在收敛于 $\vx_0$ 的点列 $\vy^{(k)}\notin\operatorname{cl}C$。由{prf:ref}`thm-02-01-separation-point`，用长度为 1 的法向量 $\va^{(k)}$ 把每个 $\vy^{(k)}$ 与 $\operatorname{cl}C$ 分离。有界序列 $\{\va^{(k)}\}$ 有一个极限点 $\va$（波尔查诺–魏尔斯特拉斯定理），$\norm{\va}=1$；在 $\va^{(k)\T}\vx<\va^{(k)\T}\vy^{(k)}$ 中沿相应的子序列取极限，即得对所有 $\vx\in C$ 有 $\va^\T\vx\le\va^\T\vx_0$。
:::

从几何上看，支撑超平面在 $\vx_0$ 处与 $C$ 接触，且 $C$ 完全位于它的一侧（{numref}`fig-02-01-separation`(b)）。在边界光滑处，它就是切平面；在角点处，支撑超平面有许多个。在最优化中，若 $\vx^*$ 是 $\vc^\T\vx$ 在 $C$ 上的极小点，则过 $\vx^*$ 的水平集 $\{\vx:\vc^\T\vx=\vc^\T\vx^*\}$ 是 $C$ 在 $\vx^*$ 处的支撑超平面，其中 $\va=-\vc$；{prf:ref}`thm-02-03-optimality`的最优性条件把这一图景推广到非线性目标函数。

## 本章回顾

### 小结

- 开集、闭集、有界集和紧集都通过球来定义（{prf:ref}`def-02-01-open-closed`、{prf:ref}`def-02-01-bounded-compact`）；一个集合是闭集当且仅当它包含其中收敛序列的极限（{prf:ref}`thm-02-01-closed-sequences`）。序列的极限点是其收敛子序列的极限（{prf:ref}`ex-02-01-limit-points`）。
- 连续函数在非空紧集上取到最小值（{prf:ref}`thm-02-01-weierstrass`）；若它是强制的，则在闭集上也取到最小值（{prf:ref}`cor-02-01-coercive`）。
- 仿射集包含直线，凸集包含线段（{prf:ref}`def-02-01-convex-set`），凸锥包含锥组合（{prf:ref}`def-02-01-cone`）；$\conv S$ 是包含 $S$ 的最小凸集。
- 超平面、半空间、多面体 $\{\vx:\mA\vx\le\vb\}$、范数球、椭球、范数锥以及锥 $\Sym^n_+$ 都是凸的。
- 交集以及仿射映射下的像和原像保持凸性（{prf:ref}`thm-02-01-intersection`、{prf:ref}`thm-02-01-affine-maps`）；并集则不然。
- {prf:ref}`ex-02-01-seven-sets`中的七个集合说明了每一种性质（{numref}`tab-02-01-seven-sets`）。
- 极点（{prf:ref}`def-02-01-extreme-point`）生成紧凸集，线性函数在这样的集合上的最小值在某个极点处取到。
- 利用投影（{prf:ref}`thm-02-01-projection`），闭凸集外的一点可以与该集合严格分离（{prf:ref}`thm-02-01-separation-point`）；凸集的每个边界点处都有支撑超平面（{prf:ref}`thm-02-01-supporting-hyperplane`）。

### 常见错误

- **把“闭”理解为“不开”。** 大多数集合，如 $\Omega_4$ 或 $[0,1)$，既不开也不闭；$\R^n$ 则既开又闭。
- **忘记并集不保持凸性。** 交集保持凸性；并集（如 $\Omega_3$）一般不保持。
- **把每个有角的凸集都称为多面体。** 多面体是*有限个闭*半空间的交；圆盘不是多面体，不是闭集的 $\Omega_4$ 也不是。
- **把魏尔斯特拉斯定理用于无界集或开集。** 闭但无界的集合上的连续函数需要额外的论证，例如强制性（{prf:ref}`cor-02-01-coercive`）。
- **混淆边界点与极点。** 紧凸集的每个极点都是边界点，但边中间的点是边界点而不是极点。
- **混淆极限点与极限。** 一个序列可以有多个极限点（{prf:ref}`ex-02-01-limit-points`）；只有当它有界且恰有一个极限点时，它才收敛。

## 参考文献

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 2: convex
  hulls, closure and interior, Weierstrass' theorem, separation and support
  of sets, convex cones, polyhedral sets and extreme points; Appendix A,
  mathematical review.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Sections 2.1–2.2, polyhedra and extreme
  points.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Chapter 2,
  Sections 2.1–2.5: affine and convex sets, important examples, operations
  that preserve convexity, separating and supporting hyperplanes; Appendix
  A.2, topology.)
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming*
  (4th ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3>
  (Appendix B, convex sets.)
