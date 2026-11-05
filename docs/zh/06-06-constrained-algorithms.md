# 约束优化算法

## 概述

{doc}`06-05-unconstrained-algorithms`中的方法在 $\R^n$ 中自由移动。有了约束，就会出现两个问题：最速下降方向可能离开可行集，而 $f$ 的驻点可能根本不可行。本章介绍绕过这些困难的两种经典途径。第一种途径使每个迭代点都保持可行，并沿**可行下降方向**（feasible descent direction）移动；其最简单的代表 **Frank–Wolfe 方法**（Frank–Wolfe method）通过在可行多面体上求解一个*线性规划问题*（linear program, LP）来确定方向，因此每次迭代把{doc}`03-03-simplex-method`中的单纯形法与{doc}`06-04-line-search`中的线搜索结合起来。该方法还在每一步给出最优值的一个下界，但只是次线性收敛。第二种途径把约束违反量的**罚项**加到目标函数中，从而去掉约束，并在罚参数 $\rho$ 不断增大的情况下求解一系列无约束问题。极小点 $\vx(\rho)$ 从可行集外部逼近解，罚项给出{doc}`06-03-kkt-conditions`中 Karush–Kuhn–Tucker（KKT）乘子的估计，代价是病态程度不断加剧。本章最后概述障碍函数法和增广拉格朗日法，它们弥补了这一弱点的部分方面。全章用到凸集和凸函数（{doc}`02-01-convex-sets`，{doc}`02-02-convex-functions`）以及{doc}`06-02-optimality-conditions`中的可行方向和下降方向。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 解释为什么负梯度方向可能不可行，以及什么使一个方向既可行又下降；
- 建立并求解 Frank–Wolfe 方法的线性子问题，把其最优值用作最优性检验和下界，并在 $[0,1]$ 上进行线搜索；
- 叙述并证明 Frank–Wolfe 方法的 $O(1/k)$ 收敛界，并识别其在顶点之间的锯齿现象；
- 构造非线性规划问题（nonlinear program, NLP）的二次罚函数，计算极小点 $\vx(\rho)$，并跟踪它们在 $\rho\to\infty$ 时的变化；
- 由罚函数法的解估计 KKT 乘子，并解释大罚参数造成的病态；
- 描述对数障碍函数和增广拉格朗日函数，以及它们与二次罚函数的区别。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\vy^{(k)}$ | 第 $k$ 次迭代时 Frank–Wolfe 线性子问题的解 |
| $\Gamma(\vx)$ | Frank–Wolfe 间隙，$\max_{\vy\in\Omega}\nabla f(\vx)^\T(\vx-\vy)\ge0$ |
| $D$, $L_f$ | $\Omega$ 的直径；$\nabla f$ 的 Lipschitz 常数 |
| $\psi(\vx)$ | 罚函数（不可行性的度量），$\psi(\vx)=0$ 当且仅当 $\vx\in\Omega$ |
| $\rho$, $\rho_k$ | 罚参数 |
| $\gamma$ | 罚参数的增长因子，$\rho_{k+1}=\gamma\rho_k$ |
| $\Phi_\rho(\vx)$ | 带罚目标函数 $f(\vx)+\rho\psi(\vx)$ |
| $\vx(\rho)$ | $\Phi_\rho$ 的一个极小点 |
| $\ell$ | 等式约束 $h_j(\vx)=0$ 的个数 |
| $\eta$ | 障碍参数 |
| $L_A(\vx,\vmu;\rho)$ | 增广拉格朗日函数 |

(sec-06-06-two-strategies)=
## 约束问题的两种策略

问题为

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)\\
\st\ & g_i(\vx)\le0,\quad i=1,\dots,m,\\
& h_j(\vx)=0,\quad j=1,\dots,\ell,
\end{aligned}
$$ (eq-06-06-nlp)

其可行集为 $\Omega$。算法的目标是局部极小点，或者退一步说，KKT 点（{doc}`06-03-kkt-conditions`）。有两种思路。

1. **可行方向法**（feasible-direction method）。从一个可行点出发，只沿既**可行**（对所有充分小的 $\alpha>0$ 有 $\vx+\alpha\vd\in\Omega$，{prf:ref}`def-06-02-feasible-direction`）又**下降**（$\nabla f(\vx)^\T\vd<0$）的方向 $\vd$ 移动。最速下降方向 $-\nabla f(\vx)$ 是自然的候选，但在边界上的点处，它往往指向 $\Omega$ 之外；方法必须找到一个在下降性与可行性之间折中的方向。每个迭代点都是可行的，因此方法可以随时停止并给出一个可用的解。
2. **变换方法。** 在 $f$ 上加一项惩罚不可行性（罚函数法）或阻止迭代点到达边界（障碍函数法）的项，把 {eq}`eq-06-06-nlp` 替换为一系列无约束问题，其解逼近 {eq}`eq-06-06-nlp` 的解，再用{doc}`06-05-unconstrained-algorithms`中的方法求解这些问题。

(sec-06-06-frank-wolfe)=
## Frank–Wolfe 方法

Frank 和 Wolfe（1956）针对*线性*约束问题提出了第一种思路。设

$$
\begin{aligned}
&\min_{\vx\in\Omega}f(\vx),\\
&\Omega=\{\vx\in\R^n:\mA\vx=\vb,\ \vx\ge\vzero\},
\end{aligned}
$$ (eq-06-06-fw-problem)

其中 $\Omega$ 非空有界（是一个多胞形），$f$ 是连续可微的凸函数。也允许出现不等式 $\mA\vx\le\vb$（添加松弛变量即可，{doc}`03-02-lp-geometry`）。

**由凸性得到可行方向。** 由于 $\Omega$ 是凸集，对任意 $\vx^{(k)}\in\Omega$ 和任意 $\vy\in\Omega$，两者之间的整条线段都是可行的：

$$
\begin{aligned}
\vx^{(k)}+\alpha(\vy-\vx^{(k)})&=(1-\alpha)\vx^{(k)}+\alpha\vy\\
&\in\Omega\qquad\text{对所有 }\alpha\in[0,1].
\end{aligned}
$$

因此每个满足 $\vy\in\Omega$ 的 $\vd=\vy-\vx^{(k)}$ 都是可行方向，并且步长限制在 $[0,1]$ 中。

**由线性规划得到下降方向。** 在这些方向中，我们要找方向导数 $\nabla f(\vx^{(k)})^\T(\vy-\vx^{(k)})$ 最负的那一个。这就是**线性化子问题**（linearized subproblem）

$$
\begin{aligned}
\min_{\vy}\ & \nabla f(\vx^{(k)})^\T\vy\\
\st\ & \mA\vy=\vb,\\
& \vy\ge\vzero,
\end{aligned}
$$ (eq-06-06-fw-lp)

它与 $\min_{\vy\in\Omega}\nabla f(\vx^{(k)})^\T(\vy-\vx^{(k)})$ 有相同的极小点，因为两个目标函数只相差常数 $\nabla f(\vx^{(k)})^\T\vx^{(k)}$。问题 {eq}`eq-06-06-fw-lp` 是一个 LP：它在 $\Omega$ 上最小化 $f$ 的一阶泰勒近似。由于 $\Omega$ 是非空多胞形，它在某个顶点处有最优解（{prf:ref}`thm-03-02-fundamental`），单纯形法可以求得该顶点（{prf:ref}`alg-03-03-simplex`）。记之为 $\vy^{(k)}$。最优值度量了可获得的下降量：

:::{prf:definition} Frank–Wolfe 间隙
:label: def-06-06-fw-gap

对 $\vx\in\Omega$，**Frank–Wolfe 间隙**（Frank–Wolfe gap）为

$$
\begin{aligned}
\Gamma(\vx)&=\max_{\vy\in\Omega}\nabla f(\vx)^\T(\vx-\vy)\\
&=-\min_{\vy\in\Omega}\nabla f(\vx)^\T(\vy-\vx)\ \ge0 .
\end{aligned}
$$

由于允许取 $\vy=\vx$，它是非负的。在第 $k$ 次迭代时，$\Gamma(\vx^{(k)})=-\nabla f(\vx^{(k)})^\T(\vy^{(k)}-\vx^{(k)})$。
:::

若 $\Gamma(\vx^{(k)})>0$，则 $\vd^{(k)}=\vy^{(k)}-\vx^{(k)}$ 是可行下降方向。若 $\Gamma(\vx^{(k)})=0$，则没有可行方向在一阶意义下是下降的，对于凸函数 $f$，当前点就是最优的。下面的定理把这一点精确化，并回答了为什么零间隙可以证明最优性。

:::{prf:theorem} 最优性检验与下界
:label: thm-06-06-fw-gap

设 $f$ 在凸集 $\Omega$ 上是可微凸函数，$f^*=\min_{\vx\in\Omega}f(\vx)$。对每个 $\vx\in\Omega$，

$$
f(\vx)-\Gamma(\vx)\le f^*\le f(\vx).
$$ (eq-06-06-fw-bound)

特别地，$\vx$ 是最优的当且仅当 $\Gamma(\vx)=0$，即对所有 $\vy\in\Omega$ 有 $\nabla f(\vx)^\T(\vy-\vx)\ge0$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

由凸性的一阶刻画（{prf:ref}`thm-02-02-first-order`），对所有 $\vy\in\Omega$ 有 $f(\vy)\ge f(\vx)+\nabla f(\vx)^\T(\vy-\vx)$。两边对 $\vy\in\Omega$ 取最小值，得 $f^*\ge f(\vx)+\min_{\vy\in\Omega}\nabla f(\vx)^\T(\vy-\vx)=f(\vx)-\Gamma(\vx)$。上界成立是因为 $\vx$ 可行。若 $\Gamma(\vx)=0$，两个界重合，$f(\vx)=f^*$。反之，若 $\Gamma(\vx)>0$，则达到最大值的方向 $\vy-\vx$ 是可行下降方向，因此 $f$ 沿该线段减小，$\vx$ 不是最优的。
:::

条件“对所有 $\vy\in\Omega$ 有 $\nabla f(\vx)^\T(\vy-\vx)\ge0$”是凸集上的一阶最优性条件（{prf:ref}`thm-02-03-optimality`，{prf:ref}`thm-06-02-first-order`）；对于多面体 $\Omega$，它等价于 KKT 条件（{prf:ref}`thm-06-03-kkt`），因为此时 LP {eq}`eq-06-06-fw-lp` 有最优解 $\vy=\vx$，而 LP 对偶提供了乘子。没有凸性时，$\Gamma(\vx)=0$ 仍然表明 $\vx$ 是 KKT 点。界 {eq}`eq-06-06-fw-bound` 是求解 LP 的副产品：每次迭代都给出一个**可证明的下界**（certified lower bound）$f(\vx^{(k)})-\Gamma(\vx^{(k)})$，这些下界中最好的那一个给出了当前值与最优值之差的上限。

**线搜索。** 新点是线段上最好的点，

$$
\begin{aligned}
\alpha_k&\in\argmin_{0\le\alpha\le1}f\bigl(\vx^{(k)}+\alpha\vd^{(k)}\bigr),\\
\vx^{(k+1)}&=\vx^{(k)}+\alpha_k\vd^{(k)},
\end{aligned}
$$ (eq-06-06-fw-line-search)

可以用{doc}`06-04-line-search`中的任何方法计算。对于黑塞矩阵为 $\mQ$ 的二次函数 $f$，$\varphi(\alpha)=f(\vx^{(k)})-\alpha\Gamma(\vx^{(k)})+\tfrac12\alpha^2\vd^{(k)\T}\mQ\vd^{(k)}$，精确步长为

$$
\alpha_k=\min\Bigl\{1,\ \frac{\Gamma(\vx^{(k)})}{\vd^{(k)\T}\mQ\vd^{(k)}}\Bigr\}
$$ (eq-06-06-fw-quadratic-step)

（若 $\vd^{(k)\T}\mQ\vd^{(k)}=0$，则取 $\alpha_k=1$）。上界 $1$ 是重要的：超过 $\alpha=1$ 后，点一般会离开 $\Omega$。

:::{prf:algorithm} Frank–Wolfe 方法（凸组合法）
:label: alg-06-06-frank-wolfe

**输入：** 连续可微的凸函数 $f$；非空多胞形 $\Omega$；可行点 $\vx^{(0)}\in\Omega$；容差 $\varepsilon>0$。
**输出：** 满足 $f(\vx)-f^*\le\varepsilon$ 的可行点 $\vx$。

1. **初始化。** 令 $k=0$。
2. **线性子问题。** 求解 LP $\min_{\vy\in\Omega}\nabla f(\vx^{(k)})^\T\vy$；设 $\vy^{(k)}$ 为一个最优顶点，$\vd^{(k)}=\vy^{(k)}-\vx^{(k)}$。
3. **终止检验。** 若 $\Gamma(\vx^{(k)})=-\nabla f(\vx^{(k)})^\T\vd^{(k)}\le\varepsilon$，则停止：由 {eq}`eq-06-06-fw-bound`，$f(\vx^{(k)})-f^*\le\varepsilon$。
4. **线搜索。** 求 $\alpha_k\in[0,1]$，使 $f(\vx^{(k)}+\alpha\vd^{(k)})$ 最小。
5. **更新。** 令 $\vx^{(k+1)}=\vx^{(k)}+\alpha_k\vd^{(k)}$，$k\leftarrow k+1$，转第 2 步。
:::

每个迭代点都是 $\vx^{(0)}$ 和顶点 $\vy^{(0)},\dots,\vy^{(k)}$ 的凸组合，因此该方法称为**凸组合法**（convex-combination method）；它也称为**线性近似法**或**条件梯度法**（conditional-gradient method）。可行的初始点本身就是一个 LP 的解（例如取 $\vx^{(0)}=\vy$，其中 $\vy$ 是单纯形法第一阶段返回的顶点，{doc}`03-04-two-phase-big-m`）。

:::{prf:example} 三角形上二次函数的 Frank–Wolfe 方法
:label: ex-06-06-frank-wolfe

从 $\vx^{(0)}=(0,0)$ 出发，用 Frank–Wolfe 方法求解

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)=-5x_1+x_1^2-8x_2+2x_2^2\\
\st\ & 3x_1+2x_2\le6,\\
& x_1,x_2\ge0,
\end{aligned}
$$
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

梯度为 $\nabla f(\vx)=(2x_1-5,\ 4x_2-8)$，黑塞矩阵为 $\mQ=\diag(2,4)$。可行集是以 $(0,0)$、$(2,0)$ 和 $(0,3)$ 为顶点的三角形，因此每个 LP 都可以通过比较三个顶点处的目标值来求解。无约束极小点 $(2.5,2)$ 违反 $3x_1+2x_2\le6$（$3\cdot2.5+2\cdot2=11.5$），因此该约束在最优点处是有效约束。

**第 0 次迭代。** $\nabla f(\vx^{(0)})=(-5,-8)$。三角形上的 LP $\min\,-5y_1-8y_2$ 在各顶点处的值为 $0$、$-10$、$-24$，因此 $\vy^{(0)}=(0,3)$，$\vd^{(0)}=(0,3)$，$\Gamma(\vx^{(0)})=24$。沿该线段，$\vx^{(0)}+\alpha\vd^{(0)}=(0,3\alpha)$，且

$$
\begin{aligned}
\varphi(\alpha)&=-24\alpha+18\alpha^2,\\
\varphi'(\alpha)&=-24+36\alpha=0\ \Rightarrow\ \alpha_0=\tfrac23 ,
\end{aligned}
$$

这与 {eq}`eq-06-06-fw-quadratic-step` 一致：$24/(\vd^\T\mQ\vd)=24/36$。因此 $\vx^{(1)}=(0,2)$，$f(\vx^{(1)})=-8$；下界为 $0-24=-24$。

**第 1 次迭代。** $\nabla f(\vx^{(1)})=(-5,0)$。LP $\min\,-5y_1$ 给出 $\vy^{(1)}=(2,0)$，$\vd^{(1)}=(2,-2)$，$\Gamma(\vx^{(1)})=10$。沿该线段，$\vx^{(1)}+\alpha\vd^{(1)}=(2\alpha,\ 2-2\alpha)$，且

$$
\begin{aligned}
\varphi(\alpha)&=-8-10\alpha+12\alpha^2,\\
\varphi'(\alpha)&=-10+24\alpha=0\ \Rightarrow\ \alpha_1=\tfrac{5}{12},
\end{aligned}
$$

因此 $\vx^{(2)}=(0,2)+\tfrac5{12}(2,-2)=(\tfrac56,\tfrac76)$，$f(\vx^{(2)})=-\tfrac{121}{12}=-10.0833$；下界改进为 $-8-10=-18$。*（这个直线函数有时被印成 $8+10\alpha-12\alpha^2$，导数为 $10-24\alpha$，这会使 $\alpha_1$ 成为极大点；正确的符号如上所示，步长 $\alpha_1=\tfrac5{12}$ 是极小点。）*

**第 2 次迭代。** $\nabla f(\vx^{(2)})=(-\tfrac{10}3,-\tfrac{10}3)$。各顶点处的值为 $0$、$-\tfrac{20}3$ 和 $-10$，因此 $\vy^{(2)}=(0,3)$，$\vd^{(2)}=(-\tfrac56,\tfrac{11}6)$，$\Gamma(\vx^{(2)})=\tfrac{10}3$，由 {eq}`eq-06-06-fw-quadratic-step` 得 $\alpha_2=\tfrac{10/3}{89/6}=\tfrac{20}{89}=0.2247$。因此 $\vx^{(3)}=(\tfrac{115}{178},\tfrac{281}{178})=(0.6461,\ 1.5787)$，$f(\vx^{(3)})=-10.4579$。

**第 3 次迭代。** $\nabla f(\vx^{(3)})=(-3.7079,\ -1.6854)$；此时 $(2,0)$ 是 LP 的最优解，$\Gamma=2.3596$，$\alpha_3=0.1731$，$\vx^{(4)}=(0.8804,\ 1.3055)$，$f(\vx^{(4)})=-10.6620$。

{numref}`tab-06-06-frank-wolfe`列出了后续迭代。最优解可由约束 $3x_1+2x_2\le6$ 有效时的 KKT 条件求得：由 $2x_1-5+3\lambda=0$，$4x_2-8+2\lambda=0$，$3x_1+2x_2=6$ 得 $\lambda=1\ge0$，$\vx^*=(1,\ 1.5)$，$f^*=-11.5$。

```{list-table} 在三角形上求 $\min\,-5x_1+x_1^2-8x_2+2x_2^2$ 的 Frank–Wolfe 迭代（4 位有效数字）
:name: tab-06-06-frank-wolfe
:header-rows: 1

* - $k$
  - $\vx^{(k)}$
  - $f(\vx^{(k)})$
  - $\vy^{(k)}$
  - $\Gamma(\vx^{(k)})$
  - $\alpha_k$
  - 最好的下界
* - 0
  - $(0,\ 0)$
  - $0$
  - $(0,3)$
  - $24$
  - $0.6667$
  - $-24$
* - 1
  - $(0,\ 2)$
  - $-8$
  - $(2,0)$
  - $10$
  - $0.4167$
  - $-18$
* - 2
  - $(0.8333,\ 1.167)$
  - $-10.08$
  - $(0,3)$
  - $3.333$
  - $0.2247$
  - $-13.42$
* - 3
  - $(0.6461,\ 1.579)$
  - $-10.46$
  - $(2,0)$
  - $2.360$
  - $0.1731$
  - $-12.82$
* - 4
  - $(0.8804,\ 1.305)$
  - $-10.66$
  - $(0,3)$
  - $1.856$
  - $0.1424$
  - $-12.52$
* - 5
  - $(0.7550,\ 1.547)$
  - $-10.79$
  - $(2,0)$
  - $1.541$
  - $0.1216$
  - $-12.33$
* - 6
  - $(0.9064,\ 1.359)$
  - $-10.89$
  - $(0,3)$
  - $1.322$
  - $0.1064$
  - $-12.21$
* - 7
  - $(0.8099,\ 1.533)$
  - $-10.96$
  - $(2,0)$
  - $1.160$
  - $0.09483$
  - $-12.12$
* - 8
  - $(0.9228,\ 1.388)$
  - $-11.01$
  - $(0,3)$
  - $1.036$
  - $0.08562$
  - $-12.05$
* - 9
  - $(0.8438,\ 1.526)$
  - $-11.06$
  - $(2,0)$
  - $0.9365$
  - $0.07812$
  - $-11.99$
```

从第 2 次迭代起，LP 的解在有效边的两个端点 $(0,3)$ 和 $(2,0)$ 之间交替，迭代点以不断缩小的步长呈锯齿形趋向 $\vx^*=(1,1.5)$（{numref}`fig-06-06-frank-wolfe`）。进展很慢：间隙在 $k=9$ 时降到 $1$ 以下，直到 $k=105$ 才降到 $0.1$ 以下，而 $f(\vx^{(k)})-f^*\le0.01$ 需要 $k=542$。但这些界始终是有效的：在 $k=9$ 时已知最优值位于 $[-11.99,\ -11.06]$ 中。

```{figure} figures/06-06-frank-wolfe.svg
:name: fig-06-06-frank-wolfe
:width: 560px
:class: dark-light
:alt: 可行三角形以及中心位于三角形外部的目标函数椭圆形等值线，Frank–Wolfe 迭代点的路径先沿 x2 轴上升，然后呈锯齿形趋向斜边上的最优解

{prf:ref}`ex-06-06-frank-wolfe` 中三角形上的 Frank–Wolfe 方法。$f$ 的等值线（中心 $(2.5,2)$ 位于 $\Omega$ 之外）、迭代点 $\vx^{(0)},\dots,\vx^{(9)}$ 以及边 $3x_1+2x_2=6$ 上的最优解 $\vx^*=(1,1.5)$。每一步都朝某个顶点 $\vy^{(k)}$（点线）移动，交替为 $(0,3)$ 和 $(2,0)$，因此路径呈锯齿形，步长不断缩小。
```
:::

**收敛性。** 在{prf:ref}`ex-06-06-frank-wolfe`的解答中看到的锯齿现象并不是该例的偶然现象：当最优解位于 $\Omega$ 的某个面上而不在顶点处时，每个 LP 的解都是顶点，因此步长方向永远不能径直指向 $\vx^*$。收敛速度是次线性的。

:::{prf:theorem} Frank–Wolfe 方法的收敛速度
:label: thm-06-06-fw-rate

设 $\Omega$ 是直径为 $D=\max_{\vx,\vy\in\Omega}\norm{\vx-\vy}$ 的紧凸集，$f$ 是梯度 Lipschitz 连续的凸函数，$\norm{\nabla f(\vx)-\nabla f(\vy)}\le L_f\norm{\vx-\vy}$。采用精确线搜索（或固定步长 $\alpha_k=2/(k+2)$）时，Frank–Wolfe 迭代点满足

$$
f(\vx^{(k)})-f^*\le\frac{2L_fD^2}{k+2},\qquad k\ge1 .
$$ (eq-06-06-fw-rate)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

令 $h_k=f(\vx^{(k)})-f^*$。梯度 Lipschitz 连续给出二次上界 $f(\vx+\alpha\vd)\le f(\vx)+\alpha\nabla f(\vx)^\T\vd+\tfrac12L_f\alpha^2\norm{\vd}^2$。取 $\vd=\vd^{(k)}$，由 {eq}`eq-06-06-fw-bound` 有 $\nabla f(\vx^{(k)})^\T\vd^{(k)}=-\Gamma(\vx^{(k)})\le-h_k$，且 $\norm{\vd^{(k)}}\le D$：

$$
h_{k+1}\le(1-\alpha)h_k+\tfrac12L_f D^2\alpha^2\qquad\text{对每个 }\alpha\in[0,1],
$$

这是因为精确线搜索至少与任何固定的 $\alpha\in[0,1]$ 一样好。对 $k=0$ 取 $\alpha=1$：$h_1\le\tfrac12L_fD^2\le\tfrac23L_fD^2$。若 $h_k\le2L_fD^2/(k+2)$，取 $\alpha=2/(k+2)$：

$$
\begin{aligned}
h_{k+1}&\le\frac{k}{k+2}\cdot\frac{2L_fD^2}{k+2}+\frac{2L_fD^2}{(k+2)^2}\\
&=\frac{2L_fD^2(k+1)}{(k+2)^2}\le\frac{2L_fD^2}{k+3},
\end{aligned}
$$

因为 $(k+1)(k+3)\le(k+2)^2$。
:::

在{prf:ref}`ex-06-06-frank-wolfe`中，$L_f=4$（$\mQ$ 的最大特征值），$D=\sqrt{13}$（从 $(2,0)$ 到 $(0,3)$ 的边），因此定理保证 $f(\vx^{(k)})-f^*\le104/(k+2)$；观察到的误差约为 $5/k$，远好于这个界，但同为 $1/k$ 量级。带**远离步**（away step，可以从凸组合中已有的顶点离开）的变体，在 $f$ 强凸时在多胞形上线性收敛。

**何时首选 Frank–Wolfe 方法？** 它的优势在于子问题是 LP，而且往往是非常特殊的 LP。在**弧费用为凸函数的网络流问题**中，例如道路网络上的交通分配（Sheffi, 1985），流量守恒约束上的 LP {eq}`eq-06-06-fw-lp` 分解为若干**最短路问题**，以边际弧费用 $\nabla f(\vx^{(k)})$ 作为长度，其解把全部需求加载到最短路上（“全有全无”分配）。于是每次迭代只需要一次最短路计算（例如{prf:ref}`alg-04-02-dijkstra`），完全不需要矩阵，这使该方法适用于有成千上万条路段的网络（{prf:ref}`ex-08-09-fw-network-flow`）。稀疏的迭代点（少数几个顶点的凸组合）在机器学习中也很有吸引力。

(sec-06-06-other-feasible-direction)=
### 其他可行方向法

Frank–Wolfe 方法按一阶意义下最好的极点来选择方向，这浪费了“当前点已经接近解”这一信息。若干可行方向法更巧妙地利用了约束（Bazaraa et al., 2006, Ch. 10；Luenberger & Ye, 2016, Ch. 12）。

- **Zoutendijk 方法**（Zoutendijk's method）在每次迭代中求解一个小规模 LP，以得到方向 $\vd$：在保持对有效约束可行以及满足诸如 $-1\le d_j\le1$ 的归一化条件下，使下降率最大（Zoutendijk, 1960）；这样方向就不再被迫指向顶点。
- **梯度投影法**（gradient projection method；Rosen, 1960）把 $-\nabla f(\vx)$ 投影到有效线性约束的子空间上，并沿投影方向移动；当投影为零时，由乘子估计值的符号决定放松哪个约束。
- **既约梯度法**（reduced-gradient method；Wolfe, 1963；后推广到非线性约束）像单纯形法（{doc}`03-03-simplex-method`）那样把变量划分为基变量和非基变量，通过等式约束消去基变量，并对非基变量执行最速下降步或拟牛顿步。
- **序列二次规划**（sequential quadratic programming, SQP；Wilson, 1963）用利用拉格朗日函数二阶信息的二次子问题代替线性子问题 {eq}`eq-06-06-fw-lp`；它是应用于 KKT 条件的牛顿法，超线性收敛（Nocedal & Wright, 2006, Ch. 18）。

它们都保留了 Frank–Wolfe 方法的原则，即在每次迭代中求解一个较简单的约束子问题，并且在没有有效约束时都化为{doc}`06-05-unconstrained-algorithms`中的方法。

(sec-06-06-penalty)=
## 罚函数法

第二种策略把约束移入目标函数。先考虑只有单个不等式约束的问题：在 $g(\vx)\le0$ 的约束下求 $\min f(\vx)$。其思想是*允许*违反约束，但要为之付出代价：当 $g(\vx)\le0$ 时不收费，当 $g(\vx)>0$ 时收取随 $g(\vx)$ 增长的费用。最简单的这种费用 $\rho\max\{0,g(\vx)\}$ 给出**精确罚函数**（exact penalty function）问题

$$
\min_{\vx}\ f(\vx)+\rho\max\{0,g(\vx)\}.
$$ (eq-06-06-exact-penalty)

之所以称为精确的，是因为若 $\vx^*$ 是约束问题的严格局部极小点，并以乘子 $\lambda^*$ 满足 KKT 条件（{prf:ref}`thm-06-03-kkt`），则对每个 $\rho>\lambda^*$，$\vx^*$ 也是 {eq}`eq-06-06-exact-penalty` 的局部极小点。其缺点是在 $g(\vx)=0$ 处有折点，而解通常恰好位于那里，因此不能应用{doc}`06-05-unconstrained-algorithms`中的光滑方法。经典的补救办法是把违反量平方。

:::{prf:definition} 二次罚函数
:label: def-06-06-quadratic-penalty

对于问题 {eq}`eq-06-06-nlp`，**二次罚函数**（quadratic penalty function）为

$$
\psi(\vx)=\sum_{i=1}^{m}\bigl[\max\{0,g_i(\vx)\}\bigr]^2+\sum_{j=1}^{\ell}h_j(\vx)^2 ,
$$ (eq-06-06-quadratic-penalty)

参数为 $\rho>0$ 的**带罚目标函数**（penalized objective）为 $\Phi_\rho(\vx)=f(\vx)+\rho\psi(\vx)$。
:::

罚函数 $\psi$ 连续、非负，并且恰在 $\Omega$ 上为零。若 $g_i$ 和 $h_j$ 连续可微，则 $\psi$ 也连续可微：函数 $t\mapsto\max\{0,t\}^2$ 的导数为 $2\max\{0,t\}$，它在 $t=0$ 处连续。更一般地，可以使用指数 $r\ge1$ 的 $\psi(\vx)=\sum_i\max\{0,g_i(\vx)\}^r+\sum_j\abs{h_j(\vx)}^r$：$r=1$ 给出精确（非光滑）罚函数，而每个 $r>1$ 都给出连续可微的 $\psi$。指数不必是偶数，因为最大值和绝对值已经使每一项非负。二次罚函数 $r=2$ 是通常的选择。它可以追溯到 Courant（1943）；把它系统地用作一系列无约束问题则归功于 Fiacco 和 McCormick（1968）。

:::{prf:example} 一维问题的二次罚函数
:label: ex-06-06-penalty-1d

用二次罚函数法求解问题：在 $x\ge2$ 的约束下求 $\min x$。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

把约束写成 $g(x)=2-x\le0$。罚函数为 $\psi(x)=\max\{0,2-x\}^2$，当 $x\ge2$ 时等于 $0$，当 $x<2$ 时等于 $(2-x)^2$，且

$$
\Phi_\rho(x)=x+\rho\max\{0,2-x\}^2 .
$$

在 $x\ge2$ 上 $\Phi_\rho(x)=x$ 递增，因此极小点必须满足 $x\le2$。当 $x<2$ 时，由 $\Phi_\rho'(x)=1-2\rho(2-x)=0$ 得

$$
x(\rho)=2-\frac{1}{2\rho},
$$

而 $\Phi_\rho''(x)=2\rho>0$ 确认这是极小点。*（导数有时被印成 $1-2\rho(2-x)^2$；那里不应有平方，只有正确的导数 $1-2\rho(2-x)$ 才给出 $x(\rho)$ 的公式。）*

```{list-table} 在 $x\ge2$ 的约束下求 $\min x$ 的二次罚函数法
:name: tab-06-06-penalty-1d
:header-rows: 1

* - $\rho$
  - $x(\rho)$
  - $\psi(x(\rho))$
  - $\rho\,\psi(x(\rho))$
  - $\Phi_\rho(x(\rho))$
  - $2\rho\max\{0,g(x(\rho))\}$
* - 1
  - 1.5
  - 0.25
  - 0.25
  - 1.75
  - 1
* - 10
  - 1.95
  - $2.5\times10^{-3}$
  - 0.025
  - 1.975
  - 1
* - 100
  - 1.995
  - $2.5\times10^{-5}$
  - $2.5\times10^{-3}$
  - 1.9975
  - 1
* - 1000
  - 1.9995
  - $2.5\times10^{-7}$
  - $2.5\times10^{-4}$
  - 1.99975
  - 1
```

当 $\rho\to\infty$ 时，$x(\rho)\to2=x^*$，且是**从可行集外部**逼近：每个 $x(\rho)$ 都稍有不可行，违反量为 $1/(2\rho)$。带罚目标函数值 $\Phi_\rho(x(\rho))=2-1/(4\rho)$ 递增趋于 $f^*=2$，并且是 $f^*$ 的下界。最后一列是常数，等于该问题的 KKT 乘子：由 $1-\lambda=0$ 得 $\lambda^*=1$。{numref}`fig-06-06-penalty-family`画出了这些带罚函数。

```{figure} figures/06-06-penalty-family.svg
:name: fig-06-06-penalty-family
:width: 520px
:class: dark-light
:alt: 若干条凸曲线 x 加 rho 乘以违反量的平方，rho 依次增大，在可行一侧都与直线 f(x) = x 重合，其极小点向右移动趋向 x = 2

{prf:ref}`ex-06-06-penalty-1d` 中 $\rho=0.5,1,2,5$ 时的带罚目标函数 $\Phi_\rho(x)=x+\rho\max\{0,2-x\}^2$。在可行一侧 $x\ge2$（阴影）上，它们都与 $f(x)=x$ 重合。极小点 $x(\rho)=2-1/(2\rho)=1$、$1.5$、$1.75$、$1.9$（圆点）位于可行集之外，随着 $\rho$ 增大逼近 $x^*=2$，同时曲线变得更陡。
```
:::

对于{prf:ref}`ex-06-06-penalty-1d`的问题，精确罚函数 {eq}`eq-06-06-exact-penalty` $x+\rho\max\{0,2-x\}$ 的表现不同：当 $\rho<1$ 时它无下界（$x\to-\infty$ 时 $(1-\rho)x+2\rho\to-\infty$），当 $\rho=1$ 时每个 $x\le2$ 都是极小点，而对每个 $\rho>1=\lambda^*$，其唯一的极小点恰为 $x^*=2$。一个有限的 $\rho$ 就足够了，代价是非光滑性。

:::{prf:algorithm} 二次罚函数法
:label: alg-06-06-penalty

**输入：** 问题 {eq}`eq-06-06-nlp`；初始参数 $\rho_0>0$；增长因子 $\gamma>1$（例如 $\gamma=10$）；初始点 $\vx^{(-1)}$；容差 $\varepsilon>0$。
**输出：** 近似解和乘子估计值。

1. **初始化。** 令 $k=0$。
2. **无约束最小化。** 从 $\vx^{(k-1)}$ 出发，用{doc}`06-05-unconstrained-algorithms`中的方法计算 $\Phi_{\rho_k}(\vx)=f(\vx)+\rho_k\psi(\vx)$ 的极小点 $\vx^{(k)}=\vx(\rho_k)$。
3. **终止检验。** 若 $\rho_k\psi(\vx^{(k)})\le\varepsilon$（或约束违反量低于容差），则停止。
4. **增大罚参数。** 令 $\rho_{k+1}=\gamma\rho_k$，$k\leftarrow k+1$，转第 2 步。
:::

每次最小化都从上一次的解出发（**热启动**，warm start）是至关重要的：$\vx(\rho_k)$ 接近 $\vx(\rho_{k+1})$，而在 $\rho$ 很大的问题上冷启动是很困难的。

**参数的选择。** 初始值 $\rho_0$ 应使罚项与目标函数大小相当，例如在初始点处取 $\rho_0\psi(\vx^{(-1)})\approx\abs{f(\vx^{(-1)})}$。大的增长因子（$\gamma=10$）只需较少的轮数，但每次热启动不那么准确；小的增长因子（$\gamma=2$）需要更多轮，但每轮的最小化更容易。由于有效约束的违反量约为 $\lambda_i^*/(2\rho)$（{ref}`sec-06-06-multipliers`），乘子为 $1$ 量级时，要使违反量达到 $10^{-6}$ 需要 $\rho\approx10^{6}$，而这大约正是下文所述的病态在双精度运算中开始损失有效数字的地方。每个子问题也不必精确求解：随 $k$ 减小的梯度容差就足以保证收敛。基本的收敛结果只需要连续性。

:::{prf:theorem} 罚函数法的收敛性
:label: thm-06-06-penalty-convergence

设 $f$、$g_i$ 和 $h_j$ 连续，{eq}`eq-06-06-nlp` 有最优解且最优值为 $f^*$，$0<\rho_k<\rho_{k+1}\to\infty$，且每个 $\vx^{(k)}$ 都是 $\Phi_{\rho_k}$ 的全局极小点。则

1. $\Phi_{\rho_k}(\vx^{(k)})\le\Phi_{\rho_{k+1}}(\vx^{(k+1)})\le f^*$，$\psi(\vx^{(k+1)})\le\psi(\vx^{(k)})$，且 $f(\vx^{(k)})\le f(\vx^{(k+1)})$；
2. $\{\vx^{(k)}\}$ 的每个极限点都是 {eq}`eq-06-06-nlp` 的最优解。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

(1) 记 $\vx_k=\vx^{(k)}$，$\psi_k=\psi(\vx_k)$，$f_k=f(\vx_k)$。对可行的最优解 $\vx^*$，$\Phi_{\rho_k}(\vx_k)\le\Phi_{\rho_k}(\vx^*)=f^*$。由于 $\rho_{k+1}>\rho_k$ 且 $\psi\ge0$，$\Phi_{\rho_k}(\vx_k)\le\Phi_{\rho_k}(\vx_{k+1})\le\Phi_{\rho_{k+1}}(\vx_{k+1})$。把最优性不等式 $f_k+\rho_k\psi_k\le f_{k+1}+\rho_k\psi_{k+1}$ 和 $f_{k+1}+\rho_{k+1}\psi_{k+1}\le f_k+\rho_{k+1}\psi_k$ 相加，得 $(\rho_{k+1}-\rho_k)(\psi_{k+1}-\psi_k)\le0$，因此 $\psi_{k+1}\le\psi_k$，再由第一个不等式得 $f_k\le f_{k+1}$。(2) 设 $\bar\vx$ 是子序列 $\vx_k$（$k\in\mathcal K$）的极限。由 $f_k+\rho_k\psi_k\le f^*$ 以及沿 $\mathcal K$ 有 $f_k\to f(\bar\vx)$ 可知，$\rho_k\psi_k$ 保持有界，因此 $\psi_k\to0$，由连续性 $\psi(\bar\vx)=0$：$\bar\vx$ 可行。此外 $f(\bar\vx)=\lim f_k\le f^*$，所以 $\bar\vx$ 是最优的。
:::

该定理假设取 $\Phi_{\rho_k}$ 的全局极小点；若取局部极小点，同样的论证在孤立局部解的附近也适用。第 (1) 部分说明函数值 $\Phi_{\rho_k}(\vx^{(k)})$ 是递增趋于 $f^*$ 的下界，正如{numref}`tab-06-06-penalty-1d`所示。

(sec-06-06-multipliers)=
### 乘子估计

在数据可微的情况下，在 $\Phi_\rho$ 的无约束极小点处，

$$
\begin{aligned}
&\nabla f(\vx(\rho))+\sum_{i=1}^{m}2\rho\max\{0,g_i(\vx(\rho))\}\nabla g_i(\vx(\rho))\\
&\quad+\sum_{j=1}^{\ell}2\rho\,h_j(\vx(\rho))\nabla h_j(\vx(\rho))=\vzero .
\end{aligned}
$$ (eq-06-06-penalty-stationarity)

与 KKT 平稳性条件 $\nabla f+\sum_i\lambda_i\nabla g_i+\sum_j\mu_j\nabla h_j=\vzero$（{prf:ref}`thm-06-03-kkt`）比较可知，量

$$
\begin{aligned}
\lambda_i(\rho)&=2\rho\max\{0,g_i(\vx(\rho))\}\ge0,\\
\mu_j(\rho)&=2\rho\,h_j(\vx(\rho))
\end{aligned}
$$ (eq-06-06-multiplier-estimates)

起着乘子的作用，并且在约束规范下，当 $\rho\to\infty$ 时它们收敛到 KKT 乘子 $\lambda_i^*$、$\mu_j^*$（Nocedal & Wright, 2006, Thm. 17.2）。一旦 $\vx(\rho)$ 接近 $\vx^*$，非有效约束就给出 $\lambda_i(\rho)=0$，这就是互补松弛性。在{prf:ref}`ex-06-06-penalty-1d`中，对每个 $\rho$ 估计值都恰为 $1$。这些估计还解释了 $\vx(\rho)$ *为什么*不可行：非零乘子需要大小为 $\lambda_i^*/(2\rho)$ 的违反量。

(sec-06-06-ill-conditioning)=
### 病态性

大 $\rho$ 的代价体现在数值上。在 $\vx(\rho)$ 附近，$\Phi_\rho$ 的黑塞矩阵近似为

$$
\begin{aligned}
\nabla^2\Phi_\rho&\approx\nabla^2_{\vx\vx}L(\vx,\vlambda,\vmu)\\
&\quad+2\rho\sum_{i\in I}\nabla g_i\nabla g_i^\T+2\rho\sum_{j}\nabla h_j\nabla h_j^\T ,
\end{aligned}
$$

即拉格朗日函数的黑塞矩阵加上在有效约束梯度方向上大小为 $\rho$ 的一项。一些特征值像 $\rho$ 一样增长，另一些保持有界，因此条件数随 $\rho$ 线性增长。在{prf:ref}`ex-06-06-penalty-1d`中，$\Phi_\rho$ 的曲率 $2\rho$ 无界增长；在{prf:ref}`ex-08-09-penalty-parabola`的二元问题中，$\nabla^2\Phi_\rho$ 的条件数从 $\rho=0.1$ 时的 $3.7$ 增长到 $\rho=10^4$ 时的约 $1.2\times10^4$。由{prf:ref}`thm-06-05-sd-rate`，最速下降法在这样的函数上慢得毫无希望，即使牛顿法也受舍入误差之害。这就是为什么要配合热启动逐步增大罚参数，也是后两节中的方法得以发展的原因。

(sec-06-06-barrier)=
## 障碍函数法

**障碍函数法**（barrier method，又称内点法）从可行集*内部*逼近解。对于内部 $\{\vx:g_i(\vx)<0\}$ 非空的不等式约束，它在内部最小化

$$
f(\vx)-\eta\sum_{i=1}^{m}\ln\bigl(-g_i(\vx)\bigr),\qquad \eta>0,
$$ (eq-06-06-log-barrier)

对数障碍函数在边界处趋于 $+\infty$，因此迭代点保持严格可行，而障碍参数（barrier parameter）$\eta$ 减小到零。对于在 $x\ge2$ 的约束下求 $\min x$ 的问题，障碍函数 $x-\eta\ln(x-2)$ 的极小点为 $x(\eta)=2+\eta$，它从可行一侧逼近 $x^*=2$；量 $\eta/(-g_i(\vx(\eta)))$（这里为 $\eta/\eta=1$）给出乘子的估计。障碍函数法与罚函数法有同样的病态问题，但其现代形式——把牛顿法应用于摄动 KKT 条件的**原始–对偶内点法**（primal–dual interior-point method）——是求解线性规划、二次规划和凸规划问题的最有效算法之一（Nocedal & Wright, 2006, Ch. 19；Fiacco & McCormick, 1968）。

(sec-06-06-augmented-lagrangian)=
## 增广拉格朗日法

二次罚函数法需要 $\rho\to\infty$，因为乘子信息是由违反量携带的。**增广拉格朗日函数**（augmented Lagrangian）显式地提供乘子。对于等式约束，它为

$$
L_A(\vx,\vmu;\rho)=f(\vx)+\sum_{j=1}^{\ell}\mu_jh_j(\vx)+\frac{\rho}{2}\sum_{j=1}^{\ell}h_j(\vx)^2 .
$$ (eq-06-06-augmented-lagrangian)

每次外层迭代在固定的 $\vmu$ 和 $\rho$ 下关于 $\vx$ 最小化 $L_A$，然后按 $\mu_j\leftarrow\mu_j+\rho\,h_j(\vx)$ 修正乘子，这是 {eq}`eq-06-06-multiplier-estimates` 的类比。若 $\vmu$ 是正确的乘子向量，并且二阶充分条件（{prf:ref}`thm-06-03-second-order`）在 $\vx^*$ 处成立，则对每个充分大但*有限*的 $\rho$，$\vx^*$ 都是 $L_A$ 的局部极小点（Nocedal & Wright, 2006, Thm. 17.5），从而避免了病态。对于在 $h(x)=2-x=0$ 的约束下求 $\min x$ 的问题，取 $\rho=1$，$\mu_0=0$：最小化 $x+\mu(2-x)+\tfrac12(2-x)^2$ 得 $x=2-(1-\mu)/\rho=1$；修正给出 $\mu_1=0+1\cdot(2-1)=1$，下一次最小化恰好给出 $x=2$，而 $\mu=1$ 正是正确的乘子。不等式约束可以通过松弛变量或修改后的修正公式来处理。该方法（也称为乘子法；Hestenes, 1969；Powell, 1969）是许多大规模 NLP 程序的基础（Nocedal & Wright, 2006, Ch. 17；Bazaraa et al., 2006, Sec. 9.3）。

## 本章回顾

### 小结

- 约束问题或者用可行方向法求解，或者通过变换为无约束问题来求解。
- Frank–Wolfe 方法（{prf:ref}`alg-06-06-frank-wolfe`）求解 LP {eq}`eq-06-06-fw-lp` 得到顶点 $\vy^{(k)}$，以 $\alpha\in[0,1]$ 沿 $\vy^{(k)}-\vx^{(k)}$ 移动，并在间隙 $\Gamma$（{prf:ref}`def-06-06-fw-gap`）很小时停止；间隙给出下界 {eq}`eq-06-06-fw-bound`（{prf:ref}`thm-06-06-fw-gap`）。
- 其收敛速度为 $O(1/k)$（{prf:ref}`thm-06-06-fw-rate`）；在{prf:ref}`ex-06-06-frank-wolfe`中，迭代点在两个顶点之间呈锯齿形趋向 $\vx^*=(1,1.5)$，$f^*=-11.5$（{numref}`fig-06-06-frank-wolfe`）。在网络流问题中，其 LP 是最短路问题。
- 二次罚函数 $\psi$（{prf:ref}`def-06-06-quadratic-penalty`）和方法{prf:ref}`alg-06-06-penalty`产生从外部逼近解的极小点 $\vx(\rho)$（{prf:ref}`thm-06-06-penalty-convergence`，{numref}`fig-06-06-penalty-family`）；$2\rho\max\{0,g_i\}$ 和 $2\rho h_j$ 给出 KKT 乘子的估计 {eq}`eq-06-06-multiplier-estimates`。
- 大的 $\rho$ 使 $\Phi_\rho$ 病态；障碍函数法 {eq}`eq-06-06-log-barrier` 从内部逼近，而增广拉格朗日函数 {eq}`eq-06-06-augmented-lagrangian` 用有限的 $\rho$ 即可达到解。

### 常见错误

- **把线段的方向写反。** Frank–Wolfe 方法的更新为 $\vx^{(k)}+\alpha(\vy^{(k)}-\vx^{(k)})=(1-\alpha)\vx^{(k)}+\alpha\vy^{(k)}$；权重 $\alpha$ 属于新的顶点。
- **让步长超过 1。** 在 Frank–Wolfe 方法中，若沿直线的无约束极小点 $\alpha>1$，必须截断为 $1$；超过 $\vy^{(k)}$ 后的点一般是不可行的。
- **使用错误点处的间隙。** 终止检验使用*当前*点处的 $\nabla f(\vx^{(k)})$；用上一次迭代的梯度得到的检验毫无意义。
- **期望 Frank–Wolfe 方法恰好停止。** 除非最优解是顶点，否则间隙只是趋于零；应在 $\Gamma(\vx^{(k)})\le\varepsilon$ 时停止，并报告下界 $f(\vx^{(k)})-\Gamma(\vx^{(k)})$。
- **把罚函数法的解当作可行解。** $\vx(\rho)$ 对有效约束的违反量约为 $\lambda_i^*/(2\rho)$；需要最后进行可行性检查或投影。
- **一开始就取极大的 $\rho$。** 这会使第一个子问题病态；应逐步增大 $\rho$ 并使用热启动。
- **把进展缓慢误认为发散。** 一个罚函数序列在几轮之后 $\rho_k\psi(\vx^{(k)})$ 仍高于容差，并不意味着失败；它只是尚未满足终止检验（{prf:ref}`ex-08-09-penalty-parabola`）。


## 参考文献

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 9: exterior
  penalty, exact penalty, augmented Lagrangian and barrier methods;
  Chapter 10: methods of feasible directions.)
- Courant, R. (1943). Variational methods for the solution of problems of
  equilibrium and vibrations. *Bulletin of the American Mathematical Society,
  49*(1), 1–23. <https://doi.org/10.1090/S0002-9904-1943-07818-4>
- Fiacco, A. V., & McCormick, G. P. (1968). *Nonlinear programming: Sequential
  unconstrained minimization techniques*. Wiley.
- Frank, M., & Wolfe, P. (1956). An algorithm for quadratic programming.
  *Naval Research Logistics Quarterly, 3*(1–2), 95–110.
  <https://doi.org/10.1002/nav.3800030109>
- Hestenes, M. R. (1969). Multiplier and gradient methods. *Journal of
  Optimization Theory and Applications, 4*(5), 303–320.
  <https://doi.org/10.1007/BF00927673>
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming* (4th
  ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3> (Chapter 12:
  primal methods, feasible directions; Chapter 13: penalty and barrier
  methods.)
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Chapter 17: quadratic
  penalty, exact penalty and augmented Lagrangian methods; Chapter 18:
  sequential quadratic programming; Chapter 19: interior-point methods.)
- Powell, M. J. D. (1969). A method for nonlinear constraints in minimization
  problems. In R. Fletcher (Ed.), *Optimization* (pp. 283–298). Academic
  Press.
- Rosen, J. B. (1960). The gradient projection method for nonlinear
  programming. Part I. Linear constraints. *Journal of the Society for
  Industrial and Applied Mathematics, 8*(1), 181–217.
- Sheffi, Y. (1985). *Urban transportation networks: Equilibrium analysis with
  mathematical programming methods*. Prentice-Hall. (Chapters 4–5: the
  convex-combination (Frank–Wolfe) algorithm and its application to traffic
  assignment.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 21: constrained nonlinear algorithms.)
- Wilson, R. B. (1963). *A simplicial algorithm for concave programming*
  [Doctoral dissertation, Harvard University].
- Wolfe, P. (1963). Methods of nonlinear programming. In R. L. Graves &
  P. Wolfe (Eds.), *Recent advances in mathematical programming*
  (pp. 67–86). McGraw-Hill.
- Zoutendijk, G. (1960). *Methods of feasible directions: A study in linear
  and nonlinear programming*. Elsevier.
