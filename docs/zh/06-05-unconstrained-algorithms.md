# 无约束优化算法

## 概述

本章把{doc}`06-04-line-search`中的一维搜索发展为最小化 $n$ 元光滑函数的方法，即求解 $\min_{\vx\in\R^n}f(\vx)$。这些方法都是**下降法**（descent method）：在当前点选择一个使 $f$ 减小的方向，再沿该方向选择步长。**最速下降法**（steepest descent method）使用负梯度方向；它简单且全局收敛，但其相邻步相互正交，路径在狭长的山谷中呈锯齿形，而且其线性收敛速度随黑塞矩阵条件数的增大而变差。**牛顿法**（Newton's method）利用黑塞矩阵构造二次模型；它一步即可求得凸二次函数的极小点，并在黑塞矩阵正定的极小点附近二次收敛，但需要二阶导数，并且在黑塞矩阵不正定时需要保护措施。介于两者之间的是**拟牛顿法**（quasi-Newton method，BFGS，以 Broyden、Fletcher、Goldfarb 和 Shanno 命名）和**共轭梯度法**（conjugate-gradient method），它们只使用梯度，却仍能超线性收敛，或在二次函数上有限步终止。所有方法都在同一个病态二次函数 $f(\vx)=(x_1-2)^2+10(x_2-2)^2$ 上运行。本章依赖于{doc}`06-02-optimality-conditions`中的最优性条件、{doc}`01-01-linear-algebra`中的特征值与有定性，以及{doc}`01-02-matrix-calculus`中的梯度、黑塞矩阵和泰勒展开。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 叙述一般的下降算法，并验证 $\mD\succ\vzero$ 时方向 $\vd=-\mD\nabla f(\vx)$ 是下降方向；
- 在二次函数上用闭式步长执行带精确线搜索的最速下降法，并解释锯齿形路径；
- 用黑塞矩阵的条件数界定最速下降法的收敛速度；
- 推导 $\R^n$ 中的牛顿法，证明它一步即可求得凸二次函数的极小点，并叙述其局部二次收敛性；
- 在黑塞矩阵不正定时修正牛顿方向；
- 描述 BFGS 方法和共轭梯度法背后的思想，并从每次迭代的计算量、存储量和收敛速度方面比较各种方法。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\vg^{(k)}$ | 迭代点处的梯度，$\vg^{(k)}=\nabla f(\vx^{(k)})$ |
| $\mH(\vx)$ | 黑塞矩阵 $\nabla^2f(\vx)$；$\mH^{(k)}=\mH(\vx^{(k)})$ |
| $\mD_k$ | 对梯度进行尺度变换的正定矩阵，$\vd^{(k)}=-\mD_k\vg^{(k)}$ |
| $\mQ$, $\vb$ | 二次函数 $f(\vx)=\tfrac12\vx^\T\mQ\vx-\vb^\T\vx$（$\mQ\succ\vzero$）的黑塞矩阵和线性项 |
| $\kappa$ | 条件数 $\lambda_{\max}(\mQ)/\lambda_{\min}(\mQ)$ |
| $\vs^{(k)}$, $\vy^{(k)}$ | 步长向量与梯度变化量，$\vs^{(k)}=\vx^{(k+1)}-\vx^{(k)}$，$\vy^{(k)}=\vg^{(k+1)}-\vg^{(k)}$（拟牛顿法） |
| $\theta_k$ | $\vd^{(k)}$ 与 $-\nabla f(\vx^{(k)})$ 之间的夹角 |
| $\gamma_k$ | BFGS 标量 $1/\vy^{(k)\T}\vs^{(k)}$ |
| $\beta_k$ | 共轭梯度系数 |
| $\Delta_k$ | 信赖域半径 |
| $\bar{\mH}_k$, $\nu_k$ | 黑塞矩阵的正定修正及构造它所用的平移量，$\bar{\mH}_k=\nabla^2f(\vx^{(k)})+\nu_k\mI$ |

(sec-06-05-descent-methods)=
## 下降法

对于 $f$ 连续可微的无约束问题 $\min_{\vx\in\R^n}f(\vx)$，本章的每种方法都是下面这一框架的特例（{doc}`06-04-line-search`，{eq}`eq-06-04-update`）。

:::{prf:algorithm} 下降法
:label: alg-06-05-descent

**输入：** 初始点 $\vx^{(0)}$；容差 $\varepsilon>0$。
**输出：** 近似驻点。

1. **初始化。** 令 $k=0$。
2. **终止检验。** 若 $\norm{\nabla f(\vx^{(k)})}\le\varepsilon$，则停止。
3. **方向。** 选择下降方向 $\vd^{(k)}$，$\nabla f(\vx^{(k)})^\T\vd^{(k)}<0$。
4. **步长。** 对 $\varphi(\alpha)=f(\vx^{(k)}+\alpha\vd^{(k)})$ 进行线搜索（精确线搜索，或满足 Wolfe 条件），选择 $\alpha_k>0$。
5. **更新。** 令 $\vx^{(k+1)}=\vx^{(k)}+\alpha_k\vd^{(k)}$，$k\leftarrow k+1$，转第 2 步。
:::

用正定矩阵对负梯度进行尺度变换，可以得到一大族方向。

:::{prf:proposition} 尺度化梯度方向是下降方向
:label: prop-06-05-scaled-gradient

若 $\nabla f(\vx)\ne\vzero$ 且 $\mD\succ\vzero$，则 $\vd=-\mD\nabla f(\vx)$ 是 $\vx$ 处的下降方向。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

由于 $\mD$ 正定且 $\nabla f(\vx)\ne\vzero$，有 $\nabla f(\vx)^\T\vd=-\nabla f(\vx)^\T\mD\nabla f(\vx)<0$。
:::

取 $\mD_k=\mI$ 得到最速下降法，取 $\mD_k=[\nabla^2f(\vx^{(k)})]^{-1}$（当它正定时）得到牛顿法，而由梯度差构造的矩阵 $\mD_k$ 给出拟牛顿法。

**全局收敛性。** 仅有下降性并不能保证收敛：方向可能变得几乎与梯度正交，或者步长可能过短。两个温和的要求可以排除这两种情况。设 $\theta_k$ 为 $\vd^{(k)}$ 与 $-\nabla f(\vx^{(k)})$ 之间的夹角。若 $f$ 有下界，其梯度 Lipschitz 连续，步长满足 Wolfe 条件（{eq}`eq-06-04-armijo`，{eq}`eq-06-04-curvature`），且对所有 $k$ 有 $\cos\theta_k\ge\delta>0$，则 $\nabla f(\vx^{(k)})\to\vzero$（Zoutendijk 定理；Zoutendijk, 1970；Nocedal & Wright, 2006, Thm. 3.2）。对于 $\vd^{(k)}=-\mD_k\vg^{(k)}$，当矩阵 $\mD_k$ 的条件数有界时，夹角条件成立。因此，若迭代点保持有界，则每个极限点都是驻点。与任何基于一阶信息的方法一样，这只保证得到驻点，而不保证得到极小点；在实践中，下降性使得收敛到鞍点或极大点的情况不大可能发生。

(sec-06-05-steepest-descent)=
## 最速下降法

在所有单位长度的方向中，哪一个使 $f$ 减小得最快？$f$ 在 $\vx$ 处沿 $\vd$ 的变化率是方向导数 $\nabla f(\vx)^\T\vd$，由柯西–施瓦茨不等式，

$$
\begin{aligned}
\nabla f(\vx)^\T\vd&\ge-\norm{\nabla f(\vx)}\,\norm{\vd}\\
&=-\norm{\nabla f(\vx)}\qquad\text{对 }\norm{\vd}=1,
\end{aligned}
$$

当 $\vd=-\nabla f(\vx)/\norm{\nabla f(\vx)}$ 时等号成立。因此负梯度方向是**最速下降方向**（steepest-descent direction）；从几何上看，它垂直于 $f$ 过 $\vx$ 的水平集（{prf:ref}`thm-01-02-gradient-level-set`）。

:::{prf:algorithm} 最速下降法（柯西法；Cauchy, 1847）
:label: alg-06-05-steepest-descent

**输入：** 初始点 $\vx^{(0)}$；容差 $\varepsilon>0$。
**输出：** 近似驻点。

1. **初始化。** 令 $k=0$。
2. **终止检验。** 计算 $\vg^{(k)}=\nabla f(\vx^{(k)})$。若 $\norm{\vg^{(k)}}\le\varepsilon$，则停止。
3. **方向。** 令 $\vd^{(k)}=-\vg^{(k)}$。
4. **线搜索。** 求在 $\alpha\ge0$ 上最小化 $f(\vx^{(k)}+\alpha\vd^{(k)})$ 的 $\alpha_k$（或满足 Wolfe 条件的 $\alpha_k$）。
5. **更新。** 令 $\vx^{(k+1)}=\vx^{(k)}+\alpha_k\vd^{(k)}$，$k\leftarrow k+1$，转第 2 步。
:::

采用精确线搜索时，该方法在{ref}`sec-06-05-descent-methods`的意义下全局收敛，而且这里的论证足够简短，可以完整给出。

:::{prf:proposition} 最速下降法的极限点是驻点
:label: prop-06-05-sd-convergence

设 $f$ 连续可微，带精确线搜索的最速下降法产生的序列 $\{\vx^{(k)}\}$ 有极限点 $\bar\vx$。则 $\nabla f(\bar\vx)=\vzero$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

假设 $\nabla f(\bar\vx)\ne\vzero$。由于 $-\nabla f(\bar\vx)$ 是下降方向，存在 $\bar\alpha>0$，使得对某个 $\delta>0$ 有 $f(\bar\vx-\bar\alpha\nabla f(\bar\vx))=f(\bar\vx)-2\delta$。由 $f$ 和 $\nabla f$ 的连续性，每个充分接近 $\bar\vx$ 的 $\vx$ 都满足 $f(\vx-\bar\alpha\nabla f(\vx))\le f(\bar\vx)-\delta$。设沿子序列 $\mathcal K$ 有 $\vx^{(k)}\to\bar\vx$。对充分大的 $k\in\mathcal K$，精确线搜索给出 $f(\vx^{(k+1)})\le f(\vx^{(k)}-\bar\alpha\vg^{(k)})\le f(\bar\vx)-\delta$。但函数值 $f(\vx^{(k)})$ 单调递减，并沿 $\mathcal K$ 收敛到 $f(\bar\vx)$，因此对所有 $j$ 有 $f(\vx^{(j)})\ge f(\bar\vx)$，矛盾。
:::

特别地，若水平集 $\{\vx:f(\vx)\le f(\vx^{(0)})\}$ 有界，则序列有极限点，并且它们全是驻点。该命题没有涉及速度；速度是本节其余部分的主题。

**二次函数的精确步长。** 对于严格凸二次函数

$$
\begin{aligned}
f(\vx)&=\tfrac12\vx^\T\mQ\vx-\vb^\T\vx+c, \qquad \mQ\succ\vzero,\\
\nabla f(\vx)&=\mQ\vx-\vb,
\end{aligned}
$$ (eq-06-05-quadratic)

沿 $\vd=-\vg$ 的直线函数是抛物线 $\varphi(\alpha)=f(\vx)-\alpha\vg^\T\vg+\tfrac12\alpha^2\vg^\T\mQ\vg$。令 $\varphi'(\alpha)=-\vg^\T\vg+\alpha\vg^\T\mQ\vg=0$，得精确步长

$$
\alpha_k=\frac{\vg^{(k)\T}\vg^{(k)}}{\vg^{(k)\T}\mQ\vg^{(k)}} .
$$ (eq-06-05-sd-step)

由瑞利商的界（{prf:ref}`prop-01-01-rayleigh`），$1/\lambda_{\max}(\mQ)\le\alpha_k\le1/\lambda_{\min}(\mQ)$。

**正交的步。** 采用精确线搜索时，$\varphi'(\alpha_k)=0$，即 $\nabla f(\vx^{(k+1)})^\T\vd^{(k)}=0$。对于最速下降法，$\vd^{(k)}=-\vg^{(k)}$，$\vd^{(k+1)}=-\vg^{(k+1)}$，因此

$$
\vd^{(k+1)\T}\vd^{(k)}=\vg^{(k+1)\T}\vg^{(k)}=0 :
$$ (eq-06-05-orthogonal)

**相邻两步相互正交**，无论 $f$ 是什么函数。新方向垂直于旧方向，而旧方向与新点处的水平集相切。在狭长的山谷中，这迫使路径以短小的正交步反复横穿山谷，而不是沿着谷底前进：这就是{numref}`fig-06-05-zigzag`中的**锯齿现象**（zigzagging）。

:::{prf:example} 病态二次函数上的最速下降法
:label: ex-06-05-steepest-descent

从 $\vx^{(0)}=(-4,-3)$ 出发，用带精确线搜索的最速下降法最小化 $f(\vx)=(x_1-2)^2+10(x_2-2)^2$。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

极小点为 $\vx^*=(2,2)$，$f^*=0$。写成 {eq}`eq-06-05-quadratic` 的形式，$\mQ=\diag(2,20)$，$\vb=(4,40)$，$\nabla f(\vx)=(2x_1-4,\ 20x_2-40)$。记 $u_k=4-2x_1^{(k)}$ 和 $v_k=40-20x_2^{(k)}$ 为 $\vd^{(k)}=-\vg^{(k)}$ 的分量。直线函数为

$$
\varphi(\alpha)=\bigl(x_1^{(k)}+\alpha u_k-2\bigr)^2+10\bigl(x_2^{(k)}+\alpha v_k-2\bigr)^2,
$$

由 $\varphi'(\alpha)=0$ 得（这就是 {eq}`eq-06-05-sd-step`）

$$
\alpha_k=\frac{u_k^2+v_k^2}{2u_k^2+20v_k^2}\in\Bigl[\frac1{20},\frac12\Bigr].
$$

在 $\vx^{(0)}$ 处：$\vg^{(0)}=(-12,-100)$，$\vd^{(0)}=(12,100)$，$\alpha_0=(144+10\,000)/(288+200\,000)=0.05065$，$\vx^{(1)}=(-4,-3)+0.05065\,(12,100)=(-3.3922,\ 2.0647)$，$f(\vx^{(1)})=29.12$，而 $f(\vx^{(0)})=286$。{numref}`tab-06-05-steepest-descent`列出后续迭代。

```{list-table} 从 $(-4,-3)$ 出发对 $f(\vx)=(x_1-2)^2+10(x_2-2)^2$ 用最速下降法
:name: tab-06-05-steepest-descent
:header-rows: 1

* - $k$
  - $\vx^{(k)}$
  - $\vd^{(k)}=-\nabla f(\vx^{(k)})$
  - $\alpha_k$
  - $f(\vx^{(k)})$
* - 0
  - $(-4.0000,\ -3.0000)$
  - $(12.0000,\ 100.0000)$
  - 0.05065
  - 286.0
* - 1
  - $(-3.3922,\ 2.0647)$
  - $(10.7845,\ -1.2941)$
  - 0.4434
  - 29.12
* - 2
  - $(1.3891,\ 1.4909)$
  - $(1.2217,\ 10.1811)$
  - 0.05065
  - 2.965
* - 3
  - $(1.4510,\ 2.0066)$
  - $(1.0980,\ -0.1318)$
  - 0.4434
  - 0.3018
* - 4
  - $(1.9378,\ 1.9482)$
  - $(0.1244,\ 1.0366)$
  - 0.05065
  - 0.03073
* - 5
  - $(1.9441,\ 2.0007)$
  - $(0.1118,\ -0.0134)$
  - 0.4434
  - $3.129\times10^{-3}$
* - 6
  - $(1.9937,\ 1.9947)$
  - $(0.0127,\ 0.1055)$
  - 0.05065
  - $3.185\times10^{-4}$
* - 7
  - $(1.9943,\ 2.0001)$
  - $(0.0114,\ -0.0014)$
  - 0.4434
  - $3.243\times10^{-5}$
* - 8
  - $(1.9994,\ 1.9995)$
  - $(0.0013,\ 0.0107)$
  - 0.05065
  - $3.302\times10^{-6}$
```

有三个特点十分突出。

- **锯齿现象。** 方向在“主要沿 $x_2$”（$k$ 为偶数，短步长 $\alpha_k=0.05065$）和“主要沿 $x_1$”（$k$ 为奇数，长步长 $\alpha_k=0.4434$）之间交替，并且相邻方向正交，例如 $(12,100)^\T(10.7845,-1.2941)=129.4-129.4=0$（{numref}`fig-06-05-zigzag`）。
- **线性收敛。** 每一步都把 $f(\vx^{(k)})-f^*$ 乘以同一个因子 $0.1018$（$29.12/286=2.965/29.12=\dots$）：每次迭代得到一位十进制有效数字。该比值是常数，是因为在二元情形下梯度在两个相互垂直的方向之间交替，而缩减因子 $1-(\vg^\T\vg)^2/\bigl((\vg^\T\mQ\vg)(\vg^\T\mQ^{-1}\vg)\bigr)$（在下面{prf:ref}`thm-06-05-sd-rate`的证明中导出）对向量 $\vg$ 及其垂直向量取相同的值。
- **代价。** 尽管极小点离初始点只有 $7.8$ 的距离，并且每一步沿其直线都是最优的，检验 $\norm{\nabla f}\le10^{-6}$ 直到 $k=17$ 才得到满足。
:::

```{figure} figures/06-05-zigzag.svg
:name: fig-06-05-zigzag
:width: 680px
:class: dark-light
:alt: 一个狭长二次函数的两幅等值线图；左图显示最速下降法的锯齿形路径、牛顿法的单步以及共轭梯度法到达极小点的两步路径；右图显示从最坏初始点出发的缓慢锯齿形路径

$f(\vx)=(x_1-2)^2+10(x_2-2)^2$ 的等值线，极小点 $\vx^*=(2,2)$。(a) 从 $\vx^{(0)}=(-4,-3)$ 出发：带精确线搜索的最速下降法（蓝色，正交的锯齿形步，{prf:ref}`ex-06-05-steepest-descent`）、牛顿法（橙色，一步，{prf:ref}`ex-06-05-newton-quadratic`）以及共轭梯度法 / BFGS 的路径（绿色，两步；其第一步就是最速下降法的第一步）。(b) 从最坏初始点 $(-8,1)$ 出发的最速下降法，$x_2$ 轴经过拉伸：每一步只把 $f$ 缩小为 {prf:ref}`thm-06-05-sd-rate` 中的因子 $((\kappa-1)/(\kappa+1))^2=0.6694$ 倍。
```

**一般情况下最速下降法有多快？** 对于二次函数，答案只取决于 $\mQ$ 的特征值。

:::{prf:theorem} 二次函数上最速下降法的收敛速度
:label: thm-06-05-sd-rate

设 $f$ 是 {eq}`eq-06-05-quadratic` 中的二次函数，$\mQ\succ\vzero$，极小点 $\vx^*=\mQ^{-1}\vb$，条件数 $\kappa=\lambda_{\max}(\mQ)/\lambda_{\min}(\mQ)$。带精确线搜索的最速下降法满足

$$
f(\vx^{(k+1)})-f(\vx^*)\le\Bigl(\frac{\kappa-1}{\kappa+1}\Bigr)^2
\bigl(f(\vx^{(k)})-f(\vx^*)\bigr).
$$ (eq-06-05-sd-rate)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

令 $E(\vx)=f(\vx)-f(\vx^*)=\tfrac12(\vx-\vx^*)^\T\mQ(\vx-\vx^*)$，$\vg=\mQ(\vx-\vx^*)$，于是 $E(\vx)=\tfrac12\vg^\T\mQ^{-1}\vg$。采用步长 {eq}`eq-06-05-sd-step` 时，$\varphi(\alpha_k)=f(\vx)-\tfrac12(\vg^\T\vg)^2/(\vg^\T\mQ\vg)$，因此

$$
\frac{E(\vx^{(k+1)})}{E(\vx^{(k)})}
=1-\frac{(\vg^\T\vg)^2}{(\vg^\T\mQ\vg)(\vg^\T\mQ^{-1}\vg)} .
$$

Kantorovich 不等式（Luenberger & Ye, 2016, Ch. 8）指出，对 $\mQ\succ\vzero$ 和每个 $\vg\ne\vzero$ 有

$$
\begin{aligned}
\frac{(\vg^\T\vg)^2}{(\vg^\T\mQ\vg)(\vg^\T\mQ^{-1}\vg)}
&\ge\frac{4\lambda_{\min}(\mQ)\,\lambda_{\max}(\mQ)}{\bigl(\lambda_{\min}(\mQ)+\lambda_{\max}(\mQ)\bigr)^2}\\
&=\frac{4\kappa}{(\kappa+1)^2} .
\end{aligned}
$$

代入得 $E(\vx^{(k+1)})/E(\vx^{(k)})\le\bigl((\kappa-1)/(\kappa+1)\bigr)^2$，这就是 {eq}`eq-06-05-sd-rate`。
:::

这个界是可以达到的。对于 $\mQ=\diag(2,20)$，$\kappa=10$，因子为 $(9/11)^2=0.6694$：从初始点 $(-8,1)$ 出发，其误差 $(-10,-1)$ 与 $(1/\lambda_1(\mQ),1/\lambda_2(\mQ))=(\tfrac12,\tfrac1{20})$ 成比例，步长始终为 $\alpha_k=2/(\lambda_1(\mQ)+\lambda_2(\mQ))=1/11$，$f$ 每步只缩小为 $0.6694$ 倍：$110$ → $73.64$ → $49.29$ → $33.00$ → $\dots$（{numref}`fig-06-05-zigzag`b）。{prf:ref}`ex-06-05-steepest-descent`的初始点 $(-4,-3)$ 恰好有利得多（因子 $0.1018$）。当 $\kappa=1$（等值线为圆）时，最速下降法一步即可到达 $\vx^*$；当 $\kappa=1000$ 时，最坏情形的因子为 $0.996$，$f(\vx^{(k)})-f^*$ 每增加一位十进制有效数字大约需要 $575$ 次迭代。

**固定步长。** 许多应用（包括机器学习）使用常数步长的最速下降法 $\vx^{(k+1)}=\vx^{(k)}-\alpha\nabla f(\vx^{(k)})$，不做线搜索。此时在二次函数 {eq}`eq-06-05-quadratic` 上，误差按线性规律演化，$\vx^{(k+1)}-\vx^*=(\mI-\alpha\mQ)(\vx^{(k)}-\vx^*)$，而 $\mI-\alpha\mQ$ 的特征值为 $1-\alpha\lambda_i(\mQ)$。迭代收敛当且仅当对所有 $i$ 有 $\abs{1-\alpha\lambda_i(\mQ)}<1$，即 $0<\alpha<2/\lambda_{\max}(\mQ)$。最佳常数步长 $\alpha=2/(\lambda_{\min}(\mQ)+\lambda_{\max}(\mQ))$ 平衡了两个极端特征值，每步把误差范数缩小为 $(\kappa-1)/(\kappa+1)$ 倍。对于 $\mQ=\diag(2,20)$，步长必须满足 $\alpha<0.1$（取 $\alpha=0.1$ 时，$x_2$ 方向的误差只是无休止地改变符号），最佳步长为 $\alpha=\tfrac1{11}$，因子为 $\tfrac9{11}=0.818$；这恰好就是从最坏初始点出发时精确线搜索在每次迭代中所选的步长。因此，好的常数步长需要知道特征值，而线搜索隐含地提供了这一信息。

对于一般的光滑函数 $f$，在满足 $\nabla^2f(\vx^*)\succ\vzero$ 的极小点 $\vx^*$ 附近，同样的界渐近地成立，其中 $\kappa$ 是 $\nabla^2f(\vx^*)$ 的条件数（Luenberger & Ye, 2016, Ch. 8）。因此，最速下降法可靠，但在尺度不良的问题上很慢。

(sec-06-05-scaling)=
### 尺度变换

条件数取决于度量变量所用的单位，因此变量替换可能有所帮助。用可逆矩阵 $\mL$ 作替换 $\vx=\mL\vy$，并对 $\tilde f(\vy)=f(\mL\vy)$ 应用最速下降法。由于 $\nabla\tilde f(\vy)=\mL^\T\nabla f(\vx)$（{prf:ref}`thm-01-02-chain-rule`），$\vy$ 中的一步 $-\alpha\nabla\tilde f(\vy)$ 就是 $\vx$ 中的一步 $-\alpha\mL\mL^\T\nabla f(\vx)$：**尺度化最速下降法**（scaled steepest descent）就是取 $\vd=-\mD\nabla f$、$\mD=\mL\mL^\T\succ\vzero$ 的下降法（{prf:ref}`prop-06-05-scaled-gradient`）。新变量下的黑塞矩阵为 $\mL^\T\nabla^2f\,\mL$，$\mL$ 的最佳选择使它接近 $\mI$ 的倍数。

对于{prf:ref}`ex-06-05-steepest-descent`中的二次函数，以 $1/\sqrt{10}$ 为单位度量第二个变量，即 $y_1=x_1-2$，$y_2=\sqrt{10}(x_2-2)$，就把 $f$ 变为 $y_1^2+y_2^2$，其等值线是圆（$\kappa=1$）。相应的 $\mD=\diag(1,\tfrac1{10})$ 在 $\vx^{(0)}$ 处给出方向 $\vd=-\diag(1,\tfrac1{10})(-12,-100)=(12,10)$，沿该方向的精确步长 $\alpha=\tfrac12$ 一次迭代即到达 $(2,2)$，与{prf:ref}`ex-06-05-newton-quadratic`中牛顿步到达的点相同。在实践中，简单的**对角尺度化**（diagonal scaling，把每个变量除以其典型量级）代价很低，并且往往能消除由单位混用造成的大部分病态，例如吨与千克、美元与百万美元。取 $\mD$ 等于黑塞矩阵的逆是理想的尺度变换，由此导出牛顿法。

(sec-06-05-newton)=
## 牛顿法

牛顿法在当前点附近用 $f$ 的二阶泰勒模型（{prf:ref}`thm-01-02-taylor`）

$$
\begin{aligned}
m_k(\vx)&=f(\vx^{(k)})+\nabla f(\vx^{(k)})^\T(\vx-\vx^{(k)})\\
&\quad+\tfrac12(\vx-\vx^{(k)})^\T\nabla^2f(\vx^{(k)})(\vx-\vx^{(k)}),
\end{aligned}
$$ (eq-06-05-newton-model)

代替 $f$，并移动到该模型的驻点。令 $\nabla m_k(\vx)=\nabla f(\vx^{(k)})+\nabla^2f(\vx^{(k)})(\vx-\vx^{(k)})=\vzero$，得到**牛顿方程**（Newton equations）和更新公式

$$
\begin{aligned}
\nabla^2f(\vx^{(k)})\,\vd^{(k)}&=-\nabla f(\vx^{(k)}),\\
\vx^{(k+1)}&=\vx^{(k)}+\vd^{(k)} .
\end{aligned}
$$ (eq-06-05-newton-step)

在{prf:ref}`alg-06-05-descent`的框架中，这就是方向 $\vd^{(k)}=-[\nabla^2f(\vx^{(k)})]^{-1}\nabla f(\vx^{(k)})$ 加上步长 $\alpha_k=1$；它是一维方法 {eq}`eq-06-04-newton-1d` 的推广。当 $\nabla^2f(\vx^{(k)})\succ\vzero$ 时，模型是凸二次函数，$\vx^{(k+1)}$ 是它的极小点，并且由{prf:ref}`prop-06-05-scaled-gradient`，$\vd^{(k)}$ 是下降方向。在计算中从不显式求逆：而是求解线性方程组 {eq}`eq-06-05-newton-step`，当黑塞矩阵正定时用 Cholesky 分解求解（约 $n^3/3$ 次运算）。

:::{prf:algorithm} 牛顿法
:label: alg-06-05-newton

**输入：** 二阶连续可微函数 $f$；初始点 $\vx^{(0)}$；容差 $\varepsilon>0$。
**输出：** 近似驻点。

1. **初始化。** 令 $k=0$。
2. **终止检验。** 若 $\norm{\nabla f(\vx^{(k)})}\le\varepsilon$，则停止。
3. **牛顿方向。** 求解 $\nabla^2f(\vx^{(k)})\vd^{(k)}=-\nabla f(\vx^{(k)})$，得 $\vd^{(k)}$。
4. **更新。** 令 $\vx^{(k+1)}=\vx^{(k)}+\vd^{(k)}$，$k\leftarrow k+1$，转第 2 步。
:::

:::{prf:proposition} 牛顿法一步求得凸二次函数的极小点
:label: prop-06-05-newton-quadratic

对于 $\mQ\succ\vzero$ 的二次函数 {eq}`eq-06-05-quadratic`，从任意初始点出发，牛顿法一步即可到达极小点 $\vx^*=\mQ^{-1}\vb$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

二次函数的泰勒模型 {eq}`eq-06-05-newton-model` 就是它本身，因此模型的极小点就是 $\vx^*$。具体地，$\vx^{(1)}=\vx^{(0)}-\mQ^{-1}(\mQ\vx^{(0)}-\vb)=\mQ^{-1}\vb$。
:::

:::{prf:example} 对同一个二次函数用牛顿法
:label: ex-06-05-newton-quadratic

从 $\vx^{(0)}=(-4,-3)$ 出发，对 $f(\vx)=(x_1-2)^2+10(x_2-2)^2$ 执行一次牛顿步。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

黑塞矩阵是常数，$\nabla^2f=\diag(2,20)$，且 $\nabla f(\vx^{(0)})=(-12,-100)$。牛顿方程 $\diag(2,20)\,\vd^{(0)}=(12,100)$ 给出 $\vd^{(0)}=(6,5)$，从而

$$
\begin{aligned}
\vx^{(1)}&=(-4,-3)+(6,5)=(2,2)=\vx^*,\\
f(\vx^{(1)})&=0 .
\end{aligned}
$$

牛顿方向 $(6,5)$ 径直指向极小点，而最速下降方向 $(12,100)$ 几乎沿着等值线密集的 $x_2$ 轴方向（{numref}`fig-06-05-zigzag`a）。
:::

如{prf:ref}`ex-06-05-newton-quadratic`所示，用黑塞矩阵的逆（这里为 $\mH^{-1}=\diag(\tfrac12,\tfrac1{20})$）对梯度进行尺度变换，恰好消除了条件数 $\kappa=10$ 造成的畸变。牛顿法不受任何线性变量替换的影响（它具有**仿射不变性**，affine invariance），因此黑塞矩阵的条件数在其中不起作用。

对于非二次函数，模型只是近似，牛顿法只有在极小点附近才快。

:::{prf:theorem} 牛顿法的局部二次收敛性
:label: thm-06-05-newton-quadratic

设 $f$ 二阶连续可微，其黑塞矩阵在 $\vx^*$ 附近 Lipschitz 连续，且 $\nabla f(\vx^*)=\vzero$，$\nabla^2f(\vx^*)\succ\vzero$。则存在 $\delta>0$，使得对每个满足 $\norm{\vx^{(0)}-\vx^*}<\delta$ 的 $\vx^{(0)}$，牛顿迭代 {eq}`eq-06-05-newton-step` 有定义，收敛到 $\vx^*$，并且对某个常数 $C$ 有 $\norm{\vx^{(k+1)}-\vx^*}\le C\norm{\vx^{(k)}-\vx^*}^2$。梯度范数 $\norm{\nabla f(\vx^{(k)})}$ 也二次收敛到零。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

（证明思路；Nocedal & Wright, 2006, Thm. 3.5。）记 $\mH_k=\nabla^2f(\vx^{(k)})$，在 $\vx^*$ 附近它可逆，且 $\norm{\mH_k^{-1}}\le2\norm{\nabla^2f(\vx^*)^{-1}}$。由于 $\nabla f(\vx^*)=\vzero$，

$$
\vx^{(k+1)}-\vx^*=\mH_k^{-1}\Bigl[\mH_k(\vx^{(k)}-\vx^*)-\bigl(\nabla f(\vx^{(k)})-\nabla f(\vx^*)\bigr)\Bigr].
$$

由积分形式的泰勒定理，方括号内的项等于 $\int_0^1\bigl(\mH_k-\nabla^2f(\vx^*+t(\vx^{(k)}-\vx^*))\bigr)(\vx^{(k)}-\vx^*)\,dt$，若 $L$ 是黑塞矩阵的 Lipschitz 常数，则其范数至多为 $\tfrac12L\norm{\vx^{(k)}-\vx^*}^2$。因此 $\norm{\vx^{(k+1)}-\vx^*}\le L\norm{\nabla^2f(\vx^*)^{-1}}\,\norm{\vx^{(k)}-\vx^*}^2$，再用{prf:ref}`thm-06-04-newton-quadratic`中的论证即可完成证明。
:::

两个假设都很重要。当 $\nabla^2f(\vx^*)$ 奇异时，收敛可能降为线性：对于 $f(\vx)=(x_1-2)^4+(x_1-2x_2)^2$，其在极小点 $(2,1)$ 处的黑塞矩阵是奇异的，每次牛顿步只把误差缩小为 $\tfrac23$ 倍（{prf:ref}`ex-08-09-newton-singular-hessian`）。而在远离 $\vx^*$ 处，牛顿步可能以一维方法（{doc}`06-04-line-search`）的所有方式失效：黑塞矩阵可能不定，步长可能极大，或者 $f$ 可能增大。

(sec-06-05-newton-modifications)=
### 使牛顿法稳健

两种修正把牛顿法变成全局收敛的算法，同时保持其快速的局部收敛速度。

**线搜索（阻尼牛顿法）。** 把 $\vd^{(k)}$ 用作方向，从 $\alpha=1$ 开始回溯来选择 $\alpha_k$（{prf:ref}`alg-06-04-backtracking`）。在解附近，单位步长满足 Armijo 条件（$c_1<\tfrac12$）而被接受，因此二次收敛速度得以保持；远离解时，较短的步长防止 $f$ 增大。

**黑塞矩阵修正。** 若 $\nabla^2f(\vx^{(k)})$ 不正定，牛顿方向可能指向上坡或指向鞍点。此时用 $\vd^{(k)}=-\bar{\mH}_k^{-1}\nabla f(\vx^{(k)})$ 代替它，其中 $\bar{\mH}_k\succ\vzero$ 接近黑塞矩阵：

- **平移：** $\bar{\mH}_k=\nabla^2f(\vx^{(k)})+\nu_k\mI$，其中 $\nu_k\ge0$ 恰好大到使 $\lambda_{\min}(\bar{\mH}_k)\ge\delta>0$（Levenberg–Marquardt 思想；Levenberg, 1944；Marquardt, 1963；$\nu_k$ 很大时，这一步变成短的最速下降步）；
- **特征值修正：** 利用谱分解 $\nabla^2f(\vx^{(k)})=\mV\diag\bigl(\lambda_i(\nabla^2f(\vx^{(k)}))\bigr)\mV^\T$（{prf:ref}`thm-01-01-spectral`），取 $\bar{\mH}_k=\mV\diag\bigl(\max\{\abs{\lambda_i(\nabla^2f(\vx^{(k)}))},\delta\}\bigr)\mV^\T$，它把负曲率翻转为正曲率；
- **修正 Cholesky 分解：** 执行 Cholesky 分解，每当某个主元过小时就增大相应的对角元，这样大约以一次分解的代价得到 $\bar{\mH}_k\succ\vzero$（实践中的首选方法）。

作为例示，取 $f(\vx)=\tfrac14x_1^4-\tfrac12x_1^2+x_2^2$，它在 $\vzero$ 处有鞍点，极小点为 $(\pm1,0)$，考虑点 $\vx=(0.2,1)$。在该点 $\nabla f=(x_1^3-x_1,\ 2x_2)=(-0.192,\ 2)$，$\nabla^2f=\diag(3x_1^2-1,\ 2)=\diag(-0.88,\ 2)$ 是不定的。纯牛顿步 $\vd=(-0.2182,\ -1)$ 在这里仍是下降方向（$\nabla f^\T\vd=-1.958$），但在 $x_1$ 坐标上它*朝着*鞍点移动：完整的一步落在 $(-0.0182,0)$，紧挨着 $\vzero$。特征值修正把 $-0.88$ 换成 $0.88$，得到 $\vd=(0.2182,\ -1)$，它远离鞍点；配合回溯，修正后的迭代在五次迭代内到达极小点 $(1,0)$，$f=-0.25$，精确到六位数字，且每次都取 $\alpha_k=1$。

**信赖域。** 使牛顿法全局化的另一种方式颠倒了两个决策的顺序：它不是先选方向再选步长，而是先确定二次模型 {eq}`eq-06-05-newton-model` 在多大范围内可信，即半径 $\Delta_k$，然后在球 $\norm{\vx-\vx^{(k)}}\le\Delta_k$ 上最小化该模型。$f$ 的实际下降量与模型预测下降量之比决定是否接受这一步，以及是扩大还是缩小 $\Delta_k$。即使黑塞矩阵不定，子问题也有解（球是紧集），而且其解具有 $-(\nabla^2f+\nu\mI)^{-1}\nabla f$（$\nu\ge0$）的形式（除去一种退化的“困难情形”），因此上面的平移策略自动出现，其中 $\nu$ 是球约束的乘子，由半径决定。信赖域（trust region）方法与线搜索方法同样被广泛使用（Nocedal & Wright, 2006, Ch. 4）。

(sec-06-05-quasi-newton)=
## 拟牛顿法

牛顿法每次迭代需要二阶导数和一次线性方程组求解；最速下降法两者都不需要，但很慢。**拟牛顿法**使用方向 $\vd^{(k)}=-\mD_k\vg^{(k)}$，其中 $\mD_k\succ\vzero$ 是黑塞矩阵*逆*的近似，在相邻迭代之间只利用梯度进行修正。修正公式的设计使 $\mD_{k+1}$ 重现沿最后一步观察到的曲率。由于 $\nabla f(\vx^{(k+1)})-\nabla f(\vx^{(k)})\approx\nabla^2f\,(\vx^{(k+1)}-\vx^{(k)})$，它要求满足**割线条件**（secant condition）

$$
\begin{aligned}
\mD_{k+1}\vy^{(k)}&=\vs^{(k)},\\
\vs^{(k)}&=\vx^{(k+1)}-\vx^{(k)},\\
\vy^{(k)}&=\vg^{(k+1)}-\vg^{(k)} .
\end{aligned}
$$ (eq-06-05-secant)

在满足 {eq}`eq-06-05-secant` 的众多矩阵中，**BFGS** 公式（Broyden, 1970；Fletcher, 1970；Goldfarb, 1970；Shanno, 1970）选取在某种加权范数下最接近 $\mD_k$ 的那一个：

$$
\begin{aligned}
\mD_{k+1}&=\Bigl(\mI-\gamma_k\vs^{(k)}\vy^{(k)\T}\Bigr)\mD_k\Bigl(\mI-\gamma_k\vy^{(k)}\vs^{(k)\T}\Bigr)\\
&\quad+\gamma_k\vs^{(k)}\vs^{(k)\T},\\
\gamma_k&=\frac1{\vy^{(k)\T}\vs^{(k)}} .
\end{aligned}
$$ (eq-06-05-bfgs)

若 $\mD_k\succ\vzero$ 且**曲率条件**（curvature condition）$\vy^{(k)\T}\vs^{(k)}>0$ 成立，则 $\mD_{k+1}\succ\vzero$；满足 Wolfe 条件 {eq}`eq-06-04-curvature` 的线搜索保证 $\vy^{(k)\T}\vs^{(k)}>0$，这就是拟牛顿程序使用 Wolfe 线搜索的原因。修正需要 $O(n^2)$ 次运算，不需要黑塞矩阵，也不需要求解线性方程组，并且 BFGS 在黑塞矩阵正定的极小点附近**超线性**收敛（Nocedal & Wright, 2006, Sec. 6.1 和 6.4）。较早的 DFP 公式（Davidon–Fletcher–Powell；Davidon, 1959；Fletcher & Powell, 1963）以同样方式修正，只是在黑塞矩阵近似中交换了 $\vs$ 和 $\vy$ 的角色；有限内存 BFGS（L-BFGS）只存储最近几对 $(\vs^{(k)},\vy^{(k)})$，能够处理有数百万个变量的问题。

在{prf:ref}`ex-06-05-steepest-descent`的二次函数上，取 $\mD_0=\mI$ 并采用精确线搜索的 BFGS 先走最速下降步，到达 $\vx^{(1)}=(-3.3922,\ 2.0647)$，并把矩阵修正为

$$
\mD_1=\begin{pmatrix}1.0122&-0.006146\\-0.006146&0.05007\end{pmatrix}.
$$

第二个方向 $-\mD_1\vg^{(1)}=(10.9240,\ -0.1311)$ 配合 $\alpha_1=0.4936$ 恰好落在 $\vx^*=(2,2)$ 上，再下一次修正给出 $\mD_2=\diag(\tfrac12,\tfrac1{20})=\mQ^{-1}$。一般地，对 $n$ 元凸二次函数采用精确线搜索时，BFGS 至多 $n$ 步终止，并最终得到 $\mD_n=\mQ^{-1}$。

(sec-06-05-conjugate-gradient)=
## 共轭梯度法

对于规模非常大的问题，即使存储一个 $n\times n$ 矩阵也负担过重。**共轭梯度法**（CG）只需要少数几个向量。若对 $i\ne j$ 有 $\vd^{(i)\T}\mQ\vd^{(j)}=0$，则称这些方向是 **$\mQ$-共轭**的（共轭方向，conjugate directions）。对于二次函数 {eq}`eq-06-05-quadratic`，沿 $n$ 个相互 $\mQ$-共轭的方向进行精确线搜索，至多 $n$ 步即可到达极小点，因为在由共轭方向确定的坐标系中，二次函数分解为 $n$ 个一维问题。CG 利用梯度在运行过程中生成这样的方向：

$$
\begin{aligned}
\vd^{(0)}&=-\vg^{(0)},\\
\vd^{(k+1)}&=-\vg^{(k+1)}+\beta_k\vd^{(k)},\\
\beta_k&=\frac{\vg^{(k+1)\T}\vg^{(k+1)}}{\vg^{(k)\T}\vg^{(k)}}
\end{aligned}
$$ (eq-06-05-cg)

（Fletcher–Reeves 公式；Hestenes & Stiefel, 1952 为线性方程组提出了 CG，Fletcher & Reeves, 1964 把它推广到非二次函数），其中 $\alpha_k$ 由精确线搜索得到（对于二次函数，$\alpha_k=-\vg^{(k)\T}\vd^{(k)}/\vd^{(k)\T}\mQ\vd^{(k)}$）。新方向是用旧方向的倍数修正过的最速下降方向；这一修正恰好消除了锯齿现象。

在{prf:ref}`ex-06-05-steepest-descent`的二次函数上，CG 的第一步就是到达 $(-3.3922,\ 2.0647)$ 的最速下降步。随后 $\beta_0=\norm{\vg^{(1)}}^2/\norm{\vg^{(0)}}^2=0.01163$，新方向 $\vd^{(1)}=(10.7845,\ -1.2941)+0.01163\,(12,100)=(10.9240,\ -0.1311)$ 与 $\vd^{(0)}$ $\mQ$-共轭，精确步长 $\alpha_1=0.4936$ 到达 $\vx^*=(2,2)$：两个变量用两步，而最速下降法要十七步（{numref}`fig-06-05-zigzag`a）。在采用精确线搜索的二次函数上，CG 和从 $\mD_0=\mI$ 开始的 BFGS 产生相同的迭代点，正如本例所示。对于非二次函数，CG 每 $n$ 次迭代重新开始一次（或者使用 Polak–Ribière 公式（Polak & Ribière, 1969）$\beta_k=\vg^{(k+1)\T}(\vg^{(k+1)}-\vg^{(k)})/\vg^{(k)\T}\vg^{(k)}$），并视问题不同而线性或超线性收敛。

(sec-06-05-practice)=
## 实践中的导数与终止

本章所有方法都需要梯度，牛顿法还需要黑塞矩阵。当它们无法手工推导时，常见的来源有三种。

- **有限差分。** 前向差分 $\partial f/\partial x_j\approx\bigl(f(\vx+h\ve_j)-f(\vx)\bigr)/h$ 每求一次梯度需要额外计算 $n$ 个函数值。其截断误差为 $h$ 量级，舍入误差为 $\epsilon_{\text{mach}}/h$ 量级，因此 $h\approx\sqrt{\epsilon_{\text{mach}}}\,\max\{1,\abs{x_j}\}\approx10^{-8}\max\{1,\abs{x_j}\}$ 使两者平衡，所得梯度的有效数字位数约为函数值的一半。中心差分 $\bigl(f(\vx+h\ve_j)-f(\vx-h\ve_j)\bigr)/(2h)$（central difference）更精确，但需要 $2n$ 个函数值。
- **自动微分。** 对计算 $f$ 的代码应用链式法则的程序，能以 $f$ 计算代价的一个小倍数（反向模式）给出精确到舍入误差的梯度，且与 $n$ 无关。建模语言和现代数值库通常都这样做（{doc}`09-06-code-nlp`）。
- **拟牛顿修正** 用梯度差代替黑塞矩阵（{ref}`sec-06-05-quasi-newton`），这就是为什么使用精确黑塞矩阵的牛顿法主要用于二阶导数计算代价低的情形。

**终止。** {prf:ref}`alg-06-05-descent`中的检验 $\norm{\nabla f(\vx^{(k)})}\le\varepsilon$ 依赖于 $f$ 和 $\vx$ 的单位。实用程序采用相对形式，例如 $\norm{\nabla f(\vx^{(k)})}\le\varepsilon\max\{1,\abs{f(\vx^{(k)})}\}$，并结合对步长的检验和迭代次数上限，同时报告是哪一项检验结束了运行（{doc}`06-04-line-search`，{ref}`sec-06-04-descent-step`）。因达到迭代次数上限而结束的运行——最速下降法在病态问题上很容易如此——并没有找到解。

(sec-06-05-comparison)=
## 方法的比较

```{list-table} $n$ 元光滑函数 $f$ 的无约束下降法
:name: tab-06-05-comparison
:header-rows: 1

* - 方法
  - 方向
  - 每次迭代的工作量
  - 存储量
  - 局部收敛速度
* - 最速下降法
  - $-\vg^{(k)}$
  - 一次梯度
  - $O(n)$
  - 线性，因子至多为 $((\kappa-1)/(\kappa+1))^2$
* - 共轭梯度法
  - $-\vg^{(k)}+\beta_{k-1}\vd^{(k-1)}$
  - 一次梯度
  - $O(n)$
  - 在二次函数上有限步终止（$\le n$ 步）；否则为线性到超线性
* - 拟牛顿法（BFGS）
  - $-\mD_k\vg^{(k)}$
  - 一次梯度，$O(n^2)$ 的修正
  - $O(n^2)$
  - 超线性
* - 牛顿法
  - $-[\nabla^2f]^{-1}\vg^{(k)}$
  - 梯度、黑塞矩阵、$O(n^3)$ 的求解
  - $O(n^2)$
  - 二次
```

选择哪种方法取决于有哪些信息可用、哪些计算代价高。若二阶导数计算代价低且 $n$ 适中，带线搜索和黑塞矩阵修正的牛顿法很难被超越。若只有梯度可用，BFGS 是标准选择，而对于 $n$ 非常大的问题则用 L-BFGS 或 CG。最速下降法作为构件，以及作为解释其他方法时的参照基准，仍然很重要。

## 本章回顾

### 小结

- 下降法（{prf:ref}`alg-06-05-descent`）沿方向 $\vd=-\mD\nabla f$（$\mD\succ\vzero$）移动，这些方向都是下降方向（{prf:ref}`prop-06-05-scaled-gradient`）；采用 Wolfe 步长且方向与梯度不太接近正交时，它们使 $\nabla f$ 趋于零。
- 最速下降法（{prf:ref}`alg-06-05-steepest-descent`）使用 $\vd=-\nabla f$；在二次函数上精确步长为 {eq}`eq-06-05-sd-step`；相邻步相互正交 {eq}`eq-06-05-orthogonal`，呈锯齿形（{numref}`fig-06-05-zigzag`）。
- 采用精确线搜索时，每个极限点都是驻点（{prf:ref}`prop-06-05-sd-convergence`）；对变量作尺度变换（{ref}`sec-06-05-scaling`）会改变 $\kappa$，从而改变速度。
- 其收敛是线性的，因子至多为 $((\kappa-1)/(\kappa+1))^2$（{prf:ref}`thm-06-05-sd-rate`）；对 $(x_1-2)^2+10(x_2-2)^2$ 从 $(-4,-3)$ 出发时，每步的因子为 $0.1018$（{prf:ref}`ex-06-05-steepest-descent`）。
- 牛顿法 {eq}`eq-06-05-newton-step` 一步即可求得凸二次函数的极小点（{prf:ref}`prop-06-05-newton-quadratic`，{prf:ref}`ex-06-05-newton-quadratic`），并在黑塞矩阵正定的极小点附近二次收敛（{prf:ref}`thm-06-05-newton-quadratic`）；线搜索和黑塞矩阵修正使它稳健。
- 拟牛顿法（BFGS，{eq}`eq-06-05-bfgs`）和共轭梯度法（{eq}`eq-06-05-cg`）只使用梯度；两者采用精确线搜索时，都至多 $n$ 步即可求得 $n$ 元凸二次函数的极小点（{numref}`tab-06-05-comparison`）。

### 常见错误

- **以为最速下降方向指向极小点。** 负梯度垂直于等值线，而不是指向 $\vx^*$；在狭长的等值线上，它几乎与指向 $\vx^*$ 的方向正交。
- **在最速下降法中使用单位步长。** $\vx-\nabla f(\vx)$ 没有理由是一个好的点；步长必须由线搜索得到（这里 $\alpha_k$ 在 $0.05065$ 和 $0.4434$ 之间交替）。
- **对黑塞矩阵求逆。** 应求解 $\nabla^2f\,\vd=-\nabla f$（Cholesky 分解）；不要构造 $[\nabla^2f]^{-1}$。
- **从远处开始运行纯牛顿法。** 没有线搜索和正定的（修正）黑塞矩阵，牛顿法可能使 $f$ 增大、发散，或收敛到鞍点或极大点。
- **终止检验中不用范数。** 对向量而言 $\abs{\nabla f(\vx)}\le\varepsilon$ 没有意义；应使用 $\norm{\nabla f(\vx)}\le\varepsilon$。


## 参考文献

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Sections 8.6–8.8:
  steepest descent, Newton's method and its modifications, conjugate
  directions and quasi-Newton methods.)
- Broyden, C. G. (1970). The convergence of a class of double-rank
  minimization algorithms 1. General considerations. *Journal of the
  Institute of Mathematics and Its Applications, 6*(1), 76–90.
  <https://doi.org/10.1093/imamat/6.1.76>
- Cauchy, A. (1847). Méthode générale pour la résolution des systèmes
  d'équations simultanées. *Comptes Rendus de l'Académie des Sciences, 25*,
  536–538.
- Davidon, W. C. (1959). *Variable metric method for minimization* (Report
  ANL-5990). Argonne National Laboratory.
- Fletcher, R. (1970). A new approach to variable metric algorithms. *The
  Computer Journal, 13*(3), 317–322. <https://doi.org/10.1093/comjnl/13.3.317>
- Fletcher, R., & Powell, M. J. D. (1963). A rapidly convergent descent method
  for minimization. *The Computer Journal, 6*(2), 163–168.
  <https://doi.org/10.1093/comjnl/6.2.163>
- Fletcher, R., & Reeves, C. M. (1964). Function minimization by conjugate
  gradients. *The Computer Journal, 7*(2), 149–154.
  <https://doi.org/10.1093/comjnl/7.2.149>
- Goldfarb, D. (1970). A family of variable-metric methods derived by
  variational means. *Mathematics of Computation, 24*(109), 23–26.
- Hestenes, M. R., & Stiefel, E. (1952). Methods of conjugate gradients for
  solving linear systems. *Journal of Research of the National Bureau of
  Standards, 49*(6), 409–436. <https://doi.org/10.6028/jres.049.044>
- Levenberg, K. (1944). A method for the solution of certain nonlinear
  problems in least squares. *Quarterly of Applied Mathematics, 2*(2),
  164–168. <https://doi.org/10.1090/qam/10666>
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming* (4th
  ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3> (Chapter 8:
  steepest descent, its rate and the Kantorovich inequality, Newton's method;
  Chapter 9: conjugate-direction methods; Chapter 10: quasi-Newton methods.)
- Marquardt, D. W. (1963). An algorithm for least-squares estimation of
  nonlinear parameters. *Journal of the Society for Industrial and Applied
  Mathematics, 11*(2), 431–441. <https://doi.org/10.1137/0111030>
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Chapters 2–3:
  line-search methods, Zoutendijk's theorem, rates of steepest descent and
  Newton's method, Hessian modifications; Chapter 4: trust-region methods;
  Chapter 5: conjugate gradients; Chapter 6: quasi-Newton methods and BFGS;
  Chapter 8: finite differences and automatic differentiation.)
- Polak, E., & Ribière, G. (1969). Note sur la convergence de méthodes de
  directions conjuguées. *Revue Française d'Informatique et de Recherche
  Opérationnelle, Série Rouge, 3*(16), 35–43.
- Shanno, D. F. (1970). Conditioning of quasi-Newton methods for function
  minimization. *Mathematics of Computation, 24*(111), 647–656.
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 21: nonlinear programming algorithms, gradient method.)
- Zoutendijk, G. (1970). Nonlinear programming, computational methods. In
  J. Abadie (Ed.), *Integer and nonlinear programming* (pp. 37–86).
  North-Holland.
