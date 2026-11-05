# 线性化技术

## 概述

许多模型*几乎*是线性的：目标函数中出现两个决策的乘积，费用取决于绝对偏差，要极小化若干线性表达式中的最坏情形，或者目标是诸如“每美元效益”之类的比值。本章汇集了一些标准技术，它们用线性约束、连续辅助变量以及（在不可避免时）0–1 变量（binary variable）把这类项精确地改写，或者给出紧的界。两个 0–1 变量之积就是{doc}`07-01-logical-constraints`中的逻辑“与”；一个 0–1 变量与一个有界连续变量之积可以用四个不等式精确线性化（linearization）；两个连续变量之积则被 **McCormick 包络**（McCormick envelope）所包围，本章由到边界的距离的四个乘积推导出它。线性函数的绝对值、最大值和最小值出现在模型的“凸的一侧”时是线性的，否则就需要 0–1 变量；本章只陈述一次这条规则，再把它应用于若干情形。最后，**Charnes–Cooper 变换**（Charnes–Cooper transformation）把两个线性函数之比转化为线性规划（linear program, LP）。每种技术都配有一个基础设施算例，其解均经过数值验证。除{doc}`07-01-logical-constraints`外，本章还用到线性规划建模（{doc}`03-01-lp-model-building`）和凸函数的上图（{prf:ref}`def-02-02-epigraph`）。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 对 0–1 变量之积以及 0–1 变量与有界连续变量之积进行线性化，并知道何时可以去掉模型的一侧；
- 推导箱体上 $xy$ 的 McCormick 不等式，说明它们何时是精确的，并通过分割加以收紧；
- 为线性函数的绝对值和最大值/最小值建模，判断用线性规划是否足够、还是需要 0–1 变量；
- 把线性分式规划转化为线性规划并还原其解；
- 通过代入 0–1 变量的整数值来检验每一种线性化。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\omega$, $\omega_{ij}$ | 代替乘积的辅助变量（$\omega=xy$，$\omega_{ij}=\delta_i\delta_j$） |
| $\underline x$, $\overline x$ | 变量 $x$ 的下界和上界 |
| $\eta$, $\eta_i$ | 代替绝对值或最小值/最大值的辅助变量 |
| $\vx^+$, $\vx^-$ | 正部与负部，$\vx=\vx^+-\vx^-$，$\vx^\pm\ge\vzero$ |
| $\vr_k^\T\vx+r_{k0}$ | 最大值或最小值中的第 $k$ 个仿射函数，$k=1,\dots,K$ |
| $\vp^\T\vx+p_0$, $\vq^\T\vx+q_0$ | 线性分式目标函数的分子和分母 |
| $\sigma$, $\vy$ | Charnes–Cooper 变量，$\sigma=1/(\vq^\T\vx+q_0)$，$\vy=\sigma\vx$ |

如果对原变量的每一个取值，新约束所允许的辅助变量值恰好就是非线性项所定义的值，则称该**线性化**是*精确的*；如果允许的值更多，则称它是一个*松弛*。本章中 0–1 变量之积和绝对值的线性化都是精确的；连续变量之积的 McCormick 包络则是一个松弛。

(sec-07-02-binary-products)=
## 0–1 变量之积

设 $\delta_1,\delta_2\in\{0,1\}$。它们的乘积为 $1$ 当且仅当二者都为 $1$，因此 $\omega=\delta_1\delta_2$ 就是逻辑“与”，可用 {eq}`eq-07-01-and` 线性化：

$$
\begin{aligned}
\omega&\le\delta_1,\\
\omega&\le\delta_2,\\
\omega&\ge\delta_1+\delta_2-1,\\
\omega&\ge0 .
\end{aligned}
$$ (eq-07-02-bin-bin)

代入四个数对：$(0,0)$、$(0,1)$ 和 $(1,0)$ 给出 $0\le\omega\le0$，$(1,1)$ 给出 $1\le\omega\le1$。注意 $\omega$ 不必声明为 0–1 变量：只要 $\delta_1,\delta_2$ 取整数，约束就会迫使 $\omega$ 取正确的整数值。$k$ 个 0–1 变量之积也可同样线性化：对每个因子有 $\omega\le\delta_i$，并且 $\omega\ge\sum_i\delta_i-(k-1)$。由于 0–1 变量满足 $\delta^2=\delta$，**0–1 二次规划**（quadratic 0–1 program）$\max\ \sum_i c_i\delta_i+\sum_{i<j}c_{ij}\delta_i\delta_j$ 就变成一个混合整数线性规划（mixed-integer linear program, MILP），每个非零的 $c_{ij}$ 对应一个 $\omega_{ij}$（Glover & Woolsey, 1974; Glover, 1975）。

**单侧形式。** 当 $\omega$ 只出现在目标函数中时，{eq}`eq-07-02-bin-bin` 的一半往往是多余的。在 $c_{ij}>0$ 的极大化问题中，优化会把 $\omega_{ij}$ 往上推，因此只需要上界 $\omega_{ij}\le\delta_i$、$\omega_{ij}\le\delta_j$；当 $c_{ij}<0$ 时，优化会把 $\omega_{ij}$ 往下压，因此只需要 $\omega_{ij}\ge\delta_i+\delta_j-1$ 和 $\omega_{ij}\ge0$。去掉错误的一侧是常见错误，下面的例子说明了这一点。

:::{prf:example} 具有协同效应的项目组合
:label: ex-07-02-synergy

四个铁路项目竞争 1.05 亿美元的预算。项目 $i$ 的成本为 $a_i=60,20,25,15$，单独实施的效益为 $b_i=80,18,24,16$（单位均为百万美元）。有些项目两两之间存在相互作用：项目 1 和 2 一起实施可增加 2 000 万美元效益（共用车站），1 和 3 增加 1 500 万美元，3 和 4 增加 600 万美元，而 2 和 3 服务同一批乘客，一起实施会*减少* 2 000 万美元。选择总效益最大的项目组合，并检查如果对每一对都只保留乘积线性化中的上界 $\omega_{ij}\le\delta_i$、$\omega_{ij}\le\delta_j$ 会出什么问题。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*模型。* 若选择项目 $i$，则 $\delta_i=1$；令 $\omega_{ij}=\delta_i\delta_j$，相互作用的项目对集合为 $P=\{(1,2),(1,3),(2,3),(3,4)\}$，

$$
\begin{aligned}
\max\ & 80\delta_1+18\delta_2+24\delta_3+16\delta_4\\
& +20\omega_{12}+15\omega_{13}-20\omega_{23}+6\omega_{34}\\
\st\ & 60\delta_1+20\delta_2+25\delta_3+15\delta_4\le105,\\
& \omega_{ij}\le\delta_i, \qquad (i,j)\in P,\\
& \omega_{ij}\le\delta_j, \qquad (i,j)\in P,\\
& \omega_{ij}\ge\delta_i+\delta_j-1, \qquad (i,j)\in P,\\
& \omega_{ij}\ge0, \qquad (i,j)\in P,\\
& \vdelta\in\{0,1\}^4 .
\end{aligned}
$$

*求解。* 枚举 $16$ 个项目组合（或求解该 MILP）得到最优组合 $\{1,3,4\}$：成本 $60+25+15=100$，效益 $80+24+16+15+6=141$，即 1.41 亿美元。次优组合 $\{1,2,3\}$ 的成本恰为 $105$，价值 $80+18+24+20+15-20=137$。

*错误的单侧模型。* 若对所有项目对都只保留上界 $\omega_{ij}\le\delta_i,\delta_j$，则当 $\delta_2=\delta_3=1$ 时没有任何约束阻止 $\omega_{23}=0$：惩罚消失了，求解器报告 $\{1,2,3\}$，“效益”为 $157$。对三个正协同项保留上界，对负协同项保留 $\omega_{23}\ge\delta_2+\delta_3-1$，即可用更少的约束行得到正确的最优值 $141$；其线性规划松弛的界为 $145.5$。
:::

(sec-07-02-binary-continuous)=
## 0–1 变量与连续变量之积

设 $\delta\in\{0,1\}$，$\underline x\le x\le\overline x$，且 $\underline x\ge0$。乘积 $\omega=\delta x$ 在 $\delta=0$ 时等于 $0$，在 $\delta=1$ 时等于 $x$。两种情形都由下式刻画：

$$
\begin{aligned}
\underline x\,\delta&\le\omega\le\overline x\,\delta,\\
x-\overline x(1-\delta)&\le\omega\le x-\underline x(1-\delta) .
\end{aligned}
$$ (eq-07-02-bin-cont)

当 $\delta=0$ 时，第一行给出 $\omega=0$，第二行成为 $x-\overline x\le0\le x-\underline x$，它总成立。当 $\delta=1$ 时，第二行给出 $\omega=x$，第一行成为 $\underline x\le x\le\overline x$。这两对不等式都是采用最紧常数（$x$ 的界）的大 $M$ 约束；它们与箱体 $[0,1]\times[\underline x,\overline x]$ 上 {eq}`eq-07-02-mccormick` 的 McCormick 不等式一致，这就解释了它们为什么是精确的（Glover, 1975）。

在常见的 $\underline x=0$ 情形下，第一行中的下界是 $\omega\ge0$，切不可遗漏。*（若没有 $\omega\ge0$，组合 $\delta=0$、$\omega=x-\overline x<0$ 就是可行的，而一个极大化 $-\omega$ 的模型会利用这一点。）* 当 $\underline x>0$ 时，第二行右边的不等式可以换成较弱的 $\omega\le x$，对 0–1 值的 $\delta$ 它仍是精确的，但给出的松弛较弱。

与 0–1 变量之积一样，当 $\omega$ 的目标函数系数的符号已知时，单侧形式就足够了：优化希望取大的 $\omega$ 只需要上界 $\omega\le\overline x\delta$ 和 $\omega\le x-\underline x(1-\delta)$。

:::{prf:example} 带数量折扣的供货合同
:label: ex-07-02-contract

某道路管理机构每年需要 $800$ 吨沥青结合料。供应商 1 以每吨 80 美元的价格供货，最多 $500$ 吨。供应商 2 最多供货 $600$ 吨，价格为每吨 90 美元；若该机构签订一份费用为 6000 美元的年度合同，则价格为每吨 70 美元。极小化年度费用。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*模型。* 设 $x_1,x_2$ 为购买的吨数；若签订合同，则 $\delta=1$。费用 $80x_1+90x_2-20\,\delta x_2+6000\,\delta$ 中含有乘积 $\omega=\delta x_2$，其中 $x_2\in[0,600]$：

$$
\begin{aligned}
\min\ & 80x_1+90x_2-20\omega+6000\delta\\
\st\ & x_1+x_2=800,\\
& 0\le x_1\le500,\\
& 0\le x_2\le600,\\
& 0\le\omega\le600\,\delta,\\
& x_2-600(1-\delta)\le\omega\le x_2,\\
& \delta\in\{0,1\} .
\end{aligned}
$$

*求解。* 不签合同时，从较便宜的供应商 1 购买 $500$ 吨，从供应商 2 购买 $300$ 吨：$40\,000+27\,000=\$67\,000$。签订合同后供应商 2 更便宜，于是从它那里购买 $600$ 吨，从供应商 1 购买 $200$ 吨：$16\,000+42\,000+6000=\$64\,000$。MILP 给出 $\delta^*=1$，$x_1^*=200$，$x_2^*=\omega^*=600$，费用为 64 000 美元。只要从供应商 2 购买至少 $6000/20=300$ 吨，合同就是划算的。由于 $\omega$ 的费用系数为负（希望它取大值），仅用 $\omega\le600\delta$ 和 $\omega\le x_2$ 这两行就已经能得到这个答案。
:::

**一般整数变量。** 有界整数变量 $x\in\{0,1,\dots,U\}$ 可以写成二进制形式

$$
\begin{aligned}
x&=\sum_{k=0}^{K-1}2^k\delta_k,\qquad \delta_k\in\{0,1\},\\
K&=\lceil\log_2(U+1)\rceil,
\end{aligned}
$$ (eq-07-02-binary-expansion)

当 $U+1$ 不是 2 的幂时，还要加上约束行 $x\le U$。于是乘积 $xy$（其中 $y\in[0,\overline y]$ 为有界连续变量）就变成 $\sum_k2^k\omega_k$，其中 $\omega_k=\delta_ky$，每一项都用 {eq}`eq-07-02-bin-cont` 线性化：用 $K$ 个新变量和 $4K$ 行代替一个非线性项。例如，公交车数量 $x\in\{0,\dots,12\}$ 乘以每辆车的连续载客量，需要 $K=4$ 个 0–1 变量（$2^4=16>12$）。同样的展开可把两个整数变量之积化为若干“0–1 变量×整数变量”乘积之和，从而化为线性形式（Glover, 1975）。

(sec-07-02-mccormick)=
## McCormick 包络

现在设两个因子都是连续变量，$\omega=xy$，其中 $x\in[\underline x,\overline x]$，$y\in[\underline y,\overline y]$。集合 $\{(x,y,\omega):\omega=xy\}$ 是一个鞍面，不是凸集，因此任何线性规划都无法描述它。用线性不等式最多只能把它包围起来。箱体中的每一点都满足四个显然的不等式：

$$
\begin{aligned}
(x-\underline x)(y-\underline y)&\ge0,\\
(\overline x-x)(\overline y-y)&\ge0,\\
(\overline x-x)(y-\underline y)&\ge0,\\
(x-\underline x)(\overline y-y)&\ge0 .
\end{aligned}
$$

把每个乘积展开并用 $\omega$ 代替 $xy$，就得到 **McCormick 不等式**（McCormick inequalities）（McCormick, 1976）：

$$
\begin{aligned}
\omega&\ge\underline y\,x+\underline x\,y-\underline x\,\underline y,\\
\omega&\ge\overline y\,x+\overline x\,y-\overline x\,\overline y,\\
\omega&\le\underline y\,x+\overline x\,y-\overline x\,\underline y,\\
\omega&\le\overline y\,x+\underline x\,y-\underline x\,\overline y .
\end{aligned}
$$ (eq-07-02-mccormick)

前两个是 $xy$ 的**下估计函数**（under-estimator），后两个是**上估计函数**（over-estimator）。它们的最大值和最小值分别是 $xy$ 在箱体上的凸包络和凹包络（Al-Khayyal & Falk, 1983），即位于 $xy$ 下方的最大凸函数和位于其上方的最小凹函数；在箱体上有效的、关于 $(x,y,\omega)$ 的任何其他线性不等式都不能切得更深。

**何时是精确的。** 若 $x$ 取在某个界上，比如 $x=\underline x$，则第一个和第四个不等式给出 $\omega\ge\underline x\,y$ 和 $\omega\le\underline x\,y$，于是 $\omega=xy$。在 $x=\overline x$、$y=\underline y$ 和 $y=\overline y$ 处同样如此。因此包络在箱体的边界上是精确的，特别地，只要有一个因子是 0–1 变量就是精确的：取 $x=\delta\in[0,1]$，{eq}`eq-07-02-mccormick` 恰好就是 {eq}`eq-07-02-bin-cont` 和 {eq}`eq-07-02-bin-bin`。

**有多宽松。** 在箱体内部，包络是一个松弛。沿着经过 $(\underline x,\underline y)$ 和 $(\overline x,\overline y)$ 的对角线，乘积是凸的，上估计函数位于其上方；沿着另一条对角线，乘积是凹的，下估计函数位于其下方（{numref}`fig-07-02-mccormick`）。任一侧的最大间隙都是 $(\overline x-\underline x)(\overline y-\underline y)/4$，在箱体中心取得。把一个区间的宽度减半，间隙也减半：这就是**分割**（partitioning）的思想。把 $[\underline x,\overline x]$ 分成 $K$ 段，用 0–1 变量选出包含 $x$ 的那一段（如{doc}`07-01-logical-constraints`中那样），并用所选区段的界写出包络；间隙便按 $1/K$ 的速度缩小。

```{figure} figures/07-02-mccormick.svg
:name: fig-07-02-mccormick
:width: 640px
:class: dark-light
:alt: 乘积 xy 沿箱体两条对角线的取值，以及 McCormick 下估计函数和上估计函数

{prf:ref}`ex-07-02-mccormick` 中箱体 $[1,4]\times[1,3]$ 上的乘积 $\omega=xy$，分别沿主对角线（左图，$xy$ 在其上是凸的）和另一条对角线（右图，$xy$ 在其上是凹的）。阴影带位于两个 McCormick 下估计函数中较大者与两个上估计函数中较小者之间；它在箱体的角点处与 $xy$ 相接，在中心处最宽，宽度为 $(3\cdot2)/4=1.5$。
```

:::{prf:example} 双线性目标函数的定界
:label: ex-07-02-mccormick

在 $2x+3y\le12$、$1\le x\le4$、$1\le y\le3$ 的约束下极大化 $\omega=xy$。（可以把 $x$ 看作价格，把 $y$ 看作由预算线联系起来的数量。）求出精确最优解和 McCormick 松弛给出的上界，并说明把 $[1,4]$ 分成 $K$ 个等长区段时该界如何改进。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*精确最优解。* 在最优解处预算约束是有效约束，$y=(12-2x)/3$，而 $xy=x(12-2x)/3$ 在 $x=3$ 处取得最大值：$\vx^*=(3,2)$，$xy=6$。

*McCormick 松弛。* 对 $\underline x=1$、$\overline x=4$、$\underline y=1$、$\overline y=3$ 用 {eq}`eq-07-02-mccormick` 将 $xy$ 替换为 $\omega$，得到线性规划

$$
\begin{aligned}
\max\ & \omega\\
\st\ & \omega\ge x+y-1,\\
& \omega\ge3x+4y-12,\\
& \omega\le x+4y-4,\\
& \omega\le3x+y-3,\\
& 2x+3y\le12,\\
& 1\le x\le4,\\
& 1\le y\le3 .
\end{aligned}
$$

其最优解位于两个上估计函数与预算线的交点：由 $x+4y-4=3x+y-3$ 得 $3y=2x+1$，再结合 $2x+3y=12$ 解得 $x=\tfrac{11}4$，$y=\tfrac{13}6$，$\omega=\tfrac{89}{12}\approx7.417$。松弛把真实最大值高估了 $1.417$，而该点处真实的乘积为 $\tfrac{11}4\cdot\tfrac{13}6=\tfrac{143}{24}\approx5.958$。

*分割。* 把 $[1,4]$ 分成 $K$ 个等长区段，并取各区段上界中的最好者（这正是每段配一个 0–1 变量的 MILP 所计算的），得到

| $K$ | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| $\max xy$ 的上界 | 7.417 | 6.722 | 6.367 | 6.185 | 6.093 |
| 与 6 的间隙 | 1.417 | 0.722 | 0.367 | 0.185 | 0.093 |

（精确值为 $\tfrac{89}{12}$、$\tfrac{121}{18}$、$\tfrac{191}{30}$、$\tfrac{167}{27}$、$\tfrac{1243}{204}$。）正如所预言的，$K$ 每加倍一次，间隙就减半。
:::

(sec-07-02-convex-side)=
## 绝对值、最大值与最小值：凸的一侧

以下各项技术都基于一个观察。当凸函数 $\chi(\vx)$ 是若干仿射函数的最大值时，

$$
\begin{aligned}
&\chi(\vx)=\max_{k=1,\dots,K}\bigl(\vr_k^\T\vx+r_{k0}\bigr)\le\eta\\
\iff\ &\vr_k^\T\vx+r_{k0}\le\eta\quad\text{对所有 }k,
\end{aligned}
$$ (eq-07-02-max-le)

就可以用线性约束从**上方**界定它，因为最大值不超过 $\eta$ 当且仅当每一项都不超过 $\eta$。任何凸分段线性函数都具有这种形式，而 $\abs{x}=\max\{x,-x\}$ 是最简单的情形。同理，凹分段线性函数（若干仿射函数的最小值）可以从**下方**界定：$\min_k(\cdot)\ge\eta$ 当且仅当每一项都 $\ge\eta$。更一般地，$\min_af(a,\vx)\ge0$ 等价于对所有 $a$ 有 $f(a,\vx)\ge0$，$\max_af(a,\vx)\le0$ 等价于对所有 $a$ 有 $f(a,\vx)\le0$。

相反方向的陈述 $\max_k(\cdot)\ge\eta$ 或 $\min_k(\cdot)\le\eta$ 说的是*至少有一项*满足该界。它们是析取条件，需要 0–1 变量（$K$ 个中至少 $1$ 个，{eq}`eq-07-01-k-of-m`）：

$$
\begin{aligned}
\vr_k^\T\vx+r_{k0}&\ge\eta-M_k(1-\delta_k),\quad k=1,\dots,K,\\
\sum_{k=1}^K\delta_k&\ge1,\\
\delta_k&\in\{0,1\},\quad k=1,\dots,K .
\end{aligned}
$$ (eq-07-02-max-ge)

在目标函数中，这条规则表述为：**极小化凸**（或极大化凹）分段线性函数是线性规划问题，借助**上图**（epigraph）变量 $z$ 即可实现：

$$
\begin{gathered}
\begin{aligned}
\min_{\vx}\ & \max_k\bigl(\vr_k^\T\vx+r_{k0}\bigr)\\
\st\ & \mA\vx\ge\vb
\end{aligned}\\
\Updownarrow\\
\begin{aligned}
\min_{\vx,z}\ & z\\
\st\ & z\ge\vr_k^\T\vx+r_{k0},\qquad k=1,\dots,K,\\
& \mA\vx\ge\vb .
\end{aligned}
\end{gathered}
$$ (eq-07-02-minmax)

在最优解处，$z$ 落在最大的那一项上，因为任何更大的 $z$ 都还能减小。**极大化凸**函数（或极小化凹函数）则不是线性规划问题：在 $\max_{\vx}\max_k(\cdot)$ 中，上图约束行 $z\ge\vr_k^\T\vx+r_{k0}$ 毫无用处（优化会把 $z$ 推向无穷大），而亚图（hypograph）条件 $z\le\max_k(\cdot)$ 正是析取条件 {eq}`eq-07-02-max-ge`。可以使用 0–1 变量，

$$
\begin{aligned}
\max_{\vx,z,\vdelta}\ & z\\
\st\ & z\le\vr_k^\T\vx+r_{k0}+M_k(1-\delta_k),\qquad k=1,\dots,K,\\
& \textstyle\sum_k\delta_k=1,\\
& \mA\vx\ge\vb,\\
& \delta_k\in\{0,1\},\qquad k=1,\dots,K,
\end{aligned}
$$

也可以简单地求解 $K$ 个线性规划 $\max\{\vr_k^\T\vx+r_{k0}:\mA\vx\ge\vb\}$ 并保留最好的一个；当 $K$ 较小时，后者往往更快。

**精确地定义最小值。** 有时一个变量必须*等于*某个最小值 $\eta=\min\{x_1,x_2\}$，例如当通过某个交汇点的流量是两个容量中的较小者、而这两个容量本身又是决策变量时。约束行 $\eta\le x_1$、$\eta\le x_2$ 只给出 $\eta\le\min$；要得到等式，还需要析取条件“$\eta\ge x_1$ 或 $\eta\ge x_2$”：

$$
\begin{aligned}
\eta&\le x_1,\\
\eta&\le x_2,\\
\eta&\ge x_1-M(1-\delta),\\
\eta&\ge x_2-M\delta,\qquad \delta\in\{0,1\},
\end{aligned}
$$ (eq-07-02-eta-min)

其中 $M\ge\max\abs{x_1-x_2}$（在变量的界上取最大）。若 $x_1\le x_2$，则可取 $\delta=1$，并迫使 $\eta=x_1$；否则 $\delta=0$，$\eta=x_2$。对于 $\eta=\max\{x_1,x_2\}$，将上述方法应用于 $-\eta=\min\{-x_1,-x_2\}$；对于 $n$ 项的 $\eta=\min_ix_i$，每一项配一个 0–1 变量，并令 $\sum_i\delta_i=1$（或 $\ge1$）。

(sec-07-02-abs)=
### 绝对值

根据上述规则：

* $\abs{f(\vx)}\le c$（$f$ 为仿射函数，$c\ge0$）即 $-c\le f(\vx)\le c$，是两个线性约束行；
* $\abs{f(\vx)}\ge c$ 是析取条件“$f(\vx)\ge c$ 或 $f(\vx)\le-c$”，需要一个 0–1 变量：$f(\vx)\ge c-M\delta$ 和 $f(\vx)\le-c+M(1-\delta)$；
* 对常数 $c\ge0$，$\abs{f(\vx)}=c$ 表示 $f(\vx)=c$ 或 $f(\vx)=-c$，一个 0–1 变量即可精确表达：$f(\vx)=c\,(2\delta-1)$。*（把 $\abs{f(\vx)}=c$ 写成 $f(\vx)=\max\{-c,c\}$ 是错误的：这个最大值就是 $c$。）* 若 $c$ 本身是有界变量，则乘积 $\delta c$ 用 {eq}`eq-07-02-bin-cont` 线性化。

在目标函数中，考虑

$$
\begin{aligned}
\min_{\vx}\ & \sum_{i=1}^nc_i\abs{x_i}\qquad (c_i\ge0)\\
\st\ & \mA\vx\ge\vb .
\end{aligned}
$$ (eq-07-02-abs-obj)

它有两种等价的线性规划形式。**上图形式**（epigraph form）：通过 $\eta_i\ge x_i$ 和 $\eta_i\ge-x_i$ 引入 $\eta_i\ge\abs{x_i}$，并极小化 $\sum_ic_i\eta_i$。**拆分形式**（split form）：令 $\vx=\vx^+-\vx^-$，$\vx^+,\vx^-\ge\vzero$，在 $\mA\vx^+-\mA\vx^-\ge\vb$ 的约束下极小化 $\sum_ic_i(x_i^++x_i^-)$。在拆分形式中，没有什么能禁止 $x_i^+>0$ 与 $x_i^->0$ 同时成立，此时 $x_i^++x_i^->\abs{x_i}$。但当 $c_i>0$ 时这样的解绝不是最优的：把两部分同时减去 $\min\{x_i^+,x_i^-\}$，$\vx$ 和约束都保持不变，而费用减小了。因此在最优解处 $x_i^+x_i^-=0$，$x_i^++x_i^-=\abs{x_i}$（Bertsimas & Tsitsiklis, 1997, §1.3）。

这两种论证都需要 $c_i\ge0$。若某个 $c_i<0$，模型就是在*极大化* $\abs{x_i}$，这是非凸的方向：在上图形式中 $\eta_i$ 可以无限增大，在拆分形式中 $x_i^+$ 和 $x_i^-$ 一起增大，无论哪种线性规划都会无界或给出错误结果。对每个这样的变量配一个 0–1 变量即可恢复其含义：$x_i^+\le U_i\delta_i$，$x_i^-\le U_i(1-\delta_i)$，其中 $U_i$ 是 $\abs{x_i}$ 的一个界。

下例中的三个距离函数说明了这条规则（见其解答中的{numref}`fig-07-02-abs-max`）：绝对值的加权和与绝对值的最大值是凸函数（极小化它们是线性规划问题），而绝对值的最小值不是凹函数（极大化它需要 0–1 变量）。

:::{prf:example} 走廊上的设施选址
:label: ex-07-02-corridor

五个城镇沿一条 $50$ 千米长的公路走廊分布，位置为 $a_i=0,12,20,34,50$ 千米，人口为 $w_i=8,5,12,4,6$ 千人。要在位置 $x\in[0,50]$ 处设置一个设施。分别求下列设施的最佳位置：(a) 使居民总距离 $\sum_iw_i\abs{x-a_i}$ 最小的养护站；(b) 使到最远城镇的距离最小的应急站；(c) 使到最近城镇的距离最大的垃圾填埋场。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

```{figure} figures/07-02-abs-max.svg
:name: fig-07-02-abs-max
:width: 680px
:class: dark-light
:alt: 走廊上选址位置的三个分段线性函数：加权距离和、最大距离与最小距离

{prf:ref}`ex-07-02-corridor` 中五个城镇的距离函数（坐标轴上标出了城镇位置）。左图：加权总距离 $\sum_iw_i\abs{x-a_i}$，凸函数，在 $x=20$ 处取最小值。中图：最大距离 $\max_i\abs{x-a_i}$，凸函数，在 $x=25$ 处取最小值。右图：最小距离 $\min_i\abs{x-a_i}$，既非凸也非凹，在 $x=42$ 处取最大值。
```

**(a) 养护站。** 令 $\eta_i\ge\abs{x-a_i}$，

$$
\begin{aligned}
\min\ & \sum_{i=1}^5w_i\eta_i\\
\st\ & \eta_i\ge x-a_i, \qquad i=1,\dots,5,\\
& \eta_i\ge a_i-x, \qquad i=1,\dots,5,\\
& 0\le x\le50 .
\end{aligned}
$$

线性规划给出 $x^*=20$ 千米，目标函数值为 $436$ 千人·千米。验证：$\sum_iw_i\abs{x-a_i}$ 的斜率等于（左侧人口）减去（右侧人口）；在 $20$ 的紧左侧它为 $8+5-(12+4+6)=-9$，在紧右侧为 $8+5+12-(4+6)=+15$，因此 $x=20$ 是极小点，即一个**加权中位数**（weighted median）（累计人口 $8,13,25$ 在第三个城镇处首次超过 $35$ 的一半）。该处的费用为 $8\cdot20+5\cdot8+0+4\cdot14+6\cdot30=436$；在 $x=12$ 处为 $508$，在 $x=34$ 处为 $646$。

**(b) 应急站。** 在对所有 $i$ 都有 $z\ge x-a_i$ 和 $z\ge a_i-x$ 的约束下极小化 $z$（{eq}`eq-07-02-minmax`）。只有两端的城镇起作用，线性规划给出 $x^*=25$ 千米，$z^*=25$ 千米，即 $0$ 与 $50$ 的中点。

**(c) 垃圾填埋场。** 极大化 $\min_i\abs{x-a_i}$ 要求对每个 $i$ 有 $z\le\abs{x-a_i}$，而其中每一个都是析取条件（“城镇位于设施左侧至少 $z$ 处，或位于其右侧至少 $z$ 处”）。令 $\delta_i=1$ 表示“设施位于城镇 $i$ 的左侧”，取 $M=100$ 千米（走廊长度的两倍可以界定两侧），

$$
\begin{aligned}
\max\ & z\\
\st\ & z\le x-a_i+M\delta_i, \qquad i=1,\dots,5,\\
& z\le a_i-x+M(1-\delta_i), \qquad i=1,\dots,5,\\
& 0\le x\le50,\\
& \delta_i\in\{0,1\}, \qquad i=1,\dots,5 .
\end{aligned}
$$

MILP 给出 $x^*=42$ 千米，$z^*=8$ 千米，$\delta_5=1$，其余 0–1 变量为 $0$：填埋场位于 $34$ 千米与 $50$ 千米处两个城镇的正中间，那里是最大的空隙。（各空隙为 $12$、$8$、$14$ 和 $16$ 千米；走廊的端点 $0$ 和 $50$ 处都有城镇，因此没有更好的端点。）该模型的线性规划松弛毫无用处：它报告 $z=50$，因为分数值的 $\delta_i$ 松弛了每一行。
:::

(sec-07-02-fractional)=
## 线性分式目标函数

**线性分式规划**（linear-fractional program）优化的是仿射函数之比，例如每美元效益或每位乘客的成本：

$$
\begin{aligned}
\max_{\vx}\ & \frac{\vp^\T\vx+p_0}{\vq^\T\vx+q_0}\\
\st\ & \mA\vx\le\vb,\\
& \vx\ge\vzero,
\end{aligned}
$$ (eq-07-02-fractional)

其中分母在可行集 $\Omega=\{\vx:\mA\vx\le\vb,\vx\ge\vzero\}$ 上为正，并假定 $\Omega$ 非空且有界。目标函数既非凸也非凹，但它是拟线性的，**Charnes 和 Cooper (1962)** 证明了通过变量代换可以把问题化为线性规划。令

$$
\begin{aligned}
\sigma&=\frac1{\vq^\T\vx+q_0}>0,\\
\vy&=\sigma\vx .
\end{aligned}
$$

把分子、分母和约束都乘以 $\sigma$，得到

$$
\begin{aligned}
\max_{\vy,\sigma}\ & \vp^\T\vy+p_0\,\sigma\\
\st\ & \mA\vy-\vb\,\sigma\le\vzero,\\
& \vq^\T\vy+q_0\,\sigma=1,\\
& \vy\ge\vzero,\\
& \sigma\ge0 .
\end{aligned}
$$ (eq-07-02-charnes-cooper)

:::{prf:theorem} Charnes–Cooper 变换
:label: thm-07-02-charnes-cooper

若 $\Omega$ 非空有界，且在 $\Omega$ 上 $\vq^\T\vx+q_0>0$，则 {eq}`eq-07-02-fractional` 与 {eq}`eq-07-02-charnes-cooper` 的最优值相同，并且对 {eq}`eq-07-02-charnes-cooper` 的任一最优解 $(\vy^*,\sigma^*)$，$\vx^*=\vy^*/\sigma^*$ 是 {eq}`eq-07-02-fractional` 的最优解。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

对 $\vx\in\Omega$，数对 $\sigma=1/(\vq^\T\vx+q_0)$、$\vy=\sigma\vx$ 是 {eq}`eq-07-02-charnes-cooper` 的可行解，其目标函数值等于 $\vx$ 处的比值。反之，设 $(\vy,\sigma)$ 可行。若 $\sigma=0$，则 $\mA\vy\le\vzero$，$\vy\ge\vzero$ 且 $\vq^\T\vy=1$，故 $\vy\ne\vzero$，并且对每个 $\vx\in\Omega$ 和 $\lambda\ge0$ 有 $\vx+\lambda\vy\in\Omega$，与有界性矛盾。因此 $\sigma>0$，而 $\vx=\vy/\sigma$ 满足 $\mA\vx\le\vb$、$\vx\ge\vzero$ 和 $\vq^\T\vx+q_0=1/\sigma$，故其比值等于 $\vp^\T\vy+p_0\sigma$。于是两个问题取得的函数值集合相同，最优解一一对应。
:::

对极小化问题，同样的变换也适用。若分母在 $\Omega$ 上为负，则先把分子和分母都乘以 $-1$。

:::{prf:example} 最佳效益–费用比
:label: ex-07-02-fractional

某机构把养护班组分配到路面重铺（$x_1$ 个班组·周）和排水修复（$x_2$ 个班组·周）两项工作。效益为 $5x_1+3x_2$，费用为 $2x_1+2x_2+4$（单位为百万美元，含 400 万美元的固定调遣费用）。可用的班组·周至多为 $12$，路面重铺至多能吸纳 $8$，排水修复至少需要 $2$。极大化效益–费用比。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*变换。* 令 $\sigma=1/(2x_1+2x_2+4)$，$\vy=\sigma\vx$，则线性规划 {eq}`eq-07-02-charnes-cooper` 为

$$
\begin{aligned}
\max\ & 5y_1+3y_2\\
\st\ & y_1+y_2-12\sigma\le0,\\
& y_1-8\sigma\le0,\\
& -y_2+2\sigma\le0,\\
& 2y_1+2y_2+4\sigma=1,\\
& y_1,y_2,\sigma\ge0 .
\end{aligned}
$$

*求解。* 线性规划的最优解为 $\sigma^*=\tfrac1{24}$，$\vy^*=(\tfrac13,\tfrac1{12})$，最优值为 $\tfrac{23}{12}\approx1.917$，于是 $\vx^*=\vy^*/\sigma^*=(8,2)$ 班组·周。验证：$(8,2)$ 处的比值为 $(40+6)/(16+4+4)=46/24=\tfrac{23}{12}$。$\Omega$ 的各顶点处的比值为 $(0,2)$：$\tfrac34$；$(0,12)$：$\tfrac97$；$(8,2)$：$\tfrac{23}{12}$；$(8,4)$：$\tfrac{13}7\approx1.857$，因此 $(8,2)$ 确实最好。注意，极大化*净*效益 $(5x_1+3x_2)-(2x_1+2x_2+4)$ 会选择 $(8,4)$，净效益为 2 400 万美元：比值与差值是不同的准则，比值偏好规模较小的方案。
:::

(sec-07-02-toolbox)=
## 各种技术一览

{numref}`tab-07-02-toolbox`汇总了本章的各种技术。“精确”是指线性模型与非线性模型的解相同。

```{list-table} 线性化技术
:name: tab-07-02-toolbox
:header-rows: 1
:widths: 30 40 30

* - 项
  - 线性模型
  - 是否精确？0–1 变量？
* - $\delta_1\delta_2$，0–1 变量
  - {eq}`eq-07-02-bin-bin`
  - 精确；不增加
* - $\delta x$, $x\in[\underline x,\overline x]$
  - {eq}`eq-07-02-bin-cont`
  - 精确；不增加
* - $x\cdot y$，$x$ 为整数，$y$ 有界
  - 二进制展开 {eq}`eq-07-02-binary-expansion` + {eq}`eq-07-02-bin-cont`
  - 精确；$\lceil\log_2(U+1)\rceil$ 个
* - $xy$，两者都是箱体上的连续变量
  - McCormick {eq}`eq-07-02-mccormick`
  - 松弛；仅在分割时需要 0–1 变量
* - $\max_k(\cdot)\le\eta$, $\abs{f}\le c$, $\min\max$, $\max\min$
  - 每项一行 {eq}`eq-07-02-max-le`, {eq}`eq-07-02-minmax`
  - 精确；不需要
* - $\max_k(\cdot)\ge\eta$, $\abs{f}\ge c$, $\max\max$, $\min\min$, $\eta=\min$
  - 析取条件 {eq}`eq-07-02-max-ge`, {eq}`eq-07-02-eta-min`
  - 精确；每项一个（两项时只需一个）
* - $\abs{f}=c$
  - $f=c(2\delta-1)$
  - 精确；一个
* - 极小化 $\sum c_i\abs{x_i}$，$c_i\ge0$
  - 上图形式或拆分形式
  - 精确；不需要
* - 仿射函数之比
  - Charnes–Cooper {eq}`eq-07-02-charnes-cooper`
  - 精确；不需要
```

## 本章回顾

### 小结

- 0–1 变量之积可用 {eq}`eq-07-02-bin-bin` 精确线性化；当乘积只出现在目标函数中时，其系数的符号决定需要哪一半约束。
- 0–1 变量与有界连续变量之积可用四个不等式 {eq}`eq-07-02-bin-cont` 精确线性化。
- McCormick 不等式 {eq}`eq-07-02-mccormick` 给出 $xy$ 在箱体上的凸包络和凹包络；它们在箱体边界上是精确的，间隙至多为 $(\overline x-\underline x)(\overline y-\underline y)/4$，分割可使间隙缩小（{numref}`fig-07-02-mccormick`）。
- 从上方界定凸分段线性函数（或从下方界定凹分段线性函数）是线性的（{eq}`eq-07-02-max-le`、{eq}`eq-07-02-minmax`）；相反的方向是析取条件，需要 0–1 变量（{eq}`eq-07-02-max-ge`、{eq}`eq-07-02-eta-min`、{numref}`fig-07-02-abs-max`）。
- 在被极小化且权重非负的目标函数中，绝对值用上图形式或拆分形式表达都是线性的；$\abs{f}=c$ 即 $f=c(2\delta-1)$。
- 仿射函数之比通过 Charnes–Cooper 变换 {eq}`eq-07-02-charnes-cooper` 化为线性规划（{prf:ref}`thm-07-02-charnes-cooper`）。

### 常见错误

- **去掉了乘积线性化中错误的一半。** 单侧形式只对其所针对的目标函数系数符号有效（{prf:ref}`ex-07-02-synergy`）。
- **在 0–1 变量与连续变量之积中遗漏 $\omega\ge0$**（或 $\omega\ge\underline x\delta$）。
- **把 McCormick 包络当作精确的。** 它只在箱体边界上是精确的；在内部，线性规划的最优值会高估乘积（{prf:ref}`ex-07-02-mccormick`）。
- **在错误的一侧使用上图技巧。** $\min\max$ 和 $\max\min$ 是线性规划问题；$\max\max$、$\min\min$、$\max\abs{\cdot}$ 以及约束 $\abs{f}\ge c$ 都需要 0–1 变量。
- **费用为负时仍拆分 $x=x^+-x^-$。** “两部分中有一部分为零”的论证依赖于 $c_i>0$。
- **把 Charnes–Cooper 变换中的 $\vy$ 当作 $\vx$。** 原问题的解是 $\vx=\vy/\sigma$。

## 参考文献

- Al-Khayyal, F. A., & Falk, J. E. (1983). Jointly constrained biconvex
  programming. *Mathematics of Operations Research, 8*(2), 273–286.
  <https://doi.org/10.1287/moor.8.2.273>
- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Section 11.4:
  linear fractional programming.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Section 1.3: piecewise linear convex
  objective functions and problems involving absolute values.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Sections
  3.1.7 and 3.2.3: epigraphs and pointwise maxima; Section 4.3.2:
  linear-fractional programming.)
- Charnes, A., & Cooper, W. W. (1962). Programming with linear fractional
  functionals. *Naval Research Logistics Quarterly, 9*(3–4), 181–186.
- Glover, F. (1975). Improved linear integer programming formulations of
  nonlinear integer problems. *Management Science, 22*(4), 455–460.
  <https://doi.org/10.1287/mnsc.22.4.455>
- Glover, F., & Woolsey, E. (1974). Converting the 0–1 polynomial programming
  problem to a 0–1 linear program. *Operations Research, 22*(1), 180–182.
  <https://doi.org/10.1287/opre.22.1.180>
- McCormick, G. P. (1976). Computability of global solutions to factorable
  nonconvex programs: Part I — Convex underestimating problems.
  *Mathematical Programming, 10*(1), 147–175.
  <https://doi.org/10.1007/BF01580665>
- Williams, H. P. (2013). *Model building in mathematical programming*
  (5th ed.). Wiley. (Chapters 9–10: logical conditions, products of
  variables and alternative formulations of integer models.)
