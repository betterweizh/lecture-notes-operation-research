# 逻辑约束与大 M 法

## 概述

线性规划（linear program, LP）描述的是可以连续调节规模的决策：生产 3.4 吨、运送 12.5 车的货物、建造 0.3 座桥。许多规划问题并非如此。一座桥要么建，要么不建；一座工厂一旦开设，就要全额支付其固定费用；一台水泵要么完全不运行，要么以不低于其最低流量的速率运行；两条技术规则中*至少有一条*必须成立。本章说明如何用**0–1 变量**（binary variable）$\delta_i\in\{0,1\}$ 和线性不等式表达这类“是/否”决策与逻辑条件，从而把模型变为**混合整数线性规划**（mixed-integer linear program, MILP）。核心工具是**大 $M$ 法**（big-$M$ method）：0–1 变量通过在约束右端项上加一个常数 $M$ 来“打开”或“关闭”该约束。我们将推导 $M$ 至少要取多大（以及为什么不应取得更大），然后在此基础上建立标准的逻辑建模模式：择一约束、若–则约束、“$n$ 个约束中有 $k$ 个成立”的条件、固定费用、半连续变量以及 0–1 变量之间的逻辑关系。择一约束的几何意义解释了为什么某些正确的模型形式优于另一些。本章最后给出两个基础设施模型：一个带先后次序规则的项目组合问题和一个设施开设问题。本章以线性规划建模（{doc}`03-01-lp-model-building`）、多面体几何（{doc}`03-02-lp-geometry`）和凸包（{prf:ref}`def-02-01-convex-hull`）为基础。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 解释 0–1 变量和常数 $M$ 如何打开或关闭一个不等式，并能由变量的界计算出有效且紧的 $M$；
- 解释为什么不必要地大的 $M$ 会削弱线性规划松弛并引起数值困难；
- 用 0–1 变量表达择一条件、若–则条件、“$m$ 个约束中至少/至多/恰有 $k$ 个成立”的条件、固定费用条件和半连续条件；
- 把关于“是/否”决策的逻辑陈述转化为线性不等式，并通过枚举加以验证；
- 以凸包为标尺，按线性规划松弛的强弱比较同一条件的不同模型形式；
- 为带有逻辑附加条件的项目选择问题和设施开设问题建立模型。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\delta_i\in\{0,1\}$ | 0–1 变量（指示变量） |
| $\Omega_0$ | 模型中所有*其他*约束所允许的点集（通常是一个箱体 $\vl\le\vx\le\vu$ 或有界多面体） |
| $g_i(\vx)\le 0$ | 由 0–1 变量打开或关闭的约束 |
| $M$, $M_i$ | 大 $M$ 常数：$g_i(\vx)$ 在 $\Omega_0$ 上的上界 |
| $\underline g$, $\underline g_i$ | $g(\vx)$、$g_i(\vx)$ 在 $\Omega_0$ 上的下界 |
| $\varepsilon$ | 用于表达严格不等式的小正容差 |
| $l_j$, $u_j$ | 变量 $x_j$ 的下界和上界 |
| $F$, $F_i$ | 固定费用（启动费用或开设费用），以货币为单位 |
| $\conv S$ | 集合 $S$ 的凸包 |
| $\omega$ | 定义为其他 0–1 变量的逻辑函数（与、或）的 0–1 变量 |

(sec-07-01-milp)=
## 0–1 变量与混合整数规划

**混合整数线性规划**的目标函数和约束都是线性的，并且部分变量限定取整数值：

$$
\begin{aligned}
\min_{\vx,\vdelta}\ & \vc^\T\vx+\vq^\T\vdelta \\
\st\ & \mA\vx+\mD\vdelta\le\vb, \\
& \vx\ge\vzero, \\
& \vdelta\in\{0,1\}^p .
\end{aligned}
$$ (eq-07-01-milp)

其中 $\vdelta=(\delta_1,\dots,\delta_p)$ 汇集了全部 0–1 变量，$\vq$ 是它们的目标函数系数。其可行集 $\Omega$ 是有限个多面体的并集（$\vdelta$ 的每个取值对应一个多面体），一般不是凸集。去掉整数性要求，即把 $\delta_i\in\{0,1\}$ 换成 $0\le\delta_i\le1$，就得到**线性规划松弛**（LP relaxation）。由于其可行集包含 $\Omega$，极小化问题的线性规划松弛的最优值不超过 MILP 的最优值：它是一个**下界**（对极大化问题则是上界）。分支定界法（branch and bound）求解器在搜索树的每个节点上都计算这样的界，因此松弛越接近 MILP，搜索就越快。整数解相同的两种模型形式的表现可能差别很大，本章的很大一部分内容就是讨论如何写出更强的那一种。

0–1 变量有两种自然的含义：表示一个**决策**（“建桥”，$\delta=1$），或者作为某个约束是否被强制执行的**指示变量**（indicator）。大 $M$ 法把这两种含义联系起来。

(sec-07-01-big-m)=
## 大 M 法

设 $\Omega_0$ 是模型中除某一个约束 $g(\vx)\le0$ 之外的所有约束所允许的点集，并设 $g$ 在 $\Omega_0$ 上有已知的上界 $M$：

$$
g(\vx)\le M\quad\text{对所有 }\vx\in\Omega_0 .
$$ (eq-07-01-valid-m)

于是约束 $g(\vx)\le M$ 在 $\Omega_0$ 上是**冗余的**：它不切掉任何点。引入一个 0–1 变量 $\delta$，把 $g(\vx)\le0$ 替换为

$$
\begin{gathered}
g(\vx)\le M(1-\delta),\\
\delta\in\{0,1\},\quad \vx\in\Omega_0 .
\end{gathered}
$$ (eq-07-01-bigm-relax)

当 $\delta=1$ 时右端项为 $0$，约束被强制执行；当 $\delta=0$ 时约束变为 $g(\vx)\le M$，$\Omega_0$ 中的每一点都满足它，因此约束被**松弛**（它可能成立，也可能不成立）。用文字表述就是 $\delta=1\Rightarrow g(\vx)\le0$。这个蕴含关系只有一个方向：$\delta=0$ 并不强制 $g(\vx)>0$。

方向相反的约束 $g(\vx)\ge0$ 先改写为 $-g(\vx)\le0$：

$$
\begin{aligned}
-g(\vx)&\le M(1-\delta),\\
M&\ge\max_{\vx\in\Omega_0}\bigl(-g(\vx)\bigr) .
\end{aligned}
$$ (eq-07-01-bigm-ge)

**等式**约束 $g(\vx)=0$ 不能只用一个大 $M$ 项来关闭，因为右端项上的任何常数都不能把等式变成冗余约束。应把它拆成两个不等式 $g(\vx)\le0$ 和 $-g(\vx)\le0$，用*同一个* 0–1 变量、各自的常数分别松弛：$g(\vx)\le M^+(1-\delta)$ 和 $-g(\vx)\le M^-(1-\delta)$，其中 $M^+$ 和 $M^-$ 分别是 $g$ 和 $-g$ 在 $\Omega_0$ 上的上界。

(sec-07-01-choosing-m)=
### M 的选取

条件 {eq}`eq-07-01-valid-m` 说明了什么是**有效的** $M$：任何不小于 $\max_{\vx\in\Omega_0}g(\vx)$ 的数。取值过小是一个错误，而不仅仅是一个弱点：若对某个 $\bar\vx\in\Omega_0$ 有 $M<g(\bar\vx)$，那么即使约束本应被松弛，点 $\bar\vx$ 也会被切掉，模型可能因此丢失真正的最优解。

对于箱体 $\vl\le\vx\le\vu$ 上的线性约束 $g(\vx)=\va^\T\vx-b$，最大值可以逐项求得：当 $a_j>0$ 时 $a_jx_j$ 在 $u_j$ 处最大，当 $a_j<0$ 时在 $l_j$ 处最大。因此最小的有效常数为

$$
M=\sum_{j:\,a_j>0}a_ju_j+\sum_{j:\,a_j<0}a_jl_j-b .
$$ (eq-07-01-m-from-bounds)

当 $\Omega_0$ 是一般的多面体时，$\max_{\Omega_0}g$ 本身就是一个线性规划问题；箱体界 {eq}`eq-07-01-m-from-bounds` 仍然有效，但可能偏大。若某个变量没有有限的界，则不存在有效的 $M$。这时要么根据问题的物理背景推出界（任何工厂的发货量都不超过总需求，任何流量都不超过管道的容量），要么把条件作为**指示约束**（indicator constraint）交给求解器：目前大多数 MILP 求解器都能直接接受“$\delta=1\Rightarrow\va^\T\vx\le b$”这样的表述，并通过分支来保证它成立，无需任何 $M$。

为什么不干脆取 $M=10^6$ 以求保险？原因有二。

1. **松弛弱。** 在线性规划松弛中 $\delta$ 可以取分数值。此时任何满足 $\delta\le1-g(\vx)/M$ 的 $\vx$ 都满足约束 {eq}`eq-07-01-bigm-relax`：$M$ 越大，在约束被违反的情况下 $\delta$ 仍可越接近 $1$；反过来说，为关闭该约束而需“付出”的 $\delta$ 的份额就越小。$M$ 很大时，松弛问题的表现几乎就像逻辑条件不存在一样，其界很差，分支定界法必须搜索多得多的节点。{prf:ref}`ex-07-01-either-or`和{ref}`sec-07-01-fixed-charge`中的固定费用模型用数值说明了这一点。
2. **数值误差。** 求解器把与某个整数相差在容差（通常为 $10^{-6}$）以内的值视为整数。取 $M=10^{6}$ 时，值 $\delta=1-10^{-6}$ 被当作 $1$，但此时约束 $g(\vx)\le M(1-\delta)=1$ 却允许整整一个单位的违反。大系数还会使线性规划的基矩阵变得病态。

因此规则是：**取 $M$ 为最小的有效界**，由 {eq}`eq-07-01-m-from-bounds` 或更紧的论证算出，并且对每个约束使用各自的常数 $M_i$。

(sec-07-01-reverse)=
### 强制 0–1 变量取值：反向蕴含

模型形式 {eq}`eq-07-01-bigm-relax` 表达的是 $\delta=1\Rightarrow g(\vx)\le0$。有时需要的是其逆命题：0–1 变量必须*记录*约束是否成立，即 $g(\vx)\le0\Rightarrow\delta=1$，或等价地 $\delta=0\Rightarrow g(\vx)>0$。线性规划中不能施加严格不等式，因此要使用容差 $\varepsilon>0$：$\delta=0\Rightarrow g(\vx)\ge\varepsilon$。设 $g$ 在 $\Omega_0$ 上有下界 $\underline g\le g(\vx)$（$\underline g<0$），则该条件为

$$
g(\vx)\ge\varepsilon-(\varepsilon-\underline g)\,\delta .
$$ (eq-07-01-reverse)

当 $\delta=0$ 时它给出 $g(\vx)\ge\varepsilon$；当 $\delta=1$ 时它给出 $g(\vx)\ge\underline g$，这总是成立的。{eq}`eq-07-01-bigm-relax` 与 {eq}`eq-07-01-reverse` 合在一起，使得 $\delta=1$ 当且仅当 $g(\vx)\le0$，只是满足 $0<g(\vx)<\varepsilon$ 的点被排除在外。若 $g$ 在 $\Omega_0$ 上只取整数值，则取 $\varepsilon=1$ 是精确的；否则 $\varepsilon$ 必须大于求解器的可行性容差。

(sec-07-01-either-or)=
## 择一约束

若两个约束 $g_1(\vx)\le0$ 和 $g_2(\vx)\le0$ 中至少有一个必须成立，则称它们互为**备选**约束。典型的情形有：一台机器要么先加工工件 $A$ 再加工工件 $B$，要么先加工 $B$ 再加工 $A$；一座建筑要么建在道路以西，要么建在河流以北。为每个约束配一个 0–1 变量，

$$
\begin{aligned}
g_1(\vx)&\le M_1(1-\delta_1),\\
g_2(\vx)&\le M_2(1-\delta_2),\\
\delta_1+\delta_2&\ge1,\qquad \delta_1,\delta_2\in\{0,1\} .
\end{aligned}
$$

于是至少有一个约束被强制执行。由于从不需要同时强制执行两个约束，也可以改用 $\delta_1+\delta_2=1$，进而消去 $\delta_2=1-\delta_1$：

$$
\begin{aligned}
g_1(\vx)&\le M_1(1-\delta),\\
g_2(\vx)&\le M_2\,\delta,\qquad \delta\in\{0,1\} .
\end{aligned}
$$ (eq-07-01-either-or)

$\delta=1$ 时第一个约束成立，$\delta=0$ 时第二个约束成立。若某个备选约束是由若干不等式组成的方程组（一个多面体），则它的每一行都使用同一个 0–1 变量，每一行各有自己的常数。

从几何上看，择一条件的可行集是两个集合的**并集**，通常不是凸集。任何线性规划松弛都是凸的，因此它最好也只能是该并集的**凸包**；它与凸包的接近程度取决于所取的常数。

:::{prf:example} 三角形或正方形
:label: ex-07-01-either-or

在箱体 $\Omega_0=[0,4]^2$ 上，点 $\vx=(x_1,x_2)$ 必须位于三角形 $T=\{x_1+x_2\le2\}$ 中或正方形 $S=\{x_1\ge3,\ x_2\ge3\}$ 中。极大化 $-x_1+2x_2$。用大 $M$ 约束表达该条件，求出整数最优解，并比较采用最小有效常数与对每一行都取 $M=10$ 时线性规划松弛给出的界。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*建模。* 令 $\delta=1$ 表示选择 $T$，$\delta=0$ 表示选择 $S$。三角形的约束行 $x_1+x_2-2\le0$ 在 $\delta=0$ 时被松弛；正方形的约束行 $3-x_1\le0$ 和 $3-x_2\le0$ 在 $\delta=1$ 时被松弛：

$$
\begin{aligned}
x_1+x_2-2&\le M_T(1-\delta),\\
3-x_1&\le M_1\delta,\\
3-x_2&\le M_2\delta .
\end{aligned}
$$

由 {eq}`eq-07-01-m-from-bounds`，最小的有效常数为 $M_T=4+4-2=6$ 和 $M_1=M_2=3-0=3$。

*整数解。* $\delta=1$ 时，$T$ 中的最好点是 $(0,2)$，目标函数值为 $4$；$\delta=0$ 时，$S$ 中的最好点是 $(3,4)$，目标函数值为 $5$。最优解为 $\vx^*=(3,4)$，$\delta^*=0$，最优值为 $5$。

*线性规划松弛。* {numref}`fig-07-01-either-or`给出了各线性规划松弛在 $(x_1,x_2)$ 平面上的投影。$T\cup S$ 的凸包是顶点为 $(0,0)$、$(2,0)$、$(4,3)$、$(4,4)$、$(3,4)$、$(0,2)$ 的六边形，面积为 $10$。采用紧常数时的松弛是更大的六边形 $(0,0)$、$(2,0)$、$(4,2)$、$(4,4)$、$(2,4)$、$(0,2)$，面积为 $12$，其线性规划最优值为 $6$，在 $\vx=(2,4)$、$\delta=\tfrac13$ 处取得：该点满足 $x_1+x_2=6\le2+6\cdot\tfrac23$ 和 $3-x_1=1\le3\cdot\tfrac13$，因此三分之一个“三角形决策”就让它同时逃出了两个集合。取 $M_T=M_1=M_2=10$ 时，松弛是整个箱体（面积 $16$）：线性规划最优值为 $8$，在 $(0,4)$ 处取得，例如取 $\delta=0.3$。两个 MILP 给出相同的整数最优值 $5$，不同的只是界：$6$（间隙 $20\%$）与 $8$（间隙 $60\%$）。

```{figure} figures/07-01-either-or.svg
:name: fig-07-01-either-or
:width: 620px
:class: dark-light
:alt: 箱体中一个三角形与一个正方形的并集、其凸包以及两种大 M 松弛的投影

{prf:ref}`ex-07-01-either-or` 的择一约束。左图：可行集 $T\cup S$（深色）及其凸包（虚线轮廓）。右图：采用紧常数 $M_T=6$、$M_1=M_2=3$ 时线性规划松弛的投影（阴影六边形，面积 12）和取 $M=10$ 时的投影（整个箱体，面积 16），并重复画出 $T\cup S$ 及其凸包以便比较。图中标出了 $\max\,-x_1+2x_2$ 的最优解 $\vx^*=(3,4)$ 和松弛问题的最优解 $(2,4)$。
```
:::

由{prf:ref}`ex-07-01-either-or`的解答可见，即使采用紧常数，松弛也严格大于 $T\cup S$ 的凸包。凸包本身可以通过一个规模更大的模型形式得到。令 $\vx=\vy+\vz$，要求 $\vy\in\delta\,T$ 和 $\vz\in(1-\delta)S$，即把每个集合的右端项（以及变量的界）乘以各自的 0–1 变量：$y_1+y_2\le2\delta$，$\vzero\le\vy\le4\delta\vone$，$z_1\ge3(1-\delta)$，$z_2\ge3(1-\delta)$，$\vz\le4(1-\delta)\vone$。当 $\delta\in\{0,1\}$ 时，两个副本中一个为零，另一个位于其对应的集合中；当 $\delta$ 取分数值时，点 $\vx$ 是 $T$ 中一点与 $S$ 中一点的凸组合。析取规划的这种**分解形式**（disaggregated formulation，又称凸包形式）（Balas, 1985）的松弛恰好是凸包：在{prf:ref}`ex-07-01-either-or`中它的线性规划最优值为 $5$，等于整数最优值。它需要为每个备选项复制一份 $\vx$，因此在界比模型规模更重要时使用（Williams, 2013, 第 9–10 章）。

(sec-07-01-if-then)=
## 若–则约束

**若–则**（if–then）条件的含义是：*若* $g_1(\vx)>0$，*则* $g_2(\vx)\le0$。当 $g_1(\vx)\le0$ 时，对 $g_2$ 没有任何要求。在逻辑上，$A\Rightarrow B$ 等价于 $(\text{非 }A)\text{ 或 }B$，因此该条件就是择一约束

$$
g_1(\vx)\le0\quad\text{或}\quad g_2(\vx)\le0,
$$

由 {eq}`eq-07-01-either-or` 得

$$
\begin{aligned}
g_1(\vx)&\le M_1(1-\delta),\\
g_2(\vx)&\le M_2\,\delta,\qquad \delta\in\{0,1\} .
\end{aligned}
$$ (eq-07-01-if-then)

验证其逻辑：若 $g_1(\vx)>0$，第一个不等式排除了 $\delta=1$，于是 $\delta=0$，第二个不等式强制 $g_2(\vx)\le0$。若 $g_1(\vx)\le0$，取 $\delta=1$ 可满足第一个不等式并松弛第二个不等式，于是 $g_2$ 不受限制。用集合的语言，该条件描述的是

$$
\begin{aligned}
&\{\vx: g_1(\vx)>0,\ g_2(\vx)\le0\}\cup\{\vx: g_1(\vx)\le0\}\\
&\quad=\{\vx: g_2(\vx)\le0\}\cup\{\vx: g_1(\vx)\le0\},
\end{aligned}
$$

当 $g_1,g_2$ 为仿射函数时，它是两个半平面的并集。例如取 $g_1=x_1+x_2-1$ 和 $g_2=1-x_1+x_2$（“若 $x_1+x_2>1$，则 $x_1-x_2\ge1$”），可行集是整个平面去掉满足 $x_1+x_2>1$ 且 $x_1-x_2<1$ 的开楔形区域。

:::{prf:example} 0–1 决策之间的若–则条件
:label: ex-07-01-if-then-binaries

四个“是/否”决策 $\delta_1,\dots,\delta_4$ 之间有如下联系：若 $\delta_1=1$，则 $\delta_2=\delta_3=\delta_4=0$（选择方案 1 就排除方案 2–4）。先由一般的若–则规则 {eq}`eq-07-01-if-then`，再按每对互斥方案各写一行，把该条件写成线性约束，并比较两个约束组的线性规划松弛。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*由一般规则。* 前提 $\delta_1=1$ 即 $g_1=\delta_1>0$，结论为 $g_2=\delta_2+\delta_3+\delta_4\le0$。在 $\Omega_0=\{0,1\}^4$ 上，界为 $M_1=1$ 和 $M_2=3$，引入辅助 0–1 变量 $y$ 后 {eq}`eq-07-01-if-then` 成为 $\delta_1\le1-y$，$\delta_2+\delta_3+\delta_4\le3y$。消去 $y$（这两个不等式存在某个解 $y\in\{0,1\}$，当且仅当 $\delta_2+\delta_3+\delta_4\le3(1-\delta_1)$），得到**聚合**约束

$$
\delta_2+\delta_3+\delta_4\le3(1-\delta_1) .
$$

*（前提必须写成 $\delta_1>0$ 或 $\delta_1\ge1$；0–1 变量永远不会大于 $1$，因此写成前提“$\delta_1>1$”的条件永远不会被触发。）*

*分解形式。* 同一条件也可以表述为：对每个 $j$，“$\delta_1$ 与 $\delta_j$ 中至多有一个为 1”：

$$
\delta_1+\delta_j\le1,\qquad j=2,3,4 .
$$

枚举全部 $16$ 个 0–1 向量可以验证，两个约束组接受的恰好都是满足 $\delta_1=0$ 或 $\delta_2=\delta_3=\delta_4=0$ 的向量。但它们的松弛不同：分数点 $(\delta_1,\dots,\delta_4)=(\tfrac23,1,0,0)$ 满足聚合约束（$1\le3\cdot\tfrac13$），却违反 $\delta_1+\delta_2\le1$。而分解形式的松弛中的每一点都满足聚合约束（把三个不等式相加即得）。
:::

由{prf:ref}`ex-07-01-if-then-binaries`的解答可见，分解形式绝不弱于聚合形式，因而更为可取，尽管它用三行而不是一行。

(sec-07-01-k-of-n)=
## m 个约束中至少、至多和恰有 k 个成立

设 $g_i(\vx)\le0$，$i=1,\dots,m$，是一组约束，其有效常数为 $M_i$。要求其中**至少 $k$ 个**成立，是择一情形的直接推广：

$$
\begin{aligned}
g_i(\vx)&\le M_i(1-\delta_i),\quad i=1,\dots,m,\\
\sum_{i=1}^m\delta_i&\ge k,\\
\delta_i&\in\{0,1\},\quad i=1,\dots,m .
\end{aligned}
$$ (eq-07-01-k-of-m)

（择一约束就是 $k=1$、$m=2$ 的情形。）单个常数 $M=\max_iM_i$ 对每一行都有效，但由{ref}`sec-07-01-choosing-m`，各行分别取常数更好。

“**至多 $k$ 个**成立”和“**恰有 $k$ 个**成立”性质不同。在 {eq}`eq-07-01-k-of-m` 中把 $\ge k$ 换成 $\le k$ 或 $=k$，限制的只是被*强制执行*的约束个数；其余约束只是被松弛，而不是被违反，它们照样可能成立。因此 $\sum_i\delta_i\le k$ 根本不构成限制（允许 $\vdelta=\vzero$），而 $\sum_i\delta_i=k$ 的含义仍然是“至少 $k$ 个”。要统计被满足的约束个数，还必须通过反向蕴含 {eq}`eq-07-01-reverse`，在约束成立时强制相应的 0–1 变量取 $1$：

$$
\begin{aligned}
g_i(\vx)&\le M_i(1-\delta_i),\quad i=1,\dots,m,\\
g_i(\vx)&\ge\varepsilon-(\varepsilon-\underline g_i)\,\delta_i,\quad i=1,\dots,m,\\
\sum_{i=1}^m\delta_i&\le k\quad(\text{或 }=k),\\
\delta_i&\in\{0,1\},\quad i=1,\dots,m .
\end{aligned}
$$ (eq-07-01-at-most-k)

此时 $\delta_i=1$ 当且仅当 $g_i(\vx)\le0$（在容差 $\varepsilon$ 意义下），基数约束统计的就是被满足的约束个数。“至多”类条件在实践中很少见（它们要求模型有意违反规则）；例如，当至多允许 $k$ 项服务标准被突破时就会出现这类条件。

(sec-07-01-fixed-charge)=
## 固定费用与半连续变量

**固定费用**（fixed charge）型成本包括一个启动部分 $F>0$（只要使用某项活动就须支付）和每单位的可变部分 $c\ge0$：

$$
C(x)=\begin{cases}0, & x=0,\\ F+cx, & 0<x\le u .\end{cases}
$$

$C$ 在 $0$ 处不连续，在 $[0,u]$ 上是凹函数，因此不能出现在线性规划中。引入 0–1 变量 $\delta$（“使用该活动”）后，成本变为 $F\delta+cx$，而 $x$ 与 $\delta$ 之间的联系是大 $M$ 约束 $x\le M\delta$（若 $\delta=0$，则 $x\le0$）。最小的有效常数就是界 $u$ 本身，于是该约束与容量约束合而为一：

$$
\begin{aligned}
\min_{x,\delta}\ & F\delta+cx\\
\st\ & 0\le x\le u\,\delta,\\
& \delta\in\{0,1\} .
\end{aligned}
$$ (eq-07-01-fixed-charge)

没有约束在 $x=0$ 时强制 $\delta=0$：这由极小化完成，因为 $\delta=1$ 会白白花费 $F$。

在线性规划松弛中，$\delta$ 只需达到 $x/u$，因此生产 $x$ 的松弛成本为 $(F/u+c)\,x$：固定费用被摊到容量上，仿佛它是一项可变成本。{numref}`fig-07-01-fixed-charge`针对 $F=\$500$、每吨 $c=\$20$ 和 $u=100$ 吨的情形展示了这一点。生产 $x=60$ 吨的真实成本是 $500+20\cdot60=\$1700$；松弛问题只计 $(5+20)\cdot60=\$1500$。若采用不必要地大的常数 $M=1000$，松弛问题只计 $(0.5+20)\cdot60=\$1230$：下界丢掉了固定费用的 $94\%$，而不是 $40\%$。

```{figure} figures/07-01-fixed-charge.svg
:name: fig-07-01-fixed-charge
:width: 520px
:class: dark-light
:alt: 在零点处有跳跃的固定费用成本函数及其两种线性规划松弛对应的两条直线

$F=\$500$、$c=\$20$/t、容量 $u=100$ t 时的固定费用成本 $C(x)$（实线；实心圆点表示 $C(0)=0$，空心圆点表示 $x\downarrow0$ 时 $C$ 趋近但取不到的值 $F$）。$M=u$ 时的线性规划松弛用直线 $25x$ 代替它，这条直线是 $C$ 在 $[0,u]$ 上的凸包络（虚线）；$M=1000$ 时该直线为 $20.5x$（点线）。在 $x=60$ t 处，三者的值分别为 1700 美元、1500 美元和 1230 美元。
```

若变量要么为 $0$，要么介于一个正的下界与一个上界之间，即 $x\in\{0\}\cup[l,u]$，$0<l\le u$，则称它为**半连续变量**（semi-continuous variable）。一台水泵要么关闭，要么以 $l=20$ 至 $u=100$ 升/秒的流量供水；一家供应商要么不接受订单，要么接受不少于其最小订货量的订单。令 $\delta=1$ 表示“开启”，则

$$
l\,\delta\le x\le u\,\delta,\qquad \delta\in\{0,1\} .
$$ (eq-07-01-semicontinuous)

当 $\delta=0$ 时上下界均为 $0$；当 $\delta=1$ 时上下界分别为 $l$ 和 $u$。固定费用与最低水平常常同时出现，此时 {eq}`eq-07-01-semicontinuous` 同时表达了两者。有些求解器可以直接接受半连续变量，并对析取条件 $x=0$ 或 $x\ge l$ 进行分支，无需显式的 0–1 变量。

(sec-07-01-logic)=
## 0–1 变量之间的逻辑关系

一旦每个“是/否”陈述都有了对应的 0–1 变量，陈述之间的逻辑关系就成为这些 0–1 变量之间的线性不等式。设 $\delta_A,\delta_B,\delta_C$ 分别指示陈述 $A$、$B$、$C$（它们可以是决策，也可以是通过 {eq}`eq-07-01-bigm-relax` 使之成立的约束）。以下三条规则足以涵盖一切：

* 否定：“非 $A$”即 $1-\delta_A$；
* 若干文字的析取（或）成立，当且仅当它们之和至少为 $1$；
* 若干条件的合取（与）就是把它们的约束全部列出。

任何逻辑条件都可以化为**合取范式**（conjunctive normal form），即若干“或”子句的“与”，而每个子句 $\ell_1\lor\dots\lor\ell_r$ 都化为 $\sum_{\text{肯定 }\ell}\delta+\sum_{\text{否定 }\ell}(1-\delta)\ge1$（前一个和取遍肯定文字，后一个和取遍否定文字）。例如，$A\Rightarrow B$ 是子句 $\lnot A\lor B$，即 $(1-\delta_A)+\delta_B\ge1$，亦即 $\delta_B\ge\delta_A$。{numref}`tab-07-01-logic`列出了常见情形；每一行都通过枚举全部 $8$ 个 0–1 三元组与真值表核对过。

```{list-table} 0–1 变量的逻辑条件及其线性形式
:name: tab-07-01-logic
:header-rows: 1
:widths: 40 60

* - 条件
  - 线性约束
* - 若 $A$ 则 $B$
  - $\delta_B\ge\delta_A$
* - 若非 $A$ 则 $B$（$A$ 或 $B$）
  - $\delta_A+\delta_B\ge1$
* - $A$、$B$ 至多一个成立
  - $\delta_A+\delta_B\le1$
* - $A$、$B$ 恰有一个成立
  - $\delta_A+\delta_B=1$
* - $A$ 当且仅当 $B$
  - $\delta_A=\delta_B$
* - 若（$A$ 或 $C$）则 $B$
  - $\delta_B\ge\delta_A$ 且 $\delta_B\ge\delta_C$（较弱的单行形式：$2\delta_B\ge\delta_A+\delta_C$）
* - 若（$A$ 且 $C$）则 $B$
  - $\delta_B\ge\delta_A+\delta_C-1$
* - 若 $A$ 则（$B$ 或 $C$）
  - $\delta_B+\delta_C\ge\delta_A$
* - 若 $A$ 则（$B$ 且 $C$）
  - $\delta_B\ge\delta_A$ 且 $\delta_C\ge\delta_A$（较弱的单行形式：$\delta_B+\delta_C\ge2\delta_A$）
```

*（对于“若（$A$ 或 $C$）则 $B$”，有时会给出单行约束 $\delta_B\ge\delta_A+\delta_C$；这是错误的，因为它完全禁止了 $\delta_A=\delta_C=1$：此时右端为 $2$。）* 表中作为替代给出的单行形式与多行形式的整数解相同，但松弛较弱，正如{prf:ref}`ex-07-01-if-then-binaries`所示。

还可以把一个新的 0–1 变量*定义*为其他 0–1 变量的逻辑函数。对 $\omega=\delta_A\land\delta_B$（与），

$$
\begin{aligned}
\omega&\le\delta_A,\\
\omega&\le\delta_B,\\
\omega&\ge\delta_A+\delta_B-1,
\end{aligned}
$$ (eq-07-01-and)

对 $\omega=\delta_A\lor\delta_B$（或），有 $\omega\ge\delta_A$，$\omega\ge\delta_B$，$\omega\le\delta_A+\delta_B$；对异或，有 $\omega\ge\delta_A-\delta_B$，$\omega\ge\delta_B-\delta_A$，$\omega\le\delta_A+\delta_B$，$\omega\le2-\delta_A-\delta_B$。枚举四个数对 $(\delta_A,\delta_B)$ 即可验证每个定义。“与”变量就是乘积 $\delta_A\delta_B$，这正是{doc}`07-02-linearization-techniques`的出发点。

(sec-07-01-infrastructure)=
## 基础设施中的应用

下面两个例子综合运用了本章的各种建模模式。它们的最优解由 MILP 求解器求得，并通过枚举全部候选决策加以确认。

:::{prf:example} 交通项目的选择
:label: ex-07-01-project-selection

某市有 1.2 亿美元的资本预算，可用于六个候选项目（{numref}`tab-07-01-projects`），效益以净效益的现值计量。

```{list-table} 候选项目（金额单位：百万美元）
:name: tab-07-01-projects
:header-rows: 1

* - $j$
  - 项目
  - 成本
  - 效益
* - 1
  - 跨河大桥，方案 A
  - 40
  - 70
* - 2
  - 跨河大桥，方案 B
  - 30
  - 52
* - 3
  - 环路
  - 35
  - 50
* - 4
  - 环路上的立交桥
  - 15
  - 28
* - 5
  - 有轨电车线
  - 45
  - 75
* - 6
  - 快速公交（BRT）线
  - 20
  - 33
```

规则如下：(i) 两种大桥方案至多建一种；(ii) 立交桥只能与环路一起建；(iii) 三个走廊项目（环路、有轨电车、BRT）中至少建两个；(iv) 若有轨电车线和 BRT 线都建，则必须建一座大桥来承载其中之一。哪些项目能使总效益最大？将答案与仅受预算限制时的选择以及线性规划松弛进行比较。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*建模。* 若选择项目 $j$，则令 $\delta_j=1$，

$$
\begin{aligned}
\max_{\vdelta}\ & 70\delta_1+52\delta_2+50\delta_3+28\delta_4\\
& \quad+75\delta_5+33\delta_6\\
\st\ & 40\delta_1+30\delta_2+35\delta_3+15\delta_4\\
& \quad+45\delta_5+20\delta_6\le120, \\
& \delta_1+\delta_2\le1, \\
& \delta_4\le\delta_3, \\
& \delta_3+\delta_5+\delta_6\ge2, \\
& \delta_1+\delta_2\ge\delta_5+\delta_6-1, \\
& \vdelta\in\{0,1\}^6 .
\end{aligned}
$$

规则 (ii) 是一个**先后次序约束**（precedence constraint）（“若 4 则 3”，{numref}`tab-07-01-logic`的第一行）；规则 (iii) 是“$3$ 个中至少 $2$ 个”；规则 (iv) 是“若（5 且 6）则（1 或 2）”，由“与”那一行把 $\delta_B$ 换成 $\delta_1+\delta_2$ 得到（规则 (i) 保证这个和仍是 0–1 值）。

*求解。* 在 $64$ 个 0–1 向量中，有 $12$ 个满足全部约束。其中最好的选择是大桥方案 A、环路和有轨电车，$\vdelta^*=(1,0,1,0,1,0)$，成本恰为 1.2 亿美元，效益为 1.95 亿美元。若没有规则 (i)–(iv)，仅凭预算可以选择大桥方案 A、立交桥、有轨电车和 BRT 线（效益 2.06 亿美元），但这样就会在一条并不存在的道路上建立交桥。线性规划松弛给出 2.014 亿美元，例如取 $\vdelta=(1,0,0.3,0.3,1,1)$：分数项目填满了预算，并以 $0.3+1+1\ge2$ 满足规则 (iii)。
:::

:::{prf:example} 回收设施的开设
:label: ex-07-01-facility

四个候选地点可以建设回收设施，为五个城区服务。地点 $i$ 的年固定费用为 $F_i$（千美元/年），容量为 $u_i$（千吨/年）；城区 $j$ 每年产生 $d_j$ 千吨废弃物；把废弃物从城区 $j$ 运到地点 $i$ 的费用为每吨 $c_{ij}$ 美元（{numref}`tab-07-01-facility`）。应开设哪些地点？该模型的线性规划松弛有多强？

```{list-table} 设施开设问题的数据
:name: tab-07-01-facility
:header-rows: 1

* - 地点 $i$
  - $c_{i1}$
  - $c_{i2}$
  - $c_{i3}$
  - $c_{i4}$
  - $c_{i5}$
  - $F_i$
  - $u_i$
* - 1
  - 4
  - 6
  - 9
  - 10
  - 8
  - 300
  - 90
* - 2
  - 7
  - 3
  - 5
  - 9
  - 10
  - 250
  - 70
* - 3
  - 10
  - 8
  - 4
  - 3
  - 6
  - 280
  - 80
* - 4
  - 9
  - 10
  - 7
  - 5
  - 3
  - 220
  - 60
* - $d_j$
  - 30
  - 25
  - 35
  - 20
  - 30
  -
  -
```
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*建模。* 设 $x_{ij}$（千吨/年）为从城区 $j$ 运往地点 $i$ 的废弃物量；若开设地点 $i$，则 $\delta_i=1$。每吨 1 美元的单位费用作用于一千吨就是 1 千美元，因此目标函数的单位为千美元/年：

$$
\begin{aligned}
\min_{\vx,\vdelta}\ & \sum_{i=1}^4F_i\delta_i+\sum_{i=1}^4\sum_{j=1}^5c_{ij}x_{ij}\\
\st\ & \sum_{i=1}^4x_{ij}=d_j, & j&=1,\dots,5,\\
& \sum_{j=1}^5x_{ij}\le u_i\,\delta_i, & i&=1,\dots,4,\\
& x_{ij}\ge0, & i&=1,\dots,4,\ j=1,\dots,5,\\
& \delta_i\in\{0,1\}, & i&=1,\dots,4 .
\end{aligned}
$$

容量约束行就是固定费用联系约束 {eq}`eq-07-01-fixed-charge`：关闭的地点不处理任何废弃物。

*求解。* 废弃物总量为每年 $140$ 千吨，超过任何单个地点的容量，因此至少要开设两个地点。对每一个可行的开设地点组合求解运输线性规划，得到下表中的总费用（固定费用加运输费用，千美元/年）；地点 $\{2,4\}$ 的容量只有 $130$。

| 开设地点 | {1,2} | {1,3} | {1,4} | {2,3} | {3,4} | {1,2,3} | {1,2,4} | {1,3,4} | {2,3,4} | 全部 |
|---|---|---|---|---|---|---|---|---|---|---|
| 总费用 | 1350 | 1240 | 1275 | **1200** | 1260 | 1410 | 1330 | 1360 | 1325 | 1535 |

最优方案开设地点 2 和 3：固定费用 $250+280=530$，运输费用 $670$，总计每年 120 万美元。地点 2 接收城区 1 和城区 2 的废弃物（$30$ 和 $25$）以及城区 3 的 $5$；地点 3 接收城区 3 其余的 $30$ 以及城区 4 和城区 5 的废弃物，恰好达到其容量 $80$。MILP 求解器给出同样的解。

*三种模型形式的强弱。* 线性规划松弛的界取决于 $x_{ij}$ 与 $\delta_i$ 之间的联系约束如何书写：

1. 在 $\sum_jx_{ij}\le140\,\delta_i$ 中采用有效但宽松的常数 $M=140$（废弃物总量），再加上容量约束行 $\sum_jx_{ij}\le u_i$，线性规划的界为 $751.07$；
2. 采用上面模型中的 $M=u_i$，界为 $976.79$；
3. 再加上**分解形式**的联系约束 $x_{ij}\le\min(d_j,u_i)\,\delta_i$（若地点关闭，则每一单个流量都为零），界提高到 $1168.46$。

它们与最优值 $1200$ 的间隙分别为 $37\%$、$19\%$ 和 $2.6\%$。三种模型形式的整数解相同；应当交给求解器的是第三种。
:::

## 本章回顾

### 小结

- 0–1 变量 $\delta$ 通过 $g(\vx)\le M(1-\delta)$（{eq}`eq-07-01-bigm-relax`）打开或关闭 $g(\vx)\le0$，当 $M\ge\max_{\Omega_0}g$（{eq}`eq-07-01-valid-m`）时有效；对于箱体上的线性函数 $g$，紧的取值为 {eq}`eq-07-01-m-from-bounds`。等式约束要拆成两个不等式；反向蕴含需要 {eq}`eq-07-01-reverse`。
- $M$ 越大，线性规划松弛越弱，数值误差也越多；在没有已知界时，指示约束可以避免使用 $M$。
- 择一条件和若–则条件各用一个 0–1 变量（{eq}`eq-07-01-either-or`、{eq}`eq-07-01-if-then`）；“$m$ 个中至少 $k$ 个”为 {eq}`eq-07-01-k-of-m`，而“至多/恰有 $k$ 个”需要 {eq}`eq-07-01-at-most-k`。
- 固定费用和半连续变量分别用 $x\le u\delta$ 和 $l\delta\le x\le u\delta$ 建模（{eq}`eq-07-01-fixed-charge`、{eq}`eq-07-01-semicontinuous`）；松弛把固定费用摊为每单位 $F/u$（{numref}`fig-07-01-fixed-charge`）。
- 0–1 变量之间的逻辑关系可由否定、子句与合取得到（{numref}`tab-07-01-logic`、{eq}`eq-07-01-and`）。
- 析取条件的可行集是若干多面体的并集；其凸包是理想的松弛（{numref}`fig-07-01-either-or`），可由分解形式达到；紧常数和分解形式的约束行使大 $M$ 模型向它靠近（{prf:ref}`ex-07-01-facility`）。

### 常见错误

- **$M$ 过小。** 小于 $\max_{\Omega_0}g$ 的常数会悄无声息地切掉可行点。应由变量的界按 {eq}`eq-07-01-m-from-bounds` 推出 $M$，切勿凭猜测。
- **$M$ 过大。** 它虽然有效，却会削弱松弛并招致容差误差；{prf:ref}`ex-07-01-either-or`和{prf:ref}`ex-07-01-facility`显示了界随 $M$ 增大而变差。
- **把蕴含方向读反。** $g(\vx)\le M(1-\delta)$ 仅表示 $\delta=1\Rightarrow g(\vx)\le0$。若模型需要 $g(\vx)\le0\Rightarrow\delta=1$，应加上 {eq}`eq-07-01-reverse`。
- **仅用 $\sum\delta_i\le k$ 表达“至多 $k$ 个”。** 没有反向蕴含时，它限制的是被强制执行的约束个数，而不是被满足的约束个数。
- **严格不等式。** 不能施加 $g(\vx)>0$；应使用 $g(\vx)\ge\varepsilon$，其中 $\varepsilon$ 大于求解器容差（对整数数据可取 $1$）。
- **把逻辑约束行聚合。** 用一行表达多个蕴含关系是正确的，但较弱；应优先使用分解形式的各行。

## 参考文献

- Balas, E. (1985). Disjunctive programming and a hierarchy of relaxations
  for discrete optimization problems. *SIAM Journal on Algebraic and Discrete
  Methods, 6*(3), 466–486.
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Sections 10.1–10.2: modeling with
  integer variables, fixed charges, disjunctive constraints, facility
  location, strength of formulations.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Sections 12.2–12.3: binary decision
  variables, either–or, $K$ out of $N$, fixed charges.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Section 9.1: integer-programming applications, fixed-charge and
  either–or constraints.)
- Williams, H. P. (2013). *Model building in mathematical programming*
  (5th ed.). Wiley. (Chapters 9–10: building integer programming models,
  logical conditions and zero–one variables, disjunctions, alternative
  formulations and the strength of their relaxations.)
