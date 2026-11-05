# 分段线性逼近

## 概述

拥堵道路上的行程时间随交通量急剧增长；污水处理厂的费用则比其处理能力增长得慢。这类非线性函数可以被引入线性规划（linear program, LP）和混合整数线性规划（mixed-integer linear program, MILP），方法是用在选定的**分段点**（breakpoint）处与之相等的**分段线性**（piecewise-linear）函数来代替它们。本章说明具体做法。分段线性函数可以写成其分段点的凸组合，权重为 $\alpha_i$（**凸组合形式**（convex-combination formulation），又称 $\lambda$ 形式）；也可以写成已填满的各区段之和（**增量形式**（incremental formulation），又称 $\delta$ 形式）；还可以通过选择一个区段来表达（**多重选择形式**（multiple-choice formulation））。关键的区别在于凸性：当模型极小化凸函数（或极大化凹函数）时，一个普通的线性规划就足够了；否则权重必须满足一个相邻性条件，即构成**第二类特殊有序集**（special ordered set of type 2, SOS2），这需要 0–1 变量或求解器的支持。我们将给出光滑函数的逼近误差界，把这一思想推广到可分规划和二元函数，并按规模和线性规划松弛的强弱比较各种形式。本章求解了两个例子——一个拥堵的双路径网络和两座具有规模经济的污水处理厂——并把结果与精确最优解进行比较。本章以{doc}`07-01-logical-constraints`中的 0–1 变量和大 $M$ 约束以及{doc}`07-02-linearization-techniques`中仿射函数的最大值为基础。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 由分段点写出函数的分段线性插值函数，并判断它是否为凸函数；
- 解释为什么极小化凸分段线性函数不需要 0–1 变量，而非凸的情形需要；
- 用 SOS2 权重、0–1 相邻性约束、增量形式和多重选择形式表达分段线性函数，并对它们进行比较；
- 根据插值误差界为二阶可微函数选取分段点；
- 用混合整数线性模型逼近可分规划和二元函数；
- 将近似解与精确最优解进行比较。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $f$, $f^a$ | 被逼近的函数及其分段线性逼近 |
| $a_1<a_2<\dots<a_n$ | 分段点；$[a_k,a_{k+1}]$ 为区段 $k$，$k=1,\dots,n-1$ |
| $h_k=a_{k+1}-a_k$, $s_k$ | 区段 $k$ 的长度和斜率，$s_k=\bigl(f(a_{k+1})-f(a_k)\bigr)/h_k$ |
| $\alpha_i\ge0$ | 分段点 $i$ 的权重（凸组合形式），$\sum_i\alpha_i=1$ |
| $\delta_k\in\{0,1\}$ | 与区段 $k$ 对应的 0–1 变量 |
| $\xi_k\in[0,1]$ | 区段 $k$ 已填满的比例（增量形式） |
| $\bar x_k$ | $x$ 位于区段 $k$ 中的部分（区段形式、多重选择形式） |
| $\kappa_k$ | $\abs{f''}$ 在区段 $k$ 上的界 |
| $\beta_j\in\{0,1\}$ | 选择 $y$ 区间 $[b_j,b_{j+1}]$ 的 0–1 变量（二元函数） |
| $a_i$, $b_j$ | 函数 $f(x,y)$ 在 $x$ 和 $y$ 方向上的网格坐标 |
| $\varepsilon$ | 小正容差 |

许多教材把权重记为 $\lambda_i$，因此有“$\lambda$ 形式”之名；本书把 $\lambda$ 留给拉格朗日乘子。

(sec-07-03-pwl)=
## 分段线性函数

设 $f$ 定义在 $[a_1,a_n]$ 上，选取分段点 $a_1<a_2<\dots<a_n$。**分段线性插值函数**（piecewise-linear interpolant）$f^a$ 用直线段依次连接各点 $\bigl(a_i,f(a_i)\bigr)$。点 $\bar x\in[a_k,a_{k+1}]$ 是其所在区段两端点的凸组合，$\bar x=\alpha_ka_k+\alpha_{k+1}a_{k+1}$，其中 $\alpha_k+\alpha_{k+1}=1$，$\alpha_k,\alpha_{k+1}\ge0$，并且

$$
\begin{aligned}
f^a(\bar x)&=\alpha_kf(a_k)+\alpha_{k+1}f(a_{k+1})\\
&=f(a_k)+s_k(\bar x-a_k),\\
\alpha_{k+1}&=\frac{\bar x-a_k}{h_k} .
\end{aligned}
$$ (eq-07-03-interp)

$f^a$ 连续，且在每个区段上是线性的。它是**凸函数**当且仅当其斜率递增，$s_1\le s_2\le\dots\le s_{n-1}$，此时它等于其各段的最大值，

$$
f^a(x)=\max_{k=1,\dots,n-1}\bigl(f(a_k)+s_k(x-a_k)\bigr),
$$

这正是{ref}`sec-07-02-convex-side`中所处理的形式。若 $f$ 是凸函数，则其插值函数也是凸的，并且位于 $f$ 的**上方**（凸函数的弦位于其图像上方）；若 $f$ 是凹函数，则 $f^a$ 是凹的，并且位于 $f$ 的**下方**。兼有两种曲率的函数，或者由带折扣的价目表直接给出的分段线性费用，一般既非凸也非凹。

{numref}`fig-07-03-interpolation`用{ref}`sec-07-03-examples`中两个例子的函数展示了这两种情形。

```{figure} figures/07-03-interpolation.svg
:name: fig-07-03-interpolation
:width: 660px
:class: dark-light
:alt: 一个凸的行程时间函数和一个凹的费用函数，以及它们的分段线性插值函数

左图：{prf:ref}`ex-07-03-congestion` 中路径 1 上的总行程时间 $T_1(v)$（凸函数）及其每隔 $1000$ 辆/小时取一个分段点的插值函数，插值函数位于其上方。右图：{prf:ref}`ex-07-03-economies` 中污水处理厂的凹费用函数 $C_1(q)=q^{0.6}$，每隔 $20$ 百万升/日取一个分段点；插值函数位于曲线下方，从 $0$ 到 $80$ 的虚线弦是不带相邻性条件的线性规划所采用的凸包络。
```

(sec-07-03-convex-case)=
## 凸的情形不需要 0–1 变量

对每个分段点配一个权重，写出 $x$ 和 $f^a(x)$：

$$
\begin{aligned}
x&=\sum_{i=1}^n\alpha_ia_i,\\
f^a&=\sum_{i=1}^n\alpha_if(a_i),\\
\sum_{i=1}^n\alpha_i&=1,\qquad \alpha_i\ge0 .
\end{aligned}
$$ (eq-07-03-cc)

这些条件描述的是分段点的所有凸组合，而不仅仅是两个相邻分段点的组合；数对 $(x,f^a)$ 可以是各点 $\bigl(a_i,f(a_i)\bigr)$ 的凸包中的任何一点。只有当至多两个权重为正且它们**相邻**时，才能还原出插值函数。

:::{prf:proposition} 凸的情形下相邻性自动成立
:label: prop-07-03-convex-case

设 $f^a$ 是凸函数，某线性模型*极小化*一个目标函数，其中 $f^a$ 以正系数出现，并且所有其他涉及 $x$ 的约束都只涉及 $x$（不涉及各个权重）。则带 {eq}`eq-07-03-cc` 的线性规划有一个至多两个权重为正且二者相邻的最优解；其最优值等于采用精确插值函数的模型的最优值。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

固定 $x$。在所有满足 $\sum_i\alpha_ia_i=x$、$\sum_i\alpha_i=1$、$\valpha\ge\vzero$ 的权重中，值 $\sum_i\alpha_if(a_i)$ 至少为 $f^a(x)$，因为 $f^a$ 是凸函数且在分段点处与 $f$ 相等：由 Jensen 不等式，$\sum_i\alpha_if(a_i)=\sum_i\alpha_if^a(a_i)\ge f^a\bigl(\sum_i\alpha_ia_i\bigr)$。{eq}`eq-07-03-interp` 中的两个相邻权重达到了这个下界。用它们替换任意最优权重，既不改变 $x$，也不增大目标函数值，因此采用相邻权重时模型是最优的，其最优值就是采用 $f^a$ 时的最优值。
:::

对于极大化凹函数 $f^a$，同样的结论成立。**区段形式**给出一个等价的线性规划：令 $x=a_1+\sum_k\bar x_k$，$0\le\bar x_k\le h_k$，$f^a=f(a_1)+\sum_ks_k\bar x_k$。斜率递增时，极小化会先填满较便宜的区段，因此在前一个区段填满之前绝不会使用后一个区段，不需要任何次序约束。这正是{prf:ref}`ex-03-01-piecewise-profit`和{prf:ref}`prop-03-01-concave-split`中分段线性利润背后的论证：在三百箱以内每百箱贡献 500 美元、超出部分每百箱贡献 300 美元，这在极大化问题中是凹的，因此线性规划会先按 500 美元的价格销售。

对于**非凸**的情形，线性规划会出问题。在极小化凹费用时，它倾向于用两个最外端的分段点来表示 $x$，因为二者之间的弦位于曲线下方；于是线性规划优化的是 $f^a$ 的**凸包络**（convex envelope）（{numref}`fig-07-03-interpolation`中的虚线），而不是 $f^a$ 本身。这时必须显式地施加相邻性条件。

(sec-07-03-sos)=
## 特殊有序集与凸组合形式

:::{prf:definition} 特殊有序集
:label: def-07-03-sos

由非负变量组成的有序列表 $(\alpha_1,\dots,\alpha_n)$，若其中至多一个非零，则称为**第一类特殊有序集**（special ordered set of type 1, SOS1）；若其中至多两个非零，且当有两个非零时它们在次序上是相邻的，则称为**第二类特殊有序集**（SOS2）。
:::

特殊有序集由 Beale 和 Tomlin (1970) 提出，同时提出的还有一种无需 0–1 变量即可保证其成立的分支规则：若线性规划的解违反 SOS2 条件，求解器就对适当的下标 $r$，在“对所有 $i<r$ 有 $\alpha_i=0$”与“对所有 $i>r$ 有 $\alpha_i=0$”之间进行分支。一个常见的 SOS1 是恰有一个等于 $1$ 的一组 0–1 变量，$\sum_i\delta_i=1$，$\delta_i\in\{0,1\}$：这正是下文要用到的多重选择条件。

在可以使用 SOS2 时，$f^a$ 的凸组合形式为

$$
\begin{aligned}
x&=\sum_{i=1}^n\alpha_ia_i,\\
f^a&=\sum_{i=1}^n\alpha_if(a_i),\\
\sum_{i=1}^n\alpha_i&=1,\qquad \alpha_i\ge0,\\
&(\alpha_1,\dots,\alpha_n)\ \text{为 SOS2} .
\end{aligned}
$$ (eq-07-03-sos2)

在没有求解器支持时，相邻性通过为每个区段配一个 0–1 变量 $\delta_k$（“$x$ 位于区段 $k$ 中”）来保证，并约定 $\delta_0=\delta_n=0$：

$$
\begin{aligned}
&\alpha_i\le\delta_{i-1}+\delta_i,\quad i=1,\dots,n,\\
&\sum_{k=1}^{n-1}\delta_k=1,\\
&\delta_k\in\{0,1\},\quad k=1,\dots,n-1 .
\end{aligned}
$$ (eq-07-03-lambda)

恰好选中一个区段 $k$；此时只有 $\alpha_k$ 和 $\alpha_{k+1}$ 可以为正，因为其他每个权重的两个相邻 0–1 变量都等于 $0$。在两端写出来，第一行和最后一行分别是 $\alpha_1\le\delta_1$ 和 $\alpha_n\le\delta_{n-1}$。目标函数项为 $\sum_{i=1}^n\alpha_if(a_i)$，每个分段点一项。*（二重和 $\sum_{i=1}^{n-1}\bigl[\alpha_if(a_i)+\alpha_{i+1}f(a_{i+1})\bigr]$ 把内部分段点计算了两次，是错误的。区段指示变量必须是 0–1 变量；若 $\delta_k\in[0,1]$，则约束行 {eq}`eq-07-03-lambda` 只是线性规划松弛。）*

(sec-07-03-incremental)=
## 增量形式

**增量形式**（又称 $\delta$ 形式）从左到右依次填满各区段。设 $\xi_k\in[0,1]$ 为区段 $k$ 已填满的比例：

$$
\begin{aligned}
x&=a_1+\sum_{k=1}^{n-1}h_k\,\xi_k,\\
f^a&=f(a_1)+\sum_{k=1}^{n-1}\bigl(f(a_{k+1})-f(a_k)\bigr)\xi_k,\\
\xi_{k+1}&\le\delta_k\le\xi_k,\quad k=1,\dots,n-2,\\
0&\le\xi_k\le1,\quad k=1,\dots,n-1,\\
\delta_k&\in\{0,1\},\quad k=1,\dots,n-2 .
\end{aligned}
$$ (eq-07-03-incremental)

0–1 变量 $\delta_k$ 表示“区段 $k$ 已填满”。若 $\delta_k=1$，则 $\xi_k=1$；若 $\delta_k=0$，则 $\xi_{k+1}=0$。因此各比例具有 $(1,\dots,1,\xi_k,0,\dots,0)$ 的形式，这正是区段 $k$ 上的插值函数 {eq}`eq-07-03-interp`。该形式需要 $n-2$ 个 0–1 变量。在凸的情形下，这些 0–1 变量同样可以去掉（这就是{ref}`sec-07-03-convex-case`中的区段形式）。这个名称与本书用 $\delta$ 表示 0–1 变量的习惯相冲突；称之为 $\delta$ 形式的教材把比例记为 $\delta_k$，把 0–1 变量记为 $y_k$。

(sec-07-03-multiple-choice)=
## 多重选择形式

第三种经典形式选择一个区段，并为它保留 $x$ 的一个副本 $\bar x_k$：

$$
\begin{aligned}
x&=\sum_{k=1}^{n-1}\bar x_k,\\
f^a&=\sum_{k=1}^{n-1}\Bigl(s_k\bar x_k+\bigl(f(a_k)-s_ka_k\bigr)\delta_k\Bigr),\\
a_k\delta_k&\le\bar x_k\le a_{k+1}\delta_k,\quad k=1,\dots,n-1,\\
\sum_{k=1}^{n-1}\delta_k&=1,\qquad \delta_k\in\{0,1\} .
\end{aligned}
$$ (eq-07-03-multiple-choice)

若 $\delta_k=1$，则只有 $\bar x_k$ 可以为正，它位于区段 $k$ 中，费用就是该区段上的直线。这正是{ref}`sec-07-01-either-or`中析取条件“$x$ 位于区段 1 或 2 或 $\dots$”的分解形式（凸包形式）。

(sec-07-03-comparison)=
## 各种形式的比较

{numref}`tab-07-03-compare`汇总了具有 $n$ 个分段点的单个函数的各种形式。

```{list-table} 具有 $n$ 个分段点的分段线性函数的各种形式
:name: tab-07-03-compare
:header-rows: 1

* - 形式
  - 连续变量
  - 0–1 变量
  - 说明
* - 凸组合，线性规划 {eq}`eq-07-03-cc`
  - $n$
  - $0$
  - 仅在凸的情形下精确
* - 带 SOS2 的凸组合 {eq}`eq-07-03-sos2`
  - $n$
  - $0$
  - 通过 SOS2 分支保证相邻性
* - 带 0–1 变量的凸组合 {eq}`eq-07-03-lambda`
  - $n$
  - $n-1$
  - 松弛的顶点处 $\delta$ 可能取分数值
* - 增量形式 {eq}`eq-07-03-incremental`
  - $n-1$
  - $n-2$
  - 松弛的顶点处 $\delta$ 取整数值
* - 多重选择形式 {eq}`eq-07-03-multiple-choice`
  - $n-1$
  - $n-1$
  - 松弛的顶点处 $\delta$ 取整数值
```

以下三点观察均已在下面的例子上做过数值验证。

1. **界。** 三种 0–1 形式的线性规划松弛给出相同的界：每一种都把 $f^a$ 换成了它在 $[a_1,a_n]$ 上的凸包络。在{prf:ref}`ex-07-03-economies`中，三者给出的界都是 1 784 万美元，而真实最优值为 1 989 万美元。
2. **松弛的整数性。** 增量形式和多重选择形式是*局部理想的*（locally ideal）：其线性规划松弛的每个顶点处 0–1 变量都取整数值（Vielma et al., 2010）。带 0–1 变量的凸组合形式则不然：在一项以具有五个分段点的单个函数为对象、采用 $3000$ 个随机目标函数的实验中，其松弛有 $65$ 个最优顶点的 $\delta$ 取分数值，而增量形式和多重选择形式没有一个如此。对局部理想的形式，分支定界法通常需要的节点较少。
3. **规模与便利性。** 凸组合形式最容易书写，也最容易推广到多个变量（{ref}`sec-07-03-two-variables`）；在有 SOS2 支持时，它完全不增加 0–1 变量。

(sec-07-03-error)=
## 非线性函数的逼近

逼近光滑函数 $f$ 需要选取分段点。线性插值的误差由曲率控制。

:::{prf:theorem} 插值误差
:label: thm-07-03-error

设 $f$ 在 $[a_k,a_{k+1}]$ 上二阶连续可微，且在其上 $\abs{f''(x)}\le \kappa_k$。则对该区段中的每个 $x$，有

$$
\abs{f(x)-f^a(x)}\le\frac{h_k^2}{8}\,\kappa_k .
$$ (eq-07-03-error)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

固定 $x\in(a_k,a_{k+1})$，令 $e(t)=f(t)-f^a(t)$，它在 $a_k$ 和 $a_{k+1}$ 处为零。定义 $w(t)=e(t)-e(x)\,\dfrac{(t-a_k)(t-a_{k+1})}{(x-a_k)(x-a_{k+1})}$。则 $w$ 在 $a_k$、$x$ 和 $a_{k+1}$ 处为零，因此两次应用罗尔定理可知，区段中存在 $t^*$ 使 $w''(t^*)=0$。由于在该区段上 $(f^a)''=0$，有 $w''(t^*)=f''(t^*)-2e(x)/[(x-a_k)(x-a_{k+1})]$，因此

$$
e(x)=\tfrac12f''(t^*)\,(x-a_k)(x-a_{k+1}) .
$$

因子 $\abs{(x-a_k)(x-a_{k+1})}$ 在中点处最大，等于 $h_k^2/4$。由此即得 {eq}`eq-07-03-error`。
:::

对建模的启示：

* 为保证某区段上的误差至多为 $\epsilon_f$，取 $h_k\le\sqrt{8\epsilon_f/\kappa_k}$。间距减半，误差变为原来的四分之一。
* 在 $\abs{f''}$ 大的地方密集地布置分段点，在 $f$ 近似线性的地方稀疏地布置；等间距布置会浪费分段点。
* 把范围 $[a_1,a_n]$ 限制在解实际可能取到的值上；对于陡峭的函数，宽范围末端附近的误差会主导误差界。
* 对凸函数 $f$，插值函数会高估（$f^a\ge f$），因此用 $f^a$ 进行极小化给出真实最优值的一个上界，并且当 $f$ 只出现在目标函数中时，其解对原问题是可行的；对凹函数 $f$，插值函数会低估。

(sec-07-03-separable)=
## 可分规划

若函数是若干一元函数之和，$f(\vx)=\sum_{j=1}^nf_j(x_j)$，则称它是**可分的**（separable）。**可分规划**（separable program）具有可分的目标函数和可分的约束 $g_i(\vx)=\sum_jg_{ij}(x_j)\le0$。每个一元函数都用它在自己的分段点上的插值函数代替，并使用各自的权重；结果是一个带 SOS2 条件的线性规划（Williams, 2013, 第 7 章; Bazaraa et al., 2006, §11.3）。

* 若每个 $f_j$ 都是凸函数（极小化），且每个 $g_{ij}$ 都是凸函数，则由{prf:ref}`prop-07-03-convex-case`及其关于约束的类似结论（$\le$ 约束中的凸函数 $g_{ij}$ 会被任何非相邻的组合高估，这只会使约束更难满足），相邻性条件可以去掉：可分规划作为线性规划求解。
* 否则，SOS2 条件要通过 SOS2 分支、通过 {eq}`eq-07-03-lambda` 中的 0–1 变量，或者（在历史上）通过单纯形法的**限制进基**（restricted basis entry）规则来保证，该规则只允许与某个基变量相邻的权重进基。最后这种方法只能找到逼近问题的局部最优解。

有些不可分的项可以通过变量代换化为可分的。经典的例子是乘积：令 $u_1=\tfrac12(x+y)$，$u_2=\tfrac12(x-y)$，

$$
xy=u_1^2-u_2^2,
$$

这是两个一元函数之差；然后对每个平方项进行分段线性逼近（项 $-u_2^2$ 是凹的，需要相邻性条件）。对于有界变量之积，{ref}`sec-07-02-mccormick`中的 McCormick 包络是另一种选择。

(sec-07-03-examples)=
## 两个例子

第一个例子是凸的，作为线性规划求解；第二个例子是凹的，需要 0–1 变量。两者都与精确最优解进行比较。

:::{prf:example} 考虑拥堵的系统最优路径分配
:label: ex-07-03-congestion

每小时 $D=6000$ 辆车的需求经由两条平行路径从郊区前往市中心。行程时间服从交通分配中使用的美国公路局（Bureau of Public Roads, BPR）路段阻抗函数（Sheffi, 1985），

$$
t_r(v)=t_r^0\Bigl[1+0.15\Bigl(\frac{v}{c_r}\Bigr)^4\Bigr]\ \text{分钟},
$$

其中自由流行程时间为 $t_1^0=10$ 分钟、$t_2^0=15$ 分钟，通行能力为 $c_1=3000$ 辆/小时、$c_2=4000$ 辆/小时。求使总行程时间 $T_1(v_1)+T_2(v_2)$ 最小的分配 $v_1+v_2=D$，其中 $T_r(v)=v\,t_r(v)/60$ 的单位为辆·小时/小时。先求精确最优解，再在两条路径上每隔 $1000$ 辆/小时取一个分段点，求解分段线性规划，并比较二者。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*精确最优解。* $T_r$ 是凸函数。令两条路径的导数（边际行程时间 $\mathrm d(v\,t_r)/\mathrm dv=t_r^0[1+0.75(v/c_r)^4]$）相等，得到 $v_1^*=3070.9$ 辆/小时和 $v_2^*=2929.1$ 辆/小时，两条路径上的边际时间都是 $18.235$ 分钟，实际行程时间分别为 $11.647$ 分钟和 $15.647$ 分钟，总计 $1359.97$ 辆·小时/小时。

*分段线性规划。* 在两条路径上都取分段点 $0,1000,\dots,6000$ 辆/小时，插值点处的函数值为

| $v$（辆/小时） | 0 | 1000 | 2000 | 3000 | 4000 | 5000 | 6000 |
|---|---|---|---|---|---|---|---|
| $T_1(v)$（辆·小时/小时） | 0 | 166.98 | 343.21 | 575.00 | 982.72 | 1797.84 | 3400.00 |
| $T_2(v)$（辆·小时/小时） | 0 | 250.15 | 504.69 | 785.60 | 1150.00 | 1707.76 | 2639.06 |

线性规划为

$$
\begin{aligned}
\min\ & \sum_{i=1}^7\bigl(\alpha_{1i}T_1(a_i)+\alpha_{2i}T_2(a_i)\bigr)\\
\st\ & \sum_{i=1}^7a_i(\alpha_{1i}+\alpha_{2i})=6000,\\
& \sum_{i=1}^7\alpha_{1i}=1,\\
& \sum_{i=1}^7\alpha_{2i}=1,\\
& \alpha_{ri}\ge0,\qquad r=1,2,\ i=1,\dots,7,
\end{aligned}
$$

其中 $a_i=1000(i-1)$；由于两个函数都是凸的，不需要 0–1 变量。其解为 $\alpha_{14}=\alpha_{24}=1$，即 $v_1=v_2=3000$ 辆/小时，目标函数值为 $575.00+785.60=1360.60$ 辆·小时/小时。用斜率（以每辆车的分钟数计的边际时间）验证：把流量从路径 2 移到路径 1，在路径 1 的区段 $[3000,4000]$ 上增加 $24.46$，在路径 2 的区段 $[2000,3000]$ 上节省 $16.85$；反方向移动则增加 $21.86$、节省 $13.91$。两种移动都不划算，因此该线性规划解是最优的。

*比较。* 近似分配与精确分配相差 $71$ 辆/小时，其真实总行程时间为 $1360.60$ 辆·小时/小时（这是精确值，因为两个流量都在分段点上），只比最优值多 $0.63$ 辆·小时/小时（$0.05\%$）。插值本身则要粗糙得多：由 {eq}`eq-07-03-error`，取 $\kappa_k=T_1''(6000)=\tfrac{10}{60}\cdot3\cdot6000^3/3000^4=\tfrac1{750}$，路径 1 上的误差界为 $1000^2/8\cdot\tfrac1{750}=166.7$ 辆·小时/小时，而实际最大误差为 $129.2$，出现在 $v=5523$ 附近。在最优解附近的 $[2000,4000]$ 上，最大误差为 $33.6$。目标函数在其极小点附近很平坦，因此粗糙的逼近仍能找到接近最优的分配；若每隔 $500$ 辆/小时取一个分段点，误差大约降为原来的四分之一，分别为 $36.7$ 和 $10.2$。
:::

:::{prf:example} 具有规模经济的污水处理厂
:label: ex-07-03-economies

某地区需要 $100$ 百万升/日的处理能力，由两座处理厂分担，每座的处理能力至多为 $80$ 百万升/日。以百万美元/年计的年化费用呈现规模经济：处理厂 1 为 $C_1(q)=q^{0.6}$，处理厂 2 为 $C_2(q)=0.8\,q^{0.6}+0.06\,q$（建造较便宜，但每百万升/日的处理能力每年有 6 万美元的泵送费用）。在 $q_1+q_2=100$、$0\le q_1,q_2\le80$ 的约束下极小化 $C_1(q_1)+C_2(q_2)$。求精确最优解；然后每隔 $20$ 百万升/日取一个分段点逼近两个费用函数，分别在不带和带相邻性条件的情况下求解凸组合模型。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*精确最优解。* 两个费用函数都是凹的，因此总费用沿可行线段 $q_1\in[20,80]$ 是凹的，其最小值在端点处取得：$C_1(80)+C_2(20)=13.863+6.027=19.89$，即 1 989 万美元，而 $C_1(20)+C_2(80)$ 为 2 192 万美元。最优解为 $q_1^*=80$，$q_2^*=20$。

*分段点。* 取 $a_i=0,20,40,60,80$ 百万升/日，

| $q$（百万升/日） | 0 | 20 | 40 | 60 | 80 |
|---|---|---|---|---|---|
| $C_1(q)$（百万美元/年） | 0 | 6.034 | 9.146 | 11.665 | 13.863 |
| $C_2(q)$（百万美元/年） | 0 | 6.027 | 9.717 | 12.932 | 15.890 |

斜率递减（对 $C_1$：$0.302$、$0.156$、$0.126$、$0.110$），因此插值函数是凹的。

*不带相邻性条件的线性规划。* 对两座处理厂都采用模型 {eq}`eq-07-03-cc`，并令 $\sum_ia_i(\alpha_{1i}+\alpha_{2i})=100$，得到费用 1 784 万美元，解为 $\alpha_{15}=1$（$q_1=80$）以及 $\alpha_{21}=0.75$、$\alpha_{25}=0.25$：它把 $q_2=20$ 表示为 $0$ 和 $80$ 的混合，计费 $0.25\cdot15.890=3.97$，而不是 $C_2(20)=6.03$。这些权重不相邻，所得费用是凸包络的费用。

*MILP。* 为每座处理厂加上 0–1 变量 {eq}`eq-07-03-lambda`（每座四个），得到 $\alpha_{15}=1$、$\alpha_{22}=1$，即 $q_1=80$、$q_2=20$，费用 1 989 万美元，正是精确最优解（最优解落在分段点上，因此那里的逼近误差为零）。增量形式和多重选择形式给出相同的解。三个 MILP 的线性规划松弛给出相同的界 1 784 万美元，即不带相邻性条件的线性规划的值：$C_1(80)+20\cdot C_2(80)/80$。
:::

{numref}`fig-07-03-sos2`仅就处理厂 1 说明了{prf:ref}`ex-07-03-economies`中所见的相邻权重与非相邻权重之间的区别：点 $q=50$ 用分段点 $40$ 和 $60$ 上的 SOS2 权重 $\alpha_3=\alpha_4=0.5$ 表示，得到 $f^a(50)=10.41$（真实值为 $10.46$）；若不要求相邻性，则可用 $0$ 和 $80$ 上的权重 $\alpha_1=0.375$、$\alpha_5=0.625$ 表示，得到弦上的值 $8.66$。

```{figure} figures/07-03-sos2.svg
:name: fig-07-03-sos2
:width: 660px
:class: dark-light
:alt: 同一点用分段点权重的两种表示：一种用相邻权重落在曲线上，另一种用非相邻权重落在弦上

左图：凹费用函数 $C_1(q)=q^{0.6}$、其插值函数、从 $0$ 到 $80$ 的弦（虚线）以及 $q=50$ 的两种表示。右图：权重 $\alpha_i$。SOS2 权重（两个相邻的非零权重）再现了插值函数；非相邻权重则把该点放在曲线下方的弦上。
```

(sec-07-03-two-variables)=
## 二元函数

非线性函数 $f(x,y)$ 在网格上逼近：$x$ 方向为 $a_1<\dots<a_n$，$y$ 方向为 $b_1<\dots<b_m$，网格把定义域分成 $(n-1)(m-1)$ 个矩形。常用的方法有四种（D'Ambrosio et al., 2010）。

**一维方法。** 为每个 $y$ 区间取一个代表值 $\tilde b_j\in[b_j,b_{j+1}]$（通常取 $\tilde b_j=b_j$），用权重为 $\alpha_i$（SOS2）的凸组合形式在 $x$ 方向上逼近 $f(x,\tilde b_j)$。0–1 变量 $\beta_j$，$j=1,\dots,m-1$，用于选择 $y$ 区间：

$$
\begin{aligned}
&\textstyle\sum_j\beta_jb_j\le y\le\sum_j\beta_jb_{j+1},\\
&\textstyle\sum_{j=1}^{m-1}\beta_j=1,\qquad \beta_j\in\{0,1\},\\
&\textstyle f^a\le\sum_i\alpha_if(a_i,\tilde b_j)+M(1-\beta_j),\quad j=1,\dots,m-1,\\
&\textstyle f^a\ge\sum_i\alpha_if(a_i,\tilde b_j)-M(1-\beta_j),\quad j=1,\dots,m-1 .
\end{aligned}
$$

对选中的 $j$，两行大 $M$ 约束给出 $f^a=\sum_i\alpha_if(a_i,\tilde b_j)$；其余各行被松弛。这种逼近在 $x$ 方向上是分段线性的，但在 $y$ 方向上是分段*常数*的，因此比较粗糙。在网格值 $y=b_j$ 处，有两个区间都符合条件，模型可能选取更有利的那个值；把上面一行换成 $y\le\sum_j\beta_jb_{j+1}-\varepsilon$ 可以消除这种歧义，代价是排除了一条宽度为 $\varepsilon$ 的窄带。*（对 $\beta_j$ 的选择条件是 $\sum_j\beta_j=1$，即一个 SOS1；若和等于 $0$，则不会选中任何区间。）*

**三角形方法。** 用从 $(a_i,b_j)$ 到 $(a_{i+1},b_{j+1})$ 的对角线把每个矩形分成上三角形和下三角形，在每个三角形上利用其三个顶点进行线性插值。在网格点上配权重 $\alpha_{ij}\ge0$，

$$
\begin{aligned}
\sum_{i,j}\alpha_{ij}&=1,\\
x&=\sum_{i,j}\alpha_{ij}a_i,\\
y&=\sum_{i,j}\alpha_{ij}b_j,\\
f^a&=\sum_{i,j}\alpha_{ij}f(a_i,b_j),
\end{aligned}
$$

并且非零权重必须位于同一个三角形的顶点上。为矩形 $(i,j)$ 的上三角形和下三角形分别配 0–1 变量 $\delta^u_{ij}$ 和 $\delta^l_{ij}$（下标超出范围时取零），

$$
\begin{aligned}
&\sum_{i=1}^{n-1}\sum_{j=1}^{m-1}\bigl(\delta^u_{ij}+\delta^l_{ij}\bigr)=1,\\
&\alpha_{ij}\le\delta^u_{ij}+\delta^l_{ij}+\delta^u_{i-1,j-1}+\delta^l_{i-1,j-1}\\
&\qquad+\delta^u_{i,j-1}+\delta^l_{i-1,j},\quad i=1,\dots,n,\ j=1,\dots,m .
\end{aligned}
$$

这六项对应以 $(a_i,b_j)$ 为顶点的三角形：其右上方矩形和左下方矩形的两个三角形，其右下方矩形的上三角形，以及其左上方矩形的下三角形（{numref}`fig-07-03-triangles`）。这种逼近是连续的，并且在网格点处是精确的，但需要 $nm$ 个权重和 $2(n-1)(m-1)$ 个 0–1 变量。*（有些教材把这个权重条件称为“SOS3”；它是一种三角形相邻性条件，而不是标准的特殊有序集。）*

```{figure} figures/07-03-triangles.svg
:name: fig-07-03-triangles
:width: 420px
:class: dark-light
:alt: 一个三角剖分的网格，某一网格点周围的六个三角形以阴影表示

用从 $(a_i,b_j)$ 到 $(a_{i+1},b_{j+1})$ 的对角线对网格进行三角剖分。只有当选中的三角形是标记网格点周围六个阴影三角形之一时，该网格点的权重 $\alpha_{ij}$ 才可以为正。
```

**矩形方法。** 沿选中矩形的下边在 $x$ 方向上进行线性插值，再加上 $y$ 方向上的线性修正。令 $\gamma=(y-b_j)/(b_{j+1}-b_j)$，$\Delta(i,j)=f(a_i,b_{j+1})-f(a_i,b_j)$，

$$
\begin{aligned}
f^a(x,y)={}&\alpha f(a_i,b_j)+(1-\alpha)f(a_{i+1},b_j)\\
&+\gamma\min\bigl\{\Delta(i,j),\Delta(i+1,j)\bigr\}
\end{aligned}
$$

其中 $x=\alpha a_i+(1-\alpha)a_{i+1}$。若 $\Delta(i,j)\le\Delta(i+1,j)$，这是经过角点 $(a_i,b_j)$、$(a_{i+1},b_j)$ 和 $(a_i,b_{j+1})$ 的平面；否则是经过 $(a_i,b_j)$、$(a_{i+1},b_j)$ 和 $(a_{i+1},b_{j+1})$ 的平面。*（这里的分情形比较的是两个 $x$ 端点处的 $y$ 方向增量 $\Delta(i,j)$ 和 $\Delta(i+1,j)$。）* 在 MILP 中，0–1 变量 $\delta^x_i$ 和 $\beta_j$ 分别选择 $x$ 区间和 $y$ 区间，连续变量 $\gamma_j\le\beta_j$ 通过 $y=\sum_j(\beta_jb_j+\gamma_j(b_{j+1}-b_j))$ 承载 $y$ 方向的比例，每个矩形用两行大 $M$ 约束

$$
\begin{aligned}
f^a&\le\sum_k\alpha_kf(a_k,b_j)+\gamma_j\Delta^{\min}_{ij}+M(2-\beta_j-\delta^x_i),\\
f^a&\ge\sum_k\alpha_kf(a_k,b_j)+\gamma_j\Delta^{\min}_{ij}-M(2-\beta_j-\delta^x_i),
\end{aligned}
$$

（其中 $\Delta^{\min}_{ij}=\min\{\Delta(i,j),\Delta(i+1,j)\}$）确定选中矩形上的 $f^a$。*（松弛项在两行中的符号相反：上界中为 $+M$，下界中为 $-M$。）* 把 $\bar f_j=\sum_k\alpha_kf(a_k,b_j)$ 作为单独的变量引入可以缩短这些约束行。该方法只需要 $O(n+m)$ 个变量和 0–1 变量，但需要 $O(nm)$ 行大 $M$ 约束。

**切平面方法。** Wang 和 Lo (2010) 在每个矩形的左下角处用一阶泰勒展开把一个二元行程时间函数线性化：

$$
\begin{aligned}
f^a(x,y)&=f(a_i,b_j)+f_x(a_i,b_j)(x-a_i)+f_y(a_i,b_j)(y-b_j)\\
&=p_{ij}x+q_{ij}y+r_{ij},
\end{aligned}
$$

其中 $p_{ij}=f_x(a_i,b_j)$，$q_{ij}=f_y(a_i,b_j)$，$r_{ij}=f(a_i,b_j)-p_{ij}a_i-q_{ij}b_j$。*（常数项使用的是网格点 $(a_i,b_j)$，而不是被逼近的点 $(x,y)$。）* 这种逼近在矩形边界上不连续。矩形由 0–1 变量定位：若 $x<a_i$，则 $\chi^x_i=1$；若 $y<b_j$，则 $\chi^y_j=1$：

$$
\begin{aligned}
x&\ge a_i-M\chi^x_i, & i&=1,\dots,n,\\
x&\le a_i-\varepsilon+M(1-\chi^x_i), & i&=1,\dots,n,
\end{aligned}
$$

对 $y$ 用 $\chi^y_j$ 写出同样的约束行。由于 $a_i$ 递增，$\chi^x_i=1$ 蕴含 $\chi^x_{i+1}=1$，因此 $\zeta^x_i=\chi^x_{i+1}-\chi^x_i\in\{0,1\}$ 等于 $1$ 当且仅当 $a_i\le x<a_{i+1}$；对 $y$ 同样有 $\zeta^y_j=\chi^y_{j+1}-\chi^y_j$。*（组合 $\chi^x_i=\chi^x_{i+1}=1$ 表示 $x<a_i$；$\chi^x_i=1$、$\chi^x_{i+1}=0$ 是不可能的。）* 最后，写成两行的 $\abs{f^a-(p_{ij}x+q_{ij}y+r_{ij})}\le M(2-\zeta^x_i-\zeta^y_j)$ 确定选中矩形上的 $f^a$。固定 $\chi^x_1=0$ 和 $\chi^x_n=1$ 把 $x$ 限制在 $[a_1,a_n)$ 中。一个等价的模型写成 $f^a=\sum_{i,j}(p_{ij}x+q_{ij}y+r_{ij})\zeta^x_i\zeta^y_j$，其中 0–1 变量与连续变量之积按{ref}`sec-07-02-binary-continuous`中的方法线性化。

```{list-table} $n\times m$ 网格上各二元函数逼近形式的规模（数量级）
:name: tab-07-03-two-var
:header-rows: 1

* - 方法
  - 连续变量
  - 0–1 变量
  - 约束
  - 逼近方式
* - 一维方法
  - $n$
  - $(n-1)+(m-1)$
  - $O(n+m)$
  - $x$ 方向线性，$y$ 方向常数
* - 三角形方法
  - $nm$
  - $2(n-1)(m-1)$
  - $nm+5$
  - 连续，在三角形上线性
* - 矩形方法
  - $n+(m-1)$
  - $(n-1)+(m-1)$
  - $O(nm)$
  - 在半矩形上线性
* - 切平面方法
  - $(n-1)+(m-1)$
  - $n+m$
  - $O(nm)$
  - 一阶泰勒展开，不连续
```

{numref}`tab-07-03-two-var`比较了各形式的规模（一维方法和矩形方法的 SOS2 权重计为连续变量，其相邻性 0–1 变量计入 0–1 变量）。三角形方法最精确，规模也最大；一维方法规模最小，也最粗糙；矩形方法是一种折中。

## 本章回顾

### 小结

- 以 $a_1<\dots<a_n$ 为分段点的插值函数为 {eq}`eq-07-03-interp`；当斜率递增时它是凸函数，此时它等于若干仿射段的最大值。
- 极小化凸（或极大化凹）插值函数是线性规划问题：在最优解处权重是相邻的（{prf:ref}`prop-07-03-convex-case`）。
- 非凸的情形需要 SOS2 权重（{prf:ref}`def-07-03-sos`、{eq}`eq-07-03-sos2`）、0–1 相邻性约束 {eq}`eq-07-03-lambda`、增量形式 {eq}`eq-07-03-incremental` 或多重选择形式 {eq}`eq-07-03-multiple-choice`（{numref}`tab-07-03-compare`、{numref}`fig-07-03-sos2`）。
- 它们的线性规划松弛给出凸包络界；增量形式和多重选择形式的松弛顶点是整数的。
- 插值误差至多为 $h_k^2\kappa_k/8$（{prf:ref}`thm-07-03-error`）；应在曲率大的地方布置分段点。
- 可分规划逐个函数地进行逼近；凸的可分规划是线性规划问题。
- 二元函数用一维方法、三角形方法、矩形方法或切平面方法逼近（{numref}`tab-07-03-two-var`）。
- 在{prf:ref}`ex-07-03-congestion`中，粗糙的线性规划与精确最优值相差不到 $0.05\%$；在{prf:ref}`ex-07-03-economies`中，只有 MILP 能找到真正的最优解。

### 常见错误

- **在非凸模型中去掉相邻性条件。** 此时线性规划优化的是凸包络，其“解”混合了相距很远的分段点（{prf:ref}`ex-07-03-economies`）。
- **在凸模型中添加 0–1 变量。** 它们无害但无用：极小化凸插值函数是线性规划问题（{prf:ref}`prop-07-03-convex-case`）。
- **在凸组合形式中按区段而不是按分段点对目标函数求和**，这会把内部分段点重复计算。
- **把区段指示变量声明为连续变量。** 此时 {eq}`eq-07-03-lambda` 只是一个松弛。
- **分段点的范围过宽。** 误差界 {eq}`eq-07-03-error` 取决于每个区段上的最大曲率；应限制范围，或在 $\abs{f''}$ 大的地方加密。
- **在带大 $M$ 约束行的二元函数模型中忽略边界上的并列情形**；应使用 $\varepsilon$，或接受在网格线上取更有利的值。

## 参考文献

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Section 11.3:
  separable programming.)
- Beale, E. M. L., & Tomlin, J. A. (1970). Special facilities in a general
  mathematical programming system for nonconvex problems using ordered sets
  of variables. In J. Lawrence (Ed.), *Proceedings of the Fifth International
  Conference on Operational Research* (pp. 447–454). Tavistock.
- D'Ambrosio, C., Lodi, A., & Martello, S. (2010). Piecewise linear
  approximation of functions of two variables in MILP models. *Operations
  Research Letters, 38*(1), 39–46.
  <https://doi.org/10.1016/j.orl.2009.09.005>
- Sheffi, Y. (1985). *Urban transportation networks: Equilibrium analysis with
  mathematical programming methods*. Prentice-Hall. (Chapter 1: link
  performance functions; Chapter 3: the user-equilibrium and system-optimum
  formulations.)
- Vielma, J. P., Ahmed, S., & Nemhauser, G. (2010). Mixed-integer models for
  nonseparable piecewise-linear optimization: Unifying framework and
  extensions. *Operations Research, 58*(2), 303–315.
- Wang, D. Z. W., & Lo, H. K. (2010). Global optimum of the linearized
  network design problem with equilibrium flows. *Transportation Research
  Part B: Methodological, 44*(4), 482–492.
  <https://doi.org/10.1016/j.trb.2009.10.003>
- Williams, H. P. (2013). *Model building in mathematical programming*
  (5th ed.). Wiley. (Chapter 7: nonlinear models and separable programming;
  Section 9.3: special ordered sets of variables.)
