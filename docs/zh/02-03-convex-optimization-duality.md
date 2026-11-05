# 凸优化与拉格朗日对偶

## 概述

凸优化问题是在由凸不等式和仿射方程描述的凸集上极小化一个凸函数。这类问题具有一般问题所不具备的两个性质：每个局部极小点都是全局极小点（{prf:ref}`thm-02-03-local-global`），并且最优性可以用一个关于梯度的简单条件来验证（{prf:ref}`thm-02-03-optimality`）。本章定义凸优化问题，说明线性规划和凸二次规划都属于这一类，然后建立**拉格朗日对偶**（Lagrange duality）理论：拉格朗日函数；给出最优值 $p^*$ 下界的对偶函数 $q$；对偶问题及其最优值 $d^*$；弱对偶性 $d^*\le p^*$（{prf:ref}`thm-02-03-weak-duality`）；Slater 条件下的强对偶性 $d^*=p^*$（{prf:ref}`thm-02-03-slater`）；互补松弛性与凸问题的 KKT 条件（Karush–Kuhn–Tucker conditions）（{prf:ref}`thm-02-03-kkt-sufficient`）；以及乘子作为最优值灵敏度的解释（{prf:ref}`thm-02-03-sensitivity`）。线性规划的对偶问题由拉格朗日函数导出；理论通过电力系统的经济调度问题、最小范数问题和一个存在对偶间隙的小问题加以演示。本章以{doc}`02-01-convex-sets`中的分离超平面和支撑超平面、{doc}`02-02-convex-functions`，以及{doc}`01-03-introduction`中的一般形式和基本术语为基础。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 识别标准形式的凸优化问题，并解释为什么其等式约束必须是仿射的；
- 证明凸问题的局部极小点是全局极小点，并应用最优性条件 $\nabla f(\vx^*)^\T(\vy-\vx^*)\ge0$；
- 写出问题的拉格朗日函数和对偶函数，并对线性和二次问题求出对偶函数的闭式表达式；
- 叙述并运用弱对偶性、强对偶性、Slater 条件和互补松弛性；
- 写出凸问题的 KKT 条件，对小规模问题求解这些条件，并把乘子解释为价格；
- 由拉格朗日函数导出线性规划的对偶问题。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\mathcal D$ | 问题的定义域：$f$、所有 $g_i$ 和所有 $h_j$ 都有定义的点 |
| $\ell$ | 等式约束的个数 |
| $L(\vx,\vlambda,\vmu)$ | 拉格朗日函数 |
| $\lambda_i$, $\mu_j$ | $g_i(\vx)\le0$ 和 $h_j(\vx)=0$ 的乘子 |
| $q(\vlambda,\vmu)$ | 拉格朗日对偶函数 |
| $p^*$, $d^*$ | 原问题和对偶问题的最优值 |
| $\vg(\vx)$, $\vh(\vx)$ | 向量 $(g_1(\vx),\dots,g_m(\vx))^\T$ 和 $(h_1(\vx),\dots,h_\ell(\vx))^\T$ |
| $\mathcal G$, $(u,t)$ | 只含一个不等式约束的问题的（约束值，目标函数值）集合及其坐标 |
| $\vu$, $\vv$; $p^*(\vu,\vv)$ | $g_i\le0$ 和 $h_j=0$ 右端项的摄动；摄动后的最优值 |

(sec-02-03-convex-problems)=
## 凸优化问题

回顾{doc}`01-03-introduction`中的一般问题：

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)\\
\st\ & g_i(\vx)\le0, && i=1,\dots,m,\\
& h_j(\vx)=0, && j=1,\dots,\ell,
\end{aligned}
$$ (eq-02-03-primal)

其定义域为 $\mathcal D=\dom f\cap\bigcap_i\dom g_i\cap\bigcap_j\dom h_j$，可行集为 $\Omega$，最优值为 $p^*=\inf\{f(\vx):\vx\in\Omega\}$。

:::{prf:definition} 凸优化问题
:label: def-02-03-convex-problem

若 $f$ 和 $g_1,\dots,g_m$ 是凸函数，$h_1,\dots,h_\ell$ 是仿射函数 $h_j(\vx)=\va_j^\T\vx-b_j$，则称问题 {eq}`eq-02-03-primal` 为**凸优化问题**（convex optimization problem），也称凸规划（convex program）。此时等式约束可以写成 $\mA\vx=\vb$。
:::

凸问题的可行集是凸集：它是下水平集 $\{\vx:g_i(\vx)\le0\}$（由{prf:ref}`thm-02-02-sublevel`知其为凸集）与超平面 $\{\vx:\va_j^\T\vx=b_j\}$ 的交（{prf:ref}`thm-02-01-intersection`）。因此凸问题就是在凸集上极小化凸函数。在凸集上极大化凹函数也是凸问题，因为它等价于极小化该函数的相反数（{eq}`eq-01-03-max-min`）。

等式约束必须是仿射的，仅仅是凸的还不够。集合 $\{\vx:h(\vx)=0\}=\{\vx:h(\vx)\le0\}\cap\{\vx:-h(\vx)\le0\}$ 要使上述论证成立，需要 $h$ 和 $-h$ 都是凸函数，这意味着 $h$ 是仿射的。若 $h$ 是凸的但非线性，该集合通常不是凸集：$x_1^2+x_2^2-1=0$ 是一个圆周。*（因此，允许“凸等式约束”的定义是错误的；凸函数 $h$ 只能出现在不等式 $h(\vx)\le0$ 中。）*

:::{prf:theorem} 凸问题的局部极小点是全局极小点
:label: thm-02-03-local-global

对凸优化问题，每个局部极小点都是全局极小点。最优解集是凸集；若 $f$ 严格凸，则最优解至多有一个。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

设 $\vx^*$ 是局部极小点：对所有满足 $\norm{\vx-\vx^*}<\varepsilon$ 的可行点 $\vx$ 有 $f(\vx^*)\le f(\vx)$。假设存在可行点 $\vy$ 使 $f(\vy)<f(\vx^*)$。对 $\theta\in(0,1]$，由于 $\Omega$ 是凸集，点 $\vz_\theta=\vx^*+\theta(\vy-\vx^*)$ 是可行的，并且由 $f$ 的凸性

$$
f(\vz_\theta)\le(1-\theta)f(\vx^*)+\theta f(\vy)<f(\vx^*).
$$

当 $\theta<\varepsilon/\norm{\vy-\vx^*}$ 时，点 $\vz_\theta$ 位于 $\vx^*$ 的 $\varepsilon$ 邻域内，这与局部最优性矛盾。因此对每个可行点 $\vy$ 有 $f(\vy)\ge f(\vx^*)$。

最优解集为 $\Omega\cap\{\vx:f(\vx)\le p^*\}$，是凸集的交。若 $\vx^*\ne\vy^*$ 都是最优解且 $f$ 严格凸，则它们的中点可行且 $f(\tfrac12\vx^*+\tfrac12\vy^*)<p^*$，这不可能。
:::

这个定理是凸性之所以重要的主要实际原因。局部搜索方法只把一个点与其邻近的点比较，却能找到凸问题的全局最优解。对非凸问题则会出现{prf:ref}`def-01-03-local-global`中的情形，{ref}`sec-06-01-difficulties`讨论了其中的困难。

(sec-02-03-optimality)=
### 一个最优性条件

对可微的目标函数，一阶条件 {eq}`eq-02-02-first-order` 导出最优点的一个刻画。

:::{prf:theorem} 凸问题的最优性条件
:label: thm-02-03-optimality

设 $f$ 是可微凸函数，$\Omega$ 是凸集。点 $\vx^*\in\Omega$ 是 $f$ 在 $\Omega$ 上的极小点，当且仅当

$$
\nabla f(\vx^*)^\T(\vy-\vx^*)\ge0\qquad\text{对所有 }\vy\in\Omega .
$$ (eq-02-03-optimality)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

若 {eq}`eq-02-03-optimality` 成立，由一阶条件，对每个 $\vy\in\Omega$ 有 $f(\vy)\ge f(\vx^*)+\nabla f(\vx^*)^\T(\vy-\vx^*)\ge f(\vx^*)$。反之，若对某个 $\vy\in\Omega$ 有 $\nabla f(\vx^*)^\T(\vy-\vx^*)<0$，则 $\theta\mapsto f(\vx^*+\theta(\vy-\vx^*))$ 在 $\theta=0$ 处的导数为负，所以对充分小的 $\theta>0$ 有 $f(\vx^*+\theta(\vy-\vx^*))<f(\vx^*)$；由凸性这些点都是可行的，因而 $\vx^*$ 不是最优的。
:::

条件 {eq}`eq-02-03-optimality` 是说：没有任何可行方向 $\vy-\vx^*$ 是下降方向。从几何上看，若 $\nabla f(\vx^*)\ne\vzero$，则过 $\vx^*$、以 $-\nabla f(\vx^*)$ 为法向量的超平面在 $\vx^*$ 处支撑 $\Omega$（{prf:ref}`thm-02-01-supporting-hyperplane`）：整个可行集都位于 $f$ 在一阶意义下不减小的那一侧。有三种特殊情形值得单独写出。

- **无约束问题**（$\Omega=\R^n$）：取 $\vy=\vx^*\pm\ve_j$ 得 $\nabla f(\vx^*)=\vzero$。对凸函数 $f$，梯度为零是全局极小的充分必要条件。
- **等式约束**（$\Omega=\{\vx:\mA\vx=\vb\}$）：方向 $\vy-\vx^*$ 取遍 $\mA$ 的零空间，对每个方向的 $\pm$ 两个方向应用 {eq}`eq-02-03-optimality`，迫使 $\nabla f(\vx^*)$ 与该零空间正交，即对某个 $\vmu$ 有 $\nabla f(\vx^*)=-\mA^\T\vmu$（$\mA$ 的行空间中的向量，{doc}`01-01-linear-algebra`）。向量 $\vmu$ 就是拉格朗日乘子向量。
- **非负约束**（$\Omega=\R^n_+$）：条件变为 $\nabla f(\vx^*)\ge\vzero$，且对所有 $j$ 有 $x_j^*\,\partial f(\vx^*)/\partial x_j=0$。取正值的分量的偏导数必为零，取下界零的分量的偏导数必非负。这是互补松弛性的第一个例子。

(sec-02-03-lp-qp)=
## 线性规划与二次规划

**线性规划问题**（linear program，LP）的目标函数和约束都是仿射的；其不等式形式为

$$
\begin{aligned}
\min_{\vx}\ & \vc^\T\vx\\
\st\ & \mA\vx\le\vb ,
\end{aligned}
$$

其标准形为 $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$。仿射函数是凸函数，所以每个 LP 都是凸问题，其可行集是多面体（{prf:ref}`def-02-01-polyhedron`），每个局部最优解都是全局最优解。线性规划各章（{doc}`03-01-lp-model-building`至{doc}`03-07-sensitivity-analysis`）利用了更多的结构：标准形 LP 若有最优解，则可以在某个极点上找到最优解（{prf:ref}`thm-03-02-fundamental`）。

**二次规划**（quadratic programming，QP）问题的目标函数是二次的，约束是仿射的：

$$
\begin{aligned}
\min_{\vx}\ & \tfrac12\vx^\T\mQ\vx+\vc^\T\vx\\
\st\ & \mA\vx\le\vb,\\
& \mE\vx=\vd,
\end{aligned}
$$ (eq-02-03-qp)

其中 $\mQ$ 对称。它是凸问题当且仅当 $\mQ\succeq0$（{doc}`02-02-convex-functions`）。带线性限制的最小二乘拟合、以方差为风险的投资组合选择、以及燃料成本为二次函数的发电机组调度，都是凸 QP。若 $\mQ$ 不定，问题是非凸的，可能有多个局部极小点。下面的例子贯穿本章。

:::{prf:example} 两台发电机组的经济调度
:label: ex-02-03-dispatch

某电力系统须由两台发电机组满足 $D=475$ 兆瓦的负荷需求。机组 1 发电 $x_1$ 兆瓦，成本为 $C_1(x_1)=8x_1+0.01x_1^2$ 美元/小时；机组 2 发电 $x_2$ 兆瓦，成本为 $C_2(x_2)=10x_2+0.02x_2^2$ 美元/小时。与火电厂的情况一样，**边际成本**（marginal cost）$C_1'(x_1)=8+0.02x_1$ 和 $C_2'(x_2)=10+0.04x_2$（美元/兆瓦时）随出力增加而增大。调度问题为

$$
\begin{aligned}
\min_{\vx}\ & 8x_1+0.01x_1^2+10x_2+0.02x_2^2\\
\st\ & x_1+x_2=475 .
\end{aligned}
$$ (eq-02-03-dispatch)

它是一个凸 QP，$\mQ=\diag(0.02,0.04)\succ0$，只有一个仿射等式约束。（出力下限 $x_1,x_2\ge0$ 被略去，因为它们最终是非有效约束；{prf:ref}`ex-02-03-dispatch-capacity`中将加入出力上限。）

用{prf:ref}`thm-02-03-optimality`的最优性条件求最优调度方案，并解释其乘子。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

把约束写成 $h(\vx)=475-x_1-x_2=0$，$\nabla h=(-1,-1)^\T$。由{prf:ref}`thm-02-03-optimality`的等式约束情形，最优解对某个 $\mu$ 满足 $\nabla f(\vx^*)+\mu\nabla h(\vx^*)=\vzero$。因此

$$
\begin{aligned}
8+0.02x_1&=\mu,\\
10+0.04x_2&=\mu :
\end{aligned}
$$

**两台机组以相同的边际成本 $\mu$ 运行。** 如果一台机组的边际成本低于另一台，把少量出力转移给它就能节省费用。把 $x_1=(\mu-8)/0.02$ 和 $x_2=(\mu-10)/0.04$ 代入 $x_1+x_2=475$，得 $50\mu-400+25\mu-250=475$，所以 $\mu=15$ 美元/兆瓦时，且

$$
\begin{aligned}
\vx^*&=(350,\ 125)\ \text{兆瓦},\\
p^*&=C_1(350)+C_2(125)=4025+1562.5=5587.5\ \text{\$/h}.
\end{aligned}
$$

共同的边际成本 $\mu=15$ 美元/兆瓦时就是该系统中的电价：多满足 1 兆瓦的需求，每小时大约多花 15 美元。下面的拉格朗日理论将解释这一含义。
:::

(sec-02-03-lagrangian)=
## 拉格朗日函数与对偶函数

从现在起，不再要求问题 {eq}`eq-02-03-primal` 是凸的。拉格朗日对偶的思想是给约束标上价格，把它们移入目标函数。

:::{prf:definition} 拉格朗日函数与对偶函数
:label: def-02-03-lagrangian

{eq}`eq-02-03-primal` 的**拉格朗日函数**（Lagrangian）是函数 $L:\mathcal D\times\R^m\times\R^\ell\to\R$，

$$
\begin{aligned}
L(\vx,\vlambda,\vmu)&=f(\vx)+\sum_{i=1}^m\lambda_ig_i(\vx)+\sum_{j=1}^\ell\mu_jh_j(\vx)\\
&=f(\vx)+\vlambda^\T\vg(\vx)+\vmu^\T\vh(\vx).
\end{aligned}
$$ (eq-02-03-lagrangian)

数 $\lambda_i$ 和 $\mu_j$ 称为约束 $g_i(\vx)\le0$ 和 $h_j(\vx)=0$ 的**拉格朗日乘子**（Lagrange multiplier），也称对偶变量（dual variable）。**拉格朗日对偶函数**（Lagrange dual function）为

$$
q(\vlambda,\vmu)=\inf_{\vx\in\mathcal D}L(\vx,\vlambda,\vmu).
$$ (eq-02-03-dual-function)
:::

下确界是在不带约束的定义域上取的，并且可能等于 $-\infty$。对偶函数之所以有用，在于以下两个事实。

*它是凹函数*，无论原问题是否为凸问题。对每个固定的 $\vx$，$L(\vx,\vlambda,\vmu)$ 是 $(\vlambda,\vmu)$ 的仿射函数，而 $q$ 是这族仿射函数的逐点下确界，因此是凹函数（{prf:ref}`thm-02-02-operations`，第 3 部分）。

*只要 $\vlambda\ge\vzero$，它就给出最优值的下界。*

:::{prf:theorem} 下界性质
:label: thm-02-03-lower-bound

若 $\vlambda\ge\vzero$，则对每个 $\vmu$ 有 $q(\vlambda,\vmu)\le p^*$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

设 $\tilde\vx$ 可行。由 $g_i(\tilde\vx)\le0$ 和 $\lambda_i\ge0$ 得 $\vlambda^\T\vg(\tilde\vx)\le0$，且 $\vh(\tilde\vx)=\vzero$。因此

$$
\begin{aligned}
q(\vlambda,\vmu)&\le L(\tilde\vx,\vlambda,\vmu)\\
&=f(\tilde\vx)+\vlambda^\T\vg(\tilde\vx)+\vmu^\T\vh(\tilde\vx)\\
&\le f(\tilde\vx).
\end{aligned}
$$

对所有可行点 $\tilde\vx$ 取下确界即得 $q\le p^*$。（若问题不可行，则 $p^*=+\infty$，无需证明。）
:::

拉格朗日函数有一种经济学解释。把 $\lambda_i$ 看作违反约束 $i$ 时每单位的价格。极小化 $L$ 允许决策者违反约束，只要为此付费；同时对留有松弛的约束按同样的价格给予奖励（$g_i(\vx)<0$ 得到 $-\lambda_ig_i(\vx)>0$）。价格非负时，每个可行点要么不付费，要么得到奖励，所以 $L$ 的无约束极小值只可能低于 $f$ 的约束极小值。若价格选得恰好使约束在 $L$ 的极小点处成立，下界就最好：这正是强对偶性的内容。

(sec-02-03-dual-problem)=
## 对偶问题与弱对偶性

每个满足 $\vlambda\ge\vzero$ 的 $(\vlambda,\vmu)$ 都给出 $p^*$ 的一个下界 $q(\vlambda,\vmu)$。其中最好的下界是**拉格朗日对偶问题**（Lagrange dual problem）

$$
\begin{aligned}
\max_{\vlambda,\vmu}\ & q(\vlambda,\vmu)\\
\st\ & \vlambda\ge\vzero .
\end{aligned}
$$ (eq-02-03-dual)

的最优值，记作 $d^*$。满足 $\vlambda\ge\vzero$ 且 $q(\vlambda,\vmu)>-\infty$ 的 $(\vlambda,\vmu)$ 称为**对偶可行**（dual feasible）的，达到 $d^*$ 的对偶可行点称为**对偶最优**的，记作 $(\vlambda^*,\vmu^*)$。与 {eq}`eq-02-03-primal` 不同，问题 {eq}`eq-02-03-dual` *总是*凸问题：它在凸集上极大化凹函数。条件 $q>-\infty$ 通常在求出 $q$ 之后表现为显式约束，如{ref}`sec-02-03-lp-dual`中 LP 的对偶问题。

:::{prf:theorem} 弱对偶性
:label: thm-02-03-weak-duality

对每个问题 {eq}`eq-02-03-primal`，无论是否为凸问题，都有 $d^*\le p^*$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

在{prf:ref}`thm-02-03-lower-bound`中对 $\vlambda\ge\vzero$ 和 $\vmu$ 取上确界即可。
:::

差 $p^*-d^*\ge0$ 称为**对偶间隙**（duality gap）。弱对偶性有两个直接用途。对偶可行点给出一个有保证的下界：若可行点 $\tilde\vx$ 和对偶可行点 $(\vlambda,\vmu)$ 满足 $f(\tilde\vx)-q(\vlambda,\vmu)\le\varepsilon$，则 $\tilde\vx$ 在 $\varepsilon$ 的误差内是最优的，许多算法正是以此作为终止检验。另外，若对偶问题无界，即 $d^*=+\infty$，则原问题不可行。

:::{prf:example} 调度问题的对偶函数
:label: ex-02-03-dispatch-dual

取 $h(\vx)=475-x_1-x_2$，计算 {eq}`eq-02-03-dispatch` 的对偶函数，并求解对偶问题。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

拉格朗日函数可以拆成每台机组各一项：

$$
\begin{aligned}
L(\vx,\mu)
&=8x_1+0.01x_1^2+10x_2+0.02x_2^2+\mu\,(475-x_1-x_2)\\
&=\bigl(0.01x_1^2-(\mu-8)x_1\bigr)+\bigl(0.02x_2^2-(\mu-10)x_2\bigr)+475\mu .
\end{aligned}
$$

每个方括号都是开口向上的抛物线。第一个在 $x_1=(\mu-8)/0.02$ 处取极小值 $-(\mu-8)^2/0.04$，第二个在 $x_2=(\mu-10)/0.04$ 处取极小值 $-(\mu-10)^2/0.08$。因此

$$
q(\mu)=475\mu-\frac{(\mu-8)^2}{0.04}-\frac{(\mu-10)^2}{0.08},
$$ (eq-02-03-dispatch-q)

这是一条开口向下的抛物线，对所有 $\mu$ 取有限值。若以价格 $\mu$ 收购每台机组的出力，每台机组都会选择使其边际成本等于 $\mu$ 的出力；$q(\mu)$ 就是这些出力的总成本，加上供给相对于需求的缺额 $475-x_1-x_2$ 按价格 $\mu$ 计算的价值（供给超过需求时为负）。

等式约束的乘子是自由的，所以对偶问题是 $q$ 的无约束极大化。令 $q'(\mu)=475-(\mu-8)/0.02-(\mu-10)/0.04=0$，得 $75\mu=1125$，所以 $\mu^*=15$ 美元/兆瓦时，且

$$
d^*=q(15)=7125-1225-312.5=5587.5\ \text{\$/h}=p^* .
$$

{numref}`fig-02-03-dispatch-dual`画出了 $q$。$q$ 的每个值都是最优成本的下界（例如 $q(10)=4650$ 和 $q(20)=4650$ 美元/小时），最好的下界等于 $p^*$，极大点 $\mu^*$ 重现了{prf:ref}`ex-02-03-dispatch`中相等的边际成本。在 $\mu^*$ 处，两台机组各自选择的出力 $x_1=350$、$x_2=125$ 恰好满足需求。

```{figure} figures/02-03-dispatch-dual.svg
:name: fig-02-03-dispatch-dual
:width: 520px
:class: dark-light
:alt: 关于 mu 的开口向下的抛物线 q，在 mu 等于 15 处取最大值 5587.5，与水平线 p-star 等于 5587.5 相切

调度问题的对偶函数 {eq}`eq-02-03-dispatch-q`。每个值 $q(\mu)$ 都是最优成本 $p^*=5587.5$ 美元/小时（虚线）的下界；在价格 $\mu^*=15$ 美元/兆瓦时处取得的最大值等于 $p^*$，因此不存在对偶间隙。
```
:::

(sec-02-03-strong-duality)=
## 强对偶性与 Slater 条件

**强对偶性**（strong duality）是指 $d^*=p^*$。它一般并不成立（见{ref}`sec-02-03-duality-gap`），但对大多数凸问题成立。保证强对偶性成立的条件称为**约束规范**（constraint qualification），其中最简单的是 Slater 给出的条件。

:::{prf:theorem} Slater 条件
:label: thm-02-03-slater

考虑凸问题 {eq}`eq-02-03-primal`，其中 $\vh(\vx)=\mA\vx-\vb$。设它是**严格可行**（strictly feasible）的：存在 $\mathcal D$ 的相对内部中的点 $\tilde\vx$，使得

$$
\begin{aligned}
g_i(\tilde\vx)&<0\quad(i=1,\dots,m),\\
\mA\tilde\vx&=\vb .
\end{aligned}
$$

则 $d^*=p^*$；并且若 $p^*$ 有限，则对偶最优值可以达到：存在对偶最优解 $(\vlambda^*,\vmu^*)$。对仿射的约束 $g_i$，不需要严格不等式。
:::

证明用一个超平面分离两个凸集（{prf:ref}`thm-02-01-separating-hyperplane`）；Boyd 和 Vandenberghe（2004，§5.3.2）给出了完整证明。其思想可以从下面的几何图景中看出。对本书来说有两个重要推论：

- 对约束全为**仿射**的凸问题（如 LP 或 QP），Slater 条件退化为可行性。若这样的问题可行且 $p^*$ 有限，则强对偶性成立，且对偶最优值可以达到。对 LP，这就是强对偶定理（{prf:ref}`thm-03-05-strong-duality`）；
- 当定义域 $\mathcal D=\R^n$ 时，相对内部就是整个 $\R^n$，条件只要求存在一个严格满足非线性不等式的点。

调度问题 {eq}`eq-02-03-dispatch` 只有一个仿射约束且是可行的，这解释了{prf:ref}`ex-02-03-dispatch-dual`中对偶间隙为零的结果。

(sec-02-03-geometry)=
### 几何解释

对只含一个不等式约束、不含等式约束的问题，考虑（约束值，目标函数值）这一对数所取值的集合

$$
\mathcal G=\{(g(\vx),\ f(\vx)):\vx\in\mathcal D\}\subseteq\R^2 ,
$$

其坐标记为 $(u,t)$。最优值是 $\mathcal G$ 中满足 $u\le0$ 的点的最小 $t$：$p^*=\inf\{t:(u,t)\in\mathcal G,\ u\le0\}$。对偶函数为

$$
q(\lambda)=\inf\{t+\lambda u:(u,t)\in\mathcal G\} .
$$

对固定的 $\lambda\ge0$，直线 $t+\lambda u=\text{常数}$ 的斜率为 $-\lambda$，下确界由仍与 $\mathcal G$ 接触的最低的这种直线取得，即 $\mathcal G$ 的一条支撑线。它与轴 $u=0$ 交于高度 $q(\lambda)$。由于这条直线斜率非正，且位于整个 $\mathcal G$ 之下，特别是位于满足 $u\le0$ 的点之下，所以其截距至多为 $p^*$；这就是弱对偶性的图示。对偶问题倾斜这条支撑线，使截距尽可能升高。当 $\mathcal G$（更确切地说，是 $\mathcal G$ 连同其右上方的全部点构成的集合）是凸集时，最好的支撑线经过点 $(0,p^*)$，且 $d^*=p^*$；满足 $u<0$ 的 Slater 点使这条直线不会变成竖直的。{numref}`fig-02-03-duality-geometry`同时给出了一个凸的情形和一个有对偶间隙的非凸情形。

```{figure} figures/02-03-duality-geometry.svg
:name: fig-02-03-duality-geometry
:width: 680px
:class: dark-light
:alt: 左图：(u, t) 平面中的抛物线形曲线 G，一条斜率为负 3 的支撑线在 u 等于 0 处与之相切，另有一条较平缓、截距较低的直线。右图：两个孤立点及过这两点的直线，其截距 0.5 低于最优值 1

一个不等式约束下对偶性的几何解释。(a) 带 {prf:ref}`ex-02-03-dispatch-capacity` 中容量限制的调度问题，化为只含变量 $x_1$ 的问题：$\mathcal G$ 是黑色曲线 $(x_1-300,\ C_1(x_1)+C_2(475-x_1))$。斜率为 $-\lambda^*=-3$ 的支撑线与 $u=0$ 交于 $q(3)=5662.5=p^*$；斜率为 $-1$ 的直线与之交于更低处 $q(1)=5629.17$。(b) {prf:ref}`ex-02-03-duality-gap` 中只能开或关的水泵：$\mathcal G$ 只有两个点，最好的支撑线与 $u=0$ 交于 $d^*=0.5<p^*=1$。灰色带标出可行的一侧 $u\le0$。
```

(sec-02-03-complementary-slackness)=
## 互补松弛性

设强对偶性成立，$\vx^*$ 是原问题最优解，$(\vlambda^*,\vmu^*)$ 是对偶最优解。则

$$
\begin{aligned}
f(\vx^*)=q(\vlambda^*,\vmu^*)
&=\inf_{\vx\in\mathcal D}\Bigl(f(\vx)+\sum_i\lambda_i^*g_i(\vx)+\sum_j\mu_j^*h_j(\vx)\Bigr)\\
&\le f(\vx^*)+\sum_i\lambda_i^*g_i(\vx^*)+\sum_j\mu_j^*h_j(\vx^*)\\
&\le f(\vx^*).
\end{aligned}
$$

最后一步用到 $\lambda_i^*\ge0$、$g_i(\vx^*)\le0$ 和 $h_j(\vx^*)=0$。第一个表达式与最后一个表达式相等，所以两个不等式都取等号。由此得到两个结论。

1. $\vx^*$ 在 $\mathcal D$ 上极小化拉格朗日函数 $L(\vx,\vlambda^*,\vmu^*)$。最优乘子把约束问题变成一个无约束问题，而 $\vx^*$ 是后者的解之一。
2. $\sum_i\lambda_i^*g_i(\vx^*)=0$。每一项都 $\le0$，所以每一项都为零：

$$
\lambda_i^*\,g_i(\vx^*)=0,\qquad i=1,\dots,m .
$$ (eq-02-03-complementary-slackness)

条件 {eq}`eq-02-03-complementary-slackness` 称为**互补松弛性**（complementary slackness）：$\lambda_i^*>0$ 蕴含 $g_i(\vx^*)=0$（价格为正的约束是有效约束，或称紧约束）；$g_i(\vx^*)<0$ 蕴含 $\lambda_i^*=0$（有松弛的约束价格为零）。从经济上说，没有用完的资源在边际上没有价值。在 LP 中，这就是互补松弛定理（{prf:ref}`thm-03-05-complementary-slackness`）。

(sec-02-03-kkt)=
## 凸问题的 KKT 条件

现在设 $f$、$g_i$ 和 $h_j$ 都可微。在强对偶性下，最优的 $\vx^*$ 极小化 $L(\cdot,\vlambda^*,\vmu^*)$，因此它关于 $\vx$ 的梯度在 $\vx^*$ 处为零。把这一条件与可行性和互补松弛性合在一起，就得到 **KKT 条件**（Karush–Kuhn–Tucker conditions）：

$$
\begin{aligned}
&\text{原始可行性：} && g_i(\vx)\le0 && \text{对所有 }i,\\
& && h_j(\vx)=0 && \text{对所有 }j,\\
&\text{对偶可行性：} && \lambda_i\ge0 && \text{对所有 }i,\\
&\text{互补松弛性：} && \lambda_ig_i(\vx)=0 && \text{对所有 }i,\\
&\text{平稳性条件：} && \nabla f(\vx)+\sum_{i=1}^m\lambda_i\nabla g_i(\vx)+\sum_{j=1}^\ell\mu_j\nabla h_j(\vx)=\vzero .
\end{aligned}
$$ (eq-02-03-kkt)

刚才的论证表明：**对任何函数可微且强对偶性成立的问题，每个原问题最优解 $\vx^*$ 和对偶最优解 $(\vlambda^*,\vmu^*)$ 都满足 KKT 条件。** 对凸问题，这些条件也是充分的。

:::{prf:theorem} KKT 条件对凸问题是充分的
:label: thm-02-03-kkt-sufficient

设 {eq}`eq-02-03-primal` 是凸问题，$f$ 和 $g_i$ 可微。若 $(\tilde\vx,\tilde\vlambda,\tilde\vmu)$ 满足 {eq}`eq-02-03-kkt`，则 $\tilde\vx$ 是原问题最优解，$(\tilde\vlambda,\tilde\vmu)$ 是对偶最优解，且 $p^*=d^*$。此外，若 Slater 条件成立，则可行点 $\vx^*$ 是最优解，当且仅当存在乘子 $(\vlambda^*,\vmu^*)$，使其与 $\vx^*$ 一起满足 {eq}`eq-02-03-kkt`。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

由于 $\tilde\lambda_i\ge0$，函数 $L(\vx,\tilde\vlambda,\tilde\vmu)=f(\vx)+\sum_i\tilde\lambda_ig_i(\vx)+ \sum_j\tilde\mu_jh_j(\vx)$ 关于 $\vx$ 是凸的：它是凸函数的非负组合加上仿射函数（$\tilde\mu_j$ 的符号无关紧要，因为 $h_j$ 是仿射的）。平稳性条件表明它的梯度在 $\tilde\vx$ 处为零，所以由 {eq}`eq-02-02-first-order`，$\tilde\vx$ 是 $L(\cdot,\tilde\vlambda,\tilde\vmu)$ 的全局极小点。因此，在最后一步利用互补松弛性和 $\vh(\tilde\vx)=\vzero$，得

$$
\begin{aligned}
q(\tilde\vlambda,\tilde\vmu)&=L(\tilde\vx,\tilde\vlambda,\tilde\vmu)\\
&=f(\tilde\vx)+\sum_i\tilde\lambda_ig_i(\tilde\vx)+\sum_j\tilde\mu_jh_j(\tilde\vx)\\
&=f(\tilde\vx).
\end{aligned}
$$

所以可行点 $\tilde\vx$ 与对偶可行点 $(\tilde\vlambda,\tilde\vmu)$ 的目标函数值相等。由弱对偶性 $f(\tilde\vx)\ge p^*\ge d^*\ge q(\tilde\vlambda,\tilde\vmu)=f(\tilde\vx)$，所以它们全部相等：$\tilde\vx$ 是原问题最优解，$(\tilde\vlambda,\tilde\vmu)$ 是对偶最优解。

在 Slater 条件下，强对偶性成立且对偶最优值可以达到（{prf:ref}`thm-02-03-slater`），因此由定理之前的论证，每个最优解 $\vx^*$ 都与某个对偶最优解一起满足 KKT 条件。
:::

KKT 条件推广了{ref}`sec-02-03-optimality`中得到的条件 $\nabla f(\vx^*)=\vzero$（无约束）和 $\nabla f(\vx^*)=-\mA^\T\vmu$（等式约束）。对凸问题，它们把最优化归结为求解一个由方程和不等式组成的系统。对非凸问题，KKT 条件在约束规范下是必要的，但不再是充分的；{doc}`06-03-kkt-conditions`由 LP 对偶和 Farkas 引理导出这些条件（{prf:ref}`thm-06-03-kkt`），讨论约束规范，并补充二阶条件。这些条件由 Karush（1939）以及 Kuhn 和 Tucker（1951）导出。

:::{prf:example} 带容量限制的调度
:label: ex-02-03-dispatch-capacity

在{prf:ref}`ex-02-03-dispatch`中，设机组 1 的出力至多为 $300$ 兆瓦。用 KKT 条件求解该问题，并解释乘子的含义。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

问题为

$$
\begin{aligned}
\min_{\vx}\ & 8x_1+0.01x_1^2+10x_2+0.02x_2^2\\
\st\ & g(\vx)=x_1-300\le0,\\
& h(\vx)=475-x_1-x_2=0 .
\end{aligned}
$$

它是凸问题，并且 Slater 条件成立，例如取 $\tilde\vx=(250,225)$：$g(\tilde\vx)=-50<0$，$h(\tilde\vx)=0$。所以 KKT 条件是充分必要的。记 $g$ 的乘子为 $\lambda$，$h$ 的乘子为 $\mu$，平稳性条件为

$$
\begin{aligned}
8+0.02x_1+\lambda-\mu&=0,\\
10+0.04x_2-\mu&=0,
\end{aligned}
$$

此外还有 $x_1+x_2=475$、$x_1\le300$、$\lambda\ge0$ 和 $\lambda(x_1-300)=0$。

*情形 $\lambda=0$。* 条件与{prf:ref}`ex-02-03-dispatch`相同，解为 $x_1=350$，违反 $x_1\le300$。这种情形不给出 KKT 点。

*情形 $x_1=300$。* 此时 $x_2=175$，由第二个方程得 $\mu=10+0.04\cdot175=17$，由第一个方程得 $\lambda=\mu-8-0.02\cdot300=17-14=3\ge0$。所有条件都成立。

因此 $\vx^*=(300,175)$ 兆瓦，$\lambda^*=3$，$\mu^*=17$，且

$$
p^*=C_1(300)+C_2(175)=3300+2362.5=5662.5\ \text{\$/h},
$$

比没有容量限制时每小时多 75 美元。对偶函数为 $q(\lambda,\mu)=475\mu-300\lambda-(\mu-\lambda-8)^2/0.04-(\mu-10)^2/0.08$（计算方法同{prf:ref}`ex-02-03-dispatch-dual`），并且 $q(3,17)=8075-900-900-612.5=5662.5=p^*$，与定理的保证一致。

乘子就是价格。$\mu^*=17$ 美元/兆瓦时是需求的边际成本：它等于机组 2 的边际成本，而机组 2 现在是唯一能够响应需求变化的机组。$\lambda^*=3$ 美元/兆瓦时是机组 1 处容量的价值：它是价格 $17$ 与机组 1 的边际成本 $C_1'(300)=14$ 美元/兆瓦时之差；把容量上限从 300 兆瓦提高到 301 兆瓦，成本每小时约降低 3 美元（精确值为 2.97 美元）。由于容量约束是紧约束，互补松弛性成立。
:::

(sec-02-03-sensitivity)=
## 乘子作为灵敏度

乘子的价格解释可以精确化。对 {eq}`eq-02-03-primal` 的约束右端项作摄动：

$$
\begin{aligned}
g_i(\vx)&\le u_i\quad(i=1,\dots,m),\\
h_j(\vx)&=v_j\quad(j=1,\dots,\ell),
\end{aligned}
$$

记 $p^*(\vu,\vv)$ 为摄动后问题的最优值，于是 $p^*(\vzero,\vzero)=p^*$。正的 $u_i$ 放松约束 $i$，负的 $u_i$ 收紧约束 $i$。

:::{prf:theorem} 全局灵敏度不等式
:label: thm-02-03-sensitivity

设未摄动问题的强对偶性成立，且对偶最优解 $(\vlambda^*,\vmu^*)$ 可以达到。则对所有 $\vu$ 和 $\vv$

$$
p^*(\vu,\vv)\ge p^*-\vlambda^{*\T}\vu-\vmu^{*\T}\vv .
$$ (eq-02-03-sensitivity)

此外，若 $p^*(\vu,\vv)$ 在 $(\vzero,\vzero)$ 处可微，则在该处 $\partial p^*/\partial u_i=-\lambda_i^*$，$\partial p^*/\partial v_j=-\mu_j^*$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

设 $\vx$ 对摄动后的问题可行。由于 $\lambda_i^*\ge0$，有 $\lambda_i^*g_i(\vx)\le\lambda_i^*u_i$，且 $\mu_j^*h_j(\vx)=\mu_j^*v_j$，所以

$$
\begin{aligned}
p^*&=q(\vlambda^*,\vmu^*)\le f(\vx)+\vlambda^{*\T}\vg(\vx)+\vmu^{*\T}\vh(\vx)\\
&\le f(\vx)+\vlambda^{*\T}\vu+\vmu^{*\T}\vv .
\end{aligned}
$$

对所有这样的 $\vx$ 取极小即得 {eq}`eq-02-03-sensitivity`。对于导数，取 $\vu=t\ve_i$、$\vv=\vzero$：由不等式得，当 $t>0$ 时 $(p^*(t\ve_i,\vzero)-p^*)/t\ge-\lambda_i^*$，当 $t<0$ 时该商 $\le-\lambda_i^*$，所以若导数存在，则等于 $-\lambda_i^*$。对 $v_j$ 的论证相同。
:::

该不等式是全局成立的，不只适用于小摄动，并且可以从两个方向来理解。若 $\lambda_i^*$ 很大，收紧约束 $i$（$u_i<0$）会使最优成本至少增加 $\lambda_i^*\abs{u_i}$：这个约束代价高昂。若 $\lambda_i^*$ 很小，在一阶意义下放松它（$u_i>0$）只能使成本略微降低。$\lambda_i^*=0$ 的约束在边际上不值得放松，这又是互补松弛性。

在{prf:ref}`ex-02-03-dispatch-capacity`中，容量约束 $x_1-300\le u$ 允许机组 1 出力 $300+u$ 兆瓦，定理预测成本变化约为 $-\lambda^*u=-3u$ 美元/小时；$u=1$ 时的精确变化为每小时 $-2.97$ 美元。需求约束 $475-x_1-x_2=v$ 对应于需求 $475-v$，所以需求增加 1 MW（$v=-1$）使成本每小时约增加 $\mu^*=17$ 美元（精确值为 17.02 美元）。LP 的乘子具有相同的含义，在那里称为影子价格（{prf:ref}`thm-03-05-shadow-price`）；{doc}`03-07-sensitivity-analysis`计算使影子价格保持有效的摄动范围（{ref}`sec-03-07-rhs`）。

(sec-02-03-lp-dual)=
## 线性规划的对偶

利用拉格朗日函数，几行就可以写出任何 LP 的对偶问题，而不必使用{doc}`03-05-lp-duality`中的规则表（{ref}`sec-03-05-dual-rules`）。

:::{prf:example} 由拉格朗日函数导出 LP 的对偶问题
:label: ex-02-03-lp-dual

导出 LP

$$
\begin{aligned}
\min_{\vx}\ & \vc^\T\vx\\
\st\ & \mA\vx\ge\vb,\\
& \vx\ge\vzero ,
\end{aligned}
$$ (eq-02-03-lp-primal)

的对偶问题，其中 $\mA\in\R^{m\times n}$。再对标准形 $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$ 做同样的推导。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*对称形式。* 把约束写成 $\vb-\mA\vx\le\vzero$，乘子为 $\vpi\in\R^m$（沿用线性规划各章的记号），以及 $-\vx\le\vzero$，乘子为 $\vlambda\in\R^n$。拉格朗日函数为

$$
\begin{aligned}
L(\vx,\vpi,\vlambda)&=\vc^\T\vx+\vpi^\T(\vb-\mA\vx)-\vlambda^\T\vx\\
&=\vb^\T\vpi+(\vc-\mA^\T\vpi-\vlambda)^\T\vx .
\end{aligned}
$$

它是 $\vx$ 的仿射函数。仿射函数在 $\R^n$ 上有下界，仅当其线性部分为零，所以

$$
q(\vpi,\vlambda)=
\begin{cases}
\vb^\T\vpi, & \mA^\T\vpi+\vlambda=\vc,\\
-\infty, & \text{其他情形。}
\end{cases}
$$

对偶问题是在 $\vpi\ge\vzero$、$\vlambda\ge\vzero$ 上极大化 $q$。把隐含约束写成显式约束，并消去 $\vlambda=\vc-\mA^\T\vpi\ge\vzero$，得

$$
\begin{aligned}
\max_{\vpi}\ & \vb^\T\vpi\\
\st\ & \mA^\T\vpi\le\vc,\\
& \vpi\ge\vzero ,
\end{aligned}
$$ (eq-02-03-lp-dual)

这就是{doc}`03-05-lp-duality`中对称形式的原始–对偶问题对（{prf:ref}`def-03-05-canonical-pair`）。被消去的乘子 $\lambda_j=c_j-\mA_{\cdot j}^\T\vpi$ 正是单纯形法中的**检验数**（reduced cost）（{eq}`eq-03-03-reduced-cost`），而互补松弛性 {eq}`eq-02-03-complementary-slackness` 变为

$$
\begin{aligned}
\pi_i\,(\mA_{i\cdot}\vx-b_i)&=0,\quad i=1,\dots,m,\\
(c_j-\mA_{\cdot j}^\T\vpi)\,x_j&=0,\quad j=1,\dots,n :
\end{aligned}
$$

即有松弛的约束的对偶价格为零，取正值的变量的检验数为零。

*标准形。* 取 $\vh(\vx)=\vb-\mA\vx=\vzero$，乘子为 $\vpi$（现在符号不限），以及 $-\vx\le\vzero$，乘子 $\vlambda\ge\vzero$。同样的计算给出：若 $\mA^\T\vpi+\vlambda=\vc$ 则 $q=\vb^\T\vpi$，否则为 $-\infty$。所以对偶问题为 $\max\{\vb^\T\vpi:\mA^\T\vpi\le\vc\}$，$\vpi$ 为自由变量。
:::

对{prf:ref}`ex-02-03-lp-dual`的解答中导出的对偶问题对 {eq}`eq-02-03-lp-primal`–{eq}`eq-02-03-lp-dual`，所有约束都是仿射的，因此 Slater 条件退化为可行性（{ref}`sec-02-03-strong-duality`）：若原问题可行且有界，则 $d^*=p^*$，且两个最优值都能达到。弱对偶性（对可行的 $\vx$ 和 $\vpi$ 有 $\vb^\T\vpi\le\vc^\T\vx$）也可以直接看出：$\vc^\T\vx\ge(\mA^\T\vpi)^\T\vx=\vpi^\T\mA\vx\ge\vpi^\T\vb$，其中第一个不等式用到 $\vx\ge\vzero$，第二个用到 $\vpi\ge\vzero$。这正是{prf:ref}`ex-01-03-water-supply`中对供水模型手工求得的下界。

(sec-02-03-least-norm)=
## 线性方程组的最小范数解

:::{prf:example} 最小范数解
:label: ex-02-03-least-norm

欠定线性方程组 $\mA\vx=\vb$（$\mA\in\R^{m\times n}$ 行满秩，秩为 $m<n$）有无穷多个解。通过对偶问题求出其中欧氏范数最小的解，

$$
\begin{aligned}
\min_{\vx}\ & \vx^\T\vx\\
\st\ & \mA\vx=\vb ,
\end{aligned}
$$

并对

$$
\begin{aligned}
\mA&=\begin{pmatrix}1&1&1\\1&2&3\end{pmatrix},\\
\vb&=\begin{pmatrix}6\\14\end{pmatrix}.
\end{aligned}
$$

求出其数值。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*对偶函数。* 取 $\vh(\vx)=\mA\vx-\vb$，乘子 $\vmu\in\R^m$，拉格朗日函数 $L(\vx,\vmu)=\vx^\T\vx+\vmu^\T(\mA\vx-\vb)$ 是 $\vx$ 的严格凸二次函数。它的梯度 $2\vx+\mA^\T\vmu$ 在 $\vx=-\tfrac12\mA^\T\vmu$ 处为零，代入得

$$
\begin{aligned}
q(\vmu)&=\tfrac14\vmu^\T\mA\mA^\T\vmu-\tfrac12\vmu^\T\mA\mA^\T\vmu-\vb^\T\vmu\\
&=-\tfrac14\vmu^\T\mA\mA^\T\vmu-\vb^\T\vmu ,
\end{aligned}
$$

这是一个凹二次函数（由于 $\mA$ 行满秩，$\mA\mA^\T\succ0$）。

*对偶问题。* 在 $\vmu\in\R^m$ 上极大化 $q$：梯度 $-\tfrac12\mA\mA^\T\vmu-\vb$ 在 $\vmu^*=-2(\mA\mA^\T)^{-1}\vb$ 处为零，且

$$
d^*=q(\vmu^*)=\vb^\T(\mA\mA^\T)^{-1}\vb .
$$

约束是仿射的且可行，因此强对偶性成立。原问题的解是 $L(\cdot,\vmu^*)$ 的极小点：

$$
\vx^*=-\tfrac12\mA^\T\vmu^*=\mA^\T(\mA\mA^\T)^{-1}\vb ,
$$ (eq-02-03-least-norm)

确实有 $\norm{\vx^*}^2=\vb^\T(\mA\mA^\T)^{-1}\mA\mA^\T(\mA\mA^\T)^{-1}\vb=d^*$。

*数值。* 这里 $\mA\mA^\T=\begin{pmatrix}3&6\\6&14\end{pmatrix}$，行列式为 $42-36=6$，所以

$$
(\mA\mA^\T)^{-1}\vb=\frac16\begin{pmatrix}14&-6\\-6&3\end{pmatrix}
\begin{pmatrix}6\\14\end{pmatrix}
=\frac16\begin{pmatrix}0\\6\end{pmatrix}
=\begin{pmatrix}0\\1\end{pmatrix}.
$$

因此 $\vx^*=\mA^\T(0,1)^\T=(1,2,3)$，$\vmu^*=(0,-2)$，$p^*=1+4+9=14=d^*=\vb^\T(0,1)^\T$。*验算：* $1+2+3=6$，$1+4+9=14$，所以 $\mA\vx^*=\vb$。其他每个解都具有 $\vx^*+\vz$ 的形式，其中 $\mA\vz=\vzero$；由于 $\vx^*$ 位于 $\mA$ 的行空间中，当 $\vz\ne\vzero$ 时 $\norm{\vx^*+\vz}^2=\norm{\vx^*}^2+\norm{\vz}^2>14$。
:::

同样的计算适用于任何带线性不等式约束的凸 QP。对 $\min\{\vx^\T\mP\vx:\mA\vx\le\vb\}$（$\mP\succ0$），拉格朗日函数 $\vx^\T\mP\vx+\vlambda^\T(\mA\vx-\vb)$ 在 $\vx=-\tfrac12\mP^{-1}\mA^\T\vlambda$ 处取极小，且

$$
q(\vlambda)=-\tfrac14\vlambda^\T\mA\mP^{-1}\mA^\T\vlambda-\vb^\T\vlambda ,
$$

需在 $\vlambda\ge\vzero$ 上极大化。对偶问题仍是 QP，每个约束对应一个变量；当约束个数少于变量个数时，这一点很有吸引力。

(sec-02-03-duality-gap)=
## 强对偶性不成立的情形

强对偶性可能以两种方式不成立：问题不是凸的；或者问题是凸的，但不满足任何约束规范。

:::{prf:example} 对偶间隙：只能开或关的水泵
:label: ex-02-03-duality-gap

某泵站有一台水泵，它要么关闭（$x=0$），要么满负荷运行（$x=1$，以其容量为单位）；运行费用为 $1$（以其每日运行费用为单位）。泵站必须至少输送其容量的一半。比较下面问题的最优值与对偶最优值：

$$
\begin{aligned}
\min_{x}\ & x\\
\st\ & g(x)=0.5-x\le0,\\
& x\in\mathcal D=\{0,1\}.
\end{aligned}
$$
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

*原问题。* 只有 $x=1$ 可行，所以 $p^*=1$。

*对偶问题。* 拉格朗日函数 $L(x,\lambda)=x+\lambda(0.5-x)$ 在 $\mathcal D$ 的两个点处取值：

$$
q(\lambda)=\min\{L(0,\lambda),\ L(1,\lambda)\}
=\min\{0.5\lambda,\ 1-0.5\lambda\}.
$$

第一项关于 $\lambda$ 递增，第二项关于 $\lambda$ 递减；两者在 $\lambda=1$ 处相等，所以 $d^*=q(1)=0.5$。对偶间隙为 $p^*-d^*=0.5$。

该问题不是凸问题，因为它的定义域 $\{0,1\}$ 不是凸集。在{numref}`fig-02-03-duality-geometry`(b) 的图景中，$\mathcal G=\{(0.5,0),(-0.5,1)\}$ 只有两个点，没有一条经过 $(0,1)$ 的直线能支撑它。对偶值 $0.5$ 是*松弛*问题的最优值，在松弛问题中 $x$ 可以取 $[0,1]$（即 $\mathcal D$ 的凸包）中的任何值：对偶问题“看不到”整数性。这让我们初步看到，为什么整数规划（{doc}`07-01-logical-constraints`）比凸问题更难，以及为什么它们的松弛问题能给出界。
:::

没有严格可行点的**凸**问题仍可能满足强对偶性，但对偶最优值可能达不到。对 $\min\{x:x^2\le0\}$，唯一的可行点是 $x=0$，所以 $p^*=0$，而 Slater 条件不成立。对偶函数为：当 $\lambda>0$ 时 $q(\lambda)=\inf_x(x+\lambda x^2)=-1/(4\lambda)$，当 $\lambda=0$ 时为 $-\infty$；其上确界为 $d^*=0=p^*$，当 $\lambda\to\infty$ 时趋近但永远达不到。不存在 KKT 点，因为平稳性条件 $1+2\lambda x=0$ 在 $x=0$ 处不成立。真正存在对偶间隙的凸问题也是有的，但需要特殊的定义域：对定义在 $\mathcal D=\{\vx:x_2>0\}$ 上的 $\min\{e^{-x_1}:x_1^2/x_2\le0\}$，可行点满足 $x_1=0$，所以 $p^*=1$，而 $d^*=0$（Boyd & Vandenberghe, 2004, Exercise 5.21）。

## 本章回顾

### 小结

- 凸问题是在凸约束 $g_i\le0$ 和仿射约束 $h_j=0$ 下极小化凸函数 $f$（{prf:ref}`def-02-03-convex-problem`）；其局部极小点是全局极小点（{prf:ref}`thm-02-03-local-global`），并且 $\vx^*$ 是最优解当且仅当对所有可行点 $\vy$ 有 $\nabla f(\vx^*)^\T(\vy-\vx^*)\ge0$（{prf:ref}`thm-02-03-optimality`）。LP 以及 $\mQ\succeq0$ 的 QP 都是凸问题。
- 拉格朗日函数 {eq}`eq-02-03-lagrangian` 为约束定价；对偶函数 $q$（{prf:ref}`def-02-03-lagrangian`）是凹函数，并且当 $\vlambda\ge\vzero$ 时是 $p^*$ 的下界（{prf:ref}`thm-02-03-lower-bound`）。
- 对偶问题 {eq}`eq-02-03-dual` 总是凸问题，且 $d^*\le p^*$（弱对偶性，{prf:ref}`thm-02-03-weak-duality`）。对有严格可行点的凸问题，$d^*=p^*$ 且对偶最优值可以达到（Slater，{prf:ref}`thm-02-03-slater`）；仿射约束只需可行。
- 在强对偶性下，$\vx^*$ 极小化 $L(\cdot,\vlambda^*,\vmu^*)$，且 $\lambda_i^*g_i(\vx^*)=0$（{eq}`eq-02-03-complementary-slackness`）；加上可微性，就得到 KKT 条件 {eq}`eq-02-03-kkt`，它们对凸问题是充分的，在 Slater 条件下也是必要的（{prf:ref}`thm-02-03-kkt-sufficient`）。
- 乘子是价格，也是灵敏度：把约束 $i$ 摄动 $u_i$，最优值的变化至少为 $-\lambda_i^*u_i$（{prf:ref}`thm-02-03-sensitivity`）；在经济调度中，$\mu^*$ 是需求的边际成本，$\lambda^*$ 是额外容量的价值（{prf:ref}`ex-02-03-dispatch-dual`、{prf:ref}`ex-02-03-dispatch-capacity`）。
- LP 的对偶问题可由拉格朗日函数导出（{prf:ref}`ex-02-03-lp-dual`）；最小范数解为 $\mA^\T(\mA\mA^\T)^{-1}\vb$（{prf:ref}`ex-02-03-least-norm`）；没有凸性时可能出现对偶间隙（{prf:ref}`ex-02-03-duality-gap`、{numref}`fig-02-03-duality-geometry`）。

### 常见错误

- **在凸问题中允许非线性等式约束。** $h_j$ 必须是仿射的。
- **忘记符号约束 $\vlambda\ge\vzero$。** 不等式 $g_i\le0$ 的乘子非负；等式的乘子是自由的。写成 $g_i\ge0$ 的约束必须先变号。
- **在可行集上极小化拉格朗日函数。** 对偶函数是在定义域 $\mathcal D$ 上、不带约束地取极小，正因如此它才容易计算。
- **不检验约束规范就假定强对偶性成立。** 弱对偶性总是成立；强对偶性则需要例如凸性加 Slater 条件。
- **把非凸问题的 KKT 点当作全局极小点。** KKT 条件的充分性（{prf:ref}`thm-02-03-kkt-sufficient`）用到了凸性。
- **把互补松弛性的方向弄反。** $\lambda_i^*>0$ 迫使约束成为有效约束，但有效约束的 $\lambda_i^*$ 可能为 $0$。

## 参考文献

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 4, Sections
  4.2–4.3: KKT conditions and their sufficiency for convex problems;
  Chapter 6, Sections 6.1–6.2: Lagrangian duality, duality theorems and
  saddle-point optimality.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Section 4.1, motivation of LP duality
  through Lagrange multipliers.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Chapter 4,
  Sections 4.1–4.4: convex optimization problems, LP and QP; Chapter 5,
  Sections 5.1–5.6: the Lagrange dual function and problem, geometric
  interpretation, optimality conditions, perturbation and sensitivity
  analysis; Exercise 5.21, a convex problem without strong duality.)
- Karush, W. (1939). *Minima of functions of several variables with
  inequalities as side constraints* [Master's thesis, University of
  Chicago].
- Kuhn, H. W., & Tucker, A. W. (1951). Nonlinear programming. In J. Neyman
  (Ed.), *Proceedings of the Second Berkeley Symposium on Mathematical
  Statistics and Probability* (pp. 481–492). University of California
  Press.
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming*
  (4th ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3>
  (Chapter 14, duality and dual methods.)
