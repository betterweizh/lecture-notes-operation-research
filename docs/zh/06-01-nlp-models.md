# 非线性规划模型

## 概述

**非线性规划问题**（nonlinear program, NLP）是目标函数或约束为决策变量的非线性函数的最优化问题。工程、经济和统计中的许多模型都属于这一类：用最小二乘法使模型拟合数据、在拥堵网络上选择流量、确定结构尺寸、为产品定价等。本章建立一般模型（{eq}`eq-06-01-nlp`），把一个最小二乘估计问题表述为 NLP 并进行数值求解（{prf:ref}`ex-06-01-consumer-regression`），比较 NLP 与线性规划（linear programming, LP），并介绍后续各章分别处理的几类 NLP：无约束问题、线性约束问题、二次规划和凸规划。随后定义局部极小点与全局极小点（{prf:ref}`def-06-01-local-global`），回顾极小点在何种条件下必定存在，并列出使 NLP 比 LP 更难的若干困难。刻画极小点的最优性条件见{doc}`06-02-optimality-conditions`和{doc}`06-03-kkt-conditions`；算法见{doc}`06-04-line-search`至{doc}`06-06-constrained-algorithms`。本章用到{doc}`01-01-linear-algebra`和{doc}`01-02-matrix-calculus`中的梯度、黑塞矩阵和二次型，以及{doc}`02-01-convex-sets`和{doc}`02-02-convex-functions`中的紧集、凸集和凸函数。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 把最优化问题写成含不等式约束 $g_i(\vx)\le0$ 和等式约束 $h_j(\vx)=0$ 的 NLP 标准形式；
- 把参数估计问题表述为带约束的最小二乘 NLP，并解释其解；
- 说明为什么与线性规划不同，NLP 的最优解不一定位于可行集的顶点；
- 识别无约束问题、线性约束问题、二次规划和凸规划，并把二次目标函数写成 $\vc^\T\vx+\tfrac12\vx^\T\mQ\vx$ 的形式，其中 $\mQ$ 对称；
- 区分局部极小点与全局极小点，并用魏尔斯特拉斯定理证明极小点存在；
- 说出非线性规划的主要困难，并在小例子中识别它们。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $m$, $\ell$ | NLP 的不等式约束个数与等式约束个数 |
| $\vg(\vx)$, $\vh(\vx)$ | 向量 $(g_1(\vx),\dots,g_m(\vx))^\T$ 与 $(h_1(\vx),\dots,h_\ell(\vx))^\T$ |
| $B(\vx^*,\varepsilon)$ | 开球 $\{\vx:\norm{\vx-\vx^*}<\varepsilon\}$ |
| $t$, $E_t$ | {prf:ref}`ex-06-01-consumer-regression`中的月份序号与月预算（美元） |
| $p_{ti}$, $Q_{ti}$ | 第 $t$ 月产品 $i$ 的价格（美元/单位）与购买量（单位） |
| $q_i$, $\beta_i$ | 产品 $i$ 的必需购买量（单位）与预算份额（待估参数） |
| $r_{ti}$ | 第 $t$ 月产品 $i$ 的残差（美元） |

(sec-06-01-general-form)=
## 一般的非线性规划问题

NLP 选取 $n$ 个实决策变量 $\vx=(x_1,\dots,x_n)^\T$，在不等式约束和等式约束下使目标函数最小：

$$
\begin{aligned}
\min_{\vx\in\R^n}\ & f(\vx)\\
\st\ & g_i(\vx)\le0, && i=1,\dots,m,\\
& h_j(\vx)=0, && j=1,\dots,\ell.
\end{aligned}
$$ (eq-06-01-nlp)

其中函数 $f,g_1,\dots,g_m,h_1,\dots,h_\ell:\R^n\to\R$ 是给定的。用向量记号，约束写作 $\vg(\vx)\le\vzero$ 和 $\vh(\vx)=\vzero$，问题即为 $\min\{f(\vx):\vx\in\Omega\}$，其中**可行集**（feasible set）为

$$
\Omega=\{\vx\in\R^n:\ \vg(\vx)\le\vzero,\ \vh(\vx)=\vzero\}.
$$ (eq-06-01-feasible-set)

$\Omega$ 中的点称为**可行解**（feasible solution）；目标函数值最小的可行点称为**最优解**（optimal solution）或**极小点**（minimum），记作 $\vx^*$。在关于 NLP 的这一部分中，我们始终假设 $f$、$g_i$ 和 $h_j$ **连续可微**（continuously differentiable），在用到二阶导数时假设它们二阶连续可微。正是光滑性保证了梯度存在，从而{doc}`06-02-optimality-conditions`中的最优性条件可以使用。

形式 {eq}`eq-06-01-nlp` 经过简单改写即可涵盖各种变体：

- 最大化问题 $\max f(\vx)$ 与 $\min\,-f(\vx)$ 有相同的最优解，最优值满足 $\max f=-\min(-f)$；
- 约束 $g(\vx)\ge b$ 改写为 $b-g(\vx)\le0$，约束 $g(\vx)\le b$ 改写为 $g(\vx)-b\le0$；
- 等式 $h(\vx)=0$ 等价于一对不等式 $h(\vx)\le0$，$-h(\vx)\le0$，{doc}`06-03-kkt-conditions`用这一技巧推导等式约束的乘子；
- 符号约束 $x_j\ge0$ 即不等式 $-x_j\le0$。

不等式采用“$\le0$”的约定，这关系到{doc}`06-03-kkt-conditions`中乘子的符号；这一部分的每个模型都先化为这种形式。

在 NLP 各章中目标函数始终记为 $f$，而 $z$ 保留给线性规划的目标函数值。下面的例子说明同一个目标函数在三个不同可行集上的表现。

:::{prf:example} 一个目标函数，三个可行集
:label: ex-06-01-three-feasible-sets

设 $f(\vx)=(x_1-2)^2+(x_2-2)^2$，即到点 $(2,2)$ 的距离的平方。考虑三个问题：

1. 在 $\Omega_1=\R^2$ 上求 $\min f(\vx)$（无约束）；
2. $\min f(\vx)$，满足 $x_1+x_2-5\le0$，$x_1\ge0$，$x_2\ge0$；
3. $\min f(\vx)$，满足 $(x_1-1)^2+(x_2-1)^2-1\le0$，$x_1+x_2\ge1$，$x_1\ge0$，$x_2\ge0$。

求每个问题的最优解和最优值，并指出哪些约束是紧约束。
:::

::::{dropdown} 解答
:class-container: solution
:animate: fade-in

*问题 1。* $f$ 处处满足 $f\ge0$，且 $f(2,2)=0$，因此 $\vx^*=(2,2)$，$f^*=0$。

*问题 2。* 可行集是以 $(0,0)$、$(5,0)$、$(0,5)$ 为顶点的三角形。点 $(2,2)$ 满足 $2+2=4<5$，位于三角形内部，因此无约束极小点是可行的，从而是最优的：$\vx^*=(2,2)$，$f^*=0$。没有紧约束。

*问题 3。* 第一个约束是以 $(1,1)$ 为圆心、半径为 1 的闭圆盘；它与坐标轴相切于 $(1,0)$ 和 $(0,1)$，因此符号约束是冗余的。第二个约束切去圆盘位于直线 $x_1+x_2=1$ 下方的一小角。点 $(2,2)$ 位于圆盘之外，到圆心的距离为 $\sqrt2$。在圆盘上使到 $(2,2)$ 的距离最小，得到的是从圆心指向 $(2,2)$ 的线段与圆盘边界的交点：

$$
\begin{aligned}
\vx^*&=(1,1)+\tfrac{1}{\sqrt2}(1,1)\\
&=\Bigl(1+\tfrac{1}{\sqrt2},\ 1+\tfrac{1}{\sqrt2}\Bigr)\approx(1.7071,\ 1.7071),
\end{aligned}
$$

它满足 $x_1+x_2\approx3.414\ge1$。其目标函数值为 $f^*=2\bigl(1-\tfrac1{\sqrt2}\bigr)^2=3-2\sqrt2\approx0.1716$。此时最优解位于弯曲的边界上，只有圆盘约束是紧约束。数值求解器（SLSQP，一个序列二次规划程序；Kraft, 1988）给出同一个点。

```{figure} figures/06-01-three-feasible-sets.svg
:name: fig-06-01-three-feasible-sets
:width: 680px
:class: dark-light
:alt: 三幅子图，显示到点 (2,2) 的距离平方的圆形等值线；第二幅图显示包含 (2,2) 的三角形，第三幅图显示以 (1,1) 为圆心、切去一角的圆盘，最优解位于圆盘边界上

$f(\vx)=(x_1-2)^2+(x_2-2)^2$ 在 {prf:ref}`ex-06-01-three-feasible-sets` 的三个可行集上的等值线。在问题 1 和问题 2 中，最优解（橙色）是无约束极小点 $(2,2)$，位于 $\Omega$ 的内部；在问题 3 中，最优解是最小的那条等值线与圆盘相切的点。三个最优解都不是顶点。
```
::::

(sec-06-01-regression)=
## 建立 NLP 模型：最小二乘估计

统计估计是 NLP 最常见的来源之一：选取模型的未知参数，使模型预测值与观测数据之间的某种失配度量最小。当预测值非线性地依赖于参数，且参数还须满足附加条件时，就得到一个带约束的 NLP。

:::{prf:example} 估计消费者的支出模型
:label: ex-06-01-consumer-regression

某家庭每月将 $E$ 美元的预算用于购买三种产品 $i=1,2,3$，价格为 $p_i$（美元/单位）。一个经典的消费者需求模型（线性支出系统）假设：家庭先购买每种产品的必需量 $q_i\ge0$，再把剩余的可自由支配金额 $E-\sum_j p_jq_j$ 按固定份额 $\beta_i\ge0$ 分配给各产品，其中 $\beta_1+\beta_2+\beta_3=1$。份额 $\beta_i$ 是该家庭把可自由支配金额用于产品 $i$ 的边际支出倾向。模型预测的产品 $i$ 支出为

$$
p_iq_i+\beta_i\Bigl(E-\sum_{j=1}^3p_jq_j\Bigr),\qquad i=1,2,3.
$$ (eq-06-01-les)

参数 $\vq=(q_1,q_2,q_3)^\T$ 和 $\vbeta=(\beta_1,\beta_2,\beta_3)^\T$ 需要根据四个月的记录来估计：

```{list-table} 第 $t=1,\dots,4$ 月的月预算 $E_t$（美元）、价格 $p_{ti}$（美元/单位）和购买量 $Q_{ti}$（单位）。
:name: tab-06-01-purchase-records
:header-rows: 1

* - 月份 $t$
  - 1
  - 2
  - 3
  - 4
* - $E_t$（美元）
  - 1500
  - 2000
  - 1500
  - 1200
* - $p_{t1}$（美元/单位）
  - 2.5
  - 2.7
  - 2.8
  - 2.6
* - $Q_{t1}$（单位）
  - 198
  - 200
  - 100
  - 152
* - $p_{t2}$（美元/单位）
  - 4.0
  - 3.5
  - 4.0
  - 3.6
* - $Q_{t2}$（单位）
  - 124
  - 200
  - 80
  - 115
* - $p_{t3}$（美元/单位）
  - 5.0
  - 4.9
  - 5.6
  - 4.5
* - $Q_{t3}$（单位）
  - 100
  - 110
  - 60
  - 65
```

把这个估计问题表述为 NLP，求解之，并解释估计结果。
:::

::::{dropdown} 解答
:class-container: solution
:animate: fade-in

**决策变量。** 六个模型参数 $\beta_1,\beta_2,\beta_3,q_1,q_2,q_3$。数据 $E_t$、$p_{ti}$ 和 $Q_{ti}$ 虽然在统计意义上是“变量”，但在这个最优化问题中是常数。

**目标函数。** 第 $t$ 月产品 $i$ 的实际支出为 $p_{ti}Q_{ti}$，而 {eq}`eq-06-01-les` 根据该月的预算和价格对它作出预测。残差为

$$
r_{ti}(\vbeta,\vq)=p_{ti}Q_{ti}-p_{ti}q_i-\beta_i\Bigl(E_t-\sum_{j=1}^3p_{tj}q_j\Bigr)
$$ (eq-06-01-residual)

其中 $t=1,\dots,4$，$i=1,2,3$，单位为美元。最小二乘法使十二个残差的平方和最小：

$$
f(\vbeta,\vq)=\sum_{t=1}^4\sum_{i=1}^3r_{ti}(\vbeta,\vq)^2 .
$$ (eq-06-01-sse)

**约束。** 份额非负且和为 1；必需购买量非负，并且不能超过任何一个月的实际购买量，即对所有 $t$ 有 $q_i\le Q_{ti}$。利用残差 {eq}`eq-06-01-residual`，完整模型为

$$
\begin{aligned}
\min_{\vbeta,\vq}\ & f(\vbeta,\vq)=\sum_{t=1}^4\sum_{i=1}^3r_{ti}(\vbeta,\vq)^2\\
\st\ & \beta_1+\beta_2+\beta_3=1,\\
& q_i\le Q_{ti}, && i=1,2,3,\ t=1,\dots,4,\\
& \vbeta\ge\vzero,\\
& \vq\ge\vzero .
\end{aligned}
$$ (eq-06-01-regression-nlp)

写成 {eq}`eq-06-01-nlp` 的形式，它有 $n=6$ 个变量、一个等式约束 $h_1=\beta_1+\beta_2+\beta_3-1$ 和 $12+6=18$ 个不等式约束。每个 $q_i$ 只有最紧的上界起作用，即 $q_i\le\min_tQ_{ti}$，也就是 $q_1\le100$、$q_2\le80$、$q_3\le60$ 单位；其余九个上界都是冗余的。

**为什么问题是非线性且非凸的。** 残差中含有乘积 $\beta_iq_j$，因此它关于 $(\vbeta,\vq)$ 是双线性的，$f$ 是四次多项式。固定 $\vbeta$ 时，残差关于 $\vq$ 是仿射的，$f$ 是 $\vq$ 的凸二次函数；固定 $\vq$ 时，$f$ 是 $\vbeta$ 的凸二次函数。但是就两组变量整体而言，$f$ 不是凸的：它在下文求得的解处的黑塞矩阵有一个负特征值（$\approx-4.55$）。约束是线性的，因此这是一个目标函数非凸的**线性约束** NLP，求解器只能保证得到局部极小点。

**求解。** 用序列二次规划方法（SLSQP；Kraft, 1988；见{doc}`06-06-constrained-algorithms`及{doc}`09-06-code-nlp`中的代码）从 1500 个随机可行点出发求解。凡是停止于满足{doc}`06-03-kkt-conditions`一阶最优性条件的点的运行，都停在同一个点上，且没有一次运行得到更小的目标函数值：

$$
\begin{aligned}
\vbeta^*&\approx(0.2474,\ 0.4635,\ 0.2891),\\
\vq^*&\approx(89.85,\ 16.67,\ 37.32)\ \text{单位},\\
f^*&\approx1.3547\times10^{5}\ \$^2 .
\end{aligned}
$$

在 $(\vbeta^*,\vq^*)$ 处所有界约束都是非有效约束；只有等式 $\sum_i\beta_i=1$ 是紧的。$f$ 在平面 $\sum_i\beta_i=1$ 上的限制在该点的黑塞矩阵正定，因此由{doc}`06-03-kkt-conditions`的二阶条件，该点是严格局部极小点。

**解释。** 据估计，该家庭每月在这三种产品上的必需支出约为 $\sum_jp_{tj}q^*_j\approx\$462$–$\$527$（金额随价格变化），其余部分按 $25\,\%:46\,\%:29\,\%$ 的比例支出。{numref}`fig-06-01-regression-fit`比较了拟合支出与观测值。均方根残差为 $\sqrt{f^*/12}\approx\$106$，失配集中在第 3 月，该月三个残差都介于 $-\$154$ 与 $-\$212$ 之间。

原因在于模型本身，而不在于求解器。把 {eq}`eq-06-01-les` 对 $i$ 求和并利用 $\sum_i\beta_i=1$ 可知，第 $t$ 月的预测总支出总是恰好等于 $E_t$。而观测到的总支出 $\sum_ip_{ti}Q_{ti}$ 分别为 1491 美元、1779 美元、936 美元和 1101.70 美元，因此第 3 月该家庭在这些产品上只花掉了 1500 美元预算中的 936 美元，任何参数选择都无法再现这一点。如果用观测总支出代替 $E_t$，同样的方法给出 $f^*\approx1.241\times10^4\ \$^2$，小了一个数量级，此时 $\vbeta^*\approx(0.306,0.416,0.277)$，$\vq^*\approx(27.08,\,0,\,19.57)$：现在界约束 $q_2\ge0$ 是紧的。选择哪些数据进入模型，本身就是建模的一部分。

```{figure} figures/06-01-regression-fit.svg
:name: fig-06-01-regression-fit
:width: 520px
:class: dark-light
:alt: 三种产品四个月的预测月支出对观测月支出的散点图，附 45 度线；第 3 月的三个点远在直线上方

在 {prf:ref}`ex-06-01-consumer-regression` 的最小二乘估计下，十二个产品–月份组合的观测支出 $p_{ti}Q_{ti}$ 与拟合支出的对比。位于虚线上的点被精确拟合。第 3 月的三个点（已标注）远在直线上方：模型迫使预测支出之和等于预算 $E_3=\$1500$，而实际只花了 936 美元。
```
::::

(sec-06-01-lp-vs-nlp)=
## 线性规划与非线性规划

线性规划是 {eq}`eq-06-01-nlp` 中 $f$、每个 $g_i$ 和每个 $h_j$ 都是仿射函数的特殊情形。比较这两类问题可以说明为什么 NLP 需要不同的理论和不同的算法。

```{list-table} 线性规划与非线性规划的比较。
:name: tab-06-01-lp-vs-nlp
:header-rows: 1

* - 性质
  - 线性规划
  - 非线性规划
* - 可行集
  - 多面体（凸的，顶点个数有限）
  - 由光滑的不等式和方程描述的任意集合；可能非凸或不连通
* - 最优解的位置
  - 若存在，则可在顶点（基本可行解）处取得
  - 任何位置：内部、弯曲的边界上或顶点处
* - 局部与全局
  - 每个局部极小点都是全局极小点
  - 除非问题是凸的，局部极小点不一定是全局极小点
* - 存在性
  - 非空可行集上有界的目标函数总能达到最优值
  - 有界的目标函数可能达不到其下确界
* - 最优性检验
  - 检验数 $\bar c_j\ge0$，可精确检验
  - Karush–Kuhn–Tucker（KKT）梯度条件，它是必要条件，但一般不是充分条件
* - 算法
  - 单纯形法：有限次主元变换，可以使用精确算术
  - 迭代方法，在极限意义下收敛到局部极小点或驻点
* - 对偶性
  - 只要最优解存在，强对偶性总成立
  - 对凸问题在约束规范下成立
```

第二行是最明显的差别。线性规划基本定理（{prf:ref}`thm-03-02-fundamental`）表明，线性目标函数在多面体上的最小值在某个顶点处取得，这正是单纯形法可以只在顶点中搜索的原因。非线性目标函数的等值线是弯曲的，与可行集相切的最小等值线可能在任何地方与之相切：{prf:ref}`ex-06-01-three-feasible-sets`中有两个内部最优解和一个位于弯曲边界上的最优解。因此，NLP 算法在可行点构成的连续统中移动，由导数引导。

(sec-06-01-classes)=
## 非线性规划的分类

算法和理论都要针对 $f$ 和 $\Omega$ 的结构来设计。下面几类问题在全书中反复出现。

**无约束问题**（unconstrained problem）。$\Omega=\R^n$：完全没有约束，例如 $\min x_1^2+x_2^2-1$ 或{prf:ref}`ex-06-01-three-feasible-sets`的问题 1。约束在最优解处都是非有效约束的问题，在局部上表现得像无约束问题。其最优性条件最简单（梯度为零），可以使用{doc}`06-05-unconstrained-algorithms`中的方法。

**约束问题**（constrained problem）。$\Omega\subsetneq\R^n$。重要的子类如下。

- **界约束问题**（bound-constrained problem）：只有 $l_j\le x_j\le u_j$。
- **线性约束问题**（linearly constrained problem）：$\Omega=\{\vx:\mA\vx\le\vb,\ \mA'\vx=\vb'\}$ 是多面体，只有 $f$ 是非线性的，如{prf:ref}`ex-06-01-consumer-regression`。{doc}`06-06-constrained-algorithms`中的 Frank–Wolfe 方法利用这一结构，在每次迭代中求解一个 LP。
- **二次规划**（quadratic programming, QP）：多面体上的二次目标函数，

  $$
  \begin{aligned}
  \min_{\vx}\ & f(\vx)=\vc^\T\vx+\tfrac12\vx^\T\mQ\vx\\
  \st\ & \mA\vx=\vb,\\
  & \vx\ge\vzero,
  \end{aligned}
  $$ (eq-06-01-qp)

  其中 $\vc\in\R^n$，$\mQ\in\Sym^n$ 是**对称**矩阵。要求对称不失一般性：由于 $\vx^\T\mQ\vx=\vx^\T\tfrac12(\mQ+\mQ^\T)\vx$，任何矩阵都可以换成它的对称部分而不改变 $f$。当 $\mQ$ 对称时，$\nabla f(\vx)=\vc+\mQ\vx$，$\nabla^2f=\mQ$（{numref}`tab-01-02-vector-derivatives`），对角元 $Q_{jj}$ 是 $x_j^2$ 系数的两倍，而非对角元 $Q_{jk}=Q_{kj}$ 都等于 $x_jx_k$ 的系数。目标函数为凸函数当且仅当 $\mQ\succeq\vzero$，此时该 QP 是凸规划。
- **线性规划**：{eq}`eq-06-01-qp` 中 $\mQ=\vzero$。LP 是第三部分{doc}`03-01-lp-model-building`至{doc}`03-07-sensitivity-analysis`处理的特殊情形。
- **凸规划**（convex program）：$f$ 是凸函数且 $\Omega$ 是凸集。对于 {eq}`eq-06-01-nlp`，当 $f$ 和每个 $g_i$ 都是凸函数且每个 $h_j$ 都是**仿射**函数时，这一条件成立。（凸的不等式约束函数给出凸的下水平集 $\{g_i\le0\}$；而非线性等式 $h_j(\vx)=0$ 一般不定义凸集，即使 $h_j$ 是凸函数也是如此，圆周 $x_1^2+x_2^2=1$ 就是一例。）凸规划是{doc}`02-03-convex-optimization-duality`的主题。

要把二次目标函数写成 QP 形式，只需读出各项系数。对于

$$
f(\vx)=x_1+3x_2+x_1^2-5x_1x_2+2x_2^2
$$

由线性部分得 $\vc=(1,3)^\T$，由二次部分得 $Q_{11}=2\cdot1=2$，$Q_{22}=2\cdot2=4$，$Q_{12}=Q_{21}=-5$：

$$
\begin{aligned}
\mQ&=\begin{pmatrix}2&-5\\-5&4\end{pmatrix},\\
\tfrac12\vx^\T\mQ\vx&=x_1^2-5x_1x_2+2x_2^2 .
\end{aligned}
$$

其特征值为 $\lambda_{1,2}(\mQ)=3\mp\sqrt{26}\approx-2.099,\ 8.099$，因此 $\mQ$ 不定，这个 QP 不是凸的。*（此例的另一种较早写法使用非对称矩阵 $\bigl(\begin{smallmatrix}2&-4\\-2&2\end{smallmatrix}\bigr)$，它给出 $\tfrac12\vx^\T\mQ\vx=x_1^2-3x_1x_2+x_2^2$，并不能还原原目标函数。）*

(sec-06-01-local-global)=
## 局部极小点与全局极小点

:::{prf:definition} 全局极小点与局部极小点
:label: def-06-01-local-global

设 $\Omega\subseteq\R^n$，$f:\Omega\to\R$，$\vx^*\in\Omega$。

1. 若对所有 $\vx\in\Omega$ 有 $f(\vx)\ge f(\vx^*)$，则称 $\vx^*$ 为 $f$ 在 $\Omega$ 上的**全局极小点**（global minimum）。
2. 若存在 $\varepsilon>0$，使得对所有满足 $\norm{\vx-\vx^*}<\varepsilon$ 的 $\vx\in\Omega$，即对所有 $\vx\in\Omega\cap B(\vx^*,\varepsilon)$，都有 $f(\vx)\ge f(\vx^*)$，则称 $\vx^*$ 为**局部极小点**（local minimum）。
3. 若对所有这样的 $\vx\ne\vx^*$ 上述不等式严格成立，则称该局部（全局）极小点是**严格的**。

把不等号反向即得极大点的定义。
:::

每个全局极小点都是局部极小点；反之一般不成立。局部极小点只需优于它的可行邻点，而半径 $\varepsilon$ 可能非常小。注意局部极小点是相对于 $\Omega$ 定义的：$\Omega$ 边界上的点即使在离开 $\Omega$ 的方向上 $f$ 会继续减小，也可以是局部极小点。

一个一元函数就能同时展示所有情形。考虑

$$
f(x)=(x^2)^x=e^{2x\ln\abs{x}},\qquad -5\le x\le1.7,
$$ (eq-06-01-x-power-2x)

并令 $f(0)=1$（即 $x\to0$ 时 $f$ 的极限）。当 $x\ne0$ 时，$f'(x)=2f(x)(\ln\abs{x}+1)$，它只在 $x=\pm1/e$ 处为零。

- 在 $[-5,-1/e)$ 上有 $\ln\abs{x}>-1$，因此 $f'>0$，$f$ 递增；在 $(-1/e,0)$ 和 $(0,1/e)$ 上 $f$ 递减；在 $(1/e,1.7]$ 上 $f$ 递增。
- $x=-1/e\approx-0.368$ 是**局部极大点**，$f=e^{2/e}\approx2.087$。
- $x=1/e\approx0.368$ 是**局部极小点**，$f=e^{-2/e}\approx0.4791$，但不是全局极小点。
- **全局极小点**是左端点 $x=-5$，$f(-5)=25^{-5}=1.024\times10^{-7}$；全局极大点是右端点，$f(1.7)\approx6.075$。

由此得到两点启示。第一，全局极小点不是驻点：$f'(-5)\approx5.3\times10^{-7}>0$，$x=-5$ 之所以是极小点，只是因为约束 $x\ge-5$ 阻止了函数值继续减小。因此约束问题的最优性条件必须允许梯度非零（{doc}`06-02-optimality-conditions`）。第二，只利用局部信息的算法若从 $x=0.5$ 附近出发，将停在局部极小点 $1/e$。{numref}`fig-06-01-x-power-2x`还标出了 $f$ 的凸区间（在 $[-5,-0.808)$ 和 $(0,1.7]$ 上 $f''>0$）和凹区间（在 $(-0.808,0)$ 上 $f''<0$）；正因为 $f$ 在整个区间上不是凸的，才会出现不是全局极小点的局部极小点。*（由凸变凹的转折点在 $x\approx-0.808$ 处，而不是有时图中所画的 $x\approx-1.1$ 处。）*

```{figure} figures/06-01-x-power-2x.svg
:name: fig-06-01-x-power-2x
:width: 560px
:class: dark-light
:alt: 函数 f(x) = (x²)^x 在区间 -5 到 1.7 上的图像，全局极小点在 x = -5，局部极大点在 -1/e，局部极小点在 1/e，全局极大点在 1.7；阴影带表示凸的部分和凹的部分

函数 {eq}`eq-06-01-x-power-2x` 在 $[-5,1.7]$ 上的图像。局部极小点 $1/e$ 不是全局极小点；全局极小点和全局极大点位于可行区间的端点。阴影部分：$f$ 为凹函数的区间 $(-0.808,0)$。
```

对于凸规划，局部极小点与全局极小点的区别消失了。

:::{prf:theorem} 凸规划的局部极小点是全局极小点
:label: thm-06-01-convex-local-global

若 $f$ 在凸集 $\Omega$ 上是凸函数，则 $f$ 在 $\Omega$ 上的每个局部极小点都是全局极小点，并且全局极小点的集合是凸集。若 $f$ 严格凸，则至多有一个极小点。
:::

证明很简短：如果局部极小点 $\vx^*$ 有一个更好的点 $\vy\in\Omega$，那么由 $\Omega$ 的凸性，点 $\vx^*+\theta(\vy-\vx^*)$（$\theta\in(0,1]$）都是可行的；由 $f$ 的凸性，它们都优于 $\vx^*$，其中包括任意接近 $\vx^*$ 的点。完整证明见{prf:ref}`thm-02-03-local-global`，基于梯度的版本见{prf:ref}`thm-06-02-convex-sufficient`。这个定理是值得花力气识别凸性（{doc}`02-02-convex-functions`）的主要原因：对于凸规划，任何能找到局部极小点的方法都能找到全局极小点。

(sec-06-01-existence)=
## 极小点的存在性

即使目标函数有界，问题也可能没有最优解（例子见{prf:ref}`ex-08-08-difficulties`）。下面两个充分条件涵盖了大多数模型。

:::{prf:theorem} 极小点的存在性
:label: thm-06-01-existence

设 $\Omega\subseteq\R^n$ 非空，$f$ 在 $\Omega$ 上连续。

1. （魏尔斯特拉斯）若 $\Omega$ 是有界闭集，则 $f$ 在 $\Omega$ 上取得全局极小值和全局极大值。
2. 若 $\Omega$ 是闭集且 $f$ 是**强制函数**（coercive function），即当 $\vx\in\Omega$ 且 $\norm{\vx}\to\infty$ 时 $f(\vx)\to\infty$，则 $f$ 在 $\Omega$ 上取得全局极小值。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

第 1 部分即魏尔斯特拉斯定理（{prf:ref}`thm-02-01-weierstrass`）。对于第 2 部分，任取 $\vx^0\in\Omega$。由强制性，存在半径 $R$，使得对每个满足 $\norm{\vx}>R$ 的 $\vx\in\Omega$ 都有 $f(\vx)>f(\vx^0)$。集合 $\Omega_R=\{\vx\in\Omega:\norm{\vx}\le R\}$ 是有界闭集且包含 $\vx^0$，因此由第 1 部分，$f$ 在 $\Omega_R$ 上有极小点 $\vx^*$。于是对 $\Omega_R$ 之外的所有可行点 $\vx$ 有 $f(\vx^*)\le f(\vx^0)<f(\vx)$，所以 $\vx^*$ 是 $\Omega$ 上的全局极小点。
:::

由连续函数通过“$\le$”和“$=$”定义的约束给出闭集，因此对于 {eq}`eq-06-01-nlp`，闭性通常是自动满足的；而像 $x_1+x_2<2$ 这样的严格不等式则不然。在{prf:ref}`ex-06-01-three-feasible-sets`中，问题 2 和问题 3 的可行集是紧集，问题 1 的目标函数是强制函数，因此无需任何计算就知道三个问题都有最优解。在{prf:ref}`ex-06-01-consumer-regression`中，可行集是紧集（$0\le\beta_i\le1$，$0\le q_i\le\min_tQ_{ti}$），因此全局极小点存在；困难只在于找到它。

(sec-06-01-difficulties)=
## 非线性规划的困难

{numref}`tab-06-01-lp-vs-nlp`中与 LP 的比较可以归纳为下面一系列会出问题的地方，它们各自影响着后续各章中方法的设计。

1. **局部极小点不一定是全局极小点。** 基于导数的算法找到的是驻点和局部极小点；要证明全局最优性，需要凸性或穷举搜索（{numref}`fig-06-01-x-power-2x`，{prf:ref}`ex-06-01-consumer-regression`）。
2. **最优解不一定在顶点处。** 可行集可能根本没有顶点，而在有限个候选点中枚举（单纯形法背后的思想）行不通（{prf:ref}`ex-06-01-three-feasible-sets`）。
3. **算法是迭代的、数值的。** 算法产生的点列只在极限意义下收敛，必须用容差来终止。尺度不良使这一点变得棘手：在{prf:ref}`ex-06-01-consumer-regression`中，$f$ 对份额的偏导数是 $10^5$ 量级，而对数量的偏导数是 1 的量级，许多次求解器运行都过早停止。舍入误差和病态问题将结合{doc}`06-04-line-search`和{doc}`06-05-unconstrained-algorithms`中的算法讨论。
4. **有界的目标函数可能没有最优解。** 当 $\Omega$ 无界或不闭时，$f$ 的下确界可能只能逼近而无法达到（{prf:ref}`ex-08-08-difficulties`）；{prf:ref}`thm-06-01-existence`给出了排除这种情况的条件。
5. **可行集可能非凸，甚至不连通。** 像 $g_1(x)\ge g_2(x)$ 这样介于两条曲线之间的约束可能把 $\Omega$ 分成互不相连的几块（{prf:ref}`ex-08-08-difficulties`），在 $\Omega$ 内连续移动的搜索无法从一块到达另一块。
6. **最优性条件只是必要条件。** {doc}`06-02-optimality-conditions`和{doc}`06-03-kkt-conditions`中的条件筛选出候选点，但在非凸问题中，候选点可能是鞍点或局部极大点；在退化情形下，除非满足约束规范（constraint qualification），极小点甚至可能不满足一阶条件。

凸的情形消除了困难 1、5 和 6 的大部分，因此接下来的两章对每个结果都先给出一般形式，再给出凸规划的形式。

## 本章回顾

### 小结

- NLP 在 $g_i(\vx)\le0$ 和 $h_j(\vx)=0$ 的约束下最小化 $f(\vx)$（{eq}`eq-06-01-nlp`）；最大化问题和“$\ge$”约束通过改变符号进行转化。
- 带附加条件的最小二乘估计是典型的 NLP（{prf:ref}`ex-06-01-consumer-regression`）；该例中目标函数非凸，其解通过数值方法求得，而拟合结果暴露了模型本身的一个特点（{numref}`fig-06-01-regression-fit`）。
- 与 LP 不同，NLP 的最优解可能位于 $\Omega$ 的内部或弯曲的边界上（{prf:ref}`ex-06-01-three-feasible-sets`，{numref}`tab-06-01-lp-vs-nlp`）。
- 主要类别：无约束问题、线性约束问题、二次规划（{eq}`eq-06-01-qp`，$\mQ$ 对称）、线性规划和凸规划（$f$ 与 $g_i$ 为凸函数，$h_j$ 为仿射函数）。
- {prf:ref}`def-06-01-local-global`区分了局部极小点与全局极小点；对于凸规划，每个局部极小点都是全局极小点（{prf:ref}`thm-06-01-convex-local-global`）。
- 非空紧集上的连续函数、或闭集上的强制函数必有极小点（{prf:ref}`thm-06-01-existence`）。
- NLP 的困难：局部最优解、非顶点最优解、迭代数值计算、达不到的下确界、非凸或不连通的可行集，以及只是必要条件的最优性条件。

### 常见错误

- **期望解位于顶点。** 只检查 $\Omega$ 的角点会漏掉内部或弯曲边界上的最优解，例如{prf:ref}`ex-06-01-three-feasible-sets`中的那些。
- **忘记因子 $\tfrac12$ 或 $\mQ$ 的对称性。** 在 $\vc^\T\vx+\tfrac12\vx^\T\mQ\vx$ 中，$\mQ$ 的对角元是 $x_j^2$ 系数的两倍，而每对非对角元各自等于 $x_jx_k$ 的完整系数。
- **把含非线性等式约束的问题称为凸问题。** 即使 $f$ 和 $g_i$ 都是凸函数，等式约束也必须是仿射的。
- **把求解器返回的第一个局部极小点当作非凸问题“唯一”的最优解。** 应使用多个初始点，或者证明问题是凸的。
- **把上确界当作最大值。** 每个有界函数都有 $\sup f$；但只有在{prf:ref}`thm-06-01-existence`之类的条件下才存在最大值点。


## 参考文献

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 1,
  introduction and examples of NLP models; Section 2.3, Weierstrass'
  theorem; Section 3.4, minima of convex functions.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Sections
  4.1–4.4, optimization problems, convex problems, quadratic programs.)
- Kraft, D. (1988). *A software package for sequential quadratic
  programming* (Report DFVLR-FB 88-28). Deutsche Forschungs- und
  Versuchsanstalt für Luft- und Raumfahrt.
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming*
  (4th ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3> (Chapter 7,
  local and global minima, first-order necessary conditions and examples of
  unconstrained problems.)
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Chapter 1,
  classes of optimization problems; Section 2.1, local and global
  solutions.)
- Sheffi, Y. (1985). *Urban transportation networks: Equilibrium analysis
  with mathematical programming methods*. Prentice-Hall. (Chapter 2, basic
  concepts in minimization problems.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 20, classical optimization theory.)
