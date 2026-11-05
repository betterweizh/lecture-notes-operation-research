# 线性规划的几何

## 概述

本章把线性规划问题（linear program，LP）看作几何对象。只有两个变量时，可行域可以画出来，**图解法**（graphical method）（{prf:ref}`alg-03-02-graphical`）通过在可行域上平移一条目标函数值恒定的直线来求最优解；把它应用于{doc}`03-01-lp-model-building`中的产品组合问题（{prf:ref}`ex-03-02-graphical`）可以看到：最优解位于一个角点上，一些约束是紧的，另一些约束有松弛，并且对分数最优解取整与求解整数问题并不是一回事。LP 的四种可能结果（唯一最优解、多重最优解、无界、不可行）见{numref}`fig-03-02-graphical-cases`。为了处理两个以上的变量，本章转向代数。每个 LP 都可以写成**标准形**（standard form）（{prf:ref}`prop-03-02-standard-form`），其可行集是一个多面体，而该多面体的角点恰好是方程组的**基本可行解**（basic feasible solution，BFS）（{prf:ref}`thm-03-02-extreme-bfs`）。**线性规划基本定理**（fundamental theorem of linear programming）（{prf:ref}`thm-03-02-fundamental`）指出，有最优解的 LP 必有最优基本可行解，因此搜索范围可以限制在有限个点上。最后，用基来表示 LP，得到**约化形式**（reduced form）{eq}`eq-03-02-reduced-form` 和第一个最优性检验（{prf:ref}`thm-03-02-optimality-test`），这是{doc}`03-03-simplex-method`中单纯形法的出发点。本章用到{doc}`01-01-linear-algebra`中的线性无关、秩和逆矩阵，以及{doc}`02-01-convex-sets`中的凸集、多面体和极点。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 用图解法求解两变量 LP，并识别唯一最优解、多重最优解、无界和不可行等情形；
- 计算松弛量和剩余量，并找出紧约束；
- 把任何 LP 化为标准形；
- 解释为什么 LP 的可行集是凸多面体；
- 计算 $\mA\vx=\vb$ 的基本解，判断其是否可行，并识别退化的基本解；
- 证明极点与基本可行解是一致的，并叙述和证明 LP 基本定理；
- 关于一个基把 LP 写成约化形式，并应用检验数最优性检验。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $P$ | LP 的可行集（多面体），标准形下 $P=\{\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$ |
| $s_i$ | 约束 $i$ 的松弛量（$\le$ 行）或剩余量（$\ge$ 行） |
| $\mB$, $\mN$ | 基矩阵（$\mA$ 的 $m$ 个线性无关的列）以及其余列构成的矩阵 |
| $\vx_B$, $\vx_N$; $\vc_B$, $\vc_N$ | 基变量和非基变量；它们的目标函数系数 |
| $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ | $x_j$ 的检验数，其中 $\vpi^\T=\vc_B^\T\mB^{-1}$ |
| $J(\vx)$ | $\vx$ 的支撑集，即指标集 $\{j:x_j>0\}$ |

(sec-03-02-graphical)=
## 图解法

含两个变量 $x_1,x_2$ 的 LP 可以通过作图求解。每个线性不等式 $a_1x_1+a_2x_2\le b$ 在直线 $a_1x_1+a_2x_2=b$ 的一侧成立，因此可行域是所有这些半平面重叠的多边形（可能无界）。具有给定目标函数值 $z=c_1x_1+c_2x_2$ 的点构成一条直线，称为**等值线**（contour，即等利润线或等成本线）；改变 $z$ 会使这条直线平行移动，$z$ 增大时它沿向量 $(c_1,c_2)$ 的方向移动。

:::{prf:algorithm} 两变量 LP 的图解法
:label: alg-03-02-graphical

**输入：** 变量为 $x_1,x_2$ 的 LP。

1. **画出约束。** 对每个约束，画出它以等式成立的直线，并确定它在哪一侧成立（检验一个点，例如原点）。
2. **画出可行域的阴影**，即所有这些半平面与非负象限的交。若它是空集，则停止：LP 不可行。
3. **画一条等值线** $c_1x_1+c_2x_2=z_0$，$z_0$ 取某个方便的值。
4. **平移等值线**：使其沿改进方向（最大化时为 $+(c_1,c_2)$，最小化时为 $-(c_1,c_2)$）平行移动，直到它仍与可行域接触的最远位置。若它可以无限移动，则停止：LP 无界。
5. **读出最优解。** 最后接触到的点就是最优解。通过求解在该处相交的两条直线的方程得到其坐标，并计算 $z$。
:::

:::{prf:example} 产品组合问题的图解
:label: ex-03-02-graphical

求解{prf:ref}`ex-03-01-product-mix`中的产品组合 LP {eq}`eq-03-01-product-mix`：

$$
\begin{aligned}
\max\ & z=500x_1+450x_2 \\
\st\ & 6x_1+5x_2\le60 \\
& 10x_1+20x_2\le150 \\
& x_1\le8 \\
& x_1,x_2\ge0,
\end{aligned}
$$

其中 $x_1,x_2$ 的单位为百箱/周。求出最优解处每个约束的松弛量，并把 LP 最优解与以整百箱计的最好方案进行比较。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**可行域。** 机器约束线 $6x_1+5x_2=60$ 与坐标轴交于 $(10,0)$ 和 $(0,12)$；存储约束线 $10x_1+20x_2=150$ 交于 $(15,0)$ 和 $(0,7.5)$；需求约束线为 $x_1=8$。三个约束在原点处都成立，所以可行域位于每条直线靠原点的一侧。它是一个五边形，顶点为

$$
\begin{gathered}
(0,0),\quad (8,0),\quad (8,\tfrac{12}{5}),\\
(\tfrac{45}{7},\tfrac{30}{7}),\quad (0,\tfrac{15}{2}),
\end{gathered}
$$

如{numref}`fig-03-02-product-mix-region`所示。顶点 $(\tfrac{45}{7},\tfrac{30}{7})$ 是机器约束线与存储约束线的交点：把 $6x_1+5x_2=60$ 乘以 4 再减去 $10x_1+20x_2=150$，得 $14x_1=90$，所以 $x_1=\tfrac{45}{7}$，$x_2=(60-6\cdot\tfrac{45}{7})/5=\tfrac{30}{7}$。直线 $x_1=8$ 与 $6x_1+5x_2=60$ 交于 $(8,\tfrac{12}{5})$。

**等值线。** 等值线 $500x_1+450x_2=z$ 的斜率为 $-\tfrac{10}{9}$。等值线 $z=2250$ 穿过可行域；$z=4000$ 经过顶点 $(8,0)$。继续移动，等值线最后接触到的可行域中的点是 $(\tfrac{45}{7},\tfrac{30}{7})$。等值线的斜率 $-\tfrac{10}{9}$ 介于在该处相交的两条直线的斜率 $-\tfrac65$（机器）和 $-\tfrac12$（存储）之间，这就是为什么最后接触到的是这个角点而不是其他角点。因此

$$
\begin{aligned}
x_1^*&=\tfrac{45}{7}=6\tfrac37,\\
x_2^*&=\tfrac{30}{7}=4\tfrac27,\\
z^*&=500\cdot\tfrac{45}{7}+450\cdot\tfrac{30}{7}\\
&=\tfrac{36000}{7}=5142\tfrac67 .
\end{aligned}
$$

作为验算，其他顶点处的目标函数值分别为 0、4000、5080 和 3375，都更小。

**松弛量。** 在最优解处，机器使用了 $6\cdot\tfrac{45}{7}+5\cdot \tfrac{30}{7}=60$ 小时，仓库存放了 $10\cdot\tfrac{45}{7}+20\cdot \tfrac{30}{7}=150$ 立方米：两个约束都是紧约束，松弛量为零。需求约束的松弛量为 $8-\tfrac{45}{7}=\tfrac{11}{7}=1\tfrac47$ 百箱，即每周约还有 157 箱果汁杯可以售出。

*（本例的一个早期版本把松弛量写成每周 $1\tfrac47$ 箱；由于 $x_1$ 以百箱为单位，松弛量应为每周 $1\tfrac47$ 百箱，即 $157\tfrac17$ 箱。）*

**整百箱。** 若方案必须以整百箱计，把最优解向下取整得到 $(6,4)$，它是可行的，$z=4800$。但它不是最好的整数方案：可行域中满足 $z\ge4800$ 的整数点有 $(8,2)$（$z=4900$）、$(7,3)$（$z=4850$）和 $(6,4)$（$z=4800$），所以整数最优解是 $(8,2)$，在 $x_1$ 方向上离 LP 最优解很远。对 LP 解取整与求解整数规划不是一回事。

```{figure} figures/03-02-product-mix-region.svg
:name: fig-03-02-product-mix-region
:width: 540px
:class: dark-light
:alt: 由三条约束线和坐标轴围成的五边形可行域，以及三条平行的虚线目标函数等值线；最高的一条只在顶点 (45/7, 30/7) 处与可行域接触

产品组合 LP 的可行域。虚线等值线 $500x_1+450x_2=z$（$z=2250$、$4000$ 和 $\tfrac{36000}{7}$）随 $z$ 增大向右上方移动；最后一条只在最优顶点 $\vx^*=(\tfrac{45}{7},\tfrac{30}{7})$ 处与可行域接触。
```
:::

{prf:ref}`ex-03-02-graphical`的解答用到了松弛量和紧约束等术语，它们具有一般意义。

:::{prf:definition} 松弛量、剩余量与紧约束
:label: def-03-02-slack

在可行点 $\vx$ 处，约束 $\mA_{i\cdot}\vx\le b_i$ 的**松弛量**（slack）为 $s_i=b_i-\mA_{i\cdot}\vx\ge0$，即资源的未使用量；约束 $\mA_{i\cdot}\vx\ge b_i$ 的**剩余量**（surplus）为 $s_i=\mA_{i\cdot}\vx-b_i\ge0$，即超出需求的量。在 $\vx$ 处松弛量或剩余量为零的约束称为在 $\vx$ 处的**紧约束**（binding constraint，又称有效约束）；否则称为**非紧约束**（non-binding constraint）。
:::

紧约束代表稀缺资源：放松它可以使目标函数改进（在{doc}`03-05-lp-duality`中，每单位的改进量就是影子价格，{ref}`sec-03-05-shadow-prices`）。非紧约束可以在其松弛量范围内收紧而不改变最优解。

(sec-03-02-graphical-cases)=
### 四种可能的结果

产品组合问题有唯一最优解。另外还可能出现三种结果，{numref}`fig-03-02-graphical-cases`对每一种都给出了一个小规模 LP。

- **唯一最优解。** 对 $\max\ x_1+x_2$，约束为 $x_1+2x_2\le6$、$2x_1+x_2\le6$、$\vx\ge\vzero$，最后一条等值线只在 $(2,2)$ 处与可行域接触，$z^*=4$。
- **多重最优解**（alternative optima）。在同一可行域上取目标函数 $\max\ x_1+2x_2$，等值线平行于边 $x_1+2x_2=6$，最后一条等值线 $z^*=6$ 包含从 $(0,3)$ 到 $(2,2)$ 的整条边。每个点 $\theta(0,3)+(1-\theta)(2,2)$（$0\le\theta\le1$）都是最优的。最优解集仍是可行域的一个面，并且仍然包含角点。
- **无界 LP。** 对 $\max\ x_1+x_2$，约束为 $-x_1+x_2\le2$、$x_1-2x_2\le2$、$\vx\ge\vzero$，可行域无界，目标函数也无界：从 $(0,0)$ 出发，点 $t\,(2,1)$（$t\ge0$）都是可行的（$-2t+t\le0\le2$，$2t-2t=0\le2$），其目标函数值 $3t$ 无限增大。仅仅*可行域*无界并不会使 LP 无界：在同一可行域上 $\min\ x_1+x_2$ 有最优解 $(0,0)$。
- **不可行 LP。** 约束 $x_1+x_2\le2$ 和 $x_1+2x_2\ge6$ 在 $\vx\ge\vzero$ 下没有公共点，因为对非负的 $\vx$ 有 $x_1+2x_2\le2(x_1+x_2)\le4$。

```{figure} figures/03-02-graphical-cases.svg
:name: fig-03-02-graphical-cases
:width: 640px
:class: dark-light
:alt: 四个小图：一个只有一个最优角点的可行域；同一可行域，最优解为一条边；一个无界可行域，带有目标函数沿之增大的箭头；两个不重叠的半平面

LP 的四种结果。(a) 唯一最优解 $(2,2)$。(b) 多重最优解：等值线 $x_1+2x_2=6$ 包含一整条边。(c) 无界：目标函数 $x_1+x_2$ 沿从原点出发、方向为 $(2,1)$ 的射线增大。(d) 不可行：半平面 $x_1+x_2\le2$ 与 $x_1+2x_2\ge6$ 在非负象限中不相交。
```

在每种有最优解的情形中，可行域的某个角点都是最优的。这一观察是单纯形法的基础，本章其余部分将对任意多个变量证明这一点（{prf:ref}`thm-03-02-fundamental`）。证明需要两样东西：角点的代数描述（在 $\R^n$ 中无法靠作图得到），以及书写所有 LP 的统一方式。

(sec-03-02-standard-form)=
## 标准形

:::{prf:definition} 标准形
:label: def-03-02-standard-form

若 LP 写成

$$
\begin{aligned}
\min_{\vx}\ & z=\vc^\T\vx \\
\st\ & \mA\vx=\vb, \\
& \vx\ge\vzero,
\end{aligned}
$$ (eq-03-02-standard-form)

其中 $\mA\in\R^{m\times n}$，$\vb\in\R^m$，$\vc\in\R^n$，即所有约束都是方程且所有变量都非负，则称该 LP 为**标准形**（standard form）。当需要构造初始基时（{doc}`03-03-simplex-method`），还要求 $\vb\ge\vzero$。
:::

下面的理论都针对最小化问题来叙述；约束相同的最大化问题可以用 {eq}`eq-03-01-max-min` 处理，或者在最优性检验中把不等号反向。本章的贯穿示例都是最大化问题。

每个 LP 都可以化为标准形，因此理论和单纯形法只需处理这一种形式。

:::{prf:proposition} 化为标准形
:label: prop-03-02-standard-form

每个 LP 都可以通过以下运算改写为标准形 LP；这些运算既不改变最优方案的集合（在撤销变量代换之后），也不改变最优值（第 1 步中的符号除外）。

1. **目标函数方向。** 把 $\max\ \vc^\T\vx$ 替换为 $\min\ (-\vc)^\T\vx$。
2. **$\le$ 约束。** 把 $\mA_{i\cdot}\vx\le b_i$ 替换为 $\mA_{i\cdot}\vx+s_i=b_i$，其中 $s_i\ge0$ 是新的**松弛变量**（slack variable）。
3. **$\ge$ 约束。** 把 $\mA_{i\cdot}\vx\ge b_i$ 替换为 $\mA_{i\cdot}\vx-s_i=b_i$，其中 $s_i\ge0$ 是新的**剩余变量**（surplus variable）。
4. **非正变量。** 若 $x_j\le0$，代换 $x_j=-x_j'$，$x_j'\ge0$。
5. **自由变量**（free variable）。若 $x_j$ 没有符号限制，代换 $x_j=x_j^+-x_j^-$，$x_j^+,x_j^-\ge0$。
6. **界。** 若 $x_j\ge l_j$，代换 $x_j=l_j+x_j'$，$x_j'\ge0$；上界 $x_j\le u_j$ 是普通的 $\le$ 约束（第 2 步）。
7. **负的右端项。** 把 $b_i<0$ 的方程乘以 $-1$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

第 1 步就是 {eq}`eq-03-01-max-min`。在第 2 步中，映射 $\vx\mapsto(\vx,\ s_i=b_i-\mA_{i\cdot}\vx)$ 是满足该不等式的点与满足该方程且 $s_i\ge0$ 的点 $(\vx,s_i)$ 之间的一一对应；目标函数不含 $s_i$，因此在对应点处取值相同。第 3 步与此相同，只是 $s_i=\mA_{i\cdot}\vx-b_i$。第 4 步和第 6 步是变量代换，分别把 $\{x_j\le0\}$ 和 $\{x_j\ge l_j\}$ 一一映射到 $\{x_j'\ge0\}$ 上，并且在第 6 步中目标函数只改变常数 $c_jl_j$。在第 5 步中，每个实数 $x_j$ 都可以写成 $x_j^+-x_j^-$，$x_j^\pm\ge0$（例如 $x_j^+=\max\{x_j,0\}$，$x_j^-=\max\{-x_j,0\}$），而每个这样的数对都给出一个实数 $x_j$；该映射是满射但不是单射，并且目标函数只依赖于差，所以最优值相同，每个最优的 $(x_j^+,x_j^-)$ 都给出一个最优的 $x_j$。第 7 步不改变方程的解集。
:::

每种转换都是一种建模手段，而不是近似。请注意两个推论。松弛变量或剩余变量在每个点处都具有{prf:ref}`def-03-02-slack`中的含义，所以它直接指出哪些约束是紧的。另外，一个自由变量会变成互为相反数的两列；这样的两列线性相关，因此（由下面的结果）在基本可行解处 $x_j^+$、$x_j^-$ 中至多有一个为正，于是拆分恰好还原出 $x_j^+=\max\{x_j,0\}$，$x_j^-=\max\{-x_j,0\}$。

对产品组合问题，把松弛变量记为 $x_3,x_4,x_5$（未使用的机器小时数、未使用的存储空间立方米数、以百箱计的未满足需求），其标准形为

$$
\begin{aligned}
\max\ & z=500x_1+450x_2 \\
\st\ & 6x_1+5x_2+x_3=60 \\
& 10x_1+20x_2+x_4=150 \\
& x_1+x_5=8 \\
& x_1,\dots,x_5\ge0,
\end{aligned}
$$ (eq-03-02-product-mix-standard)

于是

$$
\begin{gathered}
\mA=\begin{pmatrix}6&5&1&0&0\\10&20&0&1&0\\1&0&0&0&1\end{pmatrix},\\
\vb=\begin{pmatrix}60\\150\\8\end{pmatrix},\\
\vc=(500,\,450,\,0,\,0,\,0)^\T .
\end{gathered}
$$

（这里保留了最大化形式；在{doc}`03-03-simplex-method`的单纯形法中，最大化的规则与最小化的规则并列给出。）方案 $x_1=\tfrac{45}{7}$、$x_2=\tfrac{30}{7}$ 对应于 $\vx=(\tfrac{45}{7},\tfrac{30}{7},0,0,\tfrac{11}{7})$：为零的松弛变量 $x_3=x_4=0$ 标示出紧约束。

(sec-03-02-polyhedra)=
## 半空间与多面体

对非零的 $\va\in\R^n$ 和 $b\in\R$，集合 $\{\vx:\va^\T\vx=b\}$ 是一个**超平面**（hyperplane），集合 $\{\vx:\va^\T\vx\ge b\}$ 和 $\{\vx:\va^\T\vx\le b\}$ 是由它分出的两个闭**半空间**（half-space）。在平面上，超平面是一条直线，半空间是位于其一侧的闭半平面（{numref}`fig-03-02-half-space`）。LP 的每个线性不等式约束都描述一个半空间，每个非负条件 $x_j\ge0$ 都是 $\va=\ve_j$、$b=0$ 的半空间。

```{figure} figures/03-02-half-space.svg
:name: fig-03-02-half-space
:width: 420px
:class: dark-light
:alt: 经过 (0,1) 和 (1,0) 的直线 x1 + x2 = 1，其右上方区域涂有阴影，并标注 x1 + x2 至少为 1

闭半空间 $x_1+x_2\ge1$（阴影部分）及其边界，即经过 $(0,1)$ 和 $(1,0)$ 的超平面 $x_1+x_2=1$。法向量 $(1,1)$ 指向半空间内部。
```

等式约束又如何呢？方程 $\va^\T\vx=b$ 成立当且仅当 $\va^\T\vx\le b$ 和 $\va^\T\vx\ge b$ 都成立，所以它是两个半空间的交。因此，无论 $\le$、$\ge$ 和 $=$ 行如何搭配，每个 LP 的可行集都是有限个闭半空间的交，即一个**多面体**（polyhedron）（{prf:ref}`def-02-01-polyhedron`）；有界的多面体称为**多胞形**（polytope）。

:::{prf:proposition} LP 的可行集是凸多面体
:label: prop-03-02-convex-polyhedron

LP 的可行集是闭的凸多面体。特别地，$P=\{\vx\in\R^n:\mA\vx=\vb,\ \vx\ge\vzero\}$ 是凸集：若 $\vy,\vw\in P$ 且 $0\le\theta\le1$，则 $\theta\vy+(1-\theta)\vw\in P$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

如上所述，它是有限个闭半空间的交。也可以直接验证：$\mA(\theta\vy+(1-\theta)\vw)=\theta\vb+(1-\theta)\vb=\vb$，并且由于两项都非负，$\theta\vy+(1-\theta)\vw\ge\vzero$。闭性成立，因为每个半空间都是一个闭区间在某个连续线性函数下的原像。
:::

凸性使 LP 易于处理。凸集中不是最优的点，总可以沿着一条保持可行的线段向更好的点移动而得到改进，而线性目标函数沿线段线性变化；因此 LP 的*局部*最优解就是全局最优解（{prf:ref}`thm-02-03-local-global`）。可行域的角点由极点的概念来刻画，极点已在{prf:ref}`def-02-01-extreme-point`中定义，这里重述如下。

:::{prf:definition} 极点
:label: def-03-02-extreme-point

若凸集 $C$ 中的点 $\vx$ 不能写成 $\vx=\theta\vy+(1-\theta)\vw$，其中 $\vy,\vw\in C$，$\vy\ne\vw$，$0<\theta<1$，则称 $\vx$ 为 $C$ 的**极点**（extreme point，又称角点或顶点）；也就是说，$\vx$ 不严格位于 $C$ 中任何线段的内部。
:::

产品组合可行域的极点就是{prf:ref}`ex-03-02-graphical`的解答中列出的五个顶点。边上的点（例如 $(8,1)=\tfrac12(8,0)+\tfrac12(8,2)$）或内点都不是极点。要在没有图的情况下找出 $\R^n$ 中的极点，必须用代数方法来描述它们。

(sec-03-02-basic-solutions)=
## 基本解

从现在起，LP 都是标准形 {eq}`eq-03-02-standard-form`，其可行集为 $P=\{\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$。单纯形法各章始终采用以下假设。

:::{prf:assumption} 行满秩
:label: asm-03-02-full-rank

$\mA\in\R^{m\times n}$ 的秩为 $m$，因此它的 $m$ 行线性无关，且 $m\le n$。
:::

这个假设不失一般性。若 $\mA$ 的某一行是其他行的线性组合，则要么方程组 $\mA\vx=\vb$ 不相容（LP 不可行），要么 $\vb$ 也满足同样的组合关系，此时该行是可以删去的**冗余**（redundant）方程。（单纯形法的第一阶段会检测出这样的行，{prf:ref}`thm-03-04-phase-one`。）

方程 $\mA\vx=\vb$ 表明 $\vb$ 是 $\mA$ 的各列以 $x_j$ 为系数的组合：

$$
x_1\mA_{\cdot1}+x_2\mA_{\cdot2}+\dots+x_n\mA_{\cdot n}=\vb .
$$ (eq-03-02-column-view)

一个可行点就是把 $\vb$ 写成各列的*非负*组合的一种方式。最简单的这类表示使用尽可能少的列。

:::{prf:definition} 基、基本解、基本可行解
:label: def-03-02-basic-solution

在{prf:ref}`asm-03-02-full-rank`下：

1. $\mA$ 的 $m$ 个线性无关的列构成的集合称为一个**基**（basis）。由这些列构成的 $m\times m$ 矩阵 $\mB$（**基矩阵**，basis matrix）是非奇异的。其余各列构成矩阵 $\mN$。
2. $\mB$ 中各列对应的变量称为**基变量**（basic variable）$\vx_B$；其余变量称为**非基变量**（nonbasic variable）$\vx_N$。
3. 该基的**基本解**（basic solution）是 $\mA\vx=\vb$ 满足 $\vx_N=\vzero$ 的唯一解，即
   $$
   \begin{aligned}
   \vx_B&=\mB^{-1}\vb,\\
   \vx_N&=\vzero .
   \end{aligned}
   $$
4. 满足 $\vx_B\ge\vzero$ 的基本解称为**基本可行解**（basic feasible solution，BFS）。
:::

计算基本解分三步：选取 $m$ 列，检验它们线性无关（$\det\mB\ne0$），求解 $\mB\vx_B=\vb$。选取列的方式至多有 $\binom{n}{m}$ 种，因此基本解至多有 $\binom{n}{m}$ 个。

:::{prf:example} 枚举一个小方程组的基本解
:label: ex-03-02-basic-solutions

求方程组
$$
\begin{gathered}
2x_1+3x_2+x_3=5,\\
x_1+2x_2-x_3=3,\\
x_1,x_2,x_3\ge0,
\end{gathered}
$$

的所有基本解，并判断哪些是可行的。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

这里 $m=2$，$n=3$，$\mA=\begin{pmatrix}2&3&1\\1&2&-1\end{pmatrix}$，$\vb=\begin{pmatrix}5\\3\end{pmatrix}$。两行线性无关（哪一行都不是另一行的倍数），所以{prf:ref}`asm-03-02-full-rank`成立。候选基共有 $\binom32=3$ 个。

- $\{x_1,x_2\}$：$\mB=\begin{pmatrix}2&3\\1&2\end{pmatrix}$，$\det\mB=1$，$\mB^{-1}=\begin{pmatrix}2&-3\\-1&2\end{pmatrix}$，$\vx_B=\mB^{-1}\vb=(1,1)^\T$。基本解 $(1,1,0)$ 是可行的。
- $\{x_1,x_3\}$：$\mB=\begin{pmatrix}2&1\\1&-1\end{pmatrix}$，$\det\mB=-3$，$\mB^{-1}=\begin{pmatrix}\tfrac13&\tfrac13\\ \tfrac13&-\tfrac23\end{pmatrix}$，$\vx_B=(\tfrac83,-\tfrac13)^\T$。基本解 $(\tfrac83,0,-\tfrac13)$ 不可行（$x_3<0$）。
- $\{x_2,x_3\}$：$\mB=\begin{pmatrix}3&1\\2&-1\end{pmatrix}$，$\det\mB=-5$，$\mB^{-1}=\begin{pmatrix}\tfrac15&\tfrac15\\ \tfrac25&-\tfrac35\end{pmatrix}$，$\vx_B=(\tfrac85,\tfrac15)^\T$。基本解 $(0,\tfrac85,\tfrac15)$ 是可行的。

每个解都满足两个方程，例如 $2\cdot0+3\cdot\tfrac85+\tfrac15=5$ 和 $0+2\cdot\tfrac85-\tfrac15=3$。该方程组有两个基本可行解 $(1,1,0)$ 和 $(0,\tfrac85,\tfrac15)$。可行集是连接它们的线段；{prf:ref}`ex-08-01-feasible-segment`画出了这条线段。
:::

对标准形下的产品组合问题 {eq}`eq-03-02-product-mix-standard`（$m=3$，$n=5$），共有 $\binom53=10$ 个候选基。{numref}`tab-03-02-product-mix-basic`列出了这些基，{numref}`fig-03-02-product-mix-basic-solutions`在 $(x_1,x_2)$ 平面上标出了每个基本解：每个基本解都是五条直线 $6x_1+5x_2=60$、$10x_1+20x_2=150$、$x_1=8$、$x_1=0$、$x_2=0$ 中某两条的交点，因为两个非基变量被置为零，而每个变量为零就意味着“位于相应的直线上”。

```{list-table} 标准形下产品组合问题的十个候选基。非基变量为零；破折号表示奇异的基。
:name: tab-03-02-product-mix-basic
:header-rows: 1

* - 基变量
  - $(x_1,x_2,x_3,x_4,x_5)$
  - 是否可行
  - $z$
* - $x_1,x_2,x_3$
  - $(8,\tfrac72,-\tfrac{11}{2},0,0)$
  - 否
  - 5575
* - $x_1,x_2,x_4$
  - $(8,\tfrac{12}{5},0,22,0)$
  - 是
  - 5080
* - $x_1,x_2,x_5$
  - $(\tfrac{45}{7},\tfrac{30}{7},0,0,\tfrac{11}{7})$
  - 是
  - $\tfrac{36000}{7}$
* - $x_1,x_3,x_4$
  - $(8,0,12,70,0)$
  - 是
  - 4000
* - $x_1,x_3,x_5$
  - $(15,0,-30,0,-7)$
  - 否
  - 7500
* - $x_1,x_4,x_5$
  - $(10,0,0,50,-2)$
  - 否
  - 5000
* - $x_2,x_3,x_4$
  - —
  - 奇异
  -
* - $x_2,x_3,x_5$
  - $(0,\tfrac{15}{2},\tfrac{45}{2},0,8)$
  - 是
  - 3375
* - $x_2,x_4,x_5$
  - $(0,12,0,-90,8)$
  - 否
  - 5400
* - $x_3,x_4,x_5$
  - $(0,0,60,150,8)$
  - 是
  - 0
```

```{figure} figures/03-02-product-mix-basic-solutions.svg
:name: fig-03-02-product-mix-basic-solutions
:width: 540px
:class: dark-light
:alt: 产品组合可行域，其角点处有五个实心点，即基本可行解；可行域外有四个空心点，位于 (8, 3.5)、(10, 0)、(15, 0) 和 (0, 12)，即不可行的基本解

在 $(x_1,x_2)$ 平面上画出的产品组合问题的基本解。实心点是五个基本可行解，它们恰好是可行域的角点；空心点是四个不可行的基本解，即违反某个其他约束的约束线交点，例如 $(8,\tfrac72)$，它需要 65.5 机器小时。
```

表格证实了图中所显示的结论。五个基本可行解就是可行域的五个角点，其中最大的 $z$ 值 $\tfrac{36000}{7}$ 就是用图解法求得的最优值。不可行的基本解可能有更大的 $z$（在 $(15,0)$ 处为 7500），这就是必须检验可行性的原因。基 $\{x_2,x_3,x_4\}$ 是奇异的：它的各列在第三行都是零，而事实上直线 $x_1=0$ 与 $x_1=8$ 平行，永不相交。

基本解中的基变量可以为零。这种看似次要的情形对单纯形法很重要。

:::{prf:definition} 退化基本解
:label: def-03-02-degenerate

若至少有一个基变量为零，即正分量少于 $m$ 个，则称基本解是**退化**的（degenerate）。否则称它是**非退化**的（nondegenerate）。
:::

退化的 BFS 属于不止一个基：可以把为零的基变量与一个合适的非基变量交换而不改变该点。从几何上看，在 $\R^2$ 中，围成可行域的直线中有两条以上经过该角点。在{doc}`03-03-simplex-method`中，退化会导致改变基但不改变点的主元变换，在极少数情况下还会导致循环（cycling）（{ref}`sec-03-03-degeneracy`）。{prf:ref}`ex-08-01-bfs-test`中含有一个退化的 BFS。

(sec-03-02-independence-test)=
### 用主元变换检验线性无关性

要判断一组列是否构成基，或者判断某个给定解所用的列是否线性无关，就需要检验线性无关性。除了计算行列式（只适用于方阵）或秩之外，一个简单的主元变换检验法还能在线性相关时给出相关关系。

:::{prf:algorithm} 用主元变换检验线性无关性
:label: alg-03-02-independence

**输入：** 向量 $\vv_1,\dots,\vv_k\in\R^m$。

1. **建立**表 $[\,\mI_k\mid\mV\,]$，其中 $\mV$ 的第 $r$ 行为 $\vv_r^\T$（左边部分是 $k\times k$ 单位矩阵）。
2. **依次考察**第 $r=1,\dots,k$ 行。
   - 若第 $r$ 行的右边部分为零，则**停止**：这些向量线性相关，第 $r$ 行的左边部分 $(\gamma_1,\dots,\gamma_k)$ 给出相关关系 $\sum_i\gamma_i\vv_i=\vzero$。
   - 否则，在第 $r$ 行右边部分中选一个非零元素作为主元，用它除第 $r$ 行，再从其他各行中减去第 $r$ 行的适当倍数，使主元列的其余元素都变为零。
3. **结论。** 若所有 $k$ 行都完成了主元变换，则这些向量线性无关。
:::

为什么右边部分为零就意味着线性相关？表中的每一行都是初始各行 $[\ve_r^\T\mid\vv_r^\T]$ 的组合。等于 $\sum_i\gamma_i[\ve_i^\T\mid\vv_i^\T]$ 的行，其左边部分为 $(\gamma_1,\dots,\gamma_k)$，右边部分为 $\sum_i\gamma_i\vv_i^\T$，并且行变换保持这种结构。所以右边部分为零的行表明 $\sum_i\gamma_i\vv_i=\vzero$，而其左边部分不为零，因为左边的分块在行变换下始终非奇异。因此这些向量线性相关，且相关系数恰好就是这些数。反之，若 $k$ 次主元变换都成功，则右边部分已被化为 $k$ 个线性无关的行，而由于左边分块非奇异，原来的向量张成同一个 $k$ 维空间。{prf:ref}`ex-08-01-independence-pivoting`把这一检验应用于 $\R^5$ 中的四个向量，得到 $-\vv_1-\vv_2-\vv_3+\vv_4=\vzero$。

(sec-03-02-extreme-bfs)=
## 极点就是基本可行解

产品组合问题的表和图表明，它的角点与基本可行解是一致的。这对每个 LP 都成立。

:::{prf:theorem} 极点与基本可行解
:label: thm-03-02-extreme-bfs

设 $P=\{\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$，$\mA$ 行满秩。点 $\vx\in P$ 是 $P$ 的极点，当且仅当它是基本可行解。此外，$\vx\in P$ 是 BFS，当且仅当其正分量对应的列 $\{\mA_{\cdot j}:j\in J(\vx)\}$ 线性无关。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

*BFS $\Rightarrow$ 极点。* 设 $\vx$ 是基 $\mB$ 的 BFS，并设 $\vx=\theta\vy+(1-\theta)\vw$，其中 $\vy,\vw\in P$，$0<\theta<1$。对每个非基指标 $j$，由 $0=x_j=\theta y_j+(1-\theta)w_j$ 且 $y_j,w_j\ge0$ 得 $y_j=w_j=0$。于是 $\mB\vy_B=\mA\vy=\vb$，所以 $\vy_B=\mB^{-1}\vb=\vx_B$；同理 $\vw_B=\vx_B$。因此 $\vy=\vw=\vx$，$\vx$ 是极点。

*正分量的列线性无关 $\Rightarrow$ BFS。* 令 $J=J(\vx)$，并设列 $\mA_{\cdot j}$（$j\in J$）线性无关。则 $\abs{J}\le m$，并且由于 $\mA$ 的秩为 $m$，可以再从 $\mA$ 中补充 $m-\abs{J}$ 列，把这些列扩充为一个基 $\mB$。点 $\vx$ 满足 $\mA\vx=\vb$，并且在 $J$ 之外、从而在基之外的分量都为零；由{prf:ref}`def-03-02-basic-solution`第 3 条中的唯一性，它就是 $\mB$ 的基本解，而且是可行的。

*极点 $\Rightarrow$ 正分量的列线性无关。* 设列 $\mA_{\cdot j}$（$j\in J$）线性相关：存在 $\vd\ne\vzero$，当 $j\notin J$ 时 $d_j=0$，且 $\mA\vd=\sum_{j\in J}d_j\mA_{\cdot j}=\vzero$。由于对 $j\in J$ 有 $x_j>0$，当 $\varepsilon>0$ 充分小时，点 $\vx+\varepsilon\vd$ 和 $\vx-\varepsilon\vd$ 都非负，并且都满足 $\mA(\vx\pm\varepsilon\vd)=\vb$。所以两者都属于 $P$，它们互不相同，而 $\vx=\tfrac12(\vx+\varepsilon\vd)+\tfrac12(\vx-\varepsilon\vd)$ 不是极点。

这三个蕴含关系证明了两个结论：BFS 的正分量对应的列线性无关（它们都在 $\mB$ 的列之中），而正分量对应的列线性无关的点是 BFS。
:::

这个定理把一个几何性质转化为一个可以检验的代数性质：可行点是角点，当且仅当它所用的列线性无关。由此立即得到两个推论。

:::{prf:corollary} 极点个数有限
:label: cor-03-02-finite

多面体 $P$ 至多有 $\binom{n}{m}$ 个极点。非退化的 BFS 恰好对应一个基；退化的 BFS 可能对应多个基。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

每个极点都是至多 $\binom nm$ 个基中至少一个基的基本解，而一个基只有一个基本解。非退化的 BFS 有 $m$ 个正分量，它们对应的列必定构成其基；对退化的 BFS，正分量对应的列可以用不止一种方式扩充为基。
:::

(sec-03-02-fundamental-theorem)=
## 线性规划基本定理

从图上看，有最优解的 LP 必有一个最优角点。其代数表述就是下面的定理，它的证明沿着保持约束成立的方向移动，把任意可行点或最优点化为基本的点。

:::{prf:theorem} 线性规划基本定理
:label: thm-03-02-fundamental

考虑 LP {eq}`eq-03-02-standard-form`，其中 $\mA$ 行满秩。

1. 若该 LP 有可行解，则它有基本可行解。
2. 若该 LP 有最优解，则它有最优基本可行解。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

我们证明第 2 部分；第 1 部分是去掉目标函数的同一论证（取 $\vc=\vzero$，此时每个可行点都是最优的）。

在所有最优解中选取正分量个数最少的一个 $\vx$，令 $J=J(\vx)$。若列 $\mA_{\cdot j}$（$j\in J$）线性无关，则由{prf:ref}`thm-03-02-extreme-bfs`知 $\vx$ 是 BFS，证毕。设它们线性相关：存在 $\vd\ne\vzero$，当 $j\notin J$ 时 $d_j=0$，且 $\mA\vd=\vzero$。必要时用 $-\vd$ 代替 $\vd$，可以假定 $\vd$ 至少有一个负分量。

与{prf:ref}`thm-03-02-extreme-bfs`的证明一样，对每个 $\abs{\varepsilon}$ 充分小的 $\varepsilon$，$\vx+\varepsilon\vd$ 都是可行的。其目标函数值为 $\vc^\T\vx+\varepsilon\,\vc^\T\vd$。若 $\vc^\T\vd\ne0$，取一个与 $\vc^\T\vd$ 符号相反的小 $\varepsilon$，就会得到一个目标函数值比最优解 $\vx$ 更小的可行点，这不可能。因此 $\vc^\T\vd=0$，只要 $\vx+\varepsilon\vd$ 可行，它就是最优的。

现在把 $\varepsilon$ 从 0 开始增大，直到非负性所允许的最大值：

$$
\varepsilon^*=\min\left\{\frac{x_j}{-d_j}:\ d_j<0\right\}>0 .
$$

点 $\vx'=\vx+\varepsilon^*\vd$ 满足 $\mA\vx'=\vb$，由 $\varepsilon^*$ 的取法它是非负的，并且由于 $\vc^\T\vx'=\vc^\T\vx$，它是最优的。它在 $J$ 之外的分量仍为零，而达到最小值的分量 $j$ 已变为零。所以 $\vx'$ 是正分量比 $\vx$ 少的最优解，与 $\vx$ 的取法矛盾。因此 $J$ 对应的列线性无关，$\vx$ 是最优 BFS。
:::

结合{prf:ref}`thm-03-02-extreme-bfs`，这个定理表明：*若 LP 有最优解，则其可行集的某个极点是最优的。* 它并没有说每个最优解都是极点：存在多重最优解时（{numref}`fig-03-02-graphical-cases`(b)），最优边内部的点也是最优的，但这条边的端点是最优角点。

这个定理还说明了最优解何时存在。由第 1 部分，非空的标准形多面体有 BFS；而{doc}`03-03-simplex-method`将证明，从一个 BFS 出发，采用防循环规则的单纯形法要么找到最优 BFS，要么找到一条使目标函数无限改进的射线（{prf:ref}`thm-03-03-bland`）。所以每个 LP 恰好以以下三种方式之一结束：**不可行**、**无界**，或者有**最优 BFS**。特别地，与一般最优化问题不同，可行且目标函数在可行集上有界的 LP 总能达到其最优值。

**为什么不直接枚举？** 这个定理给出了一个有限算法：计算所有基本解，舍去不可行的，保留最好的。对产品组合问题，这意味着求解十个小方程组（{numref}`tab-03-02-product-mix-basic`），{doc}`08-01-examples-lp-modeling`就用这种方法求解了几个 LP（{prf:ref}`ex-08-01-basic-solutions-2d`、{prf:ref}`ex-08-01-enumerate-bases`）。但候选基的个数会爆炸式增长：有 $m=25$ 个方程、$n=50$ 个变量的 LP 有 $\binom{50}{25}\approx1.26\times10^{14}$ 个候选基。单纯形法则从一个 BFS 移到一个相邻的、更好的 BFS，它需要一种检验，无需与其他所有 BFS 比较就能判断一个 BFS 是否最优。约化形式提供了这种检验。

(sec-03-02-reduced-form)=
## 基与约化形式

固定一个基。重新排列各列后，记 $\mA=[\,\mB\ \mN\,]$，$\vx=(\vx_B,\vx_N)$，$\vc=(\vc_B,\vc_N)$。约束和目标函数变为

$$
\begin{aligned}
\mB\vx_B+\mN\vx_N&=\vb,\\
z&=\vc_B^\T\vx_B+\vc_N^\T\vx_N .
\end{aligned}
$$

由于 $\mB$ 可逆，可以从约束中解出基变量 $\vx_B=\mB^{-1}\vb-\mB^{-1}\mN\vx_N$，代入 $z$ 得

$$
z=\vc_B^\T\mB^{-1}\vb+\bigl(\vc_N^\T-\vc_B^\T\mB^{-1}\mN\bigr)\vx_N .
$$

因此，$\mA\vx=\vb$ 的每个解都由其非基部分 $\vx_N$ 确定，LP 等价于只含非基变量的**约化形式**（reduced form）：

$$
\begin{aligned}
\min_{\vx_N}\ & z=\vc_B^\T\mB^{-1}\vb+\bar\vc_N^\T\vx_N \\
\st\ & \mB^{-1}\vb-\mB^{-1}\mN\vx_N\ge\vzero, \\
& \vx_N\ge\vzero,
\end{aligned}
$$ (eq-03-02-reduced-form)

其中第一个约束表示 $\vx_B\ge\vzero$，且

$$
\begin{gathered}
\bar\vc_N^\T=\vc_N^\T-\vc_B^\T\mB^{-1}\mN,\\
\text{即}\quad \bar c_j=c_j-\vpi^\T\mA_{\cdot j},\\
\text{其中}\quad\vpi^\T=\vc_B^\T\mB^{-1}.
\end{gathered}
$$ (eq-03-02-reduced-costs)

数 $\bar c_j$ 称为 $x_j$ 的**检验数**（reduced cost，又称相对费用）：当非基变量 $x_j$ 从零开始增加、基变量随之调整以保持 $\mA\vx=\vb$ 时，$z$ 的变化率。对基变量，同一公式给出 $\bar c_j=c_j-\vc_B^\T\mB^{-1}\mB\ve_j=0$。向量 $\vpi$ 将在{doc}`03-05-lp-duality`中作为单纯形乘子或影子价格向量再次出现（{prf:ref}`thm-03-05-shadow-price`）。在基本解处（$\vx_N=\vzero$），目标函数值为 $z=\vc_B^\T\mB^{-1}\vb=\vc_B^\T\vx_B$。

:::{prf:theorem} 基本可行解的最优性检验
:label: thm-03-02-optimality-test

设 $\vx$ 是基 $\mB$ 的 BFS，其检验数由 {eq}`eq-03-02-reduced-costs` 给出。

1. 对**最小化**问题，若对所有非基指标 $j$ 有 $\bar c_j\ge0$，则 $\vx$ 是最优的。
2. 对**最大化**问题，若对所有非基指标 $j$ 有 $\bar c_j\le0$，则 $\vx$ 是最优的。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

设 $\vy$ 是任一可行解。它满足 $\mA\vy=\vb$，所以由约化形式的推导，其目标函数值为 $\vc^\T\vy=\vc_B^\T\mB^{-1}\vb+\bar\vc_N^\T\vy_N$；而由于 $\vx_N=\vzero$，$\vc^\T\vx=\vc_B^\T\mB^{-1}\vb$。因此

$$
\vc^\T\vy-\vc^\T\vx=\sum_{j\ \text{非基}}\bar c_j\,y_j .
$$

每个 $y_j\ge0$。在情形 1 中每个 $\bar c_j\ge0$，所以和非负，$\vc^\T\vy\ge\vc^\T\vx$；在情形 2 中每个 $\bar c_j\le0$，所以 $\vc^\T\vy\le\vc^\T\vx$。由于 $\vy$ 是任意的，$\vx$ 是最优的。
:::

这个条件是**充分**的。对非退化的 BFS，它也是必要的：在最小化问题中，若某个 $\bar c_j<0$ 且所有基变量都为正，则稍微增大 $x_j$ 仍能保持 $\vx_B>\vzero$ 并使 $z$ 减小，所以 $\vx$ 不是最优的。对退化的 BFS，它不是必要的：在 $\min\{-x_1:\ x_1+x_2=0,\ x_1,x_2\ge0\}$ 中，唯一的可行点 $(0,0)$ 是最优的，但对基 $\{x_2\}$，$x_1$ 的检验数为 $\bar c_1=-1<0$。同一退化点的另一个基能通过检验（这里是 $\{x_1\}$，$\bar c_2=1$），而单纯形法会找到它。

:::{prf:example} 约化形式与基的最优性
:label: ex-03-02-reduced-form

考虑

$$
\begin{aligned}
\max\ & z=3x_1+x_2+3x_3 \\
\st\ & 2x_1+x_2+x_3+x_4=5 \\
& x_1+2x_2+3x_3+x_5=5 \\
& 2x_1+2x_2+x_3+x_6=6 \\
& x_1,\dots,x_6\ge0,
\end{aligned}
$$

其中 $x_4,x_5,x_6$ 是松弛变量。取由列 $\mA_{\cdot1},\mA_{\cdot3},\mA_{\cdot6}$ 构成的基。把 LP 写成约化形式，求出 BFS，并判断它是否最优。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**分块。** 取 $\vx_B=(x_1,x_3,x_6)$，$\vx_N=(x_2,x_4,x_5)$，

$$
\begin{gathered}
\mB=\begin{pmatrix}2&1&0\\1&3&0\\2&1&1\end{pmatrix},\\
\mN=\begin{pmatrix}1&1&0\\2&0&1\\2&0&0\end{pmatrix},\\
\vc_B=\begin{pmatrix}3\\3\\0\end{pmatrix},\\
\vc_N=\begin{pmatrix}1\\0\\0\end{pmatrix}.
\end{gathered}
$$

$\det\mB=5\ne0$，所以这些列构成一个基，且

$$
\begin{gathered}
\mB^{-1}=\begin{pmatrix}\tfrac35&-\tfrac15&0\\-\tfrac15&\tfrac25&0\\
-1&0&1\end{pmatrix},\\
\mB^{-1}\vb=\begin{pmatrix}2\\1\\1\end{pmatrix},\\
\mB^{-1}\mN=\begin{pmatrix}\tfrac15&\tfrac35&-\tfrac15\\
\tfrac35&-\tfrac15&\tfrac25\\1&-1&0\end{pmatrix}.
\end{gathered}
$$

**BFS。** $\vx_B=\mB^{-1}\vb=(2,1,1)\ge\vzero$，所以基本解 $\vx=(2,0,1,0,0,1)$ 是可行的，$z=3\cdot2+3\cdot1+0\cdot1=9$。它满足三个方程：$4+1+0=5$，$2+3+0=5$，$4+1+1=6$。

**检验数。** $\vpi^\T=\vc_B^\T\mB^{-1}=(\tfrac65,\tfrac35,0)$，所以

$$
\begin{gathered}
\bar c_2=1-\vpi^\T\mA_{\cdot2}=1-\bigl(\tfrac65+\tfrac65+0\bigr)=-\tfrac75,\\
\bar c_4=0-\tfrac65=-\tfrac65,\\
\bar c_5=0-\tfrac35=-\tfrac35 .
\end{gathered}
$$

**约化形式。**

$$
\begin{aligned}
\max\ & z=9-\tfrac75x_2-\tfrac65x_4-\tfrac35x_5 \\
\st\ & x_1=2-\tfrac15x_2-\tfrac35x_4+\tfrac15x_5\ge0 \\
& x_3=1-\tfrac35x_2+\tfrac15x_4-\tfrac25x_5\ge0 \\
& x_6=1-x_2+x_4\ge0 \\
& x_2,x_4,x_5\ge0 .
\end{aligned}
$$

**最优性。** 三个检验数都为负，所以由{prf:ref}`thm-03-02-optimality-test`（最大化）知该 BFS 是最优的：任何可行点都有 $x_2,x_4,x_5\ge0$，而其中任何一个每增加一个单位都会使 $z$ 降到 9 以下。最优解为 $x_1^*=2$、$x_3^*=1$，前两个约束是紧约束，第三个约束的松弛量为 $x_6^*=1$，$z^*=9$，与求解器的结果一致。在该方程组的 $\binom63=20$ 个候选基中，19 个是非奇异的，其中 8 个基本解是可行的；检验数检验无需考察其他基本解就识别出了最优的那个。

*（本例的一个早期版本在重新整理第一个和第三个约束时，把 $x_2+x_3$ 和 $2x_2+x_5$ 作为非基部分；对基 $\{x_1,x_3,x_6\}$，非基部分应为 $x_2+x_4$ 和 $2x_2$，并且非负条件的指标范围是 $1,\dots,6$。）*
:::

当检验不通过时，约化形式指明了下一步该怎么做。在最大化问题中，$\bar c_j>0$ 的非基变量从零开始增加时以变化率 $\bar c_j$ 改进 $z$；它可以一直增加，直到某个基变量（按 {eq}`eq-03-02-reduced-form` 的约束行变化）降为零。该基变量出基，$x_j$ 进基，从而得到一个目标函数值更好（或在退化时相等）的新 BFS。重复这一步骤就是{doc}`03-03-simplex-method`中的**单纯形法**（simplex method）；从几何上看，它沿着使目标函数改进的边，从 $P$ 的一个角点走到相邻的角点。

## 本章回顾

### 小结

- 两变量 LP 可以用图解法求解：画出可行域并平移目标函数等值线（{prf:ref}`alg-03-02-graphical`）；产品组合问题的最优解位于角点 $(\tfrac{45}{7},\tfrac{30}{7})$，$z^*=\tfrac{36000}{7}$（{prf:ref}`ex-03-02-graphical`、{numref}`fig-03-02-product-mix-region`）。松弛量、剩余量和紧约束描述了资源的使用情况（{prf:ref}`def-03-02-slack`）。
- LP 要么有唯一最优解，要么有多重最优解，要么无界，要么不可行（{numref}`fig-03-02-graphical-cases`）。
- 每个 LP 都可以化为标准形 $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$（{prf:ref}`prop-03-02-standard-form`）；其可行集是凸多面体（{prf:ref}`prop-03-02-convex-polyhedron`）。
- 基是 $m$ 个线性无关的列构成的集合；其基本解为 $\vx_B=\mB^{-1}\vb$、$\vx_N=\vzero$，若 $\vx_B\ge\vzero$ 则为 BFS（{prf:ref}`def-03-02-basic-solution`）；有基变量为零的 BFS 是退化的（{prf:ref}`def-03-02-degenerate`）。
- 极点恰好就是 BFS，并且至多有 $\binom nm$ 个（{prf:ref}`thm-03-02-extreme-bfs`、{prf:ref}`cor-03-02-finite`）。
- 若 LP 有最优解，则它有最优 BFS（{prf:ref}`thm-03-02-fundamental`）。
- 在约化形式 {eq}`eq-03-02-reduced-form` 中，目标函数为 $\vc_B^\T\mB^{-1}\vb+\bar\vc_N^\T\vx_N$，检验数为 $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$；满足 $\bar\vc_N\ge\vzero$（最小化）或 $\bar\vc_N\le\vzero$（最大化）的 BFS 是最优的（{prf:ref}`thm-03-02-optimality-test`、{prf:ref}`ex-03-02-reduced-form`）。

### 常见错误

- 从图上读取最优解时不精确计算两条紧约束线的交点；图只能给出是哪个角点，而不能给出其坐标。
- 因为可行域无界就断定 LP 无界；目标函数还必须沿某条射线改进。
- 忘记检验 $\vx_B\ge\vzero$：基本解可能位于可行域之外，并且目标函数值比最优值更好。
- 选取 $m$ 列时不检验它们是否线性无关；奇异的列组不是基。
- 在退化的 BFS 处把检验数条件当作必要条件，或者对最大化问题使用最小化的规则（$\bar c_j\ge0$）。
- 变量经过缩放（百箱、千桶）时，用错误的单位报告松弛量。

## 参考文献

- Bazaraa, M. S., Jarvis, J. J., & Sherali, H. D. (2010). *Linear
  programming and network flows* (4th ed.). Wiley. (Ch. 2, §§2.4–2.6:
  convex sets, polyhedral sets, extreme points; ch. 3, §§3.1–3.4: basic
  feasible solutions, key to the simplex method.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Ch. 1, §1.4: graphical
  representation; ch. 2, §§2.1–2.6: polyhedra, extreme points, vertices and
  basic feasible solutions, standard form, degeneracy, existence and
  optimality of extreme points.)
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Ch. 7: the geometry of linear programs.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Ch. 3, §3.1: graphical solution;
  ch. 5, §5.1: foundations of the simplex method.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Ch. 2, §2.2: graphical LP solution; ch. 3, §§3.1–3.2: LP in
  equation form, from graphical to algebraic solution; ch. 7: from extreme
  points to basic solutions.)
