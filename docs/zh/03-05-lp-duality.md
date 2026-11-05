# 对偶理论

## 概述

每个线性规划问题都伴随着由同样数据构成的另一个线性规划问题，即它的**对偶问题**（dual problem）。原来的（原始）问题选择各项活动的水平，对偶问题则为资源或需求定价，两个问题互为界：对极小化问题而言，任何可行的价格向量都给出费用的一个下界，而任何可行的计划都给出价格价值的一个上界。本章从定价论证出发推导对偶问题，给出写出任意线性规划问题的对偶问题的规则，并证明几个主要定理：**弱对偶性**（weak duality，即上述界）、**强对偶性**（strong duality，两个界在最优处相遇，用单纯形法证明）和**互补松弛性**（complementary slackness，哪些约束与变量可以同时为正）。本章说明单纯形法的最优性检验就是对偶可行性的检验，说明如何从最终单纯形表中读出最优对偶值，以及为什么这些值是**影子价格**（shadow price），即右端项改变时最优目标函数值的变化率。本章最后讨论 Farkas 引理（Farkas' lemma），这是由对偶性导出的一个择一定理，在非线性规划的 Karush–Kuhn–Tucker（KKT）条件中还会再次用到。本章以{doc}`03-02-lp-geometry`中的基与基本定理，以及{doc}`03-03-simplex-method`和{doc}`03-04-two-phase-big-m`中的单纯形表与检验数为基础；{ref}`sec-03-05-dual-rules`中的拉格朗日观点在{doc}`02-03-convex-optimization-duality`中针对一般凸问题展开。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 把线性规划问题的对偶问题解释为一个定价问题，并借助规则表写出任意线性规划问题的对偶问题；
- 证明弱对偶性，并用它给最优值定界和证明最优性；
- 叙述强对偶性，从最优单纯形基出发证明它，并列出原问题与对偶问题可能出现的状态组合；
- 叙述并证明互补松弛性，并用它由原问题的解求出对偶问题的解，反之亦然；
- 判断一个基是否原始可行、对偶可行或最优，并从最终单纯形表中读出最优对偶值；
- 把对偶值解释为影子价格和隐含成本，并说明这种解释何时成立；
- 由线性规划对偶理论推导 Farkas 引理。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\vpi\in\R^m$, $\pi_i$ | 对偶变量，每个原始约束对应一个（价格） |
| $w=\vb^\T\vpi$ | 对偶目标函数值；$w^*$ 为其最优值 |
| $\mA_{i\cdot}$, $\mA_{\cdot j}$ | $\mA$ 的第 $i$ 行和第 $j$ 列；在图中 $\va_i=\mA_{i\cdot}^\T$ 是约束 $i$ 的法向量 |
| $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ | $x_j$ 的检验数；对对偶可行的 $\vpi$，它也是**对偶约束** $j$ **的松弛量** |
| $s_i=\mA_{i\cdot}\vx-b_i$ | 原始约束 $i$ 的剩余量（在对称形式中，对可行的 $\vx$ 有 $\ge0$） |
| $t_j$ | 对偶单纯形表中对偶约束 $j$ 的松弛变量 |
| $z^*(\vb)$ | 作为右端项函数的最优值 |

(sec-03-05-diet-pricing)=
## 每个线性规划问题背后的定价问题

对偶的思想在一个小型费用极小化问题中最容易看清。

:::{prf:example} 营养配餐与药片制造商的定价
:label: ex-03-05-diet

一个家庭希望用六种食物配出最便宜的组合，至少提供 9 个单位的维生素 A 和 19 个单位的维生素 B。{numref}`tab-03-05-diet-data`给出了每种食物每千克的维生素含量及其价格。

```{list-table} 营养配餐问题的数据
:name: tab-03-05-diet-data
:header-rows: 1
:stub-columns: 1

* - 食物 $j$
  - 1
  - 2
  - 3
  - 4
  - 5
  - 6
  - 需求量（单位）
* - 维生素 A（单位/千克）
  - 1
  - 0
  - 2
  - 2
  - 1
  - 2
  - 9
* - 维生素 B（单位/千克）
  - 0
  - 1
  - 3
  - 1
  - 3
  - 2
  - 19
* - 价格 $c_j$（美分/千克）
  - 35
  - 30
  - 60
  - 50
  - 27
  - 22
  -
```

一家制造商以药片形式出售维生素 A 和维生素 B，希望这个家庭从食物转向药片。它把维生素 A 定价为每单位 $\pi_1$ 美分，维生素 B 定价为每单位 $\pi_2$ 美分。试建立家庭购买最便宜配餐的问题，以及制造商在没有任何食物比药片更便宜地提供维生素的条件下选择价格以使收入最大的问题。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

设 $x_j$ 为购买食物 $j$ 的千克数，家庭要求解

$$
\begin{aligned}
\min_{\vx}\ & z=35x_1+30x_2+60x_3+50x_4+27x_5+22x_6\\
\st\ & x_1+2x_3+2x_4+x_5+2x_6\ge 9\\
& x_2+3x_3+x_4+3x_5+2x_6\ge 19\\
& x_1,\dots,x_6\ge 0,
\end{aligned}
$$ (eq-03-05-diet-primal)

即 $\min\{\vc^\T\vx:\mA\vx\ge\vb,\ \vx\ge\vzero\}$，其中 $\mA=\begin{pmatrix}1&0&2&2&1&2\\0&1&3&1&3&2\end{pmatrix}$，$\vb=(9,19)^\T$，$\vc=(35,30,60,50,27,22)^\T$。

给定药片价格 $\pi_1$ 和 $\pi_2$，家庭会把药片与每种食物作比较：一千克食物 3 提供 2 个单位的 A 和 3 个单位的 B，按药片价格值 $2\pi_1+3\pi_2$ 美分，而它的价格是 60 美分。若 $2\pi_1+3\pi_2>60$，则食物 3 是这些维生素的更便宜来源，家庭会继续购买它。只有当没有任何食物比药片更便宜，即对每种食物 $j$ 都有 $\pi_1a_{1j}+\pi_2a_{2j}\le c_j$ 时，价格才是**有竞争力的**。转向药片的家庭恰好购买其需求量，即 9 个单位的 A 和 19 个单位的 B，所以制造商在有竞争力的价格中极大化其收入 $9\pi_1+19\pi_2$：

$$
\begin{aligned}
\max_{\vpi}\ & w=9\pi_1+19\pi_2\\
\st\ & \pi_1\le 35\\
& \pi_2\le 30\\
& 2\pi_1+3\pi_2\le 60\\
& 2\pi_1+\pi_2\le 50\\
& \pi_1+3\pi_2\le 27\\
& 2\pi_1+2\pi_2\le 22\\
& \pi_1,\pi_2\ge 0,
\end{aligned}
$$ (eq-03-05-diet-dual)

即 $\max\{\vb^\T\vpi:\mA^\T\vpi\le\vc,\ \vpi\ge\vzero\}$。

这两个问题都将在本章后面求解。家庭最便宜的配餐是食物 5 $x_5=5$ 千克、食物 6 $x_6=2$ 千克，费用 $z^*=179$ 美分；制造商的最优价格是每单位 $\vpi^*=(3,8)$ 美分，收入 $w^*=9\cdot3+19\cdot8=179$ 美分。
:::

{prf:ref}`ex-03-05-diet`中的定价问题 {eq}`eq-03-05-diet-dual` 是营养配餐问题 {eq}`eq-03-05-diet-primal` 的**对偶问题**。这一对问题的若干特征对每个原始–对偶问题对（primal–dual pair）都成立：

1. 一个问题的约束矩阵是另一个问题约束矩阵的转置；
2. 一个问题的右端项是另一个问题的目标函数系数；
3. 每个原始约束对应一个对偶变量（需求量的单位价格），每个原始变量对应一个对偶约束（该食物的价格上限）；
4. 带 $\ge$ 约束的极小化问题与带 $\le$ 约束的极大化问题配对；
5. 制造商的收入永远不会超过最便宜配餐的费用。任何满足需求的配餐 $\vx$ 和任何有竞争力的价格 $\vpi$ 都满足 $9\pi_1+19\pi_2\le\sum_j c_jx_j$，因为配餐至少含有所需的维生素，而每千克食物的价格至少等于其维生素含量按药片价格计算的价值。下面的定理将使这一点严格化，并证明两个最优值实际上相等。

(sec-03-05-canonical-pair)=
## 对称形式的对偶问题对

:::{prf:definition} 对称形式的原始–对偶问题对
:label: def-03-05-canonical-pair

设 $\mA\in\R^{m\times n}$，$\vb\in\R^m$，$\vc\in\R^n$。**原问题**（primal problem）(P) 及其**对偶问题** (D) 为

$$
\begin{aligned}
\text{(P)}\quad\min_{\vx}\ & z=\vc^\T\vx\\
\st\ & \mA\vx\ge\vb\\
& \vx\ge\vzero\\[1ex]
\text{(D)}\quad\max_{\vpi}\ & w=\vb^\T\vpi\\
\st\ & \mA^\T\vpi\le\vc\\
& \vpi\ge\vzero
\end{aligned}
$$ (eq-03-05-canonical-pair)

对偶问题对每个原始约束 $\mA_{i\cdot}\vx\ge b_i$ 有一个变量 $\pi_i$，对每个原始变量 $x_j$ 有一个约束 $\vpi^\T\mA_{\cdot j}\le c_j$。从右往左读，同一个定义给出了带 $\le$ 约束的极大化问题的对偶问题：$\max\{\vc^\T\vx:\mA\vx\le\vb,\vx\ge\vzero\}$ 的对偶问题为 $\min\{\vb^\T\vpi:\mA^\T\vpi\ge\vc,\vpi\ge\vzero\}$。
:::

对称形式的问题对是一个定义，但其他任何线性规划问题的对偶问题都可以由它导出，因为每个线性规划问题都可以改写成对称形式。对于{doc}`03-02-lp-geometry`中的标准形，把每个方程 $\mA\vx=\vb$ 写成两个不等式 $\mA\vx\ge\vb$ 和 $-\mA\vx\ge-\vb$。对这两组约束分别引入对偶变量 $\vpi^+\ge\vzero$ 和 $\vpi^-\ge\vzero$，对称形式的对偶问题为 $\max\{\vb^\T(\vpi^+-\vpi^-):\mA^\T(\vpi^+-\vpi^-)\le\vc\}$，而起作用的只是差 $\vpi=\vpi^+-\vpi^-$，它可以取任意符号。因此，标准形问题 (P) 的对偶问题 (D) 含有自由变量：

$$
\begin{aligned}
\text{(P)}\quad\min_{\vx}\ & \vc^\T\vx\\
\st\ & \mA\vx=\vb\\
& \vx\ge\vzero\\[1ex]
\text{(D)}\quad\max_{\vpi}\ & \vb^\T\vpi\\
\st\ & \mA^\T\vpi\le\vc\\
& \vpi\ \text{自由。}
\end{aligned}
$$ (eq-03-05-standard-pair)

一个线性规划问题可以用许多等价的方式改写，所以自然要问：对偶问题是否依赖于原问题的写法？答案是否定的。

:::{prf:proposition} 等价问题的对偶问题等价
:label: prop-03-05-equivalent-duals

若 (P$'$) 是由 (P) 通过引入松弛变量或剩余变量、把自由变量写成两个非负变量之差、把某个约束乘以 $-1$ 或把一个方程换成两个不等式而得到的，则 (P$'$) 的对偶问题 (D$'$) 与 (P) 的对偶问题 (D) 等价：两者有相同的可行价格（相差同一个变量替换）和相同的目标函数值。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

我们证明最常用的情形：(P) 是对称形式的原问题，(P$'$) 加入了剩余变量，即 $\min\{\vc^\T\vx+\vzero^\T\vs:\mA\vx-\mI\vs=\vb,\ \vx\ge\vzero,\ \vs\ge\vzero\}$。(P$'$) 是约束矩阵为 $(\mA\ \ -\mI)$、费用向量为 $(\vc,\vzero)$ 的标准形，所以由 {eq}`eq-03-05-standard-pair`，其对偶问题为

$$
\begin{aligned}
\max_{\vpi}\ & \vb^\T\vpi\\
\st\ & \mA^\T\vpi\le\vc\\
& -\mI\vpi\le\vzero\\
& \vpi\ \text{自由}.
\end{aligned}
$$

约束 $-\mI\vpi\le\vzero$ 就是 $\vpi\ge\vzero$，所以 (D$'$) 就是 (D)。其他变换可以用同样的方法验证：拆分 $x_j=x_j^+-x_j^-$ 给出两个对偶约束 $\vpi^\T\mA_{\cdot j}\le c_j$ 和 $-\vpi^\T\mA_{\cdot j}\le -c_j$，即一个方程；把约束乘以 $-1$ 改变其对偶变量的符号；拆分方程的情形上面已经讨论过。
:::

该命题还说明对偶问题在这种改写的意义下是唯一的，所以我们可以说*该*对偶问题。把这一构造应用两次就回到原问题。

:::{prf:theorem} 对偶问题的对偶问题是原问题
:label: thm-03-05-dual-of-dual

{eq}`eq-03-05-canonical-pair` 中 (D) 的对偶问题是 (P)。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

把 (D) 写成对称形式的极小化问题：$\min\{-\vb^\T\vpi:-\mA^\T\vpi\ge-\vc,\ \vpi\ge\vzero\}$。其约束矩阵为 $-\mA^\T$，右端项为 $-\vc$，费用向量为 $-\vb$。由{prf:ref}`def-03-05-canonical-pair`，它的以 $\vx$ 为变量的对偶问题为 $\max\{-\vc^\T\vx:(-\mA^\T)^\T\vx\le-\vb,\ \vx\ge\vzero\}$，即 $\max\{-\vc^\T\vx:\mA\vx\ge\vb,\ \vx\ge\vzero\}$，这与 $\min\{\vc^\T\vx:\mA\vx\ge\vb,\ \vx\ge\vzero\}$ 相同。
:::

因此两个问题地位对等。哪一个称为原问题只取决于从哪里出发，下面的每个命题都有一个角色互换的镜像版本。

(sec-03-05-dual-rules)=
## 写出任意线性规划问题对偶问题的规则

把{prf:ref}`prop-03-05-equivalent-duals`应用于每一类约束和变量，就得到{numref}`tab-03-05-dual-rules`中的规则。当约束的方向与目标函数推动的方向一致时，称该约束具有**自然类型**（natural type）：在极小化问题中为 $\ge$（需求限制了费用能降多低），在极大化问题中为 $\le$。自然类型的约束对应非负的对偶变量，方程对应自由的对偶变量，相反类型的约束对应非正的对偶变量。对称地，非负的原始变量给出对偶问题中自然类型的对偶约束，自由变量给出方程，非正变量给出相反类型的约束。

```{list-table} 写出对偶问题的规则。原问题为极小化时从左往右读，原问题为极大化时从右往左读。
:name: tab-03-05-dual-rules
:header-rows: 1

* - 极小化问题
  - 极大化问题
* - 目标函数系数 $\vc$
  - 右端项 $\vc$
* - 右端项 $\vb$
  - 目标函数系数 $\vb$
* - 约束 $i$：$\mA_{i\cdot}\vx\ge b_i$
  - 变量 $\pi_i\ge0$
* - 约束 $i$：$\mA_{i\cdot}\vx\le b_i$
  - 变量 $\pi_i\le0$
* - 约束 $i$：$\mA_{i\cdot}\vx=b_i$
  - 变量 $\pi_i$ 自由
* - 变量 $x_j\ge0$
  - 约束 $j$：$\vpi^\T\mA_{\cdot j}\le c_j$
* - 变量 $x_j\le0$
  - 约束 $j$：$\vpi^\T\mA_{\cdot j}\ge c_j$
* - 变量 $x_j$ 自由
  - 约束 $j$：$\vpi^\T\mA_{\cdot j}=c_j$
```

### 规则为何成立：拉格朗日界

这些规则不是随意规定的；它们恰好是使价格给出有效界的条件。考虑一个含有三种类型约束和三种符号变量的极小化问题，把约束 $i$ 乘以价格 $\pi_i$。我们希望对每个可行的 $\vx$ 有

$$
\pi_i\,(b_i-\mA_{i\cdot}\vx)\le 0 ,
$$

这样把这些项加到费用上只会使其减小。对于 $\ge$ 约束，$b_i-\mA_{i\cdot}\vx\le0$，所以需要 $\pi_i\ge0$；对于 $\le$ 约束，需要 $\pi_i\le0$；对于方程，无论 $\pi_i$ 的符号如何，该项都为零。在这样的价格下，每个可行的 $\vx$ 满足

$$
\begin{aligned}
\vc^\T\vx\ \ge\ L(\vx,\vpi)&=\vc^\T\vx+\vpi^\T(\vb-\mA\vx)\\
&=\vb^\T\vpi+(\vc-\mA^\T\vpi)^\T\vx .
\end{aligned}
$$ (eq-03-05-lagrangian)

右端仍依赖于 $\vx$。它对每个可行的 $\vx$ 都至少为 $\vb^\T\vpi$，当且仅当每一项 $(c_j-\vpi^\T\mA_{\cdot j})x_j$ 对 $x_j$ 的每个允许值都非负：$x_j\ge0$ 时 $c_j-\vpi^\T\mA_{\cdot j}\ge0$，$x_j\le0$ 时 $\le0$，$x_j$ 自由时 $=0$。这些正是{numref}`tab-03-05-dual-rules`中的对偶约束，在这些约束下 $\vc^\T\vx\ge\vb^\T\vpi$。对偶问题要求的就是这种下界中最好（最大）的一个。

函数 $L$ 是该线性规划问题的**拉格朗日函数**（Lagrangian），而在 $\vx$ 的符号限制下取下确界所得的 $q(\vpi)=\inf_{\vx}L(\vx,\vpi)$，当 $\vpi$ 满足对偶约束时等于 $\vb^\T\vpi$，否则等于 $-\infty$。因此，对偶线性规划问题就是{ref}`sec-02-03-dual-problem`中的拉格朗日对偶问题 $\max_{\vpi}q(\vpi)$ 在线性函数情形下的特例。与一般凸规划的一个区别是，线性规划的强对偶性（{prf:ref}`thm-03-05-strong-duality`）不需要 Slater 条件（{prf:ref}`thm-02-03-slater`）之类的约束规范。

有两个简单的检查可以防止符号错误。第一，对偶问题的对偶问题必须回到原问题（{prf:ref}`thm-03-05-dual-of-dual`）。第二，类型错误的约束总可以先乘以 $-1$ 再应用对称形式的问题对；其对偶变量随之以相反的符号出现，这就是 $\pi_i\le0$ 规则。包括一个含全部三类约束和两个自由变量的原问题、以及一个分块形式的原问题在内的例题，见{doc}`08-03-examples-duality`。

(sec-03-05-weak-duality)=
## 弱对偶性

:::{prf:theorem} 弱对偶性
:label: thm-03-05-weak-duality

若 $\vx$ 对 {eq}`eq-03-05-canonical-pair` 中的 (P) 可行，$\vpi$ 对 (D) 可行，则

$$
\vc^\T\vx\ \ge\ \vpi^\T\mA\vx\ \ge\ \vb^\T\vpi .
$$ (eq-03-05-weak-chain)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

因为 $\vx\ge\vzero$ 且 $\mA^\T\vpi\le\vc$，$\vc^\T\vx-\vpi^\T\mA\vx=\sum_j(c_j-\vpi^\T\mA_{\cdot j})x_j$ 的每一项都是两个非负数之积，所以 $\vc^\T\vx\ge\vpi^\T\mA\vx$。因为 $\vpi\ge\vzero$ 且 $\mA\vx\ge\vb$，$\vpi^\T\mA\vx-\vb^\T\vpi=\sum_i\pi_i(\mA_{i\cdot}\vx-b_i)$ 的每一项也都是两个非负数之积，所以 $\vpi^\T\mA\vx\ge\vb^\T\vpi$。
:::

两个符号条件都是必需的：第一个不等式需要 $\vx\ge\vzero$，第二个需要 $\vpi\ge\vzero$。对于按{numref}`tab-03-05-dual-rules`的规则写出的原始–对偶问题对，同样的论证可以通过 {eq}`eq-03-05-lagrangian` 逐项进行，所以弱对偶性对每个问题对都成立，而不仅仅对对称形式成立。把证明的两部分相加，得到恒等式

$$
\begin{gathered}
\vc^\T\vx-\vb^\T\vpi=\sum_{j=1}^{n}\bar c_j\,x_j+\sum_{i=1}^{m}\pi_i\,s_i,\\
\bar c_j=c_j-\vpi^\T\mA_{\cdot j},\\
s_i=\mA_{i\cdot}\vx-b_i ,
\end{gathered}
$$ (eq-03-05-duality-gap)

它把**对偶间隙**（duality gap）$\vc^\T\vx-\vb^\T\vpi$ 写成若干乘积之和，每个乘积都是一个问题的变量与另一个问题中对应约束的松弛量（或剩余量）之积。这个恒等式在互补松弛性中还会用到。

{numref}`fig-03-05-weak-duality-line`展示了{prf:ref}`ex-03-05-diet`中营养配餐问题的弱对偶性。三个有竞争力的价格向量 $(11,0)$、$(5,6)$ 和 $(0,9)$ 分别给出 99、159 和 171 美分的收入；三个可行配餐分别给出 191、209 和 243 美分的费用（5 千克食物 6 加 3 千克食物 5；只用 9.5 千克食物 6；只用 9 千克食物 5）。每个收入都位于每个费用的左边。

```{figure} figures/03-05-weak-duality-line.svg
:name: fig-03-05-weak-duality-line
:width: 600px
:class: dark-light
:alt: 数轴：左边是对偶目标函数值 99、159 和 171，右边是原始目标函数值 191、209 和 243，两者之间是共同的最优值 179

营养配餐问题的弱对偶性（单位：美分）。可行价格的对偶值 $w=\vb^\T\vpi$（圆点）永远不超过可行配餐的原始费用 $z=\vc^\T\vx$（方块）；两者的最优值在 $z^*=w^*=179$ 处相遇（菱形）。
```

弱对偶性有几个直接推论。

:::{prf:corollary} 界与证书
:label: cor-03-05-optimality-certificate

设 (P) 和 (D) 是一个原始–对偶问题对，(P) 为极小化问题。

1. 任何原始可行解 $\vx$ 的费用 $\vc^\T\vx$ 都是 (D) 最优值的上界，任何对偶可行解 $\vpi$ 的值 $\vb^\T\vpi$ 都是 (P) 最优值的下界。
2. 若 $\vx$ 原始可行，$\vpi$ 对偶可行，且 $\vc^\T\vx=\vb^\T\vpi$，则 $\vx$ 是 (P) 的最优解，$\vpi$ 是 (D) 的最优解。
3. 若 (P) 无界（$z\to-\infty$），则 (D) 不可行。若 (D) 无界（$w\to+\infty$），则 (P) 不可行。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

第 1 部分就是从两端分别读 {eq}`eq-03-05-weak-chain`。对第 2 部分，每个原始可行解 $\vx'$ 都满足 $\vc^\T\vx'\ge\vb^\T\vpi=\vc^\T\vx$，所以 $\vx$ 是最优的；对 $\vpi$ 的论证相同。对第 3 部分，若存在对偶可行解 $\vpi$，它将给每个原始费用一个有限的下界 $\vb^\T\vpi$，这与 $z\to-\infty$ 矛盾；第二句是其镜像。
:::

第 2 部分是一种实用的**最优性证书**（optimality certificate）。声称已经求解了某个线性规划问题的人，可以通过给出一个原始解和一个对偶解来证明这一点：检查者验证原始解满足原始约束、对偶解满足对偶约束，并且两个目标函数值相等。这些检查都不需要求解任何问题。对营养配餐问题，$\vx=(0,0,0,0,5,2)$ 满足两项需求（$5+2\cdot2=9$ 和 $3\cdot5+2\cdot2=19$），$\vpi=(3,8)$ 满足全部六个价格上限（在{prf:ref}`ex-03-05-diet-cs`中验证），并且两个值都是 179 美分，所以两者都是最优的。

第 3 部分的逆命题不成立：原问题不可行并不迫使对偶问题无界，因为两个问题可能都不可行。{prf:ref}`ex-08-03-both-infeasible`就是一个例子。

(sec-03-05-strong-duality)=
## 强对偶性

弱对偶性没有回答最好的原始值与最好的对偶值之间是否还有间隙。对线性规划而言，没有间隙。

:::{prf:theorem} 强对偶性
:label: thm-03-05-strong-duality

若原始–对偶问题对中的任一问题有最优解，则另一个问题也有最优解，并且两个最优值相等：$z^*=w^*$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

由{prf:ref}`prop-03-05-equivalent-duals`，可以把原问题写成标准形 $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$，其对偶问题为 {eq}`eq-03-05-standard-pair`。还可以假设 $\mA$ 行满秩（秩为 $m$），因为删去冗余方程不改变可行集（{prf:ref}`thm-03-04-phase-one`）；删去一个方程相当于把其对偶变量固定为零，不改变对偶最优值。

设原问题有最优解。由线性规划基本定理（{prf:ref}`thm-03-02-fundamental`）和带防循环规则的单纯形法（{prf:ref}`thm-03-03-bland`），存在最优基 $\mB$：$\vx_B^*=\mB^{-1}\vb\ge\vzero$，$\vx_N^*=\vzero$，且所有检验数非负，

$$
\bar c_j=c_j-\vc_B^\T\mB^{-1}\mA_{\cdot j}\ge 0\qquad\text{对 } j=1,\dots,n .
$$

定义**单纯形乘子**（simplex multipliers）

$$
\vpi^{*\T}=\vc_B^\T\mB^{-1}.
$$ (eq-03-05-simplex-multipliers)

于是对每个 $j$ 有 $\bar c_j=c_j-\vpi^{*\T}\mA_{\cdot j}\ge0$，即 $\mA^\T\vpi^*\le\vc$，所以 $\vpi^*$ 对对偶问题（其变量是自由的）可行。它的值为

$$
\begin{aligned}
\vb^\T\vpi^*&=\vc_B^\T\mB^{-1}\vb\\
&=\vc_B^\T\vx_B^*=\vc^\T\vx^* .
\end{aligned}
$$

由{prf:ref}`cor-03-05-optimality-certificate`的第 2 部分，$\vpi^*$ 是对偶问题的最优解，且 $w^*=z^*$。

若是对偶问题有最优解，则把对偶问题写成标准形的极小化问题，再对它应用刚才的论证；由{prf:ref}`thm-03-05-dual-of-dual`，它的对偶问题就是原问题，因此原问题有最优解且最优值相同。
:::

这个定理可以追溯到 von Neumann 关于博弈的工作；关于线性规划的第一个公开发表的证明是 Gale 等（1951）给出的。这里通过最优单纯形基给出的证明属于 Dantzig（1963）。

该证明不只是证明了存在性：它由最优基构造出最优对偶解 $\vpi^{*\T}=\vc_B^\T\mB^{-1}$，并表明单纯形法的最优性检验 $\bar c\ge\vzero$ *就是*这个 $\vpi^*$ 的对偶可行性。{ref}`sec-03-05-dual-feasible-basis`将展开这一观察。

强对偶性与弱对偶性一起决定了哪些结果组合是可能的。一个线性规划问题要么不可行，要么无界，要么有最优解（线性规划基本定理），所以一个问题对在理论上有九种组合。

```{list-table} 原始–对偶问题对的可能结果（原问题为极小化）
:name: tab-03-05-status
:header-rows: 1
:stub-columns: 1

* - 原问题 \ 对偶问题
  - 有最优解
  - 无界（$w\to+\infty$）
  - 不可行
* - 有最优解
  - 可能，$z^*=w^*$
  - 不可能
  - 不可能
* - 无界（$z\to-\infty$）
  - 不可能
  - 不可能
  - 可能
* - 不可行
  - 不可能
  - 可能
  - 可能
```

第一行和第一列中的“不可能”来自强对偶性（一个问题有最优解，另一个也有）。由{prf:ref}`cor-03-05-optimality-certificate`的第 3 部分，原问题和对偶问题不可能都无界。其余三种组合都会出现：两者都有最优解（营养配餐问题）；一个无界而另一个不可行（把混合约束的例子{prf:ref}`ex-08-03-mixed-dual`中一个系数的符号改变）；两者都不可行（{prf:ref}`ex-08-03-both-infeasible`）。从表中可以读出两个有用的推论。

:::{prf:corollary} 可行性决定结果
:label: cor-03-05-feasible-pair

若一个问题对中的两个问题都可行，则两者都有最优解且最优值相等。若原问题可行而对偶问题不可行，则原问题无界；若对偶问题可行而原问题不可行，则对偶问题无界。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

若两者都可行，则都不无界（{prf:ref}`cor-03-05-optimality-certificate`的第 3 部分），所以各自都有最优解，且由{prf:ref}`thm-03-05-strong-duality`，两者的最优值相等。若原问题可行而对偶问题不可行，则原问题不可能有最优解（否则对偶问题也会有最优解），所以它无界。最后一句是其镜像。
:::

(sec-03-05-complementary-slackness)=
## 互补松弛性

间隙恒等式 {eq}`eq-03-05-duality-gap` 把 $\vc^\T\vx-\vb^\T\vpi$ 写成若干非负乘积之和。在最优解对处间隙为零，所以每个乘积都必须为零。

:::{prf:theorem} 互补松弛性
:label: thm-03-05-complementary-slackness

设 $\vx$ 对 {eq}`eq-03-05-canonical-pair` 中的 (P) 可行，$\vpi$ 对 (D) 可行。则 $\vx$ 和 $\vpi$ 都是最优解，当且仅当

$$
\begin{aligned}
\pi_i\,(\mA_{i\cdot}\vx-b_i)&=0 && i=1,\dots,m,\\
(c_j-\vpi^\T\mA_{\cdot j})\,x_j&=0 && j=1,\dots,n.
\end{aligned}
$$ (eq-03-05-cs)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

由 {eq}`eq-03-05-duality-gap`，$\vc^\T\vx-\vb^\T\vpi$ 是 {eq}`eq-03-05-cs` 中 $m+n$ 个乘积之和，对可行的 $\vx$ 和 $\vpi$，每个乘积都非负。若两者都最优，则由强对偶性 $\vc^\T\vx=\vb^\T\vpi$，所以和为零，每个乘积都为零。反之，若每个乘积都为零，则两个目标函数值相等，由{prf:ref}`cor-03-05-optimality-certificate`的第 2 部分，两个解都是最优的。
:::

用文字表述就是：在一对最优解处，只要一个问题的某个约束有正的松弛量，另一个问题中对应的变量就为零；只要某个变量为正，另一个问题中对应的约束就是紧的。下面三点说明可以防止常见的误读。

- 这些条件对**自由**变量没有任何要求，因为与之对应的约束是方程，没有松弛量。在标准形问题对 {eq}`eq-03-05-standard-pair` 中，{eq}`eq-03-05-cs` 只剩下第二行。
- 松弛量与其对应的变量可以**都**为零。互补松弛性说的是两者中至少有一个为零，而不是恰好一个为零。（对每个下标两者中都恰好有一个为正的解对称为*严格*互补的；只要两个问题都有最优解，这样的解对就存在（Goldman & Tucker, 1956），但这里不需要求出它。）
- 该定理适用于任何一对最优解，而不仅仅适用于基本最优解。当一个问题有多个最优解时，原问题的每个最优解都与对偶问题的每个最优解互补。

最常见的用法是由一个问题的解求出另一个问题的解。

:::{prf:example} 由最优配餐求药片价格
:label: ex-03-05-diet-cs

{prf:ref}`ex-03-05-diet`的最优配餐为 $\vx^*=(0,0,0,0,5,2)$ 千克。求最优价格。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

由于 $x_5^*=5>0$ 且 $x_6^*=2>0$，对偶约束 5 和 6 必须是紧的：

$$
\begin{aligned}
\pi_1+3\pi_2&=27,\\
2\pi_1+2\pi_2&=22 .
\end{aligned}
$$

用第二个方程减去第一个方程的两倍，得 $-4\pi_2=-32$，所以 $\pi_2=8$，$\pi_1=27-24=3$：每单位 $\vpi=(3,8)$ 美分。

互补松弛性只给出候选解；还必须验证它是对偶可行的。六个价格上限为

$$
\begin{aligned}
3&\le35,\\
8&\le30,\\
2\cdot3+3\cdot8=30&\le60,\\
2\cdot3+8=14&\le50,\\
3+3\cdot8=27&\le27,\\
2\cdot3+2\cdot8=22&\le22,
\end{aligned}
$$

且 $\vpi\ge\vzero$，所以 $\vpi$ 可行。{eq}`eq-03-05-cs` 的其余条件也成立：两个价格都为正，而两个维生素约束在 $\vx^*$ 处都是紧的（剩余量 $s_1=0+2\cdot0+2\cdot0+5+2\cdot2-9=0$，$s_2=0+3\cdot0+0+3\cdot5+2\cdot2-19=0$）。由{prf:ref}`thm-03-05-complementary-slackness`，$\vx^*$ 和 $\vpi^*=(3,8)$ 都是最优解，且 $z^*=w^*=27\cdot5+22\cdot2=9\cdot3+19\cdot8=179$ 美分。

如果验证不通过，$\vx^*$ 就不是最优的。这个论证也可以反过来进行：仅由 $\vpi^*=(3,8)$ 可知，两个维生素约束都必须是紧的（两个价格都为正），并且只能购买食物 5 和 6（只有它们的价格上限是紧的），所以 $x_5+2x_6=9$，$3x_5+2x_6=19$，由此得 $x_5=5$，$x_6=2$。
:::

{numref}`fig-03-05-diet-dual-region`展示了营养配餐例子的对偶问题。最优价格位于食物 5 和食物 6 的价格上限线的交点，这两种食物正是最优配餐中的食物；食物 3 和 4 的上限线不经过最优点（食物 1 和 2 的上限线在图窗之外），这些食物不被购买。

```{figure} figures/03-05-diet-dual-region.svg
:name: fig-03-05-diet-dual-region
:width: 520px
:class: dark-light
:alt: 对偶营养配餐问题在 pi1-pi2 平面上的可行域，最优价格向量 (3, 8) 位于食物 5 和食物 6 两条直线的交点

药片制造商的可行价格。直线 ($j$) 是食物 $j$ 的价格上限 $\vpi^\T\mA_{\cdot j}=c_j$。收入 $w=9\pi_1+19\pi_2$ 沿方向 $\vb$ 增大，在 $\vpi^*=(3,8)$ 处最大，该处只有食物 5 和 6 的价格上限是紧的；它们正是 $x_j^*>0$ 的食物。
```

### 几何解释

互补松弛性在原始空间中有一个直观的图像。考虑极大化问题 $\max\{\vc^\T\vx:\mA\vx\le\vb,\ \vx\ge\vzero\}$，设其最优解 $\vx^*$ 对每个 $j$ 都有 $x_j^*>0$，因而符号约束都不是有效约束。于是每个对偶约束都是紧的，$\mA^\T\vpi^*=\vc$，并且对每个在 $\vx^*$ 处非有效的约束都有 $\pi_i^*=0$。记 $\va_i=\mA_{i\cdot}^\T$ 为约束 $i$ 的外法向量，$I(\vx^*)$ 为有效约束的集合，上述结论就是

$$
\begin{aligned}
\vc&=\sum_{i\in I(\vx^*)}\pi_i^*\,\va_i,\\
\pi_i^*&\ge 0, \qquad i\in I(\vx^*) :
\end{aligned}
$$ (eq-03-05-cone-condition)

即在最优顶点处，目标函数的梯度位于有效约束的法向量所张成的锥中，而最优对偶值就是相应的权重。这是{prf:ref}`thm-06-03-kkt`中 KKT 条件的线性情形。产品组合问题展示了这一点。

:::{prf:example} 产品组合问题的对偶问题
:label: ex-03-05-product-mix-dual

在{prf:ref}`ex-03-01-product-mix`的产品组合问题中，一家成型车间生产果汁杯和鸡尾酒杯；$x_1$ 和 $x_2$ 是每周产量（单位：百箱），每百箱的利润分别为 500 美元和 450 美元。每百箱分别需要 6 小时和 5 小时机时，每周可用 60 小时；每百箱分别占用 10 立方米和 20 立方米仓储空间，可用 150 立方米；果汁杯每周至多能售出 8 百箱：

$$
\begin{aligned}
\max_{\vx}\ & z=500x_1+450x_2\\
\st\ & 6x_1+5x_2\le 60 &&\text{（机时，小时/周）}\\
& 10x_1+20x_2\le 150 &&\text{（仓储空间，立方米）}\\
& x_1\le 8 &&\text{（需求，百箱/周）}\\
& x_1,x_2\ge 0 .
\end{aligned}
$$

其最优解为 $\vx^*=(45/7,\,30/7)\approx(6.429,\,4.286)$，每周 $z^*=36000/7\approx\$5142.86$（{prf:ref}`ex-03-01-product-mix`）。写出对偶问题，用互补松弛性求出其最优解，并在 $\vx^*$ 处验证锥条件 {eq}`eq-03-05-cone-condition`。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**对偶问题。** 原问题是带 $\le$ 约束（自然类型）的极大化问题，所以由{numref}`tab-03-05-dual-rules`，对偶问题是一个极小化问题，每种资源对应一个非负变量，每种产品对应一个 $\ge$ 约束：

$$
\begin{aligned}
\min_{\vpi}\ & w=60\pi_1+150\pi_2+8\pi_3\\
\st\ & 6\pi_1+10\pi_2+\pi_3\ge 500\\
& 5\pi_1+20\pi_2\ge 450\\
& \pi_1,\pi_2,\pi_3\ge 0 .
\end{aligned}
$$

**用互补松弛性求对偶解。** $x_1^*$ 和 $x_2^*$ 都为正，所以两个对偶约束都是紧的。需求约束在 $\vx^*$ 处的松弛量为 $8-45/7=11/7>0$，所以 $\pi_3^*=0$。紧约束成为 $6\pi_1+10\pi_2=500$ 和 $5\pi_1+20\pi_2=450$。第一式的两倍减去第二式得 $7\pi_1=550$，所以

$$
\begin{gathered}
\pi_1^*=\frac{550}{7}\approx 78.57,\\
\pi_2^*=\frac{500-6\cdot\frac{550}{7}}{10}=\frac{20}{7}\approx 2.857,\\
\pi_3^*=0 .
\end{gathered}
$$

这些价格非负，且以等式满足两个对偶约束，所以它们是对偶可行的，并且 $w=60\cdot\frac{550}{7}+150\cdot\frac{20}{7}=\frac{33000+3000}{7}=\frac{36000}{7}=z^*$，这证明了两个解都是最优的。

**几何。** {numref}`fig-03-05-cs-cone`显示了最优顶点。有效约束是机时和仓储约束，法向量为 $\va_1=(6,5)$ 和 $\va_2=(10,20)$，而且确实有

$$
\begin{aligned}
\frac{550}{7}\begin{pmatrix}6\\5\end{pmatrix}+\frac{20}{7}\begin{pmatrix}10\\20\end{pmatrix}
&=\frac{1}{7}\begin{pmatrix}3300+200\\2750+400\end{pmatrix}\\
&=\begin{pmatrix}500\\450\end{pmatrix}=\vc ,
\end{aligned}
$$

所以 $\vc$ 位于 $\va_1$ 和 $\va_2$ 所张成的锥中。从 $\vx^*$ 出发的任何可行移动都有 $\va_1^\T\vd\le0$ 和 $\va_2^\T\vd\le0$，因此 $\vc^\T\vd=\pi_1^*\va_1^\T\vd+\pi_2^*\va_2^\T\vd\le0$：没有可行方向能改善目标函数。数 $550/7$ 和 $20/7$ 的经济含义是{ref}`sec-03-05-shadow-prices`的主题。

```{figure} figures/03-05-cs-cone.svg
:name: fig-03-05-cs-cone
:width: 520px
:class: dark-light
:alt: 产品组合问题的可行域及最优顶点；目标函数梯度 c 位于两个有效约束的法向量 a1 与 a2 之间

产品组合问题最优顶点处的互补松弛性。目标函数梯度 $\vc$ 位于有效约束 (1) 和 (2) 的法向量 $\va_1$、$\va_2$ 所张成的阴影锥中，权重等于最优对偶值；非有效的需求约束 (3) 的权重为零。虚线是最优目标函数等值线 $z=36000/7$。
```
:::

(sec-03-05-dual-feasible-basis)=
## 基的对偶可行性

强对偶性的证明给每个基都附加了一个对偶向量。现在我们把它用于（尚）不是最优的基。考虑标准形 $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$ 及其对偶问题 {eq}`eq-03-05-standard-pair`，设 $\mB$ 是一个基，非基列为 $\mN$。基本解为 $\vx_B=\mB^{-1}\vb$，$\vx_N=\vzero$；该基的对偶向量为 $\vpi^\T=\vc_B^\T\mB^{-1}$，它是 $\mB^\T\vpi=\vc_B$ 的唯一解，也就是把基变量对应的对偶约束写成等式所得方程组的唯一解。

:::{prf:definition} 原始可行基与对偶可行基
:label: def-03-05-dual-feasible-basis

若 $\mB^{-1}\vb\ge\vzero$，则称标准形问题的基 $\mB$ 是**原始可行**（primal feasible）的；若其乘子 $\vpi^\T=\vc_B^\T\mB^{-1}$ 满足全部对偶约束 $\mA^\T\vpi\le\vc$，则称其为**对偶可行**（dual feasible）的。由于基变量对应的对偶约束按构造以等式成立，一个基是对偶可行的当且仅当它的所有检验数都非负：

$$
\bar\vc_N^\T=\vc_N^\T-\vc_B^\T\mB^{-1}\mN\ \ge\ \vzero^\T .
$$ (eq-03-05-dual-feasible)

对极大化问题，对偶可行性意味着 $\bar\vc_N\le\vzero$。
:::

检验数 $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ 就是对偶约束 $j$ 在该基的乘子处的松弛量。因此{doc}`03-03-simplex-method`中的单纯形表一直带着对偶解：写成矩阵形式，基 $\mB$ 的单纯形表的约束行为 $(\mB^{-1}\mA\ \mid\ \mB^{-1}\vb)$，目标行为

$$
\bigl(\vc^\T-\vpi^\T\mA\ \mid\ -\vpi^\T\vb\bigr)=\bigl(\bar\vc^\T\ \mid\ -z\bigr),
$$

所以其最后一行列出了对偶松弛量，右下角存放两个基本解的共同值 $\vc_B^\T\mB^{-1}\vb=\vb^\T\vpi$ 的相反数。（把初始表变换为这张表的加边逆矩阵就是 $\mT$ 矩阵 {eq}`eq-03-06-bordered-inverse`，在{doc}`03-07-sensitivity-analysis`中还会用到。）

:::{prf:proposition} 既原始可行又对偶可行的基是最优基
:label: prop-03-05-basis-optimal

若基 $\mB$ 既是原始可行的又是对偶可行的，则 $\vx_B=\mB^{-1}\vb$，$\vx_N=\vzero$ 是原问题的最优解，$\vpi^\T=\vc_B^\T\mB^{-1}$ 是对偶问题的最优解。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

两个向量都可行，且 $\vc^\T\vx=\vc_B^\T\mB^{-1}\vb=\vpi^\T\vb$；应用{prf:ref}`cor-03-05-optimality-certificate`的第 2 部分即可。
:::

基对应的这对解自动是互补的：基变量有 $\bar c_j=0$，非基变量有 $x_j=0$，所以 {eq}`eq-03-05-cs` 中的每个乘积 $\bar c_jx_j$ 都为零。因此一个基总给出一对互补的解，唯一的问题是这两个解中哪一个可行。这是描述两种单纯形法最清楚的方式。

- **原始单纯形法**（primal simplex method）保持基原始可行（$\mB^{-1}\vb\ge\vzero$）且互补，并进行主元变换直到基也对偶可行（$\bar\vc\ge\vzero$）。负的检验数就是一个被违反的对偶约束，选它作为进基变量就是试图修复这个约束。
- **对偶单纯形法**（dual simplex method，{prf:ref}`alg-03-06-dual-simplex`）保持基对偶可行，并进行主元变换直到它也原始可行。

对给定基进行检验的例题，以及从对偶不可行的基经主元变换到达最优基的例题，见{ref}`sec-08-03-feasibility-tests`。

(sec-03-05-reading-duals)=
## 从最优单纯形表读出对偶解

单纯形法停止时，最优对偶值已经在最终单纯形表中，无须另行求解对偶问题。由 {eq}`eq-03-05-simplex-multipliers` 可以得到两种读法。

1. **由 $\mB^{-1}$ 读出。** 若初始表中含有单位矩阵（松弛变量列或人工变量列），则最终表中的相同列就是 $\mB^{-1}$，由此可算出 $\vpi^\T=\vc_B^\T\mB^{-1}$。
2. **由目标行读出。** 每一列都有 $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$。对于带符号的单位向量列，这直接给出一个 $\pi_i$，如{numref}`tab-03-05-reading-duals`所列。

```{list-table} 由最优表的目标行读出对偶值（$\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$；$\ve_i$ 是第 $i$ 个单位向量）
:name: tab-03-05-reading-duals
:header-rows: 1

* - 约束 $i$ 中的列
  - $\mA_{\cdot j}$
  - $c_j$
  - 元素 $\bar c_j$
  - 因此
* - $\le$ 行的松弛变量
  - $\ve_i$
  - $0$
  - $-\pi_i$
  - $\pi_i=-\bar c_j$
* - $\ge$ 行的剩余变量
  - $-\ve_i$
  - $0$
  - $\pi_i$
  - $\pi_i=\bar c_j$
* - 人工变量，极小化（费用 $+M$）
  - $\ve_i$
  - $M$
  - $M-\pi_i$
  - $\pi_i=M-\bar c_j$
* - 人工变量，极大化（费用 $-M$）
  - $\ve_i$
  - $-M$
  - $-M-\pi_i$
  - $\pi_i=-M-\bar c_j$
```

这些规则对极小化和极大化同样适用，因为它们都来自同一个公式 $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$。$\pi_i$ 的符号则遵循{numref}`tab-03-05-dual-rules`的规则：在产品组合问题（极大化，$\le$ 行）的最优表中，三个松弛变量列的目标行元素为 $-550/7$、$-20/7$ 和 $0$，所以 $\vpi^*=(550/7,\,20/7,\,0)$，与{prf:ref}`ex-03-05-product-mix-dual`中求得的一致。即使人工变量已经离开了问题，仍然可以读取人工变量的列：它们的元素依然满足该公式。

对营养配餐问题，原问题（标准形，剩余变量为 $x_7$ 和 $x_8$）的最优表和对偶问题（标准形，六个价格上限的松弛变量为 $t_1,\dots,t_6$）的最优表分别见{numref}`tab-03-05-diet-primal-tableau`和{numref}`tab-03-05-diet-dual-tableau`。

```{list-table} 营养配餐问题的最优表（极小化；$x_7$、$x_8$ 是维生素约束的剩余变量）
:name: tab-03-05-diet-primal-tableau
:header-rows: 1
:stub-columns: 1

* - 基变量
  - $x_1$
  - $x_2$
  - $x_3$
  - $x_4$
  - $x_5$
  - $x_6$
  - $x_7$
  - $x_8$
  - 右端项
* - $x_6$
  - $\tfrac34$
  - $-\tfrac14$
  - $\tfrac34$
  - $\tfrac54$
  - $0$
  - $1$
  - $-\tfrac34$
  - $\tfrac14$
  - $2$
* - $x_5$
  - $-\tfrac12$
  - $\tfrac12$
  - $\tfrac12$
  - $-\tfrac12$
  - $1$
  - $0$
  - $\tfrac12$
  - $-\tfrac12$
  - $5$
* - $-z$
  - $32$
  - $22$
  - $30$
  - $36$
  - $0$
  - $0$
  - $3$
  - $8$
  - $-179$
```

```{list-table} 对偶营养配餐问题的最优表（极大化；$t_j$ 是食物 $j$ 价格上限的松弛变量）
:name: tab-03-05-diet-dual-tableau
:header-rows: 1
:stub-columns: 1

* - 基变量
  - $\pi_1$
  - $\pi_2$
  - $t_1$
  - $t_2$
  - $t_3$
  - $t_4$
  - $t_5$
  - $t_6$
  - 右端项
* - $\pi_1$
  - $1$
  - $0$
  - $0$
  - $0$
  - $0$
  - $0$
  - $-\tfrac12$
  - $\tfrac34$
  - $3$
* - $\pi_2$
  - $0$
  - $1$
  - $0$
  - $0$
  - $0$
  - $0$
  - $\tfrac12$
  - $-\tfrac14$
  - $8$
* - $t_1$
  - $0$
  - $0$
  - $1$
  - $0$
  - $0$
  - $0$
  - $\tfrac12$
  - $-\tfrac34$
  - $32$
* - $t_2$
  - $0$
  - $0$
  - $0$
  - $1$
  - $0$
  - $0$
  - $-\tfrac12$
  - $\tfrac14$
  - $22$
* - $t_3$
  - $0$
  - $0$
  - $0$
  - $0$
  - $1$
  - $0$
  - $-\tfrac12$
  - $-\tfrac34$
  - $30$
* - $t_4$
  - $0$
  - $0$
  - $0$
  - $0$
  - $0$
  - $1$
  - $\tfrac12$
  - $-\tfrac54$
  - $36$
* - $-w$
  - $0$
  - $0$
  - $0$
  - $0$
  - $0$
  - $0$
  - $-5$
  - $-2$
  - $-179$
```

这两张表包含相同的信息，只是交叉排列。

- 原始表的目标行列出了对偶松弛量：食物 1–4 的检验数 $(32,22,30,36)$ 就是对偶表中 $t_1,\dots,t_4$ 的值，而食物 5 和 6 的检验数为 0，因为它们的价格上限是紧的。
- 原始目标行中剩余变量列 $x_7$、$x_8$ 的元素为 $\bar c=0-\vpi^\T(-\ve_i)=\pi_i$，即 $\vpi^*=(3,8)$（{numref}`tab-03-05-reading-duals`）。
- 对偶目标行中松弛变量列 $t_5$、$t_6$ 的元素为 $(-5,-2)=-(x_5^*,x_6^*)$：对偶问题的对偶问题是原问题（{prf:ref}`thm-03-05-dual-of-dual`），同样的读取规则给出最优配餐。$t_1,\dots,t_4$ 下以及 $\pi_1$、$\pi_2$ 下的零元素分别对应于未购买的食物和两个维生素约束的零剩余量。
- 两张表的右下角都是 $-179$。

*（此表的较早版本在目标行的 $x_4$ 列中写的是 6；食物 4 的检验数为 $50-(2\cdot3+1\cdot8)=36$，即对偶表中 $t_4$ 的值。）*

(sec-03-05-shadow-prices)=
## 影子价格与经济解释

### 对偶值作为变化率

**影子价格**是对已经最优地利用某种稀缺资源的决策者而言，多一个单位该资源（或少一个单位需求）的价值：即在边际上，右端项每增加一个单位所带来的最优目标函数值的改善量。对线性规划而言，只要最优基不变，影子价格恰好就是最优对偶值。

:::{prf:theorem} 最优对偶值是影子价格
:label: thm-03-05-shadow-price

设 $\mB$ 是 $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$ 的一个最优基，其乘子为 $\vpi^{*\T}=\vc_B^\T\mB^{-1}$。对每个满足 $\mB^{-1}(\vb+\Delta\vb)\ge\vzero$ 的改变量 $\Delta\vb$，基 $\mB$ 仍是最优的，并且

$$
z^*(\vb+\Delta\vb)=z^*(\vb)+\vpi^{*\T}\Delta\vb .
$$ (eq-03-05-shadow-price)

特别地，若 $\mB$ 非退化（$\mB^{-1}\vb>\vzero$），则上式对所有充分小的 $\Delta\vb$ 成立，并且对每个 $i$ 有 $\partial z^*/\partial b_i=\pi_i^*$。对极大化问题（$z^*$ 为其最大值），同样的结论成立。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

检验数 $\vc_N^\T-\vc_B^\T\mB^{-1}\mN$ 不依赖于 $\vb$，所以改变后 $\mB$ 仍是对偶可行的。若还有 $\mB^{-1}(\vb+\Delta\vb)\ge\vzero$，则 $\mB$ 对新的右端项是原始可行的，因而由{prf:ref}`prop-03-05-basis-optimal`是最优的，并且 $z^*(\vb+\Delta\vb)=\vc_B^\T\mB^{-1}(\vb+\Delta\vb)=z^*(\vb)+\vpi^{*\T}\Delta\vb$。若 $\mB^{-1}\vb>\vzero$，则条件 $\mB^{-1}(\vb+\Delta\vb)\ge\vzero$ 对 $\vzero$ 的某个邻域内的所有 $\Delta\vb$ 都成立，所以 $z^*$ 在该邻域内是线性的，梯度为 $\vpi^*$。
:::

只改变一个右端项 $\Delta\vb=\delta\ve_i$ 时，只要 $\mB^{-1}\vb+\delta\,\mB^{-1}\ve_i\ge\vzero$，基就保持最优；这是关于 $\delta$ 的一组线性不等式，它确定一个区间，即影子价格 $\pi_i^*$ 的**有效范围**（range of validity）。在这个区间之外，基会改变，价格也随之改变。作为 $b_i$ 的单变量函数，$z^*$ 是分段线性的，对极小化问题是凸的，对极大化问题是凹的；对于极大化问题的资源约束，这意味着每增加一个单位资源的价值不超过前一个单位的价值。在退化的最优解处，最优对偶解不一定唯一，$z^*$ 的单侧变化率可能不同；这时对极小化问题，每个最优的 $\vpi^*$ 都给出一个有效的界 $z^*(\vb+\Delta\vb)\ge z^*(\vb)+\vpi^{*\T}\Delta\vb$（由弱对偶性，因为 $\vpi^*$ 仍是对偶可行的）。关于变化范围（包括费用系数和新增约束）的系统讨论见{doc}`03-07-sensitivity-analysis`。

对产品组合问题（{prf:ref}`ex-03-05-product-mix-dual`），最优基为 $(x_1,x_2,x_5)$，其中 $x_5$ 是需求约束的松弛变量，且

$$
\begin{gathered}
\mB^{-1}=\begin{pmatrix}\tfrac27&-\tfrac1{14}&0\\[2pt]-\tfrac17&\tfrac3{35}&0\\[2pt]-\tfrac27&\tfrac1{14}&1\end{pmatrix},
\\
\mB^{-1}\vb=\begin{pmatrix}\tfrac{45}{7}\\[2pt]\tfrac{30}{7}\\[2pt]\tfrac{11}{7}\end{pmatrix}.
\end{gathered}
$$

增加 $\delta$ 小时机时，基本解变为 $\bigl(\tfrac{45}{7}+\tfrac27\delta,\ \tfrac{30}{7}-\tfrac17\delta,\ \tfrac{11}{7}-\tfrac27\delta\bigr)$，当 $-\tfrac{45}{2}\le\delta\le\tfrac{11}{2}$ 时它保持非负。因此，当每周生产时间在 37.5 小时到 65.5 小时之间时，每小时 $\pi_1^*=550/7\approx\$78.57$ 的影子价格都成立：多一个小时可使每周利润从 5142.86 美元增加到 5221.43 美元。超过 65.5 小时后，果汁杯 8 百箱的需求上限成为紧约束，基发生改变，额外的机时就毫无价值（最优解停留在 $\vx=(8,3.5)$，$z=\$5575$）。同样的计算给出仓储价格每立方米 $\pi_2^*=20/7\approx\$2.86$ 的有效范围为 128 到 240 立方米；而 $\pi_3^*=0$ 的需求上限可以降到 $45/7\approx6.43$ 百箱，才开始起作用。

### 资源的价值

对产品组合问题，强对偶性表述为

$$
\begin{aligned}
\underbrace{500x_1^*+450x_2^*}_{\text{贡献（美元/周）}}
&=\underbrace{60\,\pi_1^*}_{\text{机时}}+\underbrace{150\,\pi_2^*}_{\text{仓储}}+\underbrace{8\,\pi_3^*}_{\text{需求}}\\
&=\frac{33000}{7}+\frac{3000}{7}+0=\frac{36000}{7}.
\end{aligned}
$$

右边每一项都是一种资源的数量乘以其单位价格：60 小时，每小时 78.57 美元；150 立方米，每立方米 2.86 美元；8 百箱需求，价格 0 美元。最优利润完全由资源的价值来解释，而未用完的资源（这里是需求上限）价格为零，这就是互补松弛性。对任何可行的生产计划和任何可行的价格，弱对偶性表明利润不超过资源的价值，$\vc^\T\vx\le\vb^\T\vpi$，在最优处等号成立。

### 隐含成本与对偶约束

对偶约束也有同样直接的解释。对原问题 $\max\{\vc^\T\vx:\mA\vx\le\vb,\ \vx\ge\vzero\}$，量

$$
\sum_{i=1}^{m}a_{ij}\,\pi_i=\vpi^\T\mA_{\cdot j}
$$ (eq-03-05-imputed-cost)

是活动 $j$ 一个单位的**隐含成本**（imputed cost）：它所用的资源（资源 $i$ 用 $a_{ij}$ 个单位）按价格 $\pi_i$ 计算的价值。对偶可行性 $\vpi^\T\mA_{\cdot j}\ge c_j$ 表明，没有任何活动的收益超过它所消耗的资源的价值，而检验数 $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ 就是单位收益减去隐含成本。在一个尚非最优的基的乘子 $\vpi^\T=\vc_B^\T\mB^{-1}$ 下，某个活动有 $c_j>\vpi^\T\mA_{\cdot j}$：它的单位收益超过了所需资源的隐含成本，因此增大它是合算的，单纯形法就让它进基。在最优处，每个被采用的活动恰好收支平衡，没有任何活动在边际上有利可图。在产品组合问题中，一百箱果汁杯使用 6 小时和 10 立方米，隐含成本为 $6\cdot\frac{550}{7}+10\cdot\frac{20}{7}+1\cdot0=500$ 美元；一百箱鸡尾酒杯使用 5 小时和 20 立方米，隐含成本为 $5\cdot\frac{550}{7}+20\cdot\frac{20}{7}=450$ 美元：两者都等于各自的利润，对于产量为正的产品必然如此。

:::{prf:remark} 线性规划之外的影子价格
:label: rem-03-05-shadow-price-meaning

影子价格的思想远远超出了线性规划的范围。在公共项目的成本–收益分析中，许多影响没有市场价格：一条新的高速公路每年能为通勤者节省 500 000 小时、每年挽救五条生命并减少空气污染，其现值成本为 2.5 亿美元；只有在给这些时间和生命赋予货币价值之后，才能把它们与成本相比较，而这些价值就是时间和安全的影子价格。它们根据人们的支付意愿来估计，因而依赖于一些假设。在运营管理中，一条已经按每周 40 小时上限运行的生产线，其时间具有影子价格：管理者为多运行一小时最多应支付的金额，是这一小时所能带来的最优利润的增量，线性模型把它报告为时间约束的对偶值。在非线性规划中，KKT 条件中的拉格朗日乘子（Lagrange multiplier）起着同样的作用（{prf:ref}`thm-06-03-kkt`）；{prf:ref}`ex-08-03-nlp-multipliers`是一个例子。
:::

(sec-03-05-farkas)=
## Farkas 引理

对偶性还能回答一个完全不涉及目标函数的问题：线性方程组何时有非负解？**择一定理**（theorem of the alternative）断言两个方程组中恰有一个有解；于是第二个方程组的解就*证明*了第一个方程组无解。其中最著名的一个出自 Farkas（1902）。

:::{prf:theorem} Farkas 引理
:label: thm-03-05-farkas

设 $\mA\in\R^{m\times n}$，$\vb\in\R^m$。下列两个系统中恰有一个有解：

1. $\mA\vx=\vb$, $\vx\ge\vzero$；
2. $\mA^\T\vpi\le\vzero$, $\vb^\T\vpi>0$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

考虑线性规划问题 $\min\{\vzero^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$。由 {eq}`eq-03-05-standard-pair`，其对偶问题为 $\max\{\vb^\T\vpi:\mA^\T\vpi\le\vzero\}$，它是可行的（$\vpi=\vzero$）。

设系统 1 有解。则原问题可行且最优值为 0，由弱对偶性，每个对偶可行的 $\vpi$ 都有 $\vb^\T\vpi\le0$，所以系统 2 无解。

设系统 1 无解。则原问题不可行而对偶问题可行，所以由{prf:ref}`cor-03-05-feasible-pair`，对偶问题无界：存在 $\vpi$ 使 $\mA^\T\vpi\le\vzero$ 且 $\vb^\T\vpi>0$（事实上 $\vb^\T\vpi$ 可以任意大），系统 2 有解。
:::

从几何上看，系统 1 表示 $\vb$ 位于 $\mA$ 的各列所生成的锥 $\cone\{\mA_{\cdot1},\dots,\mA_{\cdot n}\}$ 中。系统 2 表示超平面 $\{\vy:\vpi^\T\vy=0\}$ 把 $\vb$ 与这个锥分离：所有列都位于一侧（$\vpi^\T\mA_{\cdot j}\le0$），而 $\vb$ 严格位于另一侧（$\vpi^\T\vb>0$）。因此 Farkas 引理是关于一个点与一个多面体锥的分离超平面定理（{prf:ref}`thm-02-01-separation-point`；有限生成的锥是闭的）。这里给出的证明依赖于强对偶性，而强对偶性是用单纯形法证明的；也可以从分离定理出发证明 Farkas 引理，然后反过来由 Farkas 引理推出强对偶性。通过通常的变换可以得到其他形式，例如：$\mA\vx\le\vb$, $\vx\ge\vzero$ 与 $\mA^\T\vpi\ge\vzero$, $\vpi\ge\vzero$, $\vb^\T\vpi<0$ 中恰有一个有解（对 $(\mA\ \ \mI)$ 应用该引理）。该引理是从线性规划通向非线性规划最优性条件的桥梁：在{doc}`06-03-kkt-conditions`中，它把“不存在可行的改善方向”转化为 KKT 乘子的存在性，正如 {eq}`eq-03-05-cone-condition` 在线性规划的顶点处所做的那样。

## 本章回顾

### 小结

- 每个线性规划问题都有对偶问题；对于对称形式的问题对 {eq}`eq-03-05-canonical-pair`，即 $\min\{\vc^\T\vx:\mA\vx\ge\vb,\vx\ge\vzero\}$ 与 $\max\{\vb^\T\vpi:\mA^\T\vpi\le\vc,\vpi\ge\vzero\}$，对偶问题就是{prf:ref}`ex-03-05-diet`中的定价问题。标准形问题对的对偶变量是自由的 {eq}`eq-03-05-standard-pair`。
- 任何线性规划问题的对偶问题都可以由{numref}`tab-03-05-dual-rules`得到；这些规则正是使拉格朗日函数 {eq}`eq-03-05-lagrangian` 给出下界的条件。等价问题的对偶问题等价（{prf:ref}`prop-03-05-equivalent-duals`），对偶问题的对偶问题是原问题（{prf:ref}`thm-03-05-dual-of-dual`）。
- 弱对偶性（{prf:ref}`thm-03-05-weak-duality`）：原始可行值与对偶可行值互为界；两者相等即证明了最优性（{prf:ref}`cor-03-05-optimality-certificate`）。
- 强对偶性（{prf:ref}`thm-03-05-strong-duality`）：若一个问题有最优解，则两者都有，且 $z^*=w^*$；最优对偶解为 $\vpi^{*\T}=\vc_B^\T\mB^{-1}$。问题对的可能结果列于{numref}`tab-03-05-status`。
- 互补松弛性（{prf:ref}`thm-03-05-complementary-slackness`）：可行的 $\vx$、$\vpi$ 是最优解，当且仅当 {eq}`eq-03-05-cs` 中每个变量与松弛量的乘积都为零；在最优顶点处，$\vc$ 位于有效约束法向量所张成的锥中 {eq}`eq-03-05-cone-condition`（{numref}`fig-03-05-cs-cone`）。
- 一个基是对偶可行的，当且仅当其检验数具有最优符号 {eq}`eq-03-05-dual-feasible`；既原始可行又对偶可行就意味着最优（{prf:ref}`prop-03-05-basis-optimal`）。最优对偶值可以从最终单纯形表中读出（{numref}`tab-03-05-reading-duals`）。
- 在有效范围内，最优对偶值是影子价格，$\partial z^*/\partial b_i=\pi_i^*$（{prf:ref}`thm-03-05-shadow-price`）；对偶约束把单位收益与隐含成本 {eq}`eq-03-05-imputed-cost` 相比较。
- Farkas 引理（{prf:ref}`thm-03-05-farkas`）可由对偶性导出：要么 $\mA\vx=\vb$ 有解 $\vx\ge\vzero$，要么存在向量 $\vpi$ 把 $\vb$ 与 $\mA$ 各列所生成的锥分离。

### 常见错误

- **对偶变量的符号错误。** 相反类型的约束（极小化问题中的 $\le$ 行）对应 $\pi_i\le0$，而不是 $\pi_i\ge0$；方程对应自由的 $\pi_i$。没有把握时，把该行乘以 $-1$ 再使用对称形式的问题对。
- **在弱对偶性中漏掉符号条件。** 在对称形式的问题对中，$\vpi^\T\mA\vx\ge\vpi^\T\vb$ 这一步需要 $\vpi\ge\vzero$；价格自由时，除非约束是方程，否则这一步不成立。
- **使用互补松弛性而不检验可行性。** 由互补松弛性得到的方程只给出候选解；只有当它还满足所在问题的全部约束时才是最优的（{prf:ref}`ex-03-05-diet-cs`）。
- **把互补松弛性误读为“恰有一个为零”。** 变量与其对应的松弛量可以都为零。
- **从单纯形表读对偶值时出现符号错误。** 应使用 $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$，代入松弛变量、剩余变量或人工变量的实际列和费用（{numref}`tab-03-05-reading-duals`），而不是死记某一种单纯形表约定下的规则。
- **在有效范围之外使用影子价格。** 只有当最优基保持原始可行时，$\pi_i^*$ 才能预测 $z^*$ 的变化；$b_i$ 的大幅改变通常使每单位的收益减少。
- **由不可行推出无界。** 原问题不可行并不意味着对偶问题无界；两者可能都不可行（{numref}`tab-03-05-status`）。

## 参考文献

- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Chapter 4: Sections 4.1–4.3 on the dual
  problem and the duality theorems, Section 4.4 on optimal dual variables as
  marginal costs, Section 4.6 on Farkas' lemma.)
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Chapter 6: the simplex algorithm and the duality theorem;
  Chapter 12: the price concept in linear programming.)
- Farkas, J. (1902). Theorie der einfachen Ungleichungen. *Journal für die
  reine und angewandte Mathematik, 124*, 1–27.
- Gale, D., Kuhn, H. W., & Tucker, A. W. (1951). Linear programming and the
  theory of games. In T. C. Koopmans (Ed.), *Activity analysis of production
  and allocation* (pp. 317–329). Wiley.
- Goldman, A. J., & Tucker, A. W. (1956). Theory of linear programming. In
  H. W. Kuhn & A. W. Tucker (Eds.), *Linear inequalities and related systems*
  (pp. 53–97). Princeton University Press.
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 6: Sections 6.1–6.4 on the
  essence of duality theory, its economic interpretation, primal–dual
  relationships and the dual of other primal forms.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 4: Sections 4.1–4.3 on the definition of the dual,
  primal–dual relationships and the economic interpretation of duality;
  Section 7.4: duality in matrix form.)
