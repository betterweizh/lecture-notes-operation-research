# 向量与矩阵微积分

## 概述

非线性最优化就是多元微积分。{doc}`06-02-optimality-conditions`中的一阶条件令梯度为零，{doc}`06-03-kkt-conditions`中的 Karush–Kuhn–Tucker（KKT）条件把目标函数和约束的梯度组合起来，二阶条件考察黑塞矩阵，而从{doc}`06-04-line-search`到{doc}`06-06-constrained-algorithms`的每一个算法都由泰勒展开建立目标函数的线性或二次模型。本章提供这些工具。本章确定全书的**布局约定**（layout convention）（梯度是列向量，雅可比矩阵的各行是转置后的梯度），定义偏导数、方向导数、梯度、雅可比矩阵和黑塞矩阵，并推导反复出现的导数：线性函数、二次函数、最小二乘目标函数和范数的导数（{numref}`tab-01-02-vector-derivatives`），以及迹、行列式和逆矩阵对矩阵的导数（{numref}`tab-01-02-matrix-derivatives`）。本章给出向量形式的链式法则（{prf:ref}`thm-01-02-chain-rule`）和带余项的泰勒定理（{prf:ref}`thm-01-02-taylor`），证明梯度与水平集正交并指向最速上升方向（{prf:ref}`thm-01-02-gradient-level-set`、{prf:ref}`thm-01-02-steepest-ascent`），最后给出三个应用：最小二乘的正规方程、二次函数的极小化以及拉格朗日函数的梯度。本章用到{doc}`01-01-linear-algebra`中的矩阵代数（内积、对称矩阵、二次型、范数与投影）以及一元微积分（包括中值定理）。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 按照本书的布局约定计算偏导数、梯度、雅可比矩阵、黑塞矩阵和方向导数；
- 对线性型和二次型、最小二乘目标函数以及范数求导，并使用向量形式的链式法则；
- 求迹、对数行列式和逆矩阵对矩阵的导数；
- 写出带余项的一阶和二阶泰勒展开，并判断它们对函数的逼近程度；
- 解释为什么梯度与水平集正交，以及为什么 $-\nabla f$ 是最速下降方向；
- 推导正规方程、二次函数的极小点以及拉格朗日函数的梯度。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\partial f/\partial x_j$ | $f$ 对 $x_j$ 的偏导数 |
| $\nabla f(\vx)$, $\nabla^2 f(\vx)$ | $f$ 在 $\vx$ 处的梯度（列向量）和黑塞矩阵 |
| $f'(\vx;\vd)$ | $f$ 在 $\vx$ 处沿方向 $\vd$ 的方向导数 |
| $D\vg(\vx)$ | $\vg:\R^n\to\R^m$ 在 $\vx$ 处的雅可比矩阵（$m\times n$） |
| $C^1$, $C^2$ | 连续可微，二阶连续可微 |
| $\bar\vx$ | 函数的展开点 |
| $L_c(f)$ | 水平集 $\{\vx:f(\vx)=c\}$ |
| $\mX$, $\partial f/\partial\mX$ | 矩阵变量，以及标量函数对它的导数（与 $\mX$ 形状相同） |
| $o(\norm{\vd}^k)$ | 当 $\vd\to\vzero$ 时比 $\norm{\vd}^k$ 更快趋于零的项 |
| $O(\norm{\vd}^k)$ | 当 $\vd$ 很小时以常数乘 $\norm{\vd}^k$ 为界的项 |
| $L(\vx,\vlambda,\vmu)$ | 非线性规划的拉格朗日函数 |

(sec-01-02-layout)=
## 多元函数与布局约定

本章研究**标量函数**（scalar function）$f:\R^n\to\R$，例如目标函数 $f(\vx)$；以及分量为 $g_1,\dots,g_m$ 的**向量值函数**（vector function）$\vg:\R^n\to\R^m$，例如非线性规划的约束函数。向量值函数的导数是数组，其排列方式是一种约定。不同教材的约定不同，混用约定是转置放错位置的最常见原因。本书始终采用下面的约定。

:::{prf:definition} 布局约定
:label: def-01-02-layout

1. 标量函数的**梯度**（gradient）$\nabla f(\vx)$ 是 $\R^n$ 中的**列**向量，其分量为 $\partial f/\partial x_j$。
2. 向量值函数 $\vg:\R^n\to\R^m$ 的**雅可比矩阵**（Jacobian matrix）$D\vg(\vx)$ 是元素为 $\partial g_i/\partial x_j$ 的 $m\times n$ 矩阵：第 $i$ 行是分量 $g_i$ 的梯度的转置，

   $$
   D\vg(\vx)=\begin{pmatrix}\nabla g_1(\vx)^\T\\\vdots\\\nabla g_m(\vx)^\T\end{pmatrix}.
   $$ (eq-01-02-jacobian-rows)

3. 对标量函数，$Df(\vx)=\nabla f(\vx)^\T$ 是一个行向量，而**黑塞矩阵**（Hessian matrix）是梯度的雅可比矩阵 $\nabla^2 f(\vx)=D(\nabla f)(\vx)$，是一个 $n\times n$ 矩阵。
:::

对雅可比矩阵而言，这是 Magnus and Neudecker (2019) 的“分子布局”，再配以列向量形式的梯度。其优点是一阶近似式与一元公式 $g(x+d)\approx g(x)+g'(x)d$ 的形式相同：$\vg(\vx+\vd)\approx\vg(\vx)+D\vg(\vx)\vd$，$f(\vx+\vd)\approx f(\vx)+\nabla f(\vx)^\T\vd$。有些教材（主要是统计学和机器学习方面的）使用转置的雅可比矩阵（$n\times m$）；它们的公式与本书的公式相差一个转置。

(sec-01-02-partial)=
## 偏导数、方向导数、梯度与黑塞矩阵

**偏导数**（partial derivative）$\partial f/\partial x_j(\vx)$ 是固定其他变量时 $f$ 对 $x_j$ 的普通导数，即 $f$ 沿单位向量 $\ve_j$ 的导数。偏导数只描述 $f$ 沿坐标轴的变化；一个函数可以在某点具有全部偏导数，却在该点不能被线性函数很好地近似。有用的概念是下面这个。

:::{prf:definition} 可微性与梯度
:label: def-01-02-gradient

若存在向量 $\nabla f(\vx)\in\R^n$（称为**梯度**），使得

$$
f(\vx+\vd)=f(\vx)+\nabla f(\vx)^\T\vd+o(\norm{\vd}) .
$$ (eq-01-02-differentiable)

则称 $f:\R^n\to\R$ 在 $\vx$ 处**可微**（differentiable）。取 $\vd=\alpha\ve_j$ 并令 $\alpha\to0$，可知梯度就是偏导数构成的向量 $\nabla f(\vx)=\bigl(\partial f/\partial x_1(\vx),\dots,\partial f/\partial x_n(\vx)\bigr)^\T$。若函数的偏导数在一个开集上存在且连续，则称它**连续可微**（continuously differentiable，$C^1$）；每个 $C^1$ 函数都是可微的。
:::

:::{prf:definition} 方向导数
:label: def-01-02-directional-derivative

$f$ 在 $\vx$ 处沿方向 $\vd$ 的**方向导数**（directional derivative）为

$$
f'(\vx;\vd)=\lim_{\alpha\downarrow0}\frac{f(\vx+\alpha\vd)-f(\vx)}{\alpha},
$$

即从 $\vx$ 出发沿 $\vd$ 移动时 $f$ 的初始变化率。若 $f$ 在 $\vx$ 处可微，则把 $\alpha\vd$ 代入 {eq}`eq-01-02-differentiable` 得

$$
f'(\vx;\vd)=\nabla f(\vx)^\T\vd .
$$ (eq-01-02-directional)
:::

公式 {eq}`eq-01-02-directional` 就是梯度进入最优化的方式：$\nabla f(\vx)^\T\vd$ 的符号说明了沿 $\vd$ 走一小步时 $f$ 是增加还是减少。满足 $\nabla f(\vx)^\T\vd<0$ 的方向称为**下降方向**（descent direction）（{prf:ref}`def-06-02-descent-direction`）。对某些不可微函数也可以定义方向导数，例如 $\abs{x}$ 在 $0$ 处，此时 $f'(0;d)=\abs{d}$ 关于 $d$ 不是线性的。

:::{prf:definition} 黑塞矩阵
:label: def-01-02-hessian

若 $f$ 在 $\vx$ 处二阶可微，则它的**黑塞矩阵**是由二阶偏导数构成的 $n\times n$ 矩阵：

$$
\nabla^2 f(\vx)=\left(\frac{\partial^2 f}{\partial x_i\,\partial x_j}(\vx)\right)_{i,j=1}^n .
$$

若二阶偏导数连续（$f\in C^2$），则求导次序无关紧要（施瓦茨定理），因此黑塞矩阵是**对称的**。
:::

正因为对称，{doc}`01-01-linear-algebra`中关于对称矩阵的理论（实特征值、谱定理、有定性判别法）才适用于黑塞矩阵。

:::{prf:example} 二次函数的梯度与黑塞矩阵
:label: ex-01-02-gradient-hessian-quadratic

设 $f(\vx)=2x_1^2+3x_1x_2+x_2^2$，定义在 $\R^2$ 上。

求它的梯度和黑塞矩阵，它在 $\bar\vx=(1,1)$ 处沿单位方向 $(1,-1)/\sqrt2$ 和 $(-3,4)/5$ 的方向导数，以及它的驻点的类型。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**梯度。** 固定 $x_2$，$\partial f/\partial x_1=4x_1+3x_2$；固定 $x_1$，$\partial f/\partial x_2=3x_1+2x_2$。因此

$$
\begin{aligned}
\nabla f(\vx)&=\begin{pmatrix}4x_1+3x_2\\3x_1+2x_2\end{pmatrix},\\
\nabla^2 f(\vx)&=\begin{pmatrix}4&3\\3&2\end{pmatrix}.
\end{aligned}
$$

黑塞矩阵是常数矩阵（每个二次函数都是如此），且是对称的：$\partial^2f/\partial x_1\partial x_2=\partial^2f/\partial x_2\partial x_1=3$。

**在 $\bar\vx=(1,1)$ 处的方向导数。** 这里 $\nabla f(\bar\vx)=(7,5)$。沿单位方向 $\vd=(1,-1)/\sqrt2$，$f'(\bar\vx;\vd)=(7-5)/\sqrt2=\sqrt2\approx1.414$，因此 $f$ 增加。沿单位方向 $\vd=(-3,4)/5$，$f'(\bar\vx;\vd)=(-21+20)/5=-0.2<0$，因此 $\vd$ 是下降方向，尽管它与 $-\nabla f$ 相差很远。在所有单位方向中，最大变化率为 $\norm{(7,5)}=\sqrt{74}\approx8.602$（{prf:ref}`thm-01-02-steepest-ascent`）。

**与二次型的联系。** $f(\vx)=\vx^\T\mA\vx$，其中对称矩阵 $\mA=\begin{pmatrix}2&1.5\\1.5&1\end{pmatrix}$，确实有 $\nabla f=2\mA\vx$ 和 $\nabla^2f=2\mA$（{numref}`tab-01-02-vector-derivatives`）。黑塞矩阵的行列式为 $8-9=-1<0$，因此它是不定的；其特征值为 $3\pm\sqrt{10}$，约为 $-0.1623$ 和 $6.162$。唯一的驻点是 $\vx=\vzero$（方程组 $4x_1+3x_2=0$，$3x_1+2x_2=0$ 的行列式为 $-1\ne0$），它是一个鞍点：$f$ 可分解为 $(2x_1+x_2)(x_1+x_2)$，在过原点的某些直线上为正，在另一些直线上为负（比较{numref}`fig-01-01-quadratic-forms`(c)）。
:::

(sec-01-02-jacobian)=
## 雅可比矩阵与导数表

对在 $\vx$ 处可微的向量值函数 $\vg:\R^n\to\R^m$，雅可比矩阵 {eq}`eq-01-02-jacobian-rows` 给出一阶近似

$$
\vg(\vx+\vd)=\vg(\vx)+D\vg(\vx)\,\vd+o(\norm{\vd}) .
$$ (eq-01-02-jacobian-approx)

书中反复出现两个例子。**仿射映射**（affine map）$\vg(\vx)=\mA\vx+\vb$ 在每一点都有 $D\vg(\vx)=\mA$：线性约束的梯度是常数，即 $\mA$ 的各行。**非线性约束向量**，例如 $\vg(\vx)=(x_1^2+x_2^2-1,\ x_1x_2)$，有

$$
D\vg(\vx)=\begin{pmatrix}2x_1&2x_2\\x_2&x_1\end{pmatrix},
$$

其各行为 $\nabla g_1(\vx)^\T$ 和 $\nabla g_2(\vx)^\T$；在 $\vx=(1,0)$ 处它等于 $\begin{pmatrix}2&0\\0&1\end{pmatrix}$，各行线性无关。有效约束的雅可比矩阵各行线性无关，就是{prf:ref}`thm-06-03-cq-sufficient`中的线性无关约束规范（linear independence constraint qualification, LICQ）。

{numref}`tab-01-02-vector-derivatives`中的导数在全书中都会用到。每一个都可以由定义用几行推导出来。

```{list-table} $\vx\in\R^n$ 的常见标量函数的梯度和黑塞矩阵（$\va\in\R^n$，$\mA\in\R^{n\times n}$，$\mM\in\R^{m\times n}$，$\vb\in\R^m$，$\vc\in\R^n$，$\mQ$ 对称）。
:name: tab-01-02-vector-derivatives
:header-rows: 1

* - $f(\vx)$
  - $\nabla f(\vx)$
  - $\nabla^2 f(\vx)$
* - $\va^\T\vx$
  - $\va$
  - $\vzero$
* - $\vx^\T\mA\vx$
  - $(\mA+\mA^\T)\vx$
  - $\mA+\mA^\T$
* - $\vx^\T\mQ\vx$（$\mQ$ 对称）
  - $2\mQ\vx$
  - $2\mQ$
* - $\tfrac12\vx^\T\mQ\vx+\vc^\T\vx$
  - $\mQ\vx+\vc$
  - $\mQ$
* - $\norm{\mM\vx-\vb}_2^2$
  - $2\mM^\T(\mM\vx-\vb)$
  - $2\mM^\T\mM$
* - $\norm{\vx}_2$（$\vx\ne\vzero$）
  - $\vx/\norm{\vx}_2$
  - $\bigl(\mI-\vx\vx^\T/\norm{\vx}_2^2\bigr)/\norm{\vx}_2$
* - $\vu(\vx)^\T\vv(\vx)$
  - $D\vu(\vx)^\T\vv(\vx)+D\vv(\vx)^\T\vu(\vx)$
  - —
* - $f(\mM\vy+\vb)$，作为 $\vy$ 的函数
  - $\mM^\T\nabla f(\mM\vy+\vb)$
  - $\mM^\T\nabla^2 f(\mM\vy+\vb)\mM$
```

**线性函数。** $\va^\T\vx=\sum_ja_jx_j$，因此 $\partial/\partial x_k=a_k$。

**二次型。** $\vx^\T\mA\vx=\sum_i\sum_ja_{ij}x_ix_j$。变量 $x_k$ 出现在 $i=k$ 的项和 $j=k$ 的项中，因此

$$
\begin{aligned}
\frac{\partial}{\partial x_k}\vx^\T\mA\vx&=\sum_ja_{kj}x_j+\sum_ia_{ik}x_i\\
&=(\mA\vx)_k+(\mA^\T\vx)_k ,
\end{aligned}
$$

由此得 $\nabla(\vx^\T\mA\vx)=(\mA+\mA^\T)\vx$。对线性函数 $(\mA+\mA^\T)\vx$ 再求一次导，得到常数黑塞矩阵 $\mA+\mA^\T$。对对称矩阵 $\mQ$，它们变为 $2\mQ\vx$ 和 $2\mQ$；$\tfrac12\vx^\T\mQ\vx$ 中的因子 $\tfrac12$ 正是为了消去这个 2。在{prf:ref}`ex-01-02-gradient-hessian-quadratic`中，非对称矩阵 $\begin{pmatrix}2&3\\0&1\end{pmatrix}$ 表示同一个 $f$，而 $(\mA+\mA^\T)\vx$ 给出相同的梯度。

**最小二乘。** 展开得 $\norm{\mM\vx-\vb}_2^2=\vx^\T\mM^\T\mM\vx-2\vb^\T\mM\vx+\vb^\T\vb$，这是一个以对称矩阵 $\mM^\T\mM$ 为矩阵的二次型加上一个 $\va=\mM^\T\vb$ 的线性项。由表中前面两行得到梯度 $2\mM^\T\mM\vx-2\mM^\T\vb$ 和黑塞矩阵 $2\mM^\T\mM\succeq0$。

**欧氏范数。** $\norm{\vx}_2=(\vx^\T\vx)^{1/2}$ 是 $s=\vx^\T\vx$ 的函数，由链式法则（{ref}`sec-01-02-chain-rule`）得 $\tfrac12 s^{-1/2}\cdot2\vx=\vx/\norm{\vx}_2$：一个背离原点的单位向量。范数在 $\vzero$ 处不可微（它在那里的表现与 $\abs{x}$ 类似），这就是最优化模型尽可能使用光滑的 $\norm{\vx}_2^2$ 的原因。它的黑塞矩阵满足 $\nabla^2\norm{\vx}_2\,\vx=\vzero$：沿从原点出发的射线，范数线性增长，没有曲率。

**乘积法则。** 对 $\vu,\vv:\R^n\to\R^m$，内积 $\vu^\T\vv=\sum_iu_iv_i$ 的梯度为 $\sum_i(v_i\nabla u_i+u_i\nabla v_i)$，由 {eq}`eq-01-02-jacobian-rows`，这就是表中相应的一项。

(sec-01-02-chain-rule)=
## 链式法则

:::{prf:theorem} 链式法则
:label: thm-01-02-chain-rule

设 $\vu:\R^n\to\R^m$ 在 $\vx$ 处可微，$\vv:\R^m\to\R^p$ 在 $\vu(\vx)$ 处可微。则 $\vv\circ\vu$ 在 $\vx$ 处可微，且

$$
D(\vv\circ\vu)(\vx)=D\vv\bigl(\vu(\vx)\bigr)\,D\vu(\vx).
$$ (eq-01-02-chain-rule)

对标量外函数 $f:\R^m\to\R$，上式成为 $\nabla(f\circ\vu)(\vx)=D\vu(\vx)^\T\,\nabla f\bigl(\vu(\vx)\bigr)$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

对 $\vu$ 应用 {eq}`eq-01-02-jacobian-approx`，得 $\vu(\vx+\vd)=\vu(\vx)+\ve$，其中 $\ve=D\vu(\vx)\vd+o(\norm{\vd})$，并且对很小的 $\vd$ 和某个常数 $c$ 有 $\norm{\ve}\le c\norm{\vd}$。对 $\vv$ 在 $\vu(\vx)$ 处以增量 $\ve$ 应用 {eq}`eq-01-02-jacobian-approx`，得 $\vv(\vu(\vx+\vd))=\vv(\vu(\vx))+D\vv(\vu(\vx))\ve+o(\norm{\ve}) =\vv(\vu(\vx))+D\vv(\vu(\vx))D\vu(\vx)\vd+o(\norm{\vd})$。标量形式是行向量 $Df(\vu(\vx))D\vu(\vx)$ 的转置。
:::

雅可比矩阵相乘的次序与函数复合的次序相同，与一元情形完全一样；正是{prf:ref}`def-01-02-layout`的布局约定使这一点成立。下面三种特殊情形会反复用到。

1. **仿射变量替换。** 对 $\vu(\vy)=\mM\vy+\vb$，$D\vu=\mM$，因此 $\nabla_{\vy}f(\mM\vy+\vb)=\mM^\T\nabla f(\mM\vy+\vb)$；再求一次导得到黑塞矩阵 $\mM^\T\nabla^2 f\,\mM$（{numref}`tab-01-02-vector-derivatives`的最后一行）。
2. **在直线上的限制。** 对固定的 $\vx$ 和 $\vd$，一元函数 $\varphi(\alpha)=f(\vx+\alpha\vd)$ 满足

   $$
   \begin{aligned}
   \varphi'(\alpha)&=\nabla f(\vx+\alpha\vd)^\T\vd,\\
   \varphi''(\alpha)&=\vd^\T\nabla^2 f(\vx+\alpha\vd)\,\vd .
   \end{aligned}
   $$ (eq-01-02-line-restriction)

   线搜索（{doc}`06-04-line-search`、{eq}`eq-06-04-phi`）极小化的正是这个 $\varphi$，而 {eq}`eq-01-02-line-restriction` 把它的一维导数与 $\nabla f$ 和 $\nabla^2 f$ 联系起来。
3. **与光滑标量函数的复合。** 对“log-sum-exp”函数 $f(\vx)=\ln\sum_{i=1}^m\exp(\va_i^\T\vx+b_i)$，写 $f=h\circ\vu$，其中 $\vu(\vx)=\mM\vx+\vb$（$\mM$ 的各行为 $\va_i^\T$），$h(\vy)=\ln\sum_i e^{y_i}$，其梯度的分量为 $e^{y_i}/\sum_ke^{y_k}$。由链式法则得 $\nabla f(\vx)=\mM^\T\vp$，其中在 $\vy=\mM\vx+\vb$ 处 $p_i=e^{y_i}/\sum_ke^{y_k}$：这是 $\va_i$ 的加权组合，权重为正且和为一。

(sec-01-02-matrix-derivatives)=
## 对矩阵的导数

最优化和统计学中的某些函数以矩阵 $\mX\in\R^{m\times n}$ 为自变量：协方差矩阵的对数行列式、加权乘积的迹、残差的 Frobenius 范数。标量函数 $f(\mX)$ 的导数逐个元素定义为 $\bigl(\partial f/\partial\mX\bigr)_{ij}=\partial f/\partial x_{ij}$，它与 $\mX$ 形状相同。利用矩阵内积 $\langle\mG,\mE\rangle=\tr(\mG^\T\mE)=\sum_{i,j} g_{ij}e_{ij}$，{eq}`eq-01-02-differentiable` 的类比是

$$
f(\mX+\mE)=f(\mX)+\tr\Bigl(\bigl(\tfrac{\partial f}{\partial\mX}\bigr)^\T\mE\Bigr)+o(\norm{\mE}) .
$$ (eq-01-02-matrix-gradient)

这给出了一个实用的方法：展开 $f(\mX+\mE)$，保留关于 $\mE$ 的线性项，利用循环法则 $\tr(\mP\mQ)=\tr(\mQ\mP)$ 把它们化为 $\tr(\mG^\T\mE)$ 的形式，即可读出 $\partial f/\partial\mX=\mG$。{numref}`tab-01-02-matrix-derivatives`列出了最常用的结果（Petersen & Pedersen, 2012, Section 2; Magnus & Neudecker, 2019, Chapters 8–9）。

```{list-table} 矩阵 $\mX$ 的标量函数的导数（各元素视为相互独立；$\mA$、$\va$、$\vb$ 为常数且维数相容；最后三行中 $\mX$ 为可逆方阵，对数一行中还要求 $\det\mX>0$）。
:name: tab-01-02-matrix-derivatives
:header-rows: 1

* - $f(\mX)$
  - $\partial f/\partial\mX$
* - $\tr(\mA\mX)$
  - $\mA^\T$
* - $\tr(\mX^\T\mA)$
  - $\mA$
* - $\va^\T\mX\vb$
  - $\va\vb^\T$
* - $\tr(\mX^\T\mA\mX)$
  - $(\mA+\mA^\T)\mX$
* - $\norm{\mX}_F^2=\tr(\mX^\T\mX)$
  - $2\mX$
* - $\ln\det\mX$
  - $\mX^{-\T}$
* - $\det\mX$
  - $(\det\mX)\,\mX^{-\T}$
* - $\va^\T\mX^{-1}\vb$
  - $-\mX^{-\T}\va\vb^\T\mX^{-\T}$
```

这里推导其中两项。

**乘积的迹。** $\tr(\mA\mX)=\sum_i\sum_ka_{ik}x_{ki}$。元素 $x_{kl}$ 只在 $i=l$ 的项中出现一次，系数为 $a_{lk}$。因此 $\partial\tr(\mA\mX)/\partial x_{kl}=a_{lk}$，即 $\partial\tr(\mA\mX)/\partial\mX=\mA^\T$。用上述方法也得到同样的结果，因为 $\tr(\mA(\mX+\mE)) =\tr(\mA\mX)+\tr\bigl((\mA^\T)^\T\mE\bigr)$。由于 $\va^\T\mX\vb=\tr(\vb\va^\T\mX)$，表中第三行是 $\mA=\vb\va^\T$ 的特例。

**对数行列式。** 对可逆矩阵 $\mX$ 和很小的 $\mE$，

$$
\begin{aligned}
\det(\mX+\mE)&=\det\mX\cdot\det(\mI+\mX^{-1}\mE)\\
&=\det\mX\cdot\bigl(1+\tr(\mX^{-1}\mE)+O(\norm{\mE}^2)\bigr),
\end{aligned}
$$

这是因为在 $\det(\mI+\mF)$ 的展开式中，对角元素之积贡献 $1+\sum_if_{ii}+O(\norm{\mF}^2)$，而其他每一项都至少含有 $\mF$ 的两个非对角元素。取对数，由 $\ln(1+t)=t+O(t^2)$ 得 $\ln\det(\mX+\mE)=\ln\det\mX+\tr(\mX^{-1}\mE) +O(\norm{\mE}^2)$。写出 $\tr(\mX^{-1}\mE)=\tr\bigl((\mX^{-\T})^\T\mE\bigr)$ 并与 {eq}`eq-01-02-matrix-gradient` 比较，得 $\partial\ln\det\mX/\partial\mX=\mX^{-\T}$；对对称矩阵 $\mX$，它就是 $\mX^{-1}$。再由链式法则即得下一行中 $\det\mX=\exp(\ln\det\mX)$ 的导数。

最后一行来自逆矩阵的一阶展开 $(\mX+\mE)^{-1}=\mX^{-1}-\mX^{-1}\mE\mX^{-1}+O(\norm{\mE}^2)$，乘以 $\mX+\mE$ 即可验证：乘积为 $\mI-\mE\mX^{-1}\mE\mX^{-1}$。

(sec-01-02-taylor)=
## 泰勒定理

泰勒展开在点 $\bar\vx$ 附近用多项式代替函数：关于步长 $\vd=\vx-\bar\vx$ 的一阶（线性）模型 $f(\bar\vx)+\nabla f(\bar\vx)^\T\vd$ 和二阶（二次）模型 $f(\bar\vx)+\nabla f(\bar\vx)^\T\vd+\tfrac12\vd^\T\nabla^2f(\bar\vx)\vd$。泰勒定理（Taylor's theorem）说明误差有多大。

:::{prf:theorem} 泰勒定理
:label: thm-01-02-taylor

设 $f$ 定义在一个包含从 $\bar\vx$ 到 $\bar\vx+\vd$ 的线段的开集上。

1. 若 $f\in C^1$，则存在 $t\in(0,1)$，使得

   $$
   f(\bar\vx+\vd)=f(\bar\vx)+\nabla f(\bar\vx+t\vd)^\T\vd ,
   $$ (eq-01-02-taylor-1)

   且 $f(\bar\vx+\vd)=f(\bar\vx)+\nabla f(\bar\vx)^\T\vd+o(\norm{\vd})$。
2. 若 $f\in C^2$，则存在 $t\in(0,1)$，使得

   $$
   f(\bar\vx+\vd)=f(\bar\vx)+\nabla f(\bar\vx)^\T\vd+\tfrac12\,\vd^\T\nabla^2 f(\bar\vx+t\vd)\,\vd ,
   $$ (eq-01-02-taylor-2)

   且 $f(\bar\vx+\vd)=f(\bar\vx)+\nabla f(\bar\vx)^\T\vd +\tfrac12\vd^\T\nabla^2f(\bar\vx)\vd+o(\norm{\vd}^2)$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

令 $\varphi(s)=f(\bar\vx+s\vd)$，$s\in[0,1]$。由 {eq}`eq-01-02-line-restriction`，$\varphi'(s)=\nabla f(\bar\vx+s\vd)^\T\vd$；当 $f\in C^2$ 时，$\varphi''(s)=\vd^\T\nabla^2f(\bar\vx+s\vd)\vd$。由一元中值定理，存在 $t\in(0,1)$ 使 $\varphi(1)=\varphi(0)+\varphi'(t)$，这就是 {eq}`eq-01-02-taylor-1`；由带拉格朗日余项的泰勒公式得 $\varphi(1)=\varphi(0)+\varphi'(0)+\tfrac12 \varphi''(t)$，这就是 {eq}`eq-01-02-taylor-2`。对于 $o$ 形式，把中间点换成 $\bar\vx$；误差分别为 $\bigl(\nabla f(\bar\vx+t\vd)-\nabla f(\bar\vx)\bigr)^\T\vd$ 和 $\tfrac12\vd^\T\bigl(\nabla^2f(\bar\vx+t\vd)-\nabla^2f(\bar\vx)\bigr)\vd$，由柯西–施瓦茨不等式以及 $\nabla f$（$\nabla^2f$）的连续性，它们分别是 $o(\norm{\vd})$ 和 $o(\norm{\vd}^2)$。
:::

注意各项在何处取值：在精确形式 {eq}`eq-01-02-taylor-2` 中，梯度在展开点处取值，只有最高阶导数在未知的中间点 $\bar\vx+t\vd$ 处取值。精确形式是证明的主力：例如，若 $\nabla f(\bar\vx)=\vzero$ 且 $\nabla^2f$ 在整条线段上半正定，则由 {eq}`eq-01-02-taylor-2` 立即得到 $f(\bar\vx+\vd)\ge f(\bar\vx)$（{prf:ref}`thm-02-02-second-order`、{prf:ref}`thm-06-02-second-order-sufficient`）。$o$ 形式是算法的基础：最速下降法在一个球上极小化线性模型，牛顿法极小化二次模型（{prf:ref}`alg-06-05-steepest-descent`、{eq}`eq-06-05-newton-model`）。

:::{prf:example} 一元函数的线性模型与二次模型
:label: ex-01-02-taylor-models

设 $f(x)=e^x-2x$，$\bar x=0$。

求 $f$ 在 $\bar x$ 处的线性泰勒模型 $m_1$ 和二次泰勒模型 $m_2$，在 $x=0.5$ 和 $x=2$ 处把它们与 $f$ 作比较，并比较 $m_2$ 的极小点与 $f$ 的极小点。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

此时 $f(0)=1$，$f'(x)=e^x-2$，$f'(0)=-1$，$f''(x)=e^x$，$f''(0)=1$。两个模型为

$$
\begin{aligned}
m_1(x)&=1-x,\\
m_2(x)&=1-x+\tfrac12x^2 .
\end{aligned}
$$

在 $x=0.5$ 处：$f(0.5)=e^{0.5}-1\approx0.6487$，$m_1(0.5)=0.5$（误差 $0.1487$），$m_2(0.5)=0.625$（误差 $0.0237$）。这里二次模型的精度是线性模型的六倍。$m_2$ 的误差是三阶余项 $f'''(\xi)\,0.5^3/6=e^{\xi}/48$，其中 $\xi\in(0,0.5)$，因此它必定介于 $1/48\approx0.0208$ 和 $e^{0.5}/48\approx0.0343$ 之间，实际也确实如此。精确形式 {eq}`eq-01-02-taylor-2` 在 $t\approx0.3475$ 时成立（即 $f(0.5)=1-0.5+\tfrac12e^{0.5t}\,0.25$），{eq}`eq-01-02-taylor-1` 在 $t\approx0.5208$ 时成立。

{numref}`fig-01-02-taylor-models`表明，两个模型都只在 $\bar x$ 附近才准确。二次模型的极小点是 $x=1$，对真正的极小点 $\ln2\approx0.6931$（此处 $f\approx0.6137$）而言，这是一个合理的初始估计，而走到这一点的一步就是牛顿法的一步。远离展开点时模型失效：$f(2)\approx3.389$，而 $m_2(2)=1$。

```{figure} figures/01-02-taylor-models.svg
:name: fig-01-02-taylor-models
:width: 520px
:class: dark-light
:alt: 凸曲线 e 的 x 次方减 2x，及其在 0 处的切线和在 0 处的密切抛物线；抛物线在 x 等于 1 处取最小值，曲线在 ln 2 处取最小值

$f(x)=e^x-2x$ 在 $\bar x=0$ 处的一阶和二阶泰勒模型。切线（虚线）与函数值和斜率吻合，抛物线还与曲率吻合，并在更大范围内贴近 $f$。二次模型的极小点 $x=1$ 近似于 $f$ 的极小点 $\ln2$。
```
:::

(sec-01-02-gradient-geometry)=
## 梯度的几何意义

$f$ 在高度 $c$ 处的**水平集**（level set）为 $L_c(f)=\{\vx:f(\vx)=c\}$。在 $\R^2$ 中，水平集就是地图上的等高线；在 $\R^3$ 中，水平集是曲面。下面两个事实使梯度成为下降法的核心对象。

:::{prf:theorem} 梯度与水平集正交
:label: thm-01-02-gradient-level-set

设 $f\in C^1$，$\bar\vx\in L_c(f)$，并设 $\vgamma(s)$，$s\in(-\varepsilon, \varepsilon)$，是 $L_c(f)$ 中满足 $\vgamma(0)=\bar\vx$ 的可微曲线。则 $\nabla f(\bar\vx)^\T\vgamma'(0)=0$：梯度与水平集的每个切方向都正交。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

对所有 $s$ 有 $f(\vgamma(s))=c$。用链式法则求导得 $\nabla f(\vgamma(s))^\T\vgamma'(s)=0$；令 $s=0$ 即可。
:::

:::{prf:theorem} 最速上升与最速下降
:label: thm-01-02-steepest-ascent

设 $f$ 在 $\vx$ 处可微，且 $\nabla f(\vx)\ne\vzero$。

1. 在所有单位向量 $\vd$ 中，方向导数 $f'(\vx;\vd)=\nabla f(\vx)^\T\vd$ 在 $\vd=\nabla f(\vx)/\norm{\nabla f(\vx)}$ 处最大，其值为 $\norm{\nabla f(\vx)}$；在 $\vd=-\nabla f(\vx)/\norm{\nabla f(\vx)}$ 处最小，其值为 $-\norm{\nabla f(\vx)}$。
2. 若 $\nabla f(\vx)^\T\vd<0$，则对所有充分小的 $\alpha>0$ 有 $f(\vx+\alpha\vd)<f(\vx)$。特别地，$-\nabla f(\vx)$ 是下降方向，$\nabla f(\vx)$ 是上升方向。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

(1) 由柯西–施瓦茨不等式，当 $\norm{\vd}=1$ 时 $\abs{\nabla f(\vx)^\T\vd}\le\norm{\nabla f(\vx)}$，等号成立当且仅当 $\vd$ 与梯度平行。(2) 由 {eq}`eq-01-02-differentiable`，$f(\vx+\alpha\vd)-f(\vx)=\alpha\bigl(\nabla f(\vx)^\T\vd+o(\alpha)/\alpha\bigr)$，而当 $\alpha>0$ 很小时括号内为负。
:::

这两个定理相互契合：沿水平集移动时，$f$ 在一阶意义下不变，因此变化最快的方向必定与水平集垂直。有两点需要注意。第一，“最速”是一个局部的说法：沿 $-\nabla f$ 走一大步可能走过头而使 $f$ 增加，这就是下降法需要步长规则的原因（{doc}`06-04-line-search`）。第二，梯度垂直于 $\vx$ 空间中的水平集，而不是垂直于 $\R^{n+1}$ 中 $f$ 的图像。

:::{prf:example} 线性函数的等值线
:label: ex-01-02-linear-level-lines

设 $f(\vx)=x_1+x_2$。

描述它的梯度和水平集，验证{prf:ref}`thm-01-02-gradient-level-set`，并求从 $\bar\vx=(2,1)$ 出发走长度为 $0.5$ 的一步最多能使 $f$ 增加多少。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

它的梯度 $\nabla f=(1,1)$ 在每一点都相同。水平集是一族平行直线 $x_1+x_2=c$，方向为 $(1,-1)$，而 $(1,1)^\T(1,-1)=0$ 验证了{prf:ref}`thm-01-02-gradient-level-set`。在 $\bar\vx=(2,1)$ 处函数值为 $f=3$；从 $\bar\vx$ 出发沿 $\nabla f$ 移动会依次穿过水平 $3.5$、$4$、…，沿 $-\nabla f$ 移动会依次穿过 $2.5$、$2$、…（{numref}`fig-01-02-linear-level-lines`）。沿单位梯度方向 $(1,1)/\sqrt2$ 走长度为 $0.5$ 的一步，使 $f$ 增加 $0.5\sqrt2\approx0.7071$；同样长度的步不可能做得更好，因为最大变化率为 $\norm{\nabla f}=\sqrt2$。

```{figure} figures/01-02-linear-level-lines.svg
:name: fig-01-02-linear-level-lines
:width: 460px
:class: dark-light
:alt: x1 加 x2 在取值 0 到 3 时的一族平行直线等值线，直线 f 等于 3 上的点 (2, 1)，以及从该点出发的两个箭头：指向右上方的梯度 (1, 1) 和指向左下方的负梯度，二者都与这些直线垂直

$f(\vx)=x_1+x_2$ 的等值线 $f=c$（标注了 $c=0,\dots,3$）。在 $(2,1)$ 处，梯度 $(1,1)$ 指向更高的水平，$-\nabla f$ 指向更低的水平；二者都与等值线垂直。
```
:::

对线性函数而言，沿负梯度方向 $f$ 永远以同样的速率减小，这就是在无界可行集上线性目标函数可能无界的原因（{ref}`sec-03-02-graphical-cases`）。

:::{prf:example} 抛物面的等值线
:label: ex-01-02-paraboloid-contours

设 $f(\vx)=(x_1-2)^2+(x_2-2)^2$，即到 $(2,2)$ 的距离的平方。

描述它的图像和水平集，计算它的梯度和黑塞矩阵，并在 $\bar\vx=(1.4,1.4)$ 处验证梯度与水平集的正交性以及上升与下降性质。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

它的图像是一个在 $(2,2)$ 处取最小值 $0$ 的抛物面，其水平集 $L_c(f)$（$c>0$）是以 $(2,2)$ 为圆心、半径为 $\sqrt c$ 的圆（{numref}`fig-01-02-paraboloid-contours`）。梯度和黑塞矩阵为

$$
\begin{aligned}
\nabla f(\vx)&=\begin{pmatrix}2(x_1-2)\\2(x_2-2)\end{pmatrix},\\
\nabla^2f(\vx)&=2\mI .
\end{aligned}
$$

在 $\bar\vx=(1.4,1.4)$ 处：$f(\bar\vx)=0.36+0.36=0.72$，因此 $\bar\vx$ 位于半径为 $\sqrt{0.72}\approx0.8485$ 的圆上，且 $\nabla f(\bar\vx)=(-1.2,-1.2)$。

- **正交性。** 圆在 $\bar\vx$ 处的切线方向为 $(1,-1)$，而 $(-1.2,-1.2)^\T(1,-1)=0$。
- **上升。** 梯度背离圆心，指向更大的圆，变化率为 $\norm{\nabla f(\bar\vx)}=1.2\sqrt2\approx1.697$；沿切线方向的变化率为 0，沿 $\ve_1$ 的变化率为 $-1.2$。
- **下降。** $-\nabla f(\bar\vx)=(1.2,1.2)$ 直指极小点：$\bar\vx-\tfrac12\nabla f(\bar\vx)=(2,2)$。之所以如此，是因为等值线是圆（$\nabla^2f=2\mI$）；对于狭长的椭圆等值线，$-\nabla f$ 会偏离中心，这正是{doc}`06-05-unconstrained-algorithms`中研究的最速下降法的弱点（{prf:ref}`thm-06-05-sd-rate`）。

```{figure} figures/01-02-paraboloid-contours.svg
:name: fig-01-02-paraboloid-contours
:width: 720px
:class: dark-light
:alt: 左图：定义在 1 到 3 的正方形上的碗形曲面，曲面上画有一条圆形水平曲线，以及它在底面上的投影。右图：围绕 (2, 2) 的同心圆等值线，(1.4, 1.4) 处背离圆心的梯度、指向圆心的负梯度，以及一条点线切线

(a) $f(\vx)=(x_1-2)^2+(x_2-2)^2$ 的图像，曲面上画出了水平 $f=0.72$，其投影是 $\vx$ 平面中的一个圆。(b) 等值线 $f=0.1,0.25,0.5,0.72,1,1.5$。在 $\bar\vx=(1.4,1.4)$ 处，梯度 $(-1.2,-1.2)$（按其长度的三分之一画出）与点线切线垂直并指向更高的水平；$-\nabla f$ 指向极小点 $(2,2)$。
```
:::

(sec-01-02-applications)=
## 应用

### 最小二乘与正规方程

用一个关于参数 $\vx\in\R^n$ 为线性的模型去拟合 $m>n$ 个观测值，会得到一个超定方程组 $\mM\vx\approx\vb$，其中 $\mM$ 的第 $i$ 行是第 $i$ 个观测的回归变量，$b_i$ 是其测量值。**最小二乘**（least squares）估计极小化 $f(\vx)=\norm{\mM\vx-\vb}_2^2$。由{numref}`tab-01-02-vector-derivatives`，

$$
\nabla f(\vx)=2\mM^\T(\mM\vx-\vb)=\vzero
\quad\Longleftrightarrow\quad
\mM^\T\mM\,\vx=\mM^\T\vb ,
$$ (eq-01-02-normal-equations)

这就是{doc}`01-01-linear-algebra`中从几何角度得到的**正规方程** {eq}`eq-01-01-normal-equations`。黑塞矩阵 $2\mM^\T\mM$ 是半正定的，因此 $f$ 是凸函数，{eq}`eq-01-02-normal-equations` 的每个解都是全局极小点（{prf:ref}`thm-02-02-first-order`）；若 $\mM$ 列满秩，则 $\mM^\T\mM\succ0$，极小点唯一。残差 $\vr=\vb-\mM\vx^*$ 满足 $\mM^\T\vr=\vzero$：它与 $\mM$ 的每一列都正交。

:::{prf:example} 用直线拟合交通流量计数
:label: ex-01-02-least-squares

一台交通计数器记录了某条道路在早高峰开始后 $t=0,1,2,3$ 小时的流量 $y$（单位：千辆/小时）：$y=1,3,4,4$。用最小二乘法拟合直线 $y\approx x_1+x_2t$。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

数据矩阵由一列 1 和一列时间组成，$\vb=(1,3,4,4)$：

$$
\begin{aligned}
\mM&=\begin{pmatrix}1&0\\1&1\\1&2\\1&3\end{pmatrix},\\
\mM^\T\mM&=\begin{pmatrix}4&6\\6&14\end{pmatrix},\\
\mM^\T\vb&=\begin{pmatrix}12\\23\end{pmatrix}.
\end{aligned}
$$

正规方程 $4x_1+6x_2=12$ 和 $6x_1+14x_2=23$（行列式为 $56-36=20\ne0$）给出 $x_2=1$，$x_1=1.5$。拟合的直线为 $y=1.5+1.0\,t$：截距为每小时 1500 辆，并且每过一小时流量增加每小时 1000 辆。残差 $\vr=\vb-\mM\vx^* =(-0.5,0.5,0.5,-0.5)$ 之和为零，且与时间列正交（$0\cdot(-0.5)+1\cdot0.5+2\cdot0.5+3\cdot(-0.5)=0$），正如 $\mM^\T\vr=\vzero$ 所要求的；最小平方和为 $1.0$。外推到 $t=4$ 小时得到 $5.5$ 千辆/小时，但应谨慎看待这个结果：最后两个计数表明流量已趋于平稳，而直线无法表示这一点。
:::

### 二次函数的极小化

设 $f(\vx)=\tfrac12\vx^\T\mQ\vx+\vc^\T\vx$，其中 $\mQ$ 对称。则 $\nabla f(\vx)=\mQ\vx+\vc$，$\nabla^2f=\mQ$，$f$ 的形状由 $\mQ$ 的有定性决定（{prf:ref}`def-01-01-definiteness`）。

- **$\mQ\succ0$。** 唯一的驻点是 $\vx^*=-\mQ^{-1}\vc$，配方得

  $$
  \begin{aligned}
  f(\vx)&=f(\vx^*)+\tfrac12(\vx-\vx^*)^\T\mQ(\vx-\vx^*),\\
  f(\vx^*)&=-\tfrac12\vc^\T\mQ^{-1}\vc ,
  \end{aligned}
  $$ (eq-01-02-quadratic-minimum)

  因此 $\vx^*$ 是唯一的全局极小点。对 $\mQ=\begin{pmatrix}2&1\\1&2\end{pmatrix}$ 和 $\vc=(-3,0)$，求解 $\mQ\vx=(3,0)$ 得 $\vx^*=(2,-1)$，$f(\vx^*)=\tfrac12\vc^\T\vx^*=-3$。
- **$\mQ\succeq0$ 且奇异。** 极小点存在当且仅当 $\mQ\vx=-\vc$ 有解，即 $\vc\in\col(\mQ)$；这时全部极小点构成一个仿射集。对 $\mQ=\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$ 和 $\vc=(1,-1)$，每个 $\vx=(-1,0)+s(1,1)$ 都是极小点，最优值为 $-0.5$。对与 $\col(\mQ)$ 正交的 $\vc=(1,1)$，$f(s,s)=2s$，因此 $f$ 无下界（例如 $f(-10,-10)=-20$）。
- **$\mQ$ 不定。** 沿负特征值 $\lambda_i(\mQ)<0$ 的特征向量 $\vv$，当 $s\to\pm\infty$ 时 $f(s\vv)=\tfrac12\lambda_i(\mQ)s^2\norm{\vv}^2+s\,\vc^\T\vv\to-\infty$；驻点（若存在）是鞍点。

二次规划（{ref}`sec-06-01-classes`）和牛顿法（{ref}`sec-06-05-newton`）都依赖于这三种情形。

### 拉格朗日函数的梯度

对非线性规划 $\min f(\vx)$，满足 $g_i(\vx)\le0$（$i=1,\dots,m$）和 $h_j(\vx)=0$（$j=1,\dots,\ell$），其**拉格朗日函数**（Lagrangian）为

$$
\begin{aligned}
L(\vx,\vlambda,\vmu)&=f(\vx)+\sum_{i=1}^m\lambda_ig_i(\vx)+\sum_{j=1}^\ell\mu_jh_j(\vx)\\
&=f(\vx)+\vlambda^\T\vg(\vx)+\vmu^\T\vh(\vx).
\end{aligned}
$$

它关于 $\vx$ 的梯度可由导数表和 {eq}`eq-01-02-jacobian-rows` 得到：

$$
\begin{aligned}
\nabla_{\vx}L(\vx,\vlambda,\vmu)&=\nabla f(\vx)+\sum_i\lambda_i\nabla g_i(\vx)+\sum_j\mu_j\nabla h_j(\vx)\\
&=\nabla f(\vx)+D\vg(\vx)^\T\vlambda+D\vh(\vx)^\T\vmu .
\end{aligned}
$$ (eq-01-02-lagrangian-gradient)

令其为零，就是 KKT 条件中的平稳性条件（{eq}`eq-02-03-kkt`、{prf:ref}`thm-06-03-kkt`）。对线性约束 $\vg(\vx)=\mA\vx-\vb$，雅可比矩阵为 $\mA$，条件成为 $\nabla f(\vx)+\mA^\T\vlambda=\vzero$。

:::{prf:example} 用拉格朗日函数求平面上的最近点
:label: ex-01-02-min-norm-lagrangian

求平面 $x_1+2x_2+2x_3=9$ 上距原点最近的点，即在约束 $h(\vx)=\va^\T\vx-9=0$ 下极小化 $f(\vx)=\tfrac12\norm{\vx}_2^2$，其中 $\va=(1,2,2)$。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

拉格朗日函数为 $L(\vx,\mu)=\tfrac12\vx^\T\vx+\mu(\va^\T\vx-9)$，由 {eq}`eq-01-02-lagrangian-gradient` 得 $\nabla_{\vx}L=\vx+\mu\va$。令其为零得 $\vx=-\mu\va$；代入约束得 $-\mu\,\va^\T\va=-9\mu=9$，因此 $\mu=-1$，$\vx^*=\va=(1,2,2)$。验证：$1+4+4=9$。最优值为 $f(\vx^*)=4.5$，原点到平面的距离为 $\norm{\vx^*}=3=9/\norm{\va}$。

$\vx^*$ 是全局极小点这一点可以直接看出。对任意可行点 $\vx$，$\va^\T(\vx-\vx^*)=0$，因此 $\vx^{*\T}(\vx-\vx^*)=0$，由勾股定理得 $\norm{\vx}^2=\norm{\vx^*}^2+\norm{\vx-\vx^*}^2\ge\norm{\vx^*}^2$。例如，可行点 $(9,0,0)$ 的函数值为 $f=40.5$。
:::

对含多个等式约束 $\mA\vx=\vb$（$\mA$ 行满秩）的情形，做与{prf:ref}`ex-01-02-min-norm-lagrangian`的解答中同样的计算，即得最小范数解 $\vx^*=\mA^\T(\mA\mA^\T)^{-1}\vb$。

## 本章回顾

### 小结

- 梯度是列向量，雅可比矩阵以转置后的梯度为行，黑塞矩阵是梯度的雅可比矩阵（{prf:ref}`def-01-02-layout`、{eq}`eq-01-02-jacobian-rows`）。
- 可微意味着有良好的线性近似 {eq}`eq-01-02-differentiable`；此时 $f'(\vx;\vd)=\nabla f(\vx)^\T\vd$ {eq}`eq-01-02-directional`，且 $C^2$ 函数的黑塞矩阵是对称的（{prf:ref}`def-01-02-hessian`）。
- 标准的梯度与黑塞矩阵：$\va^\T\vx$、$\vx^\T\mA\vx$、$\norm{\mM\vx-\vb}_2^2$、$\norm{\vx}_2$（{numref}`tab-01-02-vector-derivatives`）；迹、对数行列式和逆矩阵对矩阵的导数（{numref}`tab-01-02-matrix-derivatives`、{eq}`eq-01-02-matrix-gradient`）。
- 链式法则按复合的次序把雅可比矩阵相乘（{prf:ref}`thm-01-02-chain-rule`）；限制在直线上即得 {eq}`eq-01-02-line-restriction`。
- 泰勒定理给出带中间点的精确形式和 $o$ 形式（{prf:ref}`thm-01-02-taylor`）；{prf:ref}`ex-01-02-taylor-models`和{numref}`fig-01-02-taylor-models`比较了这些模型。
- 梯度与水平集正交（{prf:ref}`thm-01-02-gradient-level-set`），并给出最速上升方向；$-\nabla f$ 是最速下降方向（{prf:ref}`thm-01-02-steepest-ascent`、{numref}`fig-01-02-linear-level-lines`、{numref}`fig-01-02-paraboloid-contours`）。
- 令梯度为零，可得正规方程 {eq}`eq-01-02-normal-equations`、正定二次函数的极小点 {eq}`eq-01-02-quadratic-minimum` 以及拉格朗日函数的平稳性条件 {eq}`eq-01-02-lagrangian-gradient`。

### 常见错误

- **混用布局约定。** 本书中 $\nabla f$ 是列向量，$\vg:\R^n\to\R^m$ 的雅可比矩阵是 $m\times n$ 的。$\mA\vx$ 的雅可比矩阵是 $\mA$，但 $f(\mA\vx)$ 的梯度是 $\mA^\T\nabla f(\mA\vx)$；忘掉这个转置会得到长度错误的向量。
- **对非对称矩阵 $\mA$ 使用 $2\mA\vx$。** $\vx^\T\mA\vx$ 的梯度是 $(\mA+\mA^\T)\vx$；只有当 $\mA$ 对称时它才等于 $2\mA\vx$。
- **弄错泰勒定理中各项的取值点。** 在精确的二阶形式 {eq}`eq-01-02-taylor-2` 中，$f$ 和 $\nabla f$ 在展开点 $\bar\vx$ 处取值，只有黑塞矩阵在中间点 $\bar\vx+t\vd$ 处取值；因子 $\tfrac12$ 不能丢掉，而且 $\nabla f(\bar\vx)$ 是 $f$ 的梯度在 $\bar\vx$ 处的值，而不是对 $\bar\vx$ 求导。
- **把“沿梯度方向 $f$ 增加”理解为全局性的结论。** 它只对充分小的步长成立；一大步可能走过头。
- **以为 $-\nabla f$ 指向极小点。** 只有当等值线是圆时才如此；一般而言，它只是垂直于过当前点的等值线。
- **在 $\vzero$ 处对 $\norm{\vx}_2$ 求导。** 欧氏范数在原点处不可微；需要光滑函数时，应使用 $\norm{\vx}_2^2$。
- **把对称矩阵变量当作无约束的。** {numref}`tab-01-02-matrix-derivatives`把 $\mX$ 的所有元素视为相互独立。若 $\mX$ 限定为对称矩阵，则非对角元素的导数必须考虑 $x_{ij}=x_{ji}$。

## 参考文献

- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Appendix
  A.4: derivatives, gradient, chain rule, second derivative; Appendix
  A.5.5: Schur complement; Section 6.1: least-squares approximation.)
- Magnus, J. R., & Neudecker, H. (2019). *Matrix differential calculus with
  applications in statistics and econometrics* (3rd ed.). Wiley.
  <https://doi.org/10.1002/9781119541219> (Chapters 5–6: differentials,
  Jacobians, the chain rule and Taylor's theorem; Chapters 8–9: first
  differentials of the trace, the determinant and the inverse, and the
  layout of Jacobian matrices.)
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Section 2.1:
  Taylor's theorem; Appendix A.2: elements of analysis, directional
  derivatives.)
- Petersen, K. B., & Pedersen, M. S. (2012). *The matrix cookbook* (Version
  November 15, 2012). Technical University of Denmark. (Section 2:
  derivatives of determinants, inverses, traces and quadratic forms.)
- Strang, G. (2016). *Introduction to linear algebra* (5th ed.).
  Wellesley-Cambridge Press. (Section 4.3: least-squares approximations.)
