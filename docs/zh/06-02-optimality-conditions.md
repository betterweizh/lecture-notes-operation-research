# 最优性条件

## 概述

不把一个点与所有其他可行点逐一比较，怎样才能识别非线性规划问题的极小点？本章用导数回答这个问题。理论围绕两类方向展开：**可行方向**（feasible direction），沿它移动一小段不会离开 $\Omega$；**下降方向**（descent direction），沿它 $f$ 减小。在局部极小点处，没有一个方向能同时是这两类方向（{prf:ref}`thm-06-02-first-order`）。由这一原理导出一阶必要条件（内点处梯度为零，凸集上的变分不等式）、二阶必要条件（{prf:ref}`thm-06-02-second-order-necessary`）以及无约束问题的二阶充分条件（{prf:ref}`thm-06-02-second-order-sufficient`）。对于凸规划，一阶条件也是全局极小点的充分条件（{prf:ref}`thm-06-02-convex-sufficient`）。本章的条件用方向来表述；{doc}`06-03-kkt-conditions`把它们转化为关于约束函数的条件，即 Karush–Kuhn–Tucker（KKT）条件。本章以{doc}`06-01-nlp-models`中的非线性规划（nonlinear programming, NLP）模型和局部极小点、{doc}`01-02-matrix-calculus`中的梯度、黑塞矩阵和泰勒定理，以及{doc}`02-01-convex-sets`和{doc}`02-02-convex-functions`中的凸集和凸函数为基础。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 判断一个方向在多面体或凸集的某点处是否可行，以及它是否是 $f$ 的下降方向；
- 叙述并证明一阶必要条件：对所有可行方向有 $\nabla f(\vx^*)^\T\vd\ge0$，并由此推出内点处 $\nabla f(\vx^*)=\vzero$；
- 叙述并证明无约束极小点的二阶必要条件和二阶充分条件；
- 利用黑塞矩阵把驻点分类为局部极小点、局部极大点或鞍点；
- 证明对于凸规划，满足一阶条件的点是全局极小点。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\vd$ | $\R^n$ 中的方向；$d_j$ 为其分量 |
| $\alpha$, $\bar\alpha$ | 沿 $\vd$ 的步长；使 $\vx+\alpha\vd$ 保持可行的步长上限 |
| $D(\vx)$ | $\Omega$ 在 $\vx$ 处的可行方向集 |
| $\ve_j$ | $\R^n$ 的第 $j$ 个单位向量 |
| $\vartheta\in(0,1)$ | 泰勒定理中中间点的位置 |
| $\sigma$ | {prf:ref}`thm-06-02-second-order-sufficient`证明中的最小特征值 $\lambda_{\min}(\nabla^2f(\vx^*))$ |

本章始终假设 $f$ 在包含 $\Omega$ 的某个开集上连续可微，在出现黑塞矩阵时假设它二阶连续可微。问题为 $\min\{f(\vx):\vx\in\Omega\}$，其中 $\Omega\subseteq\R^n$；用约束函数描述 $\Omega$ 要到{doc}`06-03-kkt-conditions`才需要。

(sec-06-02-feasible-directions)=
## 可行方向

从一个可行点出发，只能沿某些方向移动而不离开 $\Omega$。下面的定义把这一点精确化；它要求从 $\vx$ 到 $\vx+\bar\alpha\vd$ 的整条线段都可行，而不只是其端点可行。

:::{prf:definition} 可行方向
:label: def-06-02-feasible-direction

设 $\vx\in\Omega$。若存在 $\bar\alpha>0$，使得

$$
\vx+\alpha\vd\in\Omega\qquad\text{对所有 }0\le\alpha\le\bar\alpha .
$$ (eq-06-02-feasible-direction)

则称向量 $\vd\in\R^n$ 为 $\Omega$ 在 $\vx$ 处的**可行方向**。$\vx$ 处的可行方向集记为 $D(\vx)$。
:::

若 $\vd\in D(\vx)$，则对每个 $t>0$ 有 $t\vd\in D(\vx)$（取 $\bar\alpha/t$ 即可），因此 $D(\vx)$ 是一个锥，起作用的只是 $\vd$ 的方向而非其长度。下面三种情形涵盖了大多数应用。

:::{prf:proposition} 三种常见情形下的可行方向
:label: prop-06-02-feasible-directions

1. 若 $\vx$ 是 $\Omega$ 的内点，则每个 $\vd\in\R^n$ 都是 $\vx$ 处的可行方向。
2. 若 $\Omega$ 是凸集且 $\vx,\vy\in\Omega$，则 $\vd=\vy-\vx$ 是 $\vx$ 处的可行方向（取 $\bar\alpha=1$）。
3. 若 $\Omega=\{\vx:\va_i^\T\vx\le b_i,\ i=1,\dots,m\}$ 是多面体，则 $\vd$ 是 $\vx\in\Omega$ 处的可行方向，当且仅当对每个满足 $\va_i^\T\vx=b_i$ 的 $i$ 有 $\va_i^\T\vd\le0$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

依次讨论三种情形。

1. 存在某个球 $B(\vx,\varepsilon)$ 含于 $\Omega$。对 $\vd\ne\vzero$ 取 $\bar\alpha=\varepsilon/(2\norm{\vd})$；则当 $0\le\alpha\le\bar\alpha$ 时 $\norm{\alpha\vd}\le\varepsilon/2<\varepsilon$，因此 $\vx+\alpha\vd\in B(\vx,\varepsilon)\subseteq\Omega$。
2. 当 $0\le\alpha\le1$ 时，$\vx+\alpha(\vy-\vx)=(1-\alpha)\vx+\alpha\vy$ 是 $\Omega$ 中两点的凸组合。
3. 对于满足 $\va_i^\T\vx=b_i$ 的约束（*有效*约束），$\va_i^\T(\vx+\alpha\vd)=b_i+\alpha\va_i^\T\vd$，它对某个 $\alpha>0$ 满足 $\le b_i$ 当且仅当 $\va_i^\T\vd\le0$。满足 $\va_i^\T\vx<b_i$ 的约束有松弛量 $s_i=b_i-\va_i^\T\vx>0$：若 $\va_i^\T\vd>0$，则当 $\alpha\le s_i/\va_i^\T\vd$ 时它仍满足；否则对所有 $\alpha\ge0$ 都满足。取 $\bar\alpha$ 为这有限个上界中的最小者，即得结论。
:::

第 3 部分说明，在多面体的点处，只有有效约束（active constraint）限制方向，并且是通过一个线性检验来限制的。对非线性约束的同样检验是{doc}`06-03-kkt-conditions`的主题，那里的情况更为微妙。

:::{prf:remark} 弯曲集合上的可行方向
:label: rem-06-02-curved

在弯曲的边界上，可行方向可能比人们预想的少得多。在单位圆盘 $\{\vx:x_1^2+x_2^2\le1\}$ 的点 $(1,0)$ 处，方向 $\vd=(d_1,d_2)$ 可行当且仅当 $d_1<0$ 或 $\vd=\vzero$：沿切线方向 $(0,1)$ 移动会立即离开圆盘，因为 $1+\alpha^2>1$。在由方程定义的单位*圆周* $\{\vx:x_1^2+x_2^2=1\}$ 上，任意点处唯一的可行方向是 $\vd=\vzero$，因为圆周上不含任何线段。此时{prf:ref}`thm-06-02-first-order`平凡地成立，不提供任何信息。对于这样的集合，必须沿可行*曲线*而不是直线移动；{doc}`06-03-kkt-conditions`借助切锥做到这一点，得到的条件对弯曲约束同样有效。本章的结果对内点、多面体以及（通过{prf:ref}`prop-06-02-feasible-directions`）凸集是精确的。
:::

(sec-06-02-descent-directions)=
## 下降方向

:::{prf:definition} 下降方向
:label: def-06-02-descent-direction

若 $\nabla f(\vx)^\T\vd<0$，则称向量 $\vd$ 为 $f$ 在 $\vx$ 处的**下降方向**。$\vx\in\Omega$ 处的**可行下降方向**（feasible descent direction）是既可行又下降的方向。
:::

数 $\nabla f(\vx)^\T\vd$ 是 $f$ 在 $\vx$ 处沿方向 $\vd$ 的方向导数（{prf:ref}`def-01-02-directional-derivative`）。它的符号有几何意义：$\nabla f(\vx)^\T\vd<0$ 当且仅当 $\vd$ 与梯度的夹角为钝角，即 $\vd$ 指向水平集 $\{\vy:f(\vy)=f(\vx)\}$ 的切平面上 $f$ 较小的那一侧。下面的引理说明了这一名称的合理性。

:::{prf:lemma} 沿下降方向 $f$ 减小
:label: lem-06-02-descent

若 $\nabla f(\vx)^\T\vd<0$，则存在 $\bar\alpha>0$，使得对所有 $0<\alpha\le\bar\alpha$ 有 $f(\vx+\alpha\vd)<f(\vx)$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

令 $\varphi(\alpha)=f(\vx+\alpha\vd)$。由链式法则，$\varphi'(0)=\nabla f(\vx)^\T\vd<0$，由导数的定义，

$$
\lim_{\alpha\downarrow0}\frac{f(\vx+\alpha\vd)-f(\vx)}{\alpha}=\varphi'(0)<0 .
$$

因此对所有充分小的 $\alpha>0$，该差商为负，这就是所要证明的结论。
:::

逆命题不成立：即使 $\nabla f(\vx)^\T\vd=0$，$f$ 也可能沿 $\vd$ 减小（对于 $f(x)=-x^2$，在 $x=0$ 处每个 $d\ne0$ 都使 $f$ 减小）。满足 $\nabla f(\vx)^\T\vd=0$ 的方向在一阶意义下无法判定，需要用二阶导数来判断（{ref}`sec-06-02-second-order`）。在单位长度的方向中，由柯西–施瓦茨不等式 $\nabla f(\vx)^\T\vd\ge-\norm{\nabla f(\vx)}$ 可知，$\vd=-\nabla f(\vx)/\norm{\nabla f(\vx)}$ 给出最负的方向导数；这一**最速下降方向**（steepest-descent direction）是{doc}`06-05-unconstrained-algorithms`中方法的基础。

:::{prf:example} 三角形上的可行方向与下降方向
:label: ex-06-02-directions

设 $\Omega=\{\vx\in\R^2:x_1+x_2\le2,\ x_1\ge0,\ x_2\ge0\}$ 是以 $(0,0)$、$(2,0)$、$(0,2)$ 为顶点的三角形，$f(\vx)=(x_1-1)^2+(x_2-2)^2$。

(a) 在点 $(1,0.5)$ 和顶点 $(2,0)$ 处，哪些方向是可行的？$\vd=(-1,0.5)$ 和 $\vd=(1,1.5)$ 在 $(2,0)$ 处是否可行？

(b) 判断 $(1,0.5)$ 处的方向 $\vd_1=(-0.4,1)$、$\vd_2=(0,1)$、$\vd_3=(0.45,1)$、$\vd_4=(-1,-0.8)$ 以及 $(2,0)$ 处的方向 $\vd_5=(-1,0.5)$、$\vd_6=(0.6,1)$、$\vd_7=(1,0.13)$ 是否为下降方向、是否可行。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in


**(a)** 点 $(1,0.5)$ 严格满足全部三个约束。它到三条边界线的距离分别为 $1$（到 $x_1=0$）、$0.5$（到 $x_2=0$）和 $\abs{1+0.5-2}/\sqrt2\approx0.354$（到 $x_1+x_2=2$），因此以它为中心、半径为 $0.354$ 的球含于 $\Omega$。该点是内点，由{prf:ref}`prop-06-02-feasible-directions`，那里每个方向都是可行的：无论方向如何，步长小于 $0.354/\norm{\vd}$ 的一步都留在球内。

在 $(2,0)$ 处，约束 $x_1+x_2\le2$ 和 $x_2\ge0$（即 $-x_2\le0$）是有效约束，而 $x_1\ge0$ 不是。由{prf:ref}`prop-06-02-feasible-directions`的第 3 部分，可行方向恰为

$$
D(2,0)=\{\vd:\ d_1+d_2\le0,\ -d_2\le0\}.
$$

对于 $\vd=(-1,0.5)$：$d_1+d_2=-0.5\le0$ 且 $d_2=0.5\ge0$，因此它是可行的；事实上，当 $0\le\alpha\le2$ 时 $(2,0)+\alpha(-1,0.5)\in\Omega$，取 $\alpha=1$ 即到达 $(1,0.5)$。对于 $\vd=(1,1.5)$：$d_1+d_2=2.5>0$，因此任何一步都会穿过边 $x_1+x_2=2$ 离开 $\Omega$；它不可行。

**(b)** 梯度为 $\nabla f(\vx)=(2(x_1-1),\,2(x_2-2))^\T$，因此 $\nabla f(1,0.5)=(0,-3)^\T$，$\nabla f(2,0)=(2,-4)^\T$。

在 $(1,0.5)$ 处，$\nabla f^\T\vd=-3d_2$：一个方向是下降方向当且仅当它指向上方，即 $d_2>0$。因此 $\vd_1,\vd_2,\vd_3$（每个都有 $\nabla f^\T\vd=-3$）是下降方向，而 $\vd_4$（$\nabla f^\T\vd_4=2.4$）是上升方向。由于该点是内点，四个方向都可行。

在 $(2,0)$ 处：

```{list-table} 顶点 $(2,0)$ 处的方向。
:name: tab-06-02-directions-vertex
:header-rows: 1

* - $\vd$
  - $\nabla f(2,0)^\T\vd$
  - $d_1+d_2\le0$，$d_2\ge0$？
  - 结论
* - $\vd_5=(-1,0.5)$
  - $-2-2=-4$
  - 是
  - 可行下降方向
* - $\vd_6=(0.6,1)$
  - $1.2-4=-2.8$
  - 否（$d_1+d_2=1.6$）
  - 下降，但不可行
* - $\vd_7=(1,0.13)$
  - $2-0.52=1.48$
  - 否（$d_1+d_2=1.13$）
  - 上升且不可行
```

只有 $\vd_5$ 是可行下降方向。既然存在可行下降方向，$(2,0)$ 就不是局部极小点：由{prf:ref}`lem-06-02-descent`，沿 $\vd_5$ 走一小步仍然可行且使 $f$ 减小。在 $(1,0.5)$ 处也是如此。那么极小点在哪里？它是 $\Omega$ 中离 $(1,2)$ 最近的点，即 $\vx^*=(0.5,1.5)$，$f(\vx^*)=0.5$；{ref}`sec-06-02-convex`将证明这一点。{numref}`fig-06-02-directions`画出了这两个点以及 $f$ 的等值线。

```{figure} figures/06-02-directions.svg
:name: fig-06-02-directions
:width: 680px
:class: dark-light
:alt: 两幅子图，显示三角形 x1 + x2 <= 2 以及以 (1,2) 为中心的 f 的圆形等值线；左图为 (1,0.5) 处的方向 d1 至 d4 和梯度，右图为顶点 (2,0) 处的方向 d5 至 d7 和梯度

{prf:ref}`ex-06-02-directions` 中的方向。左：在内点 $(1,0.5)$ 处每个方向都可行，$\vd_1,\vd_2,\vd_3$ 与 $\nabla f=(0,-3)$ 成钝角（它们指向等值线 $f=2.25$ 的内侧）。右：在顶点 $(2,0)$ 处，只有 $\vd_5$ 既指向 $\Omega$ 内部又与 $\nabla f=(2,-4)$ 成钝角；$\vd_6$ 是离开 $\Omega$ 的下降方向。
```
:::

(sec-06-02-first-order)=
## 一阶必要条件

局部极小点处不存在可行下降方向。用方向导数来表述，这就是基本的一阶条件。

:::{prf:theorem} 一阶必要条件
:label: thm-06-02-first-order

若 $\vx^*$ 是 $f$ 在 $\Omega$ 上的局部极小点，则

$$
\nabla f(\vx^*)^\T\vd\ge0\qquad\text{对每个可行方向 }\vd\in D(\vx^*).
$$ (eq-06-02-first-order)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

假设某个 $\vd\in D(\vx^*)$ 满足 $\nabla f(\vx^*)^\T\vd<0$。设 $\varepsilon>0$ 是局部极小点定义（{prf:ref}`def-06-01-local-global`）中的半径。由 {eq}`eq-06-02-feasible-direction` 和{prf:ref}`lem-06-02-descent`，存在 $\bar\alpha_1,\bar\alpha_2>0$，使得当 $\alpha\le\bar\alpha_1$ 时 $\vx^*+\alpha\vd\in\Omega$，当 $0<\alpha\le\bar\alpha_2$ 时 $f(\vx^*+\alpha\vd)<f(\vx^*)$。选取满足 $\alpha\norm{\vd}<\varepsilon$ 的 $0<\alpha\le\min\{\bar\alpha_1,\bar\alpha_2\}$。则 $\vx^*+\alpha\vd$ 是与 $\vx^*$ 距离小于 $\varepsilon$ 且目标函数值更小的可行点，与局部极小性矛盾。
:::

从几何上看，{eq}`eq-06-02-first-order` 说的是：在局部极小点处，每个可行方向与梯度的夹角至多为 $90^\circ$；下降最快的方向 $-\nabla f(\vx^*)$ 指向 $\Omega$ 之外，或者与所有可行方向正交。

在内点处所有方向都可行，条件就变成熟悉的形式。

:::{prf:corollary} 内部极小点处的平稳性
:label: cor-06-02-stationary

若 $\vx^*$ 是 $f$ 在 $\Omega$ 上的局部极小点，且是 $\Omega$ 的内点（特别地，若 $\Omega=\R^n$），则

$$
\nabla f(\vx^*)=\vzero .
$$ (eq-06-02-stationary)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

由{prf:ref}`prop-06-02-feasible-directions`，$\ve_j$ 和 $-\ve_j$ 都是可行方向，由 {eq}`eq-06-02-first-order` 得 $\nabla f(\vx^*)^\T\ve_j=\partial f(\vx^*)/\partial x_j\ge0$ 以及 $\nabla f(\vx^*)^\T(-\ve_j)=-\partial f(\vx^*)/\partial x_j\ge0$。因此对 $j=1,\dots,n$ 有 $\partial f(\vx^*)/\partial x_j=0$。
:::

满足 $\nabla f(\vx)=\vzero$ 的点称为 $f$ 的**驻点**（stationary point，又称临界点）。驻点可能是局部极小点、局部极大点，或者两者都不是；两者都不是的驻点称为**鞍点**（saddle point）。{eq}`eq-06-02-stationary` 是含 $n$ 个未知量的 $n$ 个方程组成的方程组，求解它是寻找无约束极小点的第一步。在{prf:ref}`ex-06-01-three-feasible-sets`的问题 1 中，$\nabla f=(2(x_1-2),2(x_2-2))$ 只在 $(2,2)$ 处为零，这正是那里求得的极小点。

当 $\Omega$ 是凸集时，每个 $\vy-\vx^*$（$\vy\in\Omega$）都是可行方向，由 {eq}`eq-06-02-first-order` 推出**变分不等式**（variational inequality）

$$
\nabla f(\vx^*)^\T(\vy-\vx^*)\ge0\qquad\text{对所有 }\vy\in\Omega .
$$ (eq-06-02-variational-inequality)

两个特殊情形说明了它的作用。

- *区间。* 对于一元问题 $\min\{f(x):a\le x\le b\}$，{eq}`eq-06-02-variational-inequality` 说的是：若 $a<x^*<b$，则 $f'(x^*)=0$；若 $x^*=a$，则 $f'(a)\ge0$；若 $x^*=b$，则 $f'(b)\le0$。{doc}`06-01-nlp-models`中的函数 $(x^2)^x$ 在 $[-5,1.7]$ 上的全局极小点是左端点，且 $f'(-5)\approx5.3\times10^{-7}>0$，与该条件一致。
- *非负卦限。* 对于 $\Omega=\{\vx:\vx\ge\vzero\}$，方向 $\ve_j$ 总是可行的，而当 $x^*_j>0$ 时 $-\ve_j$ 也可行。因此局部极小点满足

  $$
  \begin{aligned}
  &\frac{\partial f(\vx^*)}{\partial x_j}\ge0, && j=1,\dots,n,\\
  &x_j^*\ge0, && j=1,\dots,n,\\
  &x_j^*\,\frac{\partial f(\vx^*)}{\partial x_j}=0, && j=1,\dots,n .
  \end{aligned}
  $$ (eq-06-02-orthant)

  每个变量要么处于其界上且偏导数非负，要么为正且偏导数为零。这种*互补性*将在{doc}`06-03-kkt-conditions`的 KKT 条件中再次出现。

(sec-06-02-second-order)=
## 二阶条件

梯度为零并不能保证是极小点。函数 $f(\vx)=x_1^2-x_2^2$ 满足 $\nabla f(\vzero)=\vzero$，但它沿 $x_1$ 轴增大、沿 $x_2$ 轴减小：原点是鞍点，形如山口。一阶导数无法区分鞍点和极小点，因为它们只描述切平面，而在两种情况下切平面都是水平的。由黑塞矩阵描述的曲率则可以区分。

:::{prf:theorem} 二阶必要条件
:label: thm-06-02-second-order-necessary

设 $f$ 二阶连续可微，$\vx^*$ 是 $f$ 在 $\Omega$ 上的局部极小点。则对每个可行方向 $\vd\in D(\vx^*)$：

1. $\nabla f(\vx^*)^\T\vd\ge0$；
2. 若 $\nabla f(\vx^*)^\T\vd=0$，则 $\vd^\T\nabla^2f(\vx^*)\vd\ge0$。

特别地，若 $\vx^*$ 是 $\Omega$ 的内点（或 $\Omega=\R^n$），则 $\nabla f(\vx^*)=\vzero$ 且 $\nabla^2f(\vx^*)$ 半正定。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

第 1 部分即{prf:ref}`thm-06-02-first-order`。对于第 2 部分，设 $\vd\in D(\vx^*)$ 满足 $\nabla f(\vx^*)^\T\vd=0$，并假设 $\vd^\T\nabla^2f(\vx^*)\vd<0$。由于 $\nabla^2f$ 连续，存在 $\alpha_1>0$，使得当 $0\le\alpha\le\alpha_1$ 时 $\vd^\T\nabla^2f(\vx^*+\alpha\vd)\vd<0$。由带拉格朗日余项的泰勒定理（{prf:ref}`thm-01-02-taylor`），对每个这样的 $\alpha$，存在 $\vartheta\in(0,1)$ 使得

$$
\begin{aligned}
f(\vx^*+\alpha\vd)
&=f(\vx^*)+\alpha\nabla f(\vx^*)^\T\vd+\tfrac12\alpha^2\vd^\T\nabla^2f(\vx^*+\vartheta\alpha\vd)\vd\\
&=f(\vx^*)+\tfrac12\alpha^2\vd^\T\nabla^2f(\vx^*+\vartheta\alpha\vd)\vd ,
\end{aligned}
$$

其中一阶项因 $\nabla f(\vx^*)^\T\vd=0$ 而消失。由于 $\vartheta\alpha\le\alpha_1$，最后一项为负，因此对所有充分小的 $\alpha>0$ 有 $f(\vx^*+\alpha\vd)<f(\vx^*)$。这些点是可行的（当 $\alpha\le\bar\alpha$ 时）且任意接近 $\vx^*$，矛盾。

在内点处，由{prf:ref}`cor-06-02-stationary`有 $\nabla f(\vx^*)=\vzero$，每个 $\vd$ 都可行且满足 $\nabla f(\vx^*)^\T\vd=0$，第 2 部分说明对所有 $\vd\in\R^n$ 有 $\vd^\T\nabla^2f(\vx^*)\vd\ge0$，这就是半正定性。
:::

对称矩阵的半正定性等价于特征值非负，即 $\lambda_i(\nabla^2f(\vx^*))\ge0$（{prf:ref}`thm-01-01-eigen-test`）。因此，黑塞矩阵有负特征值的驻点不是局部极小点；如果它还有正特征值，那么它也不是局部极大点，而是鞍点。

把半正定加强为正定，就得到一个充分条件。

:::{prf:theorem} 二阶充分条件（无约束）
:label: thm-06-02-second-order-sufficient

设 $f$ 在 $\R^n$ 上二阶连续可微。若

$$
\begin{aligned}
\nabla f(\vx^*)&=\vzero,\\
\nabla^2f(\vx^*)&\succ\vzero ,
\end{aligned}
$$ (eq-06-02-sosc)

则 $\vx^*$ 是 $f$ 的严格局部极小点。当 $\vx^*$ 是 $\Omega$ 的内点时，这一结论对 $\min_{\vx\in\Omega}f$ 同样成立。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

令 $\sigma=\lambda_{\min}(\nabla^2f(\vx^*))>0$，于是对所有 $\vd$ 有 $\vd^\T\nabla^2f(\vx^*)\vd\ge\sigma\norm{\vd}^2$。由 $\nabla^2f$ 的连续性，存在 $\varepsilon>0$，使得对所有 $\vy\in B(\vx^*,\varepsilon)$，谱范数 $\norm{\nabla^2f(\vy)-\nabla^2f(\vx^*)}\le\sigma/2$；对这样的 $\vy$，

$$
\begin{aligned}
\vd^\T\nabla^2f(\vy)\vd
&\ge\vd^\T\nabla^2f(\vx^*)\vd-\norm{\nabla^2f(\vy)-\nabla^2f(\vx^*)}\,\norm{\vd}^2\\
&\ge\tfrac{\sigma}{2}\norm{\vd}^2 .
\end{aligned}
$$

现在设 $\vx\in B(\vx^*,\varepsilon)$，$\vx\ne\vx^*$，$\vd=\vx-\vx^*$。由泰勒定理，存在 $\vartheta\in(0,1)$ 使得

$$
\begin{aligned}
f(\vx)&=f(\vx^*)+\nabla f(\vx^*)^\T\vd+\tfrac12\vd^\T\nabla^2f(\vx^*+\vartheta\vd)\vd\\
&\ge f(\vx^*)+\tfrac{\sigma}{4}\norm{\vx-\vx^*}^2>f(\vx^*),
\end{aligned}
$$

这是因为 $\nabla f(\vx^*)=\vzero$ 且 $\vx^*+\vartheta\vd\in B(\vx^*,\varepsilon)$。
:::

必要条件（半正定）与充分条件（正定）之间存在空隙。在黑塞矩阵奇异且半定的驻点处，二阶导数什么也判定不了：$x^4$、$-x^4$ 和 $x^3$ 都满足 $f'(0)=f''(0)=0$，而原点分别是极小点、极大点和两者都不是。反过来，充分条件也不是必要的：$x^4$ 在 $0$ 处有严格极小点，但 $f''(0)=0$。

对于约束问题，{prf:ref}`thm-06-02-second-order-necessary`的条件 2 只涉及 $f$ 在一阶意义下平坦的可行方向；满足 $\nabla f(\vx^*)^\T\vd>0$ 的方向已经使 $f$ 增大，无需检验曲率。当约束本身是弯曲的时，可行方向无法反映边界的曲率，正确的二阶条件要用拉格朗日函数的黑塞矩阵，见{doc}`06-03-kkt-conditions`。

:::{prf:example} 一个鞍点及其约束形式
:label: ex-06-02-saddle

设 $f(\vx)=x_1^2+x_1x_2-x_1-x_2$。

(a) 求 $f$ 在 $\R^2$ 上的驻点，并判断它们是否为局部极小点。

(b) 现在在 $x_1\ge1$，$x_2\ge0$ 的约束下最小化 $f$。在 $\vx^*=(1,0)$ 处检验一阶和二阶必要条件，并判断它是否为极小点。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** 梯度和黑塞矩阵为

$$
\begin{aligned}
\nabla f(\vx)&=\begin{pmatrix}2x_1+x_2-1\\x_1-1\end{pmatrix},\\
\nabla^2f(\vx)&=\begin{pmatrix}2&1\\1&0\end{pmatrix}.
\end{aligned}
$$

令 $\nabla f=\vzero$ 得 $x_1=1$，进而 $x_2=1-2x_1=-1$：唯一的驻点是 $(1,-1)$，$f(1,-1)=0$。黑塞矩阵是常数矩阵，其行列式为 $2\cdot0-1\cdot1=-1<0$，因此两个特征值异号：$\lambda_{1,2}(\nabla^2f)=1\mp\sqrt2\approx-0.414,\ 2.414$。它是不定的，二阶必要条件不成立，$(1,-1)$ 不是局部极小点；对 $-f$ 作同样的论证可知它也不是局部极大点。它是鞍点。*（此例的一种较早写法在“$(1,0)$ 处”计算黑塞矩阵；驻点实为 $(1,-1)$，但由于黑塞矩阵是常数，结论相同。）*

鞍点也可以直接看出。因式分解得

$$
\begin{aligned}
f(\vx)&=x_1(x_1+x_2)-(x_1+x_2)\\
&=(x_1-1)(x_1+x_2),
\end{aligned}
$$

因此在相交于 $(1,-1)$ 的两条直线 $x_1=1$ 和 $x_1+x_2=0$ 上 $f=0$，且 $f$ 穿过每条直线时都变号。沿 $\vd=(1,0)$ 有 $f(1+\alpha,-1)=\alpha^2>0$，而沿 $\vd=(1,-2)$ 有 $f(1+\alpha,-1-2\alpha)=\alpha\cdot(-\alpha)=-\alpha^2<0$。这与曲率 $\vd^\T\nabla^2f\vd=2$ 和 $-2$ 一致（对于二次函数，在驻点处恰有 $f(\vx^*+\alpha\vd)=f(\vx^*)+\tfrac12\alpha^2\vd^\T\nabla^2f\vd$）。

**(b)** 在 $\vx^*=(1,0)$ 处，两个约束 $x_1\ge1$ 和 $x_2\ge0$ 都是有效约束，又因为它们是线性的，由{prf:ref}`prop-06-02-feasible-directions`得可行方向为满足 $d_1\ge0$，$d_2\ge0$ 的 $\vd=(d_1,d_2)$。梯度为 $\nabla f(1,0)=(2+0-1,\,1-1)^\T=(1,0)^\T\ne\vzero$。

- 一阶：对每个可行的 $\vd$ 有 $\nabla f(\vx^*)^\T\vd=d_1\ge0$。条件成立。
- 二阶：$\nabla f(\vx^*)^\T\vd=0$ 意味着 $d_1=0$，于是 $\vd=(0,d_2)$，$d_2\ge0$，且 $\vd^\T\nabla^2f\vd=2d_1^2+2d_1d_2=0\ge0$。条件成立。

两个必要条件都成立，但它们并不能证明最优性（黑塞矩阵不定，而且{prf:ref}`thm-06-02-second-order-sufficient`本来就不适用于边界点）。因式分解解决了这个问题：在 $\Omega$ 上 $x_1-1\ge0$ 且 $x_1+x_2\ge1>0$，因此 $f(\vx)\ge0=f(1,0)$。所以 $(1,0)$ 是全局极小点。它不是唯一的：对每个 $x_2\ge0$ 有 $f(1,x_2)=0$，因此整条射线 $\{(1,x_2):x_2\ge0\}$ 都由全局极小点组成。这说明了为什么沿可行方向 $(0,1)$ 的曲率恰好为零：$f$ 沿该方向是常数。

```{figure} figures/06-02-saddle.svg
:name: fig-06-02-saddle
:width: 520px
:class: dark-light
:alt: f = (x1 - 1)(x1 + x2) 的等值线图，鞍点 (1,-1) 周围为双曲线形水平曲线，零水平线 x1 = 1 和 x1 + x2 = 0，以及阴影表示的可行象限 x1 >= 1, x2 >= 0，其左边缘是一条由极小点组成的射线

{prf:ref}`ex-06-02-saddle` 中 $f(\vx)=(x_1-1)(x_1+x_2)$ 的等值线。零水平集由两条虚线组成，它们相交于鞍点 $(1,-1)$；在左右两个扇形区域内 $f>0$，在上下两个扇形区域内 $f<0$。在阴影所示的可行集 $x_1\ge1$，$x_2\ge0$ 上，最小值 $0$ 在整条边 $x_1=1$（粗线）上取得。
```
:::

(sec-06-02-procedure)=
## 条件的使用

对于极小点为内点的问题，这些定理可以组合成一套常规步骤；稍加修改也适用于简单的界约束。

:::{prf:algorithm} 求无约束极小点并分类
:label: alg-06-02-classify

**输入：** 开集 $U\subseteq\R^n$（通常 $U=\R^n$）上的二阶连续可微函数 $f$。

**输出：** $f$ 在 $U$ 中的局部极小点，以及（在可能时）全局极小点。

1. **驻点。** 在 $U$ 中求解 $\nabla f(\vx)=\vzero$（{prf:ref}`cor-06-02-stationary`）。每个局部极小点都在这些解之中。
2. **分类。** 在每个驻点处计算 $\nabla^2f$ 及其特征值（或顺序主子式）：正定 $\Rightarrow$ 严格局部极小点（{prf:ref}`thm-06-02-second-order-sufficient`）；负定 $\Rightarrow$ 严格局部极大点；不定 $\Rightarrow$ 鞍点（对 $f$ 和 $-f$ 应用{prf:ref}`thm-06-02-second-order-necessary`）；奇异半定 $\Rightarrow$ 无法判定，需直接考察 $f$。
3. **全局极小点。** 若 $f$ 是凸函数，则任何驻点都是全局极小点（{prf:ref}`thm-06-02-convex-sufficient`）。否则比较各局部极小点处的函数值，并确认极小点确实存在，例如利用强制性（{prf:ref}`thm-06-01-existence`）；在有界区域上还要与边界上的值比较。
:::

第 2 步依赖于极大点的二阶条件，它们可以通过对 $-f$ 应用上述定理得到：在局部极大点处 $\nabla f=\vzero$ 且 $\nabla^2f\preceq\vzero$；而 $\nabla f=\vzero$ 且 $\nabla^2f\prec\vzero$ 给出严格局部极大点。第 3 步不能省略：一个函数可以只有一个驻点且它是严格局部极小点，但函数仍然无下界（例如 $x_1^2+x_2^2(1-x_1)^3$，这是一个经典的二元反例），因此没有其他信息时，局部检验永远不能证明全局极小性。

对于{numref}`fig-06-01-x-power-2x`中区间 $[-5,1.7]$ 上的函数 $f(x)=(x^2)^x$，这套步骤如下进行。对 $f'(x)=2f(x)(\ln\abs{x}+1)$ 再求一次导数得

$$
f''(x)=4f(x)\bigl(\ln\abs{x}+1\bigr)^2+\frac{2f(x)}{x},\qquad x\ne0 .
$$

在驻点 $x=\pm1/e$ 处第一项为零，因此 $f''(1/e)=2e\cdot e^{-2/e}\approx2.605>0$（严格局部极小点），$f''(-1/e)=-2e\cdot e^{2/e}\approx-11.35<0$（严格局部极大点）。$f$ 不可微的点 $x=0$，以及 $f'$ 不必为零的端点 $-5$ 和 $1.7$，需要单独考察：{eq}`eq-06-02-variational-inequality` 的端点条件在 $x=-5$ 处成立（$f'(-5)>0$），比较函数值，$f(-5)\approx1.0\times10^{-7}$ 小于 $f(1/e)\approx0.479$。因此 $1/e$ 处的局部极小点不是全局极小点，与{doc}`06-01-nlp-models`中的结论一致。

对于线性约束问题，同样的方案可以借助{prf:ref}`thm-06-02-second-order-necessary`的条件用于边界点：先在可行方向锥上检验 $\nabla f(\vx^*)^\T\vd\ge0$，再在其中满足 $\nabla f(\vx^*)^\T\vd=0$ 的部分上检验曲率。若对每个这样的非零 $\vd$ 曲率都*为正*，则该点是严格局部极小点；这是{doc}`06-03-kkt-conditions`中所证二阶充分条件（{prf:ref}`thm-06-03-second-order`）在线性约束情形下的特例。在{prf:ref}`ex-06-02-saddle`的 (b) 中，沿临界方向 $(0,1)$ 的曲率为零，因此这一检验无法得出结论，这与那里求得的由非严格极小点组成的射线是一致的。

数值算法把同样的条件用作终止准则。无约束最小化方法在 $\norm{\nabla f(\vx^{(k)})}$ 小于容差时停止，约束问题的方法则在{doc}`06-03-kkt-conditions`的 KKT 条件近似成立时停止。这样的点只是近似驻点；它是否为极小点，要像上面那样由二阶条件或凸性来判定（{doc}`06-04-line-search`，{doc}`06-05-unconstrained-algorithms`）。

(sec-06-02-convex)=
## 凸问题：一阶条件是充分的

对于凸规划，一阶必要条件也是充分条件，并且它保证的是*全局*极小点。

:::{prf:theorem} 凸规划的最优性
:label: thm-06-02-convex-sufficient

设 $\Omega$ 是非空凸集，$f$ 在包含 $\Omega$ 的某个开集上连续可微，且在 $\Omega$ 上是凸函数。点 $\vx^*\in\Omega$ 是 $f$ 在 $\Omega$ 上的全局极小点，当且仅当

$$
\nabla f(\vx^*)^\T(\vy-\vx^*)\ge0\qquad\text{对所有 }\vy\in\Omega .
$$ (eq-06-02-convex-condition)

特别地，对于 $\Omega=\R^n$（或 $\vx^*$ 为内点），$\vx^*$ 是全局极小点当且仅当 $\nabla f(\vx^*)=\vzero$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

*必要性。* 全局极小点是局部极小点，且每个 $\vy-\vx^*$ 都是可行方向（{prf:ref}`prop-06-02-feasible-directions`），因此 {eq}`eq-06-02-convex-condition` 由{prf:ref}`thm-06-02-first-order`推出。

*充分性。* 假设 {eq}`eq-06-02-convex-condition` 成立，但某个 $\vy\in\Omega$ 满足 $f(\vy)<f(\vx^*)$。可微凸函数位于其切平面之上（{prf:ref}`thm-02-02-first-order`）：

$$
f(\vy)\ge f(\vx^*)+\nabla f(\vx^*)^\T(\vy-\vx^*)\ge f(\vx^*),
$$

矛盾。$\Omega=\R^n$ 的情形与{prf:ref}`cor-06-02-stationary`一样，取 $\vy=\vx^*\pm\ve_j$ 即得。
:::

结合以下事实：严格凸函数至多有一个极小点，以及在凸集上黑塞矩阵半正定的二阶可微函数在该集合上是凸的（{prf:ref}`thm-02-02-second-order`），这个定理给出了处理凸问题的完整方法：找到一个满足一阶条件的点；它就是全局极小点，并且当 $f$ 严格凸时它是唯一的。对于 $\mQ\succ\vzero$ 的二次函数 $f(\vx)=\vc^\T\vx+\tfrac12\vx^\T\mQ\vx$，其在 $\R^n$ 上唯一的全局极小点是 $\mQ\vx=-\vc$ 的解。

{prf:ref}`ex-06-02-directions`的目标函数是凸的，可行集是多面体，因此定理适用。在 $\vx^*=(0.5,1.5)$ 处，$\nabla f(\vx^*)=(2(0.5-1),\,2(1.5-2))^\T=(-1,-1)^\T$，对每个 $\vy\in\Omega$ 有

$$
\begin{aligned}
\nabla f(\vx^*)^\T(\vy-\vx^*)&=-(y_1-0.5)-(y_2-1.5)\\
&=2-(y_1+y_2)\ge0,
\end{aligned}
$$

这是因为在 $\Omega$ 上 $y_1+y_2\le2$。因此 $(0.5,1.5)$ 是全局极小点，$f^*=0.25+0.25=0.5$，又因为 $f$ 严格凸，它是唯一的。从几何上看，$-\nabla f(\vx^*)=(1,1)$ 垂直于边 $x_1+x_2=2$ 并指向 $\Omega$ 之外；{doc}`06-03-kkt-conditions`把这一观察发展为 KKT 条件。

## 条件汇总

```{list-table} 本章的最优性条件。“必要”条件在每个局部极小点处成立；“充分”条件保证得到局部极小点。
:name: tab-06-02-conditions
:header-rows: 1

* -
  - 无约束（$\Omega=\R^n$ 或 $\vx^*$ 为内点）
  - 有约束（可行方向 $\vd\in D(\vx^*)$）
* - 一阶必要条件
  - $\nabla f(\vx^*)=\vzero$
  - 对所有 $\vd\in D(\vx^*)$ 有 $\nabla f(\vx^*)^\T\vd\ge0$；在约束规范下为 KKT 条件（{doc}`06-03-kkt-conditions`）
* - 二阶必要条件
  - $\nabla f(\vx^*)=\vzero$ 且 $\nabla^2f(\vx^*)\succeq\vzero$
  - 另外：$\nabla f(\vx^*)^\T\vd=0\Rightarrow\vd^\T\nabla^2f(\vx^*)\vd\ge0$；二阶 KKT 条件（{doc}`06-03-kkt-conditions`）
* - 二阶充分条件
  - $\nabla f(\vx^*)=\vzero$ 且 $\nabla^2f(\vx^*)\succ\vzero$（严格局部极小点）
  - 二阶 KKT 充分条件（{doc}`06-03-kkt-conditions`）
* - $f$ 与 $\Omega$ 均为凸
  - $\nabla f(\vx^*)=\vzero$ 是全局极小点的充要条件
  - 对所有 $\vy\in\Omega$ 有 $\nabla f(\vx^*)^\T(\vy-\vx^*)\ge0$ 是全局极小点的充要条件
```

## 本章回顾

### 小结

- 可行方向使整条线段留在 $\Omega$ 中（{prf:ref}`def-06-02-feasible-direction`）；在内点处所有方向都可行，在多面体的点处，可行方向恰为对有效约束满足 $\va_i^\T\vd\le0$ 的方向（{prf:ref}`prop-06-02-feasible-directions`）。
- 下降方向满足 $\nabla f(\vx)^\T\vd<0$，步长足够小时使 $f$ 减小（{prf:ref}`lem-06-02-descent`）；$-\nabla f$ 是最速的下降方向。
- 一阶必要条件：不存在可行下降方向，即 $\nabla f(\vx^*)^\T\vd\ge0$（{prf:ref}`thm-06-02-first-order`）；在内点处 $\nabla f(\vx^*)=\vzero$（{prf:ref}`cor-06-02-stationary`）；在凸集上为变分不等式 {eq}`eq-06-02-variational-inequality`。
- 二阶必要条件：沿斜率为零的可行方向曲率非负；在内部极小点处黑塞矩阵半正定（{prf:ref}`thm-06-02-second-order-necessary`）。
- 二阶充分条件：梯度为零且黑塞矩阵正定，则为严格局部极小点（{prf:ref}`thm-06-02-second-order-sufficient`）。
- {prf:ref}`alg-06-02-classify`把这些条件组合起来：求解 $\nabla f=\vzero$，用黑塞矩阵分类，再利用凸性、强制性或函数值比较确定全局最优性。
- 对于凸规划，一阶条件是全局极小点的充要条件（{prf:ref}`thm-06-02-convex-sufficient`）。
- {numref}`tab-06-02-conditions`汇总了这些条件。

### 常见错误

- **到 $\nabla f(\vx^*)=\vzero$ 就停止。** 驻点可能是鞍点（{prf:ref}`ex-06-02-saddle`）或极大点。应检验黑塞矩阵或凸性。
- **要求边界点处 $\nabla f(\vx^*)=\vzero$。** 在约束极小点处梯度通常不为零；只要求沿可行方向的方向导数非负。
- **由黑塞矩阵半正定断定为极小点。** 半正定是必要的，正定是充分的；介于两者之间时，仅凭二阶导数得不出任何结论。
- **用子式检验非对称矩阵的有定性。** 二阶连续可微函数的黑塞矩阵是对称的；顺序主子式只能判定对称矩阵的*有定性*（{prf:ref}`thm-01-01-sylvester`），而顺序主子式非负并不意味着半正定（{prf:ref}`thm-01-01-psd-minors`）。
- **把“必要”当作“充分”。** {prf:ref}`ex-06-02-saddle`中的点 $(1,0)$ 满足两个必要条件；仅凭这一点不能说它是极小点。它的最优性需要另外的论证。


## 参考文献

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Section 3.4,
  minima of convex functions; Section 4.1, optimality conditions for
  unconstrained problems; Section 4.2, feasible and improving directions.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Section
  4.2.3, optimality criterion for differentiable objectives.)
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming*
  (4th ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3>
  (Sections 7.1–7.5, first- and second-order conditions, convex and concave
  functions, minimization of convex functions.)
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Section 2.1,
  recognizing a local minimum: Theorems 2.2–2.5.)
- Sheffi, Y. (1985). *Urban transportation networks: Equilibrium analysis
  with mathematical programming methods*. Prentice-Hall. (Chapter 2, basic
  concepts in minimization problems.)
