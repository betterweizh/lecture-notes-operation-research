# 最优性条件与 KKT 条件

本章讲解非线性规划（nonlinear programming, NLP）前三章{doc}`06-01-nlp-models`、{doc}`06-02-optimality-conditions`和{doc}`06-03-kkt-conditions`的例题与习题。

各题按主题分组：非线性规划问题的建模与分类、局部最优与全局最优、可行方向与下降方向、一阶与二阶条件、约束规范，以及 Karush–Kuhn–Tucker 条件（KKT 条件）。在写出 KKT 条件之前，每个问题都先化为 {eq}`eq-06-01-nlp` 的形式 $\min f$, $g_i\le0$, $h_j=0$。每道题的解答都先说明所用的结论，再手工完成计算并检验答案；所有最优解还用通用求解器（`scipy.optimize.minimize`，见{doc}`09-06-code-nlp`）做了数值验证。

(sec-08-08-formulation)=
## 建模与分类

:::{prf:example} 三角形上的凸二次规划
:label: ex-08-08-convex-qp

考虑

$$
\begin{aligned}
\min\ & f(\vx)=2x_1+x_2+x_1^2+2x_1x_2+x_2^2\\
\st\ & x_1+x_2\le2,\\
& x_1,\ x_2\ge0 .
\end{aligned}
$$

(a) 将 $f$ 写成 $\vc^\T\vx+\tfrac12\vx^\T\mQ\vx$ 的形式（$\mQ$ 对称），并证明该问题是凸规划（convex program）。
(b) 求最优解，并用 KKT 条件加以验证。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** 由线性部分得 $\vc=(2,1)^\T$。由二次部分 $x_1^2+2x_1x_2+x_2^2$ 得 $Q_{11}=2\cdot1$、$Q_{22}=2\cdot1$ 以及 $Q_{12}=Q_{21}=2$：

$$
\begin{aligned}
\mQ&=\begin{pmatrix}2&2\\2&2\end{pmatrix},\\
\lambda_{1,2}(\mQ)&=0,\ 4 .
\end{aligned}
$$

两个特征值都非负，故 $\mQ\succeq\vzero$，$f$ 是凸函数（但不是严格凸的，因为 $\mQ$ 奇异）。可行集是多面体，因而是凸集。所以该问题是凸二次规划（quadratic programming, QP；{ref}`sec-06-01-classes`）。

**(b)** 将 $f$ 写成 $f(\vx)=(x_1+x_2)^2+2x_1+x_2$。在可行集上每一项都非负，故 $f(\vx)\ge0=f(0,0)$，且 $f(\vx)=0$ 迫使 $x_1=x_2=0$。唯一的最优解是 $\vx^*=(0,0)$，$f^*=0$。

用 $g_1=x_1+x_2-2$、$g_2=-x_1$、$g_3=-x_2$ 做 KKT 检验：在原点处 $g_1=-2<0$，故 $\lambda_1=0$，而 $g_2$ 与 $g_3$ 是有效约束（active constraint）。梯度为 $\nabla f(\vx)=(2+2x_1+2x_2,\ 1+2x_1+2x_2)^\T$，因此 $\nabla f(\vzero)=(2,1)^\T$，由平稳性条件（stationarity）

$$
\begin{pmatrix}2\\1\end{pmatrix}+\lambda_2\begin{pmatrix}-1\\0\end{pmatrix}+\lambda_3\begin{pmatrix}0\\-1\end{pmatrix}=\vzero
$$

得 $\lambda_2=2\ge0$，$\lambda_3=1\ge0$。KKT 条件成立，由{prf:ref}`thm-06-03-convex-kkt`，原点是全局极小点（global minimum）。注意，这里的最优解恰好是顶点，但这只是因为 $\nabla f(\vzero)$ 指向可行集内部；起决定作用的是梯度，而不是顶点结构。（作为对比，凸函数 $f$ 在该三角形上的*最大值*在顶点 $(2,0)$ 处取得，$f=8$；凸函数在多胞形上的最大值总在某个顶点处取得。）
:::

(sec-08-08-local-global)=
## 局部最优与全局最优

:::{prf:example} 达不到的上确界与不连通的可行集
:label: ex-08-08-difficulties

(a) 设 $f(x)=1-1/x$，$\Omega=\{x:x\ge0.5\}$。$\max_{x\in\Omega}f(x)$ 是多少？{prf:ref}`thm-06-01-existence`的哪一条假设不成立？

(b) 对 $x>0$，令 $g_1(x)=(x-2)^2/x$，$g_2(x)=2.5-0.3x$。描述集合 $\Omega=\{x>0:g_1(x)\ge g_2(x)\}$，并说明它为什么会给搜索方法带来困难。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** $f'(x)=1/x^2>0$，故 $f$ 在 $\Omega$ 上单调递增。对每个 $x\in\Omega$ 都有 $f(x)<1$，且当 $x\to\infty$ 时 $f(x)\to1$（例如 $f(100)=0.99$）。因此 $\sup_{\Omega}f=1$，但没有任何可行点能达到它：**该问题没有最大值。** 最小值存在，在端点 $x=0.5$ 处取得，$f(0.5)=-1$。魏尔斯特拉斯定理（Weierstrass theorem）不适用，因为 $\Omega$ 是闭集但无界，且 $-f$ 在 $\Omega$ 上不是强制函数（$-f\to-1$，而不是 $+\infty$）。仅凭目标函数值有界，并不能保证最优解存在。

**(b)** 将 $g_1\ge g_2$ 两边乘以 $x>0$：

$$
(x-2)^2\ge x(2.5-0.3x)\iff 1.3x^2-6.5x+4\ge0 .
$$

$1.3x^2-6.5x+4$ 的根为 $x=\bigl(6.5\mp\sqrt{6.5^2-4\cdot1.3\cdot4}\bigr)/2.6=(6.5\mp\sqrt{21.45})/2.6$，即 $x\approx0.7187$ 与 $x\approx4.2813$。该二次函数在两根之外非负，因此

$$
\Omega=(0,\ 0.7187]\ \cup\ [4.2813,\ \infty),
$$

是两个互不相连的区间。（在 $x=2$ 处检验：$g_1=0<g_2=1.9$，故 $2\notin\Omega$。）像{doc}`06-05-unconstrained-algorithms`和{doc}`06-06-constrained-algorithms`中的下降法那样在可行点之间连续移动的方法，会停留在出发点所在的那一段，永远看不到另一段。要找到全局最优解，必须分别搜索每个连通部分；而在高维情形下，这些连通部分甚至无法事先知道。{numref}`fig-08-08-difficulties`展示了这两种情形。

```{figure} figures/08-08-difficulties.svg
:name: fig-08-08-difficulties
:width: 680px
:class: dark-light
:alt: 左图：x >= 0.5 上递增的曲线 1 - 1/x 逼近高度为 1 的虚线渐近线；右图：一条凸曲线与一条递减直线相交两次，x 轴上曲线位于直线上方的两段被突出显示

左图：$x\ge0.5$ 上的 $f(x)=1-1/x$ 逼近 $1$ 但永远达不到（{prf:ref}`ex-08-08-difficulties` (a)）。右图：$\{x:g_1(x)\ge g_2(x)\}$ 由突出显示的两个区间组成（{prf:ref}`ex-08-08-difficulties` (b)）。
```
:::

(sec-08-08-directions)=
## 可行方向与下降方向

:::{prf:example} 顶点处的可行下降方向
:label: ex-08-08-descent-cone

考虑

$$
\begin{aligned}
\min\ & f(\vx)=(x_1-6)^2+(x_2-7)^2\\
\st\ & x_1+x_2\le5,\\
& x_1-x_2\le1,\\
& x_1\ge0 .
\end{aligned}
$$

(a) 给出可行点 $(1,1)$ 处的两个下降方向（descent direction）。
(b) 求可行点 $(3,2)$ 处的三个可行下降方向（feasible descent direction），并说明它们为何符合要求。
(c) 求最优解及其乘子。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

可行集是以 $(0,5)$、$(3,2)$、$(0,-1)$ 为顶点的三角形，且 $\nabla f(\vx)=(2(x_1-6),\ 2(x_2-7))^\T$。

**(a)** 在 $(1,1)$ 处：$1+1=2<5$，$1-1=0<1$，$1>0$，故该点是内点，任何方向都是可行方向（feasible direction）（{prf:ref}`prop-06-02-feasible-directions`）。由 $\nabla f(1,1)=(-10,-12)^\T$，当 $-10d_1-12d_2<0$ 时方向是下降方向。例如

$$
\begin{aligned}
\vd&=(1,2):&\nabla f^\T\vd&=-10-24=-34<0,\\
\vd&=(1,0):&\nabla f^\T\vd&=-10<0 .
\end{aligned}
$$

两者都是可行下降方向。

**(b)** 在 $(3,2)$ 处，约束 $x_1+x_2\le5$ 与 $x_1-x_2\le1$ 是有效约束，而 $x_1\ge0$ 不是。所有约束都是线性的，故由{prf:ref}`prop-06-02-feasible-directions`，可行方向为

$$
\begin{aligned}
&d_1+d_2\le0,\\
&d_1-d_2\le0
\end{aligned}
\quad\Longleftrightarrow\quad d_1\le-\abs{d_2} .
$$

梯度为 $\nabla f(3,2)=(-6,-10)^\T$，当 $-6d_1-10d_2<0$，即 $d_2>-0.6\,d_1$ 时，$\vd$ 是下降方向。由于可行性迫使 $d_1\le0$，两个条件合起来为

$$
0.6\abs{d_1}<d_2\le\abs{d_1},
$$

因此，在相差一个正倍数的意义下，可行下降方向为 $\vd=(-1,s)$，$0.6<s\le1$。取其中三个：

```{list-table} $(3,2)$ 处的可行下降方向。
:header-rows: 1

* - $\vd$
  - $d_1+d_2$
  - $d_1-d_2$
  - $\nabla f(3,2)^\T\vd$
* - $(-3,3)$
  - $0$
  - $-6$
  - $18-30=-12$
* - $(-3,2)$
  - $-1$
  - $-5$
  - $18-20=-2$
* - $(-5,4)$
  - $-1$
  - $-9$
  - $30-40=-10$
```

第一个方向沿边 $x_1+x_2=5$ 指向 $(0,5)$；另外两个指向三角形内部。这个锥很窄（{numref}`fig-08-08-descent-cone`）：它位于边方向 $(-1,1)$ 与同梯度正交的方向 $(-1,0.6)$ 之间。*（本解答的早期版本还列出了 $(-2,2)$，但它与 $(-3,3)$ 是同一个方向；该版本写下降条件时下标有误，把第二个分量写成了 $d_{22}$；并且用 $-3/\alpha$ 去界定第一个分量，而约束 $x_1\ge0$ 在 $(3,2)$ 处是非有效约束，并不限制方向。）*

**(c)** 由于 $(1,1)$ 和 $(3,2)$ 处都存在可行下降方向，两者都不是最优解。目标函数是到点 $(6,7)$ 的距离的平方，而该点位于边 $x_1+x_2=5$ 的另一侧；将它投影到这条直线上得 $(6,7)-\tfrac{13-5}{2}(1,1)=(2,3)$，它满足 $2-3=-1\le1$ 与 $2\ge0$。KKT 检验：只有 $g_1=x_1+x_2-5$ 是有效约束，且 $\nabla f(2,3)=(-8,-8)^\T=-8\nabla g_1$，故 $\lambda_1=8\ge0$，$\lambda_2=\lambda_3=0$。该问题是凸的，所以 $\vx^*=(2,3)$，$f^*=16+16=32$ 是全局极小点（{prf:ref}`thm-06-03-convex-kkt`）。

```{figure} figures/08-08-descent-cone.svg
:name: fig-08-08-descent-cone
:width: 520px
:class: dark-light
:alt: 以 (0,5)、(3,2)、(0,-1) 为顶点的三角形；(3,2) 处的梯度指向左下方，负梯度指向右上方；(3,2) 处一个狭窄的阴影锥表示可行下降方向，它介于指向 (0,5) 的边与同梯度正交的直线之间

{prf:ref}`ex-08-08-descent-cone` 中顶点 $(3,2)$ 处的可行下降方向。可行方向指向三角形内部；下降方向与 $\nabla f(3,2)=(-6,-10)$ 成钝角。两者的交是介于边方向 $(-1,1)$ 与同梯度正交的方向 $(-1,0.6)$（点线）之间的狭窄阴影锥。最优解 $\vx^*=(2,3)$ 位于这条边上。
```
:::

(sec-08-08-first-second)=
## 一阶与二阶条件

:::{prf:example} 三元二次函数
:label: ex-08-08-quadratic-three

设

$$
f(\vx)=3x_1^2+2x_1x_2+6x_2^2-x_2x_3+5x_3^2-4x_1-2x_2-6x_3+9 .
$$

(a) 求 $f$ 的驻点（stationary point）。(b) 证明它是局部极小点（local minimum）。(c) 证明它是全局极小点，并给出最小值。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** 令梯度为零，

$$
\nabla f(\vx)=\begin{pmatrix}6x_1+2x_2-4\\2x_1+12x_2-x_3-2\\-x_2+10x_3-6\end{pmatrix}=\vzero .
$$

由第一个和第三个方程得 $x_1=(2-x_2)/3$ 与 $x_3=(6+x_2)/10$。代入第二个方程并乘以 30：

$$
20(2-x_2)+360x_2-3(6+x_2)-60=0\iff337x_2=38 .
$$

因此

$$
\vx^*=\Bigl(\tfrac{212}{337},\ \tfrac{38}{337},\ \tfrac{206}{337}\Bigr)\approx(0.6291,\ 0.1128,\ 0.6113).
$$

**(b)** 黑塞矩阵（Hessian matrix）为常数矩阵，

$$
\nabla^2f=\begin{pmatrix}6&2&0\\2&12&-1\\0&-1&10\end{pmatrix},
$$

其顺序主子式为 $6>0$、$6\cdot12-2\cdot2=68>0$ 以及 $\det\nabla^2f=6(120-1)-2(20-0)+0=674>0$。顺序主子式全为正的对称矩阵是正定的（{prf:ref}`thm-01-01-sylvester`）；其特征值为 $5.376$、$9.704$ 和 $12.92$。由{prf:ref}`thm-06-02-second-order-sufficient`，$\vx^*$ 是严格局部极小点。

**(c)** 由于处处有 $\nabla^2f\succ\vzero$，$f$ 在 $\R^3$ 上严格凸，由{prf:ref}`thm-06-02-convex-sufficient`，该驻点是唯一的全局极小点。记 $f(\vx)=\tfrac12\vx^\T\mH\vx-\vb^\T\vx+9$，其中 $\mH=\nabla^2f$，$\vb=(4,2,6)^\T$，且 $\mH\vx^*=\vb$，则最小值为

$$
\begin{aligned}
f^*&=9-\tfrac12\vb^\T\vx^*=9-\frac{4\cdot212+2\cdot38+6\cdot206}{2\cdot337}\\
&=9-\frac{1080}{337}=\frac{1953}{337}\approx5.795 .
\end{aligned}
$$

从原点出发的拟牛顿求解器给出相同的点和函数值。
:::

:::{prf:example} 蛋托曲面
:label: ex-08-08-egg-crate

$[-4,4]^2$ 上的曲面 $f(\vx)=\sin x_1\sin x_2$ 凸起与凹陷交替出现。证明原点是驻点但不是局部极小点，并对该正方形内的所有驻点进行分类。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

梯度与黑塞矩阵为

$$
\begin{aligned}
\nabla f(\vx)&=\begin{pmatrix}\cos x_1\sin x_2\\ \sin x_1\cos x_2\end{pmatrix},\\
\nabla^2f(\vx)&=\begin{pmatrix}-\sin x_1\sin x_2&\cos x_1\cos x_2\\ \cos x_1\cos x_2&-\sin x_1\sin x_2\end{pmatrix}.
\end{aligned}
$$

**原点。** $\nabla f(\vzero)=\vzero$，而 $\nabla^2f(\vzero)=\bigl(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\bigr)$ 的特征值为 $\pm1$：它是不定的，故二阶必要条件（second-order necessary condition）不成立（{prf:ref}`thm-06-02-second-order-necessary`），原点是鞍点。也可以直接看出：对较小的 $t\ne0$，$f(t,t)=\sin^2t>0$，$f(t,-t)=-\sin^2t<0$，而 $f(\vzero)=0$。仅凭一阶条件无法把这个点与极小点区分开来。

**全部驻点。** $\cos x_1\sin x_2=0$ 且 $\sin x_1\cos x_2=0$。若 $\cos x_1=0$，则 $\sin x_1=\pm1$，第二个方程迫使 $\cos x_2=0$；若 $\sin x_2=0$，则 $\cos x_2=\pm1$，从而 $\sin x_1=0$。在 $[-4,4]$ 中，$\cos$ 的零点为 $\pm\pi/2$，$\sin$ 的零点为 $0,\pm\pi$（$\pm3\pi/2\approx\pm4.71$ 在区间之外）。所以共有 $4+9=13$ 个驻点：

- $(\pm\pi/2,\pm\pi/2)$：$\nabla^2f=-\sin x_1\sin x_2\,\mI$。在 $\sin x_1\sin x_2=1$ 处，即 $(\pi/2,\pi/2)$ 和 $(-\pi/2,-\pi/2)$，黑塞矩阵为 $-\mI\prec\vzero$：它们是严格局部极大点，$f=1$。在 $\sin x_1\sin x_2=-1$ 处，黑塞矩阵为 $\mI\succ\vzero$：它们是严格局部极小点，$f=-1$。
- $x_1,x_2\in\{-\pi,0,\pi\}$ 的九个点：黑塞矩阵的对角元为零，非对角元为 $\pm1$，故特征值为 $\pm1$：这九个点都是鞍点。
:::

:::{prf:example} 边界处的一阶条件
:label: ex-08-08-bound-first-order

在约束 $0\le x\le1$ 下极小化 $f(x)=(x-2)^2$，并证明一阶必要条件（first-order necessary condition）在极小点处成立。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

$f$ 的驻点 $x=2$ 不可行。在 $[0,1]$ 上 $f'(x)=2(x-2)<0$，故 $f$ 单调递减，极小点在右端点：$x^*=1$，$f^*=1$（比较 $f(0)=4$）。

在 $x^*=1$ 处，可行方向为 $d\le0$（只能向左移动），并且

$$
f'(x^*)\,d=-2d\ge0\qquad\text{对所有 }d\le0,
$$

所以虽然 $f'(1)\ne0$，{eq}`eq-06-02-first-order` 仍然成立。写成 KKT 形式，取 $g_1=-x\le0$（非有效约束，$\lambda_1=0$）与 $g_2=x-1\le0$（有效约束）：由 $f'(1)+\lambda_2=0$ 得 $\lambda_2=2\ge0$。该问题是凸的，因此这个 KKT 点（KKT point）就是全局极小点。在另一个端点处条件不成立：在 $x=0$ 处可行方向为 $d\ge0$，而对 $d>0$ 有 $f'(0)d=-4d<0$，所以 $0$ 不是极小点。
:::

(sec-08-08-cq)=
## 约束规范

:::{prf:example} 不存在 KKT 乘子的极小点
:label: ex-08-08-cusp

考虑

$$
\begin{aligned}
\min\ & f(\vx)=-x_1\\
\st\ & g_1(\vx)=x_2-(1-x_1)^3\le0,\\
& g_2(\vx)=-x_2\le0 .
\end{aligned}
$$

(a) 证明 $\vx^*=(1,0)$ 是唯一的最优解。
(b) 证明在 $\vx^*$ 处不存在 KKT 乘子。
(c) 哪个约束规范（constraint qualification）不成立？比较 $L(\vx^*)$ 与 $T(\vx^*)$。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** 可行性要求 $0\le x_2\le(1-x_1)^3$，因此 $(1-x_1)^3\ge0$，$x_1\le1$。于是在 $\Omega$ 上 $f=-x_1\ge-1$，等号仅在 $x_1=1$ 时成立，而此时 $0\le x_2\le0$ 迫使 $x_2=0$。唯一的最优解是 $\vx^*=(1,0)$，$f^*=-1$。可行集是一个指向右方的尖角，其尖端位于 $\vx^*$。

**(b)** 两个约束在 $\vx^*$ 处都是有效约束。梯度为 $\nabla f=(-1,0)^\T$、$\nabla g_1(\vx)=(3(1-x_1)^2,\ 1)^\T$，故 $\nabla g_1(\vx^*)=(0,1)^\T$，以及 $\nabla g_2=(0,-1)^\T$。平稳性条件

$$
\begin{pmatrix}-1\\0\end{pmatrix}+\lambda_1\begin{pmatrix}0\\1\end{pmatrix}+\lambda_2\begin{pmatrix}0\\-1\end{pmatrix}=\vzero
$$

要求第一个分量满足 $-1=0$：这不可能。全局极小点不是 KKT 点。（从 $(0.5,0.05)$ 出发运行数值求解器，它会逼近 $(1,0)$，但通常伴有警告，因为它的终止检验是基于 KKT 条件的。）

**(c)** 有效约束的梯度 $(0,1)$ 与 $(0,-1)$ 线性相关，故线性无关约束规范（linear independence constraint qualification, LICQ）不成立；约束不是线性的，$g_1$ 也不是凸函数，因此{prf:ref}`thm-06-03-cq-sufficient`中其他两个充分条件也都不适用。线性化可行方向锥（cone of linearized feasible directions）为

$$
L(\vx^*)=\{\vd:\ d_2\le0,\ -d_2\le0\}=\{(d_1,0):d_1\in\R\},
$$

即整个 $x_1$ 轴，而 $\vx^*$ 附近的可行集位于 $\vx^*$ 的左侧，故 $T(\vx^*)=\{(d_1,0):d_1\le0\}$。方向 $\vd=(1,0)$ 属于 $L(\vx^*)$ 但不属于 $T(\vx^*)$，且 $\nabla f(\vx^*)^\T\vd=-1<0$：线性化给出了一个实际上并不存在的改进方向，于是 Farkas 引理否定了乘子的存在。Fritz John 条件（{prf:ref}`rem-06-03-history`）则成立，其中 $\lambda_0=0$，$\lambda_1=\lambda_2=1$。这个例子出自 Kuhn 和 Tucker（1951）。
:::

(sec-08-08-kkt)=
## KKT 条件

:::{prf:example} 写出 KKT 方程组
:label: ex-08-08-kkt-formulation

(a) 写出下列问题的 KKT 条件：

$$
\begin{aligned}
\min\ & f(\vx)=4x_1^2+4x_1x_2+2x_2^2+3x_1+e^{x_1+2x_2}\\
\st\ & x_1^2+x_2^2\le10,\\
& -e^{-x_1^2}\le0.5,\\
& x_1\ge0 .
\end{aligned}
$$

(b) 求解这些条件。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** 化为 $g_i\le0$ 的形式：$g_1=x_1^2+x_2^2-10$，$g_2=-e^{-x_1^2}-0.5$，$g_3=-x_1$。梯度为

$$
\begin{aligned}
\nabla f&=\begin{pmatrix}8x_1+4x_2+3+e^{x_1+2x_2}\\4x_1+4x_2+2e^{x_1+2x_2}\end{pmatrix},\\
\nabla g_1&=\begin{pmatrix}2x_1\\2x_2\end{pmatrix},\\
\nabla g_2&=\begin{pmatrix}2x_1e^{-x_1^2}\\0\end{pmatrix},\\
\nabla g_3&=\begin{pmatrix}-1\\0\end{pmatrix}.
\end{aligned}
$$

KKT 条件为

$$
\begin{aligned}
8x_1+4x_2+3+e^{x_1+2x_2}+2\lambda_1x_1+2\lambda_2x_1e^{-x_1^2}-\lambda_3&=0,\\
4x_1+4x_2+2e^{x_1+2x_2}+2\lambda_1x_2&=0,
\end{aligned}
$$

$$
\begin{aligned}
\lambda_1(x_1^2+x_2^2-10)&=0,\\
\lambda_2\bigl(-e^{-x_1^2}-0.5\bigr)&=0,\\
\lambda_3x_1&=0,
\end{aligned}
$$

$$
\begin{aligned}
x_1^2+x_2^2&\le10,\\
-e^{-x_1^2}&\le0.5,\\
x_1&\ge0,\\
\lambda_1,\lambda_2,\lambda_3&\ge0 .
\end{aligned}
$$

由于 $-e^{-x_1^2}\in[-1,0)$，第二个约束对每个 $x_1$ 都严格成立：它永远不是有效约束，$\lambda_2=0$，可以去掉。*（本解答的早期版本把第二个平稳性方程中的 $2e^{x_1+2x_2}$ 印成了 $e^{x_1+2x_2}$，并把 $\lambda_1$ 印成了 $\lambda_2$；$e^{x_1+2x_2}$ 对 $x_2$ 的偏导数带有因子 2，而 $x_2$ 只出现在 $g_1$ 中。）*

**(b)** *凸性。* $f$ 的二次部分的黑塞矩阵为 $\bigl(\begin{smallmatrix}8&4\\4&4\end{smallmatrix}\bigr)$，特征值为 $1.528$ 与 $10.47$，而 $e^{x_1+2x_2}$ 是仿射函数的指数，故 $f$ 严格凸。去掉 $g_2$ 后各约束都是凸的，且 $(1,0)$ 是 Slater 点。因此任何 KKT 点都是唯一的全局极小点，求解时找到第一个 KKT 点即可停止。圆盘 $x_1^2+x_2^2\le10$ 很大，所以先尝试 $g_1$ 非有效（$\lambda_1=0$）的情形。

*情形 $x_1>0$，$\lambda_1=0$。* 此时 $\lambda_3=0$，记 $E=e^{x_1+2x_2}>0$，两个方程为 $8x_1+4x_2+3+E=0$ 与 $4x_1+4x_2+2E=0$。由第二个方程得 $x_2=-x_1-E/2$；代入得 $4x_1+3-E=0$，故 $E=4x_1+3$，$x_2=-3x_1-1.5$。但这样一来 $E=e^{x_1+2x_2}=e^{-5x_1-3}\le e^{-3}\approx0.0498$，与 $E=4x_1+3\ge3$ 矛盾。此类情形没有 KKT 点。

*情形 $x_1=0$，$\lambda_1=0$。* 第二个方程变为 $4x_2+2e^{2x_2}=0$，即 $e^{2x_2}=-2x_2$。令 $\zeta=-2x_2$，则方程为 $\zeta e^{\zeta}=1$，其解为欧米加常数（omega constant）$\zeta\approx0.567143$；故 $x_2=-\zeta/2\approx-0.28357$。由第一个方程得

$$
\begin{aligned}
\lambda_3&=4x_2+3+e^{2x_2}=4x_2+3-2x_2\\
&=3+2x_2=3-\zeta\approx2.4329>0,
\end{aligned}
$$

并且 $g_1=x_2^2-10<0$ 与 $\lambda_1=0$ 相容。所有 KKT 条件都成立：

$$
\begin{aligned}
\vx^*&\approx(0,\ -0.2836),\\
\vlambda^*&\approx(0,\ 0,\ 2.4329),\\
f^*&=2x_2^2+e^{2x_2}=\tfrac12\zeta^2+\zeta\approx0.7280 .
\end{aligned}
$$

界约束 $x_1\ge0$ 是唯一的紧约束：没有它，$f$ 在 $x_1$ 取负值时还会继续减小。由于问题是凸的且已经找到一个 KKT 点，$g_1$ 有效（$x_1^2+x_2^2=10$）的情形无需再考察：它们不可能产生第二个极小点。数值求解器验证了该点及其函数值。
:::

:::{prf:example} 最大化问题：求解与唯一性
:label: ex-08-08-kkt-max

考虑

$$
\begin{aligned}
\max\ & 14x_1-x_1^2+6x_2-x_2^2\\
\st\ & x_1+x_2\le2,\\
& x_1+2x_2\le3 .
\end{aligned}
$$

(a) 用 KKT 条件求最优解。
(b) 证明最优解唯一。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** 化为最小化问题：$\min f(\vx)=x_1^2-14x_1+x_2^2-6x_2$，约束为 $g_1=x_1+x_2-2\le0$，$g_2=x_1+2x_2-3\le0$。KKT 条件为

$$
\begin{aligned}
2x_1-14+\lambda_1+\lambda_2&=0,\\
2x_2-6+\lambda_1+2\lambda_2&=0,
\end{aligned}
$$

以及 $\lambda_1g_1=0$、$\lambda_2g_2=0$、$g_1,g_2\le0$ 和 $\lambda_1,\lambda_2\ge0$。

无约束极大点 $(7,3)$ 违反两个约束，故至少有一个约束是有效的。三种情形：

- $S=\{1\}$（$\lambda_2=0$）：$x_1=7-\lambda_1/2$，$x_2=3-\lambda_1/2$，由 $x_1+x_2=10-\lambda_1=2$ 得 $\lambda_1=8$，$\vx=(3,-1)$。检验：$g_2=3-2-3=-2\le0$，$\lambda_1=8\ge0$。**是 KKT 点。**
- $S=\{2\}$（$\lambda_1=0$）：$x_1=7-\lambda_2/2$，$x_2=3-\lambda_2$，由 $x_1+2x_2=13-2.5\lambda_2=3$ 得 $\lambda_2=4$，$\vx=(5,-1)$；但 $g_1=5-1-2=2>0$。不可行。
- $S=\{1,2\}$：$\vx=(1,1)$；此时由 $\lambda_1+\lambda_2=12$ 与 $\lambda_1+2\lambda_2=4$ 得 $\lambda_2=-8<0$。舍去。

唯一的 KKT 点是 $\vx^*=(3,-1)$，$\vlambda^*=(8,0)$，最大值为 $14\cdot3-9+6\cdot(-1)-1=26$。由于 $f$ 是凸函数且约束是线性的，它就是全局最优解（{prf:ref}`thm-06-03-convex-kkt`）。在原来的最大化形式中，$\nabla(14x_1-x_1^2+6x_2-x_2^2)=(8,8)^\T=8\nabla g_1$：目标函数的梯度沿有效约束的法向指向可行集外部，正如{numref}`tab-06-03-signs`第三行所述。

**(b)** 设 $\vx^*$ 与 $\vy^*$ 都是最优解。由凸可行集上的一阶条件（{eq}`eq-06-02-variational-inequality`），

$$
\begin{aligned}
\nabla f(\vx^*)^\T(\vy^*-\vx^*)&\ge0,\\
\nabla f(\vy^*)^\T(\vx^*-\vy^*)&\ge0 .
\end{aligned}
$$

两式相加得 $\bigl(\nabla f(\vx^*)-\nabla f(\vy^*)\bigr)^\T(\vx^*-\vy^*)\le0$。对这里的 $f$，$\nabla f(\vx)-\nabla f(\vy)=2(\vx-\vy)$，故该不等式即 $2\norm{\vx^*-\vy^*}^2\le0$，从而 $\vx^*=\vy^*$。（对一般的 $f$，由中值定理，左端可写成 $(\vx^*-\vy^*)^\T\nabla^2f(\vxi)(\vx^*-\vy^*)$，其中 $\vxi$ 是线段上某一点，再由黑塞矩阵的正定性可得同样的结论；这正是{prf:ref}`thm-06-01-convex-local-global`中的严格凸性论证。）
:::

:::{prf:example} 线性规划的 KKT 条件
:label: ex-08-08-kkt-lp

考虑线性规划（linear program, LP）$\min\ \vc^\T\vx$，满足 $\mA\vx=\vb$，$\vx\ge\vzero$，其中 $\mA\in\R^{m\times n}$。

(a) 写出它的 KKT 条件。
(b) 证明这些条件就是线性规划的最优性条件，特别是互补松弛性（complementary slackness）。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** 约束为 $h_i(\vx)=\mA_{i\cdot}\vx-b_i=0$（$i=1,\dots,m$，乘子 $\mu_i$ 符号不限）与 $g_j(\vx)=-x_j\le0$（$j=1,\dots,n$，乘子 $\lambda_j\ge0$）。由 $\nabla h_i=\mA_{i\cdot}^\T$ 与 $\nabla g_j=-\ve_j$，KKT 条件为

$$
\begin{aligned}
&\vc+\mA^\T\vmu-\vlambda=\vzero,\\
&\mA\vx=\vb,\\
&\vx\ge\vzero,\\
&\vlambda\ge\vzero,\\
&\lambda_jx_j=0\quad(j=1,\dots,n).
\end{aligned}
$$

由于约束是线性的，这些条件在每个最优解处都是必要的（{prf:ref}`thm-06-03-cq-sufficient`）；又由于线性规划是凸的，它们也是充分的（{prf:ref}`thm-06-03-convex-kkt`）。

**(b)** 令 $\vpi=-\vmu$。由平稳性条件得 $\vlambda=\vc-\mA^\T\vpi$，即检验数（reduced cost）$\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ 构成的向量。KKT 条件变为

$$
\begin{aligned}
&\mA\vx=\vb,\ \vx\ge\vzero && \text{（原始可行性）},\\
&\mA^\T\vpi\le\vc && \text{（对偶可行性）},\\
&x_j\,(c_j-\vpi^\T\mA_{\cdot j})=0,\ j=1,\dots,n && \text{（互补松弛性）},
\end{aligned}
$$

其中 $\mA^\T\vpi\le\vc$ 是对偶线性规划 $\max\{\vb^\T\vpi:\mA^\T\vpi\le\vc\}$ 的约束。这些恰好就是{doc}`03-05-lp-duality`中的最优性条件（{prf:ref}`thm-03-05-complementary-slackness`）：一个变量只有在其检验数为零时才可能取正值，而等式约束的乘子就是对偶变量（单纯形乘子）。将平稳性条件乘以 $\vx$ 还可以看出对偶间隙（duality gap）为零：$\vc^\T\vx=\vpi^\T\mA\vx+\vlambda^\T\vx=\vb^\T\vpi+0$。

*数值验证。* 对于 $\min\,-3x_1-5x_2$，满足 $x_1+x_3=4$，$2x_2+x_4=12$，$3x_1+2x_2+x_5=18$，$\vx\ge\vzero$，线性规划求解器给出 $\vx^*=(2,6,2,0,0)$ 与 $\vpi=(0,-1.5,-1)$。检验数为 $\vc-\mA^\T\vpi=(0,0,0,1.5,1)\ge\vzero$，每个乘积 $x_j\bar c_j$ 都为零，且 $\vc^\T\vx^*=-36=\vb^\T\vpi$。*（本解答的早期版本把可行性条件而不是 $\lambda_jx_j=0$ 称为互补松弛条件。）*
:::

:::{prf:example} 作为灵敏度的乘子
:label: ex-08-08-multiplier-sensitivity

考虑 $\min f(\vx)=x_1^2+2x_2^2+2x_1x_2-2x_1-4x_2$，约束为 $-x_1-x_2+2\le0$。

(a) 用 KKT 条件求解该问题，并求乘子 $\lambda^*$。
(b) 将约束放松为 $-x_1-x_2+2\le\delta$（$\delta>0$ 很小），重新求解，并证明 $\lambda^*$ 度量了最优值的变化率。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** 由 $g=-x_1-x_2+2$ 与 $\nabla g=(-1,-1)^\T$，平稳性条件为

$$
\begin{aligned}
2x_1+2x_2-2&=\lambda,\\
2x_1+4x_2-4&=\lambda .
\end{aligned}
$$

无约束驻点（$\lambda=0$）是 $(0,1)$，在该点 $g=1>0$：不可行。所以约束是有效的，$x_1+x_2=2$，由第一个方程得 $\lambda=2\cdot2-2=2$。第二个方程变为 $x_1+2x_2=3$，故 $x_2=1$，$x_1=1$。于是

$$
\begin{aligned}
\vx^*&=(1,1),\\
\lambda^*&=2\ge0,\\
f^*&=1+2+2-2-4=-1 .
\end{aligned}
$$

黑塞矩阵 $\bigl(\begin{smallmatrix}2&2\\2&4\end{smallmatrix}\bigr)$ 的特征值为 $0.764$ 与 $5.236$，故 $f$ 严格凸，$\vx^*$ 是唯一的全局极小点。

**(b)** 放松后的约束为 $x_1+x_2\ge2-\delta$。对 $0\le\delta<1$，同样的推理给出 $x_1+x_2=2-\delta$，$\lambda(\delta)=2(2-\delta)-2=2-2\delta>0$，$x_1+2x_2=3-\delta$，因此

$$
\begin{aligned}
\vx^*(\delta)&=(1-\delta,\ 1),\\
f^*(\delta)&=(1-\delta)^2+2+2(1-\delta)-2(1-\delta)-4\\
&=\delta^2-2\delta-1 .
\end{aligned}
$$

于是

$$
\frac{\mathrm df^*}{\mathrm d\delta}\Big|_{\delta=0}=-2=-\lambda^*,
$$

这与 {eq}`eq-06-03-sensitivity` 一致：把约束放松 $\delta$，最优值约降低 $\lambda^*\delta$。例如，$\delta=0.1$ 时 $f^*=-1.19$，而线性预测值为 $-1-2(0.1)=-1.2$；差商 $\bigl(f^*(\delta)-f^*\bigr)/(-\delta)=2-\delta$ 趋于 $\lambda^*=2$。当 $\delta=1$ 时乘子降为 $0$，无约束极小点 $(0,1)$ 满足放松后的约束，再继续放松也无济于事。*（本解答的早期版本把有效约束写成 $-x_1-x_2=2$；应为 $x_1+x_2=2$。）*
:::

:::{prf:example} 没有有效约束的 KKT 点
:label: ex-08-08-interior-kkt

用 KKT 条件求解

$$
\begin{aligned}
\min\ & f(\vx)=4(x_1-2)^2+3(x_2-4)^2\\
\st\ & -x_1-x_2+5\le0,\\
& -x_1+1\le0,\\
& -x_2+2\le0,
\end{aligned}
$$

并证明解是唯一的。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

KKT 条件为

$$
\begin{aligned}
8(x_1-2)-\lambda_1-\lambda_2&=0,\\
6(x_2-4)-\lambda_1-\lambda_3&=0,
\end{aligned}
$$

以及 $\lambda_1(-x_1-x_2+5)=0$、$\lambda_2(1-x_1)=0$、$\lambda_3(2-x_2)=0$、可行性条件和 $\lambda_i\ge0$。

先尝试没有有效约束的情形 $\vlambda=\vzero$：此时 $x_1=2$，$x_2=4$。可行性：$-2-4+5=-1\le0$，$1-2=-1\le0$，$2-4=-2\le0$。所有条件都成立，故 $\vx^*=(2,4)$，$\vlambda^*=\vzero$，$f^*=0$。所有约束都是非有效的：无约束极小点恰好是可行的。

唯一性：$\nabla^2f=\bigl(\begin{smallmatrix}8&0\\0&6\end{smallmatrix}\bigr)\succ\vzero$，故 $f$ 严格凸；可行集是凸集，而严格凸函数在凸集上至多有一个极小点（{prf:ref}`thm-06-01-convex-local-global`）。*（本解答的早期版本给出的黑塞矩阵为 $\diag(12,6)$；由系数 $4(x_1-2)^2$ 得到的应是 $8$。）*
:::

## 参考文献

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 4, Fritz John
  and KKT conditions; Chapter 5, the cone of tangents and constraint
  qualifications.)
- John, F. (1948). Extremum problems with inequalities as subsidiary
  conditions. In *Studies and essays presented to R. Courant on his 60th
  birthday* (pp. 187–204). Interscience.
- Karush, W. (1939). *Minima of functions of several variables with
  inequalities as side constraints* [Master's thesis, University of Chicago].
- Kuhn, H. W., & Tucker, A. W. (1951). Nonlinear programming. In J. Neyman
  (Ed.), *Proceedings of the Second Berkeley Symposium on Mathematical
  Statistics and Probability* (pp. 481–492). University of California
  Press.
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming*
  (4th ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3>
  (Chapters 7 and 11, optimality conditions; Section 11.7, sensitivity.)
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Section 2.1 and
  Chapter 12, optimality conditions.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 20, classical optimization theory and the KKT
  conditions.)
