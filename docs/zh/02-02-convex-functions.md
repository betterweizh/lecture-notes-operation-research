# 凸函数

## 概述

当一个最优化问题的目标函数和不等式约束都是凸函数时，它在确切的意义下是“容易”的：此时每个局部极小点都是全局极小点，并且可以用高效的算法求出极小点（{doc}`02-03-convex-optimization-duality`）。本章定义凸函数、严格凸函数和凹函数（{prf:ref}`def-02-02-convex-function`），通过上图把它们与{doc}`02-01-convex-sets`中的凸集联系起来（{prf:ref}`thm-02-02-epigraph`），并导出两个实用的凸性判别法：一阶条件——凸函数位于其切平面之上（{prf:ref}`thm-02-02-first-order`）；二阶条件——凸函数的黑塞矩阵半正定（{prf:ref}`thm-02-02-second-order`）。随后列出常见的凸函数、保凸运算、下水平集的凸性、Jensen 不等式，并简要介绍拟凸函数。本章用到{doc}`01-01-linear-algebra`中的有定性判别法，以及{doc}`01-02-matrix-calculus`中的梯度、黑塞矩阵和泰勒定理。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 叙述凸函数、严格凸函数、凹函数和严格凹函数的定义，并解释为什么定义域必须是凸集；
- 利用定义、上图或一阶、二阶条件证明一个函数是凸函数；
- 通过黑塞矩阵的特征值或主子式判别二次函数的凸性；
- 利用保凸运算由简单的凸函数构造凸函数，并识别常见的例子（范数、指数函数、对数函数、最大值函数、log-sum-exp 函数）；
- 利用下水平集的凸性证明可行集是凸集；
- 应用 Jensen 不等式并识别拟凸函数。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\dom f$ | $f$ 的定义域（$f$ 取有限值的集合） |
| $\theta$, $\theta_i$ | 凸组合的权重，$\theta\in[0,1]$ |
| $\tilde f$ | $f$ 的扩展值延拓，在 $\dom f$ 之外取值 $+\infty$ |
| $\iota_C$ | 集合 $C$ 的示性函数：在 $C$ 上为 $0$，在 $C$ 之外为 $+\infty$ |
| $\operatorname{epi} f$ | 上图 $\{(\vx,t):\vx\in\dom f,\ f(\vx)\le t\}$ |
| $S_\gamma$ | 下水平集 $\{\vx\in\dom f: f(\vx)\le\gamma\}$ |
| $\operatorname{lse}(\vx)$ | log-sum-exp 函数 $\log\sum_{i}e^{x_i}$ |

本章中的对数均为自然对数。

(sec-02-02-definition)=
## 凸函数与凹函数

:::{prf:definition} 凸函数与凹函数
:label: def-02-02-convex-function

设 $C\subseteq\R^n$ 是非空凸集，$f:C\to\R$。若对所有 $\vx,\vy\in C$ 和所有 $\theta\in[0,1]$ 有

$$
f\bigl(\theta\vx+(1-\theta)\vy\bigr)\le\theta f(\vx)+(1-\theta)f(\vy).
$$ (eq-02-02-convex)

则称 $f$ 在 $C$ 上是**凸函数**（convex function）。进一步，

- 若当 $\vx\ne\vy$ 且 $\theta\in(0,1)$ 时不等式 {eq}`eq-02-02-convex` 严格成立，则称 $f$ 在 $C$ 上**严格凸**（strictly convex）；
- 若 $-f$ 是凸（严格凸）函数，即 {eq}`eq-02-02-convex` 在 $\ge$ 下成立（当 $\vx\ne\vy$、$\theta\in(0,1)$ 时在 $>$ 下成立），则称 $f$ 在 $C$ 上是**凹函数**（concave function）（**严格凹**，strictly concave）。
:::

*（严格性只能对不同的两点以及 $\theta\in(0,1)$ 提出要求：当 $\theta=0$、$\theta=1$ 或 $\vx=\vy$ 时，{eq}`eq-02-02-convex` 两边相等，因此若定义要求对所有 $\theta\in[0,1]$ 都有 $<$，就不会有任何函数是严格凸的。）*

{eq}`eq-02-02-convex` 的右端是从 $(\vx,f(\vx))$ 到 $(\vy,f(\vy))$ 的**弦**（chord）在点 $\theta\vx+(1-\theta)\vy$ 处的高度。因此，凸性意味着每条弦都位于图像之上或与图像重合（{numref}`fig-02-02-chord-tangent`(a)）；严格凸意味着弦在两个端点之间严格位于图像之上。凹函数的每条弦都位于图像之下或与图像重合。

定义域必须是凸集，因为 {eq}`eq-02-02-convex` 要在点 $\theta\vx+(1-\theta)\vy$ 处计算 $f$ 的值，而对任意两点和任意 $\theta$，这个点都必须属于定义域。因此，“$1/x$ 是凸函数”这样的说法是不完整的：$1/x$ 在 $(0,\infty)$ 上是凸的，在 $(-\infty,0)$ 上是凹的；而在它的整个定义域 $\R\setminus\{0\}$ 上，由于该集合不是凸集，定义不适用。

```{figure} figures/02-02-chord-tangent.svg
:name: fig-02-02-chord-tangent
:width: 640px
:class: dark-light
:alt: 左图：一条凸曲线，两点 x 与 y 之间的弦位于曲线之上。右图：同一条曲线，在 x 处的切线位于曲线之下，并标出在点 y 处的差距

一元凸函数的两种刻画。(a) $(x,f(x))$ 与 $(y,f(y))$ 之间的弦位于图像之上：在 $\theta x+(1-\theta)y$ 处，弦的高度 $\theta f(x)+(1-\theta)f(y)$ 不小于 $f(\theta x+(1-\theta)y)$。(b) $x$ 处的切线位于图像之下：对每个 $y$ 有 $f(y)\ge f(x)+f'(x)(y-x)$。
```

仅凭定义就可以回答以下三个问题。

*线性函数都是凸函数吗？* 是的。对仿射函数 $f(\vx)=\va^\T\vx+b$，{eq}`eq-02-02-convex` 两边相等，因为 $\va^\T(\theta\vx+(1-\theta)\vy)+b=\theta(\va^\T\vx+b)+(1-\theta)(\va^\T\vy+b)$。所以仿射函数既是凸函数又是凹函数，但既不是严格凸的，也不是严格凹的。反之，在 $\R^n$ 上既凸又凹的函数必是仿射函数。

*抛物线都是凸的吗？* 不是。对 $f(x)=ax^2+bx+c$，简单计算可得，对所有 $x,y$ 和 $\theta$ 有

$$
\theta f(x)+(1-\theta)f(y)-f\bigl(\theta x+(1-\theta)y\bigr)
=a\,\theta(1-\theta)(x-y)^2 .
$$ (eq-02-02-parabola-gap)

与任何仿射函数一样，一次项和常数项相互抵消；对二次项有 $\theta x^2+(1-\theta)y^2-(\theta x+(1-\theta)y)^2=\theta(1-\theta)(x-y)^2$。差值 {eq}`eq-02-02-parabola-gap` 对所有 $x,y,\theta$ 非负，当且仅当 $a\ge0$。所以抛物线是凸的当且仅当其首项系数非负，是严格凸的当且仅当 $a>0$，是凹的当且仅当 $a\le0$。抛物线 $-x^2$ 是严格凹的。

*一般情况下如何验证凸性？* 对所有点对逐一检验 {eq}`eq-02-02-convex` 很少可行。下面各节给出三种更好的工具：导数判别法（{ref}`sec-02-02-first-order`、{ref}`sec-02-02-second-order`）、已知凸函数的目录（{ref}`sec-02-02-examples`）以及组合它们的规则（{ref}`sec-02-02-operations`）。

一个有用的方法是化归为一元情形：函数 $f$ 是凸函数，当且仅当对所有 $\vx\in\dom f$ 和所有方向 $\vv$，它在每条直线上的限制（restriction to a line）$\varphi(t)=f(\vx+t\vv)$（定义在 $\{t:\vx+t\vv\in\dom f\}$ 上）都是凸函数。这是因为直线上两点的凸组合仍是同一线段上的点，所以 $f$ 满足 {eq}`eq-02-02-convex` 与所有限制 $\varphi$ 满足 {eq}`eq-02-02-convex` 是同一个命题。下面许多证明都先对 $n=1$ 进行，再用这种方式推广到 $\R^n$。

(sec-02-02-epigraph)=
## 扩展值与上图

在定义域之外令凸函数取值 $+\infty$，从而把它定义在整个 $\R^n$ 上，常常很方便。凸函数 $f$ 的**扩展值延拓**（extended-value extension）为

$$
\tilde f(\vx)=
\begin{cases}
f(\vx), & \vx\in\dom f,\\
+\infty, & \vx\notin\dom f .
\end{cases}
$$

约定运算规则：当 $\theta>0$ 时 $\theta\cdot\infty=\infty$，且对每个 $a$ 有 $a\le\infty$。于是不等式 {eq}`eq-02-02-convex` 对 $\tilde f$ 在所有 $\vx,\vy\in\R^n$ 和 $\theta\in(0,1)$ 处成立：若 $\vx$ 或 $\vy$ 在定义域之外，右端为 $+\infty$。定义域可以由 $\dom f=\{\vx:\tilde f(\vx)<\infty\}$ 还原。一个特例是凸集 $C$ 的**示性函数**（indicator function）$\iota_C$，它在 $C$ 上取 $0$，在 $C$ 之外取 $+\infty$。在 $C$ 上极小化 $f$ 等价于在 $\R^n$ 上极小化 $f+\iota_C$，因此约束可以写进目标函数。此后我们不加说明地使用这种延拓。

凸函数与凸集之间的桥梁是位于图像之上或图像上的点构成的集合。

:::{prf:definition} 上图
:label: def-02-02-epigraph

$f:\dom f\to\R$ 的**上图**（epigraph）是 $\R^{n+1}$ 的子集

$$
\operatorname{epi} f=\{(\vx,t):\vx\in\dom f,\ t\in\R,\ f(\vx)\le t\}.
$$
:::

:::{prf:theorem} 凸函数的上图是凸集
:label: thm-02-02-epigraph

定义域为凸集的函数 $f$ 是凸函数，当且仅当 $\operatorname{epi} f$ 是凸集。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

设 $f$ 是凸函数，$(\vx,s),(\vy,t)\in\operatorname{epi} f$，$\theta\in[0,1]$。则 $f(\theta\vx+(1-\theta)\vy)\le\theta f(\vx)+(1-\theta)f(\vy)\le\theta s+(1-\theta)t$，所以组合点 $(\theta\vx+(1-\theta)\vy,\ \theta s+(1-\theta)t)$ 属于上图。

反之，设上图是凸集，$\vx,\vy\in\dom f$。点 $(\vx,f(\vx))$ 和 $(\vy,f(\vy))$ 属于上图，从而 $(\theta\vx+(1-\theta)\vy,\ \theta f(\vx)+(1-\theta)f(\vy))$ 也属于上图，这正是 {eq}`eq-02-02-convex`。
:::

这个定理把关于函数的命题转化为关于集合的命题。例如，凸函数的逐点最大值是凸函数，因为它的上图是这些函数上图的交（{prf:ref}`thm-02-01-intersection`）；下一节的切平面性质就是把支撑超平面定理（{prf:ref}`thm-02-01-supporting-hyperplane`）应用于上图的结果。{numref}`fig-02-02-epigraph-sublevel`(a) 画出了一个上图。

(sec-02-02-first-order)=
## 一阶条件

此后所说的可微函数，是指其梯度在包含其（凸）定义域的某个开集上存在的函数；梯度 $\nabla f$ 和黑塞矩阵 $\nabla^2f$ 见{doc}`01-02-matrix-calculus`。

:::{prf:theorem} 一阶条件
:label: thm-02-02-first-order

设 $C$ 是凸集，$f$ 在 $C$ 上可微。则 $f$ 在 $C$ 上是凸函数，当且仅当

$$
f(\vy)\ge f(\vx)+\nabla f(\vx)^\T(\vy-\vx)\qquad\text{对所有 }\vx,\vy\in C .
$$ (eq-02-02-first-order)

$f$ 在 $C$ 上严格凸，当且仅当 {eq}`eq-02-02-first-order` 对所有 $\vx\ne\vy$ 以严格不等式成立。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

*凸 ⇒ {eq}`eq-02-02-first-order`。* 设 $\vx,\vy\in C$，$\theta\in(0,1]$。对点 $\vx+\theta(\vy-\vx) =\theta\vy+(1-\theta)\vx$ 利用凸性，得

$$
f\bigl(\vx+\theta(\vy-\vx)\bigr)\le f(\vx)+\theta\bigl(f(\vy)-f(\vx)\bigr).
$$

两边减去 $f(\vx)$ 再除以 $\theta$，得

$$
\frac{f\bigl(\vx+\theta(\vy-\vx)\bigr)-f(\vx)}{\theta}\le f(\vy)-f(\vx).
$$

当 $\theta\downarrow0$ 时，左端趋于方向导数 $\nabla f(\vx)^\T(\vy-\vx)$，由此得到 {eq}`eq-02-02-first-order`。

*{eq}`eq-02-02-first-order` ⇒ 凸。* 令 $\vz=\theta\vx+(1-\theta)\vy$。在 $\vz$ 处两次应用 {eq}`eq-02-02-first-order`：

$$
\begin{aligned}
f(\vx)&\ge f(\vz)+\nabla f(\vz)^\T(\vx-\vz),\\
f(\vy)&\ge f(\vz)+\nabla f(\vz)^\T(\vy-\vz).
\end{aligned}
$$

第一个不等式乘以 $\theta$，第二个乘以 $1-\theta$，然后相加。由于 $\theta(\vx-\vz)+(1-\theta)(\vy-\vz)=\vzero$，梯度项相互抵消，剩下的就是 $\theta f(\vx)+(1-\theta)f(\vy)\ge f(\vz)$。

*严格情形。* 若对 $\vx\ne\vy$ 严格不等式成立，则对 $\theta\in(0,1)$ 用同样的论证即得严格凸性。反之，设 $f$ 严格凸，$\vx\ne\vy$，$\vz=\tfrac12(\vx+\vy)$。由 $\vx$ 处的非严格条件和严格凸性，

$$
\begin{aligned}
&f(\vx)+\tfrac12\nabla f(\vx)^\T(\vy-\vx)\\
&=f(\vx)+\nabla f(\vx)^\T(\vz-\vx)\\
&\le f(\vz)\\
&<\tfrac12 f(\vx)+\tfrac12 f(\vy),
\end{aligned}
$$

两边乘以 $2$ 并移项即得 $f(\vy)>f(\vx)+\nabla f(\vx)^\T(\vy-\vx)$。
:::

{eq}`eq-02-02-first-order` 的右端是 $f$ 在 $\vx$ 处的一阶泰勒近似，它是 $\vy$ 的仿射函数。定理表明，对凸函数而言，这一局部的线性信息是一个*全局*下估计函数（{numref}`fig-02-02-chord-tangent`(b)）。用上图的语言说，$\R^{n+1}$ 中的超平面 $t=f(\vx)+\nabla f(\vx)^\T(\vy-\vx)$ 在 $(\vx,f(\vx))$ 处支撑 $\operatorname{epi} f$。

最重要的推论与驻点有关。若 $\nabla f(\vx^*)=\vzero$，则由 {eq}`eq-02-02-first-order` 得，对所有 $\vy\in C$ 有 $f(\vy)\ge f(\vx^*)$：**凸函数的驻点（stationary point）是全局极小点。** 对非凸函数而言，驻点可能是局部极小点、局部极大点或鞍点（{doc}`06-02-optimality-conditions`）。这一命题的约束版本是{prf:ref}`thm-02-03-optimality`中的最优性条件。

将 {eq}`eq-02-02-first-order` 交换 $\vx$ 与 $\vy$ 的角色后再用一次，两式相加得

$$
\bigl(\nabla f(\vy)-\nabla f(\vx)\bigr)^\T(\vy-\vx)\ge0
\qquad\text{对所有 }\vx,\vy\in C ,
$$ (eq-02-02-monotone-gradient)

所以凸函数的梯度是**单调**的（monotone）：当 $n=1$ 时，导数 $f'$ 单调不减。

(sec-02-02-second-order)=
## 二阶条件

对二阶可微函数，曲率决定凸性。

:::{prf:theorem} 二阶条件
:label: thm-02-02-second-order

设 $C$ 是凸集，$f$ 在包含 $C$ 的某个开集上二阶连续可微。

1. 若对所有 $\vx\in C$ 有 $\nabla^2f(\vx)\succeq0$，则 $f$ 在 $C$ 上是凸函数。
2. 若对所有 $\vx\in C$ 有 $\nabla^2f(\vx)\succ0$，则 $f$ 在 $C$ 上严格凸。
3. 若 $C$ 是开集且 $f$ 在 $C$ 上是凸函数，则对所有 $\vx\in C$ 有 $\nabla^2f(\vx)\succeq0$。

特别地，在开凸集（例如 $\R^n$）上，$f$ 是凸函数当且仅当它的黑塞矩阵处处半正定。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

*第 1、2 部分。* 设 $\vx,\vy\in C$，$\vx\ne\vy$。由带中值型余项的泰勒定理（{prf:ref}`thm-01-02-taylor`），在线段 $[\vx,\vy]\subseteq C$ 上存在点 $\vz$，使得

$$
f(\vy)=f(\vx)+\nabla f(\vx)^\T(\vy-\vx)
+\tfrac12(\vy-\vx)^\T\nabla^2f(\vz)(\vy-\vx).
$$

若 $\nabla^2f(\vz)\succeq0$，则最后一项非负，于是 {eq}`eq-02-02-first-order` 成立，由{prf:ref}`thm-02-02-first-order`知 $f$ 是凸函数。若 $\nabla^2f(\vz)\succ0$，则最后一项为正，$f$ 严格凸。

*第 3 部分。* 假设对某个 $\vx\in C$ 和某个 $\vv$ 有 $\vv^\T\nabla^2f(\vx)\vv<0$。由于 $C$ 是开集，对所有充分小的 $t>0$ 有 $\vx+t\vv\in C$，由泰勒定理和黑塞矩阵的连续性，当 $t$ 足够小时

$$
\begin{aligned}
f(\vx+t\vv)&=f(\vx)+t\nabla f(\vx)^\T\vv+\tfrac{t^2}{2}\vv^\T\nabla^2f(\vz_t)\vv\\
&<f(\vx)+t\nabla f(\vx)^\T\vv
\end{aligned}
$$

这是因为 $\vz_t\to\vx$。这与 {eq}`eq-02-02-first-order` 矛盾，所以 $\nabla^2f(\vx)\succeq0$。
:::

下面两条注可以避免对该定理的常见误读。

:::{prf:remark} 正定性不是严格凸的必要条件
:label: rem-02-02-strict-hessian

第 2 部分的逆命题不成立。函数 $f(x)=x^4$ 在 $\R$ 上严格凸（它严格满足 {eq}`eq-02-02-first-order`，因为当 $y\ne x$ 时 $y^4-x^4-4x^3(y-x)=(y-x)^2(y^2+2xy+3x^2)>0$），但 $f''(0)=0$。严格凸函数的黑塞矩阵可以在孤立点处奇异。
:::

:::{prf:remark} 黑塞矩阵判别法需要开集
:label: rem-02-02-hessian-open

第 3 部分要求定义域是开集，或者至少内部非空。在没有内点的集合上，黑塞矩阵判别法只是充分条件。例如，$f(\vx)=x_1^2-x_2^2$ 在直线 $C=\{\vx:x_2=0\}$ 上是凸函数（在该直线上它等于 $x_1^2$），尽管 $\nabla^2f=\diag(2,-2)$ 处处不定。*（这一判别法常被不加此条件地表述为“在 $\Omega$ 上是凸函数当且仅当黑塞矩阵在 $\Omega$ 上半正定”；一般情况下只有“当”的方向成立。）*
:::

对一元函数，定理成为：$f$ 在开区间上是凸函数，当且仅当在该区间上 $f''\ge0$。这就是大家熟悉的“图像向上弯曲”的说法。

### 黑塞矩阵的判别

判断 $\nabla^2f(\vx)\succeq0$ 是否成立是一个关于对称矩阵的问题，已在{doc}`01-01-linear-algebra`中解决（{prf:ref}`thm-01-01-eigen-test`、{prf:ref}`thm-01-01-sylvester`、{prf:ref}`thm-01-01-psd-minors`）。对对称矩阵 $\mH$：

- $\mH\succeq0$ 当且仅当所有特征值 $\lambda_i(\mH)\ge0$；$\mH\succ0$ 当且仅当所有 $\lambda_i(\mH)>0$；
- $\mH\succ0$ 当且仅当所有**顺序**主子式（leading principal minor）为正（西尔维斯特判据）；
- $\mH\succeq0$ 当且仅当**所有**主子式（principal minor，即删去相同的行和列后得到的子矩阵的行列式）非负。

后两个判别法不能互换。*顺序*主子式非负并不能推出半正定：$\mH=\diag(0,-1)$ 的顺序主子式为 $0$ 和 $0$，但它有特征值 $-1$。两个子式判别法都假定 $\mH$ 是对称矩阵，而二阶连续可微函数的黑塞矩阵总是对称的。*（省略对称性条件，或只用顺序主子式判别半正定性的子式判别法是错误的。）*

:::{prf:example} 一个凸但不严格凸的二次函数
:label: ex-02-02-hessian-test

$f(\vx)=x_1^2-2x_1x_2+x_2^2$ 在 $\R^2$ 上是凸函数吗？是严格凸的吗？
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

梯度和黑塞矩阵为

$$
\begin{aligned}
\nabla f(\vx)&=\begin{pmatrix}2x_1-2x_2\\-2x_1+2x_2\end{pmatrix},\\
\nabla^2f(\vx)&=\begin{pmatrix}2&-2\\-2&2\end{pmatrix}.
\end{aligned}
$$

黑塞矩阵是常数矩阵。它的特征多项式为 $(2-\lambda)^2-4=\lambda(\lambda-4)$，所以特征值为 $0$ 和 $4$。等价地，它的主子式为 $2$、$2$（对角元素）和 $\det\nabla^2f=4-4=0$，全部非负。因此在开集 $\R^2$ 上 $\nabla^2f\succeq0$，由{prf:ref}`thm-02-02-second-order`知 **$f$ 在 $\R^2$ 上是凸函数**。

黑塞矩阵不是正定的（它有特征值 $0$），所以定理的第 2 部分不适用，而由{prf:ref}`rem-02-02-strict-hessian`，仅凭这一点还不能判定问题的答案。这里的答案是否定的：特征值 $0$ 对应的特征向量是 $(1,1)$，$f$ 沿这个方向完全不弯曲。事实上 $f(\vx)=(x_1-x_2)^2$ 在整条直线 $x_1=x_2$ 上为零。取 $\vx=(0,0)$、$\vy=(1,1)$、$\theta=\tfrac12$，{eq}`eq-02-02-convex` 两边都等于 $0$，严格不等式不成立，所以 **$f$ 是凸函数但不是严格凸的**。由于有特征值 $4$，它也不是凹函数。它的极小点构成整条直线 $x_1=x_2$，这说明凸函数可以有无穷多个全局极小点。
:::

(sec-02-02-examples)=
## 凸函数的例子

下面这些函数是大多数凸模型的基本构件。每个结论都可以用{prf:ref}`thm-02-02-second-order`或定义来验证。

**仿射函数** $\va^\T\vx+b$ 既是凸函数又是凹函数。

**二次函数** $f(\vx)=\tfrac12\vx^\T\mQ\vx+\vc^\T\vx+d$（$\mQ$ 对称）满足 $\nabla^2f=\mQ$。它是凸函数当且仅当 $\mQ\succeq0$，是严格凸的当且仅当 $\mQ\succ0$（对二次函数而言，若 $\vv\ne\vzero$ 使 $\vv^\T\mQ\vv=0$，则 $f$ 沿 $\vv$ 是仿射的，如{prf:ref}`ex-02-02-hessian-test`）。最小二乘目标函数 $\norm{\mA\vx-\vb}^2$ 的黑塞矩阵为 $2\mA^\T\mA\succeq0$，因而是凸函数。

**范数。** 每个范数都是凸函数：由三角不等式和齐次性，$\norm{\theta\vx+(1-\theta)\vy}\le\theta\norm{\vx}+(1-\theta)\norm{\vy}$。没有范数是严格凸的，因为沿射线 $\norm{t\vx}=t\norm{\vx}$ 是线性的。欧氏范数的平方 $\norm{\vx}^2=\vx^\T\vx$ 是严格凸的（$\mQ=2\mI$）。范数在 $\vzero$ 处不可微，这表明凸性并不要求可微性。

**一元函数。**

- 对每个 $a\in\R$，$e^{ax}$ 在 $\R$ 上是凸函数（二阶导数 $a^2e^{ax}\ge0$）。
- $(0,\infty)$ 上的 $x^p$：当 $p\ge1$ 或 $p\le0$ 时是凸函数，当 $0\le p\le1$ 时是凹函数（二阶导数 $p(p-1)x^{p-2}$）。
- $\log x$ 在 $(0,\infty)$ 上是凹函数，$-\log x$ 是凸函数（$\log x$ 的二阶导数为 $-1/x^2$）。
- $x\log x$ 在 $(0,\infty)$ 上是凸函数（二阶导数 $1/x$）。
- 当 $p\ge1$ 时，$\abs{x}^p$ 在 $\R$ 上是凸函数。

**仿射函数的最大值。** $f(\vx)=\max_{i=1,\dots,k}(\va_i^\T\vx+b_i)$ 是凸函数，因为它的上图是 $k$ 个半空间 $\{(\vx,t):\va_i^\T\vx+b_i\le t\}$ 的交。这类分段线性凸函数可以描述边际费率递增的成本，并且可以用线性规划求其极小（{ref}`sec-07-03-convex-case`）。

**log-sum-exp 函数。** 函数 $\operatorname{lse}(\vx)=\log\bigl(e^{x_1}+\dots+e^{x_n}\bigr)$ 是最大值函数的光滑近似：

$$
\max_i x_i\le\operatorname{lse}(\vx)\le\max_i x_i+\log n ,
$$ (eq-02-02-lse-bounds)

因为 $e^{\max_i x_i}\le\sum_ie^{x_i}\le n\,e^{\max_i x_i}$。对 $\vx=(1,2,3)$，$\operatorname{lse}(\vx)=3.408$，介于 $3$ 与 $3+\log3=4.099$ 之间。它的黑塞矩阵为 $\diag(\vp)-\vp\vp^\T$，其中 $p_i=e^{x_i}/\sum_ke^{x_k}$，于是 $\vp\ge\vzero$ 且 $\vone^\T\vp=1$。对任意 $\vv$，由柯西–施瓦茨不等式，

$$
\vv^\T\nabla^2\operatorname{lse}(\vx)\vv
=\sum_ip_iv_i^2-\Bigl(\sum_ip_iv_i\Bigr)^2\ge0 ,
$$

这是因为 $\bigl(\sum_ip_iv_i\bigr)^2=\bigl(\sum_i\sqrt{p_i}\,\sqrt{p_i}v_i\bigr)^2\le \sum_ip_i\sum_ip_iv_i^2$。因此 log-sum-exp 函数是凸函数。它出现在 Logit 路径选择模型以及出行需求的离散选择模型中。

(sec-02-02-operations)=
## 保凸运算

与集合的情形一样（{ref}`sec-02-01-operations`），证明凸性的通常做法是用已知的凸“部件”组装出所需的函数。

:::{prf:theorem} 保凸运算
:label: thm-02-02-operations

1. **非负加权和。** 若 $f_1,\dots,f_k$ 是凸函数且 $w_1,\dots,w_k\ge0$，则 $\sum_iw_if_i$ 是凸函数。
2. **与仿射映射复合。** 若 $f$ 是凸函数，则 $g(\vx)=f(\mA\vx+\vb)$ 是凸函数。
3. **逐点最大值与上确界。** 若对指标集 $\mathcal I$（有限或无限）中的每个 $i$，$f_i$ 都是凸函数，则 $g(\vx)=\sup_{i\in\mathcal I}f_i(\vx)$ 是凸函数。
4. **与标量函数复合。** 设 $g:\R^n\to\R$，$h:\R\to\R$。若 $h$ 是单调不减的凸函数且 $g$ 是凸函数，或 $h$ 是单调不增的凸函数且 $g$ 是凹函数，则 $h(g(\vx))$ 是凸函数。
5. **部分极小化。** 若 $f(\vx,\vy)$ 关于 $(\vx,\vy)$ 联合凸，$C$ 是非空凸集，则 $g(\vx)=\inf_{\vy\in C}f(\vx,\vy)$ 是凸函数（在其取有限值处）。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

第 1、2 部分可由 {eq}`eq-02-02-convex` 直接得到；在第 2 部分中，映射 $\vx\mapsto\mA\vx+\vb$ 保持凸组合。对第 3 部分，$\operatorname{epi} g=\bigcap_i\operatorname{epi} f_i$ 是凸集的交（{prf:ref}`thm-02-01-intersection`）。对第 4 部分的第一种情形，由 $g$ 的凸性和 $h$ 的单调性得 $h(g(\theta\vx+(1-\theta)\vy))\le h(\theta g(\vx)+(1-\theta)g(\vy))$，再由 $h$ 的凸性，上式右端不超过 $\theta h(g(\vx))+(1-\theta)h(g(\vy))$；第二种情形的证明相同，只需把关于 $g$ 的不等式反向。对第 5 部分，任给 $\varepsilon>0$，选取 $\vy_1,\vy_2\in C$，使 $f(\vx_1,\vy_1)\le g(\vx_1)+\varepsilon$，$f(\vx_2,\vy_2)\le g(\vx_2)+\varepsilon$。由于 $\theta\vy_1+(1-\theta)\vy_2\in C$，

$$
\begin{aligned}
&g\bigl(\theta\vx_1+(1-\theta)\vx_2\bigr)\\
&\le f\bigl(\theta\vx_1+(1-\theta)\vx_2,\ \theta\vy_1+(1-\theta)\vy_2\bigr)\\
&\le\theta g(\vx_1)+(1-\theta)g(\vx_2)+\varepsilon ,
\end{aligned}
$$

而 $\varepsilon$ 是任意的。
:::

一些推论：

- $\norm{\mA\vx-\vb}$、$\norm{\mA\vx-\vb}^2$ 和 $\max_i(\va_i^\T\vx-b_i)$ 都是凸函数（第 2、3 部分）；
- 当 $g$ 是凸函数时 $e^{g(\vx)}$ 是凸函数；当 $g$ 是正的凹函数时 $1/g(\vx)$ 和 $-\log g(\vx)$ 是凸函数（第 4 部分，分别取 $h(u)=e^u$、$h(u)=1/u$ 和 $h(u)=-\log u$，$u>0$）；
- 点到凸集的距离 $\inf_{\vy\in C}\norm{\vx-\vy}$ 是凸函数（第 5 部分）；
- 仿射函数的逐点**下确界**是凹函数，例如对偶函数 {eq}`eq-02-03-dual-function`（对其相反数应用第 3 部分）。

凸函数的乘积和复合一般不是凸函数：$x$ 和 $x^2$ 都是凸函数，但 $x\cdot x^2=x^3$ 不是；$h(u)=-u$ 是凸函数，而 $h(x^2)=-x^2$ 是凹函数。

**一个基础设施的例子。** 在交通分配中，路段上的行程时间是路段流量 $x$（辆/小时）的函数。一种广泛使用的形式是美国公路局（US Bureau of Public Roads，BPR）的路段阻抗函数（link performance function），

$$
t(x)=t_0\Bigl(1+0.15\bigl(x/c\bigr)^4\Bigr),
$$ (eq-02-02-bpr)

其中 $t_0$ 为自由流行程时间，$c$ 为实际通行能力。其二阶导数 $1.8\,t_0x^2/c^4$ 非负，所以当 $x\ge0$ 时 $t$ 是凸函数：每增加一辆车对其他车辆造成的延误都比前一辆车更大。路段上所有车辆的总行程时间 $x\,t(x)=t_0x+0.15\,t_0x^5/c^4$ 是一个线性函数与一个二阶导数为 $3\,t_0x^3/c^4\ge0$ 的函数之和，因此当 $x\ge0$ 时是凸函数；而路网的总行程时间是各路段总行程时间之和，因而是路段流量的凸函数（第 1 部分）。这就是系统最优交通分配是凸优化问题的原因（Sheffi, 1985）。

(sec-02-02-sublevel)=
## 下水平集

:::{prf:definition} 下水平集
:label: def-02-02-sublevel-set

$f$ 在水平 $\gamma\in\R$ 上的**下水平集**（sublevel set）是 $S_\gamma=\{\vx\in\dom f: f(\vx)\le\gamma\}$。
:::

:::{prf:theorem} 凸函数的下水平集是凸集
:label: thm-02-02-sublevel

若 $f$ 是凸函数，则每个下水平集 $S_\gamma$ 都是凸集。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

若 $f(\vx)\le\gamma$ 且 $f(\vy)\le\gamma$，则 $f(\theta\vx+(1-\theta)\vy)\le\theta\gamma+(1-\theta)\gamma=\gamma$。
:::

{numref}`fig-02-02-epigraph-sublevel`(b) 画出了一个凸二次函数的下水平集：它们是实心椭圆，与{prf:ref}`def-02-01-ellipsoid`一致。这个定理说明了为什么凸优化中的约束写成 $g_i(\vx)\le0$ 且 $g_i$ 为凸函数的形式：集合 $\{\vx:g_i(\vx)\le0\}$ 是 $g_i$ 的下水平集 $S_0$，因而是凸集；由{prf:ref}`thm-02-01-intersection`，可行集 $\bigcap_i\{\vx:g_i(\vx)\le0\}$ 也是凸集。类似地，当 $g$ 是凹函数时 $\{\vx:g(\vx)\ge0\}$ 是凸集。而对非线性凸函数 $g$，*等式* $g(\vx)=0$ 给出的不是凸集：$x_1^2+x_2^2-1=0$ 是一个圆周。

```{figure} figures/02-02-epigraph-sublevel.svg
:name: fig-02-02-epigraph-sublevel
:width: 640px
:class: dark-light
:alt: 左图：一条凸曲线，其上方区域作为上图涂上阴影，一条位于水平 gamma 的水平线，以及曲线位于该线之下的区间。右图：一个凸二次函数的嵌套椭圆等值线，其中一个下水平集涂上阴影

(a) 一元凸函数的上图（阴影部分）是凸集。下水平集 $S_\gamma$ 是图像位于水平 $\gamma$ 之下或与之相交的 $x$ 构成的区间。(b) 凸二次函数 $f(\vx)=(x_1-1)^2+2(x_2-1)^2-\tfrac12x_1x_2$ 的等值线；下水平集 $S_1$ 涂有阴影。所有下水平集都是相互嵌套的凸椭圆。
```

{prf:ref}`thm-02-02-sublevel`的逆命题不成立。函数 $\sqrt{\abs{x}}$ 不是凸函数（它在 $[0,\infty)$ 上是凹的），但它的下水平集 $\{x:\abs{x}\le\gamma^2\}$ 都是区间。所有下水平集都是凸集的函数构成一个更大的函数类，即拟凸函数（{ref}`sec-02-02-quasiconvex`）。

(sec-02-02-jensen)=
## Jensen 不等式

定义 {eq}`eq-02-02-convex` 涉及两个点。正如凸集包含其点的所有凸组合（{prf:ref}`thm-02-01-convex-combinations`）一样，它可以推广到任意有限个点。

:::{prf:theorem} Jensen 不等式
:label: thm-02-02-jensen

设 $f$ 是凸函数，$\vx_1,\dots,\vx_k\in\dom f$，$\theta_i\ge0$ 且 $\sum_{i=1}^k\theta_i=1$。则

$$
f\Bigl(\sum_{i=1}^k\theta_i\vx_i\Bigr)\le\sum_{i=1}^k\theta_if(\vx_i).
$$ (eq-02-02-jensen)

若 $X$ 是取值于 $\dom f$ 且期望有限的随机向量，则 $f(\mathrm E\,X)\le\mathrm E\,f(X)$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

对 $\operatorname{epi} f$ 应用{prf:ref}`thm-02-01-convex-combinations`，由{prf:ref}`thm-02-02-epigraph`知它是凸集：点 $(\vx_i,f(\vx_i))$ 属于上图，从而它们的凸组合 $\bigl(\sum_i\theta_i\vx_i,\ \sum_i\theta_if(\vx_i)\bigr)$ 也属于上图，这正是 {eq}`eq-02-02-jensen`。概率形式是同一命题在概率为 $\theta_i$ 的离散分布下的表述，对一般分布可通过取极限得到。
:::

许多经典不等式都源于 Jensen 不等式。例如，对凹函数 $\log x$ 应用它可得算术–几何平均不等式：对 $a,b>0$ 有 $\sqrt{ab}\le\tfrac12(a+b)$，因为 $\log\tfrac12(a+b)\ge\tfrac12\log a+\tfrac12\log b$。在规划中它提醒我们：在*平均*需求处计算的凸成本会低估*平均*成本。

:::{prf:example} 流量波动下的行程时间
:label: ex-02-02-jensen-travel-time

某路段的自由流行程时间为 $t_0=10$ 分钟，通行能力为 $c=2{,}000$ 辆/小时，行程时间由 {eq}`eq-02-02-bpr` 给出。早高峰期间，一半日子的流量为 $1{,}600$ 辆/小时，另一半日子为 $2{,}400$ 辆/小时，因此平均流量为 $2{,}000$ 辆/小时。比较平均流量下的行程时间与平均行程时间。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

在平均流量下，$t(2000)=10\,(1+0.15\cdot1^4)=11.5$ 分钟。在两类日子里，

$$
\begin{aligned}
t(1600)&=10\,\bigl(1+0.15\cdot0.8^4\bigr)=10\,(1+0.15\cdot0.4096)=10.61\ \text{分钟},\\
t(2400)&=10\,\bigl(1+0.15\cdot1.2^4\bigr)=10\,(1+0.15\cdot2.0736)=13.11\ \text{分钟},
\end{aligned}
$$

所以平均行程时间为 $\tfrac12(10.61+13.11)=11.86$ 分钟。正如 Jensen 不等式对凸函数 $t$ 所预言的，平均行程时间超过了平均流量下的行程时间，这里超出 $0.36$ 分钟，即 $3.2\,\%$。用平均流量做规划会低估平均延误，并且误差随流量的波动程度和 $t$ 的曲率增大而增大。
:::

(sec-02-02-quasiconvex)=
## 拟凸函数

有些有用的函数不是凸函数，但具有凸函数下水平集为凸集这一性质。

:::{prf:definition} 拟凸函数
:label: def-02-02-quasiconvex

设 $C$ 是凸集，$f:C\to\R$。若

$$
f\bigl(\theta\vx+(1-\theta)\vy\bigr)\le\max\{f(\vx),f(\vy)\}
\qquad\text{对所有 }\vx,\vy\in C,\ \theta\in[0,1],
$$ (eq-02-02-quasiconvex)

则称 $f$ 是**拟凸函数**（quasiconvex function）；若当 $\vx\ne\vy$ 且 $\theta\in(0,1)$ 时不等式严格成立，则称 $f$ **严格拟凸**（strictly quasiconvex）。若 $-f$ 是（严格）拟凸函数，则称 $f$ 是（严格）**拟凹函数**（quasiconcave function）。
:::

函数是拟凸函数，当且仅当它的所有下水平集都是凸集：若 $f(\vx),f(\vy)\le\gamma$，则由 {eq}`eq-02-02-quasiconvex` 得 $f(\theta\vx+(1-\theta)\vy)\le\gamma$；反之，取 $\gamma=\max\{f(\vx),f(\vy)\}$ 即得 {eq}`eq-02-02-quasiconvex`。每个凸函数都是拟凸函数，因为 {eq}`eq-02-02-convex` 中的弦值至多为 $\max\{f(\vx),f(\vy)\}$；每个严格凸函数都是严格拟凸的。每个一元单调函数既是拟凸函数又是拟凹函数；$\sqrt{\abs{x}}$ 和 $\log x$（在 $x>0$ 上）是拟凸函数，但不是凸函数。

在一元情形，区间上的严格拟凸函数是**单峰**的（unimodal）：它在极小点之前递减，之后递增，且不需要任何可微性。黄金分割法等线搜索方法所依赖的只有这一性质（{prf:ref}`def-06-04-unimodal`、{prf:ref}`alg-06-04-golden-section`）。有两点需要注意。拟凸函数可能有平坦部分，因此除非函数严格拟凸，否则局部极小点不一定是全局极小点。此外，拟凸函数之和不一定是拟凸函数：$x^3$ 和 $-3x$ 都是单调的，因而是拟凸函数，但 $x^3-3x$ 的下水平集 $\{x:x^3-3x\le0\}=(-\infty,-\sqrt3]\cup[0,\sqrt3]$ 不是区间。

## 本章回顾

### 小结

- 凸集上的函数是凸函数，当且仅当每条弦都位于其图像之上（{prf:ref}`def-02-02-convex-function`、{numref}`fig-02-02-chord-tangent`），也当且仅当其上图是凸集（{prf:ref}`thm-02-02-epigraph`）。仿射函数既凸又凹；抛物线是凸的当且仅当其首项系数非负（{eq}`eq-02-02-parabola-gap`）。
- 可微函数是凸函数当且仅当它位于其切平面之上（{prf:ref}`thm-02-02-first-order`）；因此凸函数的驻点是全局极小点。
- 开凸集上的二阶连续可微函数是凸函数，当且仅当其黑塞矩阵半正定；黑塞矩阵正定蕴含严格凸，但反之不成立（{prf:ref}`thm-02-02-second-order`、{prf:ref}`rem-02-02-strict-hessian`、{prf:ref}`rem-02-02-hessian-open`）。函数 $x_1^2-2x_1x_2+x_2^2$ 是凸函数但不是严格凸的（{prf:ref}`ex-02-02-hessian-test`）。
- 范数、凸二次函数、指数函数、$-\log x$、仿射函数的最大值和 log-sum-exp 函数都是凸函数；非负加权和、与仿射映射复合、逐点最大值、单调复合和部分极小化都保持凸性（{prf:ref}`thm-02-02-operations`）。
- 凸函数的下水平集是凸集（{prf:ref}`thm-02-02-sublevel`），所以 $g_i$ 为凸函数的约束 $g_i(\vx)\le0$ 定义的可行集是凸集。
- Jensen 不等式（{prf:ref}`thm-02-02-jensen`）把凸性推广到平均值；对于凸的行程时间函数，平均行程时间超过平均流量下的行程时间（{prf:ref}`ex-02-02-jensen-travel-time`）。
- 拟凸函数（{prf:ref}`def-02-02-quasiconvex`）恰好是下水平集均为凸集的函数；一元严格拟凸函数是单峰的。

### 常见错误

- **忘记定义域。** 凸性是函数*在某个凸集上*的性质；$1/x$ 只在 $x>0$ 上是凸函数。
- **要求在所有 $\theta\in[0,1]$ 处严格成立。** 严格凸性比较的是不同的两点，且 $\theta\in(0,1)$；在端点处 {eq}`eq-02-02-convex` 两边相等。
- **把“$\nabla^2f\succ0$”当作严格凸的必要条件。** 它只是充分条件；$x^4$ 严格凸，但 $f''(0)=0$。
- **用顺序主子式判别半正定性。** 顺序主子式判定的是*正定性*；半正定性需要考察所有主子式或特征值（$\diag(0,-1)$ 就是反例）。
- **在没有内点的集合上把黑塞矩阵判别法当作“当且仅当”使用**（{prf:ref}`rem-02-02-hessian-open`）。
- **在“凸”模型中写非线性等式约束。** 只有仿射等式才能保持可行集的凸性。
- **以为凸函数的乘积或复合仍是凸函数。** 应使用{prf:ref}`thm-02-02-operations`中的复合规则。

## 参考文献

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 3, Sections
  3.1–3.3 and 3.5: convex functions, subgradients, differentiable convex
  functions, generalizations of convex functions.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Chapter 3,
  Sections 3.1–3.4: basic properties and examples, operations that preserve
  convexity, quasiconvex functions.)
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming*
  (4th ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3>
  (Section 7.4, convex and concave functions.)
- Sheffi, Y. (1985). *Urban transportation networks: Equilibrium analysis
  with mathematical programming methods*. Prentice-Hall. (Chapter 3: the
  user-equilibrium and system-optimum programs and their convexity, with
  convex link performance functions.)
