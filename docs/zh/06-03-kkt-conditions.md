# Karush–Kuhn–Tucker 条件

## 概述

{doc}`06-02-optimality-conditions`通过方向刻画了局部极小点：在极小点处，没有一个可行方向是下降方向。要使用这一准则，必须知道可行方向，而对于由非线性不等式和方程描述的集合，可行方向并不显然。本章用约束函数代替方向。把在 $\vx^*$ 处有效的约束线性化，得到一个容易写出的方向锥 $L(\vx^*)$；**约束规范**（constraint qualification）保证它描述的是 $\Omega$ 真实的局部几何。在这一假设下，Farkas 引理（Farkas' lemma，等价地，线性规划（linear programming, LP）对偶）把“不存在线性化的可行下降方向”转化为乘子 $\lambda_i\ge0$ 和 $\mu_j$ 的存在性，它们满足 $\nabla f(\vx^*)+\sum_i\lambda_i\nabla g_i(\vx^*)+\sum_j\mu_j\nabla h_j(\vx^*)=\vzero$：这就是 **Karush–Kuhn–Tucker（KKT）条件**（Karush–Kuhn–Tucker conditions，{prf:ref}`thm-06-03-kkt`）。接着，本章从几何上并通过拉格朗日函数解释乘子，确定最大化问题和“$\ge$”约束的符号约定，用有效集分情形讨论求解一个 KKT 方程组（{prf:ref}`ex-06-03-kkt-example`），给出临界锥上的二阶条件，并证明对于凸规划，KKT 点是全局极小点（{prf:ref}`thm-06-03-convex-kkt`）。除前两章外，本章还用到{doc}`03-05-lp-duality`中的 Farkas 引理和 LP 对偶、{doc}`02-02-convex-functions`中的凸函数，并从另一个角度再次遇到{doc}`02-03-convex-optimization-duality`中的拉格朗日函数。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 确定可行点处的有效集 $I(\vx)$ 和线性化可行方向锥；
- 解释为什么线性化锥可能与真实的局部几何不同，定义切锥，并检验常用的约束规范：线性约束、线性无关约束规范（LICQ）和 Slater 条件；
- 由 Farkas 引理推导 KKT 条件，并用拉格朗日函数 $L(\vx,\vlambda,\vmu)$ 写出含“$\le$”“$\ge$”和“$=$”约束的最小化问题与最大化问题的 KKT 条件；
- 通过对有效集分情形讨论求解小规模 KKT 方程组，并验证其解；
- 应用临界锥上的二阶必要条件和二阶充分条件；
- 利用凸性断定 KKT 点是全局极小点。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $I(\vx)$ | 不等式约束在 $\vx$ 处的有效集 $\{i:g_i(\vx)=0\}$ |
| $L(\vx^*)$ | $\vx^*$ 处的线性化可行方向锥（{eq}`eq-06-03-linearised-cone`） |
| $T(\vx^*)$ | $\Omega$ 在 $\vx^*$ 处的切锥（{prf:ref}`def-06-03-tangent-cone`） |
| $\lambda_i$, $\vlambda$ | $g_i(\vx)\le0$ 的乘子；$\vlambda\in\R^m$ |
| $\mu_j$, $\vmu$ | $h_j(\vx)=0$ 的乘子；$\vmu\in\R^\ell$ |
| $L(\vx,\vlambda,\vmu)$ | 拉格朗日函数 $f(\vx)+\vlambda^\T\vg(\vx)+\vmu^\T\vh(\vx)$ |
| $\nabla_{\vx}L$, $\nabla^2_{\vx\vx}L$ | 拉格朗日函数关于 $\vx$ 的梯度和黑塞矩阵 |
| $C(\vx^*,\vlambda^*)$ | 临界锥（{prf:ref}`def-06-03-critical-cone`） |

字母 $L$ 既用于只有一个自变量的锥 $L(\vx^*)$，也用于有三个自变量的拉格朗日函数 $L(\vx,\vlambda,\vmu)$；由自变量的个数即可区分。本章始终讨论问题 {eq}`eq-06-01-nlp`，

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)\\
\st\ & g_i(\vx)\le0, && i=1,\dots,m,\\
& h_j(\vx)=0, && j=1,\dots,\ell ,
\end{aligned}
$$

其中 $f$、$g_i$、$h_j$ 连续可微。

(sec-06-03-active)=
## 有效约束与线性化方向

在可行点 $\vx$ 处，若 $g_i(\vx)=0$，则称不等式约束为**有效约束**（active constraint，也称紧约束）；若 $g_i(\vx)<0$，则称其为**非有效约束**（inactive constraint）。**有效集**（active set）为

$$
I(\vx)=\{i\in\{1,\dots,m\}:\ g_i(\vx)=0\}.
$$ (eq-06-03-active-set)

等式约束在每个可行点处都以等式成立，算作有效约束。非有效约束有松弛量，由连续性，它在 $\vx$ 的某个邻域内保持满足；它不限制短距离的移动。只有有效约束决定 $\vx$ 附近可行集的形状。

:::{prf:example} 三角形的有效约束
:label: ex-06-03-triangle-active

设 $\Omega=\{\vx\in\R^2:x_1+x_2\le2,\ x_1\ge0,\ x_2\ge0\}$，写成 $g_1(\vx)=x_1+x_2-2\le0$，$g_2(\vx)=-x_1\le0$，$g_3(\vx)=-x_2\le0$，梯度为 $\nabla g_1=(1,1)^\T$，$\nabla g_2=(-1,0)^\T$，$\nabla g_3=(0,-1)^\T$。求顶点 $(2,0)$ 和 $(0,0)$ 处的有效集和可行方向。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

由于约束是线性的，$g_i(\vx^*+\alpha\vd)=g_i(\vx^*)+\alpha\nabla g_i^\T\vd$ 精确成立；若对所有充分小的 $\alpha>0$ 和所有 $i$ 它都 $\le0$，则方向 $\vd$ 是可行的。

```{list-table} 三角形的约束 $g_1=x_1+x_2-2$、$g_2=-x_1$、$g_3=-x_2$ 在顶点 $(2,0)$ 和 $(0,0)$ 处的情况：$g_i$ 的值（为 $0$ 时有效），以及每个约束对可行方向施加的条件（“对 $\vd$”），即对每个有效的 $g_i$ 有 $\nabla g_i^\T\vd\le0$。
:name: tab-06-03-triangle-active
:header-rows: 1

* - $g_i$
  - $g_i(2,0)$
  - 对 $\vd$
  - $g_i(0,0)$
  - 对 $\vd$
* - $g_1$
  - $0$
  - $d_1+d_2\le0$
  - $-2$
  - 无
* - $g_2$
  - $-2$
  - 无
  - $0$
  - $-d_1\le0$
* - $g_3$
  - $0$
  - $-d_2\le0$
  - $0$
  - $-d_2\le0$
```

因此 $I(2,0)=\{1,3\}$，那里的可行方向为 $\{\vd:d_1+d_2\le0,\ d_2\ge0\}$；$I(0,0)=\{2,3\}$，可行方向为 $\{\vd:\vd\ge\vzero\}$。
:::

如{prf:ref}`ex-06-03-triangle-active`的解答所示，在两个顶点处，一个方向可行当且仅当它与每个有效约束的梯度的夹角至少为 $90^\circ$。梯度 $\nabla g_i$ 指向 $\Omega$ 之外，即 $g_i$ 增大、趋于违反约束的方向。

对于非线性约束，这一夹角条件仍然是必要的。

:::{prf:proposition} 可行方向满足线性化约束
:label: prop-06-03-linearised

若 $\vd$ 是 $\Omega$ 在 $\vx^*\in\Omega$ 处的可行方向，则

$$
\begin{aligned}
\nabla g_i(\vx^*)^\T\vd&\le0, && i\in I(\vx^*),\\
\nabla h_j(\vx^*)^\T\vd&=0, && j=1,\dots,\ell.
\end{aligned}
$$
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

对于 $i\in I(\vx^*)$ 和 $0<\alpha\le\bar\alpha$，有 $g_i(\vx^*+\alpha\vd)\le0=g_i(\vx^*)$，因此差商 $\bigl(g_i(\vx^*+\alpha\vd)-g_i(\vx^*)\bigr)/\alpha$ 的值 $\le0$；它在 $\alpha\downarrow0$ 时的极限是方向导数 $\nabla g_i(\vx^*)^\T\vd$，所以后者 $\le0$。对于 $h_j$，差商恒为 $0$。
:::

该命题所允许的方向构成**线性化可行方向锥**（cone of linearized feasible directions）

$$
\begin{aligned}
L(\vx^*)=\bigl\{\vd\in\R^n:\ &\nabla g_i(\vx^*)^\T\vd\le0\ \ (i\in I(\vx^*)),\\
&\nabla h_j(\vx^*)^\T\vd=0\ \ (j=1,\dots,\ell)\bigr\}.
\end{aligned}
$$ (eq-06-03-linearised-cone)

它是一个多面体锥，由关于 $\vd$ 的线性不等式定义，仅凭梯度就能写出。该命题说明每个可行方向都在 $L(\vx^*)$ 中。问题在于，反过来 $L(\vx^*)$ 是否如实地描述了可以移动的方向。

(sec-06-03-cq)=
## 约束规范

下面的例子表明，即使对于非常简单的集合，答案也可能是否定的。

:::{prf:example} 三角形的三次描述
:label: ex-06-03-cubic

设

$$
\begin{aligned}
g_1(\vx)&=(x_1+x_2-1)^3\le0,\\
g_2(\vx)&=-x_1\le0,\\
g_3(\vx)&=-x_2\le0,
\end{aligned}
$$

于是 $\Omega=\{\vx\in\R^2:g_1(\vx)\le0,\ g_2(\vx)\le0,\ g_3(\vx)\le0\}$ 是以 $(0,0)$、$(1,0)$、$(0,1)$ 为顶点的三角形；设 $\vx^*=(0.5,0.5)$。(a) 求 $I(\vx^*)$ 和 $L(\vx^*)$。(b) $\vd=(2,1)$ 是可行方向吗？(c) 在 $\Omega$ 上最小化 $f(\vx)=(x_1-1)^2+(x_2-1)^2$，并检验在极小点处是否存在 {eq}`eq-06-03-kkt-stationarity` 中那样的乘子。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** $g_1(\vx^*)=0$，$g_2(\vx^*)=g_3(\vx^*)=-0.5<0$，因此 $I(\vx^*)=\{1\}$。梯度 $\nabla g_1(\vx)=3(x_1+x_2-1)^2(1,1)^\T$ 在 $\vx^*$ 处为零，所以 {eq}`eq-06-03-linearised-cone` 中唯一的条件成为 $0\le0$，$L(\vx^*)=\R^2$：每个方向都能通过线性化检验。

**(b)** 沿 $\vd=(2,1)$，$x_1+x_2-1=0.5+2\alpha+0.5+\alpha-1=3\alpha$，因此对每个 $\alpha>0$ 有 $g_1(\vx^*+\alpha\vd)=(3\alpha)^3=27\alpha^3>0$。该方向立即离开 $\Omega$：尽管它属于 $L(\vx^*)$，但它不可行。$g_1$ 的线性化是零函数，它不包含 $\Omega$ 位于直线 $x_1+x_2=1$ 哪一侧的任何信息。真正的可行方向是 $\{\vd:d_1+d_2\le0\}$。

**(c)** 极小点是三角形中离 $(1,1)$ 最近的点，即从 $(1,1)$ 向边 $x_1+x_2=1$ 所作垂线的垂足 $\vx^*=(0.5,0.5)$，$f^*=0.5$（求解器证实了这一点）。在该点 $\nabla f(\vx^*)=(-1,-1)^\T$。由于 $g_2$ 和 $g_3$ 是非有效约束，它们的乘子必须为零，平稳性条件 $\nabla f(\vx^*)+\lambda_1\nabla g_1(\vx^*)=\vzero$ 化为 $(-1,-1)^\T+\lambda_1(0,0)^\T=\vzero$，它无解。该极小点不是 KKT 点。
:::

如果用线性约束 $x_1+x_2-1\le0$ 描述{prf:ref}`ex-06-03-cubic`的三角形，则 $\nabla g_1=(1,1)^\T$，$L(\vx^*)=\{d_1+d_2\le0\}$ 恰好就是可行方向集，且 $\lambda_1=1$ 满足 $(-1,-1)^\T+1\cdot(1,1)^\T=\vzero$。集合相同，只是代数描述不同。线性化锥依赖于描述方式，而真实的局部几何只依赖于集合本身。{numref}`fig-06-03-cubic-cq`显示了这两种效应。

```{figure} figures/06-03-cubic-cq.svg
:name: fig-06-03-cubic-cq
:width: 680px
:class: dark-light
:alt: 左图为一个三角形，斜边上有点 (0.5,0.5)，箭头 d = (2,1) 离开三角形；右图为单位圆盘，点 (1,0) 处有向外的梯度 (2,0) 和竖直的切向箭头 (0,1)，该箭头离开圆盘，而圆周向内弯曲

线性化锥可能产生误导的两种情形。左（{prf:ref}`ex-06-03-cubic`）：取 $g_1=(x_1+x_2-1)^3$ 时，梯度在 $\vx^*=(0.5,0.5)$ 处为零，因此每个方向（包括 $\vd=(2,1)$）都能通过线性化检验，但 $\vd$ 指向 $\Omega$ 之外。右：对于圆盘 $x_1^2+x_2^2\le1$ 在 $(1,0)$ 处，梯度 $(2,0)$ 非零，切方向 $(0,1)$ 能通过检验，但沿它的每一步都离开圆盘；只有随边界弯曲的曲线才保持可行。
```

第二个例子表明，问题不只出在退化的梯度上。对于单位圆盘 $\Omega=\{\vx:g(\vx)=x_1^2+x_2^2-1\le0\}$ 在 $\vx^*=(1,0)$ 处，$\nabla g(\vx^*)=(2,0)^\T\ne\vzero$，$L(\vx^*)=\{\vd:d_1\le0\}$。方向 $\vd=(0,1)$ 属于 $L(\vx^*)$，但对所有 $\alpha>0$ 有 $g(1,\alpha)=\alpha^2>0$，因此它不是可行方向：在严格弯曲的边界的点处，没有一个切方向是可行的。尽管如此，仍然可以沿着圆周 $\vx(t)=(\cos t,\sin t)$“沿”$(0,1)$ 移动，该曲线在 $t=0$ 处的速度为 $(0,1)$。正确的研究对象是可行*曲线*，而不是可行*直线*。

:::{prf:definition} 切锥
:label: def-06-03-tangent-cone

若存在点 $\vx^{(k)}\in\Omega$ 和数 $t_k>0$，使得 $\vx^{(k)}\to\vx^*$，$t_k\to0$，且

$$
\frac{\vx^{(k)}-\vx^*}{t_k}\to\vd\qquad(k\to\infty).
$$

则称向量 $\vd$ 为 $\Omega$ 在 $\vx^*\in\Omega$ 处的**切方向**（tangent direction）。切方向的集合 $T(\vx^*)$ 称为 $\Omega$ 在 $\vx^*$ 处的**切锥**（tangent cone）。
:::

每个可行方向都是切方向（取 $\vx^{(k)}=\vx^*+t_k\vd$），过 $\vx^*$ 的光滑可行曲线的每个速度向量也是切方向。对于圆盘在 $(1,0)$ 处，$T(\vx^*)=\{d_1\le0\}=L(\vx^*)$；对于三角形的三次描述，$T(\vx^*)=\{d_1+d_2\le0\}\ne L(\vx^*)=\R^2$。切锥只依赖于集合 $\Omega$。

:::{prf:lemma} 切锥与一阶条件
:label: lem-06-03-tangent

设 $\vx^*\in\Omega$。

1. $T(\vx^*)\subseteq L(\vx^*)$。
2. 若 $\vx^*$ 是 $f$ 在 $\Omega$ 上的局部极小点，则对所有 $\vd\in T(\vx^*)$ 有 $\nabla f(\vx^*)^\T\vd\ge0$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

设 $\vd\in T(\vx^*)$，序列 $\vx^{(k)}$、$t_k$ 如定义所述，记 $\vx^{(k)}=\vx^*+t_k\vd^{(k)}$，其中 $\vd^{(k)}\to\vd$。对于可微函数 $\phi$，

$$
\phi(\vx^{(k)})=\phi(\vx^*)+t_k\nabla\phi(\vx^*)^\T\vd^{(k)}+o(t_k).
$$

1. 对于 $i\in I(\vx^*)$，取 $\phi=g_i$：由 $g_i(\vx^{(k)})\le0=g_i(\vx^*)$ 得 $t_k\nabla g_i(\vx^*)^\T\vd^{(k)}+o(t_k)\le0$；除以 $t_k$ 并令 $k\to\infty$，得 $\nabla g_i(\vx^*)^\T\vd\le0$。对于 $h_j$，用等式作同样的论证得 $\nabla h_j(\vx^*)^\T\vd=0$。
2. 取 $\phi=f$：当 $k$ 充分大时，点 $\vx^{(k)}$ 接近 $\vx^*$，因此 $f(\vx^{(k)})\ge f(\vx^*)$，同样除以 $t_k$ 得 $\nabla f(\vx^*)^\T\vd\ge0$。
:::

第 2 部分是把{prf:ref}`thm-06-02-first-order`的一阶条件从可行方向推广到切方向。要把它变成关于梯度的命题，需要第 1 部分中的包含关系成为等式。

:::{prf:definition} 约束规范
:label: def-06-03-cq

若 $T(\vx^*)=L(\vx^*)$，则称约束在 $\vx^*\in\Omega$ 处满足 **Abadie 约束规范**（Abadie's constraint qualification；Abadie, 1967）。更一般地，**约束规范**（CQ）是指约束函数在 $\vx^*$ 处满足的任何这样的条件：在该条件下，对每个光滑目标函数，只要 $\vx^*$ 是其局部极小点，就存在 KKT 乘子。
:::

*（约束规范有时被表述为“对有效的 $i$ 满足 $\nabla g_i(\vx^*)^\T\vd\le0$ 的每个 $\vd$ 都是可行方向”。正如圆盘所示，这种表述对每个严格弯曲的约束都不成立，即使梯度线性无关也是如此；切锥才是正确的对象。）*

直接检验 $T(\vx^*)=L(\vx^*)$ 很少是方便的。下面的条件容易检验，并且每一个都蕴含它。

:::{prf:theorem} 约束规范的充分条件
:label: thm-06-03-cq-sufficient

下列每一条都蕴含 $T(\vx^*)=L(\vx^*)$。

1. **线性约束**：所有 $g_i$ 和 $h_j$ 都是仿射函数。
2. **线性无关约束规范（LICQ）**：梯度 $\{\nabla g_i(\vx^*):i\in I(\vx^*)\}\cup\{\nabla h_j(\vx^*):j=1,\dots,\ell\}$ 线性无关。
3. **Slater 条件**（Slater, 1950）：每个 $g_i$ 是凸函数，每个 $h_j$ 是仿射函数，且存在点 $\hat{\vx}$，使得对所有 $i$ 有 $g_i(\hat{\vx})<0$ 且 $\vh(\hat{\vx})=\vzero$。（只有非线性的 $g_i$ 才需要严格不等式。）
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

我们在每种情形下证明 $L(\vx^*)\subseteq T(\vx^*)$；反向包含关系就是{prf:ref}`lem-06-03-tangent`的第 1 部分。

1. 对于仿射约束，$g_i(\vx^*+\alpha\vd)=g_i(\vx^*)+\alpha\nabla g_i^\T\vd$ 和 $h_j(\vx^*+\alpha\vd)=\alpha\nabla h_j^\T\vd$ 精确成立。若 $\vd\in L(\vx^*)$，则有效不等式约束和等式约束对所有 $\alpha\ge0$ 成立，而非有效约束对充分小的 $\alpha$ 成立，与{prf:ref}`prop-06-02-feasible-directions`相同。因此每个 $\vd\in L(\vx^*)$ 都是可行方向，从而是切方向：$L(\vx^*)\subseteq T(\vx^*)$，再由{prf:ref}`lem-06-03-tangent`的第 1 部分得等式。
2. 给定 $\vd\in L(\vx^*)$，把有效约束组成的方程组（视为等式）应用隐函数定理，可以得到一条速度为 $\vd$ 的光滑可行曲线；细节见 Nocedal and Wright（2006，Lemma 12.2）。
3. 令 $\vd^0=\hat{\vx}-\vx^*$。对于有效的 $g_i$，由凸函数的梯度不等式（{prf:ref}`thm-02-02-first-order`）得 $\nabla g_i(\vx^*)^\T\vd^0\le g_i(\hat{\vx})-g_i(\vx^*)=g_i(\hat{\vx})$，它对非线性的 $g_i$ 为 $<0$，对仿射的 $g_i$ 为 $\le0$；又因为 $h_j$ 是仿射的，$\nabla h_j(\vx^*)^\T\vd^0=h_j(\hat{\vx})-h_j(\vx^*)=0$。因此，给定 $\vd\in L(\vx^*)$ 和 $\varepsilon>0$，方向 $\vd^\varepsilon=\vd+\varepsilon\vd^0$ 对有效的非线性 $g_i$ 满足 $\nabla g_i(\vx^*)^\T\vd^\varepsilon<0$，所以这些函数沿该方向从值 $0$ 开始减小（{prf:ref}`lem-06-02-descent`）；有效的仿射 $g_i$ 和 $h_j$ 分别改变 $\alpha\nabla g_i(\vx^*)^\T\vd^\varepsilon\le0$ 和 $\alpha\nabla h_j(\vx^*)^\T\vd^\varepsilon=0$；非有效的 $g_i$ 在小步长下保持为负。因此 $\vd^\varepsilon$ 是可行方向，属于 $T(\vx^*)$。当 $\varepsilon\downarrow0$ 时 $\vd^\varepsilon\to\vd$，又因为切锥是闭的（用对角线方法选取序列可知，切方向的极限仍是切方向），所以 $\vd\in T(\vx^*)$。另见 Bazaraa et al.（2006，Section 5.2）。
:::

满足 LICQ 的点称为**正则点**（regular point）。LICQ 是实践中最常用的条件；它在{prf:ref}`ex-06-03-cubic`中不成立，因为 $\nabla g_1(\vx^*)=\vzero$，并且只要有效约束多于 $n$ 个，它就不成立。线性约束完全无需检验，这涵盖了 LP、二次规划（quadratic programming, QP）以及{prf:ref}`ex-06-03-kkt-example`这样的问题。Slater 条件是整个问题的性质而不是某个点的性质，是凸规划的标准假设（{prf:ref}`thm-02-03-slater`）。还存在更弱的条件，例如 Mangasarian–Fromovitz 条件（Mangasarian and Fromovitz, 1967：存在 $\vd$，使得对 $i\in I(\vx^*)$ 有 $\nabla g_i(\vx^*)^\T\vd<0$，且 $\nabla h_j(\vx^*)^\T\vd=0$，其中各 $\nabla h_j(\vx^*)$ 线性无关）；见 Bazaraa et al.（2006，Chapter 5）。

(sec-06-03-kkt)=
## KKT 条件

设 $\vx^*$ 是局部极小点，且约束规范成立，从而 $T(\vx^*)=L(\vx^*)$。由{prf:ref}`lem-06-03-tangent`，

$$
\nabla f(\vx^*)^\T\vd\ge0\qquad\text{对所有 }\vd\in L(\vx^*).
$$ (eq-06-03-no-linear-descent)

{eq}`eq-06-03-no-linear-descent` 中的一切关于 $\vd$ 都是线性的。它说的是线性规划问题

$$
\begin{aligned}
\min_{\vd\in\R^n}\ & \nabla f(\vx^*)^\T\vd\\
\st\ & \nabla g_i(\vx^*)^\T\vd\le0, && i\in I(\vx^*),\\
& \nabla h_j(\vx^*)^\T\vd=0, && j=1,\dots,\ell,
\end{aligned}
$$ (eq-06-03-direction-lp)

的最优值为 $0$，在 $\vd^*=\vzero$ 处取得。（这个 LP 是齐次的：$\vd=\vzero$ 可行且目标值为 $0$；若某个可行的 $\vd$ 的目标值为负，把它乘以 $t\to\infty$ 会使目标值趋于 $-\infty$。因此其最优值要么是 $0$，要么是 $-\infty$，而 {eq}`eq-06-03-no-linear-descent` 说明它是 $0$。）具有有限最优值的 LP 的对偶问题是可行的（{prf:ref}`thm-03-05-strong-duality`）。对每个“$\le$”行取乘子 $\lambda_i\ge0$，对每个等式行取自由乘子 $\mu_j$，{eq}`eq-06-03-direction-lp` 的对偶问题为

$$
\begin{aligned}
\max\ & 0\\
\st\ & \sum_{i\in I(\vx^*)}\lambda_i\nabla g_i(\vx^*)+\sum_{j=1}^{\ell}\mu_j\nabla h_j(\vx^*)\\
& \qquad=-\nabla f(\vx^*),\\
& \lambda_i\ge0, && i\in I(\vx^*),\\
& \mu_j\ \text{无符号限制}, && j=1,\dots,\ell.
\end{aligned}
$$

（自由的原变量 $d_1,\dots,d_n$ 在对偶问题中给出等式约束，每个坐标一个。）对偶问题可行恰好就是乘子存在。同样的结论也可以直接由 Farkas 引理的一个变形得到（Farkas, 1902；{prf:ref}`thm-03-05-farkas`，其中把 $\vmu=\vmu^+-\vmu^-$ 拆成非负部分）：在关于未知量 $\vd$ 的方程组

$$
\begin{aligned}
\text{(I)}\quad&\nabla f(\vx^*)^\T\vd<0,\\
&\vd\in L(\vx^*)
\end{aligned}
$$

与关于未知量 $\lambda_i$、$\mu_j$ 的方程组

$$
\begin{aligned}
\text{(II)}\quad&\nabla f(\vx^*)+\sum_{i\in I(\vx^*)}\lambda_i\nabla g_i(\vx^*)\\
&\qquad+\sum_{j=1}^{\ell}\mu_j\nabla h_j(\vx^*)=\vzero,\\
&\lambda_i\ge0, \qquad i\in I(\vx^*)
\end{aligned}
$$

中，恰有一个有解，而 {eq}`eq-06-03-no-linear-descent` 说明 (I) 无解。最后，对非有效约束令 $\lambda_i=0$，就可以让求和遍及所有 $i$，代价是增加条件 $\lambda_ig_i(\vx^*)=0$。

:::{prf:theorem} Karush–Kuhn–Tucker 必要条件
:label: thm-06-03-kkt

设 $\vx^*$ 是 {eq}`eq-06-01-nlp` 的局部极小点，且在该点约束规范成立（例如线性约束、LICQ 或 Slater 条件）。则存在乘子 $\vlambda^*\in\R^m$ 和 $\vmu^*\in\R^\ell$，使得

$$
\nabla f(\vx^*)+\sum_{i=1}^m\lambda_i^*\nabla g_i(\vx^*)+\sum_{j=1}^{\ell}\mu_j^*\nabla h_j(\vx^*)=\vzero ,
$$ (eq-06-03-kkt-stationarity)

$$
\begin{aligned}
g_i(\vx^*)&\le0, && i=1,\dots,m,\\
h_j(\vx^*)&=0, && j=1,\dots,\ell,
\end{aligned}
$$ (eq-06-03-kkt-primal)

$$
\lambda_i^*\ge0,\quad i=1,\dots,m,
$$ (eq-06-03-kkt-dual)

$$
\lambda_i^*\,g_i(\vx^*)=0,\quad i=1,\dots,m .
$$ (eq-06-03-kkt-cs)

在 LICQ 下乘子是唯一的。
:::

这四组条件分别称为**平稳性条件**（stationarity，或 KKT 方程）、**原始可行性**（primal feasibility）、**对偶可行性**（dual feasibility，即不等式乘子非负）和**互补松弛性**（complementary slackness）。若对点 $\vx^*$ 存在满足 {eq}`eq-06-03-kkt-stationarity` 至 {eq}`eq-06-03-kkt-cs` 的乘子，则称 $\vx^*$ 为 **KKT 点**（KKT point）。LICQ 下的唯一性成立，是因为两个不同的乘子向量将给出 $-\nabla f(\vx^*)$ 用线性无关的有效约束梯度表示的两种表达式。

**等式乘子是自由的。** 等式 $h_j=0$ 也可以写成两个不等式 $h_j\le0$ 和 $-h_j\le0$，乘子分别为 $\mu_j^+\ge0$ 和 $\mu_j^-\ge0$。它们对平稳性条件的贡献为 $\mu_j^+\nabla h_j-\mu_j^-\nabla h_j=(\mu_j^+-\mu_j^-)\nabla h_j$，因此起作用的只是差 $\mu_j=\mu_j^+-\mu_j^-$，它可正可负。这就是条件对 $\mu_j$ 不加符号限制的原因。（上面的推导直接处理等式，因为拆分后的一对不等式具有线性相关的梯度 $\pm\nabla h_j$，永远不满足 LICQ。）

**该定理是必要条件，不是充分条件。** KKT 点可能是局部极大点或鞍点；要断定是极小点，需要二阶条件或凸性（{ref}`sec-06-03-second-order`，{ref}`sec-06-03-convex`）。而没有约束规范时，极小点可能根本不是 KKT 点，如{prf:ref}`ex-06-03-cubic`和{prf:ref}`ex-08-08-cusp`中经典的尖点。

:::{prf:remark} 历史与 Fritz John 条件
:label: rem-06-03-history

这些条件由 Karush（1939）在一篇未发表的硕士论文中导出，Kuhn 和 Tucker（1951）独立地得到了它们，并通过其论文使之广为人知；只有等式约束时，它们就化为拉格朗日乘子法则。没有任何约束规范时，局部极小点总满足较弱的 **Fritz John 条件**（Fritz John conditions；John, 1948）：存在不全为零的 $\lambda_0\ge0$、$\vlambda\ge\vzero$、$\vmu$，使得 $\lambda_0\nabla f+\sum_i\lambda_i\nabla g_i+\sum_j\mu_j\nabla h_j=\vzero$，并且互补松弛性成立。当 $\lambda_0=0$ 时，这些条件对 $f$ 一无所言；当 $\lambda_0>0$ 时，除以 $\lambda_0$ 即得 KKT 乘子。约束规范保证存在 $\lambda_0>0$ 的乘子，而 Mangasarian–Fromovitz 条件恰好是使*每个* Fritz John 乘子向量都满足 $\lambda_0>0$ 的条件（Bazaraa et al., 2006，Sections 4.2–4.3）。
:::

### 拉格朗日函数

{eq}`eq-06-01-nlp` 的**拉格朗日函数**（Lagrangian）为

$$
\begin{aligned}
L(\vx,\vlambda,\vmu)&=f(\vx)+\sum_{i=1}^m\lambda_ig_i(\vx)+\sum_{j=1}^{\ell}\mu_jh_j(\vx)\\
&=f(\vx)+\vlambda^\T\vg(\vx)+\vmu^\T\vh(\vx).
\end{aligned}
$$ (eq-06-03-lagrangian)

它的偏导数简洁地重现了 KKT 方程组：

$$
\begin{aligned}
&\nabla_{\vx}L(\vx^*,\vlambda^*,\vmu^*)=\vzero,\\
&\frac{\partial L}{\partial\mu_j}=h_j(\vx^*)=0,\\
&\frac{\partial L}{\partial\lambda_i}=g_i(\vx^*)\le0,\\
&\lambda_i^*\ge0,\\
&\lambda_i^*\,\frac{\partial L}{\partial\lambda_i}=0 .
\end{aligned}
$$

平稳性条件说明 $\vx^*$ 是拉格朗日函数关于 $\vx$ 的驻点：约束以价格 $\lambda_i$、$\mu_j$ 被“计价”进目标函数，在这些价格下，约束极小点看起来就像一个无约束驻点。对于凸问题，拉格朗日对偶（{prf:ref}`thm-02-03-kkt-sufficient`）使这一点精确化，其中 KKT 乘子是对偶问题的最优解。

### 几何解释

只保留有效的不等式约束时，平稳性条件为

$$
\begin{aligned}
-\nabla f(\vx^*)&=\sum_{i\in I(\vx^*)}\lambda_i^*\nabla g_i(\vx^*),\\
\lambda_i^*&\ge0, \qquad i\in I(\vx^*).
\end{aligned}
$$ (eq-06-03-kkt-cone)

最速下降方向 $-\nabla f(\vx^*)$ 位于由有效约束的外法向量 $\nabla g_i(\vx^*)$ 生成的锥中。如果它指向别处，就会存在某个方向，既使 $f$ 减小，又保持线性化的有效约束成立；而 Farkas 引理说的正是这种可能性被排除了。可以想象一个小球沿 $f$ 的等值线滚下山坡，最后靠在 $\Omega$ 的墙壁上静止下来：“力”$-\nabla f$ 被它所接触的墙壁提供的非负反作用力 $\lambda_i\nabla g_i$ 所平衡，而它没有接触的墙壁不施加力（对非有效的 $i$，$\lambda_i=0$）。{numref}`fig-06-03-kkt-geometry`显示了有两个弯曲有效约束的情形。

```{figure} figures/06-03-kkt-geometry.svg
:name: fig-06-03-kkt-geometry
:width: 560px
:class: dark-light
:alt: 两条曲线在原点相交，围成透镜形可行域，左侧有第三条非有效直线，f 的椭圆形等值线中心在右侧，原点处画出两个有效约束的梯度，以及位于它们之间、指向右方的 f 的负梯度

KKT 条件的几何意义。可行集位于曲线 $g_1(\vx)=x_2+0.7x_1+0.3x_1^2=0$ 与 $g_2(\vx)=-x_2+0.7x_1+0.3x_1^2=0$ 之间，并位于直线 $g_3=0$ 的右侧；目标函数为 $f(\vx)=(x_1-1.5)^2+2x_2^2$。在 $\vx^*=(0,0)$ 处，$-\nabla f=(3,0)$ 位于由 $\nabla g_1=(0.7,1)$ 和 $\nabla g_2=(0.7,-1)$ 张成的锥中：$-\nabla f=\tfrac{15}{7}\nabla g_1+\tfrac{15}{7}\nabla g_2$。约束 $g_3$ 是非有效约束，其乘子为 $0$。
```

在{numref}`fig-06-03-kkt-geometry`中，$f$ 在透镜形区域上的极小点为 $\vx^*=(0,0)$，$f^*=2.25$，$\nabla f(\vx^*)=(-3,0)^\T$。求解 $(3,0)=\lambda_1(0.7,1)+\lambda_2(0.7,-1)$，由第二个分量得 $\lambda_1=\lambda_2$，由第一个分量得 $1.4\lambda_1=3$，因此 $\lambda_1=\lambda_2=15/7\approx2.143$，两者都为正，正如定理所要求的。指向透镜内部的每个方向 $\vd$ 与 $\nabla g_1$ 和 $\nabla g_2$ 的夹角都至少为 $90^\circ$，从而与它们的非负组合 $-\nabla f$ 的夹角也至少为 $90^\circ$；所以 $\nabla f(\vx^*)^\T\vd\ge0$。

同样的图景也解释了{doc}`06-02-optimality-conditions`结尾的结果：在那里 $\vx^*=(0.5,1.5)$ 在三角形 $x_1+x_2\le2$，$\vx\ge\vzero$ 上最小化 $(x_1-1)^2+(x_2-2)^2$：此时 $-\nabla f(\vx^*)=(1,1)^\T=1\cdot\nabla g_1$，$\lambda_1=1$，而符号约束是非有效约束。

### 互补松弛性与乘子的意义

互补松弛性 {eq}`eq-06-03-kkt-cs` 说明，对每个不等式，$\lambda_i^*$ 与 $-g_i(\vx^*)$ 中至少有一个为零：有松弛量的约束乘子为零，而正的乘子迫使约束为有效约束。两者可以同时为零（**退化**约束，即有效但乘子为零）；若这种情况从不发生，则称**严格互补性**（strict complementarity）成立。同样的结构在 LP 中表现为原变量与对偶松弛变量之间的互补松弛性（{prf:ref}`thm-03-05-complementary-slackness`），事实上，LP 的 KKT 条件恰好就是原始可行性、对偶可行性和互补松弛性（{prf:ref}`ex-08-08-kkt-lp`）。

乘子还度量最优值对约束的反应。若把第 $i$ 个约束放宽为 $g_i(\vx)\le b_i$，其中 $b_i>0$ 很小，则在 LICQ、{prf:ref}`thm-06-03-second-order`的二阶充分条件以及严格互补性成立时，摄动问题在 $\vx^*$ 附近有局部极小点 $\vx^*(b_i)$，其值 $f^*(b_i)=f(\vx^*(b_i))$ 可微，且

$$
\frac{\mathrm d f^*}{\mathrm d b_i}\Big|_{b_i=0}=-\lambda_i^* .
$$ (eq-06-03-sensitivity)

大的乘子表明放宽该约束收益显著；零乘子表明该约束的微小变化在一阶意义下无关紧要。这是 LP 中影子价格（{prf:ref}`thm-03-05-shadow-price`，{doc}`03-07-sensitivity-analysis`）在非线性情形下的对应物；{prf:ref}`ex-08-08-multiplier-sensitivity`在一个例子上验证了 {eq}`eq-06-03-sensitivity`（Luenberger and Ye, 2016，Section 11.7）。

(sec-06-03-signs)=
## 符号约定

{prf:ref}`thm-06-03-kkt`是针对约束为 $g_i\le0$ 的最小化问题叙述的。其他形式可以通过改变符号化为这种形式，结果值得记住，因为这里的错误很常见。

```{list-table} 常见问题形式的 KKT 条件。每一行中 $\lambda_i\ge0$，$\mu_j$ 自由，$h_j(\vx)=0$，互补松弛性为 $\lambda_ig_i(\vx)=0$。
:name: tab-06-03-signs
:header-rows: 1

* - 问题
  - 约束
  - 平稳性条件
* - $\min f$
  - $g_i(\vx)\le0$
  - $\nabla f+\sum_i\lambda_i\nabla g_i+\sum_j\mu_j\nabla h_j=\vzero$
* - $\min f$
  - $g_i(\vx)\ge0$
  - $\nabla f-\sum_i\lambda_i\nabla g_i+\sum_j\mu_j\nabla h_j=\vzero$
* - $\max f$
  - $g_i(\vx)\le0$
  - $\nabla f-\sum_i\lambda_i\nabla g_i+\sum_j\mu_j\nabla h_j=\vzero$
* - $\max f$
  - $g_i(\vx)\ge0$
  - $\nabla f+\sum_i\lambda_i\nabla g_i+\sum_j\mu_j\nabla h_j=\vzero$
```

例如，在 $g_i\le0$ 约束下求 $\max f$，就是在 $g_i\le0$ 约束下求 $\min(-f)$，平稳性条件为 $-\nabla f+\sum_i\lambda_i\nabla g_i=\vzero$：在约束极大点处，$+\nabla f$（增长最快的方向）位于外法向量生成的锥中。$\mu_j$ 的符号无关紧要，因为 $\mu_j$ 是自由的。一条简单的规则是：在最优点处，目标函数的梯度必须沿人们希望移动的方向“顶住”有效约束，而乘子度量它顶得有多用力。

有些教材把约束写成 $c_i(\vx)\ge0$，把拉格朗日函数写成 $f-\sum_i\lambda_ic_i$（Nocedal and Wright, 2006），这就是表中的第二行。所有约定给出相同的点，乘子值至多相差一个符号；重要的是始终如一地使用同一种约定。

**变量的符号约束。** 对于在 $g_i\le0$ 和 $\vx\ge\vzero$ 约束下求 $\min f$ 的问题，界约束 $-x_j\le0$ 的乘子可以消去：由平稳性条件，它们就是 $\nabla f+\sum_i\lambda_i\nabla g_i$ 的各分量，于是它们的符号条件和互补松弛条件变为

$$
\begin{aligned}
&\nabla f(\vx^*)+\sum_i\lambda_i\nabla g_i(\vx^*)\ge\vzero,\\
&\vx^*\ge\vzero,\\
&(\vx^*)^\T\Bigl(\nabla f(\vx^*)+\sum_i\lambda_i\nabla g_i(\vx^*)\Bigr)=0,
\end{aligned}
$$

这推广了 {eq}`eq-06-02-orthant`。

(sec-06-03-example)=
## 用有效集求解 KKT 条件

KKT 条件是一个由方程和不等式组成的方程组，而互补松弛条件使它具有组合性质：对每个 $i$，要么 $\lambda_i=0$，要么 $g_i(\vx)=0$。对于小规模问题，可以按照哪些约束有效来分情形讨论。

:::{prf:algorithm} KKT 条件的有效集分情形讨论
:label: alg-06-03-active-set

**输入：** $f$，$g_1,\dots,g_m$，$h_1,\dots,h_\ell$（可微）。

**输出：** 所有 KKT 点及其乘子。

1. **无约束候选点。** 求解 $\nabla f(\vx)+\sum_j\mu_j\nabla h_j(\vx)=\vzero$，$\vh(\vx)=\vzero$（所有 $\lambda_i=0$）。保留满足 $\vg(\vx)\le\vzero$ 的解。
2. **选择有效集。** 对集合 $S\subseteq\{1,\dots,m\}$，令 $i\in S$ 时 $g_i(\vx)=0$，$i\notin S$ 时 $\lambda_i=0$。
3. **求解。** 联立平稳性条件 {eq}`eq-06-03-kkt-stationarity`、$g_i(\vx)=0$（$i\in S$）和 $\vh(\vx)=\vzero$，求解 $\vx$、$\lambda_i$（$i\in S$）和 $\vmu$。
4. **检验。** 若对 $i\notin S$ 有 $g_i(\vx)\le0$，对 $i\in S$ 有 $\lambda_i\ge0$，则接受该解；否则舍弃。
5. **重复。** 对未被推理排除的集合 $S$ 重复第 2–4 步；在找到的 KKT 点中，全局极小点（若它是 KKT 点）具有最小的 $f$。
:::

共有 $2^m$ 个有效集，因此在实践中要利用符号论证、无约束极小点违反某些约束这一观察或几何直观来剪除一些情形。对于凸问题，可以在找到第一个 KKT 点时停止（{ref}`sec-06-03-convex`）。下面的例子同时展示了剪枝推理和完全枚举。

:::{prf:example} 求解含五个约束的 KKT 方程组
:label: ex-06-03-kkt-example

用 KKT 条件求

$$
\begin{aligned}
\min\ & f(\vx)=x_1^2+x_2^2+x_3^2\\
\st\ & g_1(\vx)=2x_1+x_2-5\le0,\\
& g_2(\vx)=x_1+x_3-2\le0,\\
& g_3(\vx)=1-x_1\le0,\\
& g_4(\vx)=2-x_2\le0,\\
& g_5(\vx)=-x_3\le0,
\end{aligned}
$$

的极小点。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in


**KKT 方程组。** 梯度为 $\nabla f=(2x_1,2x_2,2x_3)^\T$，$\nabla g_1=(2,1,0)^\T$，$\nabla g_2=(1,0,1)^\T$，$\nabla g_3=(-1,0,0)^\T$，$\nabla g_4=(0,-1,0)^\T$，$\nabla g_5=(0,0,-1)^\T$。平稳性条件 {eq}`eq-06-03-kkt-stationarity` 按分量写出为

$$
\begin{aligned}
2x_1+2\lambda_1+\lambda_2-\lambda_3&=0,\\
2x_2+\lambda_1-\lambda_4&=0,\\
2x_3+\lambda_2-\lambda_5&=0;
\end{aligned}
$$ (eq-06-03-example-stationarity)

互补松弛条件为

$$
\begin{aligned}
\lambda_1(2x_1+x_2-5)&=0,\\
\lambda_2(x_1+x_3-2)&=0,\\
\lambda_3(1-x_1)&=0,\\
\lambda_4(2-x_2)&=0,\\
\lambda_5x_3&=0;
\end{aligned}
$$

此外还有 $g_i(\vx)\le0$ 和 $\lambda_i\ge0$，$i=1,\dots,5$。*（第二个互补松弛条件涉及 $g_2=x_1+x_3-2$；它有时被误印为 $\lambda_2(x_1+x_2-2)=0$。）*

**分情形讨论。** 无约束极小点 $\vx=\vzero$ 违反 $x_1\ge1$ 和 $x_2\ge2$，因此某些约束是有效的。与其尝试全部 $2^5=32$ 个有效集，不如利用符号。

1. *约束 4 有效。* {eq}`eq-06-03-example-stationarity` 的第二个方程给出 $\lambda_4=2x_2+\lambda_1$。可行性要求 $x_2\ge2$，又 $\lambda_1\ge0$，因此 $\lambda_4\ge4>0$。由互补松弛性 $g_4=0$：$x_2=2$，进而 $\lambda_4=4+\lambda_1$。
2. *约束 3 有效。* 利用 $x_1\ge1$，第一个方程给出 $\lambda_3=2x_1+2\lambda_1+\lambda_2\ge2\cdot1+0+0=2>0$。因此 $g_3=0$：$x_1=1$，且 $\lambda_3=2+2\lambda_1+\lambda_2$。
3. *约束 1 非有效。* 由 $x_1=1$，$x_2=2$ 得 $g_1=2+2-5=-1<0$，因此 $\lambda_1=0$。从而 $\lambda_4=4$，$\lambda_3=2+\lambda_2$。
4. *约束 2 和约束 5。* 由 $x_1=1$ 得 $g_2=x_3-1$，第三个方程给出 $\lambda_5=2x_3+\lambda_2$。若 $\lambda_2>0$，互补松弛性迫使 $x_3=1$，于是 $\lambda_5=2+\lambda_2>0$ 又迫使 $x_3=0$：矛盾。所以 $\lambda_2=0$，$\lambda_5=2x_3$；此时 $\lambda_5x_3=2x_3^2=0$ 给出 $x_3=0$ 和 $\lambda_5=0$。

因此 KKT 方程组恰有一个解：

$$
\begin{aligned}
\vx^*&=(1,2,0),\\
\vlambda^*&=(0,0,2,4,0),\\
f(\vx^*)&=1+4+0=5 .
\end{aligned}
$$

**完全枚举。** 对全部 32 个子集 $S$ 运行{prf:ref}`alg-06-03-active-set`也得到同样的结果。由于所有约束都是线性的且 $f$ 是二次函数，每种情形下第 3 步都是一个线性方程组。{numref}`tab-06-03-active-set-cases`列出了一些有代表性的情形。

```{list-table} 五约束例题有效集分析中的部分情形。“不可行”表示 $S$ 之外的某个约束被违反；“不相容”表示方程 $g_i=0$，$i\in S$ 没有公共解。
:name: tab-06-03-active-set-cases
:header-rows: 1

* - 有效集 $S$
  - 第 3 步的解
  - 结果
* - $\varnothing$
  - $\vx=(0,0,0)$
  - 不可行（$g_3,g_4>0$）
* - $\{3\}$
  - $\vx=(1,0,0)$, $\lambda_3=2$
  - 不可行（$g_4>0$）
* - $\{4\}$
  - $\vx=(0,2,0)$, $\lambda_4=4$
  - 不可行（$g_3>0$）
* - $\{1,4\}$
  - $\vx=(1.5,2,0)$, $\lambda_1=-1.5$, $\lambda_4=2.5$
  - $\lambda_1<0$
* - $\{2,4\}$
  - $\vx=(1,2,1)$, $\lambda_2=-2$, $\lambda_4=4$
  - $\lambda_2<0$
* - $\{3,4\}$
  - $\vx=(1,2,0)$, $\lambda_3=2$, $\lambda_4=4$
  - **KKT 点**
* - $\{3,4,5\}$
  - $\vx=(1,2,0)$, $\lambda_3=2$, $\lambda_4=4$, $\lambda_5=0$
  - **KKT 点**（同一个点）
* - $\{1,3,4\}$
  - $x_1=1$, $x_2=2$，但 $2x_1+x_2=4\ne5$
  - 不相容
```

在 32 个集合中，14 个给出不可行点，8 个给出负乘子，8 个不相容，2 个给出 KKT 点 $(1,2,0)$。两个被接受的集合的区别仅在于是否把退化约束 $g_5$（有效且 $\lambda_5=0$）算作有效约束。

**验证。** 在 $\vx^*$ 处：$\vg(\vx^*)=(-1,-1,0,0,0)\le\vzero$；平稳性条件 $(2,4,0)+2(-1,0,0)+4(0,-1,0)=(0,0,0)$；互补松弛性逐项成立。有效集为 $I(\vx^*)=\{3,4,5\}$，其梯度 $-\ve_1,-\ve_2,-\ve_3$ 线性无关，因此 LICQ 成立，乘子唯一。数值求解器（SLSQP，一个序列二次规划程序）给出同一个点，$f=5.000$。

**全局最优性。** 目标函数严格凸且约束是线性的，因此由{prf:ref}`thm-06-03-convex-kkt`，该 KKT 点是唯一的全局极小点。从几何上看，$\vx^*$ 是多面体中离原点最近的点：$x_1\ge1$、$x_2\ge2$ 和 $x_3\ge0$ 不能同时改进，而这三个界约束的角点 $(1,2,0)$ 恰好满足另外两个约束。

**解释。** $\lambda_3^*=2$ 和 $\lambda_4^*=4$ 是界约束 $x_1\ge1$ 和 $x_2\ge2$ 的边际成本：把第一个下界提高到 $1+\delta$，最优值提高到 $(1+\delta)^2+4$，其在 $\delta=0$ 处的导数为 $2$；把第二个下界提高到 $2+\delta$，最优值为 $1+(2+\delta)^2$，导数为 $4$，这与 {eq}`eq-06-03-sensitivity` 一致。约束 $x_3\ge0$ 是有效约束，但没有代价（$\lambda_5^*=0$），因为即使没有它，$x_3=0$ 也是最优的。
:::

(sec-06-03-second-order)=
## 二阶条件

在 KKT 点处一阶项相互平衡，与{doc}`06-02-optimality-conditions`一样，由曲率来区分极小点和鞍点。对于弯曲的约束需要作两处修改：曲率用**拉格朗日函数**的黑塞矩阵来度量，它包含了按乘子加权的约束曲率；并且只在一阶项无法作出判定的方向上检验曲率。

:::{prf:definition} 临界锥
:label: def-06-03-critical-cone

设 $(\vx^*,\vlambda^*,\vmu^*)$ 满足 KKT 条件。**临界锥**（critical cone）为

$$
\begin{aligned}
C(\vx^*,\vlambda^*)=\bigl\{\vd\in L(\vx^*):\ &\nabla g_i(\vx^*)^\T\vd=0\\
&\text{对所有满足 }\lambda_i^*>0\text{ 的 }i\in I(\vx^*)\bigr\}.
\end{aligned}
$$ (eq-06-03-critical-cone)
:::

对于 $\vd\in L(\vx^*)$，由平稳性条件得 $\nabla f(\vx^*)^\T\vd=-\sum_{i\in I}\lambda_i^*\nabla g_i(\vx^*)^\T\vd\ge0$（等式约束项因 $\nabla h_j(\vx^*)^\T\vd=0$ 而消失），等号成立当且仅当对每个满足 $\lambda_i^*>0$ 的 $i$ 有 $\nabla g_i^\T\vd=0$。因此临界锥由满足 $\nabla f(\vx^*)^\T\vd=0$ 的线性化可行方向组成：这是{prf:ref}`thm-06-02-second-order-necessary`条件 2 的类比。

:::{prf:theorem} 约束问题的二阶条件
:label: thm-06-03-second-order

设 $f$、$g_i$、$h_j$ 二阶连续可微，记 $\nabla^2_{\vx\vx}L^*=\nabla^2f(\vx^*)+\sum_i\lambda_i^*\nabla^2g_i(\vx^*)+\sum_j\mu_j^*\nabla^2h_j(\vx^*)$。

1. （*必要性*）若 $\vx^*$ 是 LICQ 成立的局部极小点，$(\vlambda^*,\vmu^*)$ 是其（唯一的）乘子，则对所有 $\vd\in C(\vx^*,\vlambda^*)$ 有 $\vd^\T\nabla^2_{\vx\vx}L^*\vd\ge0$。
2. （*充分性*）若 $\vx^*$ 可行，$(\vx^*,\vlambda^*,\vmu^*)$ 满足 KKT 条件，且

   $$
   \vd^\T\nabla^2_{\vx\vx}L(\vx^*,\vlambda^*,\vmu^*)\,\vd>0\qquad\text{对所有 }\vd\in C(\vx^*,\vlambda^*),\ \vd\ne\vzero,
   $$ (eq-06-03-sosc)

   则 $\vx^*$ 是严格局部极小点。这里不需要约束规范。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

我们证明第 2 部分；第 1 部分见 Nocedal and Wright（2006，Theorem 12.5）。假设 $\vx^*$ 不是严格局部极小点。则存在可行点 $\vx^{(k)}\ne\vx^*$，使得 $\vx^{(k)}\to\vx^*$ 且 $f(\vx^{(k)})\le f(\vx^*)$。令 $t_k=\norm{\vx^{(k)}-\vx^*}$，$\vd^{(k)}=(\vx^{(k)}-\vx^*)/t_k$；取子序列，可设 $\vd^{(k)}\to\vd$，$\norm{\vd}=1$，并且由{prf:ref}`lem-06-03-tangent`有 $\vd\in T(\vx^*)\subseteq L(\vx^*)$。

*$\vd$ 是临界方向。* 由 $f(\vx^{(k)})-f(\vx^*)\le0$，除以 $t_k$ 并令 $k\to\infty$ 得 $\nabla f(\vx^*)^\T\vd\le0$。又因为在 $L(\vx^*)$ 上 $\nabla f(\vx^*)^\T\vd\ge0$（已在定理之前证明），所以 $\nabla f(\vx^*)^\T\vd=0$，$\vd\in C(\vx^*,\vlambda^*)$。

*曲率。* 记 $L(\vx)$ 为 $L(\vx,\vlambda^*,\vmu^*)$。由于 $\lambda_i^*\ge0$，$g_i(\vx^{(k)})\le0$，$h_j(\vx^{(k)})=0$，再利用 $\vx^*$ 处的互补松弛性，

$$
L(\vx^{(k)})\le f(\vx^{(k)})\le f(\vx^*)=L(\vx^*).
$$

在 $\nabla_{\vx}L=\vzero$ 的点 $\vx^*$ 处对 $L$ 应用泰勒定理，得 $L(\vx^{(k)})=L(\vx^*)+\tfrac12t_k^2(\vd^{(k)})^\T\nabla^2_{\vx\vx}L(\vy^{(k)})\vd^{(k)}$，其中 $\vy^{(k)}$ 位于 $\vx^*$ 与 $\vx^{(k)}$ 之间。因此 $(\vd^{(k)})^\T\nabla^2_{\vx\vx}L(\vy^{(k)})\vd^{(k)}\le0$，取极限得：对一个非零的临界方向 $\vd$ 有 $\vd^\T\nabla^2_{\vx\vx}L^*\vd\le0$，与 {eq}`eq-06-03-sosc` 矛盾。
:::

在无约束情形下 $L=f$，$C=\R^n$，该定理化为{prf:ref}`thm-06-02-second-order-necessary`和{prf:ref}`thm-06-02-second-order-sufficient`。只有线性约束时 $\nabla^2_{\vx\vx}L=\nabla^2f$；在{prf:ref}`ex-06-03-kkt-example`中，$\nabla^2f=2\mI\succ\vzero$，充分条件显然成立。当约束弯曲时，约束的曲率就起作用了，如下例所示。

**例：圆周上的线性目标函数。** 在 $h(\vx)=x_1^2+x_2^2-2=0$ 的约束下最小化 $f(\vx)=x_1+x_2$。平稳性条件 $(1,1)+\mu(2x_1,2x_2)=\vzero$ 给出 $x_1=x_2=-1/(2\mu)$，约束给出 $2/(4\mu^2)=2$，因此 $\mu=\pm\tfrac12$。两个 KKT 点为

- $\vx=(-1,-1)$，$\mu=\tfrac12$，$f=-2$；
- $\vx=(1,1)$，$\mu=-\tfrac12$，$f=2$。

LICQ 在两点处都成立（$\nabla h=(2x_1,2x_2)\ne\vzero$）。$f$ 的黑塞矩阵为零，因此仅对 $f$ 进行检验无法区分这两个点。拉格朗日函数的黑塞矩阵为 $\nabla^2_{\vx\vx}L=\nabla^2f+\mu\nabla^2h=2\mu\mI$，临界锥是切线 $\{\vd:\nabla h^\T\vd=0\}=\{t(1,-1)\}$。对于 $\vd=(1,-1)$，$\vd^\T\nabla^2_{\vx\vx}L\,\vd=4\mu$。在 $(-1,-1)$ 处它等于 $2>0$：充分条件成立，该点是严格局部极小点（实际上是全局极小点）。在 $(1,1)$ 处它等于 $-2<0$：必要条件不成立，该点不是局部极小点；它是 $f$ 在圆周上的最大值点。

(sec-06-03-convex)=
## 凸规划的 KKT 条件

对于凸规划，KKT 条件既是必要的也是充分的，而且结论是全局性的。

:::{prf:theorem} 凸规划的 KKT 点是全局极小点
:label: thm-06-03-convex-kkt

设 $f$ 和 $g_1,\dots,g_m$ 是连续可微的凸函数，$h_1,\dots,h_\ell$ 是仿射函数。若 $(\vx^*,\vlambda^*,\vmu^*)$ 满足 KKT 条件 {eq}`eq-06-03-kkt-stationarity`–{eq}`eq-06-03-kkt-cs`，则 $\vx^*$ 是全局极小点。反之，若 Slater 条件成立，则每个全局极小点都是 KKT 点。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

由于 $\lambda_i^*\ge0$ 且 $h_j$ 是仿射的，$\vx\mapsto L(\vx,\vlambda^*,\vmu^*)$ 是凸函数，而平稳性条件说明它的梯度在 $\vx^*$ 处为零。由凸函数的梯度不等式（{prf:ref}`thm-02-02-first-order`），$\vx^*$ 在 $\R^n$ 上最小化 $L(\cdot,\vlambda^*,\vmu^*)$。对任意可行点 $\vx$，

$$
\begin{aligned}
f(\vx)&\ge f(\vx)+\sum_i\lambda_i^*g_i(\vx)+\sum_j\mu_j^*h_j(\vx)\\
&=L(\vx,\vlambda^*,\vmu^*)\\
&\ge L(\vx^*,\vlambda^*,\vmu^*)=f(\vx^*),
\end{aligned}
$$

其中第一个不等式利用了 $\lambda_i^*g_i(\vx)\le0$ 和 $h_j(\vx)=0$，最后一个等式利用了互补松弛性和 $\vh(\vx^*)=\vzero$。逆命题就是以 Slater 条件作为约束规范的{prf:ref}`thm-06-03-kkt`。
:::

充分性不需要约束规范，而等式约束必须是仿射的（即使 $h_j$ 是凸函数，非线性等式也会使 $\Omega$ 非凸）。*（在这个定理中要求 $\vx^*$ 是“正则点”，或允许 $h_j$ 为非线性凸函数的表述，前者要求过强，后者是错误的。）*因此，对于满足 Slater 条件的凸规划，KKT 条件是全局最优性的**充要条件**，并且乘子是拉格朗日对偶问题的最优解（{doc}`02-03-convex-optimization-duality`）。{doc}`08-08-examples-nlp-optimality`中的大多数问题都属于这种情形，这也是{doc}`06-06-constrained-algorithms`中的求解方法对凸问题值得信赖的原因。

## 本章回顾

### 小结

- 只有有效约束 $I(\vx^*)$ 限制局部移动；可行方向满足线性化条件，属于 $L(\vx^*)$（{prf:ref}`prop-06-03-linearised`，{eq}`eq-06-03-linearised-cone`）。
- 切锥 $T(\vx^*)$（{prf:ref}`def-06-03-tangent-cone`）描述真实的局部几何；总有 $T\subseteq L$，并且在局部极小点处，在 $T$ 上有 $\nabla f^\T\vd\ge0$（{prf:ref}`lem-06-03-tangent`）。
- 约束规范保证 $T=L$；线性约束、LICQ 和 Slater 条件都是充分条件（{prf:ref}`thm-06-03-cq-sufficient`）。三次约束（{prf:ref}`ex-06-03-cubic`）说明了没有约束规范时会出什么问题。
- 然后由 Farkas 引理（或对 {eq}`eq-06-03-direction-lp` 应用 LP 对偶）得到 KKT 条件：平稳性条件、原始可行性与对偶可行性、互补松弛性（{prf:ref}`thm-06-03-kkt`），用拉格朗日函数 {eq}`eq-06-03-lagrangian` 可简写为 $\nabla_{\vx}L=\vzero$。
- 从几何上看，$-\nabla f(\vx^*)$ 位于有效约束梯度生成的锥中（{numref}`fig-06-03-kkt-geometry`）；乘子度量最优值的灵敏度（{eq}`eq-06-03-sensitivity`）。
- $\max$ 和 $\ge$ 的符号约定汇总于{numref}`tab-06-03-signs`。
- 有效集分情形讨论（{prf:ref}`alg-06-03-active-set`）可以求解小规模 KKT 方程组；{prf:ref}`ex-06-03-kkt-example`有唯一的 KKT 点 $(1,2,0)$，$\vlambda^*=(0,0,2,4,0)$，$f^*=5$。
- 二阶条件在临界锥上使用 $\nabla^2_{\vx\vx}L$（{prf:ref}`thm-06-03-second-order`）。
- 对于凸规划，每个 KKT 点都是全局极小点；在 Slater 条件下，每个极小点都是 KKT 点（{prf:ref}`thm-06-03-convex-kkt`）。

### 常见错误

- **符号错误。** 混用{numref}`tab-06-03-signs`中的约定，特别是对于最大化问题或“$\ge$”约束，会得到符号错误的乘子，从而否定真正的最优解。应先把问题化为 $\min f$，$g_i\le0$ 的形式。
- **忘记互补松弛性**，或对有效约束强加 $\lambda_i>0$。有效约束的乘子可以是 $\lambda_i=0$（{prf:ref}`ex-06-03-kkt-example`中的 $g_5$）。
- **限制等式乘子的符号。** $\mu_j$ 是自由的。
- **把 KKT 条件当作充分条件。** 没有凸性时，KKT 点可能是极大点（圆周上的点 $(1,1)$）或鞍点。应使用以*拉格朗日函数*（而不是 $f$）的黑塞矩阵表述的二阶条件。
- **忽视约束规范。** 在有效约束梯度线性相关的点处（例如{prf:ref}`ex-06-03-cubic`中梯度为零），极小点可能没有乘子；应单独搜索这样的点。
- **只通过可行方向检验 $T=L$。** 在弯曲的边界上，切方向不是可行方向；约束规范是关于切锥的。


## 参考文献

- Abadie, J. (1967). On the Kuhn–Tucker theorem. In J. Abadie (Ed.),
  *Nonlinear programming* (pp. 19–36). North-Holland.
- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Sections 4.2–4.4,
  Fritz John and KKT conditions, second-order conditions; Sections
  5.1–5.3, the cone of tangents and constraint qualifications.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Section 4.6, Farkas' lemma.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Section
  5.5.3, KKT optimality conditions; Section 5.6, perturbation and
  sensitivity.)
- Farkas, J. (1902). Theorie der einfachen Ungleichungen. *Journal für die
  reine und angewandte Mathematik, 124*, 1–27.
- John, F. (1948). Extremum problems with inequalities as subsidiary
  conditions. In *Studies and essays presented to R. Courant on his 60th
  birthday* (pp. 187–204). Interscience.
- Karush, W. (1939). *Minima of functions of several variables with
  inequalities as side constraints* [Master's thesis, University of
  Chicago].
- Kuhn, H. W., & Tucker, A. W. (1951). Nonlinear programming. In J. Neyman
  (Ed.), *Proceedings of the Second Berkeley Symposium on Mathematical
  Statistics and Probability* (pp. 481–492). University of California
  Press.
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming*
  (4th ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3>
  (Sections 11.1–11.8, constrained minimization conditions, sensitivity and
  inequality constraints.)
- Mangasarian, O. L., & Fromovitz, S. (1967). The Fritz John necessary
  optimality conditions in the presence of equality and inequality
  constraints. *Journal of Mathematical Analysis and Applications, 17*(1),
  37–47.
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Sections
  12.1–12.6, tangent cone, constraint qualifications, first- and
  second-order optimality conditions; Section 12.8, sensitivity.)
- Slater, M. (1950). *Lagrange multipliers revisited: A contribution to
  non-linear programming* (Cowles Commission Discussion Paper: Mathematics
  403). Cowles Commission for Research in Economics.
