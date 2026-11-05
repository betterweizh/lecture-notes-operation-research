# 线性代数基础

## 概述

本书中几乎每一种方法本质上都是线性代数。线性规划是一个带符号约束的线性方程组 $\mA\vx=\vb$；单纯形法通过高斯–若尔当主元变换在它的解之间移动；修正单纯形法把基矩阵的逆存储为一系列主元变换矩阵的乘积；非线性规划的最优性条件则是关于梯度、黑塞矩阵和正定矩阵的命题。本章汇集了后续各章要用到的线性代数知识。先介绍向量、矩阵、内积与线性无关，再讨论用消元法求解 $\mA\vx=\vb$，并特别关注**主元变换**（pivot operation）及其**主元变换矩阵**（pivot matrix）（{prf:ref}`alg-01-01-pivot`、{prf:ref}`prop-01-01-pivot-matrix`、{prf:ref}`prop-01-01-basis-exchange`）。四个基本子空间回答了方程组何时无解、有唯一解或有无穷多解（{prf:ref}`thm-01-01-existence-uniqueness`）。随后讨论逆矩阵、分块矩阵与行列式。本章后半部分讨论方阵：特征值、对称矩阵的谱定理（{prf:ref}`thm-01-01-spectral`）、二次型以及（半）正定性的判别法（{prf:ref}`thm-01-01-eigen-test`、{prf:ref}`thm-01-01-sylvester`、{prf:ref}`thm-01-01-psd-minors`）。本章最后讨论范数、正交性与投影。这是一个参考性的章节：既可以先读，也可以在后面的章节需要某个结论时再来查阅。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 计算 $\R^n$ 中的内积、长度、距离和夹角，并判断一组向量是否线性无关；
- 用高斯消元法和高斯–若尔当消元法求解 $\mA\vx=\vb$，并把每一步消元和每一次主元变换写成左乘一个初等矩阵；
- 在基矩阵的某一列被替换时，以乘积形式更新基矩阵的逆；
- 求矩阵四个基本子空间的基，并判断 $\mA\vx=\vb$ 无解、有唯一解还是有无穷多解；
- 求矩阵（包括分块矩阵）的逆，并运用行列式的主要性质；
- 计算特征值与特征向量，将矩阵对角化，并对对称矩阵应用谱定理；
- 用特征值判别法和子式判别法判断一个二次型是正定、半正定、负定、半负定还是不定的；
- 计算向量范数、矩阵范数以及正交投影。
:::

## 符号说明

| 符号 | 含义 |
|---|---|
| $\ve_i$ | 第 $i$ 个单位向量（第 $i$ 个分量为 1，其余为 0） |
| $\mI$ ($\mI_m$) | （$m$ 阶）单位矩阵 |
| $\mA_{\cdot j}$, $\mA_{i\cdot}$ | $\mA$ 的第 $j$ 列和第 $i$ 行 |
| $\col(\mA)$, $\row(\mA)$ | $\mA$ 的列空间和行空间 |
| $\nullsp(\mA)$, $\nullsp(\mA^\T)$ | $\mA$ 的零空间和左零空间 |
| $\rank\mA$ | $\mA$ 的秩 |
| $\mE$ | 初等矩阵；主元变换矩阵（eta 矩阵） |
| $\mL$, $\mU$ | 下三角因子和上三角因子，$\mA=\mL\mU$ |
| $\bar\mA$, $\bar a_{rs}$ | 一张表及其主元 |
| $\mB$, $\mN$; $\vx_B$, $\vx_N$ | 线性规划的基矩阵和非基矩阵；$\vx$ 的相应部分 |
| $\vpi$ | 单纯形乘子，$\vpi^\T=\vc_B^\T\mB^{-1}$ |
| $\mT$ | 基矩阵的加边逆矩阵 |
| $\det\mA$, $\tr\mA$ | 方阵的行列式和迹 |
| $\lambda_i(\mA)$ | $\mA$ 的第 $i$ 个特征值；$\lambda_{\min}(\mA)$、$\lambda_{\max}(\mA)$ 为对称矩阵 $\mA$ 的最小和最大特征值 |
| $\mV$, $\mD$ | 特征向量矩阵和特征值对角矩阵 |
| $\succ$, $\succeq$ | 正定，半正定 |
| $\Sym^n$, $\Sym^n_+$, $\Sym^n_{++}$ | $n$ 阶对称矩阵、半正定矩阵和正定矩阵的集合 |
| $D_k$ | 第 $k$ 个顺序主子式 |
| $\norm{\cdot}_1$, $\norm{\cdot}_2$, $\norm{\cdot}_\infty$ | 向量范数及其诱导的矩阵范数；$\norm{\cdot}=\norm{\cdot}_2$ |
| $\norm{\cdot}_F$ | Frobenius 范数 |
| $\kappa(\mA)$ | 条件数 $\norm{\mA}\norm{\mA^{-1}}$ |
| $\mP$ | 投影矩阵 |

(sec-01-01-vectors)=
## 向量与矩阵

**向量**（vector）$\vx\in\R^n$ 是由 $n$ 个实数组成的有序数组，这些实数 $x_1,\dots,x_n$ 称为它的**分量**（component）。向量均指列向量；转置 $\vx^\T=(x_1,\dots,x_n)$ 是把同一组数写成一行。为节省篇幅，正文中常把列向量写成 $\vx=(x_1,\dots,x_n)$。长度相同的向量按分量相加，与标量 $\alpha\in\R$ 相乘也按分量进行：$(\vx+\vy)_j=x_j+y_j$，$(\alpha\vx)_j=\alpha x_j$。零向量记为 $\vzero$，全一向量记为 $\vone$，$\ve_i$ 是第 $i$ 个单位向量。对向量而言，$\vx\ge\vy$ 表示对每个 $j$ 都有 $x_j\ge y_j$。

**矩阵**（matrix）$\mA\in\R^{m\times n}$ 是由 $m$ 行 $n$ 列元素 $a_{ij}$ 组成的矩形数表。它的第 $j$ 列为 $\mA_{\cdot j}\in\R^m$，第 $i$ 行为 $\mA_{i\cdot}$，是一个长度为 $n$ 的行向量。转置 $\mA^\T\in\R^{n\times m}$ 的元素为 $(\mA^\T)_{ij}=a_{ji}$。满足 $\mA^\T=\mA$ 的方阵称为**对称矩阵**（symmetric matrix）。若方阵对角线以外的元素全为零，则称为**对角矩阵**（diagonal matrix，记为 $\diag(d_1,\dots,d_n)$）；若对角线以下的元素全为零，则称为**上三角矩阵**（upper triangular matrix）；若对角线以上的元素全为零，则称为**下三角矩阵**（lower triangular matrix）。**迹**（trace）$\tr\mA=\sum_i a_{ii}$ 是对角元素之和。

$\mA\in\R^{m\times n}$ 与 $\mB\in\R^{n\times p}$ 的乘积是一个 $m\times p$ 矩阵，其元素为

$$
(\mA\mB)_{ik}=\sum_{j=1}^n a_{ij}b_{jk}.
$$ (eq-01-01-matrix-product)

只有当 $\mA$ 的列数等于 $\mB$ 的行数时，乘积才有定义。式 {eq}`eq-01-01-matrix-product` 有三种理解方式，会反复用到。

1. **逐个元素。** $(\mA\mB)_{ik}$ 是 $\mA$ 的第 $i$ 行与 $\mB$ 的第 $k$ 列的内积。
2. **逐列。** $\mA\vx=\sum_{j}x_j\mA_{\cdot j}$ 是 $\mA$ 的各列以 $x_j$ 为权重的组合，$\mA\mB$ 的第 $k$ 列是 $\mA\mB_{\cdot k}$。因此，方程组 $\mA\vx=\vb$ 就是寻找一组权重，使 $\mA$ 的各列组合成 $\vb$；这正是{doc}`03-02-lp-geometry`中所用的线性规划的列视角（{eq}`eq-03-02-column-view`）。
3. **逐行。** $\mA\mB$ 的第 $i$ 行是 $\mA_{i\cdot}\mB$，即 $\mB$ 的各行的组合。因此，行变换就是左乘（{ref}`sec-01-01-elimination`）。

矩阵乘法满足结合律 $(\mA\mB)\mC=\mA(\mB\mC)$ 和分配律，但不满足交换律：一般地 $\mA\mB\ne\mB\mA$，即使两个乘积都存在也是如此。乘积的转置要颠倒次序：$(\mA\mB)^\T=\mB^\T\mA^\T$。单位矩阵 $\mI$ 满足 $\mA\mI=\mI\mA=\mA$。迹在循环移位下不变：当两个乘积都是方阵时，$\tr(\mA\mB)=\tr(\mB\mA)$；{doc}`01-02-matrix-calculus`利用这一事实对矩阵求导。

(sec-01-01-inner-products)=
## 内积、长度与夹角

$\vx,\vy\in\R^n$ 的**内积**（inner product，又称点积）是数

$$
\vx^\T\vy=\sum_{j=1}^n x_jy_j=x_1y_1+\dots+x_ny_n.
$$ (eq-01-01-inner-product)

内积是对称的（$\vx^\T\vy=\vy^\T\vx$），且对每个变元都是线性的。$\vx$ 的**欧氏长度**（Euclidean length）为 $\norm{\vx}=\sqrt{\vx^\T\vx}$，$\vx$ 与 $\vy$ 之间的**距离**（distance）是它们之差的长度：

$$
\norm{\vx-\vy}=\sqrt{(x_1-y_1)^2+\dots+(x_n-y_n)^2}.
$$ (eq-01-01-distance)

:::{prf:theorem} 柯西–施瓦茨不等式
:label: thm-01-01-cauchy-schwarz

对所有 $\vx,\vy\in\R^n$，有 $\abs{\vx^\T\vy}\le\norm{\vx}\,\norm{\vy}$，等号成立当且仅当其中一个向量是另一个向量的标量倍。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

若 $\vy=\vzero$，则两边均为零。否则，关于 $t$ 的二次多项式 $\norm{\vx-t\vy}^2=\norm{\vx}^2-2t\,\vx^\T\vy+t^2\norm{\vy}^2$ 对每个实数 $t$ 都非负，因此其判别式不为正：$4(\vx^\T\vy)^2-4\norm{\vx}^2\norm{\vy}^2\le0$。等号成立当且仅当该多项式有实根 $t$，即 $\vx=t\vy$。
:::

柯西–施瓦茨不等式保证了对非零向量，$\vx^\T\vy/(\norm{\vx}\norm{\vy})$ 位于 $[-1,1]$ 中，因此它是某个角的余弦。非零向量 $\vx$ 与 $\vy$ 之间的**夹角**（angle）$\varphi\in[0,\pi]$ 定义为

$$
\cos\varphi=\frac{\vx^\T\vy}{\norm{\vx}\,\norm{\vy}}.
$$ (eq-01-01-angle)

在 $\R^2$ 中，这与初等几何一致：若 $\vx$ 与 $\vy$ 和 $x_1$ 轴的夹角分别为 $\alpha$ 和 $\beta$，则 $\vx=\norm{\vx}(\cos\alpha,\sin\alpha)$，$\vy=\norm{\vy}(\cos\beta,\sin\beta)$，且 $\vx^\T\vy=\norm{\vx}\norm{\vy}(\cos\alpha\cos\beta+\sin\alpha\sin\beta) =\norm{\vx}\norm{\vy}\cos(\alpha-\beta)$。由于式 {eq}`eq-01-01-angle` 中的分母为正，仅凭内积的符号就能判定夹角的类型：

- $\vx^\T\vy>0$：夹角为**锐角**（$\varphi<\pi/2$）；
- $\vx^\T\vy=0$：夹角为**直角**；$\vx$ 与 $\vy$ **正交**（orthogonal），记为 $\vx\perp\vy$；
- $\vx^\T\vy<0$：夹角为**钝角**（$\varphi>\pi/2$）。

这种三分法在最优化中随处可见：若方向 $\vd$ 与梯度成钝角，即 $\nabla f(\vx)^\T\vd<0$，则 $\vd$ 是 $f$ 在 $\vx$ 处的下降方向（{prf:ref}`thm-01-02-steepest-ascent`、{prf:ref}`def-06-02-descent-direction`）。

例如，$\vx=(3,4)$ 与 $\vy=(4,-1)$ 满足 $\vx^\T\vy=12-4=8$，$\norm{\vx}=5$，$\norm{\vy}=\sqrt{17}\approx4.123$，因此 $\cos\varphi=8/(5\sqrt{17})\approx0.3881$，$\varphi\approx67.17^\circ$，是锐角。它们之间的距离为 $\norm{\vx-\vy}=\norm{(-1,5)}=\sqrt{26}\approx5.099$。

(sec-01-01-span)=
## 线性组合、张成与线性无关

向量 $\vx_1,\dots,\vx_k\in\R^n$（带下标的粗体：它们是向量，而不是分量）的**线性组合**（linear combination）是以实数 $\alpha_i$ 为权重的向量 $\alpha_1\vx_1+\dots+\alpha_k\vx_k$。所有线性组合构成的集合称为它们**张成**（span）的空间，记为 $\operatorname{span}\{\vx_1,\dots,\vx_k\}$。当这些向量是矩阵 $\mX=[\vx_1\ \cdots\ \vx_k]$ 的各列时，线性组合就是 $\mX\valpha$，张成的空间就是所有 $\mX\valpha$ 构成的集合。

```{figure} figures/01-01-span-combination.svg
:name: fig-01-01-span-combination
:width: 680px
:class: dark-light
:alt: 左图：平面中两个线性无关的向量 u 和 v，按平行四边形法则作出的组合 1.5u 加 v，以及沿 u 和 v 方向、覆盖整个平面的斜网格线。右图：两个向量 u 和 v 等于负二分之一 u，它们位于同一条过原点的直线上，这条直线就是它们张成的空间

(a) 线性无关的向量 $\vu$ 和 $\vv$ 张成 $\R^2$：每个点都是一个组合 $\alpha\vu+\beta\vv$，可在斜网格上读出；橙色向量 $1.5\vu+\vv$ 由平行四边形法则作出。(b) 线性相关的一对向量 $\vu$、$\vv=-0.5\vu$ 只张成沿 $\vu$ 方向、过原点的直线。
```

:::{prf:definition} 线性无关
:label: def-01-01-independence

若

$$
\alpha_1\vx_1+\dots+\alpha_k\vx_k=\vzero\quad\Longrightarrow\quad
\alpha_1=\dots=\alpha_k=0,
$$

则称向量 $\vx_1,\dots,\vx_k\in\R^n$ **线性无关**（linearly independent），即唯一能得到零向量的组合是平凡组合。否则称它们**线性相关**（linearly dependent）：存在某个非平凡组合等于零，且该组合中权重非零的任一向量都是其余向量的组合。
:::

{numref}`fig-01-01-span-combination`给出了 $\R^2$ 中的两种情形。线性无关容易定义，但需要一个判别方法。把这些向量写成 $\mX$ 的各列，定义说的就是 $\mX\valpha=\vzero$ 只有解 $\valpha=\vzero$，而这个齐次方程组可以用消元法求解。

:::{prf:proposition} 线性无关的判别
:label: prop-01-01-independence-test

设 $\mX=[\vx_1\ \cdots\ \vx_k]\in\R^{n\times k}$。以下命题等价：

1. $\vx_1,\dots,\vx_k$ 线性无关；
2. $\mX\valpha=\vzero$ 只有解 $\valpha=\vzero$；
3. 对 $\mX$ 消元后每一列都有主元，即 $\rank\mX=k$（{ref}`sec-01-01-basis-rank`）。

特别地，$\R^n$ 中多于 $n$ 个的向量必定线性相关；$\R^n$ 中 $n$ 个向量线性无关当且仅当 $\det\mX\ne0$（{ref}`sec-01-01-determinants`）。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

(1)⇔(2) 只是把定义写成矩阵形式。消元不改变 $\mX\valpha=\vzero$ 的解集。没有主元的列对应一个自由变量，令它为 $1$ 即得一个非零解；若每一列都有主元，则回代迫使每个 $\alpha_i=0$。这就证明了 (2)⇔(3)。主元至多有 $n$ 个（每行一个），因此 $k>n$ 列不可能每列都有主元。当 $k=n$ 时，每一列都有主元意味着 $\mX$ 可逆，由{prf:ref}`thm-01-01-invertible`，这等价于 $\det\mX\ne0$。
:::

:::{prf:example} 用消元法检验线性无关
:label: ex-01-01-independence

$\vx_1=(1,0,2)$、$\vx_2=(0,1,1)$ 和 $\vx_3=(2,3,7)$ 是否线性无关？
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

求解 $\mX\valpha=\vzero$，其中 $\mX=[\vx_1\ \vx_2\ \vx_3]$。$\mX$ 的各行为 $(1,0,2)$、$(0,1,3)$ 和 $(2,1,7)$。第 3 行减去第 1 行的 2 倍得 $(0,1,3)$；新的第 3 行再减去第 2 行得 $(0,0,0)$。只出现两个主元，因此 $\rank\mX=2<3$，这些向量线性相关。第三个变量是自由变量；取 $\alpha_3=1$，由第 2 行得 $\alpha_2=-3$，由第 1 行得 $\alpha_1=-2$。验证：$-2\vx_1-3\vx_2+\vx_3=(-2+0+2,\ 0-3+3,\ -4-3+7)=\vzero$。因此 $\vx_3=2\vx_1+3\vx_2$ 位于 $\vx_1$ 和 $\vx_2$ 张成的平面内，且 $\det\mX=0$。
:::

(sec-01-01-basis-rank)=
## 子空间、基、维数与秩

$\R^n$ 的**子空间**（subspace）是对加法和数乘封闭的非空集合 $S$：$\vx,\vy\in S$ 且 $\alpha\in\R$ 蕴含 $\vx+\vy\in S$ 和 $\alpha\vx\in S$。每个子空间都包含 $\vzero$。过原点的直线和平面是 $\R^3$ 的子空间；不过原点的直线则不是（它是仿射集，{prf:ref}`def-02-01-affine-set`）。任意一组向量张成的空间都是子空间。

子空间 $S$ 的**基**（basis）是张成 $S$ 的一组线性无关的向量。这时 $S$ 中的每个向量都能以唯一的方式表示为基向量的组合：若 $\sum\alpha_i\vx_i=\sum\beta_i\vx_i$，则 $\sum(\alpha_i-\beta_i)\vx_i=\vzero$，由线性无关得 $\alpha_i=\beta_i$。$S$ 的所有基所含向量的个数都相同，这个数称为**维数**（dimension）$\dim S$。单位向量 $\ve_1,\dots,\ve_n$ 构成 $\R^n$ 的标准基，因此 $\dim\R^n=n$。

$\mA\in\R^{m\times n}$ 的**列空间**（column space）$\col(\mA)\subseteq\R^m$ 是其各列张成的空间，即所有 $\mA\vx$ 构成的集合。**行空间**（row space）$\row(\mA)=\col(\mA^\T) \subseteq\R^n$ 是其各行张成的空间。

:::{prf:theorem} 行秩等于列秩
:label: thm-01-01-rank

对任意矩阵 $\mA$，有 $\dim\col(\mA)=\dim\row(\mA)$。这个共同的数称为 $\mA$ 的**秩**（rank），记为 $\rank\mA$。它等于消元所产生的主元个数，且 $\rank\mA\le\min\{m,n\}$。
:::

证明的思路是：行变换既不改变行空间，也不改变各列之间的线性关系。在消元所得的阶梯形中（{ref}`sec-01-01-elimination`），非零行线性无关并张成行空间；主元列线性无关，其余每一列都是主元列的组合。两者的个数都等于主元个数，并且*原*矩阵的主元列构成其列空间的一组基（Strang, 2016, Section 3.5）。若 $\rank\mA=m$，则称矩阵 $\mA\in\R^{m\times n}$ **行满秩**（full row rank）；若 $\rank\mA=n$，则称其**列满秩**（full column rank）。标准形线性规划通常假设行满秩，即 $\rank\mA=m$：没有冗余方程（{prf:ref}`asm-03-02-full-rank`、{doc}`03-04-two-phase-big-m`）。

(sec-01-01-elimination)=
## 用消元法求解线性方程组

含 $n$ 个未知数的 $m$ 个线性方程构成的方程组写为 $\mA\vx=\vb$，其中 $\mA\in\R^{m\times n}$，$\vb\in\R^m$。全部信息都包含在**增广矩阵**（augmented matrix）$[\mA\mid\vb]$ 中。下列三种**初等行变换**（elementary row operations）不改变方程组的解集：

1. 把某一行的倍数加到另一行上；
2. 用一个非零标量乘某一行；
3. 交换两行。

每种变换都是可逆的，因此既不会丢失解，也不会增加解。每种变换也都相当于左乘一个**初等矩阵**（elementary matrix）$\mE$，即对单位矩阵施以同一变换所得的矩阵。例如，对三行的情形，“第 2 行减去第 1 行的 $\ell$ 倍”就是

$$
\begin{aligned}
\mE_{21}&=\begin{pmatrix}1&0&0\\-\ell&1&0\\0&0&1\end{pmatrix},\\
\mE_{21}^{-1}&=\begin{pmatrix}1&0&0\\\ell&1&0\\0&0&1\end{pmatrix}.
\end{aligned}
$$ (eq-01-01-elementary-matrix)

它的逆是把这一行加回去。用 $\gamma\ne0$ 乘第 $r$ 行对应 $\diag(1,\dots,\gamma,\dots,1)$，交换两行对应一个**置换矩阵**（permutation matrix），即交换了两行的单位矩阵。

**高斯消元法**（Gaussian elimination）利用第一种变换，逐列在对角线下方制造零。用来消去某一列的元素称为**主元**（pivot），它必须非零。

:::{prf:algorithm} 高斯消元法与回代
:label: alg-01-01-gaussian-elimination

**输入：** $\mA\in\R^{n\times n}$，$\vb\in\R^n$。
**输出：** $\mA\vx=\vb$ 的解 $\vx$，或报告 $\mA$ 是奇异矩阵。

1. **前向消元。** 对 $k=1,\dots,n-1$：若 $a_{kk}=0$，则把第 $k$ 行与某个满足 $a_{ik}\ne0$ 的第 $i>k$ 行交换（若不存在这样的行，则报告 $\mA$ 奇异）。对每个 $i>k$ 的行，计算乘数 $\ell_{ik}=a_{ik}/a_{kk}$，并从第 $i$ 行减去第 $k$ 行的 $\ell_{ik}$ 倍（$\mA$ 和 $\vb$ 中同时进行）。
2. **检查。** 若 $a_{nn}=0$，则报告 $\mA$ 奇异。此时方程组变为 $\mU\vx=\vc$，其中 $\mU$ 是上三角矩阵，主元 $u_{11},\dots,u_{nn}$ 均非零。
3. **回代。** 对 $i=n,n-1,\dots,1$：令 $x_i=\bigl(c_i-\sum_{j>i}u_{ij}x_j\bigr)/u_{ii}$。
:::

乘数并不丢弃：它们构成分解 $\mA=\mL\mU$ 的下三角因子 $\mL$（在不需要交换行时）。在浮点运算中，**部分选主元**（partial pivoting）通过交换行，使主元成为所在列中绝对值最大的元素；这使乘数的绝对值不超过 1，从而控制舍入误差。

:::{prf:example} 消元、初等矩阵与 $\mA=\mL\mU$
:label: ex-01-01-elimination-lu

求解 $\mA\vx=\vb$，其中

$$
\begin{aligned}
\mA&=\begin{pmatrix}2&1&1\\4&3&3\\8&7&9\end{pmatrix},\\
\vb&=\begin{pmatrix}3\\7\\19\end{pmatrix}.
\end{aligned}
$$

把各消元步骤写成初等矩阵，求出 $\mA=\mL\mU$ 的因子 $\mL$ 和 $\mU$，并计算 $\det\mA$。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**前向消元。** 第一个主元是 $a_{11}=2$。乘数 $\ell_{21}=4/2=2$ 和 $\ell_{31}=8/2=4$ 给出行 $(0,1,1\mid1)$ 和 $(0,3,5\mid7)$。第二个主元是 $1$；乘数 $\ell_{32}=3$ 把第 3 行变为 $(0,0,2\mid4)$。写成矩阵形式为 $\mE_{32}\mE_{31}\mE_{21}[\mA\mid\vb]=[\mU\mid\vc]$，其中

$$
\begin{aligned}
\mU&=\begin{pmatrix}2&1&1\\0&1&1\\0&0&2\end{pmatrix},\\
\vc&=\begin{pmatrix}3\\1\\4\end{pmatrix}.
\end{aligned}
$$

**回代。** $x_3=4/2=2$，$x_2=(1-x_3)/1=-1$，$x_1=(3-x_2-x_3)/2=1$，因此 $\vx=(1,-1,2)$。验证：$\mA\vx=(2-1+2,\ 4-3+6,\ 8-7+18)=(3,7,19)=\vb$。

**因子 $\mL$。** 按相反次序撤销各步骤，得 $\mA=\mE_{21}^{-1}\mE_{31}^{-1}\mE_{32}^{-1}\mU=\mL\mU$，而这些逆矩阵的乘积恰好把各乘数收集到对角线下方：

$$
\mL=\begin{pmatrix}1&0&0\\2&1&0\\4&3&1\end{pmatrix}.
$$

主元 $2,1,2$ 构成 $\mU$ 的对角线，且 $\det\mA=2\cdot1\cdot2=4$（{ref}`sec-01-01-determinants`）。存储了 $\mA=\mL\mU$ 之后，每增加一个右端项，只需求解两个三角方程组。
:::

**高斯–若尔当消元法**（Gauss–Jordan elimination）在前向消元之后继续进行：把每个主元行除以其主元，并同时消去每个主元*上方*的元素。所得结果是**简化行阶梯形**（reduced row echelon form, RREF），其中每个主元列都是单位向量。主元列对应的变量是**基**变量（或主元变量），其余的是**自由**变量；RREF 把每个基变量用自由变量表示出来。这正是线性规划的典式（{prf:ref}`def-03-03-canonical-form`）：基变量是单位列对应的变量，令自由（非基）变量为零，即可读出一个基本解。

(sec-01-01-pivot)=
## 主元变换与主元变换矩阵

单纯形法、修正单纯形法和对偶单纯形法（{doc}`03-03-simplex-method`、{doc}`03-06-revised-dual-simplex`）都建立在同一种运算之上：对表中选定的某个元素做一次高斯–若尔当主元变换。

:::{prf:algorithm} 以 $(r,s)$ 元素为主元的主元变换
:label: alg-01-01-pivot

**输入：** 一张表 $\bar\mA\in\R^{m\times q}$（例如 $[\mA\mid\vb]$），主元行 $r$ 和主元列 $s$，**主元**（pivot element）$\bar a_{rs}\ne0$。
**输出：** 第 $s$ 列变为单位向量 $\ve_r$ 的表。

1. **主元行归一化。** 把第 $r$ 行除以 $\bar a_{rs}$：对所有 $j$，$\bar a'_{rj}=\bar a_{rj}/\bar a_{rs}$。
2. **消去第 $s$ 列的其余元素。** 对每个 $i\ne r$ 的行，减去新的第 $r$ 行的 $\bar a_{is}$ 倍：对所有 $j$，$\bar a'_{ij}=\bar a_{ij}-\bar a_{is}\,\bar a_{rj}/\bar a_{rs}$。
:::

第 2 步的更新就是“矩形法则”：新元素 $\bar a'_{ij}$ 等于旧元素减去主元行和主元列中与它构成矩形的两个元素之积再除以主元。一次主元变换是一系列初等行变换（一次数乘和 $m-1$ 次倍加），因此相当于左乘这些初等矩阵的乘积。

:::{prf:proposition} 主元变换矩阵
:label: prop-01-01-pivot-matrix

设 $\vw=\bar\mA_{\cdot s}$ 为主元列，$w_r\ne0$。{prf:ref}`alg-01-01-pivot`的主元变换把 $\bar\mA$ 变为 $\mE\bar\mA$，其中**主元变换矩阵**（pivot matrix，又称 **eta 矩阵**，eta matrix）$\mE$ 是把单位矩阵的第 $r$ 列换成向量 $\veta$ 所得的矩阵：

$$
\begin{aligned}
\mE&=\bigl[\ve_1\ \cdots\ \ve_{r-1}\ \ \veta\ \ \ve_{r+1}\ \cdots\ \ve_m\bigr],\\
\eta_r&=\frac{1}{w_r},\\
\eta_i&=-\frac{w_i}{w_r}\quad (i\ne r).
\end{aligned}
$$ (eq-01-01-eta-matrix)

它的逆是把单位矩阵的第 $r$ 列换成 $\vw$ 本身所得的矩阵。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

$\mE\bar\mA$ 的第 $i$ 行为 $\sum_k e_{ik}\bar\mA_{k\cdot}$。$\mE$ 只有第 $r$ 列与单位矩阵不同，因此 $\mE\bar\mA$ 的第 $r$ 行为 $\eta_r\bar\mA_{r\cdot}=\bar\mA_{r\cdot}/w_r$，第 $i\ne r$ 行为 $\bar\mA_{i\cdot}+\eta_i\bar\mA_{r\cdot}=\bar\mA_{i\cdot}-(w_i/w_r)\bar\mA_{r\cdot}$；这正是主元变换的两步。关于逆矩阵，设 $\mF$ 是把单位矩阵的第 $r$ 列换成 $\vw$ 所得的矩阵。则 $\mE\mF$ 在第 $r$ 列以外与 $\mE$ 相同，而其第 $r$ 列为 $\mE\vw$；由于主元变换把第 $s$ 列（即 $\vw$）变为 $\ve_r$，故 $\mE\vw=\ve_r$。因此 $\mE\mF=\mI$。
:::

$\mE$ 只需存储一列，这就是修正单纯形法保存一列 eta 向量而不保存完整逆矩阵的原因。下一个结论说明了为什么不能以*零*元素为主元，以及基矩阵的某一列被替换时其逆如何变化。

:::{prf:proposition} 替换基矩阵的一列
:label: prop-01-01-basis-exchange

设 $\mB\in\R^{m\times m}$ 可逆，$\mB'$ 是把 $\mB$ 的第 $r$ 列换成向量 $\va\in\R^m$ 所得的矩阵。令 $\vw=\mB^{-1}\va$。则

1. $\mB'$ 可逆当且仅当 $w_r\ne0$，且 $\det\mB'=w_r\det\mB$；
2. 若 $w_r\ne0$，则 $\mB'^{-1}=\mE\mB^{-1}$，其中 $\mE$ 是由 $\vw$ 和 $r$ 构造的主元变换矩阵 {eq}`eq-01-01-eta-matrix`。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

由于 $\mB\vw=\va$ 且 $\mB\ve_i=\mB_{\cdot i}$，矩阵 $\mB'$ 等于 $\mB\mF$，其中 $\mF$ 是把单位矩阵的第 $r$ 列换成 $\vw$ 所得的矩阵。将 $\det\mF$ 按第 $r$ 行展开得 $\det\mF=w_r$（$\mF$ 只在第 $r$ 列与单位矩阵不同），因此 $\det\mB'=\det\mB\cdot w_r$。这就证明了 (1)。由{prf:ref}`prop-01-01-pivot-matrix`，$\mF^{-1}=\mE$，因此 $\mB'^{-1}=\mF^{-1}\mB^{-1}=\mE\mB^{-1}$。
:::

从初始基 $\mB_0=\mI$（松弛变量列）出发反复应用该命题，就得到**逆矩阵的乘积形式**（product form of the inverse）：

$$
\mB_k^{-1}=\mE_k\mE_{k-1}\cdots\mE_1,
$$ (eq-01-01-product-form)

{doc}`03-06-revised-dual-simplex`用它来组织计算（{eq}`eq-03-06-product-form`）。{prf:ref}`prop-01-01-basis-exchange`的第 1 部分从线性代数的角度说明了：在单纯形法中，变量只能从进基列元素非零的行出基，否则新的“基”矩阵将是奇异的。{doc}`03-03-simplex-method`中的比值检验 {eq}`eq-03-03-ratio-test` 在这些行中选出一行，同时保持 $\mB^{-1}\vb\ge\vzero$。

:::{prf:example} 两次主元变换与逆矩阵的乘积形式
:label: ex-01-01-pivot-product-form

考虑方程组 $\mA\vx=\vb$，其中

$$
[\mA\mid\vb]=\left(\begin{array}{ccccc|c}
1&2&1&0&0&8\\
3&1&0&1&0&9\\
1&3&0&0&1&12
\end{array}\right).
$$

$x_3,x_4,x_5$ 的列构成单位矩阵，因此它们构成初始基 $\mB_0=\mI$，对应的基本解为 $(x_3,x_4,x_5)=(8,9,12)$。

先让 $x_1$ 在第 2 行进基，再让 $x_2$ 在第 1 行进基。求 eta 矩阵 $\mE_1$、$\mE_2$ 以及新基矩阵的逆 $\mB_2^{-1}=\mE_2\mE_1$，并验证结果。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**主元变换 1：$x_1$ 在第 2 行进基。** 主元列为 $\vw=(1,3,1)$，主元为 $3$。（比值 $8/1$、$9/3$、$12/1$ 中第 2 行的比值 $3$ 最小，因此右端项（right-hand side, RHS）保持非负；见比值检验 {eq}`eq-03-03-ratio-test`。）由 {eq}`eq-01-01-eta-matrix`，

$$
\mE_1=\begin{pmatrix}1&-\tfrac13&0\\0&\tfrac13&0\\0&-\tfrac13&1\end{pmatrix},
$$

$\mE_1[\mA\mid\vb]$ 就是{numref}`tab-01-01-pivot-tableaux`中的表（中间一块）：基为 $(x_3,x_1,x_5)$，取值为 $(5,3,9)$。

**主元变换 2：$x_2$ 在第 1 行进基。** $x_2$ 当前的列为 $\vw=(\tfrac53,\tfrac13,\tfrac83)$；比值 $5/\tfrac53=3$、$3/\tfrac13=9$、$9/\tfrac83=\tfrac{27}{8}$ 选出第 1 行，主元为 $\tfrac53$。于是

$$
\begin{aligned}
\mE_2&=\begin{pmatrix}\tfrac35&0&0\\-\tfrac15&1&0\\-\tfrac85&0&1\end{pmatrix},\\
\mB_2^{-1}&=\mE_2\mE_1=\begin{pmatrix}\tfrac35&-\tfrac15&0\\-\tfrac15&\tfrac25&0\\-\tfrac85&\tfrac15&1\end{pmatrix}.
\end{aligned}
$$

```{list-table} 每次主元变换前后的表 $[\mA\mid\vb]$。主元以粗体表示。
:name: tab-01-01-pivot-tableaux
:header-rows: 1

* - 基变量
  - $x_1$
  - $x_2$
  - $x_3$
  - $x_4$
  - $x_5$
  - RHS
* - $x_3$
  - 1
  - 2
  - 1
  - 0
  - 0
  - 8
* - $x_4$
  - **3**
  - 1
  - 0
  - 1
  - 0
  - 9
* - $x_5$
  - 1
  - 3
  - 0
  - 0
  - 1
  - 12
* - $x_3$
  - 0
  - $\mathbf{\tfrac53}$
  - 1
  - $-\tfrac13$
  - 0
  - 5
* - $x_1$
  - 1
  - $\tfrac13$
  - 0
  - $\tfrac13$
  - 0
  - 3
* - $x_5$
  - 0
  - $\tfrac83$
  - 0
  - $-\tfrac13$
  - 1
  - 9
* - $x_2$
  - 0
  - 1
  - $\tfrac35$
  - $-\tfrac15$
  - 0
  - 3
* - $x_1$
  - 1
  - 0
  - $-\tfrac15$
  - $\tfrac25$
  - 0
  - 2
* - $x_5$
  - 0
  - 0
  - $-\tfrac85$
  - $\tfrac15$
  - 1
  - 1
```

**验证。** 新基矩阵按行的次序收集 $x_2,x_1,x_5$ 的原始列：

$$
\mB_2=\begin{pmatrix}2&1&0\\1&3&0\\3&1&1\end{pmatrix},
$$

确实有 $\mB_2^{-1}\mB_2=\mI$，且 $\mB_2^{-1}\vb=(3,2,1)$，正是最后一块中 $x_2,x_1,x_5$ 的取值。由{prf:ref}`prop-01-01-basis-exchange`，$\det\mB_2=\det\mI\cdot3\cdot\tfrac53=5$。
:::

如{prf:ref}`ex-01-01-pivot-product-form`的解答所示，$\mB_2^{-1}$ 出现在{numref}`tab-01-01-pivot-tableaux`最后一块中 $x_3,x_4,x_5$ 的列下，而这些列在开始时是单位矩阵。每次主元变换都把整张表乘以相应的 $\mE$，因此开始时为 $\mI$ 的那些列最终变为 $\mE_2\mE_1\mI=\mB_2^{-1}$。从松弛变量列中读出 $\mB^{-1}$（以及随之得到的单纯形乘子）是{doc}`03-05-lp-duality`（{ref}`sec-03-05-reading-duals`）和{doc}`03-07-sensitivity-analysis`中反复使用的技巧。

(sec-01-01-four-subspaces)=
## 四个基本子空间与 $\mA\vx=\vb$ 的解

除了 $\col(\mA)$ 和 $\row(\mA)$ 之外，$\mA\in\R^{m\times n}$ 还有两个子空间：**零空间**（null space）$\nullsp(\mA)=\{\vx\in\R^n:\mA\vx=\vzero\}$ 和**左零空间**（left null space）$\nullsp(\mA^\T)=\{\vy\in\R^m:\mA^\T\vy=\vzero\}$。记 $r=\rank\mA$，它们的维数为

$$
\begin{aligned}
\dim\row(\mA)&=\dim\col(\mA)=r,\\
\dim\nullsp(\mA)&=n-r,\\
\dim\nullsp(\mA^\T)&=m-r.
\end{aligned}
$$ (eq-01-01-rank-nullity)

中间的等式就是**秩–零化度定理**（rank–nullity theorem）：消元后剩下 $r$ 个主元变量和 $n-r$ 个自由变量；令某个自由变量为 1、其余自由变量为 0，就得到零空间的一个基向量。左零空间是 $\mA^\T$ 的零空间，因此同样的计数在把 $n$ 换成 $m$ 后成立。

这些子空间成对正交。若 $\mA\vx=\vzero$，则 $\mA$ 的每一行与 $\vx$ 的内积都为零，因此在 $\R^n$ 中 $\nullsp(\mA)$ 与 $\row(\mA)$ 正交。同理，在 $\R^m$ 中 $\nullsp(\mA^\T)$ 与 $\col(\mA)$ 正交。由于维数之和分别为 $n$ 和 $m$，每一对子空间都把所在空间分解开：$\R^n$ 中每个向量都是一个行空间向量与一个零空间向量之和，$\R^m$ 中每个向量都是一个列空间向量与一个左零空间向量之和（Strang, 2016, Sections 4.1–4.2）。

:::{prf:theorem} 解的存在性与唯一性
:label: thm-01-01-existence-uniqueness

设 $\mA\in\R^{m\times n}$，$\vb\in\R^m$。

1. **存在性。** $\mA\vx=\vb$ 有解，当且仅当 $\vb\in\col(\mA)$，当且仅当 $\rank[\mA\mid\vb]=\rank\mA$，当且仅当对每个 $\vy\in\nullsp(\mA^\T)$ 都有 $\vy^\T\vb=0$。
2. **结构。** 若 $\vx_p$ 是一个解，则全部解的集合为 $\vx_p+\nullsp(\mA)=\{\vx_p+\vz:\mA\vz=\vzero\}$。
3. **唯一性。** 若解存在，则解唯一当且仅当 $\nullsp(\mA)=\{\vzero\}$，即 $\rank\mA=n$。否则有无穷多个解。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

(1) $\mA\vx$ 是各列的组合，因此解存在当且仅当 $\vb$ 是各列的组合。把 $\vb$ 添加到 $\mA$ 之后秩不增加，当且仅当 $\vb$ 线性依赖于各列。最后一个条件成立是因为 $\col(\mA)$ 恰由与 $\nullsp(\mA^\T)$ 正交的向量组成（即上面的正交对）。(2) 若 $\mA\vx_p=\vb$ 且 $\mA\vx=\vb$，则 $\mA(\vx-\vx_p)=\vzero$；反之 $\mA(\vx_p+\vz)=\vb+\vzero$。(3) 由 (2) 推出；零空间中的非零向量 $\vz$ 给出无穷多个解 $\vx_p+t\vz$，$t\in\R$。
:::

第 1 部分写成“要么 $\mA\vx=\vb$ 有解，要么存在某个 $\vy$ 满足 $\mA^\T\vy=\vzero$ 且 $\vy^\T\vb\ne0$，二者恰有一个成立”，就是 **Farkas 引理**（Farkas' lemma）的等式形式；需要 $\vx\ge\vzero$ 的不等式形式在{doc}`03-05-lp-duality`中证明（{prf:ref}`thm-03-05-farkas`）。第 2 部分说明了为什么当 $n>m$ 时，标准形线性规划可以从 $\mA\vx=\vb$ 的整个仿射集形式的解中进行选择。

对于**方阵**，这些情形归结为两种。

:::{prf:theorem} 可逆矩阵
:label: thm-01-01-invertible

对 $\mA\in\R^{n\times n}$，以下命题等价：

1. $\mA$ **可逆**（invertible，即非奇异）：存在矩阵 $\mA^{-1}$ 使 $\mA\mA^{-1}=\mA^{-1}\mA=\mI$；
2. 对每个 $\vb\in\R^n$，$\mA\vx=\vb$ 恰有一个解；
3. 仅当 $\vx=\vzero$ 时 $\mA\vx=\vzero$；
4. $\mA$ 的各列（等价地，各行）线性无关；
5. $\rank\mA=n$，即消元得到 $n$ 个非零主元；
6. $\det\mA\ne0$；
7. $0$ 不是 $\mA$ 的特征值。
:::

不可逆的方阵称为**奇异矩阵**（singular matrix）。对奇异矩阵 $\mA$，$\rank\mA<n$，零空间含有非零向量，而 $\mA\vx=\vb$ **要么无解，要么有无穷多个解**：当 $\vb\notin\col(\mA)$ 时无解，否则有无穷多个解。究竟出现哪种情形，由右端项决定，而不仅仅取决于 $\mA$。

:::{prf:example} 秩为 2 的矩阵的四个子空间
:label: ex-01-01-four-subspaces

设

$$
\mA=\begin{pmatrix}1&2&0&1\\2&4&1&4\\3&6&1&5\end{pmatrix}\in\R^{3\times4}.
$$

求它的四个基本子空间各自的一组基和维数。$\vb$ 取何值时 $\mA\vx=\vb$ 有解？分别对 $\vb=(1,3,4)$ 和 $\vb=(1,0,0)$ 求解。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

第 2 行减去第 1 行的 2 倍、第 3 行减去第 1 行的 3 倍，得到行 $(0,0,1,2)$ 和 $(0,0,1,2)$；第 3 行再减去新的第 2 行得零。RREF 为

$$
\begin{pmatrix}1&2&0&1\\0&0&1&2\\0&0&0&0\end{pmatrix},
$$

主元位于第 1 列和第 3 列，因此 $\rank\mA=2$；$x_2$ 和 $x_4$ 是自由变量。

- **行空间**（维数 2）：由 RREF 的非零行 $(1,2,0,1)$ 和 $(0,0,1,2)$ 张成。
- **列空间**（维数 2）：由 $\mA$ 本身的主元列 $(1,2,3)$ 和 $(0,1,1)$ 张成。第 2 列是第 1 列的两倍，第 4 列是第 1 列加上第 3 列的两倍。
- **零空间**（维数 $4-2=2$）：由 RREF 得 $x_1=-2x_2-x_4$，$x_3=-2x_4$。分别令 $(x_2,x_4)=(1,0)$ 和 $(0,1)$，得到基 $(-2,1,0,0)$ 和 $(-1,0,-2,1)$。正如预期，它们都与行空间的两个基向量正交。
- **左零空间**（维数 $3-2=1$）：$\mA$ 的第 3 行减去第 2 行再减去第 1 行等于零，因此 $\vy=(1,1,-1)$ 满足 $\mA^\T\vy=\vzero$。

由{prf:ref}`thm-01-01-existence-uniqueness`，$\mA\vx=\vb$ 有解当且仅当 $\vy^\T\vb=b_1+b_2-b_3=0$。对 $\vb=(1,3,4)$，条件成立；令自由变量为零得特解 $\vx_p=(1,0,1,0)$（验证：$\mA\vx_p=(1,\,2+1,\,3+1)=\vb$），通解为 $\vx_p+\alpha(-2,1,0,0)+\beta(-1,0,-2,1)$，其中 $\alpha,\beta\in\R$。对 $\vb=(1,0,0)$，$b_1+b_2-b_3=1\ne0$，无解。
:::

(sec-01-01-inverse)=
## 逆矩阵

当 $\mA$ 可逆时，$\mA^{-1}$ 是唯一的，且 $\vx=\mA^{-1}\vb$ 是 $\mA\vx=\vb$ 的解。在计算中，很少仅仅为了求解一个方程组而去求 $\mA^{-1}$；消元法更省计算量，也更精确。不过逆矩阵作为一个对象很重要：$\mB^{-1}$ 出现在矩阵形式单纯形法的每一个公式中。对同阶可逆矩阵，有以下常用规则：

$$
\begin{aligned}
(\mA\mB)^{-1}&=\mB^{-1}\mA^{-1},\\
(\mA^\T)^{-1}&=(\mA^{-1})^\T=:\mA^{-\T},\\
(\alpha\mA)^{-1}&=\alpha^{-1}\mA^{-1}\ (\alpha\ne0).
\end{aligned}
$$ (eq-01-01-inverse-rules)

对满足 $ad-bc\ne0$ 的 $2\times2$ 矩阵，

$$
\begin{pmatrix}a&b\\c&d\end{pmatrix}^{-1}
=\frac{1}{ad-bc}\begin{pmatrix}d&-b\\-c&a\end{pmatrix}.
$$ (eq-01-01-inverse-2x2)

一般地，$\mA^{-1}$ 可以通过**对 $[\mA\mid\mI]$ 做高斯–若尔当消元**求得：把 $\mA$ 变为 $\mI$ 的行变换构成一个满足 $\mE\mA=\mI$ 的矩阵 $\mE$，因此 $\mE=\mA^{-1}$，而同样的变换把右边的块 $\mI$ 变为 $\mE\mI=\mA^{-1}$。这与{prf:ref}`ex-01-01-pivot-product-form`中 $\mB^{-1}$ 出现在松弛变量列下方是同一个现象。对{prf:ref}`ex-01-01-elimination-lu`中的矩阵，在对角线上做三次主元变换得

$$
\mA^{-1}=\begin{pmatrix}\tfrac32&-\tfrac12&0\\-3&\tfrac52&-\tfrac12\\1&-\tfrac32&\tfrac12\end{pmatrix},
$$

且 $\mA^{-1}\vb=(\tfrac92-\tfrac72,\ -9+\tfrac{35}2-\tfrac{19}2,\ 3-\tfrac{21}2+\tfrac{19}2)=(1,-1,2)$，与该例中求得的解一致。

(sec-01-01-block)=
## 分块矩阵

矩阵可以用横线和竖线切分成若干**块**（block，即子矩阵）。只要切分方式相容（左因子的列切分与右因子的行切分一致），分块矩阵就可以像普通矩阵一样相乘，只是用块代替元素：

$$
\begin{pmatrix}\mA_{11}&\mA_{12}\\\mA_{21}&\mA_{22}\end{pmatrix}
\begin{pmatrix}\vx_1\\\vx_2\end{pmatrix}
=\begin{pmatrix}\mA_{11}\vx_1+\mA_{12}\vx_2\\\mA_{21}\vx_1+\mA_{22}\vx_2\end{pmatrix}.
$$ (eq-01-01-block-product)

注意每个块乘积中因子的次序，因为块之间不能交换。

**线性规划的基分块。** 重排各列后，写 $\mA=[\mB\ \ \mN]$，其中 $\mB$ 是 $m\times m$ 的可逆基矩阵，并写 $\vx=(\vx_B,\vx_N)$，$\vc=(\vc_B,\vc_N)$。由 {eq}`eq-01-01-block-product`，$\mA\vx=\vb$ 成为 $\mB\vx_B+\mN\vx_N=\vb$。左乘 $\mB^{-1}$ 得 $\vx_B=\mB^{-1}\vb-\mB^{-1}\mN\vx_N$，代入 $z=\vc_B^\T\vx_B+\vc_N^\T\vx_N$ 得

$$
z=\vc_B^\T\mB^{-1}\vb+\bigl(\vc_N^\T-\vc_B^\T\mB^{-1}\mN\bigr)\vx_N .
$$ (eq-01-01-basis-partition)

这两步分块运算就是典式、检验数和单纯形乘子 $\vpi^\T=\vc_B^\T\mB^{-1}$ 背后的全部代数，{doc}`03-03-simplex-method`将展开讨论。

**分块矩阵的逆。** 对角块可逆的分块三角矩阵可以逐块求逆。对于分块为

$$
\mM=\begin{pmatrix}\mA_{11}&\mA_{12}\\\mA_{21}&\mA_{22}\end{pmatrix}
$$

的方阵，若对角块为方阵且 $\mA_{11}$ 可逆，则 $\mA_{11}$ 的 **Schur 补**（Schur complement）为 $\mS=\mA_{22}-\mA_{21}\mA_{11}^{-1}\mA_{12}$。

:::{prf:proposition} 分块矩阵的逆
:label: prop-01-01-block-inverse

若 $\mA_{11}$ 和 $\mS=\mA_{22}-\mA_{21}\mA_{11}^{-1}\mA_{12}$ 都可逆，则 $\mM$ 可逆，$\det\mM=\det\mA_{11}\cdot\det\mS$，且

$$
\mM^{-1}=\begin{pmatrix}
\mA_{11}^{-1}+\mA_{11}^{-1}\mA_{12}\mS^{-1}\mA_{21}\mA_{11}^{-1} & -\mA_{11}^{-1}\mA_{12}\mS^{-1}\\
-\mS^{-1}\mA_{21}\mA_{11}^{-1} & \mS^{-1}
\end{pmatrix}.
$$ (eq-01-01-block-inverse)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

分块消元给出分解

$$
\mM=\begin{pmatrix}\mI&\vzero\\\mA_{21}\mA_{11}^{-1}&\mI\end{pmatrix}
\begin{pmatrix}\mA_{11}&\vzero\\\vzero&\mS\end{pmatrix}
\begin{pmatrix}\mI&\mA_{11}^{-1}\mA_{12}\\\vzero&\mI\end{pmatrix},
$$

将其乘开即可验证。外侧两个因子是对角块为单位矩阵的分块三角矩阵，因此它们的行列式为 1，其逆只需改变非对角块的符号。按相反次序对三个因子求逆再相乘，即得 {eq}`eq-01-01-block-inverse`；取行列式即得 $\det\mM=\det\mA_{11}\det\mS$。
:::

本书中最重要的特例是修正单纯形法和灵敏度分析中的**加边基矩阵**（bordered basis matrix）（{eq}`eq-03-06-bordered-inverse`、{doc}`03-07-sensitivity-analysis`）。取 $\mA_{11}=\mB$，$\mA_{12}=\vzero$，$\mA_{21}=\vc_B^\T$，$\mA_{22}=1$，则 Schur 补为 $\mS=1$，{eq}`eq-01-01-block-inverse` 化为

$$
\begin{aligned}
\mT&=\begin{pmatrix}\mB&\vzero\\\vc_B^\T&1\end{pmatrix}^{-1}\\
&=\begin{pmatrix}\mB^{-1}&\vzero\\-\vc_B^\T\mB^{-1}&1\end{pmatrix}\\
&=\begin{pmatrix}\mB^{-1}&\vzero\\-\vpi^\T&1\end{pmatrix}.
\end{aligned}
$$ (eq-01-01-bordered-inverse)

这样，一个矩阵同时包含了 $\mB^{-1}$ 和单纯形乘子。

(sec-01-01-determinants)=
## 行列式

方阵的**行列式**（determinant）$\det\mA$ 是一个数，它非零当且仅当 $\mA$ 可逆。行列式由三条性质刻画：$\det\mI=1$；交换两行改变符号；$\det\mA$ 对每一行分别是线性的（其余各行固定时）。其余性质都由此推出（Strang, 2016, Section 5.1）：

1. 有两行相等或有一行为零的矩阵，其行列式为 0；
2. 把某一行的倍数加到另一行上不改变 $\det\mA$，用 $\gamma$ 乘某一行使 $\det\mA$ 乘以 $\gamma$；
3. 三角矩阵的行列式等于其对角元素之积；
4. $\det(\mA\mB)=\det\mA\cdot\det\mB$，$\det\mA^\T=\det\mA$；因此 $\det\mA^{-1}=1/\det\mA$，且关于行的所有性质对列也成立。

由性质 2 和 3，消元法可以计算行列式：$\det\mA=\pm u_{11}\cdots u_{nn}$，即各主元之积，当行交换的次数为奇数时取负号。实际中正是这样计算行列式的；这也证明了 $\det\mA\ne0$ 当且仅当全部 $n$ 个主元都非零，即{prf:ref}`thm-01-01-invertible`的条件 5。

对于小矩阵，**代数余子式展开**（cofactor expansion）很方便。设 $M_{ij}$ 为删去 $\mA$ 的第 $i$ 行和第 $j$ 列所得的 $(n-1)\times(n-1)$ 子矩阵的行列式（称为**子式**，minor）；数 $C_{ij}=(-1)^{i+j}M_{ij}$ 称为 $a_{ij}$ 的**代数余子式**（cofactor）。按任一行 $i$ 展开得

$$
\det\mA=\sum_{j=1}^n a_{ij}C_{ij}.
$$ (eq-01-01-cofactor-expansion)

当 $n=2$ 时，这就是 $ad-bc$。对矩阵 $\begin{pmatrix}2&-1&0\\-1&2&-1\\0&-1&2\end{pmatrix}$，按第 1 行展开得 $2(4-1)-(-1)\bigl((-1)\cdot2-(-1)\cdot0\bigr)+0=6-2=4$。代数余子式还给出逆矩阵的显式表达式

$$
\mA^{-1}=\frac{1}{\det\mA}\,\mC^\T,
$$ (eq-01-01-adjugate)

其中 $\mC=(C_{ij})$ 是代数余子式矩阵；以及**克拉默法则**（Cramer's rule）$x_j=\det\mA^{(j)}/\det\mA$，其中 $\mA^{(j)}$ 是把 $\mA$ 的第 $j$ 列换成 $\vb$ 所得的矩阵。两者对大规模计算都毫无用处，但在证明中很有价值。例如，若 $\mB$ 的元素都是整数且 $\det\mB=\pm1$，则由 {eq}`eq-01-01-adjugate` 知 $\mB^{-1}$ 的元素都是整数，因此对每个整数向量 $\vb$，$\mB^{-1}\vb$ 都是整数向量；这正是网络流问题整数性的根源（{prf:ref}`thm-04-04-integrality`）。

(sec-01-01-eigen)=
## 特征值与特征向量

大多数向量在乘以方阵 $\mA$ 后会改变方向。例外的那些向量只被拉伸。

:::{prf:definition} 特征值与特征向量
:label: def-01-01-eigen

若

$$
\mA\vv=\lambda\vv .
$$ (eq-01-01-eigen)

则称数 $\lambda$（实数或复数）为 $\mA\in\R^{n\times n}$ 的**特征值**（eigenvalue），$\vv\ne\vzero$ 为相应的**特征向量**（eigenvector）。本书其他地方把 $\mA$ 的特征值连同其矩阵一起写成 $\lambda_1(\mA),\dots,\lambda_n(\mA)$，以免与拉格朗日乘子混淆。
:::

式 {eq}`eq-01-01-eigen` 说的是 $(\mA-\lambda\mI)\vv=\vzero$ 且 $\vv\ne\vzero$，因此 $\mA-\lambda\mI$ 必须是奇异的：

$$
p_{\mA}(\lambda)=\det(\mA-\lambda\mI)=0 .
$$ (eq-01-01-characteristic)

**特征多项式**（characteristic polynomial）$p_{\mA}$ 的次数为 $n$，因此按重数计算，$\mA$ 恰有 $n$ 个特征值，它们可能是复数。对每个特征值，特征向量（加上 $\vzero$）构成 $\mA-\lambda\mI$ 的零空间。比较 $p_{\mA}$ 与其因式分解 $\prod_i(\lambda_i(\mA)-\lambda)$ 的系数，得到两个验算式：

$$
\begin{aligned}
\sum_{i=1}^n\lambda_i(\mA)&=\tr\mA,\\
\prod_{i=1}^n\lambda_i(\mA)&=\det\mA .
\end{aligned}
$$ (eq-01-01-trace-det)

特别地，$\mA$ 奇异当且仅当 $0$ 是它的特征值。三角矩阵的特征值就是其对角元素。

**对角化。** 设 $\mA$ 有 $n$ 个线性无关的特征向量 $\vv_1,\dots,\vv_n$，对应的特征值为 $\lambda_1,\dots,\lambda_n$。把它们作为 $\mV$ 的各列，并把特征值放入 $\mD=\diag(\lambda_1,\dots,\lambda_n)$。$n$ 个方程 $\mA\vv_i=\lambda_i\vv_i$ 合起来就是 $\mA\mV=\mV\mD$，而 $\mV$ 可逆，因此

$$
\begin{aligned}
\mA&=\mV\mD\mV^{-1},\\
\mA^k&=\mV\mD^k\mV^{-1}\quad (k=1,2,\dots).
\end{aligned}
$$ (eq-01-01-diagonalisation)

在坐标 $\vy=\mV^{-1}\vx$ 下，矩阵的作用就是对每个坐标分别做伸缩。属于不同特征值的特征向量线性无关，因此有 $n$ 个互不相同特征值的矩阵可以对角化。并非每个矩阵都能对角化：$\begin{pmatrix}1&1\\0&1\end{pmatrix}$ 有二重特征值 $1$，但 $\mA-\mI$ 的秩为 1，因此只有一个线性无关的特征向量。对称矩阵绝不会出现这种情况，这就是下一节的内容。

(sec-01-01-symmetric)=
## 对称矩阵与谱定理

二阶连续可微函数的黑塞矩阵是对称的（{doc}`01-02-matrix-calculus`），二次型的矩阵也是对称的。对这类矩阵，特征值理论再完美不过。

:::{prf:theorem} 谱定理
:label: thm-01-01-spectral

设 $\mA\in\R^{n\times n}$ 是对称矩阵。则 $\mA$ 的所有特征值都是实数，并且存在由 $\mA$ 的特征向量构成的 $\R^n$ 的一组标准正交基。等价地，存在正交矩阵 $\mV$（$\mV^\T\mV=\mI$）和实对角矩阵 $\mD$，使得

$$
\mA=\mV\mD\mV^\T=\sum_{i=1}^n\lambda_i(\mA)\,\vv_i\vv_i^\T .
$$ (eq-01-01-spectral)
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

我们证明支撑该定理的两个事实；特征值有重复时完整的标准正交特征基的存在性可对 $n$ 用归纳法证明（Strang, 2016, Section 6.4）。*特征值为实数：* 设 $\mA\vv=\lambda\vv$，其中 $\vv\ne\vzero$ 为复向量，$\bar{\vv}$ 为其复共轭。则 $\bar{\vv}^\T\mA\vv=\lambda\,\bar{\vv}^\T\vv$。取共轭并转置，利用 $\bar{\mA}=\mA=\mA^\T$，得 $\bar{\vv}^\T\mA\vv=\bar\lambda\,\bar{\vv}^\T\vv$。由于 $\bar{\vv}^\T\vv=\sum\abs{v_j}^2>0$，故 $\lambda=\bar\lambda$ 是实数。*特征向量正交：* 若 $\mA\vv_1=\lambda_1\vv_1$，$\mA\vv_2=\lambda_2\vv_2$，且 $\lambda_1\ne\lambda_2$，则 $\lambda_1\vv_2^\T\vv_1=\vv_2^\T\mA\vv_1=(\mA\vv_2)^\T\vv_1=\lambda_2\vv_2^\T\vv_1$，因此 $(\lambda_1-\lambda_2)\vv_2^\T\vv_1=0$，即 $\vv_1\perp\vv_2$。
:::

{eq}`eq-01-01-spectral` 中的和式把 $\mA$ 写成向其各特征向量方向的秩一投影 $\vv_i\vv_i^\T$ 的组合。每个特征值是沿某个轴的伸缩因子，而这些轴相互垂直：可逆对称矩阵把单位球面映为一个椭球面（{prf:ref}`def-02-01-ellipsoid`），其沿各特征向量方向的半轴长为 $\abs{\lambda_i(\mA)}$。特征值还界定了 $\vx^\T\mA\vx$ 与 $\vx^\T\vx$ 之比。

:::{prf:proposition} 瑞利商的界
:label: prop-01-01-rayleigh

对对称矩阵 $\mA$ 和每个 $\vx\ne\vzero$，

$$
\lambda_{\min}(\mA)\le\frac{\vx^\T\mA\vx}{\vx^\T\vx}\le\lambda_{\max}(\mA),
$$ (eq-01-01-rayleigh)

当 $\vx$ 是 $\lambda_{\min}(\mA)$（$\lambda_{\max}(\mA)$）的特征向量时，左（右）边的等号成立。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

令 $\vy=\mV^\T\vx$，由 {eq}`eq-01-01-spectral` 得 $\vx^\T\mA\vx=\vy^\T\mD\vy=\sum_i\lambda_i(\mA)y_i^2$；又因 $\mV$ 是正交矩阵，$\vx^\T\vx=\vy^\T\mV^\T\mV\vy=\vy^\T\vy=\sum_iy_i^2$。以 $y_i^2/\sum_ky_k^2$ 为权重的 $\lambda_i(\mA)$ 的加权平均介于最小值与最大值之间。
:::

:::{prf:example} $2\times2$ 矩阵的谱分解
:label: ex-01-01-spectral-2x2

设 $\mA=\begin{pmatrix}2&1\\1&2\end{pmatrix}$。

求它的特征值和标准正交特征向量，写出它的谱分解，并利用谱分解计算 $\mA^5$，把 $\vx^\T\mA\vx$ 写成平方项的加权和，以及描述单位圆的像。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

其特征多项式为 $(2-\lambda)^2-1=(\lambda-1)(\lambda-3)$，因此 $\lambda_1(\mA)=3$，$\lambda_2(\mA)=1$（验算：$3+1=4=\tr\mA$，$3\cdot1=3=\det\mA$）。$\mA-3\mI$ 的零空间由 $(1,1)$ 张成，$\mA-\mI$ 的零空间由 $(1,-1)$ 张成；正如谱定理所保证的，它们相互正交。单位化后，特征值 3 对应 $\vv_1=(1,1)/\sqrt2$，特征值 1 对应 $\vv_2=(1,-1)/\sqrt2$，且

$$
\mA=3\vv_1\vv_1^\T+1\,\vv_2\vv_2^\T
=\begin{pmatrix}1.5&1.5\\1.5&1.5\end{pmatrix}
+\begin{pmatrix}0.5&-0.5\\-0.5&0.5\end{pmatrix}.
$$

现在求幂就很直接了：$\mA^5=3^5\vv_1\vv_1^\T+1^5\vv_2\vv_2^\T =\begin{pmatrix}122&121\\121&122\end{pmatrix}$。二次型化为 $\vx^\T\mA\vx=3y_1^2+y_2^2$，其中 $y_i=\vv_i^\T\vx$，这是权重为正的平方和。{numref}`fig-01-01-circle-ellipse`给出了几何图景：单位圆被映为一个椭圆，沿 $\vv_1$ 和 $\vv_2$ 的半轴长分别为 3 和 1。向量 $\vx=(1,0)$ 不是特征向量；它的像 $\mA\vx=(2,1)$ 向 $\vv_1$ 方向偏转，其瑞利商 $\vx^\T\mA\vx/\vx^\T\vx=2$ 介于 1 和 3 之间。

```{figure} figures/01-01-circle-ellipse.svg
:name: fig-01-01-circle-ellipse
:width: 480px
:class: dark-light
:alt: 单位圆及其在行为 2 1 和 1 2 的矩阵作用下的像，即一个倾斜 45 度的椭圆；特征向量 v1 沿长轴被拉伸 3 倍，特征向量 v2 在短轴上保持不变，x 等于 e1 被映为 Ax 等于 (2, 1)

对称矩阵 $\mA=\begin{pmatrix}2&1\\1&2\end{pmatrix}$ 把单位圆（灰色）映为一个椭圆（蓝色）。特征向量只被拉伸：沿长轴 $\mA\vv_1=3\vv_1$，沿短轴 $\mA\vv_2=\vv_2$。一般的向量 $\vx$ 还会被转动：$\mA\vx=(2,1)$。
```
:::

(sec-01-01-quadratic-forms)=
## 二次型与有定性

**二次型**（quadratic form）是以方阵 $\mA$ 定义的函数 $Q(\vx)=\vx^\T\mA\vx=\sum_i\sum_j a_{ij}x_ix_j$。起作用的只是 $\mA$ 的对称部分：由于 $\vx^\T\mA\vx$ 是一个数，它等于自身的转置 $\vx^\T\mA^\T\vx$，因此

$$
\vx^\T\mA\vx=\vx^\T\Bigl(\tfrac12(\mA+\mA^\T)\Bigr)\vx .
$$ (eq-01-01-symmetric-part)

例如，对 $\mA=\begin{pmatrix}2&3\\0&1\end{pmatrix}$ 和 $\begin{pmatrix}2&1.5\\1.5&1\end{pmatrix}$，$2x_1^2+3x_1x_2+x_2^2$ 都等于 $\vx^\T\mA\vx$。按照惯例，以后二次型都用对称矩阵表示。

:::{prf:definition} 有定性
:label: def-01-01-definiteness

对称矩阵 $\mA\in\Sym^n$（$n\times n$ 对称矩阵的集合）及其二次型称为

- **正定**（positive definite）的，记为 $\mA\succ0$，若对每个 $\vx\ne\vzero$ 都有 $\vx^\T\mA\vx>0$；
- **半正定**（positive semidefinite）的，记为 $\mA\succeq0$，若对每个 $\vx$ 都有 $\vx^\T\mA\vx\ge0$；
- **负定**（negative definite）的（**半负定**的，negative semidefinite），若 $-\mA$ 是正定（半正定）的；
- **不定**（indefinite）的，若 $\vx^\T\mA\vx$ 既取正值又取负值。

$\Sym^n_{++}$ 和 $\Sym^n_+$ 分别表示 $n$ 阶正定矩阵和半正定矩阵的集合。
:::

这些类别决定了二次函数的形状；并且通过黑塞矩阵和泰勒定理（{doc}`01-02-matrix-calculus`），它们决定了驻点是极小点、极大点还是鞍点（{prf:ref}`thm-06-02-second-order-sufficient`），以及函数是否为凸函数（{prf:ref}`thm-02-02-second-order`）。{numref}`fig-01-01-quadratic-forms`给出了 $\R^2$ 中的三种典型图形。

```{figure} figures/01-01-quadratic-forms.svg
:name: fig-01-01-quadratic-forms
:width: 720px
:class: dark-light
:alt: 三幅等值线图。左图：围绕原点的一族嵌套椭圆。中图：一族平行直线，其中一条过原点的绿色粗线是二次型为零的地方。右图：双曲线，一对相对的区域中为蓝色，另一对中为橙色虚线，二者由两条过原点的绿线隔开

$\vx^\T\mA\vx$ 的等值线（水平 0.5 到 8 为蓝色，负水平为橙色虚线，零水平为绿色）。(a) $\mA=\begin{pmatrix}2&1\\1&2\end{pmatrix}\succ0$：闭合的椭圆，唯一的零点是 $\vx=\vzero$。(b) $\mA=\begin{pmatrix}1&-1\\-1&1\end{pmatrix} \succeq0$，奇异：$\vx^\T\mA\vx=(x_1-x_2)^2$ 在整条直线 $x_1=x_2$（即 $\mA$ 的零空间）上为零。(c) $\mA=\begin{pmatrix}1&2\\2&1\end{pmatrix}$，不定：鞍形，在一对扇形区域中为正，在另一对中为负。
```

:::{prf:theorem} 特征值判别法
:label: thm-01-01-eigen-test

对称矩阵 $\mA$ 正定当且仅当它的所有特征值都为正，半正定当且仅当所有特征值都非负。它是不定的当且仅当它既有正特征值又有负特征值。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

由谱定理，$\vx^\T\mA\vx=\sum_i\lambda_i(\mA)y_i^2$，其中 $\vy=\mV^\T\vx$，且 $\vy\ne\vzero$ 当且仅当 $\vx\ne\vzero$。若所有 $\lambda_i(\mA)>0$，则对 $\vy\ne\vzero$ 该和式为正。若某个 $\lambda_k(\mA)\le0$，则取特征向量 $\vx=\vv_k$ 得 $\vx^\T\mA\vx=\lambda_k(\mA)\le0$。半正定和不定的情形可同样证明。
:::

当 $n>2$ 时，手工计算特征值很繁琐，而子矩阵的行列式提供了另一种方法。$\mA$ 的**主子矩阵**（principal submatrix）保留行指标与列指标属于同一指标集 $J\subseteq\{1,\dots,n\}$ 的行和列；它的行列式称为**主子式**（principal minor）。**顺序**主子矩阵保留 $J=\{1,\dots,k\}$，它们的行列式就是**顺序主子式**（leading principal minor）$D_1,\dots,D_n$。

:::{prf:theorem} 西尔维斯特判据
:label: thm-01-01-sylvester

**对称**矩阵 $\mA$ 正定当且仅当它的所有顺序主子式都为正，即对 $k=1,\dots,n$ 有 $D_k>0$。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

*必要性。* 若 $\mA\succ0$，则每个顺序主子矩阵 $\mA_k$ 也是正定的（取第 $k$ 个位置之后全为零的 $\vx$），因此它的特征值都为正，而 $D_k=\det\mA_k$ 作为它们的乘积也为正。*充分性。* 进行不交换行的高斯消元。前 $k-1$ 步在左上角 $k\times k$ 块内只做倍加变换，不改变该块的行列式，因此这些步骤之后该块是以 $d_1,\dots,d_k$ 为主元的上三角矩阵，且 $D_k=d_1d_2\cdots d_k$。若所有 $D_k>0$，则每个主元 $d_k=D_k/D_{k-1}$（其中 $D_0=1$）都为正，因此消元不会中断；再由对称性得 $\mA=\mL\mD\mL^\T$，其中 $\mL$ 是对角元为 1 的下三角矩阵，$\mD=\diag(d_1,\dots,d_n)$。由于 $\mL^\T$ 可逆，对 $\vx\ne\vzero$ 有 $\vx^\T\mA\vx=\sum_kd_k(\mL^\T\vx)_k^2>0$。
:::

对于**半正定性**，只看顺序主子式是不够的，必须检验所有主子式。

:::{prf:theorem} 半正定性的主子式判别法
:label: thm-01-01-psd-minors

对称矩阵 $\mA$ 半正定当且仅当它的**所有**主子式（对每个指标集 $J$，而不仅仅是顺序主子式）都非负。它负定当且仅当对 $k=1,\dots,n$ 有 $(-1)^kD_k>0$，即顺序主子式从 $D_1<0$ 开始正负交替。
:::

负定的结论就是把西尔维斯特判据应用于 $-\mA$（$-\mA$ 的 $k$ 阶子式等于 $\mA$ 相应子式的 $(-1)^k$ 倍）。半正定的结论见 Strang (2016, Section 6.5) 或 Magnus and Neudecker (2019, Chapter 1)。对角矩阵 $\diag(0,-1)$ 说明了为什么需要额外的子式：它的顺序主子式 $D_1=0$ 和 $D_2=0$ 都非负，但 $\ve_2^\T\diag(0,-1)\ve_2=-1<0$。$J=\{2\}$ 对应的主子式，即元素 $-1$，揭示了问题所在。

{prf:ref}`thm-01-01-sylvester`证明中的分解也是实际中最好的判别方法。当 $\mA\succ0$ 时，各主元为正，且 $\mA=\mL\mD\mL^\T=(\mL\mD^{1/2})(\mL\mD^{1/2})^\T$。这就是 **Cholesky 分解**（Cholesky factorization）$\mA=\mG\mG^\T$，其中 $\mG$ 是对角元为正的下三角矩阵。它存在当且仅当 $\mA\succ0$；其计算量约为 $n^3/3$ 次算术运算（是高斯消元的一半，远少于计算特征值）；如果分解过程中遇到非正的主元，就证明了 $\mA$ 不是正定的。

:::{prf:example} $3\times3$ 矩阵的顺序主子式
:label: ex-01-01-leading-minors

对称 $3\times3$ 矩阵 $\mA=(a_{ij})$ 的三个顺序主子式为

$$
\begin{aligned}
D_1&=a_{11},\\
D_2&=\det\begin{pmatrix}a_{11}&a_{12}\\a_{21}&a_{22}\end{pmatrix}=a_{11}a_{22}-a_{12}a_{21},\\
D_3&=\det\mA .
\end{aligned}
$$

由{prf:ref}`thm-01-01-sylvester`，$\mA\succ0$ 当且仅当 $D_1>0$、$D_2>0$ 且 $D_3>0$。

用它们判断 $\mA=\begin{pmatrix}2&-1&0\\-1&2&-1\\0&-1&2\end{pmatrix}$ 和 $\mA=\begin{pmatrix}1&2&0\\2&1&0\\0&0&1\end{pmatrix}$ 是否正定或半正定。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**一个正定的例子。** 对 $\mA=\begin{pmatrix}2&-1&0\\-1&2&-1\\0&-1&2\end{pmatrix}$，$D_1=2$，$D_2=4-1=3$，且由{ref}`sec-01-01-determinants`中的代数余子式展开，$D_3=4$。它们都为正，因此 $\mA\succ0$。特征值判别法给出同样的结论：特征值为 $2-\sqrt2\approx0.5858$、$2$ 和 $2+\sqrt2\approx3.414$，其和为 $6=\tr\mA$，其积为 $4=D_3$。消元的主元为 $D_1=2$、$D_2/D_1=\tfrac32$ 和 $D_3/D_2=\tfrac43$，Cholesky 因子为

$$
\begin{aligned}
\mG&=\begin{pmatrix}\sqrt2&0&0\\-\tfrac{1}{\sqrt2}&\sqrt{\tfrac32}&0\\0&-\sqrt{\tfrac23}&\sqrt{\tfrac43}\end{pmatrix}\\
&\approx\begin{pmatrix}1.414&0&0\\-0.7071&1.225&0\\0&-0.8165&1.155\end{pmatrix},
\end{aligned}
$$

其对角元素的平方就是各主元。

**一个不定的例子。** 对 $\mA=\begin{pmatrix}1&2&0\\2&1&0\\0&0&1\end{pmatrix}$，$D_1=1>0$，但 $D_2=1-4=-3<0$，因此 $\mA$ 不是正定的；由于半正定矩阵也不可能有负的 $2\times2$ 顺序主子式，$\mA$ 也不是半正定的。特征值 $-1$、$1$、$3$ 证实了 $\mA$ 是不定的：$\vx=(1,-1,0)$ 给出 $\vx^\T\mA\vx=1-4+1=-2$。
:::

子式判别法只对对称矩阵有效，而且条件“所有顺序主子式 $\ge0$”并不蕴含半正定性。例如，$\begin{pmatrix}1&-3\\0&1\end{pmatrix}$ 的顺序主子式为 $1$ 和 $1$，但 $\vx=(1,1)$ 给出 $\vx^\T\mA\vx=-1$；$\diag(0,-1)$ 的顺序主子式为 $0,0$，却不是半正定的。

(sec-01-01-norms)=
## 范数

$\R^n$ 上的**范数**（norm）是满足下列条件的函数 $\norm{\cdot}$：对所有 $\vx,\vy\in\R^n$ 和 $\alpha\in\R$，$\norm{\vx}\ge0$，且仅当 $\vx=\vzero$ 时等号成立；$\norm{\alpha\vx}=\abs{\alpha}\norm{\vx}$；以及**三角不等式**（triangle inequality）$\norm{\vx+\vy}\le\norm{\vx}+\norm{\vy}$。本书使用三种范数：

$$
\begin{aligned}
\norm{\vx}_1&=\sum_{j=1}^n\abs{x_j},\\
\norm{\vx}_2&=\Bigl(\sum_{j=1}^nx_j^2\Bigr)^{1/2}=\sqrt{\vx^\T\vx},\\
\norm{\vx}_\infty&=\max_{j}\abs{x_j}.
\end{aligned}
$$ (eq-01-01-vector-norms)

$\norm{\cdot}_2$ 的三角不等式由柯西–施瓦茨不等式推出：$\norm{\vx+\vy}_2^2=\norm{\vx}^2+2\vx^\T\vy+\norm{\vy}^2\le (\norm{\vx}+\norm{\vy})^2$。对 $\vx=(3,-4,0)$，$\norm{\vx}_1=7$，$\norm{\vx}_2=5$，$\norm{\vx}_\infty=4$。在 $\R^n$ 中，所有范数在相差常数倍的意义下等价：

$$
\norm{\vx}_\infty\le\norm{\vx}_2\le\norm{\vx}_1\le\sqrt n\,\norm{\vx}_2\le n\norm{\vx}_\infty ,
$$ (eq-01-01-norm-equivalence)

因此序列的收敛性（{prf:ref}`def-02-01-limit-point`）与所用的范数无关。{numref}`fig-01-01-unit-balls`画出了它们的单位球。1-范数和 $\infty$-范数的单位球是多面体，这就是含这些范数的问题能写成线性规划的原因（{ref}`sec-07-02-abs`）。

```{figure} figures/01-01-unit-balls.svg
:name: fig-01-01-unit-balls
:width: 400px
:class: dark-light
:alt: 以原点为中心、相互嵌套的菱形、圆和正方形，分别是 1-范数、2-范数和无穷范数的单位球

$\R^2$ 中单位球 $\{\vx:\norm{\vx}\le1\}$ 的边界：$\norm{\cdot}_1$ 的菱形位于 $\norm{\cdot}_2$ 的圆内，而圆又位于 $\norm{\cdot}_\infty$ 的正方形内，正如 {eq}`eq-01-01-norm-equivalence` 所述。
```

**诱导矩阵范数。** 向量范数诱导出矩阵上的一个范数，即最大伸缩因子：

$$
\norm{\mA}=\max_{\vx\ne\vzero}\frac{\norm{\mA\vx}}{\norm{\vx}} .
$$ (eq-01-01-induced-norm)

由定义，$\norm{\mA\vx}\le\norm{\mA}\norm{\vx}$，且 $\norm{\mA\mB}\le\norm{\mA}\norm{\mB}$。三种诱导范数都有闭式表达式：

$$
\begin{aligned}
\norm{\mA}_1&=\max_j\sum_i\abs{a_{ij}},\\
\norm{\mA}_\infty&=\max_i\sum_j\abs{a_{ij}},\\
\norm{\mA}_2&=\sqrt{\lambda_{\max}(\mA^\T\mA)} ,
\end{aligned}
$$ (eq-01-01-matrix-norms)

分别是最大列绝对值和、最大行绝对值和，以及**谱范数**（spectral norm）。最后一个由{prf:ref}`prop-01-01-rayleigh`推出：$\norm{\mA\vx}_2^2/\norm{\vx}_2^2=\vx^\T(\mA^\T\mA)\vx/\vx^\T\vx$，其最大值为 $\lambda_{\max}(\mA^\T\mA)$。对对称矩阵 $\mA$，$\norm{\mA}_2=\max_i\abs{\lambda_i(\mA)}$。一种非诱导但常常很方便的范数是 **Frobenius 范数**（Frobenius norm）$\norm{\mA}_F=\bigl(\sum_{i,j}a_{ij}^2\bigr)^{1/2}=\sqrt{\tr(\mA^\T\mA)}$。

对 $\mA=\begin{pmatrix}1&-2\\3&4\end{pmatrix}$：各列绝对值之和为 $4$ 和 $6$，因此 $\norm{\mA}_1=6$；各行绝对值之和为 $3$ 和 $7$，因此 $\norm{\mA}_\infty=7$；$\mA^\T\mA=\begin{pmatrix}10&10\\10&20\end{pmatrix}$ 的特征值为 $15\pm\sqrt{125}$，因此 $\norm{\mA}_2=\sqrt{15+\sqrt{125}} \approx5.117$；$\norm{\mA}_F=\sqrt{30}\approx5.477$。**条件数**（condition number）$\kappa(\mA)=\norm{\mA}\,\norm{\mA^{-1}}$ 衡量 $\vb$ 中的相对误差在 $\mA\vx=\vb$ 的解中可能被放大的程度；这里 $\kappa_2(\mA)\approx2.618$，是一个良态矩阵。黑塞矩阵的条件数很大时，梯度法会变慢（{prf:ref}`thm-06-05-sd-rate`）。

(sec-01-01-orthogonality)=
## 正交性与投影

若对 $i\ne j$ 有 $\vq_i^\T\vq_j=0$，且 $\norm{\vq_i}=1$，则称向量 $\vq_1,\dots,\vq_k$ **标准正交**（orthonormal）。标准正交向量线性无关（把 $\sum\alpha_i\vq_i=\vzero$ 与 $\vq_j$ 作内积即得 $\alpha_j=0$）。各列标准正交的方阵 $\mV$ 称为**正交矩阵**（orthogonal matrix）：$\mV^\T\mV=\mI$，因此 $\mV^{-1}=\mV^\T$。正交矩阵保持长度和夹角不变：$\norm{\mV\vx}_2=\norm{\vx}_2$，$(\mV\vx)^\T(\mV\vy)=\vx^\T\vy$；旋转和反射都是正交矩阵的例子，谱定理中的特征向量矩阵也是。子空间 $S$ 的**正交补**（orthogonal complement）为 $S^\perp=\{\vy:\vy^\T\vx=0\ \text{对所有 }\vx\in S\}$；用这一语言，{ref}`sec-01-01-four-subspaces`说的就是 $\nullsp(\mA)=\row(\mA)^\perp$ 以及 $\nullsp(\mA^\T)=\col(\mA)^\perp$。

**到直线上的投影。** 在过 $\vzero$、方向为 $\va\ne\vzero$ 的直线上，距 $\vb$ 最近的点 $\vp$ 可以通过要求误差与直线正交求得。令 $\vp=t\va$，条件 $\va^\T(\vb-t\va)=0$ 给出 $t=\va^\T\vb/\va^\T\va$，因此

$$
\begin{aligned}
\vp&=\frac{\va^\T\vb}{\va^\T\va}\,\va=\mP\vb,\\
\mP&=\frac{\va\va^\T}{\va^\T\va}.
\end{aligned}
$$ (eq-01-01-projection-line)

它是最近点，因为对直线上的任意其他点 $s\va$，由勾股定理（右边两个向量正交）有 $\norm{\vb-s\va}^2=\norm{\vb-\vp}^2+\norm{\vp-s\va}^2$。

:::{prf:proposition} 到列空间上的投影
:label: prop-01-01-projection

设 $\mA\in\R^{m\times n}$ 列满秩（秩为 $n$）。则 $\mA^\T\mA$ 可逆，且 $\col(\mA)$ 中距 $\vb\in\R^m$ 最近的点为 $\vp=\mA\hat{\vx}$，其中 $\hat{\vx}$ 是**正规方程**（normal equations）

$$
\mA^\T\mA\,\hat{\vx}=\mA^\T\vb .
$$ (eq-01-01-normal-equations)

的解。投影矩阵 $\mP=\mA(\mA^\T\mA)^{-1}\mA^\T$ 是对称的，满足 $\mP^2=\mP$，且误差 $\vb-\vp$ 位于 $\nullsp(\mA^\T)$ 中。
:::

:::{dropdown} 证明
:class-container: proof
:animate: fade-in

若 $\mA^\T\mA\vx=\vzero$，则 $\vx^\T\mA^\T\mA\vx=\norm{\mA\vx}^2=0$，故 $\mA\vx=\vzero$，再由列满秩得 $\vx=\vzero$；因此 $\mA^\T\mA$ 可逆。误差 $\vb-\mA\hat{\vx}$ 必须与 $\mA$ 的每一列正交，即 $\mA^\T(\vb-\mA\hat{\vx})=\vzero$，这就是 {eq}`eq-01-01-normal-equations`。对任意其他的 $\mA\vx$，由勾股定理得 $\norm{\vb-\mA\vx}^2=\norm{\vb-\mA\hat{\vx}}^2+\norm{\mA(\hat{\vx}-\vx)}^2$，因此 $\mA\hat{\vx}$ 最近。$\mP^2=\mA(\mA^\T\mA)^{-1}\mA^\T\mA (\mA^\T\mA)^{-1}\mA^\T=\mP$ 和对称性可直接验证。
:::

{doc}`01-02-matrix-calculus`通过令 $\norm{\mA\vx-\vb}_2^2$ 的梯度为零推导出同样的正规方程；这就是回归中的最小二乘法（{ref}`sec-06-01-regression`）。矩阵 $\mA^\T\mA$ 总是半正定的，因为 $\vx^\T\mA^\T\mA\vx=\norm{\mA\vx}^2\ge0$；在列满秩时它是正定的。

:::{prf:example} 把向量投影到直线上
:label: ex-01-01-projection-line

把 $\vb=(1,2)$ 投影到过 $\va=(2,1)$ 的直线上：求投影 $\vp$、误差 $\vr=\vb-\vp$、$\vb$ 到直线的距离以及投影矩阵 $\mP$。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

这里 $\va^\T\vb=4$，$\va^\T\va=5$，因此 $\vp=\tfrac45\va=(1.6,0.8)$。误差 $\vr=\vb-\vp=(-0.6,1.2)$ 与直线正交：$\va^\T\vr=-1.2+1.2=0$。$\vb$ 到直线的距离为 $\norm{\vr}=\sqrt{0.36+1.44}=\sqrt{1.8}\approx1.342$。投影矩阵为

$$
\mP=\frac{1}{5}\begin{pmatrix}4&2\\2&1\end{pmatrix}
=\begin{pmatrix}0.8&0.4\\0.4&0.2\end{pmatrix},
$$

确实 $\mP\vb=(0.8+0.8,\ 0.4+0.4)=(1.6,0.8)$，且 $\mP^2=\mP$（投影两次与投影一次相同）。{numref}`fig-01-01-projection-line`画出了 $\vp$ 处的直角。

```{figure} figures/01-01-projection-line.svg
:name: fig-01-01-projection-line
:width: 460px
:class: dark-light
:alt: 从原点出发的向量 a 等于 (2, 1) 和 b 等于 (1, 2)，b 在过 a 的直线上的投影 p，以及从 p 指向 b 的虚线误差向量 r，它与直线成直角

把 $\vb=(1,2)$ 正交投影到 $\va=(2,1)$ 张成的直线上：$\vp=0.8\,\va=(1.6,0.8)$，误差 $\vr=\vb-\vp$ 垂直于该直线（直角标记），因此 $\vp$ 是直线上距 $\vb$ 最近的点。
```
:::

对一列向量反复做到直线上的投影，并减去已经求得的分量，就能把任意一组基变为标准正交基（**格拉姆–施密特正交化**，Gram–Schmidt process）；写成矩阵形式，就是数值稳定的最小二乘求解器所用的分解 $\mA=\mQ\mathbf R$（Strang, 2016, Section 4.4）。

## 本章回顾

### 小结

- 乘积 $\mA\vx$ 是 $\mA$ 各列的组合；行变换就是左乘（{eq}`eq-01-01-matrix-product`、{eq}`eq-01-01-elementary-matrix`）。
- 内积给出长度、距离和夹角；它的符号区分锐角、直角和钝角（{eq}`eq-01-01-angle`、{prf:ref}`thm-01-01-cauchy-schwarz`）。
- 若 $\mX\valpha=\vzero$ 迫使 $\valpha=\vzero$，则这些向量线性无关，可用消元法检验（{prf:ref}`prop-01-01-independence-test`）；秩等于主元个数，也等于行秩和列秩（{prf:ref}`thm-01-01-rank`）。
- 以 $(r,s)$ 元素为主元的主元变换就是左乘一个 eta 矩阵（{prf:ref}`prop-01-01-pivot-matrix`）；替换基矩阵的一列需要非零主元，并按 $\mB'^{-1}=\mE\mB^{-1}$ 更新逆矩阵，由此得到乘积形式 {eq}`eq-01-01-product-form`（{prf:ref}`prop-01-01-basis-exchange`、{prf:ref}`ex-01-01-pivot-product-form`）。
- $\mA\vx=\vb$ 有解当且仅当 $\vb\in\col(\mA)$；全部解构成 $\vx_p+\nullsp(\mA)$；系数矩阵奇异的方程组要么无解，要么有无穷多个解（{prf:ref}`thm-01-01-existence-uniqueness`、{prf:ref}`thm-01-01-invertible`）。
- 分块消元给出线性规划的基分块 {eq}`eq-01-01-basis-partition`、分块矩阵的逆 {eq}`eq-01-01-block-inverse` 和加边逆矩阵 {eq}`eq-01-01-bordered-inverse`。
- 行列式等于各主元之积，它为零当且仅当矩阵奇异；代数余子式给出伴随矩阵公式 {eq}`eq-01-01-adjugate`。
- 对称矩阵的特征值都是实数，且有标准正交特征基 $\mA=\mV\mD\mV^\T$（{prf:ref}`thm-01-01-spectral`），其瑞利商介于最小和最大特征值之间（{prf:ref}`prop-01-01-rayleigh`）。
- 有定性（{prf:ref}`def-01-01-definiteness`）可以用特征值判别（{prf:ref}`thm-01-01-eigen-test`），正定性可以用顺序主子式判别（{prf:ref}`thm-01-01-sylvester`），半正定性则要用全部主子式判别（{prf:ref}`thm-01-01-psd-minors`）；{numref}`fig-01-01-quadratic-forms`给出了三种形状。
- 1-范数、2-范数、$\infty$-范数及其诱导矩阵范数 {eq}`eq-01-01-matrix-norms`；投影由正规方程 {eq}`eq-01-01-normal-equations` 求得（{prf:ref}`prop-01-01-projection`）。

### 常见错误

- **对非对称矩阵使用子式判别法。** 西尔维斯特判据和特征值判别法都是关于对称矩阵的结论。对非对称矩阵 $\mA$，应改为检验其对称部分 $\tfrac12(\mA+\mA^\T)$（{eq}`eq-01-01-symmetric-part`）：$\begin{pmatrix}1&-3\\0&1\end{pmatrix}$ 的顺序主子式都为正，但作为二次型它是不定的。
- **用顺序主子式判别半正定性。** “所有顺序主子式 $\ge0$”并不蕴含 $\mA\succeq0$；$\diag(0,-1)$ 就是反例。半正定性要求**所有**主子式 $\ge0$（{prf:ref}`thm-01-01-psd-minors`），或所有特征值 $\ge0$。
- **弄错负定性的判别条件。** 负定矩阵的顺序主子式并不全为负，而是正负交替：$D_1<0$，$D_2>0$，$D_3<0,\dots$
- **以零元素为主元。** 主元必须非零；零主元会使新的基矩阵奇异（{prf:ref}`prop-01-01-basis-exchange`）。在单纯形法中，主元还必须为正，由比值检验选出。
- **把逆和转置的乘积次序弄反。** $(\mA\mB)^{-1}=\mB^{-1}\mA^{-1}$，$(\mA\mB)^\T=\mB^\T\mA^\T$，而不是反过来；{eq}`eq-01-01-product-form` 中的 $\mE_2\mE_1$ 把最近一次主元变换写在最前面。
- **由方程组是方的就断定解唯一。** 系数矩阵奇异的方程组要么无解，要么有无穷多个解，绝不会恰有一个解；务必先检验 $\det\mA\ne0$ 或主元的个数。
- **把 $\vx^\T\vy$ 当作矩阵。** $\vx^\T\vy$ 是一个数，而 $\vx\vy^\T$ 是一个秩为一的 $n\times n$ 矩阵。

## 参考文献

- Bazaraa, M. S., Jarvis, J. J., & Sherali, H. D. (2010). *Linear programming
  and network flows* (4th ed.). Wiley. (Sections 2.1–2.3: vectors, matrices
  and simultaneous linear equations.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Appendix
  A.1, norms; Appendix A.5, linear algebra: range and null space, symmetric
  eigenvalue decomposition, Schur complement.)
- Magnus, J. R., & Neudecker, H. (2019). *Matrix differential calculus with
  applications in statistics and econometrics* (3rd ed.). Wiley.
  <https://doi.org/10.1002/9781119541219> (Chapter 1: basic properties of
  vectors and matrices, including eigenvalues, the spectral theorem,
  definite matrices and partitioned matrices.)
- Petersen, K. B., & Pedersen, M. S. (2012). *The matrix cookbook* (Version
  November 15, 2012). Technical University of Denmark. (Sections 1, 3.1–3.2
  and 9.1: determinants, inverses, block matrices, eigenvalues, norms.)
- Strang, G. (2016). *Introduction to linear algebra* (5th ed.).
  Wellesley-Cambridge Press. (Chapters 1–2: vectors and elimination;
  Chapter 3: the four subspaces and rank; Chapter 4: orthogonality and
  projections; Chapter 5: determinants; Chapter 6: eigenvalues, symmetric
  and positive definite matrices.)
