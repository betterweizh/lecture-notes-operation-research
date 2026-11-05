# 对偶

本章讲解{doc}`03-05-lp-duality`的例题与习题。

题目按主题分组：写出对偶问题、弱对偶性与强对偶性、互补松弛性、检验一个基的原始可行性与对偶可行性、从单纯形表和 $\mB^{-1}$ 中读出对偶值，以及影子价格（shadow price）。记号与{doc}`03-05-lp-duality`相同；所有数字都用精确算术重新计算过，并用 HiGHS 求解器验证。

(sec-08-03-writing-duals)=
## 写出对偶问题

:::{prf:example} 三个对偶问题
:label: ex-08-03-a4-q1

写出下列每个线性规划问题（linear program）的对偶问题（dual problem），并求解每一对中的两个问题。

**(a)** $\max\ z=66x_1-22x_2$，满足 $-x_1+x_2\le-2$，$2x_1+3x_2\le5$，$x_1,x_2\ge0$。

**(b)** $\min\ z=6x_1+3x_2$，满足 $6x_1-3x_2+x_3\ge25$，$3x_1+4x_2+x_3\ge55$，$x_1,x_2,x_3\ge0$。

**(c)** $\max\ z=x_1+x_2$，满足 $2x_1+x_2=5$，$3x_1-x_2=6$，$x_1,x_2$ 为自由变量。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

我们应用{numref}`tab-03-05-dual-rules`；对最大化的原问题，从右往左读这张表。

**(a)** 两个约束都是最大化问题中的 $\le$ 行（自然类型（natural type）），因此两个对偶变量（dual variable）都非负；两个原始变量都非负，因此两个对偶约束都是 $\ge$ 约束：

$$
\begin{aligned}
\min_{\vpi}\ & w=-2\pi_1+5\pi_2\\
\st\ & -\pi_1+2\pi_2\ge 66\\
& \pi_1+3\pi_2\ge -22\\
& \pi_1,\pi_2\ge0 .
\end{aligned}
$$

另一种等价形式是先把原问题的第一行乘以 $-1$（得到 $x_1-x_2\ge2$，一个相反类型的约束）；这时它的对偶变量 $p_1=-\pi_1$ 非正，对偶问题为 $\min\{2p_1+5p_2: p_1+2p_2\ge66,\ -p_1+3p_2\ge-22,\ p_1\le0,\ p_2\ge0\}$。

*求解。* 原问题（primal problem）的最优解为 $\vx^*=(2.5,\,0)$，$z^*=165$：第二个约束是紧约束（$2\cdot2.5=5$），第一个约束的松弛量为 $-2-(-2.5)=0.5$。由互补松弛性（complementary slackness）得 $\pi_1^*=0$；又因为 $x_1^*>0$，第一个对偶约束是紧的：$2\pi_2=66$，故 $\vpi^*=(0,\,33)$。它是对偶可行的（$0+99\ge-22$），且 $w^*=5\cdot33=165=z^*$，这同时证明了两者的最优性。

**(b)** 约束是最小化问题中的 $\ge$ 行（自然类型），变量都非负；$x_3$ 的费用为 0：

$$
\begin{aligned}
\max_{\vpi}\ & w=25\pi_1+55\pi_2\\
\st\ & 6\pi_1+3\pi_2\le 6\\
& -3\pi_1+4\pi_2\le 3\\
& \pi_1+\pi_2\le 0\\
& \pi_1,\pi_2\ge0 .
\end{aligned}
$$

第三个约束连同 $\vpi\ge\vzero$ 迫使 $\vpi^*=(0,0)$，因此 $w^*=0$。在原问题中，$x_3$ 不花费任何费用，并以系数 1 出现在两行中，因此 $\vx=(0,0,55)$ 是费用为 0 的可行解；而由于所有费用都非负，没有哪个可行点的费用更低：$z^*=0=w^*$。（每个满足 $x_3\ge55$ 的 $\vx=(0,0,x_3)$ 都是最优解。）

**(c)** 两个约束都是等式，因此两个对偶变量都是自由变量；两个原始变量都是自由变量，因此两个对偶约束都是等式：

$$
\begin{aligned}
\min_{\vpi}\ & w=5\pi_1+6\pi_2\\
\st\ & 2\pi_1+3\pi_2=1\\
& \pi_1-\pi_2=1\\
& \pi_1,\pi_2\ \text{自由}.
\end{aligned}
$$

每个问题都恰有一个可行点。原问题：$5x_1=11$，故 $\vx=(2.2,\,0.6)$，$z=2.8$。对偶问题：由 $\pi_1=1+\pi_2$ 得 $5\pi_2=-1$，故 $\vpi=(0.8,\,-0.2)$，$w=4-1.2=2.8$。两个值相等，正如强对偶性（strong duality）所要求的。
:::

:::{prf:example} 分块形式问题的对偶
:label: ex-08-03-block-dual

把一个最小化问题的变量分成非负块 $\vx_1\in\R^{n_1}$ 和自由块 $\vx_2\in\R^{n_2}$，把约束分成 $m_1$ 个等式和 $m_2$ 个 $\ge$ 型不等式：

$$
\begin{aligned}
\min_{\vx_1,\vx_2}\ & \vc_1^\T\vx_1+\vc_2^\T\vx_2\\
\st\ & \mA_{11}\vx_1+\mA_{12}\vx_2=\vb_1\\
& \mA_{21}\vx_1+\mA_{22}\vx_2\ge\vb_2\\
& \vx_1\ge\vzero,\quad \vx_2\ \text{自由},
\end{aligned}
$$

其中 $\mA_{kl}\in\R^{m_k\times n_l}$。写出它的对偶问题，并直接对这一对问题证明弱对偶性（weak duality）。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**对偶问题。** 给等式行配以自由的对偶块 $\vpi_1\in\R^{m_1}$，给 $\ge$ 行配以非负的对偶块 $\vpi_2\in\R^{m_2}$。$\vx_1$ 的列为 $\begin{pmatrix}\mA_{11}\\\mA_{21}\end{pmatrix}$，给出 $\le$ 约束（$\vx_1\ge\vzero$）；$\vx_2$ 的列给出等式（$\vx_2$ 自由）：

$$
\begin{aligned}
\max_{\vpi_1,\vpi_2}\ & \vb_1^\T\vpi_1+\vb_2^\T\vpi_2\\
\st\ & \mA_{11}^\T\vpi_1+\mA_{21}^\T\vpi_2\le\vc_1\\
& \mA_{12}^\T\vpi_1+\mA_{22}^\T\vpi_2=\vc_2\\
& \vpi_1\ \text{自由},\quad \vpi_2\ge\vzero .
\end{aligned}
$$

第一个对偶约束中出现的是 $\mA_{21}$，即在不等式行中与 $\vx_1$ 相乘的块。*（这对问题的早期版本在该约束中把 $\mA_{21}$ 写成了 $\mA_{12}$；$\vx_1$ 的列块的转置是 $(\mA_{11}^\T\ \ \mA_{21}^\T)$。）*

**弱对偶性。** 设 $(\vx_1,\vx_2)$ 与 $(\vpi_1,\vpi_2)$ 都可行。则

$$
\begin{aligned}
\vc_1^\T\vx_1+\vc_2^\T\vx_2
&\ge(\mA_{11}^\T\vpi_1+\mA_{21}^\T\vpi_2)^\T\vx_1\\
&\quad+(\mA_{12}^\T\vpi_1+\mA_{22}^\T\vpi_2)^\T\vx_2\\
&=\vpi_1^\T(\mA_{11}\vx_1+\mA_{12}\vx_2)\\
&\quad+\vpi_2^\T(\mA_{21}\vx_1+\mA_{22}\vx_2)\\
&\ge\vpi_1^\T\vb_1+\vpi_2^\T\vb_2 .
\end{aligned}
$$

第一个不等式用到了 $\vx_1\ge\vzero$ 与 $\le$ 型对偶约束（$\vx_2$ 部分是等式，与 $\vx_2$ 的符号无关）。在最后一行中，等式块恰好给出 $\vpi_1^\T\vb_1$，与 $\vpi_1$ 的符号无关；而由于 $\vpi_2\ge\vzero$，不等式块至少给出 $\vpi_2^\T\vb_2$。四条符号约定恰好是每一步所需要的。当 $n_2=m_1=0$ 时，这对问题就是对称形式的问题对 {eq}`eq-03-05-canonical-pair`；当 $n_2=m_2=0$ 时，它就是标准形问题对 {eq}`eq-03-05-standard-pair`。
:::

:::{prf:example} 人工变量对应的对偶约束
:label: ex-08-03-a4-q2

考虑

$$
\begin{aligned}
\max_{\vx}\ & z=5x_1+12x_2+4x_3\\
\st\ & x_1+2x_2+x_3\le 10\\
& 2x_1-x_2+3x_3=8\\
& x_1,x_2,x_3\ge0 .
\end{aligned}
$$

为了启动大 M 法（big-M method），在第一行加入松弛变量 $x_4$，在第二行加入目标函数系数为 $-M$ 的人工变量 $R$（{ref}`sec-03-04-big-m`）。证明人工变量不改变对偶问题：一旦 $M$ 足够大，它所产生的对偶约束就是冗余的。这里 $M$ 必须多大？
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**原问题的对偶。** 第一行是最大化问题中的 $\le$ 行（$\pi_1\ge0$），第二行是等式（$\pi_2$ 自由）：

$$
\begin{aligned}
\min_{\vpi}\ & w=10\pi_1+8\pi_2\\
\st\ & \pi_1+2\pi_2\ge5\\
& 2\pi_1-\pi_2\ge12\\
& \pi_1+3\pi_2\ge4\\
& \pi_1\ge0,\quad \pi_2\ \text{自由}.
\end{aligned}
$$

**大 M 问题的对偶。** 写成等式形式为 $\max\{5x_1+12x_2+4x_3+0x_4-MR:\ x_1+2x_2+x_3+x_4=10,\ 2x_1-x_2+3x_3+R=8,\ \vx,R\ge0\}$。两行都是等式，因此两个对偶变量都是自由变量，而每个变量都给出一个 $\ge$ 约束。$x_1,x_2,x_3$ 的列给出上面的三个约束；费用为 0 的 $x_4$ 的列 $(1,0)^\T$ 给出 $\pi_1\ge0$，这恢复了原对偶问题中的符号限制；费用为 $-M$ 的 $R$ 的列 $(0,1)^\T$ 给出

$$
\pi_2\ge-M .
$$

因此两个对偶问题只相差自由变量 $\pi_2$ 的这个下界。当 $M\to\infty$ 时，它不排除任何点，这就是它冗余的含义。

**多大才算足够大。** 原问题的最优解为 $\vx^*=(26/5,\,12/5,\,0)$，$z^*=274/5=54.8$；对偶最优解为 $\vpi^*=(29/5,\,-2/5)$：它使前两个对偶约束取等号（$5.8-0.8=5$，$11.6+0.4=12$；$x_1^*,x_2^*>0$），第三个约束有松弛（$5.8-1.2=4.6\ge4$；$x_3^*=0$），且 $w=58-3.2=54.8=z^*$。当 $-2/5\ge-M$，即 $M\ge2/5$ 时，附加约束不会切掉这个点；对每个这样的 $M$，大 M 问题都有相同的最优解，且 $R=0$。在大 M 单纯形表中，$R$ 的列给出 $\bar c_R=-M-\pi_2=-M+2/5$（{numref}`tab-03-05-reading-duals`），当 $M=100$ 时它等于 $-498/5$，并且恰好在 $M\ge2/5$ 时非正（对最大化问题即为最优）。当 $M<2/5$ 时罚项太弱：例如取 $M=0.1$，大 M 问题的最优解会让 $R=13$ 留在基中（$\vx=(0,5,0)$，$z-MR=58.7$），而相应的对偶最优解落在界 $\pi_2=-0.1$ 上。
:::

:::{prf:example} 变量有界的问题，其对偶问题总是可行的
:label: ex-08-03-a4-q10

证明问题

$$
\begin{aligned}
\max_{\vx}\ & \vc^\T\vx\\
\st\ & \mA\vx\le\vb\\
& \vl\le\vx\le\vu ,
\end{aligned}
$$

（界向量满足 $\vzero<\vl\le\vu$）的对偶问题总有可行解。这对原问题说明了什么？
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

把界看作约束，并把 $\vx$ 看作自由变量（界本身已经蕴含 $\vx\ge\vl>\vzero$）。约束行为 $\mA\vx\le\vb$（价格 $\vp\ge\vzero$）、$\mI\vx\le\vu$（价格 $\vq\ge\vzero$）和 $-\mI\vx\le-\vl$（价格 $\vr\ge\vzero$），都是最大化问题中的 $\le$ 行。由于 $\vx$ 自由，每个对偶约束都是等式：

$$
\begin{aligned}
\min_{\vp,\vq,\vr}\ & w=\vb^\T\vp+\vu^\T\vq-\vl^\T\vr\\
\st\ & \mA^\T\vp+\vq-\vr=\vc\\
& \vp,\vq,\vr\ge\vzero .
\end{aligned}
$$

一个可行点是 $\vp=\vzero$，$q_j=\max(c_j,0)$，$r_j=\max(-c_j,0)$：此时对每个 $j$ 都有 $q_j-r_j=c_j$，且三个向量都非负。（若把 $\vx\ge\vzero$ 保留为符号限制，对偶约束变为 $\mA^\T\vp+\vq-\vr\ge\vc$，同一个点仍然可行。）

*（对偶问题是最小化问题；本题解答的早期版本写的是“最大化”。）*

对原问题的推论可由{numref}`tab-03-05-status`得出：由于对偶问题可行，原问题不可能无界——这一点也很明显，因为 $\vx$ 位于一个长方体中。要么原问题不可行（约束行 $\mA\vx\le\vb$ 可能排除整个长方体），此时对偶问题无界；要么原问题可行，此时两个问题都有最优解且最优值相等。这个对偶可行点还给出界 $z^*\le\sum_j\bigl(u_j\max(c_j,0)-\ell_j\max(-c_j,0)\bigr)$：即每个变量都取其最有利的界。
:::

(sec-08-03-weak-strong)=
## 弱对偶性与强对偶性

:::{prf:example} 原问题与对偶问题都不可行
:label: ex-08-03-both-infeasible

**(a)** 写出 $\max\{x_1+3x_2:\ x_1+2x_2=1,\ x_1+2x_2=2,\ x_1,x_2\ \text{自由}\}$ 的对偶问题，并证明两个问题都没有可行解。

**(b)** 若加上符号限制 $x_1,x_2\ge0$，会有什么变化？
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** 两个等式给出自由的对偶变量；两个自由的原始变量给出对偶等式：

$$
\begin{aligned}
\min_{\vpi}\ & w=\pi_1+2\pi_2\\
\st\ & \pi_1+\pi_2=1\\
& 2\pi_1+2\pi_2=3\\
& \vpi\ \text{自由}.
\end{aligned}
$$

原问题要求 $x_1+2x_2$ 既等于 1 又等于 2，这不可能。对偶问题要求 $\pi_1+\pi_2$ 既等于 1 又等于 $3/2$，这同样不可能。这就是{numref}`tab-03-05-status`中“两者都不可行”的情形，它说明原问题不可行并不意味着对偶问题无界。

**(b)** 原问题仍然不可行。对偶约束变为 $\pi_1+\pi_2\ge1$ 和 $2\pi_1+2\pi_2\ge3$，例如每个 $\vpi=(t,\ 3/2-t)$ 都满足它们；当 $t\to\infty$ 时，沿这些点有 $w=t+3-2t=3-t\to-\infty$。此时对偶问题无界，属于表中“原问题不可行、对偶问题无界”的情形。
:::

:::{prf:example} 分别由原问题与对偶问题求最优值
:label: ex-08-03-a4-q9

某个最大化问题有两个变量 $x_1,x_2$ 和三个 $\le$ 约束，松弛变量为 $x_3,x_4,x_5$。它的最优基由 $x_1$、$x_2$ 和 $x_3$ 的列组成，且

$$
\begin{gathered}
\mB^{-1}=\begin{pmatrix}0&-1&1\\0&1&0\\1&1&-1\end{pmatrix},\\
\vx_B=(x_1,x_2,x_3)^\T=(1,3,1)^\T,\\
\vpi^*=(0,3,2)^\T .
\end{gathered}
$$

用两种方法计算最优目标函数值：一次通过原问题，一次通过对偶问题。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

对 $\mB^{-1}$ 求逆，得到基矩阵

$$
\mB=\begin{pmatrix}1&0&1\\0&1&0\\1&1&0\end{pmatrix},
$$

它的第三列 $(1,0,0)^\T$ 确实是第 1 行松弛变量的列，前两列则是 $x_1$ 和 $x_2$ 的约束列。

**对偶方法。** 右端项（right-hand side, RHS）为 $\vb=\mB\vx_B=(1+1,\ 3,\ 1+3)^\T=(2,3,4)^\T$，因此 $w^*=\vb^\T\vpi^*=0\cdot2+3\cdot3+2\cdot4=17$。

**原始方法。** 由 $\vpi^{*\T}=\vc_B^\T\mB^{-1}$ 得 $\vc_B^\T=\vpi^{*\T}\mB=(0,3,2)\mB=(2,\ 5,\ 0)$，因此 $c_1=2$，$c_2=5$，而松弛变量的费用为 0，这正是应有的结果。于是 $z^*=\vc_B^\T\vx_B=2\cdot1+5\cdot3+0\cdot1=17$。

两种方法都得到 17，正如强对偶性所要求的。这些数据描述的是问题 $\max\{2x_1+5x_2:\ x_1\le2,\ x_2\le3,\ x_1+x_2\le4,\ \vx\ge\vzero\}$，其最优解 $(1,3)$ 在第一行留下 1 个单位的松弛（因此 $\pi_1^*=0$），并使另外两行成为紧约束。
:::

(sec-08-03-complementary-slackness)=
## 互补松弛性

:::{prf:example} 含三种约束的原问题
:label: ex-08-03-mixed-dual

考虑

$$
\begin{aligned}
\min_{\vx}\ & z=5x_1-6x_2+7x_3+x_4\\
\st\ & x_1+2x_2-x_3-x_4=-7\\
& 6x_1-3x_2+x_3-7x_4\ge14\\
& -2.8x_1-17x_2+4x_3+2x_4\le-3\\
& x_1,x_2\ge0,\quad x_3,x_4\ \text{自由}.
\end{aligned}
$$

**(a)** 写出对偶问题。

**(b)** 有人声称点 $\vx^*=\bigl(0,\ \tfrac{131}{35},\ \tfrac{554}{35},\ -\tfrac{47}{35}\bigr)$ 是最优解。利用互补松弛性求出一个对偶解，并验证这一说法。

**(c)** 证明：若第二个约束中 $x_4$ 的系数是 $+7$ 而不是 $-7$，则原问题无界，并由此得出对偶问题的情况。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** 第一行是等式（$\pi_1$ 自由），第二行是自然类型的 $\ge$ 行（$\pi_2\ge0$），第三行是最小化问题中的 $\le$ 行（$\pi_3\le0$）。变量 $x_1,x_2\ge0$ 给出 $\le$ 约束，自由变量 $x_3,x_4$ 给出等式：

$$
\begin{aligned}
\max_{\vpi}\ & w=-7\pi_1+14\pi_2-3\pi_3\\
\st\ & \pi_1+6\pi_2-2.8\pi_3\le5 &&(x_1)\\
& 2\pi_1-3\pi_2-17\pi_3\le-6 &&(x_2)\\
& -\pi_1+\pi_2+4\pi_3=7 &&(x_3)\\
& -\pi_1-7\pi_2+2\pi_3=1 &&(x_4)\\
& \pi_1\ \text{自由}\\
& \pi_2\ge0\\
& \pi_3\le0 .
\end{aligned}
$$

若先把原问题的第三行乘以 $-1$（使之成为自然类型的 $\ge$ 行），得到的是同一个对偶问题，只是 $\pi_3'=-\pi_3\ge0$，且系数 $-2.8$、$-17$、$+4$、$+2$、$-3$ 换成 $+2.8$、$+17$、$-4$、$-2$、$+3$。

**(b)** 首先，$\vx^*$ 是可行的：$x_1^*=0$，$x_2^*>0$，三行的值为

$$
\begin{gathered}
\tfrac{262-554+47}{35}=-7,\\
\tfrac{-393+554+329}{35}=14,\\
\tfrac{-2227+2216-94}{35}=-3,
\end{gathered}
$$

因此三个约束都是紧约束；其费用为 $z=\tfrac{-786+3878-47}{35}=\tfrac{3045}{35}=87$。互补松弛性要求 $x_2$ 的对偶约束是紧的（$x_2^*>0$）；$x_3$ 和 $x_4$ 的约束本来就是等式，而原问题的各行都是紧的，因此对 $\pi_2,\pi_3$ 没有附加条件。解方程组

$$
\begin{aligned}
2\pi_1-3\pi_2-17\pi_3&=-6\\
-\pi_1+\pi_2+4\pi_3&=7\\
-\pi_1-7\pi_2+2\pi_3&=1
\end{aligned}
$$

得 $\vpi=(-10,\ 1,\ -1)$。符号正确（$\pi_2=1\ge0$，$\pi_3=-1\le0$），剩下的对偶约束也成立：$-10+6+2.8=-1.2\le5$（有松弛，与 $x_1^*=0$ 相符）。因此 $\vpi$ 是对偶可行的，且与 $\vx^*$ 互补；其值为 $w=70+14+3=87=z$。故 $\vx^*$ 与 $\vpi^*=(-10,1,-1)$ 都是最优解，$z^*=w^*=87$。

**(c)** 当第二行为 $+7x_4$ 时，点 $(0,0,-10,17)$ 是可行的（各行的值为 $10-17=-7$，$-10+119=109\ge14$，$-40+34=-6\le-3$），费用为 $-53$。沿方向 $\vd=(0,0,-1,1)$，各行的值分别改变 $1-1=0$、$-1+7=6\ge0$ 和 $-4+2=-2\le0$，因此每个点 $(0,0,-10-t,\ 17+t)$，$t\ge0$ 都保持可行，而 $t$ 每增加一个单位费用改变 $\vc^\T\vd=-7+1=-6$。费用趋于 $-\infty$：原问题无界，由{prf:ref}`cor-03-05-optimality-certificate`的第 3 部分，对偶问题（此时最后一个约束为 $-\pi_1+7\pi_2+2\pi_3=1$）不可行。

*（本题的早期版本在模型中写的是 $+7$，而在系数表中写的是 $-7$。只有 $-7$ 的版本有最优解；取 $+7$ 时，这对问题落入{numref}`tab-03-05-status`中“无界/不可行”的格子。）*
:::

(sec-08-03-feasibility-tests)=
## 对偶可行性与最优性检验

下面两个例子使用同一个标准形问题 $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$，其中

$$
\begin{gathered}
\mA=\begin{pmatrix}1&0&0&-1&2&4\\0&1&0&1&-2&6\\0&0&1&2&3&-7\end{pmatrix},\\
\vb=\begin{pmatrix}7\\19\\21\end{pmatrix},\\
\vc=(3,\ -4,\ -8,\ 16,\ -17,\ 19)^\T .
\end{gathered}
$$ (eq-08-03-six-variable-lp)

:::{prf:example} 这个基是最优的吗？
:label: ex-08-03-basis-optimal

对问题 {eq}`eq-08-03-six-variable-lp`，取基变量为 $(x_6,x_2,x_3)$ 的基。它是原始可行的吗？是对偶可行的吗？由此能得出什么结论？
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

基矩阵及其逆为

$$
\begin{gathered}
\mB=\begin{pmatrix}4&0&0\\6&1&0\\-7&0&1\end{pmatrix},\\
\mB^{-1}=\begin{pmatrix}\tfrac14&0&0\\-\tfrac32&1&0\\\tfrac74&0&1\end{pmatrix}.
\end{gathered}
$$

**原始可行性。** $\vx_B=\mB^{-1}\vb=\bigl(\tfrac74,\ -\tfrac{21}{2}+19,\ \tfrac{49}{4}+21\bigr)=\bigl(\tfrac74,\ \tfrac{17}{2},\ \tfrac{133}{4}\bigr)\ge\vzero$，因此该基是原始可行（primal feasible）的：$x_6=\tfrac74$，$x_2=\tfrac{17}{2}$，$x_3=\tfrac{133}{4}$，其余变量为零。

**对偶可行性。** 取 $\vc_B=(19,-4,-8)$，

$$
\begin{aligned}
\vpi^\T=\vc_B^\T\mB^{-1}&=\bigl(\tfrac{19}{4}+6-14,\ -4,\ -8\bigr)\\
&=\bigl(-\tfrac{13}{4},\ -4,\ -8\bigr).
\end{aligned}
$$

基变量的对偶约束 $\vpi^\T\mA_{\cdot j}\le c_j$ 按构造都取等号；对非基变量，

$$
\begin{aligned}
x_1:&\quad \pi_1=-\tfrac{13}{4}\le3,\\
x_4:&\quad -\pi_1+\pi_2+2\pi_3=-\tfrac{67}{4}\le16,\\
x_5:&\quad 2\pi_1-2\pi_2+3\pi_3=-\tfrac{45}{2}\le-17.
\end{aligned}
$$

检验数（reduced cost）就是这些差：$\bar c_1=3+\tfrac{13}{4}=\tfrac{25}{4}$，$\bar c_4=16+\tfrac{67}{4}=\tfrac{131}{4}$，$\bar c_5=-17+\tfrac{45}{2}=\tfrac{11}{2}$。所有检验数 $\bar\vc=\bigl(\tfrac{25}{4},0,0,\tfrac{131}{4},\tfrac{11}{2},0\bigr)$ 都非负，因此该基是对偶可行（dual feasible）的。

**结论。** 由{prf:ref}`prop-03-05-basis-optimal`，该基是最优的：$\vx^*=\bigl(0,\tfrac{17}{2},\tfrac{133}{4},0,0,\tfrac74\bigr)$ 是原问题的最优解，$\vpi^*=\bigl(-\tfrac{13}{4},-4,-8\bigr)$ 是对偶问题的最优解，且

$$
\begin{aligned}
z^*&=19\cdot\tfrac74-4\cdot\tfrac{17}{2}-8\cdot\tfrac{133}{4}\\
&=\tfrac{133-136-1064}{4}=-\tfrac{1067}{4}=-266.75\\
&=\vb^\T\vpi^*=-\tfrac{91}{4}-76-168 .
\end{aligned}
$$
:::

:::{prf:example} 检验数就是对偶松弛变量
:label: ex-08-03-dual-slacks

对问题 {eq}`eq-08-03-six-variable-lp`，取基 $(x_1,x_2,x_3)$。计算它的对偶向量以及全部六个对偶约束的松弛量。哪些对偶约束被违反？单纯形法下一步会做什么？
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

这里 $\mB=\mI$，因此 $\vx_B=\vb=(7,19,21)\ge\vzero$（原始可行），且 $\vpi^\T=\vc_B^\T=(3,-4,-8)$。把对偶约束写成“松弛量 $\ge0$”的形式，它们及其取值为：

```{list-table} $\vpi=(3,-4,-8)$ 处的对偶松弛变量
:header-rows: 1

* - 变量
  - 对偶约束
  - 松弛量 $\bar c_j$
* - $x_1$
  - $3-\pi_1\ge0$
  - $0$
* - $x_2$
  - $-4-\pi_2\ge0$
  - $0$
* - $x_3$
  - $-8-\pi_3\ge0$
  - $0$
* - $x_4$
  - $16+\pi_1-\pi_2-2\pi_3\ge0$
  - $39$
* - $x_5$
  - $-17-2\pi_1+2\pi_2-3\pi_3\ge0$
  - $-7$
* - $x_6$
  - $19-4\pi_1-6\pi_2+7\pi_3\ge0$
  - $-25$
```

这些松弛量恰好就是该基的单纯形表中的检验数，该表的目标行为 $(0,0,0,39,-7,-25\mid 223)$，且 $z=\vc_B^\T\vb=21-76-168=-223$。$x_5$ 和 $x_6$ 的对偶约束被违反，因此该基不是对偶可行的；又因为它是非退化的，所以它不是最优的：增大 $x_5$ 或 $x_6$ 都会降低费用。

单纯形法让 $x_6$ 进基（检验数最负）。它的列为 $(4,6,-7)$，比值检验给出 $\min\{7/4,\ 19/6\}=7/4$，在 $x_1$ 所在行取得，因此 $x_1$ 出基。新的基 $(x_6,x_2,x_3)$ 就是{prf:ref}`ex-08-03-basis-optimal`中的最优基，费用下降 $25\cdot\tfrac74=43.75$，从 $-223$ 降到 $-266.75$。
:::

:::{prf:example} 用对偶问题否定最优性
:label: ex-08-03-a4-q7

利用对偶问题证明：对于

$$
\begin{aligned}
\max_{\vx}\ & z=x_1+2x_2+4x_3-3x_4\\
\st\ & x_1+x_2+x_3=4\\
& x_1+4x_2+x_4=8\\
& \vx\ge\vzero ,
\end{aligned}
$$

基变量为 $(x_1,x_2)$ 的基本解不是最优的，并求出两个问题的最优解。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**对偶问题。** 两行都是等式（$\vpi$ 自由），变量都非负，因此对这个最大化问题，其对偶问题为

$$
\begin{aligned}
\min_{\vpi}\ & w=4\pi_1+8\pi_2\\
\st\ & \pi_1+\pi_2\ge1 &&(x_1)\\
& \pi_1+4\pi_2\ge2 &&(x_2)\\
& \pi_1\ge4 &&(x_3)\\
& \pi_2\ge-3 &&(x_4)\\
& \pi_1,\pi_2\ \text{自由}.
\end{aligned}
$$

**基 $(x_1,x_2)$。** 由 $\mB=\begin{pmatrix}1&1\\1&4\end{pmatrix}$，$\mB^{-1}=\tfrac13\begin{pmatrix}4&-1\\-1&1\end{pmatrix}$，基本解为 $\vx_B=\mB^{-1}\vb=\bigl(\tfrac83,\tfrac43\bigr)>\vzero$，$z=\tfrac83+\tfrac83=\tfrac{16}{3}$，单纯形乘子为 $\vpi^\T=\vc_B^\T\mB^{-1}=(1,2)\cdot\tfrac13\begin{pmatrix}4&-1\\-1&1\end{pmatrix}=\bigl(\tfrac23,\tfrac13\bigr)$。如果该基是最优的，那么这个 $\vpi$ 就是对偶问题的最优解（{prf:ref}`thm-03-05-strong-duality`），因而是可行的。但它违反了 $x_3$ 的对偶约束 $\pi_1\ge4$（检验数 $\bar c_3=4-\tfrac23=\tfrac{10}{3}>0$，在最大化问题中这是一个有利可图的变量）。由于该基本解是非退化的，$x_3$ 可以以正的水平进基，因此不仅这个基不是最优的，这个解本身也不是最优的。

**最优解。** 让 $x_3$ 进基（列 $\mB^{-1}\mA_{\cdot3}=(\tfrac43,-\tfrac13)$），$x_1$ 出基（比值 $\tfrac83/\tfrac43=2$）。新的基 $(x_3,x_2)$ 给出 $\vx^*=(0,2,2,0)$，$z^*=4+8=12$，乘子为 $\vpi^*=(4,\ -\tfrac12)$。它们满足所有对偶约束（$3.5\ge1$，$2\ge2$，$4\ge4$，$-0.5\ge-3$），且 $w=16-4=12=z^*$，因此两者都是最优解。
:::

:::{prf:example} 由基的逆矩阵求单纯形表
:label: ex-08-03-a4-q5

对问题

$$
\begin{aligned}
\min_{\vx}\ & z=2x_1+x_2\\
\st\ & 3x_1+x_2-x_3=3\\
& 4x_1+3x_2-x_4=6\\
& x_1+2x_2+x_5=3\\
& x_1,\dots,x_5\ge0,
\end{aligned}
$$

基变量为 $(x_1,x_2,x_5)$ 的基满足

$$
\mB^{-1}=\begin{pmatrix}\tfrac35&-\tfrac15&0\\-\tfrac45&\tfrac35&0\\1&-1&1\end{pmatrix}.
$$

计算这个基的完整单纯形表，并检验它的可行性和最优性。最优解是唯一的吗？
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**单纯形表。** 约束行为 $\mB^{-1}\mA$ 和 $\mB^{-1}\vb$；取 $\vc_B=(2,1,0)$，单纯形乘子为 $\vpi^\T=\vc_B^\T\mB^{-1}=\bigl(\tfrac65-\tfrac45,\ -\tfrac25+\tfrac35,\ 0\bigr)=\bigl(\tfrac25,\tfrac15,0\bigr)$，由它得到目标行 $\vc^\T-\vpi^\T\mA$ 和 $-z=-\vpi^\T\vb$：

```{list-table} 基 $(x_1,x_2,x_5)$ 的单纯形表
:header-rows: 1
:stub-columns: 1

* - 基变量
  - $x_1$
  - $x_2$
  - $x_3$
  - $x_4$
  - $x_5$
  - 右端项
* - $x_1$
  - $1$
  - $0$
  - $-\tfrac35$
  - $\tfrac15$
  - $0$
  - $\tfrac35$
* - $x_2$
  - $0$
  - $1$
  - $\tfrac45$
  - $-\tfrac35$
  - $0$
  - $\tfrac65$
* - $x_5$
  - $0$
  - $0$
  - $-1$
  - $1$
  - $1$
  - $0$
* - $-z$
  - $0$
  - $0$
  - $\tfrac25$
  - $\tfrac15$
  - $0$
  - $-\tfrac{12}{5}$
```

作为验证，$\mB^{-1}$ 乘以基变量的列 $(3,4,1)$、$(1,3,2)$、$(0,0,1)$ 得到单位向量，而剩余变量 $x_3$、$x_4$ 的列恰好是 $\mB^{-1}$ 前两列的 $-1$ 倍。

**可行性与最优性。** 右端项 $(\tfrac35,\tfrac65,0)$ 非负，因此该基是原始可行的；检验数 $(0,0,\tfrac25,\tfrac15,0)$ 非负，因此它是对偶可行的。所以该基是最优的，$\vx^*=(\tfrac35,\tfrac65,0,0,0)$，$z^*=\tfrac65+\tfrac65=\tfrac{12}{5}$，而 $\vpi=(\tfrac25,\tfrac15,0)$ 是一个对偶最优解。该基是**退化**的：基变量 $x_5$ 为零（第三个约束 $x_1+2x_2\le3$ 在 $\vx^*$ 处是紧的，尽管它的对偶值为 0）。

**唯一性。** 原问题的最优解是唯一的，因为每个非基变量的检验数都严格为正。对偶最优解则不唯一。一次退化的主元变换让 $x_3$ 替换 $x_5$ 进基（以元素 $-1$ 为主元，由于右端项为 0，这是允许的），得到同一点 $\vx^*$ 处的基 $(x_1,x_2,x_3)$，其目标行为 $(0,0,0,\tfrac35,\tfrac25\mid-\tfrac{12}{5})$，同样是对偶可行的，乘子为 $\vpi=(0,\tfrac35,-\tfrac25)$。线段

$$
\vpi=\bigl(\tfrac25+t,\ \tfrac15-t,\ t\bigr),\qquad -\tfrac25\le t\le0,
$$

上的每一点都使 $x_1,x_2$ 的两个对偶约束取等号，并满足符号条件（$\pi_1,\pi_2\ge0$；对 $\le$ 行有 $\pi_3\le0$），其值为 $w=\tfrac{12}{5}$，因此它们都是最优的。正是原问题最优解的退化性使得对偶最优解可以不止一个。相应地，第三个约束有两个单侧影子价格：放松它（增大右端项 3）完全不会降低费用，而收紧它会使费用以每单位 $\tfrac25$ 的速率增加（例如当右端项为 2.9 时 $z^*=2.44$）。
:::

(sec-08-03-tableau-duals)=
## 从单纯形表和 $\mB^{-1}$ 读出对偶值

:::{prf:example} 由大 M 单纯形表求对偶解（一）
:label: ex-08-03-a4-q3

问题

$$
\begin{aligned}
\max_{\vx}\ & z=5x_1+2x_2+3x_3\\
\st\ & x_1+5x_2+2x_3=30\\
& x_1-5x_2-6x_3\le40\\
& \vx\ge\vzero
\end{aligned}
$$

用大 M 法求解，取 $M=100$，在第一行引入人工变量 $x_4$，在第二行引入松弛变量 $x_5$，二者构成初始基。最优单纯形表为

```{list-table}
:header-rows: 1
:stub-columns: 1

* - 基变量
  - $x_1$
  - $x_2$
  - $x_3$
  - $x_4$
  - $x_5$
  - 右端项
* - $x_1$
  - $1$
  - $5$
  - $2$
  - $1$
  - $0$
  - $30$
* - $x_5$
  - $0$
  - $-10$
  - $-8$
  - $-1$
  - $1$
  - $10$
* - $-z$
  - $0$
  - $-23$
  - $-7$
  - $-105$
  - $0$
  - $-150$
```

写出对偶问题，并用两种方法求出它的最优解：由单纯形表求，以及用图解法求。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**对偶问题。** 第一行是等式（$\pi_1$ 自由），第二行是自然类型的 $\le$ 行（$\pi_2\ge0$）：

$$
\begin{aligned}
\min_{\vpi}\ & w=30\pi_1+40\pi_2\\
\st\ & \pi_1+\pi_2\ge5 &&(x_1)\\
& 5\pi_1-5\pi_2\ge2 &&(x_2)\\
& 2\pi_1-6\pi_2\ge3 &&(x_3)\\
& \pi_1\ \text{自由},\quad \pi_2\ge0 .
\end{aligned}
$$

**由单纯形表求。** 初始基为 $(x_4,x_5)$，因此最终表中 $x_4$ 和 $x_5$ 的列就是 $\mB^{-1}=\begin{pmatrix}1&0\\-1&1\end{pmatrix}$；对 $(x_1,x_5)$ 取 $\vc_B=(5,0)$，得 $\vpi^{*\T}=\vc_B^\T\mB^{-1}=(5,0)$。目标行也直接给出同样的结果（{numref}`tab-03-05-reading-duals`）：人工变量 $x_4$ 的费用为 $-M=-100$，列为 $\ve_1$，因此 $-105=-100-\pi_1$，$\pi_1=5$；松弛变量 $x_5$ 对应的元素为 $0=-\pi_2$，因此 $\pi_2=0$。

**图解法。** {numref}`fig-08-03-big-m-dual-region-1`画出了对偶可行域。其中 $x_2$ 的约束是冗余的。可行域有两个顶点：$(5,0)$，$w=150$；以及 $(\tfrac{33}{8},\tfrac78)$（$x_1$ 与 $x_3$ 的约束直线的交点），$w=\tfrac{635}{4}=158.75$；可行域沿方向 $(3,1)$ 无界，而沿该方向 $w$ 增大。最小值在 $\vpi^*=(5,0)$ 处取得，$w^*=150$。

**验证。** 从单纯形表读出原问题的最优解为 $\vx^*=(30,0,0)$，$z^*=150=w^*$。互补松弛性成立：$x_1^*>0$，且 $x_1$ 的约束是紧的（$5+0=5$）；原问题第二行的松弛量为 $40-30=10$，而 $\pi_2^*=0$。

```{figure} figures/08-03-big-m-dual-region-1.svg
:name: fig-08-03-big-m-dual-region-1
:width: 500px
:class: dark-light
:alt: pi1-pi2 平面上由 x1 和 x3 的约束直线围成的对偶可行域，最优点在 (5, 0)

{prf:ref}`ex-08-03-a4-q3` 的对偶问题的可行域。直线 ($x_j$) 是 $x_j$ 的对偶约束；虚线是等值线 $w=150$，它只在 $\vpi^*=(5,0)$ 处与可行域接触。
```
:::

:::{prf:example} 由大 M 单纯形表求对偶解（二）
:label: ex-08-03-a4-q4

问题

$$
\begin{aligned}
\max_{\vx}\ & z=x_1+5x_2+3x_3\\
\st\ & x_1+2x_2+x_3=3\\
& 2x_1-x_2=4\\
& \vx\ge\vzero
\end{aligned}
$$

用大 M 法求解，取 $M=100$，从第一行的 $x_3$ 和第二行的人工变量 $x_4$ 出发。最优单纯形表为

```{list-table}
:header-rows: 1
:stub-columns: 1

* - 基变量
  - $x_1$
  - $x_2$
  - $x_3$
  - $x_4$
  - 右端项
* - $x_3$
  - $0$
  - $\tfrac52$
  - $1$
  - $-\tfrac12$
  - $1$
* - $x_1$
  - $1$
  - $-\tfrac12$
  - $0$
  - $\tfrac12$
  - $2$
* - $-z$
  - $0$
  - $-2$
  - $0$
  - $-99$
  - $-5$
```

写出对偶问题，并用两种方法求出它的最优解。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**对偶问题。** 两行都是等式，因此两个对偶变量都是自由变量：

$$
\begin{aligned}
\min_{\vpi}\ & w=3\pi_1+4\pi_2\\
\st\ & \pi_1+2\pi_2\ge1 &&(x_1)\\
& 2\pi_1-\pi_2\ge5 &&(x_2)\\
& \pi_1\ge3 &&(x_3)\\
& \pi_1,\pi_2\ \text{自由}.
\end{aligned}
$$

（$\pi_1$ 原则上是自由的，但 $x_3$ 的约束使它至少为 3。）

**由单纯形表求。** 初始基为 $(x_3,x_4)$，因此最终表中 $x_3$ 和 $x_4$ 的列就是 $\mB^{-1}=\begin{pmatrix}1&-\tfrac12\\0&\tfrac12\end{pmatrix}$；对 $(x_3,x_1)$ 取 $\vc_B=(3,1)$，得 $\vpi^{*\T}=\vc_B^\T\mB^{-1}=\bigl(3,\ -\tfrac32+\tfrac12\bigr)=(3,\ -1)$。若改为读目标行：$x_3$ 的列为 $\ve_1$，费用为 3，因此 $0=3-\pi_1$，$\pi_1=3$；人工变量 $x_4$ 的列为 $\ve_2$，费用为 $-100$，因此 $-99=-100-\pi_2$，$\pi_2=-1$。

**图解法。** 在{numref}`fig-08-03-big-m-dual-region-2`中，可行域位于 $\pi_1=3$ 的右侧、$x_1$ 与 $x_2$ 的约束直线之间。它的顶点为 $(3,-1)$，$w=5$，以及 $(3,1)$，$w=13$；在可行域上 $w$ 随 $\pi_1$ 增大而增大（在下边界上 $w=5+(\pi_1-3)$）。最优解为 $\vpi^*=(3,-1)$，$w^*=5$。

**验证。** 原问题的最优解为 $\vx^*=(2,0,1)$：$2+1=3$，$4=4$，且 $z^*=2+3=5=w^*$。

*（这张单纯形表的早期版本把 $x_1$ 的列写成 $(1,0)$，把右下角元素写成 $+5$。基变量的列必须是在自己所在行取 1 的单位向量，因此 $x_1$ 的列是 $(0,1)$，而 $-z$ 行的右下角元素是 $-z^*=-5$。）*

```{figure} figures/08-03-big-m-dual-region-2.svg
:name: fig-08-03-big-m-dual-region-2
:width: 500px
:class: dark-light
:alt: 含自由变量的对偶可行域，位于 pi1 = 3 的右侧，最优点在 (3, -1)

{prf:ref}`ex-08-03-a4-q4` 的对偶问题的可行域。两个对偶变量都是自由变量，因此可行域延伸到 $\pi_2=0$ 以下；等值线 $w=5$（虚线）在 $\vpi^*=(3,-1)$ 处与可行域接触。
```
:::

:::{prf:example} 重构单纯形表
:label: ex-08-03-a4-q6

对问题

$$
\begin{aligned}
\max_{\vx}\ & z=5x_1+2x_2+3x_3\\
\st\ & x_1+5x_2+2x_3\le b_1\\
& x_1-5x_2-6x_3\le b_2\\
& \vx\ge\vzero ,
\end{aligned}
$$

引入松弛变量 $x_4$ 和 $x_5$ 后，对 $b_1$ 和 $b_2$ 的某组取值，最优单纯形表为

```{list-table}
:header-rows: 1
:stub-columns: 1

* - 基变量
  - $x_1$
  - $x_2$
  - $x_3$
  - $x_4$
  - $x_5$
  - 右端项
* - $x_1$
  - $1$
  - $b$
  - $2$
  - $1$
  - $0$
  - $3$
* - $x_5$
  - $0$
  - $c$
  - $-8$
  - $-1$
  - $1$
  - $1$
* - $-z$
  - $0$
  - $a$
  - $-7$
  - $d$
  - $e$
  - $-15$
```

求 **(a)** $b_1$ 和 $b_2$，**(b)** 对偶最优解，**(c)** 元素 $a,b,c,d,e$。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

松弛变量列 $x_4,x_5$ 构成初始的单位矩阵，因此它们在最终表中的列给出基 $(x_1,x_5)$ 的 $\mB^{-1}=\begin{pmatrix}1&0\\-1&1\end{pmatrix}$。

**(a)** $\mB^{-1}\vb=(3,1)$ 意味着 $b_1=3$ 且 $-b_1+b_2=1$，因此 $b_1=3$，$b_2=4$。

**(b)** $\vpi^{*\T}=\vc_B^\T\mB^{-1}=(5,0)\mB^{-1}=(5,0)$。验证：$\vb^\T\vpi^*=15=z^*$，与右下角元素 $-15$ 相符。

**(c)** $x_2$ 的列为 $\mA_{\cdot2}=(5,-5)$，因此 $\mB^{-1}\mA_{\cdot2}=(5,\ -5-5)=(5,-10)$：$b=5$，$c=-10$。它的检验数为 $a=c_2-\vpi^\T\mA_{\cdot2}=2-25=-23$。松弛变量 $x_4$ 的 $d=0-\pi_1=-5$；$x_5$ 是基变量，因此 $e=0$。给定的元素与此一致：$\bar c_3=3-5\cdot2=-7$，$\mB^{-1}\mA_{\cdot3}=(2,\ -6-2)=(2,-8)$。所有检验数都 $\le0$，因此对这个最大化问题该表是最优的，$\vx^*=(3,0,0)$，$z^*=15$。

*（本题解答的早期版本给出 $d=-7$ 和 $e=-5$：实际上 $d$ 是第一行对偶值的相反数，即 $-5$，而基变量的检验数为 0。本题的松弛变量是 $x_4$ 和 $x_5$，而不是早期题目叙述中所说的 $x_3$ 和 $x_4$。）*
:::

(sec-08-03-shadow-prices)=
## 影子价格与经济解释

:::{prf:example} 电缆生产中的影子价格
:label: ex-08-03-a4-q8

某厂生产四种电缆。每单位电缆依次经过四道工序：接线、焊接、套管和检验。各道工序所需的时间、每天可用的时间以及每单位的收入如下：

```{list-table}
:header-rows: 1
:stub-columns: 1

* - 电缆类型
  - 接线（分钟/单位）
  - 焊接（分钟/单位）
  - 套管（分钟/单位）
  - 检验（分钟/单位）
  - 收入（美元/单位）
* - 1
  - 10.5
  - 20.4
  - 3.2
  - 5
  - 9.40
* - 2
  - 9.3
  - 24.6
  - 2.5
  - 5
  - 10.80
* - 3
  - 11.6
  - 17.7
  - 3.6
  - 5
  - 8.75
* - 4
  - 8.2
  - 26.5
  - 5.5
  - 5
  - 7.80
* - 能力（分钟/日）
  - 4800
  - 9600
  - 4700
  - 4500
  -
```

一份合同要求每种电缆每天至少生产 100 单位。

**(a)** 建立使每日收入最大的线性规划模型，并求出最优生产计划。

**(b)** 根据对偶价格，你是否建议增加某道工序的每日能力？

**(c)** 最低产量要求对厂家是有利还是不利？用对偶价格加以解释。

**(d)** 如果焊接能力增加 10%，能否保证每增加一分钟获得对偶价格所给出的收入？
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** 设 $x_j$ 为第 $j$ 类电缆的日产量（单位/日）：

$$
\begin{aligned}
\max_{\vx}\ & z=9.4x_1+10.8x_2+8.75x_3+7.8x_4\\
\st\ & 10.5x_1+9.3x_2+11.6x_3+8.2x_4\le4800 &&\text{（接线）}\\
& 20.4x_1+24.6x_2+17.7x_3+26.5x_4\le9600 &&\text{（焊接）}\\
& 3.2x_1+2.5x_2+3.6x_3+5.5x_4\le4700 &&\text{（套管）}\\
& 5x_1+5x_2+5x_3+5x_4\le4500 &&\text{（检验）}\\
& x_j\ge100 && j=1,\dots,4 .
\end{aligned}
$$

最优计划是第 1、2、4 类电缆只生产最低产量，并把剩余的全部焊接时间用于第 3 类：

$$
\begin{aligned}
x_1^*&=x_2^*=x_4^*=100,\\
x_3^*&=\frac{9600-100\,(20.4+24.6+26.5)}{17.7}=\frac{2450}{17.7}\\
&=\frac{24500}{177}\approx138.42,
\end{aligned}
$$

每天 $z^*=\tfrac{709975}{177}\approx\$4011.16$。焊接能力被完全用尽；接线有 $394.35$ 分钟的松弛，套管有 $3081.69$ 分钟，检验有 $2307.91$ 分钟。（该计划是分数解；可以把它理解为平均日产率，或把 $x_3$ 向下取整为 138 单位，这仍然是可行的。）

**(b)** 由互补松弛性，三道有松弛的工序的对偶价格为 0。焊接的价格由基变量 $x_3$ 的对偶约束得出，该约束必须是紧的：$17.7\,\pi_2=8.75$，因此

$$
\pi_2^*=\frac{8.75}{17.7}=\frac{175}{354}\approx\$0.4944\ \text{/分钟}.
$$

只有焊接能力值得增加：每增加一分钟约增加 0.49 美元的收入（只要基保持最优，见 (d)），因此只要每增加一分钟的成本低于这个数，增加焊接能力就是合算的。额外的接线、套管或检验时间毫无价值。

**(c)** 下界 $x_j\ge100$ 的对偶值就是检验数 $\bar c_j=c_j-\pi_2^*a_{2j}$，其中 $a_{2j}$ 是第 $j$ 类电缆的焊接时间：

$$
\begin{aligned}
\bar c_1&=9.4-0.4944\cdot20.4\approx-0.6847,\\
\bar c_2&\approx-1.3610,\\
\bar c_3&=0,\\
\bar c_4&\approx-5.3003\qquad\text{（美元/单位）}.
\end{aligned}
$$

第 1、2、4 类电缆的收入低于它们所占用焊接时间的价值（对第 4 类，7.80 美元对 $26.5\cdot0.4944=\$13.10$）。因此它们的最低产量要求每增加一个单位都会损失收入：把第 1、2 或 4 类的最低产量降低一个单位，会使 $z^*$ 分别增加 0.68 美元、1.36 美元或 5.30 美元。最低产量要求是不利的。没有这些要求时，最优解为 $\vx\approx(0,\ 218.63,\ 238.51,\ 0)$，每天 $z\approx\$4448.20$，因此这份合同每天约使收入损失 437.04 美元。

**(d)** 只要基保持最优，即只要 $x_3=(b_2-7150)/17.7$ 使所有约束都满足，价格 $\pi_2^*$ 就有效。记 $b_2$ 为焊接能力：

- $x_3\ge100$ 要求 $b_2\ge7150+1770=8920$ 分钟；
- 接线约束 $10.5\cdot100+9.3\cdot100+8.2\cdot100+11.6\,x_3\le4800$ 要求 $x_3\le\tfrac{5000}{29}\approx172.41$，即 $b_2\le\tfrac{295850}{29}\approx10201.72$ 分钟；
- 套管和检验约束允许 $x_3$ 分别达到约 994 和 600 单位，不构成限制。

因此每分钟 0.4944 美元的价格在 $8920\le b_2\le10201.72$ 分钟范围内有效，最多允许增加 $601.72$ 分钟，即 6.27%。增加 10%，即增加到 10 560 分钟，超出了这个范围，所以答案是**不能**：对整个增量而言，这个价格得不到保证。事实上，超过 10 201.72 分钟后接线约束也成为紧约束，焊接的价格降到每分钟约 0.3636 美元；在 10 560 分钟时最优收入为 4438.89 美元，增加了 427.73 美元，而不是价格所预测的 $960\cdot0.4944=\$474.58$。

*（本题解答的早期版本把该范围的上端写为 20 101.72 分钟；正确值为 10 201.72 分钟。最大增幅为 6.27%。）*
:::

:::{prf:example} 作为影子价格的拉格朗日乘子
:label: ex-08-03-nlp-multipliers

非线性问题中也有影子价格，即拉格朗日乘子（Lagrange multiplier）。考虑

$$
\begin{aligned}
\max_{\vx}\ & f(\vx)=(x_1-2)^2+(x_2-2)^2\\
\st\ & g_1(\vx)=(x_1-1)^2+(x_2-1)^2-1\le0\\
& g_2(\vx)=x_1+x_2-1\le0 .
\end{aligned}
$$

**(a)** 求最优解。**(b)** 求乘子 $\lambda_1,\lambda_2\ge0$，使 $\nabla f(\vx^*)=\lambda_1\nabla g_1(\vx^*)+\lambda_2\nabla g_2(\vx^*)$ 且 $\lambda_ig_i(\vx^*)=0$。**(c)** 证明：当第一个约束放松为 $g_1(\vx)\le\varepsilon$ 时，$\lambda_1$ 就是最优值的增长率。
:::

:::{dropdown} 解答
:class-container: solution
:animate: fade-in

**(a)** $f$ 是到点 $(2,2)$ 的距离的平方，第一个约束是以 $(1,1)$ 为圆心、半径为 1 的圆盘。圆盘上离 $(2,2)$ 最远的点位于经过 $(2,2)$ 和圆心的直线上、圆心的另一侧：$\vx^*=(1,1)-\tfrac{1}{\sqrt2}(1,1)=\bigl(1-\tfrac{1}{\sqrt2},\ 1-\tfrac{1}{\sqrt2}\bigr)\approx(0.2929,\ 0.2929)$。它满足第二个约束 $x_1+x_2=2-\sqrt2\approx0.586\le1$，因此对更小的可行集它也是最优的。它到 $(2,2)$ 的距离为 $\sqrt2+1$，因此 $f^*=(1+\sqrt2)^2=3+2\sqrt2\approx5.828$。

**(b)** 第二个约束不是有效约束（$g_2(\vx^*)=-0.414<0$），因此 $\lambda_2=0$。在 $\vx^*$ 处的梯度为 $\nabla f=2(\vx^*-(2,2))=-2\bigl(1+\tfrac{1}{\sqrt2}\bigr)(1,1)$ 和 $\nabla g_1=2(\vx^*-(1,1))=-\tfrac{2}{\sqrt2}(1,1)$，因此

$$
\lambda_1=\frac{1+1/\sqrt2}{1/\sqrt2}=\sqrt2+1\approx2.414 .
$$

两个条件在 $\lambda_1,\lambda_2\ge0$ 下都成立。（这就是{prf:ref}`thm-06-03-kkt`中针对最大化问题写出的 KKT 条件（Karush–Kuhn–Tucker conditions）。由于 $f$ 是凸函数而我们在对它求最大值，这些条件在这里是必要的但不是充分的；最优性已在 (a) 中用几何方法证明。）

**(c)** 放松为 $g_1\le\varepsilon$ 会使圆盘的半径增大到 $\sqrt{1+\varepsilon}$，最远点沿同一条直线向外移动，而当 $\varepsilon$ 很小时第二个约束仍然是非有效约束，因此

$$
\begin{gathered}
f^*(\varepsilon)=\bigl(\sqrt2+\sqrt{1+\varepsilon}\bigr)^2,\\
\frac{\mathrm d f^*}{\mathrm d\varepsilon}\Big|_{\varepsilon=0}
=2\bigl(\sqrt2+1\bigr)\cdot\frac12=\sqrt2+1=\lambda_1 .
\end{gathered}
$$

乘子就是该约束的影子价格：把约束放松 $\varepsilon=\pm0.001$ 会使 $f^*$ 改变约 $\pm0.00241$。非有效约束的影子价格为 $\lambda_2=0$，与线性规划中一样。

*（本问题 KKT 条件的早期表述中，第二个约束的互补条件写成了含 $g_1$ 的形式；它应为 $\lambda_2g_2(\vx^*)=0$。）*
:::

## 参考文献

- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Sections 4.2–4.4: the dual problem,
  duality theorems, dual variables as marginal costs.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Sections 6.1–6.4 and 7.1–7.2: duality
  theory and its use in sensitivity analysis.)
- Karush, W. (1939). *Minima of functions of several variables with
  inequalities as side constraints* [Master's thesis, University of
  Chicago].
- Kuhn, H. W., & Tucker, A. W. (1951). Nonlinear programming. In J. Neyman
  (Ed.), *Proceedings of the Second Berkeley Symposium on Mathematical
  Statistics and Probability* (pp. 481–492). University of California Press.
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 4: Sections 4.1–4.3 and the problem sets of Sections 4.1,
  4.2 and 4.5, from which the cable-production problem, the big-M dual
  problems and the bounded-variable dual are adapted; Section 7.4: duality in
  matrix form.)
