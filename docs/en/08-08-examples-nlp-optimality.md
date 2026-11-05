# Optimality Conditions and KKT

This chapter works through the examples and problems of the first three
chapters on nonlinear programming: {doc}`06-01-nlp-models`,
{doc}`06-02-optimality-conditions` and {doc}`06-03-kkt-conditions`.

The items are grouped by topic: formulating and classifying nonlinear
programs, local and global optima, feasible and descent directions, first-
and second-order conditions, constraint qualifications, and the Karush–Kuhn–Tucker (KKT)
conditions. Every problem is brought to the form $\min f$, $g_i\le0$,
$h_j=0$ of {eq}`eq-06-01-nlp` before the KKT conditions are written. Each
solution states the result it uses, carries out the computation by hand and
checks the answer; every optimum was also confirmed numerically with a
general-purpose solver (`scipy.optimize.minimize`, see {doc}`09-06-code-nlp`).

(sec-08-08-formulation)=
## Formulation and classification

:::{prf:example} A convex quadratic program on a triangle
:label: ex-08-08-convex-qp

Consider

$$
\begin{aligned}
\min\ & f(\vx)=2x_1+x_2+x_1^2+2x_1x_2+x_2^2\\
\st\ & x_1+x_2\le2,\\
& x_1,\ x_2\ge0 .
\end{aligned}
$$

(a) Write $f$ in the form $\vc^\T\vx+\tfrac12\vx^\T\mQ\vx$ with a symmetric
$\mQ$ and show that the problem is a convex program.
(b) Find the optimal solution and verify it with the KKT conditions.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** The linear part gives $\vc=(2,1)^\T$. The quadratic part
$x_1^2+2x_1x_2+x_2^2$ gives $Q_{11}=2\cdot1$, $Q_{22}=2\cdot1$ and
$Q_{12}=Q_{21}=2$:

$$
\begin{aligned}
\mQ&=\begin{pmatrix}2&2\\2&2\end{pmatrix},\\
\lambda_{1,2}(\mQ)&=0,\ 4 .
\end{aligned}
$$

Both eigenvalues are nonnegative, so $\mQ\succeq\vzero$ and $f$ is convex
(not strictly, since $\mQ$ is singular). The feasible set is a polyhedron,
hence convex. The problem is a convex quadratic program (QP;
{ref}`sec-06-01-classes`).

**(b)** Write $f(\vx)=(x_1+x_2)^2+2x_1+x_2$. On the feasible set every term
is nonnegative, so $f(\vx)\ge0=f(0,0)$, and $f(\vx)=0$ forces
$x_1=x_2=0$. The unique optimum is $\vx^*=(0,0)$ with $f^*=0$.

KKT check with $g_1=x_1+x_2-2$, $g_2=-x_1$, $g_3=-x_2$: at the origin
$g_1=-2<0$, so $\lambda_1=0$, while $g_2$ and $g_3$ are active. The gradient
is $\nabla f(\vx)=(2+2x_1+2x_2,\ 1+2x_1+2x_2)^\T$, so
$\nabla f(\vzero)=(2,1)^\T$, and stationarity

$$
\begin{pmatrix}2\\1\end{pmatrix}+\lambda_2\begin{pmatrix}-1\\0\end{pmatrix}+\lambda_3\begin{pmatrix}0\\-1\end{pmatrix}=\vzero
$$

gives $\lambda_2=2\ge0$, $\lambda_3=1\ge0$. The KKT conditions hold, and by
{prf:ref}`thm-06-03-convex-kkt` the origin is a global minimum. Note that the
optimum is a vertex here, but only because $\nabla f(\vzero)$ points into
the feasible set; the gradient decides, not the vertex structure. (For
contrast, the *maximum* of the convex $f$ over the triangle is at the vertex
$(2,0)$ with $f=8$; a convex function attains its maximum over a polytope at
a vertex.)
:::

(sec-08-08-local-global)=
## Local and global optima

:::{prf:example} An unattained supremum and a disconnected feasible set
:label: ex-08-08-difficulties

(a) Let $f(x)=1-1/x$ on $\Omega=\{x:x\ge0.5\}$. What is
$\max_{x\in\Omega}f(x)$? Which hypothesis of {prf:ref}`thm-06-01-existence`
fails?

(b) For $x>0$ let $g_1(x)=(x-2)^2/x$ and $g_2(x)=2.5-0.3x$. Describe the
set $\Omega=\{x>0:g_1(x)\ge g_2(x)\}$ and explain why it is a difficulty for
search methods.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** $f'(x)=1/x^2>0$, so $f$ is increasing on $\Omega$. For every
$x\in\Omega$, $f(x)<1$, and $f(x)\to1$ as $x\to\infty$ (for example
$f(100)=0.99$). Hence $\sup_{\Omega}f=1$, but no feasible point attains
it: **the problem has no maximum.** The minimum exists, at the end point
$x=0.5$ with $f(0.5)=-1$. Weierstrass' theorem does not apply because
$\Omega$ is closed but unbounded, and $-f$ is not coercive on $\Omega$
($-f\to-1$, not $+\infty$). Bounded objective values alone do not
guarantee an optimum.

**(b)** Multiplying $g_1\ge g_2$ by $x>0$:

$$
(x-2)^2\ge x(2.5-0.3x)\iff 1.3x^2-6.5x+4\ge0 .
$$

The roots of $1.3x^2-6.5x+4$ are
$x=\bigl(6.5\mp\sqrt{6.5^2-4\cdot1.3\cdot4}\bigr)/2.6=(6.5\mp\sqrt{21.45})/2.6$,
that is $x\approx0.7187$ and $x\approx4.2813$. The quadratic is nonnegative
outside the roots, so

$$
\Omega=(0,\ 0.7187]\ \cup\ [4.2813,\ \infty),
$$

two separate intervals. (Check at $x=2$: $g_1=0<g_2=1.9$, so $2\notin\Omega$.)
A method that moves continuously through feasible points, as the descent
methods of {doc}`06-05-unconstrained-algorithms` and
{doc}`06-06-constrained-algorithms` do, stays in the piece where it starts and
never sees the other one. To find a global optimum each component must be
searched separately, and in higher dimensions the components are not even
known in advance. {numref}`fig-08-08-difficulties` shows both situations.

```{figure} figures/08-08-difficulties.svg
:name: fig-08-08-difficulties
:width: 680px
:class: dark-light
:alt: Left, the increasing curve 1 - 1/x for x >= 0.5 approaching the dashed asymptote at height 1; right, a convex curve and a decreasing line crossing twice, with the two parts of the x-axis where the curve lies above the line highlighted

Left: $f(x)=1-1/x$ on $x\ge0.5$ approaches $1$ but never reaches it
({prf:ref}`ex-08-08-difficulties` (a)). Right: $\{x:g_1(x)\ge g_2(x)\}$
consists of the two highlighted intervals ({prf:ref}`ex-08-08-difficulties`
(b)).
```
:::

(sec-08-08-directions)=
## Feasible and descent directions

:::{prf:example} Feasible descent directions at a vertex
:label: ex-08-08-descent-cone

Consider

$$
\begin{aligned}
\min\ & f(\vx)=(x_1-6)^2+(x_2-7)^2\\
\st\ & x_1+x_2\le5,\\
& x_1-x_2\le1,\\
& x_1\ge0 .
\end{aligned}
$$

(a) Give two descent directions at the feasible point $(1,1)$.
(b) Find three feasible descent directions at the feasible point $(3,2)$
and explain why they qualify.
(c) Find the optimal solution and its multipliers.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The feasible set is the triangle with vertices $(0,5)$, $(3,2)$ and
$(0,-1)$, and $\nabla f(\vx)=(2(x_1-6),\ 2(x_2-7))^\T$.

**(a)** At $(1,1)$: $1+1=2<5$, $1-1=0<1$ and $1>0$, so the point is
interior and every direction is feasible
({prf:ref}`prop-06-02-feasible-directions`). With
$\nabla f(1,1)=(-10,-12)^\T$, a direction is a descent direction when
$-10d_1-12d_2<0$. For example

$$
\begin{aligned}
\vd&=(1,2):&\nabla f^\T\vd&=-10-24=-34<0,\\
\vd&=(1,0):&\nabla f^\T\vd&=-10<0 .
\end{aligned}
$$

Both are feasible descent directions.

**(b)** At $(3,2)$ the constraints $x_1+x_2\le5$ and $x_1-x_2\le1$ are
active and $x_1\ge0$ is not. All constraints are linear, so by
{prf:ref}`prop-06-02-feasible-directions` the feasible directions are

$$
\begin{aligned}
&d_1+d_2\le0,\\
&d_1-d_2\le0
\end{aligned}
\quad\Longleftrightarrow\quad d_1\le-\abs{d_2} .
$$

The gradient is $\nabla f(3,2)=(-6,-10)^\T$, and $\vd$ is a descent
direction when $-6d_1-10d_2<0$, that is $d_2>-0.6\,d_1$. Since feasibility
forces $d_1\le0$, both conditions together read

$$
0.6\abs{d_1}<d_2\le\abs{d_1},
$$

so, up to positive scaling, the feasible descent directions are
$\vd=(-1,s)$ with $0.6<s\le1$. Three of them:

```{list-table} Feasible descent directions at $(3,2)$.
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

The first points along the edge $x_1+x_2=5$ towards $(0,5)$; the other two
point into the triangle. The cone is narrow
({numref}`fig-08-08-descent-cone`): it lies between the edge direction
$(-1,1)$ and the direction $(-1,0.6)$ orthogonal to the gradient.
*(An earlier version of this solution also lists $(-2,2)$, which is the
same direction as $(-3,3)$; it writes the descent condition with a
misprinted index, $d_{22}$ for the second component; and it bounds the
first component by $-3/\alpha$, although the constraint $x_1\ge0$ is
inactive at $(3,2)$ and does not restrict directions.)*

**(c)** Since feasible descent directions exist at $(1,1)$ and $(3,2)$,
neither is optimal. The objective is the squared distance to $(6,7)$, which
lies beyond the edge $x_1+x_2=5$; projecting it onto that line gives
$(6,7)-\tfrac{13-5}{2}(1,1)=(2,3)$, which satisfies $2-3=-1\le1$ and
$2\ge0$. KKT check: only $g_1=x_1+x_2-5$ is active, and
$\nabla f(2,3)=(-8,-8)^\T=-8\nabla g_1$, so $\lambda_1=8\ge0$ and
$\lambda_2=\lambda_3=0$. The problem is convex, so $\vx^*=(2,3)$ with
$f^*=16+16=32$ is the global minimum ({prf:ref}`thm-06-03-convex-kkt`).

```{figure} figures/08-08-descent-cone.svg
:name: fig-08-08-descent-cone
:width: 520px
:class: dark-light
:alt: The triangle with vertices (0,5), (3,2), (0,-1), the gradient at (3,2) pointing down-left, its negative pointing up-right, and a narrow shaded cone of feasible descent directions at (3,2) between the edge towards (0,5) and the line orthogonal to the gradient

Feasible descent directions at the vertex $(3,2)$ of
{prf:ref}`ex-08-08-descent-cone`. Feasible directions point into the triangle;
descent directions make an obtuse angle with $\nabla f(3,2)=(-6,-10)$. Their
intersection is the narrow shaded cone between the edge direction $(-1,1)$
and the direction $(-1,0.6)$ orthogonal to the gradient (dotted). The
optimum $\vx^*=(2,3)$ lies on the edge.
```
:::

(sec-08-08-first-second)=
## First- and second-order conditions

:::{prf:example} A quadratic function of three variables
:label: ex-08-08-quadratic-three

Let

$$
f(\vx)=3x_1^2+2x_1x_2+6x_2^2-x_2x_3+5x_3^2-4x_1-2x_2-6x_3+9 .
$$

(a) Find the stationary point of $f$. (b) Show that it is a local minimum.
(c) Show that it is the global minimum, and give the minimum value.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** Setting the gradient to zero,

$$
\nabla f(\vx)=\begin{pmatrix}6x_1+2x_2-4\\2x_1+12x_2-x_3-2\\-x_2+10x_3-6\end{pmatrix}=\vzero .
$$

The first and third equations give $x_1=(2-x_2)/3$ and $x_3=(6+x_2)/10$.
Substituting into the second and multiplying by 30:

$$
20(2-x_2)+360x_2-3(6+x_2)-60=0\iff337x_2=38 .
$$

Hence

$$
\vx^*=\Bigl(\tfrac{212}{337},\ \tfrac{38}{337},\ \tfrac{206}{337}\Bigr)\approx(0.6291,\ 0.1128,\ 0.6113).
$$

**(b)** The Hessian is constant,

$$
\nabla^2f=\begin{pmatrix}6&2&0\\2&12&-1\\0&-1&10\end{pmatrix},
$$

with leading principal minors $6>0$, $6\cdot12-2\cdot2=68>0$ and
$\det\nabla^2f=6(120-1)-2(20-0)+0=674>0$. A symmetric matrix with positive
leading principal minors is positive definite ({prf:ref}`thm-01-01-sylvester`);
its eigenvalues are $5.376$, $9.704$ and $12.92$. By
{prf:ref}`thm-06-02-second-order-sufficient`, $\vx^*$ is a strict local
minimum.

**(c)** Since $\nabla^2f\succ\vzero$ everywhere, $f$ is strictly convex on
$\R^3$, and by {prf:ref}`thm-06-02-convex-sufficient` the stationary point is
the unique global minimum. With $f(\vx)=\tfrac12\vx^\T\mH\vx-\vb^\T\vx+9$,
$\mH=\nabla^2f$, $\vb=(4,2,6)^\T$ and $\mH\vx^*=\vb$, the minimum value is

$$
\begin{aligned}
f^*&=9-\tfrac12\vb^\T\vx^*=9-\frac{4\cdot212+2\cdot38+6\cdot206}{2\cdot337}\\
&=9-\frac{1080}{337}=\frac{1953}{337}\approx5.795 .
\end{aligned}
$$

A quasi-Newton solver started at the origin returns the same point and
value.
:::

:::{prf:example} The egg-crate surface
:label: ex-08-08-egg-crate

The surface $f(\vx)=\sin x_1\sin x_2$ on $[-4,4]^2$ has alternating humps
and hollows. Show that the origin is a stationary point but not a local
minimum, and classify all stationary points in the square.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The gradient and Hessian are

$$
\begin{aligned}
\nabla f(\vx)&=\begin{pmatrix}\cos x_1\sin x_2\\ \sin x_1\cos x_2\end{pmatrix},\\
\nabla^2f(\vx)&=\begin{pmatrix}-\sin x_1\sin x_2&\cos x_1\cos x_2\\ \cos x_1\cos x_2&-\sin x_1\sin x_2\end{pmatrix}.
\end{aligned}
$$

**The origin.** $\nabla f(\vzero)=\vzero$, and
$\nabla^2f(\vzero)=\bigl(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\bigr)$
has eigenvalues $\pm1$: it is indefinite, so the second-order necessary
condition fails ({prf:ref}`thm-06-02-second-order-necessary`) and the origin
is a saddle point. Directly: $f(t,t)=\sin^2t>0$ and $f(t,-t)=-\sin^2t<0$
for small $t\ne0$, while $f(\vzero)=0$. The first-order condition alone
cannot distinguish this point from a minimum.

**All stationary points.** $\cos x_1\sin x_2=0$ and $\sin x_1\cos x_2=0$.
If $\cos x_1=0$ then $\sin x_1=\pm1$ and the second equation forces
$\cos x_2=0$; if instead $\sin x_2=0$ then $\cos x_2=\pm1$ and
$\sin x_1=0$. In $[-4,4]$ the zeros of $\cos$ are $\pm\pi/2$ and those of
$\sin$ are $0,\pm\pi$ ($\pm3\pi/2\approx\pm4.71$ lie outside). So there
are $4+9=13$ stationary points:

- $(\pm\pi/2,\pm\pi/2)$: $\nabla^2f=-\sin x_1\sin x_2\,\mI$. Where
  $\sin x_1\sin x_2=1$, that is at $(\pi/2,\pi/2)$ and $(-\pi/2,-\pi/2)$,
  the Hessian is $-\mI\prec\vzero$: strict local maxima with $f=1$. Where
  $\sin x_1\sin x_2=-1$ the Hessian is $\mI\succ\vzero$: strict local minima
  with $f=-1$.
- The nine points with $x_1,x_2\in\{-\pi,0,\pi\}$: the diagonal of the
  Hessian is zero and the off-diagonal entry is $\pm1$, so the eigenvalues
  are $\pm1$: all nine are saddle points.
:::

:::{prf:example} First-order condition at a bound
:label: ex-08-08-bound-first-order

Minimize $f(x)=(x-2)^2$ subject to $0\le x\le1$ and show that the
first-order necessary condition holds at the minimum.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The stationary point $x=2$ of $f$ is infeasible. On $[0,1]$,
$f'(x)=2(x-2)<0$, so $f$ is decreasing and the minimum is at the right end:
$x^*=1$, $f^*=1$ (compare $f(0)=4$).

At $x^*=1$ the feasible directions are $d\le0$ (one can only move left),
and

$$
f'(x^*)\,d=-2d\ge0\qquad\text{for all }d\le0,
$$

so {eq}`eq-06-02-first-order` holds although $f'(1)\ne0$. In KKT form, with
$g_1=-x\le0$ (inactive, $\lambda_1=0$) and $g_2=x-1\le0$ (active):
$f'(1)+\lambda_2=0$ gives $\lambda_2=2\ge0$. The problem is convex, so the
KKT point is the global minimum. At the other end point the condition
fails: at $x=0$ the feasible directions are $d\ge0$ and $f'(0)d=-4d<0$ for
$d>0$, so $0$ is not a minimum.
:::

(sec-08-08-cq)=
## Constraint qualifications

:::{prf:example} A minimum without KKT multipliers
:label: ex-08-08-cusp

Consider

$$
\begin{aligned}
\min\ & f(\vx)=-x_1\\
\st\ & g_1(\vx)=x_2-(1-x_1)^3\le0,\\
& g_2(\vx)=-x_2\le0 .
\end{aligned}
$$

(a) Show that $\vx^*=(1,0)$ is the unique optimal solution.
(b) Show that no KKT multipliers exist at $\vx^*$.
(c) Which constraint qualification fails? Compare $L(\vx^*)$ with
$T(\vx^*)$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** Feasibility requires $0\le x_2\le(1-x_1)^3$, hence $(1-x_1)^3\ge0$
and $x_1\le1$. So $f=-x_1\ge-1$ on $\Omega$, with equality only for
$x_1=1$, where $0\le x_2\le0$ forces $x_2=0$. The unique optimum is
$\vx^*=(1,0)$, $f^*=-1$. The feasible set is a cusp pointing to the right
with its tip at $\vx^*$.

**(b)** Both constraints are active at $\vx^*$. The gradients are
$\nabla f=(-1,0)^\T$, $\nabla g_1(\vx)=(3(1-x_1)^2,\ 1)^\T$, so
$\nabla g_1(\vx^*)=(0,1)^\T$, and $\nabla g_2=(0,-1)^\T$. Stationarity

$$
\begin{pmatrix}-1\\0\end{pmatrix}+\lambda_1\begin{pmatrix}0\\1\end{pmatrix}+\lambda_2\begin{pmatrix}0\\-1\end{pmatrix}=\vzero
$$

requires $-1=0$ in the first component: impossible. The global minimum is
not a KKT point. (A numerical solver run from $(0.5,0.05)$ approaches
$(1,0)$, typically with warnings, because its stopping test is based on
the KKT conditions.)

**(c)** The active gradients $(0,1)$ and $(0,-1)$ are linearly dependent,
so the linear independence constraint qualification (LICQ) fails; the constraints are not linear, and $g_1$ is not convex, so
neither of the other sufficient conditions of
{prf:ref}`thm-06-03-cq-sufficient` applies. The linearized cone is

$$
L(\vx^*)=\{\vd:\ d_2\le0,\ -d_2\le0\}=\{(d_1,0):d_1\in\R\},
$$

the whole $x_1$-axis, whereas the feasible set near $\vx^*$ lies to the left
of $\vx^*$, so $T(\vx^*)=\{(d_1,0):d_1\le0\}$. The direction $\vd=(1,0)$
belongs to $L(\vx^*)$ but not to $T(\vx^*)$, and
$\nabla f(\vx^*)^\T\vd=-1<0$: the linearization suggests an improving
direction that does not exist, and Farkas' lemma then denies the existence
of multipliers. The Fritz John conditions ({prf:ref}`rem-06-03-history`) do
hold, with $\lambda_0=0$ and $\lambda_1=\lambda_2=1$. This example is due to
Kuhn and Tucker (1951).
:::

(sec-08-08-kkt)=
## KKT conditions

:::{prf:example} Writing down a KKT system
:label: ex-08-08-kkt-formulation

(a) Formulate the KKT conditions of

$$
\begin{aligned}
\min\ & f(\vx)=4x_1^2+4x_1x_2+2x_2^2+3x_1+e^{x_1+2x_2}\\
\st\ & x_1^2+x_2^2\le10,\\
& -e^{-x_1^2}\le0.5,\\
& x_1\ge0 .
\end{aligned}
$$

(b) Solve them.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** In the form $g_i\le0$: $g_1=x_1^2+x_2^2-10$,
$g_2=-e^{-x_1^2}-0.5$, $g_3=-x_1$. The gradients are

$$
\begin{aligned}
\nabla f&=\begin{pmatrix}8x_1+4x_2+3+e^{x_1+2x_2}\\4x_1+4x_2+2e^{x_1+2x_2}\end{pmatrix},\\
\nabla g_1&=\begin{pmatrix}2x_1\\2x_2\end{pmatrix},\\
\nabla g_2&=\begin{pmatrix}2x_1e^{-x_1^2}\\0\end{pmatrix},\\
\nabla g_3&=\begin{pmatrix}-1\\0\end{pmatrix}.
\end{aligned}
$$

The KKT conditions are

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

Since $-e^{-x_1^2}\in[-1,0)$, the second constraint holds strictly for
every $x_1$: it is never active, $\lambda_2=0$, and it can be dropped.
*(An earlier version of this solution prints the second stationarity
equation with $e^{x_1+2x_2}$ instead of $2e^{x_1+2x_2}$ and with $\lambda_2$
instead of $\lambda_1$; the partial derivative of $e^{x_1+2x_2}$ with respect
to $x_2$ carries the factor 2, and $x_2$ appears only in $g_1$.)*

**(b)** *Convexity.* The quadratic part of $f$ has Hessian
$\bigl(\begin{smallmatrix}8&4\\4&4\end{smallmatrix}\bigr)$ with eigenvalues
$1.528$ and $10.47$, and $e^{x_1+2x_2}$ is the exponential of an affine
function, so $f$ is strictly convex. After dropping $g_2$ the constraints
are convex, and $(1,0)$ is a Slater point. Any KKT point is therefore the
unique global minimum, and the search can stop at the first KKT point
found. The disc $x_1^2+x_2^2\le10$ is large, so try the cases with $g_1$
inactive ($\lambda_1=0$) first.

*Case $x_1>0$, $\lambda_1=0$.* Then $\lambda_3=0$, and with
$E=e^{x_1+2x_2}>0$ the two equations are $8x_1+4x_2+3+E=0$ and
$4x_1+4x_2+2E=0$. The second gives $x_2=-x_1-E/2$; substituting,
$4x_1+3-E=0$, so $E=4x_1+3$ and $x_2=-3x_1-1.5$. But then
$E=e^{x_1+2x_2}=e^{-5x_1-3}\le e^{-3}\approx0.0498$, which contradicts
$E=4x_1+3\ge3$. No KKT point of this kind.

*Case $x_1=0$, $\lambda_1=0$.* The second equation becomes
$4x_2+2e^{2x_2}=0$, that is $e^{2x_2}=-2x_2$. With $\zeta=-2x_2$ this is
$\zeta e^{\zeta}=1$, whose solution is the omega constant $\zeta\approx0.567143$; so
$x_2=-\zeta/2\approx-0.28357$. The first equation gives

$$
\begin{aligned}
\lambda_3&=4x_2+3+e^{2x_2}=4x_2+3-2x_2\\
&=3+2x_2=3-\zeta\approx2.4329>0,
\end{aligned}
$$

and $g_1=x_2^2-10<0$ is consistent with $\lambda_1=0$. All KKT conditions
hold:

$$
\begin{aligned}
\vx^*&\approx(0,\ -0.2836),\\
\vlambda^*&\approx(0,\ 0,\ 2.4329),\\
f^*&=2x_2^2+e^{2x_2}=\tfrac12\zeta^2+\zeta\approx0.7280 .
\end{aligned}
$$

The bound $x_1\ge0$ is the only binding constraint: without it, $f$ would
decrease further for negative $x_1$. Because the problem is convex and a
KKT point has been found, the cases with $g_1$ active ($x_1^2+x_2^2=10$)
need not be examined: they cannot produce a second minimum. A numerical
solver confirms the point and value.
:::

:::{prf:example} A maximization problem: solution and uniqueness
:label: ex-08-08-kkt-max

Consider

$$
\begin{aligned}
\max\ & 14x_1-x_1^2+6x_2-x_2^2\\
\st\ & x_1+x_2\le2,\\
& x_1+2x_2\le3 .
\end{aligned}
$$

(a) Find the optimal solution with the KKT conditions.
(b) Show that the optimal solution is unique.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** Convert to minimization: $\min f(\vx)=x_1^2-14x_1+x_2^2-6x_2$
with $g_1=x_1+x_2-2\le0$, $g_2=x_1+2x_2-3\le0$. The KKT conditions are

$$
\begin{aligned}
2x_1-14+\lambda_1+\lambda_2&=0,\\
2x_2-6+\lambda_1+2\lambda_2&=0,
\end{aligned}
$$

together with $\lambda_1g_1=0$, $\lambda_2g_2=0$, $g_1,g_2\le0$ and
$\lambda_1,\lambda_2\ge0$.

The unconstrained maximizer $(7,3)$ violates both constraints, so at least
one is active. The three cases:

- $S=\{1\}$ ($\lambda_2=0$): $x_1=7-\lambda_1/2$, $x_2=3-\lambda_1/2$, and
  $x_1+x_2=10-\lambda_1=2$ gives $\lambda_1=8$, $\vx=(3,-1)$. Check:
  $g_2=3-2-3=-2\le0$ and $\lambda_1=8\ge0$. **KKT point.**
- $S=\{2\}$ ($\lambda_1=0$): $x_1=7-\lambda_2/2$, $x_2=3-\lambda_2$, and
  $x_1+2x_2=13-2.5\lambda_2=3$ gives $\lambda_2=4$, $\vx=(5,-1)$; but
  $g_1=5-1-2=2>0$. Infeasible.
- $S=\{1,2\}$: $\vx=(1,1)$; then $\lambda_1+\lambda_2=12$ and
  $\lambda_1+2\lambda_2=4$ give $\lambda_2=-8<0$. Rejected.

The only KKT point is $\vx^*=(3,-1)$ with $\vlambda^*=(8,0)$, and the
maximum value is $14\cdot3-9+6\cdot(-1)-1=26$. Since $f$ is convex and the
constraints are linear, it is the global optimum
({prf:ref}`thm-06-03-convex-kkt`). In the original maximization form,
$\nabla(14x_1-x_1^2+6x_2-x_2^2)=(8,8)^\T=8\nabla g_1$: the gradient of the
objective points out of the feasible set along the normal of the active
constraint, as the third row of {numref}`tab-06-03-signs` says.

**(b)** Suppose $\vx^*$ and $\vy^*$ are both optimal. By the first-order
condition on the convex feasible set ({eq}`eq-06-02-variational-inequality`),

$$
\begin{aligned}
\nabla f(\vx^*)^\T(\vy^*-\vx^*)&\ge0,\\
\nabla f(\vy^*)^\T(\vx^*-\vy^*)&\ge0 .
\end{aligned}
$$

Adding, $\bigl(\nabla f(\vx^*)-\nabla f(\vy^*)\bigr)^\T(\vx^*-\vy^*)\le0$.
For this $f$, $\nabla f(\vx)-\nabla f(\vy)=2(\vx-\vy)$, so the inequality
says $2\norm{\vx^*-\vy^*}^2\le0$, hence $\vx^*=\vy^*$. (For a general
$f$ the mean value theorem turns the left-hand side into
$(\vx^*-\vy^*)^\T\nabla^2f(\vxi)(\vx^*-\vy^*)$ for some $\vxi$ on the
segment, and positive definiteness of the Hessian gives the same
conclusion; this is the strict convexity argument of
{prf:ref}`thm-06-01-convex-local-global`.)
:::

:::{prf:example} KKT conditions of a linear program
:label: ex-08-08-kkt-lp

Consider the linear program (LP) $\min\ \vc^\T\vx$ s.t. $\mA\vx=\vb$,
$\vx\ge\vzero$, with $\mA\in\R^{m\times n}$.

(a) Write down its KKT conditions.
(b) Show that they are the LP optimality conditions, in particular
complementary slackness.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** The constraints are $h_i(\vx)=\mA_{i\cdot}\vx-b_i=0$
($i=1,\dots,m$, multipliers $\mu_i$ free) and $g_j(\vx)=-x_j\le0$
($j=1,\dots,n$, multipliers $\lambda_j\ge0$). With $\nabla h_i=\mA_{i\cdot}^\T$
and $\nabla g_j=-\ve_j$, the KKT conditions are

$$
\begin{aligned}
&\vc+\mA^\T\vmu-\vlambda=\vzero,\\
&\mA\vx=\vb,\\
&\vx\ge\vzero,\\
&\vlambda\ge\vzero,\\
&\lambda_jx_j=0\quad(j=1,\dots,n).
\end{aligned}
$$

Since the constraints are linear, these conditions are necessary at every
optimum ({prf:ref}`thm-06-03-cq-sufficient`), and since an LP is convex they
are also sufficient ({prf:ref}`thm-06-03-convex-kkt`).

**(b)** Put $\vpi=-\vmu$. Stationarity gives $\vlambda=\vc-\mA^\T\vpi$, the
vector of reduced costs $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$. The KKT
conditions become

$$
\begin{aligned}
&\mA\vx=\vb,\ \vx\ge\vzero && \text{(primal feasibility)},\\
&\mA^\T\vpi\le\vc && \text{(dual feasibility)},\\
&x_j\,(c_j-\vpi^\T\mA_{\cdot j})=0,\ j=1,\dots,n && \text{(complementary slackness)},
\end{aligned}
$$

where $\mA^\T\vpi\le\vc$ is the constraint of the dual LP
$\max\{\vb^\T\vpi:\mA^\T\vpi\le\vc\}$. These are exactly the optimality
conditions of {doc}`03-05-lp-duality` ({prf:ref}`thm-03-05-complementary-slackness`): a variable can be positive only if its
reduced cost is zero, and the multipliers of the equality constraints are
the dual variables (simplex multipliers). Multiplying stationarity by
$\vx$ also shows that the duality gap vanishes:
$\vc^\T\vx=\vpi^\T\mA\vx+\vlambda^\T\vx=\vb^\T\vpi+0$.

*Numerical check.* For $\min\,-3x_1-5x_2$ s.t. $x_1+x_3=4$,
$2x_2+x_4=12$, $3x_1+2x_2+x_5=18$, $\vx\ge\vzero$, an LP solver returns
$\vx^*=(2,6,2,0,0)$ and $\vpi=(0,-1.5,-1)$. The reduced costs are
$\vc-\mA^\T\vpi=(0,0,0,1.5,1)\ge\vzero$, every product $x_j\bar c_j$ is
zero, and $\vc^\T\vx^*=-36=\vb^\T\vpi$.
*(An earlier version of this solution names the feasibility conditions,
rather than $\lambda_jx_j=0$, as the complementary slackness condition.)*
:::

:::{prf:example} The multiplier as a sensitivity
:label: ex-08-08-multiplier-sensitivity

Consider $\min f(\vx)=x_1^2+2x_2^2+2x_1x_2-2x_1-4x_2$ subject to
$-x_1-x_2+2\le0$.

(a) Solve the problem with the KKT conditions and find the multiplier
$\lambda^*$.
(b) Relax the constraint to $-x_1-x_2+2\le\delta$ with a small $\delta>0$,
solve again, and show that $\lambda^*$ measures the rate at which the
optimal value changes.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** With $g=-x_1-x_2+2$ and $\nabla g=(-1,-1)^\T$, stationarity reads

$$
\begin{aligned}
2x_1+2x_2-2&=\lambda,\\
2x_1+4x_2-4&=\lambda .
\end{aligned}
$$

The unconstrained stationary point ($\lambda=0$) is $(0,1)$, where
$g=1>0$: infeasible. So the constraint is active, $x_1+x_2=2$, and the
first equation gives $\lambda=2\cdot2-2=2$. The second becomes
$x_1+2x_2=3$, so $x_2=1$, $x_1=1$. Thus

$$
\begin{aligned}
\vx^*&=(1,1),\\
\lambda^*&=2\ge0,\\
f^*&=1+2+2-2-4=-1 .
\end{aligned}
$$

The Hessian $\bigl(\begin{smallmatrix}2&2\\2&4\end{smallmatrix}\bigr)$ has
eigenvalues $0.764$ and $5.236$, so $f$ is strictly convex and $\vx^*$ is
the unique global minimum.

**(b)** The relaxed constraint is $x_1+x_2\ge2-\delta$. For $0\le\delta<1$
the same reasoning gives $x_1+x_2=2-\delta$,
$\lambda(\delta)=2(2-\delta)-2=2-2\delta>0$, $x_1+2x_2=3-\delta$, hence

$$
\begin{aligned}
\vx^*(\delta)&=(1-\delta,\ 1),\\
f^*(\delta)&=(1-\delta)^2+2+2(1-\delta)-2(1-\delta)-4\\
&=\delta^2-2\delta-1 .
\end{aligned}
$$

Then

$$
\frac{\mathrm df^*}{\mathrm d\delta}\Big|_{\delta=0}=-2=-\lambda^*,
$$

in agreement with {eq}`eq-06-03-sensitivity`: relaxing the constraint by
$\delta$ lowers the optimal value by about $\lambda^*\delta$. For example,
$\delta=0.1$ gives $f^*=-1.19$, against the linear prediction
$-1-2(0.1)=-1.2$; the difference quotient
$\bigl(f^*(\delta)-f^*\bigr)/(-\delta)=2-\delta$ tends to $\lambda^*=2$.
At $\delta=1$ the multiplier reaches $0$ and the relaxed constraint is
satisfied by the unconstrained minimum $(0,1)$, so further relaxation no
longer helps. *(An earlier version of this solution writes the active
constraint as $-x_1-x_2=2$; it is $x_1+x_2=2$.)*
:::

:::{prf:example} A KKT point with no active constraint
:label: ex-08-08-interior-kkt

Use the KKT conditions to solve

$$
\begin{aligned}
\min\ & f(\vx)=4(x_1-2)^2+3(x_2-4)^2\\
\st\ & -x_1-x_2+5\le0,\\
& -x_1+1\le0,\\
& -x_2+2\le0,
\end{aligned}
$$

and show that the solution is unique.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The KKT conditions are

$$
\begin{aligned}
8(x_1-2)-\lambda_1-\lambda_2&=0,\\
6(x_2-4)-\lambda_1-\lambda_3&=0,
\end{aligned}
$$

with $\lambda_1(-x_1-x_2+5)=0$, $\lambda_2(1-x_1)=0$, $\lambda_3(2-x_2)=0$,
feasibility and $\lambda_i\ge0$.

Try the case with no active constraint, $\vlambda=\vzero$: then
$x_1=2$, $x_2=4$. Feasibility: $-2-4+5=-1\le0$, $1-2=-1\le0$,
$2-4=-2\le0$. All conditions hold, so $\vx^*=(2,4)$ with
$\vlambda^*=\vzero$ and $f^*=0$. The constraints are all inactive: the
unconstrained minimum happens to be feasible.

Uniqueness: $\nabla^2f=\bigl(\begin{smallmatrix}8&0\\0&6\end{smallmatrix}\bigr)\succ\vzero$,
so $f$ is strictly convex, the feasible set is convex, and a strictly convex
function has at most one minimum over a convex set
({prf:ref}`thm-06-01-convex-local-global`). *(An earlier version of this
solution gives the Hessian as $\diag(12,6)$; the coefficient
$4(x_1-2)^2$ gives $8$.)*
:::

## References

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
