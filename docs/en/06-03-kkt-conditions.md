# Karush–Kuhn–Tucker Conditions

## Overview

{doc}`06-02-optimality-conditions` characterized local minima through
directions: at a minimum no feasible direction is a descent direction. To
use this criterion one must know the feasible directions, and for a set
described by nonlinear inequalities and equations they are not obvious.
This chapter replaces directions by the constraint functions. Linearizing
the constraints that are active at $\vx^*$ gives a cone of directions
$L(\vx^*)$ that is easy to write down; a **constraint qualification**
guarantees that it describes the true local geometry of $\Omega$. Under
that assumption, Farkas' lemma (equivalently, linear programming (LP) duality) turns "no
linearized feasible descent direction" into the existence of multipliers
$\lambda_i\ge0$ and $\mu_j$ with
$\nabla f(\vx^*)+\sum_i\lambda_i\nabla g_i(\vx^*)+\sum_j\mu_j\nabla h_j(\vx^*)=\vzero$:
the **Karush–Kuhn–Tucker (KKT) conditions** ({prf:ref}`thm-06-03-kkt`). The
chapter then interprets the multipliers geometrically and through the
Lagrangian, fixes the sign conventions for maximization and "$\ge$"
constraints, solves a KKT system by active-set case analysis
({prf:ref}`ex-06-03-kkt-example`), states the second-order conditions on the
critical cone, and shows that for convex programs a KKT point is a global
minimum ({prf:ref}`thm-06-03-convex-kkt`). Besides the two preceding
chapters, it uses Farkas' lemma and LP duality from {doc}`03-05-lp-duality`,
convex functions from {doc}`02-02-convex-functions`, and meets the
Lagrangian of {doc}`02-03-convex-optimization-duality` again from a
different side.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- determine the active set $I(\vx)$ and the cone of linearized feasible
  directions at a feasible point;
- explain why the linearized cone can differ from the true local geometry,
  define the tangent cone, and check the common constraint qualifications:
  linear constraints, the linear independence constraint qualification
  (LICQ) and Slater's condition;
- derive the KKT conditions from Farkas' lemma and write them with the
  Lagrangian $L(\vx,\vlambda,\vmu)$, for minimization and maximization
  problems with "$\le$", "$\ge$" and "$=$" constraints;
- solve a small KKT system by case analysis over active sets and verify
  the solution;
- apply the second-order necessary and sufficient conditions on the
  critical cone; and
- use convexity to conclude that a KKT point is a global minimum.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $I(\vx)$ | active set $\{i:g_i(\vx)=0\}$ of the inequality constraints at $\vx$ |
| $L(\vx^*)$ | cone of linearized feasible directions at $\vx^*$ ({eq}`eq-06-03-linearised-cone`) |
| $T(\vx^*)$ | tangent cone of $\Omega$ at $\vx^*$ ({prf:ref}`def-06-03-tangent-cone`) |
| $\lambda_i$, $\vlambda$ | multiplier of $g_i(\vx)\le0$; $\vlambda\in\R^m$ |
| $\mu_j$, $\vmu$ | multiplier of $h_j(\vx)=0$; $\vmu\in\R^\ell$ |
| $L(\vx,\vlambda,\vmu)$ | Lagrangian $f(\vx)+\vlambda^\T\vg(\vx)+\vmu^\T\vh(\vx)$ |
| $\nabla_{\vx}L$, $\nabla^2_{\vx\vx}L$ | gradient and Hessian of the Lagrangian with respect to $\vx$ |
| $C(\vx^*,\vlambda^*)$ | critical cone ({prf:ref}`def-06-03-critical-cone`) |

The letter $L$ is used both for the cone $L(\vx^*)$, which has one
argument, and for the Lagrangian $L(\vx,\vlambda,\vmu)$, which has three;
the argument list tells them apart. Throughout, the problem is
{eq}`eq-06-01-nlp`,

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)\\
\st\ & g_i(\vx)\le0, && i=1,\dots,m,\\
& h_j(\vx)=0, && j=1,\dots,\ell ,
\end{aligned}
$$

with $f$, $g_i$, $h_j$ continuously differentiable.

(sec-06-03-active)=
## Active constraints and linearized directions

At a feasible point $\vx$, an inequality constraint is **active** (or
binding) if $g_i(\vx)=0$ and **inactive** if $g_i(\vx)<0$. The **active
set** is

$$
I(\vx)=\{i\in\{1,\dots,m\}:\ g_i(\vx)=0\}.
$$ (eq-06-03-active-set)

Equality constraints hold with equality at every feasible point and count
as active. An inactive constraint has slack, and by continuity it stays
satisfied in a neighborhood of $\vx$; it does not restrict short moves.
Only active constraints shape the feasible set near $\vx$.

:::{prf:example} Active constraints of a triangle
:label: ex-06-03-triangle-active

Let $\Omega=\{\vx\in\R^2:x_1+x_2\le2,\ x_1\ge0,\ x_2\ge0\}$, written as
$g_1(\vx)=x_1+x_2-2\le0$, $g_2(\vx)=-x_1\le0$, $g_3(\vx)=-x_2\le0$, with
gradients $\nabla g_1=(1,1)^\T$, $\nabla g_2=(-1,0)^\T$,
$\nabla g_3=(0,-1)^\T$. Find the active sets and the feasible directions at
the vertices $(2,0)$ and $(0,0)$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Because the constraints are linear, $g_i(\vx^*+\alpha\vd)=g_i(\vx^*)+\alpha\nabla g_i^\T\vd$
exactly, and a direction $\vd$ is feasible if this is $\le0$ for all small
$\alpha>0$ and all $i$.

```{list-table} Constraints $g_1=x_1+x_2-2$, $g_2=-x_1$, $g_3=-x_2$ of the triangle at the vertices $(2,0)$ and $(0,0)$: value of $g_i$ (active when $0$), and the condition each puts on a feasible direction ("on $\vd$"), namely $\nabla g_i^\T\vd\le0$ for every active $g_i$.
:name: tab-06-03-triangle-active
:header-rows: 1

* - $g_i$
  - $g_i(2,0)$
  - on $\vd$
  - $g_i(0,0)$
  - on $\vd$
* - $g_1$
  - $0$
  - $d_1+d_2\le0$
  - $-2$
  - none
* - $g_2$
  - $-2$
  - none
  - $0$
  - $-d_1\le0$
* - $g_3$
  - $0$
  - $-d_2\le0$
  - $0$
  - $-d_2\le0$
```

So $I(2,0)=\{1,3\}$ and the feasible directions there are
$\{\vd:d_1+d_2\le0,\ d_2\ge0\}$; $I(0,0)=\{2,3\}$ and the feasible
directions are $\{\vd:\vd\ge\vzero\}$.
:::

As the solution of {prf:ref}`ex-06-03-triangle-active` shows, at both
vertices a direction is feasible exactly when it makes an angle of at least
$90^\circ$ with the gradient of every active constraint. The gradients $\nabla g_i$ point out of
$\Omega$, in the direction in which $g_i$ grows towards violation.

For nonlinear constraints the angle condition is still necessary.

:::{prf:proposition} Feasible directions satisfy the linearized constraints
:label: prop-06-03-linearised

If $\vd$ is a feasible direction of $\Omega$ at $\vx^*\in\Omega$, then

$$
\begin{aligned}
\nabla g_i(\vx^*)^\T\vd&\le0, && i\in I(\vx^*),\\
\nabla h_j(\vx^*)^\T\vd&=0, && j=1,\dots,\ell.
\end{aligned}
$$
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

For $i\in I(\vx^*)$ and $0<\alpha\le\bar\alpha$,
$g_i(\vx^*+\alpha\vd)\le0=g_i(\vx^*)$, so the difference quotient
$\bigl(g_i(\vx^*+\alpha\vd)-g_i(\vx^*)\bigr)/\alpha$ is $\le0$; its limit as
$\alpha\downarrow0$ is the directional derivative
$\nabla g_i(\vx^*)^\T\vd$, which is therefore $\le0$. For $h_j$ the
difference quotient is identically $0$.
:::

The directions allowed by the proposition form the **cone of linearized
feasible directions**

$$
\begin{aligned}
L(\vx^*)=\bigl\{\vd\in\R^n:\ &\nabla g_i(\vx^*)^\T\vd\le0\ \ (i\in I(\vx^*)),\\
&\nabla h_j(\vx^*)^\T\vd=0\ \ (j=1,\dots,\ell)\bigr\}.
\end{aligned}
$$ (eq-06-03-linearised-cone)

It is a polyhedral cone, defined by linear inequalities in $\vd$ that can be
written down from the gradients alone. The proposition says that every
feasible direction lies in $L(\vx^*)$. The question is whether, conversely,
$L(\vx^*)$ is a faithful description of the directions in which one can
move.

(sec-06-03-cq)=
## Constraint qualifications

The next example shows that the answer can be no, even for a very simple
set.

:::{prf:example} A cubic description of a triangle
:label: ex-06-03-cubic

Let

$$
\begin{aligned}
g_1(\vx)&=(x_1+x_2-1)^3\le0,\\
g_2(\vx)&=-x_1\le0,\\
g_3(\vx)&=-x_2\le0,
\end{aligned}
$$

so that $\Omega=\{\vx\in\R^2:g_1(\vx)\le0,\ g_2(\vx)\le0,\ g_3(\vx)\le0\}$ is
the triangle with vertices $(0,0)$, $(1,0)$, $(0,1)$, and
$\vx^*=(0.5,0.5)$. (a) Find $I(\vx^*)$ and $L(\vx^*)$. (b) Is $\vd=(2,1)$ a
feasible direction? (c) Minimize $f(\vx)=(x_1-1)^2+(x_2-1)^2$ over
$\Omega$ and check whether multipliers as in {eq}`eq-06-03-kkt-stationarity`
exist at the minimum.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** $g_1(\vx^*)=0$, $g_2(\vx^*)=g_3(\vx^*)=-0.5<0$, so
$I(\vx^*)=\{1\}$. The gradient
$\nabla g_1(\vx)=3(x_1+x_2-1)^2(1,1)^\T$ vanishes at $\vx^*$, so the only
condition in {eq}`eq-06-03-linearised-cone` reads $0\le0$, and
$L(\vx^*)=\R^2$: every direction passes the linearized test.

**(b)** Along $\vd=(2,1)$,
$x_1+x_2-1=0.5+2\alpha+0.5+\alpha-1=3\alpha$, so
$g_1(\vx^*+\alpha\vd)=(3\alpha)^3=27\alpha^3>0$ for every $\alpha>0$. The
direction leaves $\Omega$ immediately: it is not feasible, although it lies
in $L(\vx^*)$. The linearization of $g_1$ is the zero function, which
carries no information about the side of the line $x_1+x_2=1$ on which
$\Omega$ lies. The true feasible directions are $\{\vd:d_1+d_2\le0\}$.

**(c)** The minimum is the point of the triangle closest to $(1,1)$, the
foot of the perpendicular on the edge $x_1+x_2=1$, which is $\vx^*=(0.5,0.5)$
with $f^*=0.5$ (a solver confirms this). There
$\nabla f(\vx^*)=(-1,-1)^\T$. Since $g_2$ and $g_3$ are inactive, their
multipliers must vanish, and the stationarity condition
$\nabla f(\vx^*)+\lambda_1\nabla g_1(\vx^*)=\vzero$ reduces to
$(-1,-1)^\T+\lambda_1(0,0)^\T=\vzero$, which has no solution. The minimum
is not a KKT point.
:::

If the triangle of {prf:ref}`ex-06-03-cubic` is described by the linear
constraint $x_1+x_2-1\le0$, then $\nabla g_1=(1,1)^\T$, $L(\vx^*)=\{d_1+d_2\le0\}$ is
exactly the set of feasible directions, and $\lambda_1=1$ satisfies
$(-1,-1)^\T+1\cdot(1,1)^\T=\vzero$. The set is the same; only its
algebraic description differs. The linearized cone depends on the
description, while the true local geometry depends only on the set.
{numref}`fig-06-03-cubic-cq` shows both effects.

```{figure} figures/06-03-cubic-cq.svg
:name: fig-06-03-cubic-cq
:width: 680px
:class: dark-light
:alt: Left, a triangle with the point (0.5,0.5) on its hypotenuse and an arrow d = (2,1) leaving the triangle; right, the unit disc with the point (1,0), its outward gradient (2,0) and a vertical tangent arrow (0,1) that leaves the disc while the circle bends away

Two ways in which the linearized cone can mislead. Left
({prf:ref}`ex-06-03-cubic`): with $g_1=(x_1+x_2-1)^3$ the gradient vanishes at
$\vx^*=(0.5,0.5)$, so every direction, including $\vd=(2,1)$, passes the
linearized test, yet $\vd$ points out of $\Omega$. Right: for the disc
$x_1^2+x_2^2\le1$ at $(1,0)$ the gradient $(2,0)$ is nonzero and the
tangent direction $(0,1)$ passes the test, but every step along it leaves
the disc; only curves that bend with the boundary stay feasible.
```

A second example shows that the trouble is not only with degenerate
gradients. For the unit disc $\Omega=\{\vx:g(\vx)=x_1^2+x_2^2-1\le0\}$ at
$\vx^*=(1,0)$, $\nabla g(\vx^*)=(2,0)^\T\ne\vzero$ and
$L(\vx^*)=\{\vd:d_1\le0\}$. The direction $\vd=(0,1)$ lies in $L(\vx^*)$, but
$g(1,\alpha)=\alpha^2>0$ for all $\alpha>0$, so it is not a feasible
direction: at a point of a strictly curved boundary no tangent direction is
feasible. Nevertheless one can move "along" $(0,1)$ by following the circle,
$\vx(t)=(\cos t,\sin t)$, whose velocity at $t=0$ is $(0,1)$. Feasible
*curves* rather than feasible *straight lines* are the right objects.

:::{prf:definition} Tangent cone
:label: def-06-03-tangent-cone

A vector $\vd$ is a **tangent direction** of $\Omega$ at $\vx^*\in\Omega$ if
there are points $\vx^{(k)}\in\Omega$ and numbers $t_k>0$ with
$\vx^{(k)}\to\vx^*$, $t_k\to0$ and

$$
\frac{\vx^{(k)}-\vx^*}{t_k}\to\vd\qquad(k\to\infty).
$$

The set $T(\vx^*)$ of tangent directions is the **tangent cone** of
$\Omega$ at $\vx^*$.
:::

Every feasible direction is a tangent direction (take
$\vx^{(k)}=\vx^*+t_k\vd$), and so is every velocity vector of a smooth
feasible curve through $\vx^*$. For the disc at $(1,0)$,
$T(\vx^*)=\{d_1\le0\}=L(\vx^*)$; for the cubic description of the triangle,
$T(\vx^*)=\{d_1+d_2\le0\}\ne L(\vx^*)=\R^2$. The tangent cone depends only
on the set $\Omega$.

:::{prf:lemma} Tangent cone and first-order condition
:label: lem-06-03-tangent

Let $\vx^*\in\Omega$.

1. $T(\vx^*)\subseteq L(\vx^*)$.
2. If $\vx^*$ is a local minimum of $f$ over $\Omega$, then
   $\nabla f(\vx^*)^\T\vd\ge0$ for all $\vd\in T(\vx^*)$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\vd\in T(\vx^*)$ with sequences $\vx^{(k)}$, $t_k$ as in the
definition, and write $\vx^{(k)}=\vx^*+t_k\vd^{(k)}$ with
$\vd^{(k)}\to\vd$. For a differentiable function $\phi$,

$$
\phi(\vx^{(k)})=\phi(\vx^*)+t_k\nabla\phi(\vx^*)^\T\vd^{(k)}+o(t_k).
$$

1. For $i\in I(\vx^*)$, take $\phi=g_i$: $g_i(\vx^{(k)})\le0=g_i(\vx^*)$
   gives $t_k\nabla g_i(\vx^*)^\T\vd^{(k)}+o(t_k)\le0$; dividing by $t_k$
   and letting $k\to\infty$ gives $\nabla g_i(\vx^*)^\T\vd\le0$. For $h_j$
   the same argument with equality gives $\nabla h_j(\vx^*)^\T\vd=0$.
2. Take $\phi=f$: for large $k$ the point $\vx^{(k)}$ is close to $\vx^*$,
   so $f(\vx^{(k)})\ge f(\vx^*)$, and the same division gives
   $\nabla f(\vx^*)^\T\vd\ge0$.
:::

Part 2 is the first-order condition of {prf:ref}`thm-06-02-first-order`
extended from feasible directions to tangent directions. To turn it into a
statement about gradients, one needs the inclusion of part 1 to be an
equality.

:::{prf:definition} Constraint qualification
:label: def-06-03-cq

The constraints satisfy **Abadie's constraint qualification** (Abadie, 1967) at
$\vx^*\in\Omega$ if $T(\vx^*)=L(\vx^*)$. More generally, a **constraint
qualification** (CQ) is any condition on the constraint functions at
$\vx^*$ under which every local minimum of every smooth objective at
$\vx^*$ admits KKT multipliers.
:::

*(The constraint qualification is sometimes stated as "every $\vd$ with
$\nabla g_i(\vx^*)^\T\vd\le0$ for the active $i$ is a feasible direction".
That version fails for every strictly curved constraint, as the disc shows,
even when the gradients are linearly independent; the tangent cone is the
correct object.)*

Checking $T(\vx^*)=L(\vx^*)$ directly is rarely convenient. The following
conditions are easy to check and each implies it.

:::{prf:theorem} Sufficient conditions for a constraint qualification
:label: thm-06-03-cq-sufficient

Each of the following implies $T(\vx^*)=L(\vx^*)$.

1. **Linear constraints**: all $g_i$ and $h_j$ are affine functions.
2. **Linear independence constraint qualification (LICQ)**: the gradients
   $\{\nabla g_i(\vx^*):i\in I(\vx^*)\}\cup\{\nabla h_j(\vx^*):j=1,\dots,\ell\}$
   are linearly independent.
3. **Slater's condition** (Slater, 1950): every $g_i$ is convex, every $h_j$ is affine,
   and there is a point $\hat{\vx}$ with $g_i(\hat{\vx})<0$ for all $i$ and
   $\vh(\hat{\vx})=\vzero$. (Strict inequality is needed only for the
   nonlinear $g_i$.)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

We show $L(\vx^*)\subseteq T(\vx^*)$ in each case; the reverse inclusion
is part 1 of {prf:ref}`lem-06-03-tangent`.

1. For affine constraints, $g_i(\vx^*+\alpha\vd)=g_i(\vx^*)+\alpha\nabla g_i^\T\vd$
   and $h_j(\vx^*+\alpha\vd)=\alpha\nabla h_j^\T\vd$ exactly. If
   $\vd\in L(\vx^*)$, the active inequalities and the equalities hold for
   all $\alpha\ge0$, and the inactive ones hold for small $\alpha$, as in
   {prf:ref}`prop-06-02-feasible-directions`. So every $\vd\in L(\vx^*)$ is a
   feasible direction, hence a tangent direction:
   $L(\vx^*)\subseteq T(\vx^*)$, and part 1 of {prf:ref}`lem-06-03-tangent`
   gives equality.
2. Given $\vd\in L(\vx^*)$, the implicit function theorem, applied to the
   system of active constraints made into equations, produces a smooth
   feasible curve with velocity $\vd$; see Nocedal and Wright (2006,
   Lemma 12.2) for the details.
3. Let $\vd^0=\hat{\vx}-\vx^*$. For an active $g_i$, the gradient
   inequality of convex functions ({prf:ref}`thm-02-02-first-order`) gives
   $\nabla g_i(\vx^*)^\T\vd^0\le g_i(\hat{\vx})-g_i(\vx^*)=g_i(\hat{\vx})$,
   which is $<0$ for the nonlinear $g_i$ and $\le0$ for the affine ones;
   and $\nabla h_j(\vx^*)^\T\vd^0=h_j(\hat{\vx})-h_j(\vx^*)=0$ because $h_j$
   is affine. Given $\vd\in L(\vx^*)$ and $\varepsilon>0$, the direction
   $\vd^\varepsilon=\vd+\varepsilon\vd^0$ therefore satisfies
   $\nabla g_i(\vx^*)^\T\vd^\varepsilon<0$ for the active nonlinear $g_i$,
   so these decrease from the value $0$ along it
   ({prf:ref}`lem-06-02-descent`); the active affine $g_i$ and the $h_j$ change
   by $\alpha\nabla g_i(\vx^*)^\T\vd^\varepsilon\le0$ and
   $\alpha\nabla h_j(\vx^*)^\T\vd^\varepsilon=0$; and the inactive $g_i$
   stay negative for small steps. Hence $\vd^\varepsilon$ is a feasible
   direction and lies in $T(\vx^*)$. As $\varepsilon\downarrow0$,
   $\vd^\varepsilon\to\vd$, and since the tangent cone is closed (a
   diagonal choice of sequences shows that a limit of tangent directions is
   a tangent direction), $\vd\in T(\vx^*)$. See also Bazaraa et al. (2006,
   Section 5.2).
:::

A point at which LICQ holds is called a **regular point**. LICQ is the
condition used most often in practice; it fails in {prf:ref}`ex-06-03-cubic`
because $\nabla g_1(\vx^*)=\vzero$, and it fails whenever more than $n$
constraints are active. Linear constraints need no check at all, which
covers LP, quadratic programs (QP) and problems such as {prf:ref}`ex-06-03-kkt-example`.
Slater's condition is a property of the whole problem rather than of a
point, and is the standard assumption for convex programs
({prf:ref}`thm-02-03-slater`). Weaker conditions exist, such as
the Mangasarian–Fromovitz condition (Mangasarian and Fromovitz, 1967:
there is a $\vd$ with
$\nabla g_i(\vx^*)^\T\vd<0$ for $i\in I(\vx^*)$ and
$\nabla h_j(\vx^*)^\T\vd=0$, with the $\nabla h_j(\vx^*)$ linearly
independent); see Bazaraa et al. (2006, Chapter 5).

(sec-06-03-kkt)=
## The KKT conditions

Suppose $\vx^*$ is a local minimum and a constraint qualification holds,
so that $T(\vx^*)=L(\vx^*)$. By {prf:ref}`lem-06-03-tangent`,

$$
\nabla f(\vx^*)^\T\vd\ge0\qquad\text{for all }\vd\in L(\vx^*).
$$ (eq-06-03-no-linear-descent)

Everything in {eq}`eq-06-03-no-linear-descent` is linear in $\vd$. It says that
the linear program

$$
\begin{aligned}
\min_{\vd\in\R^n}\ & \nabla f(\vx^*)^\T\vd\\
\st\ & \nabla g_i(\vx^*)^\T\vd\le0, && i\in I(\vx^*),\\
& \nabla h_j(\vx^*)^\T\vd=0, && j=1,\dots,\ell,
\end{aligned}
$$ (eq-06-03-direction-lp)

has optimal value $0$, attained at $\vd^*=\vzero$. (The LP is homogeneous:
$\vd=\vzero$ is feasible with value $0$, and if some feasible $\vd$ had a
negative value, multiplying it by $t\to\infty$ would make the value tend to
$-\infty$. Hence its optimal value is either $0$ or $-\infty$, and
{eq}`eq-06-03-no-linear-descent` says it is $0$.) The dual of an LP with a
finite optimum is feasible ({prf:ref}`thm-03-05-strong-duality`). With a multiplier
$\lambda_i\ge0$ for each "$\le$" row and a free multiplier $\mu_j$ for each
equality row, the dual of {eq}`eq-06-03-direction-lp` is

$$
\begin{aligned}
\max\ & 0\\
\st\ & \sum_{i\in I(\vx^*)}\lambda_i\nabla g_i(\vx^*)+\sum_{j=1}^{\ell}\mu_j\nabla h_j(\vx^*)\\
& \qquad=-\nabla f(\vx^*),\\
& \lambda_i\ge0, && i\in I(\vx^*),\\
& \mu_j\ \text{free}, && j=1,\dots,\ell.
\end{aligned}
$$

(The free primal variables $d_1,\dots,d_n$ give equality constraints in the
dual, one per coordinate.) Its feasibility is exactly the existence of
multipliers. The same conclusion follows directly from a variant of Farkas' lemma
(Farkas, 1902; {prf:ref}`thm-03-05-farkas`, with $\vmu=\vmu^+-\vmu^-$ split into
nonnegative parts): of the two systems

$$
\begin{aligned}
\text{(I)}\quad&\nabla f(\vx^*)^\T\vd<0,\\
&\vd\in L(\vx^*)
\end{aligned}
$$

in the unknown $\vd$, and

$$
\begin{aligned}
\text{(II)}\quad&\nabla f(\vx^*)+\sum_{i\in I(\vx^*)}\lambda_i\nabla g_i(\vx^*)\\
&\qquad+\sum_{j=1}^{\ell}\mu_j\nabla h_j(\vx^*)=\vzero,\\
&\lambda_i\ge0, \qquad i\in I(\vx^*)
\end{aligned}
$$

in the unknowns $\lambda_i$, $\mu_j$,
exactly one has a solution, and {eq}`eq-06-03-no-linear-descent` says that (I)
has none. Finally, setting $\lambda_i=0$ for the inactive constraints lets
the sum run over all $i$, at the price of the condition
$\lambda_ig_i(\vx^*)=0$.

:::{prf:theorem} Karush–Kuhn–Tucker necessary conditions
:label: thm-06-03-kkt

Let $\vx^*$ be a local minimum of {eq}`eq-06-01-nlp` at which a constraint
qualification holds (for example linear constraints, LICQ or Slater's
condition). Then there exist multipliers $\vlambda^*\in\R^m$ and
$\vmu^*\in\R^\ell$ such that

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

Under LICQ the multipliers are unique.
:::

The four groups are called **stationarity** (or the KKT equations),
**primal feasibility**, **dual feasibility** (nonnegativity of the
inequality multipliers) and **complementary slackness**. A point $\vx^*$
for which multipliers satisfying {eq}`eq-06-03-kkt-stationarity` to
{eq}`eq-06-03-kkt-cs` exist is a **KKT point**. Uniqueness under LICQ holds
because two multiplier vectors would give two representations of
$-\nabla f(\vx^*)$ by the linearly independent active gradients.

**Equality multipliers are free.** An equality $h_j=0$ can also be written
as the two inequalities $h_j\le0$ and $-h_j\le0$, with multipliers
$\mu_j^+\ge0$ and $\mu_j^-\ge0$. They contribute
$\mu_j^+\nabla h_j-\mu_j^-\nabla h_j=(\mu_j^+-\mu_j^-)\nabla h_j$ to the
stationarity condition, so only the difference $\mu_j=\mu_j^+-\mu_j^-$
matters, and it can have either sign. This is why the conditions put no
sign restriction on $\mu_j$. (The derivation above treats equalities
directly, because the split pair has dependent gradients $\pm\nabla h_j$
and never satisfies LICQ.)

**The theorem is necessary, not sufficient.** A KKT point can be a local
maximum or a saddle point; second-order conditions or convexity are needed
to conclude a minimum ({ref}`sec-06-03-second-order`, {ref}`sec-06-03-convex`).
And without a constraint qualification a minimum need not be a KKT point at
all, as in {prf:ref}`ex-06-03-cubic` and the classical cusp of
{prf:ref}`ex-08-08-cusp`.

:::{prf:remark} History and the Fritz John conditions
:label: rem-06-03-history

The conditions were derived by Karush (1939) in an unpublished master's
thesis and independently by Kuhn and Tucker (1951), whose paper made them
widely known; for equality constraints alone they reduce to Lagrange's
multiplier rule. Without any constraint qualification, a local minimum
always satisfies the weaker **Fritz John conditions** (John, 1948): there are
$\lambda_0\ge0$, $\vlambda\ge\vzero$, $\vmu$, not all zero, with
$\lambda_0\nabla f+\sum_i\lambda_i\nabla g_i+\sum_j\mu_j\nabla h_j=\vzero$
and complementary slackness. In the case $\lambda_0=0$ the conditions say
nothing about $f$; when $\lambda_0>0$, dividing by $\lambda_0$ gives KKT
multipliers. A constraint qualification guarantees that multipliers with
$\lambda_0>0$ exist, and the Mangasarian–Fromovitz condition is precisely
the condition under which *every* Fritz John multiplier vector has
$\lambda_0>0$ (Bazaraa et al., 2006, Sections 4.2–4.3).
:::

### The Lagrangian

The **Lagrangian** of {eq}`eq-06-01-nlp` is

$$
\begin{aligned}
L(\vx,\vlambda,\vmu)&=f(\vx)+\sum_{i=1}^m\lambda_ig_i(\vx)+\sum_{j=1}^{\ell}\mu_jh_j(\vx)\\
&=f(\vx)+\vlambda^\T\vg(\vx)+\vmu^\T\vh(\vx).
\end{aligned}
$$ (eq-06-03-lagrangian)

Its partial derivatives reproduce the KKT system compactly:

$$
\begin{aligned}
&\nabla_{\vx}L(\vx^*,\vlambda^*,\vmu^*)=\vzero,\\
&\frac{\partial L}{\partial\mu_j}=h_j(\vx^*)=0,\\
&\frac{\partial L}{\partial\lambda_i}=g_i(\vx^*)\le0,\\
&\lambda_i^*\ge0,\\
&\lambda_i^*\,\frac{\partial L}{\partial\lambda_i}=0 .
\end{aligned}
$$

Stationarity says that $\vx^*$ is a stationary point of the Lagrangian in
$\vx$: the constraints are "priced" into the objective with prices
$\lambda_i$, $\mu_j$, and at those prices the constrained minimum looks like
an unconstrained stationary point. For convex problems this is made precise
by Lagrangian duality ({prf:ref}`thm-02-03-kkt-sufficient`), in which the
KKT multipliers are optimal solutions of the dual problem.

### Geometry

With only the active inequalities present, stationarity reads

$$
\begin{aligned}
-\nabla f(\vx^*)&=\sum_{i\in I(\vx^*)}\lambda_i^*\nabla g_i(\vx^*),\\
\lambda_i^*&\ge0, \qquad i\in I(\vx^*).
\end{aligned}
$$ (eq-06-03-kkt-cone)

The steepest descent direction $-\nabla f(\vx^*)$ lies in the cone generated
by the outward normals $\nabla g_i(\vx^*)$ of the active constraints. If it
pointed anywhere else, some direction would decrease $f$ while keeping the
linearized active constraints satisfied, and Farkas' lemma says that
precisely this alternative is excluded. One can picture a ball rolling
downhill on the contours of $f$ and coming to rest against the walls of
$\Omega$: the "force" $-\nabla f$ is balanced by nonnegative reaction forces
$\lambda_i\nabla g_i$ from the walls it touches, and walls it does not touch
exert no force ($\lambda_i=0$ for inactive $i$).
{numref}`fig-06-03-kkt-geometry` shows a case with two curved active
constraints.

```{figure} figures/06-03-kkt-geometry.svg
:name: fig-06-03-kkt-geometry
:width: 560px
:class: dark-light
:alt: A lens-shaped feasible region between two curves meeting at the origin, a third inactive line on the left, elliptical contours of f centered to the right, the gradients of the two active constraints at the origin and the negative gradient of f pointing to the right between them

Geometry of the KKT conditions. The feasible set lies between the curves
$g_1(\vx)=x_2+0.7x_1+0.3x_1^2=0$ and $g_2(\vx)=-x_2+0.7x_1+0.3x_1^2=0$ and to
the right of the line $g_3=0$; the objective is
$f(\vx)=(x_1-1.5)^2+2x_2^2$. At $\vx^*=(0,0)$, $-\nabla f=(3,0)$ lies in the
cone spanned by $\nabla g_1=(0.7,1)$ and $\nabla g_2=(0.7,-1)$:
$-\nabla f=\tfrac{15}{7}\nabla g_1+\tfrac{15}{7}\nabla g_2$. The constraint
$g_3$ is inactive and has multiplier $0$.
```

In {numref}`fig-06-03-kkt-geometry` the minimum of $f$ over the lens is at
$\vx^*=(0,0)$ with $f^*=2.25$, and $\nabla f(\vx^*)=(-3,0)^\T$. Solving
$(3,0)=\lambda_1(0.7,1)+\lambda_2(0.7,-1)$ gives $\lambda_1=\lambda_2$ from
the second component and $1.4\lambda_1=3$ from the first, so
$\lambda_1=\lambda_2=15/7\approx2.143$, both positive as the theorem
requires. Every direction $\vd$ into the lens makes an angle of at least
$90^\circ$ with both $\nabla g_1$ and $\nabla g_2$, hence with
$-\nabla f$, which is a nonnegative combination of them; so
$\nabla f(\vx^*)^\T\vd\ge0$.

The same picture explains the end of {doc}`06-02-optimality-conditions`, where
$\vx^*=(0.5,1.5)$ minimized $(x_1-1)^2+(x_2-2)^2$ over the triangle
$x_1+x_2\le2$, $\vx\ge\vzero$: there $-\nabla f(\vx^*)=(1,1)^\T=1\cdot\nabla g_1$
with $\lambda_1=1$, and the sign constraints are inactive.

### Complementary slackness and the meaning of the multipliers

Complementary slackness {eq}`eq-06-03-kkt-cs` says that for each inequality at
least one of $\lambda_i^*$ and $-g_i(\vx^*)$ is zero: a constraint with slack
has a zero multiplier, and a positive multiplier forces the constraint to be
active. Both may vanish together (a **degenerate** constraint, active with
zero multiplier); when this never happens, **strict complementarity**
holds. The same structure appears in LP as complementary slackness between
primal variables and dual slacks ({prf:ref}`thm-03-05-complementary-slackness`), and indeed the KKT
conditions of an LP are exactly primal feasibility, dual feasibility and
complementary slackness ({prf:ref}`ex-08-08-kkt-lp`).

The multipliers also measure how the optimal value reacts to the
constraints. If the $i$-th constraint is relaxed to $g_i(\vx)\le b_i$ with a
small $b_i>0$, then under LICQ, the second-order sufficient condition of
{prf:ref}`thm-06-03-second-order`
and strict complementarity, the perturbed problem has a local minimum
$\vx^*(b_i)$ near $\vx^*$ whose value $f^*(b_i)=f(\vx^*(b_i))$ is
differentiable, and

$$
\frac{\mathrm d f^*}{\mathrm d b_i}\Big|_{b_i=0}=-\lambda_i^* .
$$ (eq-06-03-sensitivity)

A large multiplier marks a constraint whose relaxation pays off; a zero
multiplier marks one whose small changes do not matter to first order. This
is the nonlinear counterpart of the shadow prices of LP
({prf:ref}`thm-03-05-shadow-price`, {doc}`03-07-sensitivity-analysis`); {prf:ref}`ex-08-08-multiplier-sensitivity`
checks {eq}`eq-06-03-sensitivity` on an example (Luenberger and Ye, 2016,
Section 11.7).

(sec-06-03-signs)=
## Sign conventions

{prf:ref}`thm-06-03-kkt` is stated for a minimization problem with
constraints $g_i\le0$. Other forms are reduced to it by changing signs, and
the result is worth remembering because mistakes here are common.

```{list-table} KKT conditions for the common problem forms. In every row $\lambda_i\ge0$, $\mu_j$ is free, $h_j(\vx)=0$, and complementary slackness is $\lambda_ig_i(\vx)=0$.
:name: tab-06-03-signs
:header-rows: 1

* - Problem
  - Constraints
  - Stationarity
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

For example, $\max f$ subject to $g_i\le0$ is $\min(-f)$ subject to
$g_i\le0$, with stationarity $-\nabla f+\sum_i\lambda_i\nabla g_i=\vzero$:
at a constrained maximum, $+\nabla f$ (the direction of fastest increase)
lies in the cone of the outward normals. The sign of $\mu_j$ is immaterial
because $\mu_j$ is free. A simple rule: at an optimum, the gradient of the
objective must point "against" the active constraints in the direction one
would like to move, and the multipliers measure how hard it pushes.

Some texts write constraints as $c_i(\vx)\ge0$ and the Lagrangian as
$f-\sum_i\lambda_ic_i$ (Nocedal and Wright, 2006), which is the second row
of the table. All conventions give the same points and the same multiplier
values up to sign; what matters is to use one consistently.

**Sign restrictions on variables.** For $\min f$ subject to $g_i\le0$ and
$\vx\ge\vzero$, the multipliers of the bounds $-x_j\le0$ can be eliminated:
stationarity gives them as the components of
$\nabla f+\sum_i\lambda_i\nabla g_i$, and their sign and complementary
slackness conditions become

$$
\begin{aligned}
&\nabla f(\vx^*)+\sum_i\lambda_i\nabla g_i(\vx^*)\ge\vzero,\\
&\vx^*\ge\vzero,\\
&(\vx^*)^\T\Bigl(\nabla f(\vx^*)+\sum_i\lambda_i\nabla g_i(\vx^*)\Bigr)=0,
\end{aligned}
$$

which generalizes {eq}`eq-06-02-orthant`.

(sec-06-03-example)=
## Solving the KKT conditions by active sets

The KKT conditions are a system of equations and inequalities, and the
complementary slackness conditions make it combinatorial: for each $i$,
either $\lambda_i=0$ or $g_i(\vx)=0$. For small problems one can split into
cases according to which constraints are active.

:::{prf:algorithm} Active-set case analysis of the KKT conditions
:label: alg-06-03-active-set

**Input** $f$, $g_1,\dots,g_m$, $h_1,\dots,h_\ell$ (differentiable).

**Output** all KKT points with their multipliers.

1. **Unconstrained candidate.** Solve $\nabla f(\vx)+\sum_j\mu_j\nabla h_j(\vx)=\vzero$,
   $\vh(\vx)=\vzero$ (all $\lambda_i=0$). Keep solutions with
   $\vg(\vx)\le\vzero$.
2. **Choose an active set.** For a set $S\subseteq\{1,\dots,m\}$, impose
   $g_i(\vx)=0$ for $i\in S$ and $\lambda_i=0$ for $i\notin S$.
3. **Solve.** Solve stationarity {eq}`eq-06-03-kkt-stationarity` together with
   $g_i(\vx)=0$ ($i\in S$) and $\vh(\vx)=\vzero$ for $\vx$, $\lambda_i$
   ($i\in S$) and $\vmu$.
4. **Check.** Accept a solution if $g_i(\vx)\le0$ for $i\notin S$ and
   $\lambda_i\ge0$ for $i\in S$; otherwise discard it.
5. **Repeat** steps 2–4 for the sets $S$ not excluded by reasoning; among
   the KKT points found, the global minimum (if it is a KKT point) has the
   smallest $f$.
:::

There are $2^m$ active sets, so in practice the cases are pruned by sign
arguments, by the observation that the unconstrained minimum violates some
constraints, or by geometry. For a convex problem one may stop at the first
KKT point found ({ref}`sec-06-03-convex`). The following example shows both
the pruned reasoning and the full enumeration.

:::{prf:example} Solving a KKT system with five constraints
:label: ex-06-03-kkt-example

Find the minimum of

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

using the KKT conditions.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in


**KKT system.** The gradients are $\nabla f=(2x_1,2x_2,2x_3)^\T$,
$\nabla g_1=(2,1,0)^\T$, $\nabla g_2=(1,0,1)^\T$, $\nabla g_3=(-1,0,0)^\T$,
$\nabla g_4=(0,-1,0)^\T$, $\nabla g_5=(0,0,-1)^\T$. Stationarity
{eq}`eq-06-03-kkt-stationarity`, component by component:

$$
\begin{aligned}
2x_1+2\lambda_1+\lambda_2-\lambda_3&=0,\\
2x_2+\lambda_1-\lambda_4&=0,\\
2x_3+\lambda_2-\lambda_5&=0;
\end{aligned}
$$ (eq-06-03-example-stationarity)

complementary slackness:

$$
\begin{aligned}
\lambda_1(2x_1+x_2-5)&=0,\\
\lambda_2(x_1+x_3-2)&=0,\\
\lambda_3(1-x_1)&=0,\\
\lambda_4(2-x_2)&=0,\\
\lambda_5x_3&=0;
\end{aligned}
$$

plus $g_i(\vx)\le0$ and $\lambda_i\ge0$ for $i=1,\dots,5$.
*(The second complementary slackness condition involves $g_2=x_1+x_3-2$; it
is sometimes misprinted as $\lambda_2(x_1+x_2-2)=0$.)*

**Case analysis.** The unconstrained minimum $\vx=\vzero$ violates
$x_1\ge1$ and $x_2\ge2$, so some constraints are active. Rather than trying
all $2^5=32$ active sets, use signs.

1. *Constraint 4 is active.* The second equation of
   {eq}`eq-06-03-example-stationarity` gives $\lambda_4=2x_2+\lambda_1$.
   Feasibility requires $x_2\ge2$, and $\lambda_1\ge0$, so
   $\lambda_4\ge4>0$. By complementary slackness $g_4=0$: $x_2=2$, and then
   $\lambda_4=4+\lambda_1$.
2. *Constraint 3 is active.* The first equation gives
   $\lambda_3=2x_1+2\lambda_1+\lambda_2\ge2\cdot1+0+0=2>0$, using $x_1\ge1$.
   Hence $g_3=0$: $x_1=1$, and $\lambda_3=2+2\lambda_1+\lambda_2$.
3. *Constraint 1 is inactive.* With $x_1=1$, $x_2=2$:
   $g_1=2+2-5=-1<0$, so $\lambda_1=0$. Consequently $\lambda_4=4$ and
   $\lambda_3=2+\lambda_2$.
4. *Constraints 2 and 5.* With $x_1=1$, $g_2=x_3-1$ and the third equation
   gives $\lambda_5=2x_3+\lambda_2$. If $\lambda_2>0$, complementary
   slackness forces $x_3=1$, then $\lambda_5=2+\lambda_2>0$ forces
   $x_3=0$: a contradiction. So $\lambda_2=0$ and $\lambda_5=2x_3$; now
   $\lambda_5x_3=2x_3^2=0$ gives $x_3=0$ and $\lambda_5=0$.

The KKT system therefore has exactly one solution:

$$
\begin{aligned}
\vx^*&=(1,2,0),\\
\vlambda^*&=(0,0,2,4,0),\\
f(\vx^*)&=1+4+0=5 .
\end{aligned}
$$

**Full enumeration.** The same result is obtained by running
{prf:ref}`alg-06-03-active-set` over all 32 subsets $S$. Since all constraints
are linear and $f$ is quadratic, step 3 is a linear system in each case.
{numref}`tab-06-03-active-set-cases` lists representative cases.

```{list-table} Selected cases of the active-set analysis of the five-constraint example. "Infeasible" means that a constraint outside $S$ is violated; "inconsistent" that the equations $g_i=0$, $i\in S$, have no common solution.
:name: tab-06-03-active-set-cases
:header-rows: 1

* - Active set $S$
  - Solution of step 3
  - Outcome
* - $\varnothing$
  - $\vx=(0,0,0)$
  - infeasible ($g_3,g_4>0$)
* - $\{3\}$
  - $\vx=(1,0,0)$, $\lambda_3=2$
  - infeasible ($g_4>0$)
* - $\{4\}$
  - $\vx=(0,2,0)$, $\lambda_4=4$
  - infeasible ($g_3>0$)
* - $\{1,4\}$
  - $\vx=(1.5,2,0)$, $\lambda_1=-1.5$, $\lambda_4=2.5$
  - $\lambda_1<0$
* - $\{2,4\}$
  - $\vx=(1,2,1)$, $\lambda_2=-2$, $\lambda_4=4$
  - $\lambda_2<0$
* - $\{3,4\}$
  - $\vx=(1,2,0)$, $\lambda_3=2$, $\lambda_4=4$
  - **KKT point**
* - $\{3,4,5\}$
  - $\vx=(1,2,0)$, $\lambda_3=2$, $\lambda_4=4$, $\lambda_5=0$
  - **KKT point** (the same one)
* - $\{1,3,4\}$
  - $x_1=1$, $x_2=2$ but $2x_1+x_2=4\ne5$
  - inconsistent
```

Of the 32 sets, 14 give infeasible points, 8 give a negative multiplier,
8 are inconsistent, and 2 give the KKT point $(1,2,0)$. The two
accepted sets differ only in whether the degenerate constraint $g_5$
(active with $\lambda_5=0$) is counted as active.

**Verification.** At $\vx^*$:
$\vg(\vx^*)=(-1,-1,0,0,0)\le\vzero$; stationarity
$(2,4,0)+2(-1,0,0)+4(0,-1,0)=(0,0,0)$; complementary slackness holds term by
term. The active set is $I(\vx^*)=\{3,4,5\}$ with gradients
$-\ve_1,-\ve_2,-\ve_3$, which are linearly independent, so LICQ holds and the
multipliers are unique. A numerical solver (SLSQP, a sequential quadratic programming
code) returns the same point
with $f=5.000$.

**Global optimality.** The objective is strictly convex and the constraints
are linear, so by {prf:ref}`thm-06-03-convex-kkt` the KKT point is the unique
global minimum. Geometrically, $\vx^*$ is the point of the polyhedron
closest to the origin: $x_1\ge1$, $x_2\ge2$ and $x_3\ge0$ cannot be
improved simultaneously, and the corner $(1,2,0)$ of these three bounds
happens to satisfy the other two constraints.

**Interpretation.** $\lambda_3^*=2$ and $\lambda_4^*=4$ are the marginal
costs of the bounds $x_1\ge1$ and $x_2\ge2$: raising the first bound to
$1+\delta$ raises the optimal value to $(1+\delta)^2+4$, whose derivative at
$\delta=0$ is $2$, and raising the second to $2+\delta$ gives
$1+(2+\delta)^2$, with derivative $4$, in agreement with
{eq}`eq-06-03-sensitivity`. The constraint $x_3\ge0$ is active but costs
nothing ($\lambda_5^*=0$), because $x_3=0$ would be optimal even without it.
:::

(sec-06-03-second-order)=
## Second-order conditions

At a KKT point the first-order terms balance, and as in
{doc}`06-02-optimality-conditions` curvature decides between minimum and
saddle. Two changes are needed for curved constraints: the curvature is
measured by the Hessian of the **Lagrangian**, which includes the curvature
of the constraints weighted by their multipliers, and it is tested only on
the directions along which the first-order terms give no decision.

:::{prf:definition} Critical cone
:label: def-06-03-critical-cone

Let $(\vx^*,\vlambda^*,\vmu^*)$ satisfy the KKT conditions. The **critical
cone** is

$$
\begin{aligned}
C(\vx^*,\vlambda^*)=\bigl\{\vd\in L(\vx^*):\ &\nabla g_i(\vx^*)^\T\vd=0\\
&\text{for all }i\in I(\vx^*)\text{ with }\lambda_i^*>0\bigr\}.
\end{aligned}
$$ (eq-06-03-critical-cone)
:::

For $\vd\in L(\vx^*)$, stationarity gives (the equality terms drop out
because $\nabla h_j(\vx^*)^\T\vd=0$)
$\nabla f(\vx^*)^\T\vd=-\sum_{i\in I}\lambda_i^*\nabla g_i(\vx^*)^\T\vd\ge0$,
with equality exactly when $\nabla g_i^\T\vd=0$ for every $i$ with
$\lambda_i^*>0$. So the critical cone consists of the linearized feasible
directions with $\nabla f(\vx^*)^\T\vd=0$: the analog of condition 2 of
{prf:ref}`thm-06-02-second-order-necessary`.

:::{prf:theorem} Second-order conditions for constrained problems
:label: thm-06-03-second-order

Let $f$, $g_i$, $h_j$ be twice continuously differentiable and write
$\nabla^2_{\vx\vx}L^*=\nabla^2f(\vx^*)+\sum_i\lambda_i^*\nabla^2g_i(\vx^*)+\sum_j\mu_j^*\nabla^2h_j(\vx^*)$.

1. (*Necessary*) If $\vx^*$ is a local minimum at which LICQ holds and
   $(\vlambda^*,\vmu^*)$ are its (unique) multipliers, then
   $\vd^\T\nabla^2_{\vx\vx}L^*\vd\ge0$ for all $\vd\in C(\vx^*,\vlambda^*)$.
2. (*Sufficient*) If $\vx^*$ is feasible, $(\vx^*,\vlambda^*,\vmu^*)$
   satisfies the KKT conditions, and

   $$
   \vd^\T\nabla^2_{\vx\vx}L(\vx^*,\vlambda^*,\vmu^*)\,\vd>0\qquad\text{for all }\vd\in C(\vx^*,\vlambda^*),\ \vd\ne\vzero,
   $$ (eq-06-03-sosc)

   then $\vx^*$ is a strict local minimum. No constraint qualification is
   needed.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

We prove part 2; for part 1 see Nocedal and Wright (2006, Theorem 12.5).
Suppose $\vx^*$ is not a strict local minimum. Then there are feasible
$\vx^{(k)}\ne\vx^*$ with $\vx^{(k)}\to\vx^*$ and $f(\vx^{(k)})\le f(\vx^*)$.
Let $t_k=\norm{\vx^{(k)}-\vx^*}$ and $\vd^{(k)}=(\vx^{(k)}-\vx^*)/t_k$;
passing to a subsequence, $\vd^{(k)}\to\vd$ with $\norm{\vd}=1$, and
$\vd\in T(\vx^*)\subseteq L(\vx^*)$ by {prf:ref}`lem-06-03-tangent`.

*$\vd$ is critical.* From $f(\vx^{(k)})-f(\vx^*)\le0$, dividing by $t_k$ and
letting $k\to\infty$ gives $\nabla f(\vx^*)^\T\vd\le0$. Since also
$\nabla f(\vx^*)^\T\vd\ge0$ on $L(\vx^*)$ (shown before the theorem),
$\nabla f(\vx^*)^\T\vd=0$, and $\vd\in C(\vx^*,\vlambda^*)$.

*Curvature.* Write $L(\vx)$ for $L(\vx,\vlambda^*,\vmu^*)$. Because
$\lambda_i^*\ge0$, $g_i(\vx^{(k)})\le0$ and $h_j(\vx^{(k)})=0$, and by
complementary slackness at $\vx^*$,

$$
L(\vx^{(k)})\le f(\vx^{(k)})\le f(\vx^*)=L(\vx^*).
$$

Taylor's theorem for $L$ at $\vx^*$, where $\nabla_{\vx}L=\vzero$, gives
$L(\vx^{(k)})=L(\vx^*)+\tfrac12t_k^2(\vd^{(k)})^\T\nabla^2_{\vx\vx}L(\vy^{(k)})\vd^{(k)}$
for some $\vy^{(k)}$ between $\vx^*$ and $\vx^{(k)}$. Hence
$(\vd^{(k)})^\T\nabla^2_{\vx\vx}L(\vy^{(k)})\vd^{(k)}\le0$, and in the limit
$\vd^\T\nabla^2_{\vx\vx}L^*\vd\le0$ for a nonzero critical $\vd$,
contradicting {eq}`eq-06-03-sosc`.
:::

In the unconstrained case $L=f$ and $C=\R^n$, and the theorem reduces to
{prf:ref}`thm-06-02-second-order-necessary` and
{prf:ref}`thm-06-02-second-order-sufficient`. With only linear constraints,
$\nabla^2_{\vx\vx}L=\nabla^2f$; in {prf:ref}`ex-06-03-kkt-example`,
$\nabla^2f=2\mI\succ\vzero$ and the sufficient condition holds trivially.
The constraint curvature matters when the constraints are curved, as the
following example shows.

**Example: a linear objective on a circle.** Minimize $f(\vx)=x_1+x_2$
subject to $h(\vx)=x_1^2+x_2^2-2=0$. Stationarity
$(1,1)+\mu(2x_1,2x_2)=\vzero$ gives $x_1=x_2=-1/(2\mu)$, and the constraint
gives $2/(4\mu^2)=2$, so $\mu=\pm\tfrac12$. The two KKT points are

- $\vx=(-1,-1)$ with $\mu=\tfrac12$ and $f=-2$;
- $\vx=(1,1)$ with $\mu=-\tfrac12$ and $f=2$.

LICQ holds at both ($\nabla h=(2x_1,2x_2)\ne\vzero$). The Hessian of $f$ is
zero, so a test on $f$ alone could not tell the points apart. The Hessian of
the Lagrangian is $\nabla^2_{\vx\vx}L=\nabla^2f+\mu\nabla^2h=2\mu\mI$, and the
critical cone is the tangent line $\{\vd:\nabla h^\T\vd=0\}=\{t(1,-1)\}$.
For $\vd=(1,-1)$, $\vd^\T\nabla^2_{\vx\vx}L\,\vd=4\mu$. At $(-1,-1)$ this is
$2>0$: the sufficient condition holds and the point is a strict local
minimum (in fact the global one). At $(1,1)$ it is $-2<0$: the necessary
condition fails, and the point is not a local minimum; it is the maximum of
$f$ on the circle.

(sec-06-03-convex)=
## KKT conditions for convex programs

For convex programs the KKT conditions are sufficient as well as
necessary, and the conclusion is global.

:::{prf:theorem} KKT points of convex programs are global minima
:label: thm-06-03-convex-kkt

Let $f$ and $g_1,\dots,g_m$ be convex and continuously differentiable and
$h_1,\dots,h_\ell$ affine. If $(\vx^*,\vlambda^*,\vmu^*)$ satisfies the KKT
conditions {eq}`eq-06-03-kkt-stationarity`–{eq}`eq-06-03-kkt-cs`, then $\vx^*$ is a
global minimum. Conversely, if Slater's condition holds, every global
minimum is a KKT point.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Since $\lambda_i^*\ge0$ and $h_j$ is affine, $\vx\mapsto L(\vx,\vlambda^*,\vmu^*)$
is a convex function, and stationarity says that its gradient vanishes at
$\vx^*$. By the gradient inequality of convex functions
({prf:ref}`thm-02-02-first-order`), $\vx^*$ minimizes $L(\cdot,\vlambda^*,\vmu^*)$
over $\R^n$. For any feasible $\vx$,

$$
\begin{aligned}
f(\vx)&\ge f(\vx)+\sum_i\lambda_i^*g_i(\vx)+\sum_j\mu_j^*h_j(\vx)\\
&=L(\vx,\vlambda^*,\vmu^*)\\
&\ge L(\vx^*,\vlambda^*,\vmu^*)=f(\vx^*),
\end{aligned}
$$

where the first inequality uses $\lambda_i^*g_i(\vx)\le0$ and $h_j(\vx)=0$,
and the last equality uses complementary slackness and $\vh(\vx^*)=\vzero$.
The converse is {prf:ref}`thm-06-03-kkt` with Slater's condition as the
constraint qualification.
:::

No constraint qualification is needed for sufficiency, and the equality
constraints must be affine (a nonlinear equality makes $\Omega$ nonconvex
even when $h_j$ is convex). *(Statements that require $\vx^*$ to be a
"regular point" or allow convex nonlinear $h_j$ in this theorem are,
respectively, stronger than necessary and incorrect.)* For convex programs
that satisfy Slater's condition the KKT conditions are therefore
**necessary and sufficient** for global optimality, and the multipliers are
optimal solutions of the Lagrangian dual
({doc}`02-03-convex-optimization-duality`). This is the setting of most
problems in {doc}`08-08-examples-nlp-optimality`, and the reason why the
solution methods of {doc}`06-06-constrained-algorithms` can be trusted for
convex problems.

## Summary and Common Mistakes

### Summary

- Only active constraints $I(\vx^*)$ restrict local moves; feasible
  directions satisfy the linearized conditions and lie in $L(\vx^*)$
  ({prf:ref}`prop-06-03-linearised`, {eq}`eq-06-03-linearised-cone`).
- The tangent cone $T(\vx^*)$ ({prf:ref}`def-06-03-tangent-cone`) describes the
  true local geometry; $T\subseteq L$ always, and at a local minimum
  $\nabla f^\T\vd\ge0$ on $T$ ({prf:ref}`lem-06-03-tangent`).
- A constraint qualification ensures $T=L$; linear constraints, LICQ and
  Slater's condition are sufficient ({prf:ref}`thm-06-03-cq-sufficient`). The
  cubic constraint ({prf:ref}`ex-06-03-cubic`) shows what goes wrong without
  one.
- Farkas' lemma (or LP duality applied to {eq}`eq-06-03-direction-lp`) then
  yields the KKT conditions: stationarity, primal and dual feasibility,
  complementary slackness ({prf:ref}`thm-06-03-kkt`), compactly
  $\nabla_{\vx}L=\vzero$ for the Lagrangian {eq}`eq-06-03-lagrangian`.
- Geometrically $-\nabla f(\vx^*)$ lies in the cone of the active constraint
  gradients ({numref}`fig-06-03-kkt-geometry`); multipliers measure the
  sensitivity of the optimal value ({eq}`eq-06-03-sensitivity`).
- Sign conventions for $\max$ and $\ge$ are collected in
  {numref}`tab-06-03-signs`.
- Active-set case analysis ({prf:ref}`alg-06-03-active-set`) solves small KKT
  systems; {prf:ref}`ex-06-03-kkt-example` has the unique KKT point $(1,2,0)$
  with $\vlambda^*=(0,0,2,4,0)$ and $f^*=5$.
- Second-order conditions use $\nabla^2_{\vx\vx}L$ on the critical cone
  ({prf:ref}`thm-06-03-second-order`).
- For convex programs every KKT point is a global minimum, and under Slater
  every minimum is a KKT point ({prf:ref}`thm-06-03-convex-kkt`).

### Common mistakes

- **Wrong signs.** Mixing the conventions of {numref}`tab-06-03-signs`, in
  particular for maximization problems or "$\ge$" constraints, gives
  multipliers of the wrong sign and rejects true optima. Bring the problem to
  the form $\min f$, $g_i\le0$ first.
- **Forgetting complementary slackness**, or imposing $\lambda_i>0$ for an
  active constraint. An active constraint may have $\lambda_i=0$
  ($g_5$ in {prf:ref}`ex-06-03-kkt-example`).
- **Restricting the sign of equality multipliers.** $\mu_j$ is free.
- **Treating KKT as sufficient.** Without convexity a KKT point can be a
  maximum (the point $(1,1)$ on the circle) or a saddle point. Use the
  second-order conditions with the Hessian of the *Lagrangian*, not of $f$.
- **Ignoring constraint qualifications.** At a point where the active
  gradients are dependent (for example a zero gradient as in
  {prf:ref}`ex-06-03-cubic`), a minimum may have no multipliers; search for
  such points separately.
- **Checking only $T=L$ via feasible directions.** On curved boundaries
  tangent directions are not feasible directions; the qualification is about
  the tangent cone.


## References

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
