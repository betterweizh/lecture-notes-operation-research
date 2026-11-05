# Convex Optimization and Lagrange Duality

## Overview

A convex optimization problem minimizes a convex function over a convex
set described by convex inequalities and affine equations. Such problems
have two properties that general problems lack: every local minimum is a
global minimum ({prf:ref}`thm-02-03-local-global`), and optimality can be
certified by a simple condition on the gradient
({prf:ref}`thm-02-03-optimality`). This chapter defines convex problems,
shows that linear and convex quadratic programs belong to the class, and
then develops **Lagrange duality**: the Lagrangian, the dual function
$q$, which gives lower bounds on the optimal value $p^*$, the dual problem
and its optimal value $d^*$, weak duality $d^*\le p^*$
({prf:ref}`thm-02-03-weak-duality`), strong duality $d^*=p^*$ under Slater's
condition ({prf:ref}`thm-02-03-slater`), complementary slackness and the
Karush–Kuhn–Tucker (KKT) conditions for convex problems
({prf:ref}`thm-02-03-kkt-sufficient`), and the interpretation of the
multipliers as sensitivities of the optimal value
({prf:ref}`thm-02-03-sensitivity`). The dual of a linear program is
derived from the Lagrangian, and the theory is worked through on an
economic-dispatch problem of an electricity system, on a least-norm
problem and on a small problem with a duality gap.
It builds on the separating and supporting hyperplanes of
{doc}`02-01-convex-sets`, on {doc}`02-02-convex-functions`, and on the
general form and basic terms of {doc}`01-03-introduction`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- recognize a convex optimization problem in standard form and explain
  why its equality constraints must be affine;
- prove that local minima of convex problems are global and apply the
  optimality condition $\nabla f(\vx^*)^\T(\vy-\vx^*)\ge0$;
- form the Lagrangian and the dual function of a problem and evaluate the
  dual function in closed form for linear and quadratic problems;
- state and use weak duality, strong duality, Slater's condition and
  complementary slackness;
- write the KKT conditions of a convex problem, solve them for small
  problems, and interpret the multipliers as prices; and
- derive the dual of a linear program from the Lagrangian.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\mathcal D$ | domain of a problem: the points where $f$, all $g_i$ and all $h_j$ are defined |
| $\ell$ | number of equality constraints |
| $L(\vx,\vlambda,\vmu)$ | Lagrangian |
| $\lambda_i$, $\mu_j$ | multipliers of $g_i(\vx)\le0$ and $h_j(\vx)=0$ |
| $q(\vlambda,\vmu)$ | Lagrange dual function |
| $p^*$, $d^*$ | optimal values of the primal and of the dual problem |
| $\vg(\vx)$, $\vh(\vx)$ | vectors $(g_1(\vx),\dots,g_m(\vx))^\T$ and $(h_1(\vx),\dots,h_\ell(\vx))^\T$ |
| $\mathcal G$, $(u,t)$ | set of (constraint, objective) values of a problem with one inequality, and its coordinates |
| $\vu$, $\vv$; $p^*(\vu,\vv)$ | perturbations of the right-hand sides of $g_i\le0$ and $h_j=0$; perturbed optimal value |

(sec-02-03-convex-problems)=
## Convex optimization problems

Recall the general problem of {doc}`01-03-introduction`:

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)\\
\st\ & g_i(\vx)\le0, && i=1,\dots,m,\\
& h_j(\vx)=0, && j=1,\dots,\ell,
\end{aligned}
$$ (eq-02-03-primal)

with domain $\mathcal D=\dom f\cap\bigcap_i\dom g_i\cap\bigcap_j\dom h_j$,
feasible set $\Omega$ and optimal value $p^*=\inf\{f(\vx):\vx\in\Omega\}$.

:::{prf:definition} Convex optimization problem
:label: def-02-03-convex-problem

Problem {eq}`eq-02-03-primal` is a **convex optimization problem** (convex
program) if $f$ and $g_1,\dots,g_m$ are convex functions and
$h_1,\dots,h_\ell$ are affine, $h_j(\vx)=\va_j^\T\vx-b_j$. The equality
constraints can then be written $\mA\vx=\vb$.
:::

The feasible set of a convex problem is convex: it is the intersection of
the sublevel sets $\{\vx:g_i(\vx)\le0\}$, which are convex by
{prf:ref}`thm-02-02-sublevel`, and the hyperplanes $\{\vx:\va_j^\T\vx=b_j\}$
({prf:ref}`thm-02-01-intersection`). A convex problem therefore minimizes a
convex function over a convex set. Maximizing a concave function over a
convex set is also a convex problem, because it is equivalent to
minimizing its negative ({eq}`eq-01-03-max-min`).

The equality constraints must be affine, not merely convex. The set
$\{\vx:h(\vx)=0\}=\{\vx:h(\vx)\le0\}\cap\{\vx:-h(\vx)\le0\}$ needs both $h$
and $-h$ to be convex for this argument, which means that $h$ is affine.
For a convex but nonlinear $h$ the set is usually not convex:
$x_1^2+x_2^2-1=0$ is a circle. *(Definitions that allow "convex equality
constraints" are incorrect for this reason; a convex $h$ may appear only in
an inequality $h(\vx)\le0$.)*

:::{prf:theorem} Local minima of convex problems are global
:label: thm-02-03-local-global

For a convex optimization problem, every local minimum is a global
minimum. The set of optimal solutions is convex, and if $f$ is strictly
convex there is at most one optimal solution.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\vx^*$ be a local minimum: $f(\vx^*)\le f(\vx)$ for all feasible $\vx$
with $\norm{\vx-\vx^*}<\varepsilon$. Suppose a feasible $\vy$ had
$f(\vy)<f(\vx^*)$. For $\theta\in(0,1]$ the point
$\vz_\theta=\vx^*+\theta(\vy-\vx^*)$ is feasible, because $\Omega$ is
convex, and by convexity of $f$

$$
f(\vz_\theta)\le(1-\theta)f(\vx^*)+\theta f(\vy)<f(\vx^*).
$$

For $\theta<\varepsilon/\norm{\vy-\vx^*}$ the point $\vz_\theta$ lies in the
$\varepsilon$-neighborhood of $\vx^*$, which contradicts local
optimality. Hence $f(\vy)\ge f(\vx^*)$ for every feasible $\vy$.

The optimal set is $\Omega\cap\{\vx:f(\vx)\le p^*\}$, an intersection of
convex sets. If $\vx^*\ne\vy^*$ were both optimal and $f$ strictly convex,
their midpoint would be feasible with
$f(\tfrac12\vx^*+\tfrac12\vy^*)<p^*$, which is impossible.
:::

This theorem is the main practical reason for the importance of
convexity. A local search method, which only compares a point with its
neighbors, finds a global optimum of a convex problem. For nonconvex
problems the situation of {prf:ref}`def-01-03-local-global` arises, and
{ref}`sec-06-01-difficulties` discusses the difficulties.

(sec-02-03-optimality)=
### An optimality condition

For a differentiable objective, the first-order condition
{eq}`eq-02-02-first-order` leads to a characterization of the optimal points.

:::{prf:theorem} Optimality condition for convex problems
:label: thm-02-03-optimality

Let $f$ be convex and differentiable and $\Omega$ convex. A point
$\vx^*\in\Omega$ minimizes $f$ over $\Omega$ if and only if

$$
\nabla f(\vx^*)^\T(\vy-\vx^*)\ge0\qquad\text{for all }\vy\in\Omega .
$$ (eq-02-03-optimality)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

If {eq}`eq-02-03-optimality` holds, the first-order condition gives
$f(\vy)\ge f(\vx^*)+\nabla f(\vx^*)^\T(\vy-\vx^*)\ge f(\vx^*)$ for every
$\vy\in\Omega$. Conversely, if $\nabla f(\vx^*)^\T(\vy-\vx^*)<0$ for some
$\vy\in\Omega$, the derivative of $\theta\mapsto f(\vx^*+\theta(\vy-\vx^*))$
at $\theta=0$ is negative, so $f(\vx^*+\theta(\vy-\vx^*))<f(\vx^*)$ for
small $\theta>0$; these points are feasible by convexity, and $\vx^*$ is
not optimal.
:::

Condition {eq}`eq-02-03-optimality` says that no feasible direction
$\vy-\vx^*$ is a descent direction. Geometrically, if
$\nabla f(\vx^*)\ne\vzero$, the hyperplane through $\vx^*$ with normal
$-\nabla f(\vx^*)$ supports $\Omega$ at $\vx^*$
({prf:ref}`thm-02-01-supporting-hyperplane`): the whole feasible set lies on
the side into which $f$ does not decrease to first order. Three special
cases are worth stating.

- **Unconstrained problems** ($\Omega=\R^n$): taking $\vy=\vx^*\pm\ve_j$
  gives $\nabla f(\vx^*)=\vzero$. For convex $f$, stationarity is
  necessary and sufficient for a global minimum.
- **Equality constraints** ($\Omega=\{\vx:\mA\vx=\vb\}$): the directions
  $\vy-\vx^*$ range over the null space of $\mA$, and
  {eq}`eq-02-03-optimality` for $\pm$ each direction forces
  $\nabla f(\vx^*)$ to be orthogonal to that null space, that is,
  $\nabla f(\vx^*)=-\mA^\T\vmu$ for some $\vmu$ (a vector in the row space
  of $\mA$, {doc}`01-01-linear-algebra`). The vector $\vmu$ is a vector of
  Lagrange multipliers.
- **Nonnegativity** ($\Omega=\R^n_+$): the condition becomes
  $\nabla f(\vx^*)\ge\vzero$ and $x_j^*\,\partial f(\vx^*)/\partial x_j=0$
  for all $j$. A component that is positive must have a zero partial
  derivative, and a component at its bound zero must have a nonnegative
  one. This is the first instance of complementary slackness.

(sec-02-03-lp-qp)=
## Linear and quadratic programs

A **linear program** (LP) has an affine objective and affine constraints; in
inequality form,

$$
\begin{aligned}
\min_{\vx}\ & \vc^\T\vx\\
\st\ & \mA\vx\le\vb ,
\end{aligned}
$$

and in standard form $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$. Affine
functions are convex, so every LP is a convex problem, its feasible set is
a polyhedron ({prf:ref}`def-02-01-polyhedron`), and every local optimum is
global. The LP chapters {doc}`03-01-lp-model-building` to
{doc}`03-07-sensitivity-analysis` exploit the additional structure: an
optimal solution of an LP in standard form, if one exists, can be found
at an extreme point ({prf:ref}`thm-03-02-fundamental`).

A **quadratic program** (QP) has a quadratic objective and affine
constraints,

$$
\begin{aligned}
\min_{\vx}\ & \tfrac12\vx^\T\mQ\vx+\vc^\T\vx\\
\st\ & \mA\vx\le\vb,\\
& \mE\vx=\vd,
\end{aligned}
$$ (eq-02-03-qp)

with symmetric $\mQ$. It is a convex problem if and only if $\mQ\succeq0$
({doc}`02-02-convex-functions`). Least-squares fitting with linear
restrictions, portfolio selection with variance as risk, and the dispatch
of generators with quadratic fuel costs are convex QPs. With an indefinite
$\mQ$ the problem is nonconvex and can have several local minima. The
following example is used throughout the chapter.

:::{prf:example} Economic dispatch of two generators
:label: ex-02-03-dispatch

An electricity system must supply a demand of $D=475$ MW from two
generators. Generator 1 produces $x_1$ MW at a cost of
$C_1(x_1)=8x_1+0.01x_1^2$ \$/h, generator 2 produces $x_2$ MW at
$C_2(x_2)=10x_2+0.02x_2^2$ \$/h. The **marginal costs**
$C_1'(x_1)=8+0.02x_1$ and $C_2'(x_2)=10+0.04x_2$ (\$/MWh) increase with
output, as they do for thermal plants. The dispatch problem is

$$
\begin{aligned}
\min_{\vx}\ & 8x_1+0.01x_1^2+10x_2+0.02x_2^2\\
\st\ & x_1+x_2=475 .
\end{aligned}
$$ (eq-02-03-dispatch)

It is a convex QP with $\mQ=\diag(0.02,0.04)\succ0$ and one affine
equality. (Output limits $x_1,x_2\ge0$ are omitted because they turn out
to be inactive; an upper limit is added in
{prf:ref}`ex-02-03-dispatch-capacity`.)

Find the optimal dispatch with the optimality condition of
{prf:ref}`thm-02-03-optimality`, and interpret its multiplier.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Write the constraint as $h(\vx)=475-x_1-x_2=0$, with
$\nabla h=(-1,-1)^\T$. By the equality case of
{prf:ref}`thm-02-03-optimality`, the optimum satisfies
$\nabla f(\vx^*)+\mu\nabla h(\vx^*)=\vzero$ for some $\mu$. Hence

$$
\begin{aligned}
8+0.02x_1&=\mu,\\
10+0.04x_2&=\mu :
\end{aligned}
$$

**both generators run at the same marginal cost $\mu$.** If one generator
had a lower marginal cost than the other, shifting a little output to it
would save money. Substituting $x_1=(\mu-8)/0.02$ and
$x_2=(\mu-10)/0.04$ into $x_1+x_2=475$ gives $50\mu-400+25\mu-250=475$,
so $\mu=15$ \$/MWh and

$$
\begin{aligned}
\vx^*&=(350,\ 125)\ \text{MW},\\
p^*&=C_1(350)+C_2(125)=4025+1562.5=5587.5\ \text{\$/h}.
\end{aligned}
$$

The common marginal cost $\mu=15$ \$/MWh is the price of electricity in
this system: serving one more MW of demand costs about \$15 per hour. The
Lagrangian theory below explains this interpretation.
:::

(sec-02-03-lagrangian)=
## The Lagrangian and the dual function

From now on, problem {eq}`eq-02-03-primal` is not required to be convex.
The idea of Lagrange duality is to move the constraints into the
objective with prices attached.

:::{prf:definition} Lagrangian and dual function
:label: def-02-03-lagrangian

The **Lagrangian** of {eq}`eq-02-03-primal` is the function
$L:\mathcal D\times\R^m\times\R^\ell\to\R$,

$$
\begin{aligned}
L(\vx,\vlambda,\vmu)&=f(\vx)+\sum_{i=1}^m\lambda_ig_i(\vx)+\sum_{j=1}^\ell\mu_jh_j(\vx)\\
&=f(\vx)+\vlambda^\T\vg(\vx)+\vmu^\T\vh(\vx).
\end{aligned}
$$ (eq-02-03-lagrangian)

The numbers $\lambda_i$ and $\mu_j$ are the **Lagrange multipliers** (or
dual variables) of the constraints $g_i(\vx)\le0$ and $h_j(\vx)=0$. The
**Lagrange dual function** is

$$
q(\vlambda,\vmu)=\inf_{\vx\in\mathcal D}L(\vx,\vlambda,\vmu).
$$ (eq-02-03-dual-function)
:::

The infimum is taken over the domain without the constraints, and it may
be $-\infty$. Two facts make the dual function useful.

*It is concave*, whether or not the primal problem is convex. For each
fixed $\vx$, $L(\vx,\vlambda,\vmu)$ is an affine function of
$(\vlambda,\vmu)$, and $q$ is the pointwise infimum of this family of
affine functions, hence concave ({prf:ref}`thm-02-02-operations`, part 3).

*It bounds the optimal value from below* whenever $\vlambda\ge\vzero$.

:::{prf:theorem} Lower-bound property
:label: thm-02-03-lower-bound

If $\vlambda\ge\vzero$, then $q(\vlambda,\vmu)\le p^*$ for every $\vmu$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\tilde\vx$ be feasible. Then $g_i(\tilde\vx)\le0$ and
$\lambda_i\ge0$ give $\vlambda^\T\vg(\tilde\vx)\le0$, and
$\vh(\tilde\vx)=\vzero$. Hence

$$
\begin{aligned}
q(\vlambda,\vmu)&\le L(\tilde\vx,\vlambda,\vmu)\\
&=f(\tilde\vx)+\vlambda^\T\vg(\tilde\vx)+\vmu^\T\vh(\tilde\vx)\\
&\le f(\tilde\vx).
\end{aligned}
$$

Taking the infimum over all feasible $\tilde\vx$ gives $q\le p^*$. (If the
problem is infeasible, $p^*=+\infty$ and there is nothing to prove.)
:::

The Lagrangian has an economic reading. Think of $\lambda_i$ as a price
per unit of violation of constraint $i$. Minimizing $L$ lets the decision
maker break the constraints, provided that they pay for it, and rewards
them at the same price for slack ($g_i(\vx)<0$ earns
$-\lambda_ig_i(\vx)>0$). With nonnegative prices, every feasible point
pays nothing or earns a reward, so the unconstrained minimum of $L$ can
only be lower than the constrained minimum of $f$. The bound is best when
the prices are chosen so that, at the minimizer of $L$, the constraints
hold: this is the content of strong duality.

(sec-02-03-dual-problem)=
## The dual problem and weak duality

Each $(\vlambda,\vmu)$ with $\vlambda\ge\vzero$ gives a lower bound
$q(\vlambda,\vmu)$ on $p^*$. The best such bound is the optimal value of the
**Lagrange dual problem**

$$
\begin{aligned}
\max_{\vlambda,\vmu}\ & q(\vlambda,\vmu)\\
\st\ & \vlambda\ge\vzero .
\end{aligned}
$$ (eq-02-03-dual)

Its optimal value is denoted $d^*$. A pair $(\vlambda,\vmu)$ with
$\vlambda\ge\vzero$ and $q(\vlambda,\vmu)>-\infty$ is **dual feasible**,
and a dual feasible pair attaining $d^*$ is **dual optimal**, written
$(\vlambda^*,\vmu^*)$. In contrast to {eq}`eq-02-03-primal`, problem
{eq}`eq-02-03-dual` is *always* a convex problem: it maximizes a concave
function over a convex set. The condition $q>-\infty$ usually appears as
explicit constraints once $q$ is computed, as in the LP dual of
{ref}`sec-02-03-lp-dual`.

:::{prf:theorem} Weak duality
:label: thm-02-03-weak-duality

For every problem {eq}`eq-02-03-primal`, convex or not, $d^*\le p^*$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Take the supremum over $\vlambda\ge\vzero$ and $\vmu$ in
{prf:ref}`thm-02-03-lower-bound`.
:::

The difference $p^*-d^*\ge0$ is the **duality gap**. Weak duality has two
immediate uses. A dual feasible point certifies a lower bound: if a
feasible $\tilde\vx$ and a dual feasible $(\vlambda,\vmu)$ satisfy
$f(\tilde\vx)-q(\vlambda,\vmu)\le\varepsilon$, then $\tilde\vx$ is
optimal to within $\varepsilon$, and many algorithms stop on this test.
And if the dual is unbounded, $d^*=+\infty$, the primal is infeasible.

:::{prf:example} Dual function of the dispatch problem
:label: ex-02-03-dispatch-dual

Compute the dual function of {eq}`eq-02-03-dispatch` with
$h(\vx)=475-x_1-x_2$, and solve the dual problem.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The Lagrangian separates into one term per generator:

$$
\begin{aligned}
L(\vx,\mu)
&=8x_1+0.01x_1^2+10x_2+0.02x_2^2+\mu\,(475-x_1-x_2)\\
&=\bigl(0.01x_1^2-(\mu-8)x_1\bigr)+\bigl(0.02x_2^2-(\mu-10)x_2\bigr)+475\mu .
\end{aligned}
$$

Each bracket is a convex parabola. The first is minimized at
$x_1=(\mu-8)/0.02$ with minimum $-(\mu-8)^2/0.04$, the second at
$x_2=(\mu-10)/0.04$ with minimum $-(\mu-10)^2/0.08$. Hence

$$
q(\mu)=475\mu-\frac{(\mu-8)^2}{0.04}-\frac{(\mu-10)^2}{0.08},
$$ (eq-02-03-dispatch-q)

a concave parabola, finite for all $\mu$. Each generator, offered the
price $\mu$ for its output, chooses the output at which its marginal cost
equals $\mu$, and $q(\mu)$ is the total cost of these outputs plus the
value, at the price $\mu$, of the shortfall $475-x_1-x_2$ of supply
below demand (negative when supply exceeds demand).

The equality constraint has a free multiplier, so the dual problem is the
unconstrained maximization of $q$. Setting
$q'(\mu)=475-(\mu-8)/0.02-(\mu-10)/0.04=0$ gives $75\mu=1125$, so
$\mu^*=15$ \$/MWh and

$$
d^*=q(15)=7125-1225-312.5=5587.5\ \text{\$/h}=p^* .
$$

{numref}`fig-02-03-dispatch-dual` plots $q$. Every value of $q$ is a lower
bound on the optimal cost (for example $q(10)=4650$ and $q(20)=4650$ \$/h),
the best bound equals $p^*$, and the maximizer $\mu^*$ reproduces the
equal marginal cost of {prf:ref}`ex-02-03-dispatch`. At $\mu^*$ the
generators' own choices $x_1=350$, $x_2=125$ meet the demand exactly.

```{figure} figures/02-03-dispatch-dual.svg
:name: fig-02-03-dispatch-dual
:width: 520px
:class: dark-light
:alt: A concave parabola q of mu with its maximum 5587.5 at mu equal to 15, touching the horizontal line p-star equal to 5587.5

The dual function {eq}`eq-02-03-dispatch-q` of the dispatch problem. Every
value $q(\mu)$ is a lower bound on the optimal cost $p^*=5587.5$ \$/h
(dashed); the maximum, at the price $\mu^*=15$ \$/MWh, equals $p^*$, so
there is no duality gap.
```
:::

(sec-02-03-strong-duality)=
## Strong duality and Slater's condition

**Strong duality** means $d^*=p^*$. It does not hold in general (see
{ref}`sec-02-03-duality-gap`), but it holds for most convex problems.
Conditions that guarantee it are called **constraint qualifications**; the
simplest is due to Slater.

:::{prf:theorem} Slater's condition
:label: thm-02-03-slater

Consider a convex problem {eq}`eq-02-03-primal` with $\vh(\vx)=\mA\vx-\vb$.
Suppose it is **strictly feasible**: there is a point $\tilde\vx$ in the
relative interior of $\mathcal D$ with

$$
\begin{aligned}
g_i(\tilde\vx)&<0\quad(i=1,\dots,m),\\
\mA\tilde\vx&=\vb .
\end{aligned}
$$

Then $d^*=p^*$, and if $p^*$ is finite the dual optimum is attained: there
is a dual optimal $(\vlambda^*,\vmu^*)$. The strict inequality is not
needed for constraints $g_i$ that are affine.
:::

The proof separates two convex sets by a hyperplane
({prf:ref}`thm-02-01-separating-hyperplane`); Boyd and Vandenberghe (2004,
§5.3.2) give it in full. Its idea is visible in the geometric picture
below. Two consequences are important for this book:

- for a convex problem whose constraints are all **affine**, such as an LP
  or a QP, Slater's condition reduces to feasibility. If such a problem is
  feasible and $p^*$ is finite, strong duality holds and the dual optimum
  is attained. For LP this is the strong duality theorem
  {prf:ref}`thm-03-05-strong-duality`;
- in the domain $\mathcal D=\R^n$, the relative interior is all of $\R^n$,
  and the condition just asks for a point that satisfies the nonlinear
  inequalities strictly.

The dispatch problem {eq}`eq-02-03-dispatch` has only an affine constraint
and is feasible, which explains the zero gap found in
{prf:ref}`ex-02-03-dispatch-dual`.

(sec-02-03-geometry)=
### A geometric view

For a problem with one inequality constraint and no equality constraint,
consider the set of values that the pair (constraint, objective) takes,

$$
\mathcal G=\{(g(\vx),\ f(\vx)):\vx\in\mathcal D\}\subseteq\R^2 ,
$$

with coordinates $(u,t)$. The optimal value is the lowest $t$ among the
points of $\mathcal G$ with $u\le0$:
$p^*=\inf\{t:(u,t)\in\mathcal G,\ u\le0\}$. The dual function is

$$
q(\lambda)=\inf\{t+\lambda u:(u,t)\in\mathcal G\} .
$$

For fixed $\lambda\ge0$ the lines $t+\lambda u=\text{const}$ have slope
$-\lambda$, and the infimum is attained by the lowest such line that still
touches $\mathcal G$: a supporting line of $\mathcal G$. It crosses the
axis $u=0$ at height $q(\lambda)$. Because the line has nonpositive slope
and lies below all of $\mathcal G$, in particular below the points with
$u\le0$, its intercept is at most $p^*$; this is weak duality in a
picture. The dual problem tilts the supporting line to raise the
intercept as much as possible. When $\mathcal G$ (more precisely, the set
$\mathcal G$ plus everything above and to the right of it) is convex,
the best supporting line passes through the point $(0,p^*)$ and $d^*=p^*$;
Slater's point, with $u<0$, prevents the line from becoming vertical.
{numref}`fig-02-03-duality-geometry` shows both a convex case and a
nonconvex case with a gap.

```{figure} figures/02-03-duality-geometry.svg
:name: fig-02-03-duality-geometry
:width: 680px
:class: dark-light
:alt: Left, a parabola-shaped curve G in the (u, t) plane with a supporting line of slope minus 3 touching it at u equal to 0, and a less steep line with a lower intercept. Right, two isolated points and the line through them, whose intercept 0.5 lies below the optimal value 1

Geometric view of duality with one inequality constraint. (a) The dispatch
problem with the capacity limit of {prf:ref}`ex-02-03-dispatch-capacity`,
reduced to the variable $x_1$: $\mathcal G$ is the black curve
$(x_1-300,\ C_1(x_1)+C_2(475-x_1))$. The supporting line with
slope $-\lambda^*=-3$ meets $u=0$ at $q(3)=5662.5=p^*$; the line with
slope $-1$ meets it lower, at $q(1)=5629.17$. (b) The on/off pump of
{prf:ref}`ex-02-03-duality-gap`: $\mathcal G$ has two points, and the best
supporting line meets $u=0$ at $d^*=0.5<p^*=1$. The gray band marks the
feasible side $u\le0$.
```

(sec-02-03-complementary-slackness)=
## Complementary slackness

Suppose strong duality holds, $\vx^*$ is primal optimal and
$(\vlambda^*,\vmu^*)$ is dual optimal. Then

$$
\begin{aligned}
f(\vx^*)=q(\vlambda^*,\vmu^*)
&=\inf_{\vx\in\mathcal D}\Bigl(f(\vx)+\sum_i\lambda_i^*g_i(\vx)+\sum_j\mu_j^*h_j(\vx)\Bigr)\\
&\le f(\vx^*)+\sum_i\lambda_i^*g_i(\vx^*)+\sum_j\mu_j^*h_j(\vx^*)\\
&\le f(\vx^*).
\end{aligned}
$$

The last step uses $\lambda_i^*\ge0$, $g_i(\vx^*)\le0$ and
$h_j(\vx^*)=0$. The first and last expressions are equal, so both
inequalities are equalities. Two conclusions follow.

1. $\vx^*$ minimizes the Lagrangian $L(\vx,\vlambda^*,\vmu^*)$ over
   $\mathcal D$. The optimal multipliers turn the constrained problem
   into an unconstrained one that has $\vx^*$ among its solutions.
2. $\sum_i\lambda_i^*g_i(\vx^*)=0$. Each term is $\le0$, so each is zero:

$$
\lambda_i^*\,g_i(\vx^*)=0,\qquad i=1,\dots,m .
$$ (eq-02-03-complementary-slackness)

Condition {eq}`eq-02-03-complementary-slackness` is **complementary
slackness**: $\lambda_i^*>0$ implies $g_i(\vx^*)=0$ (a constraint with a
positive price is active, or binding), and $g_i(\vx^*)<0$ implies
$\lambda_i^*=0$ (a constraint with slack has price zero). Economically, a
resource that is not used up has no value at the margin. In LP this is
the complementary slackness theorem
{prf:ref}`thm-03-05-complementary-slackness`.

(sec-02-03-kkt)=
## KKT conditions for convex problems

Now let $f$, $g_i$ and $h_j$ be differentiable. Under strong duality, the
optimal $\vx^*$ minimizes $L(\cdot,\vlambda^*,\vmu^*)$, so its gradient
with respect to $\vx$ vanishes at $\vx^*$. Collecting this with
feasibility and complementary slackness gives the **Karush–Kuhn–Tucker
(KKT) conditions**:

$$
\begin{aligned}
&\text{primal feasibility:} && g_i(\vx)\le0 && \text{for all }i,\\
& && h_j(\vx)=0 && \text{for all }j,\\
&\text{dual feasibility:} && \lambda_i\ge0 && \text{for all }i,\\
&\text{complementary slackness:} && \lambda_ig_i(\vx)=0 && \text{for all }i,\\
&\text{stationarity:} && \nabla f(\vx)+\sum_{i=1}^m\lambda_i\nabla g_i(\vx)+\sum_{j=1}^\ell\mu_j\nabla h_j(\vx)=\vzero .
\end{aligned}
$$ (eq-02-03-kkt)

The argument just given shows: **for any problem with differentiable
functions in which strong duality holds, every primal optimal $\vx^*$ and
dual optimal $(\vlambda^*,\vmu^*)$ satisfy the KKT conditions.** For
convex problems the conditions are also sufficient.

:::{prf:theorem} KKT conditions are sufficient for convex problems
:label: thm-02-03-kkt-sufficient

Let {eq}`eq-02-03-primal` be a convex problem with differentiable $f$ and
$g_i$. If $(\tilde\vx,\tilde\vlambda,\tilde\vmu)$ satisfies
{eq}`eq-02-03-kkt`, then $\tilde\vx$ is primal optimal,
$(\tilde\vlambda,\tilde\vmu)$ is dual optimal, and $p^*=d^*$. If, in
addition, Slater's condition holds, then a feasible $\vx^*$ is optimal if
and only if there are multipliers $(\vlambda^*,\vmu^*)$ that together with
$\vx^*$ satisfy {eq}`eq-02-03-kkt`.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Since $\tilde\lambda_i\ge0$, the function
$L(\vx,\tilde\vlambda,\tilde\vmu)=f(\vx)+\sum_i\tilde\lambda_ig_i(\vx)+
\sum_j\tilde\mu_jh_j(\vx)$ is convex in $\vx$: a nonnegative combination
of convex functions plus affine functions (the sign of $\tilde\mu_j$ does
not matter because $h_j$ is affine). Stationarity says that its gradient
vanishes at $\tilde\vx$, so by {eq}`eq-02-02-first-order` $\tilde\vx$ is a
global minimizer of $L(\cdot,\tilde\vlambda,\tilde\vmu)$. Hence, using
complementary slackness and $\vh(\tilde\vx)=\vzero$ in the last step,

$$
\begin{aligned}
q(\tilde\vlambda,\tilde\vmu)&=L(\tilde\vx,\tilde\vlambda,\tilde\vmu)\\
&=f(\tilde\vx)+\sum_i\tilde\lambda_ig_i(\tilde\vx)+\sum_j\tilde\mu_jh_j(\tilde\vx)\\
&=f(\tilde\vx).
\end{aligned}
$$

So the feasible point $\tilde\vx$ and the dual feasible point
$(\tilde\vlambda,\tilde\vmu)$ have equal objective values. By weak duality
$f(\tilde\vx)\ge p^*\ge d^*\ge q(\tilde\vlambda,\tilde\vmu)=f(\tilde\vx)$,
so all are equal: $\tilde\vx$ is primal optimal and
$(\tilde\vlambda,\tilde\vmu)$ dual optimal.

Under Slater's condition, strong duality holds and the dual optimum is
attained ({prf:ref}`thm-02-03-slater`), so by the argument before the
theorem every optimal $\vx^*$ satisfies the KKT conditions with a dual
optimal pair.
:::

The KKT conditions generalize the conditions $\nabla f(\vx^*)=\vzero$
(no constraints) and $\nabla f(\vx^*)=-\mA^\T\vmu$ (equality constraints)
found in {ref}`sec-02-03-optimality`. For convex problems they reduce
optimization to solving a system of equations and inequalities. For
nonconvex problems the KKT conditions are necessary under a constraint
qualification, but no longer sufficient; {doc}`06-03-kkt-conditions` derives
them ({prf:ref}`thm-06-03-kkt`) from LP duality and Farkas' lemma, discusses constraint
qualifications and adds second-order conditions. The conditions were
derived by Karush (1939) and by Kuhn and Tucker (1951).

:::{prf:example} Dispatch with a capacity limit
:label: ex-02-03-dispatch-capacity

In {prf:ref}`ex-02-03-dispatch`, suppose generator 1 can produce at most
$300$ MW. Solve the problem with the KKT conditions, and interpret the
multipliers.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The problem is

$$
\begin{aligned}
\min_{\vx}\ & 8x_1+0.01x_1^2+10x_2+0.02x_2^2\\
\st\ & g(\vx)=x_1-300\le0,\\
& h(\vx)=475-x_1-x_2=0 .
\end{aligned}
$$

It is convex, and Slater's condition holds, for example with
$\tilde\vx=(250,225)$: $g(\tilde\vx)=-50<0$ and $h(\tilde\vx)=0$. So the
KKT conditions are necessary and sufficient. With multipliers $\lambda$
for $g$ and $\mu$ for $h$, stationarity reads

$$
\begin{aligned}
8+0.02x_1+\lambda-\mu&=0,\\
10+0.04x_2-\mu&=0,
\end{aligned}
$$

together with $x_1+x_2=475$, $x_1\le300$, $\lambda\ge0$ and
$\lambda(x_1-300)=0$.

*Case $\lambda=0$.* The conditions are those of {prf:ref}`ex-02-03-dispatch`,
with solution $x_1=350$, which violates $x_1\le300$. This case gives no
KKT point.

*Case $x_1=300$.* Then $x_2=175$, the second equation gives
$\mu=10+0.04\cdot175=17$, and the first gives
$\lambda=\mu-8-0.02\cdot300=17-14=3\ge0$. All conditions hold.

Hence $\vx^*=(300,175)$ MW, $\lambda^*=3$, $\mu^*=17$, and

$$
p^*=C_1(300)+C_2(175)=3300+2362.5=5662.5\ \text{\$/h},
$$

which is \$75 per hour more than without the limit. The dual function is
$q(\lambda,\mu)=475\mu-300\lambda-(\mu-\lambda-8)^2/0.04-(\mu-10)^2/0.08$
(computed as in {prf:ref}`ex-02-03-dispatch-dual`), and
$q(3,17)=8075-900-900-612.5=5662.5=p^*$, as the theorem guarantees.

The multipliers are prices. $\mu^*=17$ \$/MWh is the marginal cost of
demand: it equals the marginal cost of generator 2, which is now the only
generator that can respond. $\lambda^*=3$ \$/MWh is the value of capacity
at generator 1: it is the gap between the price $17$ and generator 1's
marginal cost $C_1'(300)=14$ \$/MWh, and raising the limit from 300 to
301 MW would lower the cost by about \$3 per hour (exactly by \$2.97).
Complementary slackness holds because the capacity constraint is binding.
:::

(sec-02-03-sensitivity)=
## Multipliers as sensitivities

The price interpretation of the multipliers can be made precise. Perturb
the right-hand sides of the constraints of {eq}`eq-02-03-primal`,

$$
\begin{aligned}
g_i(\vx)&\le u_i\quad(i=1,\dots,m),\\
h_j(\vx)&=v_j\quad(j=1,\dots,\ell),
\end{aligned}
$$

and let $p^*(\vu,\vv)$ be the optimal value of the perturbed problem, so
that $p^*(\vzero,\vzero)=p^*$. A positive $u_i$ relaxes constraint $i$, a
negative $u_i$ tightens it.

:::{prf:theorem} Global sensitivity inequality
:label: thm-02-03-sensitivity

Suppose strong duality holds for the unperturbed problem and the dual
optimum $(\vlambda^*,\vmu^*)$ is attained. Then for all $\vu$ and $\vv$

$$
p^*(\vu,\vv)\ge p^*-\vlambda^{*\T}\vu-\vmu^{*\T}\vv .
$$ (eq-02-03-sensitivity)

If, in addition, $p^*(\vu,\vv)$ is differentiable at $(\vzero,\vzero)$,
then $\partial p^*/\partial u_i=-\lambda_i^*$ and
$\partial p^*/\partial v_j=-\mu_j^*$ there.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\vx$ be feasible for the perturbed problem. Then
$\lambda_i^*g_i(\vx)\le\lambda_i^*u_i$ because $\lambda_i^*\ge0$, and
$\mu_j^*h_j(\vx)=\mu_j^*v_j$, so

$$
\begin{aligned}
p^*&=q(\vlambda^*,\vmu^*)\le f(\vx)+\vlambda^{*\T}\vg(\vx)+\vmu^{*\T}\vh(\vx)\\
&\le f(\vx)+\vlambda^{*\T}\vu+\vmu^{*\T}\vv .
\end{aligned}
$$

Minimizing over all such $\vx$ gives {eq}`eq-02-03-sensitivity`. For the
derivatives, take $\vu=t\ve_i$, $\vv=\vzero$: the inequality gives
$(p^*(t\ve_i,\vzero)-p^*)/t\ge-\lambda_i^*$ for $t>0$ and $\le-\lambda_i^*$
for $t<0$, so the derivative, if it exists, equals $-\lambda_i^*$. The
same argument applies to $v_j$.
:::

The inequality holds globally, not only for small perturbations, and it
can be read in two directions. If $\lambda_i^*$ is large, tightening
constraint $i$ ($u_i<0$) increases the optimal cost by at least
$\lambda_i^*\abs{u_i}$: the constraint is expensive. If $\lambda_i^*$ is
small, relaxing it ($u_i>0$) can lower the cost only a little to first
order. A constraint with $\lambda_i^*=0$ is not worth relaxing at the
margin, which is complementary slackness once more.

In {prf:ref}`ex-02-03-dispatch-capacity`, the capacity constraint
$x_1-300\le u$ allows $300+u$ MW at generator 1, and the theorem predicts
a cost change of about $-\lambda^*u=-3u$ \$/h; the exact change for
$u=1$ is $-\$2.97$ per hour. The demand constraint
$475-x_1-x_2=v$ corresponds to the demand $475-v$, so a demand increase
of one MW ($v=-1$) raises the cost by about $\mu^*=17$ \$/h (exactly
\$17.02 per hour). The multipliers of an LP have the same meaning; there
they are called shadow prices ({prf:ref}`thm-03-05-shadow-price`), and
{doc}`03-07-sensitivity-analysis` computes the ranges of perturbation over
which they stay valid ({ref}`sec-03-07-rhs`).

(sec-02-03-lp-dual)=
## Duality for linear programs

The Lagrangian gives the dual of every LP in a few lines, without the
tables of rules used in {doc}`03-05-lp-duality` ({ref}`sec-03-05-dual-rules`).

:::{prf:example} The LP dual from the Lagrangian
:label: ex-02-03-lp-dual

Derive the dual of the LP

$$
\begin{aligned}
\min_{\vx}\ & \vc^\T\vx\\
\st\ & \mA\vx\ge\vb,\\
& \vx\ge\vzero ,
\end{aligned}
$$ (eq-02-03-lp-primal)

with $\mA\in\R^{m\times n}$. Then do the same for the standard form
$\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Canonical form.* Write the constraints as $\vb-\mA\vx\le\vzero$, with
multipliers $\vpi\in\R^m$ (the notation of the LP chapters), and
$-\vx\le\vzero$, with multipliers $\vlambda\in\R^n$. The Lagrangian is

$$
\begin{aligned}
L(\vx,\vpi,\vlambda)&=\vc^\T\vx+\vpi^\T(\vb-\mA\vx)-\vlambda^\T\vx\\
&=\vb^\T\vpi+(\vc-\mA^\T\vpi-\vlambda)^\T\vx .
\end{aligned}
$$

It is affine in $\vx$. An affine function is bounded below on $\R^n$ only
if its linear part is zero, so

$$
q(\vpi,\vlambda)=
\begin{cases}
\vb^\T\vpi, & \mA^\T\vpi+\vlambda=\vc,\\
-\infty, & \text{otherwise.}
\end{cases}
$$

The dual problem maximizes $q$ over $\vpi\ge\vzero$, $\vlambda\ge\vzero$.
Making the implicit constraint explicit and eliminating
$\vlambda=\vc-\mA^\T\vpi\ge\vzero$ gives

$$
\begin{aligned}
\max_{\vpi}\ & \vb^\T\vpi\\
\st\ & \mA^\T\vpi\le\vc,\\
& \vpi\ge\vzero ,
\end{aligned}
$$ (eq-02-03-lp-dual)

the canonical dual pair of {doc}`03-05-lp-duality`
({prf:ref}`def-03-05-canonical-pair`). The eliminated
multipliers $\lambda_j=c_j-\mA_{\cdot j}^\T\vpi$ are the **reduced costs**
of the simplex method ({eq}`eq-03-03-reduced-cost`), and complementary
slackness {eq}`eq-02-03-complementary-slackness` becomes

$$
\begin{aligned}
\pi_i\,(\mA_{i\cdot}\vx-b_i)&=0,\quad i=1,\dots,m,\\
(c_j-\mA_{\cdot j}^\T\vpi)\,x_j&=0,\quad j=1,\dots,n :
\end{aligned}
$$

a constraint with slack has a zero dual price, and a positive variable has
a zero reduced cost.

*Standard form.* With $\vh(\vx)=\vb-\mA\vx=\vzero$ and multipliers
$\vpi$ (now free in sign) and $-\vx\le\vzero$ with $\vlambda\ge\vzero$, the
same calculation gives $q=\vb^\T\vpi$ if $\mA^\T\vpi+\vlambda=\vc$ and
$-\infty$ otherwise, so the dual is
$\max\{\vb^\T\vpi:\mA^\T\vpi\le\vc\}$ with $\vpi$ free.
:::

For the dual pair {eq}`eq-02-03-lp-primal`–{eq}`eq-02-03-lp-dual` derived
in the solution of {prf:ref}`ex-02-03-lp-dual`, all constraints are affine,
so Slater's condition reduces to feasibility
({ref}`sec-02-03-strong-duality`): if the primal is feasible and bounded,
$d^*=p^*$ and both optima are attained. Weak duality,
$\vb^\T\vpi\le\vc^\T\vx$ for feasible $\vx$ and $\vpi$, can also be seen
directly: $\vc^\T\vx\ge(\mA^\T\vpi)^\T\vx=\vpi^\T\mA\vx\ge\vpi^\T\vb$, using
$\vx\ge\vzero$ in the first and $\vpi\ge\vzero$ in the second inequality.
This is the bound found by hand for the water-supply model of
{prf:ref}`ex-01-03-water-supply`.

(sec-02-03-least-norm)=
## The least-norm solution of a linear system

:::{prf:example} Least-norm solution
:label: ex-02-03-least-norm

An underdetermined linear system $\mA\vx=\vb$, with $\mA\in\R^{m\times n}$
of full row rank $m<n$, has infinitely many solutions. Find the one of
smallest Euclidean norm,

$$
\begin{aligned}
\min_{\vx}\ & \vx^\T\vx\\
\st\ & \mA\vx=\vb ,
\end{aligned}
$$

through the dual problem, and evaluate it for

$$
\begin{aligned}
\mA&=\begin{pmatrix}1&1&1\\1&2&3\end{pmatrix},\\
\vb&=\begin{pmatrix}6\\14\end{pmatrix}.
\end{aligned}
$$
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Dual function.* With $\vh(\vx)=\mA\vx-\vb$ and multipliers
$\vmu\in\R^m$, the Lagrangian $L(\vx,\vmu)=\vx^\T\vx+\vmu^\T(\mA\vx-\vb)$
is a strictly convex quadratic in $\vx$. Its gradient
$2\vx+\mA^\T\vmu$ vanishes at $\vx=-\tfrac12\mA^\T\vmu$, and substituting,

$$
\begin{aligned}
q(\vmu)&=\tfrac14\vmu^\T\mA\mA^\T\vmu-\tfrac12\vmu^\T\mA\mA^\T\vmu-\vb^\T\vmu\\
&=-\tfrac14\vmu^\T\mA\mA^\T\vmu-\vb^\T\vmu ,
\end{aligned}
$$

a concave quadratic ($\mA\mA^\T\succ0$ because $\mA$ has full row rank).

*Dual problem.* Maximizing $q$ over $\vmu\in\R^m$: the gradient
$-\tfrac12\mA\mA^\T\vmu-\vb$ vanishes at
$\vmu^*=-2(\mA\mA^\T)^{-1}\vb$, with

$$
d^*=q(\vmu^*)=\vb^\T(\mA\mA^\T)^{-1}\vb .
$$

The constraints are affine and feasible, so strong duality holds. The
primal solution is the minimizer of $L(\cdot,\vmu^*)$:

$$
\vx^*=-\tfrac12\mA^\T\vmu^*=\mA^\T(\mA\mA^\T)^{-1}\vb ,
$$ (eq-02-03-least-norm)

and indeed
$\norm{\vx^*}^2=\vb^\T(\mA\mA^\T)^{-1}\mA\mA^\T(\mA\mA^\T)^{-1}\vb=d^*$.

*Numbers.* Here
$\mA\mA^\T=\begin{pmatrix}3&6\\6&14\end{pmatrix}$ with determinant
$42-36=6$, so

$$
(\mA\mA^\T)^{-1}\vb=\frac16\begin{pmatrix}14&-6\\-6&3\end{pmatrix}
\begin{pmatrix}6\\14\end{pmatrix}
=\frac16\begin{pmatrix}0\\6\end{pmatrix}
=\begin{pmatrix}0\\1\end{pmatrix}.
$$

Hence $\vx^*=\mA^\T(0,1)^\T=(1,2,3)$, $\vmu^*=(0,-2)$ and
$p^*=1+4+9=14=d^*=\vb^\T(0,1)^\T$. *Check:* $1+2+3=6$ and $1+4+9=14$, so
$\mA\vx^*=\vb$. Every other solution has the form $\vx^*+\vz$ with
$\mA\vz=\vzero$, and since $\vx^*$ lies in the row space of $\mA$,
$\norm{\vx^*+\vz}^2=\norm{\vx^*}^2+\norm{\vz}^2>14$ for $\vz\ne\vzero$.
:::

The same calculation applies to any convex QP with linear inequality
constraints. For $\min\{\vx^\T\mP\vx:\mA\vx\le\vb\}$ with $\mP\succ0$, the
Lagrangian $\vx^\T\mP\vx+\vlambda^\T(\mA\vx-\vb)$ is minimized at
$\vx=-\tfrac12\mP^{-1}\mA^\T\vlambda$, and

$$
q(\vlambda)=-\tfrac14\vlambda^\T\mA\mP^{-1}\mA^\T\vlambda-\vb^\T\vlambda ,
$$

to be maximized over $\vlambda\ge\vzero$. The dual is again a QP, with one
variable per constraint, which is attractive when there are fewer
constraints than variables.

(sec-02-03-duality-gap)=
## When strong duality fails

Strong duality can fail in two ways: the problem is not convex, or it is
convex but violates every constraint qualification.

:::{prf:example} A duality gap: an on/off pump
:label: ex-02-03-duality-gap

A pumping station has one pump that is either off ($x=0$) or on at full
output ($x=1$, in units of its capacity); running it costs $1$ (in units
of its daily running cost). The station must deliver at least half its
capacity. Compare the optimal value with the dual optimal value of

$$
\begin{aligned}
\min_{x}\ & x\\
\st\ & g(x)=0.5-x\le0,\\
& x\in\mathcal D=\{0,1\}.
\end{aligned}
$$
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Primal.* Only $x=1$ is feasible, so $p^*=1$.

*Dual.* The Lagrangian $L(x,\lambda)=x+\lambda(0.5-x)$ is evaluated at the
two points of $\mathcal D$:

$$
q(\lambda)=\min\{L(0,\lambda),\ L(1,\lambda)\}
=\min\{0.5\lambda,\ 1-0.5\lambda\}.
$$

The first term increases and the second decreases in $\lambda$; they are
equal at $\lambda=1$, so $d^*=q(1)=0.5$. The duality gap is
$p^*-d^*=0.5$.

The problem is not convex, because its domain $\{0,1\}$ is not a convex
set. In the picture of {numref}`fig-02-03-duality-geometry`(b),
$\mathcal G=\{(0.5,0),(-0.5,1)\}$ has two points, and no line through
$(0,1)$ supports it. The dual value $0.5$ is the optimal value of the
*relaxed* problem in which $x$ may take any value in $[0,1]$, the convex
hull of $\mathcal D$: the dual cannot "see" the integrality. This is a
first glimpse of why integer programs ({doc}`07-01-logical-constraints`) are
harder than convex ones, and why their relaxations give bounds.
:::

A **convex** problem without a strictly feasible point can still have
strong duality, but the dual optimum may not be attained. For
$\min\{x:x^2\le0\}$ the only feasible point is $x=0$, so $p^*=0$, and
Slater's condition fails. The dual function is
$q(\lambda)=\inf_x(x+\lambda x^2)=-1/(4\lambda)$ for $\lambda>0$ and
$-\infty$ for $\lambda=0$; its supremum is $d^*=0=p^*$, approached as
$\lambda\to\infty$ but never attained. No KKT point exists, since
stationarity $1+2\lambda x=0$ fails at $x=0$. Convex problems with a
genuine gap exist as well, but they need unusual domains: for
$\min\{e^{-x_1}:x_1^2/x_2\le0\}$ on $\mathcal D=\{\vx:x_2>0\}$, the
feasible points have $x_1=0$, so $p^*=1$, whereas $d^*=0$ (Boyd &
Vandenberghe, 2004, Exercise 5.21).

## Summary and Common Mistakes

### Summary

- A convex problem minimizes a convex $f$ subject to convex $g_i\le0$ and
  affine $h_j=0$ ({prf:ref}`def-02-03-convex-problem`); its local minima are
  global ({prf:ref}`thm-02-03-local-global`), and $\vx^*$ is optimal if and
  only if $\nabla f(\vx^*)^\T(\vy-\vx^*)\ge0$ for all feasible $\vy$
  ({prf:ref}`thm-02-03-optimality`). LPs and QPs with $\mQ\succeq0$ are convex.
- The Lagrangian {eq}`eq-02-03-lagrangian` prices the constraints; the dual
  function $q$ ({prf:ref}`def-02-03-lagrangian`) is concave and a lower bound
  on $p^*$ for $\vlambda\ge\vzero$ ({prf:ref}`thm-02-03-lower-bound`).
- The dual problem {eq}`eq-02-03-dual` is always convex, and $d^*\le p^*$
  (weak duality, {prf:ref}`thm-02-03-weak-duality`). For convex problems with
  a strictly feasible point, $d^*=p^*$ and the dual optimum is attained
  (Slater, {prf:ref}`thm-02-03-slater`); affine constraints need only
  feasibility.
- Under strong duality, $\vx^*$ minimizes $L(\cdot,\vlambda^*,\vmu^*)$ and
  $\lambda_i^*g_i(\vx^*)=0$ ({eq}`eq-02-03-complementary-slackness`); with
  differentiability these give the KKT conditions {eq}`eq-02-03-kkt`, which
  are sufficient for convex problems and, under Slater, necessary
  ({prf:ref}`thm-02-03-kkt-sufficient`).
- Multipliers are prices and sensitivities: perturbing constraint $i$ by
  $u_i$ changes the optimal value by at least $-\lambda_i^*u_i$
  ({prf:ref}`thm-02-03-sensitivity`); in economic dispatch, $\mu^*$ is the
  marginal cost of demand and $\lambda^*$ the value of extra capacity
  ({prf:ref}`ex-02-03-dispatch-dual`, {prf:ref}`ex-02-03-dispatch-capacity`).
- The LP dual follows from the Lagrangian ({prf:ref}`ex-02-03-lp-dual`); the
  least-norm solution is $\mA^\T(\mA\mA^\T)^{-1}\vb$
  ({prf:ref}`ex-02-03-least-norm`); without convexity a duality gap can occur
  ({prf:ref}`ex-02-03-duality-gap`, {numref}`fig-02-03-duality-geometry`).

### Common mistakes

- **Allowing nonlinear equality constraints in a convex problem.** The
  $h_j$ must be affine.
- **Forgetting the sign constraint $\vlambda\ge\vzero$.** Multipliers of
  inequalities $g_i\le0$ are nonnegative; multipliers of equalities are
  free. A constraint written as $g_i\ge0$ must first be turned around.
- **Minimizing the Lagrangian over the feasible set.** The dual function
  minimizes over the domain $\mathcal D$, without the constraints; that is
  what makes it easy to evaluate.
- **Assuming strong duality without checking a constraint qualification.**
  Weak duality always holds; strong duality needs, for example, convexity
  and Slater's condition.
- **Treating KKT points of nonconvex problems as global minima.**
  Sufficiency of the KKT conditions ({prf:ref}`thm-02-03-kkt-sufficient`) uses
  convexity.
- **Mixing up complementary slackness.** $\lambda_i^*>0$ forces the
  constraint to be active, but an active constraint may have
  $\lambda_i^*=0$.

## References

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 4, Sections
  4.2–4.3: KKT conditions and their sufficiency for convex problems;
  Chapter 6, Sections 6.1–6.2: Lagrangian duality, duality theorems and
  saddle-point optimality.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Section 4.1, motivation of LP duality
  through Lagrange multipliers.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Chapter 4,
  Sections 4.1–4.4: convex optimization problems, LP and QP; Chapter 5,
  Sections 5.1–5.6: the Lagrange dual function and problem, geometric
  interpretation, optimality conditions, perturbation and sensitivity
  analysis; Exercise 5.21, a convex problem without strong duality.)
- Karush, W. (1939). *Minima of functions of several variables with
  inequalities as side constraints* [Master's thesis, University of
  Chicago].
- Kuhn, H. W., & Tucker, A. W. (1951). Nonlinear programming. In J. Neyman
  (Ed.), *Proceedings of the Second Berkeley Symposium on Mathematical
  Statistics and Probability* (pp. 481–492). University of California
  Press.
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming*
  (4th ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3>
  (Chapter 14, duality and dual methods.)
