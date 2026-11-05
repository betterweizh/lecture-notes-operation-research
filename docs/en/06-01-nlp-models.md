# Nonlinear Programming Models

## Overview

A **nonlinear program** (NLP) is an optimization problem whose objective or
constraints are nonlinear functions of the decision variables. Many models
in engineering, economics and statistics are of this kind: fitting a model
to data by least squares, choosing flows on a congested network, sizing a
structure, pricing a product. This chapter sets up the general model
({eq}`eq-06-01-nlp`), formulates a least-squares estimation problem as an NLP
and solves it numerically ({prf:ref}`ex-06-01-consumer-regression`), compares
NLP with linear programming (LP), and introduces the classes of NLP that later
chapters treat separately: unconstrained problems, linearly constrained
problems, quadratic programs and convex programs. It then defines local and
global minima ({prf:ref}`def-06-01-local-global`), recalls when a minimum is
guaranteed to exist, and lists the difficulties that make NLP harder than
LP. The optimality conditions that characterize minima follow in
{doc}`06-02-optimality-conditions` and {doc}`06-03-kkt-conditions`; algorithms
follow in {doc}`06-04-line-search` to {doc}`06-06-constrained-algorithms`.
The chapter uses gradients, Hessians and quadratic forms from
{doc}`01-01-linear-algebra` and {doc}`01-02-matrix-calculus`, and compact
sets, convex sets and convex functions from {doc}`02-01-convex-sets` and
{doc}`02-02-convex-functions`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- write an optimization problem in the standard NLP form with inequality
  constraints $g_i(\vx)\le0$ and equality constraints $h_j(\vx)=0$;
- formulate a parameter-estimation problem as a constrained least-squares
  NLP and interpret its solution;
- explain why the optimum of an NLP need not lie at a vertex of the
  feasible set, unlike in linear programming;
- recognize unconstrained, linearly constrained, quadratic and convex
  programs, and write a quadratic objective as
  $\vc^\T\vx+\tfrac12\vx^\T\mQ\vx$ with a symmetric $\mQ$;
- distinguish local from global minima and use Weierstrass' theorem to
  show that a minimum exists; and
- name the main difficulties of nonlinear programming and recognize them in
  small examples.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $m$, $\ell$ | number of inequality and of equality constraints of an NLP |
| $\vg(\vx)$, $\vh(\vx)$ | vectors $(g_1(\vx),\dots,g_m(\vx))^\T$ and $(h_1(\vx),\dots,h_\ell(\vx))^\T$ |
| $B(\vx^*,\varepsilon)$ | open ball $\{\vx:\norm{\vx-\vx^*}<\varepsilon\}$ |
| $t$, $E_t$ | month index and monthly budget (\$) in {prf:ref}`ex-06-01-consumer-regression` |
| $p_{ti}$, $Q_{ti}$ | price (\$/unit) and purchased quantity (units) of product $i$ in month $t$ |
| $q_i$, $\beta_i$ | committed quantity (units) and budget share of product $i$ (parameters to estimate) |
| $r_{ti}$ | residual of product $i$ in month $t$ (\$) |

(sec-06-01-general-form)=
## The general nonlinear program

An NLP chooses $n$ real decision variables $\vx=(x_1,\dots,x_n)^\T$ so as
to minimize an objective function subject to inequality and equality
constraints:

$$
\begin{aligned}
\min_{\vx\in\R^n}\ & f(\vx)\\
\st\ & g_i(\vx)\le0, && i=1,\dots,m,\\
& h_j(\vx)=0, && j=1,\dots,\ell.
\end{aligned}
$$ (eq-06-01-nlp)

The functions $f,g_1,\dots,g_m,h_1,\dots,h_\ell:\R^n\to\R$ are given. In
vector notation the constraints read $\vg(\vx)\le\vzero$ and
$\vh(\vx)=\vzero$, and the problem is $\min\{f(\vx):\vx\in\Omega\}$ with the
**feasible set**

$$
\Omega=\{\vx\in\R^n:\ \vg(\vx)\le\vzero,\ \vh(\vx)=\vzero\}.
$$ (eq-06-01-feasible-set)

A point of $\Omega$ is a **feasible solution**; a feasible point with the
smallest objective value is an **optimal solution** or **minimum**, written
$\vx^*$. Throughout this part on NLP we assume that $f$, $g_i$ and $h_j$ are
**continuously differentiable**, and twice continuously differentiable when
second derivatives are used. Smoothness is what makes gradients, and hence
the optimality conditions of {doc}`06-02-optimality-conditions`, available.

The form {eq}`eq-06-01-nlp` covers every variant by simple rewriting:

- a maximization problem $\max f(\vx)$ has the same optimal solutions as
  $\min\,-f(\vx)$, and the optimal values satisfy $\max f=-\min(-f)$;
- a constraint $g(\vx)\ge b$ becomes $b-g(\vx)\le0$, and $g(\vx)\le b$
  becomes $g(\vx)-b\le0$;
- an equality $h(\vx)=0$ is equivalent to the pair $h(\vx)\le0$,
  $-h(\vx)\le0$, a device used in {doc}`06-03-kkt-conditions` to derive the
  multipliers of equality constraints;
- a sign restriction $x_j\ge0$ is the inequality $-x_j\le0$.

The "$\le0$" convention for inequalities matters for the signs of the
multipliers in {doc}`06-03-kkt-conditions`; every model in this part is first
brought to this form.

The objective is called $f$ throughout the NLP chapters, while $z$ is kept
for the objective value of a linear program. The following example shows
how the same objective behaves on three different feasible sets.

:::{prf:example} One objective, three feasible sets
:label: ex-06-01-three-feasible-sets

Let $f(\vx)=(x_1-2)^2+(x_2-2)^2$, the squared distance from the point
$(2,2)$. Consider three problems:

1. $\min f(\vx)$ over $\Omega_1=\R^2$ (no constraints);
2. $\min f(\vx)$ s.t. $x_1+x_2-5\le0$, $x_1\ge0$, $x_2\ge0$;
3. $\min f(\vx)$ s.t. $(x_1-1)^2+(x_2-1)^2-1\le0$, $x_1+x_2\ge1$,
   $x_1\ge0$, $x_2\ge0$.

Find the optimal solution and the optimal value of each problem, and say
which constraints are binding.
:::

::::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Problem 1.* $f\ge0$ everywhere and $f(2,2)=0$, so $\vx^*=(2,2)$ with
$f^*=0$.

*Problem 2.* The feasible set is the triangle with vertices $(0,0)$,
$(5,0)$, $(0,5)$. The point $(2,2)$ satisfies $2+2=4<5$ and lies inside
it, so the unconstrained minimum is feasible and therefore optimal:
$\vx^*=(2,2)$, $f^*=0$. No constraint is binding.

*Problem 3.* The first constraint is the closed disc of radius 1 centered at
$(1,1)$; it touches the axes at $(1,0)$ and $(0,1)$, so the sign
constraints are redundant. The second constraint cuts off the small corner
of the disc below the line $x_1+x_2=1$. The point $(2,2)$ lies outside the
disc, at distance $\sqrt2$ from its center. Minimizing the distance to
$(2,2)$ over a disc gives the point of the disc on the segment from the
center towards $(2,2)$:

$$
\begin{aligned}
\vx^*&=(1,1)+\tfrac{1}{\sqrt2}(1,1)\\
&=\Bigl(1+\tfrac{1}{\sqrt2},\ 1+\tfrac{1}{\sqrt2}\Bigr)\approx(1.7071,\ 1.7071),
\end{aligned}
$$

which satisfies $x_1+x_2\approx3.414\ge1$. Its objective value is
$f^*=2\bigl(1-\tfrac1{\sqrt2}\bigr)^2=3-2\sqrt2\approx0.1716$. Here the
optimum lies on the curved boundary, and only the disc constraint is
binding. A numerical solver (SLSQP, a sequential quadratic programming
code; Kraft, 1988) returns the same point.

```{figure} figures/06-01-three-feasible-sets.svg
:name: fig-06-01-three-feasible-sets
:width: 680px
:class: dark-light
:alt: Three panels with circular contours of the squared distance to the point (2,2); the second panel shows a triangle containing (2,2), the third a disc centered at (1,1) with a corner cut off, and the optimum on the disc boundary

Contours of $f(\vx)=(x_1-2)^2+(x_2-2)^2$ over the three feasible sets of
{prf:ref}`ex-06-01-three-feasible-sets`. In problems 1 and 2 the optimum
(orange) is the unconstrained minimizer $(2,2)$, in the interior of
$\Omega$; in problem 3 it is the point where the smallest contour touches
the disc. None of the three optima is a vertex.
```
::::

(sec-06-01-regression)=
## Formulating an NLP: least-squares estimation

Statistical estimation is one of the most common sources of NLPs: the
unknown parameters of a model are chosen to minimize a measure of misfit
between the model's predictions and observed data. When the predictions
depend nonlinearly on the parameters, and the parameters must satisfy
side conditions, the result is a constrained NLP.

:::{prf:example} Estimating a consumer's spending model
:label: ex-06-01-consumer-regression

A household spends a monthly budget of $E$ dollars on three products
$i=1,2,3$ with prices $p_i$ (\$/unit). A classical model of consumer
demand (a linear expenditure system) assumes that the household first buys
a committed quantity $q_i\ge0$ of each product and then divides the
remaining, uncommitted money $E-\sum_j p_jq_j$ among the products in fixed
shares $\beta_i\ge0$ with $\beta_1+\beta_2+\beta_3=1$. The share $\beta_i$
is the household's marginal propensity to spend uncommitted money on
product $i$. The model predicts the spending on product $i$ as

$$
p_iq_i+\beta_i\Bigl(E-\sum_{j=1}^3p_jq_j\Bigr),\qquad i=1,2,3.
$$ (eq-06-01-les)

The parameters $\vq=(q_1,q_2,q_3)^\T$ and
$\vbeta=(\beta_1,\beta_2,\beta_3)^\T$ are to be estimated from
four months of records:

```{list-table} Monthly budget $E_t$ (\$), prices $p_{ti}$ (\$/unit) and purchased quantities $Q_{ti}$ (units) in months $t=1,\dots,4$.
:name: tab-06-01-purchase-records
:header-rows: 1

* - Month $t$
  - 1
  - 2
  - 3
  - 4
* - $E_t$ (\$)
  - 1500
  - 2000
  - 1500
  - 1200
* - $p_{t1}$ (\$/unit)
  - 2.5
  - 2.7
  - 2.8
  - 2.6
* - $Q_{t1}$ (units)
  - 198
  - 200
  - 100
  - 152
* - $p_{t2}$ (\$/unit)
  - 4.0
  - 3.5
  - 4.0
  - 3.6
* - $Q_{t2}$ (units)
  - 124
  - 200
  - 80
  - 115
* - $p_{t3}$ (\$/unit)
  - 5.0
  - 4.9
  - 5.6
  - 4.5
* - $Q_{t3}$ (units)
  - 100
  - 110
  - 60
  - 65
```

Formulate the estimation problem as an NLP, solve it, and interpret the
estimates.
:::

::::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Decision variables.** The six model parameters
$\beta_1,\beta_2,\beta_3,q_1,q_2,q_3$. The data $E_t$, $p_{ti}$ and $Q_{ti}$
are constants of the optimization problem, even though they are "variables"
in the statistical sense.

**Objective.** The observed spending on product $i$ in month $t$ is
$p_{ti}Q_{ti}$, and {eq}`eq-06-01-les` predicts it from that month's budget
and prices. The residual is

$$
r_{ti}(\vbeta,\vq)=p_{ti}Q_{ti}-p_{ti}q_i-\beta_i\Bigl(E_t-\sum_{j=1}^3p_{tj}q_j\Bigr)
$$ (eq-06-01-residual)

for $t=1,\dots,4$ and $i=1,2,3$, in dollars. Least squares minimizes the sum of the twelve squared residuals,

$$
f(\vbeta,\vq)=\sum_{t=1}^4\sum_{i=1}^3r_{ti}(\vbeta,\vq)^2 .
$$ (eq-06-01-sse)

**Constraints.** The shares are nonnegative and sum to one; committed
quantities are nonnegative and cannot exceed what was actually bought in
any month, $q_i\le Q_{ti}$ for all $t$. With the residuals
{eq}`eq-06-01-residual`, the complete model is

$$
\begin{aligned}
\min_{\vbeta,\vq}\ & f(\vbeta,\vq)=\sum_{t=1}^4\sum_{i=1}^3r_{ti}(\vbeta,\vq)^2\\
\st\ & \beta_1+\beta_2+\beta_3=1,\\
& q_i\le Q_{ti}, && i=1,2,3,\ t=1,\dots,4,\\
& \vbeta\ge\vzero,\\
& \vq\ge\vzero .
\end{aligned}
$$ (eq-06-01-regression-nlp)

In the form {eq}`eq-06-01-nlp` it has $n=6$ variables, one equality constraint
$h_1=\beta_1+\beta_2+\beta_3-1$, and $12+6=18$ inequalities. Only the
tightest upper bound on each $q_i$ matters, $q_i\le\min_tQ_{ti}$, that is
$q_1\le100$, $q_2\le80$, $q_3\le60$ units; the other nine bounds are
redundant.

**Why the problem is nonlinear and not convex.** The residual contains the
products $\beta_iq_j$, so it is bilinear in $(\vbeta,\vq)$ and $f$ is a
polynomial of degree four. For fixed $\vbeta$ the residuals are affine in
$\vq$ and $f$ is a convex quadratic in $\vq$; for fixed $\vq$ it is a convex
quadratic in $\vbeta$. Jointly, however, $f$ is not convex: its Hessian at
the solution found below has a negative eigenvalue ($\approx-4.55$). The
constraints are linear, so this is a **linearly constrained** NLP with a
nonconvex objective, and a solver can only promise a local minimum.

**Solution.** A sequential quadratic programming method (SLSQP; Kraft, 1988; see
{doc}`06-06-constrained-algorithms` and the code in {doc}`09-06-code-nlp`) was
started from 1500 random feasible points. Every run that stopped at a point
satisfying the first-order optimality conditions of
{doc}`06-03-kkt-conditions` stopped at the same point, and no run found a
smaller objective value:

$$
\begin{aligned}
\vbeta^*&\approx(0.2474,\ 0.4635,\ 0.2891),\\
\vq^*&\approx(89.85,\ 16.67,\ 37.32)\ \text{units},\\
f^*&\approx1.3547\times10^{5}\ \$^2 .
\end{aligned}
$$

All bounds are inactive at $(\vbeta^*,\vq^*)$; only the equality
$\sum_i\beta_i=1$ binds. The Hessian of $f$ restricted to the plane
$\sum_i\beta_i=1$ is positive definite there, so the point is a strict local
minimum by the second-order conditions of {doc}`06-03-kkt-conditions`.

**Interpretation.** The household is estimated to commit about
$\sum_jp_{tj}q^*_j\approx\$462$–$\$527$ a month to the three products (the
amount varies with prices) and to spend the rest in the ratio
$25\,\%:46\,\%:29\,\%$. The fitted spending is compared with the
observations in {numref}`fig-06-01-regression-fit`. The root-mean-square
residual is $\sqrt{f^*/12}\approx\$106$, and the misfit is concentrated in
month 3, where all three residuals are between $-\$154$ and $-\$212$.

The reason is a feature of the model rather than of the solver. Summing
{eq}`eq-06-01-les` over $i$ and using $\sum_i\beta_i=1$ shows that the
predicted total spending in month $t$ is always exactly $E_t$. The
observed totals $\sum_ip_{ti}Q_{ti}$ are \$1491, \$1779, \$936 and
\$1101.70, so in month 3 the household spent only \$936 of its \$1500 budget
on these products, and no choice of parameters can reproduce that. If $E_t$
is replaced by the observed total, the same procedure gives
$f^*\approx1.241\times10^4\ \$^2$, an order of magnitude smaller, with
$\vbeta^*\approx(0.306,0.416,0.277)$ and $\vq^*\approx(27.08,\,0,\,19.57)$:
now the bound $q_2\ge0$ is binding. Choosing the data that enter a model is
part of formulating it.

```{figure} figures/06-01-regression-fit.svg
:name: fig-06-01-regression-fit
:width: 520px
:class: dark-light
:alt: Scatter plot of predicted against observed monthly spending for three products over four months, with the 45-degree line; the three month-3 points lie far above the line

Observed spending $p_{ti}Q_{ti}$ against fitted spending for the twelve
product–month pairs of {prf:ref}`ex-06-01-consumer-regression` at the
least-squares estimate. Points on the dashed line are fitted exactly. The
three points of month 3 (labeled) lie far above the line: the model
forces predicted spending to add up to the budget $E_3=\$1500$, while only
\$936 was spent.
```
::::

(sec-06-01-lp-vs-nlp)=
## Linear versus nonlinear programs

A linear program is the special case of {eq}`eq-06-01-nlp` in which $f$, every
$g_i$ and every $h_j$ are affine. Comparing the two classes explains why NLP
needs different theory and different algorithms.

```{list-table} Linear and nonlinear programs compared.
:name: tab-06-01-lp-vs-nlp
:header-rows: 1

* - Property
  - Linear program
  - Nonlinear program
* - Feasible set
  - a polyhedron (convex, finitely many vertices)
  - any set described by smooth inequalities and equations; may be
    nonconvex or disconnected
* - Where an optimum lies
  - if one exists, at a vertex (basic feasible solution)
  - anywhere: in the interior, on a curved boundary, or at a vertex
* - Local versus global
  - every local minimum is global
  - local minima need not be global unless the problem is convex
* - Existence
  - a bounded objective on a nonempty feasible set always attains its
    optimum
  - a bounded objective may fail to attain its infimum
* - Optimality test
  - reduced costs $\bar c_j\ge0$, checked exactly
  - Karush–Kuhn–Tucker (KKT) gradient conditions that are necessary but in general not
    sufficient
* - Algorithms
  - simplex method: finitely many pivots, exact arithmetic possible
  - iterative methods that converge in the limit, to a local minimum or
    stationary point
* - Duality
  - strong duality always holds when an optimum exists
  - holds for convex problems under a constraint qualification
```

The second row is the most visible difference. The fundamental theorem of
LP ({prf:ref}`thm-03-02-fundamental`) shows that a linear objective attains its minimum
over a polyhedron at a vertex, which is why the simplex method may restrict
its search to vertices. A nonlinear objective has curved contours, and the
smallest contour that touches the feasible set may touch it anywhere:
{prf:ref}`ex-06-01-three-feasible-sets` has two interior optima and one on a
curved boundary. NLP algorithms therefore move through the continuum of
feasible points, guided by derivatives.

(sec-06-01-classes)=
## Classes of nonlinear programs

Algorithms and theory are tailored to the structure of $f$ and $\Omega$.
The following classes recur throughout the book.

**Unconstrained problems.** $\Omega=\R^n$: no constraints at all, for
instance $\min x_1^2+x_2^2-1$ or problem 1 of
{prf:ref}`ex-06-01-three-feasible-sets`. Problems with constraints that are
inactive at the optimum behave locally like unconstrained ones. The
optimality conditions are the simplest (a zero gradient), and the methods of
{doc}`06-05-unconstrained-algorithms` apply.

**Constrained problems.** $\Omega\subsetneq\R^n$. Important subclasses are
the following.

- **Bound-constrained problems**: only $l_j\le x_j\le u_j$.
- **Linearly constrained problems**: $\Omega=\{\vx:\mA\vx\le\vb,\
  \mA'\vx=\vb'\}$ is a polyhedron and only $f$ is nonlinear, as in
  {prf:ref}`ex-06-01-consumer-regression`. The Frank–Wolfe method of
  {doc}`06-06-constrained-algorithms` exploits this structure by solving an LP
  in every iteration.
- **Quadratic programs** (QP): a quadratic objective over a polyhedron,

  $$
  \begin{aligned}
  \min_{\vx}\ & f(\vx)=\vc^\T\vx+\tfrac12\vx^\T\mQ\vx\\
  \st\ & \mA\vx=\vb,\\
  & \vx\ge\vzero,
  \end{aligned}
  $$ (eq-06-01-qp)

  with $\vc\in\R^n$ and a **symmetric** matrix $\mQ\in\Sym^n$. Symmetry
  costs nothing: since $\vx^\T\mQ\vx=\vx^\T\tfrac12(\mQ+\mQ^\T)\vx$, any
  matrix may be replaced by its symmetric part without changing $f$. With
  $\mQ$ symmetric, $\nabla f(\vx)=\vc+\mQ\vx$ and $\nabla^2f=\mQ$
  ({numref}`tab-01-02-vector-derivatives`), and the diagonal entry $Q_{jj}$ is twice the
  coefficient of $x_j^2$ while the off-diagonal entries $Q_{jk}=Q_{kj}$ both
  equal the coefficient of $x_jx_k$. The objective is convex exactly when
  $\mQ\succeq\vzero$, and the QP is then a convex program.
- **Linear programs**: $\mQ=\vzero$ in {eq}`eq-06-01-qp`. LP is the special
  case that Part III treats in {doc}`03-01-lp-model-building` to
  {doc}`03-07-sensitivity-analysis`.
- **Convex programs**: $f$ is a convex function and $\Omega$ is a convex
  set. Under {eq}`eq-06-01-nlp` this holds when $f$ and every $g_i$ are convex
  and every $h_j$ is **affine**. (Convex inequality functions give convex
  sublevel sets $\{g_i\le0\}$; a nonlinear equality $h_j(\vx)=0$ generally
  does not define a convex set, even when $h_j$ is convex, as the circle
  $x_1^2+x_2^2=1$ shows.) Convex programs are the subject of
  {doc}`02-03-convex-optimization-duality`.

To write a quadratic objective in QP form, read off the coefficients. For

$$
f(\vx)=x_1+3x_2+x_1^2-5x_1x_2+2x_2^2
$$

the linear part gives $\vc=(1,3)^\T$, and the quadratic part gives
$Q_{11}=2\cdot1=2$, $Q_{22}=2\cdot2=4$ and $Q_{12}=Q_{21}=-5$:

$$
\begin{aligned}
\mQ&=\begin{pmatrix}2&-5\\-5&4\end{pmatrix},\\
\tfrac12\vx^\T\mQ\vx&=x_1^2-5x_1x_2+2x_2^2 .
\end{aligned}
$$

Its eigenvalues are $\lambda_{1,2}(\mQ)=3\mp\sqrt{26}\approx-2.099,\ 8.099$,
so $\mQ$ is indefinite and this QP is not convex. *(An earlier version of
this example uses the non-symmetric matrix
$\bigl(\begin{smallmatrix}2&-4\\-2&2\end{smallmatrix}\bigr)$, which gives
$\tfrac12\vx^\T\mQ\vx=x_1^2-3x_1x_2+x_2^2$ and does not reproduce the
objective.)*

(sec-06-01-local-global)=
## Local and global minima

:::{prf:definition} Global and local minima
:label: def-06-01-local-global

Let $\Omega\subseteq\R^n$ and $f:\Omega\to\R$, and let $\vx^*\in\Omega$.

1. $\vx^*$ is a **global minimum** of $f$ over $\Omega$ if
   $f(\vx)\ge f(\vx^*)$ for all $\vx\in\Omega$.
2. $\vx^*$ is a **local minimum** if there is an $\varepsilon>0$ such that
   $f(\vx)\ge f(\vx^*)$ for all $\vx\in\Omega$ with
   $\norm{\vx-\vx^*}<\varepsilon$, that is, for all
   $\vx\in\Omega\cap B(\vx^*,\varepsilon)$.
3. A local (global) minimum is **strict** if the inequality is strict for
   all such $\vx\ne\vx^*$.

Maxima are defined with the inequalities reversed.
:::

Every global minimum is a local minimum; the converse fails in general.
A local minimum only has to beat its feasible neighbors, and the radius
$\varepsilon$ may be tiny. Note that a local minimum is defined relative to
$\Omega$: a point on the boundary of $\Omega$ can be a local minimum even
though $f$ decreases if one leaves $\Omega$.

A one-dimensional function shows all the cases at once. Consider

$$
f(x)=(x^2)^x=e^{2x\ln\abs{x}},\qquad -5\le x\le1.7,
$$ (eq-06-01-x-power-2x)

with $f(0)=1$ (the limit of $f$ as $x\to0$). For $x\ne0$,
$f'(x)=2f(x)(\ln\abs{x}+1)$, which vanishes only at $x=\pm1/e$.

- On $[-5,-1/e)$ we have $\ln\abs{x}>-1$, so $f'>0$ and $f$ increases;
  on $(-1/e,0)$ and $(0,1/e)$ it decreases; on $(1/e,1.7]$ it increases.
- $x=-1/e\approx-0.368$ is a **local maximum** with
  $f=e^{2/e}\approx2.087$.
- $x=1/e\approx0.368$ is a **local minimum** with $f=e^{-2/e}\approx0.4791$,
  but not a global one.
- The **global minimum** is the left end point $x=-5$, with
  $f(-5)=25^{-5}=1.024\times10^{-7}$; the global maximum is the right end
  point, $f(1.7)\approx6.075$.

Two lessons follow. First, the global minimum is not a stationary point:
$f'(-5)\approx5.3\times10^{-7}>0$, and $x=-5$ is a minimum only because the
constraint $x\ge-5$ stops further decrease. Optimality conditions for
constrained problems must therefore allow a nonzero gradient
({doc}`06-02-optimality-conditions`). Second, an algorithm that only looks at
local information and is started near $x=0.5$ will stop at the local
minimum $1/e$. {numref}`fig-06-01-x-power-2x` also marks where $f$ is convex
($f''>0$ on $[-5,-0.808)$ and on $(0,1.7]$) and concave
($f''<0$ on $(-0.808,0)$); local minima that are not global can occur only
because $f$ is not convex on the whole interval.
*(The switch from convex to concave is at $x\approx-0.808$, not at
$x\approx-1.1$ as it is sometimes drawn.)*

```{figure} figures/06-01-x-power-2x.svg
:name: fig-06-01-x-power-2x
:width: 560px
:class: dark-light
:alt: Graph of f(x) = (x squared) to the power x on the interval from -5 to 1.7, with the global minimum at x = -5, a local maximum at -1/e, a local minimum at 1/e and the global maximum at 1.7; shaded bands show convex and concave parts

The function {eq}`eq-06-01-x-power-2x` on $[-5,1.7]$. The local minimum
$1/e$ is not global; the global minimum and maximum are at the end points
of the feasible interval. Shaded: the interval $(-0.808,0)$ where $f$ is
concave.
```

The difference between local and global minima disappears for convex
programs.

:::{prf:theorem} Local minima of convex programs are global
:label: thm-06-01-convex-local-global

If $f$ is convex on the convex set $\Omega$, every local minimum of $f$ over
$\Omega$ is a global minimum, and the set of global minima is convex. If
$f$ is strictly convex, there is at most one minimum.
:::

The proof is short: if a local minimum $\vx^*$ had a better point
$\vy\in\Omega$, the points $\vx^*+\theta(\vy-\vx^*)$, $\theta\in(0,1]$,
would be feasible by convexity of $\Omega$ and better than $\vx^*$ by
convexity of $f$, including points arbitrarily close to $\vx^*$. It is given
in full as {prf:ref}`thm-02-03-local-global`, and a version based on
gradients is {prf:ref}`thm-06-02-convex-sufficient`. The theorem is the
main reason why recognizing convexity ({doc}`02-02-convex-functions`) is worth
the effort: for a convex program, any method that finds a local minimum
finds a global one.

(sec-06-01-existence)=
## Existence of a minimum

A problem may have no optimal solution even when its objective is bounded
(an example is in {prf:ref}`ex-08-08-difficulties`). Two sufficient conditions
cover most models.

:::{prf:theorem} Existence of minima
:label: thm-06-01-existence

Let $\Omega\subseteq\R^n$ be nonempty and $f$ continuous on $\Omega$.

1. (Weierstrass) If $\Omega$ is closed and bounded, $f$ attains a global
   minimum and a global maximum on $\Omega$.
2. If $\Omega$ is closed and $f$ is **coercive**, that is
   $f(\vx)\to\infty$ whenever $\norm{\vx}\to\infty$ with $\vx\in\Omega$,
   then $f$ attains a global minimum on $\Omega$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Part 1 is Weierstrass' theorem ({prf:ref}`thm-02-01-weierstrass`). For part
2, fix any $\vx^0\in\Omega$. By coercivity there is a radius $R$ with
$f(\vx)>f(\vx^0)$ for every $\vx\in\Omega$ with $\norm{\vx}>R$. The set
$\Omega_R=\{\vx\in\Omega:\norm{\vx}\le R\}$ is closed, bounded and contains
$\vx^0$, so by part 1 $f$ has a minimum $\vx^*$ on $\Omega_R$. Then
$f(\vx^*)\le f(\vx^0)<f(\vx)$ for all feasible $\vx$ outside $\Omega_R$,
so $\vx^*$ is a global minimum over $\Omega$.
:::

Constraints defined by continuous functions with "$\le$" and "$=$" give
closed sets, so closedness is usually automatic for {eq}`eq-06-01-nlp`; strict
inequalities such as $x_1+x_2<2$ do not. In
{prf:ref}`ex-06-01-three-feasible-sets`, problems 2 and 3 have compact feasible
sets, and in problem 1 the objective is coercive, so all three have optimal
solutions before any calculation. In {prf:ref}`ex-06-01-consumer-regression`
the feasible set is compact ($0\le\beta_i\le1$, $0\le q_i\le\min_tQ_{ti}$),
so a global minimum exists; the difficulty is only to find it.

(sec-06-01-difficulties)=
## Difficulties of nonlinear programming

The comparison with LP in {numref}`tab-06-01-lp-vs-nlp` can be summarized as a
list of things that go wrong, each of which shapes the methods in later
chapters.

1. **A local minimum need not be global.** Algorithms based on derivatives
   find stationary points and local minima; certifying global optimality
   requires convexity or an exhaustive search
   ({numref}`fig-06-01-x-power-2x`, {prf:ref}`ex-06-01-consumer-regression`).
2. **The optimum need not be at a vertex.** There may be no vertices at
   all, and enumeration of finitely many candidate points, the idea behind
   the simplex method, is not available
   ({prf:ref}`ex-06-01-three-feasible-sets`).
3. **Algorithms are iterative and numerical.** They produce a sequence of
   points that converges only in the limit and must be stopped by a
   tolerance. Poor scaling makes this delicate: in
   {prf:ref}`ex-06-01-consumer-regression` the partial derivatives of $f$ with
   respect to the shares are of order $10^5$ while those with respect to the
   quantities are of order one, and many solver runs stop early. Rounding
   errors and ill-conditioning are discussed with the algorithms in
   {doc}`06-04-line-search` and {doc}`06-05-unconstrained-algorithms`.
4. **A bounded objective may have no optimum.** The infimum of $f$ may be
   approached but not attained when $\Omega$ is unbounded or not closed
   ({prf:ref}`ex-08-08-difficulties`); {prf:ref}`thm-06-01-existence` gives
   conditions that exclude this.
5. **The feasible set may be nonconvex or even disconnected.** A
   constraint such as $g_1(x)\ge g_2(x)$ between two curves can split
   $\Omega$ into separate pieces ({prf:ref}`ex-08-08-difficulties`), and a
   search that moves continuously inside $\Omega$ cannot pass from one
   piece to another.
6. **Optimality conditions are only necessary.** The conditions of
   {doc}`06-02-optimality-conditions` and {doc}`06-03-kkt-conditions` single out
   candidate points, but in a nonconvex problem a candidate can be a saddle
   point or a local maximum, and in degenerate cases a minimum can even fail
   to satisfy the first-order conditions unless a constraint qualification
   holds.

The convex case removes difficulties 1, 5 and most of 6, which is why the
next two chapters state each result first in general and then for convex
programs.

## Summary and Common Mistakes

### Summary

- An NLP minimizes $f(\vx)$ subject to $g_i(\vx)\le0$ and $h_j(\vx)=0$
  ({eq}`eq-06-01-nlp`); maximization and "$\ge$" constraints are converted by
  changing signs.
- Least-squares estimation with side conditions is a typical NLP
  ({prf:ref}`ex-06-01-consumer-regression`); its objective there is nonconvex,
  its solution was found numerically, and the fit exposed a feature of the
  model ({numref}`fig-06-01-regression-fit`).
- Unlike LP, the optimum of an NLP may lie in the interior or on a curved
  boundary of $\Omega$ ({prf:ref}`ex-06-01-three-feasible-sets`,
  {numref}`tab-06-01-lp-vs-nlp`).
- Main classes: unconstrained, linearly constrained, quadratic
  ({eq}`eq-06-01-qp`, with symmetric $\mQ$), linear and convex programs (convex
  $f$ and $g_i$, affine $h_j$).
- Local and global minima are distinguished in
  {prf:ref}`def-06-01-local-global`; for convex programs every local minimum is
  global ({prf:ref}`thm-06-01-convex-local-global`).
- A minimum exists for a continuous function on a nonempty compact set, or
  for a coercive function on a closed set ({prf:ref}`thm-06-01-existence`).
- Difficulties of NLP: local optima, non-vertex optima, iterative numerics,
  unattained infima, nonconvex or disconnected feasible sets, and optimality
  conditions that are only necessary.

### Common mistakes

- **Expecting a vertex solution.** Checking only the corners of $\Omega$
  misses interior and curved-boundary optima such as those in
  {prf:ref}`ex-06-01-three-feasible-sets`.
- **Forgetting the factor $\tfrac12$ or the symmetry of $\mQ$.** In
  $\vc^\T\vx+\tfrac12\vx^\T\mQ\vx$ the diagonal of $\mQ$ holds twice the
  coefficients of $x_j^2$ and each off-diagonal pair holds the full
  coefficient of $x_jx_k$.
- **Calling a problem with nonlinear equality constraints convex.** Even
  with convex $f$ and $g_i$, the equalities must be affine.
- **Reporting the first local minimum a solver returns as "the" optimum** of
  a nonconvex problem. Use several starting points, or prove convexity.
- **Treating a supremum as a maximum.** $\sup f$ exists for every bounded
  function; a maximizer exists only under conditions such as
  {prf:ref}`thm-06-01-existence`.


## References

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 1,
  introduction and examples of NLP models; Section 2.3, Weierstrass'
  theorem; Section 3.4, minima of convex functions.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Sections
  4.1–4.4, optimization problems, convex problems, quadratic programs.)
- Kraft, D. (1988). *A software package for sequential quadratic
  programming* (Report DFVLR-FB 88-28). Deutsche Forschungs- und
  Versuchsanstalt für Luft- und Raumfahrt.
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming*
  (4th ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3> (Chapter 7,
  local and global minima, first-order necessary conditions and examples of
  unconstrained problems.)
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Chapter 1,
  classes of optimization problems; Section 2.1, local and global
  solutions.)
- Sheffi, Y. (1985). *Urban transportation networks: Equilibrium analysis
  with mathematical programming methods*. Prentice-Hall. (Chapter 2, basic
  concepts in minimization problems.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 20, classical optimization theory.)
