# Convex Functions

## Overview

An optimization problem is easy, in a precise sense, when its objective
and its inequality constraints are convex functions: then every local
minimum is a global one and the minima can be found by efficient
algorithms ({doc}`02-03-convex-optimization-duality`). This chapter defines
convex, strictly convex and concave functions
({prf:ref}`def-02-02-convex-function`), relates them to the convex sets of
{doc}`02-01-convex-sets` through the epigraph
({prf:ref}`thm-02-02-epigraph`), and derives the two practical tests of
convexity: the first-order condition, which says that a convex function
lies above its tangent planes ({prf:ref}`thm-02-02-first-order`), and the
second-order condition, which says that its Hessian is positive
semidefinite ({prf:ref}`thm-02-02-second-order`). It then lists the standard
convex functions, the operations that preserve convexity, the convexity of
sublevel sets, Jensen's inequality and, briefly, quasiconvex functions.
It uses the definiteness tests of {doc}`01-01-linear-algebra` and the
gradients, Hessians and Taylor's theorem of {doc}`01-02-matrix-calculus`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- state the definitions of convex, strictly convex, concave and strictly
  concave functions and explain why the domain must be a convex set;
- prove that a function is convex from the definition, through its
  epigraph, or with the first- or second-order condition;
- test the convexity of a quadratic function through the eigenvalues or
  the principal minors of its Hessian;
- build convex functions from simple ones with operations that preserve
  convexity, and recognize the standard examples (norms, exponentials,
  logarithms, maxima, log-sum-exp);
- use the convexity of sublevel sets to show that a feasible set is convex;
  and
- apply Jensen's inequality and recognize quasiconvex functions.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\dom f$ | domain of $f$ (the set on which $f$ is finite) |
| $\theta$, $\theta_i$ | weights of a convex combination, $\theta\in[0,1]$ |
| $\tilde f$ | extended-value extension of $f$, equal to $+\infty$ outside $\dom f$ |
| $\iota_C$ | indicator function of a set $C$: $0$ on $C$, $+\infty$ outside |
| $\operatorname{epi} f$ | epigraph $\{(\vx,t):\vx\in\dom f,\ f(\vx)\le t\}$ |
| $S_\gamma$ | sublevel set $\{\vx\in\dom f: f(\vx)\le\gamma\}$ |
| $\operatorname{lse}(\vx)$ | log-sum-exp function $\log\sum_{i}e^{x_i}$ |

Logarithms are natural logarithms.

(sec-02-02-definition)=
## Convex and concave functions

:::{prf:definition} Convex and concave functions
:label: def-02-02-convex-function

Let $C\subseteq\R^n$ be a nonempty convex set and $f:C\to\R$. The
function $f$ is **convex** on $C$ if for all $\vx,\vy\in C$ and all
$\theta\in[0,1]$

$$
f\bigl(\theta\vx+(1-\theta)\vy\bigr)\le\theta f(\vx)+(1-\theta)f(\vy).
$$ (eq-02-02-convex)

It is

- **strictly convex** on $C$ if the inequality {eq}`eq-02-02-convex` is
  strict whenever $\vx\ne\vy$ and $\theta\in(0,1)$;
- **concave** (**strictly concave**) on $C$ if $-f$ is convex
  (strictly convex), that is, if {eq}`eq-02-02-convex` holds with $\ge$
  (with $>$ for $\vx\ne\vy$, $\theta\in(0,1)$).
:::

*(Strictness can only be required for distinct points and for
$\theta\in(0,1)$: for $\theta=0$, $\theta=1$ or $\vx=\vy$ both sides of
{eq}`eq-02-02-convex` are equal, so a definition that asks for $<$ for all
$\theta\in[0,1]$ would make no function strictly convex.)*

The right-hand side of {eq}`eq-02-02-convex` is the height of the **chord**
from $(\vx,f(\vx))$ to $(\vy,f(\vy))$ above the point
$\theta\vx+(1-\theta)\vy$. Convexity therefore means that every chord lies
on or above the graph ({numref}`fig-02-02-chord-tangent`(a)); strict
convexity means that it lies strictly above the graph between its end
points. A concave function has every chord on or below its graph.

The domain must be convex because {eq}`eq-02-02-convex` evaluates $f$ at the
point $\theta\vx+(1-\theta)\vy$, which must lie in the domain for every
pair of points and every $\theta$. A statement such as "$1/x$ is convex"
is therefore incomplete: $1/x$ is convex on $(0,\infty)$ and concave on
$(-\infty,0)$, and on its whole domain $\R\setminus\{0\}$, which is not
convex, the definition does not apply.

```{figure} figures/02-02-chord-tangent.svg
:name: fig-02-02-chord-tangent
:width: 640px
:class: dark-light
:alt: Left, a convex curve with a chord between two points x and y lying above the curve. Right, the same curve with its tangent line at x lying below the curve, and the gap at a point y

The two characterizations of a convex function of one variable. (a) The
chord between $(x,f(x))$ and $(y,f(y))$ lies above the graph: at
$\theta x+(1-\theta)y$ its height $\theta f(x)+(1-\theta)f(y)$ is at
least $f(\theta x+(1-\theta)y)$. (b) The tangent at $x$ lies below the
graph: $f(y)\ge f(x)+f'(x)(y-x)$ for every $y$.
```

Three questions can be answered from the definition alone.

*Is every linear function convex?* Yes. For an affine function
$f(\vx)=\va^\T\vx+b$, both sides of {eq}`eq-02-02-convex` are equal, because
$\va^\T(\theta\vx+(1-\theta)\vy)+b=\theta(\va^\T\vx+b)+(1-\theta)(\va^\T\vy+b)$.
An affine function is therefore both convex and concave, but neither
strictly convex nor strictly concave. Conversely, a function that is both
convex and concave on $\R^n$ is affine.

*Is every parabola convex?* No. For $f(x)=ax^2+bx+c$ a short calculation
gives, for all $x,y$ and $\theta$,

$$
\theta f(x)+(1-\theta)f(y)-f\bigl(\theta x+(1-\theta)y\bigr)
=a\,\theta(1-\theta)(x-y)^2 .
$$ (eq-02-02-parabola-gap)

The linear and constant terms cancel, as for any affine function, and for
the quadratic term
$\theta x^2+(1-\theta)y^2-(\theta x+(1-\theta)y)^2=\theta(1-\theta)(x-y)^2$.
The gap {eq}`eq-02-02-parabola-gap` is nonnegative for all $x,y,\theta$
exactly when $a\ge0$. So a parabola is convex if and only if its leading
coefficient is nonnegative, strictly convex if and only if $a>0$, and
concave if and only if $a\le0$. The parabola $-x^2$ is strictly concave.

*How can convexity be verified in general?* Checking
{eq}`eq-02-02-convex` for all pairs of points is rarely practical. The
sections below give three better tools: the derivative tests
({ref}`sec-02-02-first-order`, {ref}`sec-02-02-second-order`), a catalog of
known convex functions ({ref}`sec-02-02-examples`) and rules that combine
them ({ref}`sec-02-02-operations`).

A useful reduction to one variable is the following. A function $f$ is
convex if and only if its restriction to every line,
$\varphi(t)=f(\vx+t\vv)$ on $\{t:\vx+t\vv\in\dom f\}$, is convex for all
$\vx\in\dom f$ and all directions $\vv$. This holds because a convex
combination of two points on a line is a point of the same segment, so
{eq}`eq-02-02-convex` for $f$ and for all restrictions $\varphi$ are the same
statements. Many proofs below are done for $n=1$ and then lifted to
$\R^n$ in this way.

(sec-02-02-epigraph)=
## Extended values and the epigraph

It is often convenient to define a convex function on all of $\R^n$ by
giving it the value $+\infty$ outside its domain. The **extended-value
extension** of a convex $f$ is

$$
\tilde f(\vx)=
\begin{cases}
f(\vx), & \vx\in\dom f,\\
+\infty, & \vx\notin\dom f .
\end{cases}
$$

With the arithmetic rules $\theta\cdot\infty=\infty$ for $\theta>0$ and
$a\le\infty$ for every $a$, inequality {eq}`eq-02-02-convex` holds for
$\tilde f$ at all $\vx,\vy\in\R^n$ and $\theta\in(0,1)$: if $\vx$ or $\vy$
lies outside the domain, the right-hand side is $+\infty$. The domain is
recovered as $\dom f=\{\vx:\tilde f(\vx)<\infty\}$. A special case is the
**indicator function** $\iota_C$ of a convex set $C$, which is $0$ on $C$
and $+\infty$ outside. Minimizing $f$ over $C$ is the same as minimizing
$f+\iota_C$ over $\R^n$, so constraints can be written into the objective.
From now on the extension is used without comment.

The link between convex functions and convex sets is the set of points on
or above the graph.

:::{prf:definition} Epigraph
:label: def-02-02-epigraph

The **epigraph** of $f:\dom f\to\R$ is the subset of $\R^{n+1}$

$$
\operatorname{epi} f=\{(\vx,t):\vx\in\dom f,\ t\in\R,\ f(\vx)\le t\}.
$$
:::

:::{prf:theorem} Convex functions have convex epigraphs
:label: thm-02-02-epigraph

A function $f$ with a convex domain is convex if and only if
$\operatorname{epi} f$ is a convex set.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $f$ be convex, $(\vx,s),(\vy,t)\in\operatorname{epi} f$ and
$\theta\in[0,1]$. Then
$f(\theta\vx+(1-\theta)\vy)\le\theta f(\vx)+(1-\theta)f(\vy)\le\theta
s+(1-\theta)t$, so the combination
$(\theta\vx+(1-\theta)\vy,\ \theta s+(1-\theta)t)$ lies in the epigraph.

Conversely, let the epigraph be convex and $\vx,\vy\in\dom f$. The points
$(\vx,f(\vx))$ and $(\vy,f(\vy))$ belong to it, hence so does
$(\theta\vx+(1-\theta)\vy,\ \theta f(\vx)+(1-\theta)f(\vy))$, which is
{eq}`eq-02-02-convex`.
:::

The theorem turns statements about functions into statements about sets.
For instance, a pointwise maximum of convex functions is convex because
its epigraph is the intersection of their epigraphs
({prf:ref}`thm-02-01-intersection`), and the tangent-plane property of the
next section is the supporting hyperplane theorem
({prf:ref}`thm-02-01-supporting-hyperplane`) applied to the epigraph.
{numref}`fig-02-02-epigraph-sublevel`(a) shows an epigraph.

(sec-02-02-first-order)=
## The first-order condition

From here on, a differentiable function is one whose gradient exists on an
open set containing its (convex) domain; {doc}`01-02-matrix-calculus`
explains the gradient $\nabla f$ and the Hessian $\nabla^2f$.

:::{prf:theorem} First-order condition
:label: thm-02-02-first-order

Let $C$ be convex and $f$ differentiable on $C$. Then $f$ is convex on $C$
if and only if

$$
f(\vy)\ge f(\vx)+\nabla f(\vx)^\T(\vy-\vx)\qquad\text{for all }\vx,\vy\in C .
$$ (eq-02-02-first-order)

$f$ is strictly convex on $C$ if and only if {eq}`eq-02-02-first-order` holds
with strict inequality for all $\vx\ne\vy$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

*Convex ⇒ {eq}`eq-02-02-first-order`.* Let $\vx,\vy\in C$ and
$\theta\in(0,1]$. By convexity, with the point $\vx+\theta(\vy-\vx)
=\theta\vy+(1-\theta)\vx$,

$$
f\bigl(\vx+\theta(\vy-\vx)\bigr)\le f(\vx)+\theta\bigl(f(\vy)-f(\vx)\bigr).
$$

Subtracting $f(\vx)$ and dividing by $\theta$,

$$
\frac{f\bigl(\vx+\theta(\vy-\vx)\bigr)-f(\vx)}{\theta}\le f(\vy)-f(\vx).
$$

As $\theta\downarrow0$ the left-hand side tends to the directional
derivative $\nabla f(\vx)^\T(\vy-\vx)$, which gives {eq}`eq-02-02-first-order`.

*{eq}`eq-02-02-first-order` ⇒ convex.* Let $\vz=\theta\vx+(1-\theta)\vy$.
Apply {eq}`eq-02-02-first-order` at $\vz$ twice:

$$
\begin{aligned}
f(\vx)&\ge f(\vz)+\nabla f(\vz)^\T(\vx-\vz),\\
f(\vy)&\ge f(\vz)+\nabla f(\vz)^\T(\vy-\vz).
\end{aligned}
$$

Multiply the first inequality by $\theta$, the second by $1-\theta$ and
add. The gradient terms cancel because
$\theta(\vx-\vz)+(1-\theta)(\vy-\vz)=\vzero$, and what remains is
$\theta f(\vx)+(1-\theta)f(\vy)\ge f(\vz)$.

*Strict version.* If the strict inequality holds for $\vx\ne\vy$, the same
argument with $\theta\in(0,1)$ gives strict convexity. Conversely, let $f$
be strictly convex, $\vx\ne\vy$ and $\vz=\tfrac12(\vx+\vy)$. By the
nonstrict condition at $\vx$ and by strict convexity,

$$
\begin{aligned}
&f(\vx)+\tfrac12\nabla f(\vx)^\T(\vy-\vx)\\
&=f(\vx)+\nabla f(\vx)^\T(\vz-\vx)\\
&\le f(\vz)\\
&<\tfrac12 f(\vx)+\tfrac12 f(\vy),
\end{aligned}
$$

and multiplying by $2$ and rearranging gives
$f(\vy)>f(\vx)+\nabla f(\vx)^\T(\vy-\vx)$.
:::

The right-hand side of {eq}`eq-02-02-first-order` is the first-order Taylor
approximation of $f$ at $\vx$, an affine function of $\vy$. The theorem
says that for a convex function this local, linear information is a
*global* underestimator ({numref}`fig-02-02-chord-tangent`(b)). In terms of the
epigraph, the hyperplane $t=f(\vx)+\nabla f(\vx)^\T(\vy-\vx)$ in
$\R^{n+1}$ supports $\operatorname{epi} f$ at $(\vx,f(\vx))$.

The most important consequence concerns stationary points. If
$\nabla f(\vx^*)=\vzero$, then {eq}`eq-02-02-first-order` gives
$f(\vy)\ge f(\vx^*)$ for all $\vy\in C$: **a stationary point of a convex
function is a global minimizer.** For a nonconvex function a stationary
point may be a local minimum, a local maximum or a saddle point
({doc}`06-02-optimality-conditions`). The constrained version of this
statement is the optimality condition of
{prf:ref}`thm-02-03-optimality`.

Applied twice, with the roles of $\vx$ and $\vy$ exchanged, and added,
{eq}`eq-02-02-first-order` gives

$$
\bigl(\nabla f(\vy)-\nabla f(\vx)\bigr)^\T(\vy-\vx)\ge0
\qquad\text{for all }\vx,\vy\in C ,
$$ (eq-02-02-monotone-gradient)

so the gradient of a convex function is **monotone**: for $n=1$, the
derivative $f'$ is nondecreasing.

(sec-02-02-second-order)=
## The second-order condition

For twice differentiable functions the curvature decides convexity.

:::{prf:theorem} Second-order condition
:label: thm-02-02-second-order

Let $C$ be a convex set and $f$ twice continuously differentiable on an
open set containing $C$.

1. If $\nabla^2f(\vx)\succeq0$ for all $\vx\in C$, then $f$ is convex on
   $C$.
2. If $\nabla^2f(\vx)\succ0$ for all $\vx\in C$, then $f$ is strictly
   convex on $C$.
3. If $C$ is open and $f$ is convex on $C$, then $\nabla^2f(\vx)\succeq0$
   for all $\vx\in C$.

In particular, on an open convex set (for example $\R^n$), $f$ is convex if
and only if its Hessian is positive semidefinite everywhere.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

*Parts 1 and 2.* Let $\vx,\vy\in C$, $\vx\ne\vy$. By Taylor's theorem with
the mean-value remainder ({prf:ref}`thm-01-02-taylor`) there is a point
$\vz$ on the segment $[\vx,\vy]\subseteq C$ with

$$
f(\vy)=f(\vx)+\nabla f(\vx)^\T(\vy-\vx)
+\tfrac12(\vy-\vx)^\T\nabla^2f(\vz)(\vy-\vx).
$$

If $\nabla^2f(\vz)\succeq0$ the last term is nonnegative, so
{eq}`eq-02-02-first-order` holds and $f$ is convex by
{prf:ref}`thm-02-02-first-order`. If $\nabla^2f(\vz)\succ0$ the last term is
positive and $f$ is strictly convex.

*Part 3.* Suppose $\vv^\T\nabla^2f(\vx)\vv<0$ for some $\vx\in C$ and some
$\vv$. Because $C$ is open, $\vx+t\vv\in C$ for all small $t>0$, and by
Taylor's theorem and continuity of the Hessian

$$
\begin{aligned}
f(\vx+t\vv)&=f(\vx)+t\nabla f(\vx)^\T\vv+\tfrac{t^2}{2}\vv^\T\nabla^2f(\vz_t)\vv\\
&<f(\vx)+t\nabla f(\vx)^\T\vv
\end{aligned}
$$

for $t$ small enough, since $\vz_t\to\vx$. This contradicts
{eq}`eq-02-02-first-order`, so $\nabla^2f(\vx)\succeq0$.
:::

Two remarks prevent common misreadings of the theorem.

:::{prf:remark} Positive definiteness is not necessary for strict convexity
:label: rem-02-02-strict-hessian

Part 2 has no converse. The function $f(x)=x^4$ is strictly convex on $\R$
(it satisfies {eq}`eq-02-02-first-order` strictly, because
$y^4-x^4-4x^3(y-x)=(y-x)^2(y^2+2xy+3x^2)>0$ for $y\ne x$), yet
$f''(0)=0$. A strictly convex function can have a singular Hessian at
isolated points.
:::

:::{prf:remark} The Hessian test needs an open set
:label: rem-02-02-hessian-open

Part 3 needs an open domain, or at least one with nonempty interior. On a
set without interior the Hessian test is only sufficient. For example,
$f(\vx)=x_1^2-x_2^2$ is convex on the line $C=\{\vx:x_2=0\}$, where it
equals $x_1^2$, although $\nabla^2f=\diag(2,-2)$ is indefinite everywhere.
*(The test is often stated as "convex on $\Omega$ if and only if the
Hessian is positive semidefinite on $\Omega$" without this condition; only
the "if" direction holds in general.)*
:::

For a function of one variable the theorem reads: $f$ is convex on an
open interval if and only if $f''\ge0$ there. This is the familiar
statement that the graph "curves upward".

### Testing the Hessian

Deciding whether $\nabla^2f(\vx)\succeq0$ is a question about a symmetric
matrix, answered in {doc}`01-01-linear-algebra` ({prf:ref}`thm-01-01-eigen-test`,
{prf:ref}`thm-01-01-sylvester`, {prf:ref}`thm-01-01-psd-minors`). For a symmetric $\mH$:

- $\mH\succeq0$ if and only if all eigenvalues $\lambda_i(\mH)\ge0$, and
  $\mH\succ0$ if and only if all $\lambda_i(\mH)>0$;
- $\mH\succ0$ if and only if all **leading** principal minors are
  positive (Sylvester's criterion);
- $\mH\succeq0$ if and only if **all** principal minors (the determinants
  of the submatrices obtained by deleting the same rows and columns) are
  nonnegative.

The last two tests are not interchangeable. Nonnegative *leading*
principal minors do not imply positive semidefiniteness:
$\mH=\diag(0,-1)$ has leading minors $0$ and $0$, but the eigenvalue $-1$.
Both minor tests also assume that $\mH$ is symmetric, which a Hessian of a
twice continuously differentiable function always is. *(A version of the
minor test that omits symmetry or uses only leading minors for
semidefiniteness is incorrect.)*

:::{prf:example} A convex quadratic that is not strictly convex
:label: ex-02-02-hessian-test

Is $f(\vx)=x_1^2-2x_1x_2+x_2^2$ convex on $\R^2$? Is it strictly convex?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The gradient and Hessian are

$$
\begin{aligned}
\nabla f(\vx)&=\begin{pmatrix}2x_1-2x_2\\-2x_1+2x_2\end{pmatrix},\\
\nabla^2f(\vx)&=\begin{pmatrix}2&-2\\-2&2\end{pmatrix}.
\end{aligned}
$$

The Hessian is constant. Its characteristic polynomial is
$(2-\lambda)^2-4=\lambda(\lambda-4)$, so its eigenvalues are $0$ and $4$.
Equivalently, its principal minors are $2$, $2$ (the diagonal entries) and
$\det\nabla^2f=4-4=0$, all nonnegative. Hence $\nabla^2f\succeq0$ on the
open set $\R^2$, and by {prf:ref}`thm-02-02-second-order` **$f$ is convex on
$\R^2$**.

The Hessian is not positive definite (it has the eigenvalue $0$), so part 2
of the theorem does not apply, and by {prf:ref}`rem-02-02-strict-hessian` this
alone would not decide the question. Here the answer is negative: the
eigenvector of the eigenvalue $0$ is $(1,1)$, and along this direction $f$
does not curve at all. Indeed $f(\vx)=(x_1-x_2)^2$ is zero on the whole
line $x_1=x_2$. For $\vx=(0,0)$, $\vy=(1,1)$ and $\theta=\tfrac12$ both
sides of {eq}`eq-02-02-convex` equal $0$, so the strict inequality fails and
**$f$ is convex but not strictly convex**. It is not concave either,
because of the eigenvalue $4$. Its minimizers form the whole line
$x_1=x_2$, which illustrates that a convex function can have infinitely
many global minimizers.
:::

(sec-02-02-examples)=
## Examples of convex functions

The following functions are the building blocks of most convex models.
Each claim can be checked with {prf:ref}`thm-02-02-second-order` or with the
definition.

**Affine functions** $\va^\T\vx+b$ are convex and concave.

**Quadratic functions** $f(\vx)=\tfrac12\vx^\T\mQ\vx+\vc^\T\vx+d$ with
symmetric $\mQ$ have $\nabla^2f=\mQ$. They are convex if and only if
$\mQ\succeq0$, and strictly convex if and only if $\mQ\succ0$ (for a
quadratic, $\vv^\T\mQ\vv=0$ with $\vv\ne\vzero$ makes $f$ affine along
$\vv$, as in {prf:ref}`ex-02-02-hessian-test`). Least-squares objectives
$\norm{\mA\vx-\vb}^2$, with Hessian $2\mA^\T\mA\succeq0$, are convex.

**Norms.** Every norm is convex: by the triangle inequality and
homogeneity,
$\norm{\theta\vx+(1-\theta)\vy}\le\theta\norm{\vx}+(1-\theta)\norm{\vy}$.
No norm is strictly convex, because along a ray $\norm{t\vx}=t\norm{\vx}$
is linear. The squared Euclidean norm $\norm{\vx}^2=\vx^\T\vx$ is strictly
convex ($\mQ=2\mI$). Norms are not differentiable at $\vzero$, which shows
that convexity does not require differentiability.

**Functions of one variable.**

- $e^{ax}$ is convex on $\R$ for every $a\in\R$ (second derivative
  $a^2e^{ax}\ge0$).
- $x^p$ on $(0,\infty)$ is convex for $p\ge1$ or $p\le0$ and concave for
  $0\le p\le1$ (second derivative $p(p-1)x^{p-2}$).
- $\log x$ is concave on $(0,\infty)$ and $-\log x$ is convex (second
  derivative of $\log x$ is $-1/x^2$).
- $x\log x$ is convex on $(0,\infty)$ (second derivative $1/x$).
- $\abs{x}^p$ is convex on $\R$ for $p\ge1$.

**Maximum of affine functions.** $f(\vx)=\max_{i=1,\dots,k}(\va_i^\T\vx+b_i)$
is convex, since its epigraph is the intersection of the $k$ half-spaces
$\{(\vx,t):\va_i^\T\vx+b_i\le t\}$. Such piecewise-linear convex functions
model costs with increasing marginal rates and can be minimized by linear programming
({ref}`sec-07-03-convex-case`).

**Log-sum-exp.** The function
$\operatorname{lse}(\vx)=\log\bigl(e^{x_1}+\dots+e^{x_n}\bigr)$ is a smooth
approximation of the maximum:

$$
\max_i x_i\le\operatorname{lse}(\vx)\le\max_i x_i+\log n ,
$$ (eq-02-02-lse-bounds)

because $e^{\max_i x_i}\le\sum_ie^{x_i}\le n\,e^{\max_i x_i}$. For
$\vx=(1,2,3)$, $\operatorname{lse}(\vx)=3.408$, between $3$ and
$3+\log3=4.099$. Its Hessian is $\diag(\vp)-\vp\vp^\T$ with
$p_i=e^{x_i}/\sum_ke^{x_k}$, so $\vp\ge\vzero$ and $\vone^\T\vp=1$. For any
$\vv$, by the Cauchy–Schwarz inequality,

$$
\vv^\T\nabla^2\operatorname{lse}(\vx)\vv
=\sum_ip_iv_i^2-\Bigl(\sum_ip_iv_i\Bigr)^2\ge0 ,
$$

since $\bigl(\sum_ip_iv_i\bigr)^2=\bigl(\sum_i\sqrt{p_i}\,\sqrt{p_i}v_i\bigr)^2\le
\sum_ip_i\sum_ip_iv_i^2$. Hence log-sum-exp is convex. It appears in
logit route-choice and discrete-choice models of travel demand.

(sec-02-02-operations)=
## Operations that preserve convexity

As with sets ({ref}`sec-02-01-operations`), convexity is usually proved by
assembling a function from known convex pieces.

:::{prf:theorem} Operations that preserve convexity
:label: thm-02-02-operations

1. **Nonnegative weighted sums.** If $f_1,\dots,f_k$ are convex and
   $w_1,\dots,w_k\ge0$, then $\sum_iw_if_i$ is convex.
2. **Affine composition.** If $f$ is convex, then $g(\vx)=f(\mA\vx+\vb)$ is
   convex.
3. **Pointwise maximum and supremum.** If $f_i$ is convex for every $i$ in
   an index set $\mathcal I$ (finite or infinite), then
   $g(\vx)=\sup_{i\in\mathcal I}f_i(\vx)$ is convex.
4. **Scalar composition.** Let $g:\R^n\to\R$ and $h:\R\to\R$. Then
   $h(g(\vx))$ is convex if $h$ is convex and nondecreasing and $g$ is
   convex, or if $h$ is convex and nonincreasing and $g$ is concave.
5. **Partial minimization.** If $f(\vx,\vy)$ is jointly convex in
   $(\vx,\vy)$ and $C$ is a nonempty convex set, then
   $g(\vx)=\inf_{\vy\in C}f(\vx,\vy)$ is convex (where it is finite).
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Parts 1 and 2 follow directly from {eq}`eq-02-02-convex`; in part 2 the map
$\vx\mapsto\mA\vx+\vb$ preserves convex combinations. For part 3,
$\operatorname{epi} g=\bigcap_i\operatorname{epi} f_i$ is an intersection
of convex sets ({prf:ref}`thm-02-01-intersection`). For the first case of
part 4, convexity of $g$ and monotonicity of $h$ give
$h(g(\theta\vx+(1-\theta)\vy))\le h(\theta g(\vx)+(1-\theta)g(\vy))$, and
convexity of $h$ bounds this by $\theta h(g(\vx))+(1-\theta)h(g(\vy))$;
the second case is the same with the inequality for $g$ reversed. For
part 5, given $\varepsilon>0$ choose $\vy_1,\vy_2\in C$ with
$f(\vx_1,\vy_1)\le g(\vx_1)+\varepsilon$ and
$f(\vx_2,\vy_2)\le g(\vx_2)+\varepsilon$. Since
$\theta\vy_1+(1-\theta)\vy_2\in C$,

$$
\begin{aligned}
&g\bigl(\theta\vx_1+(1-\theta)\vx_2\bigr)\\
&\le f\bigl(\theta\vx_1+(1-\theta)\vx_2,\ \theta\vy_1+(1-\theta)\vy_2\bigr)\\
&\le\theta g(\vx_1)+(1-\theta)g(\vx_2)+\varepsilon ,
\end{aligned}
$$

and $\varepsilon$ is arbitrary.
:::

Some consequences:

- $\norm{\mA\vx-\vb}$, $\norm{\mA\vx-\vb}^2$ and
  $\max_i(\va_i^\T\vx-b_i)$ are convex (parts 2 and 3);
- $e^{g(\vx)}$ is convex when $g$ is convex, and $1/g(\vx)$ and
  $-\log g(\vx)$ are convex when $g$ is concave and positive (part 4 with
  $h(u)=e^u$, $h(u)=1/u$ and $h(u)=-\log u$ on $u>0$);
- the distance $\inf_{\vy\in C}\norm{\vx-\vy}$ from a point to a convex set
  is convex (part 5);
- a pointwise **infimum** of affine functions, such as the dual function
  {eq}`eq-02-03-dual-function`, is concave (part 3 applied to the
  negatives).

The product and the composition of convex functions are not convex in
general: $x$ and $x^2$ are convex, but $x\cdot x^2=x^3$ is not, and
$h(u)=-u$ is convex while $h(x^2)=-x^2$ is concave.

**An infrastructure example.** In traffic assignment the travel time on a
road link is a function of the link flow $x$ (veh/h). A widely used form is
the link performance function of the US Bureau of Public Roads (BPR),

$$
t(x)=t_0\Bigl(1+0.15\bigl(x/c\bigr)^4\Bigr),
$$ (eq-02-02-bpr)

with free-flow travel time $t_0$ and practical capacity $c$. Its second
derivative $1.8\,t_0x^2/c^4$ is nonnegative, so $t$ is convex for $x\ge0$:
each additional vehicle delays the others more than the previous one did.
The total travel time of all vehicles on the link,
$x\,t(x)=t_0x+0.15\,t_0x^5/c^4$, is the sum of a linear function and a
function with second derivative $3\,t_0x^3/c^4\ge0$, hence convex for
$x\ge0$, and the total travel time of a network, a sum over its links, is
a convex function of the link flows (part 1). This is why the system-optimal traffic assignment is a
convex optimization problem (Sheffi, 1985).

(sec-02-02-sublevel)=
## Sublevel sets

:::{prf:definition} Sublevel set
:label: def-02-02-sublevel-set

The **sublevel set** of $f$ at level $\gamma\in\R$ is
$S_\gamma=\{\vx\in\dom f: f(\vx)\le\gamma\}$.
:::

:::{prf:theorem} Sublevel sets of convex functions are convex
:label: thm-02-02-sublevel

If $f$ is convex, every sublevel set $S_\gamma$ is convex.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

If $f(\vx)\le\gamma$ and $f(\vy)\le\gamma$, then
$f(\theta\vx+(1-\theta)\vy)\le\theta\gamma+(1-\theta)\gamma=\gamma$.
:::

{numref}`fig-02-02-epigraph-sublevel`(b) shows sublevel sets of a convex
quadratic: they are filled ellipses, in line with
{prf:ref}`def-02-01-ellipsoid`. The theorem is the reason why constraints in
convex optimization are written as $g_i(\vx)\le0$ with convex $g_i$: the
set $\{\vx:g_i(\vx)\le0\}$ is the sublevel set $S_0$ of $g_i$, hence
convex, and the feasible set $\bigcap_i\{\vx:g_i(\vx)\le0\}$ is convex by
{prf:ref}`thm-02-01-intersection`. Similarly, $\{\vx:g(\vx)\ge0\}$ is convex
when $g$ is concave. An *equality* $g(\vx)=0$ with a nonlinear convex $g$
does not give a convex set: $x_1^2+x_2^2-1=0$ is a circle.

```{figure} figures/02-02-epigraph-sublevel.svg
:name: fig-02-02-epigraph-sublevel
:width: 640px
:class: dark-light
:alt: Left, a convex curve with the region above it shaded as the epigraph, a horizontal line at level gamma and the interval where the curve is below the line. Right, nested elliptic contour lines of a convex quadratic with one sublevel set shaded

(a) The epigraph of a convex function of one variable (shaded) is a convex
set. The sublevel set $S_\gamma$ is the interval of $x$ where the graph is
at or below the level $\gamma$. (b) Contour lines of the convex quadratic
$f(\vx)=(x_1-1)^2+2(x_2-1)^2-\tfrac12x_1x_2$; the sublevel set $S_1$ is
shaded. All sublevel sets are nested convex ellipses.
```

The converse of {prf:ref}`thm-02-02-sublevel` is false. The function
$\sqrt{\abs{x}}$ is not convex (it is concave on $[0,\infty)$), but its
sublevel sets $\{x:\abs{x}\le\gamma^2\}$ are intervals. Functions whose
sublevel sets are all convex form the larger class of quasiconvex
functions ({ref}`sec-02-02-quasiconvex`).

(sec-02-02-jensen)=
## Jensen's inequality

The definition {eq}`eq-02-02-convex` involves two points. It extends to any
finite number of points, exactly as convex sets contain all convex
combinations ({prf:ref}`thm-02-01-convex-combinations`).

:::{prf:theorem} Jensen's inequality
:label: thm-02-02-jensen

Let $f$ be convex, $\vx_1,\dots,\vx_k\in\dom f$ and $\theta_i\ge0$ with
$\sum_{i=1}^k\theta_i=1$. Then

$$
f\Bigl(\sum_{i=1}^k\theta_i\vx_i\Bigr)\le\sum_{i=1}^k\theta_if(\vx_i).
$$ (eq-02-02-jensen)

If $X$ is a random vector with values in $\dom f$ and finite expectation,
then $f(\mathrm E\,X)\le\mathrm E\,f(X)$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Apply {prf:ref}`thm-02-01-convex-combinations` to $\operatorname{epi} f$,
which is convex by {prf:ref}`thm-02-02-epigraph`: the points
$(\vx_i,f(\vx_i))$ lie in the epigraph, hence so does their convex
combination $\bigl(\sum_i\theta_i\vx_i,\ \sum_i\theta_if(\vx_i)\bigr)$,
which is {eq}`eq-02-02-jensen`. The probabilistic form is the same statement
for a discrete distribution with probabilities $\theta_i$, and follows for
general distributions by a limiting argument.
:::

Jensen's inequality is the source of many classical inequalities. With the
concave function $\log x$, for instance, it gives the
arithmetic–geometric mean inequality
$\sqrt{ab}\le\tfrac12(a+b)$ for $a,b>0$, since
$\log\tfrac12(a+b)\ge\tfrac12\log a+\tfrac12\log b$. In planning, it
warns that a convex cost evaluated at the *average* demand underestimates
the *average* cost.

:::{prf:example} Travel time under fluctuating flow
:label: ex-02-02-jensen-travel-time

A road link has free-flow travel time $t_0=10$ min and capacity
$c=2{,}000$ veh/h, with the travel time {eq}`eq-02-02-bpr`. During the
morning peak the flow is $1{,}600$ veh/h on half of the days and
$2{,}400$ veh/h on the other half, so the average flow is $2{,}000$ veh/h.
Compare the travel time at the average flow with the average travel time.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

At the average flow, $t(2000)=10\,(1+0.15\cdot1^4)=11.5$ min. On the two
kinds of days,

$$
\begin{aligned}
t(1600)&=10\,\bigl(1+0.15\cdot0.8^4\bigr)=10\,(1+0.15\cdot0.4096)=10.61\ \text{min},\\
t(2400)&=10\,\bigl(1+0.15\cdot1.2^4\bigr)=10\,(1+0.15\cdot2.0736)=13.11\ \text{min},
\end{aligned}
$$

so the average travel time is $\tfrac12(10.61+13.11)=11.86$ min. As
Jensen's inequality predicts for the convex function $t$, the average
travel time exceeds the travel time at the average flow, here by
$0.36$ min or $3.2\,\%$. A plan that uses average flows underestimates
the average delay, and the error grows with the variability of the flow
and with the curvature of $t$.
:::

(sec-02-02-quasiconvex)=
## Quasiconvex functions

Some useful functions are not convex but share the convexity of their
sublevel sets.

:::{prf:definition} Quasiconvex function
:label: def-02-02-quasiconvex

Let $C$ be convex and $f:C\to\R$. The function $f$ is **quasiconvex** if

$$
f\bigl(\theta\vx+(1-\theta)\vy\bigr)\le\max\{f(\vx),f(\vy)\}
\qquad\text{for all }\vx,\vy\in C,\ \theta\in[0,1],
$$ (eq-02-02-quasiconvex)

and **strictly quasiconvex** if the inequality is strict whenever
$\vx\ne\vy$ and $\theta\in(0,1)$. The function is (strictly)
**quasiconcave** if $-f$ is (strictly) quasiconvex.
:::

A function is quasiconvex if and only if all its sublevel sets are convex:
if $f(\vx),f(\vy)\le\gamma$, then {eq}`eq-02-02-quasiconvex` gives
$f(\theta\vx+(1-\theta)\vy)\le\gamma$, and conversely taking
$\gamma=\max\{f(\vx),f(\vy)\}$ gives {eq}`eq-02-02-quasiconvex`. Every
convex function is quasiconvex, since the chord value in
{eq}`eq-02-02-convex` is at most $\max\{f(\vx),f(\vy)\}$, and every strictly
convex function is strictly quasiconvex. Every monotone function of one
variable is both quasiconvex and quasiconcave, and $\sqrt{\abs{x}}$ and
$\log x$ (on $x>0$) are quasiconvex without being convex.

In one variable, a strictly quasiconvex function on an interval is
**unimodal**: it decreases up to its minimizer and increases after it, and
no differentiability is needed. Line-search methods such as golden-section
search rely only on this property ({prf:ref}`def-06-04-unimodal`,
{prf:ref}`alg-06-04-golden-section`). Two warnings
apply. A quasiconvex function may have flat parts, so a local minimum need
not be global unless the function is strictly quasiconvex. And sums of
quasiconvex functions need not be quasiconvex: $x^3$ and $-3x$ are
monotone, hence quasiconvex, but $x^3-3x$ has the sublevel set
$\{x:x^3-3x\le0\}=(-\infty,-\sqrt3]\cup[0,\sqrt3]$, which is not an
interval.

## Summary and Common Mistakes

### Summary

- A function on a convex set is convex if every chord lies above its graph
  ({prf:ref}`def-02-02-convex-function`, {numref}`fig-02-02-chord-tangent`), and
  if and only if its epigraph is convex ({prf:ref}`thm-02-02-epigraph`).
  Affine functions are convex and concave; a parabola is convex if and only
  if its leading coefficient is nonnegative ({eq}`eq-02-02-parabola-gap`).
- A differentiable function is convex if and only if it lies above its
  tangent planes ({prf:ref}`thm-02-02-first-order`); hence a stationary point
  of a convex function is a global minimizer.
- A twice continuously differentiable function on an open convex set is
  convex if and only if its Hessian is positive semidefinite; a positive
  definite Hessian implies strict convexity, but not conversely
  ({prf:ref}`thm-02-02-second-order`, {prf:ref}`rem-02-02-strict-hessian`,
  {prf:ref}`rem-02-02-hessian-open`). The function
  $x_1^2-2x_1x_2+x_2^2$ is convex but not strictly convex
  ({prf:ref}`ex-02-02-hessian-test`).
- Norms, convex quadratics, exponentials, $-\log x$, maxima of affine
  functions and log-sum-exp are convex; nonnegative sums, affine
  compositions, pointwise maxima, monotone compositions and partial
  minimization preserve convexity ({prf:ref}`thm-02-02-operations`).
- Sublevel sets of convex functions are convex
  ({prf:ref}`thm-02-02-sublevel`), so constraints $g_i(\vx)\le0$ with convex
  $g_i$ define convex feasible sets.
- Jensen's inequality ({prf:ref}`thm-02-02-jensen`) extends convexity to
  averages; for a convex travel-time function, the average travel time
  exceeds the travel time at the average flow
  ({prf:ref}`ex-02-02-jensen-travel-time`).
- Quasiconvex functions ({prf:ref}`def-02-02-quasiconvex`) are exactly the
  functions with convex sublevel sets; strictly quasiconvex functions of
  one variable are unimodal.

### Common mistakes

- **Forgetting the domain.** Convexity is a property of a function *on a
  convex set*; $1/x$ is convex on $x>0$ only.
- **Requiring strictness at all $\theta\in[0,1]$.** Strict convexity
  compares distinct points with $\theta\in(0,1)$; at the end points both
  sides of {eq}`eq-02-02-convex` are equal.
- **Reading "$\nabla^2f\succ0$" as necessary for strict convexity.** It is
  sufficient only; $x^4$ is strictly convex with $f''(0)=0$.
- **Using leading principal minors to test semidefiniteness.** Leading
  minors decide *definiteness*; semidefiniteness needs all principal
  minors or the eigenvalues ($\diag(0,-1)$ is a counterexample).
- **Applying the Hessian test as "if and only if" on sets without
  interior** ({prf:ref}`rem-02-02-hessian-open`).
- **Writing nonlinear equality constraints in a "convex" model.** Only
  affine equalities keep the feasible set convex.
- **Assuming that products or compositions of convex functions are
  convex.** Use the composition rules of {prf:ref}`thm-02-02-operations`.

## References

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
