# Revised and Dual Simplex Methods

## Overview

The tableau method of {doc}`03-03-simplex-method` recomputes every entry of the
tableau at every pivot, although one iteration only needs the reduced costs,
one column and the right-hand side. The **revised simplex method** keeps the
original data $\mA,\vb,\vc$ untouched and carries only the inverse of the
current basis matrix, $\mB^{-1}$; everything else is computed from it on
demand. This chapter writes a simplex iteration in matrix form, shows that a
pivot multiplies $\mB^{-1}$ on the left by an **elementary (eta) matrix**, so
that $\mB^{-1}$ is a product of such matrices (the **product form of the
inverse**), and states the revised algorithm, its two-phase version with the
bordered inverse $\mT$, and its operation counts. The second half turns the
simplex method around: the **dual simplex method** starts from a basis whose
reduced costs are optimal but whose basic solution is infeasible, and pivots
toward feasibility while keeping the reduced costs optimal. It is the tool
for re-optimizing after a constraint is added or a right-hand side changes,
which is why {doc}`03-07-sensitivity-analysis` relies on it. Three examples are
worked in full: a revised simplex run in product form, a two-phase revised
run, and a dual simplex pivot. The chapter assumes block matrices and
elementary matrices from {doc}`01-01-linear-algebra`, the two-phase method of
{doc}`03-04-two-phase-big-m`, and the dual linear program (LP) and dual
feasibility of
{doc}`03-05-lp-duality`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- write the quantities of a simplex iteration ($\bar\vb$, $\vpi$, $\bar c_j$,
  $\bar\mA_{\cdot s}$) in terms of $\mB^{-1}$ and the original data;
- build the eta matrix of a pivot and update $\mB^{-1}$ in product form;
- run the revised simplex method, including its two-phase version with the
  bordered inverse $\mT$, by hand;
- compare the work per iteration of the tableau and the revised methods;
- derive the leaving and entering rules of the dual simplex method and run
  it, in tableau and in revised form; and
- re-optimize after adding a constraint, and start the dual simplex method
  with an artificial bounding constraint when no dual feasible basis is at
  hand.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\mB_k$ | basis matrix in iteration $k$ |
| $\mE_k$ | eta (elementary, pivot) matrix of pivot $k$: $\mB_{k+1}^{-1}=\mE_k\mB_k^{-1}$ |
| $\veta$ | eta column, the column of $\mE_k$ that differs from the identity |
| $\bar\mA_{\cdot s}=\mB^{-1}\mA_{\cdot s}$ | updated column of the entering variable $x_s$ |
| $\bar\mA_{r\cdot}$, $\bar a_{rj}$ | row $r$ of $\mB^{-1}\mA$ (pivot row of the dual simplex method) and its entries |
| $\ve_r$ | $r$-th unit vector |
| $\hat\mB$ | bordered basis matrix $\begin{pmatrix}\mB&\vzero\\\vc_B^\T&1\end{pmatrix}$; $\mT=\hat\mB^{-1}$ |
| $\theta_D$ | step of the dual simplex method (dual ratio) |

(sec-03-06-why-revised)=
## Why a revised method

Consider the standard-form LP

$$
\begin{aligned}
\min_{\vx}\ & \vc^\T\vx \\
\st\ & \mA\vx=\vb,\\
& \vx\ge\vzero,
\end{aligned}
$$ (eq-03-06-lp)

with $\mA\in\R^{m\times n}$ of rank $m$. For a basis $\mB$ (the columns of
the basic variables $\vx_B$) and the nonbasic matrix $\mN$, the simplex
tableau of {doc}`03-03-simplex-method` holds exactly the following quantities:

$$
\begin{aligned}
\bar\vb &= \mB^{-1}\vb && \text{values of the basic variables,}\\
\vpi^\T &= \vc_B^\T\mB^{-1} && \text{simplex multipliers (dual solution),}\\
\bar c_j &= c_j-\vpi^\T\mA_{\cdot j} && \text{reduced cost of } x_j,\\
\bar\mA_{\cdot j} &= \mB^{-1}\mA_{\cdot j} && \text{updated column of } x_j,\\
-z &= -\vc_B^\T\bar\vb = -\vpi^\T\vb && \text{objective entry.}
\end{aligned}
$$ (eq-03-06-matrix-quantities)

One iteration of the method (for a minimization) uses only part of them:

1. the reduced costs $\bar c_j$ of the nonbasic variables, to test
   optimality and pick an entering variable $x_s$ with $\bar c_s<0$;
2. the single column $\bar\mA_{\cdot s}$, for the ratio test
   $\theta=\min\{\bar b_i/\bar a_{is}:\bar a_{is}>0\}$ and the unboundedness
   test $\bar\mA_{\cdot s}\le\vzero$;
3. the right-hand side $\bar\vb$.

The other $n-m-1$ updated columns are computed and thrown away at every
pivot. Two observations make this waste large in practice. First,
computational experience shows that the simplex method usually reaches an
optimum in a number of pivots of the order of $m$ (figures between about $m$
and $3m$ are commonly reported), so most of the $n$ columns never enter the
basis when $n\gg m$. Second, real constraint matrices are **sparse** (most
$a_{ij}$ are zero), and the tableau destroys sparsity: after a few pivots
$\mB^{-1}\mA$ is dense even when $\mA$ is not.

The **revised simplex method** therefore stores the original data, which
never change, and a representation of $\mB^{-1}$, and evaluates
{eq}`eq-03-06-matrix-quantities` only where an iteration needs it. All that
remains is to update $\mB^{-1}$ cheaply after a pivot.

(sec-03-06-product-form)=
## Elementary matrices and the product form of the inverse

A pivot replaces one basic column: the entering column $\mA_{\cdot s}$ takes
the place of the column of the leaving variable, which sits in position $r$
of the basis. The new inverse is the old one multiplied on the left by a
matrix that differs from the identity in one column.

:::{prf:definition} Eta matrix of a pivot
:label: def-03-06-eta

Let $\bar\mA_{\cdot s}=\mB^{-1}\mA_{\cdot s}$ with pivot element
$\bar a_{rs}\ne0$. The **eta matrix** (pivot matrix, elementary matrix)
$\mE$ of the pivot on $(r,s)$ is the $m\times m$ identity matrix with
column $r$ replaced by the **eta column**

$$
\veta=\left(-\frac{\bar a_{1s}}{\bar a_{rs}},\ \dots,\
-\frac{\bar a_{r-1,s}}{\bar a_{rs}},\ \frac{1}{\bar a_{rs}},\
-\frac{\bar a_{r+1,s}}{\bar a_{rs}},\ \dots,\
-\frac{\bar a_{ms}}{\bar a_{rs}}\right)^\T .
$$ (eq-03-06-eta)
:::

Multiplying a vector $\vv$ on the left by $\mE$ performs exactly the
Gauss–Jordan row operations of a tableau pivot: row $r$ is divided by
$\bar a_{rs}$, and $\bar a_{is}/\bar a_{rs}$ times the old row $r$ is
subtracted from every other row $i$,

$$
\begin{gathered}
(\mE\vv)_r=\frac{v_r}{\bar a_{rs}},\\
(\mE\vv)_i=v_i-\frac{\bar a_{is}}{\bar a_{rs}}\,v_r\quad(i\ne r).
\end{gathered}
$$ (eq-03-06-eta-action)

In particular $\mE\bar\mA_{\cdot s}=\ve_r$: the pivot turns the entering
column into a unit column. (This is the pivot matrix {eq}`eq-01-01-eta-matrix`
of {prf:ref}`prop-01-01-pivot-matrix`, with $\vw=\bar\mA_{\cdot s}$; as an
elementary-matrix product it is $m-1$ row additions and one row scaling.)

:::{prf:theorem} Update of the basis inverse
:label: thm-03-06-update

Let $\mB'$ be obtained from $\mB$ by replacing its $r$-th column by
$\mA_{\cdot s}$, where $\bar a_{rs}\ne0$. Then $\mB'$ is nonsingular and

$$
\mB'^{-1}=\mE\,\mB^{-1},
$$ (eq-03-06-update)

with $\mE$ the eta matrix {eq}`eq-03-06-eta` of the pivot on $(r,s)$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Write $\mF$ for the identity with column $r$ replaced by
$\bar\mA_{\cdot s}$. Column $i\ne r$ of $\mB\mF$ is $\mB\ve_i$, the $i$-th
column of $\mB$, and column $r$ is $\mB\bar\mA_{\cdot s}=\mA_{\cdot s}$.
Hence $\mB'=\mB\mF$. The matrix $\mF$ is the identity except for column $r$,
so $\det\mF=\bar a_{rs}\ne0$ and $\mB'$ is nonsingular. By
{eq}`eq-03-06-eta-action`, $\mE\ve_i=\ve_i$ for $i\ne r$ and
$\mE\bar\mA_{\cdot s}=\ve_r$, so $\mE\mF=\mI$, that is $\mE=\mF^{-1}$.
Therefore $\mB'^{-1}=\mF^{-1}\mB^{-1}=\mE\mB^{-1}$.
:::

If the method starts from a slack basis, $\mB_1=\mI$, then after $k-1$
pivots

$$
\mB_k^{-1}=\mE_{k-1}\mE_{k-2}\cdots\mE_1 .
$$ (eq-03-06-product-form)

This is the **product form of the inverse** (Dantzig & Orchard-Hays, 1954),
already met as {eq}`eq-01-01-product-form` (there with the pivots counted from
$\mE_1$ at the first basis change). An eta matrix
is stored as its position $r$ and its column $\veta$, that is
$m+1$ numbers, often fewer because $\bar\mA_{\cdot s}$ is sparse. The
products in {eq}`eq-03-06-matrix-quantities` are evaluated without ever forming
$\mB^{-1}$:

* $\bar\mA_{\cdot s}=\mE_{k-1}(\cdots(\mE_1\mA_{\cdot s}))$, applying the
  etas from the first to the last (each application costs about $m$
  multiplications by {eq}`eq-03-06-eta-action`);
* $\vpi^\T=(\cdots(\vc_B^\T\mE_{k-1})\cdots)\mE_1$, applying them from the
  last to the first; multiplying a row vector $\vu^\T$ by $\mE$ changes
  only its $r$-th entry, which becomes $\vu^\T\veta$.

The list of etas grows by one per pivot, so the cost of these products grows
too, and rounding errors accumulate. Codes therefore **reinvert**
periodically: they factorize the current $\mB$ afresh (today as an LU
factorization rather than as a product of etas) and restart the list. Small
hand examples, like the one below, simply carry $\mB^{-1}$ explicitly and
update it with {eq}`eq-03-06-update`.

(sec-03-06-revised-algorithm)=
## The revised simplex algorithm

:::{prf:algorithm} Revised simplex method (minimization)
:label: alg-03-06-revised-simplex

**Input:** the data $\mA,\vb,\vc$ of {eq}`eq-03-06-lp` and a feasible basis
$\mB$ with $\mB^{-1}$ (or its eta factors) and $\vx_B=\mB^{-1}\vb\ge\vzero$.
**Output:** an optimal basis, or a proof of unboundedness.

1. **Price.** Compute the multipliers $\vpi^\T=\vc_B^\T\mB^{-1}$ and the
   reduced costs $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ of the nonbasic
   variables. If $\bar c_j\ge0$ for all of them, stop: $\vx_B=\mB^{-1}\vb$,
   $\vx_N=\vzero$ is optimal and $z^*=\vpi^\T\vb$.
2. **Choose the entering variable.** Pick $x_s$ with $\bar c_s<0$ (for
   example the most negative) and compute
   $\bar\mA_{\cdot s}=\mB^{-1}\mA_{\cdot s}$.
3. **Ratio test.** If $\bar\mA_{\cdot s}\le\vzero$, stop: the LP is
   unbounded along the ray that increases $x_s$. Otherwise compute
   $\theta=\min\{\bar b_i/\bar a_{is}:\bar a_{is}>0\}$, with
   $\bar\vb=\vx_B$, and let $r$ be a row attaining the minimum. The basic
   variable of row $r$ leaves.
4. **Update.** Form the eta matrix $\mE$ of the pivot on $(r,s)$ by
   {eq}`eq-03-06-eta`, set $\mB^{-1}\leftarrow\mE\mB^{-1}$ and
   $\vx_B\leftarrow\mE\vx_B$ (the new values are $\theta$ for $x_s$ in row
   $r$ and $\bar b_i-\theta\bar a_{is}$ elsewhere), replace the basic
   variable of row $r$ by $x_s$, and go to Step 1.
:::

```{figure} figures/03-06-revised-iteration.svg
:name: fig-03-06-revised-iteration
:width: 560px
:class: dark-light
:alt: Cycle of four boxes: price, choose the entering variable, ratio test, update the inverse, with exits for an optimal and an unbounded problem

One iteration of the revised simplex method
({prf:ref}`alg-03-06-revised-simplex`). Only $\mB^{-1}$ (or its eta factors)
is carried from one iteration to the next; $\mA$, $\vb$ and $\vc$ are never
changed.
```

{numref}`fig-03-06-revised-iteration` shows the four steps as a cycle. For a
maximization the optimality test reads $\bar c_j\le0$ and the entering
variable has $\bar c_s>0$; the other steps are unchanged. The method visits
the same bases as the tableau method when it uses the same entering and
leaving rules, so finiteness, degeneracy and Bland's rule
({prf:ref}`alg-03-03-bland`) carry over unchanged. By hand, a working table with the columns
Basic $\mid\mB^{-1}\mid\vx_B\mid\bar\mA_{\cdot s}\mid$ ratio replaces the
full tableau.

:::{prf:example} A revised simplex run in product form
:label: ex-03-06-revised-run

Solve

$$
\begin{aligned}
\min_{\vx}\ & -3x_1-x_2-3x_3 \\
\st\ & 2x_1+x_2+x_3+x_4=2,\\
& x_1+2x_2+3x_3+x_5=5,\\
& 2x_1+2x_2+x_3+x_6=6,\\
& x_1,\dots,x_6\ge0
\end{aligned}
$$

by the revised simplex method from the slack basis $(x_4,x_5,x_6)$, writing
every $\mB_k^{-1}$ as a product of eta matrices.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The data are
$\mA_{\cdot1}=(2,1,2)^\T$, $\mA_{\cdot2}=(1,2,2)^\T$,
$\mA_{\cdot3}=(1,3,1)^\T$, the slack columns $\ve_1,\ve_2,\ve_3$,
$\vb=(2,5,6)^\T$ and $\vc=(-3,-1,-3,0,0,0)^\T$. All arithmetic below is
exact.

**Iteration 1.** $\mB_1=\mI$, so $\mB_1^{-1}=\mI$, $\vx_B=(2,5,6)^\T$,
$z=0$ and $\vpi^\T=(0,0,0)$. The reduced costs are
$(\bar c_1,\bar c_2,\bar c_3)=(-3,-1,-3)$. Any variable with a negative
reduced cost may enter; we let $x_2$ enter because its pivot will be 1 and
keeps the hand arithmetic in integers (the most-negative rule would pick
$x_1$; see the remark after this solution).

| Basic | $\mB_1^{-1}$ | | | $\vx_B$ | $\bar\mA_{\cdot2}$ | ratio |
|---|---|---|---|---|---|---|
| $x_4$ | 1 | 0 | 0 | 2 | **1** | 2 |
| $x_5$ | 0 | 1 | 0 | 5 | 2 | 5/2 |
| $x_6$ | 0 | 0 | 1 | 6 | 2 | 3 |

The minimum ratio is 2 in row 1, so $x_4$ leaves. By {eq}`eq-03-06-eta` the eta
column is $(1/1,-2/1,-2/1)^\T$ in position 1:

$$
\begin{gathered}
\mE_1=\begin{pmatrix}1&0&0\\-2&1&0\\-2&0&1\end{pmatrix},\\
\mB_2^{-1}=\mE_1\mB_1^{-1}=\mE_1 .
\end{gathered}
$$

**Iteration 2.** Basis $(x_2,x_5,x_6)$, $\vx_B=\mB_2^{-1}\vb=(2,1,2)^\T$,
$z=-2$. Pricing: $\vpi^\T=\vc_B^\T\mB_2^{-1}=(-1,0,0)\mB_2^{-1}=(-1,0,0)$,
so $\bar c_1=-3-(-1)(2)=-1$, $\bar c_3=-3-(-1)(1)=-2$ and
$\bar c_4=0-(-1)(1)=1$. The most negative is $\bar c_3$; $x_3$ enters with
$\bar\mA_{\cdot3}=\mB_2^{-1}(1,3,1)^\T=(1,1,-1)^\T$.

| Basic | $\mB_2^{-1}$ | | | $\vx_B$ | $\bar\mA_{\cdot3}$ | ratio |
|---|---|---|---|---|---|---|
| $x_2$ | 1 | 0 | 0 | 2 | 1 | 2 |
| $x_5$ | −2 | 1 | 0 | 1 | **1** | 1 |
| $x_6$ | −2 | 0 | 1 | 2 | −1 | — |

Row 2 wins, $x_5$ leaves, and the eta column $(-1,1,1)^\T$ sits in
position 2:

$$
\begin{gathered}
\mE_2=\begin{pmatrix}1&-1&0\\0&1&0\\0&1&1\end{pmatrix},\\
\mB_3^{-1}=\mE_2\mE_1=\begin{pmatrix}3&-1&0\\-2&1&0\\-4&1&1\end{pmatrix}.
\end{gathered}
$$

**Iteration 3.** Basis $(x_2,x_3,x_6)$, $\vx_B=(1,1,3)^\T$, $z=-4$.
Pricing: $\vpi^\T=(-1,-3,0)\mB_3^{-1}=(3,-2,0)$, so $\bar c_1=-3-(6-2)=-7$,
$\bar c_4=0-3=-3$ and $\bar c_5=0-(-2)=2$. Now $x_1$ enters with
$\bar\mA_{\cdot1}=\mB_3^{-1}(2,1,2)^\T=(5,-3,-5)^\T$.

| Basic | $\mB_3^{-1}$ | | | $\vx_B$ | $\bar\mA_{\cdot1}$ | ratio |
|---|---|---|---|---|---|---|
| $x_2$ | 3 | −1 | 0 | 1 | **5** | 1/5 |
| $x_3$ | −2 | 1 | 0 | 1 | −3 | — |
| $x_6$ | −4 | 1 | 1 | 3 | −5 | — |

Only row 1 has a positive entry, so $x_2$ leaves with $\theta=1/5$. The eta
column is $(1/5,\,3/5,\,5/5)^\T$ in position 1:

$$
\begin{gathered}
\mE_3=\begin{pmatrix}1/5&0&0\\3/5&1&0\\1&0&1\end{pmatrix},\\
\mB_4^{-1}=\mE_3\mE_2\mE_1=\begin{pmatrix}3/5&-1/5&0\\-1/5&2/5&0\\-1&0&1\end{pmatrix}.
\end{gathered}
$$

**Iteration 4.** Basis $(x_1,x_3,x_6)$, $\vx_B=\mB_4^{-1}\vb=(1/5,8/5,4)^\T$
and $z=-3\cdot\tfrac15-3\cdot\tfrac85=-\tfrac{27}{5}$. Pricing:
$\vpi^\T=(-3,-3,0)\mB_4^{-1}=(-6/5,-3/5,0)$ and

$$
\begin{gathered}
\bar c_2=-1-\left(-\tfrac65-\tfrac65\right)=\tfrac75,\\
\bar c_4=\tfrac65,\\
\bar c_5=\tfrac35 .
\end{gathered}
$$

All reduced costs are nonnegative, so the basis is optimal:

$$
\begin{gathered}
\vx^*=\left(\tfrac15,\,0,\,\tfrac85,\,0,\,0,\,4\right)^\T,\\
z^*=-\tfrac{27}{5}=-5.4 .
\end{gathered}
$$

**Checks.** The product $\mB_4^{-1}\mB_4$ with
$\mB_4=(\mA_{\cdot1},\mA_{\cdot3},\mA_{\cdot6})$ is the identity, and the
dual objective $\vpi^\T\vb=-\tfrac65\cdot2-\tfrac35\cdot5=-\tfrac{27}{5}$
equals $z^*$, as strong duality requires ({prf:ref}`thm-03-05-strong-duality`). A
solver (HiGHS) returns the same solution. For comparison,
{numref}`tab-03-06-final-tableau` is the full final tableau
$\mB_4^{-1}\mA$ that the tableau method would have produced; the revised
method never formed it, yet its reduced-cost row and right-hand side are
exactly what Step 1 computed.

```{table} Final tableau of the revised simplex run, computed as $\mB_4^{-1}\mA$ and $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ for comparison.
:name: tab-03-06-final-tableau

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|---|---|---|---|---|---|---|---|
| $x_1$ | 1 | 1/5 | 0 | 3/5 | −1/5 | 0 | 1/5 |
| $x_3$ | 0 | 3/5 | 1 | −1/5 | 2/5 | 0 | 8/5 |
| $x_6$ | 0 | 1 | 0 | −1 | 0 | 1 | 4 |
| $-z$ | 0 | 7/5 | 0 | 6/5 | 3/5 | 0 | 27/5 |
```

*(Earlier versions of this example print the third row of $\mE_3$ as
$(-1,0,1)$, the third row of $\mB_4^{-1}$ as $(-1,1,1)$, two entries of the
final tableau as $-1/5$ and $(0,3,0,-1,1,1)$, and the solution with five
components. The values above are the correct ones: the third row of $\mE_3$
is $-\bar a_{31}/\bar a_{11}=-(-5)/5=1$, and with the misprinted
$\mB_4^{-1}$ the value of $x_6$ would be 9 instead of 4.)*
:::

The entering choice in iteration 1 of the solution of
{prf:ref}`ex-03-06-revised-run` deliberately broke the most-negative rule. With
Dantzig's rule ($x_1$ first, ties broken by the smaller index) the same LP is
solved in two pivots, $x_1$ replacing $x_4$ and then $x_3$ replacing $x_5$,
ending in the same basis $(x_1,x_3,x_6)$. Neither rule is best on every
problem.

(sec-03-06-operation-counts)=
## Operation counts

Count multiplications per iteration, the dominant cost, for dense data.

* **Tableau method.** A pivot updates all $(m+1)(n+1)$ entries of the
  tableau, about one multiplication each: roughly $mn$.
* **Revised method with explicit $\mB^{-1}$.** Pricing needs
  $\vpi^\T=\vc_B^\T\mB^{-1}$ ($m^2$) and $\vpi^\T\mA_{\cdot j}$ for the
  $n-m$ nonbasic columns ($m(n-m)$); the entering column
  $\mB^{-1}\mA_{\cdot s}$ costs $m^2$; the ratio test $m$; the update
  $\mE\mB^{-1}$ by {eq}`eq-03-06-eta-action` about $m^2$. In total about
  $3m^2+m(n-m)$.

For a dense matrix the two are of the same order: with $m=3$, $n=6$
({prf:ref}`ex-03-06-revised-run`) about 28 against 36 multiplications, with $m=100$,
$n=1000$ about $1.0\times10^5$ against $1.2\times10^5$. The revised method
wins for three reasons that the dense count hides.

1. **Sparsity.** Pricing works on the original $\mA$. If only a small
   fraction of its entries are nonzero, pricing costs only that fraction of
   $m(n-m)$. With $m=1000$, $n=10\,000$ and 1% nonzeros, one tableau pivot
   costs about
   $10^7$ multiplications, a revised iteration about $3.1\times10^6$, and
   with a sparse factorization of $\mB$ the $m^2$ terms shrink as well.
2. **Partial pricing.** Step 1 may stop at the first few attractive columns
   instead of pricing all $n-m$, which the tableau method cannot exploit
   because it updates every column anyway.
3. **Memory and accuracy.** The revised method stores $\mA$ (sparse) plus
   $\mB^{-1}$ or its factors, not a dense $m\times n$ tableau, and every
   quantity is recomputed from the original data, so rounding errors do not
   accumulate in the columns.

For these reasons all production simplex codes are revised codes. The
tableau remains the clearest way to learn and to do small examples by hand.

(sec-03-06-two-phase-revised)=
## Two-phase revised simplex and the bordered inverse

When no slack basis is available, Phase I of the two-phase method
({prf:ref}`alg-03-04-two-phase`) minimizes the sum $\phi$ of artificial variables $R_i$. In revised form it
is convenient to keep the original objective as an extra constraint row,

$$
\vc^\T\vx+(-z)=0,
$$ (eq-03-06-cost-row)

where $-z$ is treated as a free variable that is always basic. The basis of
this enlarged system is the **bordered basis matrix**

$$
\begin{aligned}
\hat\mB&=\begin{pmatrix}\mB&\vzero\\\vc_B^\T&1\end{pmatrix},\\
\mT=\hat\mB^{-1}&=\begin{pmatrix}\mB^{-1}&\vzero\\-\vc_B^\T\mB^{-1}&1\end{pmatrix}
=\begin{pmatrix}\mB^{-1}&\vzero\\-\vpi^\T&1\end{pmatrix}.
\end{aligned}
$$ (eq-03-06-bordered-inverse)

The formula for $\mT$ is checked by multiplication:
$\hat\mB\mT=\begin{pmatrix}\mB\mB^{-1}&\vzero\\\vc_B^\T\mB^{-1}-\vpi^\T&1\end{pmatrix}=\mI$.
Applied to the extended column of a variable and to the extended right-hand
side, $\mT$ returns a whole column of the tableau with its reduced cost:

$$
\begin{aligned}
\mT\begin{pmatrix}\mA_{\cdot j}\\c_j\end{pmatrix}
&=\begin{pmatrix}\bar\mA_{\cdot j}\\\bar c_j\end{pmatrix},\\
\mT\begin{pmatrix}\vb\\0\end{pmatrix}
&=\begin{pmatrix}\bar\vb\\-z\end{pmatrix}.
\end{aligned}
$$ (eq-03-06-T-column)

So the last row of $\mT$ holds $-\vpi^\T$, and the last entry of
$\mT(\vb;0)$ is the objective entry $-z$ of the tableau. The matrix $\mT$
(the **T-matrix**) is updated by $(m+1)\times(m+1)$ eta matrices exactly as
$\mB^{-1}$ is; the extra row is just another row of the pivot.

In Phase I the basis consists of the artificial variables and $-z$. Pricing
uses the Phase I costs: 1 for each artificial variable and 0 for everything
else, *including the column of $-z$*. Phase I ends when $\phi=0$. If, as
usual, no artificial variable is then basic, the basic columns are columns
of the original $\mA$ and their values satisfy $\mA\vx=\vb$, $\vx\ge\vzero$:
they form a basic feasible solution of the original LP. (An artificial
variable that is still basic at value 0 is pivoted out, or its row is
dropped as redundant, as in {prf:ref}`thm-03-04-phase-one`.) Because $-z$ was
kept in the basis throughout, the current $\mT$ *is* the bordered inverse of
this basis: its upper-left block is $\mB^{-1}$ and its last row is
$(-\vpi^\T,1)$ for the original costs. Phase II starts at once, without
recomputing anything.

:::{prf:example} Two-phase revised simplex
:label: ex-03-06-two-phase-revised

Solve

$$
\begin{aligned}
\min_{\vx}\ & 3x_1-x_2-7x_3+3x_4+x_5 \\
\st\ & 5x_1-4x_2+13x_3-2x_4+x_5=20,\\
& x_1-x_2+5x_3-x_4+x_5=8,\\
& x_1,\dots,x_5\ge0
\end{aligned}
$$

by the two-phase revised simplex method with the bordered inverse $\mT$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Set-up.** Add artificial variables $x_6,x_7$ to the two rows and the cost
row {eq}`eq-03-06-cost-row`. Phase I minimizes $\phi=x_6+x_7$. The extended
columns $(\mA_{\cdot j};c_j)$ are

| | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | $-z$ | RHS |
|---|---|---|---|---|---|---|---|---|---|
| row 1 | 5 | −4 | 13 | −2 | 1 | 1 | 0 | 0 | 20 |
| row 2 | 1 | −1 | 5 | −1 | 1 | 0 | 1 | 0 | 8 |
| cost row | 3 | −1 | −7 | 3 | 1 | 0 | 0 | 1 | 0 |
| Phase I cost | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | |

The starting basis is $(x_6,x_7,-z)$ with $\mT_1=\mI$.

**Phase I, iteration 1.** The Phase I multipliers are
$(1,1,0)\mT_1=(1,1,0)$ (the third entry is the Phase I cost 0 of $-z$), and
the Phase I reduced costs $0-(1,1,0)(\mA_{\cdot j};c_j)$ of $x_1,\dots,x_5$
are $(-6,\,5,\,-18,\,3,\,-2)$. The most negative is $-18$: $x_3$ enters,
with extended column $\mT_1(13,5,-7)^\T=(13,5,-7)^\T$.

| Basic | $\mT_1$ | | | value | $\bar\mA_{\cdot3}$ | ratio |
|---|---|---|---|---|---|---|
| $x_6$ | 1 | 0 | 0 | 20 | **13** | 20/13 |
| $x_7$ | 0 | 1 | 0 | 8 | 5 | 8/5 |
| $-z$ | 0 | 0 | 1 | 0 | −7 | — |

The row of $-z$ never takes part in the ratio test ($-z$ is free). Since
$20/13<8/5$, $x_6$ leaves and

$$
\begin{gathered}
\mE_1=\begin{pmatrix}1/13&0&0\\-5/13&1&0\\7/13&0&1\end{pmatrix},\\
\mT_2=\mE_1 .
\end{gathered}
$$

**Phase I, iteration 2.** Values $\mT_2(20,8,0)^\T=(20/13,\,4/13,\,140/13)^\T$,
so $\phi=4/13$. Multipliers $(0,1,0)\mT_2=(-5/13,\,1,\,0)$; reduced costs of
$x_1,x_2,x_4,x_5,x_6$: $(12/13,\,-7/13,\,3/13,\,-8/13,\,18/13)$. The most
negative is $-8/13$: $x_5$ enters with column $\mT_2(1,1,1)^\T=(1/13,\,8/13,\,20/13)^\T$.

| Basic | $\mT_2$ | | | value | $\bar\mA_{\cdot5}$ | ratio |
|---|---|---|---|---|---|---|
| $x_3$ | 1/13 | 0 | 0 | 20/13 | 1/13 | 20 |
| $x_7$ | −5/13 | 1 | 0 | 4/13 | **8/13** | 1/2 |
| $-z$ | 7/13 | 0 | 1 | 140/13 | 20/13 | — |

$x_7$ leaves, and

$$
\begin{gathered}
\mE_2=\begin{pmatrix}1&-1/8&0\\0&13/8&0\\0&-5/2&1\end{pmatrix},\\
\mT_3=\mE_2\mT_2=\begin{pmatrix}1/8&-1/8&0\\-5/8&13/8&0\\3/2&-5/2&1\end{pmatrix}.
\end{gathered}
$$

**End of Phase I.** The values are $\mT_3(20,8,0)^\T=(3/2,\,1/2,\,10)^\T$ for
$(x_3,x_5,-z)$, so $\phi=0$ and both artificial variables are nonbasic: the
Phase I reduced costs are now $(0,0,0,0,0,1,1)$ and $(x_3,x_5)=(3/2,1/2)$ is
a basic feasible solution of the original LP, with $z=-10$. The answer to
"what is $\mB^{-1}$ for Phase II?" can be read off $\mT_3$: its upper-left
block is $\mB^{-1}$ of the basis $(x_3,x_5)$ and its last row is
$-\vpi^\T$ for the original costs,

$$
\begin{gathered}
\mB^{-1}=\begin{pmatrix}1/8&-1/8\\-5/8&13/8\end{pmatrix},\\
\vpi^\T=\left(-\tfrac32,\ \tfrac52\right).
\end{gathered}
$$

The artificial columns are dropped from now on.

**Phase II, iteration 3.** Reduced costs of $x_1,x_2,x_4$, with
$\mA_{\cdot1}=(5,1)^\T$, $\mA_{\cdot2}=(-4,-1)^\T$, $\mA_{\cdot4}=(-2,-1)^\T$:

$$
\begin{aligned}
(\bar c_1,\bar c_2,\bar c_4)&=(3,-1,3)-\left(-\tfrac32,\tfrac52\right)
\begin{pmatrix}5&-4&-2\\1&-1&-1\end{pmatrix}\\
&=\left(8,\,-\tfrac92,\,\tfrac52\right).
\end{aligned}
$$

$x_2$ enters with $\mT_3(-4,-1,-1)^\T=(-3/8,\,7/8,\,-9/2)^\T$; the last entry
is $\bar c_2$, as {eq}`eq-03-06-T-column` promises.

| Basic | $\mB^{-1}$ | | value | $\bar\mA_{\cdot2}$ | ratio |
|---|---|---|---|---|---|
| $x_3$ | 1/8 | −1/8 | 3/2 | −3/8 | — |
| $x_5$ | −5/8 | 13/8 | 1/2 | **7/8** | 4/7 |

$x_5$ leaves with $\theta=4/7$. On the $2\times2$ block the eta matrix is
$\begin{pmatrix}1&3/7\\0&8/7\end{pmatrix}$; on the bordered $3\times3$
inverse it is the same with the extra row $(0,\,36/7,\,1)$, and

$$
\begin{gathered}
\mT_4=\begin{pmatrix}-1/7&4/7&0\\-5/7&13/7&0\\-12/7&41/7&1\end{pmatrix},\\
\mT_4\begin{pmatrix}20\\8\\0\end{pmatrix}=\begin{pmatrix}12/7\\4/7\\88/7\end{pmatrix}.
\end{gathered}
$$

**Phase II, iteration 4.** Basis $(x_3,x_2)$ with values $(12/7,\,4/7)$ and
$-z=88/7$. From the last row of $\mT_4$, $\vpi^\T=(12/7,\,-41/7)$, and

$$
(\bar c_1,\bar c_4,\bar c_5)=\left(\tfrac27,\ \tfrac47,\ \tfrac{36}{7}\right)\ge0 .
$$

The basis is optimal:

$$
\begin{gathered}
\vx^*=\left(0,\,\tfrac47,\,\tfrac{12}{7},\,0,\,0\right)^\T,\\
z^*=-\tfrac{88}{7}\approx-12.57 .
\end{gathered}
$$

Check: $5\cdot0-4\cdot\tfrac47+13\cdot\tfrac{12}{7}=\tfrac{140}{7}=20$ and
$-\tfrac47+5\cdot\tfrac{12}{7}=8$; the objective is
$-\tfrac47-7\cdot\tfrac{12}{7}=-\tfrac{88}{7}$, and
$\vpi^\T\vb=\tfrac{12}{7}\cdot20-\tfrac{41}{7}\cdot8=-\tfrac{88}{7}$. HiGHS
returns the same optimum.

*(Earlier versions of this example list the starting basic variables as
$x_5,x_6,-z$ instead of $x_6,x_7,-z$, use the cost $-7$ and the column
$(2,-1)^\T$ for $x_4$ in the Phase II pricing, and stop without stating the
optimal solution. The cost of $x_4$ is 3 and its column is $(-2,-1)^\T$;
the reduced costs $(8,-9/2,5/2)$ are correct for these data.)*
:::

(sec-03-06-dual-simplex)=
## The dual simplex method

A basis is **primal feasible** when $\bar\vb=\mB^{-1}\vb\ge\vzero$ and
**dual feasible** when all reduced costs are optimal, $\bar c_j\ge0$ for a
minimization. Dual feasibility means that $\vpi^\T=\vc_B^\T\mB^{-1}$
satisfies $\vpi^\T\mA_{\cdot j}\le c_j$ for every $j$, that is, $\vpi$ is a
feasible solution of the dual LP $\max\{\vb^\T\vpi:\mA^\T\vpi\le\vc\}$
({prf:ref}`def-03-05-dual-feasible-basis`). A basis that is both is optimal
({prf:ref}`prop-03-05-basis-optimal`). The primal simplex
method keeps primal feasibility and works toward dual feasibility; the
**dual simplex method** (Lemke, 1954) does the opposite:

* it starts from a basis that is dual feasible but primal infeasible
  (some $\bar b_r<0$);
* it pivots so that the basis stays dual feasible;
* it stops when the basis becomes primal feasible, hence optimal, or when a
  row proves that the LP has no feasible solution.

Such starting bases arise naturally: a $\ge$ constraint with nonnegative
costs gives a slack basis with negative right-hand sides
({prf:ref}`ex-08-04-dual-infeasible`), and after a right-hand side changes or a
constraint is added, the old optimal basis keeps its reduced costs but may
lose feasibility ({doc}`03-07-sensitivity-analysis`).

### Deriving the pivot rules

Take a dual feasible tableau with $\bar b_r<0$ in row $r$, whose basic
variable $x_{B(r)}$ must leave (it is negative). With $\mathcal B$ the index
set of the basic variables, row $r$ reads

$$
x_{B(r)}+\sum_{j\notin\mathcal B}\bar a_{rj}x_j=\bar b_r<0 .
$$ (eq-03-06-pivot-row)

**Infeasibility test.** If $\bar a_{rj}\ge0$ for every nonbasic $j$, the
left-hand side of {eq}`eq-03-06-pivot-row` is nonnegative for every
$\vx\ge\vzero$, while the right-hand side is negative. No feasible $\vx$
exists.

**Entering variable.** Otherwise pivot on some $\bar a_{rs}<0$; a negative
pivot is needed so that the new value $\bar b_r/\bar a_{rs}$ of $x_s$ is
positive. The pivot subtracts $\bar c_s/\bar a_{rs}$ times row $r$ from the
$-z$ row. Writing $\theta_D=-\bar c_s/\bar a_{rs}\ge0$, the new reduced
costs are

$$
\bar c_j'=\bar c_j+\theta_D\,\bar a_{rj}\qquad(j\notin\mathcal B),
$$ (eq-03-06-dual-update)

and the leaving variable gets $\bar c'=\theta_D$. For $\bar a_{rj}\ge0$
the reduced cost cannot fall. For $\bar a_{rj}<0$ it stays nonnegative
exactly when $\theta_D\le-\bar c_j/\bar a_{rj}$. Dual feasibility is kept if
and only if $\theta_D$ is the **dual ratio**

$$
\theta_D=\min\left\{\frac{\bar c_j}{-\bar a_{rj}}:\ j\notin\mathcal B,\ \bar a_{rj}<0\right\},
$$ (eq-03-06-dual-ratio)

and the entering variable $x_s$ is one that attains it. (A larger step
would make the reduced cost of the minimizing variable negative.)

**Objective.** The objective entry changes from $-z$ to
$-z+\theta_D\bar b_r$, so $z$ increases by $\theta_D\abs{\bar b_r}\ge0$.
Since $z=\vc_B^\T\bar\vb=\vpi^\T\vb$ is the dual objective of the dual
feasible $\vpi$, the method increases the dual objective monotonically: it
is the primal simplex method applied to the dual LP, carried out on the
primal tableau. If every pivot has $\theta_D>0$ (dual nondegeneracy), no
basis repeats and the method stops after finitely many pivots.

For a maximization, dual feasibility is $\bar c_j\le0$, and the same
argument gives the ratio $\min\{\abs{\bar c_j/\bar a_{rj}}:\bar a_{rj}<0\}$;
the absolute-value form covers both cases.

:::{prf:algorithm} Dual simplex method (minimization)
:label: alg-03-06-dual-simplex

**Input:** a basis $\mB$ that is dual feasible, $\bar c_j\ge0$ for all
nonbasic $j$.
**Output:** an optimal basis, or a proof that the LP is infeasible.

1. **Optimality test.** If $\bar\vb=\mB^{-1}\vb\ge\vzero$, stop: the basis
   is optimal.
2. **Leaving variable.** Choose a row $r$ with $\bar b_r<0$ (usually the
   most negative). Its basic variable leaves.
3. **Infeasibility test.** Compute row $r$ of $\mB^{-1}\mA$,
   $\bar a_{rj}=(\ve_r^\T\mB^{-1})\mA_{\cdot j}$ for nonbasic $j$. If all
   $\bar a_{rj}\ge0$, stop: the LP has no feasible solution.
4. **Entering variable.** Choose $x_s$ attaining the dual ratio
   {eq}`eq-03-06-dual-ratio`.
5. **Pivot.** Pivot on $\bar a_{rs}$ (in tableau form), or in revised form
   compute $\bar\mA_{\cdot s}=\mB^{-1}\mA_{\cdot s}$ and set
   $\mB^{-1}\leftarrow\mE\mB^{-1}$ with the eta matrix of the pivot on
   $(r,s)$. Go to Step 1.
:::

In revised form, Step 3 needs only the row vector $\ve_r^\T\mB^{-1}$, row
$r$ of the inverse, and one inner product per nonbasic column; Step 5
needs one column. The multipliers $\vpi$ are updated or recomputed as in
{prf:ref}`alg-03-06-revised-simplex` to obtain the reduced costs. A fully
worked revised run is {prf:ref}`ex-08-04-a5-q4-revised-dual`.

:::{prf:example} One dual simplex pivot
:label: ex-03-06-dual-simplex

A minimization problem has the canonical tableau below, with basic
variables $x_1,x_2,x_3$ and objective value 0. Re-optimize it.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|---|---|---|---|---|---|---|---|
| $x_1$ | 1 | 0 | 0 | 4 | −5 | 7 | 8 |
| $x_2$ | 0 | 1 | 0 | −2 | 4 | −2 | −2 |
| $x_3$ | 0 | 0 | 1 | 1 | −3 | 2 | 2 |
| $-z$ | 0 | 0 | 0 | 1 | 3 | 2 | 0 |
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

All reduced costs $(1,3,2)$ of the nonbasic variables are nonnegative, so
the basis is dual feasible, but $x_2=-2<0$: the basis is primal infeasible
and the dual simplex method applies.

*Leaving variable.* Row 2 is the only row with a negative right-hand side,
so $x_2$ leaves.

*Entering variable.* The negative entries of row 2 are $\bar a_{24}=-2$ and
$\bar a_{26}=-2$. Adding $\theta_D$ times row 2 to the $-z$ row changes
the reduced costs of $x_4$ and $x_6$ to $1-2\theta_D$ and $2-2\theta_D$;
both stay nonnegative only for $\theta_D\le\min\{1/2,\,1\}=1/2$. This is
the dual ratio {eq}`eq-03-06-dual-ratio`: $1/2$ for $x_4$ and $2/2=1$ for
$x_6$. So $x_4$ enters, and the pivot is $\bar a_{24}=-2$.

*Pivot.*

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|---|---|---|---|---|---|---|---|
| $x_1$ | 1 | 2 | 0 | 0 | 3 | 3 | 4 |
| $x_4$ | 0 | −1/2 | 0 | 1 | −2 | 1 | 1 |
| $x_3$ | 0 | 1/2 | 1 | 0 | −1 | 1 | 1 |
| $-z$ | 0 | 1/2 | 0 | 0 | 5 | 1 | −1 |

The right-hand side is now nonnegative and the reduced costs
$(1/2,5,1)$ are still nonnegative: the tableau is optimal with
$x_1=4$, $x_3=1$, $x_4=1$, the other variables 0, and $z^*=1$. The
objective rose from 0 to 1, which is $\theta_D\abs{\bar b_2}=\tfrac12\cdot2$.
Had $x_6$ entered instead (ratio 1), its pivot would have made the reduced
cost of $x_4$ equal to $1-2\cdot1=-1<0$, destroying dual feasibility;
this is why the rule takes negative entries only and the minimum ratio.
Solving the tableau's system with HiGHS gives the same optimum.
:::

(sec-03-06-adding-constraint)=
## Adding a constraint and re-optimizing

Suppose an LP has been solved and a new constraint $\va^\T\vx\le\beta$ must
be added (a new restriction, a cut, or a constraint that was left out to
save work, see {prf:ref}`ex-08-04-a6-q4-secondary`). If the optimal $\vx^*$
satisfies it, $\vx^*$ remains optimal, because adding a constraint can only
shrink the feasible set. Otherwise:

1. **Append the row.** Add a slack variable $x_{n+1}\ge0$, so that
   $\va^\T\vx+x_{n+1}=\beta$, and make $x_{n+1}$ basic in the new row.
2. **Restore canonical form.** The new row contains basic variables of the
   old tableau. Subtract $a_j$ times the row of each basic $x_j$ from it.
   The result is the row of $x_{n+1}$, with right-hand side
   $\beta-\va^\T\vx^*<0$.
3. **Re-optimize.** The reduced costs are unchanged (the slack has cost 0
   and the old rows are unchanged), so the enlarged basis is dual feasible
   and primal infeasible: apply {prf:ref}`alg-03-06-dual-simplex`, starting
   with the new row as pivot row.

In matrix terms the new basis matrix and its inverse are

$$
\begin{gathered}
\mB'=\begin{pmatrix}\mB&\vzero\\\va_B^\T&1\end{pmatrix},\\
\mB'^{-1}=\begin{pmatrix}\mB^{-1}&\vzero\\-\va_B^\T\mB^{-1}&1\end{pmatrix},
\end{gathered}
$$ (eq-03-06-new-row)

with the same bordered structure as {eq}`eq-03-06-bordered-inverse`; the new
row of $\mB'^{-1}\mA$ is $\va^\T-\va_B^\T\mB^{-1}\mA$, which is Step 2. An
equality constraint is handled as two inequalities, or by adding an
artificial variable to its row. Usually a few dual simplex pivots suffice,
far fewer than a new solve. {prf:ref}`ex-03-07-new-constraint` applies the
procedure to the running example of {doc}`03-07-sensitivity-analysis`.

(sec-03-06-artificial-constraint)=
## The artificial-constraint technique

The dual simplex method needs a dual feasible starting basis. When the
obvious slack basis is neither primal nor dual feasible, for example a
maximization with positive costs and $\ge$ constraints, a dual feasible
basis can be manufactured with one extra constraint.

Let $Q$ be the set of nonbasic variables with a wrong-signed reduced cost
($\bar c_j<0$ in a minimization, $\bar c_j>0$ in a maximization), or simply
all nonbasic variables. Add the **artificial constraint**

$$
\sum_{j\in Q}x_j+x_{n+1}=M,
$$ (eq-03-06-artificial-constraint)

where $M$ is a constant so large that the constraint cuts off no relevant
feasible point, and $x_{n+1}$ is its slack, basic with value $M$. Pivot in
this row on the column $x_s$ with the most wrong-signed reduced cost
($\bar c_s=\min_{j\in Q}\bar c_j<0$ for a minimization). The pivot element
is 1, and the new reduced costs are

$$
\begin{gathered}
\bar c_j'=\bar c_j-\bar c_s\ \ (j\in Q),\\
\bar c_{n+1}'=-\bar c_s>0,
\end{gathered}
$$

while those outside $Q$ do not change. Since $\bar c_s$ is the smallest,
every $\bar c_j'\ge0$: the basis is dual feasible. Now run the dual simplex
method, carrying $M$ as a symbol and comparing expressions $aM+b$
by $a$ first and $b$ second, as with the big-M method of
{ref}`sec-03-04-big-m`.

At the end, either the slack $x_{n+1}$ is basic and positive, the artificial
constraint is not binding, and the optimum is the optimum of the original
LP; or the artificial constraint is binding. In the second case, if the
optimal value depends on $M$ (it grows without bound with $M$), the original
LP is unbounded; if it does not, the reduced cost of $x_{n+1}$ is zero, the
value is the optimal value of the original LP, and a primal pivot that brings
$x_{n+1}$ into the basis gives an optimal solution free of $M$. If the dual simplex method stops with
an infeasible row, the original LP is infeasible. The technique is worked
in {prf:ref}`ex-08-04-a5-q5-artificial-constraint`.

(sec-03-06-primal-vs-dual)=
## Primal versus dual simplex

```{table} The primal and the dual simplex methods side by side (minimization).
:name: tab-03-06-primal-vs-dual

| | Primal simplex | Dual simplex |
|---|---|---|
| Starting basis | primal feasible, $\bar\vb\ge\vzero$ | dual feasible, $\bar c_j\ge0$ |
| Kept at every pivot | primal feasibility | dual feasibility |
| Goal | dual feasibility | primal feasibility |
| Chosen first | entering variable ($\bar c_s<0$) | leaving variable ($\bar b_r<0$) |
| Ratio test | over $\bar a_{is}>0$ in the entering column | over $\bar a_{rj}<0$ in the leaving row |
| Objective $z$ | does not increase | does not decrease |
| Failure test | column $\bar\mA_{\cdot s}\le\vzero$: unbounded | row $\bar\mA_{r\cdot}\ge\vzero$ with $\bar b_r<0$: infeasible |
| Typical use | solving from a slack or Phase I basis | re-optimizing after changes of $\vb$ or new constraints |
```

The two failure tests are dual to each other: an unbounded primal means an
infeasible dual and vice versa ({numref}`tab-03-05-status`). In sensitivity
analysis the choice is mechanical: a change that destroys dual feasibility
(a cost, a new variable) is repaired by primal simplex pivots; a change that
destroys primal feasibility (a right-hand side, a new constraint) is
repaired by dual simplex pivots.

## Summary and Common Mistakes

### Summary

- A simplex iteration needs only $\vpi^\T=\vc_B^\T\mB^{-1}$, the reduced
  costs, one column $\mB^{-1}\mA_{\cdot s}$ and $\bar\vb=\mB^{-1}\vb$
  ({eq}`eq-03-06-matrix-quantities`); the revised method computes these from
  the original data ({prf:ref}`alg-03-06-revised-simplex`).
- A pivot multiplies the inverse on the left by an eta matrix
  ({prf:ref}`def-03-06-eta`, {prf:ref}`thm-03-06-update`), so
  $\mB_k^{-1}=\mE_{k-1}\cdots\mE_1$ ({eq}`eq-03-06-product-form`).
- Per iteration the revised method costs about $3m^2+m(n-m)$ multiplications
  for dense data against about $mn$ for the tableau; sparsity, partial
  pricing and accuracy make it the method of choice in practice.
- The bordered inverse $\mT$ ({eq}`eq-03-06-bordered-inverse`) carries
  $\mB^{-1}$ and $-\vpi^\T$ together and returns whole tableau columns with
  their reduced costs ({eq}`eq-03-06-T-column`); kept through Phase I, it
  starts Phase II without extra work ({prf:ref}`ex-03-06-two-phase-revised`).
- The dual simplex method keeps $\bar c\ge\vzero$, removes a negative basic
  variable, and enters the variable with the smallest dual ratio
  ({eq}`eq-03-06-dual-ratio`, {prf:ref}`alg-03-06-dual-simplex`); a pivot row
  with no negative entry proves infeasibility.
- New constraints are added by bordering the basis ({eq}`eq-03-06-new-row`)
  and re-optimizing with the dual simplex method; an artificial bounding
  constraint ({eq}`eq-03-06-artificial-constraint`) creates a dual feasible
  start when none is available.

### Common mistakes

- Updating $\mB^{-1}$ by multiplying on the right, $\mB^{-1}\mE$; the eta
  matrix multiplies on the left, $\mE\mB^{-1}$.
- Forgetting the sign of the eta column: off-pivot entries are
  $-\bar a_{is}/\bar a_{rs}$, the pivot entry is $1/\bar a_{rs}$.
- Including the row of $-z$ (or of $\phi$) in the ratio test of the
  bordered method; $-z$ is a free variable.
- Pricing $-z$ with cost 1 in Phase I; its Phase I cost is 0.
- In the dual simplex method, taking a positive $\bar a_{rj}$ in the ratio
  test, or the largest ratio: either destroys dual feasibility.
- Starting the dual simplex method from a basis that is not dual feasible;
  use the artificial-constraint technique first.
- Concluding infeasibility from a row with $\bar b_r<0$ that still has a
  negative entry; only a row with all entries nonnegative proves it.

## References

- Bazaraa, M. S., Jarvis, J. J., & Sherali, H. D. (2010). *Linear programming
  and network flows* (4th ed.). Wiley. (Section 5.1: the revised simplex
  method; Sections 6.4 and 6.6: the dual simplex method and the
  artificial-constraint technique.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Section 3.3: implementations of the
  simplex method, revised simplex and operation counts; Section 4.5: the
  dual simplex method.)
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Chapters 8–9: pivoting and inverses in product form, and
  the simplex method using multipliers, that is, the revised simplex method.)
- Dantzig, G. B., & Orchard-Hays, W. (1954). The product form for the
  inverse in the simplex method. *Mathematical Tables and Other Aids to
  Computation, 8*(46), 64–67.
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 5: the matrix form of the
  simplex method and the revised simplex method; Section 8.1: the dual
  simplex method.)
- Lemke, C. E. (1954). The dual method of solving the linear programming
  problem. *Naval Research Logistics Quarterly, 1*(1), 36–47.
  <https://doi.org/10.1002/nav.3800010107>
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Section 4.4: the dual simplex method; Section 7.2: the revised
  simplex method.)
