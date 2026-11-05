# Artificial Variables: Two-Phase and Big-M Methods

## Overview

The simplex method of {doc}`03-03-simplex-method` starts from a linear
program (LP) in
canonical form with respect to a feasible basis. When every constraint is
of the form $\le$ with a nonnegative right-hand side, the slack variables
provide such a basis. Constraints of the form $\ge$ or $=$ do not: a surplus
variable has coefficient $-1$, and an equation has no extra variable at all.
This chapter shows how to start the simplex method anyway. **Artificial
variables** are added to the rows without a basic variable, which gives an
identity basis for an enlarged problem; the work is then to drive them to
zero. The **two-phase method** does this in a separate **Phase I**, which
minimizes the sum $\phi$ of the artificial variables. Its outcome either
proves that the LP is infeasible, or delivers a feasible basis of the
original problem from which **Phase II** optimizes the true objective. The
**big-M method** does both in one run by charging each artificial variable a
very large cost $M$. The chapter derives both methods, settles what can
happen at the end of Phase I (artificial variables left in the basis at zero
level, redundant constraints), solves one example by the two-phase method
and another by both methods, and compares them.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- explain why artificial variables are needed, where they are needed, and
  why a solution of the original LP corresponds to artificial variables at
  zero;
- set up and price out the Phase I problem $\min\phi=\sum_iR_i$, and run
  Phase I while carrying the original objective row along;
- interpret the end of Phase I: infeasibility, an artificial variable
  basic at zero level, a redundant constraint, variables that must stay at
  zero;
- run the big-M method with $M$ kept symbolic, comparing coefficients of
  the form $aM+b$; and
- compare the two methods and avoid their common pitfalls.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $R_i$ | artificial variable added to constraint $i$; $\vr$ the vector of artificial variables |
| $\mathcal R$ | index set of the constraints that receive an artificial variable |
| $\phi$ | Phase I objective $\phi=\sum_{i\in\mathcal R}R_i$, with optimal value $\phi^*$ |
| $-\phi$ | label of the Phase I objective row of a tableau (reduced costs of $\phi$; right-hand side $-\phi$) |
| $M$ | big-M penalty, a symbolic number larger than any number it is compared with |
| $aM+b$ | an objective-row entry of the big-M method, stored as the pair $(a,b)$ |

(sec-03-04-artificial)=
## Why artificial variables are needed

Bring the LP to standard form $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$
with $\vb\ge\vzero$; rows with a negative right-hand side are multiplied by
$-1$ first. The simplex method needs a basis $\mB$ with
$\mB^{-1}\vb\ge\vzero$, and the easiest one to invert is an identity matrix.
Look at each constraint row:

- A $\le$ row has a slack variable with coefficient $+1$ that appears in no
  other row. It can be the basic variable of the row, with value
  $b_i\ge0$.
- A $\ge$ row has a surplus variable with coefficient $-1$. Making it basic
  would give it the value $-b_i<0$ (for $b_i>0$), which is infeasible.
- An equality row has no extra variable.

More generally, a row needs nothing extra if some variable has a positive
coefficient in it and a zero coefficient in every other row; after the row is
divided by that coefficient, the variable can be the row's basic variable (if
it has a nonzero cost, the objective row is priced out as in
{eq}`eq-03-04-phase-one-row` below). To every remaining
row $i\in\mathcal R$ we add an **artificial variable** $R_i\ge0$ with
coefficient $+1$:

$$
\begin{gathered}
\sum_{j=1}^n a_{ij}x_j+R_i=b_i\quad(i\in\mathcal R),\\
\sum_{j=1}^n a_{ij}x_j=b_i\quad(i\notin\mathcal R).
\end{gathered}
$$ (eq-03-04-augmented)

The basic variables of the rows outside $\mathcal R$ (usually their slack
variables) and the artificial variables form an identity basis. Its basic
feasible solution (BFS) gives
each of these basic variables the right-hand side $b_i$ of its row, so in
particular $\vr=\vb_{\mathcal R}$, and all other variables are zero.

The artificial variables are not part of the model: they measure by how
much the starting point violates a constraint. The enlarged system
{eq}`eq-03-04-augmented` and the original system are related as follows.

:::{prf:lemma} Artificial variables at zero
:label: lem-03-04-artificial

A vector $\vx$ is feasible for the original constraints
$\mA\vx=\vb,\ \vx\ge\vzero$ if and only if $(\vx,\vr)=(\vx,\vzero)$ is
feasible for {eq}`eq-03-04-augmented`. Moreover a BFS of
{eq}`eq-03-04-augmented` with $\vr=\vzero$ whose basis contains no artificial
variable is a BFS of the original constraints.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

With $\vr=\vzero$ the two systems coincide. For the second statement, the
basis columns of such a BFS are columns of $\mA$, linearly independent, and
the nonbasic variables are zero.
:::

So the original LP is feasible exactly when the enlarged system has a
solution with **all** artificial variables equal to zero. Since $R_i\ge0$,
all of them are zero if and only if their sum is zero. This suggests
minimizing that sum.

(sec-03-04-phase-one)=
## Phase I

The **Phase I problem** is

$$
\begin{aligned}
\min_{\vx,\vr}\ & \phi=\sum_{i\in\mathcal R}R_i\\
\st\ & \textstyle\sum_{j}a_{ij}x_j+R_i=b_i,\qquad i\in\mathcal R,\\
& \textstyle\sum_{j}a_{ij}x_j=b_i,\qquad i\notin\mathcal R,\\
& \vx,\ \vr\ge\vzero .
\end{aligned}
$$ (eq-03-04-phase-one)

Its constraints are those of {eq}`eq-03-04-augmented`. It always has an
optimal solution: the starting BFS makes it feasible, and $\phi\ge0$ bounds
it below, so the simplex method (with an anticycling rule such as Bland's,
{prf:ref}`thm-03-03-bland`, if needed) ends with an optimal basis.

**Pricing out.** The objective $\phi$ has cost $1$ on the artificial
variables, which are basic, so the row $(0,\dots,0,1,\dots,1\mid0)$ is not
yet in canonical form. Subtracting every artificial row from it makes the
artificial columns zero:

$$
\begin{gathered}
\bar c^{\,\phi}_j=-\sum_{i\in\mathcal R}a_{ij}\quad(j=1,\dots,n),\\
\text{right-hand side}=-\sum_{i\in\mathcal R}b_i=-\phi .
\end{gathered}
$$ (eq-03-04-phase-one-row)

The entries $\bar c^{\,\phi}_j$ are the Phase I reduced costs; a variable
with a large total coefficient in the artificial rows has a very negative
entry and is attractive, because raising it reduces several violations at
once.

**Carrying the objective along.** The original objective row ($-z$) is kept
in the tableau during Phase I and is updated by every pivot like any other
row, but it plays no part in choosing pivots. When Phase I ends it is
already in canonical form for the basis found, so Phase II can start at once.
Some texts write Phase I as $\max(-\phi)$, with a "$-w$" row of the
opposite sign; the pivots are the same.

(sec-03-04-end-phase-one)=
### The end of Phase I

:::{prf:theorem} Outcomes of Phase I
:label: thm-03-04-phase-one

Let Phase I end with an optimal basis and optimal value $\phi^*$.

1. If $\phi^*>0$, the original LP is **infeasible**.
2. If $\phi^*=0$ and no artificial variable is basic, the final basis is a
   feasible basis of the original LP.
3. If $\phi^*=0$ and an artificial variable $R_k$ is basic (necessarily at
   value $0$) in row $i$, then either
   (a) some non-artificial variable $x_j$ has $\bar a_{ij}\neq0$, and a pivot
   on $\bar a_{ij}$ removes $R_k$ from the basis without changing the
   solution; or
   (b) $\bar a_{ij}=0$ for every non-artificial $x_j$, and the original
   constraint $k$ is **redundant**: it is a linear combination of the other
   constraints and can be deleted together with row $i$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

1. If the original LP had a feasible $\vx$, then $(\vx,\vzero)$ would be
   feasible for {eq}`eq-03-04-phase-one` with $\phi=0<\phi^*$, contradicting
   optimality.
2. All artificial variables are nonbasic, hence zero; apply
   {prf:ref}`lem-03-04-artificial`.
3. (a) Row $i$ has $\bar b_i=R_k=0$, so a pivot on $\bar a_{ij}$ has step
   length $\theta=0$ whatever the sign of $\bar a_{ij}$: no basic variable
   changes value, the solution stays feasible and $\phi$ stays $0$, and the
   new basis is again a basis because $\bar a_{ij}\neq 0$. (b) Row $i$ of
   the tableau, restricted to the original columns, reads
   $0=0$. It equals $\mB^{-1}$-row $i$ times $[\mA\mid\vb]$, a linear
   combination of the original rows in which row $k$ has a nonzero
   coefficient (the coefficient of $R_k$ in row $i$ is $1$). Hence row $k$
   of $[\mA\mid\vb]$ is a combination of the other rows.
:::

Case 3(a) is repeated until no artificial variable is basic or only rows of
type 3(b) remain; deleting those rows leaves an $\mA$ of full row rank.
In case 3(a) any nonzero entry will do; negative pivots are allowed here
because $\theta=0$.

A useful refinement concerns the nonbasic variables. At the end of Phase I
the $-\phi$ row says that, for every solution of the constraints,

$$
\phi=\phi^*+\sum_{j\ \text{nonbasic}}\bar c^{\,\phi}_jx_j .
$$

If $\phi^*=0$, every feasible solution of the original problem has $\phi=0$,
so every nonbasic variable with $\bar c^{\,\phi}_j>0$ must be **zero in
every feasible solution** of the original LP. Such variables, and all
artificial variables, can be deleted before Phase II.
{prf:ref}`ex-08-02-a3-q10` uses this to see that Phase II is not even needed.

(sec-03-04-phase-two)=
## Phase II and the two-phase algorithm

In **Phase II** the artificial columns and the $-\phi$ row are deleted, and
the simplex method continues from the final basis of Phase I with the
original objective row, which is already in canonical form.

:::{prf:algorithm} Two-phase simplex method
:label: alg-03-04-two-phase

**Input:** an LP in standard form with $\vb\ge\vzero$.
**Output:** an optimal BFS, or a proof of infeasibility, or an unbounded ray.

1. **Artificial variables.** Add $R_i$ to every row $i$ that has no
   variable able to serve as its basic variable; the slack variables of the
   other rows and the $R_i$ form the starting basis.
2. **Phase I set-up.** Write the $-z$ row of the original objective and the
   $-\phi$ row of $\phi=\sum_iR_i$ below the constraint rows. Price out the
   $-\phi$ row by {eq}`eq-03-04-phase-one-row` (and the $-z$ row, if a basic
   variable has nonzero cost).
3. **Phase I.** Minimize $\phi$ by {prf:ref}`alg-03-03-simplex`, choosing pivots
   from the $-\phi$ row and updating the $-z$ row with every pivot.
4. **Feasibility test.** If $\phi^*>0$, stop: the LP is infeasible.
5. **Clean-up.** Pivot every artificial variable that is still basic (at
   zero level) out of the basis on a nonzero entry of a non-artificial
   column; delete the rows in which no such entry exists (redundant
   constraints).
6. **Phase II.** Delete the artificial columns and the $-\phi$ row
   (optionally also the nonbasic variables with positive $-\phi$ entries).
   Apply {prf:ref}`alg-03-03-simplex` with the $-z$ row.
:::

:::{prf:example} Two-phase method with a degenerate artificial variable
:label: ex-03-04-two-phase

Solve

$$
\begin{aligned}
\max_{\vx}\ & z=-3x_1+x_3\\
\st\ & x_1+x_2+x_3+x_4=4,\\
& -2x_1+x_2-x_3=1,\\
& 3x_2+x_3+x_4=9,\\
& x_1,x_2,x_3,x_4\ge0 .
\end{aligned}
$$

by the two-phase method.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Set-up.** All three constraints are equations. The only candidate,
$x_4$, occurs in two rows and cannot be the basic variable of either, so no
row has a ready basic variable, and three artificial variables
$R_1,R_2,R_3$ are needed. The Phase I objective is $\phi=R_1+R_2+R_3$, whose priced-out row
{eq}`eq-03-04-phase-one-row` is $-(1-2+0,\ 1+1+3,\ 1-1+1,\ 1+0+1)=(1,-5,-1,-2)$
with right-hand side $-(4+1+9)=-14$. The starting BFS violates the
constraints by a total of $\phi=14$.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $R_1$ | $R_2$ | $R_3$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $R_1$ | $1$ | $1$ | $1$ | $1$ | $1$ | $0$ | $0$ | $4$ |
| $R_2$ | $-2$ | $\boldsymbol{1}$ | $-1$ | $0$ | $0$ | $1$ | $0$ | $1$ |
| $R_3$ | $0$ | $3$ | $1$ | $1$ | $0$ | $0$ | $1$ | $9$ |
| $-z$ | $-3$ | $0$ | $1$ | $0$ | $0$ | $0$ | $0$ | $0$ |
| $-\phi$ | $1$ | $-5$ | $-1$ | $-2$ | $0$ | $0$ | $0$ | $-14$ |

**Phase I, iteration 1.** The most negative entry of the $-\phi$ row is
$-5$, so $x_2$ enters. The ratios are $4/1=4$, $1/1=1$ and $9/3=3$; $R_2$
leaves, and $\phi$ falls by $1\cdot5$ to $9$.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $R_1$ | $R_2$ | $R_3$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $R_1$ | $\boldsymbol{3}$ | $0$ | $2$ | $1$ | $1$ | $-1$ | $0$ | $3$ |
| $x_2$ | $-2$ | $1$ | $-1$ | $0$ | $0$ | $1$ | $0$ | $1$ |
| $R_3$ | $6$ | $0$ | $4$ | $1$ | $0$ | $-3$ | $1$ | $6$ |
| $-z$ | $-3$ | $0$ | $1$ | $0$ | $0$ | $0$ | $0$ | $0$ |
| $-\phi$ | $-9$ | $0$ | $-6$ | $-2$ | $0$ | $5$ | $0$ | $-9$ |

**Phase I, iteration 2.** Now $x_1$ enters ($-9$). The ratios $3/3=1$ and
$6/6=1$ **tie**. By the smallest-index convention $R_1$ leaves, and $R_3$
stays basic with the value $6-6\cdot1=0$: the new BFS is degenerate.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $R_1$ | $R_2$ | $R_3$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $\tfrac{2}{3}$ | $\tfrac{1}{3}$ | $\tfrac{1}{3}$ | $-\tfrac{1}{3}$ | $0$ | $1$ |
| $x_2$ | $0$ | $1$ | $\tfrac{1}{3}$ | $\tfrac{2}{3}$ | $\tfrac{2}{3}$ | $\tfrac{1}{3}$ | $0$ | $3$ |
| $R_3$ | $0$ | $0$ | $0$ | $\boldsymbol{-1}$ | $-2$ | $-1$ | $1$ | $0$ |
| $-z$ | $0$ | $0$ | $3$ | $1$ | $1$ | $-1$ | $0$ | $3$ |
| $-\phi$ | $0$ | $0$ | $0$ | $1$ | $3$ | $2$ | $0$ | $0$ |

All entries of the $-\phi$ row are $\ge0$ and $\phi^*=0$: the LP is
feasible. The BFS $(x_1,x_2,x_3,x_4)=(1,3,0,0)$ satisfies all three
equations ($1+3=4$, $-2+3=1$, $9=9$).

**Clean-up.** The artificial $R_3$ is still basic, at zero level. Its row
has the nonzero entry $-1$ under $x_4$, so case 3(a) of
{prf:ref}`thm-03-04-phase-one` applies: pivot on $-1$. The step length is
$0/(-1)$, nothing changes value, and $x_4$ replaces $R_3$ at level $0$:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $R_1$ | $R_2$ | $R_3$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $\tfrac{2}{3}$ | $0$ | $-\tfrac{1}{3}$ | $-\tfrac{2}{3}$ | $\tfrac{1}{3}$ | $1$ |
| $x_2$ | $0$ | $1$ | $\tfrac{1}{3}$ | $0$ | $-\tfrac{2}{3}$ | $-\tfrac{1}{3}$ | $\tfrac{2}{3}$ | $3$ |
| $x_4$ | $0$ | $0$ | $0$ | $1$ | $2$ | $1$ | $-1$ | $0$ |
| $-z$ | $0$ | $0$ | $3$ | $0$ | $-1$ | $-2$ | $1$ | $3$ |
| $-\phi$ | $0$ | $0$ | $0$ | $0$ | $1$ | $1$ | $1$ | $0$ |

The degenerate row had a meaning. Before the clean-up it read
$-x_4-2R_1-R_2+R_3=0$, and with $\vr=\vzero$ it says $x_4=0$. Indeed,
(third constraint) $-2\times$(first) $-$ (second) gives $-x_4=0$: the
constraints force $x_4=0$, which is also what the positive entry $1$ of
$x_4$ in the $-\phi$ row of the previous tableau announced. No constraint is
redundant (the matrix $\mA$ has rank 3), so case 3(b) does not occur.

**Phase II.** Delete the artificial columns and the $-\phi$ row. The $-z$
row, carried along, reads $(0,0,3,0\mid3)$: the current value is $z=-3$
(from $x_1=1$), and $\bar c_3=3>0$, so for this maximization $x_3$ enters.
The ratios are $1/\tfrac23=\tfrac32$ (row $x_1$) and $3/\tfrac13=9$ (row
$x_2$); the third row has $\bar a=0$. So $x_1$ leaves with
$\theta=\tfrac32$, and $z$ grows by $\tfrac32\cdot3=\tfrac92$.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $\boldsymbol{\tfrac{2}{3}}$ | $0$ | $1$ |
| $x_2$ | $0$ | $1$ | $\tfrac{1}{3}$ | $0$ | $3$ |
| $x_4$ | $0$ | $0$ | $0$ | $1$ | $0$ |
| $-z$ | $0$ | $0$ | $3$ | $0$ | $3$ |

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $\tfrac{3}{2}$ | $0$ | $1$ | $0$ | $\tfrac{3}{2}$ |
| $x_2$ | $-\tfrac{1}{2}$ | $1$ | $0$ | $0$ | $\tfrac{5}{2}$ |
| $x_4$ | $0$ | $0$ | $0$ | $1$ | $0$ |
| $-z$ | $-\tfrac{9}{2}$ | $0$ | $0$ | $0$ | $-\tfrac{3}{2}$ |

All reduced costs are $\le0$, so the last tableau is optimal:
$\vx^*=(0,\tfrac52,\tfrac32,0)$ with $z^*=\tfrac32$. Check:
$0+\tfrac52+\tfrac32+0=4$, $0+\tfrac52-\tfrac32=1$,
$3\cdot\tfrac52+\tfrac32+0=9$, and $-3\cdot0+\tfrac32=\tfrac32$; an LP
solver gives the same optimum. Altogether the method needed two Phase I
pivots, one clean-up pivot and one Phase II pivot.
:::

If the clean-up pivot is skipped and $R_3$ is simply kept in the basis during
Phase II, it must not be allowed to become positive. In
{prf:ref}`ex-03-04-two-phase` the $x_3$ column has a zero in the $R_3$ row, so
nothing would go wrong; in general one must either pivot it out first or delete
the variables with positive $-\phi$ entries, as in
{ref}`sec-03-04-end-phase-one`.

(sec-03-04-big-m)=
## The big-M method

The big-M method (Charnes, 1952) combines the two phases into one LP. Every
artificial variable gets a large cost in the original objective,

$$
\begin{aligned}
\min_{\vx,\vr}\ & z_M=\vc^\T\vx+M\sum_{i\in\mathcal R}R_i\\
\st\ & \textstyle\sum_{j}a_{ij}x_j+R_i=b_i,\qquad i\in\mathcal R,\\
& \textstyle\sum_{j}a_{ij}x_j=b_i,\qquad i\notin\mathcal R,\\
& \vx,\ \vr\ge\vzero,
\end{aligned}
$$ (eq-03-04-big-m)

that is, the constraints {eq}`eq-03-04-augmented`, where $M>0$ is "very
large". For a maximization problem the penalty is
$-M\sum_iR_i$. As long as an artificial variable is positive, the penalty
outweighs anything the true objective can gain, so the simplex method
drives the artificial variables out first and then optimizes
$\vc^\T\vx$.

**Keeping $M$ symbolic.** Instead of choosing a number, treat $M$ as larger
than every number it is compared with. Every objective-row entry then has
the form $aM+b$ and is stored as the pair $(a,b)$; pivots act on both parts
separately, and comparisons are lexicographic:

$$
\begin{aligned}
aM+b<a'M+b'\iff{}& a<a'\\
&\text{or}\ (a=a'\ \text{and}\ b<b') .
\end{aligned}
$$ (eq-03-04-compare-m)

So $aM+b$ is negative when $a<0$, or when $a=0$ and $b<0$. As with Phase I,
the objective row must first be **priced out**: the artificial variables
are basic but have cost $M$, so $M$ times every artificial row is
subtracted from the row $(\vc^\T, M,\dots,M\mid0)$.

:::{prf:example} Big-M and two-phase methods on the same LP
:label: ex-03-04-big-m

Solve

$$
\begin{aligned}
\min_{\vx}\ & z=4x_1+x_2\\
\st\ & 3x_1+x_2=3,\\
& 4x_1+3x_2\ge6,\\
& x_1+2x_2\le4,\\
& x_1,x_2\ge0 .
\end{aligned}
$$

by the big-M method and by the two-phase method, and compare the two runs.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Set-up.** With the surplus variable $x_3$ and the slack variable $x_4$,
the constraints are $3x_1+x_2=3$, $4x_1+3x_2-x_3=6$ and $x_1+2x_2+x_4=4$.
The third row has its basic variable $x_4$; the first two receive $R_1$ and
$R_2$. The starting BFS is $R_1=3$, $R_2=6$, $x_4=4$. *(An earlier
version of this example gives $x_4=6$ for this starting solution; the third
constraint gives $x_4=4$.)* The objective row
$(4,1,0,0,M,M\mid0)$ is priced out by subtracting $M$ times the rows of
$R_1$ and $R_2$:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $R_1$ | $R_2$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $R_1$ | $\boldsymbol{3}$ | $1$ | $0$ | $0$ | $1$ | $0$ | $3$ |
| $R_2$ | $4$ | $3$ | $-1$ | $0$ | $0$ | $1$ | $6$ |
| $x_4$ | $1$ | $2$ | $0$ | $1$ | $0$ | $0$ | $4$ |
| $-z$ | $-7M+4$ | $-4M+1$ | $M$ | $0$ | $0$ | $0$ | $-9M$ |

The current value is $z_M=9M$: the artificial variables add up to 9.

**Iteration 1.** By {eq}`eq-03-04-compare-m` the most negative entry is
$-7M+4$, so $x_1$ enters. The ratios are $3/3=1$, $6/4=\tfrac32$ and
$4/1=4$; $R_1$ leaves. The value becomes $z_M=2M+4$.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $R_1$ | $R_2$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $\tfrac{1}{3}$ | $0$ | $0$ | $\tfrac{1}{3}$ | $0$ | $1$ |
| $R_2$ | $0$ | $\boldsymbol{\tfrac{5}{3}}$ | $-1$ | $0$ | $-\tfrac{4}{3}$ | $1$ | $2$ |
| $x_4$ | $0$ | $\tfrac{5}{3}$ | $0$ | $1$ | $-\tfrac{1}{3}$ | $0$ | $3$ |
| $-z$ | $0$ | $-\tfrac{5}{3}M-\tfrac{1}{3}$ | $M$ | $0$ | $\tfrac{7}{3}M-\tfrac{4}{3}$ | $0$ | $-2M-4$ |

**Iteration 2.** Only $x_2$ has a negative entry. The ratios are
$1/\tfrac13=3$, $2/\tfrac53=\tfrac65$ and $3/\tfrac53=\tfrac95$; $R_2$
leaves. Both artificial variables are now nonbasic, and the $M$-parts of
all entries of the non-artificial columns vanish.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $R_1$ | $R_2$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $\tfrac{1}{5}$ | $0$ | $\tfrac{3}{5}$ | $-\tfrac{1}{5}$ | $\tfrac{3}{5}$ |
| $x_2$ | $0$ | $1$ | $-\tfrac{3}{5}$ | $0$ | $-\tfrac{4}{5}$ | $\tfrac{3}{5}$ | $\tfrac{6}{5}$ |
| $x_4$ | $0$ | $0$ | $\boldsymbol{1}$ | $1$ | $1$ | $-1$ | $1$ |
| $-z$ | $0$ | $0$ | $-\tfrac{1}{5}$ | $0$ | $M-\tfrac{8}{5}$ | $M+\tfrac{1}{5}$ | $-\tfrac{18}{5}$ |

The point $(\tfrac35,\tfrac65)$, with $z=\tfrac{18}{5}$, is feasible for the
original LP but not optimal.

**Iteration 3.** $x_3$ enters ($-\tfrac15$). Its column has the positive
entries $\tfrac15$ (row $x_1$) and $1$ (row $x_4$), with ratios
$\tfrac35/\tfrac15=3$ and $1/1=1$, so $x_4$ leaves with $\theta=1$ and $z$
falls by $\tfrac15$.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $R_1$ | $R_2$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $0$ | $-\tfrac{1}{5}$ | $\tfrac{2}{5}$ | $0$ | $\tfrac{2}{5}$ |
| $x_2$ | $0$ | $1$ | $0$ | $\tfrac{3}{5}$ | $-\tfrac{1}{5}$ | $0$ | $\tfrac{9}{5}$ |
| $x_3$ | $0$ | $0$ | $1$ | $1$ | $1$ | $-1$ | $1$ |
| $-z$ | $0$ | $0$ | $0$ | $\tfrac{1}{5}$ | $M-\tfrac{7}{5}$ | $M$ | $-\tfrac{17}{5}$ |

All entries are $\ge0$ for large $M$, so the tableau is optimal:
$x_1^*=\tfrac25$, $x_2^*=\tfrac95$, surplus $x_3^*=1$, $x_4^*=0$, and
$z^*=\tfrac{17}{5}$. Check: $3\cdot\tfrac25+\tfrac95=3$,
$4\cdot\tfrac25+3\cdot\tfrac95=7\ge6$, $\tfrac25+2\cdot\tfrac95=4\le4$, and
$4\cdot\tfrac25+\tfrac95=\tfrac{17}{5}$.

**The same LP by the two-phase method.** Phase I minimizes
$\phi=R_1+R_2$; its priced-out row is $(-7,-4,1,0,0,0\mid-9)$, exactly the
$M$-part of the first big-M row, and the original row $(4,1,0,0,0,0\mid0)$
is its constant part. The Phase I pivots are therefore the same as the
first two big-M pivots ($x_1$ for $R_1$, then $x_2$ for $R_2$), and Phase I
ends with $\phi^*=0$ at

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $R_1$ | $R_2$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $\tfrac{1}{5}$ | $0$ | $\tfrac{3}{5}$ | $-\tfrac{1}{5}$ | $\tfrac{3}{5}$ |
| $x_2$ | $0$ | $1$ | $-\tfrac{3}{5}$ | $0$ | $-\tfrac{4}{5}$ | $\tfrac{3}{5}$ | $\tfrac{6}{5}$ |
| $x_4$ | $0$ | $0$ | $1$ | $1$ | $1$ | $-1$ | $1$ |
| $-z$ | $0$ | $0$ | $-\tfrac{1}{5}$ | $0$ | $-\tfrac{8}{5}$ | $\tfrac{1}{5}$ | $-\tfrac{18}{5}$ |
| $-\phi$ | $0$ | $0$ | $0$ | $0$ | $1$ | $1$ | $0$ |

After the artificial columns and the $-\phi$ row are deleted, Phase II makes
the same third pivot ($x_3$ for $x_4$) and reaches the same optimum
$(\tfrac25,\tfrac95)$, $z^*=\tfrac{17}{5}$.

```{figure} figures/03-04-phase-one-path.svg
:name: fig-03-04-phase-one-path
:width: 520px
:class: dark-light
:alt: The feasible set of the example is a short segment of the line 3x1+x2=3; the simplex path goes from the origin to (1,0) outside the feasible set, then to (3/5,6/5) and along the segment to (2/5,9/5)

The path of both methods in {prf:ref}`ex-03-04-big-m`, drawn in the
$(x_1,x_2)$-plane. The feasible set is the thick segment of the line
$3x_1+x_2=3$ between $(\tfrac25,\tfrac95)$ and $(\tfrac35,\tfrac65)$.
The first two iterates, $(0,0)$ and $(1,0)$, are BFSs of the enlarged
problem with $R_1+R_2>0$ and lie outside it (hollow markers); the end of
Phase I, $(\tfrac35,\tfrac65)$, is the first feasible point, and Phase II
moves along the segment to the optimum $(\tfrac25,\tfrac95)$, where the
contour $4x_1+x_2=\tfrac{17}{5}$ (dashed) touches it.
```
:::

### How large must $M$ be?

With $M$ symbolic, the method is exactly as reliable as the two-phase method.
With a number for $M$ it is reliable only if the number is large enough. In the
solution of {prf:ref}`ex-03-04-big-m` the final reduced cost of $R_1$ is
$M-\tfrac75$; with $M<\tfrac75$ the final basis would not be optimal for
{eq}`eq-03-04-big-m`. Indeed, solving {eq}`eq-03-04-big-m` with $M=1$ gives
$\vx=(0,2)$ with $R_1=1$ and $z_M=3<\tfrac{17}{5}$: the solver prefers to
violate the first constraint and pay the penalty. Conversely, the final tableau
is optimal for every numerical $M\ge\tfrac75$, and for $M>\tfrac75$ it is the
only optimal basis (a solver confirms this for $M=1.5,2,10$). The threshold
$\tfrac75$ is the absolute value of the optimal dual price of the first
constraint ({ref}`sec-03-05-reading-duals`): the penalty for violating a
constraint must exceed the value of violating it. In general that value is
unknown in advance, which is the main weakness of a numerical $M$; a huge $M$,
on the other hand, mixes numbers of very different magnitude in one row and
causes round-off errors in floating-point arithmetic.

**Infeasibility.** If the big-M method with symbolic $M$ ends with an optimal tableau in which
some artificial variable is still positive (the objective value has a
positive $M$-part), the original LP is infeasible; this is the analog of
$\phi^*>0$. If it stops with an unbounded column while some artificial
variable is still positive, the original LP is either infeasible or
unbounded, and a Phase I computation is needed to tell which; this is one
more reason to prefer the two-phase method.

(sec-03-04-comparison)=
## Comparison of the two methods

```{list-table} Two-phase and big-M methods compared
:name: tab-03-04-comparison
:header-rows: 1
:stub-columns: 1

* - 
  - Two-phase method
  - Big-M method
* - Objective while artificial variables are positive
  - $\phi=\sum_iR_i$ only
  - $\vc^\T\vx+M\sum_iR_i$
* - Objective rows carried
  - two ($-z$ and $-\phi$)
  - one, with entries $aM+b$
* - Choice of a parameter
  - none
  - $M$ (symbolic, or a number that must be large enough)
* - Infeasibility test
  - $\phi^*>0$
  - an artificial variable positive at the optimum
* - Artificial variable basic at zero level
  - clean-up pivots or deletion of a redundant row before Phase II
  - remains in the basis with penalty $M$; harmless as long as it stays zero
* - Pivots
  - the same as big-M with symbolic $M$ while artificials are positive
  - the same
* - Use in software
  - standard (Phase I is also used to find feasible points in general)
  - mainly in hand computation and as an idea (penalty functions)
```

With $M$ kept symbolic, the $M$-parts of the big-M objective row are exactly
the $-\phi$ row and the constant parts are the $-z$ row, and the
lexicographic comparison {eq}`eq-03-04-compare-m` lets the $M$-part decide
first, as Phase I does. While some artificial variable can still be reduced,
both methods make the same pivots (up to ties), as
{prf:ref}`ex-03-04-big-m` shows. The differences are bookkeeping and numerics.
The two-phase method needs no parameter and makes the end of Phase I an
explicit checkpoint at which infeasibility, redundant constraints and
forced-zero variables are detected. The big-M idea of replacing hard
constraints by heavy penalties returns in the penalty methods of
{doc}`06-06-constrained-algorithms` ({ref}`sec-06-06-penalty`) and in the big-M
constraints of {doc}`07-01-logical-constraints` ({ref}`sec-07-01-big-m`; a different
use of the same letter).

## Summary and Common Mistakes

### Summary

- Rows of type $\ge$ or $=$ have no natural basic variable; artificial
  variables $R_i$ create an identity basis for the enlarged system
  {eq}`eq-03-04-augmented`, whose solutions with $\vr=\vzero$ are exactly the
  feasible solutions of the LP ({prf:ref}`lem-03-04-artificial`).
- Phase I minimizes $\phi=\sum_iR_i$ ({eq}`eq-03-04-phase-one`); its row is
  priced out by {eq}`eq-03-04-phase-one-row`, and the original objective row is
  carried along.
- At the end of Phase I: $\phi^*>0$ proves infeasibility; $\phi^*=0$ gives a
  feasible basis after artificial variables at zero level are pivoted out
  or redundant rows deleted ({prf:ref}`thm-03-04-phase-one`); nonbasic
  variables with positive $-\phi$ entries are zero in every feasible
  solution.
- Phase II continues with the carried $-z$ row
  ({prf:ref}`alg-03-04-two-phase`); {prf:ref}`ex-03-04-two-phase` reaches
  $\vx^*=(0,\tfrac52,\tfrac32,0)$, $z^*=\tfrac32$.
- The big-M method {eq}`eq-03-04-big-m` penalizes artificial variables by $M$
  and compares entries $aM+b$ lexicographically ({eq}`eq-03-04-compare-m`); on
  {prf:ref}`ex-03-04-big-m` it takes the same path as the two-phase method
  ({numref}`fig-03-04-phase-one-path`) to $(\tfrac25,\tfrac95)$, $z^*=\tfrac{17}{5}$.
- A numerical $M$ must exceed the dual prices of the penalized rows
  ($M>\tfrac75$ in the example); the two-phase method avoids the choice
  ({numref}`tab-03-04-comparison`).

### Common mistakes

- **Forgetting to price out** the $-\phi$ row or the big-M row. A row with
  nonzero entries under basic variables does not hold reduced costs, and
  the pivots chosen from it are wrong.
- **Adding artificial variables where none are needed.** Rows with a slack
  variable, or with any variable that occurs only in that row with a
  positive coefficient, have a basic variable already. Fewer artificial
  variables mean fewer Phase I pivots ({prf:ref}`ex-08-02-a3-q11` needs only
  one).
- **Negative right-hand sides.** Multiply such rows by $-1$ before deciding
  which rows need artificial variables; otherwise the starting basis is
  infeasible.
- **Ending Phase I too early.** Phase I ends when the $-\phi$ row passes the
  optimality test, not when the first artificial variable leaves.
- **Declaring infeasibility when $\phi^*=0$ with an artificial variable in
  the basis.** That artificial variable is zero; the LP is feasible
  (case 3 of {prf:ref}`thm-03-04-phase-one`).
- **Letting a zero-level artificial variable become positive in Phase II**,
  for example by entering a variable with a nonzero entry in its row
  without first pivoting it out.
- **Choosing a numerical $M$ that is too small** (wrong optimum or false
  infeasibility) or absurdly large (round-off).
- **Restarting Phase II with the original costs** instead of the $-z$ row
  carried through Phase I: the original row is not in canonical form for
  the new basis.

## References

- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Section 3.5: finding an initial basic
  feasible solution, the two-phase method, driving artificial variables out
  of the basis, redundant equality constraints, the big-M method.)
- Charnes, A. (1952). Optimality and degeneracy in linear programming.
  *Econometrica, 20*(2), 160–170. <https://doi.org/10.2307/1907845>
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Chapter 5: Phase I and Phase II of the simplex method.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Section 4.6: adapting to other model
  forms, the big-M and two-phase methods.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Section 3.4: artificial starting solution, the M-method and the
  two-phase method; Chapter 3, special cases of the simplex method: infeasible
  solution.)
