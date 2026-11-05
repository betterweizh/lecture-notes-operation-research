# Duality Theory

## Overview

Every linear program comes with a second linear program built from the same
data, its **dual**. Where the original (primal) problem chooses activity
levels, the dual chooses prices for the resources or requirements, and the two
problems bound each other: any feasible price vector gives a lower bound on
the cost of a minimization problem, and any feasible plan gives an upper bound
on the value of the prices. This chapter derives the dual from a pricing
argument, gives the rules for writing the dual of any linear program, and
proves the main theorems: **weak duality** (the bound), **strong duality**
(the bounds meet at the optimum, proved with the simplex method), and
**complementary slackness** (which constraints and variables can be positive
together). It shows that the optimality test of the simplex method is a test
of dual feasibility, how to read the optimal dual values from a final tableau,
and why these values are **shadow prices**: rates at which the optimal
objective changes when a right-hand side changes. The chapter closes with
Farkas' lemma, a theorem of the alternative that follows from duality and is
used again for the Karush–Kuhn–Tucker (KKT) conditions of nonlinear
programming. It builds on
the bases and the fundamental theorem of {doc}`03-02-lp-geometry` and on the
tableaux and reduced costs of {doc}`03-03-simplex-method` and
{doc}`03-04-two-phase-big-m`; the Lagrangian view of
{ref}`sec-03-05-dual-rules` is developed for general convex problems in
{doc}`02-03-convex-optimization-duality`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- explain the dual of a linear program as a pricing problem and write the
  dual of any linear program with the table of rules;
- prove weak duality and use it to bound optimal values and to certify
  optimality;
- state strong duality, prove it from an optimal simplex basis, and list
  which primal/dual status combinations can occur;
- state and prove complementary slackness and use it to recover a dual
  solution from a primal one and conversely;
- decide whether a basis is primal feasible, dual feasible or optimal, and
  read the optimal dual values from a final tableau;
- interpret dual values as shadow prices and imputed costs, and state when
  this interpretation is valid; and
- derive Farkas' lemma from linear programming duality.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\vpi\in\R^m$, $\pi_i$ | dual variables, one per primal constraint (prices) |
| $w=\vb^\T\vpi$ | dual objective value; $w^*$ its optimal value |
| $\mA_{i\cdot}$, $\mA_{\cdot j}$ | row $i$ and column $j$ of $\mA$; $\va_i=\mA_{i\cdot}^\T$ is the normal of constraint $i$ in figures |
| $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ | reduced cost of $x_j$; for a dual feasible $\vpi$ also the **slack of dual constraint** $j$ |
| $s_i=\mA_{i\cdot}\vx-b_i$ | surplus of primal constraint $i$ ($\ge0$ for a feasible $\vx$ in the canonical form) |
| $t_j$ | slack variable of dual constraint $j$ in a dual tableau |
| $z^*(\vb)$ | optimal value as a function of the right-hand side |

(sec-03-05-diet-pricing)=
## A pricing problem behind every linear program

The idea of duality is easiest to see in a small cost-minimization problem.

:::{prf:example} A diet and the pill-maker's prices
:label: ex-03-05-diet

A household wants the cheapest mix of six foods that supplies at least 9 units
of vitamin A and 19 units of vitamin B. {numref}`tab-03-05-diet-data` gives the
vitamin content of one kilogram of each food and its price.

```{list-table} Data of the diet problem
:name: tab-03-05-diet-data
:header-rows: 1
:stub-columns: 1

* - Food $j$
  - 1
  - 2
  - 3
  - 4
  - 5
  - 6
  - Requirement (units)
* - Vitamin A (units/kg)
  - 1
  - 0
  - 2
  - 2
  - 1
  - 2
  - 9
* - Vitamin B (units/kg)
  - 0
  - 1
  - 3
  - 1
  - 3
  - 2
  - 19
* - Price $c_j$ (cents/kg)
  - 35
  - 30
  - 60
  - 50
  - 27
  - 22
  -
```

A manufacturer offers vitamin A and vitamin B as pills and wants the
household to switch from foods to pills. It sets a price of $\pi_1$ cents per
unit of vitamin A and $\pi_2$ cents per unit of vitamin B. Formulate the
household's problem of buying the cheapest diet, and the manufacturer's
problem of choosing prices that maximize its revenue while no food is a
cheaper source of the vitamins than the pills.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

With $x_j$ the kilograms of food $j$ bought, the household solves

$$
\begin{aligned}
\min_{\vx}\ & z=35x_1+30x_2+60x_3+50x_4+27x_5+22x_6\\
\st\ & x_1+2x_3+2x_4+x_5+2x_6\ge 9\\
& x_2+3x_3+x_4+3x_5+2x_6\ge 19\\
& x_1,\dots,x_6\ge 0,
\end{aligned}
$$ (eq-03-05-diet-primal)

that is, $\min\{\vc^\T\vx:\mA\vx\ge\vb,\ \vx\ge\vzero\}$ with
$\mA=\begin{pmatrix}1&0&2&2&1&2\\0&1&3&1&3&2\end{pmatrix}$, $\vb=(9,19)^\T$
and $\vc=(35,30,60,50,27,22)^\T$.

Given pill prices $\pi_1$ and $\pi_2$, the household will compare the pills
with each food: one kilogram of food 3 delivers 2 units of A and 3 units of B,
worth $2\pi_1+3\pi_2$ cents at pill prices, and costs 60 cents. If
$2\pi_1+3\pi_2>60$, food 3 is a cheaper source of these vitamins and the
household keeps buying it. The prices are **competitive** only if no food is a
cheaper source than the pills, that is, if $\pi_1a_{1j}+\pi_2a_{2j}\le c_j$ for
every food $j$. A household that switches buys exactly its requirements, 9
units of A and 19 units of B, so the manufacturer maximizes its revenue
$9\pi_1+19\pi_2$ over competitive prices:

$$
\begin{aligned}
\max_{\vpi}\ & w=9\pi_1+19\pi_2\\
\st\ & \pi_1\le 35\\
& \pi_2\le 30\\
& 2\pi_1+3\pi_2\le 60\\
& 2\pi_1+\pi_2\le 50\\
& \pi_1+3\pi_2\le 27\\
& 2\pi_1+2\pi_2\le 22\\
& \pi_1,\pi_2\ge 0,
\end{aligned}
$$ (eq-03-05-diet-dual)

that is, $\max\{\vb^\T\vpi:\mA^\T\vpi\le\vc,\ \vpi\ge\vzero\}$.

Both problems are solved later in this chapter. The household's cheapest diet
is $x_5=5$ kg of food 5 and $x_6=2$ kg of food 6 at $z^*=179$ cents, and the
manufacturer's best prices are $\vpi^*=(3,8)$ cents per unit, with revenue
$w^*=9\cdot3+19\cdot8=179$ cents.
:::

The pricing problem {eq}`eq-03-05-diet-dual` of {prf:ref}`ex-03-05-diet` is
the **dual** of the diet problem {eq}`eq-03-05-diet-primal`.
Several features of the pair hold for every primal–dual pair:

1. the constraint matrix of one problem is the transpose of the other's;
2. the right-hand sides of one problem are the objective coefficients of the
   other;
3. each primal constraint has a dual variable (a price per unit of the
   requirement) and each primal variable has a dual constraint (a price
   limit for that food);
4. a minimization problem with $\ge$ constraints is paired with a maximization
   problem with $\le$ constraints; and
5. the manufacturer can never earn more than the cheapest diet costs. Any diet
   $\vx$ that meets the requirements and any competitive prices $\vpi$ satisfy
   $9\pi_1+19\pi_2\le\sum_j c_jx_j$, because the diet contains at least the
   required vitamins and each kilogram of food is worth at least its vitamin
   content at pill prices. The theorems below make this precise and show that
   the two optimal values are in fact equal.

(sec-03-05-canonical-pair)=
## The canonical dual pair

:::{prf:definition} Canonical primal–dual pair
:label: def-03-05-canonical-pair

Let $\mA\in\R^{m\times n}$, $\vb\in\R^m$ and $\vc\in\R^n$. The **primal**
problem (P) and its **dual** problem (D) are

$$
\begin{aligned}
\text{(P)}\quad\min_{\vx}\ & z=\vc^\T\vx\\
\st\ & \mA\vx\ge\vb\\
& \vx\ge\vzero\\[1ex]
\text{(D)}\quad\max_{\vpi}\ & w=\vb^\T\vpi\\
\st\ & \mA^\T\vpi\le\vc\\
& \vpi\ge\vzero
\end{aligned}
$$ (eq-03-05-canonical-pair)

The dual has one variable $\pi_i$ for each primal constraint
$\mA_{i\cdot}\vx\ge b_i$ and one constraint $\vpi^\T\mA_{\cdot j}\le c_j$ for
each primal variable $x_j$. Read from right to left, the same definition
gives the dual of a maximization problem with $\le$ constraints:
$\max\{\vc^\T\vx:\mA\vx\le\vb,\vx\ge\vzero\}$ has the dual
$\min\{\vb^\T\vpi:\mA^\T\vpi\ge\vc,\vpi\ge\vzero\}$.
:::

The canonical pair is a definition, but the dual of every other linear
program follows from it, because every linear program can be rewritten in the
canonical form. For the standard form of {doc}`03-02-lp-geometry`, write each
equation $\mA\vx=\vb$ as the two inequalities $\mA\vx\ge\vb$ and
$-\mA\vx\ge-\vb$. With dual variables $\vpi^+\ge\vzero$ and
$\vpi^-\ge\vzero$ for the two blocks, the canonical dual is
$\max\{\vb^\T(\vpi^+-\vpi^-):\mA^\T(\vpi^+-\vpi^-)\le\vc\}$, and only the
difference $\vpi=\vpi^+-\vpi^-$ matters, which can take any sign. Hence the
standard-form problem (P) has the dual (D) with free variables:

$$
\begin{aligned}
\text{(P)}\quad\min_{\vx}\ & \vc^\T\vx\\
\st\ & \mA\vx=\vb\\
& \vx\ge\vzero\\[1ex]
\text{(D)}\quad\max_{\vpi}\ & \vb^\T\vpi\\
\st\ & \mA^\T\vpi\le\vc\\
& \vpi\ \text{free.}
\end{aligned}
$$ (eq-03-05-standard-pair)

A linear program can be rewritten in many equivalent ways, so it is natural to
ask whether the dual depends on the way the primal is written. It does not.

:::{prf:proposition} Duals of equivalent programs are equivalent
:label: prop-03-05-equivalent-duals

If (P$'$) is obtained from (P) by adding slack or surplus variables, by
writing a free variable as a difference of two nonnegative variables, by
multiplying a constraint by $-1$ or by replacing an equation by two
inequalities, then the dual (D$'$) of (P$'$) is equivalent to the dual (D)
of (P): the two have the same feasible prices (up to the same change of
variables) and the same objective values.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

We show the case used most often: (P) is the canonical primal and (P$'$)
adds surplus variables, $\min\{\vc^\T\vx+\vzero^\T\vs:\mA\vx-\mI\vs=\vb,\
\vx\ge\vzero,\ \vs\ge\vzero\}$. (P$'$) is in standard form with constraint
matrix $(\mA\ \ -\mI)$ and cost vector $(\vc,\vzero)$, so by
{eq}`eq-03-05-standard-pair` its dual is

$$
\begin{aligned}
\max_{\vpi}\ & \vb^\T\vpi\\
\st\ & \mA^\T\vpi\le\vc\\
& -\mI\vpi\le\vzero\\
& \vpi\ \text{free}.
\end{aligned}
$$

The constraint $-\mI\vpi\le\vzero$ says $\vpi\ge\vzero$, so (D$'$) is (D). The
other transformations are checked the same way: splitting $x_j=x_j^+-x_j^-$
gives the two dual constraints $\vpi^\T\mA_{\cdot j}\le c_j$ and
$-\vpi^\T\mA_{\cdot j}\le -c_j$, that is, one equation; multiplying a
constraint by $-1$ changes the sign of its dual variable; splitting an
equation was done above.
:::

The proposition also says that the dual is unique up to such rewriting, so we
may speak of *the* dual. Applying the construction twice returns the primal.

:::{prf:theorem} The dual of the dual is the primal
:label: thm-03-05-dual-of-dual

The dual of (D) in {eq}`eq-03-05-canonical-pair` is (P).
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Write (D) as a canonical minimization problem:
$\min\{-\vb^\T\vpi:-\mA^\T\vpi\ge-\vc,\ \vpi\ge\vzero\}$. Its constraint
matrix is $-\mA^\T$, its right-hand side $-\vc$ and its cost vector $-\vb$.
By {prf:ref}`def-03-05-canonical-pair` its dual, with variables $\vx$, is
$\max\{-\vc^\T\vx:(-\mA^\T)^\T\vx\le-\vb,\ \vx\ge\vzero\}$, which is
$\max\{-\vc^\T\vx:\mA\vx\ge\vb,\ \vx\ge\vzero\}$, the same as
$\min\{\vc^\T\vx:\mA\vx\ge\vb,\ \vx\ge\vzero\}$.
:::

The two problems are therefore on an equal footing. Which one is called the
primal is a matter of where we start, and every statement below has a mirror
image with the roles exchanged.

(sec-03-05-dual-rules)=
## Rules for writing the dual of any linear program

Applying {prf:ref}`prop-03-05-equivalent-duals` to every kind of constraint and
variable gives the rules in {numref}`tab-03-05-dual-rules`. A constraint has the
**natural type** when it points the way the objective pushes: $\ge$ in a
minimization problem (the requirement limits how low the cost can go) and
$\le$ in a maximization problem. A natural-type constraint has a nonnegative
dual variable, an equation has a free one and a constraint of the opposite
type a nonpositive one. Symmetrically, a nonnegative primal variable gives a
dual constraint of the natural type for the dual, a free variable gives an
equation and a nonpositive variable a constraint of the opposite type.

```{list-table} Rules for writing the dual. Read left to right for a minimization primal and right to left for a maximization primal.
:name: tab-03-05-dual-rules
:header-rows: 1

* - Minimization problem
  - Maximization problem
* - objective coefficients $\vc$
  - right-hand sides $\vc$
* - right-hand sides $\vb$
  - objective coefficients $\vb$
* - constraint $i$: $\mA_{i\cdot}\vx\ge b_i$
  - variable $\pi_i\ge0$
* - constraint $i$: $\mA_{i\cdot}\vx\le b_i$
  - variable $\pi_i\le0$
* - constraint $i$: $\mA_{i\cdot}\vx=b_i$
  - variable $\pi_i$ free
* - variable $x_j\ge0$
  - constraint $j$: $\vpi^\T\mA_{\cdot j}\le c_j$
* - variable $x_j\le0$
  - constraint $j$: $\vpi^\T\mA_{\cdot j}\ge c_j$
* - variable $x_j$ free
  - constraint $j$: $\vpi^\T\mA_{\cdot j}=c_j$
```

### Why the rules hold: the Lagrangian bound

The rules are not arbitrary; they are exactly the conditions under which the
prices give a valid bound. Take a minimization problem with constraints of the
three types and variables of the three signs, and multiply constraint $i$ by a
price $\pi_i$. For every feasible $\vx$ we want

$$
\pi_i\,(b_i-\mA_{i\cdot}\vx)\le 0 ,
$$

so that adding these terms to the cost can only lower it. For a $\ge$
constraint, $b_i-\mA_{i\cdot}\vx\le0$, so we need $\pi_i\ge0$; for a $\le$
constraint we need $\pi_i\le0$; for an equation the term is zero whatever the
sign of $\pi_i$. With such prices, every feasible $\vx$ satisfies

$$
\begin{aligned}
\vc^\T\vx\ \ge\ L(\vx,\vpi)&=\vc^\T\vx+\vpi^\T(\vb-\mA\vx)\\
&=\vb^\T\vpi+(\vc-\mA^\T\vpi)^\T\vx .
\end{aligned}
$$ (eq-03-05-lagrangian)

The right-hand side still depends on $\vx$. It is at least $\vb^\T\vpi$ for
every feasible $\vx$ exactly when each term $(c_j-\vpi^\T\mA_{\cdot j})x_j$ is
nonnegative for every allowed value of $x_j$: $c_j-\vpi^\T\mA_{\cdot j}\ge0$
when $x_j\ge0$, $\le0$ when $x_j\le0$, and $=0$ when $x_j$ is free. These are
the dual constraints of {numref}`tab-03-05-dual-rules`, and under them
$\vc^\T\vx\ge\vb^\T\vpi$. The dual problem asks for the best (largest) such
lower bound.

The function $L$ is the **Lagrangian** of the linear program, and
$q(\vpi)=\inf_{\vx}L(\vx,\vpi)$, taken over the sign restrictions on $\vx$,
equals $\vb^\T\vpi$ when $\vpi$ satisfies the dual constraints and $-\infty$
otherwise. The dual linear program is therefore the Lagrange dual problem
$\max_{\vpi}q(\vpi)$ of {ref}`sec-02-03-dual-problem`, specialized to
linear functions. A difference from general convex programs is that strong
duality for linear programs ({prf:ref}`thm-03-05-strong-duality`) needs no
constraint qualification such as Slater's condition
({prf:ref}`thm-02-03-slater`).

Two small checks guard against sign errors. First, the dual of the dual must
give back the primal ({prf:ref}`thm-03-05-dual-of-dual`). Second, a constraint of
the wrong type can always be multiplied by $-1$ before applying the canonical
pair; its dual variable then appears with the opposite sign, which is the
$\pi_i\le0$ rule. Worked examples, including a primal with all three kinds of
constraint and two free variables and a primal in block form, are in
{doc}`08-03-examples-duality`.

(sec-03-05-weak-duality)=
## Weak duality

:::{prf:theorem} Weak duality
:label: thm-03-05-weak-duality

If $\vx$ is feasible for (P) and $\vpi$ is feasible for (D) in
{eq}`eq-03-05-canonical-pair`, then

$$
\vc^\T\vx\ \ge\ \vpi^\T\mA\vx\ \ge\ \vb^\T\vpi .
$$ (eq-03-05-weak-chain)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Because $\vx\ge\vzero$ and $\mA^\T\vpi\le\vc$, each term of
$\vc^\T\vx-\vpi^\T\mA\vx=\sum_j(c_j-\vpi^\T\mA_{\cdot j})x_j$ is a product of
two nonnegative numbers, so $\vc^\T\vx\ge\vpi^\T\mA\vx$. Because
$\vpi\ge\vzero$ and $\mA\vx\ge\vb$, each term of
$\vpi^\T\mA\vx-\vb^\T\vpi=\sum_i\pi_i(\mA_{i\cdot}\vx-b_i)$ is also a product
of two nonnegative numbers, so $\vpi^\T\mA\vx\ge\vb^\T\vpi$.
:::

Both sign conditions are needed: $\vx\ge\vzero$ for the first inequality and
$\vpi\ge\vzero$ for the second. In a primal–dual pair written with the rules of
{numref}`tab-03-05-dual-rules`, the same argument works term by term through
{eq}`eq-03-05-lagrangian`, so weak duality holds for every pair, not only for the
canonical one. Adding the two parts of the proof gives the identity

$$
\begin{gathered}
\vc^\T\vx-\vb^\T\vpi=\sum_{j=1}^{n}\bar c_j\,x_j+\sum_{i=1}^{m}\pi_i\,s_i,\\
\bar c_j=c_j-\vpi^\T\mA_{\cdot j},\\
s_i=\mA_{i\cdot}\vx-b_i ,
\end{gathered}
$$ (eq-03-05-duality-gap)

which writes the **duality gap** $\vc^\T\vx-\vb^\T\vpi$ as a sum of products
of a variable of one problem and the slack (or surplus) of the matching
constraint of the other. It is used again for complementary slackness.

{numref}`fig-03-05-weak-duality-line` shows weak duality for the diet problem of
{prf:ref}`ex-03-05-diet`. Three competitive price vectors, $(11,0)$, $(5,6)$ and
$(0,9)$, give revenues of 99, 159 and 171 cents; three feasible diets give
costs of 191, 209 and 243 cents (5 kg of food 6 with 3 kg of food 5; 9.5 kg of
food 6 alone; 9 kg of food 5 alone). Every revenue lies to the left of every
cost.

```{figure} figures/03-05-weak-duality-line.svg
:name: fig-03-05-weak-duality-line
:width: 600px
:class: dark-light
:alt: Number line with dual objective values 99, 159 and 171 on the left, primal values 191, 209 and 243 on the right, and the common optimum 179 between them

Weak duality for the diet problem (cents). Dual values $w=\vb^\T\vpi$ of
feasible prices (circles) never exceed primal costs $z=\vc^\T\vx$ of feasible
diets (squares); the best values of both meet at $z^*=w^*=179$ (diamond).
```

Weak duality has several immediate consequences.

:::{prf:corollary} Bounds and certificates
:label: cor-03-05-optimality-certificate

Let (P) and (D) be a primal–dual pair, (P) a minimization problem.

1. The cost $\vc^\T\vx$ of any primal feasible $\vx$ is an upper bound on the
   optimal value of (D), and the value $\vb^\T\vpi$ of any dual feasible
   $\vpi$ is a lower bound on the optimal value of (P).
2. If $\vx$ is primal feasible, $\vpi$ is dual feasible and
   $\vc^\T\vx=\vb^\T\vpi$, then $\vx$ is optimal for (P) and $\vpi$ is optimal
   for (D).
3. If (P) is unbounded ($z\to-\infty$), then (D) is infeasible. If (D) is
   unbounded ($w\to+\infty$), then (P) is infeasible.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Part 1 is {eq}`eq-03-05-weak-chain` read from either end. For part 2, every
primal feasible $\vx'$ has $\vc^\T\vx'\ge\vb^\T\vpi=\vc^\T\vx$, so $\vx$ is
optimal; the argument for $\vpi$ is the same. For part 3, a dual feasible
$\vpi$ would give the finite lower bound $\vb^\T\vpi$ on every primal cost,
which contradicts $z\to-\infty$; the second sentence is the mirror image.
:::

Part 2 is a practical **optimality certificate**. Someone who claims to have
solved a linear program can prove it by handing over a primal solution and a
dual solution: the checker verifies that the primal solution satisfies the
primal constraints, that the dual solution satisfies the dual constraints, and
that the two objective values are equal. None of these checks requires solving
anything. For the diet problem, $\vx=(0,0,0,0,5,2)$ meets both requirements
($5+2\cdot2=9$ and $3\cdot5+2\cdot2=19$), $\vpi=(3,8)$ satisfies all six price
limits (checked in {prf:ref}`ex-03-05-diet-cs`), and both values are 179 cents,
so both are optimal.

The converse of part 3 does not hold: an infeasible primal does not force an
unbounded dual, because both problems can be infeasible. An example is
{prf:ref}`ex-08-03-both-infeasible`.

(sec-03-05-strong-duality)=
## Strong duality

Weak duality leaves open whether a gap remains between the best primal and the
best dual value. For linear programs it does not.

:::{prf:theorem} Strong duality
:label: thm-03-05-strong-duality

If either problem of a primal–dual pair has an optimal solution, so does the
other, and the two optimal values are equal: $z^*=w^*$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

By {prf:ref}`prop-03-05-equivalent-duals` we may write the primal in standard
form, $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$, whose dual is
{eq}`eq-03-05-standard-pair`. We may also assume that $\mA$ has full row rank $m$,
since redundant equations can be deleted without changing the feasible set
({prf:ref}`thm-03-04-phase-one`); deleting an equation fixes its dual variable at
zero and does not change the dual optimal value.

Suppose the primal has an optimal solution. By the fundamental theorem of
linear programming ({prf:ref}`thm-03-02-fundamental`) and the simplex method
with an anticycling rule ({prf:ref}`thm-03-03-bland`), there is an optimal basis
$\mB$: $\vx_B^*=\mB^{-1}\vb\ge\vzero$, $\vx_N^*=\vzero$, and all reduced
costs are nonnegative,

$$
\bar c_j=c_j-\vc_B^\T\mB^{-1}\mA_{\cdot j}\ge 0\qquad\text{for } j=1,\dots,n .
$$

Define the **simplex multipliers**

$$
\vpi^{*\T}=\vc_B^\T\mB^{-1}.
$$ (eq-03-05-simplex-multipliers)

Then $\bar c_j=c_j-\vpi^{*\T}\mA_{\cdot j}\ge0$ for every $j$, that is,
$\mA^\T\vpi^*\le\vc$, so $\vpi^*$ is feasible for the dual (whose variables
are free). Its value is

$$
\begin{aligned}
\vb^\T\vpi^*&=\vc_B^\T\mB^{-1}\vb\\
&=\vc_B^\T\vx_B^*=\vc^\T\vx^* .
\end{aligned}
$$

By part 2 of {prf:ref}`cor-03-05-optimality-certificate`, $\vpi^*$ is optimal for
the dual and $w^*=z^*$.

If instead the dual has an optimal solution, apply the argument just given to
the dual, written as a minimization problem in standard form; by
{prf:ref}`thm-03-05-dual-of-dual` its dual is the primal, which therefore has an
optimal solution with the same value.
:::

The theorem goes back to von Neumann's work on games; the first published
proof for linear programs is that of Gale et al. (1951). The proof given
here, through an optimal simplex basis, is the one of Dantzig (1963).

The proof does more than show existence: it constructs the optimal dual
solution from the optimal basis, $\vpi^{*\T}=\vc_B^\T\mB^{-1}$, and shows that
the simplex optimality test $\bar c\ge\vzero$ *is* the dual feasibility of
this $\vpi^*$. {ref}`sec-03-05-dual-feasible-basis` develops this observation.

Strong duality together with weak duality decides which combinations of
outcomes are possible. A linear program is either infeasible, unbounded, or
has an optimal solution (fundamental theorem of linear programming), so a pair
has nine conceivable combinations.

```{list-table} Possible outcomes of a primal–dual pair (primal: minimization)
:name: tab-03-05-status
:header-rows: 1
:stub-columns: 1

* - Primal \ Dual
  - optimal solution
  - unbounded ($w\to+\infty$)
  - infeasible
* - optimal solution
  - possible, $z^*=w^*$
  - impossible
  - impossible
* - unbounded ($z\to-\infty$)
  - impossible
  - impossible
  - possible
* - infeasible
  - impossible
  - possible
  - possible
```

The entries "impossible" in the first row and column follow from strong
duality (if one problem has an optimum, so does the other). Primal and dual
cannot both be unbounded, by part 3 of
{prf:ref}`cor-03-05-optimality-certificate`. The three remaining combinations
occur: both optimal (the diet problem); one unbounded and the other
infeasible (the mixed-constraint example {prf:ref}`ex-08-03-mixed-dual` with one
sign of a coefficient changed); and both infeasible
({prf:ref}`ex-08-03-both-infeasible`). Two useful consequences can be read off
the table.

:::{prf:corollary} Feasibility decides the outcome
:label: cor-03-05-feasible-pair

If both problems of a pair are feasible, both have optimal solutions with
equal values. If the primal is feasible and the dual infeasible, the primal is
unbounded; if the dual is feasible and the primal infeasible, the dual is
unbounded.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

If both are feasible, neither is unbounded (part 3 of
{prf:ref}`cor-03-05-optimality-certificate`), so each has an optimal solution, and
the values agree by {prf:ref}`thm-03-05-strong-duality`. If the primal is feasible
and the dual infeasible, the primal cannot have an optimal solution (the dual
would then have one), so it is unbounded. The last sentence is the mirror
image.
:::

(sec-03-05-complementary-slackness)=
## Complementary slackness

The gap identity {eq}`eq-03-05-duality-gap` writes $\vc^\T\vx-\vb^\T\vpi$ as a
sum of nonnegative products. At an optimal pair the gap is zero, so every
product must vanish.

:::{prf:theorem} Complementary slackness
:label: thm-03-05-complementary-slackness

Let $\vx$ be feasible for (P) and $\vpi$ feasible for (D) in
{eq}`eq-03-05-canonical-pair`. Then $\vx$ and $\vpi$ are both optimal if and only
if

$$
\begin{aligned}
\pi_i\,(\mA_{i\cdot}\vx-b_i)&=0 && i=1,\dots,m,\\
(c_j-\vpi^\T\mA_{\cdot j})\,x_j&=0 && j=1,\dots,n.
\end{aligned}
$$ (eq-03-05-cs)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

By {eq}`eq-03-05-duality-gap`, $\vc^\T\vx-\vb^\T\vpi$ is the sum of the
$m+n$ products in {eq}`eq-03-05-cs`, each of which is nonnegative for feasible
$\vx$ and $\vpi$. If both are optimal, strong duality gives
$\vc^\T\vx=\vb^\T\vpi$, so the sum is zero and every product is zero.
Conversely, if every product is zero, the two objective values are equal and
both solutions are optimal by part 2 of
{prf:ref}`cor-03-05-optimality-certificate`.
:::

In words: whenever a constraint of one problem has positive slack at an
optimal pair, the matching variable of the other problem is zero; whenever a
variable is positive, the matching constraint of the other problem is tight.
Three remarks prevent common misreadings.

- The conditions say nothing about a **free** variable, because its matching
  constraint is an equation and has no slack. In the standard-form pair
  {eq}`eq-03-05-standard-pair` only the second line of {eq}`eq-03-05-cs` remains.
- A slack and its matching variable may **both** be zero. Complementary
  slackness says that at least one of the two is zero, not exactly one. (A
  pair in which exactly one of the two is positive for every index is called
  *strictly* complementary; whenever both problems have optimal solutions,
  such a pair exists (Goldman & Tucker, 1956), but finding it is not needed
  here.)
- The theorem applies to any pair of optimal solutions, not only to basic
  ones. When a problem has several optimal solutions, every optimal solution
  of the primal is complementary to every optimal solution of the dual.

The most common use is to recover the solution of one problem from the
solution of the other.

:::{prf:example} Recovering the pill prices from the optimal diet
:label: ex-03-05-diet-cs

The optimal diet of {prf:ref}`ex-03-05-diet` is $\vx^*=(0,0,0,0,5,2)$ kg. Find
the optimal prices.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Since $x_5^*=5>0$ and $x_6^*=2>0$, dual constraints 5 and 6 must be tight:

$$
\begin{aligned}
\pi_1+3\pi_2&=27,\\
2\pi_1+2\pi_2&=22 .
\end{aligned}
$$

Subtracting twice the first equation from the second gives $-4\pi_2=-32$, so
$\pi_2=8$ and $\pi_1=27-24=3$: $\vpi=(3,8)$ cents per unit.

Complementary slackness only gives the candidate; it must still be checked to
be dual feasible. The six price limits are

$$
\begin{aligned}
3&\le35,\\
8&\le30,\\
2\cdot3+3\cdot8=30&\le60,\\
2\cdot3+8=14&\le50,\\
3+3\cdot8=27&\le27,\\
2\cdot3+2\cdot8=22&\le22,
\end{aligned}
$$

and $\vpi\ge\vzero$, so $\vpi$ is feasible. The remaining conditions of
{eq}`eq-03-05-cs` hold as well: both prices are positive, and both vitamin
constraints are tight at $\vx^*$ (surpluses
$s_1=0+2\cdot0+2\cdot0+5+2\cdot2-9=0$ and
$s_2=0+3\cdot0+0+3\cdot5+2\cdot2-19=0$). By
{prf:ref}`thm-03-05-complementary-slackness`, $\vx^*$ and $\vpi^*=(3,8)$ are
optimal, with $z^*=w^*=27\cdot5+22\cdot2=9\cdot3+19\cdot8=179$ cents.

Had the check failed, $\vx^*$ would not have been optimal. The argument also
runs in the other direction: from $\vpi^*=(3,8)$ alone, both vitamin
constraints must be tight (both prices are positive), and only foods 5 and 6
can be bought (only their price limits are tight), so
$x_5+2x_6=9$ and $3x_5+2x_6=19$, which gives $x_5=5$, $x_6=2$.
:::

{numref}`fig-03-05-diet-dual-region` shows the dual problem of the diet
example.
The optimal prices lie at the vertex where the price limits of foods 5 and 6
cross, the two foods in the optimal diet; the limits of foods 3 and 4 do not
touch the optimum (and those of foods 1 and 2 lie outside the window), and
these foods are not bought.

```{figure} figures/03-05-diet-dual-region.svg
:name: fig-03-05-diet-dual-region
:width: 520px
:class: dark-light
:alt: Feasible region of the dual diet problem in the pi1-pi2 plane with the optimal price vector (3, 8) at the intersection of the lines of foods 5 and 6

Feasible prices of the pill manufacturer. Line ($j$) is the price limit
$\vpi^\T\mA_{\cdot j}=c_j$ of food $j$. The revenue $w=9\pi_1+19\pi_2$ grows in
the direction $\vb$ and is largest at $\vpi^*=(3,8)$, where only the limits of
foods 5 and 6 are tight; these are the foods with $x_j^*>0$.
```

### A geometric reading

Complementary slackness has a picture in the primal space. Consider a
maximization problem $\max\{\vc^\T\vx:\mA\vx\le\vb,\ \vx\ge\vzero\}$ whose
optimum $\vx^*$ has $x_j^*>0$ for every $j$, so that the sign constraints are
not active. Then every dual constraint is tight, $\mA^\T\vpi^*=\vc$, and
$\pi_i^*=0$ for every constraint that is not active at $\vx^*$. Writing
$\va_i=\mA_{i\cdot}^\T$ for the outward normal of constraint $i$ and $I(\vx^*)$
for the active constraints, this says

$$
\begin{aligned}
\vc&=\sum_{i\in I(\vx^*)}\pi_i^*\,\va_i,\\
\pi_i^*&\ge 0, \qquad i\in I(\vx^*) :
\end{aligned}
$$ (eq-03-05-cone-condition)

at an optimal vertex, the objective gradient lies in the cone spanned by the
normals of the active constraints, and the optimal dual values are the
weights. This is the linear case of the KKT conditions of
{prf:ref}`thm-06-03-kkt`. The product-mix problem shows it.

:::{prf:example} The dual of the product-mix problem
:label: ex-03-05-product-mix-dual

In the product-mix problem of {prf:ref}`ex-03-01-product-mix`, a molder makes
juice glasses and cocktail glasses; $x_1$ and $x_2$ are the weekly outputs in
hundreds of cases, and the contributions are \$500 and \$450 per hundred cases.
Machine time is 6 and 5 hours per hundred cases with 60 hours available per
week; storage takes 10 and 20 m³ per hundred cases with 150 m³ available; and
at most 8 hundred cases of juice glasses can be sold per week:

$$
\begin{aligned}
\max_{\vx}\ & z=500x_1+450x_2\\
\st\ & 6x_1+5x_2\le 60 &&\text{(machine time, h/week)}\\
& 10x_1+20x_2\le 150 &&\text{(storage, m³)}\\
& x_1\le 8 &&\text{(demand, hundred cases/week)}\\
& x_1,x_2\ge 0 .
\end{aligned}
$$

Its optimum is $\vx^*=(45/7,\,30/7)\approx(6.429,\,4.286)$ with
$z^*=36000/7\approx\$5142.86$ per week ({prf:ref}`ex-03-01-product-mix`).
Write the dual problem, find its optimal solution by complementary
slackness, and check the cone condition {eq}`eq-03-05-cone-condition` at
$\vx^*$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**The dual.** The primal is a maximization problem with $\le$ constraints,
the natural type, so by {numref}`tab-03-05-dual-rules` the dual is a minimization
problem with one nonnegative variable per resource and one $\ge$ constraint per
product:

$$
\begin{aligned}
\min_{\vpi}\ & w=60\pi_1+150\pi_2+8\pi_3\\
\st\ & 6\pi_1+10\pi_2+\pi_3\ge 500\\
& 5\pi_1+20\pi_2\ge 450\\
& \pi_1,\pi_2,\pi_3\ge 0 .
\end{aligned}
$$

**Dual solution by complementary slackness.** Both $x_1^*$ and $x_2^*$ are
positive, so both dual constraints are tight. The demand constraint has slack
$8-45/7=11/7>0$ at $\vx^*$, so $\pi_3^*=0$. The tight constraints become
$6\pi_1+10\pi_2=500$ and $5\pi_1+20\pi_2=450$. Twice the first minus the
second gives $7\pi_1=550$, so

$$
\begin{gathered}
\pi_1^*=\frac{550}{7}\approx 78.57,\\
\pi_2^*=\frac{500-6\cdot\frac{550}{7}}{10}=\frac{20}{7}\approx 2.857,\\
\pi_3^*=0 .
\end{gathered}
$$

These prices are nonnegative and satisfy both dual constraints with equality,
so they are dual feasible, and
$w=60\cdot\frac{550}{7}+150\cdot\frac{20}{7}=\frac{33000+3000}{7}=\frac{36000}{7}=z^*$,
which certifies both solutions as optimal.

**Geometry.** {numref}`fig-03-05-cs-cone` shows the optimal vertex. The active
constraints are machine time and storage, with normals $\va_1=(6,5)$ and
$\va_2=(10,20)$, and indeed

$$
\begin{aligned}
\frac{550}{7}\begin{pmatrix}6\\5\end{pmatrix}+\frac{20}{7}\begin{pmatrix}10\\20\end{pmatrix}
&=\frac{1}{7}\begin{pmatrix}3300+200\\2750+400\end{pmatrix}\\
&=\begin{pmatrix}500\\450\end{pmatrix}=\vc ,
\end{aligned}
$$

so $\vc$ lies in the cone spanned by $\va_1$ and $\va_2$. Any feasible move
away from $\vx^*$ has $\va_1^\T\vd\le0$ and $\va_2^\T\vd\le0$, hence
$\vc^\T\vd=\pi_1^*\va_1^\T\vd+\pi_2^*\va_2^\T\vd\le0$: no feasible direction
improves the objective. The economic meaning of the numbers $550/7$ and $20/7$
is the subject of {ref}`sec-03-05-shadow-prices`.

```{figure} figures/03-05-cs-cone.svg
:name: fig-03-05-cs-cone
:width: 520px
:class: dark-light
:alt: Product-mix feasible region with the optimal vertex; the objective gradient c lies between the normals a1 and a2 of the two active constraints

Complementary slackness at the optimal vertex of the product-mix problem. The
objective gradient $\vc$ lies in the shaded cone spanned by the normals
$\va_1$, $\va_2$ of the active constraints (1) and (2), with weights equal to
the optimal dual values; the inactive demand constraint (3) has weight zero.
The dashed line is the optimal objective contour $z=36000/7$.
```
:::

(sec-03-05-dual-feasible-basis)=
## Dual feasibility of a basis

The proof of strong duality attached a dual vector to every basis. We now
use it for bases that are not (yet) optimal. Consider the standard form
$\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$ with its dual
{eq}`eq-03-05-standard-pair`, and let $\mB$ be a basis, with nonbasic columns
$\mN$. The basic solution is $\vx_B=\mB^{-1}\vb$, $\vx_N=\vzero$; the dual
vector of the basis is $\vpi^\T=\vc_B^\T\mB^{-1}$, the unique solution of
$\mB^\T\vpi=\vc_B$, that is, of the dual constraints of the basic variables
written as equations.

:::{prf:definition} Primal and dual feasible bases
:label: def-03-05-dual-feasible-basis

A basis $\mB$ of the standard-form problem is **primal feasible** if
$\mB^{-1}\vb\ge\vzero$, and **dual feasible** if its multipliers
$\vpi^\T=\vc_B^\T\mB^{-1}$ satisfy all dual constraints, $\mA^\T\vpi\le\vc$.
Since the dual constraints of the basic variables hold with equality by
construction, a basis is dual feasible if and only if all its reduced costs
are nonnegative:

$$
\bar\vc_N^\T=\vc_N^\T-\vc_B^\T\mB^{-1}\mN\ \ge\ \vzero^\T .
$$ (eq-03-05-dual-feasible)

For a maximization problem, dual feasibility means $\bar\vc_N\le\vzero$.
:::

The reduced cost $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ is the slack of dual
constraint $j$ at the multipliers of the basis. The tableau of
{doc}`03-03-simplex-method` therefore carries the dual solution along: in
matrix form, the tableau of basis $\mB$ has the constraint rows
$(\mB^{-1}\mA\ \mid\ \mB^{-1}\vb)$ and the objective row

$$
\bigl(\vc^\T-\vpi^\T\mA\ \mid\ -\vpi^\T\vb\bigr)=\bigl(\bar\vc^\T\ \mid\ -z\bigr),
$$

so its last row lists the dual slacks and its corner holds minus the common
value $\vc_B^\T\mB^{-1}\vb=\vb^\T\vpi$ of the two basic solutions. (The
bordered inverse that maps the initial tableau to this one is the $\mT$-matrix
{eq}`eq-03-06-bordered-inverse`, used again in {doc}`03-07-sensitivity-analysis`.)

:::{prf:proposition} A basis that is primal and dual feasible is optimal
:label: prop-03-05-basis-optimal

If a basis $\mB$ is both primal and dual feasible, then
$\vx_B=\mB^{-1}\vb$, $\vx_N=\vzero$ is optimal for the primal and
$\vpi^\T=\vc_B^\T\mB^{-1}$ is optimal for the dual.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Both vectors are feasible, and
$\vc^\T\vx=\vc_B^\T\mB^{-1}\vb=\vpi^\T\vb$; apply part 2 of
{prf:ref}`cor-03-05-optimality-certificate`.
:::

The basic pair is automatically complementary: basic variables have
$\bar c_j=0$ and nonbasic variables have $x_j=0$, so every product
$\bar c_jx_j$ in {eq}`eq-03-05-cs` vanishes. A basis therefore always gives a
complementary pair, and the only question is which of the two solutions is
feasible. This is the cleanest way to describe the two simplex methods.

- The **primal simplex method** keeps the basis primal feasible
  ($\mB^{-1}\vb\ge\vzero$) and complementary, and pivots until the basis is
  also dual feasible ($\bar\vc\ge\vzero$). A negative reduced cost is a
  violated dual constraint, and choosing it as the entering variable is an
  attempt to repair it.
- The **dual simplex method** ({prf:ref}`alg-03-06-dual-simplex`) keeps the
  basis dual feasible and pivots until it is also primal feasible.

A worked test of a given basis, and a pivot from a dual infeasible basis to
the optimal one, are in {ref}`sec-08-03-feasibility-tests`.

(sec-03-05-reading-duals)=
## Reading the dual solution from an optimal tableau

When the simplex method stops, the optimal dual values are already in the final
tableau; there is no need to solve the dual separately. Two ways of reading
them follow from {eq}`eq-03-05-simplex-multipliers`.

1. **From $\mB^{-1}$.** If the initial tableau contained an identity matrix
   (slack or artificial columns), the same columns of the final tableau
   contain $\mB^{-1}$, and $\vpi^\T=\vc_B^\T\mB^{-1}$ can be computed from it.
2. **From the objective row.** Every column has $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$.
   For a column that is a signed unit vector, this gives one $\pi_i$
   directly, as listed in {numref}`tab-03-05-reading-duals`.

```{list-table} Dual values from the objective row of an optimal tableau ($\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$; $\ve_i$ is the $i$th unit vector)
:name: tab-03-05-reading-duals
:header-rows: 1

* - Column in constraint $i$
  - $\mA_{\cdot j}$
  - $c_j$
  - Entry $\bar c_j$
  - Hence
* - slack of a $\le$ row
  - $\ve_i$
  - $0$
  - $-\pi_i$
  - $\pi_i=-\bar c_j$
* - surplus of a $\ge$ row
  - $-\ve_i$
  - $0$
  - $\pi_i$
  - $\pi_i=\bar c_j$
* - artificial, minimization (cost $+M$)
  - $\ve_i$
  - $M$
  - $M-\pi_i$
  - $\pi_i=M-\bar c_j$
* - artificial, maximization (cost $-M$)
  - $\ve_i$
  - $-M$
  - $-M-\pi_i$
  - $\pi_i=-M-\bar c_j$
```

The rules hold for minimization and maximization alike, because they come
from the single formula $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$. The sign of
$\pi_i$ then follows the rules of {numref}`tab-03-05-dual-rules`: in the
product-mix problem (maximization, $\le$ rows) the objective-row entries of
the three slack columns in the optimal tableau are $-550/7$, $-20/7$ and $0$, so
$\vpi^*=(550/7,\,20/7,\,0)$, as found in {prf:ref}`ex-03-05-product-mix-dual`.
Artificial columns can be read even though their variables have left the
problem: their entries still satisfy the formula.

For the diet problem, the optimal tableau of the primal (in standard form with
surplus variables $x_7$ and $x_8$) and of the dual (in standard form with slack
variables $t_1,\dots,t_6$ for the six price limits) are shown in
{numref}`tab-03-05-diet-primal-tableau` and {numref}`tab-03-05-diet-dual-tableau`.

```{list-table} Optimal tableau of the diet problem (minimization; $x_7$, $x_8$ are the surplus variables of the vitamin constraints)
:name: tab-03-05-diet-primal-tableau
:header-rows: 1
:stub-columns: 1

* - Basic
  - $x_1$
  - $x_2$
  - $x_3$
  - $x_4$
  - $x_5$
  - $x_6$
  - $x_7$
  - $x_8$
  - RHS
* - $x_6$
  - $\tfrac34$
  - $-\tfrac14$
  - $\tfrac34$
  - $\tfrac54$
  - $0$
  - $1$
  - $-\tfrac34$
  - $\tfrac14$
  - $2$
* - $x_5$
  - $-\tfrac12$
  - $\tfrac12$
  - $\tfrac12$
  - $-\tfrac12$
  - $1$
  - $0$
  - $\tfrac12$
  - $-\tfrac12$
  - $5$
* - $-z$
  - $32$
  - $22$
  - $30$
  - $36$
  - $0$
  - $0$
  - $3$
  - $8$
  - $-179$
```

```{list-table} Optimal tableau of the dual diet problem (maximization; $t_j$ is the slack of the price limit of food $j$)
:name: tab-03-05-diet-dual-tableau
:header-rows: 1
:stub-columns: 1

* - Basic
  - $\pi_1$
  - $\pi_2$
  - $t_1$
  - $t_2$
  - $t_3$
  - $t_4$
  - $t_5$
  - $t_6$
  - RHS
* - $\pi_1$
  - $1$
  - $0$
  - $0$
  - $0$
  - $0$
  - $0$
  - $-\tfrac12$
  - $\tfrac34$
  - $3$
* - $\pi_2$
  - $0$
  - $1$
  - $0$
  - $0$
  - $0$
  - $0$
  - $\tfrac12$
  - $-\tfrac14$
  - $8$
* - $t_1$
  - $0$
  - $0$
  - $1$
  - $0$
  - $0$
  - $0$
  - $\tfrac12$
  - $-\tfrac34$
  - $32$
* - $t_2$
  - $0$
  - $0$
  - $0$
  - $1$
  - $0$
  - $0$
  - $-\tfrac12$
  - $\tfrac14$
  - $22$
* - $t_3$
  - $0$
  - $0$
  - $0$
  - $0$
  - $1$
  - $0$
  - $-\tfrac12$
  - $-\tfrac34$
  - $30$
* - $t_4$
  - $0$
  - $0$
  - $0$
  - $0$
  - $0$
  - $1$
  - $\tfrac12$
  - $-\tfrac54$
  - $36$
* - $-w$
  - $0$
  - $0$
  - $0$
  - $0$
  - $0$
  - $0$
  - $-5$
  - $-2$
  - $-179$
```

The two tableaux contain the same information, arranged crosswise.

- The objective row of the primal tableau lists the dual slacks: the reduced
  costs $(32,22,30,36)$ of foods 1–4 are the values of $t_1,\dots,t_4$ in the
  dual tableau, and foods 5 and 6 have reduced cost 0 because their price
  limits are tight.
- The surplus columns $x_7$, $x_8$ of the primal objective row hold
  $\bar c=0-\vpi^\T(-\ve_i)=\pi_i$, that is, $\vpi^*=(3,8)$
  ({numref}`tab-03-05-reading-duals`).
- The slack columns $t_5$, $t_6$ of the dual objective row hold
  $(-5,-2)=-(x_5^*,x_6^*)$: the dual of the dual is the primal
  ({prf:ref}`thm-03-05-dual-of-dual`), and the same reading rule returns the
  optimal diet. The zero entries under $t_1,\dots,t_4$ and under $\pi_1$,
  $\pi_2$ correspond to the foods not bought and to the zero surpluses of the
  two vitamin constraints.
- Both corners hold $-179$.

*(An earlier version of this tableau shows 6 in the $x_4$ column of the
objective row; the reduced cost of food 4 is
$50-(2\cdot3+1\cdot8)=36$, the value of $t_4$ in the dual tableau.)*

(sec-03-05-shadow-prices)=
## Shadow prices and economic interpretation

### Dual values as rates of change

A **shadow price** is the value of one more unit of a scarce resource (or one
unit less of a requirement) to a decision maker who already uses it optimally:
the improvement in the optimal objective per unit of extra right-hand side,
measured at the margin. For linear programs the shadow prices are exactly the
optimal dual values, as long as the optimal basis does not change.

:::{prf:theorem} Optimal dual values are shadow prices
:label: thm-03-05-shadow-price

Let $\mB$ be an optimal basis of $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$
with multipliers $\vpi^{*\T}=\vc_B^\T\mB^{-1}$. For every change
$\Delta\vb$ with $\mB^{-1}(\vb+\Delta\vb)\ge\vzero$, the basis $\mB$ stays
optimal and

$$
z^*(\vb+\Delta\vb)=z^*(\vb)+\vpi^{*\T}\Delta\vb .
$$ (eq-03-05-shadow-price)

In particular, if $\mB$ is nondegenerate ($\mB^{-1}\vb>\vzero$), this holds for
all sufficiently small $\Delta\vb$, and
$\partial z^*/\partial b_i=\pi_i^*$ for every $i$. The same statements hold
for a maximization problem, with $z^*$ its maximal value.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

The reduced costs $\vc_N^\T-\vc_B^\T\mB^{-1}\mN$ do not depend on $\vb$, so
$\mB$ remains dual feasible after the change. If also
$\mB^{-1}(\vb+\Delta\vb)\ge\vzero$, $\mB$ is primal feasible for the new
right-hand side, hence optimal by {prf:ref}`prop-03-05-basis-optimal`, and
$z^*(\vb+\Delta\vb)=\vc_B^\T\mB^{-1}(\vb+\Delta\vb)=z^*(\vb)+\vpi^{*\T}\Delta\vb$.
If $\mB^{-1}\vb>\vzero$, the condition $\mB^{-1}(\vb+\Delta\vb)\ge\vzero$ holds
for all $\Delta\vb$ in a neighborhood of $\vzero$, so $z^*$ is linear there
with gradient $\vpi^*$.
:::

Changing one right-hand side, $\Delta\vb=\delta\ve_i$, keeps the basis optimal
as long as $\mB^{-1}\vb+\delta\,\mB^{-1}\ve_i\ge\vzero$, a set of linear
inequalities in $\delta$ that defines an interval, the **range of
validity** of the shadow price $\pi_i^*$. Outside it the basis changes and so
does the price. As a function of $b_i$ alone, $z^*$ is piecewise linear,
convex for a minimization problem and concave for a maximization problem; for
a resource constraint of a maximization problem this means that each extra
unit of the resource is worth no more than the one before. At a
degenerate optimum the optimal dual solution need not be unique, and the
one-sided rates of change of $z^*$ can differ; each optimal $\vpi^*$ then
gives a valid bound $z^*(\vb+\Delta\vb)\ge z^*(\vb)+\vpi^{*\T}\Delta\vb$ for a
minimization problem (by weak duality, since $\vpi^*$ stays dual feasible). The
systematic treatment of ranges, including cost coefficients and new
constraints, is in {doc}`03-07-sensitivity-analysis`.

For the product-mix problem ({prf:ref}`ex-03-05-product-mix-dual`) the optimal
basis is $(x_1,x_2,x_5)$, where $x_5$ is the slack of the demand constraint,
with

$$
\begin{gathered}
\mB^{-1}=\begin{pmatrix}\tfrac27&-\tfrac1{14}&0\\[2pt]-\tfrac17&\tfrac3{35}&0\\[2pt]-\tfrac27&\tfrac1{14}&1\end{pmatrix},
\\
\mB^{-1}\vb=\begin{pmatrix}\tfrac{45}{7}\\[2pt]\tfrac{30}{7}\\[2pt]\tfrac{11}{7}\end{pmatrix}.
\end{gathered}
$$

Adding $\delta$ hours of machine time changes the basic solution to
$\bigl(\tfrac{45}{7}+\tfrac27\delta,\ \tfrac{30}{7}-\tfrac17\delta,\
\tfrac{11}{7}-\tfrac27\delta\bigr)$, which stays nonnegative for
$-\tfrac{45}{2}\le\delta\le\tfrac{11}{2}$. The shadow price
$\pi_1^*=550/7\approx\$78.57$ per hour therefore holds for production time
between 37.5 and 65.5 hours per week: one extra hour raises the weekly
contribution from \$5142.86 to \$5221.43. Above 65.5 hours the demand limit
of 8 hundred cases of juice glasses becomes binding, the basis changes, and
extra machine time is worth nothing (the optimum stays at
$\vx=(8,3.5)$, $z=\$5575$). The same computation gives the range 128 to
240 m³ for the storage price $\pi_2^*=20/7\approx\$2.86$ per m³, and
the demand limit, with $\pi_3^*=0$, can fall to $45/7\approx6.43$ hundred
cases before it matters.

### The value of the resources

Strong duality reads, for the product-mix problem,

$$
\begin{aligned}
\underbrace{500x_1^*+450x_2^*}_{\text{contribution, \$/week}}
&=\underbrace{60\,\pi_1^*}_{\text{time}}+\underbrace{150\,\pi_2^*}_{\text{storage}}+\underbrace{8\,\pi_3^*}_{\text{demand}}\\
&=\frac{33000}{7}+\frac{3000}{7}+0=\frac{36000}{7}.
\end{aligned}
$$

Each term on the right is a quantity of a resource times its price per unit:
60 hours at \$78.57 per hour, 150 m³ at \$2.86 per m³ and 8 hundred cases
of demand at \$0. The optimal contribution is fully accounted for by the value
of the resources, and a resource that is not used up (here, the demand
limit) has price zero, which is complementary slackness. For any feasible
production plan and any feasible prices, weak duality says that the
contribution does not exceed the value of the resources,
$\vc^\T\vx\le\vb^\T\vpi$, with equality at the optimum.

### Imputed costs and the dual constraints

The dual constraints have an equally direct reading. For the primal
$\max\{\vc^\T\vx:\mA\vx\le\vb,\ \vx\ge\vzero\}$, the quantity

$$
\sum_{i=1}^{m}a_{ij}\,\pi_i=\vpi^\T\mA_{\cdot j}
$$ (eq-03-05-imputed-cost)

is the **imputed cost** of one unit of activity $j$: the resources it uses,
$a_{ij}$ units of resource $i$, valued at the prices $\pi_i$. Dual
feasibility $\vpi^\T\mA_{\cdot j}\ge c_j$ says that no activity earns more than
the value of the resources it consumes, and the reduced cost
$\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ is the unit revenue minus the imputed
cost. At the multipliers $\vpi^\T=\vc_B^\T\mB^{-1}$ of a basis that is not yet
optimal, some activity has $c_j>\vpi^\T\mA_{\cdot j}$: its unit revenue exceeds
the imputed cost of the resources it needs, so it pays to increase it, and the
simplex method brings it into the basis. At the optimum every activity in use
just breaks even and no activity is profitable at the margin. In the product
mix, a hundred cases of juice glasses use 6 hours and 10 m³, with imputed cost
$6\cdot\frac{550}{7}+10\cdot\frac{20}{7}+1\cdot0=500$ dollars, and a hundred
cases of cocktail glasses use 5 hours and 20 m³, with imputed cost
$5\cdot\frac{550}{7}+20\cdot\frac{20}{7}=450$ dollars: both equal the
contributions, as they must for products made in positive amounts.

:::{prf:remark} Shadow prices outside linear programming
:label: rem-03-05-shadow-price-meaning

The idea of a shadow price is used well beyond linear programming. In
cost–benefit analysis of public projects, many effects have no market price:
a new freeway that would save commuters 500 000 hours a year and five lives a
year, and reduce air pollution, at a present cost of \$250 million, can only be compared
with its cost after the hours and the lives are given monetary values, and
these values are shadow prices of time and of safety. They are estimated from
what people are willing to pay, and so they depend on assumptions. In
operations, a production line that already runs at its 40-hour weekly limit
has a shadow price for its time: the most the manager should pay for one
extra hour of operation is the increase in optimal profit that the hour makes
possible, which a linear model reports as the dual value of the time
constraint. In nonlinear programming, the Lagrange multipliers of the KKT
conditions play the same role ({prf:ref}`thm-06-03-kkt`); an example is
{prf:ref}`ex-08-03-nlp-multipliers`.
:::

(sec-03-05-farkas)=
## Farkas' lemma

Duality also answers a question that does not mention an objective at all:
when does a system of linear equations have a nonnegative solution? A
**theorem of the alternative** states that exactly one of two systems has a
solution; a solution of the second system then *certifies* that the first has
none. The best known one is due to Farkas (1902).

:::{prf:theorem} Farkas' lemma
:label: thm-03-05-farkas

Let $\mA\in\R^{m\times n}$ and $\vb\in\R^m$. Exactly one of the following
systems has a solution:

1. $\mA\vx=\vb$, $\vx\ge\vzero$;
2. $\mA^\T\vpi\le\vzero$, $\vb^\T\vpi>0$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Consider the linear program $\min\{\vzero^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$.
By {eq}`eq-03-05-standard-pair` its dual is
$\max\{\vb^\T\vpi:\mA^\T\vpi\le\vzero\}$, which is feasible ($\vpi=\vzero$).

Suppose system 1 has a solution. Then the primal is feasible with optimal value
0, and by weak duality every dual feasible $\vpi$ has $\vb^\T\vpi\le0$, so
system 2 has no solution.

Suppose system 1 has no solution. Then the primal is infeasible while the dual
is feasible, so by {prf:ref}`cor-03-05-feasible-pair` the dual is unbounded: there
is a $\vpi$ with $\mA^\T\vpi\le\vzero$ and $\vb^\T\vpi>0$ (indeed with
$\vb^\T\vpi$ as large as we like), and system 2 has a solution.
:::

Geometrically, system 1 says that $\vb$ lies in the cone generated by the
columns of $\mA$, $\cone\{\mA_{\cdot1},\dots,\mA_{\cdot n}\}$. System 2 says
that the hyperplane $\{\vy:\vpi^\T\vy=0\}$ separates $\vb$ from this cone: all
columns lie on one side ($\vpi^\T\mA_{\cdot j}\le0$) and $\vb$ strictly on the
other ($\vpi^\T\vb>0$). Farkas' lemma is thus a separating-hyperplane theorem
for a point and a polyhedral cone ({prf:ref}`thm-02-01-separation-point`;
a finitely generated cone is closed). The proof given
here leans on strong duality, which was proved with the simplex method; a proof
from the separation theorem is also possible, and then strong duality can be
derived from Farkas' lemma instead. Other forms follow by the usual
transformations, for example: exactly one of $\mA\vx\le\vb$, $\vx\ge\vzero$
and $\mA^\T\vpi\ge\vzero$, $\vpi\ge\vzero$, $\vb^\T\vpi<0$ has a solution
(apply the lemma to $(\mA\ \ \mI)$). The lemma is the step from linear
programming to the optimality conditions of nonlinear programming: in
{doc}`06-03-kkt-conditions` it turns "no improving feasible direction" into the
existence of KKT multipliers, exactly as {eq}`eq-03-05-cone-condition` did at a
vertex of a linear program.

## Summary and Common Mistakes

### Summary

- Every linear program has a dual; for the canonical pair
  {eq}`eq-03-05-canonical-pair`, $\min\{\vc^\T\vx:\mA\vx\ge\vb,\vx\ge\vzero\}$ and
  $\max\{\vb^\T\vpi:\mA^\T\vpi\le\vc,\vpi\ge\vzero\}$, the dual is the
  pricing problem of {prf:ref}`ex-03-05-diet`. The standard-form pair has free
  dual variables {eq}`eq-03-05-standard-pair`.
- The dual of any linear program follows from {numref}`tab-03-05-dual-rules`;
  the rules are the conditions under which the Lagrangian
  {eq}`eq-03-05-lagrangian` gives a lower bound. Duals of equivalent programs are
  equivalent ({prf:ref}`prop-03-05-equivalent-duals`) and the dual of the dual is
  the primal ({prf:ref}`thm-03-05-dual-of-dual`).
- Weak duality ({prf:ref}`thm-03-05-weak-duality`): feasible primal and dual
  values bound each other; equal values certify optimality
  ({prf:ref}`cor-03-05-optimality-certificate`).
- Strong duality ({prf:ref}`thm-03-05-strong-duality`): if one problem has an
  optimum, both do and $z^*=w^*$; the optimal dual is
  $\vpi^{*\T}=\vc_B^\T\mB^{-1}$. The possible outcomes of a pair are listed
  in {numref}`tab-03-05-status`.
- Complementary slackness ({prf:ref}`thm-03-05-complementary-slackness`):
  feasible $\vx$, $\vpi$ are optimal if and only if every variable–slack product in
  {eq}`eq-03-05-cs` is zero; at an optimal vertex $\vc$ lies in the cone of the
  active normals {eq}`eq-03-05-cone-condition` ({numref}`fig-03-05-cs-cone`).
- A basis is dual feasible if and only if its reduced costs have the optimal sign
  {eq}`eq-03-05-dual-feasible`; primal and dual feasible means optimal
  ({prf:ref}`prop-03-05-basis-optimal`). Optimal dual values can be read from
  the final tableau ({numref}`tab-03-05-reading-duals`).
- Optimal dual values are shadow prices, $\partial z^*/\partial b_i=\pi_i^*$,
  within a range of validity ({prf:ref}`thm-03-05-shadow-price`); dual
  constraints compare unit revenue with imputed cost {eq}`eq-03-05-imputed-cost`.
- Farkas' lemma ({prf:ref}`thm-03-05-farkas`) follows from duality: either
  $\mA\vx=\vb$ has a solution $\vx\ge\vzero$, or a vector $\vpi$ separates
  $\vb$ from the cone of the columns of $\mA$.

### Common mistakes

- **Wrong sign of a dual variable.** A constraint of the opposite type (a
  $\le$ row in a minimization problem) has $\pi_i\le0$, not $\pi_i\ge0$; an
  equation has a free $\pi_i$. When unsure, multiply the row by $-1$ and use
  the canonical pair.
- **Dropping the sign condition in weak duality.** The step
  $\vpi^\T\mA\vx\ge\vpi^\T\vb$ needs $\vpi\ge\vzero$ in the canonical pair;
  with free prices it fails unless the constraints are equations.
- **Using complementary slackness without checking feasibility.** The
  equations obtained from complementary slackness give a candidate; it is
  optimal only if it also satisfies all constraints of its problem
  ({prf:ref}`ex-03-05-diet-cs`).
- **Reading "exactly one is zero" into complementary slackness.** A variable
  and its matching slack may both be zero.
- **Sign errors when reading duals from a tableau.** Use
  $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ with the actual column and cost of the
  slack, surplus or artificial variable ({numref}`tab-03-05-reading-duals`)
  rather than memorized rules for one tableau convention.
- **Using a shadow price outside its range.** $\pi_i^*$ predicts the change
  of $z^*$ only while the optimal basis stays primal feasible; a large change
  of $b_i$ usually earns less per unit.
- **Concluding unboundedness from infeasibility.** An infeasible primal does
  not imply an unbounded dual; both may be infeasible
  ({numref}`tab-03-05-status`).

## References

- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Chapter 4: Sections 4.1–4.3 on the dual
  problem and the duality theorems, Section 4.4 on optimal dual variables as
  marginal costs, Section 4.6 on Farkas' lemma.)
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Chapter 6: the simplex algorithm and the duality theorem;
  Chapter 12: the price concept in linear programming.)
- Farkas, J. (1902). Theorie der einfachen Ungleichungen. *Journal für die
  reine und angewandte Mathematik, 124*, 1–27.
- Gale, D., Kuhn, H. W., & Tucker, A. W. (1951). Linear programming and the
  theory of games. In T. C. Koopmans (Ed.), *Activity analysis of production
  and allocation* (pp. 317–329). Wiley.
- Goldman, A. J., & Tucker, A. W. (1956). Theory of linear programming. In
  H. W. Kuhn & A. W. Tucker (Eds.), *Linear inequalities and related systems*
  (pp. 53–97). Princeton University Press.
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 6: Sections 6.1–6.4 on the
  essence of duality theory, its economic interpretation, primal–dual
  relationships and the dual of other primal forms.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 4: Sections 4.1–4.3 on the definition of the dual,
  primal–dual relationships and the economic interpretation of duality;
  Section 7.4: duality in matrix form.)
