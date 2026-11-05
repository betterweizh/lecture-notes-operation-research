# Linearization Techniques

## Overview

Many models are *almost* linear: a product of two decisions appears in the
objective, the cost depends on an absolute deviation, the worst case of
several linear expressions is minimized, or the target is a ratio such as
benefit per dollar. This chapter collects the standard techniques that
rewrite such terms exactly, or bound them tightly, with linear constraints,
continuous auxiliary variables and, where unavoidable, binary variables. A
product of two binaries is the logical AND of {doc}`07-01-logical-constraints`;
a product of a binary and a bounded continuous variable is linearized
exactly by four inequalities; a product of two continuous variables is
enclosed by the **McCormick envelope**, derived here from four products of
bound distances. Absolute values and maxima and minima of linear functions
are linear when they appear on the "convex side" of the model and need
binaries otherwise; the chapter states this rule once and applies it to
several cases. Finally, the **Charnes–Cooper transformation** turns a ratio
of two linear functions into a linear program (LP). Each technique is shown on an
infrastructure example whose solution is verified numerically. Besides
{doc}`07-01-logical-constraints`, the chapter uses LP modeling
({doc}`03-01-lp-model-building`) and epigraphs of convex functions
({prf:ref}`def-02-02-epigraph`).

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- linearize products of binary variables, and of a binary and a bounded
  continuous variable, and know when one side of the formulation can be
  dropped;
- derive the McCormick inequalities for $xy$ on a box, state when they are
  exact, and tighten them by partitioning;
- model absolute values and maxima/minima of linear functions, deciding
  whether an LP suffices or binaries are needed;
- transform a linear-fractional program into an LP and recover its
  solution; and
- check every linearization by substituting the integer values of its
  binaries.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\omega$, $\omega_{ij}$ | auxiliary variable that replaces a product ($\omega=xy$, $\omega_{ij}=\delta_i\delta_j$) |
| $\underline x$, $\overline x$ | lower and upper bound of a variable $x$ |
| $\eta$, $\eta_i$ | auxiliary variable that replaces an absolute value or a minimum/maximum |
| $\vx^+$, $\vx^-$ | positive and negative parts, $\vx=\vx^+-\vx^-$, $\vx^\pm\ge\vzero$ |
| $\vr_k^\T\vx+r_{k0}$ | the $k$-th affine function inside a maximum or minimum, $k=1,\dots,K$ |
| $\vp^\T\vx+p_0$, $\vq^\T\vx+q_0$ | numerator and denominator of a linear-fractional objective |
| $\sigma$, $\vy$ | Charnes–Cooper variables, $\sigma=1/(\vq^\T\vx+q_0)$, $\vy=\sigma\vx$ |

A **linearization** is *exact* if, for every value of the original
variables, the new constraints admit exactly the auxiliary values that the
nonlinear term defines; it is a *relaxation* if they admit more. The
binary products and the absolute values of this chapter are exact; the
McCormick envelope of a continuous product is a relaxation.

(sec-07-02-binary-products)=
## Products of binary variables

Let $\delta_1,\delta_2\in\{0,1\}$. Their product is $1$ exactly when both are
$1$, so $\omega=\delta_1\delta_2$ is the logical AND, and it is linearized by
{eq}`eq-07-01-and`:

$$
\begin{aligned}
\omega&\le\delta_1,\\
\omega&\le\delta_2,\\
\omega&\ge\delta_1+\delta_2-1,\\
\omega&\ge0 .
\end{aligned}
$$ (eq-07-02-bin-bin)

Substituting the four pairs: $(0,0)$, $(0,1)$ and $(1,0)$ give
$0\le\omega\le0$, and $(1,1)$ gives $1\le\omega\le1$. Note that $\omega$
need not be declared binary: the constraints force the right integer value
whenever $\delta_1,\delta_2$ are integer. A product of $k$ binaries is
linearized the same way, with $\omega\le\delta_i$ for each factor and
$\omega\ge\sum_i\delta_i-(k-1)$. Since $\delta^2=\delta$ for a binary, a
**quadratic 0–1 program** $\max\ \sum_i c_i\delta_i+\sum_{i<j}c_{ij}\delta_i\delta_j$
becomes a mixed-integer linear program (MILP) with one $\omega_{ij}$ per nonzero $c_{ij}$ (Glover & Woolsey, 1974;
Glover, 1975).

**One-sided formulations.** When $\omega$ appears only in the objective, half
of {eq}`eq-07-02-bin-bin` is often redundant. In a maximization with
$c_{ij}>0$ the optimizer pushes $\omega_{ij}$ up, so only the upper bounds
$\omega_{ij}\le\delta_i$, $\omega_{ij}\le\delta_j$ are needed; with
$c_{ij}<0$ it pushes $\omega_{ij}$ down, so only
$\omega_{ij}\ge\delta_i+\delta_j-1$ and $\omega_{ij}\ge0$ are needed. Dropping
the wrong side is a common error, as the next example shows.

:::{prf:example} Project portfolio with synergies
:label: ex-07-02-synergy

Four rail projects compete for a budget of \$105 million. Project $i$ costs
$a_i=60,20,25,15$ and has a stand-alone benefit $b_i=80,18,24,16$ (all in
\$ million). Some pairs interact: projects 1 and 2 together add \$20 million
(shared stations), 1 and 3 add \$15 million, 3 and 4 add \$6 million, while 2
and 3 serve the same passengers and together *lose* \$20 million. Choose
the portfolio of greatest total benefit, and check what goes wrong if only
the upper bounds $\omega_{ij}\le\delta_i$, $\omega_{ij}\le\delta_j$ of the
product linearization are kept for every pair.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Model.* With $\delta_i=1$ if project $i$ is chosen,
$\omega_{ij}=\delta_i\delta_j$ and the set of interacting pairs
$P=\{(1,2),(1,3),(2,3),(3,4)\}$,

$$
\begin{aligned}
\max\ & 80\delta_1+18\delta_2+24\delta_3+16\delta_4\\
& +20\omega_{12}+15\omega_{13}-20\omega_{23}+6\omega_{34}\\
\st\ & 60\delta_1+20\delta_2+25\delta_3+15\delta_4\le105,\\
& \omega_{ij}\le\delta_i, \qquad (i,j)\in P,\\
& \omega_{ij}\le\delta_j, \qquad (i,j)\in P,\\
& \omega_{ij}\ge\delta_i+\delta_j-1, \qquad (i,j)\in P,\\
& \omega_{ij}\ge0, \qquad (i,j)\in P,\\
& \vdelta\in\{0,1\}^4 .
\end{aligned}
$$

*Solution.* Enumerating the $16$ portfolios (or solving the MILP) gives the
optimum $\{1,3,4\}$: cost $60+25+15=100$, benefit
$80+24+16+15+6=\$141$ million. The runner-up $\{1,2,3\}$ costs exactly $105$
and is worth $80+18+24+20+15-20=137$.

*The wrong one-sided model.* If only the upper bounds
$\omega_{ij}\le\delta_i,\delta_j$ are kept for all pairs, nothing stops
$\omega_{23}=0$ when $\delta_2=\delta_3=1$: the penalty disappears, and the
solver reports $\{1,2,3\}$ with a "benefit" of $157$. Keeping the upper
bounds for the three positive synergies and $\omega_{23}\ge\delta_2+\delta_3-1$
for the negative one gives the correct optimum, $141$, with fewer rows; its
LP relaxation bound is $145.5$.
:::

(sec-07-02-binary-continuous)=
## Product of a binary and a continuous variable

Let $\delta\in\{0,1\}$ and $\underline x\le x\le\overline x$ with
$\underline x\ge0$. The product $\omega=\delta x$ equals $0$ when $\delta=0$
and $x$ when $\delta=1$. Both cases are captured by

$$
\begin{aligned}
\underline x\,\delta&\le\omega\le\overline x\,\delta,\\
x-\overline x(1-\delta)&\le\omega\le x-\underline x(1-\delta) .
\end{aligned}
$$ (eq-07-02-bin-cont)

For $\delta=0$ the first line gives $\omega=0$ and the second line reads
$x-\overline x\le0\le x-\underline x$, which holds. For $\delta=1$ the second
line gives $\omega=x$ and the first line reads $\underline x\le x\le\overline x$.
The pairs of inequalities are big-$M$ constraints with the tightest
constants (the bounds of $x$); they coincide with the McCormick inequalities
of {eq}`eq-07-02-mccormick` for the box
$[0,1]\times[\underline x,\overline x]$, which explains why they are exact
(Glover, 1975).

In the frequent case $\underline x=0$ the lower bound in the first line is
$\omega\ge0$; it must not be forgotten. *(Without $\omega\ge0$, the
combination $\delta=0$, $\omega=x-\overline x<0$ would be feasible, and a
model that maximizes $-\omega$ would exploit it.)* With $\underline x>0$ the
right-hand inequality of the second line may be replaced by the weaker
$\omega\le x$, which is still exact for binary $\delta$ but gives a weaker
relaxation.

As with binary products, a one-sided version suffices when the sign of
$\omega$'s objective coefficient is known: an $\omega$ that the optimizer wants
large needs only the upper bounds $\omega\le\overline x\delta$ and
$\omega\le x-\underline x(1-\delta)$.

:::{prf:example} Supply contract with a volume discount
:label: ex-07-02-contract

A road agency needs $800$ tonnes of asphalt binder per year. Supplier 1
sells at \$80 per tonne up to $500$ tonnes. Supplier 2 sells up to $600$
tonnes at \$90 per tonne, or at \$70 per tonne if the agency signs an annual
contract costing \$6000. Minimize the annual cost.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Model.* Let $x_1,x_2$ be the tonnes bought and $\delta=1$ if the contract
is signed. The cost $80x_1+90x_2-20\,\delta x_2+6000\,\delta$ contains the
product $\omega=\delta x_2$ with $x_2\in[0,600]$:

$$
\begin{aligned}
\min\ & 80x_1+90x_2-20\omega+6000\delta\\
\st\ & x_1+x_2=800,\\
& 0\le x_1\le500,\\
& 0\le x_2\le600,\\
& 0\le\omega\le600\,\delta,\\
& x_2-600(1-\delta)\le\omega\le x_2,\\
& \delta\in\{0,1\} .
\end{aligned}
$$

*Solution.* Without the contract, buy the cheaper $500$ tonnes from
supplier 1 and $300$ from supplier 2: $40\,000+27\,000=\$67\,000$. With the
contract supplier 2 is cheaper, so buy $600$ tonnes from it and $200$ from
supplier 1: $16\,000+42\,000+6000=\$64\,000$. The MILP returns
$\delta^*=1$, $x_1^*=200$, $x_2^*=\omega^*=600$, cost \$64 000. The contract
pays off once at least $6000/20=300$ tonnes are bought from supplier 2.
Because $\omega$ has a negative cost coefficient (it is wanted large), the
rows $\omega\le600\delta$ and $\omega\le x_2$ alone already give this
answer.
:::

**General integers.** A bounded integer variable $x\in\{0,1,\dots,U\}$ can
be written in binary form,

$$
\begin{aligned}
x&=\sum_{k=0}^{K-1}2^k\delta_k,\qquad \delta_k\in\{0,1\},\\
K&=\lceil\log_2(U+1)\rceil,
\end{aligned}
$$ (eq-07-02-binary-expansion)

with the extra row $x\le U$ when $U+1$ is not a power of two. A product
$xy$ with a bounded continuous $y\in[0,\overline y]$ then becomes
$\sum_k2^k\omega_k$ with $\omega_k=\delta_ky$, each linearized by
{eq}`eq-07-02-bin-cont`: $K$ new variables and $4K$ rows instead of one
nonlinear term. For example, a number of buses $x\in\{0,\dots,12\}$ times a
continuous load per bus needs $K=4$ binaries ($2^4=16>12$). The same
expansion makes a product of two integers a sum of binary×integer products,
and so linear (Glover, 1975).

(sec-07-02-mccormick)=
## McCormick envelopes

Now let both factors be continuous, $\omega=xy$ with
$x\in[\underline x,\overline x]$ and $y\in[\underline y,\overline y]$. The
set $\{(x,y,\omega):\omega=xy\}$ is a saddle surface and is not convex, so no
LP describes it. The best one can do with linear inequalities is to enclose
it. Every point of the box satisfies four obvious inequalities,

$$
\begin{aligned}
(x-\underline x)(y-\underline y)&\ge0,\\
(\overline x-x)(\overline y-y)&\ge0,\\
(\overline x-x)(y-\underline y)&\ge0,\\
(x-\underline x)(\overline y-y)&\ge0 .
\end{aligned}
$$

Expanding each product and replacing $xy$ by $\omega$ gives the
**McCormick inequalities** (McCormick, 1976):

$$
\begin{aligned}
\omega&\ge\underline y\,x+\underline x\,y-\underline x\,\underline y,\\
\omega&\ge\overline y\,x+\overline x\,y-\overline x\,\overline y,\\
\omega&\le\underline y\,x+\overline x\,y-\overline x\,\underline y,\\
\omega&\le\overline y\,x+\underline x\,y-\underline x\,\overline y .
\end{aligned}
$$ (eq-07-02-mccormick)

The first two are **under-estimators** of $xy$, the last two
**over-estimators**. Their maximum and minimum are the convex and the
concave envelope of $xy$ on the box (Al-Khayyal & Falk, 1983), that is, the
largest convex function below $xy$ and the smallest concave function above it; no other linear
inequalities in $(x,y,\omega)$ valid on the box can cut deeper.

**When it is exact.** If $x$ is at a bound, say $x=\underline x$, the first
and the fourth inequalities give $\omega\ge\underline x\,y$ and
$\omega\le\underline x\,y$, so $\omega=xy$. The same holds at $x=\overline x$,
$y=\underline y$ and $y=\overline y$. Thus the envelope is exact on the
boundary of the box and, in particular, whenever one factor is binary: with
$x=\delta\in[0,1]$, {eq}`eq-07-02-mccormick` is exactly {eq}`eq-07-02-bin-cont`
and {eq}`eq-07-02-bin-bin`.

**How loose it is.** Inside the box the envelope is a relaxation. On the
diagonal through $(\underline x,\underline y)$ and $(\overline x,\overline y)$
the product is convex and the over-estimators lie above it; on the other
diagonal it is concave and the under-estimators lie below it
({numref}`fig-07-02-mccormick`). The largest gap on either side is
$(\overline x-\underline x)(\overline y-\underline y)/4$, attained at the
center of the box. Halving the width of one interval halves the gap: this is
the idea of **partitioning**. Split $[\underline x,\overline x]$ into $K$
pieces, let a binary choose the piece that contains $x$ (as in
{doc}`07-01-logical-constraints`), and write the envelope with the bounds of the
chosen piece; the gap then shrinks like $1/K$.

```{figure} figures/07-02-mccormick.svg
:name: fig-07-02-mccormick
:width: 640px
:class: dark-light
:alt: The product xy along the two diagonals of a box, with the McCormick under- and over-estimators

The product $\omega=xy$ on the box $[1,4]\times[1,3]$ of
{prf:ref}`ex-07-02-mccormick`, along the main diagonal (left, where $xy$ is
convex) and the other diagonal (right, where it is concave). The shaded band
lies between the larger of the two McCormick under-estimators and the
smaller of the two over-estimators; it touches $xy$ at the corners of the
box and is widest, by $(3\cdot2)/4=1.5$, at the center.
```

:::{prf:example} Bounding a bilinear objective
:label: ex-07-02-mccormick

Maximize $\omega=xy$ subject to $2x+3y\le12$, $1\le x\le4$, $1\le y\le3$.
(Think of $x$ as a price and $y$ as a quantity that the budget line links.)
Find the exact optimum and the upper bound given by the McCormick
relaxation, and show how the bound improves when $[1,4]$ is split into $K$
equal pieces.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Exact optimum.* At the optimum the budget constraint is active,
$y=(12-2x)/3$, and $xy=x(12-2x)/3$ is maximized at $x=3$: $\vx^*=(3,2)$
with $xy=6$.

*McCormick relaxation.* Replacing $xy$ by $\omega$ with
{eq}`eq-07-02-mccormick` for $\underline x=1$, $\overline x=4$, $\underline y=1$,
$\overline y=3$ gives the LP

$$
\begin{aligned}
\max\ & \omega\\
\st\ & \omega\ge x+y-1,\\
& \omega\ge3x+4y-12,\\
& \omega\le x+4y-4,\\
& \omega\le3x+y-3,\\
& 2x+3y\le12,\\
& 1\le x\le4,\\
& 1\le y\le3 .
\end{aligned}
$$

Its optimum lies where the two over-estimators and the budget line meet:
$x+4y-4=3x+y-3$ gives $3y=2x+1$, and with $2x+3y=12$ this yields
$x=\tfrac{11}4$, $y=\tfrac{13}6$, $\omega=\tfrac{89}{12}\approx7.417$. The
relaxation overestimates the true maximum by $1.417$, and its point has
true product $\tfrac{11}4\cdot\tfrac{13}6=\tfrac{143}{24}\approx5.958$.

*Partitioning.* Splitting $[1,4]$ into $K$ equal pieces and taking the best
bound over the pieces (which is what the MILP with one binary per piece
computes) gives

| $K$ | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| bound on $\max xy$ | 7.417 | 6.722 | 6.367 | 6.185 | 6.093 |
| gap to 6 | 1.417 | 0.722 | 0.367 | 0.185 | 0.093 |

(exact values $\tfrac{89}{12}$, $\tfrac{121}{18}$, $\tfrac{191}{30}$,
$\tfrac{167}{27}$, $\tfrac{1243}{204}$). The gap halves each time $K$
doubles, as predicted.
:::

(sec-07-02-convex-side)=
## Absolute values, maxima and minima: the convex side

The next techniques rest on one observation. A convex function
$\chi(\vx)$ can be bounded from **above** by linear constraints when it is a
maximum of affine functions,

$$
\begin{aligned}
&\chi(\vx)=\max_{k=1,\dots,K}\bigl(\vr_k^\T\vx+r_{k0}\bigr)\le\eta\\
\iff\ &\vr_k^\T\vx+r_{k0}\le\eta\quad\text{for all }k,
\end{aligned}
$$ (eq-07-02-max-le)

because a maximum is at most $\eta$ exactly when every term is. Any convex
piecewise-linear function has this form, and $\abs{x}=\max\{x,-x\}$ is the
simplest case. By the same token a concave piecewise-linear function
(a minimum of affine functions) can be bounded from **below**:
$\min_k(\cdot)\ge\eta$ iff every term is $\ge\eta$. More generally,
$\min_af(a,\vx)\ge0$ is equivalent to $f(a,\vx)\ge0$ for all $a$, and
$\max_af(a,\vx)\le0$ to $f(a,\vx)\le0$ for all $a$.

The opposite statements, $\max_k(\cdot)\ge\eta$ or $\min_k(\cdot)\le\eta$,
say that *at least one* term satisfies the bound. They are disjunctions and
need binaries (at least $1$ of $K$, {eq}`eq-07-01-k-of-m`):

$$
\begin{aligned}
\vr_k^\T\vx+r_{k0}&\ge\eta-M_k(1-\delta_k),\quad k=1,\dots,K,\\
\sum_{k=1}^K\delta_k&\ge1,\\
\delta_k&\in\{0,1\},\quad k=1,\dots,K .
\end{aligned}
$$ (eq-07-02-max-ge)

In an objective the rule reads: **minimizing a convex** (or maximizing a
concave) piecewise-linear function is an LP, through the **epigraph**
variable $z$:

$$
\begin{gathered}
\begin{aligned}
\min_{\vx}\ & \max_k\bigl(\vr_k^\T\vx+r_{k0}\bigr)\\
\st\ & \mA\vx\ge\vb
\end{aligned}\\
\Updownarrow\\
\begin{aligned}
\min_{\vx,z}\ & z\\
\st\ & z\ge\vr_k^\T\vx+r_{k0},\qquad k=1,\dots,K,\\
& \mA\vx\ge\vb .
\end{aligned}
\end{gathered}
$$ (eq-07-02-minmax)

At an optimum $z$ settles on the largest term, since any larger $z$ could be
reduced. **Maximizing a convex** function (or minimizing a concave one) is
not an LP: in $\max_{\vx}\max_k(\cdot)$ the epigraph rows
$z\ge\vr_k^\T\vx+r_{k0}$ are useless (the optimizer would push $z$ to
infinity), and the hypograph condition $z\le\max_k(\cdot)$ is the
disjunction {eq}`eq-07-02-max-ge`. One can either use binaries,

$$
\begin{aligned}
\max_{\vx,z,\vdelta}\ & z\\
\st\ & z\le\vr_k^\T\vx+r_{k0}+M_k(1-\delta_k),\qquad k=1,\dots,K,\\
& \textstyle\sum_k\delta_k=1,\\
& \mA\vx\ge\vb,\\
& \delta_k\in\{0,1\},\qquad k=1,\dots,K,
\end{aligned}
$$

or simply solve the $K$ LPs $\max\{\vr_k^\T\vx+r_{k0}:\mA\vx\ge\vb\}$ and keep
the best; the latter is often faster when $K$ is small.

**Defining a minimum exactly.** Sometimes a variable must *equal* a minimum,
$\eta=\min\{x_1,x_2\}$, for example when the flow through a junction is the
smaller of two capacities that are themselves decisions. The rows
$\eta\le x_1$, $\eta\le x_2$ only give $\eta\le\min$; equality needs the
disjunction "$\eta\ge x_1$ or $\eta\ge x_2$":

$$
\begin{aligned}
\eta&\le x_1,\\
\eta&\le x_2,\\
\eta&\ge x_1-M(1-\delta),\\
\eta&\ge x_2-M\delta,\qquad \delta\in\{0,1\},
\end{aligned}
$$ (eq-07-02-eta-min)

with $M\ge\max\abs{x_1-x_2}$ over the bounds. If $x_1\le x_2$, $\delta=1$ is
possible and forces $\eta=x_1$; otherwise $\delta=0$ and $\eta=x_2$. For
$\eta=\max\{x_1,x_2\}$ apply this to $-\eta=\min\{-x_1,-x_2\}$; for
$\eta=\min_ix_i$ over $n$ terms use one binary per term with
$\sum_i\delta_i=1$ (or $\ge1$).

(sec-07-02-abs)=
### Absolute values

By the rule above:

* $\abs{f(\vx)}\le c$ (with $f$ affine and $c\ge0$) is
  $-c\le f(\vx)\le c$, two linear rows;
* $\abs{f(\vx)}\ge c$ is the disjunction $f(\vx)\ge c$ or $f(\vx)\le-c$, with
  one binary: $f(\vx)\ge c-M\delta$ and $f(\vx)\le-c+M(1-\delta)$;
* $\abs{f(\vx)}=c$ for a constant $c\ge0$ means $f(\vx)=c$ or
  $f(\vx)=-c$, which one binary expresses exactly:
  $f(\vx)=c\,(2\delta-1)$. *(Writing $\abs{f(\vx)}=c$ as
  $f(\vx)=\max\{-c,c\}$ is incorrect: that maximum is just $c$.)* If $c$ is
  itself a bounded variable, the product $\delta c$ is linearized by
  {eq}`eq-07-02-bin-cont`.

In the objective, consider

$$
\begin{aligned}
\min_{\vx}\ & \sum_{i=1}^nc_i\abs{x_i}\qquad (c_i\ge0)\\
\st\ & \mA\vx\ge\vb .
\end{aligned}
$$ (eq-07-02-abs-obj)

There are two equivalent LP formulations. **Epigraph form**: introduce
$\eta_i\ge\abs{x_i}$ through $\eta_i\ge x_i$ and $\eta_i\ge-x_i$, and minimize
$\sum_ic_i\eta_i$. **Split form**: write $\vx=\vx^+-\vx^-$ with
$\vx^+,\vx^-\ge\vzero$, and minimize $\sum_ic_i(x_i^++x_i^-)$ subject to
$\mA\vx^+-\mA\vx^-\ge\vb$. In the split form nothing forbids both
$x_i^+>0$ and $x_i^->0$, in which case $x_i^++x_i^->\abs{x_i}$. But such a
solution is never optimal when $c_i>0$: lowering both parts by
$\min\{x_i^+,x_i^-\}$ keeps $\vx$ and the constraints unchanged and reduces
the cost. Hence at an optimum $x_i^+x_i^-=0$ and $x_i^++x_i^-=\abs{x_i}$
(Bertsimas & Tsitsiklis, 1997, §1.3).

Both arguments need $c_i\ge0$. If some $c_i<0$ the model *maximizes*
$\abs{x_i}$, which is the nonconvex direction: in the epigraph form
$\eta_i$ can grow without limit, in the split form $x_i^+$ and $x_i^-$ grow
together, and either LP becomes unbounded or wrong. A binary per such
variable restores the meaning: $x_i^+\le U_i\delta_i$,
$x_i^-\le U_i(1-\delta_i)$, with $U_i$ a bound on $\abs{x_i}$.

The three distance functions of the next example illustrate the rule
({numref}`fig-07-02-abs-max`, in its solution): a weighted sum of absolute
values and a maximum of absolute values are convex (minimizing them is an
LP), while a minimum of absolute values is not concave (maximizing it needs
binaries).

:::{prf:example} Locating facilities on a corridor
:label: ex-07-02-corridor

Five towns lie along a $50$ km highway corridor at positions
$a_i=0,12,20,34,50$ km with populations $w_i=8,5,12,4,6$ thousand. A single
facility is placed at a position $x\in[0,50]$. Find the best position for
(a) a maintenance depot that minimizes the total resident–distance
$\sum_iw_i\abs{x-a_i}$; (b) an emergency station that minimizes the largest
distance to a town; (c) a landfill that maximizes the smallest distance to
a town.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

```{figure} figures/07-02-abs-max.svg
:name: fig-07-02-abs-max
:width: 680px
:class: dark-light
:alt: Three piecewise-linear functions of a location on a corridor: a weighted sum of distances, the maximum distance and the minimum distance

Distance functions for the five towns of {prf:ref}`ex-07-02-corridor` (town
positions marked on the axis). Left: total weighted distance
$\sum_iw_i\abs{x-a_i}$, convex, minimum at $x=20$. Middle: largest distance
$\max_i\abs{x-a_i}$, convex, minimum at $x=25$. Right: smallest distance
$\min_i\abs{x-a_i}$, neither convex nor concave, maximum at $x=42$.
```

**(a) Depot.** With $\eta_i\ge\abs{x-a_i}$,

$$
\begin{aligned}
\min\ & \sum_{i=1}^5w_i\eta_i\\
\st\ & \eta_i\ge x-a_i, \qquad i=1,\dots,5,\\
& \eta_i\ge a_i-x, \qquad i=1,\dots,5,\\
& 0\le x\le50 .
\end{aligned}
$$

The LP gives $x^*=20$ km with $436$ thousand resident-km. Check: the slope
of $\sum_iw_i\abs{x-a_i}$ is (population to the left) minus (population to
the right); it is $8+5-(12+4+6)=-9$ just left of $20$ and
$8+5+12-(4+6)=+15$ just right of it, so $x=20$ is the minimum, a
**weighted median** (the cumulative population $8,13,25$ first passes half
of $35$ at the third town). The cost there is
$8\cdot20+5\cdot8+0+4\cdot14+6\cdot30=436$; at $x=12$ it would be $508$ and at
$x=34$, $646$.

**(b) Station.** Minimize $z$ subject to $z\ge x-a_i$ and $z\ge a_i-x$ for
all $i$ ({eq}`eq-07-02-minmax`). Only the end towns matter, and the LP gives
$x^*=25$ km, $z^*=25$ km, the midpoint of $0$ and $50$.

**(c) Landfill.** Maximizing $\min_i\abs{x-a_i}$ requires
$z\le\abs{x-a_i}$ for every $i$, and each of these is a disjunction
("the town is at least $z$ to the left or at least $z$ to the right"). With
$\delta_i=1$ meaning "the facility lies left of town $i$" and $M=100$ km
(twice the corridor length bounds both sides),

$$
\begin{aligned}
\max\ & z\\
\st\ & z\le x-a_i+M\delta_i, \qquad i=1,\dots,5,\\
& z\le a_i-x+M(1-\delta_i), \qquad i=1,\dots,5,\\
& 0\le x\le50,\\
& \delta_i\in\{0,1\}, \qquad i=1,\dots,5 .
\end{aligned}
$$

The MILP gives $x^*=42$ km, $z^*=8$ km, with $\delta_5=1$ and the other
binaries $0$: the landfill sits midway between the towns at $34$ and $50$
km, the widest gap. (The gaps are $12$, $8$, $14$ and $16$ km; the end of
the corridor at $0$ or $50$ is a town, so no end point is better.) The LP
relaxation of this model is useless: it reports $z=50$, since fractional
$\delta_i$ relax every row.
:::

(sec-07-02-fractional)=
## Linear-fractional objectives

A **linear-fractional program** optimizes a ratio of affine functions, such
as benefit per dollar or cost per passenger:

$$
\begin{aligned}
\max_{\vx}\ & \frac{\vp^\T\vx+p_0}{\vq^\T\vx+q_0}\\
\st\ & \mA\vx\le\vb,\\
& \vx\ge\vzero,
\end{aligned}
$$ (eq-07-02-fractional)

where the denominator is positive on the feasible set
$\Omega=\{\vx:\mA\vx\le\vb,\vx\ge\vzero\}$, assumed nonempty and bounded.
The objective is neither convex nor concave, but it is quasilinear, and
**Charnes and Cooper (1962)** showed that a change of variables turns the
problem into an LP. Set

$$
\begin{aligned}
\sigma&=\frac1{\vq^\T\vx+q_0}>0,\\
\vy&=\sigma\vx .
\end{aligned}
$$

Multiplying numerator, denominator and constraints by $\sigma$ gives

$$
\begin{aligned}
\max_{\vy,\sigma}\ & \vp^\T\vy+p_0\,\sigma\\
\st\ & \mA\vy-\vb\,\sigma\le\vzero,\\
& \vq^\T\vy+q_0\,\sigma=1,\\
& \vy\ge\vzero,\\
& \sigma\ge0 .
\end{aligned}
$$ (eq-07-02-charnes-cooper)

:::{prf:theorem} Charnes–Cooper transformation
:label: thm-07-02-charnes-cooper

If $\Omega$ is nonempty and bounded and $\vq^\T\vx+q_0>0$ on $\Omega$, then
{eq}`eq-07-02-fractional` and {eq}`eq-07-02-charnes-cooper` have the same optimal
value, and $\vx^*=\vy^*/\sigma^*$ solves {eq}`eq-07-02-fractional` for any
optimal solution $(\vy^*,\sigma^*)$ of {eq}`eq-07-02-charnes-cooper`.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

For $\vx\in\Omega$ the pair $\sigma=1/(\vq^\T\vx+q_0)$, $\vy=\sigma\vx$ is
feasible for {eq}`eq-07-02-charnes-cooper` with objective value equal to the
ratio at $\vx$. Conversely, let $(\vy,\sigma)$ be feasible. If $\sigma=0$,
then $\mA\vy\le\vzero$, $\vy\ge\vzero$ and $\vq^\T\vy=1$, so $\vy\ne\vzero$
and $\vx+\lambda\vy\in\Omega$ for every $\vx\in\Omega$ and $\lambda\ge0$,
contradicting boundedness. Hence $\sigma>0$, and $\vx=\vy/\sigma$ satisfies
$\mA\vx\le\vb$, $\vx\ge\vzero$ and $\vq^\T\vx+q_0=1/\sigma$, so its ratio
equals $\vp^\T\vy+p_0\sigma$. The two problems therefore attain the same set
of values, and optimal solutions correspond.
:::

For minimization the same transformation applies. If the denominator were
negative on $\Omega$, multiply numerator and denominator by $-1$ first.

:::{prf:example} Best benefit–cost ratio
:label: ex-07-02-fractional

An agency assigns maintenance crews to pavement resurfacing ($x_1$
crew-weeks) and drainage repair ($x_2$ crew-weeks). The benefit is
$5x_1+3x_2$ and the cost $2x_1+2x_2+4$ (\$ million, including a fixed
mobilization cost of \$4 million). At most $12$ crew-weeks are available,
resurfacing can absorb at most $8$, and drainage needs at least $2$.
Maximize the benefit–cost ratio.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Transformation.* With $\sigma=1/(2x_1+2x_2+4)$ and $\vy=\sigma\vx$, the LP
{eq}`eq-07-02-charnes-cooper` is

$$
\begin{aligned}
\max\ & 5y_1+3y_2\\
\st\ & y_1+y_2-12\sigma\le0,\\
& y_1-8\sigma\le0,\\
& -y_2+2\sigma\le0,\\
& 2y_1+2y_2+4\sigma=1,\\
& y_1,y_2,\sigma\ge0 .
\end{aligned}
$$

*Solution.* The LP optimum is $\sigma^*=\tfrac1{24}$,
$\vy^*=(\tfrac13,\tfrac1{12})$, value $\tfrac{23}{12}\approx1.917$, so
$\vx^*=\vy^*/\sigma^*=(8,2)$ crew-weeks. Check: the ratio at $(8,2)$ is
$(40+6)/(16+4+4)=46/24=\tfrac{23}{12}$. The vertices of $\Omega$ give the
ratios $(0,2)$: $\tfrac34$; $(0,12)$: $\tfrac97$; $(8,2)$: $\tfrac{23}{12}$;
$(8,4)$: $\tfrac{13}7\approx1.857$, so $(8,2)$ is indeed the best. Note that
maximizing the *net* benefit $(5x_1+3x_2)-(2x_1+2x_2+4)$ would choose
$(8,4)$ with \$24 million: the ratio and the difference are different
criteria, and the ratio prefers a smaller program.
:::

(sec-07-02-toolbox)=
## The techniques at a glance

{numref}`tab-07-02-toolbox` collects the techniques of this chapter. "Exact"
means that the linear model has the same solutions as the nonlinear one.

```{list-table} Linearization techniques
:name: tab-07-02-toolbox
:header-rows: 1
:widths: 30 40 30

* - Term
  - Linear model
  - Exact? Binaries?
* - $\delta_1\delta_2$, binaries
  - {eq}`eq-07-02-bin-bin`
  - exact; none added
* - $\delta x$, $x\in[\underline x,\overline x]$
  - {eq}`eq-07-02-bin-cont`
  - exact; none added
* - $x\cdot y$, integer $x$ and bounded $y$
  - binary expansion {eq}`eq-07-02-binary-expansion` + {eq}`eq-07-02-bin-cont`
  - exact; $\lceil\log_2(U+1)\rceil$
* - $xy$, both continuous on a box
  - McCormick {eq}`eq-07-02-mccormick`
  - relaxation; binaries only for partitioning
* - $\max_k(\cdot)\le\eta$, $\abs{f}\le c$, $\min\max$, $\max\min$
  - one row per term {eq}`eq-07-02-max-le`, {eq}`eq-07-02-minmax`
  - exact; none
* - $\max_k(\cdot)\ge\eta$, $\abs{f}\ge c$, $\max\max$, $\min\min$, $\eta=\min$
  - disjunction {eq}`eq-07-02-max-ge`, {eq}`eq-07-02-eta-min`
  - exact; one per term (or one for two terms)
* - $\abs{f}=c$
  - $f=c(2\delta-1)$
  - exact; one
* - $\sum c_i\abs{x_i}$ minimized, $c_i\ge0$
  - epigraph or split form
  - exact; none
* - ratio of affine functions
  - Charnes–Cooper {eq}`eq-07-02-charnes-cooper`
  - exact; none
```

## Summary and Common Mistakes

### Summary

- Products of binaries are linearized exactly by {eq}`eq-07-02-bin-bin`; when
  the product only appears in the objective, the sign of its coefficient
  tells which half is needed.
- A binary times a bounded continuous variable is linearized exactly by the
  four inequalities {eq}`eq-07-02-bin-cont`.
- The McCormick inequalities {eq}`eq-07-02-mccormick` give the convex and
  concave envelopes of $xy$ on a box; they are exact on its boundary, their
  gap is at most $(\overline x-\underline x)(\overline y-\underline y)/4$,
  and partitioning shrinks it ({numref}`fig-07-02-mccormick`).
- Bounding a convex piecewise-linear function from above (or a concave one
  from below) is linear ({eq}`eq-07-02-max-le`, {eq}`eq-07-02-minmax`); the
  opposite direction is a disjunction needing binaries ({eq}`eq-07-02-max-ge`,
  {eq}`eq-07-02-eta-min`, {numref}`fig-07-02-abs-max`).
- Absolute values in a minimized objective with nonnegative weights are
  linear in epigraph or split form; $\abs{f}=c$ is $f=c(2\delta-1)$.
- A ratio of affine functions becomes an LP through the Charnes–Cooper
  transformation {eq}`eq-07-02-charnes-cooper` ({prf:ref}`thm-07-02-charnes-cooper`).

### Common mistakes

- **Dropping the wrong half of a product linearization.** One-sided
  formulations are valid only for the sign of objective coefficient they
  were built for ({prf:ref}`ex-07-02-synergy`).
- **Forgetting $\omega\ge0$** (or $\omega\ge\underline x\delta$) in the
  binary–continuous product.
- **Treating the McCormick envelope as exact.** It is exact only on the
  boundary of the box; inside, the LP optimum overestimates the product
  ({prf:ref}`ex-07-02-mccormick`).
- **Using the epigraph trick on the wrong side.** $\min\max$ and $\max\min$
  are LPs; $\max\max$, $\min\min$, $\max\abs{\cdot}$ and constraints
  $\abs{f}\ge c$ need binaries.
- **Splitting $x=x^+-x^-$ with a negative cost.** The argument that one part
  is zero relies on $c_i>0$.
- **Recovering $\vx$ from Charnes–Cooper as $\vy$.** The solution is
  $\vx=\vy/\sigma$.

## References

- Al-Khayyal, F. A., & Falk, J. E. (1983). Jointly constrained biconvex
  programming. *Mathematics of Operations Research, 8*(2), 273–286.
  <https://doi.org/10.1287/moor.8.2.273>
- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Section 11.4:
  linear fractional programming.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Section 1.3: piecewise linear convex
  objective functions and problems involving absolute values.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Sections
  3.1.7 and 3.2.3: epigraphs and pointwise maxima; Section 4.3.2:
  linear-fractional programming.)
- Charnes, A., & Cooper, W. W. (1962). Programming with linear fractional
  functionals. *Naval Research Logistics Quarterly, 9*(3–4), 181–186.
- Glover, F. (1975). Improved linear integer programming formulations of
  nonlinear integer problems. *Management Science, 22*(4), 455–460.
  <https://doi.org/10.1287/mnsc.22.4.455>
- Glover, F., & Woolsey, E. (1974). Converting the 0–1 polynomial programming
  problem to a 0–1 linear program. *Operations Research, 22*(1), 180–182.
  <https://doi.org/10.1287/opre.22.1.180>
- McCormick, G. P. (1976). Computability of global solutions to factorable
  nonconvex programs: Part I — Convex underestimating problems.
  *Mathematical Programming, 10*(1), 147–175.
  <https://doi.org/10.1007/BF01580665>
- Williams, H. P. (2013). *Model building in mathematical programming*
  (5th ed.). Wiley. (Chapters 9–10: logical conditions, products of
  variables and alternative formulations of integer models.)
