# The Simplex Method

This chapter works through the examples and problems of
{doc}`03-03-simplex-method` and of {doc}`03-04-two-phase-big-m`. Tableaux follow the
conventions of {doc}`03-03-simplex-method`: one row per constraint, labeled
by its basic variable, and the objective row last, labeled $-z$ (or $-\phi$
in Phase I), holding the reduced costs $\bar c_j$ and, in the right-hand-side
(RHS) column, minus the objective value. For a minimization a tableau is optimal
when all $\bar c_j\ge0$; for a maximization, when all $\bar c_j\le0$. Pivot
elements are printed in bold. All tableaux were computed with exact
fractions and every optimum was confirmed with a linear programming (LP)
solver.

(sec-08-02-standard-form)=
## Standard form and transformations

:::{prf:example} Slack and surplus variables
:label: ex-08-02-slack-surplus

In a blending model, $x_1,x_2,x_3,x_4\ge0$ are the amounts of four inputs,
and three quality requirements read

$$
\begin{aligned}
4.5x_1+5.5x_2+4.2x_3+x_4&\ge4.5,\\
4x_1+13x_2+5x_3&\ge35,\\
4x_1+7x_2+8x_3&\le60 .
\end{aligned}
$$

Write the constraints as equations with nonnegative variables and explain
the meaning of each new variable.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

A $\ge$ constraint becomes an equation by **subtracting** a nonnegative
**surplus** variable, which measures by how much the left-hand side exceeds
the lower limit; a $\le$ constraint becomes an equation by **adding** a
nonnegative **slack** variable, the unused part of the upper limit:

$$
\begin{aligned}
4.5x_1+5.5x_2+4.2x_3+x_4-x_5&=4.5,\\
4x_1+13x_2+5x_3-x_6&=35,\\
4x_1+7x_2+8x_3+x_7&=60,\\
x_1,\dots,x_7&\ge0 .
\end{aligned}
$$

Here $x_5$ and $x_6$ are surplus variables and $x_7$ is a slack variable.
A point satisfies the original inequalities exactly when the values of
$x_5,x_6,x_7$ computed from the equations are nonnegative; for example
$\vx=(1,2,1,0)$ gives $x_5=4.5+11+4.2-4.5=15.2$, $x_6=4+26+5-35=0$ (the
second requirement is met exactly) and $x_7=60-4-14-8=34$.

The slack $x_7$ has coefficient $+1$ in its row and appears nowhere else, so
it can serve as a starting basic variable. The surplus variables have
coefficient $-1$ and cannot (for $b_i>0$ they would be negative); rows 1 and
2 need artificial variables to start the simplex method
({doc}`03-04-two-phase-big-m`). *(The inputs are renamed $x_1,\dots,x_4$; a
name such as $x_B$ for an input would clash with the vector of basic
variables.)*
:::

:::{prf:example} A lower bound and a free variable
:label: ex-08-02-bounds-free

Bring the constraints

$$
\begin{aligned}
2x_1+3x_2+4x_3&=10,\\
x_1-x_2+2x_3&=8,\\
x_1&\ge0,\\
x_2&\ge3,\\
x_3&\ \text{free}
\end{aligned}
$$

to the form $\mA\vx=\vb$, $\vx\ge\vzero$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Lower bound.** Shift $x_2$ by its bound: $x_2=x_2'+3$ with $x_2'\ge0$.

**Free variable.** Write $x_3=x_3^+-x_3^-$ with $x_3^+,x_3^-\ge0$ (its
positive and negative parts).

Substituting in both equations,

$$
\begin{aligned}
2x_1+3(x_2'+3)+4(x_3^+-x_3^-)&=10,\\
x_1-(x_2'+3)+2(x_3^+-x_3^-)&=8,
\end{aligned}
$$

and moving the constants to the right:

$$
\begin{aligned}
2x_1+3x_2'+4x_3^+-4x_3^-&=1,\\
x_1-x_2'+2x_3^+-2x_3^-&=11,\\
x_1,x_2',x_3^+,x_3^-&\ge0 .
\end{aligned}
$$

Both right-hand sides are nonnegative, so no row has to be multiplied by
$-1$. A solution of the new system gives $x_2=x_2'+3\ge3$ and
$x_3=x_3^+-x_3^-$, and every solution of the old system arises in this way.
The columns of $x_3^+$ and $x_3^-$ are negatives of each other, so they are
never both basic ({prf:ref}`ex-08-02-a3-q2`).
:::

:::{prf:example} Minimizing the larger of two absolute values
:label: ex-08-02-a3-q1

Show that

$$
\min_{\vx\ge\vzero}\ \max\bigl\{\abs{x_1-x_2+3x_3},\ \abs{-x_1+3x_2-x_3}\bigr\}
$$

can be written as an equivalent LP. (Use that $\abs a\le t$ holds if and
only if $a\le t$ and $-a\le t$.)
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Introduce a variable $t$ for the value of the objective and require it to
bound both absolute values from above:

$$
\begin{aligned}
\min_{\vx,t}\ & t\\
\st\ & x_1-x_2+3x_3-t\le0,\\
& -x_1+x_2-3x_3-t\le0,\\
& -x_1+3x_2-x_3-t\le0,\\
& x_1-3x_2+x_3-t\le0,\\
& x_1,x_2,x_3\ge0 .
\end{aligned}
$$

**Equivalence.** For fixed $\vx$ the four constraints say exactly
$t\ge\abs{x_1-x_2+3x_3}$ and $t\ge\abs{-x_1+3x_2-x_3}$, so the smallest
feasible $t$ is $f(\vx)=\max\{\abs{\cdot},\abs{\cdot}\}$. Minimizing $t$ over
$(\vx,t)$ therefore minimizes $f$ over $\vx$, the optimal values agree, and
$\vx$ is optimal in one problem if and only if $(\vx,f(\vx))$ is optimal in
the other. The constraints imply $t\ge0$, so $t$ may be declared
nonnegative, and four slack variables $s_1,\dots,s_4\ge0$ turn the model
into standard form, for example $x_1-x_2+3x_3-t+s_1=0$.

**Check.** The optimal value is $0$, attained only at $\vx=\vzero$: both
expressions vanish only if $x_1-x_2+3x_3=0$ and $-x_1+3x_2-x_3=0$; adding
them gives $2x_2+2x_3=0$, so $x_2=x_3=0$ and then $x_1=0$. An LP solver
returns $\vx=\vzero$, $t=0$.
:::

:::{prf:example} The two parts of a free variable are never both basic
:label: ex-08-02-a3-q2

For

$$
\begin{aligned}
\max_{\vx}\ & z=2x_1+3x_2+5x_3\\
\st\ & -6x_1+7x_2-9x_3\ge4,\\
& x_1+x_2+4x_3=10,\\
& x_1,x_3\ge0,\quad x_2\ \text{free},
\end{aligned}
$$

write the equation form with $x_2=x_2^+-x_2^-$, $x_2^+,x_2^-\ge0$, and show
that no basic solution has both $x_2^+$ and $x_2^-$ basic.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

With a surplus variable $x_4$ in the first constraint,

$$
\begin{aligned}
\max\ & z=2x_1+3x_2^+-3x_2^-+5x_3\\
\st\ & -6x_1+7x_2^+-7x_2^--9x_3-x_4=4,\\
& x_1+x_2^+-x_2^-+4x_3=10,\\
& x_1,x_2^+,x_2^-,x_3,x_4\ge0 .
\end{aligned}
$$

The columns of $x_2^+$ and $x_2^-$ are $(7,1)^\T$ and $(-7,-1)^\T$. They
are linearly dependent (one is $-1$ times the other), so no basis matrix can
contain both: with $m=2$ rows, the basis $\{x_2^+,x_2^-\}$ would be singular,
and in a larger problem any set of columns containing both is dependent.
Hence in every basic solution at least one of $x_2^+$, $x_2^-$ is
nonbasic, that is, zero, and the other equals $\abs{x_2}$: a basic feasible solution (BFS) always splits
a free variable into its genuine positive and negative parts.

For completeness, the LP has the optimum $x_2=10$, $x_1=x_3=0$, $z^*=30$
(basis $\{x_2^+,x_4\}$ with $x_4=66$).
:::

(sec-08-02-basic-solutions)=
## Basic solutions and adjacency

:::{prf:example} Finding a BFS of the product-mix LP
:label: ex-08-02-product-mix-bfs

For the product-mix LP of {prf:ref}`ex-03-03-product-mix`, with slack
variables $x_3,x_4,x_5$:

$$
\begin{aligned}
\max\ & z=500x_1+450x_2\\
\st\ & 6x_1+5x_2+x_3=60,\\
& 10x_1+20x_2+x_4=150,\\
& x_1+x_5=8,\\
& x_1,\dots,x_5\ge0 .
\end{aligned}
$$

(a) Identify a basic feasible solution directly from the detached
coefficient tableau, and explain why it is a BFS.

(b) Use Gauss–Jordan row operations to obtain the canonical form with
respect to the basic variables $x_1,x_3,x_4$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** The columns of $x_3,x_4,x_5$ form an identity matrix and the
right-hand sides are nonnegative, so the tableau is already in canonical
form with respect to $\{x_3,x_4,x_5\}$:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $6$ | $5$ | $1$ | $0$ | $0$ | $60$ |
| $x_4$ | $10$ | $20$ | $0$ | $1$ | $0$ | $150$ |
| $x_5$ | $\boldsymbol{1}$ | $0$ | $0$ | $0$ | $1$ | $8$ |
| $-z$ | $500$ | $450$ | $0$ | $0$ | $0$ | $0$ |

Set the nonbasic variables to zero, $x_1=x_2=0$; each basic variable then
equals the right-hand side of its row: $\vx=(0,0,60,150,8)$. It is a basic
solution because the basis matrix $\mB=\mI$ is nonsingular, and feasible
because all values are nonnegative; it is the origin of the
$(x_1,x_2)$-plane, with $z=0$.

**(b)** The columns of $x_3$ and $x_4$ are already unit vectors, with their
$1$ in rows 1 and 2. Only $x_1$ must become a unit vector with its $1$ in
the remaining row 3 (the row of $x_5$). Pivot on the bold entry $1$: row 3
stays as it is, and 6, 10 and 500 times row 3 are subtracted from row 1,
row 2 and the objective row:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $0$ | $5$ | $1$ | $0$ | $-6$ | $12$ |
| $x_4$ | $0$ | $20$ | $0$ | $1$ | $-10$ | $70$ |
| $x_1$ | $1$ | $0$ | $0$ | $0$ | $1$ | $8$ |
| $-z$ | $0$ | $450$ | $0$ | $0$ | $-500$ | $-4000$ |

This is the canonical form with respect to $\{x_1,x_3,x_4\}$: BFS
$(x_1,x_3,x_4)=(8,12,70)$, $x_2=x_5=0$, $z=4000$, reduced costs
$\bar c_2=450$ and $\bar c_5=-500$. It agrees with the matrix computation
$\mB^{-1}\vb$, $\mB^{-1}\mN$ in {prf:ref}`ex-03-03-product-mix`. The basis is
feasible because $\bar{\vb}=(12,70,8)\ge\vzero$. The pivot row is forced:
a pivot on the entry 6 or 10 of the $x_1$ column would remove $x_3$ or
$x_4$ from the basis instead of $x_5$, and would also give a negative
right-hand side ($8-10$ or $8-15$).
:::

:::{prf:example} Basic feasible solutions with a free variable
:label: ex-08-02-a3-q3

Consider

$$
\begin{aligned}
\max_{\vx}\ & z=x_1+3x_2\\
\st\ & x_1+2x_2\le2,\\
& -x_1+x_2\le4,\\
& x_1\ \text{free},\quad x_2\ge0 .
\end{aligned}
$$

(a) Determine all basic feasible solutions.

(b) Find the best BFS by evaluating the objective at each.

(c) Solve the problem graphically and confirm the result of (b).
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** Write $x_1=x_1^+-x_1^-$ and add slack variables $x_3,x_4$:

$$
\begin{aligned}
x_1^+-x_1^-+2x_2+x_3&=2,\\
-x_1^++x_1^-+x_2+x_4&=4,
\end{aligned}
$$

with all five variables nonnegative. A basis has two of the five columns,
so there are $\binom52=10$ candidates. The pair $\{x_1^+,x_1^-\}$ is
singular ({prf:ref}`ex-08-02-a3-q2`). Solving $\mB\vx_B=\vb$ for the other
nine gives:

| Basis | $(x_1^+,x_1^-,x_2,x_3,x_4)$ | Feasible | $(x_1,x_2)$ | $z$ |
|:--|:-:|:-:|:-:|:-:|
| $x_1^-,x_2$ | $(0,2,2,0,0)$ | yes | $(-2,2)$ | $4$ |
| $x_2,x_4$ | $(0,0,1,0,3)$ | yes | $(0,1)$ | $3$ |
| $x_1^+,x_4$ | $(2,0,0,0,6)$ | yes | $(2,0)$ | $2$ |
| $x_3,x_4$ | $(0,0,0,2,4)$ | yes | $(0,0)$ | $0$ |
| $x_1^-,x_3$ | $(0,4,0,6,0)$ | yes | $(-4,0)$ | $-4$ |
| $x_1^+,x_2$ | $(-2,0,2,0,0)$ | no | | |
| $x_1^+,x_3$ | $(-4,0,0,6,0)$ | no | | |
| $x_1^-,x_4$ | $(0,-2,0,0,6)$ | no | | |
| $x_2,x_3$ | $(0,0,4,-6,0)$ | no | | |

So there are **five** basic feasible solutions.

**(b)** The largest objective value among them is $z=4$, at
$(x_1,x_2)=(-2,2)$.

**(c)** In the $(x_1,x_2)$-plane the feasible region is the triangle with
vertices $(-4,0)$, $(2,0)$ and $(-2,2)$ ({numref}`fig-08-02-free-variable-region`).
The contour $x_1+3x_2=4$ touches it only at $(-2,2)$, the intersection of
$x_1+2x_2=2$ and $-x_1+x_2=4$, so $\vx^*=(-2,2)$ with $z^*=4$, as in (b).
An LP solver gives the same point.

```{figure} figures/08-02-free-variable-region.svg
:name: fig-08-02-free-variable-region
:width: 480px
:class: dark-light
:alt: Triangle with vertices (-4,0), (2,0) and (-2,2); the optimum (-2,2) and two non-vertex BFS points (0,0) and (0,1) are marked

Feasible region of {prf:ref}`ex-08-02-a3-q3`. The three corner points
correspond to three of the five BFSs; the other two BFSs, $(0,0)$ and
$(0,1)$ (squares), have $x_1^+=x_1^-=0$ and lie on edges. The dashed
contour $x_1+3x_2=4$ touches the region at the optimum $(-2,2)$.
```

Two of the five BFSs, $(0,0)$ and $(0,1)$, are not corner points of the
triangle. In them both $x_1^+$ and $x_1^-$ are nonbasic, which forces
$x_1=0$; that is a vertex of the transformed region in $\R^5$, but in the
$(x_1,x_2)$-plane $x_1=0$ is not a constraint, and the points lie inside
edges. Splitting a free variable creates such extra BFSs; the optimum,
however, is always attained at a genuine corner as well.

*(An earlier version of this solution lists the BFS with $x_1=-1$,
$x_2=3$ as the best one and labels the corner of the graph $(-1,3)$. That
point violates $x_1+2x_2\le2$ (it gives $5$); the correct corner and
optimum are $(-2,2)$ with $z^*=4$, and the BFS list has five entries, not
three.)*
:::

:::{prf:example} Successive corner points of a cut cube
:label: ex-08-02-a3-q4

The feasible region of a three-variable LP is the unit cube
$0\le x_1,x_2,x_3\le1$ with the corner $(1,1,1)$ cut off by the plane
$x_1+x_2+x_3=\tfrac52$ ({numref}`fig-08-02-cut-cube`). Its ten extreme points
are $A=(0,0,0)$, $B=(1,0,0)$, $C=(0,1,0)$, $D=(0,0,1)$, $E=(1,1,0)$,
$F=(0,1,1)$, $G=(1,0,1)$, $H=(1,\tfrac12,1)$, $I=(1,1,\tfrac12)$ and
$J=(\tfrac12,1,1)$.

```{figure} figures/08-02-cut-cube.svg
:name: fig-08-02-cut-cube
:width: 440px
:class: dark-light
:alt: Unit cube with the corner (1,1,1) cut off by a triangle H-I-J; vertices A to J labeled

The solution space of {prf:ref}`ex-08-02-a3-q4` and
{prf:ref}`ex-08-02-a3-q5`: a unit cube whose corner $(1,1,1)$ is replaced by
the triangle $HIJ$. Hidden edges are dashed.
```

(a) Which of the pairs $(A,B)$, $(H,I)$, $(E,H)$ and $(A,I)$ cannot be the
corner points of two successive simplex iterations? Why?

(b) The simplex method starts at $A$, and the optimum is at $H$. Which of
the following paths are impossible for the simplex method, and why?
(i) $A\to B\to G\to H$; (ii) $A\to D\to F\to C\to A\to B\to G\to H$;
(iii) $A\to C\to I\to H$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

A simplex pivot with $\theta>0$ moves along one edge of the polytope to an
**adjacent** extreme point (one that shares an edge), and with $\theta=0$ it
stays put. Two extreme points of a 3-D polytope are adjacent when two
linearly independent constraints are active at both. The edges of the cut
cube give the adjacency lists

| Point | $A$ | $B$ | $C$ | $D$ | $E$ | $F$ | $G$ | $H$ | $I$ | $J$ |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| Adjacent to | $B,C,D$ | $A,E,G$ | $A,E,F$ | $A,F,G$ | $B,C,I$ | $C,D,J$ | $B,D,H$ | $G,I,J$ | $E,H,J$ | $F,H,I$ |

**(a)** $(A,B)$ is the edge along the $x_1$-axis, and $(H,I)$ is an edge
of the cutting triangle; both can be successive iterations. $(E,H)$ and
$(A,I)$ cannot: $E=(1,1,0)$ and $H=(1,\tfrac12,1)$ share only the
constraint $x_1\le1$, a face, not an edge; $A$ and $I$ have no active
constraint in common.

**(b)** (i) $A\to B\to G\to H$ is possible: every step follows an edge,
and an objective such as $z=5x_1+2x_2+4x_3$ increases along it
($0,5,9,10$).
(ii) Impossible. All steps follow edges, but the path returns to $A$. The
cut cube has no degenerate vertex (exactly three constraints are active at
each extreme point), so every pivot strictly increases $z$, and a point can
never be visited twice.
(iii) Impossible: $C$ and $I$ are not adjacent.
:::

:::{prf:example} The first pivot from the origin
:label: ex-08-02-a3-q5

For the solution space of {prf:ref}`ex-08-02-a3-q4`, starting at $A=\vzero$
with the slack basis, find the entering variable chosen by Dantzig's rule,
the next corner point and the improvement in $z$ for

(a) $\max z=x_1-2x_2+3x_3$; (b) $\max z=5x_1+2x_2+4x_3$;
(c) $\max z=-2x_1+7x_2+2x_3$; (d) $\max z=x_1+x_2+x_3$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

At $A$ all slack variables are basic and $\vc_B=\vzero$, so the reduced
costs are the objective coefficients. Dantzig's rule enters the variable
with the largest positive coefficient. Raising $x_j$ from $0$ moves along
the $x_j$-axis; the bound $x_j\le1$ stops it after one unit (the cutting
plane $x_1+x_2+x_3\le\tfrac52$ is not reached), so $\theta=1$ and
$\Delta z=\theta\bar c_j=c_j$.

| Objective | Entering | Next point | $\Delta z$ |
|:--|:-:|:-:|:-:|
| (a) $x_1-2x_2+3x_3$ | $x_3$ | $D$ | $3$ |
| (b) $5x_1+2x_2+4x_3$ | $x_1$ | $B$ | $5$ |
| (c) $-2x_1+7x_2+2x_3$ | $x_2$ | $C$ | $7$ |
| (d) $x_1+x_2+x_3$ | $x_1$, $x_2$ or $x_3$ (tie) | $B$, $C$ or $D$ | $1$ |

In (d) the three coefficients are equal, and any of the three variables may
enter; the usual tie-break (smallest index) takes $x_1$.
:::

:::{prf:example} A one-constraint LP solved by inspection
:label: ex-08-02-a3-q6

Solve by inspection, and justify the answer in terms of basic feasible
solutions:

$$
\begin{aligned}
\max_{\vx}\ & z=5x_1-6x_2+3x_3-5x_4+12x_5\\
\st\ & x_1+3x_2+5x_3+6x_4+3x_5\le30,\\
& x_1,\dots,x_5\ge0 .
\end{aligned}
$$
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

With a slack $x_6$, there is one constraint, so every basis has exactly
**one** basic variable. If $x_j$ is basic ($a_j>0$), the BFS is
$x_j=30/a_j$ with all other variables zero, and $z=30\,c_j/a_j$. The best
BFS uses the variable with the largest ratio $c_j/a_j$ (profit per unit of
the single resource):

| Basic variable | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ (slack) |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| Value | $30$ | $10$ | $6$ | $5$ | $10$ | $30$ |
| $z$ | $150$ | $-60$ | $18$ | $-25$ | $120$ | $0$ |

The optimum is $x_1=30$, $x_2=\dots=x_5=0$, $z^*=150$.

**Optimality check.** For the basis $\{x_1\}$, $\pi=c_1/a_1=5$ and the
reduced costs $\bar c_j=c_j-5a_j$ are $0,-21,-22,-35,-3$ for
$x_1,\dots,x_5$ and $-5$ for the slack: all $\le0$, so by the optimality
test the BFS is optimal (and unique, since all nonbasic reduced costs are
strictly negative). An LP solver agrees.
:::

(sec-08-02-iterations)=
## Simplex iterations

:::{prf:example} Reading a tableau
:label: ex-08-02-a3-q7

The following tableau occurs during the simplex method; all variables are
nonnegative. It is optimal neither for maximization nor for minimization.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | $x_8$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_8$ | $0$ | $3$ | $0$ | $-2$ | $-3$ | $-1$ | $5$ | $1$ | $12$ |
| $x_3$ | $0$ | $1$ | $1$ | $3$ | $1$ | $0$ | $3$ | $0$ | $6$ |
| $x_1$ | $1$ | $-1$ | $0$ | $0$ | $6$ | $-4$ | $0$ | $0$ | $0$ |
| $-z$ | $0$ | $5$ | $0$ | $-4$ | $1$ | $10$ | $0$ | $0$ | $-620$ |

(a) Classify the variables as basic and nonbasic and give the current
values of all variables and of $z$.

(b) Assume the problem is a maximization. Which nonbasic variables can
improve $z$? For each, find the leaving variable (if any) and the change in
$z$, without pivoting.

(c) Repeat (b) for a minimization problem.

(d) Which nonbasic variables leave $z$ unchanged when they enter?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** Basic: $x_8=12$, $x_3=6$, $x_1=0$ (degenerate). Nonbasic, hence
zero: $x_2,x_4,x_5,x_6,x_7$. The right-hand side of the $-z$ row is
$-620$, so $z=620$.

**(b)** For a maximization, a variable improves $z$ if $\bar c_j>0$:
$x_2$ ($5$), $x_5$ ($1$) and $x_6$ ($10$). The change in $z$ is
$\theta\bar c_j$ with $\theta$ from the ratio test.

- $x_2$: positive entries $3$ (row $x_8$) and $1$ (row $x_3$), ratios
  $12/3=4$ and $6/1=6$; $x_8$ leaves, $\theta=4$, $\Delta z=4\cdot5=20$.
- $x_5$: positive entries $1$ (row $x_3$) and $6$ (row $x_1$), ratios $6$
  and $0/6=0$; $x_1$ leaves, $\theta=0$, $\Delta z=0$ (a degenerate pivot).
- $x_6$: the column $(-1,0,-4)^\T$ has no positive entry; no variable
  leaves, and $z$ increases without bound ($\Delta z=+\infty$): the LP is
  unbounded.

**(c)** For a minimization, a variable improves $z$ if $\bar c_j<0$: only
$x_4$ ($-4$). Its only positive entry is $3$ (row $x_3$), $\theta=6/3=2$;
$x_3$ leaves and $\Delta z=2\cdot(-4)=-8$.

**(d)** $x_7$, because $\bar c_7=0$ (it would enter with $\theta=\min\{12/5,\,6/3\}=2$,
replacing $x_3$, and $z$ stays 620), and $x_5$, because its step length is
$\theta=0$.
:::

:::{prf:example} Most positive versus least positive entering rule
:label: ex-08-02-a3-q8

$$
\begin{aligned}
\max_{\vx}\ & z=16x_1+15x_2\\
\st\ & 40x_1+31x_2\le124,\\
& -x_1+x_2\le1,\\
& x_1\le3,\\
& x_1,x_2\ge0 .
\end{aligned}
$$

(a) Solve by the simplex method, entering the nonbasic variable with the
most positive $\bar c_j$.

(b) Solve again, always entering the variable with the least positive
$\bar c_j$.

(c) Compare the numbers of iterations. Does the most-positive rule always
need fewer iterations? What does this say about the optimality condition?

(d) The problem is turned into a minimization by multiplying $z$ by $-1$.
How does this change the iterations?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Slack variables $x_3,x_4,x_5$ give the starting basis.

**(a)** Most positive rule. $x_1$ enters ($16$); ratios $124/40=3.1$,
— , $3/1=3$; $x_5$ leaves.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $40$ | $31$ | $1$ | $0$ | $0$ | $124$ |
| $x_4$ | $-1$ | $1$ | $0$ | $1$ | $0$ | $1$ |
| $x_5$ | $\boldsymbol{1}$ | $0$ | $0$ | $0$ | $1$ | $3$ |
| $-z$ | $16$ | $15$ | $0$ | $0$ | $0$ | $0$ |

$x_2$ enters ($15$); ratios $4/31$ and $4/1$; $x_3$ leaves.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $0$ | $\boldsymbol{31}$ | $1$ | $0$ | $-40$ | $4$ |
| $x_4$ | $0$ | $1$ | $0$ | $1$ | $1$ | $4$ |
| $x_1$ | $1$ | $0$ | $0$ | $0$ | $1$ | $3$ |
| $-z$ | $0$ | $15$ | $0$ | $0$ | $-16$ | $-48$ |

$x_5$ enters ($\tfrac{104}{31}$); ratios: row $x_2$ negative, row $x_4$
$\tfrac{120/31}{71/31}=\tfrac{120}{71}$, row $x_1$ $3$; $x_4$ leaves.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_2$ | $0$ | $1$ | $\tfrac{1}{31}$ | $0$ | $-\tfrac{40}{31}$ | $\tfrac{4}{31}$ |
| $x_4$ | $0$ | $0$ | $-\tfrac{1}{31}$ | $1$ | $\boldsymbol{\tfrac{71}{31}}$ | $\tfrac{120}{31}$ |
| $x_1$ | $1$ | $0$ | $0$ | $0$ | $1$ | $3$ |
| $-z$ | $0$ | $0$ | $-\tfrac{15}{31}$ | $0$ | $\tfrac{104}{31}$ | $-\tfrac{1548}{31}$ |

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_2$ | $0$ | $1$ | $\tfrac{1}{71}$ | $\tfrac{40}{71}$ | $0$ | $\tfrac{164}{71}$ |
| $x_5$ | $0$ | $0$ | $-\tfrac{1}{71}$ | $\tfrac{31}{71}$ | $1$ | $\tfrac{120}{71}$ |
| $x_1$ | $1$ | $0$ | $\tfrac{1}{71}$ | $-\tfrac{31}{71}$ | $0$ | $\tfrac{93}{71}$ |
| $-z$ | $0$ | $0$ | $-\tfrac{31}{71}$ | $-\tfrac{104}{71}$ | $0$ | $-\tfrac{3948}{71}$ |

Optimal after **3 iterations**: $x_1^*=\tfrac{93}{71}\approx1.310$,
$x_2^*=\tfrac{164}{71}\approx2.310$, $z^*=\tfrac{3948}{71}\approx55.61$.

**(b)** Least positive rule. $x_2$ enters ($15<16$); ratios $124/31=4$ and
$1/1=1$; $x_4$ leaves.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $40$ | $31$ | $1$ | $0$ | $0$ | $124$ |
| $x_4$ | $-1$ | $\boldsymbol{1}$ | $0$ | $1$ | $0$ | $1$ |
| $x_5$ | $1$ | $0$ | $0$ | $0$ | $1$ | $3$ |
| $-z$ | $16$ | $15$ | $0$ | $0$ | $0$ | $0$ |

$x_1$ enters (the only positive entry, $31$); ratios $93/71$ and $3$;
$x_3$ leaves.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $\boldsymbol{71}$ | $0$ | $1$ | $-31$ | $0$ | $93$ |
| $x_2$ | $-1$ | $1$ | $0$ | $1$ | $0$ | $1$ |
| $x_5$ | $1$ | $0$ | $0$ | $0$ | $1$ | $3$ |
| $-z$ | $31$ | $0$ | $0$ | $-15$ | $0$ | $-15$ |

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $\tfrac{1}{71}$ | $-\tfrac{31}{71}$ | $0$ | $\tfrac{93}{71}$ |
| $x_2$ | $0$ | $1$ | $\tfrac{1}{71}$ | $\tfrac{40}{71}$ | $0$ | $\tfrac{164}{71}$ |
| $x_5$ | $0$ | $0$ | $-\tfrac{1}{71}$ | $\tfrac{31}{71}$ | $1$ | $\tfrac{120}{71}$ |
| $-z$ | $0$ | $0$ | $-\tfrac{31}{71}$ | $-\tfrac{104}{71}$ | $0$ | $-\tfrac{3948}{71}$ |

Optimal after **2 iterations**, with the same optimum.

**(c)** Here the least-positive rule needs fewer pivots (2 against 3). The
most-positive rule picks the largest improvement *per unit* of the entering
variable, but the step length $\theta$ and the number of remaining pivots
are not taken into account, so it is a heuristic and does not guarantee the
fewest iterations. The optimality condition does not depend on the entering
rule at all: a tableau of a maximization problem is optimal when all
$\bar c_j\le0$ in the $-z$ row, whichever path led to it; both rules end
with the same tableau (rows in a different order).
*(An earlier version of this solution reports 4 and 3 iterations; it counts
tableaux, not pivots. It also states the optimality condition as
"nonnegative $z$-row coefficients", which is the same condition in the
$z$-row convention, $-\bar c_j\ge0$.)*

**(d)** Minimizing $-z=-16x_1-15x_2$ negates every entry of the objective
row. The entering variable is now the one with the most negative entry,
which is the same variable as before, and the ratio tests are unchanged,
so the sequence of bases is identical: 3 pivots with the most-negative
rule, ending with the objective row $(0,0,\tfrac{31}{71},\tfrac{104}{71},0\mid\tfrac{3948}{71})$,
all entries $\ge0$, and optimal value $\min(-z)=-\tfrac{3948}{71}$.
:::

:::{prf:example} A given starting basis that needs pricing out
:label: ex-08-02-a3-q9

Solve, starting with $x_3$ and $x_4$ as basic variables:

$$
\begin{aligned}
\min_{\vx}\ & z=3x_1+2x_2+3x_3\\
\st\ & x_1+4x_2+x_3\ge14,\\
& 2x_1+x_2+x_4\ge20,\\
& x_1,x_2,x_3,x_4\ge0 .
\end{aligned}
$$
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

With surplus variables $x_5,x_6$ the constraints are
$x_1+4x_2+x_3-x_5=14$ and $2x_1+x_2+x_4-x_6=20$. The variable $x_3$ occurs
only in row 1 and $x_4$ only in row 2, both with coefficient $+1$, so
$\{x_3,x_4\}$ is a basis with the BFS $x_3=14$, $x_4=20$; no artificial
variable is needed. But $x_3$ has cost $3$, so the objective row
$(3,2,3,0,0,0\mid0)$ is not in canonical form:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $1$ | $4$ | $1$ | $0$ | $-1$ | $0$ | $14$ |
| $x_4$ | $2$ | $1$ | $0$ | $1$ | $0$ | $-1$ | $20$ |
| $-z$ | $3$ | $2$ | $3$ | $0$ | $0$ | $0$ | $0$ |

Subtract 3 times row 1 from the objective row:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $1$ | $\boldsymbol{4}$ | $1$ | $0$ | $-1$ | $0$ | $14$ |
| $x_4$ | $2$ | $1$ | $0$ | $1$ | $0$ | $-1$ | $20$ |
| $-z$ | $0$ | $-10$ | $0$ | $0$ | $3$ | $0$ | $-42$ |

The current value is $z=42$. The only negative reduced cost is
$\bar c_2=-10$; ratios $14/4=\tfrac72$ and $20/1=20$; $x_3$ leaves.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_2$ | $\tfrac{1}{4}$ | $1$ | $\tfrac{1}{4}$ | $0$ | $-\tfrac{1}{4}$ | $0$ | $\tfrac{7}{2}$ |
| $x_4$ | $\tfrac{7}{4}$ | $0$ | $-\tfrac{1}{4}$ | $1$ | $\tfrac{1}{4}$ | $-1$ | $\tfrac{33}{2}$ |
| $-z$ | $\tfrac{5}{2}$ | $0$ | $\tfrac{5}{2}$ | $0$ | $\tfrac{1}{2}$ | $0$ | $-7$ |

All reduced costs are $\ge0$: optimal, with $x_2^*=\tfrac72$,
$x_4^*=\tfrac{33}{2}$, $x_1^*=x_3^*=0$ and $z^*=7$ (check:
$4\cdot\tfrac72=14$, $\tfrac72+\tfrac{33}{2}=20$, $2\cdot\tfrac72=7$).
:::

(sec-08-02-artificial)=
## Artificial variables: two-phase and big-M

:::{prf:example} Two-phase method with free variables and mixed constraints
:label: ex-08-02-two-phase-free

Consider

$$
\begin{aligned}
\max_{\vy}\ & z=-3y_1+2y_2-y_3+4y_4\\
\st\ & y_1+y_2-4y_3+2y_4\ge4,\\
& -3y_1+y_2-2y_3\le6,\\
& y_2-y_4=-1,\\
& y_1+y_2-y_3=0,\\
& y_3,y_4\ge0,\quad y_1,y_2\ \text{free}.
\end{aligned}
$$

(a) Bring the LP to standard form with $y_1=x_1-x_2$, $y_2=x_3-x_4$,
$y_3=x_5$, $y_4=x_6$. Is the resulting tableau in canonical form? Which
variables could form a basis?

(b) Construct an LP whose optimal BFS is a BFS of the original problem.

(c) Solve the problem by the two-phase method.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** Substitute, add a surplus variable $x_7$ to the first constraint and
a slack variable $x_8$ to the second, and multiply the third constraint by
$-1$ so that its right-hand side is nonnegative:

$$
\begin{aligned}
\max\ & z=-3x_1+3x_2+2x_3-2x_4-x_5+4x_6\\
\st\ & x_1-x_2+x_3-x_4-4x_5+2x_6-x_7=4,\\
& -3x_1+3x_2+x_3-x_4-2x_5+x_8=6,\\
& -x_3+x_4+x_6=1,\\
& x_1-x_2+x_3-x_4-x_5=0,\\
& x_1,\dots,x_8\ge0 .
\end{aligned}
$$

This is not in canonical form. Only row 2 has a variable that can serve as
its basic variable (the slack $x_8$). In row 1 the surplus $x_7$ has
coefficient $-1$; row 3 contains $x_6$ with $+1$, but $x_6$ also appears in
row 1; row 4 has no candidate at all. So only $x_8$ is available, and three
more basic variables are missing.

**(b)** Add artificial variables $R_1,R_3,R_4$ to rows 1, 3 and 4. Then
$\{R_1,x_8,R_3,R_4\}$ is an identity basis with BFS
$R_1=4$, $x_8=6$, $R_3=1$, $R_4=0$. A solution of the original LP is a
solution of the enlarged system with $R_1=R_3=R_4=0$, so the **Phase I** LP

$$
\min\ \phi=R_1+R_3+R_4
$$

over the enlarged constraints has optimal value $0$ exactly when the
original LP is feasible, and its optimal BFS (after artificial variables at
zero level are removed) is a BFS of the original LP
({prf:ref}`thm-03-04-phase-one`).

**(c)** *Phase I.* Price out the $-\phi$ row by subtracting rows 1, 3 and 4
from $(0,\dots,0,1,1,1\mid0)$; it becomes
$(-2,2,-1,1,5,-3,1,0,0,0,0\mid-5)$. The $-z$ row
$(-3,3,2,-2,-1,4,0,0,\dots\mid0)$ is carried along.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | $x_8$ | $R_1$ | $R_3$ | $R_4$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $R_1$ | $1$ | $-1$ | $1$ | $-1$ | $-4$ | $2$ | $-1$ | $0$ | $1$ | $0$ | $0$ | $4$ |
| $x_8$ | $-3$ | $3$ | $1$ | $-1$ | $-2$ | $0$ | $0$ | $1$ | $0$ | $0$ | $0$ | $6$ |
| $R_3$ | $0$ | $0$ | $-1$ | $1$ | $0$ | $\boldsymbol{1}$ | $0$ | $0$ | $0$ | $1$ | $0$ | $1$ |
| $R_4$ | $1$ | $-1$ | $1$ | $-1$ | $-1$ | $0$ | $0$ | $0$ | $0$ | $0$ | $1$ | $0$ |
| $-z$ | $-3$ | $3$ | $2$ | $-2$ | $-1$ | $4$ | $0$ | $0$ | $0$ | $0$ | $0$ | $0$ |
| $-\phi$ | $-2$ | $2$ | $-1$ | $1$ | $5$ | $-3$ | $1$ | $0$ | $0$ | $0$ | $0$ | $-5$ |

$x_6$ enters ($-3$); ratios $4/2=2$ and $1/1=1$; $R_3$ leaves.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | $x_8$ | $R_1$ | $R_3$ | $R_4$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $R_1$ | $1$ | $-1$ | $3$ | $-3$ | $-4$ | $0$ | $-1$ | $0$ | $1$ | $-2$ | $0$ | $2$ |
| $x_8$ | $-3$ | $3$ | $1$ | $-1$ | $-2$ | $0$ | $0$ | $1$ | $0$ | $0$ | $0$ | $6$ |
| $x_6$ | $0$ | $0$ | $-1$ | $1$ | $0$ | $1$ | $0$ | $0$ | $0$ | $1$ | $0$ | $1$ |
| $R_4$ | $1$ | $-1$ | $\boldsymbol{1}$ | $-1$ | $-1$ | $0$ | $0$ | $0$ | $0$ | $0$ | $1$ | $0$ |
| $-z$ | $-3$ | $3$ | $6$ | $-6$ | $-1$ | $0$ | $0$ | $0$ | $0$ | $-4$ | $0$ | $-4$ |
| $-\phi$ | $-2$ | $2$ | $-4$ | $4$ | $5$ | $0$ | $1$ | $0$ | $0$ | $3$ | $0$ | $-2$ |

$x_3$ enters ($-4$); the ratios are $2/3$ (row $R_1$), $6/1$ (row $x_8$)
and $0/1=0$ (row $R_4$); $R_4$ leaves with $\theta=0$, a degenerate pivot.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | $x_8$ | $R_1$ | $R_3$ | $R_4$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $R_1$ | $-2$ | $\boldsymbol{2}$ | $0$ | $0$ | $-1$ | $0$ | $-1$ | $0$ | $1$ | $-2$ | $-3$ | $2$ |
| $x_8$ | $-4$ | $4$ | $0$ | $0$ | $-1$ | $0$ | $0$ | $1$ | $0$ | $0$ | $-1$ | $6$ |
| $x_6$ | $1$ | $-1$ | $0$ | $0$ | $-1$ | $1$ | $0$ | $0$ | $0$ | $1$ | $1$ | $1$ |
| $x_3$ | $1$ | $-1$ | $1$ | $-1$ | $-1$ | $0$ | $0$ | $0$ | $0$ | $0$ | $1$ | $0$ |
| $-z$ | $-9$ | $9$ | $0$ | $0$ | $5$ | $0$ | $0$ | $0$ | $0$ | $-4$ | $-6$ | $-4$ |
| $-\phi$ | $2$ | $-2$ | $0$ | $0$ | $1$ | $0$ | $1$ | $0$ | $0$ | $3$ | $4$ | $-2$ |

$x_2$ enters ($-2$); the ratios are $2/2=1$ (row $R_1$) and $6/4=\tfrac32$
(row $x_8$); $R_1$ leaves.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | $x_8$ | $R_1$ | $R_3$ | $R_4$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_2$ | $-1$ | $1$ | $0$ | $0$ | $-\tfrac{1}{2}$ | $0$ | $-\tfrac{1}{2}$ | $0$ | $\tfrac{1}{2}$ | $-1$ | $-\tfrac{3}{2}$ | $1$ |
| $x_8$ | $0$ | $0$ | $0$ | $0$ | $1$ | $0$ | $2$ | $1$ | $-2$ | $4$ | $5$ | $2$ |
| $x_6$ | $0$ | $0$ | $0$ | $0$ | $-\tfrac{3}{2}$ | $1$ | $-\tfrac{1}{2}$ | $0$ | $\tfrac{1}{2}$ | $0$ | $-\tfrac{1}{2}$ | $2$ |
| $x_3$ | $0$ | $0$ | $1$ | $-1$ | $-\tfrac{3}{2}$ | $0$ | $-\tfrac{1}{2}$ | $0$ | $\tfrac{1}{2}$ | $-1$ | $-\tfrac{1}{2}$ | $1$ |
| $-z$ | $0$ | $0$ | $0$ | $0$ | $\tfrac{19}{2}$ | $0$ | $\tfrac{9}{2}$ | $0$ | $-\tfrac{9}{2}$ | $5$ | $\tfrac{15}{2}$ | $-13$ |
| $-\phi$ | $0$ | $0$ | $0$ | $0$ | $0$ | $0$ | $0$ | $0$ | $1$ | $1$ | $1$ | $0$ |

All $-\phi$ entries are $\ge0$, $\phi^*=0$, and no artificial variable is
basic: the basis $\{x_2,x_8,x_6,x_3\}$ is feasible for the original LP.

*Phase II.* Delete the artificial columns and the $-\phi$ row. In the $-z$
row $x_5$ ($\tfrac{19}{2}$) and $x_7$ ($\tfrac92$) are improving; Dantzig's
rule takes $x_5$, whose column has a positive entry only in row $x_8$
(ratio $2/1=2$), so $x_8$ leaves.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | $x_8$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_2$ | $-1$ | $1$ | $0$ | $0$ | $0$ | $0$ | $\tfrac{1}{2}$ | $\tfrac{1}{2}$ | $2$ |
| $x_5$ | $0$ | $0$ | $0$ | $0$ | $1$ | $0$ | $2$ | $1$ | $2$ |
| $x_6$ | $0$ | $0$ | $0$ | $0$ | $0$ | $1$ | $\tfrac{5}{2}$ | $\tfrac{3}{2}$ | $5$ |
| $x_3$ | $0$ | $0$ | $1$ | $-1$ | $0$ | $0$ | $\tfrac{5}{2}$ | $\tfrac{3}{2}$ | $4$ |
| $-z$ | $0$ | $0$ | $0$ | $0$ | $0$ | $0$ | $-\tfrac{29}{2}$ | $-\tfrac{19}{2}$ | $-32$ |

All reduced costs are $\le0$: optimal. The optimum is $x_2=2$, $x_3=4$,
$x_5=2$, $x_6=5$ (the others zero), that is

$$
\begin{aligned}
y_1&=x_1-x_2=-2,\\
y_2&=x_3-x_4=4,\\
y_3&=2,\\
y_4&=5,\\
z^*&=32 .
\end{aligned}
$$

Check: $-2+4-8+10=4\ge4$, $6+4-4=6\le6$, $4-5=-1$, $-2+4-2=0$, and
$z=6+8-2+20=32$. An LP solver on the original $\vy$-model gives the same.

*(An earlier version of this example has a tableau with $-1$ instead of
$-4$ under $x_5$ in row 1 and $6$ instead of $1$ under $x_6$ in row 3, and
the priced-out $-w$ row is given as $(2,-2,1,-1,-2,8,-1,0\mid5)$ for
$\max w=-\phi$. Those coefficients do not match the model; with them the LP
would even be unbounded. The model above gives the row
$(2,-2,1,-1,-5,3,-1,0\mid5)$ in that convention, which is the negative of
the $-\phi$ row used here.)*
:::

:::{prf:example} Variables that must stay zero after Phase I
:label: ex-08-02-a3-q10

For

$$
\begin{aligned}
\max_{\vx}\ & z=3x_1+2x_2+3x_3\\
\st\ & 2x_1+x_2+x_3\le2,\\
& 3x_1+4x_2+2x_3\ge8,\\
& x_1,x_2,x_3\ge0,
\end{aligned}
$$

with slack $x_5$ in the first constraint, surplus $x_4$ and artificial $R$
in the second, Phase I ($\min\phi=R$) ends with the tableau

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $R$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_2$ | $2$ | $1$ | $1$ | $0$ | $1$ | $0$ | $2$ |
| $R$ | $-5$ | $0$ | $-2$ | $-1$ | $-4$ | $1$ | $0$ |
| $-\phi$ | $5$ | $0$ | $2$ | $1$ | $4$ | $0$ | $0$ |

Explain why the nonbasic variables $x_1,x_3,x_4,x_5$ can never become
positive in Phase II, conclude that their columns can be dropped before
Phase II, and find the optimal solution.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The $-\phi$ row is the equation $-\phi+5x_1+2x_3+x_4+4x_5=0$, valid for
every solution of the constraints. Hence

$$
\phi=R=5x_1+2x_3+x_4+4x_5 .
$$

Every feasible solution of the original LP has $R=0$
({prf:ref}`lem-03-04-artificial`), and since all four variables are
nonnegative with positive coefficients, $R=0$ forces
$x_1=x_3=x_4=x_5=0$. They are zero in **every** feasible solution, so they
can never become positive in Phase II and their columns (and $R$) can be
deleted. What remains of the constraints is the first row, $x_2=2$; the
second row becomes $0=0$.

The feasible set is therefore the single point $\vx=(0,2,0)$, with
$z^*=4$, and it is **not necessary** to carry out Phase II at all. Check:
$2\le2$ and $4\cdot2=8\ge8$, both constraints binding; an LP solver gives
the same.

The tableau arises after one Phase I pivot: in the starting tableau the
$-\phi$ row is $(-3,-4,-2,1,0,0\mid-8)$, $x_2$ enters, and the ratios
$2/1$ and $8/4$ tie; with $x_5$ leaving, $R$ stays basic at zero level.
*(An earlier version of this problem reads "it is necessary to carry out
Phase II"; the conclusion is the opposite.)*
:::

:::{prf:example} One artificial variable instead of two
:label: ex-08-02-a3-q11

For

$$
\begin{aligned}
\min_{\vx}\ & z=2x_1-4x_2+3x_3\\
\st\ & 5x_1-6x_2+2x_3\ge5,\\
& -x_1+3x_2+5x_3\ge8,\\
& 2x_1+5x_2-4x_3\le4,\\
& x_1,x_2,x_3\ge0,
\end{aligned}
$$

show how the constraints can be rewritten as equations that need only one
artificial variable, and solve the problem.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

With surplus variables $x_4,x_5$ and slack $x_6$:

$$
\begin{aligned}
5x_1-6x_2+2x_3-x_4&=5 &&(1)\\
-x_1+3x_2+5x_3-x_5&=8 &&(2)\\
2x_1+5x_2-4x_3+x_6&=4 &&(3)
\end{aligned}
$$

Row (3) has its basic variable $x_6$. The two $\ge$ rows would each need an
artificial variable. Instead, put one artificial variable $R$ into the
$\ge$ row with the **largest** right-hand side, row (2), and replace row (1)
by (2) $-$ (1):

$$
-6x_1+9x_2+3x_3+x_4-x_5=3 .
$$

Now $x_4$ has coefficient $+1$ in this row and appears in no other, and the
right-hand side $8-5=3$ is nonnegative because row (2) had the largest
right-hand side. The starting basis is $\{x_4,R,x_6\}$ with
$(x_4,R,x_6)=(3,8,4)$, and only one artificial variable is needed.

*Phase I* ($\min\phi=R$; the $-\phi$ row is minus row (2)):

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $R$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_4$ | $-6$ | $9$ | $\boldsymbol{3}$ | $1$ | $-1$ | $0$ | $0$ | $3$ |
| $R$ | $-1$ | $3$ | $5$ | $0$ | $-1$ | $0$ | $1$ | $8$ |
| $x_6$ | $2$ | $5$ | $-4$ | $0$ | $0$ | $1$ | $0$ | $4$ |
| $-z$ | $2$ | $-4$ | $3$ | $0$ | $0$ | $0$ | $0$ | $0$ |
| $-\phi$ | $1$ | $-3$ | $-5$ | $0$ | $1$ | $0$ | $0$ | $-8$ |

$x_3$ enters ($-5$); ratios $3/3=1$ and $8/5$; $x_4$ leaves.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $R$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $-2$ | $3$ | $1$ | $\tfrac{1}{3}$ | $-\tfrac{1}{3}$ | $0$ | $0$ | $1$ |
| $R$ | $\boldsymbol{9}$ | $-12$ | $0$ | $-\tfrac{5}{3}$ | $\tfrac{2}{3}$ | $0$ | $1$ | $3$ |
| $x_6$ | $-6$ | $17$ | $0$ | $\tfrac{4}{3}$ | $-\tfrac{4}{3}$ | $1$ | $0$ | $8$ |
| $-z$ | $8$ | $-13$ | $0$ | $-1$ | $1$ | $0$ | $0$ | $-3$ |
| $-\phi$ | $-9$ | $12$ | $0$ | $\tfrac{5}{3}$ | $-\tfrac{2}{3}$ | $0$ | $0$ | $-3$ |

$x_1$ enters ($-9$); only row $R$ has a positive entry; $R$ leaves.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $R$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $0$ | $\tfrac{1}{3}$ | $1$ | $-\tfrac{1}{27}$ | $-\tfrac{5}{27}$ | $0$ | $\tfrac{2}{9}$ | $\tfrac{5}{3}$ |
| $x_1$ | $1$ | $-\tfrac{4}{3}$ | $0$ | $-\tfrac{5}{27}$ | $\tfrac{2}{27}$ | $0$ | $\tfrac{1}{9}$ | $\tfrac{1}{3}$ |
| $x_6$ | $0$ | $9$ | $0$ | $\tfrac{2}{9}$ | $-\tfrac{8}{9}$ | $1$ | $\tfrac{2}{3}$ | $10$ |
| $-z$ | $0$ | $-\tfrac{7}{3}$ | $0$ | $\tfrac{13}{27}$ | $\tfrac{11}{27}$ | $0$ | $-\tfrac{8}{9}$ | $-\tfrac{17}{3}$ |
| $-\phi$ | $0$ | $0$ | $0$ | $0$ | $0$ | $0$ | $1$ | $0$ |

$\phi^*=0$ and $R$ is nonbasic. *Phase II:* delete $R$ and the $-\phi$
row; $x_2$ is the only variable with a negative reduced cost
($-\tfrac73$), the ratios are $\tfrac53/\tfrac13=5$ (row $x_3$) and
$\tfrac{10}{9}$ (row $x_6$), and $x_6$ leaves:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $0$ | $0$ | $1$ | $-\tfrac{11}{243}$ | $-\tfrac{37}{243}$ | $-\tfrac{1}{27}$ | $\tfrac{35}{27}$ |
| $x_1$ | $1$ | $0$ | $0$ | $-\tfrac{37}{243}$ | $-\tfrac{14}{243}$ | $\tfrac{4}{27}$ | $\tfrac{49}{27}$ |
| $x_2$ | $0$ | $1$ | $0$ | $\tfrac{2}{81}$ | $-\tfrac{8}{81}$ | $\tfrac{1}{9}$ | $\tfrac{10}{9}$ |
| $-z$ | $0$ | $0$ | $0$ | $\tfrac{131}{243}$ | $\tfrac{43}{243}$ | $\tfrac{7}{27}$ | $-\tfrac{83}{27}$ |

All reduced costs are $\ge0$, so the optimum is
$\vx^*=(\tfrac{49}{27},\tfrac{10}{9},\tfrac{35}{27})$ with
$z^*=\tfrac{83}{27}\approx3.074$, which an LP solver confirms. All three
constraints are binding ($x_4=x_5=x_6=0$); for instance
$5\cdot\tfrac{49}{27}-6\cdot\tfrac{30}{27}+2\cdot\tfrac{35}{27}=\tfrac{135}{27}=5$.

*(An earlier version of this solution puts $R$ in row (1), solves row (2)
for $x_2$ and starts from the basis $(R,x_2,x_6)$. That basis is not
feasible: it gives $R=21$, $x_2=\tfrac83$ and $x_6=-\tfrac{28}{3}<0$. The
construction above, with $R$ in the row with the largest right-hand side,
always yields a feasible start.)*
:::

(sec-08-02-special-cases)=
## Degeneracy, alternative optima, unboundedness and infeasibility

:::{prf:example} Simplex path through a degenerate vertex
:label: ex-08-02-a3-q12

In the feasible region of {numref}`fig-08-02-degenerate-vertex`, four
constraint lines pass through the corner $C$. The simplex method starts at
$A$, the optimum is at $D$, and the objective is such that $x_1$ enters
first at $A$.

```{figure} figures/08-02-degenerate-vertex.svg
:name: fig-08-02-degenerate-vertex
:width: 480px
:class: dark-light
:alt: Pentagon A-B-C-D-E with four constraint lines passing through the corner C

A pentagon $ABCDE$ in which the corner $C=(8,2)$ is degenerate: the lines
$x_1-2x_2=4$, $4x_1+x_2=34$, $2x_1-x_2=14$ and $4x_1-3x_2=26$ all pass
through it. The other corners are $A=(0,0)$, $B=(4,0)$, $D=(7,6)$ and
$E=(0,4)$; the dashed contour $2x_1+x_2=20$ touches the region at $D$.
```

(a) Identify the corner points on the simplex path to the optimum.

(b) What is the largest possible number of simplex iterations to reach the
optimum, assuming no cycling?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The region is
$\{x_1-2x_2\le4,\ 4x_1+x_2\le34,\ -2x_1+7x_2\le28,\ 2x_1-x_2\le14,\ 4x_1-3x_2\le26,\ \vx\ge\vzero\}$
with slack variables $x_3,\dots,x_7$ in this order; the last two
constraints are redundant but touch the region at $C$.

**(a)** Entering $x_1$ at $A$ moves along the $x_1$-axis to $B$. From $B$
the only adjacent corner other than $A$ is $C$, and returning to $A$ would
lower $z$. From $C$ the method reaches $D$. The path is
$A\to B\to C\to D$.

**(b)** Count the bases that describe each corner. At $A$, $B$ and $D$
exactly two constraints (counting the axes) are active, so each is described
by one basis. At $C$ four slack variables are zero: $x_3,x_4,x_6,x_7$. A
basis there contains $x_1$, $x_2$, the positive slack $x_5$ and two of the
four zero slacks, and since any two of the four lines through $C$
intersect, every choice works: $\binom42=6$ bases. Without cycling, no basis
is repeated, so the method visits at most $1+1+6+1=9$ bases, which means at
most **8 iterations** (pivots): one from $A$ to $B$, one from $B$ to $C$,
five degenerate pivots among the bases of $C$, and one from $C$ to $D$.
*(An earlier version of this solution gives 9; that is the number of bases,
that is, of tableaux, not of iterations.)*

For example, with $\max z=2x_1+x_2$ and Dantzig's rule the method needs only
three pivots: $x_1$ replaces $x_3$ ($A\to B$), then $x_2$ enters with a
three-way tie in the ratio test (rows $x_4$, $x_6$, $x_7$, ratio 2), and
after $x_4$ leaves the BFS at $C$ is degenerate:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $\tfrac{1}{9}$ | $\tfrac{2}{9}$ | $0$ | $0$ | $0$ | $8$ |
| $x_2$ | $0$ | $1$ | $-\tfrac{4}{9}$ | $\tfrac{1}{9}$ | $0$ | $0$ | $0$ | $2$ |
| $x_5$ | $0$ | $0$ | $\boldsymbol{\tfrac{10}{3}}$ | $-\tfrac{1}{3}$ | $1$ | $0$ | $0$ | $30$ |
| $x_6$ | $0$ | $0$ | $-\tfrac{2}{3}$ | $-\tfrac{1}{3}$ | $0$ | $1$ | $0$ | $0$ |
| $x_7$ | $0$ | $0$ | $-\tfrac{16}{9}$ | $-\tfrac{5}{9}$ | $0$ | $0$ | $1$ | $0$ |
| $-z$ | $0$ | $0$ | $\tfrac{2}{9}$ | $-\tfrac{5}{9}$ | $0$ | $0$ | $0$ | $-18$ |

In the next pivot $x_3$ enters ($\tfrac29$); the rows $x_6$ and $x_7$ at
zero level have negative entries in its column and impose no limit, the
ratios are $8/\tfrac19=72$ (row $x_1$) and $30/\tfrac{10}{3}=9$ (row
$x_5$), so $x_5$ leaves with $\theta=9$ and the method moves to $D=(7,6)$
with $z^*=20$.
:::

:::{prf:example} Temporary degeneracy
:label: ex-08-02-a3-q13

$$
\begin{aligned}
\max_{\vx}\ & z=3x_1+2x_2\\
\st\ & 4x_1-x_2\le4,\\
& 4x_1+3x_2\le6,\\
& 4x_1+x_2\le4,\\
& x_1,x_2\ge0 .
\end{aligned}
$$

(a) Show that the simplex iterations are temporarily degenerate.

(b) Verify the result graphically.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** Slack variables $x_3,x_4,x_5$.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $\boldsymbol{4}$ | $-1$ | $1$ | $0$ | $0$ | $4$ |
| $x_4$ | $4$ | $3$ | $0$ | $1$ | $0$ | $6$ |
| $x_5$ | $4$ | $1$ | $0$ | $0$ | $1$ | $4$ |
| $-z$ | $3$ | $2$ | $0$ | $0$ | $0$ | $0$ |

$x_1$ enters; the ratios $4/4$, $6/4$ and $4/4$ **tie** between rows $x_3$
and $x_5$. With $x_3$ leaving, $x_5$ stays basic at value $0$:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $-\tfrac{1}{4}$ | $\tfrac{1}{4}$ | $0$ | $0$ | $1$ |
| $x_4$ | $0$ | $4$ | $-1$ | $1$ | $0$ | $2$ |
| $x_5$ | $0$ | $\boldsymbol{2}$ | $-1$ | $0$ | $1$ | $0$ |
| $-z$ | $0$ | $\tfrac{11}{4}$ | $-\tfrac{3}{4}$ | $0$ | $0$ | $-3$ |

The BFS $(x_1,x_2)=(1,0)$ is degenerate ($x_5=0$). Next $x_2$ enters
($\tfrac{11}{4}$); the ratios are $2/4$ (row $x_4$) and $0/2=0$ (row
$x_5$), so $x_5$ leaves with $\theta=0$: a **degenerate pivot**, the point
and $z=3$ do not change.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $\tfrac{1}{8}$ | $0$ | $\tfrac{1}{8}$ | $1$ |
| $x_4$ | $0$ | $0$ | $\boldsymbol{1}$ | $1$ | $-2$ | $2$ |
| $x_2$ | $0$ | $1$ | $-\tfrac{1}{2}$ | $0$ | $\tfrac{1}{2}$ | $0$ |
| $-z$ | $0$ | $0$ | $\tfrac{5}{8}$ | $0$ | $-\tfrac{11}{8}$ | $-3$ |

Now $x_3$ enters ($\tfrac58$); the ratios are $1/\tfrac18=8$ (row $x_1$)
and $2/1=2$ (row $x_4$), so $x_4$ leaves with $\theta=2$ and the pivot moves
to a new point:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $0$ | $-\tfrac{1}{8}$ | $\tfrac{3}{8}$ | $\tfrac{3}{4}$ |
| $x_3$ | $0$ | $0$ | $1$ | $1$ | $-2$ | $2$ |
| $x_2$ | $0$ | $1$ | $0$ | $\tfrac{1}{2}$ | $-\tfrac{1}{2}$ | $1$ |
| $-z$ | $0$ | $0$ | $0$ | $-\tfrac{5}{8}$ | $-\tfrac{1}{8}$ | $-\tfrac{17}{4}$ |

The final tableau is optimal and nondegenerate:
$\vx^*=(\tfrac34,1)$, $z^*=\tfrac{17}{4}$. The degeneracy occurred only in
an intermediate tableau; such **temporary degeneracy** costs one pivot of
length zero but does not prevent the method from reaching the optimum.

**(b)** In {numref}`fig-08-02-temporary-degeneracy` the corner $(1,0)$ lies
on three lines: $4x_1-x_2=4$, $4x_1+x_2=4$ and $x_2=0$. Two of them define
the corner, so one basic variable at $(1,0)$ is zero, and there are three
bases for the same point. The optimum $(\tfrac34,1)$, where $4x_1+3x_2=6$
and $4x_1+x_2=4$ meet, is nondegenerate.

```{figure} figures/08-02-temporary-degeneracy.svg
:name: fig-08-02-temporary-degeneracy
:width: 480px
:class: dark-light
:alt: Quadrilateral with corners (0,0), (1,0), (3/4,1) and (0,2); three lines meet at the degenerate corner (1,0)

Feasible region of {prf:ref}`ex-08-02-a3-q13`. Three lines meet at the
degenerate corner $(1,0)$, where the method spends one pivot of length
zero; the path is $(0,0)\to(1,0)\to(1,0)\to(\tfrac34,1)$, and the dashed
contour $3x_1+2x_2=\tfrac{17}{4}$ touches the region at the optimum.
```

The constraint $4x_1-x_2\le4$ is redundant for the region (it only touches
it at $(1,0)$), which is the typical cause of degeneracy in practice.
:::

:::{prf:example} Alternative optima that are not corner points
:label: ex-08-02-a3-q14

$$
\begin{aligned}
\max_{\vx}\ & z=2x_1-x_2+3x_3\\
\st\ & x_1-x_2+5x_3\le5,\\
& 2x_1-x_2+3x_3\le20,\\
& x_1,x_2,x_3\ge0 .
\end{aligned}
$$

(a) From the optimal tableau, show that all the alternative optimal
solutions other than the optimal BFS are not corner points.

(b) Give a two-dimensional example of a solution space and objective that
produce the same situation.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** With slack variables $x_4,x_5$:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_4$ | $1$ | $-1$ | $\boldsymbol{5}$ | $1$ | $0$ | $5$ |
| $x_5$ | $2$ | $-1$ | $3$ | $0$ | $1$ | $20$ |
| $-z$ | $2$ | $-1$ | $3$ | $0$ | $0$ | $0$ |

The pivots are: $x_3$ enters, $x_4$ leaves; $x_1$ enters, $x_3$ leaves;
$x_2$ enters, $x_5$ leaves. The optimal tableau is

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $-2$ | $-1$ | $1$ | $15$ |
| $x_2$ | $0$ | $1$ | $-7$ | $-2$ | $1$ | $10$ |
| $-z$ | $0$ | $0$ | $0$ | $0$ | $-1$ | $-20$ |

with $\vx^*=(15,10,0)$ and $z^*=20$. Two nonbasic variables, $x_3$ and the
slack $x_4$, have reduced cost $0$, so entering either keeps $z=20$. Their
columns, $(-2,-7)^\T$ and $(-1,-2)^\T$, have **no positive entry**: the
ratio test finds no leaving variable, so no other BFS can be reached, and
increasing them traces rays of optimal solutions. Written out,

$$
\begin{aligned}
x_1&=15+2t+u,\\
x_2&=10+7t+2u,\\
x_3&=t
\end{aligned}
$$

(with $x_4=u$, $x_5=0$) is optimal for all $t,u\ge0$. For example $t=1$,
$u=0$ gives $(17,17,1)$: $17-17+5=5$, $34-17+3=20$, and
$z=34-17+3=20$. The optimal set is an unbounded two-dimensional face whose
only corner is $(15,10,0)$; all other optimal solutions are not corner
points. Note that the feasible region is unbounded, but the objective is
not.

*(An earlier version of this solution gives the optimal ray as
$(15+t,10+t,t)$. That is not even feasible: at $t=1$ the first constraint
gives $16-11+5=10>5$, and $z=24$. The coefficients come from the $x_3$
column, $x_1=15+2t$ and $x_2=10+7t$.)*

**(b)** For example,

$$
\begin{aligned}
\max\ & z=-3x_1+x_2\\
\st\ & 2x_1-x_2\le2,\\
& -3x_1+x_2\le3,\\
& x_1,x_2\ge0 .
\end{aligned}
$$

The feasible region is unbounded, and the objective is parallel to the
boundary $-3x_1+x_2=3$. Every point of the ray $(t,3+3t)$, $t\ge0$, is
optimal with $z^*=3$, and only its end point $(0,3)$ is a corner
({numref}`fig-08-02-optimal-ray`).

```{figure} figures/08-02-optimal-ray.svg
:name: fig-08-02-optimal-ray
:width: 480px
:class: dark-light
:alt: Unbounded feasible region whose upper-left boundary is the ray from (0,3) in direction (1,3); every point of that ray is optimal

The two-dimensional example of {prf:ref}`ex-08-02-a3-q14`(b). The objective
$z=-3x_1+x_2$ is constant along the boundary ray from $(0,3)$ in direction
$(1,3)$ (thick), so every point of the ray is optimal, but only $(0,3)$ is
a corner point.
```
:::

:::{prf:example} Recognizing an unbounded objective in a tableau
:label: ex-08-02-unbounded-ray

Show that the LP

$$
\begin{aligned}
\max_{\vx}\ & z=3x_3-x_4+20\\
\st\ & x_1-3x_3+3x_4=6,\\
& x_2-8x_3+4x_4=4,\\
& x_1,x_2,x_3,x_4\ge0
\end{aligned}
$$

is unbounded, and give a ray of feasible solutions along which $z\to\infty$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The LP is in canonical form with respect to $\{x_1,x_2\}$; the constant
$20$ appears as $-20$ on the right of the objective row ($-z+3x_3-x_4=-20$):

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $-3$ | $3$ | $6$ |
| $x_2$ | $0$ | $1$ | $-8$ | $4$ | $4$ |
| $-z$ | $0$ | $0$ | $3$ | $-1$ | $-20$ |

For a maximization, $x_3$ is improving ($\bar c_3=3>0$), and its column
$(-3,-8)^\T$ has no positive entry. Raise $x_3$ to $t\ge0$ and keep $x_4=0$;
then

$$
\begin{aligned}
x_1&=6+3t,\\
x_2&=4+8t,\\
x_3&=t,\\
x_4&=0,\\
z&=20+3t,
\end{aligned}
$$

which is feasible for every $t\ge0$ while $z\to\infty$ as $t\to\infty$.
The LP is unbounded (an LP solver reports the same).
:::

:::{prf:example} An unbounded direction found by inspection
:label: ex-08-02-a3-q15

$$
\begin{aligned}
\max_{\vx}\ & z=20x_1+5x_2+x_3\\
\st\ & 3x_1+5x_2-5x_3\le50,\\
& x_1\le10,\\
& x_1+3x_2-4x_3\le20,\\
& x_1,x_2,x_3\ge0 .
\end{aligned}
$$

(a) By inspecting the constraints, determine in which direction ($x_1$,
$x_2$ or $x_3$) the solution space is unbounded.

(b) Without further computation, what can be said about the optimal
objective value?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** The coefficients of $x_3$ are $-5$, $0$ and $-4$, all
$\le0$: increasing $x_3$ never violates a constraint. So the solution space
is unbounded in the direction of $x_3$: from any feasible point, for example
$\vzero$, the ray $(0,0,t)$, $t\ge0$, stays feasible. It is not unbounded
in the direction of $x_1$ ($x_1\le10$) or of $x_2$ alone (with $x_3$ fixed,
$5x_2\le50+5x_3$).

**(b)** The objective coefficient of $x_3$ is positive, so along the ray
$z(0,0,t)=t\to\infty$: the optimal value is unbounded, and the LP has no
optimal solution.

**Check by the simplex method.** Dantzig's rule makes three pivots ($x_1$
for $x_5$, $x_2$ for $x_6$, $x_3$ for $x_4$) and reaches

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $0$ | $0$ | $1$ | $\tfrac{3}{5}$ | $-\tfrac{4}{5}$ | $-1$ | $2$ |
| $x_1$ | $1$ | $0$ | $0$ | $0$ | $1$ | $0$ | $10$ |
| $x_2$ | $0$ | $1$ | $0$ | $\tfrac{4}{5}$ | $-\tfrac{7}{5}$ | $-1$ | $6$ |
| $-z$ | $0$ | $0$ | $0$ | $-\tfrac{23}{5}$ | $-\tfrac{61}{5}$ | $6$ | $-232$ |

where the improving variable $x_6$ ($\bar c_6=6$) has the column
$(-1,0,-1)^\T$ without positive entries: the tableau test for
unboundedness ({prf:ref}`alg-03-03-simplex`, Step 3) confirms (b).
:::

:::{prf:example} A degenerate optimum, and infeasibility
:label: ex-08-02-a3-q16

$$
\begin{aligned}
\max_{\vx}\ & z=3x_1+2x_2+3x_3\\
\st\ & 2x_1+x_2+x_3\le4,\\
& 3x_1+4x_2+2x_3\ge16,\\
& x_1,x_2,x_3\ge0 .
\end{aligned}
$$

Show by hand computation that the optimal solution can be represented by a
basis in which a slack or surplus variable is basic at zero level, and find
an optimal basis. Does the problem have an optimal solution?
What changes if the second right-hand side is $17$ instead of $16$?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Use the two-phase method with slack $x_4$, surplus $x_5$ and artificial $R$
in the second row; Phase I minimizes $\phi=R$.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $R$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_4$ | $2$ | $1$ | $1$ | $1$ | $0$ | $0$ | $4$ |
| $R$ | $3$ | $\boldsymbol{4}$ | $2$ | $0$ | $-1$ | $1$ | $16$ |
| $-z$ | $3$ | $2$ | $3$ | $0$ | $0$ | $0$ | $0$ |
| $-\phi$ | $-3$ | $-4$ | $-2$ | $0$ | $1$ | $0$ | $-16$ |

$x_2$ enters ($-4$). The ratios $4/1$ and $16/4$ **tie**.

*Tie broken in favor of $R$ leaving* (bold pivot above):

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $R$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_4$ | $\tfrac{5}{4}$ | $0$ | $\tfrac{1}{2}$ | $1$ | $\tfrac{1}{4}$ | $-\tfrac{1}{4}$ | $0$ |
| $x_2$ | $\tfrac{3}{4}$ | $1$ | $\tfrac{1}{2}$ | $0$ | $-\tfrac{1}{4}$ | $\tfrac{1}{4}$ | $4$ |
| $-z$ | $\tfrac{3}{2}$ | $0$ | $2$ | $0$ | $\tfrac{1}{2}$ | $-\tfrac{1}{2}$ | $-8$ |
| $-\phi$ | $0$ | $0$ | $0$ | $0$ | $0$ | $1$ | $0$ |

$\phi^*=0$ with no artificial variable basic, so the LP is feasible, and
the slack $x_4$ is basic **at zero level**: the BFS $\vx=(0,4,0)$ is
degenerate, with both constraints binding. In Phase II the $-z$ row still
has positive entries, so this basis does not pass the optimality test;
$x_3$ enters ($2$), and the only positive entry with a ratio is in row
$x_4$: ratio $0/\tfrac12=0$, a degenerate pivot.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $\tfrac{5}{2}$ | $0$ | $1$ | $2$ | $\tfrac{1}{2}$ | $0$ |
| $x_2$ | $-\tfrac{1}{2}$ | $1$ | $0$ | $-1$ | $-\tfrac{1}{2}$ | $4$ |
| $-z$ | $-\tfrac{7}{2}$ | $0$ | $0$ | $-4$ | $-\tfrac{1}{2}$ | $-8$ |

All reduced costs are now $\le0$. The optimal solution is
$\vx^*=(0,4,0)$ with $z^*=8$ (an LP solver agrees), and the optimal basis
$\{x_3,x_2\}$ has $x_3$ basic at zero level.

*Tie broken in favor of $x_4$ leaving:*

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $R$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_2$ | $2$ | $1$ | $1$ | $1$ | $0$ | $0$ | $4$ |
| $R$ | $-5$ | $0$ | $-2$ | $-4$ | $-1$ | $1$ | $0$ |
| $-z$ | $-1$ | $0$ | $1$ | $-2$ | $0$ | $0$ | $-8$ |
| $-\phi$ | $5$ | $0$ | $2$ | $4$ | $1$ | $0$ | $0$ |

Now Phase I ends with $\phi^*=0$ but $R$ basic at zero level. The $-\phi$
row says $R=5x_1+2x_3+4x_4+x_5$, so, as in {prf:ref}`ex-08-02-a3-q10`, the
variables $x_1,x_3,x_4,x_5$ must be zero in every feasible solution; the
feasible set is the single point $(0,4,0)$, and $z^*=8$ again.

**Answer.** Yes, the problem has an optimal solution, $\vx^*=(0,4,0)$ with
$z^*=8$, and it is degenerate: only one variable ($x_2$) is positive,
although there are two constraints. The point is described by four bases,
$\{x_2,x_1\}$, $\{x_2,x_3\}$, $\{x_2,x_4\}$ and $\{x_2,x_5\}$; in the last
two a slack ($x_4$) or a surplus ($x_5$) variable is basic at zero level.
Only $\{x_2,x_3\}$ passes the optimality test; the bases with the slack or
the surplus variable describe the same optimal point, but their tableaux
still show a positive reduced cost ($\bar c_3=2$ and $\bar c_3=1$). This
illustrates that the reduced-cost test is sufficient but not necessary at a
degenerate BFS ({prf:ref}`thm-03-03-optimality`).
*(An earlier version of this solution presents the tableau of the basis
$\{x_5,x_2\}$, with objective row $(-1,0,1,-2,0)$, as a second optimal
tableau. For a maximization the entry $+1$ under $x_3$ fails the
optimality test; the tableau represents the optimal point but not an
optimal basis.)*

**Right-hand side 17.** With $3x_1+4x_2+2x_3\ge17$ the same first pivot
($x_2$ enters, now $x_4$ leaves with ratio $4<\tfrac{17}{4}$) gives

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $R$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_2$ | $2$ | $1$ | $1$ | $1$ | $0$ | $0$ | $4$ |
| $R$ | $-5$ | $0$ | $-2$ | $-4$ | $-1$ | $1$ | $1$ |
| $-z$ | $-1$ | $0$ | $1$ | $-2$ | $0$ | $0$ | $-8$ |
| $-\phi$ | $5$ | $0$ | $2$ | $4$ | $1$ | $0$ | $-1$ |

All $-\phi$ entries are $\ge0$, so Phase I is optimal with $\phi^*=1>0$: the
LP is **infeasible** ({prf:ref}`thm-03-04-phase-one`). Indeed, the first
constraint gives $x_2\le4-2x_1-x_3$, so
$3x_1+4x_2+2x_3\le16-5x_1-2x_3\le16<17$.
:::

## References

- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Sections 3.2–3.5: simplex iterations,
  degeneracy and the initial basic feasible solution.)
- Charnes, A. (1952). Optimality and degeneracy in linear programming.
  *Econometrica, 20*(2), 160–170. <https://doi.org/10.2307/1907845>
  (The big-M method.)
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Chapter 5: the simplex method and its entering rule;
  Chapter 10: degeneracy.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Problem sets of Sections 3.1–3.5, from which several problems
  of this chapter are adapted with restated wording; the cut-cube and
  degenerate-pentagon figures are redrawn from data chosen to match them.)
