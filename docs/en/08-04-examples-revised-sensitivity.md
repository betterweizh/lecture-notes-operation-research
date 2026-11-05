# Revised and Dual Simplex, Sensitivity

This chapter works through the examples and problems of
{doc}`03-06-revised-dual-simplex` and of {doc}`03-07-sensitivity-analysis`. The
conventions are those of the two theory chapters: tableaux list the basic
variables on the left and the objective row last, labeled $-z$, holding the
reduced costs $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ and, in the
right-hand-side (RHS) column, $-z$. A minimization is optimal when all
$\bar c_j\ge0$, a maximization when all $\bar c_j\le0$. Every number was recomputed in exact rational arithmetic, and
every optimum, including the optimum after each change of data, was confirmed
by re-solving the modified linear program (LP) with the HiGHS solver.

(sec-08-04-revised)=
## Revised simplex and tableau conditions

:::{prf:example} Revised simplex from a given basis
:label: ex-08-04-a5-q3-revised

Solve the LP

$$
\begin{aligned}
\min_{\vx}\ & 7x_2+11x_3-10x_4+26x_6 \\
\st\ & x_2-x_3+x_5+x_6=3,\\
& x_2-x_3+x_4+3x_6=4,\\
& x_1+x_2-3x_3+x_4+x_5=6,\\
& x_1,\dots,x_6\ge0
\end{aligned}
$$

by the revised simplex method, starting from the basis $(x_2,x_4,x_5)$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Iteration 1.** The basis columns are $\mA_{\cdot2}=(1,1,1)^\T$,
$\mA_{\cdot4}=(0,1,1)^\T$ and $\mA_{\cdot5}=(1,0,1)^\T$, so

$$
\begin{gathered}
\mB_1=\begin{pmatrix}1&0&1\\1&1&0\\1&1&1\end{pmatrix},\\
\mB_1^{-1}=\begin{pmatrix}1&1&-1\\-1&0&1\\0&-1&1\end{pmatrix}.
\end{gathered}
$$

Then $\vx_B=\mB_1^{-1}\vb=(1,3,2)^\T$ for $(x_2,x_4,x_5)$, which is
feasible, with $z=7-30=-23$. Pricing with $\vc_B=(7,-10,0)$:

$$
\begin{gathered}
\vpi^\T=\vc_B^\T\mB_1^{-1}=(17,\ 7,\ -17),\\
(\bar c_1,\bar c_3,\bar c_6)=(17,\ -16,\ -12).
\end{gathered}
$$

For example $\bar c_3=11-\vpi^\T(-1,-1,-3)^\T=11-(-17-7+51)=-16$. The most
negative is $\bar c_3$, so $x_3$ enters, with
$\bar\mA_{\cdot3}=\mB_1^{-1}(-1,-1,-3)^\T=(1,-2,-2)^\T$. Only the first entry
is positive: the ratio is $1/1=1$ and $x_2$ leaves. The eta matrix and the
new inverse are

$$
\begin{gathered}
\mE_1=\begin{pmatrix}1&0&0\\2&1&0\\2&0&1\end{pmatrix},\\
\mB_2^{-1}=\mE_1\mB_1^{-1}=\begin{pmatrix}1&1&-1\\1&2&-1\\2&1&-1\end{pmatrix}.
\end{gathered}
$$

**Iteration 2.** Basis $(x_3,x_4,x_5)$, $\vx_B=\mB_2^{-1}\vb=(1,5,4)^\T$ and
$z=11-50=-39$. Pricing with $\vc_B=(11,-10,0)$:

$$
\begin{gathered}
\vpi^\T=(1,\ -9,\ -1),\\
(\bar c_1,\bar c_2,\bar c_6)=(1,\ 16,\ 52)\ge0 .
\end{gathered}
$$

The basis is optimal:
$\vx^*=(0,0,1,5,4,0)$ and $z^*=-39$.

**Check.** The three constraints give $-1+4=3$, $-1+5=4$ and $-3+5+4=6$,
and $\vpi^\T\vb=3-36-6=-39=z^*$. HiGHS returns the same solution.

*(An earlier version of this solution prints $\vx_B=(1,4,5)$ and
$\vx=(0,0,1,4,5,0)$, with the values of $x_4$ and $x_5$ exchanged; that point
violates the second constraint and has objective $-29$. The reduced cost
$\bar c_6=52$ is correct.)*
:::

:::{prf:example} Scaling a column and the right-hand side in the ratio test
:label: ex-08-04-a5-q1-scaling

Consider $\max\{\vc^\T\vx:\mA\vx\le\vb,\ \vx\ge\vzero\}$ with
$\vb\ge\vzero$, solved by the simplex method. At some iteration $x_j$ is the
entering variable and at least one entry of $\mB^{-1}\mA_{\cdot j}$ is
positive.

(a) The column $\mA_{\cdot j}$ is replaced by $\alpha\mA_{\cdot j}$ with
$\alpha>0$, and $x_j$ is still chosen to enter. How does the value that $x_j$
takes after the pivot change?

(b) Answer (a) when, in addition, $\vb$ is replaced by $\beta\vb$ with
$\beta>0$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The entering variable takes the value of the minimum ratio,

$$
\begin{gathered}
x_j=\theta=\min\left\{\frac{\bar b_i}{\bar a_{ij}}:\ \bar a_{ij}>0\right\},\\
\bar\vb=\mB^{-1}\vb,\\
\bar\mA_{\cdot j}=\mB^{-1}\mA_{\cdot j}.
\end{gathered}
$$

Since $x_j$ is nonbasic, $\mB$ does not contain its column and is not
affected by the scaling.

**(a)** The updated column becomes $\mB^{-1}(\alpha\mA_{\cdot j})=\alpha\bar\mA_{\cdot j}$.
Its positive entries are in the same rows, every ratio is divided by
$\alpha$, and the same row attains the minimum. Hence the new value is
$x_j'=\theta/\alpha$: **$1/\alpha$ times the old value**. (Scaling a column
by $\alpha$ is the same as measuring $x_j$ in units of $1/\alpha$.)

**(b)** Now also $\bar\vb'=\mB^{-1}(\beta\vb)=\beta\bar\vb$, so every ratio
is multiplied by $\beta/\alpha$ and again the same row wins:
$x_j'=(\beta/\alpha)\,\theta$, **$\beta/\alpha$ times the old value**.

**Check.** In $\max\ 3x_1+2x_2$ with $x_1+x_2\le4$, $x_1+3x_2\le6$, the
slack basis lets $x_1$ enter with $\theta=\min\{4,6\}=4$. With $\alpha=2$
the ratios are $4/2$ and $6/2$, so $x_1=2=4/2$; with $\alpha=2$ and
$\beta=3$ they are $12/2$ and $18/2$, so $x_1=6=\tfrac32\cdot4$.

The proviso "$x_j$ is still chosen to enter" matters: the scaling also
changes the reduced cost to $c_j-\alpha\vpi^\T\mA_{\cdot j}$, which may
change the choice of the entering variable.
:::

:::{prf:example} Conditions on a tableau with parameters
:label: ex-08-04-pq2-conditions

The following tableau arose while maximizing $z$ over nonnegative
variables; $p_1,\dots,p_6$ are unknown parameters.

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|---|---|---|---|---|---|---|---|
| $x_3$ | 0 | 2 | 1 | $p_2$ | 0 | 2 | $p_1$ |
| $x_1$ | 1 | $p_3$ | 0 | 2 | 0 | −1 | 1 |
| $x_5$ | 0 | 0 | 0 | $p_4$ | 1 | −4 | 3 |
| $-z$ | 0 | $p_5$ | 0 | $p_6$ | 0 | −2 | −17 |

State and justify conditions on the parameters under which:

(a) $x_4$ can enter the basis, but the objective value does not change;

(b) the LP has no feasible solution;

(c) the LP is unbounded.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The tableau is in canonical form, so it is equivalent to the LP: express
the basic variables through the nonbasic ones,

$$
\begin{aligned}
x_3&=p_1-2x_2-p_2x_4-2x_6,\\
x_1&=1-p_3x_2-2x_4+x_6,\\
x_5&=3-p_4x_4+4x_6,\\
z&=17+p_5x_2+p_6x_4-2x_6 .
\end{aligned}
$$

**(a)** In a maximization $x_4$ is eligible to enter when its reduced cost is
positive, $p_6>0$ (and it is the most-positive choice if also $p_6\ge p_5$).
The objective changes by $\bar c_4\,\theta$, so it stays the same exactly when
the step $\theta$ is zero: the ratio test must have a row with right-hand side
0 and a positive entry in column $x_4$. Row $x_1$ gives the ratio $1/2>0$ and
row $x_5$ gives $3/p_4>0$ (if $p_4>0$), so the zero ratio must come from row
$x_3$: $p_1=0$ and $p_2>0$. The conditions are $p_6>0$, $p_1=0$ and
$p_2>0$: a degenerate pivot.

A second case gives the same result in a different way. If $p_6=0$ and the
tableau is optimal ($p_1\ge0$, $p_5\le0$), then $x_4$ can be brought into
the basis (the entry 2 in row $x_1$ always gives a finite ratio) with
$\bar c_4\,\theta=0$ for any step: the pivot leads to an alternative optimal
basis with the same $z=17$.

**(b)** If $p_1<0$ the basic solution is infeasible. Row $x_3$ then states
$2x_2+x_3+p_2x_4+2x_6=p_1<0$. When $p_2\ge0$ every coefficient on the left is
nonnegative, so no $\vx\ge\vzero$ satisfies it (this is the infeasibility
test of the dual simplex method). Hence the conditions $p_1<0$ and
$p_2\ge0$ guarantee that the LP has no feasible solution.

**(c)** The usual test, a nonbasic column with positive reduced cost and no
positive entry, cannot be met here: every nonbasic column has the entry
$+2$ (column $x_2$ in row $x_3$, column $x_4$ in row $x_1$, column $x_6$ in row
$x_3$). Unboundedness must therefore be found from a ray. A direction
$\vd\ge\vzero$ in the nonbasic variables $(x_2,x_4,x_6)=(d_2,d_4,d_6)$ keeps
the basic variables nonnegative for every step length only if
$2d_2+p_2d_4+2d_6\le0$ (row $x_3$) and $p_3d_2+2d_4-d_6\le0$ (row $x_1$).
The first forces $d_4>0$ and $p_2\le0$; with $p_2=0$ it would also force
$d_2=d_6=0$, and then row $x_1$ would require $2d_4\le0$. So
**$p_2<0$ is necessary.** Taking $d_2=0$, $d_4=1$, $d_6=g\ge0$, the three
rows and the objective require

$$
\begin{aligned}
g&\le-\tfrac{p_2}{2} && \text{(row } x_3\text{)}\\
g&\ge2 && \text{(row } x_1\text{)}\\
g&\ge\tfrac{p_4}{4} && \text{(row } x_5\text{)}\\
p_6-2g&>0 && \text{(objective)}
\end{aligned}
$$

With $g_0=\max\{2,\,p_4/4\}$, a sufficient set of conditions is $p_1\ge0$
(for feasibility), $p_2\le-2g_0$ and $p_6>2g_0$. For example
$p=(1,-5,0,0,0,5)$: the simplex method lets $x_4$ enter (ratio $1/2$ in row $x_1$) and reaches

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|---|---|---|---|---|---|---|---|
| $x_3$ | 5/2 | 2 | 1 | 0 | 0 | −1/2 | 7/2 |
| $x_4$ | 1/2 | 0 | 0 | 1 | 0 | −1/2 | 1/2 |
| $x_5$ | 0 | 0 | 0 | 0 | 1 | −4 | 3 |
| $-z$ | −5/2 | 0 | 0 | 0 | 0 | 1/2 | −39/2 |

where $x_6$ has reduced cost $1/2>0$ and a column without positive entries:
the LP is unbounded, as HiGHS confirms. With $p_2=-3$ or $p_6=3$ instead
(the other values unchanged) the LP has a finite optimum, as HiGHS also
confirms.

*(As posed, part (c) suggests a condition readable directly from the given
tableau; with these entries there is none, and the answer requires the ray
argument above.)*
:::

(sec-08-04-dual-simplex)=
## Dual simplex and adding constraints

:::{prf:example} The dual simplex method detects infeasibility
:label: ex-08-04-dual-infeasible

Apply the dual simplex method to

$$
\begin{aligned}
\min_{\vx}\ & 8x_1+8x_2+9x_3 \\
\st\ & x_1+x_2+x_3\le1,\\
& 2x_1+4x_2+x_3\ge8,\\
& x_1-x_2-x_3\ge2,\\
& x_1,x_2,x_3\ge0 .
\end{aligned}
$$
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Multiply the two $\ge$ rows by $-1$ and add slack variables $x_4,x_5,x_6$:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|---|---|---|---|---|---|---|---|
| $x_4$ | 1 | 1 | 1 | 1 | 0 | 0 | 1 |
| $x_5$ | −2 | −4 | −1 | 0 | 1 | 0 | −8 |
| $x_6$ | −1 | 1 | 1 | 0 | 0 | 1 | −2 |
| $-z$ | 8 | 8 | 9 | 0 | 0 | 0 | 0 |

The slack basis is dual feasible (all reduced costs $8,8,9\ge0$) but
$\mB^{-1}\vb=(1,-8,-2)^\T\not\ge\vzero$.

**Pivot 1.** The most negative value is $x_5=-8$: $x_5$ leaves. Its row has
negative entries in $x_1,x_2,x_3$ with dual ratios $8/2=4$, $8/4=2$ and
$9/1=9$, so $x_2$ enters:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|---|---|---|---|---|---|---|---|
| $x_4$ | 1/2 | 0 | 3/4 | 1 | 1/4 | 0 | −1 |
| $x_2$ | 1/2 | 1 | 1/4 | 0 | −1/4 | 0 | 2 |
| $x_6$ | −3/2 | 0 | 3/4 | 0 | 1/4 | 1 | −4 |
| $-z$ | 4 | 0 | 7 | 0 | 2 | 0 | −16 |

**Pivot 2.** Now $x_6=-4$ leaves. Its only negative entry is $-3/2$ in column
$x_1$ (dual ratio $4/(3/2)=8/3$), so $x_1$ enters:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|---|---|---|---|---|---|---|---|
| $x_4$ | 0 | 0 | 1 | 1 | 1/3 | 1/3 | −7/3 |
| $x_2$ | 0 | 1 | 1/2 | 0 | −1/6 | 1/3 | 2/3 |
| $x_1$ | 1 | 0 | −1/2 | 0 | −1/6 | −2/3 | 8/3 |
| $-z$ | 0 | 0 | 9 | 0 | 8/3 | 8/3 | −80/3 |

**Pivot 3.** Now $x_4=-7/3$ must leave, but its row $(0,0,1,1,1/3,1/3)$ has
no negative entry. The row states $x_3+x_4+\tfrac13x_5+\tfrac13x_6=-\tfrac73$,
impossible for nonnegative variables: **the LP has no feasible solution.**

**Certificate.** The row of $x_4$ is $\vy^\T$ times the constraint rows, with
$\vy^\T=(1,\tfrac13,\tfrac13)$, the row of $\mB^{-1}$ (columns $x_4,x_5,x_6$).
Adding the $\le$ forms of the constraints with these weights,

$$
\begin{aligned}
&(x_1+x_2+x_3)+\tfrac13(-2x_1-4x_2-x_3)\\
&\quad+\tfrac13(-x_1+x_2+x_3)\le1-\tfrac83-\tfrac23,
\end{aligned}
$$

gives $x_3\le-\tfrac73$, which contradicts $x_3\ge0$. HiGHS also reports
infeasibility. Along the way $z$ rose from 0 to 16 to $80/3$: these are values
of the dual objective at dual feasible points, and the dual LP is unbounded,
as it must be when the primal is infeasible and the dual is feasible.

*(An earlier version of the last tableau prints the column of $x_3$ as
$(7/8,3/8,-1/4,8)$ and the entry of $x_5$ in row $x_1$ as $1/6$; the correct
values are $(1,1/2,-1/2,9)$ and $-1/6$. The conclusion is unaffected.)*
:::

:::{prf:example} The revised dual simplex method
:label: ex-08-04-a5-q4-revised-dual

The revised dual simplex method works with $\mB^{-1}$ only. Starting from the
slack basis $\mB=\mI$, each iteration (1) computes $\vx_B=\mB^{-1}\vb$ and
removes the most negative basic variable $x_r$, stopping if
$\vx_B\ge\vzero$; (2) computes the reduced costs
$\bar c_j=c_j-\vc_B^\T\mB^{-1}\mA_{\cdot j}$ and the entries
$\bar a_{rj}=(\mB^{-1}\mA_{\cdot j})_r$ of the leaving row for the nonbasic
$x_j$, and lets enter the variable that attains
$\min\{\abs{\bar c_j/\bar a_{rj}}:\bar a_{rj}<0\}$, stopping with "infeasible"
if no $\bar a_{rj}$ is negative; and (3) updates $\mB^{-1}$. Apply it to

$$
\begin{aligned}
\min_{\vx}\ & 3x_1+2x_2 \\
\st\ & 3x_1+x_2\ge3,\\
& 4x_1+3x_2\ge6,\\
& x_1+2x_2\le3,\\
& x_1,x_2\ge0 .
\end{aligned}
$$
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

In standard form, with the $\ge$ rows multiplied by $-1$:

$$
\begin{aligned}
-3x_1-x_2+x_3&=-3,\\
-4x_1-3x_2+x_4&=-6,\\
x_1+2x_2+x_5&=3 .
\end{aligned}
$$

**Iteration 1.** Basis $(x_3,x_4,x_5)$, $\mB^{-1}=\mI$,
$\vx_B=(-3,-6,3)^\T$. The most negative is $x_4$ (row 2). With $\vpi=\vzero$
the reduced costs are $\bar c_1=3$, $\bar c_2=2$, both nonnegative (dual
feasible). Row 2 of $\mB^{-1}$ is $(0,1,0)$, so $\bar a_{21}=-4$ and
$\bar a_{22}=-3$; the ratios are $3/4$ and $2/3$, and $x_2$ enters. Its
column is $\bar\mA_{\cdot2}=(-1,-3,2)^\T$, the pivot is $-3$, the eta column
is $(-1/3,\,-1/3,\,2/3)^\T$ in position 2, and

$$
\mB_2^{-1}=\begin{pmatrix}1&-1/3&0\\0&-1/3&0\\0&2/3&1\end{pmatrix}.
$$

**Iteration 2.** Basis $(x_3,x_2,x_5)$, $\vx_B=\mB_2^{-1}\vb=(-1,2,-1)^\T$.
Two values tie at $-1$; take $x_5$ (row 3). The multipliers are
$\vpi^\T=(0,2,0)\mB_2^{-1}=(0,-2/3,0)$, so $\bar c_1=3-\tfrac83=\tfrac13$ and
$\bar c_4=0+\tfrac23=\tfrac23$. Row 3 of $\mB_2^{-1}$ is $(0,2/3,1)$, giving
$\bar a_{31}=-\tfrac83+1=-\tfrac53$ and $\bar a_{34}=\tfrac23$. Only $x_1$
qualifies (ratio $\tfrac13/\tfrac53=\tfrac15$), so $x_1$ enters with
$\bar\mA_{\cdot1}=(-5/3,\,4/3,\,-5/3)^\T$. The eta column for the pivot
$-5/3$ in position 3 is $(-1,\,4/5,\,-3/5)^\T$, and

$$
\mB_3^{-1}=\begin{pmatrix}1&-1&-1\\0&1/5&4/5\\0&-2/5&-3/5\end{pmatrix}.
$$

**Iteration 3.** Basis $(x_3,x_2,x_1)$, $\vx_B=\mB_3^{-1}\vb=(0,\,6/5,\,3/5)^\T\ge\vzero$:
the basis is optimal. As a check, $\vpi^\T=(0,2,3)\mB_3^{-1}=(0,-4/5,-1/5)$
gives $\bar c_4=4/5$ and $\bar c_5=1/5$, still nonnegative. The optimum is

$$
\begin{aligned}
x_1^*&=\tfrac35,\\
x_2^*&=\tfrac65,\\
z^*&=\tfrac95+\tfrac{12}{5}=\tfrac{21}{5}=4.2 .
\end{aligned}
$$

It satisfies all three constraints with equality ($3x_1+x_2=3$,
$4x_1+3x_2=6$, $x_1+2x_2=3$), which is why the basic variable $x_3$ is 0
(a degenerate optimum); HiGHS returns the same point. Taking $x_3$ instead of $x_5$
at the tie in iteration 2 also brings in $x_1$ (ratios $1/5$ and $2$) and
ends at the same point, with basis $(x_1,x_2,x_5)$ and $x_5=0$.

*(An earlier version of this solution prints the reduced costs as
$\bar c_1=-3$, $\bar c_2=-2$ in iteration 1 and $\bar c_1=-1/3$,
$\bar c_4=-2/3$ in iteration 2. They are positive, as dual feasibility
requires; the ratios are unaffected because they use absolute values.)*
:::

:::{prf:example} Dual simplex with an artificial constraint
:label: ex-08-04-a5-q5-artificial-constraint

Consider

$$
\begin{aligned}
\max_{\vx}\ & 2x_1-x_2+x_3 \\
\st\ & 2x_1+3x_2-5x_3\ge4,\\
& -x_1+9x_2-x_3\ge3,\\
& 4x_1+6x_2+3x_3\le8,\\
& x_1,x_2,x_3\ge0 .
\end{aligned}
$$

The basis of the surplus variables $x_4,x_5$ and the slack $x_6$ has
$x_4=-4$ and $x_5=-3$, so it is infeasible, and its reduced costs violate the
optimality condition of a maximization, so the dual simplex method cannot
start from it either. Show that adding the constraint $x_1+x_3\le M$, with
$M$ large enough to cut off no relevant point, and pivoting on its row in the
column with the most positive reduced cost produces a dual feasible tableau.
Then finish with the dual simplex method.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Write the $\ge$ rows as $\le$ rows with slack variables (so $x_4,x_5$ are the
surplus variables of the original rows) and add the artificial constraint
with slack $x_7$ ({ref}`sec-03-06-artificial-constraint`):

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | RHS |
|---|---|---|---|---|---|---|---|---|
| $x_4$ | −2 | −3 | 5 | 1 | 0 | 0 | 0 | −4 |
| $x_5$ | 1 | −9 | 1 | 0 | 1 | 0 | 0 | −3 |
| $x_6$ | 4 | 6 | 3 | 0 | 0 | 1 | 0 | 8 |
| $x_7$ | 1 | 0 | 1 | 0 | 0 | 0 | 1 | $M$ |
| $-z$ | 2 | −1 | 1 | 0 | 0 | 0 | 0 | 0 |

**Making the tableau dual feasible.** The reduced costs $2$ and $1$ of $x_1$
and $x_3$ are positive. Pivot on row $x_7$ in column $x_1$, the most positive:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | RHS |
|---|---|---|---|---|---|---|---|---|
| $x_4$ | 0 | −3 | 7 | 1 | 0 | 0 | 2 | $2M-4$ |
| $x_5$ | 0 | −9 | 0 | 0 | 1 | 0 | −1 | $-3-M$ |
| $x_6$ | 0 | 6 | −1 | 0 | 0 | 1 | −4 | $8-4M$ |
| $x_1$ | 1 | 0 | 1 | 0 | 0 | 0 | 1 | $M$ |
| $-z$ | 0 | −1 | −1 | 0 | 0 | 0 | −2 | $-2M$ |

All reduced costs $(0,-1,-1,0,0,0,-2)$ are now $\le0$: the tableau is dual
feasible, as {ref}`sec-03-06-artificial-constraint` guarantees (subtracting the
largest positive reduced cost from all of them).

**Dual simplex, pivot 1.** For large $M$ the most negative right-hand side is
$8-4M$ (it is below $-3-M$ once $M>11/3$): $x_6$ leaves. Its row has negative
entries in $x_3$ (ratio $\abs{-1/-1}=1$) and $x_7$ (ratio
$\abs{-2/-4}=1/2$), so $x_7$ enters:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | RHS |
|---|---|---|---|---|---|---|---|---|
| $x_4$ | 0 | 0 | 13/2 | 1 | 0 | 1/2 | 0 | 0 |
| $x_5$ | 0 | −21/2 | 1/4 | 0 | 1 | −1/4 | 0 | −5 |
| $x_7$ | 0 | −3/2 | 1/4 | 0 | 0 | −1/4 | 1 | $M-2$ |
| $x_1$ | 1 | 3/2 | 3/4 | 0 | 0 | 1/4 | 0 | 2 |
| $-z$ | 0 | −4 | −1/2 | 0 | 0 | −1/2 | 0 | −4 |

**Pivot 2.** $x_5=-5$ leaves. Negative entries: $x_2$ (ratio
$4/\tfrac{21}{2}=\tfrac{8}{21}$) and $x_6$ (ratio
$\tfrac12/\tfrac14=2$), so $x_2$ enters:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | RHS |
|---|---|---|---|---|---|---|---|---|
| $x_4$ | 0 | 0 | 13/2 | 1 | 0 | 1/2 | 0 | 0 |
| $x_2$ | 0 | 1 | −1/42 | 0 | −2/21 | 1/42 | 0 | 10/21 |
| $x_7$ | 0 | 0 | 3/14 | 0 | −1/7 | −3/14 | 1 | $M-9/7$ |
| $x_1$ | 1 | 0 | 11/14 | 0 | 1/7 | 3/14 | 0 | 9/7 |
| $-z$ | 0 | 0 | −25/42 | 0 | −8/21 | −17/42 | 0 | −44/21 |

All right-hand sides are nonnegative (for $M\ge9/7$) and all reduced costs
are $\le0$: the tableau is optimal. The slack $x_7=M-9/7$ of the artificial
constraint is basic, so that constraint is not binding and the solution is
optimal for the original LP:

$$
\begin{aligned}
x_1^*&=\tfrac97,\\
x_2^*&=\tfrac{10}{21},\\
x_3^*&=0,\\
z^*&=\tfrac{18}{7}-\tfrac{10}{21}=\tfrac{44}{21}\approx2.095 .
\end{aligned}
$$

**Check.** $2\cdot\tfrac97+3\cdot\tfrac{10}{21}=4$,
$-\tfrac97+9\cdot\tfrac{10}{21}=3$ and
$4\cdot\tfrac97+6\cdot\tfrac{10}{21}=8$: all three constraints are binding
(the basic $x_4=0$ shows the degeneracy). HiGHS gives the same optimum.

If instead $x_5$ (value $-3-M$) leaves first, as the rule "any negative
basic variable" allows, the pivots are $x_2$ in, $x_5$ out and then $x_7$ in,
$x_6$ out, ending in the same final tableau; the intermediate objective
entries are $-2M$, $\tfrac13-\tfrac{17}{9}M$ and $-\tfrac{44}{21}$.

*(An earlier version of this solution follows the second path but carries
the entry 0 instead of $-1$ in row $x_6$, column $x_3$ from the second tableau
on, so its final column of $x_3$ reads $(7,0,0,1,-1)$; the correct column is
$(13/2,\,-1/42,\,3/14,\,11/14)$ with reduced cost $-25/42$. The optimal
solution is not affected.)*
:::

(sec-08-04-cost-ranging)=
## Cost ranging

:::{prf:example} Price changes for an electric-motor plant
:label: ex-08-04-a6-q2-electra

A plant makes four types of electric motor, each on its own assembly line
with daily capacities of 500, 500, 800 and 750 motors. Each motor needs an
electronic component: 8, 5, 4 and 6 units for types 1 to 4, and the supplier
delivers at most 8000 units a day. The motors sell for \$60, \$40, \$25 and
\$30.

(a) Find the revenue-maximizing daily production.

(b) Competition may force a lower price for type 2. By how much can it drop
without changing the production plan?

(c) All prices are cut by 25 percent. Does the plan change?

(d) Type 4 is not produced. By how much must its price rise before it is?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** With $x_i$ the daily number of type-$i$ motors:

$$
\begin{aligned}
\max_{\vx}\ & 60x_1+40x_2+25x_3+30x_4 \\
\st\ & 8x_1+5x_2+4x_3+6x_4+s_1=8000,\\
& x_1+s_2=500,\\
& x_2+s_3=500,\\
& x_3+s_4=800,\\
& x_4+s_5=750,\\
& \vx,\vs\ge\vzero .
\end{aligned}
$$

The optimal basis is $(x_1,x_2,x_3,s_4,s_5)$ with

$$
\begin{gathered}
\mB^{-1}=\begin{pmatrix}0&1&0&0&0\\0&0&1&0&0\\1/4&-2&-5/4&0&0\\-1/4&2&5/4&1&0\\0&0&0&0&1\end{pmatrix},\\
\vx_B=\mB^{-1}\vb=(500,\,500,\,375,\,425,\,750)^\T .
\end{gathered}
$$

The plan is 500 motors of type 1, 500 of type 2, 375 of type 3 and none of
type 4, with revenue \$59 375 a day. All 8000 components are used
($s_1=0$); line 3 has 425 motors of spare capacity and line 4 is idle
($s_4=425$, $s_5=750$). The multipliers are

$$
\vpi^\T=\vc_B^\T\mB^{-1}=\left(\tfrac{25}{4},\ 10,\ \tfrac{35}{4},\ 0,\ 0\right),
$$

and the reduced costs of the nonbasic $x_4,s_1,s_2,s_3$ are
$30-6\cdot\tfrac{25}{4}=-7.5$, $-6.25$, $-10$ and $-8.75$, all $\le0$, so the
plan is optimal. It is the greedy plan by revenue per component: types 2
(\$8.00 per unit of component) and 1 (\$7.50) at capacity, type 3 (\$6.25)
with the remaining 1500 units, type 4 (\$5.00) not at all. The component is
worth $\pi_1=\$6.25$ per unit, the revenue per unit of the marginal type 3.

**(b)** Let the type-2 price be $40+\delta$. Since $x_2$ is basic, $\vpi$
changes, but of the nonbasic reduced costs only that of $s_3$ (the slack of
$x_2\le500$) involves $\delta$:
$\bar c_{s_3}=-\tfrac{35}{4}-\delta\le0$ gives $\delta\ge-8.75$. The price
of type 2 can fall by up to **\$8.75**, to \$31.25, without changing the
plan. At \$31.25 type 2 earns \$6.25 per component, the same as type 3, and
the plan $(500,160,800,0)$ ties with the current one at \$55 000; below that
price HiGHS switches to it.

**(c)** Multiplying all prices by 0.75 multiplies $\vpi$ and every reduced
cost by 0.75, which keeps their signs. The plan is unchanged and the revenue
becomes $0.75\cdot59{,}375=\$44{,}531.25$. For example, the reduced cost of
$x_4$ is $c_4-1.5c_3$ in general (since $\pi_1=c_3/4$), and with the new
prices it is $0.75\cdot30-1.5\cdot0.75\cdot25=-5.625<0$.

**(d)** $\bar c_4=c_4-\vpi^\T\mA_{\cdot4}=c_4-6\cdot6.25=c_4-37.5$. Type 4
enters the plan only if $\bar c_4>0$, that is if its price rises **by more
than \$7.50**, above \$37.50. At exactly \$37.50 the plan with 250 type-4
motors instead of the 375 type-3 motors earns the same \$59 375.

The full cost ranges, from {eq}`eq-03-07-basic-range`, are $c_1\ge50$,
$c_2\ge31.25$ and $20\le c_3\le30$ (dollars per motor), and $c_4\le37.5$.

*(An earlier version of this solution reports $s_4=0$ and $s_5=0$; the
correct values are $s_4=425$ and $s_5=750$. In (b) it names the relevant
nonbasic variable $s_2$; it is $s_3$, the slack of the type-2 line. The
answers \$8.75 and \$7.50 are correct.)*
:::

:::{prf:example} New objective functions for a paint-mix model
:label: ex-08-04-a6-q5-reddy-mikks

A paint maker solves

$$
\begin{aligned}
\max_{\vx}\ & 5x_1+4x_2 \\
\st\ & 6x_1+4x_2\le24,\\
& x_1+2x_2\le6,\\
& -x_1+x_2\le1,\\
& x_2\le2,\\
& x_1,x_2\ge0,
\end{aligned}
$$

with optimum $\vx^*=(3,\,1.5)$, $z^*=21$. Check whether this solution stays
optimal for each of the objectives below and, if not, find the new optimum
from the current basis: (a) $3x_1+2x_2$; (b) $8x_1+10x_2$; (c) $2x_1+5x_2$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

With slacks $s_1,\dots,s_4$, the optimal basis is $(x_1,x_2,s_3,s_4)$:

$$
\begin{gathered}
\mB=\begin{pmatrix}6&4&0&0\\1&2&0&0\\-1&1&1&0\\0&1&0&1\end{pmatrix},\\
\mB^{-1}=\begin{pmatrix}1/4&-1/2&0&0\\-1/8&3/4&0&0\\3/8&-5/4&1&0\\1/8&-3/4&0&1\end{pmatrix},
\end{gathered}
$$

$\vx_B=(3,\,3/2,\,5/2,\,1/2)^\T$, and for $\vc=(5,4)$ the multipliers are
$\vpi^\T=(3/4,\,1/2,\,0,\,0)$. The basis is primal feasible whatever the
objective, and the nonbasic variables are the slacks $s_1,s_2$, whose
reduced costs are $-\pi_1$ and $-\pi_2$. So the basis stays optimal exactly
when $\pi_1\ge0$ and $\pi_2\ge0$. With $\vc=(c_1,c_2)$,
$\vpi^\T=\left(\tfrac{c_1}{4}-\tfrac{c_2}{8},\ \tfrac{3c_2}{4}-\tfrac{c_1}{2},\ 0,\ 0\right)$,
so the condition is $\tfrac12\le c_1/c_2\le\tfrac32$: the gradient
$\vc$ must lie in the cone spanned by the normals $(6,4)$ and $(1,2)$ of the
two binding constraints ({numref}`fig-08-04-reddy-mikks-objectives`).

**(a)** $\vc=(3,2)$: $\vpi^\T=(1/2,\,0,\,0,\,0)\ge\vzero$, so $(3,1.5)$
remains optimal with $z^*=12$. But $\pi_2=0$ for a binding constraint means
$\bar c_{s_2}=0$: the objective is parallel to $6x_1+4x_2=24$, and every
point of the edge from $(3,1.5)$ to $(4,0)$ is optimal, with $z^*=12$
(alternative optima).

**(b)** $\vc=(8,10)$: $\vpi^\T=(3/4,\,7/2,\,0,\,0)\ge\vzero$; the solution
stays optimal with $z^*=24+15=39$.

**(c)** $\vc=(2,5)$: $\vpi^\T=(-1/8,\,11/4,\,0,\,0)$, so
$\bar c_{s_1}=1/8>0$ and the basis is no longer optimal. The tableau of the
old basis with the new costs is

| Basic | $x_1$ | $x_2$ | $s_1$ | $s_2$ | $s_3$ | $s_4$ | RHS |
|---|---|---|---|---|---|---|---|
| $x_1$ | 1 | 0 | 1/4 | −1/2 | 0 | 0 | 3 |
| $x_2$ | 0 | 1 | −1/8 | 3/4 | 0 | 0 | 3/2 |
| $s_3$ | 0 | 0 | 3/8 | −5/4 | 1 | 0 | 5/2 |
| $s_4$ | 0 | 0 | 1/8 | −3/4 | 0 | 1 | 1/2 |
| $-z$ | 0 | 0 | 1/8 | −11/4 | 0 | 0 | −27/2 |

$s_1$ enters; the ratios are 12 (row $x_1$), $20/3$ (row $s_3$) and 4 (row
$s_4$), so $s_4$ leaves:

| Basic | $x_1$ | $x_2$ | $s_1$ | $s_2$ | $s_3$ | $s_4$ | RHS |
|---|---|---|---|---|---|---|---|
| $x_1$ | 1 | 0 | 0 | 1 | 0 | −2 | 2 |
| $x_2$ | 0 | 1 | 0 | 0 | 0 | 1 | 2 |
| $s_3$ | 0 | 0 | 0 | 1 | 1 | −3 | 1 |
| $s_1$ | 0 | 0 | 1 | −6 | 0 | 8 | 4 |
| $-z$ | 0 | 0 | 0 | −2 | 0 | −1 | −14 |

All reduced costs are $\le0$: the new optimum is $x_1=2$, $x_2=2$ with
$s_1=4$ and $z^*=14$. HiGHS confirms all three cases.

```{figure} figures/08-04-reddy-mikks-objectives.svg
:name: fig-08-04-reddy-mikks-objectives
:width: 520px
:class: dark-light
:alt: Feasible region of the paint-mix LP with the cones of objective directions for which two corners are optimal, and the four objective gradients

Optimality cones of the paint-mix LP. A corner is optimal for every gradient
$\vc$ inside the cone of the normals of its binding constraints (shaded).
The current gradient $(5,4)$ and gradient (b), $(8,10)$, lie inside the
cone at $(3,1.5)$; gradient (a), $(3,2)$, lies on its edge, so a whole edge of
the region is optimal; gradient (c), $(2,5)$, lies outside it (dashed arrow)
but inside the cone at $(2,2)$, the new optimum.
```

*(An earlier version of this solution writes the result of (c) as
$x_3=4$; the problem has no $x_3$, and the value 4 belongs to the slack
$s_1$. It also does not mention the alternative optima in (a).)*
:::

(sec-08-04-rhs-ranging)=
## Right-hand-side ranging and shadow prices

:::{prf:example} Resource prices in an electronics plant
:label: ex-08-04-a6-q1-hidec

A firm assembles two models of an electronic device from resistors,
capacitors and chips. Model 1 uses 2 resistors, 2 capacitors and no chip;
model 2 uses 3 resistors, 1 capacitor and 4 chips. At most 1200 resistors,
1000 capacitors and 800 chips are available, and the models sell for \$3 and
\$4. With $x_1,x_2$ the numbers produced, the model and its optimal tableau
($s_1,s_2,s_3$ the slacks) are

$$
\begin{aligned}
\max_{\vx}\ & 3x_1+4x_2 \\
\st\ & 2x_1+3x_2\le1200,\\
& 2x_1+x_2\le1000,\\
& 4x_2\le800,\\
& x_1,x_2\ge0,
\end{aligned}
$$

| Basic | $x_1$ | $x_2$ | $s_1$ | $s_2$ | $s_3$ | RHS |
|---|---|---|---|---|---|---|
| $x_1$ | 1 | 0 | −1/4 | 3/4 | 0 | 450 |
| $s_3$ | 0 | 0 | −2 | 2 | 1 | 400 |
| $x_2$ | 0 | 1 | 1/2 | −1/2 | 0 | 100 |
| $-z$ | 0 | 0 | −5/4 | −1/4 | 0 | −1750 |

(a) Which resources are scarce and which are abundant? (b) What are the dual
prices of the three resources? (c) Find the optimum if 1300 resistors are
available. (d) If only 350 chips are available, can the new optimum be found
from this tableau? (e) A supplier offers extra resistors at \$0.40 each,
provided at least 500 are bought. Should the firm accept?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The columns of $s_1,s_2,s_3$ hold $\mB^{-1}$ for the basis
$(x_1,s_3,x_2)$, and the $-z$ row holds $-\vpi^\T$ under them.

**(a)** $s_1=s_2=0$: resistors and capacitors are used up (scarce).
$s_3=400$: 400 chips are left over (abundant).

**(b)** $\vpi=(5/4,\,1/4,\,0)$: an extra resistor is worth \$1.25 of revenue,
an extra capacitor \$0.25, an extra chip nothing.

**(c)** Column 1 of $\mB^{-1}$ is $(-1/4,-2,1/2)^\T$, so by
{eq}`eq-03-07-rhs-range` the resistor range is

$$
\begin{aligned}
450-\tfrac14\delta&\ge0,\\
400-2\delta&\ge0,\\
100+\tfrac12\delta&\ge0,
\end{aligned}
$$

which hold exactly for $-200\le\delta\le200$, that is, $1000\le b_1\le1400$ (capacitors give $800\le b_2\le1200$, chips
$b_3\ge400$). With $b_1=1300$ ($\delta=100$) the basis stays optimal:
$(x_1,s_3,x_2)=(425,\,200,\,150)$, and the revenue is
$1750+1.25\cdot100=\$1875$, at $x_1=425$, $x_2=150$.

**(d)** With 350 chips ($\delta=-450$, outside the range) the basic values
become $(x_1,s_3,x_2)=(450,\,-50,\,100)$: the basis is dual feasible but
infeasible, so **one or more dual simplex pivots from this tableau give the
new optimum**; there is no need to solve again. Row $s_3$ is
$(0,0,-2,2,1\mid-50)$; its only negative entry is in $s_1$, with ratio
$\abs{(-5/4)/(-2)}=5/8$, so $s_1$ enters:

| Basic | $x_1$ | $x_2$ | $s_1$ | $s_2$ | $s_3$ | RHS |
|---|---|---|---|---|---|---|
| $x_1$ | 1 | 0 | 0 | 1/2 | −1/8 | 456.25 |
| $s_1$ | 0 | 0 | 1 | −1 | −1/2 | 25 |
| $x_2$ | 0 | 1 | 0 | 0 | 1/4 | 87.5 |
| $-z$ | 0 | 0 | 0 | −3/2 | −5/8 | −1718.75 |

This is optimal: $x_1=456.25$, $x_2=87.5$, revenue \$1718.75; 25 resistors
are now left over, and chips are worth \$0.625 each. HiGHS agrees.

**(e)** The dual price \$1.25 exceeds \$0.40, but it is valid only up to
$b_1=1400$, that is, for 200 extra resistors. With 500 extra ($b_1=1700$)
revenue rises from \$1750 to \$2000 (at $b_1\ge1400$ the chip constraint
binds and the resistor price drops to 0): a gain of \$250 for a cost of
$500\cdot0.40=\$200$. **Accept, buying exactly 500**: the net gain is \$50,
although 300 of the resistors will not be used. Valuing all 500 at \$1.25
(a gain of \$625) would overstate the benefit; that is the error of using a
shadow price outside its range.

*(An earlier version of this solution answers (d) with "the problem must be
solved as a new problem"; one dual simplex pivot from the given tableau
suffices.)*
:::

:::{prf:example} Weekly feed schedule for a broiler farm
:label: ex-08-04-a6-q3-ozark

A farm raises 20 000 broilers over 8 weeks. The feed per bird in weeks 1 to 8
is 0.26, 0.48, 0.75, 1, 1.3, 1.6, 1.9 and 2.1 lb. The feed is a mix of
limestone, corn and soybean meal:

| Ingredient | Calcium (lb/lb) | Protein (lb/lb) | Fiber (lb/lb) | Cost (\$/lb) |
|---|---|---|---|---|
| Limestone | 0.38 | 0 | 0 | 0.12 |
| Corn | 0.001 | 0.09 | 0.02 | 0.45 |
| Soybean meal | 0.002 | 0.50 | 0.08 | 1.60 |

The mix must contain between 0.8 and 1.2 percent calcium, at least 22 percent
protein and at most 5 percent fiber. Find the cheapest mix for week 1, and use
post-optimal analysis to schedule the other seven weeks.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Model.** Let $Q$ be the week's total feed (lb) and $x_1,x_2,x_3$ the pounds
of limestone, corn and soybean meal:

$$
\begin{aligned}
\min_{\vx}\ & 0.12x_1+0.45x_2+1.60x_3 \\
\st\ & x_1+x_2+x_3=Q,\\
& 0.008Q\le0.38x_1+0.001x_2+0.002x_3\le0.012Q,\\
& 0.09x_2+0.50x_3\ge0.22Q,\\
& 0.02x_2+0.08x_3\le0.05Q,\\
& \vx\ge\vzero .
\end{aligned}
$$

**Week 1.** $Q=20{,}000\cdot0.26=5200$ lb. At the optimum the calcium
**upper** limit and the protein limit are binding, together with the mix
equation; solving these three equations gives the proportions

$$
\begin{aligned}
\frac{\vx^*}{Q}&=\left(\tfrac{219}{7774},\ \tfrac{2521}{3887},\ \tfrac{2513}{7774}\right)\\
&\approx(0.028171,\ 0.648572,\ 0.323257),
\end{aligned}
$$

that is, 146.5 lb of limestone, 3372.6 lb of corn and 1680.9 lb of soybean meal,
at a cost of \$4224.74 ($315{,}799/388{,}700\approx\$0.8124$ per lb). Fiber
is 3.88 percent and calcium at its lower limit is slack. The shadow prices
are about \$2.81 per lb of required protein and $-\$0.20$ per lb of allowed
calcium (relaxing the upper calcium limit would save money).

**Other weeks.** Every right-hand side is proportional to $Q$:
$\vb(Q)=Q\vb_0$. For the optimal basis,
$\vx_B(Q)=Q\,\mB^{-1}\vb_0\ge\vzero$ for every $Q>0$, and the reduced costs do
not depend on $\vb$. So the week-1 basis stays optimal in every week, and the
optimal mix is the same proportions times $Q$; the cost is
$0.812449\,Q$ dollars.

| Week | Feed $Q$ (lb) | Limestone (lb) | Corn (lb) | Soybean meal (lb) | Cost (\$) |
|---|---|---|---|---|---|
| 1 | 5 200 | 146.5 | 3 372.6 | 1 680.9 | 4 224.74 |
| 2 | 9 600 | 270.4 | 6 226.3 | 3 103.3 | 7 799.51 |
| 3 | 15 000 | 422.6 | 9 728.6 | 4 848.9 | 12 186.74 |
| 4 | 20 000 | 563.4 | 12 971.4 | 6 465.1 | 16 248.98 |
| 5 | 26 000 | 732.4 | 16 862.9 | 8 404.7 | 21 123.68 |
| 6 | 32 000 | 901.5 | 20 754.3 | 10 344.2 | 25 998.37 |
| 7 | 38 000 | 1 070.5 | 24 645.7 | 12 283.8 | 30 873.07 |
| 8 | 42 000 | 1 183.2 | 27 240.0 | 13 576.8 | 34 122.87 |
| Total | 187 800 | 5 290.5 | 121 801.9 | 60 707.7 | 152 577.96 |

HiGHS, solving each week's LP separately, gives the same amounts and costs.

*(An earlier version of this solution gives the cost as $0.81221\,Q$,
computed from proportions rounded to three decimals; the exact optimal cost
is $0.812449\,Q$, about \$45 more over the eight weeks.)*
:::

:::{prf:example} Reconstructing a final tableau
:label: ex-08-04-pq1-reconstruct

For

$$
\begin{aligned}
\max_{\vx}\ & c_1x_1+2x_2+4x_3 \\
\st\ & x_1+2x_2\le b_1,\\
& 2x_1-x_2+5x_3\le b_2,\\
& x_1,x_2,x_3\ge0,
\end{aligned}
$$

with slacks $x_4,x_5$, the optimal tableau is

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|---|---|---|---|---|---|---|
| $x_2$ | $u_1$ | 1 | 0 | $u_3$ | $u_5$ | 2 |
| $x_3$ | $u_2$ | 0 | 1 | $u_4$ | $u_6$ | 2 |
| $-z$ | −2 | 0 | 0 | $u_7$ | $u_8$ | $u_9$ |

(a) Find $c_1$, $b_1$, $b_2$ and $u_1,\dots,u_9$. (b) Find the range of $b_2$
over which this basis stays optimal, and the change in $z^*$ when $b_2$
grows by 1. (c) Find the range of $c_2$ (now 2) over which the basis stays
optimal, and the change in $z^*$ when $c_2$ grows by 1.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** The basis consists of $x_2,x_3$ with original columns $(2,-1)^\T$
and $(0,5)^\T$:

$$
\begin{gathered}
\mB=\begin{pmatrix}2&0\\-1&5\end{pmatrix},\\
\mB^{-1}=\begin{pmatrix}1/2&0\\1/10&1/5\end{pmatrix}.
\end{gathered}
$$

The slack columns of the tableau are $\mB^{-1}$: $u_3=1/2$, $u_5=0$,
$u_4=1/10$, $u_6=1/5$. Column $x_1$ is
$\mB^{-1}(1,2)^\T=(1/2,\,1/2)^\T$, so $u_1=u_2=1/2$. From
$\mB^{-1}\vb=(2,2)^\T$, $\vb=\mB(2,2)^\T=(4,8)^\T$: $b_1=4$, $b_2=8$. The
multipliers are $\vpi^\T=(2,4)\mB^{-1}=(7/5,\,4/5)$, so the reduced costs of
the slacks are $u_7=-7/5$ and $u_8=-4/5$. For $x_1$,
$\bar c_1=c_1-\vpi^\T(1,2)^\T=c_1-3=-2$ gives $c_1=1$. Finally
$u_9=-z^*=-(2\cdot2+4\cdot2)=-12$. HiGHS solves the reconstructed LP to
$\vx^*=(0,2,2)$, $z^*=12$.

**(b)** $\vx_B(b_2)=\mB^{-1}(4,b_2)^\T=(2,\ \tfrac25+\tfrac{b_2}{5})^\T\ge\vzero$
holds for $b_2\ge-2$. One more unit of $b_2$ raises $z^*$ by the shadow
price $\pi_2=4/5$ (HiGHS: $z^*=12.8$ at $b_2=9$).

**(c)** With $c_2$ as a symbol, $\vpi^\T=(c_2,4)\mB^{-1}=(\tfrac{c_2}{2}+\tfrac25,\ \tfrac45)$
and, for this maximization, the nonbasic reduced costs must stay $\le0$:

$$
\begin{gathered}
\bar c_1=1-\left(\tfrac{c_2}{2}+\tfrac25+\tfrac85\right)=-1-\tfrac{c_2}{2}\le0,\\
\bar c_4=-\tfrac{c_2}{2}-\tfrac25\le0,\\
\bar c_5=-\tfrac45\le0 .
\end{gathered}
$$

The binding condition is $c_2\ge-4/5$. Inside the range the solution does
not move, so one more unit of $c_2$ raises $z^*$ by $x_2^*=2$ (HiGHS:
$z^*=14$ at $c_2=3$; at $c_2=-0.8$ the basis ties with the solution
$(0,0,1.6)$, which takes over below).

*(The unknowns of the tableau are renamed $u_1,\dots,u_9$ here so that they
do not clash with the right-hand sides $b_i$.)*
:::

(sec-08-04-new-activities)=
## New activities and constraints

:::{prf:example} Scaling down a nonbasic activity
:label: ex-08-04-a5-q2-scaled-activity

An LP $\max\{\vc^\T\vx:\mA\vx\le\vb,\ \vx\ge\vzero\}$ with $\vb\ge\vzero$ and
$\vc\ge\vzero$ has been solved, and $x_j$ is nonbasic at the optimum. It is
proposed to make $x_j$ attractive by cutting its resource use per unit to
$1/q$ of the present values ($q>1$); its profit per unit would then also fall
to $1/q$ of its present value. Does this make $x_j$ worth producing?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Since $x_j$ is nonbasic, the basis $\mB$ and the multipliers
$\vpi^\T=\vc_B^\T\mB^{-1}$ do not involve its data and do not change. At the
optimum, $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}\le0$ (maximization). After the
change the column is $\mA_{\cdot j}/q$ and the profit $c_j/q$, so

$$
\bar c_j'=\frac{c_j}{q}-\vpi^\T\frac{\mA_{\cdot j}}{q}=\frac{\bar c_j}{q}\le0 .
$$

**No:** the reduced cost keeps its sign, the current solution stays
optimal, and $x_j$ stays at zero. Scaling a whole column, cost included, is
only a change of the unit in which $x_j$ is measured
({prf:ref}`ex-08-04-a5-q1-scaling`), and the unit cannot make an activity
profitable.

If only the resource use were cut, $\bar c_j'=c_j-\vpi^\T\mA_{\cdot j}/q$,
which is positive, so that $x_j$ would enter, if and only if
$q>\vpi^\T\mA_{\cdot j}/c_j$ (for $c_j>0$). Since
$\vpi^\T\mA_{\cdot j}\ge c_j$, this threshold is at least 1, and a cut with
$1<q$ below the threshold is not enough.

*(An earlier version of this solution states that cutting only the resource
use makes $x_j$ profitable for any $q>1$, and writes the entering condition
as $\bar c_j<0$; in a maximization $x_j$ enters when $\bar c_j>0$, which
needs $q>\vpi^\T\mA_{\cdot j}/c_j$. Its final conclusion is correct.)*
:::

:::{prf:example} Secondary constraints
:label: ex-08-04-a6-q4-secondary

Instead of solving an LP with all its constraints, one may first leave out
the **secondary** constraints, those least likely to be binding, solve with
the remaining **primary** ones, and then check the secondary constraints one
at a time: a satisfied constraint is discarded, a violated one is added and
the LP re-optimized. Apply this procedure to

$$
\begin{aligned}
\max_{\vx}\ & 5x_1+6x_2+3x_3 \\
\st\ & 5x_1+5x_2+3x_3\le50 &&(1)\\
& x_1+x_2-x_3\le20 &&(2)\\
& 7x_1+6x_2-9x_3\le30 &&(3)\\
& 5x_1+5x_2+5x_3\le35 &&(4)\\
& 12x_1+6x_2\le90 &&(5)\\
& x_2-9x_3\le20 &&(6)\\
& x_1,x_2,x_3\ge0 .
\end{aligned}
$$
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Choosing the primary set.** Constraint (4) bounds $x_1+x_2+x_3\le7$, and
the objective rewards every variable, so (4) is the natural primary
constraint. It also makes (1) redundant: for $\vx\ge\vzero$,
$5x_1+5x_2+3x_3\le5(x_1+x_2+x_3)\le35<50$. Take (4) as primary and (1), (2),
(3), (5), (6) as secondary.

**Step 1: primary constraint only.** $\max\ 5x_1+6x_2+3x_3$ subject to
$5x_1+5x_2+5x_3\le35$ puts everything into the variable with the best
objective per unit of the constraint, $x_2$: $\vx=(0,7,0)$, $z=42$. With
slack $s_4$, the optimal tableau is

| Basic | $x_1$ | $x_2$ | $x_3$ | $s_4$ | RHS |
|---|---|---|---|---|---|
| $x_2$ | 1 | 1 | 1 | 1/5 | 7 |
| $-z$ | −1 | 0 | −3 | −6/5 | −42 |

**Step 2: check the secondary constraints at $(0,7,0)$.** The slacks of (1),
(2), (3), (5), (6) are $15$, $13$, $-12$, $48$, $13$. Only (3) is violated
($42>30$).

**Step 3: add (3).** Append $7x_1+6x_2-9x_3+s_3=30$ and subtract 6 times row
$x_2$ to restore canonical form ({ref}`sec-03-06-adding-constraint`):

| Basic | $x_1$ | $x_2$ | $x_3$ | $s_4$ | $s_3$ | RHS |
|---|---|---|---|---|---|---|
| $x_2$ | 1 | 1 | 1 | 1/5 | 0 | 7 |
| $s_3$ | 1 | 0 | −15 | −6/5 | 1 | −12 |
| $-z$ | −1 | 0 | −3 | −6/5 | 0 | −42 |

The tableau is dual feasible (maximization, all $\bar c_j\le0$) with
$s_3=-12$. Row $s_3$ has negative entries in $x_3$ (ratio $3/15=1/5$) and
$s_4$ (ratio $\tfrac65/\tfrac65=1$), so $x_3$ enters:

| Basic | $x_1$ | $x_2$ | $x_3$ | $s_4$ | $s_3$ | RHS |
|---|---|---|---|---|---|---|
| $x_2$ | 16/15 | 1 | 0 | 3/25 | 1/15 | 31/5 |
| $x_3$ | −1/15 | 0 | 1 | 2/25 | −1/15 | 4/5 |
| $-z$ | −6/5 | 0 | 0 | −24/25 | −1/5 | −198/5 |

This is optimal for constraints (3) and (4): $\vx=(0,\,6.2,\,0.8)$,
$z=39.6$.

**Step 4: check the rest at $(0,6.2,0.8)$.** The slacks of (1), (2), (5),
(6) are $16.6$, $14.6$, $52.8$ and $21$, all nonnegative. No constraint
needs to be added, so

$$
\begin{aligned}
\vx^*&=(0,\ 6.2,\ 0.8),\\
z^*&=39.6
\end{aligned}
$$

is optimal for the full LP. HiGHS, solving the full LP directly, agrees.
Choosing (1), (3) and (4) as primary gives the same answer in one solve; the
order above shows the re-optimization step, which costs one dual simplex
pivot instead of a new solve.

*(An earlier version of this solution does not state the optimal value;
it is $z^*=39.6$.)*
:::

(sec-08-04-parametric)=
## Parametric analysis

:::{prf:example} A right-hand side that moves with a parameter
:label: ex-08-04-a6-q6-parametric

Study how the optimum of

$$
\begin{aligned}
\min_{\vx}\ & 4x_1+x_2+2x_3 \\
\st\ & 3x_1+x_2+2x_3=6+6t,\\
& 4x_1+3x_2+2x_3\ge12+4t,\\
& x_1+2x_2+5x_3\le8-2t,\\
& x_1,x_2,x_3\ge0
\end{aligned}
$$

changes for $t\ge0$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

With surplus $x_4$ and slack $x_5$, the constraints read
$3x_1+x_2+2x_3=6+6t$, $4x_1+3x_2+2x_3-x_4=12+4t$ and
$x_1+2x_2+5x_3+x_5=8-2t$, that is, $\vb(t)=\vb+t\vd$ with $\vb=(6,12,8)^\T$ and
$\vd=(6,4,-2)^\T$.

**The basis at $t=0$.** The optimal basis is $(x_1,x_2,x_4)$:

$$
\begin{gathered}
\mB=\begin{pmatrix}3&1&0\\4&3&-1\\1&2&0\end{pmatrix},\\
\mB^{-1}=\begin{pmatrix}2/5&0&-1/5\\-1/5&0&3/5\\1&-1&1\end{pmatrix}.
\end{gathered}
$$

$\vpi^\T=(4,1,0)\mB^{-1}=(7/5,\,0,\,-1/5)$, and the nonbasic reduced costs
are $\bar c_3=2-(14/5-1)=1/5$ and $\bar c_5=0-(-1/5)=1/5$, both
nonnegative: the basis is dual feasible for every $t$, since $t$ appears only
in $\vb$.

**Range of $t$.**

$$
\vx_B(t)=\mB^{-1}\vb(t)=
\begin{pmatrix}x_1\\x_2\\x_4\end{pmatrix}=
\begin{pmatrix}\tfrac45+\tfrac{14}{5}t\\[2pt]\tfrac{18}{5}-\tfrac{12}{5}t\\[2pt]2\end{pmatrix}\ge\vzero ,
$$

that is, $0\le t\le\tfrac32$.

For $0\le t\le3/2$ the optimum is $x_1=0.8+2.8t$, $x_2=3.6-2.4t$,
$x_3=0$ (with surplus $x_4=2$), and

$$
\begin{aligned}
z^*(t)&=4\left(\tfrac45+\tfrac{14}{5}t\right)+\tfrac{18}{5}-\tfrac{12}{5}t\\
&=\tfrac{34}{5}+\tfrac{44}{5}t=6.8+8.8t,
\end{aligned}
$$

whose slope is $\vpi^\T\vd=\tfrac75\cdot6+0\cdot4-\tfrac15\cdot(-2)=\tfrac{44}{5}$.
At $t=3/2$: $\vx=(5,0,0)$, $z^*=20$.

**Beyond $t=3/2$.** $x_2$ turns negative. Its row of $\mB^{-1}\mA$ is
$(0,1,\tfrac{13}{5},0,\tfrac35)$, with no negative entry, so the dual simplex
method stops: **the LP is infeasible for every $t>3/2$**. Directly:
eliminating $x_1=(6+6t-x_2-2x_3)/3$ from the third constraint gives
$\tfrac53x_2+\tfrac{13}{3}x_3\le6-4t$, which has no solution with
$x_2,x_3\ge0$ when $t>3/2$. HiGHS confirms $z^*=6.8$, $15.6$ and $20$ at
$t=0$, 1 and 1.5, and infeasibility at $t=1.51$.

*(An earlier version of this solution labels the third component of the
optimal vector $x_3^*$; it is the surplus $x_4^*=2$, and $x_3^*=0$. It also
omits the optimality check and $z^*(t)$.)*
:::

## References

- Bazaraa, M. S., Jarvis, J. J., & Sherali, H. D. (2010). *Linear
  programming and network flows* (4th ed.). Wiley. (Chapter 6: the dual
  simplex method, the artificial-constraint technique, sensitivity and
  parametric analysis.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Section 3.3: revised simplex; Section
  4.5: dual simplex; Chapter 5: sensitivity analysis.)
- Dantzig, G. B., & Orchard-Hays, W. (1954). The product form for the
  inverse in the simplex method. *Mathematical Tables and Other Aids to
  Computation, 8*(46), 64–67.
- Lemke, C. E. (1954). The dual method of solving the linear programming
  problem. *Naval Research Logistics Quarterly, 1*(1), 36–47.
  <https://doi.org/10.1002/nav.3800010107>
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapters 3, 4 and 7: the electronics, electric-motor, broiler-feed,
  paint-mix and secondary-constraint problems are adapted from the
  exercises and examples of these chapters.)
