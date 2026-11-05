# Duality

This chapter works through the examples and problems of
{doc}`03-05-lp-duality`.

The problems are grouped by topic: writing dual problems, weak and strong
duality, complementary slackness, testing a basis for primal and dual
feasibility, reading dual values from tableaux and from $\mB^{-1}$, and
shadow prices. The notation is that of {doc}`03-05-lp-duality`; every number was
recomputed in exact arithmetic and checked with the HiGHS solver.

(sec-08-03-writing-duals)=
## Writing dual problems

:::{prf:example} Three duals
:label: ex-08-03-a4-q1

Write the dual of each linear program, and solve both problems of each pair.

**(a)** $\max\ z=66x_1-22x_2$ subject to $-x_1+x_2\le-2$, $2x_1+3x_2\le5$,
$x_1,x_2\ge0$.

**(b)** $\min\ z=6x_1+3x_2$ subject to $6x_1-3x_2+x_3\ge25$,
$3x_1+4x_2+x_3\ge55$, $x_1,x_2,x_3\ge0$.

**(c)** $\max\ z=x_1+x_2$ subject to $2x_1+x_2=5$, $3x_1-x_2=6$, $x_1,x_2$
free.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

We apply {numref}`tab-03-05-dual-rules`, reading it right to left for a
maximization primal.

**(a)** Both constraints are $\le$ rows of a maximization problem (the natural
type), so both dual variables are nonnegative; both primal variables are
nonnegative, so both dual constraints are $\ge$:

$$
\begin{aligned}
\min_{\vpi}\ & w=-2\pi_1+5\pi_2\\
\st\ & -\pi_1+2\pi_2\ge 66\\
& \pi_1+3\pi_2\ge -22\\
& \pi_1,\pi_2\ge0 .
\end{aligned}
$$

An equivalent form multiplies the first primal row by $-1$ first
($x_1-x_2\ge2$, a constraint of the opposite type); its dual variable
$p_1=-\pi_1$ is then nonpositive and the dual reads
$\min\{2p_1+5p_2: p_1+2p_2\ge66,\ -p_1+3p_2\ge-22,\ p_1\le0,\ p_2\ge0\}$.

*Solutions.* The primal optimum is $\vx^*=(2.5,\,0)$ with $z^*=165$: the
second constraint is tight ($2\cdot2.5=5$) and the first has slack
$-2-(-2.5)=0.5$. By complementary slackness $\pi_1^*=0$, and since
$x_1^*>0$ the first dual constraint is tight: $2\pi_2=66$, so
$\vpi^*=(0,\,33)$. It is dual feasible ($0+99\ge-22$) and
$w^*=5\cdot33=165=z^*$, which certifies both.

**(b)** The constraints are $\ge$ rows of a minimization problem (natural
type) and the variables are nonnegative; $x_3$ has cost 0:

$$
\begin{aligned}
\max_{\vpi}\ & w=25\pi_1+55\pi_2\\
\st\ & 6\pi_1+3\pi_2\le 6\\
& -3\pi_1+4\pi_2\le 3\\
& \pi_1+\pi_2\le 0\\
& \pi_1,\pi_2\ge0 .
\end{aligned}
$$

The third constraint with $\vpi\ge\vzero$ forces $\vpi^*=(0,0)$, so
$w^*=0$. In the primal, $x_3$ costs nothing and appears in both rows with
coefficient 1, so $\vx=(0,0,55)$ is feasible with cost 0, and no feasible
point costs less because all costs are nonnegative: $z^*=0=w^*$. (Every
$\vx=(0,0,x_3)$ with $x_3\ge55$ is optimal.)

**(c)** Both constraints are equations, so both dual variables are free; both
primal variables are free, so both dual constraints are equations:

$$
\begin{aligned}
\min_{\vpi}\ & w=5\pi_1+6\pi_2\\
\st\ & 2\pi_1+3\pi_2=1\\
& \pi_1-\pi_2=1\\
& \pi_1,\pi_2\ \text{free}.
\end{aligned}
$$

Each problem has exactly one feasible point. Primal:
$5x_1=11$, so $\vx=(2.2,\,0.6)$ and $z=2.8$. Dual: $\pi_1=1+\pi_2$ gives
$5\pi_2=-1$, so $\vpi=(0.8,\,-0.2)$ and $w=4-1.2=2.8$. The values agree, as
strong duality requires.
:::

:::{prf:example} The dual of a problem in block form
:label: ex-08-03-block-dual

The variables of a minimization problem are split into a nonnegative block
$\vx_1\in\R^{n_1}$ and a free block $\vx_2\in\R^{n_2}$, and the constraints
into $m_1$ equations and $m_2$ inequalities of type $\ge$:

$$
\begin{aligned}
\min_{\vx_1,\vx_2}\ & \vc_1^\T\vx_1+\vc_2^\T\vx_2\\
\st\ & \mA_{11}\vx_1+\mA_{12}\vx_2=\vb_1\\
& \mA_{21}\vx_1+\mA_{22}\vx_2\ge\vb_2\\
& \vx_1\ge\vzero,\quad \vx_2\ \text{free},
\end{aligned}
$$

where $\mA_{kl}\in\R^{m_k\times n_l}$. Write its dual, and prove weak duality
for this pair directly.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Dual.** Give the equations the free dual block $\vpi_1\in\R^{m_1}$ and the
$\ge$ rows the nonnegative block $\vpi_2\in\R^{m_2}$. The columns of $\vx_1$
are $\begin{pmatrix}\mA_{11}\\\mA_{21}\end{pmatrix}$ and give $\le$
constraints ($\vx_1\ge\vzero$); the columns of $\vx_2$ give equations
($\vx_2$ free):

$$
\begin{aligned}
\max_{\vpi_1,\vpi_2}\ & \vb_1^\T\vpi_1+\vb_2^\T\vpi_2\\
\st\ & \mA_{11}^\T\vpi_1+\mA_{21}^\T\vpi_2\le\vc_1\\
& \mA_{12}^\T\vpi_1+\mA_{22}^\T\vpi_2=\vc_2\\
& \vpi_1\ \text{free},\quad \vpi_2\ge\vzero .
\end{aligned}
$$

The first dual constraint involves $\mA_{21}$, the block that multiplies
$\vx_1$ in the inequality rows. *(An earlier version of this pair has
$\mA_{12}$ in place of $\mA_{21}$ in that constraint; the transpose of the
column block of $\vx_1$ is $(\mA_{11}^\T\ \ \mA_{21}^\T)$.)*

**Weak duality.** Let $(\vx_1,\vx_2)$ and $(\vpi_1,\vpi_2)$ be feasible. Then

$$
\begin{aligned}
\vc_1^\T\vx_1+\vc_2^\T\vx_2
&\ge(\mA_{11}^\T\vpi_1+\mA_{21}^\T\vpi_2)^\T\vx_1\\
&\quad+(\mA_{12}^\T\vpi_1+\mA_{22}^\T\vpi_2)^\T\vx_2\\
&=\vpi_1^\T(\mA_{11}\vx_1+\mA_{12}\vx_2)\\
&\quad+\vpi_2^\T(\mA_{21}\vx_1+\mA_{22}\vx_2)\\
&\ge\vpi_1^\T\vb_1+\vpi_2^\T\vb_2 .
\end{aligned}
$$

The first inequality uses $\vx_1\ge\vzero$ with the $\le$ dual constraint
(the $\vx_2$ part is an equality, whatever the sign of $\vx_2$). In the last
line, the equation block gives $\vpi_1^\T\vb_1$ exactly, whatever the sign of
$\vpi_1$, and the inequality block gives at least $\vpi_2^\T\vb_2$ because
$\vpi_2\ge\vzero$. The four sign conventions are exactly what each step needs.
With $n_2=m_1=0$ the pair is the canonical pair {eq}`eq-03-05-canonical-pair`,
and with $n_2=m_2=0$ it is the standard-form pair {eq}`eq-03-05-standard-pair`.
:::

:::{prf:example} The dual constraint of an artificial variable
:label: ex-08-03-a4-q2

Consider

$$
\begin{aligned}
\max_{\vx}\ & z=5x_1+12x_2+4x_3\\
\st\ & x_1+2x_2+x_3\le 10\\
& 2x_1-x_2+3x_3=8\\
& x_1,x_2,x_3\ge0 .
\end{aligned}
$$

To start the big-M method, a slack variable $x_4$ is added to the first
row and an artificial variable $R$ with objective coefficient $-M$ to the
second ({ref}`sec-03-04-big-m`). Show that the artificial variable does not
change the dual: the dual constraint it creates is redundant once $M$ is
large enough. How large must $M$ be here?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Dual of the original problem.** The first row is a $\le$ row of a
maximization problem ($\pi_1\ge0$) and the second an equation ($\pi_2$ free):

$$
\begin{aligned}
\min_{\vpi}\ & w=10\pi_1+8\pi_2\\
\st\ & \pi_1+2\pi_2\ge5\\
& 2\pi_1-\pi_2\ge12\\
& \pi_1+3\pi_2\ge4\\
& \pi_1\ge0,\quad \pi_2\ \text{free}.
\end{aligned}
$$

**Dual of the big-M problem.** In equation form,
$\max\{5x_1+12x_2+4x_3+0x_4-MR:\ x_1+2x_2+x_3+x_4=10,\
2x_1-x_2+3x_3+R=8,\ \vx,R\ge0\}$. Both rows are equations, so both dual
variables are free, and every variable gives a $\ge$ constraint. The columns
of $x_1,x_2,x_3$ give the three constraints above; the column $(1,0)^\T$ of
$x_4$ with cost 0 gives $\pi_1\ge0$, which restores the sign restriction of
the original dual; and the column $(0,1)^\T$ of $R$ with cost $-M$ gives

$$
\pi_2\ge-M .
$$

So the two duals differ only by this lower bound on the free variable
$\pi_2$. As $M\to\infty$ it excludes nothing, which is the sense in which it
is redundant.

**How large is large enough.** The primal optimum is
$\vx^*=(26/5,\,12/5,\,0)$ with $z^*=274/5=54.8$, and the dual optimum is
$\vpi^*=(29/5,\,-2/5)$: it satisfies the first two dual constraints with
equality ($5.8-0.8=5$, $11.6+0.4=12$; both $x_1^*,x_2^*>0$), the third with
slack ($5.8-1.2=4.6\ge4$; $x_3^*=0$), and $w=58-3.2=54.8=z^*$. The extra
constraint does not cut off this point when $-2/5\ge-M$, that is, for
$M\ge2/5$; for every such $M$ the big-M problem has the same optimum, with
$R=0$. In the big-M tableau the column of $R$ shows
$\bar c_R=-M-\pi_2=-M+2/5$ ({numref}`tab-03-05-reading-duals`), which is
$-498/5$ for $M=100$ and is nonpositive (optimal for maximization) exactly
when $M\ge2/5$. For $M<2/5$ the penalty is too weak: with $M=0.1$, for
example, the big-M optimum keeps $R=13$ in the basis
($\vx=(0,5,0)$, $z-MR=58.7$), and the corresponding dual optimum sits on the
bound $\pi_2=-0.1$.
:::

:::{prf:example} The dual of a problem with bounded variables is feasible
:label: ex-08-03-a4-q10

Show that the dual of

$$
\begin{aligned}
\max_{\vx}\ & \vc^\T\vx\\
\st\ & \mA\vx\le\vb\\
& \vl\le\vx\le\vu ,
\end{aligned}
$$

with bound vectors $\vzero<\vl\le\vu$, always has a feasible solution. What
does this say about the primal?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Treat the bounds as constraints and $\vx$ as free (the bounds already imply
$\vx\ge\vl>\vzero$). The rows are $\mA\vx\le\vb$ (prices $\vp\ge\vzero$),
$\mI\vx\le\vu$ (prices $\vq\ge\vzero$) and $-\mI\vx\le-\vl$ (prices
$\vr\ge\vzero$), all $\le$ rows of a maximization problem. With $\vx$ free,
each dual constraint is an equation:

$$
\begin{aligned}
\min_{\vp,\vq,\vr}\ & w=\vb^\T\vp+\vu^\T\vq-\vl^\T\vr\\
\st\ & \mA^\T\vp+\vq-\vr=\vc\\
& \vp,\vq,\vr\ge\vzero .
\end{aligned}
$$

A feasible point is $\vp=\vzero$, $q_j=\max(c_j,0)$, $r_j=\max(-c_j,0)$:
then $q_j-r_j=c_j$ for every $j$, and all three vectors are nonnegative. (If
$\vx\ge\vzero$ is kept as a sign restriction, the dual constraints become
$\mA^\T\vp+\vq-\vr\ge\vc$, and the same point is feasible.)

*(The dual is a minimization problem; an earlier version of this solution
says "maximize".)*

Consequences for the primal follow from {numref}`tab-03-05-status`: since the
dual is feasible, the primal cannot be unbounded, which is also clear because
$\vx$ lies in a box. Either the primal is infeasible (the rows
$\mA\vx\le\vb$ may exclude the whole box), and then the dual is unbounded; or
it is feasible, and then both have optimal solutions with equal values. The
feasible dual point also gives the bound
$z^*\le\sum_j\bigl(u_j\max(c_j,0)-\ell_j\max(-c_j,0)\bigr)$: each variable at
its most favorable bound.
:::

(sec-08-03-weak-strong)=
## Weak and strong duality

:::{prf:example} Primal and dual both infeasible
:label: ex-08-03-both-infeasible

**(a)** Write the dual of
$\max\{x_1+3x_2:\ x_1+2x_2=1,\ x_1+2x_2=2,\ x_1,x_2\ \text{free}\}$ and show
that neither problem has a feasible solution.

**(b)** What changes if the sign restrictions $x_1,x_2\ge0$ are added?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** Two equations give free dual variables; two free primal variables
give dual equations:

$$
\begin{aligned}
\min_{\vpi}\ & w=\pi_1+2\pi_2\\
\st\ & \pi_1+\pi_2=1\\
& 2\pi_1+2\pi_2=3\\
& \vpi\ \text{free}.
\end{aligned}
$$

The primal asks $x_1+2x_2$ to equal both 1 and 2, which is impossible. The
dual asks $\pi_1+\pi_2$ to equal both 1 and $3/2$, which is also impossible.
This is the "both infeasible" entry of {numref}`tab-03-05-status`, and it shows
that an infeasible primal does not imply an unbounded dual.

**(b)** The primal stays infeasible. The dual constraints become
$\pi_1+\pi_2\ge1$ and $2\pi_1+2\pi_2\ge3$, which are satisfied, for example,
by every $\vpi=(t,\ 3/2-t)$; along these points $w=t+3-2t=3-t\to-\infty$ as
$t\to\infty$. Now the dual is unbounded, the "primal infeasible, dual
unbounded" entry of the table.
:::

:::{prf:example} The optimal value from the primal and from the dual
:label: ex-08-03-a4-q9

A maximization problem has two variables $x_1,x_2$ and three $\le$
constraints with slack variables $x_3,x_4,x_5$. Its optimal basis consists of
the columns of $x_1$, $x_2$ and $x_3$, with

$$
\begin{gathered}
\mB^{-1}=\begin{pmatrix}0&-1&1\\0&1&0\\1&1&-1\end{pmatrix},\\
\vx_B=(x_1,x_2,x_3)^\T=(1,3,1)^\T,\\
\vpi^*=(0,3,2)^\T .
\end{gathered}
$$

Compute the optimal objective value in two ways, once through the primal and
once through the dual.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Inverting $\mB^{-1}$ gives the basis matrix

$$
\mB=\begin{pmatrix}1&0&1\\0&1&0\\1&1&0\end{pmatrix},
$$

whose third column $(1,0,0)^\T$ is indeed the slack column of row 1, and whose
first two columns are the constraint columns of $x_1$ and $x_2$.

**Dual way.** The right-hand side (RHS) is $\vb=\mB\vx_B=(1+1,\ 3,\ 1+3)^\T=(2,3,4)^\T$,
so $w^*=\vb^\T\vpi^*=0\cdot2+3\cdot3+2\cdot4=17$.

**Primal way.** From $\vpi^{*\T}=\vc_B^\T\mB^{-1}$ we get
$\vc_B^\T=\vpi^{*\T}\mB=(0,3,2)\mB=(2,\ 5,\ 0)$, so $c_1=2$, $c_2=5$ and the
slack has cost 0, as it should. Then
$z^*=\vc_B^\T\vx_B=2\cdot1+5\cdot3+0\cdot1=17$.

Both give 17, as strong duality requires. The data describe the problem
$\max\{2x_1+5x_2:\ x_1\le2,\ x_2\le3,\ x_1+x_2\le4,\ \vx\ge\vzero\}$, whose
optimum $(1,3)$ leaves 1 unit of slack in the first row (so $\pi_1^*=0$) and
makes the other two rows tight.
:::

(sec-08-03-complementary-slackness)=
## Complementary slackness

:::{prf:example} A primal with all three kinds of constraint
:label: ex-08-03-mixed-dual

Consider

$$
\begin{aligned}
\min_{\vx}\ & z=5x_1-6x_2+7x_3+x_4\\
\st\ & x_1+2x_2-x_3-x_4=-7\\
& 6x_1-3x_2+x_3-7x_4\ge14\\
& -2.8x_1-17x_2+4x_3+2x_4\le-3\\
& x_1,x_2\ge0,\quad x_3,x_4\ \text{free}.
\end{aligned}
$$

**(a)** Write the dual.

**(b)** The point $\vx^*=\bigl(0,\ \tfrac{131}{35},\ \tfrac{554}{35},\
-\tfrac{47}{35}\bigr)$ is claimed to be optimal. Use complementary slackness to
find a dual solution and confirm the claim.

**(c)** Show that if the coefficient of $x_4$ in the second constraint is $+7$
instead of $-7$, the primal is unbounded, and conclude what happens to the
dual.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** The first row is an equation ($\pi_1$ free), the second a natural
$\ge$ row ($\pi_2\ge0$), the third a $\le$ row of a minimization problem
($\pi_3\le0$). The variables $x_1,x_2\ge0$ give $\le$ constraints and the free
$x_3,x_4$ give equations:

$$
\begin{aligned}
\max_{\vpi}\ & w=-7\pi_1+14\pi_2-3\pi_3\\
\st\ & \pi_1+6\pi_2-2.8\pi_3\le5 &&(x_1)\\
& 2\pi_1-3\pi_2-17\pi_3\le-6 &&(x_2)\\
& -\pi_1+\pi_2+4\pi_3=7 &&(x_3)\\
& -\pi_1-7\pi_2+2\pi_3=1 &&(x_4)\\
& \pi_1\ \text{free}\\
& \pi_2\ge0\\
& \pi_3\le0 .
\end{aligned}
$$

Multiplying the third primal row by $-1$ first (to make it a natural $\ge$
row) gives the same dual with $\pi_3'=-\pi_3\ge0$ and coefficients
$+2.8$, $+17$, $-4$, $-2$, $+3$ in place of $-2.8$, $-17$, $+4$, $+2$, $-3$.

**(b)** First, $\vx^*$ is feasible: $x_1^*=0$, $x_2^*>0$, and the three rows
evaluate to

$$
\begin{gathered}
\tfrac{262-554+47}{35}=-7,\\
\tfrac{-393+554+329}{35}=14,\\
\tfrac{-2227+2216-94}{35}=-3,
\end{gathered}
$$

so all three constraints are tight; its cost is
$z=\tfrac{-786+3878-47}{35}=\tfrac{3045}{35}=87$. Complementary slackness
requires the dual constraint of $x_2$ to be tight ($x_2^*>0$); the constraints
of $x_3$ and $x_4$ are equations anyway, and the primal rows are tight, so
they impose no condition on $\pi_2,\pi_3$. Solving

$$
\begin{aligned}
2\pi_1-3\pi_2-17\pi_3&=-6\\
-\pi_1+\pi_2+4\pi_3&=7\\
-\pi_1-7\pi_2+2\pi_3&=1
\end{aligned}
$$

gives $\vpi=(-10,\ 1,\ -1)$. The signs are right ($\pi_2=1\ge0$,
$\pi_3=-1\le0$), and the remaining dual constraint holds:
$-10+6+2.8=-1.2\le5$ (with slack, matching $x_1^*=0$). So $\vpi$ is dual
feasible and complementary to $\vx^*$; its value is
$w=70+14+3=87=z$. Hence $\vx^*$ and $\vpi^*=(-10,1,-1)$ are optimal, and
$z^*=w^*=87$.

**(c)** With $+7x_4$ in the second row, the point $(0,0,-10,17)$ is
feasible (row values $10-17=-7$, $-10+119=109\ge14$, $-40+34=-6\le-3$) with
cost $-53$. Along the direction $\vd=(0,0,-1,1)$ the row values change by
$1-1=0$, $-1+7=6\ge0$ and $-4+2=-2\le0$, so every point
$(0,0,-10-t,\ 17+t)$, $t\ge0$, stays feasible, while the cost changes by
$\vc^\T\vd=-7+1=-6$ per unit of $t$. The cost tends to $-\infty$: the primal
is unbounded, and by part 3 of {prf:ref}`cor-03-05-optimality-certificate` the
dual (now with $-\pi_1+7\pi_2+2\pi_3=1$ as its last constraint) is
infeasible.

*(An earlier version of this problem has $+7$ in the model and $-7$ in its
tables of coefficients. Only the version with $-7$ has an optimal solution; with $+7$ the pair falls into the "unbounded/infeasible" cell of
{numref}`tab-03-05-status`.)*
:::

(sec-08-03-feasibility-tests)=
## Dual feasibility and optimality tests

The next two examples use the same standard-form problem
$\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$ with

$$
\begin{gathered}
\mA=\begin{pmatrix}1&0&0&-1&2&4\\0&1&0&1&-2&6\\0&0&1&2&3&-7\end{pmatrix},\\
\vb=\begin{pmatrix}7\\19\\21\end{pmatrix},\\
\vc=(3,\ -4,\ -8,\ 16,\ -17,\ 19)^\T .
\end{gathered}
$$ (eq-08-03-six-variable-lp)

:::{prf:example} Is this basis optimal?
:label: ex-08-03-basis-optimal

For problem {eq}`eq-08-03-six-variable-lp`, take the basis with basic variables
$(x_6,x_2,x_3)$. Is it primal feasible? Is it dual feasible? What do you
conclude?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The basis matrix and its inverse are

$$
\begin{gathered}
\mB=\begin{pmatrix}4&0&0\\6&1&0\\-7&0&1\end{pmatrix},\\
\mB^{-1}=\begin{pmatrix}\tfrac14&0&0\\-\tfrac32&1&0\\\tfrac74&0&1\end{pmatrix}.
\end{gathered}
$$

**Primal feasibility.**
$\vx_B=\mB^{-1}\vb=\bigl(\tfrac74,\ -\tfrac{21}{2}+19,\ \tfrac{49}{4}+21\bigr)=\bigl(\tfrac74,\ \tfrac{17}{2},\ \tfrac{133}{4}\bigr)\ge\vzero$,
so the basis is primal feasible: $x_6=\tfrac74$, $x_2=\tfrac{17}{2}$,
$x_3=\tfrac{133}{4}$, the other variables zero.

**Dual feasibility.** With $\vc_B=(19,-4,-8)$,

$$
\begin{aligned}
\vpi^\T=\vc_B^\T\mB^{-1}&=\bigl(\tfrac{19}{4}+6-14,\ -4,\ -8\bigr)\\
&=\bigl(-\tfrac{13}{4},\ -4,\ -8\bigr).
\end{aligned}
$$

The dual constraints $\vpi^\T\mA_{\cdot j}\le c_j$ of the basic variables hold
with equality by construction; for the nonbasic ones,

$$
\begin{aligned}
x_1:&\quad \pi_1=-\tfrac{13}{4}\le3,\\
x_4:&\quad -\pi_1+\pi_2+2\pi_3=-\tfrac{67}{4}\le16,\\
x_5:&\quad 2\pi_1-2\pi_2+3\pi_3=-\tfrac{45}{2}\le-17.
\end{aligned}
$$

The reduced costs are the differences $\bar c_1=3+\tfrac{13}{4}=\tfrac{25}{4}$,
$\bar c_4=16+\tfrac{67}{4}=\tfrac{131}{4}$ and $\bar c_5=-17+\tfrac{45}{2}=\tfrac{11}{2}$.
All reduced costs, $\bar\vc=\bigl(\tfrac{25}{4},0,0,\tfrac{131}{4},\tfrac{11}{2},0\bigr)$,
are nonnegative, so the basis is dual feasible.

**Conclusion.** By {prf:ref}`prop-03-05-basis-optimal`, the basis is optimal:
$\vx^*=\bigl(0,\tfrac{17}{2},\tfrac{133}{4},0,0,\tfrac74\bigr)$ is optimal for
the primal and $\vpi^*=\bigl(-\tfrac{13}{4},-4,-8\bigr)$ for the dual, with

$$
\begin{aligned}
z^*&=19\cdot\tfrac74-4\cdot\tfrac{17}{2}-8\cdot\tfrac{133}{4}\\
&=\tfrac{133-136-1064}{4}=-\tfrac{1067}{4}=-266.75\\
&=\vb^\T\vpi^*=-\tfrac{91}{4}-76-168 .
\end{aligned}
$$
:::

:::{prf:example} Reduced costs are dual slacks
:label: ex-08-03-dual-slacks

For problem {eq}`eq-08-03-six-variable-lp`, take the basis $(x_1,x_2,x_3)$.
Compute its dual vector and the slacks of all six dual constraints. Which dual
constraints are violated, and what does the simplex method do next?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Here $\mB=\mI$, so $\vx_B=\vb=(7,19,21)\ge\vzero$ (primal feasible) and
$\vpi^\T=\vc_B^\T=(3,-4,-8)$. The dual constraints, written as "slack
$\ge0$", and their values are:

```{list-table} Dual slacks at $\vpi=(3,-4,-8)$
:header-rows: 1

* - Variable
  - Dual constraint
  - Slack $\bar c_j$
* - $x_1$
  - $3-\pi_1\ge0$
  - $0$
* - $x_2$
  - $-4-\pi_2\ge0$
  - $0$
* - $x_3$
  - $-8-\pi_3\ge0$
  - $0$
* - $x_4$
  - $16+\pi_1-\pi_2-2\pi_3\ge0$
  - $39$
* - $x_5$
  - $-17-2\pi_1+2\pi_2-3\pi_3\ge0$
  - $-7$
* - $x_6$
  - $19-4\pi_1-6\pi_2+7\pi_3\ge0$
  - $-25$
```

These slacks are exactly the reduced costs of the tableau of this basis,
whose objective row is $(0,0,0,39,-7,-25\mid 223)$, with $z=\vc_B^\T\vb=21-76-168=-223$.
The dual constraints of $x_5$ and $x_6$ are violated, so the basis is not dual
feasible, and since it is nondegenerate it is not optimal: increasing $x_5$
or $x_6$ lowers the cost.

The simplex method enters $x_6$ (most negative reduced cost). Its column is
$(4,6,-7)$, and the ratio test gives $\min\{7/4,\ 19/6\}=7/4$ in the row of
$x_1$, so $x_1$ leaves. The new basis $(x_6,x_2,x_3)$ is the optimal basis of
{prf:ref}`ex-08-03-basis-optimal`, and the cost falls by
$25\cdot\tfrac74=43.75$, from $-223$ to $-266.75$.
:::

:::{prf:example} Disproving optimality with the dual
:label: ex-08-03-a4-q7

Use the dual problem to show that the basic solution with basic variables
$(x_1,x_2)$ is not optimal for

$$
\begin{aligned}
\max_{\vx}\ & z=x_1+2x_2+4x_3-3x_4\\
\st\ & x_1+x_2+x_3=4\\
& x_1+4x_2+x_4=8\\
& \vx\ge\vzero ,
\end{aligned}
$$

and find the optimal solutions of both problems.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Dual.** Both rows are equations ($\vpi$ free) and the variables are
nonnegative, so for this maximization problem the dual is

$$
\begin{aligned}
\min_{\vpi}\ & w=4\pi_1+8\pi_2\\
\st\ & \pi_1+\pi_2\ge1 &&(x_1)\\
& \pi_1+4\pi_2\ge2 &&(x_2)\\
& \pi_1\ge4 &&(x_3)\\
& \pi_2\ge-3 &&(x_4)\\
& \pi_1,\pi_2\ \text{free}.
\end{aligned}
$$

**The basis $(x_1,x_2)$.** With $\mB=\begin{pmatrix}1&1\\1&4\end{pmatrix}$,
$\mB^{-1}=\tfrac13\begin{pmatrix}4&-1\\-1&1\end{pmatrix}$, the basic solution
is $\vx_B=\mB^{-1}\vb=\bigl(\tfrac83,\tfrac43\bigr)>\vzero$ with
$z=\tfrac83+\tfrac83=\tfrac{16}{3}$, and the multipliers are
$\vpi^\T=\vc_B^\T\mB^{-1}=(1,2)\cdot\tfrac13\begin{pmatrix}4&-1\\-1&1\end{pmatrix}=\bigl(\tfrac23,\tfrac13\bigr)$.
If the basis were optimal, this $\vpi$ would be optimal for the dual
({prf:ref}`thm-03-05-strong-duality`), hence feasible. But it violates
$\pi_1\ge4$, the dual constraint of $x_3$ (reduced cost
$\bar c_3=4-\tfrac23=\tfrac{10}{3}>0$, a profitable variable in a
maximization problem). Because the basic solution is nondegenerate, $x_3$ can
enter at a positive level, so the solution itself, not only the basis, is not
optimal.

**The optimum.** Entering $x_3$ (column $\mB^{-1}\mA_{\cdot3}=(\tfrac43,-\tfrac13)$)
removes $x_1$ (ratio $\tfrac83/\tfrac43=2$). The new basis $(x_3,x_2)$ gives
$\vx^*=(0,2,2,0)$ with $z^*=4+8=12$ and multipliers
$\vpi^*=(4,\ -\tfrac12)$. These satisfy all dual constraints
($3.5\ge1$, $2\ge2$, $4\ge4$, $-0.5\ge-3$), and $w=16-4=12=z^*$, so both are
optimal.
:::

:::{prf:example} A tableau from the inverse of its basis
:label: ex-08-03-a4-q5

For the problem

$$
\begin{aligned}
\min_{\vx}\ & z=2x_1+x_2\\
\st\ & 3x_1+x_2-x_3=3\\
& 4x_1+3x_2-x_4=6\\
& x_1+2x_2+x_5=3\\
& x_1,\dots,x_5\ge0,
\end{aligned}
$$

the basis with basic variables $(x_1,x_2,x_5)$ has

$$
\mB^{-1}=\begin{pmatrix}\tfrac35&-\tfrac15&0\\-\tfrac45&\tfrac35&0\\1&-1&1\end{pmatrix}.
$$

Compute the whole tableau of this basis and check it for feasibility and
optimality. Are the optimal solutions unique?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Tableau.** The constraint rows are $\mB^{-1}\mA$ and $\mB^{-1}\vb$, and with
$\vc_B=(2,1,0)$ the multipliers are
$\vpi^\T=\vc_B^\T\mB^{-1}=\bigl(\tfrac65-\tfrac45,\ -\tfrac25+\tfrac35,\ 0\bigr)=\bigl(\tfrac25,\tfrac15,0\bigr)$,
which give the objective row $\vc^\T-\vpi^\T\mA$ and $-z=-\vpi^\T\vb$:

```{list-table} Tableau of the basis $(x_1,x_2,x_5)$
:header-rows: 1
:stub-columns: 1

* - Basic
  - $x_1$
  - $x_2$
  - $x_3$
  - $x_4$
  - $x_5$
  - RHS
* - $x_1$
  - $1$
  - $0$
  - $-\tfrac35$
  - $\tfrac15$
  - $0$
  - $\tfrac35$
* - $x_2$
  - $0$
  - $1$
  - $\tfrac45$
  - $-\tfrac35$
  - $0$
  - $\tfrac65$
* - $x_5$
  - $0$
  - $0$
  - $-1$
  - $1$
  - $1$
  - $0$
* - $-z$
  - $0$
  - $0$
  - $\tfrac25$
  - $\tfrac15$
  - $0$
  - $-\tfrac{12}{5}$
```

As a check, $\mB^{-1}$ times the basic columns $(3,4,1)$, $(1,3,2)$,
$(0,0,1)$ gives the unit vectors, and the surplus columns of $x_3$, $x_4$
reproduce $-1$ times the first two columns of $\mB^{-1}$.

**Feasibility and optimality.** The right-hand side
$(\tfrac35,\tfrac65,0)$ is nonnegative, so the basis is primal feasible; the
reduced costs $(0,0,\tfrac25,\tfrac15,0)$ are nonnegative, so it is dual
feasible. The basis is therefore optimal, with $\vx^*=(\tfrac35,\tfrac65,0,0,0)$
and $z^*=\tfrac65+\tfrac65=\tfrac{12}{5}$, and $\vpi=(\tfrac25,\tfrac15,0)$ is
an optimal dual solution. The basis is **degenerate**: the basic variable
$x_5$ is zero (the third constraint, $x_1+2x_2\le3$, is tight at
$\vx^*$ although its dual value is 0).

**Uniqueness.** The primal optimum is unique, because every nonbasic variable
has a strictly positive reduced cost. The dual optimum is not. A degenerate
pivot that brings $x_3$ into the basis in place of $x_5$ (pivot on the entry
$-1$, allowed because the right-hand side is 0) gives the basis
$(x_1,x_2,x_3)$ at the same point $\vx^*$, with objective row
$(0,0,0,\tfrac35,\tfrac25\mid-\tfrac{12}{5})$, which is also dual feasible,
and multipliers $\vpi=(0,\tfrac35,-\tfrac25)$. Every point of the segment

$$
\vpi=\bigl(\tfrac25+t,\ \tfrac15-t,\ t\bigr),\qquad -\tfrac25\le t\le0,
$$

satisfies both dual constraints of $x_1,x_2$ with equality and the sign
conditions ($\pi_1,\pi_2\ge0$; $\pi_3\le0$ for the $\le$ row) and has value
$w=\tfrac{12}{5}$, so all of them are optimal. Degeneracy of the primal
optimum is what allows several optimal dual solutions. Correspondingly the
third constraint has two one-sided shadow prices: loosening it (raising the
right-hand side 3) does not lower the cost at all, while tightening it raises
the cost at the rate $\tfrac25$ per unit (for example $z^*=2.44$ when the
right-hand side is 2.9).
:::

(sec-08-03-tableau-duals)=
## Dual values from tableaux and $\mB^{-1}$

:::{prf:example} Duals from a big-M tableau (1)
:label: ex-08-03-a4-q3

The problem

$$
\begin{aligned}
\max_{\vx}\ & z=5x_1+2x_2+3x_3\\
\st\ & x_1+5x_2+2x_3=30\\
& x_1-5x_2-6x_3\le40\\
& \vx\ge\vzero
\end{aligned}
$$

was solved by the big-M method with $M=100$, an artificial variable $x_4$ in
the first row and a slack variable $x_5$ in the second, which formed the
starting basis. The optimal tableau is

```{list-table}
:header-rows: 1
:stub-columns: 1

* - Basic
  - $x_1$
  - $x_2$
  - $x_3$
  - $x_4$
  - $x_5$
  - RHS
* - $x_1$
  - $1$
  - $5$
  - $2$
  - $1$
  - $0$
  - $30$
* - $x_5$
  - $0$
  - $-10$
  - $-8$
  - $-1$
  - $1$
  - $10$
* - $-z$
  - $0$
  - $-23$
  - $-7$
  - $-105$
  - $0$
  - $-150$
```

Write the dual problem and find its optimal solution in two ways: from the
tableau and graphically.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Dual.** The first row is an equation ($\pi_1$ free), the second a natural
$\le$ row ($\pi_2\ge0$):

$$
\begin{aligned}
\min_{\vpi}\ & w=30\pi_1+40\pi_2\\
\st\ & \pi_1+\pi_2\ge5 &&(x_1)\\
& 5\pi_1-5\pi_2\ge2 &&(x_2)\\
& 2\pi_1-6\pi_2\ge3 &&(x_3)\\
& \pi_1\ \text{free},\quad \pi_2\ge0 .
\end{aligned}
$$

**From the tableau.** The starting basis was $(x_4,x_5)$, so the final
columns of $x_4$ and $x_5$ hold
$\mB^{-1}=\begin{pmatrix}1&0\\-1&1\end{pmatrix}$, and with $\vc_B=(5,0)$ for
$(x_1,x_5)$,
$\vpi^{*\T}=\vc_B^\T\mB^{-1}=(5,0)$. The objective row gives the same
directly ({numref}`tab-03-05-reading-duals`): the artificial $x_4$ has cost
$-M=-100$ and column $\ve_1$, so $-105=-100-\pi_1$ and $\pi_1=5$; the slack
$x_5$ has entry $0=-\pi_2$, so $\pi_2=0$.

**Graphically.** {numref}`fig-08-03-big-m-dual-region-1` shows the dual region.
The constraint of $x_2$ is redundant there. The region has two vertices,
$(5,0)$ with $w=150$ and $(\tfrac{33}{8},\tfrac78)$ (where the constraints of
$x_1$ and $x_3$ cross) with $w=\tfrac{635}{4}=158.75$, and it is unbounded
along the direction $(3,1)$, in which $w$ increases. The minimum is at
$\vpi^*=(5,0)$ with $w^*=150$.

**Check.** The primal optimum read from the tableau is $\vx^*=(30,0,0)$ with
$z^*=150=w^*$. Complementary slackness holds: $x_1^*>0$ and the constraint of
$x_1$ is tight ($5+0=5$); the second primal row has slack $40-30=10$ and
$\pi_2^*=0$.

```{figure} figures/08-03-big-m-dual-region-1.svg
:name: fig-08-03-big-m-dual-region-1
:width: 500px
:class: dark-light
:alt: Dual feasible region in the pi1-pi2 plane bounded by the lines of x1 and x3, with the optimum at (5, 0)

Feasible region of the dual of {prf:ref}`ex-08-03-a4-q3`. Line ($x_j$) is the
dual constraint of $x_j$; the dashed line is the contour $w=150$, which
touches the region only at $\vpi^*=(5,0)$.
```
:::

:::{prf:example} Duals from a big-M tableau (2)
:label: ex-08-03-a4-q4

The problem

$$
\begin{aligned}
\max_{\vx}\ & z=x_1+5x_2+3x_3\\
\st\ & x_1+2x_2+x_3=3\\
& 2x_1-x_2=4\\
& \vx\ge\vzero
\end{aligned}
$$

was solved by the big-M method with $M=100$, starting from $x_3$ in the
first row and an artificial variable $x_4$ in the second. The optimal tableau
is

```{list-table}
:header-rows: 1
:stub-columns: 1

* - Basic
  - $x_1$
  - $x_2$
  - $x_3$
  - $x_4$
  - RHS
* - $x_3$
  - $0$
  - $\tfrac52$
  - $1$
  - $-\tfrac12$
  - $1$
* - $x_1$
  - $1$
  - $-\tfrac12$
  - $0$
  - $\tfrac12$
  - $2$
* - $-z$
  - $0$
  - $-2$
  - $0$
  - $-99$
  - $-5$
```

Write the dual problem and find its optimal solution in two ways.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Dual.** Both rows are equations, so both dual variables are free:

$$
\begin{aligned}
\min_{\vpi}\ & w=3\pi_1+4\pi_2\\
\st\ & \pi_1+2\pi_2\ge1 &&(x_1)\\
& 2\pi_1-\pi_2\ge5 &&(x_2)\\
& \pi_1\ge3 &&(x_3)\\
& \pi_1,\pi_2\ \text{free}.
\end{aligned}
$$

($\pi_1$ is free in principle, but the constraint of $x_3$ makes it at least
3.)

**From the tableau.** The starting basis was $(x_3,x_4)$, so the final columns
of $x_3$ and $x_4$ hold
$\mB^{-1}=\begin{pmatrix}1&-\tfrac12\\0&\tfrac12\end{pmatrix}$; with
$\vc_B=(3,1)$ for $(x_3,x_1)$,
$\vpi^{*\T}=\vc_B^\T\mB^{-1}=\bigl(3,\ -\tfrac32+\tfrac12\bigr)=(3,\ -1)$.
Reading the objective row instead: the column of $x_3$ is $\ve_1$ with cost
3, so $0=3-\pi_1$ and $\pi_1=3$; the artificial $x_4$ has column $\ve_2$ and
cost $-100$, so $-99=-100-\pi_2$ and $\pi_2=-1$.

**Graphically.** In {numref}`fig-08-03-big-m-dual-region-2` the region lies to
the right of $\pi_1=3$, between the lines of $x_1$ and $x_2$. Its vertices are
$(3,-1)$ with $w=5$ and $(3,1)$ with $w=13$, and along the region $w$ grows
as $\pi_1$ grows (on the lower edge, $w=5+(\pi_1-3)$). The optimum is
$\vpi^*=(3,-1)$ with $w^*=5$.

**Check.** The primal optimum is $\vx^*=(2,0,1)$: $2+1=3$, $4=4$, and
$z^*=2+3=5=w^*$.

*(An earlier version of this tableau shows the column of $x_1$ as $(1,0)$
and the corner entry as $+5$. A basic column must be a unit vector with a 1 in
its own row, so the column of $x_1$ is $(0,1)$, and the corner of the $-z$
row holds $-z^*=-5$.)*

```{figure} figures/08-03-big-m-dual-region-2.svg
:name: fig-08-03-big-m-dual-region-2
:width: 500px
:class: dark-light
:alt: Dual feasible region with free variables, right of pi1 = 3, with the optimum at (3, -1)

Feasible region of the dual of {prf:ref}`ex-08-03-a4-q4`. Both dual variables
are free, so the region extends below $\pi_2=0$; the contour $w=5$ (dashed)
touches it at $\vpi^*=(3,-1)$.
```
:::

:::{prf:example} Reconstructing a tableau
:label: ex-08-03-a4-q6

For the problem

$$
\begin{aligned}
\max_{\vx}\ & z=5x_1+2x_2+3x_3\\
\st\ & x_1+5x_2+2x_3\le b_1\\
& x_1-5x_2-6x_3\le b_2\\
& \vx\ge\vzero ,
\end{aligned}
$$

with slack variables $x_4$ and $x_5$, the optimal tableau for some values of
$b_1$ and $b_2$ is

```{list-table}
:header-rows: 1
:stub-columns: 1

* - Basic
  - $x_1$
  - $x_2$
  - $x_3$
  - $x_4$
  - $x_5$
  - RHS
* - $x_1$
  - $1$
  - $b$
  - $2$
  - $1$
  - $0$
  - $3$
* - $x_5$
  - $0$
  - $c$
  - $-8$
  - $-1$
  - $1$
  - $1$
* - $-z$
  - $0$
  - $a$
  - $-7$
  - $d$
  - $e$
  - $-15$
```

Find **(a)** $b_1$ and $b_2$, **(b)** the optimal dual solution and **(c)** the
entries $a,b,c,d,e$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The slack columns $x_4,x_5$ formed the initial identity, so their final
columns give $\mB^{-1}=\begin{pmatrix}1&0\\-1&1\end{pmatrix}$ for the basis
$(x_1,x_5)$.

**(a)** $\mB^{-1}\vb=(3,1)$ means $b_1=3$ and $-b_1+b_2=1$, so $b_1=3$ and
$b_2=4$.

**(b)** $\vpi^{*\T}=\vc_B^\T\mB^{-1}=(5,0)\mB^{-1}=(5,0)$. Check:
$\vb^\T\vpi^*=15=z^*$, matching the corner entry $-15$.

**(c)** Column $x_2$ is $\mA_{\cdot2}=(5,-5)$, so
$\mB^{-1}\mA_{\cdot2}=(5,\ -5-5)=(5,-10)$: $b=5$, $c=-10$. Its reduced cost
is $a=c_2-\vpi^\T\mA_{\cdot2}=2-25=-23$. The slack $x_4$ has
$d=0-\pi_1=-5$, and $x_5$ is basic, so $e=0$. The given entries agree:
$\bar c_3=3-5\cdot2=-7$ and $\mB^{-1}\mA_{\cdot3}=(2,\ -6-2)=(2,-8)$. All
reduced costs are $\le0$, so the tableau is optimal for this maximization
problem, with $\vx^*=(3,0,0)$ and $z^*=15$.

*(An earlier version of this solution gives $d=-7$ and $e=-5$: $d$ is
minus the dual value of the first row, $-5$, and a basic variable has
reduced cost 0. The slack variables of this problem are $x_4$ and $x_5$, not
$x_3$ and $x_4$ as an earlier statement of the problem says.)*
:::

(sec-08-03-shadow-prices)=
## Shadow prices and economic interpretation

:::{prf:example} Shadow prices in cable production
:label: ex-08-03-a4-q8

A manufacturer makes four types of cable. Every unit passes through four
operations in sequence: splicing, soldering, sleeving and inspection. The
time each operation takes, the daily time available and the revenue per unit
are:

```{list-table}
:header-rows: 1
:stub-columns: 1

* - Cable type
  - Splicing (min/unit)
  - Soldering (min/unit)
  - Sleeving (min/unit)
  - Inspection (min/unit)
  - Revenue (\$/unit)
* - 1
  - 10.5
  - 20.4
  - 3.2
  - 5
  - 9.40
* - 2
  - 9.3
  - 24.6
  - 2.5
  - 5
  - 10.80
* - 3
  - 11.6
  - 17.7
  - 3.6
  - 5
  - 8.75
* - 4
  - 8.2
  - 26.5
  - 5.5
  - 5
  - 7.80
* - Capacity (min/day)
  - 4800
  - 9600
  - 4700
  - 4500
  -
```

A contract requires at least 100 units of each type per day.

**(a)** Formulate a linear program that maximizes daily revenue and find the
optimal production plan.

**(b)** Using the dual prices, would you recommend increasing the daily
capacity of any operation?

**(c)** Are the minimum production levels an advantage or a disadvantage for
the manufacturer? Explain with the dual prices.

**(d)** Can the revenue per extra minute given by the dual price be
guaranteed if the soldering capacity is increased by 10%?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** Let $x_j$ be the daily output of cable type $j$ (units/day):

$$
\begin{aligned}
\max_{\vx}\ & z=9.4x_1+10.8x_2+8.75x_3+7.8x_4\\
\st\ & 10.5x_1+9.3x_2+11.6x_3+8.2x_4\le4800 &&\text{(splicing)}\\
& 20.4x_1+24.6x_2+17.7x_3+26.5x_4\le9600 &&\text{(soldering)}\\
& 3.2x_1+2.5x_2+3.6x_3+5.5x_4\le4700 &&\text{(sleeving)}\\
& 5x_1+5x_2+5x_3+5x_4\le4500 &&\text{(inspection)}\\
& x_j\ge100 && j=1,\dots,4 .
\end{aligned}
$$

The optimal plan makes the minimum of types 1, 2 and 4 and puts all remaining
soldering time into type 3:

$$
\begin{aligned}
x_1^*&=x_2^*=x_4^*=100,\\
x_3^*&=\frac{9600-100\,(20.4+24.6+26.5)}{17.7}=\frac{2450}{17.7}\\
&=\frac{24500}{177}\approx138.42,
\end{aligned}
$$

with $z^*=\tfrac{709975}{177}\approx\$4011.16$ per day. Soldering is used to
capacity; splicing has $394.35$ minutes of slack, sleeving $3081.69$ and
inspection $2307.91$. (The plan is fractional; read it as an average daily
rate, or round $x_3$ down to 138 units, which stays feasible.)

**(b)** By complementary slackness, the three operations with slack have dual
price 0. The soldering price follows from the dual constraint of the basic
variable $x_3$, which must be tight: $17.7\,\pi_2=8.75$, so

$$
\pi_2^*=\frac{8.75}{17.7}=\frac{175}{354}\approx\$0.4944\ \text{per minute}.
$$

Only soldering capacity is worth increasing: each extra minute adds about
\$0.49 of revenue (as long as the basis stays optimal, see (d)), so it pays
if an extra minute costs less than that. Extra splicing, sleeving or
inspection time is worth nothing.

**(c)** The dual values of the lower bounds $x_j\ge100$ are the reduced costs
$\bar c_j=c_j-\pi_2^*a_{2j}$, with $a_{2j}$ the soldering time of type $j$:

$$
\begin{aligned}
\bar c_1&=9.4-0.4944\cdot20.4\approx-0.6847,\\
\bar c_2&\approx-1.3610,\\
\bar c_3&=0,\\
\bar c_4&\approx-5.3003\qquad\text{(\$/unit)}.
\end{aligned}
$$

Types 1, 2 and 4 earn less than the value of the soldering time they use
(for type 4, \$7.80 against $26.5\cdot0.4944=\$13.10$). Each unit of their
required minimum therefore costs revenue: lowering the minimum of type 1, 2 or
4 by one unit would raise $z^*$ by \$0.68, \$1.36 or \$5.30. The minimums are
a disadvantage. Without them the optimum is
$\vx\approx(0,\ 218.63,\ 238.51,\ 0)$ with $z\approx\$4448.20$ per day, so the
contract costs about \$437.04 of revenue per day.

**(d)** The price $\pi_2^*$ holds while the basis stays optimal, that is,
while $x_3=(b_2-7150)/17.7$ keeps all constraints satisfied. Writing $b_2$
for the soldering capacity:

- $x_3\ge100$ needs $b_2\ge7150+1770=8920$ minutes;
- splicing, $10.5\cdot100+9.3\cdot100+8.2\cdot100+11.6\,x_3\le4800$, needs
  $x_3\le\tfrac{5000}{29}\approx172.41$, that is
  $b_2\le\tfrac{295850}{29}\approx10201.72$ minutes;
- sleeving and inspection allow $x_3$ up to about 994 and 600 units, which
  are not binding.

So the price \$0.4944 per minute is valid for $8920\le b_2\le10201.72$
minutes, which allows an increase of at most $601.72$ minutes, or 6.27%. A
10% increase, to 10 560 minutes, lies outside this range, so **no**: the
price is not guaranteed for the whole increase. Indeed, beyond 10 201.72
minutes splicing also becomes binding and the soldering price drops to about
\$0.3636 per minute; at 10 560 minutes the optimal revenue is \$4438.89, a gain
of \$427.73 rather than the $960\cdot0.4944=\$474.58$ that the price predicts.

*(An earlier version of this solution gives the upper end of the range as
20 101.72 minutes; the correct value is 10 201.72 minutes. The maximum
increase is 6.27%.)*
:::

:::{prf:example} Lagrange multipliers as shadow prices
:label: ex-08-03-nlp-multipliers

Shadow prices also exist for nonlinear problems, as Lagrange multipliers.
Consider

$$
\begin{aligned}
\max_{\vx}\ & f(\vx)=(x_1-2)^2+(x_2-2)^2\\
\st\ & g_1(\vx)=(x_1-1)^2+(x_2-1)^2-1\le0\\
& g_2(\vx)=x_1+x_2-1\le0 .
\end{aligned}
$$

**(a)** Find the optimal solution. **(b)** Find multipliers
$\lambda_1,\lambda_2\ge0$ with
$\nabla f(\vx^*)=\lambda_1\nabla g_1(\vx^*)+\lambda_2\nabla g_2(\vx^*)$ and
$\lambda_ig_i(\vx^*)=0$. **(c)** Show that $\lambda_1$ is the rate at which
the optimal value grows when the first constraint is relaxed to
$g_1(\vx)\le\varepsilon$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** $f$ is the squared distance from the point $(2,2)$, and the first
constraint is the disc of radius 1 around $(1,1)$. The point of the disc
farthest from $(2,2)$ lies on the line through $(2,2)$ and the center, on the
far side of the center:
$\vx^*=(1,1)-\tfrac{1}{\sqrt2}(1,1)=\bigl(1-\tfrac{1}{\sqrt2},\ 1-\tfrac{1}{\sqrt2}\bigr)\approx(0.2929,\ 0.2929)$.
It satisfies the second constraint, $x_1+x_2=2-\sqrt2\approx0.586\le1$, so it
is optimal for the smaller feasible set too. Its distance from $(2,2)$ is
$\sqrt2+1$, so $f^*=(1+\sqrt2)^2=3+2\sqrt2\approx5.828$.

**(b)** The second constraint is not active ($g_2(\vx^*)=-0.414<0$), so
$\lambda_2=0$. The gradients at $\vx^*$ are
$\nabla f=2(\vx^*-(2,2))=-2\bigl(1+\tfrac{1}{\sqrt2}\bigr)(1,1)$ and
$\nabla g_1=2(\vx^*-(1,1))=-\tfrac{2}{\sqrt2}(1,1)$, so

$$
\lambda_1=\frac{1+1/\sqrt2}{1/\sqrt2}=\sqrt2+1\approx2.414 .
$$

Both conditions hold with $\lambda_1,\lambda_2\ge0$. (These are the Karush–Kuhn–Tucker (KKT)
conditions of {prf:ref}`thm-06-03-kkt`, written for a maximization problem. Because
$f$ is convex and is being maximized, they are necessary but not sufficient
here; optimality was shown geometrically in (a).)

**(c)** Relaxing $g_1\le\varepsilon$ enlarges the disc to radius
$\sqrt{1+\varepsilon}$, the farthest point moves outward along the same
line, and the second constraint stays inactive for small $\varepsilon$, so

$$
\begin{gathered}
f^*(\varepsilon)=\bigl(\sqrt2+\sqrt{1+\varepsilon}\bigr)^2,\\
\frac{\mathrm d f^*}{\mathrm d\varepsilon}\Big|_{\varepsilon=0}
=2\bigl(\sqrt2+1\bigr)\cdot\frac12=\sqrt2+1=\lambda_1 .
\end{gathered}
$$

The multiplier is the shadow price of the constraint: relaxing it by
$\varepsilon=\pm0.001$ changes $f^*$ by about $\pm0.00241$. The inactive
constraint has shadow price $\lambda_2=0$, as in linear programming.

*(An earlier statement of the KKT conditions for this problem writes the
complementarity condition of the second constraint with $g_1$; it should read
$\lambda_2g_2(\vx^*)=0$.)*
:::

## References

- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Sections 4.2–4.4: the dual problem,
  duality theorems, dual variables as marginal costs.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Sections 6.1–6.4 and 7.1–7.2: duality
  theory and its use in sensitivity analysis.)
- Karush, W. (1939). *Minima of functions of several variables with
  inequalities as side constraints* [Master's thesis, University of
  Chicago].
- Kuhn, H. W., & Tucker, A. W. (1951). Nonlinear programming. In J. Neyman
  (Ed.), *Proceedings of the Second Berkeley Symposium on Mathematical
  Statistics and Probability* (pp. 481–492). University of California Press.
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 4: Sections 4.1–4.3 and the problem sets of Sections 4.1,
  4.2 and 4.5, from which the cable-production problem, the big-M dual
  problems and the bounded-variable dual are adapted; Section 7.4: duality in
  matrix form.)
