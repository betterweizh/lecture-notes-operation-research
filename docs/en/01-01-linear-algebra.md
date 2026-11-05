# Linear Algebra Essentials

## Overview

Almost every method in this book is linear algebra in disguise. A linear
program is a system of linear equations $\mA\vx=\vb$ with sign constraints,
the simplex method moves between its solutions by Gauss–Jordan pivots, the
revised simplex method stores the inverse of a basis matrix as a product of
pivot matrices, and the optimality conditions of nonlinear programming are
statements about gradients, Hessians and positive definite matrices. This
chapter collects the linear algebra that the later chapters use. It starts
with vectors, matrices, inner products and linear independence, then treats
the solution of $\mA\vx=\vb$ by elimination, with special care for the
**pivot operation** and its **pivot matrices**
({prf:ref}`alg-01-01-pivot`, {prf:ref}`prop-01-01-pivot-matrix`,
{prf:ref}`prop-01-01-basis-exchange`). The four fundamental subspaces answer
when a system has no, one or infinitely many solutions
({prf:ref}`thm-01-01-existence-uniqueness`). Inverses, block matrices and
determinants follow. The second half is about square matrices: eigenvalues,
the spectral theorem for symmetric matrices ({prf:ref}`thm-01-01-spectral`),
quadratic forms and the tests for positive (semi)definiteness
({prf:ref}`thm-01-01-eigen-test`, {prf:ref}`thm-01-01-sylvester`,
{prf:ref}`thm-01-01-psd-minors`). The chapter ends with norms, orthogonality
and projections. It is a reference chapter: it can be read first, or
consulted when a later chapter needs a result.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- compute inner products, lengths, distances and angles in $\R^n$, and
  decide whether a set of vectors is linearly independent;
- solve $\mA\vx=\vb$ by Gaussian and Gauss–Jordan elimination, and write
  each elimination step and each pivot as multiplication by an elementary
  matrix;
- update the inverse of a basis matrix when one column is exchanged, in
  product form;
- find bases of the four fundamental subspaces of a matrix and decide
  whether $\mA\vx=\vb$ has no, one or infinitely many solutions;
- invert matrices, including partitioned matrices, and use the main
  properties of determinants;
- compute eigenvalues and eigenvectors, diagonalize a matrix and apply the
  spectral theorem to symmetric matrices;
- classify a quadratic form as positive (semi)definite, negative
  (semi)definite or indefinite with the eigenvalue and minor tests; and
- compute vector and matrix norms and orthogonal projections.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\ve_i$ | $i$-th unit vector (1 in position $i$, 0 elsewhere) |
| $\mI$ ($\mI_m$) | identity matrix (of order $m$) |
| $\mA_{\cdot j}$, $\mA_{i\cdot}$ | column $j$ and row $i$ of $\mA$ |
| $\col(\mA)$, $\row(\mA)$ | column space and row space of $\mA$ |
| $\nullsp(\mA)$, $\nullsp(\mA^\T)$ | null space and left null space of $\mA$ |
| $\rank\mA$ | rank of $\mA$ |
| $\mE$ | elementary matrix; a pivot (eta) matrix |
| $\mL$, $\mU$ | lower and upper triangular factors, $\mA=\mL\mU$ |
| $\bar\mA$, $\bar a_{rs}$ | a tableau and its pivot element |
| $\mB$, $\mN$; $\vx_B$, $\vx_N$ | basis matrix and nonbasic matrix of a linear program; the matching parts of $\vx$ |
| $\vpi$ | simplex multipliers, $\vpi^\T=\vc_B^\T\mB^{-1}$ |
| $\mT$ | bordered inverse of a basis matrix |
| $\det\mA$, $\tr\mA$ | determinant and trace of a square matrix |
| $\lambda_i(\mA)$ | $i$-th eigenvalue of $\mA$; $\lambda_{\min}(\mA)$, $\lambda_{\max}(\mA)$ the extreme ones of a symmetric $\mA$ |
| $\mV$, $\mD$ | matrix of eigenvectors and diagonal matrix of eigenvalues |
| $\succ$, $\succeq$ | positive definite, positive semidefinite |
| $\Sym^n$, $\Sym^n_+$, $\Sym^n_{++}$ | symmetric, positive semidefinite and positive definite matrices of order $n$ |
| $D_k$ | $k$-th leading principal minor |
| $\norm{\cdot}_1$, $\norm{\cdot}_2$, $\norm{\cdot}_\infty$ | vector norms and the matrix norms they induce; $\norm{\cdot}=\norm{\cdot}_2$ |
| $\norm{\cdot}_F$ | Frobenius norm |
| $\kappa(\mA)$ | condition number $\norm{\mA}\norm{\mA^{-1}}$ |
| $\mP$ | projection matrix |

(sec-01-01-vectors)=
## Vectors and matrices

A **vector** $\vx\in\R^n$ is an ordered list of $n$ real numbers, its
**components** $x_1,\dots,x_n$. Vectors are columns; the transpose
$\vx^\T=(x_1,\dots,x_n)$ is the same list written as a row. In running text
a column vector is often written $\vx=(x_1,\dots,x_n)$ to save space.
Vectors of the same length are added component by component and multiplied
by a scalar $\alpha\in\R$ component by component:
$(\vx+\vy)_j=x_j+y_j$ and $(\alpha\vx)_j=\alpha x_j$. The zero vector is
$\vzero$, the all-ones vector is $\vone$, and $\ve_i$ is the $i$-th unit
vector. For vectors, $\vx\ge\vy$ means $x_j\ge y_j$ for every $j$.

A **matrix** $\mA\in\R^{m\times n}$ is a rectangular array of $m$ rows and
$n$ columns with entries $a_{ij}$. Its column $j$ is $\mA_{\cdot j}\in\R^m$
and its row $i$ is $\mA_{i\cdot}$, a row vector of length $n$. The
transpose $\mA^\T\in\R^{n\times m}$ has entries $(\mA^\T)_{ij}=a_{ji}$. A
square matrix with $\mA^\T=\mA$ is **symmetric**. A square matrix is
**diagonal** when all entries off the diagonal are zero ($\diag(d_1,\dots,d_n)$),
**upper triangular** when all entries below the diagonal are zero and
**lower triangular** when all entries above it are zero. The **trace**
$\tr\mA=\sum_i a_{ii}$ is the sum of the diagonal entries.

The product of $\mA\in\R^{m\times n}$ and $\mB\in\R^{n\times p}$ is the
$m\times p$ matrix with entries

$$
(\mA\mB)_{ik}=\sum_{j=1}^n a_{ij}b_{jk}.
$$ (eq-01-01-matrix-product)

It is defined only when the number of columns of $\mA$ equals the number of
rows of $\mB$. Three readings of {eq}`eq-01-01-matrix-product` are used again
and again.

1. **Entry by entry.** $(\mA\mB)_{ik}$ is the inner product of row $i$ of
   $\mA$ with column $k$ of $\mB$.
2. **Column by column.** $\mA\vx=\sum_{j}x_j\mA_{\cdot j}$ is a combination
   of the columns of $\mA$ with weights $x_j$, and column $k$ of $\mA\mB$ is
   $\mA\mB_{\cdot k}$. A system of equations $\mA\vx=\vb$ therefore asks for
   weights that combine the columns of $\mA$ into $\vb$; this is the
   column view of a linear program used in {doc}`03-02-lp-geometry`
   ({eq}`eq-03-02-column-view`).
3. **Row by row.** Row $i$ of $\mA\mB$ is $\mA_{i\cdot}\mB$, a combination
   of the rows of $\mB$. Row operations are therefore multiplications from
   the left ({ref}`sec-01-01-elimination`).

Matrix multiplication is associative, $(\mA\mB)\mC=\mA(\mB\mC)$, and
distributive, but not commutative: in general $\mA\mB\ne\mB\mA$, even when
both products exist. The transpose of a product reverses the order,
$(\mA\mB)^\T=\mB^\T\mA^\T$. The identity $\mI$ satisfies $\mA\mI=\mI\mA=\mA$.
The trace is invariant under cyclic shifts, $\tr(\mA\mB)=\tr(\mB\mA)$ when
both products are square, a fact {doc}`01-02-matrix-calculus` uses to
differentiate with respect to matrices.

(sec-01-01-inner-products)=
## Inner products, lengths and angles

The **inner product** (dot product) of $\vx,\vy\in\R^n$ is the number

$$
\vx^\T\vy=\sum_{j=1}^n x_jy_j=x_1y_1+\dots+x_ny_n.
$$ (eq-01-01-inner-product)

It is symmetric ($\vx^\T\vy=\vy^\T\vx$) and linear in each argument. The
**Euclidean length** of $\vx$ is $\norm{\vx}=\sqrt{\vx^\T\vx}$ and the
**distance** between $\vx$ and $\vy$ is the length of their difference,

$$
\norm{\vx-\vy}=\sqrt{(x_1-y_1)^2+\dots+(x_n-y_n)^2}.
$$ (eq-01-01-distance)

:::{prf:theorem} Cauchy–Schwarz inequality
:label: thm-01-01-cauchy-schwarz

For all $\vx,\vy\in\R^n$, $\abs{\vx^\T\vy}\le\norm{\vx}\,\norm{\vy}$, with
equality if and only if one vector is a scalar multiple of the other.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

If $\vy=\vzero$ both sides are zero. Otherwise the quadratic polynomial in $t$
$\norm{\vx-t\vy}^2=\norm{\vx}^2-2t\,\vx^\T\vy+t^2\norm{\vy}^2$ is
nonnegative for every real $t$, so its discriminant is not positive:
$4(\vx^\T\vy)^2-4\norm{\vx}^2\norm{\vy}^2\le0$. Equality holds exactly when
the polynomial has a real root $t$, that is, when $\vx=t\vy$.
:::

Cauchy–Schwarz guarantees that $\vx^\T\vy/(\norm{\vx}\norm{\vy})$ lies in
$[-1,1]$ for nonzero vectors, so it is the cosine of an angle. The **angle**
$\varphi\in[0,\pi]$ between nonzero $\vx$ and $\vy$ is defined by

$$
\cos\varphi=\frac{\vx^\T\vy}{\norm{\vx}\,\norm{\vy}}.
$$ (eq-01-01-angle)

In $\R^2$ this agrees with elementary geometry: if $\vx$ and $\vy$ make
angles $\alpha$ and $\beta$ with the $x_1$-axis, then
$\vx=\norm{\vx}(\cos\alpha,\sin\alpha)$, $\vy=\norm{\vy}(\cos\beta,\sin\beta)$
and $\vx^\T\vy=\norm{\vx}\norm{\vy}(\cos\alpha\cos\beta+\sin\alpha\sin\beta)
=\norm{\vx}\norm{\vy}\cos(\alpha-\beta)$. Since the denominator in
{eq}`eq-01-01-angle` is positive, the sign of the inner product alone decides
the type of angle:

- $\vx^\T\vy>0$: the angle is **acute** ($\varphi<\pi/2$);
- $\vx^\T\vy=0$: the angle is a **right angle**; $\vx$ and $\vy$ are
  **orthogonal**, written $\vx\perp\vy$;
- $\vx^\T\vy<0$: the angle is **obtuse** ($\varphi>\pi/2$).

This trichotomy is used constantly in optimization: a direction $\vd$ is a
descent direction of $f$ at $\vx$ when it makes an obtuse angle with the
gradient, $\nabla f(\vx)^\T\vd<0$ ({prf:ref}`thm-01-02-steepest-ascent`,
{prf:ref}`def-06-02-descent-direction`).

For example, $\vx=(3,4)$ and $\vy=(4,-1)$ have $\vx^\T\vy=12-4=8$,
$\norm{\vx}=5$ and $\norm{\vy}=\sqrt{17}\approx4.123$, so
$\cos\varphi=8/(5\sqrt{17})\approx0.3881$ and $\varphi\approx67.17^\circ$,
an acute angle. Their distance is
$\norm{\vx-\vy}=\norm{(-1,5)}=\sqrt{26}\approx5.099$.

(sec-01-01-span)=
## Linear combinations, span and linear independence

A **linear combination** of vectors $\vx_1,\dots,\vx_k\in\R^n$ (bold with a
subscript: these are vectors, not components) is a vector
$\alpha_1\vx_1+\dots+\alpha_k\vx_k$ with real weights $\alpha_i$. The set of
all linear combinations is their **span**,
$\operatorname{span}\{\vx_1,\dots,\vx_k\}$. When the vectors are the columns
of a matrix $\mX=[\vx_1\ \cdots\ \vx_k]$, a linear combination is
$\mX\valpha$, and the span is the set of all $\mX\valpha$.

```{figure} figures/01-01-span-combination.svg
:name: fig-01-01-span-combination
:width: 680px
:class: dark-light
:alt: Left, two independent vectors u and v in the plane, the combination 1.5u plus v built by the parallelogram rule, and a skewed grid of lines along u and v covering the plane. Right, two vectors u and v equal to minus one half u lying on one line through the origin, which is their span

(a) The independent vectors $\vu$ and $\vv$ span $\R^2$: every point is a
combination $\alpha\vu+\beta\vv$, read off on the skewed grid; the orange
vector $1.5\vu+\vv$ is built by the parallelogram rule. (b) The dependent
pair $\vu$, $\vv=-0.5\vu$ spans only the line through the origin in the
direction of $\vu$.
```

:::{prf:definition} Linear independence
:label: def-01-01-independence

Vectors $\vx_1,\dots,\vx_k\in\R^n$ are **linearly independent** if

$$
\alpha_1\vx_1+\dots+\alpha_k\vx_k=\vzero\quad\Longrightarrow\quad
\alpha_1=\dots=\alpha_k=0,
$$

that is, if the only combination that gives the zero vector is the trivial
one. Otherwise they are **linearly dependent**: some nontrivial combination
vanishes, and any vector with a nonzero weight in it is a combination of
the others.
:::

{numref}`fig-01-01-span-combination` shows the two cases in $\R^2$.
Independence is easy to define but needs a test. Writing the vectors as the
columns of $\mX$, the definition says that $\mX\valpha=\vzero$ has only the
solution $\valpha=\vzero$, and this homogeneous system can be solved by
elimination.

:::{prf:proposition} Testing linear independence
:label: prop-01-01-independence-test

Let $\mX=[\vx_1\ \cdots\ \vx_k]\in\R^{n\times k}$. The following are
equivalent:

1. $\vx_1,\dots,\vx_k$ are linearly independent;
2. $\mX\valpha=\vzero$ has only the solution $\valpha=\vzero$;
3. elimination on $\mX$ produces a pivot in every column, that is,
   $\rank\mX=k$ ({ref}`sec-01-01-basis-rank`).

In particular, more than $n$ vectors in $\R^n$ are always dependent, and
$n$ vectors in $\R^n$ are independent if and only if $\det\mX\ne0$
({ref}`sec-01-01-determinants`).
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

(1)⇔(2) is the definition written in matrix form. Elimination does not
change the solution set of $\mX\valpha=\vzero$. A column without a pivot
belongs to a free variable, which can be set to $1$ to give a nonzero
solution; if every column has a pivot, back substitution forces every
$\alpha_i=0$. This gives (2)⇔(3). There are at most $n$ pivots, one per
row, so $k>n$ columns cannot all have one. For $k=n$, a pivot in every
column means that $\mX$ is invertible, which is equivalent to $\det\mX\ne0$
by {prf:ref}`thm-01-01-invertible`.
:::

:::{prf:example} Checking independence by elimination
:label: ex-01-01-independence

Are $\vx_1=(1,0,2)$, $\vx_2=(0,1,1)$ and $\vx_3=(2,3,7)$ linearly
independent?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Solve $\mX\valpha=\vzero$ with $\mX=[\vx_1\ \vx_2\ \vx_3]$. The rows of
$\mX$ are $(1,0,2)$, $(0,1,3)$ and $(2,1,7)$. Subtracting twice row 1 from
row 3 gives $(0,1,3)$; subtracting row 2 from the new row 3 gives
$(0,0,0)$. Only two pivots appear, so $\rank\mX=2<3$ and the vectors are
dependent. The third variable is free; with $\alpha_3=1$, row 2 gives
$\alpha_2=-3$ and row 1 gives $\alpha_1=-2$. Check:
$-2\vx_1-3\vx_2+\vx_3=(-2+0+2,\ 0-3+3,\ -4-3+7)=\vzero$. So
$\vx_3=2\vx_1+3\vx_2$ lies in the plane spanned by $\vx_1$ and $\vx_2$, and
$\det\mX=0$.
:::

(sec-01-01-basis-rank)=
## Subspaces, bases, dimension and rank

A **subspace** of $\R^n$ is a nonempty set $S$ that is closed under
addition and scalar multiplication: $\vx,\vy\in S$ and $\alpha\in\R$ imply
$\vx+\vy\in S$ and $\alpha\vx\in S$. Every subspace contains $\vzero$.
Lines and planes through the origin are subspaces of $\R^3$; a line that
misses the origin is not (it is an affine set, {prf:ref}`def-02-01-affine-set`). The
span of any set of vectors is a subspace.

A **basis** of a subspace $S$ is a set of linearly independent vectors that
spans $S$. Every vector of $S$ is then a combination of the basis vectors
in exactly one way: if $\sum\alpha_i\vx_i=\sum\beta_i\vx_i$, then
$\sum(\alpha_i-\beta_i)\vx_i=\vzero$ and independence gives
$\alpha_i=\beta_i$. All bases of $S$ have the same number of vectors, the
**dimension** $\dim S$. The unit vectors $\ve_1,\dots,\ve_n$ form the
standard basis of $\R^n$, so $\dim\R^n=n$.

The **column space** $\col(\mA)\subseteq\R^m$ of
$\mA\in\R^{m\times n}$ is the span of its columns, that is, the set of all
$\mA\vx$. The **row space** $\row(\mA)=\col(\mA^\T)
\subseteq\R^n$ is the span of its rows.

:::{prf:theorem} Row rank equals column rank
:label: thm-01-01-rank

For every matrix $\mA$, $\dim\col(\mA)=\dim\row(\mA)$.
This common number is the **rank** of $\mA$, $\rank\mA$. It equals the
number of pivots produced by elimination, and
$\rank\mA\le\min\{m,n\}$.
:::

The idea of the proof is that row operations change neither the row space
nor the dependencies among the columns. In the echelon form produced by
elimination ({ref}`sec-01-01-elimination`), the nonzero rows are independent
and span the row space, and the pivot columns are independent and every
other column is a combination of them. Both counts equal the number of
pivots, and the pivot columns of the *original* matrix form a basis of its
column space (Strang, 2016, Section 3.5). A matrix $\mA\in\R^{m\times n}$
has **full row rank** if $\rank\mA=m$ and **full column rank** if
$\rank\mA=n$. Linear programs in standard form usually assume full row
rank, $\rank\mA=m$: no equation is redundant ({prf:ref}`asm-03-02-full-rank`,
{doc}`03-04-two-phase-big-m`).

(sec-01-01-elimination)=
## Solving linear systems by elimination

A system of $m$ linear equations in $n$ unknowns is written $\mA\vx=\vb$
with $\mA\in\R^{m\times n}$, $\vb\in\R^m$. All the information sits in the
**augmented matrix** $[\mA\mid\vb]$. Three **elementary row operations**
transform it without changing the solution set:

1. add a multiple of one row to another row;
2. multiply a row by a nonzero scalar;
3. exchange two rows.

Each operation is reversible, so no solution is lost or gained. Each is
also a multiplication from the left by an **elementary matrix** $\mE$, the
identity with the same operation applied to it. For example, with three
rows, "subtract $\ell$ times row 1 from row 2" is

$$
\begin{aligned}
\mE_{21}&=\begin{pmatrix}1&0&0\\-\ell&1&0\\0&0&1\end{pmatrix},\\
\mE_{21}^{-1}&=\begin{pmatrix}1&0&0\\\ell&1&0\\0&0&1\end{pmatrix}.
\end{aligned}
$$ (eq-01-01-elementary-matrix)

The inverse adds the row back. Scaling row $r$ by $\gamma\ne0$ is
$\diag(1,\dots,\gamma,\dots,1)$, and a row exchange is a **permutation
matrix**, the identity with two rows exchanged.

**Gaussian elimination** uses the first operation to create zeros below the
diagonal, column by column. The entry used to eliminate a column is the
**pivot**; it must be nonzero.

:::{prf:algorithm} Gaussian elimination with back substitution
:label: alg-01-01-gaussian-elimination

**Input:** $\mA\in\R^{n\times n}$, $\vb\in\R^n$.
**Output:** the solution $\vx$ of $\mA\vx=\vb$, or the message that $\mA$
is singular.

1. **Forward elimination.** For $k=1,\dots,n-1$: if $a_{kk}=0$, exchange
   row $k$ with a row $i>k$ that has $a_{ik}\ne0$ (if there is none, report
   that $\mA$ is singular). For each row $i>k$, compute the multiplier
   $\ell_{ik}=a_{ik}/a_{kk}$ and subtract $\ell_{ik}$ times row $k$ from
   row $i$ (both in $\mA$ and in $\vb$).
2. **Check.** If $a_{nn}=0$, report that $\mA$ is singular. The system is
   now $\mU\vx=\vc$ with $\mU$ upper triangular and nonzero pivots
   $u_{11},\dots,u_{nn}$.
3. **Back substitution.** For $i=n,n-1,\dots,1$: set
   $x_i=\bigl(c_i-\sum_{j>i}u_{ij}x_j\bigr)/u_{ii}$.
:::

The multipliers are not thrown away: they form the lower triangular factor
$\mL$ of the factorization $\mA=\mL\mU$ (when no row exchanges are needed).
In floating-point arithmetic, **partial pivoting** exchanges rows so that
the pivot is the entry of largest absolute value in its column; this keeps
the multipliers at most 1 in absolute value and controls rounding errors.

:::{prf:example} Elimination, elementary matrices and $\mA=\mL\mU$
:label: ex-01-01-elimination-lu

Solve $\mA\vx=\vb$ with

$$
\begin{aligned}
\mA&=\begin{pmatrix}2&1&1\\4&3&3\\8&7&9\end{pmatrix},\\
\vb&=\begin{pmatrix}3\\7\\19\end{pmatrix}.
\end{aligned}
$$

Write the elimination steps as elementary matrices, find the factors $\mL$
and $\mU$ of $\mA=\mL\mU$, and compute $\det\mA$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Forward elimination.** The first pivot is $a_{11}=2$. The multipliers
$\ell_{21}=4/2=2$ and $\ell_{31}=8/2=4$ give rows $(0,1,1\mid1)$ and
$(0,3,5\mid7)$. The second pivot is $1$; the multiplier $\ell_{32}=3$ turns
row 3 into $(0,0,2\mid4)$. In matrix form,
$\mE_{32}\mE_{31}\mE_{21}[\mA\mid\vb]=[\mU\mid\vc]$ with

$$
\begin{aligned}
\mU&=\begin{pmatrix}2&1&1\\0&1&1\\0&0&2\end{pmatrix},\\
\vc&=\begin{pmatrix}3\\1\\4\end{pmatrix}.
\end{aligned}
$$

**Back substitution.** $x_3=4/2=2$, $x_2=(1-x_3)/1=-1$,
$x_1=(3-x_2-x_3)/2=1$, so $\vx=(1,-1,2)$. Check:
$\mA\vx=(2-1+2,\ 4-3+6,\ 8-7+18)=(3,7,19)=\vb$.

**The factor $\mL$.** Undoing the steps in reverse order gives
$\mA=\mE_{21}^{-1}\mE_{31}^{-1}\mE_{32}^{-1}\mU=\mL\mU$, and the product of
the inverses simply collects the multipliers below the diagonal:

$$
\mL=\begin{pmatrix}1&0&0\\2&1&0\\4&3&1\end{pmatrix}.
$$

The pivots $2,1,2$ are the diagonal of $\mU$, and
$\det\mA=2\cdot1\cdot2=4$ ({ref}`sec-01-01-determinants`). With $\mA=\mL\mU$
stored, every further right-hand side costs only two triangular solves.
:::

**Gauss–Jordan elimination** continues after the forward phase: it divides
each pivot row by its pivot and also eliminates the entries *above* each
pivot. The result is the **reduced row echelon form** (RREF), in which each
pivot column is a unit vector. Variables of pivot columns are **basic** (or
pivot) variables, the others are **free** variables; the RREF expresses
every basic variable in terms of the free ones. This is exactly the
canonical form of a linear program ({prf:ref}`def-03-03-canonical-form`): the basic
variables are those of the identity columns, and setting the free
(nonbasic) variables to zero reads off a basic solution.

(sec-01-01-pivot)=
## The pivot operation and pivot matrices

The simplex method, the revised simplex method and the dual simplex method
({doc}`03-03-simplex-method`, {doc}`03-06-revised-dual-simplex`) are built from
one operation: a Gauss–Jordan pivot on a single chosen entry of a tableau.

:::{prf:algorithm} Pivot on entry $(r,s)$
:label: alg-01-01-pivot

**Input:** a tableau $\bar\mA\in\R^{m\times q}$ (for example $[\mA\mid\vb]$), a
pivot row $r$ and a pivot column $s$ with **pivot element** $\bar a_{rs}\ne0$.
**Output:** the tableau in which column $s$ is the unit vector $\ve_r$.

1. **Normalize the pivot row.** Divide row $r$ by $\bar a_{rs}$:
   $\bar a'_{rj}=\bar a_{rj}/\bar a_{rs}$ for all $j$.
2. **Eliminate column $s$ elsewhere.** For every row $i\ne r$, subtract
   $\bar a_{is}$ times the new row $r$: $\bar a'_{ij}=\bar a_{ij}-\bar a_{is}\,\bar a_{rj}/\bar a_{rs}$
   for all $j$.
:::

The update in step 2 is the "rectangle rule": the new entry $\bar a'_{ij}$ is
the old entry minus the product of the two entries in the pivot row and
pivot column that complete the rectangle, divided by the pivot. A pivot is
a sequence of elementary row operations (one scaling and $m-1$ row
additions), so it is multiplication from the left by their product.

:::{prf:proposition} Pivot matrices
:label: prop-01-01-pivot-matrix

Let $\vw=\bar\mA_{\cdot s}$ be the pivot column, with $w_r\ne0$. The pivot of
{prf:ref}`alg-01-01-pivot` replaces $\bar\mA$ by $\mE\bar\mA$, where the **pivot
matrix** (also called an **eta matrix**) $\mE$ is the identity with column
$r$ replaced by the vector $\veta$,

$$
\begin{aligned}
\mE&=\bigl[\ve_1\ \cdots\ \ve_{r-1}\ \ \veta\ \ \ve_{r+1}\ \cdots\ \ve_m\bigr],\\
\eta_r&=\frac{1}{w_r},\\
\eta_i&=-\frac{w_i}{w_r}\quad (i\ne r).
\end{aligned}
$$ (eq-01-01-eta-matrix)

Its inverse is the identity with column $r$ replaced by $\vw$ itself.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Row $i$ of $\mE\bar\mA$ is $\sum_k e_{ik}\bar\mA_{k\cdot}$. Only column $r$ of
$\mE$ differs from the identity, so row $r$ of $\mE\bar\mA$ is
$\eta_r\bar\mA_{r\cdot}=\bar\mA_{r\cdot}/w_r$ and row $i\ne r$ is
$\bar\mA_{i\cdot}+\eta_i\bar\mA_{r\cdot}=\bar\mA_{i\cdot}-(w_i/w_r)\bar\mA_{r\cdot}$; these
are the two steps of the pivot. For the inverse, let $\mF$ be the identity
with column $r$ replaced by $\vw$. Then $\mE\mF$ agrees with $\mE$ outside
column $r$, and its column $r$ is $\mE\vw$, which is $\ve_r$ because the
pivot turns column $s$ (that is, $\vw$) into $\ve_r$. Hence $\mE\mF=\mI$.
:::

Only one column of $\mE$ has to be stored, which is why the revised simplex
method keeps a list of eta vectors instead of a full inverse. The next
result explains why a pivot on a *zero* entry is impossible and how the
inverse of a basis matrix changes when one column is exchanged.

:::{prf:proposition} Exchanging one column of a basis
:label: prop-01-01-basis-exchange

Let $\mB\in\R^{m\times m}$ be invertible and let $\mB'$ be obtained from
$\mB$ by replacing its column $r$ with a vector $\va\in\R^m$. Put
$\vw=\mB^{-1}\va$. Then

1. $\mB'$ is invertible if and only if $w_r\ne0$, and
   $\det\mB'=w_r\det\mB$;
2. if $w_r\ne0$, then $\mB'^{-1}=\mE\mB^{-1}$, with $\mE$ the pivot matrix
   {eq}`eq-01-01-eta-matrix` built from $\vw$ and $r$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Since $\mB\vw=\va$ and $\mB\ve_i=\mB_{\cdot i}$, the matrix $\mB'$ equals
$\mB\mF$, where $\mF$ is the identity with column $r$ replaced by $\vw$.
Expanding $\det\mF$ along row $r$ gives $\det\mF=w_r$ ($\mF$ differs from
the identity only in column $r$), so $\det\mB'=\det\mB\cdot w_r$. This
proves (1). By {prf:ref}`prop-01-01-pivot-matrix`, $\mF^{-1}=\mE$, hence
$\mB'^{-1}=\mF^{-1}\mB^{-1}=\mE\mB^{-1}$.
:::

Applied repeatedly from a starting basis $\mB_0=\mI$ (slack columns), the
proposition gives the **product form of the inverse**,

$$
\mB_k^{-1}=\mE_k\mE_{k-1}\cdots\mE_1,
$$ (eq-01-01-product-form)

which {doc}`03-06-revised-dual-simplex` uses to organize the computations
({eq}`eq-03-06-product-form`).
Part 1 of {prf:ref}`prop-01-01-basis-exchange` is the linear-algebra reason why, in the simplex method, a variable
can leave the basis only from a row where the entering column has a
nonzero entry: otherwise the new "basis" matrix would be singular. The
ratio test {eq}`eq-03-03-ratio-test` of {doc}`03-03-simplex-method` chooses, among these rows, one that
also keeps $\mB^{-1}\vb\ge\vzero$.

:::{prf:example} Two pivots and the product form of the inverse
:label: ex-01-01-pivot-product-form

Take the system $\mA\vx=\vb$ with

$$
[\mA\mid\vb]=\left(\begin{array}{ccccc|c}
1&2&1&0&0&8\\
3&1&0&1&0&9\\
1&3&0&0&1&12
\end{array}\right).
$$

The columns of $x_3,x_4,x_5$ form the identity, so they are a starting
basis $\mB_0=\mI$ with basic solution $(x_3,x_4,x_5)=(8,9,12)$.

Pivot $x_1$ into the basis in row 2, then $x_2$ into the basis in row 1.
Find the eta matrices $\mE_1$, $\mE_2$ and the inverse $\mB_2^{-1}=\mE_2\mE_1$
of the new basis matrix, and check the result.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Pivot 1: $x_1$ into the basis in row 2.** The pivot column is
$\vw=(1,3,1)$ and the pivot element is $3$. (The ratios $8/1$, $9/3$,
$12/1$ make row 2 the row with the smallest ratio $3$, so the right-hand
side (RHS) stays nonnegative; see the ratio test {eq}`eq-03-03-ratio-test`.)
By {eq}`eq-01-01-eta-matrix`,

$$
\mE_1=\begin{pmatrix}1&-\tfrac13&0\\0&\tfrac13&0\\0&-\tfrac13&1\end{pmatrix},
$$

and $\mE_1[\mA\mid\vb]$ is the tableau of {numref}`tab-01-01-pivot-tableaux`
(middle block): basis $(x_3,x_1,x_5)$ with values $(5,3,9)$.

**Pivot 2: $x_2$ into the basis in row 1.** The current column of $x_2$ is
$\vw=(\tfrac53,\tfrac13,\tfrac83)$; the ratios $5/\tfrac53=3$,
$3/\tfrac13=9$, $9/\tfrac83=\tfrac{27}{8}$ select row 1, pivot
$\tfrac53$. Then

$$
\begin{aligned}
\mE_2&=\begin{pmatrix}\tfrac35&0&0\\-\tfrac15&1&0\\-\tfrac85&0&1\end{pmatrix},\\
\mB_2^{-1}&=\mE_2\mE_1=\begin{pmatrix}\tfrac35&-\tfrac15&0\\-\tfrac15&\tfrac25&0\\-\tfrac85&\tfrac15&1\end{pmatrix}.
\end{aligned}
$$

```{list-table} The tableau $[\mA\mid\vb]$ before and after each pivot. Pivot elements in bold.
:name: tab-01-01-pivot-tableaux
:header-rows: 1

* - Basic
  - $x_1$
  - $x_2$
  - $x_3$
  - $x_4$
  - $x_5$
  - RHS
* - $x_3$
  - 1
  - 2
  - 1
  - 0
  - 0
  - 8
* - $x_4$
  - **3**
  - 1
  - 0
  - 1
  - 0
  - 9
* - $x_5$
  - 1
  - 3
  - 0
  - 0
  - 1
  - 12
* - $x_3$
  - 0
  - $\mathbf{\tfrac53}$
  - 1
  - $-\tfrac13$
  - 0
  - 5
* - $x_1$
  - 1
  - $\tfrac13$
  - 0
  - $\tfrac13$
  - 0
  - 3
* - $x_5$
  - 0
  - $\tfrac83$
  - 0
  - $-\tfrac13$
  - 1
  - 9
* - $x_2$
  - 0
  - 1
  - $\tfrac35$
  - $-\tfrac15$
  - 0
  - 3
* - $x_1$
  - 1
  - 0
  - $-\tfrac15$
  - $\tfrac25$
  - 0
  - 2
* - $x_5$
  - 0
  - 0
  - $-\tfrac85$
  - $\tfrac15$
  - 1
  - 1
```

**Checks.** The new basis matrix collects the original columns of
$x_2,x_1,x_5$ in the order of the rows,

$$
\mB_2=\begin{pmatrix}2&1&0\\1&3&0\\3&1&1\end{pmatrix},
$$

and indeed $\mB_2^{-1}\mB_2=\mI$ and $\mB_2^{-1}\vb=(3,2,1)$, the values of
$x_2,x_1,x_5$ in the last block. By {prf:ref}`prop-01-01-basis-exchange`,
$\det\mB_2=\det\mI\cdot3\cdot\tfrac53=5$.
:::

As the solution of {prf:ref}`ex-01-01-pivot-product-form` shows,
$\mB_2^{-1}$ appears in the last block of {numref}`tab-01-01-pivot-tableaux`
under the columns of $x_3,x_4,x_5$, which were the identity at the start. Every pivot multiplies the whole tableau by its
$\mE$, so the columns that started as $\mI$ end as
$\mE_2\mE_1\mI=\mB_2^{-1}$. Reading $\mB^{-1}$ (and with it the simplex
multipliers) off the slack columns is a device used throughout
{doc}`03-05-lp-duality` ({ref}`sec-03-05-reading-duals`) and
{doc}`03-07-sensitivity-analysis`.

(sec-01-01-four-subspaces)=
## The four fundamental subspaces and the solutions of $\mA\vx=\vb$

Besides $\col(\mA)$ and $\row(\mA)$, two more
subspaces belong to $\mA\in\R^{m\times n}$: the **null space**
$\nullsp(\mA)=\{\vx\in\R^n:\mA\vx=\vzero\}$ and the **left null
space** $\nullsp(\mA^\T)=\{\vy\in\R^m:\mA^\T\vy=\vzero\}$.
With $r=\rank\mA$, their dimensions are

$$
\begin{aligned}
\dim\row(\mA)&=\dim\col(\mA)=r,\\
\dim\nullsp(\mA)&=n-r,\\
\dim\nullsp(\mA^\T)&=m-r.
\end{aligned}
$$ (eq-01-01-rank-nullity)

The middle equation is the **rank–nullity theorem**: elimination leaves $r$
pivot variables and $n-r$ free variables, and each free variable, set to 1
with the other free variables at 0, gives one basis vector of the null
space. The left null space is the null space of $\mA^\T$, so the same count
applies with $m$ in place of $n$.

The subspaces come in orthogonal pairs. If $\mA\vx=\vzero$, then every row
of $\mA$ has zero inner product with $\vx$, so $\nullsp(\mA)$
is orthogonal to $\row(\mA)$ in $\R^n$. Likewise
$\nullsp(\mA^\T)$ is orthogonal to $\col(\mA)$ in
$\R^m$. Because the dimensions add up to $n$ and $m$, each pair splits its
space: every vector of $\R^n$ is a row-space vector plus a null-space
vector, and every vector of $\R^m$ is a column-space vector plus a
left-null-space vector (Strang, 2016, Sections 4.1–4.2).

:::{prf:theorem} Existence and uniqueness of solutions
:label: thm-01-01-existence-uniqueness

Let $\mA\in\R^{m\times n}$ and $\vb\in\R^m$.

1. **Existence.** $\mA\vx=\vb$ has a solution if and only if
   $\vb\in\col(\mA)$, if and only if
   $\rank[\mA\mid\vb]=\rank\mA$, if and only if $\vy^\T\vb=0$ for every
   $\vy\in\nullsp(\mA^\T)$.
2. **Structure.** If $\vx_p$ is one solution, the set of all solutions is
   $\vx_p+\nullsp(\mA)=\{\vx_p+\vz:\mA\vz=\vzero\}$.
3. **Uniqueness.** A solution, if one exists, is unique if and only if
   $\nullsp(\mA)=\{\vzero\}$, that is, $\rank\mA=n$.
   Otherwise there are infinitely many solutions.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

(1) $\mA\vx$ is a combination of the columns, so a solution exists exactly
when $\vb$ is one. Appending $\vb$ to $\mA$ does not raise the rank exactly
when $\vb$ depends on the columns. The last condition holds because
$\col(\mA)$ consists of the vectors orthogonal to
$\nullsp(\mA^\T)$ (the orthogonal pair above). (2) If
$\mA\vx_p=\vb$ and $\mA\vx=\vb$, then $\mA(\vx-\vx_p)=\vzero$; conversely
$\mA(\vx_p+\vz)=\vb+\vzero$. (3) follows from (2); a nonzero $\vz$ in the
null space gives the infinitely many solutions $\vx_p+t\vz$, $t\in\R$.
:::

Part 1, written as "either $\mA\vx=\vb$ has a solution, or some $\vy$ has
$\mA^\T\vy=\vzero$ and $\vy^\T\vb\ne0$, but not both", is the
equation version of **Farkas' lemma**; the version with inequalities,
which needs $\vx\ge\vzero$, is proved in {doc}`03-05-lp-duality`
({prf:ref}`thm-03-05-farkas`). Part 2 is
why a linear program in standard form has a whole affine set of solutions
of $\mA\vx=\vb$ to choose from when $n>m$.

For a **square** matrix the cases collapse to two.

:::{prf:theorem} Invertible matrices
:label: thm-01-01-invertible

For $\mA\in\R^{n\times n}$ the following are equivalent:

1. $\mA$ is **invertible** (nonsingular): there is a matrix $\mA^{-1}$ with
   $\mA\mA^{-1}=\mA^{-1}\mA=\mI$;
2. $\mA\vx=\vb$ has exactly one solution for every $\vb\in\R^n$;
3. $\mA\vx=\vzero$ only for $\vx=\vzero$;
4. the columns (equivalently, the rows) of $\mA$ are linearly independent;
5. $\rank\mA=n$, that is, elimination finds $n$ nonzero pivots;
6. $\det\mA\ne0$;
7. $0$ is not an eigenvalue of $\mA$.
:::

A square matrix that is not invertible is **singular**. For a singular
$\mA$, $\rank\mA<n$, the null space contains a nonzero vector, and
$\mA\vx=\vb$ has **either no solution or infinitely many**: none when
$\vb\notin\col(\mA)$, infinitely many otherwise. Which case
occurs is decided by the right-hand side, not by $\mA$ alone.

:::{prf:example} The four subspaces of a rank-2 matrix
:label: ex-01-01-four-subspaces

Let

$$
\mA=\begin{pmatrix}1&2&0&1\\2&4&1&4\\3&6&1&5\end{pmatrix}\in\R^{3\times4}.
$$

Find a basis and the dimension of each of its four fundamental subspaces.
For which $\vb$ is $\mA\vx=\vb$ solvable? Solve it for $\vb=(1,3,4)$ and for
$\vb=(1,0,0)$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Subtracting 2 times row 1 from row 2 and 3 times row 1 from row 3 gives the
rows $(0,0,1,2)$ and $(0,0,1,2)$; subtracting the new row 2 from row 3
gives zero. The RREF is

$$
\begin{pmatrix}1&2&0&1\\0&0&1&2\\0&0&0&0\end{pmatrix},
$$

with pivots in columns 1 and 3, so $\rank\mA=2$; $x_2$ and $x_4$ are free.

- **Row space** (dimension 2): spanned by $(1,2,0,1)$ and $(0,0,1,2)$, the
  nonzero rows of the RREF.
- **Column space** (dimension 2): spanned by the pivot columns of $\mA$
  itself, $(1,2,3)$ and $(0,1,1)$. Column 2 is twice column 1, and column 4
  is column 1 plus twice column 3.
- **Null space** (dimension $4-2=2$): the RREF says $x_1=-2x_2-x_4$ and
  $x_3=-2x_4$. Setting $(x_2,x_4)=(1,0)$ and $(0,1)$ gives the basis
  $(-2,1,0,0)$ and $(-1,0,-2,1)$. Both are orthogonal to both row-space
  basis vectors, as they must be.
- **Left null space** (dimension $3-2=1$): row 3 minus row 2 minus row 1 of
  $\mA$ is zero, so $\vy=(1,1,-1)$ satisfies $\mA^\T\vy=\vzero$.

By {prf:ref}`thm-01-01-existence-uniqueness`, $\mA\vx=\vb$ is solvable exactly
when $\vy^\T\vb=b_1+b_2-b_3=0$. For $\vb=(1,3,4)$ the condition holds;
setting the free variables to zero gives the particular solution
$\vx_p=(1,0,1,0)$ (check: $\mA\vx_p=(1,\,2+1,\,3+1)=\vb$), and the general
solution is $\vx_p+\alpha(-2,1,0,0)+\beta(-1,0,-2,1)$ with
$\alpha,\beta\in\R$. For $\vb=(1,0,0)$, $b_1+b_2-b_3=1\ne0$, and there is
no solution.
:::

(sec-01-01-inverse)=
## The inverse matrix

When $\mA$ is invertible, $\mA^{-1}$ is unique, and $\vx=\mA^{-1}\vb$ is
the solution of $\mA\vx=\vb$. In computation one rarely forms $\mA^{-1}$
just to solve one system; elimination is cheaper and more accurate. The
inverse matters as an object, though: $\mB^{-1}$ appears in every formula
of the simplex method in matrix form. Useful rules, for invertible
matrices of the same order:

$$
\begin{aligned}
(\mA\mB)^{-1}&=\mB^{-1}\mA^{-1},\\
(\mA^\T)^{-1}&=(\mA^{-1})^\T=:\mA^{-\T},\\
(\alpha\mA)^{-1}&=\alpha^{-1}\mA^{-1}\ (\alpha\ne0).
\end{aligned}
$$ (eq-01-01-inverse-rules)

For a $2\times2$ matrix with $ad-bc\ne0$,

$$
\begin{pmatrix}a&b\\c&d\end{pmatrix}^{-1}
=\frac{1}{ad-bc}\begin{pmatrix}d&-b\\-c&a\end{pmatrix}.
$$ (eq-01-01-inverse-2x2)

In general, $\mA^{-1}$ is found by **Gauss–Jordan elimination on
$[\mA\mid\mI]$**: the row operations that turn $\mA$ into $\mI$ form a
matrix $\mE$ with $\mE\mA=\mI$, so $\mE=\mA^{-1}$, and the same operations
turn the right block $\mI$ into $\mE\mI=\mA^{-1}$. This is the same
observation as in {prf:ref}`ex-01-01-pivot-product-form`, where $\mB^{-1}$
appeared under the slack columns. For the matrix of
{prf:ref}`ex-01-01-elimination-lu`, three pivots on the diagonal give

$$
\mA^{-1}=\begin{pmatrix}\tfrac32&-\tfrac12&0\\-3&\tfrac52&-\tfrac12\\1&-\tfrac32&\tfrac12\end{pmatrix},
$$

and $\mA^{-1}\vb=(\tfrac92-\tfrac72,\ -9+\tfrac{35}2-\tfrac{19}2,\
3-\tfrac{21}2+\tfrac{19}2)=(1,-1,2)$, the solution found there.

(sec-01-01-block)=
## Partitioned matrices

A matrix can be cut by horizontal and vertical lines into **blocks**
(submatrices). Partitioned matrices multiply like ordinary matrices with
blocks in place of entries, as long as the cuts are compatible (the column
cut of the left factor matches the row cut of the right factor):

$$
\begin{pmatrix}\mA_{11}&\mA_{12}\\\mA_{21}&\mA_{22}\end{pmatrix}
\begin{pmatrix}\vx_1\\\vx_2\end{pmatrix}
=\begin{pmatrix}\mA_{11}\vx_1+\mA_{12}\vx_2\\\mA_{21}\vx_1+\mA_{22}\vx_2\end{pmatrix}.
$$ (eq-01-01-block-product)

Note the order of the factors inside each block product, which matters
because blocks do not commute.

**The basis partition of a linear program.** After reordering the
columns, write $\mA=[\mB\ \ \mN]$ with an invertible $m\times m$ basis
$\mB$ and $\vx=(\vx_B,\vx_N)$, $\vc=(\vc_B,\vc_N)$. By
{eq}`eq-01-01-block-product`, $\mA\vx=\vb$ reads $\mB\vx_B+\mN\vx_N=\vb$.
Multiplying from the left by $\mB^{-1}$ gives
$\vx_B=\mB^{-1}\vb-\mB^{-1}\mN\vx_N$, and substituting into
$z=\vc_B^\T\vx_B+\vc_N^\T\vx_N$ gives

$$
z=\vc_B^\T\mB^{-1}\vb+\bigl(\vc_N^\T-\vc_B^\T\mB^{-1}\mN\bigr)\vx_N .
$$ (eq-01-01-basis-partition)

These two block manipulations are the whole algebra behind the canonical
form, the reduced costs and the simplex multipliers
$\vpi^\T=\vc_B^\T\mB^{-1}$, which {doc}`03-03-simplex-method` develops.

**Block inverses.** A block-triangular matrix with invertible diagonal
blocks is inverted block by block. For a square matrix partitioned as

$$
\mM=\begin{pmatrix}\mA_{11}&\mA_{12}\\\mA_{21}&\mA_{22}\end{pmatrix}
$$

with square diagonal blocks and $\mA_{11}$ invertible, the **Schur
complement** of $\mA_{11}$ is $\mS=\mA_{22}-\mA_{21}\mA_{11}^{-1}\mA_{12}$.

:::{prf:proposition} Block inverse
:label: prop-01-01-block-inverse

If $\mA_{11}$ and $\mS=\mA_{22}-\mA_{21}\mA_{11}^{-1}\mA_{12}$ are
invertible, then $\mM$ is invertible,
$\det\mM=\det\mA_{11}\cdot\det\mS$, and

$$
\mM^{-1}=\begin{pmatrix}
\mA_{11}^{-1}+\mA_{11}^{-1}\mA_{12}\mS^{-1}\mA_{21}\mA_{11}^{-1} & -\mA_{11}^{-1}\mA_{12}\mS^{-1}\\
-\mS^{-1}\mA_{21}\mA_{11}^{-1} & \mS^{-1}
\end{pmatrix}.
$$ (eq-01-01-block-inverse)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Block elimination gives the factorization

$$
\mM=\begin{pmatrix}\mI&\vzero\\\mA_{21}\mA_{11}^{-1}&\mI\end{pmatrix}
\begin{pmatrix}\mA_{11}&\vzero\\\vzero&\mS\end{pmatrix}
\begin{pmatrix}\mI&\mA_{11}^{-1}\mA_{12}\\\vzero&\mI\end{pmatrix},
$$

which is checked by multiplying out. The outer factors are block
triangular with identity diagonal blocks, so their determinants are 1 and
their inverses are obtained by changing the sign of the off-diagonal block.
Inverting the three factors in reverse order and multiplying gives
{eq}`eq-01-01-block-inverse`; taking determinants gives
$\det\mM=\det\mA_{11}\det\mS$.
:::

The most important special case in this book is the **bordered basis
matrix** of the revised simplex method and of sensitivity analysis
({eq}`eq-03-06-bordered-inverse`, {doc}`03-07-sensitivity-analysis`). With
$\mA_{11}=\mB$, $\mA_{12}=\vzero$, $\mA_{21}=\vc_B^\T$ and $\mA_{22}=1$,
the Schur complement is $\mS=1$ and {eq}`eq-01-01-block-inverse`
reduces to

$$
\begin{aligned}
\mT&=\begin{pmatrix}\mB&\vzero\\\vc_B^\T&1\end{pmatrix}^{-1}\\
&=\begin{pmatrix}\mB^{-1}&\vzero\\-\vc_B^\T\mB^{-1}&1\end{pmatrix}\\
&=\begin{pmatrix}\mB^{-1}&\vzero\\-\vpi^\T&1\end{pmatrix}.
\end{aligned}
$$ (eq-01-01-bordered-inverse)

A single matrix thus carries both $\mB^{-1}$ and the simplex multipliers.

(sec-01-01-determinants)=
## Determinants

The **determinant** $\det\mA$ of a square matrix is a number that is
nonzero exactly when $\mA$ is invertible. It is characterized by three
properties: $\det\mI=1$; exchanging two rows changes the sign; and
$\det\mA$ is linear in each row separately (when the other rows are held
fixed). Everything else follows (Strang, 2016, Section 5.1):

1. a matrix with two equal rows, or a zero row, has determinant 0;
2. adding a multiple of one row to another does not change $\det\mA$, and
   scaling one row by $\gamma$ multiplies $\det\mA$ by $\gamma$;
3. the determinant of a triangular matrix is the product of its diagonal
   entries;
4. $\det(\mA\mB)=\det\mA\cdot\det\mB$ and $\det\mA^\T=\det\mA$; hence
   $\det\mA^{-1}=1/\det\mA$, and all row properties hold for columns too.

By 2 and 3, elimination computes the determinant: $\det\mA=\pm u_{11}\cdots
u_{nn}$, the product of the pivots, with the minus sign when an odd number
of row exchanges was made. This is how determinants are computed in
practice, and it proves that $\det\mA\ne0$ exactly when all $n$ pivots are
nonzero, which is condition 5 of {prf:ref}`thm-01-01-invertible`.

For small matrices, the **cofactor expansion** is convenient. Let $M_{ij}$
be the determinant of the $(n-1)\times(n-1)$ submatrix obtained by deleting
row $i$ and column $j$ of $\mA$ (a **minor**); the number
$C_{ij}=(-1)^{i+j}M_{ij}$ is the **cofactor** of $a_{ij}$. Expansion along
any row $i$ gives

$$
\det\mA=\sum_{j=1}^n a_{ij}C_{ij}.
$$ (eq-01-01-cofactor-expansion)

For $n=2$ this is $ad-bc$. For the matrix
$\begin{pmatrix}2&-1&0\\-1&2&-1\\0&-1&2\end{pmatrix}$, expansion along row 1
gives $2(4-1)-(-1)\bigl((-1)\cdot2-(-1)\cdot0\bigr)+0=6-2=4$. The
cofactors also give an explicit inverse,

$$
\mA^{-1}=\frac{1}{\det\mA}\,\mC^\T,
$$ (eq-01-01-adjugate)

where $\mC=(C_{ij})$ is the matrix of cofactors, and **Cramer's
rule** $x_j=\det\mA^{(j)}/\det\mA$, with $\mA^{(j)}$ equal to $\mA$ with
column $j$ replaced by $\vb$. Both are useless for large computations but
valuable in proofs. For instance, if $\mB$ has integer entries and
$\det\mB=\pm1$, then {eq}`eq-01-01-adjugate` shows that $\mB^{-1}$ has integer
entries, so $\mB^{-1}\vb$ is integral for every integral $\vb$; this is
the root of the integrality of network flow problems
({prf:ref}`thm-04-04-integrality`).

(sec-01-01-eigen)=
## Eigenvalues and eigenvectors

Most vectors change direction when multiplied by a square matrix $\mA$. The
exceptional ones are only stretched.

:::{prf:definition} Eigenvalue and eigenvector
:label: def-01-01-eigen

A number $\lambda$ (real or complex) is an **eigenvalue** of
$\mA\in\R^{n\times n}$ and $\vv\ne\vzero$ is an associated
**eigenvector** if

$$
\mA\vv=\lambda\vv .
$$ (eq-01-01-eigen)

Elsewhere in the book the eigenvalues of $\mA$ are written with their
matrix, $\lambda_1(\mA),\dots,\lambda_n(\mA)$, so that they are not
confused with Lagrange multipliers.
:::

Equation {eq}`eq-01-01-eigen` says $(\mA-\lambda\mI)\vv=\vzero$ with
$\vv\ne\vzero$, so $\mA-\lambda\mI$ must be singular:

$$
p_{\mA}(\lambda)=\det(\mA-\lambda\mI)=0 .
$$ (eq-01-01-characteristic)

The **characteristic polynomial** $p_{\mA}$ has degree $n$, so $\mA$ has
exactly $n$ eigenvalues counted with multiplicity, possibly complex. For
each eigenvalue, the eigenvectors (with $\vzero$ added) form the null space
of $\mA-\lambda\mI$. Comparing coefficients of $p_{\mA}$ with its
factorization $\prod_i(\lambda_i(\mA)-\lambda)$ gives two checks:

$$
\begin{aligned}
\sum_{i=1}^n\lambda_i(\mA)&=\tr\mA,\\
\prod_{i=1}^n\lambda_i(\mA)&=\det\mA .
\end{aligned}
$$ (eq-01-01-trace-det)

In particular, $\mA$ is singular if and only if $0$ is an eigenvalue. The
eigenvalues of a triangular matrix are its diagonal entries.

**Diagonalization.** Suppose $\mA$ has $n$ linearly independent
eigenvectors $\vv_1,\dots,\vv_n$ with eigenvalues $\lambda_1,\dots,\lambda_n$.
Put them as columns of $\mV$ and the eigenvalues into
$\mD=\diag(\lambda_1,\dots,\lambda_n)$. The $n$ equations
$\mA\vv_i=\lambda_i\vv_i$ read $\mA\mV=\mV\mD$, and $\mV$ is invertible, so

$$
\begin{aligned}
\mA&=\mV\mD\mV^{-1},\\
\mA^k&=\mV\mD^k\mV^{-1}\quad (k=1,2,\dots).
\end{aligned}
$$ (eq-01-01-diagonalisation)

In the coordinates $\vy=\mV^{-1}\vx$, the matrix acts by scaling each
coordinate separately. Eigenvectors belonging to distinct eigenvalues are
independent, so a matrix with $n$ distinct eigenvalues is diagonalizable.
Not every matrix is: $\begin{pmatrix}1&1\\0&1\end{pmatrix}$ has the double
eigenvalue $1$, but $\mA-\mI$ has rank 1, so there is only one independent
eigenvector. Symmetric matrices never fail in this way, which is the
content of the next section.

(sec-01-01-symmetric)=
## Symmetric matrices and the spectral theorem

Hessians of twice continuously differentiable functions are symmetric
({doc}`01-02-matrix-calculus`), and so are the matrices of quadratic forms. For
them, the eigenvalue theory is as good as it can be.

:::{prf:theorem} Spectral theorem
:label: thm-01-01-spectral

Let $\mA\in\R^{n\times n}$ be symmetric. Then all eigenvalues of $\mA$ are
real, and there is an orthonormal basis of $\R^n$ consisting of
eigenvectors of $\mA$. Equivalently, there is an orthogonal matrix $\mV$
($\mV^\T\mV=\mI$) and a real diagonal $\mD$ with

$$
\mA=\mV\mD\mV^\T=\sum_{i=1}^n\lambda_i(\mA)\,\vv_i\vv_i^\T .
$$ (eq-01-01-spectral)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

We prove the two facts that carry the theorem; the existence of a full
orthonormal eigenbasis when eigenvalues repeat is proved by induction on
$n$ (Strang, 2016, Section 6.4). *Real eigenvalues:* let
$\mA\vv=\lambda\vv$ with complex $\vv\ne\vzero$, and let $\bar{\vv}$ be the
complex conjugate. Then $\bar{\vv}^\T\mA\vv=\lambda\,\bar{\vv}^\T\vv$.
Conjugating and transposing, using $\bar{\mA}=\mA=\mA^\T$, gives
$\bar{\vv}^\T\mA\vv=\bar\lambda\,\bar{\vv}^\T\vv$. Since
$\bar{\vv}^\T\vv=\sum\abs{v_j}^2>0$, $\lambda=\bar\lambda$ is real.
*Orthogonal eigenvectors:* if $\mA\vv_1=\lambda_1\vv_1$ and
$\mA\vv_2=\lambda_2\vv_2$ with $\lambda_1\ne\lambda_2$, then
$\lambda_1\vv_2^\T\vv_1=\vv_2^\T\mA\vv_1=(\mA\vv_2)^\T\vv_1=\lambda_2\vv_2^\T\vv_1$,
so $(\lambda_1-\lambda_2)\vv_2^\T\vv_1=0$ and $\vv_1\perp\vv_2$.
:::

The sum in {eq}`eq-01-01-spectral` writes $\mA$ as a combination of the
rank-one projections $\vv_i\vv_i^\T$ onto its eigenvector directions. Each
eigenvalue is a stretch factor along an axis, and the axes are
perpendicular: an invertible symmetric matrix maps the unit sphere to an
ellipsoid ({prf:ref}`def-02-01-ellipsoid`) whose semi-axes are $\abs{\lambda_i(\mA)}$ along the eigenvectors. The
eigenvalues also bound the ratio between $\vx^\T\mA\vx$ and $\vx^\T\vx$.

:::{prf:proposition} Rayleigh quotient bounds
:label: prop-01-01-rayleigh

For symmetric $\mA$ and every $\vx\ne\vzero$,

$$
\lambda_{\min}(\mA)\le\frac{\vx^\T\mA\vx}{\vx^\T\vx}\le\lambda_{\max}(\mA),
$$ (eq-01-01-rayleigh)

with equality on the left (right) for an eigenvector of $\lambda_{\min}(\mA)$
($\lambda_{\max}(\mA)$).
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

With $\vy=\mV^\T\vx$, {eq}`eq-01-01-spectral` gives
$\vx^\T\mA\vx=\vy^\T\mD\vy=\sum_i\lambda_i(\mA)y_i^2$, and
$\vx^\T\vx=\vy^\T\mV^\T\mV\vy=\vy^\T\vy=\sum_iy_i^2$ because $\mV$ is
orthogonal. A weighted average of the $\lambda_i(\mA)$ with weights
$y_i^2/\sum_ky_k^2$ lies between the smallest and the largest.
:::

:::{prf:example} Spectral decomposition of a $2\times2$ matrix
:label: ex-01-01-spectral-2x2

Let $\mA=\begin{pmatrix}2&1\\1&2\end{pmatrix}$.

Find its eigenvalues and orthonormal eigenvectors, write its spectral
decomposition, and use it to compute $\mA^5$, to write $\vx^\T\mA\vx$ as a
weighted sum of squares and to describe the image of the unit circle.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Its characteristic
polynomial is $(2-\lambda)^2-1=(\lambda-1)(\lambda-3)$, so
$\lambda_1(\mA)=3$ and $\lambda_2(\mA)=1$ (check: $3+1=4=\tr\mA$,
$3\cdot1=3=\det\mA$). The null space of $\mA-3\mI$ is spanned by
$(1,1)$, that of $\mA-\mI$ by $(1,-1)$; these are orthogonal, as the
spectral theorem promises. Normalized, $\vv_1=(1,1)/\sqrt2$ for eigenvalue
3 and $\vv_2=(1,-1)/\sqrt2$ for eigenvalue 1, and

$$
\mA=3\vv_1\vv_1^\T+1\,\vv_2\vv_2^\T
=\begin{pmatrix}1.5&1.5\\1.5&1.5\end{pmatrix}
+\begin{pmatrix}0.5&-0.5\\-0.5&0.5\end{pmatrix}.
$$

Powers are now immediate: $\mA^5=3^5\vv_1\vv_1^\T+1^5\vv_2\vv_2^\T
=\begin{pmatrix}122&121\\121&122\end{pmatrix}$. The quadratic form becomes
$\vx^\T\mA\vx=3y_1^2+y_2^2$ with $y_i=\vv_i^\T\vx$, a sum of squares with
positive weights. {numref}`fig-01-01-circle-ellipse` shows the geometry: the
unit circle is mapped to an ellipse with semi-axes 3 and 1 along $\vv_1$
and $\vv_2$. The vector $\vx=(1,0)$ is not an eigenvector; its image
$\mA\vx=(2,1)$ is turned toward $\vv_1$, and its Rayleigh quotient
$\vx^\T\mA\vx/\vx^\T\vx=2$ lies between 1 and 3.

```{figure} figures/01-01-circle-ellipse.svg
:name: fig-01-01-circle-ellipse
:width: 480px
:class: dark-light
:alt: The unit circle and its image under the matrix with rows 2 1 and 1 2, an ellipse tilted at 45 degrees; the eigenvector v1 is stretched by 3 along the long axis, the eigenvector v2 is unchanged on the short axis, and x equal to e1 is mapped to Ax equal to (2, 1)

The symmetric matrix $\mA=\begin{pmatrix}2&1\\1&2\end{pmatrix}$ maps the
unit circle (gray) to an ellipse (blue). The eigenvectors are only
stretched: $\mA\vv_1=3\vv_1$ along the long axis, $\mA\vv_2=\vv_2$ along
the short axis. A generic vector $\vx$ is also turned: $\mA\vx=(2,1)$.
```
:::

(sec-01-01-quadratic-forms)=
## Quadratic forms and definiteness

A **quadratic form** is a function $Q(\vx)=\vx^\T\mA\vx=\sum_i\sum_j
a_{ij}x_ix_j$ with a square matrix $\mA$. Only the symmetric part of $\mA$
matters: since $\vx^\T\mA\vx$ is a number, it equals its transpose
$\vx^\T\mA^\T\vx$, and so

$$
\vx^\T\mA\vx=\vx^\T\Bigl(\tfrac12(\mA+\mA^\T)\Bigr)\vx .
$$ (eq-01-01-symmetric-part)

For example, $2x_1^2+3x_1x_2+x_2^2$ is $\vx^\T\mA\vx$ for
$\mA=\begin{pmatrix}2&3\\0&1\end{pmatrix}$ and for
$\begin{pmatrix}2&1.5\\1.5&1\end{pmatrix}$ alike. From now on quadratic
forms are written with a symmetric matrix, as is customary.

:::{prf:definition} Definiteness
:label: def-01-01-definiteness

A symmetric matrix $\mA\in\Sym^n$ (the set of symmetric $n\times n$
matrices) and its quadratic form are

- **positive definite**, $\mA\succ0$, if $\vx^\T\mA\vx>0$ for every
  $\vx\ne\vzero$;
- **positive semidefinite**, $\mA\succeq0$, if $\vx^\T\mA\vx\ge0$ for every
  $\vx$;
- **negative definite** (semidefinite) if $-\mA$ is positive definite
  (semidefinite);
- **indefinite** if $\vx^\T\mA\vx$ takes both positive and negative values.

$\Sym^n_{++}$ and $\Sym^n_+$ denote the sets of positive definite and
positive semidefinite matrices of order $n$.
:::

These classes decide the shape of a quadratic function, and through the
Hessian and Taylor's theorem ({doc}`01-02-matrix-calculus`) they decide whether
a stationary point is a minimum, a maximum or a saddle
({prf:ref}`thm-06-02-second-order-sufficient`) and whether a function is
convex ({prf:ref}`thm-02-02-second-order`). {numref}`fig-01-01-quadratic-forms` shows the
three typical pictures in $\R^2$.

```{figure} figures/01-01-quadratic-forms.svg
:name: fig-01-01-quadratic-forms
:width: 720px
:class: dark-light
:alt: Three contour plots. Left, nested ellipses around the origin. Middle, parallel straight lines with a thick green line through the origin where the form is zero. Right, hyperbolas, blue on one pair of opposite sides and dashed orange on the other, separated by two green lines through the origin

Contours of $\vx^\T\mA\vx$ (levels 0.5 to 8 in blue, negative levels in
dashed orange, the zero level in green). (a)
$\mA=\begin{pmatrix}2&1\\1&2\end{pmatrix}\succ0$: closed ellipses, the
only zero is $\vx=\vzero$. (b) $\mA=\begin{pmatrix}1&-1\\-1&1\end{pmatrix}
\succeq0$, singular: $\vx^\T\mA\vx=(x_1-x_2)^2$ vanishes on the whole line
$x_1=x_2$, the null space of $\mA$. (c)
$\mA=\begin{pmatrix}1&2\\2&1\end{pmatrix}$, indefinite: a saddle, positive
in one pair of sectors and negative in the other.
```

:::{prf:theorem} Eigenvalue test
:label: thm-01-01-eigen-test

A symmetric matrix $\mA$ is positive definite if and only if all its
eigenvalues are positive, and positive semidefinite if and only if all are
nonnegative. It is indefinite if and only if it has a positive and a
negative eigenvalue.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

By the spectral theorem, $\vx^\T\mA\vx=\sum_i\lambda_i(\mA)y_i^2$ with
$\vy=\mV^\T\vx$, and $\vy\ne\vzero$ exactly when $\vx\ne\vzero$. If all
$\lambda_i(\mA)>0$ the sum is positive for $\vy\ne\vzero$. If some
$\lambda_k(\mA)\le0$, the eigenvector $\vx=\vv_k$ gives
$\vx^\T\mA\vx=\lambda_k(\mA)\le0$. The semidefinite and indefinite cases
are proved in the same way.
:::

Computing eigenvalues by hand is tedious beyond $n=2$, and determinants of
submatrices give an alternative. A **principal submatrix** of $\mA$ keeps
the rows and columns with the same index set $J\subseteq\{1,\dots,n\}$; its
determinant is a **principal minor**. The **leading** principal submatrices
keep $J=\{1,\dots,k\}$, and their determinants are the **leading principal
minors** $D_1,\dots,D_n$.

:::{prf:theorem} Sylvester's criterion
:label: thm-01-01-sylvester

A **symmetric** matrix $\mA$ is positive definite if and only if all its
leading principal minors are positive, $D_k>0$ for $k=1,\dots,n$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

*Necessity.* If $\mA\succ0$, each leading submatrix $\mA_k$ is positive
definite too (take $\vx$ with zeros after position $k$), so its eigenvalues
are positive and $D_k=\det\mA_k$, their product, is positive.
*Sufficiency.* Run Gaussian elimination without row exchanges. The first
$k-1$ steps act inside the leading $k\times k$ block as row additions and do not change
its determinant, so after them the block is upper triangular with pivots
$d_1,\dots,d_k$ and $D_k=d_1d_2\cdots d_k$. If all $D_k>0$, every pivot
$d_k=D_k/D_{k-1}$ (with $D_0=1$) is positive, so elimination never stops,
and symmetry gives $\mA=\mL\mD\mL^\T$ with $\mL$ unit lower triangular and
$\mD=\diag(d_1,\dots,d_n)$. Then $\vx^\T\mA\vx=\sum_kd_k(\mL^\T\vx)_k^2>0$
for $\vx\ne\vzero$, since $\mL^\T$ is invertible.
:::

For **semidefiniteness**, leading minors are not enough; all principal
minors must be checked.

:::{prf:theorem} Principal-minor test for semidefiniteness
:label: thm-01-01-psd-minors

A symmetric matrix $\mA$ is positive semidefinite if and only if **all**
its principal minors (for every index set $J$, not only the leading ones)
are nonnegative. It is negative definite if and only if
$(-1)^kD_k>0$ for $k=1,\dots,n$, that is, the leading minors alternate in
sign starting with $D_1<0$.
:::

The negative-definite statement is Sylvester's criterion applied to $-\mA$
(the minors of $-\mA$ of order $k$ are $(-1)^k$ times those of $\mA$). For
the semidefinite statement see Strang (2016, Section 6.5) or Magnus and
Neudecker (2019, Chapter 1). The diagonal matrix $\diag(0,-1)$ shows why
the extra minors are needed: its leading minors are $D_1=0$ and $D_2=0$,
both nonnegative, but $\ve_2^\T\diag(0,-1)\ve_2=-1<0$. The principal minor
for $J=\{2\}$, the entry $-1$, reveals the problem.

The factorization in the proof of {prf:ref}`thm-01-01-sylvester` is also the
best practical test. When $\mA\succ0$, the pivots are positive and
$\mA=\mL\mD\mL^\T=(\mL\mD^{1/2})(\mL\mD^{1/2})^\T$. This is the
**Cholesky factorization** $\mA=\mG\mG^\T$ with $\mG$
lower triangular with positive diagonal. It exists if and only if
$\mA\succ0$, it costs about $n^3/3$ arithmetic operations (half the work
of Gaussian elimination and far less than an eigenvalue computation), and an attempt that meets a nonpositive pivot proves that
$\mA$ is not positive definite.

:::{prf:example} Leading principal minors of a $3\times3$ matrix
:label: ex-01-01-leading-minors

For a symmetric $3\times3$ matrix $\mA=(a_{ij})$ the three leading
principal minors are

$$
\begin{aligned}
D_1&=a_{11},\\
D_2&=\det\begin{pmatrix}a_{11}&a_{12}\\a_{21}&a_{22}\end{pmatrix}=a_{11}a_{22}-a_{12}a_{21},\\
D_3&=\det\mA .
\end{aligned}
$$

By {prf:ref}`thm-01-01-sylvester`, $\mA\succ0$ if and only if $D_1>0$,
$D_2>0$ and $D_3>0$.

Use them to decide whether
$\mA=\begin{pmatrix}2&-1&0\\-1&2&-1\\0&-1&2\end{pmatrix}$ and
$\mA=\begin{pmatrix}1&2&0\\2&1&0\\0&0&1\end{pmatrix}$ are positive definite
or positive semidefinite.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**A positive definite instance.** For
$\mA=\begin{pmatrix}2&-1&0\\-1&2&-1\\0&-1&2\end{pmatrix}$,
$D_1=2$, $D_2=4-1=3$ and, by the cofactor expansion of
{ref}`sec-01-01-determinants`, $D_3=4$. All are positive, so $\mA\succ0$. The
eigenvalue test agrees: the eigenvalues are $2-\sqrt2\approx0.5858$, $2$
and $2+\sqrt2\approx3.414$, with sum $6=\tr\mA$ and product $4=D_3$. The
elimination pivots are $D_1=2$, $D_2/D_1=\tfrac32$ and $D_3/D_2=\tfrac43$,
and the Cholesky factor is

$$
\begin{aligned}
\mG&=\begin{pmatrix}\sqrt2&0&0\\-\tfrac{1}{\sqrt2}&\sqrt{\tfrac32}&0\\0&-\sqrt{\tfrac23}&\sqrt{\tfrac43}\end{pmatrix}\\
&\approx\begin{pmatrix}1.414&0&0\\-0.7071&1.225&0\\0&-0.8165&1.155\end{pmatrix},
\end{aligned}
$$

whose squared diagonal entries are the pivots.

**An indefinite instance.** For
$\mA=\begin{pmatrix}1&2&0\\2&1&0\\0&0&1\end{pmatrix}$, $D_1=1>0$ but
$D_2=1-4=-3<0$, so $\mA$ is not positive definite; since a negative
$2\times2$ leading minor is impossible for a semidefinite matrix as well,
$\mA$ is not positive semidefinite either. The eigenvalues $-1$, $1$, $3$
confirm that $\mA$ is indefinite: $\vx=(1,-1,0)$ gives
$\vx^\T\mA\vx=1-4+1=-2$.
:::

The minor test is valid only for symmetric matrices, and the condition
"all leading principal minors $\ge0$" does not imply positive
semidefiniteness. For example, $\begin{pmatrix}1&-3\\0&1\end{pmatrix}$ has
leading minors $1$ and $1$, yet $\vx=(1,1)$ gives $\vx^\T\mA\vx=-1$; and
$\diag(0,-1)$ has leading minors $0,0$ but is not semidefinite.

(sec-01-01-norms)=
## Norms

A **norm** on $\R^n$ is a function $\norm{\cdot}$ with, for all
$\vx,\vy\in\R^n$ and $\alpha\in\R$: $\norm{\vx}\ge0$ with equality only for
$\vx=\vzero$; $\norm{\alpha\vx}=\abs{\alpha}\norm{\vx}$; and the
**triangle inequality** $\norm{\vx+\vy}\le\norm{\vx}+\norm{\vy}$. Three
norms are used in this book:

$$
\begin{aligned}
\norm{\vx}_1&=\sum_{j=1}^n\abs{x_j},\\
\norm{\vx}_2&=\Bigl(\sum_{j=1}^nx_j^2\Bigr)^{1/2}=\sqrt{\vx^\T\vx},\\
\norm{\vx}_\infty&=\max_{j}\abs{x_j}.
\end{aligned}
$$ (eq-01-01-vector-norms)

The triangle inequality for $\norm{\cdot}_2$ follows from Cauchy–Schwarz:
$\norm{\vx+\vy}_2^2=\norm{\vx}^2+2\vx^\T\vy+\norm{\vy}^2\le
(\norm{\vx}+\norm{\vy})^2$. For $\vx=(3,-4,0)$, $\norm{\vx}_1=7$,
$\norm{\vx}_2=5$ and $\norm{\vx}_\infty=4$. In $\R^n$ all norms are
equivalent up to constants,

$$
\norm{\vx}_\infty\le\norm{\vx}_2\le\norm{\vx}_1\le\sqrt n\,\norm{\vx}_2\le n\norm{\vx}_\infty ,
$$ (eq-01-01-norm-equivalence)

so convergence of a sequence ({prf:ref}`def-02-01-limit-point`) does not
depend on the norm used. {numref}`fig-01-01-unit-balls` shows their unit balls.
The 1- and $\infty$-norm balls are polyhedra, which is why problems with
these norms can be written as linear programs
({ref}`sec-07-02-abs`).

```{figure} figures/01-01-unit-balls.svg
:name: fig-01-01-unit-balls
:width: 400px
:class: dark-light
:alt: A diamond, a circle and a square centered at the origin, the unit balls of the 1-norm, 2-norm and infinity-norm, nested inside each other

Boundaries of the unit balls $\{\vx:\norm{\vx}\le1\}$ in $\R^2$: the
diamond of $\norm{\cdot}_1$ lies inside the circle of $\norm{\cdot}_2$,
which lies inside the square of $\norm{\cdot}_\infty$, as
{eq}`eq-01-01-norm-equivalence` says.
```

**Induced matrix norms.** A vector norm induces a norm on matrices, the
largest stretch factor:

$$
\norm{\mA}=\max_{\vx\ne\vzero}\frac{\norm{\mA\vx}}{\norm{\vx}} .
$$ (eq-01-01-induced-norm)

By definition $\norm{\mA\vx}\le\norm{\mA}\norm{\vx}$, and
$\norm{\mA\mB}\le\norm{\mA}\norm{\mB}$. The three induced norms have
closed forms:

$$
\begin{aligned}
\norm{\mA}_1&=\max_j\sum_i\abs{a_{ij}},\\
\norm{\mA}_\infty&=\max_i\sum_j\abs{a_{ij}},\\
\norm{\mA}_2&=\sqrt{\lambda_{\max}(\mA^\T\mA)} ,
\end{aligned}
$$ (eq-01-01-matrix-norms)

the largest absolute column sum, the largest absolute row sum, and the
**spectral norm**. The last follows from {prf:ref}`prop-01-01-rayleigh`:
$\norm{\mA\vx}_2^2/\norm{\vx}_2^2=\vx^\T(\mA^\T\mA)\vx/\vx^\T\vx$, whose
maximum is $\lambda_{\max}(\mA^\T\mA)$. For symmetric $\mA$,
$\norm{\mA}_2=\max_i\abs{\lambda_i(\mA)}$. A norm that is not induced but
often convenient is the **Frobenius norm**
$\norm{\mA}_F=\bigl(\sum_{i,j}a_{ij}^2\bigr)^{1/2}=\sqrt{\tr(\mA^\T\mA)}$.

For $\mA=\begin{pmatrix}1&-2\\3&4\end{pmatrix}$: the absolute column sums
are $4$ and $6$, so $\norm{\mA}_1=6$; the row sums are $3$ and $7$, so
$\norm{\mA}_\infty=7$; $\mA^\T\mA=\begin{pmatrix}10&10\\10&20\end{pmatrix}$
has eigenvalues $15\pm\sqrt{125}$, so $\norm{\mA}_2=\sqrt{15+\sqrt{125}}
\approx5.117$; and $\norm{\mA}_F=\sqrt{30}\approx5.477$. The **condition
number** $\kappa(\mA)=\norm{\mA}\,\norm{\mA^{-1}}$ measures how much
relative errors in $\vb$ can be amplified in the solution of $\mA\vx=\vb$;
here $\kappa_2(\mA)\approx2.618$, a well-conditioned matrix. Large condition
numbers of Hessians slow down gradient methods
({prf:ref}`thm-06-05-sd-rate`).

(sec-01-01-orthogonality)=
## Orthogonality and projections

Vectors $\vq_1,\dots,\vq_k$ are **orthonormal** if
$\vq_i^\T\vq_j=0$ for $i\ne j$ and $\norm{\vq_i}=1$. Orthonormal vectors are
independent (take the inner product of $\sum\alpha_i\vq_i=\vzero$ with
$\vq_j$ to get $\alpha_j=0$). A square matrix $\mV$ with orthonormal columns
is **orthogonal**: $\mV^\T\mV=\mI$, so $\mV^{-1}=\mV^\T$. Orthogonal
matrices preserve lengths and angles, $\norm{\mV\vx}_2=\norm{\vx}_2$ and
$(\mV\vx)^\T(\mV\vy)=\vx^\T\vy$; rotations and reflections are examples,
and so is the eigenvector matrix of the spectral theorem. The
**orthogonal complement** of a subspace $S$ is
$S^\perp=\{\vy:\vy^\T\vx=0\ \text{for all}\ \vx\in S\}$; in this language
{ref}`sec-01-01-four-subspaces` says
$\nullsp(\mA)=\row(\mA)^\perp$ and
$\nullsp(\mA^\T)=\col(\mA)^\perp$.

**Projection onto a line.** The point $\vp$ on the line through $\vzero$
in the direction $\va\ne\vzero$ that is closest to $\vb$ is found by
requiring the error to be orthogonal to the line. With $\vp=t\va$, the
condition $\va^\T(\vb-t\va)=0$ gives $t=\va^\T\vb/\va^\T\va$, so

$$
\begin{aligned}
\vp&=\frac{\va^\T\vb}{\va^\T\va}\,\va=\mP\vb,\\
\mP&=\frac{\va\va^\T}{\va^\T\va}.
\end{aligned}
$$ (eq-01-01-projection-line)

It is the closest point because, for any other point $s\va$ of the line,
$\norm{\vb-s\va}^2=\norm{\vb-\vp}^2+\norm{\vp-s\va}^2$ by Pythagoras (the
two vectors on the right are orthogonal).

:::{prf:proposition} Projection onto a column space
:label: prop-01-01-projection

Let $\mA\in\R^{m\times n}$ have full column rank $n$. Then $\mA^\T\mA$ is
invertible, and the point of $\col(\mA)$ closest to
$\vb\in\R^m$ is $\vp=\mA\hat{\vx}$, where $\hat{\vx}$ solves the **normal
equations**

$$
\mA^\T\mA\,\hat{\vx}=\mA^\T\vb .
$$ (eq-01-01-normal-equations)

The projection matrix $\mP=\mA(\mA^\T\mA)^{-1}\mA^\T$ is symmetric and
satisfies $\mP^2=\mP$, and the error $\vb-\vp$ lies in
$\nullsp(\mA^\T)$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

If $\mA^\T\mA\vx=\vzero$, then $\vx^\T\mA^\T\mA\vx=\norm{\mA\vx}^2=0$, so
$\mA\vx=\vzero$ and, by full column rank, $\vx=\vzero$; hence $\mA^\T\mA$
is invertible. The error $\vb-\mA\hat{\vx}$ must be orthogonal to every
column of $\mA$, that is, $\mA^\T(\vb-\mA\hat{\vx})=\vzero$, which is
{eq}`eq-01-01-normal-equations`. For any other $\mA\vx$, Pythagoras gives
$\norm{\vb-\mA\vx}^2=\norm{\vb-\mA\hat{\vx}}^2+\norm{\mA(\hat{\vx}-\vx)}^2$,
so $\mA\hat{\vx}$ is closest. $\mP^2=\mA(\mA^\T\mA)^{-1}\mA^\T\mA
(\mA^\T\mA)^{-1}\mA^\T=\mP$ and symmetry are checked directly.
:::

{doc}`01-02-matrix-calculus` derives the same normal equations by setting the
gradient of $\norm{\mA\vx-\vb}_2^2$ to zero; this is the least-squares
method of regression ({ref}`sec-06-01-regression`). The matrix $\mA^\T\mA$ is
always positive semidefinite, since $\vx^\T\mA^\T\mA\vx=\norm{\mA\vx}^2\ge0$,
and positive definite under full column rank.

:::{prf:example} Projecting a vector onto a line
:label: ex-01-01-projection-line

Project $\vb=(1,2)$ onto the line through $\va=(2,1)$: find the
projection $\vp$, the error $\vr=\vb-\vp$, the distance from $\vb$ to the
line and the projection matrix $\mP$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Here
$\va^\T\vb=4$ and $\va^\T\va=5$, so $\vp=\tfrac45\va=(1.6,0.8)$. The error
$\vr=\vb-\vp=(-0.6,1.2)$ is orthogonal to the line:
$\va^\T\vr=-1.2+1.2=0$. The distance from $\vb$ to the line is
$\norm{\vr}=\sqrt{0.36+1.44}=\sqrt{1.8}\approx1.342$. The projection matrix
is

$$
\mP=\frac{1}{5}\begin{pmatrix}4&2\\2&1\end{pmatrix}
=\begin{pmatrix}0.8&0.4\\0.4&0.2\end{pmatrix},
$$

and indeed $\mP\vb=(0.8+0.8,\ 0.4+0.4)=(1.6,0.8)$ and $\mP^2=\mP$
(projecting twice changes nothing). {numref}`fig-01-01-projection-line` shows
the right angle at $\vp$.

```{figure} figures/01-01-projection-line.svg
:name: fig-01-01-projection-line
:width: 460px
:class: dark-light
:alt: Vectors a equal to (2, 1) and b equal to (1, 2) from the origin, the projection p of b on the line through a, and the dashed error vector r from p to b meeting the line at a right angle

Orthogonal projection of $\vb=(1,2)$ onto the line spanned by
$\va=(2,1)$: $\vp=0.8\,\va=(1.6,0.8)$, and the error $\vr=\vb-\vp$ is
perpendicular to the line (right-angle mark), so $\vp$ is the closest point
of the line to $\vb$.
```
:::

Repeating the line projection on a sequence of vectors and subtracting the
components already found turns any basis into an orthonormal one
(**Gram–Schmidt process**); in matrix form this is the factorization
$\mA=\mQ\mathbf R$ used by numerically stable least-squares solvers
(Strang, 2016, Section 4.4).

## Summary and Common Mistakes

### Summary

- The product $\mA\vx$ is a combination of the columns of $\mA$; row
  operations are multiplications from the left
  ({eq}`eq-01-01-matrix-product`, {eq}`eq-01-01-elementary-matrix`).
- The inner product gives lengths, distances and angles; its sign
  distinguishes acute, right and obtuse angles ({eq}`eq-01-01-angle`,
  {prf:ref}`thm-01-01-cauchy-schwarz`).
- Vectors are independent when $\mX\valpha=\vzero$ forces $\valpha=\vzero$,
  checked by elimination ({prf:ref}`prop-01-01-independence-test`); the rank
  counts pivots and equals row rank and column rank
  ({prf:ref}`thm-01-01-rank`).
- A pivot on entry $(r,s)$ is multiplication by an eta matrix
  ({prf:ref}`prop-01-01-pivot-matrix`); exchanging one basis column needs a
  nonzero pivot and updates the inverse by $\mB'^{-1}=\mE\mB^{-1}$, which
  gives the product form {eq}`eq-01-01-product-form`
  ({prf:ref}`prop-01-01-basis-exchange`, {prf:ref}`ex-01-01-pivot-product-form`).
- $\mA\vx=\vb$ is solvable if and only if $\vb\in\col(\mA)$; solutions
  form $\vx_p+\nullsp(\mA)$; a singular square system has none
  or infinitely many ({prf:ref}`thm-01-01-existence-uniqueness`,
  {prf:ref}`thm-01-01-invertible`).
- Block elimination gives the basis partition of a linear program
  {eq}`eq-01-01-basis-partition`, the block inverse
  {eq}`eq-01-01-block-inverse` and the bordered inverse
  {eq}`eq-01-01-bordered-inverse`.
- Determinants are products of pivots and vanish exactly for singular
  matrices; cofactors give the adjugate formula {eq}`eq-01-01-adjugate`.
- Symmetric matrices have real eigenvalues and an orthonormal eigenbasis,
  $\mA=\mV\mD\mV^\T$ ({prf:ref}`thm-01-01-spectral`), and their Rayleigh
  quotients lie between the extreme eigenvalues
  ({prf:ref}`prop-01-01-rayleigh`).
- Definiteness ({prf:ref}`def-01-01-definiteness`) is tested by eigenvalues
  ({prf:ref}`thm-01-01-eigen-test`), by leading minors for positive
  definiteness ({prf:ref}`thm-01-01-sylvester`) and by all principal minors
  for semidefiniteness ({prf:ref}`thm-01-01-psd-minors`);
  {numref}`fig-01-01-quadratic-forms` shows the three shapes.
- The 1-, 2- and $\infty$-norms and their induced matrix norms
  {eq}`eq-01-01-matrix-norms`; projections solve the normal equations
  {eq}`eq-01-01-normal-equations` ({prf:ref}`prop-01-01-projection`).

### Common mistakes

- **Applying the minor test to a nonsymmetric matrix.** Sylvester's
  criterion and the eigenvalue test are statements about symmetric
  matrices. For a nonsymmetric $\mA$, test the symmetric part
  $\tfrac12(\mA+\mA^\T)$ instead ({eq}`eq-01-01-symmetric-part`):
  $\begin{pmatrix}1&-3\\0&1\end{pmatrix}$ has positive leading minors but
  is indefinite as a quadratic form.
- **Using leading minors to test semidefiniteness.** "All leading
  principal minors $\ge0$" does not imply $\mA\succeq0$; $\diag(0,-1)$ is
  a counterexample. Positive semidefiniteness needs **all** principal minors
  $\ge0$ ({prf:ref}`thm-01-01-psd-minors`), or all eigenvalues $\ge0$.
- **Getting negative definiteness wrong.** A negative definite matrix does
  not have all leading minors negative; they alternate in sign,
  $D_1<0$, $D_2>0$, $D_3<0,\dots$
- **Pivoting on a zero entry.** The pivot element must be nonzero; a zero
  pivot would make the new basis matrix singular
  ({prf:ref}`prop-01-01-basis-exchange`). In the simplex method the pivot is
  moreover positive, chosen by the ratio test.
- **Reversing the order in products of inverses and transposes.**
  $(\mA\mB)^{-1}=\mB^{-1}\mA^{-1}$ and $(\mA\mB)^\T=\mB^\T\mA^\T$, not the
  other way round; and $\mE_2\mE_1$ in {eq}`eq-01-01-product-form` lists the
  latest pivot first.
- **Concluding uniqueness from a square system.** A square system with a
  singular matrix has no solution or infinitely many, never exactly one;
  always check $\det\mA\ne0$ or the number of pivots first.
- **Treating $\vx^\T\vy$ as a matrix.** $\vx^\T\vy$ is a number, but
  $\vx\vy^\T$ is an $n\times n$ matrix of rank one.

## References

- Bazaraa, M. S., Jarvis, J. J., & Sherali, H. D. (2010). *Linear programming
  and network flows* (4th ed.). Wiley. (Sections 2.1–2.3: vectors, matrices
  and simultaneous linear equations.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Appendix
  A.1, norms; Appendix A.5, linear algebra: range and null space, symmetric
  eigenvalue decomposition, Schur complement.)
- Magnus, J. R., & Neudecker, H. (2019). *Matrix differential calculus with
  applications in statistics and econometrics* (3rd ed.). Wiley.
  <https://doi.org/10.1002/9781119541219> (Chapter 1: basic properties of
  vectors and matrices, including eigenvalues, the spectral theorem,
  definite matrices and partitioned matrices.)
- Petersen, K. B., & Pedersen, M. S. (2012). *The matrix cookbook* (Version
  November 15, 2012). Technical University of Denmark. (Sections 1, 3.1–3.2
  and 9.1: determinants, inverses, block matrices, eigenvalues, norms.)
- Strang, G. (2016). *Introduction to linear algebra* (5th ed.).
  Wellesley-Cambridge Press. (Chapters 1–2: vectors and elimination;
  Chapter 3: the four subspaces and rank; Chapter 4: orthogonality and
  projections; Chapter 5: determinants; Chapter 6: eigenvalues, symmetric
  and positive definite matrices.)
