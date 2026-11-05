# Vector and Matrix Calculus

## Overview

Nonlinear optimization is calculus in many variables. The first-order
conditions of {doc}`06-02-optimality-conditions` set a gradient to zero, the
Karush–Kuhn–Tucker (KKT) conditions of {doc}`06-03-kkt-conditions` combine the gradients of the
objective and the constraints, the second-order conditions inspect a
Hessian, and every algorithm of {doc}`06-04-line-search` to
{doc}`06-06-constrained-algorithms` builds a linear or quadratic model of the
objective from a Taylor expansion. This chapter supplies these tools. It
fixes the **layout convention** of the book (gradients are column vectors,
the rows of a Jacobian are transposed gradients), defines partial and
directional derivatives, gradients, Jacobians and Hessians, and derives the
derivatives that occur again and again: of linear and quadratic functions,
of least-squares objectives and norms ({numref}`tab-01-02-vector-derivatives`),
and of the trace, the determinant and the inverse with respect to a matrix
({numref}`tab-01-02-matrix-derivatives`). It states the chain rule in vector
form ({prf:ref}`thm-01-02-chain-rule`) and Taylor's theorem with remainder
({prf:ref}`thm-01-02-taylor`), proves that the gradient is orthogonal to level
sets and points in the direction of steepest ascent
({prf:ref}`thm-01-02-gradient-level-set`, {prf:ref}`thm-01-02-steepest-ascent`),
and closes with three applications: the normal equations of least squares,
the minimization of a quadratic function and the gradient of a Lagrangian.
It uses the matrix algebra of {doc}`01-01-linear-algebra` (inner products,
symmetric matrices, quadratic forms, norms and projections) and one-variable
calculus, including the mean value theorem.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- compute partial derivatives, gradients, Jacobians, Hessians and
  directional derivatives, using the book's layout convention;
- differentiate linear and quadratic forms, least-squares objectives and
  norms, and use the chain rule in vector form;
- differentiate the trace, the log-determinant and the inverse with respect
  to a matrix;
- write first- and second-order Taylor expansions with remainder and judge
  how well they approximate a function;
- explain why the gradient is orthogonal to level sets and why $-\nabla f$ is
  the direction of steepest descent; and
- derive the normal equations, the minimizer of a quadratic function and
  the gradient of a Lagrangian.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\partial f/\partial x_j$ | partial derivative of $f$ with respect to $x_j$ |
| $\nabla f(\vx)$, $\nabla^2 f(\vx)$ | gradient (a column vector) and Hessian of $f$ at $\vx$ |
| $f'(\vx;\vd)$ | directional derivative of $f$ at $\vx$ in direction $\vd$ |
| $D\vg(\vx)$ | Jacobian matrix of $\vg:\R^n\to\R^m$ at $\vx$ ($m\times n$) |
| $C^1$, $C^2$ | continuously differentiable, twice continuously differentiable |
| $\bar\vx$ | point at which a function is expanded |
| $L_c(f)$ | level set $\{\vx:f(\vx)=c\}$ |
| $\mX$, $\partial f/\partial\mX$ | matrix variable and the derivative of a scalar function with respect to it (same shape as $\mX$) |
| $o(\norm{\vd}^k)$ | a term that tends to zero faster than $\norm{\vd}^k$ as $\vd\to\vzero$ |
| $O(\norm{\vd}^k)$ | a term bounded by a constant times $\norm{\vd}^k$ for small $\vd$ |
| $L(\vx,\vlambda,\vmu)$ | Lagrangian of a nonlinear program |

(sec-01-02-layout)=
## Functions of several variables and the layout convention

This chapter studies **scalar functions** $f:\R^n\to\R$, such as an
objective $f(\vx)$, and **vector functions** $\vg:\R^n\to\R^m$ with
components $g_1,\dots,g_m$, such as the constraint functions of a nonlinear
program. Derivatives of vector functions are arrays, and the way they are
arranged is a convention. Texts differ, and mixing conventions is the most
common source of misplaced transposes. This book uses the following one
throughout.

:::{prf:definition} Layout convention
:label: def-01-02-layout

1. The **gradient** $\nabla f(\vx)$ of a scalar function is a **column**
   vector in $\R^n$ with components $\partial f/\partial x_j$.
2. The **Jacobian** $D\vg(\vx)$ of a vector function $\vg:\R^n\to\R^m$ is
   the $m\times n$ matrix with entries $\partial g_i/\partial x_j$: row $i$ is
   the transposed gradient of the component $g_i$,

   $$
   D\vg(\vx)=\begin{pmatrix}\nabla g_1(\vx)^\T\\\vdots\\\nabla g_m(\vx)^\T\end{pmatrix}.
   $$ (eq-01-02-jacobian-rows)

3. For a scalar function, $Df(\vx)=\nabla f(\vx)^\T$ is a row, and the
   **Hessian** is the Jacobian of the gradient, $\nabla^2 f(\vx)=D(\nabla
   f)(\vx)$, an $n\times n$ matrix.
:::

This is the "numerator layout" of Magnus and Neudecker (2019) for
Jacobians, combined with column gradients. Its advantage is that the
first-order approximations read like the one-variable formula
$g(x+d)\approx g(x)+g'(x)d$: $\vg(\vx+\vd)\approx\vg(\vx)+D\vg(\vx)\vd$ and
$f(\vx+\vd)\approx f(\vx)+\nabla f(\vx)^\T\vd$. Some texts, mainly in
statistics and machine learning, use the transposed Jacobian ($n\times m$);
their formulas differ from ours by a transpose.

(sec-01-02-partial)=
## Partial and directional derivatives, gradient and Hessian

The **partial derivative** $\partial f/\partial x_j(\vx)$ is the ordinary
derivative of $f$ with respect to $x_j$ while the other variables are held
fixed, that is, the derivative of $f$ along the unit vector $\ve_j$. Partial
derivatives describe $f$ only along the coordinate axes; a function can have
all partial derivatives at a point without being well approximated by a
linear function there. The useful notion is the following.

:::{prf:definition} Differentiability and gradient
:label: def-01-02-gradient

$f:\R^n\to\R$ is **differentiable** at $\vx$ if there is a vector
$\nabla f(\vx)\in\R^n$, the **gradient**, such that

$$
f(\vx+\vd)=f(\vx)+\nabla f(\vx)^\T\vd+o(\norm{\vd}) .
$$ (eq-01-02-differentiable)

Taking $\vd=\alpha\ve_j$ with $\alpha\to0$ shows that the gradient is the
vector of partial derivatives,
$\nabla f(\vx)=\bigl(\partial f/\partial x_1(\vx),\dots,\partial
f/\partial x_n(\vx)\bigr)^\T$. A function whose partial derivatives exist
and are continuous on an open set is **continuously differentiable**
($C^1$); every $C^1$ function is differentiable.
:::

:::{prf:definition} Directional derivative
:label: def-01-02-directional-derivative

The **directional derivative** of $f$ at $\vx$ in the direction $\vd$ is

$$
f'(\vx;\vd)=\lim_{\alpha\downarrow0}\frac{f(\vx+\alpha\vd)-f(\vx)}{\alpha},
$$

the initial rate of change of $f$ when moving from $\vx$ along $\vd$. If
$f$ is differentiable at $\vx$, then substituting $\alpha\vd$ into
{eq}`eq-01-02-differentiable` gives

$$
f'(\vx;\vd)=\nabla f(\vx)^\T\vd .
$$ (eq-01-02-directional)
:::

Formula {eq}`eq-01-02-directional` is how gradients enter optimization: the
sign of $\nabla f(\vx)^\T\vd$ says whether $f$ increases or decreases along
$\vd$ for small steps. A direction with $\nabla f(\vx)^\T\vd<0$ is a
**descent direction** ({prf:ref}`def-06-02-descent-direction`). Directional
derivatives are also defined for some nondifferentiable functions, such as
$\abs{x}$ at $0$, where $f'(0;d)=\abs{d}$ is not linear in $d$.

:::{prf:definition} Hessian
:label: def-01-02-hessian

If $f$ is twice differentiable at $\vx$, its **Hessian** is the $n\times n$
matrix of second partial derivatives,

$$
\nabla^2 f(\vx)=\left(\frac{\partial^2 f}{\partial x_i\,\partial x_j}(\vx)\right)_{i,j=1}^n .
$$

If the second partial derivatives are continuous ($f\in C^2$), the order of
differentiation does not matter (Schwarz's theorem), so the Hessian is
**symmetric**.
:::

Symmetry is why the theory of symmetric matrices in
{doc}`01-01-linear-algebra` (real eigenvalues, the spectral theorem, the
definiteness tests) applies to Hessians.

:::{prf:example} Gradient and Hessian of a quadratic function
:label: ex-01-02-gradient-hessian-quadratic

Let $f(\vx)=2x_1^2+3x_1x_2+x_2^2$ on $\R^2$.

Find its gradient and Hessian, its directional derivatives at
$\bar\vx=(1,1)$ in the unit directions $(1,-1)/\sqrt2$ and $(-3,4)/5$, and
the nature of its stationary point.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Gradient.** Holding $x_2$ fixed, $\partial f/\partial x_1=4x_1+3x_2$;
holding $x_1$ fixed, $\partial f/\partial x_2=3x_1+2x_2$. So

$$
\begin{aligned}
\nabla f(\vx)&=\begin{pmatrix}4x_1+3x_2\\3x_1+2x_2\end{pmatrix},\\
\nabla^2 f(\vx)&=\begin{pmatrix}4&3\\3&2\end{pmatrix}.
\end{aligned}
$$

The Hessian is constant (as for every quadratic function) and symmetric:
$\partial^2f/\partial x_1\partial x_2=\partial^2f/\partial x_2\partial x_1=3$.

**Directional derivatives at $\bar\vx=(1,1)$.** Here $\nabla f(\bar\vx)=(7,5)$.
In the unit direction $\vd=(1,-1)/\sqrt2$,
$f'(\bar\vx;\vd)=(7-5)/\sqrt2=\sqrt2\approx1.414$, so $f$ increases. In
the unit direction $\vd=(-3,4)/5$, $f'(\bar\vx;\vd)=(-21+20)/5=-0.2<0$, so
$\vd$ is a descent direction although it is far from $-\nabla f$. The
largest rate over all unit directions is $\norm{(7,5)}=\sqrt{74}\approx8.602$
({prf:ref}`thm-01-02-steepest-ascent`).

**Link with quadratic forms.** $f(\vx)=\vx^\T\mA\vx$ with the symmetric
$\mA=\begin{pmatrix}2&1.5\\1.5&1\end{pmatrix}$, and indeed
$\nabla f=2\mA\vx$ and $\nabla^2f=2\mA$ ({numref}`tab-01-02-vector-derivatives`).
The Hessian has determinant $8-9=-1<0$, so it is indefinite; its
eigenvalues are $3\pm\sqrt{10}$, about $-0.1623$ and $6.162$. The only
stationary point is $\vx=\vzero$ (the system $4x_1+3x_2=0$,
$3x_1+2x_2=0$ has determinant $-1\ne0$), and it is a saddle point: $f$
factors as $(2x_1+x_2)(x_1+x_2)$, which is positive on some lines through
the origin and negative on others (compare {numref}`fig-01-01-quadratic-forms`(c)).
:::

(sec-01-02-jacobian)=
## The Jacobian and the derivative table

For a vector function $\vg:\R^n\to\R^m$ differentiable at $\vx$, the
Jacobian {eq}`eq-01-02-jacobian-rows` gives the first-order approximation

$$
\vg(\vx+\vd)=\vg(\vx)+D\vg(\vx)\,\vd+o(\norm{\vd}) .
$$ (eq-01-02-jacobian-approx)

Two examples recur in the book. An **affine map** $\vg(\vx)=\mA\vx+\vb$ has
$D\vg(\vx)=\mA$ at every point: linear constraints have constant
gradients, the rows of $\mA$. A **nonlinear constraint vector** such as
$\vg(\vx)=(x_1^2+x_2^2-1,\ x_1x_2)$ has

$$
D\vg(\vx)=\begin{pmatrix}2x_1&2x_2\\x_2&x_1\end{pmatrix},
$$

whose rows are $\nabla g_1(\vx)^\T$ and $\nabla g_2(\vx)^\T$; at
$\vx=(1,0)$ it is $\begin{pmatrix}2&0\\0&1\end{pmatrix}$, with linearly
independent rows. Independence of the rows of the Jacobian of the active
constraints is the linear independence constraint qualification (LICQ)
of {prf:ref}`thm-06-03-cq-sufficient`.

The derivatives in {numref}`tab-01-02-vector-derivatives` are used throughout
the book. Each is derived from the definitions in a few lines.

```{list-table} Gradients and Hessians of common scalar functions of $\vx\in\R^n$ ($\va\in\R^n$, $\mA\in\R^{n\times n}$, $\mM\in\R^{m\times n}$, $\vb\in\R^m$, $\vc\in\R^n$, $\mQ$ symmetric).
:name: tab-01-02-vector-derivatives
:header-rows: 1

* - $f(\vx)$
  - $\nabla f(\vx)$
  - $\nabla^2 f(\vx)$
* - $\va^\T\vx$
  - $\va$
  - $\vzero$
* - $\vx^\T\mA\vx$
  - $(\mA+\mA^\T)\vx$
  - $\mA+\mA^\T$
* - $\vx^\T\mQ\vx$ ($\mQ$ symmetric)
  - $2\mQ\vx$
  - $2\mQ$
* - $\tfrac12\vx^\T\mQ\vx+\vc^\T\vx$
  - $\mQ\vx+\vc$
  - $\mQ$
* - $\norm{\mM\vx-\vb}_2^2$
  - $2\mM^\T(\mM\vx-\vb)$
  - $2\mM^\T\mM$
* - $\norm{\vx}_2$ ($\vx\ne\vzero$)
  - $\vx/\norm{\vx}_2$
  - $\bigl(\mI-\vx\vx^\T/\norm{\vx}_2^2\bigr)/\norm{\vx}_2$
* - $\vu(\vx)^\T\vv(\vx)$
  - $D\vu(\vx)^\T\vv(\vx)+D\vv(\vx)^\T\vu(\vx)$
  - —
* - $f(\mM\vy+\vb)$ as a function of $\vy$
  - $\mM^\T\nabla f(\mM\vy+\vb)$
  - $\mM^\T\nabla^2 f(\mM\vy+\vb)\mM$
```

**Linear function.** $\va^\T\vx=\sum_ja_jx_j$, so $\partial/\partial x_k=a_k$.

**Quadratic form.** $\vx^\T\mA\vx=\sum_i\sum_ja_{ij}x_ix_j$. The variable
$x_k$ appears in the terms with $i=k$ and in those with $j=k$, so

$$
\begin{aligned}
\frac{\partial}{\partial x_k}\vx^\T\mA\vx&=\sum_ja_{kj}x_j+\sum_ia_{ik}x_i\\
&=(\mA\vx)_k+(\mA^\T\vx)_k ,
\end{aligned}
$$

which gives $\nabla(\vx^\T\mA\vx)=(\mA+\mA^\T)\vx$. Differentiating the
linear function $(\mA+\mA^\T)\vx$ once more gives the constant Hessian
$\mA+\mA^\T$. For symmetric $\mQ$ these become $2\mQ\vx$ and $2\mQ$; the
factor $\tfrac12$ in $\tfrac12\vx^\T\mQ\vx$ is chosen to cancel the 2.
In {prf:ref}`ex-01-02-gradient-hessian-quadratic` the nonsymmetric matrix
$\begin{pmatrix}2&3\\0&1\end{pmatrix}$ represents the same $f$, and
$(\mA+\mA^\T)\vx$ gives the same gradient.

**Least squares.** Expanding,
$\norm{\mM\vx-\vb}_2^2=\vx^\T\mM^\T\mM\vx-2\vb^\T\mM\vx+\vb^\T\vb$, a quadratic form with symmetric matrix $\mM^\T\mM$ plus a linear term with $\va=\mM^\T\vb$. The previous two rows of
the table give the gradient $2\mM^\T\mM\vx-2\mM^\T\vb$ and
the Hessian $2\mM^\T\mM\succeq0$.

**Euclidean norm.** $\norm{\vx}_2=(\vx^\T\vx)^{1/2}$ is a function of
$s=\vx^\T\vx$, and the chain rule ({ref}`sec-01-02-chain-rule`) gives
$\tfrac12 s^{-1/2}\cdot2\vx=\vx/\norm{\vx}_2$: a unit vector pointing
away from the origin. The norm is not differentiable at $\vzero$ (it
behaves like $\abs{x}$ there), which is why optimization models use the
smooth $\norm{\vx}_2^2$ wherever possible. Its Hessian satisfies
$\nabla^2\norm{\vx}_2\,\vx=\vzero$: along a ray from the origin the norm
grows linearly, without curvature.

**Product rule.** For $\vu,\vv:\R^n\to\R^m$, the inner product
$\vu^\T\vv=\sum_iu_iv_i$ has gradient
$\sum_i(v_i\nabla u_i+u_i\nabla v_i)$, which is the entry of the table
by {eq}`eq-01-02-jacobian-rows`.

(sec-01-02-chain-rule)=
## The chain rule

:::{prf:theorem} Chain rule
:label: thm-01-02-chain-rule

Let $\vu:\R^n\to\R^m$ be differentiable at $\vx$ and $\vv:\R^m\to\R^p$ be
differentiable at $\vu(\vx)$. Then $\vv\circ\vu$ is differentiable at $\vx$
and

$$
D(\vv\circ\vu)(\vx)=D\vv\bigl(\vu(\vx)\bigr)\,D\vu(\vx).
$$ (eq-01-02-chain-rule)

For a scalar outer function $f:\R^m\to\R$ this reads
$\nabla(f\circ\vu)(\vx)=D\vu(\vx)^\T\,\nabla f\bigl(\vu(\vx)\bigr)$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

By {eq}`eq-01-02-jacobian-approx` for $\vu$,
$\vu(\vx+\vd)=\vu(\vx)+\ve$ with $\ve=D\vu(\vx)\vd+o(\norm{\vd})$, and
$\norm{\ve}\le c\norm{\vd}$ for small $\vd$ and some constant $c$.
Applying {eq}`eq-01-02-jacobian-approx` to $\vv$ at $\vu(\vx)$ with increment
$\ve$ gives $\vv(\vu(\vx+\vd))=\vv(\vu(\vx))+D\vv(\vu(\vx))\ve+o(\norm{\ve})
=\vv(\vu(\vx))+D\vv(\vu(\vx))D\vu(\vx)\vd+o(\norm{\vd})$. The scalar form
is the transpose of the row $Df(\vu(\vx))D\vu(\vx)$.
:::

The Jacobians multiply in the same order as the functions are composed,
exactly as in one variable; the layout convention of
{prf:ref}`def-01-02-layout` is what makes this work. Three special cases are
used repeatedly.

1. **Affine change of variables.** For $\vu(\vy)=\mM\vy+\vb$,
   $D\vu=\mM$, so
   $\nabla_{\vy}f(\mM\vy+\vb)=\mM^\T\nabla f(\mM\vy+\vb)$;
   differentiating once more gives the Hessian
   $\mM^\T\nabla^2 f\,\mM$ (last row of
   {numref}`tab-01-02-vector-derivatives`).
2. **Restriction to a line.** For fixed $\vx$ and $\vd$, the function
   $\varphi(\alpha)=f(\vx+\alpha\vd)$ of one variable has

   $$
   \begin{aligned}
   \varphi'(\alpha)&=\nabla f(\vx+\alpha\vd)^\T\vd,\\
   \varphi''(\alpha)&=\vd^\T\nabla^2 f(\vx+\alpha\vd)\,\vd .
   \end{aligned}
   $$ (eq-01-02-line-restriction)

   Line searches ({doc}`06-04-line-search`, {eq}`eq-06-04-phi`) minimize exactly this
   $\varphi$, and {eq}`eq-01-02-line-restriction` connects their
   one-dimensional derivatives with $\nabla f$ and $\nabla^2 f$.
3. **Composition with a smooth scalar function.** For the
   "log-sum-exp" function $f(\vx)=\ln\sum_{i=1}^m\exp(\va_i^\T\vx+b_i)$,
   write $f=h\circ\vu$ with $\vu(\vx)=\mM\vx+\vb$ ($\mM$ has rows
   $\va_i^\T$) and $h(\vy)=\ln\sum_i e^{y_i}$, whose gradient has components
   $e^{y_i}/\sum_ke^{y_k}$. The chain rule gives
   $\nabla f(\vx)=\mM^\T\vp$, where $p_i=e^{y_i}/\sum_ke^{y_k}$
   at $\vy=\mM\vx+\vb$: a weighted combination of the $\va_i$ with
   weights that are positive and sum to one.

(sec-01-02-matrix-derivatives)=
## Derivatives with respect to a matrix

Some functions in optimization and statistics take a matrix argument
$\mX\in\R^{m\times n}$: the log-determinant of a covariance matrix,
the trace of a weighted product, the Frobenius norm of a residual. The
derivative of a scalar function $f(\mX)$ is defined entry by entry,
$\bigl(\partial f/\partial\mX\bigr)_{ij}=\partial f/\partial x_{ij}$,
and has the same shape as $\mX$. With the matrix inner product
$\langle\mG,\mE\rangle=\tr(\mG^\T\mE)=\sum_{i,j}
g_{ij}e_{ij}$, the analog of {eq}`eq-01-02-differentiable` is

$$
f(\mX+\mE)=f(\mX)+\tr\Bigl(\bigl(\tfrac{\partial f}{\partial\mX}\bigr)^\T\mE\Bigr)+o(\norm{\mE}) .
$$ (eq-01-02-matrix-gradient)

This gives a practical recipe: expand $f(\mX+\mE)$, keep the
terms linear in $\mE$, bring them into the form
$\tr(\mG^\T\mE)$ with the cyclic rule $\tr(\mP\mQ)=\tr(\mQ\mP)$, and read off $\partial f/\partial\mX=\mG$. {numref}`tab-01-02-matrix-derivatives` lists the results needed most often
(Petersen & Pedersen, 2012, Section 2; Magnus & Neudecker, 2019,
Chapters 8–9).

```{list-table} Derivatives of scalar functions of a matrix $\mX$ (entries treated as independent; $\mA$, $\va$, $\vb$ constant and of compatible sizes; $\mX$ square and invertible in the last three rows, with $\det\mX>0$ for the logarithm).
:name: tab-01-02-matrix-derivatives
:header-rows: 1

* - $f(\mX)$
  - $\partial f/\partial\mX$
* - $\tr(\mA\mX)$
  - $\mA^\T$
* - $\tr(\mX^\T\mA)$
  - $\mA$
* - $\va^\T\mX\vb$
  - $\va\vb^\T$
* - $\tr(\mX^\T\mA\mX)$
  - $(\mA+\mA^\T)\mX$
* - $\norm{\mX}_F^2=\tr(\mX^\T\mX)$
  - $2\mX$
* - $\ln\det\mX$
  - $\mX^{-\T}$
* - $\det\mX$
  - $(\det\mX)\,\mX^{-\T}$
* - $\va^\T\mX^{-1}\vb$
  - $-\mX^{-\T}\va\vb^\T\mX^{-\T}$
```

Two entries are derived here.

**The trace of a product.** $\tr(\mA\mX)=\sum_i\sum_ka_{ik}x_{ki}$.
The entry $x_{kl}$ occurs once, in the term with $i=l$, with coefficient
$a_{lk}$. Hence $\partial\tr(\mA\mX)/\partial x_{kl}=a_{lk}$, and
$\partial\tr(\mA\mX)/\partial\mX=\mA^\T$. The same result
follows from the recipe, since $\tr(\mA(\mX+\mE))
=\tr(\mA\mX)+\tr\bigl((\mA^\T)^\T\mE\bigr)$. Because
$\va^\T\mX\vb=\tr(\vb\va^\T\mX)$, the third row of the table is
the special case $\mA=\vb\va^\T$.

**The log-determinant.** For invertible $\mX$ and a small $\mE$,

$$
\begin{aligned}
\det(\mX+\mE)&=\det\mX\cdot\det(\mI+\mX^{-1}\mE)\\
&=\det\mX\cdot\bigl(1+\tr(\mX^{-1}\mE)+O(\norm{\mE}^2)\bigr),
\end{aligned}
$$

because in the expansion of $\det(\mI+\mF)$ the product of the
diagonal entries contributes $1+\sum_if_{ii}+O(\norm{\mF}^2)$ and
every other term contains at least two off-diagonal entries of $\mF$.
Taking logarithms, $\ln(1+t)=t+O(t^2)$ gives
$\ln\det(\mX+\mE)=\ln\det\mX+\tr(\mX^{-1}\mE)
+O(\norm{\mE}^2)$. Writing $\tr(\mX^{-1}\mE)=\tr\bigl((\mX^{-\T})^\T\mE\bigr)$ and comparing with
{eq}`eq-01-02-matrix-gradient` gives $\partial\ln\det\mX/\partial\mX=\mX^{-\T}$, which is $\mX^{-1}$ for symmetric $\mX$. The
chain rule then gives the derivative of $\det\mX=\exp(\ln\det\mX)$ in the next row.

The last row comes from the first-order expansion of the inverse,
$(\mX+\mE)^{-1}=\mX^{-1}-\mX^{-1}\mE\mX^{-1}+O(\norm{\mE}^2)$, which is verified by multiplying with
$\mX+\mE$: the product is
$\mI-\mE\mX^{-1}\mE\mX^{-1}$.

(sec-01-02-taylor)=
## Taylor's theorem

Taylor expansions replace a function near a point $\bar\vx$ by a
polynomial: the first-order (linear) model
$f(\bar\vx)+\nabla f(\bar\vx)^\T\vd$ and the second-order (quadratic) model
$f(\bar\vx)+\nabla f(\bar\vx)^\T\vd+\tfrac12\vd^\T\nabla^2f(\bar\vx)\vd$ in
the step $\vd=\vx-\bar\vx$. The theorem says how large the error is.

:::{prf:theorem} Taylor's theorem
:label: thm-01-02-taylor

Let $f$ be defined on an open set containing the segment from $\bar\vx$ to
$\bar\vx+\vd$.

1. If $f\in C^1$, there is a $t\in(0,1)$ with

   $$
   f(\bar\vx+\vd)=f(\bar\vx)+\nabla f(\bar\vx+t\vd)^\T\vd ,
   $$ (eq-01-02-taylor-1)

   and $f(\bar\vx+\vd)=f(\bar\vx)+\nabla f(\bar\vx)^\T\vd+o(\norm{\vd})$.
2. If $f\in C^2$, there is a $t\in(0,1)$ with

   $$
   f(\bar\vx+\vd)=f(\bar\vx)+\nabla f(\bar\vx)^\T\vd+\tfrac12\,\vd^\T\nabla^2 f(\bar\vx+t\vd)\,\vd ,
   $$ (eq-01-02-taylor-2)

   and $f(\bar\vx+\vd)=f(\bar\vx)+\nabla f(\bar\vx)^\T\vd
   +\tfrac12\vd^\T\nabla^2f(\bar\vx)\vd+o(\norm{\vd}^2)$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\varphi(s)=f(\bar\vx+s\vd)$ for $s\in[0,1]$. By
{eq}`eq-01-02-line-restriction`, $\varphi'(s)=\nabla f(\bar\vx+s\vd)^\T\vd$ and,
for $f\in C^2$, $\varphi''(s)=\vd^\T\nabla^2f(\bar\vx+s\vd)\vd$. The
one-variable mean value theorem gives $\varphi(1)=\varphi(0)+\varphi'(t)$
for some $t\in(0,1)$, which is {eq}`eq-01-02-taylor-1`; Taylor's formula with
Lagrange remainder gives $\varphi(1)=\varphi(0)+\varphi'(0)+\tfrac12
\varphi''(t)$, which is {eq}`eq-01-02-taylor-2`. For the $o$-forms, replace
the intermediate point by $\bar\vx$; the error is
$\bigl(\nabla f(\bar\vx+t\vd)-\nabla f(\bar\vx)\bigr)^\T\vd$, respectively
$\tfrac12\vd^\T\bigl(\nabla^2f(\bar\vx+t\vd)-\nabla^2f(\bar\vx)\bigr)\vd$,
which by Cauchy–Schwarz and continuity of $\nabla f$ (of $\nabla^2f$) is
$o(\norm{\vd})$ (respectively $o(\norm{\vd}^2)$).
:::

Note what is evaluated where: in the exact forms the gradient in
{eq}`eq-01-02-taylor-2` sits at the expansion point and only the highest
derivative sits at the unknown intermediate point $\bar\vx+t\vd$. The
exact forms are the workhorses of proofs: if, for instance,
$\nabla f(\bar\vx)=\vzero$ and $\nabla^2f$ is positive semidefinite on the
whole segment, {eq}`eq-01-02-taylor-2` gives $f(\bar\vx+\vd)\ge f(\bar\vx)$ at
once ({prf:ref}`thm-02-02-second-order`,
{prf:ref}`thm-06-02-second-order-sufficient`). The
$o$-forms are the basis of algorithms: steepest descent minimizes the
linear model over a ball, Newton's method minimizes the quadratic model
({prf:ref}`alg-06-05-steepest-descent`, {eq}`eq-06-05-newton-model`).

:::{prf:example} Linear and quadratic models of a function of one variable
:label: ex-01-02-taylor-models

Let $f(x)=e^x-2x$ and $\bar x=0$.

Find the linear and quadratic Taylor models $m_1$ and $m_2$ of $f$ at
$\bar x$, compare them with $f$ at $x=0.5$ and $x=2$, and compare the
minimizer of $m_2$ with that of $f$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Then $f(0)=1$, $f'(x)=e^x-2$, $f'(0)=-1$,
$f''(x)=e^x$, $f''(0)=1$. The two models are

$$
\begin{aligned}
m_1(x)&=1-x,\\
m_2(x)&=1-x+\tfrac12x^2 .
\end{aligned}
$$

At $x=0.5$: $f(0.5)=e^{0.5}-1\approx0.6487$, $m_1(0.5)=0.5$ (error
$0.1487$) and $m_2(0.5)=0.625$ (error $0.0237$). The quadratic model is
six times more accurate here. The error of $m_2$ is the third-order
remainder $f'''(\xi)\,0.5^3/6=e^{\xi}/48$ for some $\xi\in(0,0.5)$, so it
must lie between $1/48\approx0.0208$ and $e^{0.5}/48\approx0.0343$, and it
does. The exact form {eq}`eq-01-02-taylor-2` holds with $t\approx0.3475$
(that is, $f(0.5)=1-0.5+\tfrac12e^{0.5t}\,0.25$), and {eq}`eq-01-02-taylor-1`
with $t\approx0.5208$.

{numref}`fig-01-02-taylor-models` shows that both models are good only near
$\bar x$. The minimizer of the quadratic model is $x=1$, a
reasonable first guess for the true minimizer $\ln2\approx0.6931$ (where
$f\approx0.6137$), and a step to it is one step of Newton's method. Far
away the model fails: $f(2)\approx3.389$ but $m_2(2)=1$.

```{figure} figures/01-02-taylor-models.svg
:name: fig-01-02-taylor-models
:width: 520px
:class: dark-light
:alt: The convex curve e to the x minus 2x with its tangent line at 0 and the osculating parabola at 0; the parabola has its minimum at x equal to 1, the curve at ln 2

First- and second-order Taylor models of $f(x)=e^x-2x$ at $\bar x=0$. The
tangent line (dashed) matches value and slope, the parabola also matches
the curvature and follows $f$ further. The minimizer $x=1$ of the
quadratic model approximates the minimizer $\ln2$ of $f$.
```
:::

(sec-01-02-gradient-geometry)=
## Geometry of the gradient

The **level set** of $f$ at height $c$ is $L_c(f)=\{\vx:f(\vx)=c\}$. In
$\R^2$ level sets are the contour lines of a map; in $\R^3$ they are
surfaces. Two facts make the gradient the central object of descent
methods.

:::{prf:theorem} The gradient is orthogonal to level sets
:label: thm-01-02-gradient-level-set

Let $f\in C^1$, $\bar\vx\in L_c(f)$, and let $\vgamma(s)$, $s\in(-\varepsilon,
\varepsilon)$, be a differentiable curve in $L_c(f)$ with
$\vgamma(0)=\bar\vx$. Then $\nabla f(\bar\vx)^\T\vgamma'(0)=0$: the
gradient is orthogonal to every tangent direction of the level set.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

$f(\vgamma(s))=c$ for all $s$. Differentiating with the chain rule,
$\nabla f(\vgamma(s))^\T\vgamma'(s)=0$; put $s=0$.
:::

:::{prf:theorem} Steepest ascent and descent
:label: thm-01-02-steepest-ascent

Let $f$ be differentiable at $\vx$ with $\nabla f(\vx)\ne\vzero$.

1. Among all unit vectors $\vd$, the directional derivative
   $f'(\vx;\vd)=\nabla f(\vx)^\T\vd$ is largest for
   $\vd=\nabla f(\vx)/\norm{\nabla f(\vx)}$, with value
   $\norm{\nabla f(\vx)}$, and smallest for
   $\vd=-\nabla f(\vx)/\norm{\nabla f(\vx)}$, with value
   $-\norm{\nabla f(\vx)}$.
2. If $\nabla f(\vx)^\T\vd<0$, then $f(\vx+\alpha\vd)<f(\vx)$ for all
   sufficiently small $\alpha>0$. In particular $-\nabla f(\vx)$ is a
   descent direction and $\nabla f(\vx)$ an ascent direction.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

(1) By Cauchy–Schwarz, $\abs{\nabla f(\vx)^\T\vd}\le\norm{\nabla f(\vx)}$
for $\norm{\vd}=1$, with equality exactly when $\vd$ is parallel to the
gradient. (2) By {eq}`eq-01-02-differentiable`,
$f(\vx+\alpha\vd)-f(\vx)=\alpha\bigl(\nabla f(\vx)^\T\vd+o(\alpha)/\alpha\bigr)$,
and the bracket is negative for small $\alpha>0$.
:::

The two theorems fit together: moving along a level set does not change
$f$ to first order, so the direction of fastest change must be
perpendicular to it. Two cautions. First, "steepest" is a local
statement: a long step along $-\nabla f$ can overshoot and increase $f$,
which is why descent methods need a step-size rule
({doc}`06-04-line-search`). Second, the gradient is perpendicular to the level
sets in the space of $\vx$, not to the graph of $f$ in $\R^{n+1}$.

:::{prf:example} Level lines of a linear function
:label: ex-01-02-linear-level-lines

Let $f(\vx)=x_1+x_2$.

Describe its gradient and its level sets, check
{prf:ref}`thm-01-02-gradient-level-set`, and find how much a step of length
$0.5$ from $\bar\vx=(2,1)$ can raise $f$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Its gradient $\nabla f=(1,1)$ is the same at every
point. The level sets are the parallel lines $x_1+x_2=c$, with direction
$(1,-1)$, and $(1,1)^\T(1,-1)=0$ confirms
{prf:ref}`thm-01-02-gradient-level-set`. At $\bar\vx=(2,1)$ the value is
$f=3$; moving from $\bar\vx$ along $\nabla f$ crosses the levels $3.5$,
$4$, … and moving along $-\nabla f$ crosses $2.5$, $2$, …
({numref}`fig-01-02-linear-level-lines`). A step of length $0.5$ along the
unit gradient $(1,1)/\sqrt2$ raises $f$ by $0.5\sqrt2\approx0.7071$; no step
of the same length does better, since the largest rate is
$\norm{\nabla f}=\sqrt2$.

```{figure} figures/01-02-linear-level-lines.svg
:name: fig-01-02-linear-level-lines
:width: 460px
:class: dark-light
:alt: Parallel straight level lines of x1 plus x2 for the values 0 to 3, a point at (2, 1) on the line f equal to 3, and two arrows from it, the gradient (1, 1) pointing up-right and its negative pointing down-left, both perpendicular to the lines

Level lines $f=c$ of $f(\vx)=x_1+x_2$ (labeled for $c=0,\dots,3$). At
$(2,1)$ the gradient $(1,1)$ points to higher levels and $-\nabla f$ to
lower levels; both are perpendicular to the level lines.
```
:::

For a linear function the opposite gradient
decreases $f$ at the same rate forever, which is why a linear objective
over an unbounded feasible set can be unbounded
({ref}`sec-03-02-graphical-cases`).

:::{prf:example} Contours of a paraboloid
:label: ex-01-02-paraboloid-contours

Let $f(\vx)=(x_1-2)^2+(x_2-2)^2$, the squared distance from $(2,2)$.

Describe its graph and level sets, compute its gradient and Hessian, and
check orthogonality to the level set, ascent and descent at
$\bar\vx=(1.4,1.4)$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Its
graph is a paraboloid with minimum value $0$ at $(2,2)$, and its level sets
$L_c(f)$, $c>0$, are circles of radius $\sqrt c$ around $(2,2)$
({numref}`fig-01-02-paraboloid-contours`). The gradient and Hessian are

$$
\begin{aligned}
\nabla f(\vx)&=\begin{pmatrix}2(x_1-2)\\2(x_2-2)\end{pmatrix},\\
\nabla^2f(\vx)&=2\mI .
\end{aligned}
$$

At $\bar\vx=(1.4,1.4)$: $f(\bar\vx)=0.36+0.36=0.72$, so $\bar\vx$ lies on
the circle of radius $\sqrt{0.72}\approx0.8485$, and
$\nabla f(\bar\vx)=(-1.2,-1.2)$.

- **Orthogonality.** The tangent to the circle at $\bar\vx$ has direction
  $(1,-1)$, and $(-1.2,-1.2)^\T(1,-1)=0$.
- **Ascent.** The gradient points away from the center, toward larger
  circles, with rate $\norm{\nabla f(\bar\vx)}=1.2\sqrt2\approx1.697$; along
  the tangent the rate is 0, along $\ve_1$ it is $-1.2$.
- **Descent.** $-\nabla f(\bar\vx)=(1.2,1.2)$ points straight at the
  minimizer: $\bar\vx-\tfrac12\nabla f(\bar\vx)=(2,2)$. This happens because
  the contours are circles ($\nabla^2f=2\mI$); for elongated elliptic
  contours $-\nabla f$ misses the center, which is the weakness of
  steepest descent studied in {doc}`06-05-unconstrained-algorithms`
  ({prf:ref}`thm-06-05-sd-rate`).

```{figure} figures/01-02-paraboloid-contours.svg
:name: fig-01-02-paraboloid-contours
:width: 720px
:class: dark-light
:alt: Left, a bowl-shaped surface over the square from 1 to 3 with a circular level curve drawn on it and its shadow on the floor. Right, concentric circular contours around (2, 2) with the gradient at (1.4, 1.4) pointing away from the center, its negative pointing to the center, and a dotted tangent line

(a) Graph of $f(\vx)=(x_1-2)^2+(x_2-2)^2$ with the level $f=0.72$ on the
surface and its projection, a circle in the $\vx$-plane. (b) Contours
$f=0.1,0.25,0.5,0.72,1,1.5$. At $\bar\vx=(1.4,1.4)$ the gradient
$(-1.2,-1.2)$ (drawn at one third of its length) is perpendicular to the
dotted tangent and points to higher levels; $-\nabla f$ points to the
minimizer $(2,2)$.
```
:::

(sec-01-02-applications)=
## Applications

### Least squares and the normal equations

Fitting a model that is linear in its parameters $\vx\in\R^n$ to $m>n$
observations leads to an overdetermined system $\mM\vx\approx\vb$,
where row $i$ of $\mM$ holds the regressors of observation $i$ and
$b_i$ its measured value. The **least-squares** estimate minimizes
$f(\vx)=\norm{\mM\vx-\vb}_2^2$. By
{numref}`tab-01-02-vector-derivatives`,

$$
\nabla f(\vx)=2\mM^\T(\mM\vx-\vb)=\vzero
\quad\Longleftrightarrow\quad
\mM^\T\mM\,\vx=\mM^\T\vb ,
$$ (eq-01-02-normal-equations)

the **normal equations** {eq}`eq-01-01-normal-equations` found geometrically in
{doc}`01-01-linear-algebra`.
The Hessian $2\mM^\T\mM$ is positive semidefinite, so $f$ is
convex and every solution of {eq}`eq-01-02-normal-equations` is a global
minimizer ({prf:ref}`thm-02-02-first-order`); if $\mM$ has full column
rank, $\mM^\T\mM\succ0$ and the minimizer is unique. The
residual $\vr=\vb-\mM\vx^*$ satisfies $\mM^\T\vr=\vzero$: it is
orthogonal to every column of $\mM$.

:::{prf:example} Fitting a straight line to traffic counts
:label: ex-01-02-least-squares

A traffic counter records the flow $y$ on a road (in thousand vehicles per
hour) at times $t=0,1,2,3$ h after the start of the morning: $y=1,3,4,4$.
Fit the line $y\approx x_1+x_2t$ by least squares.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The data matrix has a column of ones and a column of times, and
$\vb=(1,3,4,4)$:

$$
\begin{aligned}
\mM&=\begin{pmatrix}1&0\\1&1\\1&2\\1&3\end{pmatrix},\\
\mM^\T\mM&=\begin{pmatrix}4&6\\6&14\end{pmatrix},\\
\mM^\T\vb&=\begin{pmatrix}12\\23\end{pmatrix}.
\end{aligned}
$$

The normal equations $4x_1+6x_2=12$ and $6x_1+14x_2=23$ (determinant
$56-36=20\ne0$) give $x_2=1$ and $x_1=1.5$. The fitted line is
$y=1.5+1.0\,t$: an intercept of 1500 vehicles per hour and an increase of
1000 vehicles per hour every hour. The residuals $\vr=\vb-\mM\vx^*
=(-0.5,0.5,0.5,-0.5)$ sum to zero and are orthogonal to the time column
($0\cdot(-0.5)+1\cdot0.5+2\cdot0.5+3\cdot(-0.5)=0$), as
$\mM^\T\vr=\vzero$ requires; the minimal sum of squares is $1.0$.
Extrapolating to $t=4$ h gives $5.5$ thousand vehicles per hour, which
should be read with caution: the last two counts suggest the flow has
leveled off, and a straight line cannot represent that.
:::

### Minimizing a quadratic function

Let $f(\vx)=\tfrac12\vx^\T\mQ\vx+\vc^\T\vx$ with symmetric $\mQ$. Then
$\nabla f(\vx)=\mQ\vx+\vc$ and $\nabla^2f=\mQ$, and the shape of $f$ is
decided by the definiteness of $\mQ$ ({prf:ref}`def-01-01-definiteness`).

- **$\mQ\succ0$.** The unique stationary point is $\vx^*=-\mQ^{-1}\vc$, and
  completing the square gives

  $$
  \begin{aligned}
  f(\vx)&=f(\vx^*)+\tfrac12(\vx-\vx^*)^\T\mQ(\vx-\vx^*),\\
  f(\vx^*)&=-\tfrac12\vc^\T\mQ^{-1}\vc ,
  \end{aligned}
  $$ (eq-01-02-quadratic-minimum)

  so $\vx^*$ is the unique global minimizer. For
  $\mQ=\begin{pmatrix}2&1\\1&2\end{pmatrix}$ and $\vc=(-3,0)$, solving
  $\mQ\vx=(3,0)$ gives $\vx^*=(2,-1)$ and $f(\vx^*)=\tfrac12\vc^\T\vx^*=-3$.
- **$\mQ\succeq0$ singular.** A minimizer exists if and only if
  $\mQ\vx=-\vc$ is solvable, that is, $\vc\in\col(\mQ)$; then the
  minimizers form an affine set. For $\mQ=\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$
  and $\vc=(1,-1)$, every $\vx=(-1,0)+s(1,1)$ is a minimizer with value
  $-0.5$. For $\vc=(1,1)$, which is orthogonal to $\col(\mQ)$,
  $f(s,s)=2s$, so $f$ is unbounded below (for example,
  $f(-10,-10)=-20$).
- **$\mQ$ indefinite.** Along an eigenvector $\vv$ of a negative eigenvalue
  $\lambda_i(\mQ)<0$,
  $f(s\vv)=\tfrac12\lambda_i(\mQ)s^2\norm{\vv}^2+s\,\vc^\T\vv\to-\infty$
  as $s\to\pm\infty$; the
  stationary point, if there is one, is a saddle.

Quadratic programs ({ref}`sec-06-01-classes`) and Newton's method
({ref}`sec-06-05-newton`) rely on these three cases.

### The gradient of a Lagrangian

For the nonlinear program $\min f(\vx)$ subject to $g_i(\vx)\le0$
($i=1,\dots,m$) and $h_j(\vx)=0$ ($j=1,\dots,\ell$), the **Lagrangian** is

$$
\begin{aligned}
L(\vx,\vlambda,\vmu)&=f(\vx)+\sum_{i=1}^m\lambda_ig_i(\vx)+\sum_{j=1}^\ell\mu_jh_j(\vx)\\
&=f(\vx)+\vlambda^\T\vg(\vx)+\vmu^\T\vh(\vx).
\end{aligned}
$$

Its gradient with respect to $\vx$ follows from the table and
{eq}`eq-01-02-jacobian-rows`:

$$
\begin{aligned}
\nabla_{\vx}L(\vx,\vlambda,\vmu)&=\nabla f(\vx)+\sum_i\lambda_i\nabla g_i(\vx)+\sum_j\mu_j\nabla h_j(\vx)\\
&=\nabla f(\vx)+D\vg(\vx)^\T\vlambda+D\vh(\vx)^\T\vmu .
\end{aligned}
$$ (eq-01-02-lagrangian-gradient)

Setting it to zero is the stationarity part of the KKT conditions
({eq}`eq-02-03-kkt`, {prf:ref}`thm-06-03-kkt`). For
linear constraints $\vg(\vx)=\mA\vx-\vb$ the Jacobian is $\mA$, and the
condition reads $\nabla f(\vx)+\mA^\T\vlambda=\vzero$.

:::{prf:example} Closest point of a plane through the Lagrangian
:label: ex-01-02-min-norm-lagrangian

Find the point of the plane $x_1+2x_2+2x_3=9$ closest to the origin, that is,
minimize $f(\vx)=\tfrac12\norm{\vx}_2^2$ subject to $h(\vx)=\va^\T\vx-9=0$
with $\va=(1,2,2)$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The Lagrangian is $L(\vx,\mu)=\tfrac12\vx^\T\vx+\mu(\va^\T\vx-9)$, and
{eq}`eq-01-02-lagrangian-gradient` gives $\nabla_{\vx}L=\vx+\mu\va$. Setting it
to zero gives $\vx=-\mu\va$; substituting into the constraint,
$-\mu\,\va^\T\va=-9\mu=9$, so $\mu=-1$ and $\vx^*=\va=(1,2,2)$. Check:
$1+4+4=9$. The optimal value is $f(\vx^*)=4.5$, and the distance from the
origin to the plane is $\norm{\vx^*}=3=9/\norm{\va}$.

That $\vx^*$ is the global minimizer can be seen directly. For any feasible
$\vx$, $\va^\T(\vx-\vx^*)=0$, so $\vx^{*\T}(\vx-\vx^*)=0$ and by Pythagoras
$\norm{\vx}^2=\norm{\vx^*}^2+\norm{\vx-\vx^*}^2\ge\norm{\vx^*}^2$. The
feasible point $(9,0,0)$, for instance, has $f=40.5$.
:::

The same computation as in the solution of
{prf:ref}`ex-01-02-min-norm-lagrangian`, with several equality constraints
$\mA\vx=\vb$ ($\mA$ of full row rank), gives the minimum-norm solution $\vx^*=\mA^\T(\mA\mA^\T)^{-1}\vb$.

## Summary and Common Mistakes

### Summary

- Gradients are columns, Jacobians have the transposed gradients as rows,
  and the Hessian is the Jacobian of the gradient
  ({prf:ref}`def-01-02-layout`, {eq}`eq-01-02-jacobian-rows`).
- Differentiability means a good linear approximation
  {eq}`eq-01-02-differentiable`; then $f'(\vx;\vd)=\nabla f(\vx)^\T\vd$
  {eq}`eq-01-02-directional`, and the Hessian of a $C^2$ function is
  symmetric ({prf:ref}`def-01-02-hessian`).
- Standard gradients and Hessians: $\va^\T\vx$, $\vx^\T\mA\vx$,
  $\norm{\mM\vx-\vb}_2^2$, $\norm{\vx}_2$
  ({numref}`tab-01-02-vector-derivatives`); matrix derivatives of the trace,
  the log-determinant and the inverse ({numref}`tab-01-02-matrix-derivatives`,
  {eq}`eq-01-02-matrix-gradient`).
- The chain rule multiplies Jacobians in the order of composition
  ({prf:ref}`thm-01-02-chain-rule`); restricted to a line it gives
  {eq}`eq-01-02-line-restriction`.
- Taylor's theorem gives exact forms with an intermediate point and
  $o$-forms ({prf:ref}`thm-01-02-taylor`); {prf:ref}`ex-01-02-taylor-models` and
  {numref}`fig-01-02-taylor-models` compare the models.
- The gradient is orthogonal to level sets ({prf:ref}`thm-01-02-gradient-level-set`)
  and gives the direction of steepest ascent; $-\nabla f$ is the steepest
  descent direction ({prf:ref}`thm-01-02-steepest-ascent`,
  {numref}`fig-01-02-linear-level-lines`, {numref}`fig-01-02-paraboloid-contours`).
- Setting gradients to zero yields the normal equations
  {eq}`eq-01-02-normal-equations`, the minimizer of a positive definite
  quadratic {eq}`eq-01-02-quadratic-minimum` and the stationarity condition of
  the Lagrangian {eq}`eq-01-02-lagrangian-gradient`.

### Common mistakes

- **Mixing layout conventions.** In this book $\nabla f$ is a column and
  the Jacobian of $\vg:\R^n\to\R^m$ is $m\times n$. The Jacobian of
  $\mA\vx$ is $\mA$, but the gradient of $f(\mA\vx)$ is
  $\mA^\T\nabla f(\mA\vx)$; forgetting this transpose produces vectors of
  the wrong length.
- **Using $2\mA\vx$ for a nonsymmetric $\mA$.** The gradient of
  $\vx^\T\mA\vx$ is $(\mA+\mA^\T)\vx$; it equals $2\mA\vx$ only when $\mA$
  is symmetric.
- **Misplacing the arguments in Taylor's theorem.** In the exact
  second-order form {eq}`eq-01-02-taylor-2`, $f$ and $\nabla f$ are evaluated
  at the expansion point $\bar\vx$ and only the Hessian at the intermediate
  point $\bar\vx+t\vd$; the factor $\tfrac12$ must not be dropped, and
  $\nabla f(\bar\vx)$ is the gradient of $f$ evaluated at $\bar\vx$, not a
  derivative with respect to $\bar\vx$.
- **Reading "along the gradient $f$ increases" as a global statement.** It
  holds for sufficiently small steps only; a long step can overshoot.
- **Expecting $-\nabla f$ to point at the minimizer.** It does so only when
  the contours are circles; in general it is merely perpendicular to the
  contour through the current point.
- **Differentiating $\norm{\vx}_2$ at $\vzero$.** The Euclidean norm is not
  differentiable at the origin; use $\norm{\vx}_2^2$ when a smooth function
  is needed.
- **Treating symmetric matrix variables as unconstrained.**
  {numref}`tab-01-02-matrix-derivatives` treats all entries of $\mX$ as
  independent. If $\mX$ is restricted to symmetric matrices, the
  off-diagonal derivatives must account for $x_{ij}=x_{ji}$.

## References

- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Appendix
  A.4: derivatives, gradient, chain rule, second derivative; Appendix
  A.5.5: Schur complement; Section 6.1: least-squares approximation.)
- Magnus, J. R., & Neudecker, H. (2019). *Matrix differential calculus with
  applications in statistics and econometrics* (3rd ed.). Wiley.
  <https://doi.org/10.1002/9781119541219> (Chapters 5–6: differentials,
  Jacobians, the chain rule and Taylor's theorem; Chapters 8–9: first
  differentials of the trace, the determinant and the inverse, and the
  layout of Jacobian matrices.)
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Section 2.1:
  Taylor's theorem; Appendix A.2: elements of analysis, directional
  derivatives.)
- Petersen, K. B., & Pedersen, M. S. (2012). *The matrix cookbook* (Version
  November 15, 2012). Technical University of Denmark. (Section 2:
  derivatives of determinants, inverses, traces and quadratic forms.)
- Strang, G. (2016). *Introduction to linear algebra* (5th ed.).
  Wellesley-Cambridge Press. (Section 4.3: least-squares approximations.)
