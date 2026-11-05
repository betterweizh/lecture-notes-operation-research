# Notation

The book uses one notation throughout. Vectors are bold lower-case letters
($\vx$, $\vc$, $\vpi$) and matrices bold upper-case letters ($\mA$, $\mB$);
scalars and the components of vectors and matrices are italic ($x_j$,
$a_{ij}$). All vectors are column vectors: $\vx^\T$ is the transpose, so
$\vc^\T\vx$ is a scalar product, and the gradient $\nabla f(\vx)$ is a column
vector. Inequalities between vectors, such as $\vx\ge\vzero$ or
$\mA\vx\le\vb$, hold component by component; the symbols $\succeq$ and
$\succ$ are kept for the (semi)definiteness of symmetric matrices.
Eigenvalues are always written with their matrix, $\lambda_i(\mA)$, so that
$\lambda_i$ alone is free for Lagrange multipliers. The tables list the
symbols that are used in more than one chapter or that a reader is likely to
look up; symbols that belong to a single chapter or example are defined
where they are introduced (theory chapters also list them in their Notation
section). The last column names the chapter in which a symbol is
introduced.

## General and linear algebra

| Symbol | Meaning | First used in |
|---|---|---|
| $\R^n$, $\R^{m\times n}$ | real $n$-vectors, real $m\times n$ matrices | {doc}`01-01-linear-algebra` |
| $\vx$, $x_j$ | vector and its $j$-th component; $\vx_1,\dots,\vx_k$ (bold with a subscript) is a list of vectors | {doc}`01-01-linear-algebra` |
| $\vx^\T$, $\mA^\T$, $\mA^{-\T}$ | transpose; inverse of the transpose | {doc}`01-01-linear-algebra` |
| $\vzero$, $\vone$, $\ve_i$ | zero vector, all-ones vector, $i$-th unit vector | {doc}`01-01-linear-algebra` |
| $\mI$, $\mI_m$ | identity matrix (of order $m$) | {doc}`01-01-linear-algebra` |
| $a_{ij}$, $\mA_{\cdot j}$, $\mA_{i\cdot}$ | entry, column $j$ and row $i$ of $\mA$ | {doc}`01-01-linear-algebra` |
| $\col(\mA)$, $\row(\mA)$, $\nullsp(\mA)$ | column space, row space, null space | {doc}`01-01-linear-algebra` |
| $\rank\mA$, $\det\mA$, $\tr\mA$ | rank, determinant, trace | {doc}`01-01-linear-algebra` |
| $\diag(\cdot)$ | diagonal matrix with the given diagonal | {doc}`01-01-linear-algebra` |
| $\mE$, $\mE_k$; $\veta$ | elementary (pivot, eta) matrix; its eta column | {doc}`01-01-linear-algebra` |
| $\bar\mA$, $\bar a_{rs}$ | current tableau; pivot element in row $r$, column $s$ | {doc}`01-01-linear-algebra` |
| $\mL$, $\mU$ | lower and upper triangular factors, $\mA=\mL\mU$ | {doc}`01-01-linear-algebra` |
| $\lambda_i(\mA)$, $\lambda_{\min}(\mA)$, $\lambda_{\max}(\mA)$ | eigenvalues of $\mA$ | {doc}`01-01-linear-algebra` |
| $\mV$, $\mD$ | eigenvector matrix and diagonal eigenvalue matrix, $\mA=\mV\mD\mV^\T$ | {doc}`01-01-linear-algebra` |
| $\Sym^n$, $\Sym^n_+$, $\Sym^n_{++}$ | symmetric, positive semidefinite, positive definite $n\times n$ matrices | {doc}`01-01-linear-algebra` |
| $\mA\succeq\vzero$, $\mA\succ\vzero$ | $\mA$ positive semidefinite, positive definite | {doc}`01-01-linear-algebra` |
| $D_k$ | $k$-th leading principal minor | {doc}`01-01-linear-algebra` |
| $\norm{\cdot}_1$, $\norm{\cdot}_2$, $\norm{\cdot}_\infty$, $\norm{\cdot}_F$ | vector norms and induced matrix norms, Frobenius norm; $\norm{\cdot}=\norm{\cdot}_2$ | {doc}`01-01-linear-algebra` |
| $\kappa(\mA)$ | condition number $\norm{\mA}\,\norm{\mA^{-1}}$ | {doc}`01-01-linear-algebra` |
| $\mP$ | projection matrix | {doc}`01-01-linear-algebra` |
| $\abs{S}$ | number of elements of a finite set $S$ | {doc}`03-03-simplex-method` |
| $\vx^{(k)}$ | $k$-th term of a sequence; iterate $k$ of an algorithm | {doc}`02-01-convex-sets` |
| $\varepsilon$ | small positive number; tolerance of a stopping rule | {doc}`01-02-matrix-calculus` |

## Calculus

| Symbol | Meaning | First used in |
|---|---|---|
| $\nabla f(\vx)$, $\nabla^2 f(\vx)$ | gradient (a column vector) and Hessian of $f$ | {doc}`01-02-matrix-calculus` |
| $D\vg(\vx)$ | Jacobian of $\vg:\R^n\to\R^m$, an $m\times n$ matrix with rows $\nabla g_i(\vx)^\T$ | {doc}`01-02-matrix-calculus` |
| $f'(\vx;\vd)$ | directional derivative of $f$ at $\vx$ along $\vd$ | {doc}`01-02-matrix-calculus` |
| $C^1$, $C^2$ | once, twice continuously differentiable | {doc}`01-02-matrix-calculus` |
| $\bar\vx$ | expansion point of a Taylor formula; the intermediate point is $\bar\vx+t\vd$, $0<t<1$ (written $\bar\vx+\vartheta\vd$ from {doc}`06-02-optimality-conditions` on) | {doc}`01-02-matrix-calculus` |
| $o(\cdot)$, $O(\cdot)$ | Landau symbols | {doc}`01-02-matrix-calculus` |
| $\mX$, $\partial f/\partial\mX$ | matrix variable and the derivative with respect to it | {doc}`01-02-matrix-calculus` |
| $\mQ$ | symmetric matrix of a quadratic function $\tfrac12\vx^\T\mQ\vx+\vc^\T\vx$ | {doc}`01-02-matrix-calculus` |
| $\varphi(\alpha)$ | $f$ restricted to a line, $\varphi(\alpha)=f(\vx+\alpha\vd)$ | {doc}`01-02-matrix-calculus` |

## Convex analysis and duality

| Symbol | Meaning | First used in |
|---|---|---|
| $f$, $g_i$, $h_j$ | objective, inequality constraints $g_i(\vx)\le0$ ($i=1,\dots,m$), equality constraints $h_j(\vx)=0$ ($j=1,\dots,\ell$) | {doc}`01-02-matrix-calculus` |
| $\vg(\vx)$, $\vh(\vx)$ | vectors of the constraint functions | {doc}`01-02-matrix-calculus` |
| $m$, $\ell$ | numbers of inequality and of equality constraints | {doc}`01-02-matrix-calculus` |
| $\Omega$ | feasible set of an optimization problem | {doc}`01-03-introduction` |
| $\vx^*$ | optimal solution (minimizer) | {doc}`01-02-matrix-calculus` |
| $p^*$ | optimal value of a minimization problem | {doc}`01-03-introduction` |
| $C$, $K$ | a (convex) set; a cone | {doc}`02-01-convex-sets` |
| $\theta$, $\theta_i$ | weights of an affine, convex or conic combination | {doc}`02-01-convex-sets` |
| $\aff S$, $\conv S$, $\cone S$ | affine, convex and conic hull of $S$ | {doc}`02-01-convex-sets` |
| $\R^n_+$ | nonnegative orthant | {doc}`02-01-convex-sets` |
| $B(\vx,r)$ | open ball of radius $r$ about $\vx$ | {doc}`02-01-convex-sets` |
| $\intr S$ | interior of $S$ | {doc}`02-01-convex-sets` |
| $\dom f$ | domain of $f$ | {doc}`02-02-convex-functions` |
| $\operatorname{epi} f$ | epigraph of $f$ | {doc}`02-02-convex-functions` |
| $L(\vx,\vlambda,\vmu)$ | Lagrangian $f(\vx)+\vlambda^\T\vg(\vx)+\vmu^\T\vh(\vx)$ | {doc}`01-02-matrix-calculus` |
| $\lambda_i$, $\vlambda$ | multipliers of the inequality constraints $g_i\le0$ ($\lambda_i\ge0$) | {doc}`01-02-matrix-calculus` |
| $\mu_j$, $\vmu$ | multipliers of the equality constraints $h_j=0$ (free in sign) | {doc}`01-02-matrix-calculus` |
| $q(\vlambda,\vmu)$ | Lagrange dual function | {doc}`02-03-convex-optimization-duality` |
| $d^*$ | optimal value of the dual problem | {doc}`02-03-convex-optimization-duality` |

## Linear programming

| Symbol | Meaning | First used in |
|---|---|---|
| $\mA\in\R^{m\times n}$, $\vb\in\R^m$, $\vc\in\R^n$ | constraint matrix, right-hand side, cost (or profit) vector; $m$ constraints, $n$ variables | {doc}`03-01-lp-model-building` |
| $z$, $z^*$ | objective value $\vc^\T\vx$ of the (primal) LP; its optimal value | {doc}`03-01-lp-model-building` |
| $s_i$ | slack or surplus variable of constraint $i$ | {doc}`03-02-lp-geometry` |
| $x_j^+$, $x_j^-$ | positive and negative part of a free variable | {doc}`03-02-lp-geometry` |
| $P$ | feasible polyhedron $\{\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$ | {doc}`03-02-lp-geometry` |
| $\mB$, $\mN$ | basis matrix and nonbasic matrix | {doc}`01-01-linear-algebra` |
| $\vx_B$, $\vx_N$, $\vc_B$, $\vc_N$ | basic and nonbasic parts of $\vx$ and $\vc$ | {doc}`01-01-linear-algebra` |
| $\mathcal B$, $B(i)$ | index set of the basic variables; basic variable of row $i$ | {doc}`03-03-simplex-method` |
| $\bar\vb=\mB^{-1}\vb$, $\bar a_{ij}$ | updated right-hand side and tableau entries | {doc}`03-03-simplex-method` |
| $\vpi$, $\pi_i$ | simplex multipliers, dual variables, shadow prices: $\vpi^\T=\vc_B^\T\mB^{-1}$ | {doc}`01-01-linear-algebra` |
| $\bar c_j$ | reduced (relative) cost $c_j-\vpi^\T\mA_{\cdot j}$ | {doc}`03-02-lp-geometry` |
| $-z$ row | objective row of a tableau: the $\bar c_j$ and, on the right, $-z$ | {doc}`03-03-simplex-method` |
| $s$, $r$ | entering column and pivot row of a simplex iteration | {doc}`03-03-simplex-method` |
| $\theta$ | step length of a pivot, the minimum ratio $\min\{\bar b_i/\bar a_{is}:\bar a_{is}>0\}$ | {doc}`03-03-simplex-method` |
| $R_i$ | artificial variable of row $i$ | {doc}`03-04-two-phase-big-m` |
| $\phi$ | Phase I objective, the sum of the artificial variables | {doc}`03-04-two-phase-big-m` |
| $M$ | big-M penalty (a large constant) | {doc}`03-04-two-phase-big-m` |
| (P), (D) | primal and dual problem of a pair | {doc}`03-05-lp-duality` |
| $w$, $w^*$ | dual objective value $\vb^\T\vpi$; its optimal value | {doc}`03-05-lp-duality` |
| $I(\vx)$ | index set of the constraints active (binding) at $\vx$ | {doc}`03-05-lp-duality` |
| $\mT$ | bordered inverse $\begin{pmatrix}\mB&\vzero\\\vc_B^\T&1\end{pmatrix}^{-1}$ (the $T$-matrix) | {doc}`01-01-linear-algebra` |
| $\mB_k$, $\mE_k$ | basis and eta matrix of iteration $k$ of the revised simplex method | {doc}`03-06-revised-dual-simplex` |
| $\theta_D$ | dual ratio of the dual simplex method | {doc}`03-06-revised-dual-simplex` |
| $\vl$, $\vu$ | lower and upper bound vectors of the variables | {doc}`04-04-minimum-cost-flow` |

## Networks

| Symbol | Meaning | First used in |
|---|---|---|
| $G=(\mathcal N,\mathcal A)$ | network with node set $\mathcal N$ and arc (edge) set $\mathcal A$ | {doc}`04-01-networks-mst` |
| $(i,j)$ | arc from tail $i$ to head $j$; in an undirected network the edge $i$–$j$ | {doc}`04-01-networks-mst` |
| $\abs{\mathcal N}$, $\abs{\mathcal A}$ | number of nodes, number of arcs | {doc}`04-01-networks-mst` |
| $\mathcal A(i)$ | arcs leaving node $i$ (adjacency list) | {doc}`04-01-networks-mst` |
| $\deg(i)$ | degree of node $i$ | {doc}`04-01-networks-mst` |
| $c_{ij}$ | length or unit cost of arc $(i,j)$ | {doc}`04-01-networks-mst` |
| $\mathcal T$ | a (spanning) tree, as a set of edges | {doc}`04-01-networks-mst` |
| $\mathcal W$, $c(\mathcal W)$ | a cycle and its length or cost | {doc}`04-01-networks-mst` |
| $[S,\bar S]$ | cut between $S$ and $\bar S=\mathcal N\setminus S$ | {doc}`04-01-networks-mst` |
| $[\gamma(j),P(j)]$ | label of Prim's algorithm: cheapest edge cost to the tree, and its end in the tree | {doc}`04-01-networks-mst` |
| $\mA$ | node–arc incidence matrix; the column of arc $(i,j)$ is $\ve_i-\ve_j$ | {doc}`04-01-networks-mst` |
| $s$, $t$ | origin and destination (shortest paths); source and sink (flows) | {doc}`04-02-shortest-paths` |
| $V(j)$, $P(j)$ | shortest-path label: length of a shortest path from $s$ to $j$, and the predecessor of $j$ on it | {doc}`04-02-shortest-paths` |
| $S$, $T$ | sets of nodes with permanent and with temporary labels (Dijkstra) | {doc}`04-02-shortest-paths` |
| $\mD_k=(d^{(k)}_{ij})$, $\mP_k=(p^{(k)}_{ij})$ | Floyd–Warshall distance and predecessor matrices after iteration $k$ | {doc}`04-02-shortest-paths` |
| $\bar c_{ij}$ | reduced length or cost of arc $(i,j)$: $c_{ij}+V(i)-V(j)$ with distance labels, $c_{ij}-\pi_i+\pi_j$ with node potentials | {doc}`04-02-shortest-paths` |
| $\pi_i$ | node potential, the dual variable of the balance equation of node $i$ | {doc}`04-02-shortest-paths` |
| $x_{ij}$ | flow on arc $(i,j)$ | {doc}`04-03-maximum-flow` |
| $u_{ij}$ | capacity (upper bound) of arc $(i,j)$ | {doc}`04-03-maximum-flow` |
| $l_{ij}$ | lower bound of the flow on arc $(i,j)$ | {doc}`04-04-minimum-cost-flow` |
| $v$, $v^*$ | value of a flow; maximum flow value | {doc}`04-03-maximum-flow` |
| $r_{ij}$, $G(\vx)$ | residual capacity; residual network of the flow $\vx$ | {doc}`04-03-maximum-flow` |
| $\mathcal P$, $\Delta$ | augmenting path; amount of flow sent along it (or around a cycle) | {doc}`04-03-maximum-flow` |
| $u[S,\bar S]$ | capacity of an $s$–$t$ cut | {doc}`04-03-maximum-flow` |
| $b_i$ | net supply of node $i$ ($b_i>0$ supply, $b_i<0$ demand) | {doc}`04-04-minimum-cost-flow` |

## Dynamic programming

| Symbol | Meaning | First used in |
|---|---|---|
| $k$, $N$ | stage index; number of stages | {doc}`05-01-dp-fundamentals` |
| $y$, $\mathcal Y_k$ | state; set of states at stage $k$ | {doc}`05-01-dp-fundamentals` |
| $d$, $\mathcal D_k(y)$ | decision; decisions allowed in state $y$ at stage $k$ | {doc}`05-01-dp-fundamentals` |
| $t_k(y,d)$ | transition function: the next state | {doc}`05-01-dp-fundamentals` |
| $c_k(y,d)$, $g(y)$ | stage cost; terminal cost | {doc}`05-01-dp-fundamentals` |
| $V_k(y)$ | optimal value function: best value from state $y$ at stage $k$ to the end | {doc}`05-01-dp-fundamentals` |
| $P_k(y)$ | optimal policy function: the optimal decision(s) in state $y$ at stage $k$ | {doc}`05-01-dp-fundamentals` |
| $F_k(y)$ | forward value: best value from the start to state $y$ at stage $k$ | {doc}`05-01-dp-fundamentals` |
| $(k,y)$, $c^{\uparrow}(k,y)$, $c^{\downarrow}(k,y)$ | vertex of a lattice (column, level); costs of its up and down arcs | {doc}`05-01-dp-fundamentals` |
| $\mathrm U$, $\mathrm D$ | decisions "up" and "down" on a lattice | {doc}`05-01-dp-fundamentals` |
| $\mathrm B$, $\mathrm K$, $\mathrm O$, $\mathrm L$ | decisions buy, keep, overhaul, lease (equipment replacement) | {doc}`05-02-dp-equipment-replacement` |
| $i$; $c(i)$, $p$, $t(i)$, $s(i)$ | age of a machine; operating cost at age $i$, price of a new machine, trade-in value, salvage value | {doc}`05-02-dp-equipment-replacement` |
| $U_k$ | value function of the one-argument replacement model | {doc}`05-02-dp-equipment-replacement` |
| $x$, $X$ | units allocated to an activity; total units of the resource | {doc}`05-03-dp-resource-allocation` |
| $r_k(x)$ | return of activity $k$ from $x$ units | {doc}`05-03-dp-resource-allocation` |
| $V^{(k)}(j)$, $P^{(k)}(j)$ | Bellman–Ford value and predecessor using at most $k$ arcs (procedure 1) | {doc}`05-04-dp-shortest-paths` |
| $\tilde V^{(k)}(j)$, $\hat V^{(k)}(j)$ | values of procedures 2 and 3: shortest walks with at most $k-1$ decreases, reversals | {doc}`05-04-dp-shortest-paths` |
| $D^{(k)}(i,j)$, $P^{(k)}(i,j)$ | Floyd's value: shortest $i$–$j$ distance with intermediate nodes in $\{1,\dots,k\}$; predecessor of $j$ | {doc}`05-04-dp-shortest-paths` |
| $\tilde D^{(k)}(i,j)$ | Dantzig's value: the same with origins and destinations in $\{1,\dots,k\}$ | {doc}`05-04-dp-shortest-paths` |
| $\mathcal S_k$, $i_k$, $T^{(k)}(j)$ | Dijkstra's method as a DP: the $k$ nodes nearest to node 1, the node added at stage $k$, the tentative label of $j$ | {doc}`05-04-dp-shortest-paths` |
| $\bar V(i)$ | length of a shortest path from $i$ to the last node (backward recursion) | {doc}`05-04-dp-shortest-paths` |

## Nonlinear programming

| Symbol | Meaning | First used in |
|---|---|---|
| $f^*$ | value of $f$ at a local or global minimum $\vx^*$ | {doc}`06-01-nlp-models` |
| $\vd$, $\alpha$ | direction; step length along it | {doc}`06-02-optimality-conditions` |
| $D(\vx)$ | cone of feasible directions of $\Omega$ at $\vx$ | {doc}`06-02-optimality-conditions` |
| $L(\vx^*)$, $T(\vx^*)$ | cone of linearized feasible directions; tangent cone of $\Omega$ at $\vx^*$ | {doc}`06-03-kkt-conditions` |
| $\nabla_{\vx}L$, $\nabla^2_{\vx\vx}L$ | gradient and Hessian of the Lagrangian with respect to $\vx$ | {doc}`06-03-kkt-conditions` |
| $C(\vx^*,\vlambda^*)$ | critical cone | {doc}`06-03-kkt-conditions` |
| $\lambda_0$ | multiplier of $f$ in the Fritz John conditions | {doc}`06-03-kkt-conditions` |
| $\vd^{(k)}$, $\alpha_k$ | search direction and step length of iteration $k$: $\vx^{(k+1)}=\vx^{(k)}+\alpha_k\vd^{(k)}$ | {doc}`06-04-line-search` |
| $[a_k,b_k]$ | interval of uncertainty of a one-dimensional search | {doc}`06-04-line-search` |
| $\tau$ | golden-section ratio $(\sqrt5-1)/2\approx0.618$ | {doc}`06-04-line-search` |
| $F_N$ | Fibonacci numbers | {doc}`06-04-line-search` |
| $c_1$, $c_2$ | Armijo and curvature (Wolfe) constants, $0<c_1<c_2<1$ | {doc}`06-04-line-search` |
| $\gamma$ | reduction factor of backtracking | {doc}`06-04-line-search` |
| $\vg^{(k)}$ | gradient $\nabla f(\vx^{(k)})$ at iterate $k$ | {doc}`06-05-unconstrained-algorithms` |
| $\mH(\vx)$, $\bar\mH_k$ | Hessian $\nabla^2 f(\vx)$; a positive definite modification of it | {doc}`06-05-unconstrained-algorithms` |
| $\mD_k$ | positive definite scaling matrix, $\vd^{(k)}=-\mD_k\vg^{(k)}$ (inverse-Hessian approximation of quasi-Newton methods) | {doc}`06-05-unconstrained-algorithms` |
| $\vs^{(k)}$, $\vy^{(k)}$ | step and gradient change of a quasi-Newton update | {doc}`06-05-unconstrained-algorithms` |
| $\Gamma(\vx)$ | Frank–Wolfe gap | {doc}`06-06-constrained-algorithms` |
| $\psi(\vx)$, $\rho$ | penalty function and penalty parameter | {doc}`06-06-constrained-algorithms` |
| $\Phi_\rho(\vx)$, $\vx(\rho)$ | penalized objective $f+\rho\psi$ and its minimizer | {doc}`06-06-constrained-algorithms` |
| $L_A(\vx,\vmu;\rho)$ | augmented Lagrangian | {doc}`06-06-constrained-algorithms` |

## Linearization

| Symbol | Meaning | First used in |
|---|---|---|
| $\delta_i\in\{0,1\}$, $\vdelta$ | binary (indicator) variables and their vector | {doc}`07-01-logical-constraints` |
| $M$, $M_i$ | big-M constants: upper bounds of $g$, $g_i$ on $\Omega_0$ | {doc}`07-01-logical-constraints` |
| $\underline g$, $\underline g_i$ | lower bounds of $g$, $g_i$ on $\Omega_0$ | {doc}`07-01-logical-constraints` |
| $\Omega_0$ | set allowed by all the other constraints of a model | {doc}`07-01-logical-constraints` |
| $F$, $F_i$ | fixed (set-up, opening) cost | {doc}`07-01-logical-constraints` |
| $\omega$, $\omega_{ij}$ | variable that replaces a product ($\delta_A\delta_B$, $\delta x$, $xy$) | {doc}`07-01-logical-constraints` |
| $\underline x$, $\overline x$ | lower and upper bound of a scalar variable | {doc}`07-02-linearization-techniques` |
| $\vx^+$, $\vx^-$ | positive and negative parts of $\vx$ | {doc}`07-02-linearization-techniques` |
| $\eta$ | auxiliary variable that bounds an absolute value, a maximum or a minimum | {doc}`07-02-linearization-techniques` |
| $f^a$ | piecewise-linear approximation of $f$ | {doc}`07-03-piecewise-linear-approximation` |
| $a_1<\dots<a_n$ | breakpoints of a piecewise-linear approximation | {doc}`07-03-piecewise-linear-approximation` |
| $\alpha_i$ | interpolation (SOS2) weights of the breakpoints | {doc}`07-03-piecewise-linear-approximation` |
| $\xi_k$ | filled fraction of segment $k$ (incremental formulation) | {doc}`07-03-piecewise-linear-approximation` |

## Symbols with more than one meaning

A few letters are used deliberately in two roles that never meet in the same
argument. The context, and the arguments of the symbol, always tell them
apart.

| Symbol | Meanings |
|---|---|
| $L$ | the Lagrangian $L(\vx,\vlambda,\vmu)$ ({doc}`01-02-matrix-calculus`, {doc}`02-03-convex-optimization-duality`) and the cone of linearized feasible directions $L(\vx^*)$, with one argument ({doc}`06-03-kkt-conditions`); bold $\mL$ is a triangular factor ({doc}`01-01-linear-algebra`) or the Laplacian matrix of a graph ({doc}`09-04-code-networks`) |
| $\theta$ | the minimum ratio of a simplex pivot ({doc}`03-03-simplex-method`) and the weight of a convex combination ({doc}`02-01-convex-sets`, {doc}`02-02-convex-functions`) |
| $s$, $t$ | source and sink nodes ({doc}`04-02-shortest-paths`–{doc}`04-04-minimum-cost-flow`); $s$ is also the entering column of a pivot ({doc}`03-03-simplex-method`) and $s_i$ a slack variable ({doc}`03-02-lp-geometry`); $t_k(y,d)$ is a DP transition ({doc}`05-01-dp-fundamentals`) and $t(i)$ a trade-in value ({doc}`05-02-dp-equipment-replacement`) |
| $S$, $T$ | permanent and temporary label sets of Dijkstra's algorithm ({doc}`04-02-shortest-paths`); $S$ is also one side of a cut ({doc}`04-01-networks-mst`, {doc}`04-03-maximum-flow`), $T^{(k)}(j)$ a tentative label of Dijkstra's method in DP form ({doc}`05-04-dp-shortest-paths`), $T(\vx^*)$ the tangent cone ({doc}`06-03-kkt-conditions`) and bold $\mT$ the bordered inverse ({doc}`03-06-revised-dual-simplex`) |
| $P$ | predecessor $P(j)$ ({doc}`04-02-shortest-paths`) and optimal policy $P_k(y)$ ({doc}`05-01-dp-fundamentals`), the same idea in two settings; the polyhedron $P$ of {doc}`03-02-lp-geometry`; bold $\mP$ a projection matrix ({doc}`01-01-linear-algebra`) or the Floyd–Warshall predecessor matrix $\mP_k$ ({doc}`04-02-shortest-paths`) |
| $\mD_k$ | the Floyd–Warshall distance matrix ({doc}`04-02-shortest-paths`) and the scaling matrix of a descent method ({doc}`06-05-unconstrained-algorithms`) |
| $D$ | the Jacobian operator $D\vg$ ({doc}`01-02-matrix-calculus`), the leading principal minors $D_k$ ({doc}`01-01-linear-algebra`) and the cone of feasible directions $D(\vx)$ ({doc}`06-02-optimality-conditions`) |
| $\vy^{(k)}$ | the gradient change of a quasi-Newton update ({doc}`06-05-unconstrained-algorithms`) and the vertex returned by the linear subproblem of the Frank–Wolfe method ({doc}`06-06-constrained-algorithms`) |
| $\alpha$ | step length ({doc}`06-02-optimality-conditions`, {doc}`06-04-line-search`) and, with a subscript $\alpha_i$, the SOS2 weights of {doc}`07-03-piecewise-linear-approximation` |
| $\eta$ | the guaranteed return of a maximin model ({doc}`03-01-lp-model-building`), the barrier parameter ({doc}`06-06-constrained-algorithms`) and an auxiliary bounding variable ({doc}`07-02-linearization-techniques`); bold $\veta$ is the eta column ({doc}`01-01-linear-algebra`) |
| $\omega$ | the direction in which a lattice vertex is left ({doc}`05-01-dp-fundamentals`) and a product variable ({doc}`07-01-logical-constraints`, {doc}`07-02-linearization-techniques`) |
| $\gamma$ | the Prim label $\gamma(j)$ ({doc}`04-01-networks-mst`), the backtracking factor ({doc}`06-04-line-search`) and the growth factor of the penalty parameter ({doc}`06-06-constrained-algorithms`) |
| $\delta$ | binary variables $\delta_i$ ({doc}`07-01-logical-constraints`) and, without a subscript, the change of a cost or right-hand side in sensitivity analysis ({doc}`03-07-sensitivity-analysis`) |
| $\kappa$ | the condition number $\kappa(\mA)$ ({doc}`01-01-linear-algebra`, {doc}`06-05-unconstrained-algorithms`), the fractions $\kappa_r$ of the 100 percent rule ({doc}`03-07-sensitivity-analysis`) and the curvature bounds $\kappa_k$ ({doc}`07-03-piecewise-linear-approximation`) |
| $c_1$, $c_2$ | the Armijo and curvature constants ({doc}`06-04-line-search`) and the unit prices of two resources ({doc}`05-03-dp-resource-allocation`) |
