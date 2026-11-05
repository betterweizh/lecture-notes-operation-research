# Introduction to Operations Research

## Overview

Operations research (OR) is the discipline of making decisions with
mathematical models: a real decision problem is written as an
**optimization model** with decision variables, an objective and
constraints, the model is solved by an algorithm, and the answer is
checked against reality before it is used. This chapter closes Part I
(Preliminaries): after the mathematical tools of {doc}`01-01-linear-algebra`
and {doc}`01-02-matrix-calculus`, it explains what OR is and where it came from, shows where it is used, with examples from the
planning and operation of transport, water and energy infrastructure, and
describes the modeling process from problem definition to implementation.
It then fixes the vocabulary used in the rest of the book: the general
form {eq}`eq-01-03-general-form` of an optimization problem, feasible and
optimal solutions, optimal values, infeasible and unbounded problems, and
local and global optima. The chapter ends with the classes of problems the
book treats and a map of its nine parts ({numref}`fig-01-03-roadmap`) with
suggestions on how to read them.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- explain what operations research is and outline its history;
- describe the phases of an OR study, from problem definition to
  implementation;
- identify the decision variables, parameters, objective and constraints
  of a verbal decision problem and write it in the general form
  {eq}`eq-01-03-general-form`;
- use the terms feasible, infeasible, unbounded, optimal solution, optimal
  value, local and global optimum correctly; and
- name the problem classes of this book and find the parts and chapters
  that treat them.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\vx\in\R^n$ | vector of decision variables |
| $f$ | objective function |
| $g_i$, $i=1,\dots,m$ | inequality constraint functions, $g_i(\vx)\le0$ |
| $h_j$, $j=1,\dots,\ell$ | equality constraint functions, $h_j(\vx)=0$ |
| $\Omega$ | feasible set |
| $\vx^*$ | optimal solution |
| $p^*$ | optimal value of a minimization problem, $p^*=\inf\{f(\vx):\vx\in\Omega\}$ |

(sec-01-03-what-is-or)=
## What operations research is

Every organization that runs a system must decide how to use resources
that are limited: money, staff, vehicles, machines, energy, water, time.
A transit agency decides how many buses to run on each route and when; a
water utility decides how much water to take from each source and how to
schedule its pumps; an electricity system operator decides which power
plants to switch on and how much each should produce in every hour; a road
authority decides which pavements to resurface this year and which can
wait. Such decisions have three features in common. There are many
possible choices, often infinitely many. The choices are restricted by
physical, technical, financial or legal limits. And some choices are
clearly better than others with respect to a measurable goal such as cost,
travel time, reliability or emissions.

**Operations research** is the application of mathematical models and
algorithms to such decisions. The name comes from its origin in the study
of military *operations* (see {ref}`sec-01-03-history`); the terms
*management science* and *decision science* are used for the same field,
and *operational research* is the British name. Three ideas characterize
the OR approach:

1. **Modeling.** The decision problem is described by a mathematical model
   that keeps the essential trade-offs and leaves out details that do not
   change the decision. The model makes the assumptions explicit, so that
   they can be discussed and tested.
2. **Optimization.** Among all choices allowed by the model, OR looks for
   one that is best with respect to a stated objective, or at least
   provably close to the best. This distinguishes OR from simulation and
   from rules of thumb, which evaluate a given plan but do not search for
   a better one.
3. **Algorithms.** The models of practice have thousands or millions of
   variables. Solving them requires algorithms whose correctness can be
   proved and whose effort grows reasonably with the size of the problem,
   implemented in software (solvers).

The questions OR asks of a model are therefore of three kinds: *where* the
optimal solutions are (the geometry and the optimality conditions of the
problem), *how* to find them (algorithms), and *what they mean* (the
sensitivity of the answer to the data and the economic interpretation of
the results, for example as prices of scarce resources). The parts of
this book answer these three questions for each class of models in turn.

The focus of OR is on decisions about the **operation and management** of
systems, which is why infrastructure is a natural field of application.
Transport networks, water supply and drainage systems, power grids,
ports and telecommunication networks are large, expensive and shared by
many users; small improvements in how they are planned and operated save
large amounts of money, energy and time.

(sec-01-03-history)=
## A short history

Mathematical ideas that we now count as OR are older than the name.
Agner Krarup Erlang analyzed telephone traffic in Copenhagen around 1909 and
founded queueing theory in order to size telephone exchanges, an early
infrastructure planning problem. In 1939 Leonid Kantorovich formulated
production-planning problems as what we now call linear programs, and in
1941 Frank Hitchcock stated the transportation problem
(Hillier & Lieberman, 2021, Chapter 1).

The name itself dates from the late 1930s, when British scientists were
asked to study how the new radar stations could best be used within the air
defense system. During the Second World War, mixed teams of physicists,
mathematicians, biologists and engineers analyzed military operations: the
deployment of radar and anti-aircraft guns, the size of merchant convoys,
the search patterns of anti-submarine patrols. Their work was called
*operational research* in Britain and *operations research* in the United
States, where similar groups were set up in the Navy and the Army Air
Forces. The teams were valued because they treated operations as systems,
collected data, built simple quantitative models and tested their
recommendations in the field.

After the war, the methods were turned toward industry and government, and
the field acquired its mathematical core.

- In 1947 George B. Dantzig, working on planning problems of the US Air
  Force, formulated the general **linear programming** (LP) problem and
  invented the **simplex method** to solve it; in the same year John von
  Neumann sketched the duality theory of LP. Dantzig's account of these
  origins is in Dantzig (1963). Linear programming is the subject of
  {doc}`03-01-lp-model-building` to {doc}`03-07-sensitivity-analysis`.
- In 1951 Harold W. Kuhn and Albert W. Tucker published optimality
  conditions for **nonlinear programming** (Kuhn & Tucker, 1951). William
  Karush had derived the same conditions in his 1939 master's thesis
  (Karush, 1939), and they now carry the names of all three. These
  **Karush–Kuhn–Tucker (KKT) conditions**
  appear in {doc}`02-03-convex-optimization-duality` and {doc}`06-03-kkt-conditions`.
- In the early 1950s Richard Bellman developed **dynamic programming** for
  multistage decision problems, based on his principle of optimality
  (Bellman, 1957); see {doc}`05-01-dp-fundamentals`.
- Network optimization took shape in the same decade: the maximum-flow
  algorithm of Ford and Fulkerson (1956), the minimum spanning
  tree algorithms of Kruskal (1956) and Prim (1957), and Dijkstra's shortest-path
  algorithm (Dijkstra, 1959) ({doc}`04-01-networks-mst` to
  {doc}`04-04-minimum-cost-flow`).

Professional societies and journals followed (the Operations Research
Society of America was founded in 1952), and with the spread of computers
OR became part of engineering, business and public administration.
Later milestones include integer-programming methods (cutting planes and
branch and bound, around 1960), interior-point methods for LP in the 1980s,
and, since the 1990s, commercial and open-source solvers that solve
models with millions of variables on a desktop computer
({doc}`09-01-code-setup-solvers`). Today OR methods are also part of machine
learning, where training a model is an optimization problem, and of the
digital twins used to operate infrastructure in real time.

(sec-01-03-applications)=
## Where operations research is used

OR is used wherever decisions are repeated often, involve large sums, or
must coordinate many interacting parts. {numref}`tab-01-03-applications` lists
typical decisions in infrastructure systems together with the model classes
that describe them and the chapters of this book that treat those classes.

```{list-table} Typical OR decisions in infrastructure systems
:name: tab-01-03-applications
:header-rows: 1
:widths: 18 42 22 18

* - Sector
  - Decision
  - Model class
  - Chapters
* - Transport
  - Routes and timetables of buses and trains; crew rostering
  - LP, integer programming
  - {doc}`3.1 <03-01-lp-model-building>`, {doc}`7.1 <07-01-logical-constraints>`–{doc}`7.2 <07-02-linearization-techniques>`
* - Transport
  - Fastest or cheapest routes for vehicles, emergency services and freight
  - Shortest paths
  - {doc}`4.2 <04-02-shortest-paths>`, {doc}`5.4 <05-04-dp-shortest-paths>`
* - Transport
  - Traffic assignment: how traffic distributes over a congested road
    network, and which links to add or widen
  - Convex and nonlinear programming
  - {doc}`2.2 <02-02-convex-functions>`–{doc}`2.3 <02-03-convex-optimization-duality>`, {doc}`6.1 <06-01-nlp-models>`–{doc}`6.6 <06-06-constrained-algorithms>`
* - Water
  - Daily draw from reservoirs, groundwater and desalination; pump
    schedules
  - LP, nonlinear programming
  - {doc}`3.1 <03-01-lp-model-building>`, {doc}`6.1 <06-01-nlp-models>`
* - Water
  - Release policy of a reservoir over a season with uncertain inflows
  - Dynamic programming
  - {doc}`5.1 <05-01-dp-fundamentals>`–{doc}`5.3 <05-03-dp-resource-allocation>`
* - Energy
  - Output of each power plant to meet demand at least cost (economic
    dispatch); which plants to switch on (unit commitment)
  - Convex programming, mixed-integer programming
  - {doc}`2.3 <02-03-convex-optimization-duality>`, {doc}`7.1 <07-01-logical-constraints>`–{doc}`7.3 <07-03-piecewise-linear-approximation>`
* - Energy and water
  - Capacity of pipes, cables and transmission lines in a network
  - Network flows, mixed-integer programming
  - {doc}`4.3 <04-03-maximum-flow>`–{doc}`4.4 <04-04-minimum-cost-flow>`, {doc}`7.1 <07-01-logical-constraints>`
* - All sectors
  - When to repair or replace pavements, bridges, pumps, vehicles
  - Dynamic programming
  - {doc}`5.2 <05-02-dp-equipment-replacement>`
* - All sectors
  - Allocation of a limited budget to projects
  - Integer and dynamic programming
  - {doc}`5.3 <05-03-dp-resource-allocation>`, {doc}`7.1 <07-01-logical-constraints>`
```

Outside infrastructure, the same models plan production and inventories in
manufacturing, route parcels in logistics, schedule operating theaters in
hospitals, allocate portfolios in finance and assign frequencies in
telecommunications. Hillier and Lieberman (2021, Chapter 1) and Taha (2017,
Chapter 1) describe many such applications. The traffic-assignment models
mentioned in the table are treated in depth by Sheffi (1985).

(sec-01-03-modelling-process)=
## The modeling process

An OR study is more than solving an equation. It is usually described in
five phases ({numref}`fig-01-03-modelling-process`).

```{figure} figures/01-03-modelling-process.svg
:name: fig-01-03-modelling-process
:width: 680px
:class: dark-light
:alt: Five boxes from problem definition to implementation connected by arrows, with feedback arrows from validation back to formulation and to problem definition

The phases of an OR study. Validation often sends the analyst back to the
formulation or even to the definition of the problem; the process is a
cycle rather than a straight line.
```

1. **Problem definition.** Find out what is to be decided, by whom, with
   which goal and within which limits. The decision maker's first
   description is rarely precise: "reduce the operating cost of the water
   network" may mean energy cost of pumping, chemical cost of treatment,
   or both, and may or may not allow changes to the reservoirs' operating
   rules. This phase fixes the scope, the objective and the alternatives.
2. **Model formulation.** Translate the problem into mathematics: choose
   the decision variables, express the objective and the constraints as
   functions of them, and collect the data (parameters). The art lies in
   choosing the level of detail: enough to capture the trade-offs, not so
   much that the model cannot be solved or its data cannot be obtained.
   Every model rests on assumptions, such as proportionality and
   additivity in an LP ({ref}`sec-03-01-assumptions`), and these
   assumptions should be written down.
3. **Solution.** Solve the model with a suitable algorithm, usually with a
   solver. Parts II–VII of this book explain the algorithms, because
   knowing how they work tells you which models can be solved, how large they may
   be, and what the output means (for example, dual values and ranges of
   optimality, {doc}`03-05-lp-duality` and {doc}`03-07-sensitivity-analysis`).
4. **Validation.** Check whether the model's answer makes sense. Does the
   model reproduce past decisions when it is given past data? Does the
   solution change in a plausible way when the data change? Are all
   constraints that matter in reality present? A solution that looks too
   good, or a model that is unbounded or infeasible (see
   {ref}`sec-01-03-basic-terms`), usually means that the formulation has an
   error.
5. **Implementation.** Put the decision into practice, which needs the
   trust of the people who will use it, and maintain the model as the
   system changes. Many models are run routinely, for example a dispatch
   model that is solved every few minutes.

A shorter version of the same cycle, used throughout this book, has four
steps: state the decision problem; build the model; design or choose a
solution algorithm; and test the model and algorithm on examples and case
studies. The book stresses the second and third steps, model building and
algorithm design, because these are the skills that transfer from one
application to the next.

(sec-01-03-model-parts)=
## The parts of an optimization model

Every optimization model in this book has four parts.

- **Decision variables** are the quantities the decision maker controls,
  collected in a vector $\vx=(x_1,\dots,x_n)^\T\in\R^n$: amounts of water
  drawn from each source, power output of each generator, number of buses
  on each route, binary choices such as "build this link or not".
- **Parameters** (data) are the quantities the decision maker cannot change
  in the model: demands, capacities, unit costs, travel times. They are
  treated as known constants; how the optimal decision depends on them is
  the subject of sensitivity analysis.
- The **objective function** $f(\vx)$ measures the quality of a decision:
  a cost or travel time to be minimized, or a profit or throughput to be
  maximized.
- **Constraints** describe which decisions are allowed: physical limits
  (capacities, conservation of flow), technical rules (pressure limits,
  safety margins), policies and budgets.

:::{prf:example} Supplying a district with water
:label: ex-01-03-water-supply

A water utility must deliver at least 120 ML (megaliters) of treated water
per day to a district. It can draw water from a reservoir, at a pumping and
treatment cost of \$250 per ML, but the reservoir's safe yield allows at
most 80 ML/day. It can also buy water from a desalination plant at \$1 100
per ML, up to the plant's capacity of 100 ML/day. How much should it take
from each source? Formulate the decision as an optimization model, naming
its decision variables, parameters, objective and constraints, and solve it.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Decision variables.* $x_1$ = water drawn from the reservoir and $x_2$ =
water bought from the desalination plant, both in ML/day.

*Parameters.* The demand of 120 ML/day, the capacities 80 and 100 ML/day,
and the unit costs \$250 and \$1 100 per ML.

*Objective.* The daily cost $f(\vx)=250x_1+1100x_2$ (\$/day), to be
minimized.

*Constraints.* Demand must be met, $x_1+x_2\ge120$; the capacities
$x_1\le80$ and $x_2\le100$; and the amounts cannot be negative,
$x_1,x_2\ge0$. The model is

$$
\begin{aligned}
\min_{\vx}\ & 250x_1+1100x_2\\
\st\ & x_1+x_2\ge120,\\
& x_1\le80,\\
& x_2\le100,\\
& x_1,x_2\ge0.
\end{aligned}
$$ (eq-01-03-water-model)

*Optimal plan.* The reservoir is cheaper, so it should be used to its limit and
the desalination plant should cover the rest: $\vx^*=(80,40)$ with cost
$250\cdot80+1100\cdot40=\$64{,}000$ per day. The intuition can be turned
into a proof. For every feasible $\vx$, the demand constraint gives
$x_2\ge120-x_1$, so

$$
\begin{aligned}
250x_1+1100x_2
&\ge 250x_1+1100\,(120-x_1)\\
&= 132{,}000-850\,x_1\\
&\ge 132{,}000-850\cdot80\\
&= 64{,}000,
\end{aligned}
$$

where the last inequality uses $x_1\le80$. No feasible plan costs less than
\$64,000 per day, and $\vx^*$ attains this bound.
:::

The model {eq}`eq-01-03-water-model` is a linear program: the objective and
all constraints are linear in $\vx$. Small as it is, it shows the typical
structure: a cheap but limited resource, an expensive but plentiful one,
and a requirement that forces a trade-off. The argument at the end of the
solution of {prf:ref}`ex-01-03-water-supply`, which combines constraints
into a lower bound on the objective, is the germ of duality theory
({doc}`02-03-convex-optimization-duality`, {doc}`03-05-lp-duality`).

(sec-01-03-general-form)=
## The general form of an optimization problem

Write the constraints of a model as $g_i(\vx)\le0$ for $i=1,\dots,m$
(inequalities) and $h_j(\vx)=0$ for $j=1,\dots,\ell$ (equalities). The set
of decisions that satisfy all of them,

$$
\begin{aligned}
\Omega=\bigl\{\vx\in\R^n:\ & g_i(\vx)\le0,\ i=1,\dots,m,\\
& h_j(\vx)=0,\ j=1,\dots,\ell\bigr\},
\end{aligned}
$$ (eq-01-03-feasible-set)

is the **feasible set** (or feasible region). The general optimization
problem is

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)\\
\st\ & \vx\in\Omega .
\end{aligned}
$$ (eq-01-03-general-form)

Every constraint can be brought into one of the two forms: $a(\vx)\ge
b(\vx)$ becomes $b(\vx)-a(\vx)\le0$. In the water model
{eq}`eq-01-03-water-model`, $g_1(\vx)=120-x_1-x_2$, $g_2(\vx)=x_1-80$,
$g_3(\vx)=x_2-100$, $g_4(\vx)=-x_1$ and $g_5(\vx)=-x_2$, with $m=5$ and
$\ell=0$.

Maximization needs no separate theory, because

$$
\max_{\vx\in\Omega} f(\vx) = -\min_{\vx\in\Omega}\bigl(-f(\vx)\bigr),
$$ (eq-01-03-max-min)

and both problems have the same optimal solutions. This book writes the
general theory for minimization; LP chapters state the rules for both
senses where it helps.

(sec-01-03-basic-terms)=
## Basic terms

The following terms are used throughout the book.

:::{prf:definition} Feasible and optimal solutions
:label: def-01-03-feasible-optimal

Consider the problem {eq}`eq-01-03-general-form`.

1. A point $\vx\in\Omega$ is a **feasible solution**; a point
   $\vx\notin\Omega$ is **infeasible**.
2. The problem is **infeasible** if $\Omega=\emptyset$.
3. The **optimal value** is $p^*=\inf\{f(\vx):\vx\in\Omega\}$, with the
   conventions $p^*=+\infty$ if the problem is infeasible and
   $p^*=-\infty$ if $f$ takes arbitrarily small values on $\Omega$; in the
   latter case the problem is **unbounded**.
4. A feasible $\vx^*$ with $f(\vx^*)=p^*$ is an **optimal solution** (a
   minimizer). The set of all optimal solutions may be empty, a single
   point, or infinitely many points.
:::

The distinction between the optimal *solution* $\vx^*$ (a decision) and the
optimal *value* $p^*$ (a number) matters. In {prf:ref}`ex-01-03-water-supply`
the optimal solution is $\vx^*=(80,40)$ ML/day and the optimal value is
$p^*=\$64{,}000$ per day. If the desalination cost were also \$250 per ML,
the optimal value would be $\$30{,}000$ per day and every plan with
$x_1+x_2=120$, $20\le x_1\le80$ would be optimal: one value, infinitely many
optimal solutions (alternative optima, {ref}`sec-03-03-special-cases`).

Three things can go wrong, and each has a typical cause in practice.

- **Infeasibility.** If the district needed 200 ML/day, the constraints
  $x_1+x_2\ge200$, $x_1\le80$, $x_2\le100$ could not hold together, since
  the two sources supply at most 180 ML/day. The model then has no
  feasible solution. In practice an infeasible model means that the
  requirements are too strict or that the data contain an error.
- **Unboundedness.** Suppose the model also let the utility sell surplus
  water to a neighboring district at \$1 300 per ML, and the modeler
  forgot the capacity of the desalination plant. Each extra ML bought for
  \$1,100 and sold for \$1 300 would lower the cost by \$200, without
  limit, and $p^*=-\infty$. Real systems are never unbounded; an unbounded
  model has a missing constraint.
- **No attained optimum.** The infimum can be finite without being
  attained. Minimizing $f(x)=1/x$ over $\Omega=\{x\in\R:x\ge1\}$ gives
  $p^*=0$, but $f(x)>0$ for every feasible $x$, so no optimal solution
  exists. Weierstrass' theorem ({prf:ref}`thm-02-01-weierstrass`) gives
  conditions (a continuous objective on a nonempty, closed and bounded
  feasible set) under which a minimum always exists.

Finally, a solution can be best among its neighbors without being best
overall.

:::{prf:definition} Local and global minima
:label: def-01-03-local-global

A feasible point $\vx^*$ is a **global minimum** of {eq}`eq-01-03-general-form`
if $f(\vx^*)\le f(\vx)$ for all $\vx\in\Omega$. It is a **local minimum** if
there is an $\varepsilon>0$ such that $f(\vx^*)\le f(\vx)$ for all
$\vx\in\Omega$ with $\norm{\vx-\vx^*}<\varepsilon$.
:::

Every global minimum is a local one, but not conversely. If the
desalination plant had strong economies of scale, so that its cost per ML
fell steeply with the amount bought, the cost function would no longer be
linear and the model could have a local minimum that uses mostly reservoir
water and a different, global minimum that uses mostly desalinated water.
Algorithms that only improve a solution step by step can stop at the
wrong one. For **convex** problems this cannot happen: every local minimum
is global ({prf:ref}`thm-02-03-local-global`). This is the main reason
why convexity ({doc}`02-01-convex-sets`, {doc}`02-02-convex-functions`) is
studied early in this book, in Part II, right after the preliminaries.
Local and global optima of general nonlinear programs are the subject of
{doc}`06-01-nlp-models` and {doc}`06-02-optimality-conditions`.

(sec-01-03-problem-classes)=
## Problem classes in this book

The general form {eq}`eq-01-03-general-form` covers every model in this book.
What distinguishes the classes is the type of $f$, $g_i$, $h_j$ and of the
variables, and this type decides which algorithms apply.

```{list-table} Problem classes treated in this book
:name: tab-01-03-problem-classes
:header-rows: 1
:widths: 24 46 30

* - Class
  - Structure
  - Main methods (part, chapters)
* - Convex optimization
  - $f$, $g_i$ convex, $h_j$ affine
  - Convexity, Lagrange duality, KKT conditions (Part II, {doc}`2.1 <02-01-convex-sets>`–{doc}`2.3 <02-03-convex-optimization-duality>`)
* - Linear programming (LP)
  - $f$, $g_i$, $h_j$ linear (affine); continuous variables
  - Simplex, duality, sensitivity (Part III, {doc}`3.1 <03-01-lp-model-building>`–{doc}`3.7 <03-07-sensitivity-analysis>`)
* - Network optimization
  - LPs on a graph: flow conservation at nodes, arc capacities and costs
  - Spanning trees, shortest paths, maximum flow, minimum-cost flow (Part IV, {doc}`4.1 <04-01-networks-mst>`–{doc}`4.4 <04-04-minimum-cost-flow>`)
* - Dynamic programming (DP)
  - Decisions in stages; the objective adds up over the stages
  - Recursion on value functions (Part V, {doc}`5.1 <05-01-dp-fundamentals>`–{doc}`5.4 <05-04-dp-shortest-paths>`)
* - Nonlinear programming (NLP)
  - $f$, $g_i$, $h_j$ smooth, possibly nonconvex
  - Optimality conditions, line search, descent and penalty methods (Part VI, {doc}`6.1 <06-01-nlp-models>`–{doc}`6.6 <06-06-constrained-algorithms>`)
* - Mixed-integer programming (MIP)
  - Some variables restricted to integers or to $\{0,1\}$
  - Modeling with binary variables, linearization (Part VII, {doc}`7.1 <07-01-logical-constraints>`–{doc}`7.3 <07-03-piecewise-linear-approximation>`)
```

The classes overlap. An LP is a convex problem; a network flow problem is
an LP with special structure that faster algorithms exploit; a shortest-path
problem can be solved as an LP, by Dijkstra's method or by dynamic
programming; and nonlinear or logical relations can often be turned into a
mixed-integer *linear* model ({doc}`07-02-linearization-techniques`,
{doc}`07-03-piecewise-linear-approximation`), which a MIP solver can handle.

(sec-01-03-roadmap)=
## A map of the book

{numref}`fig-01-03-roadmap` shows the nine parts of the book and how they
depend on each other.

```{figure} figures/01-03-roadmap.svg
:name: fig-01-03-roadmap
:width: 760px
:class: dark-light
:alt: Block diagram of the book: a band for Part I preliminaries at the top; six boxes for Parts II to VII (convex optimization, linear programming, network optimization, dynamic programming, nonlinear programming, linearization) joined by dependency arrows; bands for Part VIII examples and problems and Part IX implementation at the bottom

The nine parts of the book and their dependencies. A solid arrow means that
the later part uses the earlier one; a dashed line joins parts that treat
the same problem from two sides (KKT conditions derived from LP duality;
shortest paths by labeling and by dynamic programming). Part VIII collects
worked examples and problems for Parts III–VI, Part IX implements the
methods in Python.
```

**Part I** (Preliminaries) collects the tools used later: matrices and
linear systems, eigenvalues and definiteness of quadratic forms
({doc}`01-01-linear-algebra`), gradients, Hessians and Taylor expansions
({doc}`01-02-matrix-calculus`), and, in this chapter, the vocabulary of
optimization models. The first two chapters can be read first or consulted
when needed; this chapter is worth reading before any later part.

**Part II** (Convex Optimization), {doc}`02-01-convex-sets` to
{doc}`02-03-convex-optimization-duality`, develops convex sets, convex
functions and Lagrange duality. It gives the language for the geometry of
LP ({doc}`03-02-lp-geometry`), for LP duality ({doc}`03-05-lp-duality`) and
for the optimality conditions of nonlinear programs
({doc}`06-02-optimality-conditions`, {doc}`06-03-kkt-conditions`).

**Part III** (Linear Programming), {doc}`03-01-lp-model-building` to
{doc}`03-07-sensitivity-analysis`, builds LP models, solves them by the
simplex method and its variants, and analyzes duality and sensitivity.
**Part IV** (Network Optimization), {doc}`04-01-networks-mst` to
{doc}`04-04-minimum-cost-flow`, uses LP theory but has its own, faster
algorithms. **Part V** (Dynamic Programming), {doc}`05-01-dp-fundamentals`
to {doc}`05-04-dp-shortest-paths`, needs little beyond Part I and can be
read at any time; its last chapter revisits {doc}`04-02-shortest-paths`.
**Part VI** (Nonlinear Programming), {doc}`06-01-nlp-models` to
{doc}`06-06-constrained-algorithms`, builds on convexity and duality.
**Part VII** (Linearization), {doc}`07-01-logical-constraints` to
{doc}`07-03-piecewise-linear-approximation`, shows how to model logical
conditions and nonlinear functions with binary variables, so that a
mixed-integer linear solver can handle them.

**Part VIII** (Examples and Problems) contains, for Parts III–VI, further
worked examples and problems with folded solutions. The theory chapters
keep only the examples that are needed to develop a method, also with
folded solutions; the rest are in Part VIII. In both places, try an example
yourself before you open its solution. **Part IX** (Implementation) shows
how to solve the models with Python and how to program the algorithms of
Parts III–VI; the code reproduces the numbers of the worked examples.

A reader interested mainly in linear models can follow Parts III and IV
({doc}`03-01-lp-model-building` to {doc}`04-04-minimum-cost-flow`) and then
Part VII ({doc}`07-01-logical-constraints` to
{doc}`07-03-piecewise-linear-approximation`). A reader interested mainly in
nonlinear models should read Part II ({doc}`02-01-convex-sets` to
{doc}`02-03-convex-optimization-duality`), then Part VI
({doc}`06-01-nlp-models` to {doc}`06-06-constrained-algorithms`), with
{doc}`03-05-lp-duality` as background for the KKT conditions. Either way,
the matching chapters of Parts VIII and IX can be read alongside each
theory chapter.

## Summary

- Operations research uses mathematical models and algorithms to make
  decisions about systems with limited resources; infrastructure systems
  (transport, water, energy) are a major field of application
  ({numref}`tab-01-03-applications`).
- Its mathematical core was formed between 1939 and 1960: linear
  programming and the simplex method (1947), the KKT conditions (1939,
  1951), dynamic programming and network algorithms (1950s).
- An OR study passes through problem definition, model formulation,
  solution, validation and implementation, usually more than once
  ({numref}`fig-01-03-modelling-process`).
- A model consists of decision variables, parameters, an objective and
  constraints; every model in the book has the general form
  {eq}`eq-01-03-general-form` with feasible set {eq}`eq-01-03-feasible-set`, and
  maximization reduces to minimization by {eq}`eq-01-03-max-min`.
- A problem can be infeasible, unbounded, have a finite optimal value that
  is not attained, or have one or many optimal solutions
  ({prf:ref}`def-01-03-feasible-optimal`); local minima need not be global
  ({prf:ref}`def-01-03-local-global`), except in convex problems.
- The book has nine parts: preliminaries (Part I); convex optimization, LP,
  network optimization, DP, nonlinear programming and linearization for
  mixed-integer models (Parts II–VII); worked examples and problems
  (Part VIII); and the implementation in Python (Part IX)
  ({numref}`tab-01-03-problem-classes`, {numref}`fig-01-03-roadmap`).

## References

- Bellman, R. (1957). *Dynamic programming*. Princeton University Press.
  (Preface and Chapter 3, the principle of optimality.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Chapter 1,
  mathematical optimization and its problem classes.)
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Chapters 1–2, the origins of linear programming.)
- Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
  *Numerische Mathematik, 1*(1), 269–271. <https://doi.org/10.1007/BF01386390>
- Ford, L. R., & Fulkerson, D. R. (1956). Maximal flow through a network.
  *Canadian Journal of Mathematics, 8*, 399–404.
  <https://doi.org/10.4153/CJM-1956-045-5>
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 1, the origins and nature of
  OR; Chapter 2, overview of the modeling approach.)
- Karush, W. (1939). *Minima of functions of several variables with
  inequalities as side constraints* [Master's thesis, University of
  Chicago].
- Kruskal, J. B. (1956). On the shortest spanning subtree of a graph and
  the traveling salesman problem. *Proceedings of the American Mathematical
  Society, 7*(1), 48–50. <https://doi.org/10.1090/S0002-9939-1956-0078686-7>
- Kuhn, H. W., & Tucker, A. W. (1951). Nonlinear programming. In J. Neyman
  (Ed.), *Proceedings of the Second Berkeley Symposium on Mathematical
  Statistics and Probability* (pp. 481–492). University of California
  Press.
- Prim, R. C. (1957). Shortest connection networks and some
  generalizations. *Bell System Technical Journal, 36*(6), 1389–1401.
  <https://doi.org/10.1002/j.1538-7305.1957.tb01515.x>
- Sheffi, Y. (1985). *Urban transportation networks: Equilibrium analysis
  with mathematical programming methods*. Prentice-Hall. (Chapter 1,
  network representation and equilibrium.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 1, what is operations research.)
- Williams, H. P. (2013). *Model building in mathematical programming*
  (5th ed.). Wiley. (Chapter 1, the concept of a model and mathematical
  programming models.)
