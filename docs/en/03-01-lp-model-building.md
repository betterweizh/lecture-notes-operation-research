# Linear Programming: Model Building

## Overview

A linear program (LP) chooses the values of a few, or many thousands of,
continuous decision variables so that a linear objective is as large or as
small as possible while linear constraints hold. Most of the work in
applying LP lies not in solving the model, which software does reliably,
but in writing it down: deciding what the variables are and in which
units, which limits matter, and what is being optimized. This chapter
teaches that craft. It defines an LP ({prf:ref}`def-03-01-lp`), sets out the
modeling steps, and builds the product-mix model
({prf:ref}`ex-03-01-product-mix`) that chapters
{doc}`03-02-lp-geometry` to {doc}`03-07-sensitivity-analysis` use as their running
example. It then shows the standard patterns: percentage specifications
and ratio constraints, which become linear after multiplying out
({prf:ref}`prop-03-01-ratio`); blending models; multi-period models held
together by balance equations {eq}`eq-03-01-inventory-balance`; objectives
that turn linear after a reformulation, including concave piecewise-linear
profits ({prf:ref}`prop-03-01-concave-split`). It ends with the four
assumptions behind every LP and with the situations in which LP is the
wrong tool. It uses the general optimization problem and the terms
feasible, optimal, bounded and unbounded from {doc}`01-03-introduction`,
and matrix–vector notation from {doc}`01-01-linear-algebra`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- state the general form of an LP and name its components;
- formulate an LP from a verbal description, defining every decision
  variable with its units;
- write percentage specifications, ratio requirements and blending rules
  as linear constraints;
- build multi-period models with inventory, capital or cash balance
  equations and sensible end-of-horizon conditions;
- turn a maximization into a minimization, and model a concave
  piecewise-linear profit (or convex cost) with split variables; and
- check a model against the proportionality, additivity, divisibility and
  certainty assumptions and decide when LP cannot be applied.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $x_j$, $j=1,\dots,n$ | level of activity $j$ (a decision variable, with units) |
| $c_j$ | contribution (profit) or cost per unit of activity $j$ |
| $a_{ij}$ | amount of resource or requirement $i$ used per unit of activity $j$ |
| $b_i$ | amount of resource $i$ available, or requirement $i$ to be met |
| $x_{ij}$ | a doubly indexed variable; its meaning (distillate $i$ in blend $j$, product $i$ in period $j$, …) is defined in each model |
| $p_t$, $d_t$, $I_t$ | production, demand and end-of-period inventory in period $t$ (multi-period models) |

(sec-03-01-what-is-lp)=
## What a linear program is

Every optimization model has three parts: **decision variables**, the
quantities the decision maker controls; an **objective function** of the
variables, to be maximized or minimized; and **constraints**, the
conditions that the variables must satisfy. The numbers that appear in
the model but are not decided (prices, capacities, demands, technical
coefficients) are the **parameters** or data.

:::{prf:definition} Linear program
:label: def-03-01-lp

A **linear program** (LP) is an optimization problem in finitely many real
variables $x_1,\dots,x_n$ whose objective is a linear function and whose
constraints are linear equations or linear inequalities ($\le$ or $\ge$,
never strict):

$$
\begin{aligned}
\max_{\vx}\ \text{or}\ \min_{\vx}\ & z=c_1x_1+c_2x_2+\dots+c_nx_n \\
\st\ & a_{i1}x_1+a_{i2}x_2+\dots+a_{in}x_n\ \{\le,=,\ge\}\ b_i, \qquad i=1,\dots,m.
\end{aligned}
$$ (eq-03-01-general-lp)

Bounds such as $x_j\ge0$ or $x_j\le u_j$ are constraints of the same kind.
:::

A point $\vx$ that satisfies all the constraints is **feasible**; the set of
feasible points is the **feasible region**. A feasible point with the best
objective value is **optimal**, and that value $z^*$ is the **optimal
value**. An LP may have no feasible point (it is **infeasible**), or its
objective may improve without limit over the feasible region (it is
**unbounded**); {doc}`03-02-lp-geometry` shows that these are the only ways an
LP can fail to have an optimal solution ({prf:ref}`thm-03-02-fundamental`).

The adjective "linear" carries the whole theory: each $x_j$ enters every
function to the first power only, multiplied by a constant, and the terms
are added. No products $x_1x_2$, no powers $x_j^2$, no $\abs{x_j}$, no
$\max\{x_1,x_2\}$ and no ratios $x_1/x_2$ may appear. Several of these can
nevertheless be modeled with linear constraints, and this chapter shows
how ({ref}`sec-03-01-specifications`, {ref}`sec-03-01-objectives`).

With vectors, a model with $\le$ constraints and nonnegative variables
reads

$$
\begin{aligned}
\max_{\vx}\ & \vc^\T\vx \\
\st\ & \mA\vx\le\vb \\
& \vx\ge\vzero,
\end{aligned}
$$ (eq-03-01-matrix-form)

where $\vc=(c_1,\dots,c_n)^\T$, $\mA=(a_{ij})\in\R^{m\times n}$ and
$\vb=(b_1,\dots,b_m)^\T$; inequalities between vectors hold component by
component. Column $\mA_{\cdot j}$ lists what one unit of activity $j$
uses of every resource; row $\mA_{i\cdot}$ lists how every activity uses
resource $i$. Every LP can be rewritten in the **standard form** $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$ used by
the simplex method ({prf:ref}`prop-03-02-standard-form`).

(sec-03-01-modelling-steps)=
## The modeling process

A verbal problem description rarely lists its variables and constraints.
A reliable way to extract them is to ask four questions, in this order.

1. **Activities.** What can the decision maker do, and at what level?
   Producing a product, shipping on a route, buying a raw material,
   depositing money in an account, starting a work shift: each activity
   gets one decision variable, its level.
2. **Resources and requirements.** What do the activities consume, and what
   must they achieve? Machine hours, storage space, budgets and raw
   materials are limited; demands, quality specifications and service
   levels must be met.
3. **Constraints.** For every resource, the total use by all activities
   cannot exceed the amount available; for every requirement, the total
   contribution of all activities must reach the required level; balance
   conditions link quantities that flow from one activity or period to the
   next.
4. **Objective.** Which single measure judges a plan: profit, cost,
   revenue, number of staff, final wealth?

Then the model is written down.

- **Define each decision variable in words, with its units**, for example
  "$x_1$ = number of hundreds of cases of juice glasses produced per
  week". The units fix every coefficient: a coefficient is "units of the
  row quantity per unit of the variable". Many modeling errors are unit
  errors, and they are found by checking that every term of a constraint
  has the units of its right-hand side.
- **Write the objective** as the sum over activities of (contribution or
  cost per unit) × (level).
- **Write the constraints**, one family at a time, each with a short name
  (capacity, demand, balance, specification). Do not forget the
  **nonnegativity constraints** $x_j\ge0$; activity levels are rarely
  allowed to be negative, and a solver does not know that unless told.
- **Check the model**: on a trivial plan (all variables zero), on a plan
  you can evaluate by hand, and on the solver's answer. A model that is
  unbounded usually lacks a constraint; one that is infeasible usually has
  a sign or unit error, or data that cannot all be satisfied.

(sec-03-01-product-mix)=
## A product-mix model

The first model is the book's running example. Chapters
{doc}`03-02-lp-geometry` to {doc}`03-07-sensitivity-analysis` solve it
graphically, with the simplex method, through its dual and with
sensitivity analysis.

:::{prf:example} Product mix of a glassware maker
:label: ex-03-01-product-mix

A small molding shop runs one injection-molding machine, on which it can
mount one of two molds. The first mold makes 200 cc juice glasses and
the second 350 cc cocktail glasses; the data per *hundred cases* of each
product are given in the table. The shop plans 60 machine hours per week
and keeps a week's output in a store of 150 m³. At most 800 cases of juice
glasses can be sold per week; cocktail glasses sell without limit. How
many cases of each glass should be made per week to maximize the total
contribution?

| | Juice glasses | Cocktail glasses | Available |
|---|---|---|---|
| Machine time (h per 100 cases) | 6 | 5 | 60 h/week |
| Storage (m³ per 100 cases) | 10 | 20 | 150 m³ |
| Contribution (\$ per 100 cases) | 500 | 450 | |
| Demand limit (cases/week) | 800 | none | |
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Activities and resources.** The two activities are producing juice
glasses and producing cocktail glasses. They use two resources, machine
time and storage space, and juice-glass sales are limited by demand. The
objective is total weekly contribution.

**Decision variables.**

- $x_1$ = juice glasses produced, in hundreds of cases per week;
- $x_2$ = cocktail glasses produced, in hundreds of cases per week.

Measuring in hundreds of cases keeps the coefficients small: \$5 per case
becomes \$500 per hundred cases, 0.1 m³ per case becomes 10 m³ per hundred
cases, and the demand limit of 800 cases becomes $x_1\le8$.

**Model.**

$$
\begin{aligned}
\max_{x_1,x_2}\ & z=500x_1+450x_2 && \text{(contribution, dollars/week)} \\
\st\ & 6x_1+5x_2\le60 && \text{(machine time, h/week)} \\
& 10x_1+20x_2\le150 && \text{(storage, m}^3\text{)} \\
& x_1\le8 && \text{(demand, hundreds of cases/week)} \\
& x_1,\ x_2\ge0. &&
\end{aligned}
$$ (eq-03-01-product-mix)

**Solution.** The optimal plan, found graphically in
{prf:ref}`ex-03-02-graphical` and by the simplex method in
{doc}`03-03-simplex-method`, is

$$
\begin{gathered}
x_1^*=\tfrac{45}{7}\approx6.429,\\
x_2^*=\tfrac{30}{7}\approx4.286,\\
z^*=\tfrac{36000}{7}\approx5142.86 ,
\end{gathered}
$$

that is, about 643 cases of juice glasses and 429 cases of cocktail
glasses per week, for a contribution of \$5142.86 per week. Both the
machine ($6\cdot\tfrac{45}{7}+5\cdot\tfrac{30}{7}=60$ h) and the store
($10\cdot\tfrac{45}{7}+20\cdot\tfrac{30}{7}=150$ m³) are fully used; the
demand limit is not reached, since $8-\tfrac{45}{7}=\tfrac{11}{7}$ hundred,
about 157 cases per week, of juice-glass demand stays unserved.
:::

As the solution of {prf:ref}`ex-03-01-product-mix` shows, the product
with the higher contribution per case, juice glasses, is *not* produced up
to its demand limit. Neither product is better on every count: per machine hour,
cocktail glasses earn more ($450/5=90$ against $500/6\approx83.3$ \$/h),
while per cubic meter of storage, juice glasses earn more ($500/10=50$
against $450/20=22.5$ \$/m³). The optimum mixes the two products so that
both scarce resources are used up. Such trade-offs are hard to see by
inspection and are exactly what LP resolves.

(sec-03-01-resource-allocation)=
### The resource-allocation model

{prf:ref}`ex-03-01-product-mix` is an instance of the general
**resource-allocation** (product-mix) model: $n$ activities compete for
$m$ limited resources,

$$
\begin{aligned}
\max_{\vx}\ & z=\sum_{j=1}^n c_jx_j \\
\st\ & \sum_{j=1}^n a_{ij}x_j\le b_i, \qquad i=1,\dots,m, \\
& x_j\ge0, \qquad j=1,\dots,n,
\end{aligned}
$$ (eq-03-01-resource-allocation)

where $c_j$ is the contribution per unit of activity $j$, $b_i$ the amount
of resource $i$ available and $a_{ij}$ the amount of resource $i$ used per
unit of activity $j$. Demand limits are resources too ($x_1\le8$ has
$a_{31}=1$, $a_{32}=0$, $b_3=8$). Much of LP's vocabulary comes from this
model: a constraint that holds with equality at the optimum is
**binding**, the unused amount of a resource is its **slack**
({doc}`03-02-lp-geometry`), and the value of one more unit of a resource is
its **shadow price** ({ref}`sec-03-05-shadow-prices`).
{doc}`08-01-examples-lp-modeling` works through product-mix models with
three products and a production ratio
({prf:ref}`ex-08-01-three-models`) and a toll-revenue model in which the
"resource" is the green time of a traffic light
({prf:ref}`ex-08-01-traffic-light`).

(sec-03-01-specifications)=
## Specifications, percentages and ratios

Many products must meet quality specifications stated as percentages: an
alloy with at least 0.45 % manganese, a fuel blend with at most 60 % of
one component, a bag of hardware with at least 25 % bolts by weight. How
such a specification is written depends on whether the total quantity is
fixed or is itself a decision.

**Fixed total.** If the total amount is known, the percentage is a fixed
amount. A foundry that must deliver exactly 1000 kg of castings with at
least 0.45 % manganese needs at least $0.0045\cdot1000=4.5$ kg of
manganese, and the specification is the linear constraint "manganese
contributed by all inputs $\ge4.5$ kg". If input $j$ contains a fraction
$q_j$ of manganese and $x_j$ kg of it is used, the constraint reads
$\sum_j q_jx_j\ge4.5$. This is how the foundry model
({prf:ref}`ex-08-01-foundry`) writes its manganese and silicon limits.

**Variable total.** If the total is a decision, the specification is a
ratio of variables. Suppose a blend is made of components with amounts
$x_1,\dots,x_k\ge0$ and component 1 may make up at most a fraction $p$
of it:

$$
\frac{x_1}{x_1+x_2+\dots+x_k}\le p .
$$ (eq-03-01-ratio-spec)

This is not linear, but it becomes linear after multiplying by the
denominator.

:::{prf:proposition} Linearizing a ratio specification
:label: prop-03-01-ratio

Let $x_1,\dots,x_k\ge0$ with $s=x_1+\dots+x_k>0$ and $0\le p\le1$. Then
{eq}`eq-03-01-ratio-spec` holds if and only if

$$
(1-p)\,x_1-p\,x_2-\dots-p\,x_k\le0 .
$$ (eq-03-01-ratio-linear)

The same holds with $\ge$ in place of $\le$ throughout.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Because $s>0$, multiplying both sides of $x_1/s\le p$ by $s$ preserves the
inequality, so {eq}`eq-03-01-ratio-spec` is equivalent to $x_1\le p\,s=
p\,x_1+p\,x_2+\dots+p\,x_k$. Moving all terms to the left gives
{eq}`eq-03-01-ratio-linear`. Every step is reversible, and with $\ge$ the
same argument applies.
:::

When the blend may be empty ($s=0$), the ratio is undefined but
{eq}`eq-03-01-ratio-linear` holds with equality at $\vx=\vzero$, which is
what one wants: an unproduced blend violates no specification. The
multiplication must never be done with a denominator that can be
*negative*, because that would reverse the inequality.

For example, "component A is at most 60 % of the deluxe blend" with
$x_{A1},x_{B1},x_{C1}$ liters of the three components becomes
$0.4x_{A1}-0.6x_{B1}-0.6x_{C1}\le0$, and "component C is at least 20 % of
it" becomes $-0.2x_{A1}-0.2x_{B1}+0.8x_{C1}\ge0$.

**Fixed proportions.** A requirement that products be made in the ratio
$3:2:5$ is a pair of equations between the variables. With $x_1,x_2,x_3$
the quantities, $x_1:x_2:x_3=3:2:5$ means $x_1/3=x_2/2=x_3/5$, that is

$$
\begin{aligned}
2x_1-3x_2&=0,\\
5x_2-2x_3&=0 .
\end{aligned}
$$

Equivalently, introduce one variable $k\ge0$ (the number of "batches" of
the mix) and set $x_1=3k$, $x_2=2k$, $x_3=5k$; the model then has one
decision instead of three ({prf:ref}`ex-08-01-three-models`).

**Counting by weight.** Specifications sometimes mix weights and counts.
If a bolt weighs ten times as much as a nut and the number of bolts in a
bag may not exceed the number of nuts, then with $x_b$, $x_n$ the weights
(lb) and $\beta$ the weight of one bolt, the counts are $x_b/\beta$ and
$10x_n/\beta$, and the constraint is $x_b\le10x_n$. The unknown unit weight
cancels, as it must ({prf:ref}`ex-08-01-hardware-bag`).

(sec-03-01-max-min)=
### Maximization and minimization

Maximizing $f$ and minimizing $-f$ over the same feasible set $\Omega$
have the same optimal solutions, and the
optimal values differ only in sign:

$$
\max_{\vx\in\Omega} f(\vx)=-\min_{\vx\in\Omega}\bigl(-f(\vx)\bigr).
$$ (eq-03-01-max-min)

Constants in the objective do not affect the optimal solution either. A
foundry that has a firm order for 1000 kg of castings at \$0.85 per kg
earns a fixed revenue of \$850, so

$$
\max\ (850-\text{cost}(\vx))\quad\Longleftrightarrow\quad\min\ \text{cost}(\vx)
$$

over the same constraints; the two models have the same optimal
inputs, and the maximum profit equals \$850 minus the minimum cost
({prf:ref}`ex-08-01-foundry`). Whether a problem is stated as a maximization
or as a minimization is therefore a matter of convenience. The simplex
method of {doc}`03-03-simplex-method` is written for one sense and handles the
other through {eq}`eq-03-01-max-min`.

(sec-03-01-blending)=
## Blending models

In a **blending** problem, raw materials (crude oils, distillates, ores,
feeds) are mixed into products that must meet specifications. The natural
variables are doubly indexed:

$$
x_{ij}=\text{amount of input } i \text{ used in product } j ,
$$

so that $\sum_j x_{ij}$ is the total use of input $i$ and $\sum_i x_{ij}$
the total amount of product $j$ made (when volumes add, as they do for
liters of miscible liquids). A blending model then has

- **supply constraints** $\sum_j x_{ij}\le S_i$ for each input $i$;
- **specification constraints**, one per specification, written with
  {prf:ref}`prop-03-01-ratio`, for example "input $i$ is at most a fraction
  $p_{ij}$ of product $j$":
  $x_{ij}\le p_{ij}\sum_{h}x_{hj}$;
- **demand constraints**, if any, on $\sum_i x_{ij}$;
- the objective "revenue minus cost",
  $\sum_j r_j\sum_i x_{ij}-\sum_i k_i\sum_j x_{ij}=\sum_{i,j}(r_j-k_i)x_{ij}$,
  with selling price $r_j$ per unit of product $j$ and purchase cost $k_i$
  per unit of input $i$.

The last identity shows how to read a blending model: every unit of input
$i$ placed in product $j$ earns the margin $r_j-k_i$, and the
specifications limit how freely the inputs can be placed. If the most
expensive product can absorb everything that is available without
violating its specifications, the optimum is to make only that product,
as happens in the distillate-blending example of
{prf:ref}`ex-08-01-distillate-blending`. When the specifications bind, the
model trades off higher prices against tighter recipes. The foundry model
({prf:ref}`ex-08-01-foundry`) and the hardware bag
({prf:ref}`ex-08-01-hardware-bag`) are blending models with a single product;
the refinery models ({prf:ref}`ex-08-01-oilco`, {prf:ref}`ex-08-01-refinery`)
are the reverse situation, in which one input splits into several products
in *fixed* proportions.

(sec-03-01-multi-period)=
## Multi-period models

Decisions in one period often affect the next: goods produced now can be
sold later, money invested now returns later, equipment busy now is free
later. A **multi-period** (dynamic) LP indexes its variables by period
$t=1,\dots,T$ and links consecutive periods by **balance equations** of
the form

$$
\begin{aligned}
&\text{(amount at the start of } t\text{)}+\text{(inflow during } t\text{)} \\
&\quad=\text{(outflow during } t\text{)}+\text{(amount carried to } t+1\text{)} .
\end{aligned}
$$

The simplest case is a single product with production $p_t$, demand $d_t$
and end-of-period inventory $I_t$ (with $I_0$ given):

$$
\begin{aligned}
I_t&=I_{t-1}+p_t-d_t, && t=1,\dots,T,\\
I_t&\ge0, && t=1,\dots,T .
\end{aligned}
$$ (eq-03-01-inventory-balance)

The condition $I_t\ge0$ says that demand must be met on time (no
backlogging): the cumulative production up to $t$ covers the cumulative
demand up to $t$. Adding capacity limits $p_t\le\bar p_t$ and the cost
$\sum_t(k_tp_t+h_tI_t)$, with unit production cost $k_t$ and holding cost
$h_t$ per unit carried for one period, gives the classical
**production–inventory model**. {numref}`fig-03-01-inventory-flow` shows its
structure: each period is a node, and the balance equation says that what
enters the node equals what leaves it.

```{figure} figures/03-01-inventory-flow.svg
:name: fig-03-01-inventory-flow
:width: 560px
:class: dark-light
:alt: Three period nodes in a row; production arrows enter each node from above, demand arrows leave each node downwards, and inventory arrows link each node to the next

Flow view of the balance equations {eq}`eq-03-01-inventory-balance` for three
periods. Production $p_t$ enters period $t$, demand $d_t$ leaves it, and the
inventory $I_t$ carries goods to period $t+1$. In every node, inflow equals
outflow.
```

The same pattern covers many applications.

- **Cash flow.** In an investment plan, "money at the start of year $t$" is
  the balance quantity: returns of earlier investments and bank deposits
  with interest flow in; new investments and deposits flow out
  ({prf:ref}`ex-08-01-project-investment`, {prf:ref}`ex-08-01-two-plans`).
- **Delayed activities.** When an activity started in period $t$ delivers
  its output in period $t+2$, the variable appears with a minus sign in the
  balance of period $t$ (as a use) and with a plus sign in the balance of
  period $t+2$ (as a supply). The capital-goods model of
  {prf:ref}`ex-08-01-capital-consumer` has two such activities.
- **Several products** sharing one capacity add one balance per product and
  one capacity constraint per period
  ({prf:ref}`ex-08-01-production-inventory`).
- **Projects** executed over several years, with yearly budgets, use the
  fraction of each project completed in each year as variables
  ({prf:ref}`ex-08-01-renewal-projects`).

**The horizon.** A model with $T$ periods ignores what happens afterwards,
and an optimizer exploits this: without a further condition it ends the
horizon with nothing in stock, nothing invested and all equipment sold, a
behavior called the **end-of-horizon effect**. A good model therefore
states what the final state is worth, either by requiring a final
inventory ($I_T\ge\bar I$) or by valuing whatever remains at the end at a
terminal price. The capital-goods example values unused capital and
consumer goods at the prices of the period after the horizon. The choice
matters: different terminal values lead to different plans in the last
periods, and they should be chosen to reflect what the remaining goods are
really worth.

(sec-03-01-objectives)=
## Objectives that become linear

Some objectives that are not linear as stated can be replaced by an
equivalent linear model. Three cases arise in this chapter's problems;
{doc}`07-02-linearization-techniques` treats them in general.

**Maximin.** To maximize the smallest of several linear functions,
$\max_{\vx}\min_k \ell_k(\vx)$, introduce a free variable $\eta$ and solve
$\max\{\eta: \eta\le \ell_k(\vx)\ \text{for all}\ k,\ \vx\in\Omega\}$. At the
optimum $\eta$ equals the smallest $\ell_k$, because $\eta$ is pushed up until one
of the constraints $\eta\le \ell_k(\vx)$ binds. A gambler who maximizes the
worst-case return of a bet is the classical instance
({prf:ref}`ex-08-01-max-min-bets`).

**Absolute values in a minimization.** To minimize $\sum_i\abs{e_i}$, where
each $e_i$ is a linear function of the variables, either bound each term
from above, $u_i\ge e_i$ and $u_i\ge-e_i$, and minimize $\sum_iu_i$; or
split $e_i=e_i^+-e_i^-$ with $e_i^\pm\ge0$ and minimize $\sum_i(e_i^++e_i^-)$.
({ref}`sec-07-02-abs` treats both in general.)
Fitting a straight line by least absolute deviations is an LP for this
reason ({prf:ref}`ex-08-01-lad-regression`). Both tricks work only when the
objective pushes the absolute values *down*; maximizing
$\abs{e_i}$ is not an LP.

**Concave profit or convex cost.** A profit whose marginal rate decreases
(or a cost whose marginal rate increases) is a piecewise-linear function
that the next section models with split variables.

(sec-03-01-piecewise)=
## Piecewise-linear profits and costs

The product-mix model assumes that each case earns the same contribution.
Suppose instead that the juice glasses earn \$500 per hundred cases for the
first 300 cases per week, but only \$300 per hundred cases beyond that,
because the extra cases must be sold at a discount. The contribution of
juice glasses is then the piecewise-linear function

$$
r(x_1)=\begin{cases}500x_1, & 0\le x_1\le3,\\ 1500+300(x_1-3), & x_1>3,\end{cases}
$$

shown in {numref}`fig-03-01-piecewise-profit`. It is concave: its slope
decreases from 500 to 300 at the breakpoint $x_1=3$.

```{figure} figures/03-01-piecewise-profit.svg
:name: fig-03-01-piecewise-profit
:width: 480px
:class: dark-light
:alt: A piecewise-linear increasing function of x1 with slope 500 up to x1 = 3 and slope 300 beyond, with the linear function 500 x1 dashed for comparison

Contribution of juice glasses with a quantity discount: slope 500 \$ per
hundred cases up to $x_1=3$, slope 300 beyond. The dashed line is the
constant-rate contribution $500x_1$ of {prf:ref}`ex-03-01-product-mix`.
```

:::{prf:example} Product mix with a concave piecewise-linear profit
:label: ex-03-01-piecewise-profit

How can the product-mix model {eq}`eq-03-01-product-mix` be changed so that it
uses the contribution $r(x_1)$ and remains an LP, and what is its optimal
plan?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Split variables.** Write the juice-glass output as the sum of the part
sold at the full rate and the part sold at the discounted rate:

$$
\begin{gathered}
x_1=x_{11}+x_{12},\\
0\le x_{11}\le3,\\
x_{12}\ge0 ,
\end{gathered}
$$

both in hundreds of cases per week. Replace $x_1$ by $x_{11}+x_{12}$
everywhere, and value each part at its own rate:

$$
\begin{aligned}
\max\ & z=500x_{11}+300x_{12}+450x_2 \\
\st\ & 6(x_{11}+x_{12})+5x_2\le60 \\
& 10(x_{11}+x_{12})+20x_2\le150 \\
& x_{11}+x_{12}\le8 \\
& 0\le x_{11}\le3 \\
& x_{12},\ x_2\ge0 .
\end{aligned}
$$ (eq-03-01-piecewise-lp)

**Solution.** The optimum (computed with the HiGHS solver; the simplex
method of {doc}`03-03-simplex-method` gives the same) is

$$
\begin{gathered}
x_{11}^*=3,\\
x_{12}^*=\tfrac{24}{7}\approx3.429,\\
x_2^*=\tfrac{30}{7}\approx4.286,\\
z^*=\tfrac{31200}{7}\approx4457.14 .
\end{gathered}
$$

The plan is the same as before, $x_1^*=3+\tfrac{24}{7}=\tfrac{45}{7}$
hundred cases of juice glasses and $\tfrac{30}{7}$ of cocktail glasses;
only the contribution falls, from \$5142.86 to \$4457.14 per week, because
$\tfrac{24}{7}$ hundred cases now earn \$300 instead of \$500. The full-rate
part is used up ($x_{11}^*=3$) before any discounted cases are counted.
:::

Nothing in the split model {eq}`eq-03-01-piecewise-lp` says that $x_{11}$ must reach its
bound 3 before $x_{12}$ becomes positive; a plan such as $x_{11}=1$,
$x_{12}=2$ is feasible but makes no sense, since it sells cases at the
discount while full-price sales are still possible. Such plans are never
optimal, and the reason is concavity.

:::{prf:proposition} Split variables for concave profits and convex costs
:label: prop-03-01-concave-split

Let a variable $x\ge0$ be split as $x=\sum_{k=1}^K y_k$ with
$0\le y_k\le\Delta_k$ ($\Delta_k$ the length of segment $k$, $\Delta_K$ possibly
infinite), and let segment $k$ earn $r_k$ per unit with
$r_1>r_2>\dots>r_K$ (a concave profit in a maximization). Then every
optimal solution of the split LP fills the segments in order: if
$y_k>0$, then $y_i=\Delta_i$ for all $i<k$. The same holds for a
minimization with costs $r_1<r_2<\dots<r_K$ (a convex cost).
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Take a feasible solution with $y_k>0$ and $y_i<\Delta_i$ for some $i<k$. Let
$\varepsilon=\min\{y_k,\ \Delta_i-y_i\}>0$ and move $\varepsilon$ units from
segment $k$ to segment $i$. The sum $x=\sum_ky_k$ is unchanged, so every
constraint in which the split variables appear only through $x$ still
holds, and the bounds on $y_i$, $y_k$ hold by the choice of
$\varepsilon$. The objective changes by $(r_i-r_k)\varepsilon>0$ in the
maximization (by $(r_i-r_k)\varepsilon<0$ in the minimization), so the
original solution was not optimal.
:::

The argument fails when the rates *increase* in a maximization (a convex
profit, as with a volume bonus) or decrease in a minimization (a concave
cost, as with economies of scale). If the juice-glass rates were \$300 for
the first 300 cases and \$500 beyond, the split LP would set $x_{11}=0$,
$x_{12}=\tfrac{45}{7}$ and report \$5142.86, crediting every case with the
bonus rate, although the true contribution of that plan is only
$900+500\cdot\tfrac{24}{7}+450\cdot\tfrac{30}{7}=\tfrac{31800}{7}\approx
\$4542.86$. Non-concave piecewise-linear functions need binary variables
that enforce the order of the segments; {doc}`07-03-piecewise-linear-approximation`
develops these formulations, and {doc}`02-02-convex-functions` explains why
convexity is the dividing line.

(sec-03-01-assumptions)=
## Assumptions of linear programming

An LP is a faithful model of a decision problem only if four assumptions
hold, at least approximately.

1. **Proportionality.** The contribution of an activity to the objective
   and its use of each resource are proportional to its level. If one
   kilogram of product $j$ needs 3 kg of raw material $i$, then $x_j$ kg
   need $3x_j$ kg, for every $x_j\ge0$; if one unit earns $c_j$, then
   $x_j$ units earn $c_jx_j$. There are no set-up costs, no quantity
   discounts and no economies of scale.
2. **Additivity.** The total use of each resource, and the total value of
   the objective, are the sums of the contributions of the individual
   activities. If product 1 uses 2 units of a resource per unit and
   product 2 uses 3, then $x_1$ and $x_2$ units together use $2x_1+3x_2$.
   There are no interactions between activities, that is, no
   cross-product terms such as $x_1x_2$.
3. **Divisibility** (continuity). Every variable may take any real value
   in its range; fractional activity levels, such as $\tfrac{45}{7}$
   hundred cases, are meaningful.
4. **Certainty.** All parameters $c_j$, $a_{ij}$, $b_i$ are known
   constants.

Proportionality and additivity together are the **linearity
assumptions**: they make the objective and the constraint functions
linear. Divisibility makes the feasible region a continuum rather than a
set of lattice points, which is what allows the geometry of
{doc}`03-02-lp-geometry` and the simplex method. Certainty makes a single
optimal plan meaningful.

(sec-03-01-when-lp-fails)=
### When LP cannot be applied

LP is used so widely because large models can be solved fast and
reliably, and because the dual solution prices the resources. Its
assumptions nevertheless fail in many situations, and then a different
model is needed.

- **Fixed costs and economies of scale** violate proportionality. Opening
  a depot, setting up a production run or building a road has a cost that
  does not depend on how much the facility is used ({ref}`sec-07-01-fixed-charge`);
  a cost that grows
  more slowly than the quantity (a concave cost) cannot be modeled with
  the split variables of {prf:ref}`prop-03-01-concave-split`. These need
  mixed-integer models ({doc}`07-01-logical-constraints`,
  {doc}`07-03-piecewise-linear-approximation`).
- **Congestion and other nonlinear physics.** The travel time on a road
  grows faster than linearly with the traffic volume, the pressure loss in
  a pipe grows with the square of the flow, revenue is price times
  quantity when the price depends on the quantity sold. Such models are
  nonlinear programs ({doc}`06-01-nlp-models`). When the nonlinearity is
  convex in the right direction, a piecewise-linear approximation keeps the
  model an LP.
- **Interactions** between activities violate additivity, for example two
  treatments whose combined effect differs from the sum of their separate
  effects; a product $x_1x_2$ then appears ({doc}`07-02-linearization-techniques`
  shows the cases that can still be linearized).
- **Indivisible decisions** violate divisibility: the number of buses on a
  route, the number of staff on a shift, whether a project is built at all.
  If the numbers are large, rounding an LP solution is often good enough;
  if they are small, rounding can be infeasible or far from optimal, and an
  integer program is needed. Rounding is *not* the same as solving the
  integer program: in the product-mix example, rounding the LP optimum
  down gives $(6,4)$ with \$4800 per week, while the best integer plan is
  $(8,2)$ with \$4900 ({prf:ref}`ex-03-02-graphical`).
- **Logical conditions** such as "either this constraint or that one",
  "if A is built then B must be built" or "at most two of five projects"
  make the feasible region nonconvex; they need binary variables
  ({doc}`07-01-logical-constraints`).
- **Uncertainty** violates certainty. Demands, prices and travel times are
  forecasts. When the uncertainty is mild, sensitivity analysis
  ({doc}`03-07-sensitivity-analysis`) tells how far the data may move before
  the optimal plan changes; when it is essential, stochastic or robust
  optimization is required. The maximin model of
  {prf:ref}`ex-08-01-max-min-bets` is a simple robust model.
- **Several conflicting objectives** (cost, emissions, equity) cannot be
  combined into one linear function without choosing weights, which is a
  value judgment that the model cannot make.

Even when an assumption fails, an LP is often a useful first model: it is
quick to build, quick to solve, and its dual and sensitivity information
show which data and assumptions matter most.

## Summary and Common Mistakes

### Summary

- An LP has a linear objective and finitely many linear equality and
  inequality constraints ({prf:ref}`def-03-01-lp`, {eq}`eq-03-01-general-lp`);
  in matrix form $\max\{\vc^\T\vx:\mA\vx\le\vb,\ \vx\ge\vzero\}$
  {eq}`eq-03-01-matrix-form`.
- Modeling starts from activities, resources, constraints and the
  objective; every variable is defined in words with its units
  ({ref}`sec-03-01-modelling-steps`).
- The product-mix LP {eq}`eq-03-01-product-mix` is the book's running
  example; its optimum is $\vx^*=(\tfrac{45}{7},\tfrac{30}{7})$ hundred cases
  per week with $z^*=\tfrac{36000}{7}\approx\$5142.86$ per week
  ({prf:ref}`ex-03-01-product-mix`).
- A ratio specification with a positive denominator is linear after
  multiplying out ({prf:ref}`prop-03-01-ratio`); fixed proportions are pairs
  of equations; max and min are interchangeable {eq}`eq-03-01-max-min`.
- Blending models use variables $x_{ij}$ for input $i$ in product $j$;
  multi-period models are held together by balance equations
  {eq}`eq-03-01-inventory-balance` ({numref}`fig-03-01-inventory-flow`) and need
  an explicit end-of-horizon condition.
- Maximin objectives and sums of absolute values in a minimization become
  linear with auxiliary variables; concave profits and convex costs become
  linear with split variables ({prf:ref}`ex-03-01-piecewise-profit`,
  {prf:ref}`prop-03-01-concave-split`).
- LP rests on proportionality, additivity, divisibility and certainty
  ({ref}`sec-03-01-assumptions`); fixed costs, congestion, interactions,
  indivisibility, logical conditions and uncertainty call for integer,
  nonlinear or stochastic models ({ref}`sec-03-01-when-lp-fails`).

### Common mistakes

- Defining a variable without units, or mixing units in one constraint
  (cases with hundreds of cases, kg with tonnes). Read every term of a
  constraint in the units of its right-hand side.
- Forgetting the nonnegativity constraints, or imposing them on a
  variable that may legitimately be negative (a net return, a free
  intercept, a deviation).
- Writing a percentage specification as $x_1\le p$ (a fraction compared
  with an amount) instead of $x_1\le p\sum_jx_j$.
- Multiplying a ratio constraint by a denominator that can be negative.
- Leaving out a balance equation in a multi-period model, or letting the
  model end with no value for what remains at the horizon.
- Using split variables for a convex profit or a concave cost: the LP then
  fills the most profitable segment first and overstates the objective.
- Rounding a fractional LP solution and calling it the integer optimum.

## References

- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Ch. 3, the linear programming model.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Ch. 3, §§3.1–3.4: prototype example,
  the LP model, assumptions of LP, additional examples.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Ch. 2, §§2.1–2.4: two-variable LP model, LP applications;
  ch. 3, §3.1: LP model in equation form.)
- Williams, H. P. (2013). *Model building in mathematical programming*
  (5th ed.). Wiley. (Chs. 1–5: building linear programming models,
  including blending, multi-period and structured models; ch. 7:
  nonlinear models and separable programming.)
