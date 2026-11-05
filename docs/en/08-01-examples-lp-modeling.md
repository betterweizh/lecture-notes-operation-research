# LP Modeling and Geometry

This chapter works through the examples and problems of the first two
chapters on linear programming (LP). The first five sections, on
{doc}`03-01-lp-model-building`, are modeling problems, grouped by the
pattern they use; the last three, on {doc}`03-02-lp-geometry`, practice the
graphical method, standard form, basic solutions and bases. Each item
states the problem; the folded solution gives the model
(decision variables with units, objective, constraints), the optimal
solution and value, and a short interpretation. All optima were computed
with the HiGHS LP solver and checked by hand or with exact fractions; the
simplex method of {doc}`03-03-simplex-method` gives the same results.

(sec-08-01-product-mix)=
## Product mix and resource allocation

:::{prf:example} Three models with a fixed production ratio
:label: ex-08-01-three-models

A manufacturer makes three models, I, II and III, of a product from two
raw materials. One unit of model I, II and III needs 2, 3 and 5 units of
raw material A (4000 units available) and 4, 2 and 7 units of raw material
B (6000 units available). A unit of model I takes twice as much labor as a
unit of model II and three times as much as a unit of model III; the whole
workforce could make 1500 units of model I if it made nothing else. At
least 200, 200 and 150 units of models I, II and III must be made, and the
market requires the three models in the ratio $3:2:5$. The selling prices
are \$30, \$20 and \$50 per unit. Formulate an LP that maximizes revenue
and solve it.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $x_1,x_2,x_3$ = units of models I, II, III produced.

**Labor.** Let $t$ be the labor time of one unit of model I; models II
and III take $t/2$ and $t/3$, and the workforce supplies $1500t$. Dividing
by $t$, the labor constraint is $x_1+\tfrac12x_2+\tfrac13x_3\le1500$
(in "model-I equivalents").

**Model.**

$$
\begin{aligned}
\max\ & z=30x_1+20x_2+50x_3 \\
\st\ & 2x_1+3x_2+5x_3\le4000 && \text{(raw material A)} \\
& 4x_1+2x_2+7x_3\le6000 && \text{(raw material B)} \\
& x_1+\tfrac12x_2+\tfrac13x_3\le1500 && \text{(labor)} \\
& 2x_1-3x_2=0 && \text{(ratio } 3:2:5\text{)} \\
& 5x_2-2x_3=0 && \text{(ratio } 3:2:5\text{)} \\
& x_1\ge200 && \text{(minimum output)} \\
& x_2\ge200 \\
& x_3\ge150 .
\end{aligned}
$$

**Solution.** The ratio equations mean $x_1=3k$, $x_2=2k$, $x_3=5k$ for
some $k\ge0$ (the number of "batches" of the mix). One batch uses
$6+6+25=37$ units of A, $12+4+35=51$ units of B and $3+1+\tfrac53=\tfrac{17}{3}$
labor units, and earns $90+40+250=\$380$. The tightest limit is raw
material A, $37k\le4000$, so $k=\tfrac{4000}{37}$ and

$$
\begin{aligned}
x_1^*&=\tfrac{12000}{37}\approx324.32,\\
x_2^*&=\tfrac{8000}{37}\approx216.22,\\
x_3^*&=\tfrac{20000}{37}\approx540.54,\\
z^*&=\tfrac{1520000}{37}\approx\$41\,081.08 .
\end{aligned}
$$

The minimum demands hold ($324.32\ge200$, $216.22\ge200$, $540.54\ge150$);
raw material B ($5513.5$ of 6000) and labor ($612.6$ of 1500) have
slack.

**Interpretation.** Once the ratio is fixed, the LP has a single degree of
freedom, and the scarcest resource, raw material A, decides the output. If
whole units are required, $k=108$ batches give $(324,216,540)$, which keeps
the exact ratio, uses 3996 units of A and earns \$41 040.
:::

:::{prf:example} Green times of a traffic light
:label: ex-08-01-traffic-light

Traffic from three highways H1, H2 and H3 waits at a signal before
entering a toll road. The tolls are \$4, \$5 and \$6 per car from H1, H2 and
H3, and the arrival rates are 550, 650 and 450 cars per hour. The signal
cycle may last at most 2.2 minutes; each highway gets a green phase of at
least 22 s, followed by a yellow phase of 10 s during which no car moves.
The toll gate can process at most 500 cars per hour of green time; for its
capacity per cycle, count the green time of a cycle of maximum length. Assume
that, while its light is green, a highway discharges cars at its arrival
rate. Find the green times that maximize the toll revenue per cycle.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $g_1,g_2,g_3$ = green time of H1, H2, H3 (s per cycle).

**Constraints.** The cycle is at most $2.2\cdot60=132$ s, of which
$3\cdot10=30$ s are yellow, so $g_1+g_2+g_3\le102$. A highway with green
time $g_i$ sends $(\nu_i/3600)g_i$ cars per cycle, with $\nu_i$ its rate in
cars per hour. The gate can process $500/3600$ cars per second during at
most 102 s of green, that is, at most $500\cdot102/3600=14.17$ cars per
cycle.

**Model.**

$$
\begin{aligned}
\max\ & z=4\cdot\tfrac{550}{3600}g_1+5\cdot\tfrac{650}{3600}g_2+6\cdot\tfrac{450}{3600}g_3
&& \text{(dollars per cycle)} \\
\st\ & \tfrac{550}{3600}g_1+\tfrac{650}{3600}g_2+\tfrac{450}{3600}g_3\le\tfrac{500}{3600}\cdot102
&& \text{(toll gate, cars per cycle)} \\
& g_1+g_2+g_3\le102 && \text{(cycle length, s)} \\
& g_1,\ g_2,\ g_3\ge22 . &&
\end{aligned}
$$

**Solution.**

$$
\begin{gathered}
g_1^*=22\ \text{s},\\
g_2^*=22\ \text{s},\\
g_3^*=\tfrac{164}{3}\approx54.67\ \text{s},\\
z^*=\tfrac{2675}{36}\approx\$74.31\ \text{per cycle}.
\end{gathered}
$$

The toll gate is binding ($14.17$ cars per cycle), the cycle constraint is
not (the cycle lasts $98.67+30=128.67$ s), and H1 and H2 get only their
minimum green. Over an hour this is about $74.31\cdot3600/128.67\approx
\$2079$.

**Interpretation.** The scarce resource is the capacity of the gate, so
each unit of it should go to the highway whose cars pay the highest toll,
H3. Without the gate limit the answer would change: the cycle constraint
would bind and the green time would go to H2, which earns the most per
second of green ($5\cdot650=3250$ against 2200 and 2700 dollars per hour of
green), giving $(22,58,22)$ s and \$82.31 per cycle.

*(An earlier version of this solution gives $g=(22,\,45.1,\,34.9)$ s and
\$80.33 "per hour". That plan sends 15.87 cars per cycle to a gate that can
take only 14.17, so it is infeasible, and the objective of the model is
revenue per cycle, not per hour.)*
:::

:::{prf:example} Planning an urban renewal project
:label: ex-08-01-land-use-plan

A city plans an urban renewal project with lower-, middle- and
upper-income housing, public housing, an elementary school and retail
units. The number of classrooms must be enough for the pupils (at most 25
pupils per room), and the retail capacity must cover the retail demand
generated by the housing and the school. The data are:

| | Lower income | Middle income | Upper income | Public housing | School room | Retail unit |
|---|---|---|---|---|---|---|
| Minimum number of units | 100 | 125 | 75 | 300 | — | 0 |
| Maximum number of units | 200 | 190 | 260 | 600 | — | 25 |
| Land per unit (acre) | 0.05 | 0.07 | 0.03 | 0.025 | 0.045 | 0.1 |
| Pupils per unit | 1.3 | 1.2 | 0.5 | 1.4 | — | — |
| Retail demand per unit | 0.023 | 0.034 | 0.046 | 0.023 | 0.034 | — |
| Annual income per unit (\$) | 7000 | 12 000 | 20 000 | 5000 | — | 15 000 |

The school may occupy at most 2 acres, and each classroom costs \$10 000
per year to run. The project uses a vacant city-owned site of 50 acres and,
if useful, part of an adjacent area with 200 condemned houses of 0.25 acre
each, which cost \$7000 each to buy and demolish. Open space, streets and
parking take 15 % of all land used. Formulate an LP that maximizes the net
income and solve it.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $x_L,x_M,x_U,x_P$ = number of lower-, middle-, upper-income
and public housing units; $x_S$ = number of classrooms; $x_R$ = number of
retail units; $y$ = number of condemned houses bought.

**Model.** We read the retail demand as a number of retail units per
housing unit (or classroom), and, since the statement asks for a single net
income, combine the annual income with the one-time demolition cost in one
objective (dollars).

$$
\begin{aligned}
\max\ & z=7000x_L+12\,000x_M+20\,000x_U+5000x_P+15\,000x_R \\
& \qquad-10\,000x_S-7000y \\
\st\ & 100\le x_L\le200 \\
& 125\le x_M\le190 \\
& 75\le x_U\le260 \\
& 300\le x_P\le600 \\
& 0\le x_R\le25 \\
& 0\le y\le200 \\
& x_S\ge0 \\
& 0.045x_S\le2 && \text{(school area)} \\
& 1.3x_L+1.2x_M+0.5x_U+1.4x_P\le25x_S && \text{(pupils)} \\
& 0.023x_L+0.034x_M+0.046x_U+0.023x_P+0.034x_S\le x_R && \text{(retail)} \\
& 0.05x_L+0.07x_M+0.03x_U+0.025x_P+0.045x_S+0.1x_R \\
& \qquad\le0.85\,(50+0.25y) && \text{(land)}
\end{aligned}
$$

**Solution.** Lower-income, middle-income and public housing stay at their
minimum numbers, $x_L^*=100$, $x_M^*=125$ and $x_P^*=300$; the retail units
are at their maximum, $x_R^*=25$; $y^*=0$; and

$$
\begin{aligned}
x_U^*&=\tfrac{264950}{1167}\approx227.04,\\
x_S^*&=\tfrac{37975}{1167}\approx32.54,
\end{aligned}
$$

with $z^*=\tfrac{3224925000}{389}\approx\$8\,290\,295.63$. The pupil and
retail constraints are binding; the plan uses 32.03 of the 42.5 usable
acres of the vacant site, so no condemned houses are bought, and the
school occupies 1.46 of its 2 acres.

**Interpretation.** Upper-income units earn the most per unit, so the plan
builds as many as the retail capacity allows: every upper-income unit
creates retail demand (0.046 units) and pupils who need classroom space,
which in turn creates further retail demand. The other housing types stay
at their minimum. In whole units, the best plan is
$(x_L,x_M,x_U,x_P,x_S,x_R,y)=(101,125,226,300,33,25,0)$ with \$8 272 000.
If the retail demand were measured in acres and compared with the
$0.1x_R\le2.5$ acres of retail space, no plan would be feasible: even the
minimum numbers of units generate 17.9 units of retail demand.

*(An earlier version of this model bounds $x_L$, instead of $x_P$, by 300
and 600, gives the school area as $0.45x_S$ instead of $0.045x_S$, and
credits the upper-income units with the coefficient of $x_M$. With these
slips the model is infeasible; the model above corrects them.)*
:::

(sec-08-01-blending)=
## Blending and percentage specifications

:::{prf:example} Cheapest inputs for castings
:label: ex-08-01-foundry

A foundry has a firm order for 1000 kg of castings that must contain at
least 0.45 % manganese and between 3.5 % and 6 % silicon; the castings sell
for \$0.85 per kg. Three pig irons are available:

| | Pig iron A | Pig iron B | Pig iron C |
|---|---|---|---|
| Manganese | 0.45 % | 0.55 % | 0.42 % |
| Silicon | 4 % | 1.3 % | 0.5 % |
| Cost (\$ per 1000 kg) | 40 | 45 | 25 |

Pure manganese can be added at \$14 per kg, and melting pig iron costs
0.9 cents per kg. Which inputs maximize the profit on the order?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $x_1,x_2,x_3$ = pig iron A, B, C used (thousands of kg);
$x_4$ = pure manganese added (kg).

**Costs.** Melting costs \$9 per 1000 kg, so a thousand kilograms of pig
iron A, B, C cost \$49, \$54 and \$34. The revenue is fixed at
$0.85\cdot1000=\$850$, so maximizing profit $850-\text{cost}$ is the same
as minimizing cost {eq}`eq-03-01-max-min`.

**Specifications.** With a fixed total of 1000 kg, the percentages are
amounts: at least 4.5 kg of manganese and between 35 and 60 kg of silicon.
One thousand kilograms of pig iron A contain 4.5 kg of manganese and 40 kg
of silicon, and so on.

**Model.**

$$
\begin{aligned}
\min\ & 49x_1+54x_2+34x_3+14x_4 && \text{(dollars)} \\
\st\ & 1000x_1+1000x_2+1000x_3+x_4=1000 && \text{(weight, kg)} \\
& 4.5x_1+5.5x_2+4.2x_3+x_4\ge4.5 && \text{(manganese, kg)} \\
& 40x_1+13x_2+5x_3\ge35 && \text{(silicon, kg)} \\
& 40x_1+13x_2+5x_3\le60 && \text{(silicon, kg)} \\
& x_1,\dots,x_4\ge0 .
\end{aligned}
$$

**Solution.** The manganese and lower silicon constraints are binding and
no pure manganese is added. Solving $x_1+x_2+x_3=1$,
$4.5x_1+5.5x_2+4.2x_3=4.5$, $40x_1+13x_2+5x_3=35$ gives

$$
\begin{gathered}
x_1^*=\tfrac{366}{431}\approx0.8492,\\
x_2^*=\tfrac{15}{431}\approx0.0348,\\
x_3^*=\tfrac{50}{431}\approx0.1160,\\
x_4^*=0,
\end{gathered}
$$

that is, about 849.2 kg of A, 34.8 kg of B and 116.0 kg of C, at a cost of
$\tfrac{20444}{431}\approx\$47.43$. The maximum profit is
$850-47.43=\$802.57$.

**Interpretation.** Cheap pig iron C is used as far as the silicon and
manganese minimums allow; A supplies most of the silicon. Pure manganese
is not worth buying: meeting one more kilogram of manganese requirement
through the pig-iron mix costs about \$13.46 (the dual value of the
manganese row), less than the \$14 price of pure manganese.

*(An earlier version of this model prints the silicon coefficient of pig
iron A as 4.0 instead of 40 in both silicon constraints. With 4.0 the
model is infeasible, because no mix would then contain 35 kg of silicon.)*
:::

:::{prf:example} Blending distillates into fuel additives
:label: ex-08-01-distillate-blending

A company blends three distillates, A, B and C, into three additives. The
specifications and selling prices of the additives, and the daily supply
and cost of the distillates, are:

| Additive | A at most | C at least | Price (\$/liter) |
|---|---|---|---|
| Deluxe | 60 % | 20 % | 1.80 |
| Standard | 15 % | 60 % | 1.50 |
| Economy | — | 50 % | 1.20 |

| Distillate | Supply (liters/day) | Cost (\$/liter) |
|---|---|---|
| A | 16 000 | 0.20 |
| B | 20 000 | 0.15 |
| C | 10 000 | 0.11 |

All additive produced can be sold. Find the blending plan that maximizes
the daily profit.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $x_{ij}$ = liters of distillate $i\in\{A,B,C\}$ used per day
in additive $j\in\{1,2,3\}$ (deluxe, standard, economy).

**Model.** By {prf:ref}`prop-03-01-ratio`, "A at most 60 % of deluxe" is
$x_{A1}\le0.6(x_{A1}+x_{B1}+x_{C1})$, and so on. With the selling prices
$r=(r_1,r_2,r_3)=(1.80,1.50,1.20)$ and the costs
$k=(k_A,k_B,k_C)=(0.20,0.15,0.11)$ dollars per liter:

$$
\begin{aligned}
\max\ & z=\sum_{j}r_j\sum_i x_{ij}-\sum_i k_i\sum_j x_{ij} \\
\st\ & x_{A1}+x_{A2}+x_{A3}\le16\,000 && \text{(supply of A)} \\
& x_{B1}+x_{B2}+x_{B3}\le20\,000 && \text{(supply of B)} \\
& x_{C1}+x_{C2}+x_{C3}\le10\,000 && \text{(supply of C)} \\
& 0.4x_{A1}-0.6x_{B1}-0.6x_{C1}\le0 && \text{(A in deluxe)} \\
& 0.85x_{A2}-0.15x_{B2}-0.15x_{C2}\le0 && \text{(A in standard)} \\
& -0.2x_{A1}-0.2x_{B1}+0.8x_{C1}\ge0 && \text{(C in deluxe)} \\
& -0.6x_{A2}-0.6x_{B2}+0.4x_{C2}\ge0 && \text{(C in standard)} \\
& -0.5x_{A3}-0.5x_{B3}+0.5x_{C3}\ge0 && \text{(C in economy)} \\
& x_{ij}\ge0 .
\end{aligned}
$$

**Solution.** All distillate goes into the deluxe additive:

$$
\begin{aligned}
x_{A1}^*&=16\,000\ \text{liters/day},\\
x_{B1}^*&=20\,000\ \text{liters/day},\\
x_{C1}^*&=10\,000\ \text{liters/day},
\end{aligned}
$$

all other $x_{ij}^*=0$, and

$$
\begin{aligned}
z^*&=1.8\cdot46\,000\\
&\quad-(0.20\cdot16\,000+0.15\cdot20\,000+0.11\cdot10\,000) \\
&=82\,800-7300\\
&=\$75\,500\ \text{per day}.
\end{aligned}
$$

**Check.** The deluxe blend contains $16\,000/46\,000=34.8\,\%$ of A
($\le60\,\%$) and $10\,000/46\,000=21.7\,\%$ of C ($\ge20\,\%$).

**Interpretation.** Each liter of distillate $i$ in additive $j$ earns the
margin $r_j-k_i$ ({ref}`sec-03-01-blending`), which is highest in the deluxe
additive for every distillate. Because the whole supply, mixed together,
already meets the deluxe specifications, nothing forces any distillate into
the cheaper additives. With a limited demand for deluxe the other additives
would be made as well, and their specifications would become binding.
:::

:::{prf:example} The cheapest handyman bag
:label: ex-08-01-hardware-bag

A hardware store fills handyman bags with screws, bolts, nuts and washers.
Screws come in 100-lb boxes at \$120, bolts in 100-lb boxes at \$175, nuts
in 80-lb boxes at \$75 and washers in 30-lb boxes at \$25. A bag must weigh
at least 1 lb and contain, by weight, at least 10 % screws and at least
25 % bolts, and at most 15 % nuts and at most 10 % washers. To balance the
bag, the number of bolts may not exceed the number of nuts, nor the number
of washers. A bolt weighs as much as 10 nuts or 50 washers. Find the
cheapest bag.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $x_s,x_b,x_n,x_w$ = weight of screws, bolts, nuts and
washers in a bag (lb); $y$ = weight of the bag (lb).

**Costs per lb.** Screws $120/100=\$1.20$, bolts $175/100=\$1.75$, nuts
$75/80=\$0.9375$, washers $25/30\approx\$0.8333$.

**Counts.** If a bolt weighs $\beta$ lb, a nut weighs $\beta/10$ and a washer
$\beta/50$, so the bag holds $x_b/\beta$ bolts, $10x_n/\beta$ nuts and $50x_w/\beta$
washers. "Bolts $\le$ nuts" is $x_b\le10x_n$ and "bolts $\le$ washers" is
$x_b\le50x_w$; $\beta$ cancels.

**Model.**

$$
\begin{aligned}
\min\ & 1.2x_s+1.75x_b+0.9375x_n+\tfrac56x_w && \text{(dollars per bag)} \\
\st\ & x_s+x_b+x_n+x_w=y && \text{(weight of the bag)} \\
& y\ge1 \\
& x_s\ge0.10y && \text{(screws)} \\
& x_b\ge0.25y && \text{(bolts)} \\
& x_n\le0.15y && \text{(nuts)} \\
& x_w\le0.10y && \text{(washers)} \\
& x_b\le10x_n && \text{(no more bolts than nuts)} \\
& x_b\le50x_w && \text{(no more bolts than washers)} \\
& x_s,x_b,x_n,x_w,y\ge0 .
\end{aligned}
$$

**Solution.** $y^*=1$ lb with $x_s^*=0.50$, $x_b^*=0.25$, $x_n^*=0.15$,
$x_w^*=0.10$ lb, at a cost of

$$
\begin{aligned}
&1.2\cdot0.5+1.75\cdot0.25+0.9375\cdot0.15+\tfrac56\cdot0.1 \\
&\quad=0.6+0.4375+0.140625+0.08\overline{3}\\
&\quad=\tfrac{1211}{960}\approx\$1.2615 .
\end{aligned}
$$

**Interpretation.** Nuts and washers are the cheapest per pound and are
used up to their maximum shares; bolts are the most expensive and are held
at their minimum share; screws fill the rest. The count constraints are
not binding ($0.25\le1.5$ and $0.25\le5$), and no bag heavier than the
minimum of 1 lb pays.

*(An earlier version of this solution reports the cost of this bag as
\$1.12; the correct cost is \$1.2615.)*
:::

:::{prf:example} Refinery with shortage penalties and storage costs
:label: ex-08-01-refinery

A refinery distils two crude oils, A and B, into regular gasoline, premium
gasoline and jet fuel. Every barrel of crude yields all three products in
fixed fractions:

| | Regular | Premium | Jet fuel | Cost (\$/bbl) | Available (bbl/day) |
|---|---|---|---|---|---|
| Crude A | 0.20 | 0.10 | 0.25 | 30 | 2500 |
| Crude B | 0.25 | 0.30 | 0.10 | 40 | 3000 |
| Demand (bbl/day) | 500 | 700 | 400 | | |
| Revenue (\$/bbl) | 50 | 70 | 120 | | |
| Storage cost of surplus (\$/bbl) | 2 | 3 | 4 | | |
| Penalty for unmet demand (\$/bbl) | 10 | 15 | 20 | | |

The refinery earns the revenue on its own output up to the demand.
Output beyond the demand is not sold now and incurs a storage cost;
demand that the refinery cannot meet is bought from outside and costs a
penalty per barrel. Formulate an LP for the daily plan and solve it.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $x_A,x_B$ = crude A and B distilled (bbl/day); for each
product $k\in\{R,P,J\}$, $s_k$ = surplus output and $u_k$ = unmet demand
(bbl/day).

**Balances.** Output minus surplus plus shortage equals demand:

$$
\begin{aligned}
0.20x_A+0.25x_B-s_R+u_R&=500 && \text{(regular)} \\
0.10x_A+0.30x_B-s_P+u_P&=700 && \text{(premium)} \\
0.25x_A+0.10x_B-s_J+u_J&=400 && \text{(jet fuel)}
\end{aligned}
$$

Own output sold is $d_k-u_k$ for product $k$.

**Model.**

$$
\begin{aligned}
\max\ & z=50(500-u_R)+70(700-u_P)+120(400-u_J)-30x_A-40x_B \\
& \qquad-(2s_R+3s_P+4s_J)-(10u_R+15u_P+20u_J) \\
\st\ & \text{the three balances above} \\
& 0\le x_A\le2500 \\
& 0\le x_B\le3000 \\
& s_k,\ u_k\ge0,\quad k\in\{R,P,J\} .
\end{aligned}
$$

**Solution.** Regular gasoline and jet fuel exactly meet demand,
$0.20x_A+0.25x_B=500$ and $0.25x_A+0.10x_B=400$, so

$$
\begin{aligned}
x_A^*&=\tfrac{20000}{17}\approx1176.47\ \text{bbl/day},\\
x_B^*&=\tfrac{18000}{17}\approx1058.82\ \text{bbl/day}.
\end{aligned}
$$

Premium output is $\tfrac{7400}{17}\approx435.29$ bbl/day, short of demand by
$u_P^*=\tfrac{4500}{17}\approx264.71$ bbl/day; there is no surplus. The
daily profit is

$$
\begin{aligned}
z^*&=\underbrace{103\,470.59}_{\text{revenue}}-\underbrace{77\,647.06}_{\text{crude}}
-\underbrace{3970.59}_{\text{penalty}}\\
&=\tfrac{371500}{17}\approx\$21\,852.94 .
\end{aligned}
$$

**Interpretation.** Making more premium would require more crude, and
every extra barrel of crude also yields regular gasoline and jet fuel
that would have to be stored. One more barrel of crude B, for instance,
costs \$40 and yields 0.30 bbl of premium, worth at most
$0.30\cdot(70+15)=\$25.50$ in revenue and avoided penalty, plus
$0.25+0.10$ bbl of surplus that costs \$0.90 to store. So the refinery
accepts the premium shortage. The dual value of the premium balance,
\$85 per bbl, is exactly the revenue plus penalty of one more barrel of
premium demand met.

*(An earlier formulation of this problem assigns each barrel of crude to a
single product, with variables such as "crude A used for regular
gasoline". That ignores that every barrel yields all three products at
once; its optimum, a loss of \$15 500 per day, does not describe the
refinery. The joint-yield model above is the correct one.)*
:::

(sec-08-01-staffing)=
## Staffing and covering

:::{prf:example} Scheduling office runners
:label: ex-08-01-errand-staffing

A university office hires students to run errands between 8:00 and 17:00.
It needs at least 2 students from 8:00 to 10:00, 4 from 10:00 to 11:00, 3
from 11:00 to 13:00 and 2 from 13:00 to 17:00. A student works 3
consecutive hours, except that students starting at 15:00 work 2 hours and
those starting at 16:00 work 1 hour. Students can start at any full hour
from 8:00 to 16:00 except at 12:00. Find a schedule with the fewest
students.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $x_i$ = number of students starting at hour $i$, where
$i=1,\dots,9$ stands for 8:00, 9:00, …, 16:00; $x_5=0$ (no start at noon).

**Covering constraints.** A student who starts at hour $i$ is present in
hours $i$, $i+1$, $i+2$ (as far as 17:00). The students present in hour
$t$ are therefore those with $t-2\le i\le t$:

$$
\begin{aligned}
\min\ & z=x_1+x_2+x_3+x_4+x_6+x_7+x_8+x_9 \\
\st\ & x_1\ge2 && \text{(8–9)} \\
& x_1+x_2\ge2 && \text{(9–10)} \\
& x_1+x_2+x_3\ge4 && \text{(10–11)} \\
& x_2+x_3+x_4\ge3 && \text{(11–12)} \\
& x_3+x_4\ge3 && \text{(12–13)} \\
& x_4+x_6\ge2 && \text{(13–14)} \\
& x_6+x_7\ge2 && \text{(14–15)} \\
& x_6+x_7+x_8\ge2 && \text{(15–16)} \\
& x_7+x_8+x_9\ge2 && \text{(16–17)} \\
& x_i\ge0 .
\end{aligned}
$$

**Solution.** The minimum is $z^*=8$ students, for example the schedule
$x_1=2$, $x_3=2$, $x_4=1$, $x_6=1$, $x_7=1$, $x_9=1$ (all other $x_i=0$),
which puts $2,2,4,3,3,2,2,2,2$ students in the nine hours, meeting every
requirement. The LP optimum is integral here; other optimal schedules
exist.

**Why 8 is a lower bound.** The hours 10–11, 13–14 and 16–17 are covered by
disjoint groups of start times, $\{1,2,3\}$, $\{4,5,6\}$ and $\{7,8,9\}$.
Adding their constraints gives $x_1+\dots+x_9\ge4+2+2=8$ for every
feasible schedule.

*(An earlier version of this solution lists $x_1=x_3=2$ and
$x_6=x_7=x_9=1$ with $z=8$. Those values add up to 7 students and leave
11:00–13:00 with 2 students and 13:00–14:00 with 1; the student starting at
11:00, $x_4=1$, is missing.)*
:::

:::{prf:example} Salespersons with two consecutive days off
:label: ex-08-01-cyclic-staffing

A department store open seven days a week needs at least 12, 18, 20, 28,
32, 40 and 40 salespersons on Monday to Sunday. Each salesperson works
five consecutive days and then has two consecutive days off, and the
off-days are staggered over the week. How many salespersons should be
hired, and how should their days off be arranged?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $x_i$ = number of salespersons whose five-day stretch starts
on day $i$ ($i=1$ Monday, …, $i=7$ Sunday); they are off on the two days
before day $i$.

**Covering constraints.** On Monday the people who started on Thursday,
Friday, Saturday, Sunday or Monday are at work, and similarly for the other
days:

$$
\begin{aligned}
\min\ & z=x_1+x_2+\dots+x_7 \\
\st\ & x_1+x_4+x_5+x_6+x_7\ge12 && \text{(Mon)} \\
& x_1+x_2+x_5+x_6+x_7\ge18 && \text{(Tue)} \\
& x_1+x_2+x_3+x_6+x_7\ge20 && \text{(Wed)} \\
& x_1+x_2+x_3+x_4+x_7\ge28 && \text{(Thu)} \\
& x_1+x_2+x_3+x_4+x_5\ge32 && \text{(Fri)} \\
& x_2+x_3+x_4+x_5+x_6\ge40 && \text{(Sat)} \\
& x_3+x_4+x_5+x_6+x_7\ge40 && \text{(Sun)} \\
& x_i\ge0 .
\end{aligned}
$$

**Solution.** The minimum is 42 salespersons, for example
$x^*=(0,2,24,0,6,8,2)$: 2 start on Tuesday (off Sunday–Monday), 24 on
Wednesday (off Monday–Tuesday), 6 on Friday (off Wednesday–Thursday), 8 on
Saturday (off Thursday–Friday) and 2 on Sunday (off Friday–Saturday). The
staff present on Monday to Sunday are 16, 18, 36, 28, 32, 40 and 40. The
LP optimum is integral; other optimal schedules exist.

**Why 42 is a lower bound.** Add the constraints for Tuesday, Thursday,
Saturday and Sunday. Every $x_i$ appears in at most three of them (and
$x_1$ in two), so $3(x_1+\dots+x_7)\ge18+28+40+40=126$ and
$z\ge42$. (These weights, $\tfrac13$ on each of the four days, are the
optimal dual solution; the argument is weak duality,
{prf:ref}`thm-03-05-weak-duality`.)
:::

(sec-08-01-multi-period)=
## Multi-period planning and investment

:::{prf:example} Capital goods and consumer goods over five periods
:label: ex-08-01-capital-consumer

A firm makes capital goods and consumer goods and plans five periods. At
the start of period 1 it holds 100 units of capital goods and 50 units of
consumer goods. One unit of capital goods, committed for two periods,
either produces 3 new units of capital goods, consuming 3 units of consumer
goods from stock at the start, or produces 10 units of consumer goods; in
both cases the unit itself is available again after the two periods
({numref}`fig-08-01-capital-activities`). Capital goods can also be left idle
for a period or sold; consumer goods can be sold or kept for the next
period. The prices per unit are:

| Period | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| Capital goods | 4 | 6 | 7 | 6 | 4 | 3 |
| Consumer goods | 0.10 | 0.15 | 0.15 | 0.15 | 0.10 | 0.10 |

Whatever the firm holds at the start of period 6, including goods that
activities started in period 4 deliver then, is valued at the period-6
prices. Find the plan that maximizes the total value.
:::

```{figure} figures/08-01-capital-activities.svg
:name: fig-08-01-capital-activities
:width: 620px
:class: dark-light
:alt: Two process boxes. Left: 1 unit of capital and 3 units of consumer goods in period t become 4 units of capital in period t plus 2. Right: 1 unit of capital in period t becomes 1 unit of capital and 10 units of consumer goods in period t plus 2

The two production activities of {prf:ref}`ex-08-01-capital-consumer`. Both
take two periods and return the unit of capital that runs them.
```

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** For each period $t=1,\dots,5$, in units of goods:

- $x_{t1}$ = capital committed to making capital (returns 4 units of
  capital in period $t+2$);
- $x_{t2}$ = capital committed to making consumer goods (returns 1 unit of
  capital and 10 units of consumer goods in period $t+2$);
- $x_{t3}$ = capital left idle (available in period $t+1$);
- $x_{t4}$ = capital sold; $x_{t5}$ = consumer goods sold;
- $x_{t6}$ = consumer goods carried to period $t+1$.

**Balances** (start of period $t$; terms with a period index below 1 are
omitted):

$$
\begin{aligned}
x_{t1}+x_{t2}+x_{t3}+x_{t4}&=[100]_{t=1}+x_{t-1,3}+4x_{t-2,1}+x_{t-2,2}
&& \text{(capital)} \\
3x_{t1}+x_{t5}+x_{t6}&=[50]_{t=1}+x_{t-1,6}+10x_{t-2,2} && \text{(consumer goods)}
\end{aligned}
$$

where $[100]_{t=1}$ is 100 in period 1 and 0 otherwise.

**Objective.** Sales in periods 1–5 plus the terminal value of what is on
hand at the start of period 6: idle capital $x_{53}$, carried consumer
goods $x_{56}$ and the output of the activities started in period 4
(activities started in period 5 deliver after the horizon and earn
nothing):

$$
\begin{aligned}
\max\ z={}&\sum_{t=1}^5\bigl(p^{\mathrm{cap}}_t x_{t4}+p^{\mathrm{con}}_t x_{t5}\bigr)
+3\bigl(x_{53}+4x_{41}+x_{42}\bigr)+0.1\bigl(x_{56}+10x_{42}\bigr),
\end{aligned}
$$

with the prices $p^{\mathrm{cap}}_t$, $p^{\mathrm{con}}_t$ of the table. All
variables are nonnegative.

**Solution.** $z^*=2457.5$ monetary units, with

$$
\begin{aligned}
x_{11}^*&=\tfrac{50}{3}\approx16.67,\\
x_{12}^*&=\tfrac{250}{3}\approx83.33,\\
x_{31}^*&=150,\\
x_{35}^*&=\tfrac{1150}{3}\approx383.33,\\
x_{54}^*&=600,
\end{aligned}
$$

and all other variables zero. Read period by period:

- Period 1: the 50 consumer goods allow $50/3$ units of capital to make
  capital; the other $250/3$ units make consumer goods.
- Period 3: $4\cdot\tfrac{50}{3}+\tfrac{250}{3}=150$ units of capital and
  $10\cdot\tfrac{250}{3}=\tfrac{2500}{3}$ consumer goods return. All 150
  units make capital, using 450 consumer goods; the remaining
  $\tfrac{1150}{3}$ are sold at 0.15.
- Period 5: $4\cdot150=600$ units of capital return and are sold at 4.

The value is $0.15\cdot\tfrac{1150}{3}+4\cdot600=57.5+2400=2457.5$.

**Interpretation.** Making capital multiplies it by four every two
periods, but it needs consumer goods, so the first period produces exactly
the consumer goods that the period-3 expansion requires. Activities
started in period 5 would deliver after the horizon, so the capital is
sold in period 5 at 4 rather than held to the start of period 6, where it
would be valued at only 3. The optimum is not unique: the consumer goods
sold in period 3 could equally be carried to period 4 and sold there at the
same price.

*(An earlier version of this example speaks of a five-period horizon
but gives prices for six periods; the model above plans periods 1–5 and
values the holdings at the start of period 6. An earlier version of the
model also prices consumer goods sold in period 4 at 0.1 instead of 0.15
and omits the consumer goods used by $x_{41}$ from one balance; neither
slip changes the optimum, since $x_{41}^*=x_{45}^*=0$.)*
:::

:::{prf:example} Production and inventory of two products
:label: ex-08-01-production-inventory

A company must deliver two products, A and B, in June, July and August:

| | June | July | August |
|---|---|---|---|
| Demand for A (units) | 500 | 5000 | 750 |
| Demand for B (units) | 1000 | 1200 | 1200 |
| Production capacity (h) | 3000 | 3500 | 3000 |

A unit of A takes 0.75 h of capacity and a unit of B takes 1 h. Demand
must be met on time, but goods made in one month may be held for a later
month at a holding cost of \$0.90 (A) and \$0.75 (B) per unit per month.
Production costs are \$30 per unit of A and \$28 per unit of B. Find the
cheapest production plan.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** For product $k\in\{A,B\}$ and month $t=1,2,3$ (June to
August): $p_{kt}$ = units produced, $I_{kt}$ = units in stock at the end of
the month ($I_{k0}=0$).

**Model.** With $d_{kt}$ the demand of the table and $H_t$ the capacity of
month $t$, $H=(3000,3500,3000)$ h, the balances
{eq}`eq-03-01-inventory-balance` with $I_{kt}\ge0$ ensure that demand is met
on time.

$$
\begin{aligned}
\min\ & z=\sum_{t=1}^3\bigl(30p_{At}+28p_{Bt}+0.90I_{At}+0.75I_{Bt}\bigr) \\
\st\ & 0.75p_{At}+p_{Bt}\le H_t,\quad t=1,2,3 && \text{(capacity, h)} \\
& I_{At}=I_{A,t-1}+p_{At}-d_{At},\quad t=1,2,3 && \text{(stock of A)} \\
& I_{Bt}=I_{B,t-1}+p_{Bt}-d_{Bt},\quad t=1,2,3 && \text{(stock of B)} \\
& p_{kt},\ I_{kt}\ge0 .
\end{aligned}
$$

**Solution.**

| | June | July | August |
|---|---|---|---|
| Production of A | 833.33 | 4666.67 | 750 |
| Production of B | 2200 | 0 | 1200 |
| Stock of A at month end | 333.33 | 0 | 0 |
| Stock of B at month end | 1200 | 0 | 0 |
| Hours used | 2825 | 3500 | 1762.5 |

The minimum cost is $z^*=\$283\,900$: production $30\cdot6250+28\cdot3400=
\$282\,700$ plus holding $0.90\cdot333.33+0.75\cdot1200=\$1200$.

**Interpretation.** July's demand needs $0.75\cdot5000+1200=4950$ h, 1450 h
more than July's capacity, so part of it must be made in June. Freeing one
hour of July capacity costs \$0.75 of holding if one unit of B is made in
June instead, but $\tfrac43\cdot0.90=\$1.20$ if $\tfrac43$ units of A are.
So all of July's B (1200 h) is made in June, and only the remaining 250 h,
that is 333.33 units of A, are made early. June keeps 175 h of slack.

*(Read as production rates of 0.75 and 1 units per hour, the data would
make June and July together need
$(500+5000)/0.75+2200\approx9533$ h against 6500 h available, and the
problem would be infeasible. We read the rate for A as 0.75 hours per
unit.)*
:::

:::{prf:example} Investing in four projects and a bank account
:label: ex-08-01-project-investment

An investor has \$10 000. Four projects are open; per dollar invested they
produce the following cash flows (\$) at the start of years 1 to 5 (a
negative entry is the investment itself):

| Project | Year 1 | Year 2 | Year 3 | Year 4 | Year 5 |
|---|---|---|---|---|---|
| 1 | −1.00 | 0.50 | 0.30 | 1.80 | 1.20 |
| 2 | −1.00 | 0.60 | 0.20 | 1.50 | 1.30 |
| 3 | 0.00 | −1.00 | 0.80 | 1.90 | 0.80 |
| 4 | −1.00 | 0.40 | 0.60 | 1.80 | 0.95 |

Money can also be kept in a bank account that pays 6.5 % per year, and all
money returned can be reinvested. Find the investment plan that maximizes
the cash available at the start of year 5.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $x_i$ = dollars invested in project $i$ (projects 1, 2, 4 at
the start of year 1, project 3 at the start of year 2); $y_t$ = dollars
deposited in the bank at the start of year $t$, $t=1,\dots,4$; $y_5$ = cash
at the start of year 5.

**Cash balances.** At the start of each year, the money available (returns
plus bank balance with interest) is either invested or deposited:

$$
\begin{aligned}
\max\ & z=y_5 \\
\st\ & x_1+x_2+x_4+y_1\le10\,000 && \text{(year 1)} \\
& x_3+y_2=0.5x_1+0.6x_2+0.4x_4+1.065y_1 && \text{(year 2)} \\
& y_3=0.3x_1+0.2x_2+0.8x_3+0.6x_4+1.065y_2 && \text{(year 3)} \\
& y_4=1.8x_1+1.5x_2+1.9x_3+1.8x_4+1.065y_3 && \text{(year 4)} \\
& y_5=1.2x_1+1.3x_2+0.8x_3+0.95x_4+1.065y_4 && \text{(year 5)} \\
& x_i,\ y_t\ge0 .
\end{aligned}
$$

**Solution.** $x_2^*=\$10\,000$, $x_3^*=\$6000$, $y_3^*=\$6800$,
$y_4^*=\$33\,642$ and

$$
\begin{aligned}
z^*=y_5^*&=1.3\cdot10\,000+0.8\cdot6000+1.065\cdot33\,642\\
&=\$53\,628.73 .
\end{aligned}
$$

**Check.** Year 2: project 2 returns \$6000, all invested in project 3.
Year 3: $0.2\cdot10\,000+0.8\cdot6000=6800$ goes to the bank. Year 4:
$1.5\cdot10\,000+1.9\cdot6000+1.065\cdot6800=33\,642$.

**Interpretation.** Project 2's early return (0.60 per dollar in year 2)
finances project 3, which returns 3.50 per dollar over the following
three years. Project 1 has a higher total return per dollar (3.80 against
3.60), but less of it arrives in year 2, when project 3 can absorb it. The
bank is used only for money that no project can take.
:::

:::{prf:example} One-year and two-year investment plans
:label: ex-08-01-two-plans

An executive can invest in two plans. Plan A returns \$1.70 one year later
for every dollar invested and accepts money every year. Plan B returns
\$3 two years later for every dollar and accepts money only for whole
two-year periods. How should \$100 000 be invested to maximize the value at
the end of year 3?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables** (in \$1000). $x_{tA}$ = amount put in plan A at the start of
year $t=1,2,3$; $x_{tB}$ = amount put in plan B at the start of year
$t=1,2$ (money put in plan B in year 3 would not return by the end of
year 3).

**Model.**

$$
\begin{aligned}
\max\ & z=3x_{2B}+1.7x_{3A} && \text{(value at the end of year 3)} \\
\st\ & x_{1A}+x_{1B}\le100 && \text{(year 1)} \\
& x_{2A}+x_{2B}=1.7x_{1A} && \text{(year 2)} \\
& x_{3A}=3x_{1B}+1.7x_{2A} && \text{(year 3)} \\
& x_{tA},\ x_{tB}\ge0 .
\end{aligned}
$$

**Solution.** $z^*=510$, that is, \$510 000 at the end of year 3 (earnings of
\$410 000). There are two optimal plans:

- A in year 1, then B: $x_{1A}=100$, $x_{2B}=170$, value $3\cdot170=510$;
- B in year 1, then A: $x_{1B}=100$, $x_{3A}=300$, value $1.7\cdot300=510$.

Both give a factor $1.7\cdot3=5.1$ over three years, and any split between
them is optimal too. Plan A three times gives only $1.7^3=4.913$.

*(An earlier version of this solution reports \$500 000 as the optimal
value; the plan it lists, $x_{1A}=100$ and $x_{2B}=170$, is worth
$3\cdot170\,000=\$510\,000$.)*
:::

:::{prf:example} Scheduling four renewal projects
:label: ex-08-01-renewal-projects

A city will carry out four housing projects over the next five years:

| Project | Years | Cost (million \$) | Annual income when complete (million \$) |
|---|---|---|---|
| 1 | 1–3 | 5.0 | 0.05 |
| 2 | 2–5 | 8.0 | 0.07 |
| 3 | 1–5 | 15.0 | 0.15 |
| 4 | 3–4 | 1.2 | 0.02 |

The budgets for years 1 to 5 are 3, 6, 7, 7 and 7 million dollars. Projects
1 and 4 must be completed within their years; projects 2 and 3 may be
completed only partly, but at least 25 % of each must be done within its
years. A part completed in a year is occupied at the end of that year and
from then on earns the same fraction of the annual income; for example,
40 % of project 1 completed in year 1 and 60 % in year 3 earn
$(4\cdot0.4+2\cdot0.6)\cdot\$50\,000$ over the five years. Ignoring the time
value of money, find the schedule that maximizes the total income over the
five years.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $x_{ij}$ = fraction of project $i$ completed in year $j$,
for the years $j$ in which project $i$ runs.

**Income.** A fraction completed in year $j$ earns income in years
$j+1,\dots,5$, that is, for $5-j$ years.

**Model** (in million dollars).

$$
\begin{aligned}
\max\ & z=0.05(4x_{11}+3x_{12}+2x_{13})+0.07(3x_{22}+2x_{23}+x_{24}) \\
& \qquad+0.15(4x_{31}+3x_{32}+2x_{33}+x_{34})+0.02(2x_{43}+x_{44}) \\
\st\ & x_{11}+x_{12}+x_{13}=1 && \text{(project 1)} \\
& 0.25\le x_{22}+x_{23}+x_{24}+x_{25}\le1 && \text{(project 2)} \\
& 0.25\le x_{31}+\dots+x_{35}\le1 && \text{(project 3)} \\
& x_{43}+x_{44}=1 && \text{(project 4)} \\
& 5x_{11}+15x_{31}\le3 && \text{(year 1)} \\
& 5x_{12}+8x_{22}+15x_{32}\le6 && \text{(year 2)} \\
& 5x_{13}+8x_{23}+15x_{33}+1.2x_{43}\le7 && \text{(year 3)} \\
& 8x_{24}+15x_{34}+1.2x_{44}\le7 && \text{(year 4)} \\
& 8x_{25}+15x_{35}\le7 && \text{(year 5)} \\
& x_{ij}\ge0 .
\end{aligned}
$$

**Solution.** $z^*=0.52375$, that is, \$523 750. One optimal schedule
completes project 1 with $x_{11}=0.6$ and $x_{12}=0.4$, project 2 with
$x_{24}=0.225$ and $x_{25}=0.025$, project 3 with $x_{32}=\tfrac{4}{15}$,
$x_{33}=\tfrac{29}{75}$ and $x_{34}=\tfrac{26}{75}$, and project 4 with
$x_{43}=1$; all other $x_{ij}=0$.

**Check.** Year budgets: $5\cdot0.6=3$; $5\cdot0.4+15\cdot\tfrac4{15}=6$;
$15\cdot\tfrac{29}{75}+1.2=7$; $8\cdot0.225+15\cdot\tfrac{26}{75}=7$;
$8\cdot0.025=0.2\le7$. Project 2 is 25 % done, project 3 100 %. Income:
$0.05\cdot3.6+0.07\cdot0.225+0.15\cdot1.92+0.02\cdot2=0.18+0.01575+0.288+0.04
=0.52375$.

**Interpretation.** Per million dollars spent, project 4 earns the most
income per year, then projects 1 and 3, then project 2; so project 2 is held
at its 25 % minimum and pushed late, and the budgets of years 1 to 4 are
used completely. The optimum is not unique: for example, year 1 can fund
project 3 instead of project 1, since both earn \$0.04 million of income
per million dollars spent in year 1.

*(An earlier version of this problem speaks of five projects but describes
four, and an earlier version of this solution gives $x_{24}=0.255$ and $x_{25}=0.045$;
with those values year 4 spends 7.23 million dollars, over its budget. The
values above give the same optimal income and are feasible.)*
:::

(sec-08-01-objectives)=
## Maximin and absolute-deviation objectives

:::{prf:example} Betting to maximize the worst-case return
:label: ex-08-01-max-min-bets

A game lets a player divide a stake among four choices; then one of three
outcomes occurs. The net return per dollar placed on each choice is:

| Outcome | Choice 1 | Choice 2 | Choice 3 | Choice 4 |
|---|---|---|---|---|
| 1 | −3 | 4 | −7 | 15 |
| 2 | 5 | −3 | 9 | 4 |
| 3 | 3 | −9 | 10 | −8 |

The player has \$1500 to play once and does not know which outcome will
occur, so she wants to maximize the smallest of the three possible
returns. Write this as a nonlinear model, turn it into an LP and solve it.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $x_j$ = dollars placed on choice $j$, $j=1,\dots,4$.

**Nonlinear model.** With $\ell_k(\vx)$ the return of outcome $k$,

$$
\begin{aligned}
\max_{\vx}\ & \min\{\ell_1(\vx),\ell_2(\vx),\ell_3(\vx)\}\\
\st\ & x_1+x_2+x_3+x_4\le1500\\
& \vx\ge\vzero,
\end{aligned}
$$

where $\ell_1=-3x_1+4x_2-7x_3+15x_4$, $\ell_2=5x_1-3x_2+9x_3+4x_4$ and
$\ell_3=3x_1-9x_2+10x_3-8x_4$. The objective is the minimum of linear
functions, which is not linear.

**LP.** Introduce a free variable $\eta$ for the guaranteed return
({ref}`sec-03-01-objectives`):

$$
\begin{aligned}
\max\ & \eta \\
\st\ & \eta\le-3x_1+4x_2-7x_3+15x_4 \\
& \eta\le5x_1-3x_2+9x_3+4x_4 \\
& \eta\le3x_1-9x_2+10x_3-8x_4 \\
& x_1+x_2+x_3+x_4\le1500 \\
& x_j\ge0,\quad j=1,\dots,4 \\
& \eta\ \text{free}.
\end{aligned}
$$

$\eta$ must be free, because the guaranteed return could be negative.

**Solution.** $x_3^*=862.50$, $x_4^*=637.50$, $x_1^*=x_2^*=0$, with returns
$\ell_1=\ell_3=\$3525$ and $\ell_2=\$10\,312.50$; the guaranteed return is
$\eta^*=\$3525$. At the optimum the two worst outcomes are balanced:
$-7x_3+15x_4=10x_3-8x_4$ with $x_3+x_4=1500$ gives $x_3=\tfrac{1725}{2}$,
$x_4=\tfrac{1275}{2}$.

**Why no plan does better.** Weight outcome 1 by 0.45 and outcome 3 by
0.55. Per dollar, the weighted return of the four choices is 0.30, −3.15,
2.35 and 2.35, never more than 2.35. Hence for every plan
$\min_k\ell_k\le0.45\ell_1+0.55\ell_3\le2.35\cdot1500=3525$. (These weights are the
dual solution of the LP, and the argument is weak duality,
{prf:ref}`thm-03-05-weak-duality`.)
:::

:::{prf:example} Fitting a line by least absolute deviations
:label: ex-08-01-lad-regression

In a ten-week typing course, the average speed of the students (words per
minute) after $x$ weeks was:

| Week $x$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| Speed $y$ (wpm) | 5 | 9 | 15 | 19 | 21 | 24 | 26 | 30 | 31 | 35 |

Find the line $\hat y=ax+b$ that minimizes the sum of the absolute
deviations $\sum_i\abs{y_i-(ax_i+b)}$. Formulate it as an LP in two ways,
explain why each is equivalent to the original problem, and solve it.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $a$ (wpm per week) and $b$ (wpm), both free; the deviation
of point $i$ is $e_i=y_i-ax_i-b$.

**First LP: bound each absolute value.**

$$
\begin{aligned}
\min\ & \sum_{i=1}^{10}u_i\\
\st\ & u_i\ge y_i-ax_i-b,\quad i=1,\dots,10,\\
& u_i\ge-(y_i-ax_i-b),\quad i=1,\dots,10,\\
& a,b\ \text{free}.
\end{aligned}
$$

The two constraints say $u_i\ge\abs{e_i}$, so $u_i\ge0$ automatically.
Since the objective pushes every $u_i$ down and $u_i$ appears nowhere else,
at an optimum $u_i=\max\{e_i,-e_i\}=\abs{e_i}$: otherwise lowering $u_i$
would stay feasible and reduce the objective.

**Second LP: split each deviation.**

$$
\begin{aligned}
\min\ & \sum_{i=1}^{10}(e_i^++e_i^-)\\
\st\ & e_i^+-e_i^-=y_i-ax_i-b,\quad i=1,\dots,10,\\
& e_i^+,e_i^-\ge0,\quad i=1,\dots,10,\\
& a,b\ \text{free}.
\end{aligned}
$$

Any $e_i$ can be written as $e_i^+-e_i^-$, and $e_i^++e_i^-\ge\abs{e_i}$
with equality exactly when one of the two parts is zero. If both were
positive at an optimum, reducing both by $\min\{e_i^+,e_i^-\}$ would keep
the constraint and lower the objective; so at an optimum
$e_i^++e_i^-=\abs{e_i}$. (At a basic feasible solution (BFS) this is automatic:
the columns of $e_i^+$ and $e_i^-$ are negatives of each other, so by
{prf:ref}`thm-03-02-extreme-bfs` they cannot both be positive.)

**Solution.**

$$
\begin{gathered}
a^*=\tfrac{20}{7}\approx2.857,\\
b^*=\tfrac{45}{7}\approx6.429,\\
\sum_i\abs{e_i}=\tfrac{81}{7}\approx11.57\ \text{wpm}.
\end{gathered}
$$

The line passes exactly through the data points of weeks 3 and 10
($3a+b=15$, $10a+b=35$); the deviations of the ten points are
$-\tfrac{30}{7},-\tfrac{22}{7},0,\tfrac87,\tfrac27,\tfrac37,-\tfrac37,
\tfrac57,-\tfrac87,0$, and their absolute values add up to $\tfrac{81}{7}$.
{numref}`fig-08-01-lad-fit` compares this line with the least-squares line
$\hat y=3.170x+4.067$, whose sum of absolute deviations is 12.51.

**Interpretation.** A least-absolute-deviation line passes through (at
least) two data points, the analog of a basic solution, and it is less
affected by the slow start in weeks 1 and 2 than the least-squares line,
which squares those large deviations.
:::

```{figure} figures/08-01-lad-fit.svg
:name: fig-08-01-lad-fit
:width: 500px
:class: dark-light
:alt: Ten data points of typing speed against week with two fitted lines; the least-absolute-deviation line passes through the points of weeks 3 and 10, the dashed least-squares line is slightly steeper

Typing speed against week, with the least-absolute-deviation line
$\hat y=\tfrac{20}{7}x+\tfrac{45}{7}$ of {prf:ref}`ex-08-01-lad-regression`
(solid; it passes through two data points, marked) and the least-squares
line (dashed).
```

(sec-08-01-graphical)=
## Graphical solution

:::{prf:example} Minimum capacity of a refinery
:label: ex-08-01-oilco

A company is building a refinery that must supply at least 14 000 bbl/day
of diesel, 30 000 bbl/day of gasoline, 10 000 bbl/day of lubricants and
8000 bbl/day of jet fuel. Its crude comes from two sources, Iraq and Dubai,
and its supply contracts require that at least 40 % of the crude come from
Iraq. A barrel of Iraqi crude yields 0.20 bbl of diesel, 0.25 bbl of
gasoline, 0.10 bbl of lubricant and 0.15 bbl of jet fuel; a barrel of Dubai
crude yields 0.10, 0.60, 0.15 and 0.10 bbl. Find the smallest refinery
capacity (bbl/day of crude) that meets the demands, using the graphical
method.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Variables.** $x_1,x_2$ = crude from Iraq and from Dubai (thousand
bbl/day).

**Model.** "At least 40 % from Iraq" is $x_1\ge0.4(x_1+x_2)$, that is
$-0.6x_1+0.4x_2\le0$ or $x_2\le1.5x_1$ ({prf:ref}`prop-03-01-ratio`).

$$
\begin{aligned}
\min\ & z=x_1+x_2 && \text{(capacity, thousand bbl/day)} \\
\st\ & -0.6x_1+0.4x_2\le0 && \text{(Iraqi share)} \\
& 0.20x_1+0.10x_2\ge14 && \text{(diesel)} \\
& 0.25x_1+0.60x_2\ge30 && \text{(gasoline)} \\
& 0.10x_1+0.15x_2\ge10 && \text{(lubricants)} \\
& 0.15x_1+0.10x_2\ge8 && \text{(jet fuel)} \\
& x_1,x_2\ge0 .
\end{aligned}
$$

**Graph.** {numref}`fig-08-01-oilco-region` shows the region. It is unbounded
(more crude always satisfies the demands), with the corners

| Corner | Lines that meet | $z=x_1+x_2$ |
|---|---|---|
| $(40,60)$ | Iraqi share, diesel | 100 |
| $(55,30)$ | diesel, lubricants | 85 |
| $(\tfrac{200}{3},\tfrac{200}{9})\approx(66.67,22.22)$ | lubricants, gasoline | 88.89 |
| $(120,0)$ | gasoline, $x_2=0$ | 120 |

The jet-fuel constraint is redundant: at every point of the region the jet
fuel output is at least 11.25 thousand bbl/day. Since the objective is
minimized, the contour $x_1+x_2=z$ is moved down and to the left; the last
point of the region it touches is $(55,30)$.

**Solution.** $x_1^*=55$, $x_2^*=30$ thousand bbl/day, and the minimum
capacity is $z^*=85\,000$ bbl/day. Diesel ($11+3=14$) and lubricants
($5.5+4.5=10$) are binding; gasoline has a surplus of 1.75 and jet fuel of
3.25 thousand bbl/day. Iraqi crude makes up $55/85=64.7\,\%$, above the
required 40 %.

**Interpretation.** The region is unbounded, but the LP is not: an
unbounded region with a minimized objective that grows in every feasible
direction still has an optimal corner ({ref}`sec-03-02-graphical-cases`).
:::

```{figure} figures/08-01-oilco-region.svg
:name: fig-08-01-oilco-region
:width: 540px
:class: dark-light
:alt: Unbounded feasible region above and to the right of four constraint lines, cut on the left by the line x2 = 1.5 x1; the dashed contour x1 + x2 = 85 touches the region at the corner (55, 30)

Feasible region of {prf:ref}`ex-08-01-oilco` (shaded, unbounded towards the
upper right) in thousand bbl/day. The dashed contour $x_1+x_2=85$ touches
the region only at $(55,30)$, where the diesel and lubricant lines meet. The
jet-fuel line lies entirely outside the region.
```

(sec-08-01-basic-solutions)=
## Standard form and basic solutions

:::{prf:example} All basic solutions of a two-variable LP
:label: ex-08-01-basic-solutions-2d

Consider

$$
\begin{aligned}
\max\ & z=2x_1+3x_2\\
\st\ & x_1+3x_2\le12\\
& 3x_1+2x_2\le12\\
& x_1,x_2\ge0 .
\end{aligned}
$$

**(a)** Write the LP in standard (equation) form.
**(b)** Find all its basic solutions and classify them as feasible or
infeasible.
**(c)** Find the best basic feasible solution by evaluating $z$.
**(d)** Solve the LP graphically and confirm that the solution of (c) is
optimal.
**(e)** Where do the infeasible basic solutions lie in the graph?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** With slack variables $s_1,s_2$:

$$
\begin{aligned}
x_1+3x_2+s_1&=12,\\
3x_1+2x_2+s_2&=12,\\
x_1,x_2,s_1,s_2&\ge0 .
\end{aligned}
$$

**(b)** With $m=2$ equations and $n=4$ variables there are $\binom42=6$
candidate bases, all nonsingular:

| Basic variables | $(x_1,x_2,s_1,s_2)$ | Feasible? | $z$ |
|---|---|---|---|
| $x_1,x_2$ | $(\tfrac{12}{7},\tfrac{24}{7},0,0)$ | yes | $\tfrac{96}{7}$ |
| $x_1,s_1$ | $(4,0,8,0)$ | yes | 8 |
| $x_1,s_2$ | $(12,0,0,-24)$ | no | — |
| $x_2,s_1$ | $(0,6,-6,0)$ | no | — |
| $x_2,s_2$ | $(0,4,0,4)$ | yes | 12 |
| $s_1,s_2$ | $(0,0,12,12)$ | yes | 0 |

For example, the basis $\{x_1,x_2\}$ solves $x_1+3x_2=12$, $3x_1+2x_2=12$:
three times the first minus the second gives $7x_2=24$.

**(c)** The best of the four basic feasible solutions is
$(x_1,x_2)=(\tfrac{12}{7},\tfrac{24}{7})$ with $z=\tfrac{96}{7}\approx13.71$.

**(d)** {numref}`fig-08-01-basic-solutions-2d` shows the region with corners
$(0,0)$, $(4,0)$, $(\tfrac{12}{7},\tfrac{24}{7})$ and $(0,4)$. The contour
$2x_1+3x_2=\tfrac{96}{7}$ has slope $-\tfrac23$, between the slopes
$-\tfrac13$ and $-\tfrac32$ of the two constraint lines, so it touches the
region only at $(\tfrac{12}{7},\tfrac{24}{7})$: the graphical optimum is the
best BFS, as {prf:ref}`thm-03-02-fundamental` guarantees.

**(e)** The infeasible basic solutions $(12,0)$ and $(0,6)$ are the points
where a constraint line crosses an axis outside the region: $(12,0)$ lies
on $x_1+3x_2=12$ but violates $3x_1+2x_2\le12$ (hence $s_2=-24$), and
$(0,6)$ lies on $3x_1+2x_2=12$ but violates $x_1+3x_2\le12$ ($s_1=-6$).
:::

```{figure} figures/08-01-basic-solutions-2d.svg
:name: fig-08-01-basic-solutions-2d
:width: 500px
:class: dark-light
:alt: Quadrilateral feasible region with filled dots at its four corners and hollow dots at (12, 0) and (0, 6) outside it; a dashed contour touches the region at (12/7, 24/7)

The LP of {prf:ref}`ex-08-01-basic-solutions-2d`. Filled dots are the four
basic feasible solutions, hollow dots the two infeasible basic solutions;
the dashed contour $2x_1+3x_2=\tfrac{96}{7}$ touches the region at the
optimum $(\tfrac{12}{7},\tfrac{24}{7})$.
```

:::{prf:example} Solving LPs by enumerating basic solutions
:label: ex-08-01-enumerate-bases

Find the optimal solution of each LP by enumerating all its basic
solutions.

**(a)**
$\max\ z=2x_1-4x_2+5x_3-6x_4$ subject to $x_1+4x_2-2x_3+8x_4\le2$,
$-x_1+2x_2+3x_3+4x_4\le1$, $x_1,\dots,x_4\ge0$.

**(b)**
$\min\ z=x_1+2x_2-3x_3-2x_4$ subject to $x_1+2x_2-3x_3+x_4=4$,
$x_1+2x_2+x_3+2x_4=4$, $x_1,\dots,x_4\ge0$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** Add slacks $x_5,x_6$. There are $\binom62=15$ candidate bases. The
pair $\{x_2,x_4\}$ is singular (its columns $(4,2)$ and $(8,4)$ are
proportional). The other 14 give only 8 distinct basic solutions, because
several are degenerate and belong to more than one basis:

| Basic solution $(x_1,\dots,x_6)$ | Bases | Feasible? | $z$ |
|---|---|---|---|
| $(0,\tfrac12,0,0,0,0)$ | $\{1,2\},\{2,3\},\{2,5\},\{2,6\}$ | yes (degenerate) | $-2$ |
| $(8,0,3,0,0,0)$ | $\{1,3\}$ | yes | **31** |
| $(0,0,0,\tfrac14,0,0)$ | $\{1,4\},\{3,4\},\{4,5\},\{4,6\}$ | yes (degenerate) | $-\tfrac32$ |
| $(-1,0,0,0,3,0)$ | $\{1,5\}$ | no | — |
| $(2,0,0,0,0,3)$ | $\{1,6\}$ | yes | 4 |
| $(0,0,\tfrac13,0,\tfrac83,0)$ | $\{3,5\}$ | yes | $\tfrac53$ |
| $(0,0,-1,0,0,4)$ | $\{3,6\}$ | no | — |
| $(0,0,0,0,2,1)$ | $\{5,6\}$ | yes | 0 |

The optimum is $\vx^*=(8,0,3,0)$ with $z^*=31$. Check: $8-6=2$ and
$-8+9=1$, so both constraints are binding.

**(b)** There are $\binom42=6$ candidate bases.

| Basic variables | Basic solution $(x_1,x_2,x_3,x_4)$ | Feasible? | $z$ |
|---|---|---|---|
| $x_1,x_2$ | singular | — | — |
| $x_1,x_3$ | $(4,0,0,0)$ | yes (degenerate) | 4 |
| $x_1,x_4$ | $(4,0,0,0)$ | yes (degenerate) | 4 |
| $x_2,x_3$ | $(0,2,0,0)$ | yes (degenerate) | 4 |
| $x_2,x_4$ | $(0,2,0,0)$ | yes (degenerate) | 4 |
| $x_3,x_4$ | $(0,0,-\tfrac47,\tfrac{16}{7})$ | no | — |

Both basic feasible solutions have $z=4$, so the optimal value is
$z^*=4$, attained at $(4,0,0,0)$ and $(0,2,0,0)$. In fact every feasible
point is optimal: subtracting the first equation from the second gives
$4x_3+x_4=0$, which with $x_3,x_4\ge0$ forces $x_3=x_4=0$. The feasible set
is the segment $x_1+2x_2=4$, $x_1,x_2\ge0$, on which
$z=x_1+2x_2=4$ is constant.
:::

:::{prf:example} The feasible set of a small system
:label: ex-08-01-feasible-segment

Describe and draw, in two dimensions, the feasible set of the system of
{prf:ref}`ex-03-02-basic-solutions`,

$$
\begin{aligned}
2x_1+3x_2+x_3&=5,\\
x_1+2x_2-x_3&=3,\\
x_1,x_2,x_3&\ge0 .
\end{aligned}
$$
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Parametrize.** Adding the two equations eliminates $x_3$:
$3x_1+5x_2=8$. With $x_2=t$ as parameter,

$$
\begin{aligned}
x_1&=\frac{8-5t}{3},\\
x_2&=t,\\
x_3&=5-2x_1-3x_2=\frac{t-1}{3} .
\end{aligned}
$$

**Nonnegativity.** $x_1\ge0$ requires $t\le\tfrac85$ and $x_3\ge0$ requires
$t\ge1$. The feasible set is therefore the segment

$$
\Bigl\{\bigl(\tfrac{8-5t}{3},\ t,\ \tfrac{t-1}{3}\bigr):\ 1\le t\le\tfrac85\Bigr\}
$$

in $\R^3$. Its end points, $t=1$ and $t=\tfrac85$, are the two basic
feasible solutions $(1,1,0)$ and $(0,\tfrac85,\tfrac15)$ found in
{prf:ref}`ex-03-02-basic-solutions`, as {prf:ref}`thm-03-02-extreme-bfs` predicts:
a segment has exactly its two end points as extreme points. An interior
point such as $t=\tfrac{13}{10}$, $\vx=(\tfrac12,\tfrac{13}{10},\tfrac1{10})$,
has three positive components and is not basic.

**Drawing.** Projected on the $(x_1,x_2)$ plane, the set is the piece of
the line $3x_1+5x_2=8$ between $(1,1)$ and $(0,\tfrac85)$
({numref}`fig-08-01-feasible-segment`). The infeasible basic solution
$(\tfrac83,0,-\tfrac13)$ projects to $(\tfrac83,0)$ on the same line, beyond
the end where $x_3$ becomes negative.
:::

```{figure} figures/08-01-feasible-segment.svg
:name: fig-08-01-feasible-segment
:width: 460px
:class: dark-light
:alt: The line 3 x1 + 5 x2 = 8 in the (x1, x2) plane with a thick segment from (0, 1.6) to (1, 1) marked as the feasible set, and a hollow dot at (8/3, 0)

Projection of the feasible set of {prf:ref}`ex-08-01-feasible-segment` on the
$(x_1,x_2)$ plane: the thick segment of $3x_1+5x_2=8$ between the two basic
feasible solutions. The hollow dot is the infeasible basic solution
$(\tfrac83,0,-\tfrac13)$; along the line beyond $(1,1)$, $x_3$ is negative.
```

:::{prf:example} Is a given point a basic feasible solution?
:label: ex-08-01-bfs-test

Consider the system $\mA\vx=\vb$, $\vx\ge\vzero$ with

$$
\begin{gathered}
\mA=\begin{pmatrix}1&0&-1&3&-1\\0&1&1&4&2\\0&0&0&-7&3\end{pmatrix},\\
\vb=\begin{pmatrix}1\\4\\0\end{pmatrix}.
\end{gathered}
$$

**(a)** Is $\vx=(2,3,1,0,0)$ a basic feasible solution?
**(b)** Is $\vx=(1,4,0,0,0)$?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The rows of $\mA$ are linearly independent (the $3\times3$ submatrix of
columns 1, 2, 4 has determinant $-7\ne0$), so
$\rank\mA=3=m$ and {prf:ref}`thm-03-02-extreme-bfs` applies: a feasible point is
a BFS exactly when the columns of its positive components are linearly
independent.

**(a)** The point is feasible: $2-1=1$, $3+1=4$, $0=0$, and $\vx\ge\vzero$.
Its positive components are $x_1,x_2,x_3$, with columns $(1,0,0)$,
$(0,1,0)$ and $(-1,1,0)$. These are dependent:
$\mA_{\cdot3}=-\mA_{\cdot1}+\mA_{\cdot2}$. So $(2,3,1,0,0)$ is **not** a BFS.
Indeed, with $\vd=(1,-1,1,0,0)$, which satisfies $\mA\vd=\vzero$,

$$
(2,3,1,0,0)=\tfrac12(1,4,0,0,0)+\tfrac12(3,2,2,0,0),
$$

and both points on the right are feasible, so the point is not extreme.

**(b)** The point is feasible ($1=1$, $4=4$, $0=0$). Its positive components
$x_1,x_2$ have the independent columns $(1,0,0)$ and $(0,1,0)$, so it
**is** a BFS. It is **degenerate**: it has only two positive components
while $m=3$ ({prf:ref}`def-03-02-degenerate`). Its basis is not unique: the
columns can be completed by $\mA_{\cdot4}$ (determinant $-7$) or by
$\mA_{\cdot5}$ (determinant 3), and both bases give $\vx_B=(1,4,0)$. Adding
$\mA_{\cdot3}$ instead gives a singular matrix, since
$\mA_{\cdot3}=-\mA_{\cdot1}+\mA_{\cdot2}$.
:::

(sec-08-01-bases)=
## Bases and linear independence

:::{prf:example} Testing linear independence by pivoting
:label: ex-08-01-independence-pivoting

Are the columns of

$$
\mA=\begin{pmatrix}0&1&1&2\\1&0&1&2\\2&1&1&4\\0&1&1&2\\1&0&1&2\end{pmatrix}
$$

linearly independent? Use {prf:ref}`alg-03-02-independence`, and if they are
dependent, find a dependence relation.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Write $\vv_j=\mA_{\cdot j}$. The tableau has one row per column of $\mA$: the
left part $\gamma_1,\dots,\gamma_4$ starts as the identity, the right part
holds $\vv_j^\T$. The first pivot is marked in bold.

**Row 1.** Pivot on the second entry of the right part (value 1):

```{list-table}
:header-rows: 1

* - $\gamma_1$
  - $\gamma_2$
  - $\gamma_3$
  - $\gamma_4$
  - right part
* - 1
  - 0
  - 0
  - 0
  - $(0,\ \mathbf1,\ 2,\ 0,\ 1)$
* - 0
  - 1
  - 0
  - 0
  - $(1,\ 0,\ 1,\ 1,\ 0)$
* - 0
  - 0
  - 1
  - 0
  - $(1,\ 1,\ 1,\ 1,\ 1)$
* - 0
  - 0
  - 0
  - 1
  - $(2,\ 2,\ 4,\ 2,\ 2)$
```

After the pivot, rows 3 and 4 become $(-1,0,1,0\mid1,0,-1,1,0)$ and
$(-2,0,0,1\mid2,0,0,2,0)$.

**Row 2.** Its right part $(1,0,1,1,0)$ is nonzero; pivot on its first
entry (value 1). Rows 3 and 4 become $(-1,-1,1,0\mid0,0,-2,0,0)$ and
$(-2,-2,0,1\mid0,0,-2,0,0)$.

**Row 3.** Its right part $(0,0,-2,0,0)$ is nonzero; pivot on $-2$. The
final tableau is

```{list-table}
:header-rows: 1

* - $\gamma_1$
  - $\gamma_2$
  - $\gamma_3$
  - $\gamma_4$
  - right part
* - 0
  - $-1$
  - 1
  - 0
  - $(0,\ 1,\ 0,\ 0,\ 1)$
* - $-\tfrac12$
  - $\tfrac12$
  - $\tfrac12$
  - 0
  - $(1,\ 0,\ 0,\ 1,\ 0)$
* - $\tfrac12$
  - $\tfrac12$
  - $-\tfrac12$
  - 0
  - $(0,\ 0,\ 1,\ 0,\ 0)$
* - $-1$
  - $-1$
  - $-1$
  - 1
  - $(0,\ 0,\ 0,\ 0,\ 0)$
```

**Row 4.** Its right part is zero, so the columns are linearly
**dependent**, and the left part gives the relation

$$
-\vv_1-\vv_2-\vv_3+\vv_4=\vzero .
$$

**Check.** $-(0,1,2,0,1)-(1,0,1,1,0)-(1,1,1,1,1)+(2,2,4,2,2)=(0,0,0,0,0)$.
The first three columns are independent (three pivots succeeded), so
$\rank\mA=3$; the zero row means that $\vv_4=\vv_1+\vv_2+\vv_3$ adds nothing
new. In LP terms, these four columns could never all be basic together.
:::

:::{prf:example} Which column sets form a basis?
:label: ex-08-01-basis-test

For the system $x_1\vp_1+x_2\vp_2+x_3\vp_3+x_4\vp_4=\vb$ with

$$
\begin{gathered}
\bigl(\vp_1\ \ \vp_2\ \ \vp_3\ \ \vp_4\bigr)=\begin{pmatrix}1&0&1&2\\2&2&4&0\\3&1&2&0\end{pmatrix},\\
\vb=\begin{pmatrix}3\\4\\2\end{pmatrix},
\end{gathered}
$$

decide whether each of the following sets is a basis:
**(a)** $\{\vp_1,\vp_2,\vp_3\}$; **(b)** $\{\vp_1,\vp_3,\vp_4\}$;
**(c)** $\{\vp_2,\vp_3,\vp_4\}$; **(d)** $\{\vp_1,\vp_2,\vp_3,\vp_4\}$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

A basis of this system is a set of $m=3$ linearly independent columns,
that is, a $3\times3$ matrix with nonzero determinant.

**(a)** $\det(\vp_1\ \vp_2\ \vp_3)=1\cdot(2\cdot2-4\cdot1)-0+1\cdot(2\cdot1-2\cdot3)=-4\ne0$:
a basis. (Its basic solution is $\vx=(0,-4,3,0)$, which is infeasible.)

**(b)** $\det(\vp_1\ \vp_3\ \vp_4)=-16\ne0$ (expand along the third
column: $2\cdot(2\cdot2-4\cdot3)$): a basis. Its basic solution is
$\vx=(0,0,1,1)$, a degenerate basic feasible solution.

**(c)** $\det(\vp_2\ \vp_3\ \vp_4)=0$: not a basis. Indeed
$\vp_3=2\vp_2+\tfrac12\vp_4$.

**(d)** Not a basis: four vectors in $\R^3$ are always linearly dependent,
and a basis of this system has exactly three columns.
:::

:::{prf:example} Linear systems with a singular matrix
:label: ex-08-01-singular-systems

Let $\mB$ be a square matrix. True or false?

**(a)** $\mB\vx=\vb$ can have a unique solution even if $\mB$ is singular.
**(b)** If $\mB$ is singular and $\vb$ is not a linear combination of the
columns of $\mB$, then $\mB\vx=\vb$ has no solution.
**(c)** If $\mB$ is singular and $\vb$ is a linear combination of the
columns of $\mB$, then $\mB\vx=\vb$ has infinitely many solutions.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

A singular $\mB$ has a nonzero null-space vector $\vd$, $\mB\vd=\vzero$.

**(a) False.** If $\vx_0$ is a solution, so is $\vx_0+t\vd$ for every real
$t$; a singular system has either no solution or infinitely many, never
exactly one.

**(b) True.** $\mB\vx$ is always a combination of the columns of $\mB$, so
if $\vb$ is not in the column space, no $\vx$ solves the system.

**(c) True.** If $\vb$ is in the column space, there is a solution
$\vx_0$, and $\vx_0+t\vd$, $t\in\R$, are infinitely many more.

This is why a basis must be nonsingular: only then does
$\vx_B=\mB^{-1}\vb$ define one basic solution ({prf:ref}`def-03-02-basic-solution`).
:::

## References

- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Ch. 1: formulation examples; ch. 2:
  basic feasible solutions.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Ch. 3, §3.4: additional LP
  examples.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Ch. 2, §§2.2 and 2.4: graphical solution and LP applications, the
  source of several problem settings restated here; ch. 3, §§3.1–3.2; ch. 7.)
- Williams, H. P. (2013). *Model building in mathematical programming*
  (5th ed.). Wiley. (Chs. 1–5: blending, multi-period and refinery
  models.)
