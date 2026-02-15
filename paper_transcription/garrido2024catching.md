![](_page_0_Picture_1.jpeg)

## **Catching-Up Algorithm with Approximate Projections for Moreau's Sweeping Processes**

**Juan Guillermo Garrido1 · Emilio Vilches<sup>2</sup>**

Received: 14 August 2023 / Accepted: 7 February 2024 / Published online: 18 March 2024 © The Author(s), under exclusive licence to Springer Science+Business Media, LLC, part of Springer Nature 2024

#### **Abstract**

In this paper, we develop an enhanced version of the catching-up algorithm for sweeping processes through an appropriate concept of approximate projection. We establish some properties of this notion of approximate projection. Then, under suitable assumptions, we show the convergence of the enhanced catching-up algorithm for prox-regular, subsmooth, and merely closed sets. Finally, we briefly discuss some efficient numerical methods for obtaining approximate projections. Our results recover classical existence results in the literature and provide new insights into the numerical simulation of sweeping processes.

**Keywords** Sweeping process · Differential inclusions · Approximate projections

**Mathematics Subject Classification** 34A60 · 49J52 · 34G25 · 49J53

#### **1 Introduction**

Given a Hilbert space *H*, Moreau's sweeping process is a first-order differential inclusion involving the normal cone to a family of closed moving sets (*C*(*t*))*t*∈[0,*<sup>T</sup>* ]. In its simplest form, it can be written as

$$
\dot{x}(t) \in -N(C(t); x(t)) \quad \text{a.e. } t \in [0, T],
$$
  
\n
$$
x(0) = x_0 \in C(0),
$$
 (SP)

Communicated by Giovanni Colombo.

B Emilio Vilches emilio.vilches@uoh.cl

> Juan Guillermo Garrido jgarrido@dim.uchile.cl

<sup>1</sup> Departamento de Ingeniería Matemática, Universidad de Chile, Santiago, Chile

<sup>2</sup> Instituto de Ciencias de la Ingeniería, Universidad de O'Higgins, Rancagua, Chile

where *N*(*C*(*t*); ·) denotes an appropriate normal cone to the sets (*C*(*t*))*t*∈[0,*<sup>T</sup>* ]. Since its introduction by J.J. Moreau in [25, 26], the sweeping process has allowed the development of various applications in contact mechanics, electrical circuits, and crowd motion, among others (see, e.g., [1, 9, 24]). Furthermore, so far, we have a well-consolidated existence theory for moving sets in the considerable class of proxregular sets.

The most prominent (and constructive) method for solving the sweeping process is the so-called *catching-up algorithm*. Developed by J.J. Moreau in [26] for convex moving sets, it consists in taking a time discretization {*t n k* } *n <sup>k</sup>*=<sup>0</sup> of the interval [0, *<sup>T</sup>* ] and defining a piecewise linear and continuous function *xn* : [0, *T* ] → *H* with nodes

$$
x_{k+1}^n := \text{proj}_{C(t_{k+1}^n)}(x_k^n) \text{ for all } k \in \{0, \dots, n-1\}.
$$

Moreover, under general assumptions, it could be proved that the sequence (*xn*) converges to the unique solution of (SP) (see, e.g., [8]).

The applicability, from the numerical point of view, of the catching-up algorithm is based on the possibility of calculating an exact formula for the projection to the moving sets. However, for the majority of sets, the projection onto a closed set is not possible to obtain exactly, and only numerical approximations can be computed. Since there are still no guarantees on the convergence of the catching-up algorithm with approximate projections, in this paper, we develop a theoretical framework for the numerical approximation of the solutions of the sweeping process using an appropriate concept of approximate projection that is consistent with the numerical methods for the computation of the projection onto a closed set.

Regarding numerical approximations of sweeping processes, we are aware of the paper [33], where the author proposes an implementable numerical method for the particular case of the intersection of the complement of convex sets, which is used to study crowd motion. Our approach follows a different path and is based on numerical optimization methods to find an approximate projection in the following sense: given a closed set *C* ⊂ *H*, ε > 0 and *x* ∈ *H*, we say that *x*¯ ∈ *C* is an *approximate projection* of *C* at *x* ∈ *H* if

$$
||x - \bar{x}||^2 < \inf_{y \in C} ||x - y||^2 + \varepsilon.
$$

We observe that the set of approximate projections is always nonempty and can be obtained through numerical optimization methods. Hence, in this paper, we study the properties of approximate projections and propose a general numerical method for the sweeping process based on approximate projections. We prove that this algorithm converges in three general cases: (i) prox-regular moving sets (without compactness assumptions), (ii) ball-compact subsmooth moving sets, and (iii) general ball-compact fixed closed sets. Hence, our results cover a wide range of existence results for the sweeping process and provide important insights into the numerical simulation of sweeping processes.

The paper is organized as follows. Section 2 provides the mathematical tools needed for the presentation of the paper and also develops the theoretical properties of approximate projections. Section 3 is devoted to presenting the proposed algorithm and its main properties. Then, in Sect. 4, we prove the convergence of the algorithm when the moving set has uniformly prox-regular values (without compactness assumptions). Next, in Sect. 5, we provide the convergence of the proposed algorithm for ball-compact subsmooth moving sets. Section 6 shows the convergence for a fixed ball-compact set. Finally, Sect. 7 discusses numerical aspects for obtaining approximate projections. The paper ends with concluding remarks.

#### **2 Preliminaries**

From now on, *H* stands for a real Hilbert space, whose norm, denoted by ·, is induced by an inner product ·, ·. The closed (resp. open) ball centered at *x* with radius *<sup>r</sup>* <sup>&</sup>gt; 0 is denoted by <sup>B</sup>[*x*,*r*] (resp. <sup>B</sup>(*x*,*r*)), and the closed unit ball is denoted by <sup>B</sup>. The notation *<sup>H</sup>*<sup>w</sup> stands for *<sup>H</sup>* equipped with the weak topology, and *xn<sup>x</sup>* denotes the weak convergence of a sequence (*xn*) to *x*. For a given set *S* ⊂ *H*, the *support* and the *distance function* of *S* of at *x* ∈ *H* are defined, respectively, as

$$
\sigma(x, S) := \sup_{z \in S} \langle x, z \rangle \text{ and } d_S(x) := \inf_{z \in S} ||x - z||.
$$

Given ρ ∈]0, +∞] and γ < 1 positive, the ρ-*enlargement* and the γρ-*enlargement* of *S* are defined, respectively, as

$$
U_{\rho}(S) = \{x \in \mathcal{H} : d_S(x) < \rho\} \text{ and } U_{\rho}^{\gamma}(S) := \{x \in \mathcal{H} : d_S(x) < \gamma \rho\}.
$$

Given *A*, *B* ⊂ *H* two sets, we define the *excess* of *A* over *B* as the quantity *e*(*A*, *B*) := sup*x*∈*<sup>A</sup> dB*(*x*). From this, we define the *Hausdorff distance* between *A* and *B* as

$$
d_H(A, B) := \max\{e(A, B), e(B, A)\}.
$$

Further properties about Hausdorff distance can be found in [3, Sec. 3.16].

A vector *h* ∈ *H* belongs to the Clarke tangent cone *T* (*S*; *x*) (see [10]); when for every sequence (*xn*) in *S* converging to *x* and every sequence of positive numbers (*tn*) converging to 0, there exists a sequence (*hn*) in *H* converging to *h* such that *xn* <sup>+</sup> *tnhn* <sup>∈</sup> *<sup>S</sup>* for all *<sup>n</sup>* <sup>∈</sup> <sup>N</sup>. This cone is closed and convex, and its negative polar *N*(*S*; *x*) is the Clarke normal cone to *S* at *x* ∈ *S*, that is,

$$
N(S; x) := \{ v \in \mathcal{H} : \langle v, h \rangle \le 0 \text{ for all } h \in T(S; x) \}.
$$

As usual, *N*(*S*; *x*) = ∅ if *x* ∈/ *S*. Through that normal cone, the *Clarke subdifferential* of a function *<sup>f</sup>* : *<sup>H</sup>* <sup>→</sup> <sup>R</sup> ∪ {+∞} is defined by

$$
\partial f(x) := \{ v \in \mathcal{H} : (v, -1) \in N (\text{epi } f, (x, f(x))) \},
$$

where epi *<sup>f</sup>* := {(*y*,*r*) <sup>∈</sup> *<sup>H</sup>* <sup>×</sup> <sup>R</sup> : *<sup>f</sup>* (*y*) <sup>≤</sup> *<sup>r</sup>*} is the epigraph of *<sup>f</sup>* . When the function *f* is finite and locally Lipschitzian around *x*, the Clarke subdifferential is characterized (see [11]) in the following simple and amenable way

## 123

$$
\partial f(x) = \left\{ v \in \mathcal{H} : \langle v, h \rangle \le f^{\circ}(x; h) \text{ for all } h \in \mathcal{H} \right\},\
$$

where

$$
f^{\circ}(x; h) := \limsup_{(t,y)\to(0^+, x)} t^{-1} [f(y+th) - f(y)],
$$

is the *generalized directional derivative* of the locally Lipschitzian function *f* at *x* in the direction *h* ∈ *H*. The function *f* ◦(*x*; ·) is in fact the support of ∂ *f* (*x*), i.e., *f* ◦(*x*; *h*) = sup*z*∈<sup>∂</sup> *<sup>f</sup>* (*x*)*h*,*z*. That characterization easily yields that the Clarke subdifferential of any locally Lipschitzian function is a set-valued map with nonempty and convex values satisfying the important property of upper semicontinuity from *H* into *H*w.

Let *<sup>f</sup>* : *<sup>H</sup>* <sup>→</sup> <sup>R</sup> ∪ {+∞} be an lsc (*lower semicontinuous*) function and *<sup>x</sup>* <sup>∈</sup> dom *<sup>f</sup>* . We say that

(i) An element ζ belongs to the *proximal subdifferential* of *f* at *x*, denoted by ∂*<sup>P</sup> f* (*x*), if there exist two non-negative numbers σ and η such that

$$
f(y) \ge f(x) + \langle \zeta, y - x \rangle - \sigma ||y - x||^2
$$
 for all  $y \in \mathbb{B}(x; \eta)$ .

(ii) An element ζ ∈ *H* belongs to the *Fréchet subdifferential* of *f* at *x*, denoted by ∂*<sup>F</sup> f* (*x*), if

$$
\liminf_{h \to 0} \frac{f(x+h) - f(x) - \langle \zeta, h \rangle}{\|h\|} \ge 0.
$$

(iii) An element ζ belongs to the *limiting subdifferential* of *f* at *x*, denoted by ∂*<sup>L</sup> f* (*x*), if there exist sequences (ζ*n*) and (*xn*) such that <sup>ζ</sup>*<sup>n</sup>* <sup>∈</sup> <sup>∂</sup>*<sup>P</sup> <sup>f</sup>* (*xn*) for all *<sup>n</sup>* <sup>∈</sup> <sup>N</sup> and *xn* → *x*, ζ*n*ζ , and *f* (*xn*) → *f* (*x*).

Through these concepts, we can define the proximal, Fréchet, and limiting normal cone of a given set *S* ⊂ *H* at *x* ∈ *S*, respectively, as

$$
N^{P}(S; x) := \partial_{P} I_{S}(x), N^{F}(C; x) := \partial_{F} I_{C}(x) \text{ and } N^{L}(S; x) := \partial_{L} I_{S}(x),
$$

where *IS* is the indicator function of *S* ⊂ *H* (recall that *IS*(*x*) = 0 if *x* ∈ *S* and *IS*(*x*) = +∞ if *x* ∈/ *S*). It is well-known that (see [7, Theorem 4.1])

$$
N^{P}(S; x) \cap \mathbb{B} = \partial_{P} d_{S}(x) \quad \text{for all } x \in S. \tag{1}
$$

The equality (see [11])

$$
N(S; x) = \overline{\text{co}}^* N^L(S; x) = \text{cl}^* (\mathbb{R}_+ \partial d_S(x)) \quad \text{for } x \in S,
$$

gives an expression of the Clarke normal cone in terms of the distance function.

#### 123

Now, we recall the concept of uniformly prox-regular sets. Introduced by Federer in the finite-dimensional case (see [17]) and later developed by Rockafellar, Poliquin, and Thibault in [30], the prox-regularity generalizes and unifies convexity and nonconvex bodies with *C*<sup>2</sup> boundary. We refer to [12, 31] for a survey.

**Definition 1** Let *S* be a closed subset of *H* and ρ ∈]0, +∞]. The set *S* is called <sup>ρ</sup>-uniformly prox-regular if for all *<sup>x</sup>* <sup>∈</sup> *<sup>S</sup>* and <sup>ζ</sup> <sup>∈</sup> *<sup>N</sup> <sup>P</sup>*(*S*; *<sup>x</sup>*) one has

$$
\langle \zeta, x' - x \rangle \le \frac{\|\zeta\|}{2\rho} \|x' - x\|^2 \text{ for all } x' \in S.
$$

It is important to emphasize that convex sets are ρ-uniformly prox-regular for any ρ > 0. The following proposition provides a characterization of uniformly proxregular sets (see, e.g., [12, 27]).

**Proposition 1** *Let S* ⊂ *H be a closed set and* ρ ∈]0, +∞]*. The following assertions are equivalent:*

- *(a) S is* ρ*-uniformly prox-regular.*
- *(b) For any positive* γ < <sup>1</sup> *the mapping* proj*<sup>S</sup> is well-defined on U*<sup>γ</sup> <sup>ρ</sup> (*S*) *and Lipschitz continuous on U*<sup>γ</sup> <sup>ρ</sup> (*S*) *with* (<sup>1</sup> <sup>−</sup> γ )−<sup>1</sup> *as a Lipschitz constant, i.e.,*

$$
\|\text{proj}_S(u_1) - \text{proj}_S(u_2)\| \le (1 - \gamma)^{-1} \|u_1 - u_2\|
$$

*for all u*1, *<sup>u</sup>*<sup>2</sup> <sup>∈</sup> *<sup>U</sup>*<sup>γ</sup> <sup>ρ</sup> (*S*)*.*

*(c) For any xi* <sup>∈</sup> *<sup>S</sup>*, v*<sup>i</sup>* <sup>∈</sup> *<sup>N</sup> <sup>P</sup>* (*S*; *xi*)*, with i* <sup>=</sup> <sup>1</sup>, <sup>2</sup>*, one has*

$$
\langle v_1 - v_2, x_1 - x_2 \rangle \ge -\frac{1}{2\rho} (\|v_1\| + \|v_2\|) \|x_1 - x_2\|^2,
$$

*that is, the set-valued mapping N <sup>P</sup>*(*S*; ·) <sup>∩</sup> <sup>B</sup> *is* <sup>1</sup>/ρ*-hypomonotone. (d) For all* <sup>γ</sup> ∈]0, <sup>1</sup>[*, for all x*, *<sup>x</sup>* <sup>∈</sup> *<sup>U</sup>*<sup>γ</sup> <sup>ρ</sup> (*S*)*, for all* ξ ∈ ∂*PdS*(*x*)*, one has*

$$
\langle \xi, x' - x \rangle \le \frac{1}{2\rho(1 - \gamma)^2} \|x' - x\|^2 + d_S(x') - d_S(x).
$$

Next, we recall the class of subsmooth sets that includes the concepts of convex and uniformly prox-regular sets (see [4] and also [31, Chapter 8] for a survey).

**Definition 2** Let *S* be a closed subset of *H*. We say that *S* is *subsmooth* at *x*<sup>0</sup> ∈ *S*, if for every ε > 0 there exists δ > 0 such that

$$
\langle \xi_2 - \xi_1, x_2 - x_1 \rangle \ge -\varepsilon \| x_2 - x_1 \| \,, \tag{2}
$$

whenever *<sup>x</sup>*1, *<sup>x</sup>*<sup>2</sup> <sup>∈</sup> <sup>B</sup>[*x*0, δ] <sup>∩</sup> *<sup>S</sup>* and <sup>ξ</sup>*<sup>i</sup>* <sup>∈</sup> *<sup>N</sup>* (*S*; *xi*) <sup>∩</sup> <sup>B</sup> for *<sup>i</sup>* ∈ {1, <sup>2</sup>}. The set *<sup>S</sup>* is said *subsmooth* if it is subsmooth at each point of *S*. We further say that *S* is *uniformly subsmooth*, if for every ε > 0 there exists δ > 0, such that (2) holds for all *x*1, *x*<sup>2</sup> ∈ *S* satisfying *x*<sup>1</sup> <sup>−</sup> *<sup>x</sup>*2 <sup>≤</sup> <sup>δ</sup> and all <sup>ξ</sup>*<sup>i</sup>* <sup>∈</sup> *<sup>N</sup>* (*S*; *xi*) <sup>∩</sup> <sup>B</sup> for *<sup>i</sup>* ∈ {1, <sup>2</sup>}.

#### 123

Let (*S*(*t*))*t*∈*<sup>I</sup>* be a family of closed sets of *H* indexed by a nonempty set *I*. The family is called *equi-uniformly subsmooth*, if for all ε > 0, there exists δ > 0 such that for all *t* ∈ *I*, inequality (2) holds for all *x*1, *x*<sup>2</sup> ∈ *S*(*t*) satisfying *x*<sup>1</sup> − *x*2 ≤ δ and all <sup>ξ</sup>*<sup>i</sup>* <sup>∈</sup> *<sup>N</sup>*(*S*(*t*); *xi*) <sup>∩</sup> <sup>B</sup> with *<sup>i</sup>* ∈ {1, <sup>2</sup>}.

Given an interval *I*, a set-valued map *F* : *I* ⇒ *H* is said to be measurable if for all open set *<sup>U</sup>* of *<sup>H</sup>*, the inverse image *<sup>F</sup>*−1(*U*) = {*<sup>t</sup>* <sup>∈</sup> *<sup>I</sup>* : *<sup>F</sup>*(*t*) <sup>∩</sup> *<sup>U</sup>* = ∅} is a Lebesgue measurable set. When *F* takes nonempty and closed values and *H* is separable, this notion is equivalent to the *L* ⊗ *B*(*H*)-measurability of the graph gph *F* := {(*t*, *x*) ∈ *I* × *H* : *x* ∈ *F*(*t*)} (see, e.g., [28, Theorem 6.2.20]).

Given a set-valued map *F* : *H* ⇒ *H*, we say *F* is upper semicontinuous from *H* into *<sup>H</sup>*<sup>w</sup> if for all weakly closed set*<sup>C</sup>* of *<sup>H</sup>*, the inverse image *<sup>F</sup>*−1(*C*)is a closed set of *<sup>H</sup>*. It is known (see, e.g., see [28, Proposition 6.1.15 (c)]) that if *F* is upper semicontinuous, then the map *x* → σ (ξ , *F*(*x*)) is upper semicontinuous for all ξ ∈ *H*. When *F* takes convex and weakly compact values, these two properties are equivalent (see [28, Proposition 6.1.17]).

A set *<sup>S</sup>* <sup>⊂</sup> *<sup>H</sup>* is said ball compact if the set *<sup>S</sup>* <sup>∩</sup> *<sup>r</sup>*<sup>B</sup> is compact for all *<sup>r</sup>* <sup>&</sup>gt; 0. The *projection* onto *S* ⊂ *H* is the (possibly empty) set

$$
Proj_S(x) := \{ z \in S : d_S(x) = ||x - z|| \}.
$$

When the projection set is a singleton, we denote it as proj*S*(*x*). For ε > 0, we define the set of *approximate projections*:

$$
\text{proj}_{S}^{E}(x) := \left\{ z \in S : \|x - z\|^2 < d_{S}^{2}(x) + \varepsilon \right\}.
$$

By definition, the above set is nonempty and open. Moreover, it satisfies similar properties as the projection map (see Proposition 2 below). The approximate projections have been considered several times in variational analysis. In particular, they were used to characterize the subdifferential of the Asplund function of a given set. Indeed, let *S* ⊂ *H* and consider the Asplund function of the set *S*

$$
\varphi_S(x) := \frac{1}{2} ||x||^2 - \frac{1}{2} d_S^2(x) \quad x \in \mathcal{H}.
$$

Then, the following formula holds (see, e.g., [21, p. 467]):

$$
\partial \varphi_S(x) = \bigcap_{\varepsilon > 0} \overline{\text{co}}(\text{proj}_S^{\varepsilon}(x)).
$$

We recall that for any set *S* ⊂ *H* and *x* ∈ *H*, where Proj*S*(*x*) = ∅, the following formula is a consequence of formula (1):

$$
x - z \in d_S(x) \partial_P d_S(z) \text{ for all } z \in \text{Proj}_S(x).
$$

The next result provides an approximate version of the above formula for any closed set *S* ⊂ *H*.

#### 123

**Lemma 1** *Let S* <sup>⊂</sup> *<sup>H</sup> be a closed set, x* <sup>∈</sup> *<sup>H</sup>, and* ε > <sup>0</sup>*. For each z* <sup>∈</sup> proj<sup>ε</sup> *<sup>S</sup>*(*x*) *there is* <sup>v</sup> <sup>∈</sup> proj<sup>ε</sup> *<sup>S</sup>*(*x*) *such that z* − v < 2 <sup>√</sup><sup>ε</sup> *and*

$$
x - z \in (4\sqrt{\varepsilon} + d_S(x))\partial_P d_S(v) + 3\sqrt{\varepsilon} \mathbb{B}.
$$

*Proof* Fix ε > 0, *<sup>x</sup>* <sup>∈</sup> *<sup>H</sup>* and *<sup>z</sup>* <sup>∈</sup> proj<sup>ε</sup> *S*(*x*). According to the Borwein-Preiss Variational Principle [6, Theorem 2.6] applied to *<sup>y</sup>* → *<sup>g</sup>*(*y*) := *<sup>x</sup>* <sup>−</sup> *<sup>y</sup>*<sup>2</sup> <sup>+</sup> *IS*(*y*), there exists <sup>v</sup> <sup>∈</sup> proj<sup>ε</sup> *<sup>S</sup>*(*x*) such that *z* − v < 2 <sup>√</sup><sup>ε</sup> and 0 <sup>∈</sup> <sup>∂</sup>*<sup>P</sup> <sup>g</sup>*(v) <sup>+</sup> <sup>2</sup> <sup>√</sup>εB. Then, by the sum rule for the proximal subdifferential (see, e.g., [11, Proposition 2.11]), we obtain that

$$
x - v \in N^P(S; v) + \sqrt{\varepsilon} \mathbb{B},
$$

which implies that *<sup>x</sup>* <sup>−</sup> *<sup>z</sup>* <sup>∈</sup> *<sup>N</sup> <sup>P</sup>*(*S*; v) <sup>+</sup> <sup>3</sup> <sup>√</sup>εB. Next, since *<sup>x</sup>* <sup>−</sup> *<sup>z</sup>* ≤ *dS*(*x*) <sup>+</sup> <sup>√</sup>ε, we obtain that

$$
x - z \in N^{P}(S; v) \cap (4\sqrt{\varepsilon} + d_{S}(x))\mathbb{B} + 3\sqrt{\varepsilon} \mathbb{B}.
$$

Finally, the result follows from formula (1) and the above inclusion.

The following proposition displays some properties of approximation projections for uniformly prox-regular sets.

**Proposition 2** *Let S* ⊂ *H be a* ρ*-uniformly prox-regular set. Then, one has:*

- *(a) Let* (*xn*) *be a sequence converging to x* ∈ *U*ρ(*S*)*. Then for any* (*zn*) *and any sequence of positive numbers* (ε*n*) *converging to* <sup>0</sup> *with zn* <sup>∈</sup> projε*<sup>n</sup> <sup>S</sup>* (*xn*) *for all <sup>n</sup>* <sup>∈</sup> <sup>N</sup>*, we have that zn* <sup>→</sup> proj*S*(*x*)*.*
- *(b) Let* γ ∈]0, 1[ *and* ε ∈]0, ε0] *where* ε<sup>0</sup> *is such that*

$$
\gamma + 4\sqrt{\varepsilon_0} \left( 1 + \gamma + \frac{1}{\rho} (1 + 4\sqrt{\varepsilon_0}) \right) = 1.
$$

*Then, for all zi* <sup>∈</sup> proj<sup>ε</sup> *<sup>S</sup>*(*xi*) *with xi* <sup>∈</sup> *<sup>U</sup>*<sup>γ</sup> <sup>ρ</sup> (*S*) *for i* ∈ {1, 2}*, we have*

$$
(1 - F) \|z_1 - z_2\|^2 \le \sqrt{\varepsilon} \|x_1 - x_2\|^2 + M \sqrt{\varepsilon} + \langle x_1 - x_2, z_1 - z_2 \rangle,
$$

*where* - := <sup>α</sup> <sup>ρ</sup> +4 √ε <sup>1</sup> <sup>+</sup> <sup>α</sup> <sup>ρ</sup> <sup>+</sup> <sup>1</sup> <sup>ρ</sup> (<sup>1</sup> <sup>+</sup> <sup>√</sup>ε) *with* α := max{*dS*(*x*1), *dS*(*x*2)} *and M is a non-negative constant only dependent on* , ρ, γ *.*

*Proof* (*a*) We observe that for all *<sup>n</sup>* <sup>∈</sup> <sup>N</sup>

$$
||z_n|| \le ||z_n - x_n|| + ||x_n|| \le d_C(x_n) + \sqrt{\varepsilon_n} + ||x_n||.
$$

#### 123

Hence, since ε*<sup>n</sup>* → 0 and *xn* → *x*, we obtain (*yn*) is bounded. On the other hand, since *x* ∈ *U*ρ(*S*), we obtain proj*S*(*x*) is well-defined and

$$
||z_n - \text{proj}_S(x)||^2 = ||z_n - x_n||^2 - ||x_n - \text{proj}_S(x)||^2
$$
  
+ 2(x - \text{proj}\_S(x), z\_n - \text{proj}\_S(x)) + 2\langle z\_n - \text{proj}\_S(x), x\_n - x \rangle  
 
$$
\leq d_S^2(x_n) + \varepsilon_n - ||x_n - \text{proj}_S(x)||^2
$$
  
+ 2(x - \text{proj}\_S(x), z\_n - \text{proj}\_S(x)) + 2\langle z\_n - \text{proj}\_S(x), x\_n - x \rangle  
 
$$
\leq \varepsilon_n + 2\langle x - \text{proj}_S(x), z_n - \text{proj}_S(x) \rangle
$$
  
+ 2\langle z\_n - \text{proj}\_S(x), x\_n - x \rangle

where we have used *zn* <sup>∈</sup> projε*<sup>n</sup> <sup>S</sup>* (*xn*) and that *<sup>d</sup>*<sup>2</sup> *<sup>S</sup>*(*xn*) ≤ *xn* <sup>−</sup> proj*S*(*x*)2. Moreover, since *<sup>x</sup>* <sup>−</sup> proj*S*(*x*) <sup>∈</sup> *<sup>N</sup> <sup>P</sup>*(*S*; proj*S*(*x*)) and *<sup>S</sup>* is <sup>ρ</sup>-uniformly prox-regular, we obtain that

$$
2\langle x - \text{proj}_S(x), z_n - \text{proj}_S(x) \rangle \le \frac{d_S(x)}{\rho} \|z_n - \text{proj}_S(x)\|^2.
$$

Therefore, by using the above inequality and rearranging terms, we obtain that

$$
||z_n - \text{proj}_S(x)||^2 \leq \frac{\rho}{\rho - d_S(x)} \left( \varepsilon_n + 2\langle z_n - \text{proj}_S(x), x_n - x \rangle \right).
$$

Finally, since *xn* → *x* and (*zn*) is bounded, we conclude that *zn* → proj*S*(*x*).

(*b*) By virtue of Lemma 1, for *i* ∈ {1, 2} there exists v*i*, *bi* ∈ *H* such that

$$
b_i \in \mathbb{B}
$$
,  $v_i \in \text{proj}_S^{\varepsilon}(x_i)$ ,  $||z_i - v_i|| < 2\sqrt{\varepsilon}$  and  $\frac{x_i - z_i - 3\sqrt{\varepsilon}b_i}{4\sqrt{\varepsilon} + d_S(x_i)} \in \partial_P d_S(v_i)$ .

The hypomonotonicity of proximal normal cone (see Proposition 1 (b)) implies that

$$
\langle \zeta_1 - \zeta_2, v_1 - v_2 \rangle \ge \frac{-1}{\rho} ||v_1 - v_2||^2
$$

where <sup>ζ</sup>*<sup>i</sup>* := *xi*−*zi*−<sup>3</sup> <sup>√</sup>ε*bi* 4 <sup>√</sup>ε+<sup>α</sup> for *<sup>i</sup>* ∈ {1, <sup>2</sup>} and <sup>α</sup> := max{*dS*(*x*1), *dS*(*x*2)}. On the one hand, we have

$$
||v_1 - v_2|| \le ||v_1 - z_1|| + ||z_1 - z_2|| + ||z_2 - v_2|| \le 4\sqrt{\varepsilon} + ||z_1 - z_2||,
$$

and for all *z* ∈ *H* and *i* ∈ {1, 2}

$$
|\langle z, v_i - z_i \rangle| \le \frac{\sqrt{\varepsilon} \|z\|^2}{2} + \frac{\|v_i - z_i\|^2}{2\sqrt{\varepsilon}} \le \frac{\sqrt{\varepsilon} \|z\|^2}{2} + 2\sqrt{\varepsilon}.
$$

### 123

On the other hand,

$$
\langle (x_1 - z_1 - 3\sqrt{\varepsilon}b_1) - (x_2 - z_2 - 3\sqrt{\varepsilon}b_2), v_1 - v_2 \rangle
$$
  
=  $3\sqrt{\varepsilon} \langle b_2 - b_1, v_1 - v_2 \rangle + \langle (x_1 - x_2) - (z_1 - z_2), v_1 - v_2 \rangle$   
=  $3\sqrt{\varepsilon} \langle b_2 - b_1, v_1 - v_2 \rangle + \langle x_1 - x_2, v_1 - z_1 \rangle + \langle x_1 - x_2, z_1 - z_2 \rangle$   
+  $\langle x_1 - x_2, z_2 - v_2 \rangle - \langle z_1 - z_2, v_1 - z_1 \rangle - ||z_1 - z_2||^2 - \langle z_1 - z_2, z_2 - v_2 \rangle$   
 $\leq 6\sqrt{\varepsilon} \langle 4\sqrt{\varepsilon} + ||z_1 - z_2|| \rangle + \sqrt{\varepsilon} ||x_1 - x_2||^2 + 8\sqrt{\varepsilon} + \langle x_1 - x_2, z_1 - z_2 \rangle$   
-  $(1 - \sqrt{\varepsilon}) ||z_1 - z_2||^2$   
 $\leq 24\varepsilon + 11\sqrt{\varepsilon} + \sqrt{\varepsilon} ||x_1 - x_2||^2 + \langle x_1 - x_2, z_1 - z_2 \rangle - (1 - 4\sqrt{\varepsilon}) ||z_1 - z_2||^2.$ 

It follows that

$$
\left[1 - \frac{\alpha}{\rho} - 4\sqrt{\varepsilon}(1 + \frac{1}{\rho}(1 + 4\sqrt{\varepsilon} + \alpha))\right] \|z_1 - z_2\|^2
$$
  

$$
\leq \sqrt{\varepsilon} \|x_1 - x_2\|^2 + \langle x_1 - x_2, z_1 - z_2 \rangle + 4(4\varepsilon + \sqrt{\varepsilon})(4\frac{\sqrt{\varepsilon}}{\rho} + \gamma) + 24\varepsilon + 11\sqrt{\varepsilon}
$$

which proves the desired inequality.

The following result provides a stability result for a family of equi-uniformly subsmooth sets. We refer to see [20, Lemma 2.7] for a similar result.

**Lemma 2** *Let C* = {*Cn*}*n*∈<sup>N</sup> ∪ {*C*} *be a family of nonempty, closed, and equi-uniformly subsmooth sets. Assume that*

$$
\lim_{n\to\infty} d_{C_n}(x) = 0, \text{ for all } x \in C.
$$

*Then, for any sequence* <sup>α</sup>*<sup>n</sup>* <sup>→</sup> <sup>α</sup> <sup>∈</sup> <sup>R</sup> *and any sequence* (*yn*) *converging to y with yn* ∈ *Cn and y* ∈ *C, one has*

$$
\limsup_{n\to\infty}\sigma(\xi,\alpha_n\partial d_{C_n}(y_n))\leq \sigma(\xi,\alpha\partial d_C(y))\text{ for all }\xi\in\mathcal{H}.
$$

*Proof* Fix <sup>ξ</sup> <sup>∈</sup> *<sup>H</sup>*. Since <sup>∂</sup>*dS*(*x*) <sup>⊂</sup> <sup>B</sup> for all *<sup>x</sup>* <sup>∈</sup> *<sup>H</sup>*, we observe that

$$
\beta := \limsup_{n \to \infty} \sigma(\xi, \alpha_n \partial d_{C_n}(y_n)) < +\infty.
$$

Let us consider a subsequence (*nk* ) such that

$$
\beta = \lim_{k \to \infty} \sigma(\xi, \alpha_{n_k} \partial d_{C_{n_k}}(y_{n_k})).
$$

Given that <sup>∂</sup>*dCnk* (*ynk* ) is weakly compact for all *<sup>k</sup>* <sup>∈</sup> <sup>N</sup>, there is <sup>v</sup>*nk* <sup>∈</sup> <sup>∂</sup>*dCnk* (*ynk* ) such that

$$
\sigma(\xi, \alpha_{n_k} \partial d_{C_{n_k}}(y_{n_k})) = \langle \xi, \alpha_{n_k} v_{n_k} \rangle \text{ for all } k \in \mathbb{N}.
$$

123

$$
\qquad \qquad \Box
$$

Moreover, the sequence (v*nk* ) is bounded. Hence, without loss of generality, we can assume that <sup>v</sup>*nk*v <sup>∈</sup> <sup>B</sup>. It follows that <sup>β</sup> = ξ, αv. By equi-uniformly subsmoothness of *C*, for any ε > 0, there is δ > 0 such that for all *D* ∈ *C* and *x*1, *x*<sup>2</sup> ∈ *D* with *x*<sup>1</sup> − *x*2 < δ, one has

$$
\langle \zeta_1 - \zeta_2, x_1 - x_2 \rangle \ge -\varepsilon \|x_1 - x_2\|,\tag{3}
$$

whenever <sup>ζ</sup>*<sup>i</sup>* <sup>∈</sup> *<sup>N</sup>*(*D*; *xi*)∩<sup>B</sup> for *<sup>i</sup>* ∈ {1, <sup>2</sup>}. Next, let *<sup>y</sup>* <sup>∈</sup> *<sup>C</sup>* such that *<sup>y</sup>* <sup>−</sup> *<sup>y</sup>* < δ/2. Then, since *dCnk* (*y* ) converges to 0, there is a sequence (*y nk* ) converging to *y* with *y nk* <sup>∈</sup> *Cnk* for all *<sup>k</sup>* <sup>∈</sup> <sup>N</sup>. Hence, there is *<sup>k</sup>*<sup>0</sup> <sup>∈</sup> <sup>N</sup> such that *y nk* − *y* < δ/2 for all *k* ≥ *k*0. On the other hand, since *yn* → *y*, then there is *k* <sup>0</sup> <sup>∈</sup> <sup>N</sup> such that *ynk* <sup>−</sup> *<sup>y</sup>* < δ/<sup>2</sup> for all *k* ≥ *k* <sup>0</sup>. Hence, if *k* ≥ max{*k*0, *k* <sup>0</sup>} =: *k*ˆ we have *ynk* − *y nk* < δ. Therefore, it follows from the fact that 0 ∈ ∂*dCnk* (*y nk* ) and inequality (3) that

$$
\langle v_{n_k}, y_{n_k} - y'_{n_k} \rangle \ge -\varepsilon \| y_{n_k} - y'_{n_k} \| \text{ for all } k \ge \hat{k}.
$$

By taking *k* → ∞, we obtain that

$$
\langle v, y - y' \rangle \ge -\varepsilon \|y - y'\| \text{ for all } y' \in C \cap \mathbb{B}(y, \delta/2),
$$

which implies that <sup>v</sup> <sup>∈</sup> *<sup>N</sup> <sup>F</sup>* (*C*; *<sup>y</sup>*). Then, by [29, Lemma 4.21],

$$
v \in N^{F}(C; y) \cap \mathbb{B} = \partial_F d_C(y) \subset \partial d_C(y).
$$

Finally, we have proved that

$$
\beta = \langle \xi, \alpha v \rangle \le \sigma(\xi, \alpha \partial d_C(y)),
$$

which ends the proof.

The following lemma is a convergence theorem for a set-valued map from a topological space into a Hilbert space.

**Lemma 3** *Let* (*E*,τ) *be a topological space and G* : *E* ⇒ *H be a set-valued map with nonempty, closed, and convex values. Consider sequences* (*xn*) ⊂ *E,* (*yn*) ⊂ *H and* (ε*n*) <sup>⊂</sup> <sup>R</sup><sup>+</sup> *such that*

*(i) xn* → *x (in E), yny (weakly in H) and* ε*<sup>n</sup>* → 0*;*

- *(ii) For all n* <sup>∈</sup> <sup>N</sup>*, yn* <sup>∈</sup> *co*(*G*(*xk* ) <sup>+</sup> <sup>ε</sup>*k*<sup>B</sup> : *<sup>k</sup>* <sup>≥</sup> *<sup>n</sup>*)*;*
- *(iii)* lim sup *n*→∞ σ (ξ , *G*(*xn*)) ≤ σ (ξ , *G*(*x*)) *for all* ξ ∈ *H.*

*Then, y* ∈ *G*(*x*)*.*

*Proof* Assume by contradiction that *<sup>y</sup>* <sup>∈</sup>/ *<sup>G</sup>*(*x*). By virtue of Hahn–Banach theorem there exists <sup>ξ</sup> <sup>∈</sup> *<sup>H</sup>* \ {0}, δ > 0 and <sup>α</sup> <sup>∈</sup> <sup>R</sup> such that

$$
\langle \xi, y' \rangle + \delta \le \alpha \le \langle \xi, y \rangle, \ \forall y' \in \mathcal{G}(x).
$$

### 123

Then, it follows that σ (ξ , *G*(*x*)) ≤ α − δ. Besides, according to (ii) we have for all *<sup>n</sup>* <sup>∈</sup> <sup>N</sup>, there is a finite set *Jn* <sup>⊂</sup> <sup>N</sup> such that for all *<sup>m</sup>* <sup>∈</sup> *Jn*, *<sup>m</sup>* <sup>≥</sup> *<sup>n</sup>* and

$$
y_n = \sum_{j \in J_n} \alpha_j (y'_j + \varepsilon_j v_j)
$$

where for all *<sup>j</sup>* <sup>∈</sup> *Jn*, <sup>α</sup>*<sup>j</sup>* <sup>≥</sup> 0, <sup>v</sup> *<sup>j</sup>* <sup>∈</sup> <sup>B</sup>, *<sup>y</sup> <sup>j</sup>* ∈ *G*(*x <sup>j</sup>*) and *<sup>j</sup>*∈*Jn* <sup>α</sup>*<sup>j</sup>* <sup>=</sup> 1. Also, there exists *<sup>N</sup>* <sup>∈</sup> <sup>N</sup> such that for all *<sup>n</sup>* <sup>≥</sup> *<sup>N</sup>*, <sup>ε</sup>*<sup>n</sup>* <sup>&</sup>lt; <sup>δ</sup> <sup>2</sup>ξ . Thus, for *<sup>n</sup>* <sup>≥</sup> *<sup>N</sup>*

$$
\langle \xi, y_n \rangle = \sum_{j \in J_n} \alpha_j \langle \xi, y'_j + \varepsilon_j v_j \rangle
$$
  
\n
$$
\leq \sum_{j \in J_n} \alpha_j \sup_{k \geq n} \sigma(\xi, \mathcal{G}(x_k)) + \sum_{j \in J_n} \alpha_j \varepsilon_j \langle \xi, v_j \rangle
$$
  
\n
$$
\leq \sup_{k \geq n} \sigma(\xi, \mathcal{G}(x_k)) + ||\xi|| \sum_{j \in J_n} \alpha_j \frac{\delta}{2||\xi||} \leq \sup_{k \geq n} \sigma(\xi, \mathcal{G}(x_k)) + \frac{\delta}{2}.
$$

Therefore, as *yny*, letting *n* → ∞ in the last inequality we obtain that

$$
\langle \xi, y \rangle \leq \limsup_{n \to \infty} \sigma(\xi, \mathcal{G}(x_n)) + \frac{\delta}{2} \leq \sigma(\xi, \mathcal{G}(x)) + \frac{\delta}{2}.
$$

Therefore, ξ , *y* ≤ α − δ/2 ≤ ξ , *y* − δ/2, which is a contradiction. The proof is then complete.

The next lemma is a technical result whose proof can be found in [23, Lemma 2.2].

**Lemma 4** *Let* (*xn*) *be a sequence of absolutely continuous functions from* [0, *T* ] *into <sup>H</sup> with xn* (0) <sup>=</sup> *<sup>x</sup><sup>n</sup>* <sup>0</sup> *. Assume that for all n* <sup>∈</sup> <sup>N</sup>

$$
\|\dot{x}_n(t)\| \le \psi(t) \quad a.e \ t \in [0, T]
$$

*where* <sup>ψ</sup> <sup>∈</sup> *<sup>L</sup>*1([0, *<sup>T</sup>* ]; <sup>R</sup>+) *and that x<sup>n</sup>* <sup>0</sup> → *x*<sup>0</sup> *as n* → ∞*. Then, there exists a subsequence xnk of* (*xn*) *and an absolutely continuous function x such that*

*(i) xnk* (*t*)*x*(*t*) *in H as k* → +∞ *for all t* ∈ [0, *T* ]*. (ii) xnkx in L*1([0, *<sup>T</sup>* ]; *<sup>H</sup>*) *as k* → +∞*. (iii) <sup>x</sup>*˙*nkx in L* ˙ <sup>1</sup> ([0, *<sup>T</sup>* ]; *<sup>H</sup>*) *as k* → +∞*.*

*(iv)* ˙*x*(*t*) ≤ ψ(*t*) *a.e. t* ∈ [0, *T* ]*.*

#### **3 Catching-Up Algorithm with Errors for Sweeping Processes**

In this section, we propose a numerical method for the existence of solutions for the sweeping process:

$$
\dot{x}(t) \in -N(C(t); x(t)) + F(t, x(t)) \quad \text{a.e. } t \in [0, T],
$$
  
\n
$$
x(0) = x_0 \in C(0),
$$
\n(4)

123

where *C* : [0, *T* ] ⇒ *H* is a set-valued map with closed values in a Hilbert space *H*, *N* (*C*(*t*); *x*) stands for the Clarke normal cone to *C*(*t*) at *x*, and *F* : [0, *T* ] × *H* ⇒ *H* is a given set-valued map with nonempty closed and convex values. Our algorithm is based on the catching-up algorithm, except that we do not ask for an exact calculation of the projections.

The proposed algorithm is given as follows. For *<sup>n</sup>* <sup>∈</sup> <sup>N</sup>∗, let (*<sup>t</sup> n <sup>k</sup>* : *k* = 0, 1,..., *n*) be a uniform partition of [0, *T* ] with uniform time step μ*<sup>n</sup>* := *T* /*n*. Let (ε*n*) be a sequence of positive numbers such that ε*n*/μ<sup>2</sup> *<sup>n</sup>* → 0. We consider a sequence of piecewise continuous linear approximations (*xn*) defined as *xn*(0) = *x*<sup>0</sup> and for any *k* ∈ {0,..., *n* − 1} and *t* ∈]*t n <sup>k</sup>* , *t n <sup>k</sup>*+1]

$$
x_n(t) = x_k^n + \frac{t - t_k^n}{\mu_n} \left( x_{k+1}^n - x_k^n - \int_{t_k^n}^{t_{k+1}^n} f(s, x_k^n) \, ds \right) + \int_{t_k^n}^t f(s, x_k^n) \, ds, \tag{5}
$$

where *x<sup>n</sup>* <sup>0</sup> = *x*<sup>0</sup> and

$$
x_{k+1}^n \in \text{proj}_{C(t_{k+1}^n)}^{\varepsilon_n} \left( x_k^n + \int_{t_k^n}^{t_{k+1}^n} f(s, x_k^n) \, \text{d}s \right) \text{ for } k \in \{0, 1, \dots, n-1\}. \tag{6}
$$

Here *f* (*t*, *x*) denotes any selection of *F*(*t*, *x*) such that *f* (·, *x*) is measurable for all *<sup>x</sup>* <sup>∈</sup> *<sup>H</sup>*. For simplicity, we consider *<sup>f</sup>* (*t*, *<sup>x</sup>*) <sup>∈</sup> proj<sup>γ</sup> *F*(*t*,*x*) (0) for some γ > 0. In Proposition 3, we prove that it is possible to obtain such measurable selection under mild assumptions.

The above algorithm is called *catching-up algorithm with approximate projections* because the projection is not necessarily exactly calculated. We will prove that the above algorithm converges for several families of algorithms as long as inclusion (6) is verified.

Let us consider functions δ*n*(·) and θ*n*(·) defined as

$$
\delta_n(t) = \begin{cases} t_k^n & \text{if } t \in [t_k^n, t_{k+1}^n[ \\ t_{n-1}^n & \text{if } t = T, \end{cases} \text{ and } \theta_n(t) = \begin{cases} t_{k+1}^n & \text{if } t \in [t_k^n, t_{k+1}^n[ \\ T & \text{if } t = T. \end{cases}
$$

In what follows, we show useful properties satisfied for the above algorithm, which will help us to prove the existence of sweeping process (4) in three cases:

- (i) The set-valued map *t* ⇒ *C*(*t*) takes uniformly prox-regular values.
- (ii) The set-valued map *t* ⇒ *C*(*t*) takes subsmooth and ball-compact values.
- (iii) *C*(*t*) ≡ *C* in [0, *T* ] and *C* is ball-compact.

Throughout this section, *F* : [0, *T* ] × *H* ⇒ *H* will be a set-valued map with nonempty, closed, and convex values. Moreover, we will consider the following conditions:

(*H<sup>F</sup>* <sup>1</sup> ) For all *t* ∈ [0, *T* ], *F*(*t*, ·) is upper semicontinuous from *H* into *H*w.

123

(*H<sup>F</sup>* <sup>2</sup> ) There exists *<sup>h</sup>* : *<sup>H</sup>* <sup>→</sup> <sup>R</sup><sup>+</sup> Lipschitz continuous (with constant *Lh* <sup>&</sup>gt; 0) such that

$$
d(0, F(t, x)) := \inf\{\|w\| : w \in F(t, x)\} \le h(x),
$$

for all *x* ∈ *H* and a.e. *t* ∈ [0, *T* ].

(*H<sup>F</sup>* <sup>3</sup> ) There is γ > 0 such that the set-valued map (*t*, *<sup>x</sup>*) <sup>⇒</sup> proj<sup>γ</sup> *F*(*t*,*x*) (0) has a selection *f* : [0, *T* ] × *H* → *H* such that *f* (·, *x*) is measurable for all *x* ∈ *H*.

The following proposition provides conditions for the feasibility of hypothesis (*H<sup>F</sup>* <sup>3</sup> ).

**Proposition 3** *Let us assume that H is a separable Hilbert space. Moreover we suppose <sup>F</sup>*(·, *<sup>x</sup>*) *is measurable for all x* <sup>∈</sup> *<sup>H</sup>; then,* (*H<sup>F</sup>* <sup>3</sup> ) *holds for all* γ > 0*.*

*Proof* Let γ > 0 and fix *<sup>x</sup>* <sup>∈</sup> *<sup>H</sup>*. Since the set-valued map *<sup>F</sup>*(·, *<sup>x</sup>*) is measurable, the map *t* → *d*(0, *F*(*t*, *x*)) is a measurable function. Let us define the set-valued map *<sup>F</sup><sup>x</sup>* : *<sup>t</sup>* <sup>⇒</sup> proj<sup>γ</sup> *F*(*t*,*x*) (0). Then,

$$
\begin{aligned} \text{gph}\,\mathcal{F}_x &= \{(t, y) \in [0, T] \times \mathcal{H} : y \in \text{proj}_{F(t, x)}^{\mathcal{V}}(0)\} \\ &= \{(t, y) \in [0, T] \times \mathcal{H} : \|y\|^2 < d(0, F(t, x))^2 + \gamma \text{ and } y \in F(t, x)\} \\ &= \text{gph}\, F(\cdot, x) \cap \{(t, y) \in [0, T] \times \mathcal{H} : \|y\|^2 < d(0, F(t, x))^2 + \gamma\}. \end{aligned}
$$

Hence, gph *F<sup>x</sup>* is a measurable set. Consequently, *F<sup>x</sup>* has a measurable selection (see [28, Theorem 6.3.20]). Denoting by *t* → *f* (*t*, *x*) such selection, we obtain the result. 

Now, we establish the main properties of the proposed algorithm.

**Theorem 1** *Assume, in addition to* (*H<sup>F</sup>* <sup>1</sup> )*,* (*H<sup>F</sup>* <sup>2</sup> ) *and* (*H<sup>F</sup>* <sup>3</sup> )*, that C* : [0, *T* ] ⇒ *H is a set-valued map with nonempty and closed values such that*

$$
d_H(C(t), C(s)) \le L_C |t - s| \text{ for all } t, s \in [0, T]. \tag{7}
$$

*Then, the sequence of functions* (*xn* : [0, *T* ] → *H*) *generated by numerical scheme* (5) *and* (6) *satisfies the following properties:*

- *(a) There are non-negative constants K*1, *<sup>K</sup>*2, *<sup>K</sup>*3, *<sup>K</sup>*4, *<sup>K</sup>*<sup>5</sup> *such that for all n* <sup>∈</sup> <sup>N</sup> *and t* ∈ [0, *T* ]*:*
	- *(i) dC*(θ*<sup>n</sup>* (*<sup>t</sup>*))(*xn*(δ*n*(*t*))+ <sup>θ</sup>*<sup>n</sup>* (*t*) <sup>δ</sup>*<sup>n</sup>* (*t*) *<sup>f</sup>* (*s*, *xn*(δ*n*(*t*)))d*s*) <sup>≤</sup> (*LC* <sup>+</sup>*h*(*x*(δ*n*(*t*)))+√γ )μ*n*.
	- *(ii) xn*(θ*n*(*t*)) − *x*0 ≤ *K*1.
	- *(iii) xn*(*t*) ≤ *K*2.
	- *(iv) xn*(θ*n*(*t*)) <sup>−</sup> *xn*(δ*n*(*t*)) ≤ *<sup>K</sup>*3μ*<sup>n</sup>* <sup>+</sup> <sup>√</sup>ε*n*.
	- *(v) xn*(*t*) − *xn*(θ*n*(*t*)) ≤ *K*4μ*<sup>n</sup>* + 2 <sup>√</sup>ε*n.*
- *(b) There exists K*<sup>5</sup> <sup>&</sup>gt; <sup>0</sup> *such that for all t* ∈ [0, *<sup>T</sup>* ] *and m*, *<sup>n</sup>* <sup>∈</sup> <sup>N</sup> *we have*

$$
d_{C(\theta_n(t))}(x_m(t)) \leq K_5\mu_m + L_C\mu_n + 2\sqrt{\varepsilon_m}.
$$

123

- *(c) There exists K*<sup>6</sup> <sup>&</sup>gt; <sup>0</sup> *such that for all n* <sup>∈</sup> <sup>N</sup> *and almost all t* ∈ [0, *<sup>T</sup>* ]*,* ˙*xn*(*t*) ≤ *K*6*.*
- *(d) For all n* <sup>∈</sup> <sup>N</sup> *and k* ∈ {0, <sup>1</sup>,..., *<sup>n</sup>* <sup>−</sup> <sup>1</sup>}*, there is* <sup>v</sup>*<sup>n</sup> <sup>k</sup>*+<sup>1</sup> <sup>∈</sup> *<sup>C</sup>*(*<sup>t</sup> n <sup>k</sup>*+1) *such that for all t* ∈]*t n <sup>k</sup>* , *t n <sup>k</sup>*+1[*:*

$$
\dot{x}_n(t) \in -\frac{\lambda_n(t)}{\mu_n} \partial_P d_{C(\theta_n(t))}(v_{k+1}^n) + f(t, x_n(\delta_n(t))) + \frac{3\sqrt{\varepsilon_n}}{\mu_n} \mathbb{B},\tag{8}
$$

*where* λ*n*(*t*) = 4 <sup>√</sup>ε*<sup>n</sup>* <sup>+</sup> (*LC* <sup>+</sup> *<sup>h</sup>*(*x*(δ*n*(*t*))) <sup>+</sup> <sup>√</sup>γ )μ*n. Moreover,* v*<sup>n</sup> <sup>k</sup>*+<sup>1</sup> <sup>−</sup> *xn*(θ*n*(*t*)) <sup>&</sup>lt; <sup>2</sup> <sup>√</sup>ε*n.*

*Proof* (*a*): Setμ*<sup>n</sup>* := *<sup>T</sup>* /*<sup>n</sup>* and let(ε*n*) be a sequence of non-negative numbers such that ε*n*/μ<sup>2</sup> *<sup>n</sup>* → 0. We define c := sup*n*∈<sup>N</sup> <sup>√</sup>ε*<sup>n</sup>* <sup>μ</sup>*<sup>n</sup>* . We denote by *Lh* the Lipschitz constant of *h*. For all *<sup>t</sup>* ∈ [0, *<sup>T</sup>* ] and *<sup>n</sup>* <sup>∈</sup> <sup>N</sup>, we define <sup>τ</sup>*n*(*t*) := *xn*(δ*n*(*t*)) <sup>+</sup> <sup>θ</sup>*<sup>n</sup>* (*t*) <sup>δ</sup>*<sup>n</sup>* (*t*) *f* (*s*, *xn*(δ*n*(*t*)))d*s*. Since *<sup>f</sup>* (*t*, *xn*(δ*n*(*t*))) <sup>∈</sup> proj<sup>γ</sup> *<sup>F</sup>*(*t*,*xn* (δ*<sup>n</sup>* (*t*)))(0) we obtain that

$$
d_{C(\theta_n(t))}(\tau_n(t)) \leq d_{C(\theta_n(t))}(x_n(\delta_n(t))) + \left\| \int_{\delta_n(t)}^{\theta_n(t)} f(s, x_n(\delta_n(t))) ds \right\|
$$
  
\n
$$
\leq L_C \mu_n + \int_{\delta_n(t)}^{\theta_n(t)} \|f(s, x_n(\delta_n(t)))\| ds
$$
  
\n
$$
\leq L_C \mu_n + \int_{\delta_n(t)}^{\theta_n(t)} (h(x_n(\delta_n(t))) + \sqrt{\gamma}) ds
$$
  
\n
$$
\leq (L_C + h(x_n(\delta_n(t))) + \sqrt{\gamma}) \mu_n,
$$

which proves (*i*). Moreover, since *xn*(θ*n*(*t*)) <sup>∈</sup> projε*<sup>n</sup> <sup>C</sup>*(θ*<sup>n</sup>* (*t*))(τ*n*(*t*)), we get that

$$
||x_n(\theta_n(t)) - \tau_n(t)|| \leq d_{C(\theta_n(t))}(\tau_n(t)) + \sqrt{\varepsilon_n}
$$
  
\n
$$
\leq (L_C + h(x_n(\delta_n(t))) + \sqrt{\gamma})\mu_n + \sqrt{\varepsilon_n},
$$
\n(9)

which yields

$$
||x_n(\theta_n(t)) - x_n(\delta_n(t))|| \le (L_C + 2h(x_n(\delta_n(t))) + 2\sqrt{\gamma})\mu_n + \sqrt{\varepsilon_n}
$$
  
\n
$$
\le (L_C + 2h(x_0) + 2\sqrt{\gamma} + 2L_h ||x_n(\delta_n(t)) - x_0||)\mu_n(10)
$$
  
\n
$$
+ \sqrt{\varepsilon_n}.
$$

Hence, for all *t* ∈ [0, *T* ]

$$
||x_n(\theta_n(t)) - x_0|| \le (1 + 2L_h\mu_n) ||x_n(\delta_n(t)) - x_0|| + (L_C + 2h(x_0) + 2\sqrt{\gamma})\mu_n + \sqrt{\varepsilon_n}.
$$

The above inequality means that for all *k* ∈ {0, 1,..., *n* − 1}:

$$
||x_{k+1}^n - x_0|| \le (1 + 2L_h\mu_n) ||x_k^n - x_0|| + (L_C + 2h(x_0) + 2\sqrt{\gamma})\mu_n + \sqrt{\varepsilon_n}.
$$

#### 123

Then, by [11, p. 183], we obtain that for all *k* ∈ {0,..., *n* − 1}

$$
||x_{k+1}^n - x_0|| \le (k+1)((L_C + 2h(x_0) + 2\sqrt{\gamma})\mu_n + \sqrt{\varepsilon_n}) \exp(2L_h(k+1)\mu_n)
$$
  
\n
$$
\le T(L_C + 2h(x_0) + \sqrt{\gamma} + \mathfrak{c}) \exp(2L_h T) =: K_1.
$$

which proves (*ii*).

(*iii*): By definition of *xn*, for *t* ∈]*t n <sup>k</sup>* , *t n <sup>k</sup>*+1] and *<sup>k</sup>* ∈ {0, <sup>1</sup> ..., *<sup>n</sup>* <sup>−</sup> <sup>1</sup>}, using (5)

$$
||x_n(t)|| \le ||x_k^n|| + ||x_{k+1}^n - \tau_n(t)|| + \int_{t_k^n}^t ||f(s, x_k^n)|| ds
$$
  
 
$$
\le K_1 + ||x_0|| + (L_C + \sqrt{\gamma} + h(x_k^n))\mu_n + \sqrt{\varepsilon_n} + (h(x_k^n) + \sqrt{\gamma})\mu_n,
$$

where we have used (9). Moreover, it is clear that for *k* ∈ {0,..., *n*}

$$
h(x_k^n) \le h(x_0) + L_h \|x_k^n - x_0\| \le h(x_0) + L_h K_1.
$$

Therefore, for all *t* ∈ [0, *T* ]

$$
||x_n(t)|| \le K_1 + ||x_0|| + (L_C + 2(h(x_0) + L_h K_1 + \sqrt{\gamma}))\mu_n + \sqrt{\varepsilon_n}
$$
  
 
$$
\le K_1 + ||x_0|| + T(L_C + 2(h(x_0) + L_h K_1 + \sqrt{\gamma}) + \mathfrak{c}) =: K_2,
$$

which proves (*iii*).

(*i*v): From (10) and (11) it is easy to see that there exists *K*<sup>3</sup> > 0 such that for all *<sup>n</sup>* <sup>∈</sup> <sup>N</sup> and *<sup>t</sup>* ∈ [0, *<sup>T</sup>* ]: *xn*(θ*n*(*t*)) <sup>−</sup> *xn*(δ*n*(*t*)) ≤ *<sup>K</sup>*3μ*<sup>n</sup>* <sup>+</sup> <sup>√</sup>ε*n*.

(v): To conclude this part, we consider *t* ∈]*t n <sup>k</sup>* , *t n <sup>k</sup>*+1] for some *<sup>k</sup>* ∈ {0, <sup>1</sup>,..., *<sup>n</sup>*−1}. Then *xn*(θ*n*(*t*)) <sup>=</sup> *<sup>x</sup><sup>n</sup> <sup>k</sup>*+<sup>1</sup> and also

$$
||x_n(\theta_n(t)) - x_n(t)|| \le ||x_{k+1}^n - x_k^n|| + ||x_{k+1}^n - \tau_n(t)|| + \int_{t_k^n}^t ||f(s, x_k^n)||ds
$$
  
\n
$$
\le K_3 \mu_n + \sqrt{\varepsilon_n} + (L_C + \sqrt{\gamma} + h(x_0) + L_h K_1)\mu_n + \sqrt{\varepsilon_n}
$$
  
\n
$$
+ \mu_n (h(x_k^n) + \sqrt{\gamma})
$$
  
\n
$$
\le (\underline{K_3 + L_C + 2(h(x_0) + L_h K_1) + 2\sqrt{\gamma}})\mu_n + 2\sqrt{\varepsilon_n},
$$
  
\n
$$
=:K_4
$$

and we conclude this first part. (*b*): Let *<sup>m</sup>*, *<sup>n</sup>* <sup>∈</sup> <sup>N</sup> and *<sup>t</sup>* ∈ [0, *<sup>T</sup>* ], then

$$
d_{C(\theta_n(t))}(x_m(t)) \leq d_{C(\theta_n(t))}(x_m(\theta_m(t))) + ||x_m(\theta_m(t)) - x_m(t)||
$$
  
\n
$$
\leq d_H(C(\theta_n(t)), C(\theta_m(t))) + K_4\mu_m + 2\sqrt{\varepsilon_m}
$$
  
\n
$$
\leq L_C|\theta_n(t) - \theta_m(t)| + K_4\mu_m + 2\sqrt{\varepsilon_m}
$$
  
\n
$$
\leq L_C(\mu_n + \mu_m) + K_4\mu_m + 2\sqrt{\varepsilon_m}
$$

#### 123

where we have used (v). Hence, by setting *K*<sup>5</sup> := *K*<sup>4</sup> + *LC* we prove (b). (*c*): Let *<sup>n</sup>* <sup>∈</sup> <sup>N</sup>, *<sup>k</sup>* ∈ {0, <sup>1</sup>,..., *<sup>n</sup>* <sup>−</sup> <sup>1</sup>} and *<sup>t</sup>* ∈]*<sup>t</sup> n <sup>k</sup>* , *t n <sup>k</sup>*+1]. Then,

$$
\|\dot{x}_n(t)\| = \left\| \frac{1}{\mu_n} \left( x_{k+1}^n - x_k^n - \int_{t_k^n}^{t_{k+1}^n} f(s, x_k^n) ds \right) + f(t, x_k^n) \right\|
$$
  
\n
$$
\leq \frac{1}{\mu_n} \|x_n(\theta_n(t)) - \tau_n(t)\| + \|f(t, x_k^n)\|
$$
  
\n
$$
\leq \frac{1}{\mu_n} ((L_C + h(x_k^n) + \sqrt{\gamma})\mu_n + \sqrt{\varepsilon_n}) + h(x_k^n) + \sqrt{\gamma}
$$
  
\n
$$
\leq \frac{\sqrt{\varepsilon_n}}{\mu_n} + L_C + 2(h(x_0) + L_h K_1 + \sqrt{\gamma})
$$
  
\n
$$
\leq c + L_C + 2(h(x_0) + L_h K_1 + \sqrt{\gamma}) =: K_6,
$$

which proves (*c*).

(*d*): Fix *k* ∈ {0, 1,..., *n* − 1} and *t* ∈]*t n <sup>k</sup>* , *t n <sup>k</sup>*+1[. Then, *<sup>x</sup><sup>n</sup> <sup>k</sup>*+<sup>1</sup> <sup>∈</sup> projε*<sup>n</sup> C*(*t n <sup>k</sup>*+1) (τ*n*(*t*)). Hence, by Lemma 1, there exists v*<sup>n</sup> <sup>k</sup>*+<sup>1</sup> <sup>∈</sup> *<sup>C</sup>*(*<sup>t</sup> n <sup>k</sup>*+1) such that *xk*+<sup>1</sup> <sup>−</sup> <sup>v</sup>*<sup>n</sup> <sup>k</sup>*+1 <sup>&</sup>lt; <sup>2</sup> <sup>√</sup>ε*<sup>n</sup>* and

$$
\tau_n(t) - x_{k+1}^n \in \alpha_n(t) \partial \rho d_{C(t_{k+1}^n)}(v_{k+1}^n) + 3\sqrt{\varepsilon_n} \mathbb{B}, \ \forall t \in ]t_k^n, t_{k+1}^n[
$$

where α*n*(*t*) = 4 <sup>√</sup>ε*<sup>n</sup>* <sup>+</sup> *dC*(θ*<sup>n</sup>* (*<sup>t</sup>*))(τ*n*(*t*)). By virtue of (*i*),

$$
\alpha_n(t) \leq 4\sqrt{\varepsilon_n} + (L_C + h(x(\delta_n(t))) + \sqrt{\gamma})\mu_n =: \lambda_n(t).
$$

Then, for all *t* ∈]*t n <sup>k</sup>* , *t n <sup>k</sup>*+1[

$$
-\mu_n(\dot{x}_n(t)-f(t,x_k^n))\in\lambda_n(t)\partial_Pd_{C(t_{k+1}^n)}(v_{k+1}^n)+3\sqrt{\varepsilon_n}\mathbb{B},
$$

which implies that *t* ∈]*t n <sup>k</sup>* , *t n <sup>k</sup>*+1[

$$
\dot{x}_n(t) \in -\frac{\lambda_n(t)}{\mu_n} \partial_P d_{C(t_{k+1}^n)}(v_{k+1}^n) + f(t, x_k^n) + \frac{3\sqrt{\varepsilon_n}}{\mu_n} \mathbb{B}.
$$

#### **4 Prox-Regular Case**

In this section, we will study the algorithm under the assumption of uniform proxregularity of the moving sets. The classical catching-up algorithm in this framework was studied in [8], where the existence of solutions for (4) was established for a set-valued map *F* taking values in a fixed compact set.

#### 123

**Theorem 2** *Suppose, in addition to the assumptions of Theorem* 1*, that C*(*t*) *is* ρ*uniformly prox-regular for all t* ∈ [0, *T* ]*, and for all r* > 0*, there exists a non-negative integrable function kr such that for all t* ∈ [0, *<sup>T</sup>* ] *and x*, *<sup>x</sup>* <sup>∈</sup> *<sup>r</sup>*<sup>B</sup> *one has*

$$
\langle y - y', x - x' \rangle \le k_r(t) \|x - x'\|^2, \ \forall y \in F(t, x), \forall y' \in F(t, x'). \tag{12}
$$

*Then, the sequence of functions* (*xn*) *generated by algorithm* (5) *and* (6) *converges uniformly to an absolutely continuous function x, which is a solution of* (4)*. Moreover, if F satisfies the following growth condition,*

$$
\sup_{y \in F(t,x)} \|y\| \le c(t) (\|x\| + 1), \forall x \in \mathcal{H}, t \in [0, T],
$$
\n(13)

*where c* <sup>∈</sup> *<sup>L</sup>*1([0, *<sup>T</sup>* ]; <sup>R</sup>+)*, then the solution x is unique.*

*Proof* Consider *<sup>m</sup>*, *<sup>n</sup>* <sup>∈</sup> <sup>N</sup> with *<sup>m</sup>* <sup>≥</sup> *<sup>n</sup>* big enough such that for all *<sup>t</sup>* ∈ [0, *<sup>T</sup>* ], *dC*(θ*<sup>n</sup>* (*<sup>t</sup>*))(*xm*(*t*)) < ρ, this can be guaranteed by Theorem 1. Then, for a.e. *t* ∈ [0, *T* ]

$$
\frac{\mathrm{d}}{\mathrm{d}t}\left(\frac{1}{2}\|x_n(t)-x_m(t)\|^2\right)=\langle\dot{x}_n(t)-\dot{x}_m(t),x_n(t)-x_m(t)\rangle.
$$

Let *t* ∈ [0, *T* ] where the above equality holds. Let *k*, *j* ∈ {0, 1,..., *n* − 1} such that *t* ∈]*t n <sup>k</sup>* , *t n <sup>k</sup>*+1] and *<sup>t</sup>* ∈]*<sup>t</sup> m <sup>j</sup>* , *t m <sup>j</sup>*+1]. On the one hand, we have that

$$
\langle \dot{x}_n(t) - \dot{x}_m(t), x_n(t) - x_m(t) \rangle = \langle \dot{x}_n(t) - \dot{x}_m(t), x_n(t) - x_{k+1}^n \rangle \n+ \langle \dot{x}_n(t) - \dot{x}_m(t), x_{k+1}^n - v_{k+1}^n \rangle \n+ \langle \dot{x}_n(t) - \dot{x}_m(t), v_{k+1}^n - v_{j+1}^m \rangle \n+ \langle \dot{x}_n(t) - \dot{x}_m(t), v_{j+1}^m - x_{j+1}^m \rangle \n+ \langle \dot{x}_n(t) - \dot{x}_m(t), x_{j+1}^m - x_m(t) \rangle \n\le 2K_6(K_4(\mu_n + \mu_m) + 4(\sqrt{\varepsilon_n} + \sqrt{\varepsilon_m})) \n+ \langle \dot{x}_n(t) - \dot{x}_m(t), v_{k+1}^n - v_{j+1}^m \rangle,
$$
\n(14)

where v*<sup>n</sup> <sup>k</sup>*+<sup>1</sup> <sup>∈</sup> *<sup>C</sup>*(*<sup>t</sup> n <sup>k</sup>*+1) and <sup>v</sup>*<sup>m</sup> <sup>j</sup>*+<sup>1</sup> <sup>∈</sup> *<sup>C</sup>*(*<sup>t</sup> m <sup>j</sup>*+1) are the given in Theorem 1. We can see that

$$
\max \left\{ d_{C(t_{k+1}^n)}(v_{j+1}^m), d_{C(t_{j+1}^m)}(v_{k+1}^n) \right\} \le d_H(C(t_{j+1}^m), C(t_{k+1}^n))
$$
  

$$
\le L_C |t_{j+1}^m - t_{k+1}^n| \le L_C(\mu_n + \mu_m).
$$

From now, *<sup>m</sup>*, *<sup>n</sup>* <sup>∈</sup> <sup>N</sup> are big enough such that *LC*(μ*<sup>n</sup>* <sup>+</sup> <sup>μ</sup>*m*) < <sup>ρ</sup> <sup>2</sup> . Moreover, as *h* is *Lh*-Lipschitz, we have that for all *<sup>p</sup>* <sup>∈</sup> <sup>N</sup>, *<sup>i</sup>* ∈ {0, <sup>1</sup>,..., *<sup>p</sup>*} and *<sup>t</sup>* ∈ [0, *<sup>T</sup>* ]

$$
|| f(t, x_i^p) || \le h(x_i^p) + \sqrt{\gamma} \le h(x_0) + L_h K_1 + \sqrt{\gamma} =: \alpha.
$$

#### 123

On the other hand, using (8) and Proposition 1 we have that

$$
\frac{1}{F} \max\{\left\langle \zeta_n - \dot{x}_n(t), v_{j+1}^m - v_{k+1}^n \right\rangle, \left\langle \zeta_m - \dot{x}_m(t), v_{k+1}^n - v_{j+1}^m \right\rangle\} \n\leq \frac{2}{\rho} \|v_{k+1}^n - v_{j+1}^m\|^2 + L_C(\mu_n + \mu_m),
$$

where <sup>ξ</sup>*n*, ξ*<sup>m</sup>* <sup>∈</sup> <sup>B</sup>, - := sup{ λ(*t*) μ : *<sup>t</sup>* ∈ [0, *<sup>T</sup>* ], <sup>∈</sup> <sup>N</sup>} and <sup>ζ</sup>*<sup>i</sup>* := *<sup>f</sup>* (*t*, *xi*(δ*i*(*t*))) <sup>+</sup> 3 <sup>√</sup>ε*<sup>i</sup>* <sup>μ</sup>*<sup>i</sup>* ξ*<sup>i</sup>* for *i* ∈ {*n*, *m*}. Therefore, we have that

$$
\langle \dot{x}_n(t) - \dot{x}_m(t), v_{k+1}^n - v_{j+1}^m \rangle
$$
  
\n
$$
= \langle \dot{x}_n(t) - \zeta_n, v_{k+1}^n - v_{j+1}^m \rangle + \langle \zeta_n - \zeta_m, v_{k+1}^n - v_{j+1}^m \rangle
$$
  
\n
$$
+ \langle \zeta_m - \dot{x}_m(t), v_{k+1}^n - v_{j+1}^m \rangle
$$
  
\n
$$
\leq 2F \left( \frac{2}{\rho} ||v_{k+1}^n - v_{j+1}^m||^2 + L_C(\mu_n + \mu_m) \right) + \langle \zeta_n - \zeta_m, v_{k+1}^n - v_{j+1}^m \rangle
$$
  
\n
$$
\leq \frac{4F}{\rho} (||x_n(t) - x_m(t)|| + 3(\sqrt{\varepsilon_n} + \sqrt{\varepsilon_m}) + K_4(\mu_n + \mu_m))^2
$$
  
\n
$$
+ 2FL_C(\mu_n + \mu_m) + \langle \zeta_n - \zeta_m, v_{k+1}^n - v_{j+1}^m \rangle.
$$

Moreover, by virtue of Theorem 1, we have max{*xn*∞, *xm*∞} ≤ *K*2. Hence, there is *<sup>k</sup>* <sup>∈</sup> *<sup>L</sup>*1([0, *<sup>T</sup>* ]; <sup>R</sup>+) satisfying (12) on *<sup>K</sup>*2B. Therefore, it follows that

$$
\langle \zeta_n - \zeta_m, v_{k+1}^n - v_{j+1}^m \rangle
$$
  
=\langle f(t, x\_n(\delta\_n(t))) - f(t, x\_m(\delta\_m(t))), x\_n(\delta\_n(t)) - x\_m(\delta\_m(t))) +  
+ \langle f(t, x\_n(\delta\_n(t))) - f(t, x\_m(\delta\_m(t))), v\_{k+1}^n - x\_{k+1}^n \rangle  
+ \langle f(t, x\_n(\delta\_n(t))) - f(t, x\_m(\delta\_m(t))), x\_{k+1}^n - x\_k^n \rangle  
+ \langle f(t, x\_n(\delta\_n(t))) - f(t, x\_m(\delta\_m(t))), x\_j^m - x\_{j+1}^m \rangle  
+ \langle f(t, x\_n(\delta\_n(t))) - f(t, x\_m(\delta\_m(t))), x\_{j+1}^m - v\_{j+1}^m \rangle  
+ \frac{3\sqrt{\varepsilon\_n}}{\mu\_n} \langle \xi\_n, v\_{k+1}^n - v\_{j+1}^m \rangle + \frac{3\sqrt{\varepsilon\_m}}{\mu\_m} \langle \xi\_m, v\_{j+1}^m - v\_{k+1}^n \rangle  
= k(t) \|x\_n(\delta\_n(t)) - x\_m(\delta\_m(t)) \|^2  
+ 2\alpha (3(\sqrt{\varepsilon\_n} + \sqrt{\varepsilon\_m}) + K\_3(\mu\_n + \mu\_m))  
+ \frac{3\sqrt{\varepsilon\_n}}{\mu\_n} \|v\_{k+1}^n - v\_{j+1}^m\| + \frac{3\sqrt{\varepsilon\_m}}{\mu\_m} \|v\_{j+1}^m - v\_{k+1}^n\|  
\le k(t) (\|x\_n(t) - x\_m(t)\| + 3(\sqrt{\varepsilon\_n} + \sqrt{\varepsilon\_m}) + (K\_3 + K\_4)(\mu\_n + \mu\_m))^2  
+ 2\alpha (3(\sqrt{\varepsilon\_n} + \sqrt{\varepsilon\_m}) + K\_3(\mu\_n + \mu\_m))  
+ 6\left(\frac{\sqrt{\varepsilon\_n}}{\mu\_n} + \frac{\sqrt{\varepsilon\_m}}{\mu\_m}\right) (\sqrt{\varepsilon\_n} + \sqrt{\varepsilon\_m} + K\_2).

123

These two inequalities and (14) yield

$$
\frac{d}{dt} ||x_n(t) - x_m(t)||^2
$$
\n
$$
\leq 4\left(\frac{4F}{\rho} + k(t)\right) ||x_n(t) - x_m(t)||^2 + 4\alpha(3(\sqrt{\varepsilon_n} + \sqrt{\varepsilon_m}) + K_3(\mu_n + \mu_m))
$$
\n
$$
+ 4FL_C(\mu_n + \mu_m) + 12\left(\frac{\sqrt{\varepsilon_n}}{\mu_n} + \frac{\sqrt{\varepsilon_m}}{\mu_m}\right)(\sqrt{\varepsilon_n} + \sqrt{\varepsilon_m} + K_2)
$$
\n
$$
+ \frac{16F}{\rho}(3(\sqrt{\varepsilon_n} + \sqrt{\varepsilon_m}) + K_4(\mu_n + \mu_m))^2
$$
\n
$$
+ 4k(t)(3(\sqrt{\varepsilon_n} + \sqrt{\varepsilon_m}) + (K_3 + K_4)(\mu_n + \mu_m))^2.
$$

Hence, using Gronwall's inequality, we have for all *t* ∈ [0, *T* ] and *n*, *m* big enough:

$$
||x_n(t) - x_m(t)||^2 \le A_{m,n} \exp\left(\frac{16F}{\rho}T + 4\int_0^T k(s)ds\right),
$$
 (15)

where

$$
A_{m,n} = 4\alpha T (3(\sqrt{\varepsilon_n} + \sqrt{\varepsilon_m}) + K_3(\mu_n + \mu_m))
$$
  
+ 4TFL<sub>C</sub>( $\mu_n + \mu_m$ ) + 12T  $\left(\frac{\sqrt{\varepsilon_n}}{\mu_n} + \frac{\sqrt{\varepsilon_m}}{\mu_m}\right) (\sqrt{\varepsilon_n} + \sqrt{\varepsilon_m} + K_2)$   
+  $\frac{16TF}{\rho} (3(\sqrt{\varepsilon_n} + \sqrt{\varepsilon_m}) + K_4(\mu_n + \mu_m))^2$   
+ 4 $||k||_1 (3(\sqrt{\varepsilon_n} + \sqrt{\varepsilon_m}) + (K_3 + K_4)(\mu_n + \mu_m))^2$ .

Since *Am*,*<sup>n</sup>* goes to 0 when *m*, *n* → ∞, it shows that (*xn*) is a Cauchy sequence in the space of continuous functions with the uniform convergence. Therefore, it converges uniformly to some continuous function *x* : [0, *T* ] → *H*. It remains to check that *x* is absolutely continuous, and it is the unique solution of (4). First of all, by Theorem 1 and Lemma 4, *x* is absolutely continuous and there is a subsequence of (*x*˙*n*) which converges weakly in *<sup>L</sup>*1([0, *<sup>T</sup>* ]; *<sup>H</sup>*) to *<sup>x</sup>*˙. So, without relabeling, we have *<sup>x</sup>*˙*nx*˙ in *<sup>L</sup>*1([0, *<sup>T</sup>* ]; *<sup>H</sup>*). On the other hand, using Theorem <sup>1</sup> and defining <sup>v</sup>*n*(*t*) := <sup>v</sup>*<sup>n</sup> <sup>k</sup>*+<sup>1</sup> for *t* ∈]*t n <sup>k</sup>* , *t n <sup>k</sup>*+1] we have

$$
\dot{x}_n(t) \in -\frac{\lambda_n(t)}{\mu_n} \partial_P d_{C(\theta_n(t))}(v_n(t)) + f(t, x_n(\delta_n(t))) + \frac{3\sqrt{\varepsilon_n}}{\mu_n} \mathbb{B}
$$
  

$$
\in -\kappa_1 \partial d_{C(\theta_n(t))}(v_n(t)) + \kappa_2 \mathbb{B} \cap F(t, x_n(\delta_n(t))) + \frac{3\sqrt{\varepsilon_n}}{\mu_n} \mathbb{B},
$$

where, by Theorem 1, κ<sup>1</sup> and κ<sup>2</sup> are non-negative numbers which do not depend of *<sup>n</sup>* <sup>∈</sup> <sup>N</sup> and *<sup>t</sup>* ∈ [0, *<sup>T</sup>* ]. We also have <sup>v</sup>*<sup>n</sup>* <sup>→</sup> *<sup>x</sup>*, <sup>θ</sup>*<sup>n</sup>* <sup>→</sup> Id[0,*<sup>T</sup>* ] and <sup>δ</sup>*<sup>n</sup>* <sup>→</sup> Id[0,*<sup>T</sup>* ] uniformly. Theorem 1 ensures that *x*(*t*) ∈ *C*(*t*) for all *t* ∈ [0, *T* ]. By Mazur's lemma, there is a

### 123

sequence (*y <sup>j</sup>*) such that for all *n*, *yn* ∈ co(*x*˙*<sup>k</sup>* : *k* ≥ *n*) and (*yn*) converges strongly to *<sup>x</sup>*˙ in *<sup>L</sup>*1([0, *<sup>T</sup>* ]; *<sup>H</sup>*). That is to say

$$
y_n(t) \in \text{co}\left(-\kappa_1 \partial d_{C(\theta_k(t))}(v_k(t)) + \kappa_2 \mathbb{B} \cap F(t, x_k(\delta_k(t))) + \frac{3\sqrt{\varepsilon_k}}{\mu_k} \mathbb{B} : k \ge n\right).
$$

Hence, there exists (*yn <sup>j</sup>*) which converges to *x*˙ almost everywhere in [0, *T* ]. Then, by virtue of Lemma 2, (*H<sup>F</sup>* <sup>1</sup> ) and Lemma 3, we obtain that

$$
\dot{x}(t) \in -\kappa_1 \partial d_{C(t)}(x(t)) + \kappa_2 \mathbb{B} \cap F(t, x(t)) \text{ for a.e. } t \in [0, T].
$$

Since ∂*dC*(*<sup>t</sup>*)(*x*(*t*)) ⊂ *N*(*C*(*t*); *x*(*t*)) for all *t* ∈ [0, *T* ], we have *x* is the solution of (4).

To end the proof, we are going to prove that (4) has a unique solution under growth condition (13). First, take any solution *x* of (4). Then, for a.e. *t* ∈ [0, *T* ] there is *f* (*t*, *x*(*t*)) ∈ *F*(*t*, *x*(*t*)) such that

$$
\mathcal{R}_x(t) := f(t, x(t)) - \dot{x}(t) \in N(C(t); x(t)).
$$
\n(16)

Take any *t* ∈]0, *T* ] satisfying (16). Suppose that *x*˙(*t*) = *f* (*t*, *x*(*t*)), then using (1) and the uniform prox-regularity of *C*(*t*) we have that

$$
\frac{\mathcal{R}_x(t)}{\|\mathcal{R}_x(t)\|} \in \partial_P d_{C(t)}(x(t)).
$$

Take any <sup>γ</sup> ∈]0, <sup>1</sup>[, by continuity there is δ > 0 such that *<sup>x</sup>*(*s*) <sup>∈</sup> *<sup>U</sup>*<sup>γ</sup> <sup>ρ</sup> (*C*(*t*)) for all *s* ∈]*t* − δ, *t* + δ[, using Proposition 1 we have

$$
\left\langle \frac{\mathcal{R}_x(t)}{\|\mathcal{R}_x(t)\|}, x(s) - x(t) \right\rangle \le \frac{1}{2\rho(1-\gamma)^2} \|x(s) - x(t)\|^2 + d_{C(t)}(x(s))
$$
  

$$
\le \frac{1}{2\rho(1-\gamma)^2} \|x(s) - x(t)\|^2 + L_C|t - s|.
$$

Dividing by *t* − *s* for *s* ∈]*t* − δ, *t*[ and taking the limit *s t*, we obtain that

$$
\left\langle \frac{\mathcal{R}_x(t)}{\|\mathcal{R}_x(t)\|}, -\dot{x}(t) \right\rangle \leq L_C \implies \|\mathcal{R}_x(t)\| \leq \|f(t, x(t))\| + L_C.
$$

When *x*˙(*t*) = *f* (*t*, *x*(*t*)), the above inequality always holds. Hence, for a.e. *t* ∈ [0, *T* ].

Now, take two solutions *x*1, *x*<sup>2</sup> of (4) with *x*1(0) = *x*2(0) = *x*0, then using the hypomonotonicity given in Proposition 1, we have

#### 123

$$
\langle \mathcal{R}_{x_1}(t) - \mathcal{R}_{x_2}(t), x_1(t) - x_2(t) \rangle \geq \frac{-1}{2\rho} (\|\mathcal{R}_{x_1}(t)\| + \|\mathcal{R}_{x_2}(t)\|) \|x_1(t) - x_2(t)\|^2.
$$

Defining *<sup>r</sup>* <sup>=</sup> max{*xi*∞ : *<sup>i</sup>* <sup>=</sup> <sup>1</sup>, <sup>2</sup>}, there is *kr* <sup>∈</sup> *<sup>L</sup>*1([0, *<sup>T</sup>* ]; <sup>R</sup>+) satisfying (12) on *r*B. Hence, by using growth condition (13), we have a.e.

$$
\frac{d}{dt}(\|x_1(t) - x_2(t)\|^2) \le \|x_1(t) - x_2(t)\|^2 [2k_r(t) + \frac{1}{\rho}(\|\mathcal{R}_{x_1}(t)\| + \|\mathcal{R}_{x_2}(t)\|)]
$$
  
\n
$$
\le \|x_1(t) - x_2(t)\|^2 [2k_r(t) + \frac{2L_C}{\rho} + 2c(t) \left(\frac{1}{\rho} + c\right)],
$$

which, by virtue of Gronwall's inequality, implies that *x*<sup>1</sup> ≡ *x*2. The result is proven. 

*Remark 1* The property required for *F* in (12) is a classical monotonicity assumption in the theory of existence of solutions for differential inclusions (see, e.g., [15, Theorem 10.5]).

*Remark 2* [Rate of convergence] In the precedent proof, we have established the following estimation:

$$
||x_n(t) - x_m(t)||^2 \le A_{m,n} \exp\left(\frac{16F}{\rho}T + 4\int_0^T k(s)ds\right)
$$

for *<sup>m</sup>*, *<sup>n</sup>* such that <sup>μ</sup>*<sup>n</sup>* <sup>+</sup> <sup>μ</sup>*<sup>m</sup>* <sup>&</sup>lt; <sup>ρ</sup> <sup>2</sup>*LC* . Hence, by letting *m* → ∞, we obtain that

$$
||x_n(t) - x(t)||^2 \le A_n \exp\left(\frac{16F}{\rho}T + 4\int_0^T k(s)ds\right) \text{ for all } n > \frac{2L_C T}{\rho},
$$

where

$$
A_n := \lim_{m \to \infty} A_{m,n} \le D\left(\sqrt{\varepsilon_n} + \mu_n + \frac{\sqrt{\varepsilon_n}}{\mu_n}\right),
$$

where *D* is a non-negative constant. Hence, the above estimation provides a rate of convergence for our scheme.

#### **5 Subsmooth Case**

In this section, we study sweeping process (4) for the class of subsmooth sets, which strictly includes the class of uniformly prox-regular sets. We now assume (*C*(*t*))*t*∈[0,*<sup>T</sup>* ] is a equi-uniformly subsmooth family. The classical catching-up algorithm was studied in [20] under this framework. In this case, we assume the ball compactness of the moving sets, required in the infinite-dimensional setting. We will see that our algorithm allows us to prove the existence of a solution, but we only ensure that a subsequence

#### 123

converges to this solution, which is expected due to the lack of uniqueness of solutions in this case.

**Theorem 3** *Suppose, in addition to assumptions of theorem* 1*, that the family* (*C*(*t*))*t*∈[0,*<sup>T</sup>* ] *is equi-uniformly subsmooth and the set C*(*t*) *are ball-compact for all t* ∈ [0, *T* ]*. Then, the sequence of continuous functions* (*xn*) *generated by algorithm* (5) *and* (6) *converges uniformly (up to a subsequence) to an absolutely continuous function x, which is a solution of* (4)*.*

*Proof* From Theorem <sup>1</sup> we have for all *<sup>n</sup>* <sup>∈</sup> <sup>N</sup> and *<sup>k</sup>* ∈ {0,..., *<sup>n</sup>* <sup>−</sup> <sup>1</sup>}, there is v*n <sup>k</sup>*+<sup>1</sup> <sup>∈</sup> *<sup>C</sup>*(*<sup>t</sup> n <sup>k</sup>*+1) such that v*<sup>n</sup> <sup>k</sup>*+<sup>1</sup> <sup>−</sup> *<sup>x</sup><sup>n</sup> <sup>k</sup>*+1 <sup>&</sup>lt; <sup>2</sup> <sup>√</sup>ε*<sup>n</sup>* and for all *<sup>t</sup>* ∈]*<sup>t</sup> n <sup>k</sup>* , *t n <sup>k</sup>*+1]:

$$
\dot{x}_n(t) \in -\frac{\lambda_n(t)}{\mu_n} \partial_P d_{C(\theta_n(t))}(v_{k+1}^n) + f(t, x_n(\delta_n(t))) + \frac{3\sqrt{\varepsilon_n}}{\mu_n} \mathbb{B},
$$

where λ*n*(*t*) = 4 <sup>√</sup>ε*<sup>n</sup>* <sup>+</sup> (*LC* <sup>+</sup> *<sup>h</sup>*(*x*(δ*n*(*t*))) <sup>+</sup> <sup>√</sup>γ )μ*n*. As *<sup>h</sup>* is *Lh*-Lipschitz it follows that

$$
\lambda_n(t) \le (4\mathfrak{c} + L_C + h(x_0) + \sqrt{\gamma} + L_h K_1)\mu_n.
$$

Defining <sup>v</sup>*n*(*t*) := <sup>v</sup>*<sup>n</sup> <sup>k</sup>*+<sup>1</sup> on ]*<sup>t</sup> n <sup>k</sup>* , *t n <sup>k</sup>*+1], then for all *<sup>n</sup>* <sup>∈</sup> <sup>N</sup> and almost all *<sup>t</sup>* ∈ [0, *<sup>T</sup>* ]

$$
\dot{x}_n(t) \in -M \partial_P d_{C(\theta_n(t))}(v_n(t)) + f(t, x_n(\delta_n(t))) + \frac{3\sqrt{\varepsilon_n}}{\mu_n} \mathbb{B}
$$
\n
$$
\in -M \partial d_{C(\theta_n(t))}(v_n(t)) + M \mathbb{B} \cap F(t, x_n(\delta_n(t))) + \frac{3\sqrt{\varepsilon_n}}{\mu_n} \mathbb{B},
$$
\n(17)

where *<sup>M</sup>* := <sup>4</sup><sup>c</sup> <sup>+</sup> *LC* <sup>+</sup> *<sup>h</sup>*(*x*0) <sup>+</sup> *LhK*<sup>1</sup> <sup>+</sup> <sup>√</sup><sup>γ</sup> . Moreover, by Theorem 1, we have

$$
d_{C(t)}(x_n(t)) \le d_{C(\theta_n(t))}(x_n(t)) + L_C \mu_n \le (K_5 + 2L_C)\mu_n + 2\sqrt{\varepsilon_n}.
$$
 (18)

for all *t* ∈ [0, *T* ].

Next, fix *<sup>t</sup>* ∈ [0, *<sup>T</sup>* ] and define *<sup>K</sup>*(*t*) := {*xn*(*t*) : *<sup>n</sup>* <sup>∈</sup> <sup>N</sup>}. We claim that *<sup>K</sup>*(*t*) is relatively compact. Indeed, let *xm*(*t*) ∈ *K*(*t*) and take *ym*(*t*) ∈ Proj*C*(*<sup>t</sup>*)(*xm*(*t*)) (the projection exists due to the ball compactness of *C*(*t*) and the boundedness of *K*(*t*)). Moreover, according to (18) and Theorem 1,

$$
||y_n(t)|| \leq d_{C(t)}(x_n(t)) + ||x_n(t)|| \leq (K_5 + 2L_C)\mu_n + 2\sqrt{\varepsilon_n} + K_2.
$$

This entails that *yn*(*t*) <sup>∈</sup> *<sup>C</sup>*(*t*) <sup>∩</sup> *<sup>R</sup>* <sup>B</sup> for all *<sup>n</sup>* <sup>∈</sup> <sup>N</sup> for some *<sup>R</sup>* <sup>&</sup>gt; 0. Thus, by the ball compactness of *C*(*t*), there exists a subsequence (*ymk* (*t*)) of (*ym*(*t*)) converging to some *y*(*t*) as *k* → +∞. Then,

#### 123

$$
||x_{m_k}(t) - y(t)|| \leq d_{C(t)}(x_{m_k}(t)) + ||y_{m_k}(t) - y(t)||
$$
  
\n
$$
\leq (K_5 + 2L_C)\mu_{m_k} + 2\sqrt{\varepsilon_{m_k}} + ||y_{m_k}(t) - y(t)||,
$$

which implies that *K*(*t*) is relatively compact. Moreover, it is not difficult to see by Theorem 1 that *K* := (*xn*) is equicontinuous. Therefore, by virtue of Theorem 1, Arzela-Ascoli's and Lemma 4, we obtain the existence of a Lipschitz function *x* and a subsequence (*x <sup>j</sup>*) of (*xn*) such that

- (i) (*x <sup>j</sup>*) converges uniformly to *x* on [0, *T* ].
- (ii) *<sup>x</sup>*˙*jx*˙ in *<sup>L</sup>*<sup>1</sup> ([0, *<sup>T</sup>* ]; *<sup>H</sup>*).
- (iii) *x <sup>j</sup>*(θ *<sup>j</sup>*(*t*)) → *x*(*t*) for all *t* ∈ [0, *T* ].
- (iv) *x <sup>j</sup>*(δ *<sup>j</sup>*(*t*)) → *x*(*t*) for all *t* ∈ [0, *T* ].
- (v) v *<sup>j</sup>*(*t*) → *x*(*t*) for all *t* ∈ [0, *T* ].

From (18) it is clear that *x*(*t*) ∈ *C*(*t*) for all *t* ∈ [0, *T* ]. By Mazur's lemma, there is a sequence (*y <sup>j</sup>*) such that for all *j*, *y <sup>j</sup>* ∈ co(*x*˙*<sup>k</sup>* : *k* ≥ *j*) and (*y <sup>j</sup>*) converges strongly to *<sup>x</sup>*˙ in *<sup>L</sup>*1([0, *<sup>T</sup>* ]; *<sup>H</sup>*). That is to say

$$
y_j(t) \in \text{co}\left(-M\partial d_{C(\theta_n(t))}(v_n(t)) + M\mathbb{B} \cap F(t, x_n(\delta_n(t))) + \frac{3\sqrt{\varepsilon_n}}{\mu_n} \mathbb{B} : n \ge j\right).
$$

On the other hand, there exists(*yn <sup>j</sup>*) which converges to *x*˙ almost everywhere in [0, *T* ]. Then, using Lemma 2, Lemma 3, and (*H<sup>F</sup>* <sup>1</sup> ), we have

$$
\dot{x}(t) \in -M \partial d_{C(t)}(x(t)) + M \mathbb{B} \cap F(t, x(t))
$$
 a.e.

Finally, since ∂*dC*(*<sup>t</sup>*)(*x*(*t*)) ⊂ *N*(*C*(*t*); *x*(*t*)) for all *t* ∈ [0, *T* ], it follows that *x* is the solution of (4).

#### **6 Fixed Set**

In this section, we consider a closed and nonempty set *C* ⊂ *H*, and we look for a solution of the particular case of (4) given by

$$
\dot{x}(t) \in -N(C; x(t)) + F(t, x(t)) \quad \text{a.e. } t \in [0, T],
$$
  
\n
$$
x(0) = x_0 \in C,
$$
\n(19)

where *F* : [0, *T* ] × *H* ⇒ *H* is a set-valued map defined as above. The existence of a solution using classical catching up was done in [34]. Now, we use similar ideas to get the existence of a solution using our proposed algorithm. We emphasize that in this case, no regularity of the set *C* is required.

**Theorem 4** *Let C* ⊂ *H be a ball-compact set and F* : [0, *T* ] × *H* ⇒ *H be a setvalued map satisfying* (*H<sup>F</sup>* <sup>1</sup> )*,* (*H<sup>F</sup>* <sup>2</sup> ) *and* (*H<sup>F</sup>* <sup>3</sup> )*. Then, for any x*<sup>0</sup> ∈ *S, the sequence of functions* (*xn*) *generated by algorithm* (6) *converges uniformly (up to a subsequence)*

## 123

*to a Lipschitz solution x of sweeping process* (19) *such that*

$$
\|\dot{x}(t)\| \le 2(h(x(t)) + \sqrt{\gamma}) \quad a.e. \ t \in [0, T]. \tag{20}
$$

*Proof* We are going to use the properties of Theorem 1, where now we have *LC* <sup>=</sup> 0. First of all, from Theorem <sup>1</sup> we have for all *<sup>n</sup>* <sup>∈</sup> <sup>N</sup> and *<sup>k</sup>* ∈ {0, <sup>1</sup>,..., *<sup>n</sup>* <sup>−</sup> <sup>1</sup>}, there is v*n <sup>k</sup>*+<sup>1</sup> <sup>∈</sup> *<sup>C</sup>* such that v*<sup>n</sup> <sup>k</sup>*+<sup>1</sup> <sup>−</sup> *<sup>x</sup><sup>n</sup> <sup>k</sup>*+1 <sup>&</sup>lt; <sup>2</sup> <sup>√</sup>ε*<sup>n</sup>* and for all *<sup>t</sup>* ∈]*<sup>t</sup> n <sup>k</sup>* , *t n <sup>k</sup>*+1]:

$$
\dot{x}_n(t) \in -\frac{\lambda_n(t)}{\mu_n} \partial_P d_C(v_{k+1}^n) + f(t, x_n(\delta_n(t))) + \frac{3\sqrt{\varepsilon_n}}{\mu_n} \mathbb{B},
$$

where λ*n*(*t*) = 4 <sup>√</sup>ε*<sup>n</sup>* <sup>+</sup> (*h*(*x*(δ*n*(*t*))) <sup>+</sup> <sup>√</sup>γ )μ*n*. Defining <sup>v</sup>*n*(*t*) := <sup>v</sup>*<sup>n</sup> <sup>k</sup>*+<sup>1</sup> on ]*<sup>t</sup> n <sup>k</sup>* , *t n <sup>k</sup>*+1], we get that for all *<sup>n</sup>* <sup>∈</sup> <sup>N</sup> and a.e. *<sup>t</sup>* ∈ [0, *<sup>T</sup>* ]

$$
\dot{x}_n(t) \in -\frac{\lambda_n(t)}{\mu_n} \partial_P d_C(v_n(t)) + f(t, x_n(\delta_n(t))) + \frac{3\sqrt{\varepsilon_n}}{\mu_n} \mathbb{B}
$$
\n
$$
\in -\frac{\lambda_n(t)}{\mu_n} \partial d_C(v_n(t)) + (h(t, x_n(\delta_n(t))) + \sqrt{\gamma}) \mathbb{B} \cap F(t, x_n(\delta_n(t))) + \frac{3\sqrt{\varepsilon_n}}{\mu_n} \mathbb{B}.
$$

Moreover, by Theorem 1, we have

$$
d_C(x_n(t)) \le K_5\mu_n + 2\sqrt{\varepsilon_n} \text{ for all } t \in [0, T].
$$

Next, fix *<sup>t</sup>* ∈ [0, *<sup>T</sup>* ] and define *<sup>K</sup>*(*t*) := {*xn*(*t*) : *<sup>n</sup>* <sup>∈</sup> <sup>N</sup>}. We claim that *<sup>K</sup>*(*t*) is relatively compact. Indeed, let *xm*(*t*) ∈ *K*(*t*) and take *ym*(*t*) ∈ Proj*C*(*xm*(*t*)) (the projection exists due to the ball compactness of *C* and the boundedness of *K*(*t*)). Moreover, according to the above inequality and Theorem 1,

$$
||y_n(t)|| \leq d_C(x_n(t)) + ||x_n(t)|| \leq K_5\mu_n + 2\sqrt{\varepsilon_n} + K_2,
$$

which entails that *yn*(*t*) <sup>∈</sup> *<sup>C</sup>* <sup>∩</sup> *<sup>R</sup>* <sup>B</sup> for all *<sup>n</sup>* <sup>∈</sup> <sup>N</sup> for some *<sup>R</sup>* <sup>&</sup>gt; 0. Thus, by the ball-compactness of *C*, there exists a subsequence (*ymk* (*t*)) of (*ym*(*t*)) converging to some *y*(*t*) as *k* → +∞. Then,

$$
||x_{m_k}(t) - y(t)|| \leq d_C(x_{m_k}(t)) + ||y_{m_k}(t) - y(t)||
$$
  
\n
$$
\leq K_5 \mu_{m_k} + 2\sqrt{\varepsilon_{m_k}} + ||y_{m_k}(t) - y(t)||,
$$

which implies that *K*(*t*) is relatively compact. Moreover, it is not difficult to see by Theorem 1 that the set *K* := (*xn*) is equicontinuous. Therefore, by virtue of Theorem 1, Arzela-Ascoli's and Lemma 4, we obtain the existence of a Lipschitz function *x* and a subsequence (*x <sup>j</sup>*) of (*xn*) such that

- (i) (*x <sup>j</sup>*) converges uniformly to *x* on [0, *T* ].
- (ii) *<sup>x</sup>*˙*jx*˙ in *<sup>L</sup>*<sup>1</sup> ([0, *<sup>T</sup>* ]; *<sup>H</sup>*).
- (iii) *x <sup>j</sup>*(θ *<sup>j</sup>*(*t*)) → *x*(*t*) for all *t* ∈ [0, *T* ].
- (iv) *x <sup>j</sup>*(δ *<sup>j</sup>*(*t*)) → *x*(*t*) for all *t* ∈ [0, *T* ].

### 123

(v) v *<sup>j</sup>*(*t*) → *x*(*t*) for all *t* ∈ [0, *T* ]. (vi) *x*(*t*) ∈ *C* for all *t* ∈ [0, *T* ].

By Mazur's lemma, there is a sequence (*y <sup>j</sup>*) such that for all *j*, *y <sup>j</sup>* ∈ co(*x*˙*<sup>k</sup>* : *k* ≥ *j*) and (*<sup>y</sup> <sup>j</sup>*) converges strongly to *<sup>x</sup>*˙ in *<sup>L</sup>*1([0, *<sup>T</sup>* ]; *<sup>H</sup>*). i.e.,

$$
y_j(t) \in \text{co}\left(-\alpha_n \partial d_C(v_n(t)) + \beta_n \mathbb{B} \cap F(t, x_n(\delta_n(t))) + \frac{3\sqrt{\varepsilon_n}}{\mu_n} \mathbb{B} : n \ge j\right),\,
$$

where α*<sup>n</sup>* := 4 <sup>√</sup>ε*<sup>n</sup>* <sup>μ</sup>*<sup>n</sup>* <sup>+</sup> *<sup>h</sup>*(*t*, *xn*(δ*n*(*t*))) <sup>+</sup> <sup>√</sup><sup>γ</sup> and <sup>β</sup>*<sup>n</sup>* := 4 <sup>√</sup>ε*<sup>n</sup>* <sup>μ</sup>*<sup>n</sup>* + *h*(*t*, *xn*(δ*n*(*t*))). On the other hand, there exists (*yn <sup>j</sup>*) which converges to *x*˙ almost everywhere in [0, *T* ]. Then, using Lemma 2, Lemma 3, and (*H<sup>F</sup>* <sup>1</sup> ), we have

$$
\dot{x}(t) \in -(h(x(t)) + \sqrt{\gamma})\partial d_C(x(t)) + (h(x(t)) + \sqrt{\gamma})\mathbb{B} \cap F(t, x(t)) \text{ for a.e. } t \in [0, T].
$$

It is clear that *x* satisfies bound (20). Finally, since ∂*dC*(*x*(*t*)) ⊂ *N*(*C*; *x*(*t*)) for all *t* ∈ [0, *T* ], we obtain that *x* is the solution of (19).

#### **7 Numerical Methods for Approximate Projections**

As stated before, in most cases, finding an explicit formula for the projection onto a closed set is not possible. Therefore, one must resort to numerical algorithms to obtain approximate projections. Several papers discuss this issue for different notions of approximate projections (see, e.g., [32]). These algorithms are called *projection oracles* and provide an approximate solution *z*¯ ∈ *H* to the following optimization problem:

$$
\min_{z \in C} \|x - z\|^2, \tag{P_x}
$$

where *C* is a given closed set and *x* ∈ *H*. Whether the approximate solution *z*¯ belongs to the set *C* or not depends on the notion of approximate projection. In our case, to implement our algorithm, we need that *z*¯ ∈ *C*. In this line, a well-known projection oracle fulfilling this property can be obtained via the celebrated Frank–Wolfe algorithm (see, e.g., [18, 22]), where a linear sub-problem of (*Px* ) is solved in each iteration. For several types of convex sets, this method has been successfully developed (see [5, 13, 22]). Besides, in [16], it was shown that an approximate solution of the linear sub-problem is enough to obtain a projection oracle.

Another important approach to obtaining approximate projections is the use of the Frank–Wolfe algorithm with separation oracles (see [14]). Roughly speaking, a separation oracle determines whether a given point belongs to a set and, in the negative case, provides a hyperplane separating the point from the set (see [19] for more details). For particular sets, it is easy to get an explicit separation oracle (see [19, p. 49]). An important example is the case of a sublevel set: let *<sup>g</sup>* : *<sup>H</sup>* <sup>→</sup> <sup>R</sup> be a continuous convex function and <sup>λ</sup> <sup>∈</sup> <sup>R</sup>. Then [*<sup>g</sup>* <sup>≤</sup> <sup>λ</sup>] := {*<sup>x</sup>* <sup>∈</sup> *<sup>H</sup>* : *<sup>g</sup>*(*x*) <sup>≤</sup> <sup>λ</sup>} has a separation oracle described as follows: to verify that any point belongs to [*g* ≤ λ] is straightforward.

### 123

When a point *x* ∈ *H* does not belong to [*g* ≤ λ], we can consider any *x*<sup>∗</sup> ∈ ∂*g*(*x*). Then, for all *y* ∈ [*g* ≤ λ],

$$
\langle x^*, x \rangle \ge g(x) - g(y) + \langle x^*, y \rangle > \langle x^*, y \rangle,
$$

where we have used that *g*(*x*)>λ ≥ *g*(*y*). Hence, the above inequality shows the existence of the desired hyperplane, which provides a separation oracle for [*g* ≤ λ]. Therefore, if *C* is the sublevel set of some convex function, we can use the algorithm proposed in [14] to get an approximate solution *<sup>z</sup>*¯ <sup>∈</sup> proj *<sup>S</sup>*(*x*). Moreover, the sublevel set enables us to consider the case

$$
C(t,x) := \bigcap_{i=1}^{m} \{x \in \mathcal{H} : g_i(t,x) \leq 0\} = \left\{x \in \mathcal{H} : g(t,x) := \max_{i=1,\dots,m} g_i(t,x) \leq 0\right\},\
$$

where for all *<sup>t</sup>* ∈ [0, *<sup>T</sup>* ], *gi*(*t*, ·): *<sup>H</sup>* <sup>→</sup> <sup>R</sup>, *<sup>i</sup>* <sup>=</sup> <sup>1</sup>,..., *<sup>m</sup>* are convex functions. We refer to [2, Proposition 5.1] for the proper assumptions on these functions to ensure the Lipschitz property of the map *t* ⇒ *C*(*t*) holds (7).

#### **8 Concluding Remarks**

In this paper, we have developed an enhanced version of the catching-up algorithm for sweeping processes through an appropriate concept of approximate projections. We provide the proposed algorithm's convergence for three frameworks: prox-regular, subsmooth, and merely closed sets. Some insights into numerical procedures to obtain approximate projections were given mainly in the convex case. Finally, the convergence of our algorithm for other notions of approximate solutions will be explored in forthcoming works.

**Acknowledgements** This work was supported by the Center for Mathematical Modeling (CMM) and ANID-Chile under BASAL funds for Center of Excellence FB210005, Fondecyt Regular 1200283, Fondecyt Regular 1240120 and Fondecyt Exploración 13220097.We also acknowledge ECOS-Anid Project PC23E11 and MathAmsud Project 23-MATH-17. The authors would like to express their sincere gratitude to the anonymous referees for their valuable discussions and insightful suggestions. Their contributions have played a vital role in enhancing the overall quality of this manuscript.

#### **References**

- 1. Acary, V., Bonnefon, O., Brogliato, B.: Nonsmooth Modeling and Simulation for Switched Circuits. Springer, Berlin (2011)
- 2. Adly, S., Haddad, T.: Well-posedness of nonconvex degenerate sweeping process via unconstrained evolution problems. Nonlinear Anal. Hybrid Syst. **36**, 100832 (2020)
- 3. Aliprantis, C.D., Border, K.C.: Infinite Dimensional Analysis, 3rd edn. Springer, Berlin (2006)
- 4. Aussel, D., Daniilidis, A., Thibault, L.: Subsmooth sets: functional characterizations and related concepts. Trans. Am. Math. Soc. **357**(4), 1275–1301 (2005)
- 5. Bomze, I.M., Rinaldi, F., Zeffiro, D.: Frank–Wolfe and friends: a journey into projection-free first-order optimization methods. 4OR **19**(3), 313–345 (2021)

123

- 6. Borwein, J.M., Preiss, D.: A smooth variational principle with applications to subdifferentiability and to differentiability of convex functions. Trans. Am. Math. Soc. **303**(2), 517–527 (1987)
- 7. Bounkhel, M., Thibault, L.: On various notions of regularity of sets in nonsmooth analysis. Nonlinear Anal. **48**(2), 223–246 (2002)
- 8. Bounkhel, M., Thibault, L.: Nonconvex sweeping process and prox-regularity in Hilbert space. J. Nonlinear Convex Anal. **6**(2), 359–374 (2005)
- 9. Brogliato, B.: Nonsmooth Mechanics (Communications in Numerical Methods in Engineering), 3rd edn. Springer, Berlin (2016)
- 10. Clarke, F.: Optimization and Nonsmooth Analysis. Wiley Intersciences, New York (1983)
- 11. Clarke, F., Ledyaev, Y., Stern, R., Wolenski, P.: Nonsmooth Analysis and Control Theory (Graduate Texts in Mathematics), vol. 178. Springer, New York (1998)
- 12. Colombo, G., Thibault, L.: Prox-regular sets and applications. In: Gao, D.Y., Motreanu, D. (eds.) Handbook of Nonconvex Analysis and Applications, pp. 99–182. Int. Press, Somerville (2010)
- 13. Combettes, C.W., Pokutta, S.: Complexity of linear minimization and projection on some sets. Oper. Res. Lett. **49**(4), 565–571 (2021)
- 14. Dadush, D., Hojny, C., Huiberts, S., Weltge, S.: A simple method for convex optimization in the oracle model. In: International Conference on Integer Programming and Combinatorial Optimization, pp. 154–167. Springer (2022)
- 15. Deimling, K.: Multivalued Differential Equations, De Gruyter Series in Nonlinear Analysis and Applications, vol. 1. Walter de Gruyter & Co., Berlin (1992)
- 16. Ding, L., Udell, M.: Frank–Wolfe style algorithms for large scale optimization. In: Large-Scale and Distributed Optimization. Lecture Notes in Mathematics, vol. 2227, pp. 215–245. Springer, Cham (2018)
- 17. Federer, H.: Curvature measures. Trans. Am. Math. Soc. **93**, 418–491 (1959)
- 18. Frank, M., Wolfe, P., et al.: An algorithm for quadratic programming. Naval Res. Logist. Q. **3**(1–2), 95–110 (1956)
- 19. Grötschel, M., Lovász, L., Schrijver, A.: Geometric Algorithms and Combinatorial Optimization, Algorithms and Combinatorics: Study and Research Texts, vol. 2. Springer, Berlin (1988)
- 20. Haddad, T., Noel, J., Thibault, L.: Perturbed sweeping process with a subsmooth set depending on the state. Linear Nonlinear Anal. **2**(1), 155–174 (2016)
- 21. Hiriart-Urruty, J.B., López, M.A., Volle, M.: The -strategy in variational analysis: illustration with the closed convexification of a function. Rev. Mat. Iberoam. **27**(2), 449–474 (2011)
- 22. Jaggi, M.: Revisiting Frank–Wolfe: projection-free sparse convex optimization. In: Dasgupta, S., McAllester, D. (eds.) Proceedings of the 30th International Conference on Machine Learning, Proceedings of Machine Learning Research, vol. 28, no. 1, pp. 427–435. PMLR, Atlanta, Georgia, USA (2013)
- 23. Jourani, A., Vilches, E.: Moreau–Yosida regularization of state-dependent sweeping processes with nonregular sets. J. Optim. Theory Appl. **173**(1), 91–116 (2017)
- 24. Maury, B., Venel, J.: Un modéle de mouvement de foule. ESAIM Proc. **18**, 143–152 (2007)
- 25. Moreau, J.J.: Rafle par un convexe variable I, expo. 15. Sém, Anal. Conv. Mont., pp. 1–43 (1971)
- 26. Moreau, J.J.: Rafle par un convexe variable II, expo. 3. Sém, Anal. Conv. Mont., pp. 1–36 (1972)
- 27. Nacry, F., Thibault, L.: Distance function associated to a prox-regular set. Set-Valued Var. Anal. **30**(2), 731–750 (2022)
- 28. Papageorgiou, N.S., Kyritsi-Yiallourou, S.T.: Handbook of Applied Analysis (Advances in Mechanics and Mathematics), vol. 19. Springer, New York (2009)
- 29. Penot, J.P.: Calculus Without Derivatives (Graduate Texts in Mathematics), vol. 266. Springer, New York (2013)
- 30. Poliquin, R.A., Rockafellar, R.T., Thibault, L.: Local differentiability of distance functions. Trans. Am. Math. Soc. **352**(11), 5231–5249 (2000)
- 31. Thibault, L.: Unilateral Variational Analysis in Banach Spaces. Part II: Special Classes of Functions and Sets. World Scientific, Singapore (2023)
- 32. Usmanova, I., Kamgarpour, M., Krause, A., Levy, K.: Fast projection onto convex smooth constraints. In: Meila, M., Zhang, T. (eds.) Proceedings of the 38th International Conference on Machine Learning, Proceedings of Machine Learning Research, vol. 139, pp. 10,476–10,486. PMLR (2021)
- 33. Venel, J.: A numerical scheme for a class of sweeping processes. Numer. Math. **118**(2), 367–400 (2011)
- 34. Vilches, E.: Existence and Lyapunov pairs for the perturbed sweeping process governed by a fixed set. Set-Valued Var. Anal. **27**(2), 569–583 (2019)

# 123

**Publisher's Note** Springer Nature remains neutral with regard to jurisdictional claims in published maps and institutional affiliations.

Springer Nature or its licensor (e.g. a society or other partner) holds exclusive rights to this article under a publishing agreement with the author(s) or other rightsholder(s); author self-archiving of the accepted manuscript version of this article is solely governed by the terms of such publishing agreement and applicable law.

123

## Terms and Conditions

 Springer Nature journal content, brought to you courtesy of Springer Nature Customer Service Center GmbH ("Springer Nature").

Springer Nature supports a reasonable amount of sharing of research papers by authors, subscribers and authorised users ("Users"), for small-scale personal, non-commercial use provided that all copyright, trade and service marks and other proprietary notices are maintained. By accessing, sharing, receiving or otherwise using the Springer Nature journal content you agree to these terms of use ("Terms"). For these purposes, Springer Nature considers academic use (by researchers and students) to be non-commercial.

These Terms are supplementary and will apply in addition to any applicable website terms and conditions, a relevant site licence or a personal subscription. These Terms will prevail over any conflict or ambiguity with regards to the relevant terms, a site licence or a personal subscription (to the extent of the conflict or ambiguity only). For Creative Commons-licensed articles, the terms of the Creative Commons license used will apply.

We collect and use personal data to provide access to the Springer Nature journal content. We may also use these personal data internally within ResearchGate and Springer Nature and as agreed share it, in an anonymised way, for purposes of tracking, analysis and reporting. We will not otherwise disclose your personal data outside the ResearchGate or the Springer Nature group of companies unless we have your permission as detailed in the Privacy Policy.

While Users may use the Springer Nature journal content for small scale, personal non-commercial use, it is important to note that Users may not:

- 1. use such content for the purpose of providing other users with access on a regular or large scale basis or as a means to circumvent access control;
- 2. use such content where to do so would be considered a criminal or statutory offence in any jurisdiction, or gives rise to civil liability, or is otherwise unlawful;
- 3. falsely or misleadingly imply or suggest endorsement, approval , sponsorship, or association unless explicitly agreed to by Springer Nature in writing;
- 4. use bots or other automated methods to access the content or redirect messages
- 5. override any security feature or exclusionary protocol; or
- 6. share the content in order to create substitute for Springer Nature products or services or a systematic database of Springer Nature journal content.

In line with the restriction against commercial use, Springer Nature does not permit the creation of a product or service that creates revenue, royalties, rent or income from our content or its inclusion as part of a paid for service or for other commercial gain. Springer Nature journal content cannot be used for inter-library loans and librarians may not upload Springer Nature journal content on a large scale into their, or any other, institutional repository.

These terms of use are reviewed regularly and may be amended at any time. Springer Nature is not obligated to publish any information or content on this website and may remove it or features or functionality at our sole discretion, at any time with or without notice. Springer Nature may revoke this licence to you at any time and remove access to any copies of the Springer Nature journal content which have been saved.

To the fullest extent permitted by law, Springer Nature makes no warranties, representations or guarantees to Users, either express or implied with respect to the Springer nature journal content and all parties disclaim and waive any implied warranties or warranties imposed by law, including merchantability or fitness for any particular purpose.

Please note that these rights do not automatically extend to content, data or other material published by Springer Nature that may be licensed from third parties.

If you would like to use or distribute our Springer Nature journal content to a wider audience or on a regular basis or in any other manner not expressly permitted by these Terms, please contact Springer Nature at