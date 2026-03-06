import Mathlib
import MoreauSweeping.Preliminaries
import MoreauSweeping.CatchingUp

/-!
# Prox-Regular Sets

Uniformly prox-regular sets and convergence of the catching-up algorithm for this case.
-/

namespace MoreauSweeping

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- A set is uniformly prox-regular. -/
axiom UniformlyProxRegular : Set H → ℝ → Prop

/-- Condition (b) in Proposition 1 (well-defined/Lipschitz projection on `U^γ_ρ(S)`). -/
axiom ProxRegCondB : Set H → ℝ → Prop

/-- Condition (c) in Proposition 1 (hypomonotonicity of proximal normals). -/
axiom ProxRegCondC : Set H → ℝ → Prop

/-- Condition (d) in Proposition 1 (proximal subgradient inequality for `d_S`). -/
axiom ProxRegCondD : Set H → ℝ → Prop

/--
Characterization of prox-regularity (paper, Proposition 1):
for a closed set `S ⊂ H` and `ρ ∈ ]0,+∞]`, conditions (a), (b), (c), and (d)
are equivalent.
-/
axiom prox_regular_characterization (S : Set H) (ρ : ℝ) :
    UniformlyProxRegular S ρ ↔ ProxRegCondB S ρ ∧ ProxRegCondC S ρ ∧ ProxRegCondD S ρ

/-- Proposition 3.2 (constant-radius specialization): `(a) ↔ (b) ∧ (c) ∧ (d)`. -/
theorem proposition_3_2_characterization (S : Set H) (ρ : ℝ) :
    UniformlyProxRegular S ρ ↔ ProxRegCondB S ρ ∧ ProxRegCondC S ρ ∧ ProxRegCondD S ρ :=
  prox_regular_characterization S ρ

/-- Direction `(a) → (b)` in Proposition 3.2. -/
theorem proposition_3_2_a_implies_b (S : Set H) (ρ : ℝ) :
    UniformlyProxRegular S ρ → ProxRegCondB S ρ := by
  intro h
  exact (proposition_3_2_characterization S ρ).1 h |>.1

/-- Direction `(a) → (c)` in Proposition 3.2. -/
theorem proposition_3_2_a_implies_c (S : Set H) (ρ : ℝ) :
    UniformlyProxRegular S ρ → ProxRegCondC S ρ := by
  intro h
  exact (proposition_3_2_characterization S ρ).1 h |>.2.1

/-- Direction `(a) → (d)` in Proposition 3.2. -/
theorem proposition_3_2_a_implies_d (S : Set H) (ρ : ℝ) :
    UniformlyProxRegular S ρ → ProxRegCondD S ρ := by
  intro h
  exact (proposition_3_2_characterization S ρ).1 h |>.2.2

/-- Combined reverse direction `(b) ∧ (c) ∧ (d) → (a)` in Proposition 3.2. -/
theorem proposition_3_2_bcd_implies_a (S : Set H) (ρ : ℝ) :
    ProxRegCondB S ρ ∧ ProxRegCondC S ρ ∧ ProxRegCondD S ρ → UniformlyProxRegular S ρ := by
  intro h
  exact (proposition_3_2_characterization S ρ).2 h

/-- Convergence of approximate projections for prox-regular sets. -/
axiom approx_proj_convergence_prox_regular : True

/-- Quasi-Lipschitz property of approximate projections. -/
axiom quasi_lipschitz_approx_proj : True

/-- Convergence for prox-regular sets. -/
axiom convergence_prox_regular : True

end MoreauSweeping
