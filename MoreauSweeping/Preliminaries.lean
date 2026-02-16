import Mathlib

/-!
# Preliminaries

Basic definitions for the catching-up algorithm for Moreau's sweeping processes.
-/

namespace MoreauSweeping

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- `approxProj S x ε` is the set of `ε`-approximate projections of `x` on `S`. -/
def approxProj (S : Set H) (x : H) (ε : ℝ) : Set H :=
  {z : H | z ∈ S ∧ ‖x - z‖ ^ 2 < (infDist x S) ^ 2 + ε}

/-- A point of a set has zero distance to the set. -/
theorem proximal_normal_characterization (S : Set H) (x : H) (hx : x ∈ S) :
    infDist x S = 0 := by
  simpa using infDist_eq_zero_of_mem hx

/-- Monotonicity in the approximation error. -/
theorem approxProj_mono (S : Set H) (x : H) {ε₁ ε₂ : ℝ} (hε : ε₁ ≤ ε₂) :
    approxProj S x ε₁ ⊆ approxProj S x ε₂ := by
  intro z hz
  refine ⟨hz.1, ?_⟩
  exact lt_of_lt_of_le hz.2 (add_le_add_left hε ((infDist x S) ^ 2))

/-- Any point of `S` is an `ε`-approximate projection of itself for `ε > 0`. -/
theorem approximate_projection_formula (S : Set H) (x : H) (hx : x ∈ S) {ε : ℝ}
    (hε : 0 < ε) :
    x ∈ approxProj S x ε := by
  refine ⟨hx, ?_⟩
  have hdist : infDist x S = 0 := infDist_eq_zero_of_mem hx
  simpa [approxProj, hdist] using hε

/-- Nonemptiness of approximate projections on nonempty sets for positive tolerance. -/
theorem approxProj_nonempty (S : Set H) (x : H) (hx : x ∈ S) {ε : ℝ} (hε : 0 < ε) :
    (approxProj S x ε).Nonempty :=
  ⟨x, approximate_projection_formula S x hx hε⟩

end MoreauSweeping
