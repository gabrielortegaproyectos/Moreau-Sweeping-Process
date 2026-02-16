import Mathlib

/-!
# Preliminaries

Basic definitions for the catching-up algorithm for Moreau's sweeping processes.
-/

namespace MoreauSweeping

variable {H : Type*} [NormedAddCommGroup H]

/-- `approxProj S x ε` is the set of `ε`-approximate projections of `x` on `S`. -/
def approxProj (S : Set H) (x : H) (ε : ℝ) : Set H :=
  {z : H | z ∈ S ∧ ∀ y ∈ S, ‖x - z‖ ^ 2 < ‖x - y‖ ^ 2 + ε}

/-- A point is always at least as close to itself as to any other point of the set. -/
theorem proximal_normal_characterization (S : Set H) (x y : H) (hx : x ∈ S) (hy : y ∈ S) :
    ‖x - x‖ ≤ ‖x - y‖ := by
  have _ : x ∈ S := hx
  have _ : y ∈ S := hy
  simpa using norm_nonneg (x - y)

/-- Monotonicity in the approximation error. -/
theorem approxProj_mono (S : Set H) (x : H) {ε₁ ε₂ : ℝ} (hε : ε₁ ≤ ε₂) :
    approxProj S x ε₁ ⊆ approxProj S x ε₂ := by
  intro z hz
  refine ⟨hz.1, ?_⟩
  intro y hy
  have hz' : ‖x - z‖ ^ 2 < ‖x - y‖ ^ 2 + ε₁ := hz.2 y hy
  linarith

/-- Any point of `S` is an `ε`-approximate projection of itself for `ε > 0`. -/
theorem approximate_projection_formula (S : Set H) (x : H) (hx : x ∈ S) {ε : ℝ}
    (hε : 0 < ε) :
    x ∈ approxProj S x ε := by
  refine ⟨hx, ?_⟩
  intro y _
  have hySq : 0 ≤ ‖x - y‖ ^ 2 := sq_nonneg ‖x - y‖
  have hεle : ε ≤ ‖x - y‖ ^ 2 + ε := by linarith
  have hlt : 0 < ‖x - y‖ ^ 2 + ε := lt_of_lt_of_le hε hεle
  simpa using hlt

/-- Nonemptiness of approximate projections on nonempty sets for positive tolerance. -/
theorem approxProj_nonempty (S : Set H) (x : H) (hx : x ∈ S) {ε : ℝ} (hε : 0 < ε) :
    (approxProj S x ε).Nonempty :=
  ⟨x, approximate_projection_formula S x hx hε⟩

end MoreauSweeping
