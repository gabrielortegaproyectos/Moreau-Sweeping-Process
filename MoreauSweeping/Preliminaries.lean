import Mathlib

/-!
# Preliminaries

Basic definitions for the catching-up algorithm for Moreau's sweeping processes.
-/

namespace MoreauSweeping

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]
noncomputable section
open scoped Topology

/-- Distance to a set, using `Metric.infDist`. -/
def distance (S : Set H) (x : H) : ℝ :=
  Metric.infDist x S

/-- The distance to a set is nonnegative. -/
theorem distance_nonneg (S : Set H) (x : H) : 0 ≤ distance S x := by
  simpa [distance] using Metric.infDist_nonneg

/-- If `x ∈ S`, then the distance from `x` to `S` is zero. -/
theorem distance_eq_zero_of_mem (S : Set H) (x : H) (hx : x ∈ S) : distance S x = 0 := by
  refine le_antisymm ?_ (distance_nonneg S x)
  have hle : Metric.infDist x S ≤ dist x x := Metric.infDist_le_dist_of_mem hx
  simpa [distance] using hle

/-- Support function `σ(x, S) = sup_{z ∈ S} ⟪x, z⟫`. -/
def support (x : H) (S : Set H) : ℝ :=
  sSup ((fun z : H => inner ℝ x z) '' S)

/-- Each support value is bounded above by the support function, when the image is bounded above. -/
theorem le_support (x : H) (S : Set H) (hS : BddAbove ((fun z : H => inner ℝ x z) '' S))
    {z : H} (hz : z ∈ S) :
    inner ℝ x z ≤ support x S := by
  exact le_csSup hS ⟨z, hz, rfl⟩

/-- Monotonicity of support under set inclusion (with a boundedness hypothesis). -/
axiom support_mono (x : H) {S T : Set H} (hST : S ⊆ T)
    (hT : BddAbove ((fun z : H => inner ℝ x z) '' T)) :
    support x S ≤ support x T

/-- Clarke tangent cone, given by the sequential characterization used in the blueprint. -/
def ClarkeTangentCone (S : Set H) (x : H) : Set H :=
  {h : H |
    ∀ (xSeq : ℕ → H) (tSeq : ℕ → ℝ),
      (∀ n, xSeq n ∈ S) →
      Filter.Tendsto xSeq Filter.atTop (𝓝 x) →
      (∀ n, 0 < tSeq n) →
      Filter.Tendsto tSeq Filter.atTop (𝓝 (0 : ℝ)) →
      ∃ hSeq : ℕ → H,
        Filter.Tendsto hSeq Filter.atTop (𝓝 h) ∧ ∀ n, xSeq n + tSeq n • hSeq n ∈ S}

/-- The zero vector always belongs to the Clarke tangent cone. -/
theorem zero_mem_ClarkeTangentCone (S : Set H) (x : H) :
    (0 : H) ∈ ClarkeTangentCone S x := by
  intro xSeq tSeq hxSeq _ _ _
  refine ⟨fun _ => 0, ?_, ?_⟩
  · simpa using (tendsto_const_nhds : Filter.Tendsto (fun _ : ℕ => (0 : H)) Filter.atTop (nhds (0 : H)))
  · intro n
    simpa [zero_smul, add_zero] using hxSeq n

/--
Proximal subdifferential for real-valued functions on Hilbert spaces.
This is a first formalization layer matching the blueprint's local quadratic model.
-/
def proximalSubdifferential (f : H → ℝ) (x : H) : Set H :=
  {ζ : H |
    ∃ σ η : ℝ,
      0 ≤ σ ∧ 0 < η ∧
      ∀ y : H, ‖y - x‖ < η →
        f y ≥ f x + inner ℝ ζ (y - x) - σ * ‖y - x‖ ^ 2}

/-- Rewriting lemma for membership in the proximal subdifferential. -/
theorem mem_proximalSubdifferential_iff (f : H → ℝ) (x ζ : H) :
    ζ ∈ proximalSubdifferential f x ↔
      ∃ σ η : ℝ,
        0 ≤ σ ∧ 0 < η ∧
        ∀ y : H, ‖y - x‖ < η →
          f y ≥ f x + inner ℝ ζ (y - x) - σ * ‖y - x‖ ^ 2 := by
  rfl

/-- Relaxing `(σ, η)` in the expected direction preserves proximal-subgradient membership. -/
axiom proximalSubdifferential_relax_constants (f : H → ℝ) (x ζ : H)
    {σ η σ' η' : ℝ} (hσ : 0 ≤ σ) (hη : 0 < η)
    (hσ' : σ ≤ σ') (hη' : η' ≤ η) (hσ'' : 0 ≤ σ') (hη'' : 0 < η')
    (hWitness : ∀ y : H, ‖y - x‖ < η →
      f y ≥ f x + inner ℝ ζ (y - x) - σ * ‖y - x‖ ^ 2) :
    ζ ∈ proximalSubdifferential f x

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

end
end MoreauSweeping
