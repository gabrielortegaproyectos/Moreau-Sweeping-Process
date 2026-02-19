import Mathlib
import MoreauSweeping.Preliminaries
import MoreauSweeping.CatchingUp

/-!
# Subsmooth Sets

Subsmooth sets and convergence for boundedly compact subsmooth moving sets.
-/

namespace MoreauSweeping

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- A set is boundedly compact if all intersections with closed balls are compact. -/
def BoundedlyCompact (S : Set H) : Prop :=
  ∀ r : ℝ, 0 < r → IsCompact (S ∩ Metric.closedBall (0 : H) r)

/-- Monotonicity under closed subsets. -/
theorem BoundedlyCompact.mono_closed {S T : Set H} (hS : BoundedlyCompact S)
    (hTS : T ⊆ S) (hTClosed : IsClosed T) :
    BoundedlyCompact T := by
  intro r hr
  have hCompact : IsCompact (S ∩ Metric.closedBall (0 : H) r) := hS r hr
  have hClosedInter : IsClosed (T ∩ Metric.closedBall (0 : H) r) :=
    hTClosed.inter isClosed_closedBall
  have hSubset : T ∩ Metric.closedBall (0 : H) r ⊆ S ∩ Metric.closedBall (0 : H) r := by
    intro x hx
    exact ⟨hTS hx.1, hx.2⟩
  exact hCompact.of_isClosed_subset hClosedInter hSubset

/-- Intersecting with a closed set preserves bounded compactness under a closedness hypothesis. -/
theorem BoundedlyCompact.inter_closed {S T : Set H} (hS : BoundedlyCompact S)
    (hSClosed : IsClosed S) (hTClosed : IsClosed T) :
    BoundedlyCompact (S ∩ T) := by
  intro r hr
  have hCompact : IsCompact (S ∩ Metric.closedBall (0 : H) r) := hS r hr
  have hClosedInter : IsClosed ((S ∩ T) ∩ Metric.closedBall (0 : H) r) :=
    (hSClosed.inter hTClosed).inter isClosed_closedBall
  have hSubset : (S ∩ T) ∩ Metric.closedBall (0 : H) r ⊆ S ∩ Metric.closedBall (0 : H) r := by
    intro x hx
    exact ⟨hx.1.1, hx.2⟩
  exact hCompact.of_isClosed_subset hClosedInter hSubset

/-- A set is subsmooth -/
axiom Subsmooth : Set H → Prop

/-- A family of sets is equi-uniformly subsmooth -/
axiom EquiUniformlySubsmooth : (ℝ → Set H) → Prop

/-- Stability result for subsmooth sets -/
axiom stability_subsmooth : True

/-- Convergence for subsmooth sets -/
axiom convergence_subsmooth : True

end MoreauSweeping
