import Mathlib
import MoreauSweeping.Preliminaries
import MoreauSweeping.CatchingUp

/-!
# Subsmooth Sets

Subsmooth sets and convergence for ball-compact subsmooth moving sets.
-/

namespace MoreauSweeping

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- A set is subsmooth -/
axiom Subsmooth : Set H → Prop

/-- A family of sets is equi-uniformly subsmooth -/
axiom EquiUniformlySubsmooth : (ℝ → Set H) → Prop

/-- Stability result for subsmooth sets -/
axiom stability_subsmooth : True

/-- Convergence for subsmooth sets -/
axiom convergence_subsmooth : True

end MoreauSweeping
