import Mathlib
import MoreauSweeping.Preliminaries

/-!
# Catching-Up Algorithm

Definition and properties of the catching-up algorithm with approximate projections.
-/

namespace MoreauSweeping

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- A discrete error schedule is admissible if every error is positive. -/
def AdmissibleError (ε : ℕ → ℝ) : Prop :=
  ∀ k : ℕ, 0 < ε k

/-- Under admissible errors, each point of the set is an approximate projection of itself. -/
theorem self_mem_approxProj_of_admissibleError (S : Set H) (x : H) (hx : x ∈ S)
    (ε : ℕ → ℝ) (hε : AdmissibleError ε) (k : ℕ) :
    x ∈ approxProj S x (ε k) :=
  approximate_projection_formula S x hx (hε k)

/-- Placeholder for algorithm properties -/
axiom algorithm_properties : True

end MoreauSweeping
