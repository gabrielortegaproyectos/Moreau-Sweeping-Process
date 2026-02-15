import Mathlib
import MoreauSweeping.Preliminaries
import MoreauSweeping.CatchingUp

/-!
# Fixed Set Case

Convergence of the catching-up algorithm for a fixed ball-compact closed set.
-/

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- Convergence for fixed sets -/
axiom convergence_fixed_set : True
