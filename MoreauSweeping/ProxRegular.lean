import Mathlib
import MoreauSweeping.Preliminaries
import MoreauSweeping.CatchingUp

/-!
# Prox-Regular Sets

Uniformly prox-regular sets and convergence of the catching-up algorithm for this case.
-/

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- A set is uniformly prox-regular -/
axiom UniformlyProxRegular : Set H → Prop

/-- Characterization of prox-regularity -/
axiom prox_regular_characterization : True

/-- Convergence for prox-regular sets -/
axiom convergence_prox_regular : True
