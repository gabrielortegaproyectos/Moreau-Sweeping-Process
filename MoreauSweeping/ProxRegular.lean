import Mathlib
import MoreauSweeping.Preliminaries
import MoreauSweeping.CatchingUp

/-!
# Prox-Regular Sets

Uniformly prox-regular sets and convergence of the catching-up algorithm for this case.
-/

namespace MoreauSweeping

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- A set is uniformly prox-regular -/
axiom UniformlyProxRegular : Set H → Prop

/-- Characterization of prox-regularity -/
axiom prox_regular_characterization : True

/-- Convergence of approximate projections for prox-regular sets -/
axiom approx_proj_convergence_prox_regular : True

/-- Quasi-Lipschitz property of approximate projections -/
axiom quasi_lipschitz_approx_proj : True

/-- Convergence for prox-regular sets -/
axiom convergence_prox_regular : True

end MoreauSweeping
