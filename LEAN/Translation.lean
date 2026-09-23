/-
  EdgeSingularity/Translation.lean

  Sparenberg's translation theorem for actuator surfaces with constant normal
  load (J. A. Sparenberg 1984; G. A. M. van Kuik, Diss. TU Eindhoven 1991,
  Appendix A), formalised in Lean 4 / Mathlib in a weak (distributional) form.

  Setting (2D actuator strip, Cartesian plane ℝ × ℝ).  A strip S carrying the
  constant normal load F is the force distribution
      ⟨f_S, φ⟩ = F ∫_S φ·n ds         (φ : test vector field)
  i.e. F times the flux of φ through S.  Two strips with the same edges are
  compared:
      S₁ : the flat strip  x = a₁,  a₂ ≤ y ≤ b₂,   normal +e_x ;
      S₂ : the "displaced" strip consisting of the three remaining sides of the
           rectangle [a₁,b₁] × [a₂,b₂], with the normal continued consistently
           (= the outward normal of the rectangle).  This is exactly the pair of
           surfaces in van Kuik's Fig. A.2 (flat surface vs. surface following the
           slipstream and closing downstream).
  Both strips have the same edges (a₁,a₂) and (a₁,b₂).

  Theorems proved:

    T1  conservative_difference :
        ⟨f_{S₂}, φ⟩ − ⟨f_{S₁}, φ⟩ = F ∫_{rect} div φ        for every C¹ field φ,
        i.e.  f_{S₂} − f_{S₁} = −∇(F·1_{rect})  in the sense of distributions:
        the two force fields differ by the gradient of a scalar (a conservative
        field, absorbed into the pressure).                       [Gauss theorem]

    T2  div_rot_eq_zero :
        div (rot ψ) = 0  for every C² scalar ψ, rot ψ := (∂_y ψ, −∂_x ψ).
                                                          [Schwarz/Clairaut]

    T3  translation_theorem :
        ⟨curl f_{S₂}, ψ⟩ = ⟨curl f_{S₁}, ψ⟩  for every C² test function ψ,
        where ⟨curl f_S, ψ⟩ := ⟨f_S, rot ψ⟩ is the weak (2D) curl.
        The vorticity-generating part ∇×f of the load is the same for both
        strips — van Kuik's (A-2)–(A-4).  Since ∇×f is the only way the force
        field enters the vorticity equation Dω/Dt = ρ⁻¹∇×f, the flow induced by
        a constant normal load depends only on F and on the position of the
        edges, not on the position of the surface.

  What is NOT formalised: the step "same ∇×f and same boundary conditions ⇒
  same velocity field" (uniqueness for the Euler equations), which van Kuik
  also only asserts.

  A. P. Schaffarczyk / Claude, September 2026.
-/
import Mathlib.MeasureTheory.Integral.DivergenceTheorem
import Mathlib.Analysis.Calculus.FDeriv.Symmetric
import Mathlib.Analysis.Calculus.ContDiff.Basic

open MeasureTheory Set

noncomputable section

namespace Translation

/-! ### Fluxes through the two strips and the weak divergence -/

/-- Flux of the vector field `φ` through the flat strip `S₁ = {a₁} × [a₂,b₂]`,
normal `+e_x`. -/
def flux₁ (a b : ℝ × ℝ) (φ : ℝ × ℝ → ℝ × ℝ) : ℝ :=
  ∫ y in a.2..b.2, (φ (a.1, y)).1

/-- Flux of `φ` through the displaced strip `S₂` = bottom ∪ right ∪ top side of the
rectangle, with the outward normal (`-e_y`, `+e_x`, `+e_y`). -/
def flux₂ (a b : ℝ × ℝ) (φ : ℝ × ℝ → ℝ × ℝ) : ℝ :=
  (∫ x in a.1..b.1, (φ (x, b.2)).2) - (∫ x in a.1..b.1, (φ (x, a.2)).2)
    + ∫ y in a.2..b.2, (φ (b.1, y)).1

/-- Action of the constant-load force distribution of strip `Sᵢ` on a test field. -/
def force₁ (F : ℝ) (a b : ℝ × ℝ) (φ : ℝ × ℝ → ℝ × ℝ) : ℝ := F * flux₁ a b φ
def force₂ (F : ℝ) (a b : ℝ × ℝ) (φ : ℝ × ℝ → ℝ × ℝ) : ℝ := F * flux₂ a b φ

/-- Divergence of a vector field on the plane. -/
def div (φ : ℝ × ℝ → ℝ × ℝ) (x : ℝ × ℝ) : ℝ :=
  fderiv ℝ (fun p => (φ p).1) x (1, 0) + fderiv ℝ (fun p => (φ p).2) x (0, 1)

/-- Rotated gradient of a scalar field, `rot ψ = (∂_y ψ, −∂_x ψ)`; `⟨f, rot ψ⟩` is the
weak 2D curl of `f` tested against `ψ`. -/
def rot (ψ : ℝ × ℝ → ℝ) (x : ℝ × ℝ) : ℝ × ℝ :=
  (fderiv ℝ ψ x (0, 1), -(fderiv ℝ ψ x (1, 0)))

/-- Weak curl of the force distribution of strip `Sᵢ`. -/
def curlForce₁ (F : ℝ) (a b : ℝ × ℝ) (ψ : ℝ × ℝ → ℝ) : ℝ := force₁ F a b (rot ψ)
def curlForce₂ (F : ℝ) (a b : ℝ × ℝ) (ψ : ℝ × ℝ → ℝ) : ℝ := force₂ F a b (rot ψ)

/-! ### T1: the two loads differ by a conservative field -/

lemma contDiff_fst_comp {φ : ℝ × ℝ → ℝ × ℝ} (hφ : ContDiff ℝ 1 φ) :
    ContDiff ℝ 1 (fun p => (φ p).1) := contDiff_fst.comp hφ

lemma contDiff_snd_comp {φ : ℝ × ℝ → ℝ × ℝ} (hφ : ContDiff ℝ 1 φ) :
    ContDiff ℝ 1 (fun p => (φ p).2) := contDiff_snd.comp hφ

/-- For a `C¹` function `g`, `x ↦ fderiv ℝ g x v` is continuous. -/
lemma continuous_fderiv_apply {g : ℝ × ℝ → ℝ} (hg : ContDiff ℝ 1 g) (v : ℝ × ℝ) :
    Continuous (fun x => fderiv ℝ g x v) := by
  have h : Continuous (fderiv ℝ g) := hg.continuous_fderiv one_ne_zero
  exact (ContinuousLinearMap.apply ℝ ℝ v).continuous.comp h

lemma continuous_div {φ : ℝ × ℝ → ℝ × ℝ} (hφ : ContDiff ℝ 1 φ) : Continuous (div φ) :=
  (continuous_fderiv_apply (contDiff_fst_comp hφ) (1, 0)).add
    (continuous_fderiv_apply (contDiff_snd_comp hφ) (0, 1))

/-- **T1 (Gauss).**  `⟨f_{S₂}, φ⟩ − ⟨f_{S₁}, φ⟩ = F ∫_{[a,b]} div φ` for every `C¹` field `φ`:
the force fields of the flat and of the displaced strip differ by the gradient of the
scalar `−F·1_{[a,b]}`, i.e. by a conservative field. -/
theorem conservative_difference (F : ℝ) (a b : ℝ × ℝ) (hab : a ≤ b)
    (φ : ℝ × ℝ → ℝ × ℝ) (hφ : ContDiff ℝ 1 φ) :
    force₂ F a b φ - force₁ F a b φ = F * ∫ x in Icc a b, div φ x := by
  have h1 := contDiff_fst_comp hφ
  have h2 := contDiff_snd_comp hφ
  have hdiv : (∫ x in Icc a b, div φ x) = flux₂ a b φ - flux₁ a b φ := by
    unfold div flux₁ flux₂
    have := integral_divergence_prod_Icc_of_hasFDerivAt_of_le
      (fun p => (φ p).1) (fun p => (φ p).2)
      (fun x => fderiv ℝ (fun p => (φ p).1) x) (fun x => fderiv ℝ (fun p => (φ p).2) x)
      a b hab h1.continuous.continuousOn h2.continuous.continuousOn
      (fun x _ => (h1.differentiable one_ne_zero x).hasFDerivAt)
      (fun x _ => (h2.differentiable one_ne_zero x).hasFDerivAt)
      ((continuous_div hφ).continuousOn.integrableOn_compact isCompact_Icc)
    linarith [this]
  unfold force₁ force₂
  rw [hdiv]; ring

/-! ### T2: div (rot ψ) = 0 -/

/-- `x ↦ fderiv ℝ ψ x v` is differentiable with derivative `w ↦ fderiv (fderiv ψ) x w v`. -/
lemma hasFDerivAt_fderiv_apply {ψ : ℝ × ℝ → ℝ} (hψ : ContDiff ℝ 2 ψ) (v x : ℝ × ℝ) :
    HasFDerivAt (fun y => fderiv ℝ ψ y v) ((fderiv ℝ (fderiv ℝ ψ) x).flip v) x := by
  have hd : ContDiff ℝ 1 (fderiv ℝ ψ) := hψ.fderiv_right (m := 1) (by norm_num)
  have hc : HasFDerivAt (fderiv ℝ ψ) (fderiv ℝ (fderiv ℝ ψ) x) x :=
    (hd.differentiable one_ne_zero x).hasFDerivAt
  have h := hc.clm_apply (hasFDerivAt_const v x)
  simpa using h

/-- Derivative of `x ↦ fderiv ℝ ψ x v` in direction `w` is the second derivative
`fderiv (fderiv ψ) x w v`. -/
lemma fderiv_fderiv_apply {ψ : ℝ × ℝ → ℝ} (hψ : ContDiff ℝ 2 ψ) (v w x : ℝ × ℝ) :
    fderiv ℝ (fun y => fderiv ℝ ψ y v) x w = fderiv ℝ (fderiv ℝ ψ) x w v := by
  rw [(hasFDerivAt_fderiv_apply hψ v x).fderiv]
  simp

lemma contDiff_fderiv_apply {ψ : ℝ × ℝ → ℝ} (hψ : ContDiff ℝ 2 ψ) (v : ℝ × ℝ) :
    ContDiff ℝ 1 (fun y => fderiv ℝ ψ y v) := by
  have hd : ContDiff ℝ 1 (fderiv ℝ ψ) := hψ.fderiv_right (m := 1) (by norm_num)
  exact hd.clm_apply contDiff_const

lemma contDiff_rot {ψ : ℝ × ℝ → ℝ} (hψ : ContDiff ℝ 2 ψ) : ContDiff ℝ 1 (rot ψ) :=
  (contDiff_fderiv_apply hψ (0, 1)).prodMk (contDiff_fderiv_apply hψ (1, 0)).neg

/-- **T2 (Schwarz).**  `div (rot ψ) = 0` for every `C²` scalar field `ψ`. -/
theorem div_rot_eq_zero {ψ : ℝ × ℝ → ℝ} (hψ : ContDiff ℝ 2 ψ) (x : ℝ × ℝ) :
    div (rot ψ) x = 0 := by
  unfold div rot
  simp only
  have hsymm : IsSymmSndFDerivAt ℝ ψ x :=
    hψ.contDiffAt.isSymmSndFDerivAt (by simp [minSmoothness_of_isRCLikeNormedField])
  have e1 : fderiv ℝ (fun p : ℝ × ℝ => fderiv ℝ ψ p (0, 1)) x (1, 0)
      = fderiv ℝ (fderiv ℝ ψ) x (1, 0) (0, 1) := fderiv_fderiv_apply hψ _ _ _
  have e2 : fderiv ℝ (fun p : ℝ × ℝ => -(fderiv ℝ ψ p (1, 0))) x (0, 1)
      = -(fderiv ℝ (fderiv ℝ ψ) x (0, 1) (1, 0)) := by
    have h : HasFDerivAt (fun p : ℝ × ℝ => -(fderiv ℝ ψ p (1, 0)))
        (-((fderiv ℝ (fderiv ℝ ψ) x).flip (1, 0))) x :=
      (hasFDerivAt_fderiv_apply hψ (1, 0) x).neg
    rw [h.fderiv]
    simp
  rw [e1, e2, hsymm.eq (1, 0) (0, 1)]
  ring

/-! ### T3: Sparenberg's translation theorem (weak form) -/

/-- **T3 (Sparenberg / van Kuik App. A).**  The weak curl of the constant-load force
distribution is the same for the flat strip `S₁` and the displaced strip `S₂` with the
same edges: `⟨curl f_{S₂}, ψ⟩ = ⟨curl f_{S₁}, ψ⟩` for every `C²` test function `ψ`. -/
theorem translation_theorem (F : ℝ) (a b : ℝ × ℝ) (hab : a ≤ b)
    (ψ : ℝ × ℝ → ℝ) (hψ : ContDiff ℝ 2 ψ) :
    curlForce₂ F a b ψ = curlForce₁ F a b ψ := by
  unfold curlForce₁ curlForce₂
  have h := conservative_difference F a b hab (rot ψ) (contDiff_rot hψ)
  have hz : (∫ x in Icc a b, div (rot ψ) x) = 0 := by
    have : (fun x => div (rot ψ) x) = fun _ => (0 : ℝ) := funext (div_rot_eq_zero hψ)
    rw [this]; simp
  rw [hz, mul_zero, sub_eq_zero] at h
  exact h

end Translation

end
