/-
  EdgeSingularity/Basic.lean

  Logarithmic singularity of the velocity at the end point of a cut vortex
  sheet (the local model of the disk edge of a uniformly loaded actuator
  disk), formalised in Lean 4 / Mathlib.

  Model (2D, local): a straight vortex sheet of constant strength γ on the
  real segment [0, L] of the complex plane.  Its complex velocity is
      w(z) = u - i v = (γ / (2π i)) ∫₀ᴸ dt / (z - t) ,      z ∉ [0, L].
  Statements proved below (for z with Im z ≠ 0, i.e. off the axis of the
  sheet):

    T1  ∫₀ᴸ (z - t)⁻¹ dt = log z - log (z - L)                (closed form)
    T2  u(z) = (γ/2π) (arg z - arg (z - L)),
        v(z) = (γ/2π) (log‖z‖ - log‖z - L‖)                  (explicit u, v)
    T3  |u(z)| ≤ |γ|                                         (tangential
                                                              velocity bounded)
    T4  v(z) → -∞ (γ > 0) resp. +∞ (γ < 0) as z → 0          (normal velocity
                                                              log-singular)
    T5  v(z) - (γ/2π) log‖z‖  is bounded near z = 0          (the singularity
                                                              is exactly log)

  T4/T5 is the statement "a vortex sheet that ends at a point with non-zero
  strength induces a logarithmically infinite normal velocity there"
  (Muskhelishvili, Singular Integral Equations, §29, for the special case of
  constant density and a straight arc).  Combined with the Bernoulli jump
  condition γ = Δh / V_t,m ≠ 0 of the uniformly loaded disk it shows that the
  slipstream boundary cannot leave the disk edge with a finite slope.

  A. P. Schaffarczyk / Claude, September 2026.
-/
import Mathlib.Analysis.SpecialFunctions.Complex.Log
import Mathlib.Analysis.SpecialFunctions.Complex.LogDeriv
import Mathlib.Analysis.Complex.RealDeriv
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Analysis.SpecialFunctions.Log.Basic

open Complex Filter Topology

noncomputable section

namespace EdgeSingularity

/-- The velocity kernel of the sheet element at `t` seen from `z`. -/
def kernel (z : ℂ) (t : ℝ) : ℂ := (z - (t : ℂ))⁻¹

/-- Points off the real axis are in the slit plane after any real shift. -/
lemma sub_ofReal_mem_slitPlane {z : ℂ} (hz : z.im ≠ 0) (t : ℝ) :
    z - (t : ℂ) ∈ slitPlane := by
  rw [mem_slitPlane_iff]
  right
  simpa using hz

lemma sub_ofReal_ne_zero {z : ℂ} (hz : z.im ≠ 0) (t : ℝ) : z - (t : ℂ) ≠ 0 := by
  intro h
  have : (z - (t : ℂ)).im = 0 := by rw [h]; simp
  simp at this
  exact hz this

/-- Derivative (in the real variable `t`) of `-log (z - t)`. -/
lemma hasDerivAt_neg_log_sub {z : ℂ} (hz : z.im ≠ 0) (t : ℝ) :
    HasDerivAt (fun s : ℝ => -Complex.log (z - (s : ℂ))) (kernel z t) t := by
  have h1 : HasDerivAt (fun w : ℂ => z - w) (-1 : ℂ) (t : ℂ) :=
    (hasDerivAt_id (t : ℂ)).const_sub z
  have h2 : HasDerivAt Complex.log (z - (t : ℂ))⁻¹ (z - (t : ℂ)) :=
    Complex.hasDerivAt_log (sub_ofReal_mem_slitPlane hz t)
  have h3 : HasDerivAt (fun w : ℂ => Complex.log (z - w)) ((z - (t : ℂ))⁻¹ * (-1)) (t : ℂ) :=
    h2.comp (t : ℂ) h1
  have h4 : HasDerivAt (fun w : ℂ => -Complex.log (z - w)) (-((z - (t : ℂ))⁻¹ * (-1))) (t : ℂ) :=
    h3.neg
  have h5 := h4.comp_ofReal
  simp only [kernel]
  convert h5 using 1
  ring

lemma continuous_kernel {z : ℂ} (hz : z.im ≠ 0) : Continuous (kernel z) := by
  unfold kernel
  exact (continuous_const.sub Complex.continuous_ofReal).inv₀ (sub_ofReal_ne_zero hz)

/-- **T1.** Closed form of the sheet integral for `z` off the real axis. -/
theorem integral_kernel {z : ℂ} (hz : z.im ≠ 0) (L : ℝ) :
    ∫ t in (0 : ℝ)..L, kernel z t = Complex.log z - Complex.log (z - L) := by
  have h := intervalIntegral.integral_eq_sub_of_hasDerivAt
    (f := fun s : ℝ => -Complex.log (z - (s : ℂ))) (f' := kernel z)
    (fun t _ => hasDerivAt_neg_log_sub hz t)
    ((continuous_kernel hz).intervalIntegrable 0 L)
  rw [h]
  simp
  ring

/-- Complex velocity `w = u - i v` of a sheet of strength `γ` on `[0, L]`:
`w = (γ / (2π i)) ∫ kernel = (γ/(2π)) · (-i) · ∫ kernel`. -/
def w (γ L : ℝ) (z : ℂ) : ℂ :=
  ((γ / (2 * Real.pi) : ℝ) : ℂ) * (-I) * ∫ t in (0 : ℝ)..L, kernel z t

/-- Tangential velocity component. -/
def u (γ L : ℝ) (z : ℂ) : ℝ := (w γ L z).re

/-- Normal velocity component (`w = u - i v`). -/
def v (γ L : ℝ) (z : ℂ) : ℝ := -(w γ L z).im

/-- **T2 (tangential part).** `u = (γ/2π) (arg z − arg (z − L))`. -/
theorem u_eq {z : ℂ} (hz : z.im ≠ 0) (γ L : ℝ) :
    u γ L z = (γ / (2 * Real.pi)) * (Complex.arg z - Complex.arg (z - L)) := by
  unfold u w
  rw [integral_kernel hz L]
  generalize γ / (2 * Real.pi) = c
  simp [Complex.mul_re, Complex.mul_im, Complex.log_im, Complex.log_re]
  try ring

/-- **T2 (normal part).** `v = (γ/2π) (log‖z‖ − log‖z − L‖)`. -/
theorem v_eq {z : ℂ} (hz : z.im ≠ 0) (γ L : ℝ) :
    v γ L z = (γ / (2 * Real.pi)) * (Real.log ‖z‖ - Real.log ‖z - L‖) := by
  unfold v w
  rw [integral_kernel hz L]
  generalize γ / (2 * Real.pi) = c
  simp [Complex.mul_re, Complex.mul_im, Complex.log_im, Complex.log_re]
  try ring

/-- **T3.** The tangential velocity stays bounded up to the end point: `|u| ≤ |γ|`. -/
theorem abs_u_le {z : ℂ} (hz : z.im ≠ 0) (γ L : ℝ) : |u γ L z| ≤ |γ| := by
  rw [u_eq hz]
  have h1 : |Complex.arg z - Complex.arg (z - L)| ≤ 2 * Real.pi := by
    have a1 := Complex.abs_arg_le_pi z
    have a2 := Complex.abs_arg_le_pi (z - L)
    calc |Complex.arg z - Complex.arg (z - L)|
        ≤ |Complex.arg z| + |Complex.arg (z - L)| := abs_sub _ _
      _ ≤ Real.pi + Real.pi := add_le_add a1 a2
      _ = 2 * Real.pi := by ring
  rw [abs_mul, abs_div]
  have hpi : 0 < 2 * Real.pi := by positivity
  calc |γ| / |2 * Real.pi| * |Complex.arg z - Complex.arg (z - L)|
      ≤ |γ| / |2 * Real.pi| * (2 * Real.pi) := by
        apply mul_le_mul_of_nonneg_left h1 (by positivity)
    _ = |γ| := by rw [abs_of_pos hpi]; field_simp

/-- The regular part of the normal velocity. -/
def vreg (γ L : ℝ) (z : ℂ) : ℝ := -(γ / (2 * Real.pi)) * Real.log ‖z - L‖

lemma v_split {z : ℂ} (hz : z.im ≠ 0) (γ L : ℝ) :
    v γ L z = (γ / (2 * Real.pi)) * Real.log ‖z‖ + vreg γ L z := by
  rw [v_eq hz]; unfold vreg; ring

/-- **T5.** The regular part is continuous at the end point (for `L ≠ 0`),
hence bounded near it: the singularity of `v` is exactly `(γ/2π) log‖z‖`. -/
theorem continuousAt_vreg (γ : ℝ) {L : ℝ} (hL : L ≠ 0) :
    ContinuousAt (vreg γ L) 0 := by
  unfold vreg
  apply ContinuousAt.mul continuousAt_const
  apply ContinuousAt.log
  · exact (continuous_id.sub continuous_const).norm.continuousAt
  · simp [hL]

/-- **T4.** For `γ > 0` the normal velocity tends to `-∞` at the end point
(approached off the axis of the sheet); i.e. it is logarithmically singular. -/
theorem tendsto_v_atBot {γ L : ℝ} (hγ : 0 < γ) (hL : L ≠ 0) :
    Tendsto (v γ L) (𝓝[{z : ℂ | z.im ≠ 0}] 0) atBot := by
  -- v = (γ/2π) log‖z‖ + vreg, with log‖z‖ → -∞ and vreg → vreg 0
  have hlog : Tendsto (fun z : ℂ => Real.log ‖z‖) (𝓝[{z : ℂ | z.im ≠ 0}] 0) atBot := by
    have h1 : Tendsto (fun z : ℂ => ‖z‖) (𝓝[{z : ℂ | z.im ≠ 0}] 0) (𝓝[>] 0) := by
      apply tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within
      · have h0 : Tendsto (fun z : ℂ => ‖z‖) (𝓝 (0 : ℂ)) (𝓝 ‖(0 : ℂ)‖) :=
          continuous_norm.tendsto (0 : ℂ)
        simpa using h0.mono_left nhdsWithin_le_nhds
      · filter_upwards [self_mem_nhdsWithin] with z hz
        have hz' : z ≠ 0 := by
          intro h
          apply hz
          rw [h]
          simp
        exact norm_pos_iff.mpr hz'
    exact Real.tendsto_log_nhdsGT_zero.comp h1
  have hc : Tendsto (vreg γ L) (𝓝[{z : ℂ | z.im ≠ 0}] 0) (𝓝 (vreg γ L 0)) :=
    (continuousAt_vreg γ hL).tendsto.mono_left nhdsWithin_le_nhds
  have hmul : Tendsto (fun z : ℂ => (γ / (2 * Real.pi)) * Real.log ‖z‖)
      (𝓝[{z : ℂ | z.im ≠ 0}] 0) atBot := by
    have hpos : 0 < γ / (2 * Real.pi) := by positivity
    exact Tendsto.const_mul_atBot hpos hlog
  have := hmul.atBot_add hc
  refine this.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with z hz
  exact (v_split hz γ L).symm

/-- **T4'.** For `γ < 0` the normal velocity tends to `+∞`. -/
theorem tendsto_v_atTop {γ L : ℝ} (hγ : γ < 0) (hL : L ≠ 0) :
    Tendsto (v γ L) (𝓝[{z : ℂ | z.im ≠ 0}] 0) atTop := by
  have h := tendsto_v_atBot (γ := -γ) (L := L) (by linarith) hL
  have hneg : ∀ z : ℂ, z.im ≠ 0 → v γ L z = -(v (-γ) L z) := by
    intro z hz
    rw [v_eq hz, v_eq hz]
    ring
  have : Tendsto (fun z : ℂ => -(v (-γ) L z)) (𝓝[{z : ℂ | z.im ≠ 0}] 0) atTop :=
    tendsto_neg_atBot_atTop.comp h
  refine this.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with z hz
  exact (hneg z hz).symm

end EdgeSingularity

end
