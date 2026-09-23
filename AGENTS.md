# containment — for agents

The why behind the code. `CLAUDE.md` has the commands; this file has the
reasoning, the traps that were actually hit, and what is and is not known.
Built 2026-09-23 from Allan's request, "a ball of plasma being constrained by
an electromagnetic field, whilst trying to decompose and expand". The spec is
`~/Projects/resolume/specs/SPEC-containment.md`. Every number below was
measured by `cttest` on this Mac (M4 Max, macOS 26.4).

## The one idea

The picture is a ball of hot plasma, held in a magnetic bottle, and the plasma
obeys 2.5-D compressible ideal MHD. The clip is the plasma's colour, carried
with its mass. Nothing is drawn:

- **the ball** is a hot, dense region, bright as ρ²√T (bremsstrahlung);
- **its expansion** is ∇p pushing out against magnetic pressure;
- **its edge and ringing** are p = B²/2 and the fast waves;
- **its decomposition** is interchange (magnetic Rayleigh–Taylor) driven by an
  effective gravity standing in for bad curvature, plus leakage through the
  cusps;
- **the fireball** is the free expansion when the coils quench.

## The solver

Units: frame height 1, μ₀ = 1, ρ₀ = B_ref = 1, so τ_A = 1. γ = 5/3.

State, four RGBA32F textures ping-ponged:

    A ( rho, rho u, rho v, rho w )
    B ( E, Bx, By, Bz )
    C ( psi, summed signal speed, rho K, flags )   K = p / rho^gamma
    D ( rho r, rho g, rho b, rho chi )             the picture; chi = the ball marker

One substep, six passes:

1. **reduce**: 4×4 blocks down to ≤ 16 texels, max summed signal speed, the
   floor flags and the entropy-branch flags.
2. **clock**: finishes the reduction and writes dt = min(C Δx / s_max, what
   is left of the frame / substeps left), the frame's time done, and c_h.
3. **predict**: MUSCL–Hancock's half step with MC-limited slopes.
4. **flux x** / **flux y**: each face once. U* is reconstructed with the MC
   limiter, GLM (Dedner eq. 42) sets the face's normal B and ψ, and **HLLD**
   (Miyoshi & Kusano 2005) with both degenerate cases handled their way.
   Resistive E = ηJ and the wall overrides are applied here.
5. **update**: flux differences, the EGLM energy term, the body force at the
   half step, ψ damping, the dual-energy switch, cooling and the floors. It
   writes the next state and each cell's signal speed.

**HLLD, not HLL.** HLLD behaved in float from the start. HLL stays as a
switch only for the negative control, which it fails: the contact goes from
4 cells wide to 6, and the waves sit 6–7 cells off.

**dt lives on the GPU.** The CPU plans the substep count from last frame's
speed (read back at the start of the next frame, so it never stalls the
frame it is in), plus 10%, capped at 48. If the plasma sped up since, each
substep still takes only its CFL step and the frame covers less time.
Simulated time runs slow, never unstable, and Diag logs it. After an Ignite
or a pellet the plan is made from a synchronous readback.

## The bottle

- **Coils**: N line currents of alternating sign on a circle of radius
  R_c ≥ 1.05× the frame's half-diagonal, normalised so |B| at radius 0.5 in
  the first gap is Field. Plus a uniform guide Bz. The closed form is written
  in `Physics.cpp` and in the shader, both marked `//= mirrored`.
- **Open**: the ghost's plasma is zero-gradient but never below the ambient
  background. Its field is B_vac(ghost) + (B − B_vac)(the cell it copies).
  ψ leaves.
- **Wall**: the plasma is mirrored with the whole velocity reversed (no
  slip). The field is B_vac(ghost) plus the mirrored cell's perturbation. The
  flux pass zeroes the mass, energy, picture, entropy and all B fluxes through
  the face. The coils' slow changes (a turn, a quench) are carried through the
  vessel: once a frame, the sources pass adds B_vac(now) − B_vac(last frame)
  to every cell, with the pressure left alone. That addition is curl-free and
  divergence-free.
- **Curvature**: g_eff radial from the ball's centre, scaled by T / T_ball
  (the curvature drift goes as T / R_c), tapered to 0 inside half a ball
  radius.
- **Drive**: 12 Fourier modes of a stream function with OU amplitudes (PCG,
  correlation time 0.6 τ_A), windowed by a Gaussian of 1.5 ball radii. The
  force is the curl of the windowed stream function, so it is
  divergence-free exactly.

**Default: Wall, guide field 1.4 × Field, Field 0.5, six poles.** The guide
field confines the ball and the cusp shapes it. Two reasons. With the cusp
alone, the ball drains through the cusps within about 1.5 τ_A (correct
physics, but no orb). And the cusp field at the frame's corners, 4–7× Field,
sets v_A and so the step: a weak cusp is what brought the default to ~22
substeps a frame.

## The light

- **Emission** on the grid: ρ²√T × the tracer colour, or × a mix towards the
  temperature ramp (Hot: deep red → violet → white; Aurora: the 557.7 nm
  oxygen green, going white-green when hot).
- **Glow**: the emission is box-reduced 2× exactly until its short side is at
  most 256. Three cascaded, normalised Gaussians follow (σ 0.012, 0.045 and
  0.15 of the frame height), each mirrored at a half sample. A symmetric
  kernel reflected that way is a symmetric operator, so it conserves light to
  rounding. The core loses the fraction the glare takes.
- **Field lines**: A_z from ∇²A = −J_z, one multigrid V-cycle a frame
  (weighted Jacobi, the coils' potential as the Dirichlet boundary), seeded
  from last frame. Contours are one pixel wide, lit by the plasma on them, and
  fade where they are closer than two pixels.
- **Composite**: exposure, then the fleet's shoulder (linear to 0.8, then an
  exponential approach to 1). Mix 0 returns the input by `texelFetch`, bit
  for bit.

## The traps

Ordered by how much time they cost.

**GLM's c_h at the summed signal speed blows up.** dt is set by the sum of
the two directions' speeds (unsplit stability). GLM's waves also cross both
directions at once, so c_h = sum puts them at a Courant number of 1.6. The
oblique Alfvén wave went to NaN in a few steps. c_h = sum / 2 fixed it.

**Plain GLM makes low-β pressure negative.** GLM moves the normal B without
touching E, so the pressure pays for the change in B²/2. At β ≈ 10⁻³ (the
background beside the coils) that is bigger than p. The vacuum field alone
went to the floors at the open edge. The fix is Dedner's EGLM energy term in
its exact per-step form: E gets back what the ψ fluxes changed in B²/2.

**A mirrored field at a wall is a kink.** Where coil field lines cross the
wall, mirroring B puts a rotational discontinuity at the face every step,
and the plasma blew up. Copying B (zero gradient) broke the vacuum
equilibrium instead: 19 M floor clamps with no ball at all. What works is
the coils' own field at the ghost plus the mirrored cell's departure from
it.

**A slip wall cannot be a perfect conductor where field crosses it.**
E_t = v_t B_n ≠ 0, and zeroing the face's B flux against it made a current
sheet that flung the plasma (|v| 80 in 50 frames). The wall is no-slip,
which is line-tied.

**Driving the wall with the coils' E cannot work either.** A rotating coil
set needs E_z = −∂A/∂t at the wall, and a no-slip, line-tied ideal wall
cannot carry it: another current sheet and a runaway. The coils' change is
added through the volume instead.

**The low-β background: the dual-energy switch.** Even with EGLM, the
scheme's truncation error in B²/2 exceeds p where β ≈ 10⁻³. There p comes
from the entropy K, advected as a passive scalar (Ryu et al.; Balsara &
Spicer). Energy decides everywhere else and resynchronises K, so shocks heat
as they should. On the default look 58% of cell-steps take the entropy branch
(the background), and floors then fire on 7×10⁻⁶ of cell-steps.

**The open boundary drains the background.** A zero-gradient ghost lets the
tenuous plasma stream out along field lines that cross the edge, leaving
β ≈ 0 holes that only the floors held up (millions of clamps, then a
runaway). The ghost now never holds less than the ambient background. Open
is still not indefinitely stable: on the default look it runs away after
about 3 τ_A as |B| climbs at the edges (see What is not done). A pure-vacuum
field ghost is worse: NaN by t = 4.5.

**Uniform g_eff stratifies the background.** A radial gravity of fixed size
also acts on the cold background, whose scale height T_bg / g is short. It
emptied the space round the ball (ρ → 0.01). Scaling g_eff by T / T_ball
keeps the force density ∝ p, continuous across an isobaric interface and
jumping at the ball's edge, where interchange is driven. A drive that stirred
the whole frame did the same; the drive is windowed about the ball.

**The per-substep cost is latency, not work.** Seven dependent passes a
substep, and a long register-heavy update (four HLLD solves per fragment),
never filled the GPU at the default grid: 0.3 ms a substep whatever the grid
size. Splitting the faces into their own passes cut Detail 128 from 7.2 to
5.1 ms. Reducing in 16×16 blocks, or finishing the reduction in one fragment,
made it 2.4× slower: a fragment's serial loop is latency too.

**GLSL `tanh` is NaN past |x| ≈ 44 on this GPU** (inf / inf). A top-hat
ball's far field is hundreds of cells out, so every cell outside it came out
NaN. The argument is clamped.

**A sampler declared and left on texture 0 is logged** as "unloadable". The
driver's unit number in the message is the sampler's, not the bound unit's.
Every declared sampler gets a real texture.

**The field-line span measured in a cusp gap is zero.** There B is radial,
so ∂A/∂r = −B_φ = 0 and a multipole's A is 0. The lines came out infinitely
close, every pixel lit, and Field Line Count was dead (the sweep found it).
The span is now the largest |A − A_centre| round the circle. Lines closer
than two pixels fade.

**Harness traps that looked like physics bugs:**
- a loaded state has no signal speeds, so dt came out infinite: a load now
  runs a zero sources pass;
- the drive advanced only per frame, so `StepForTest` never stirred:
  `--frozen`'s "the pattern moved" guard caught it;
- the last substep's floor flags are only counted by a later reduction, so
  `StepForTest` reduces once more at the end;
- **Rayleigh–Taylor seeded with velocity alone grows as cosh(γt), not e^γt.**
  Fitting ln v read γ tanh γt, 5–6% low, and it looked like a resolution
  problem because it did not converge. Fitting acosh(v/v₀)/t converges:
  1.25% → 0.44%;
- the Brio–Wu "published solution" is read off plots; see below for what
  stands in for it;
- the first 1-2-3 problem (Mach 3) sits exactly at the vacuum threshold and
  never needed a floor, so its negative control could not fail;
- `--resist`'s "second order" flow was first order: the pressure hole that
  balances the field bump stays put as B diffuses, and the gas squeezes the
  frozen Bz. At p = 10 that was 3% of the bump; at ρ = 100, p = 1000 it is
  10⁻⁵;
- the first shader mutation (in an HLLD branch Brio–Wu barely visits) passed.
  The shipped one is in the MC limiter and the check asserts the text exists.

## Every numeric check, and where its tolerance comes from

Asked of each: would it still hold on another rasteriser, at another raster?
Every physics check runs on the grid, whose size is set by Detail and not by
the output raster. Each prints the numbers it compared.

| check | bound | why that number |
| --- | --- | --- |
| `--briowu` positions | 3 cells | a captured discontinuity spreads over 2–3 cells and its steepest point is within one; the heads, found at 0.1% of the state, sit up to 3 cells inside |
| `--briowu` plateaus | 2% | overshoot of a second-order slow shock over a cell or two; the plateaus are ≥ 20 cells wide |
| `--briowu` contact | 4 cells, 10–90%, Detail 512 | HLLD resolves the contact as its own wave; HLL has none and measures 6 |
| `--alfven` speed | (k Δx)² | the phase error of a second-order scheme |
| `--alfven` rotation | 0.02 rad, correlation < −0.99 | an exact solution; a quarter turn either way is O(1) |
| `--alfven` order | ≥ 1.5 | MC clips each extremum to first order, so the global order of a smooth wave is between 1.5 and 2 |
| `--conserve` | 3 u × substeps, u = 2⁻²⁴ | each substep rounds each cell's value ~3 times by ≤ u of itself; summed over cells that is 3u of the total (the cell count cancels), linear over steps |
| `--divb` | max 0.1, rms 0.01 of \|B\|/Δx | O(1) is a monopole a cell across; GLM holds it to truncation level. Measured max 0.055 / 0.084 at Detail 128 / 256; with GLM off 0.060 / 0.21 -- **only Detail 256 tells them apart**: on a 128 grid this run's field is smooth enough that the scheme's own divergence stays small either way |
| `--balance` pressure, Bz | max ρv² / p_total + 0.5% | residual ringing is balanced by ρ dv/dt ~ ρv²/L; half a percent for the edge's own width |
| `--balance` edge | 3 cells | the edge from the marker's area against flux conservation |
| `--balance` β = 1 | inside the edge layer ± 1 cell | β = 1 lies where 2p = B², inside the layer between marker 0.9 and 0.1 |
| `--rt` rate | (kL/2) gkA/γ² + g/(k c_s²) + 2%, L = 2 cells | a diffuse interface weakens only the buoyancy term; plus compressibility and the fit |
| `--rt` stable | peak < 2× the seed | the unstable cases reach > 20× in the same time |
| `--cusp` angles | half a bin + 2Δx/r | the histogram's bin and the grid at the sampling circle |
| `--cusp` contrast | 1.5× | leaks, not ripple (an eight-pole cusp's gaps are narrow: 2.6×) |
| `--frozen` | correlation > 0.99 | two advections of one mass flux, differing only in numerical diffusion at 32 cells a wavelength |
| `--quench` front | 2Δx / Δt | the front found to a cell at each end |
| `--quench` decay | 10⁻⁶ of the field | the plugin's product of per-frame decays is exp(−Σ dt / τ_q) exactly |
| `--resist` | 1% | the 5-point Laplacian on a Gaussian 20 cells wide is (Δx/w)²/12; the compression is 10⁻⁵ |
| `--floors` | < 10⁻⁴ of cell-steps | "rarely" made a number; the default look measures 7×10⁻⁶ |
| `--still` | bit for bit | `texelFetch` of the input at Mix 0 |
| `--glow` grid | 6 × 233 u | sums of up to 233 products, two passes a stage, three stages |
| `--glow` raster | 10⁻⁴ | at a raster exactly twice the grid, bilinear resampling keeps a sum |

**The Brio–Wu reference.** The published solution exists as plots, not a
table, and there is no network here. The reference is computed by the
harness: 1-D ideal MHD in double precision, MUSCL–Hancock with minmod and the
**Rusanov** flux (a different Riemann solver from the plugin's), at 8192
cells. Its two fast-rarefaction heads are checked against the closed-form
fast speeds of the initial states (0.3195 vs 0.3208, 0.8686 vs 0.8684). Its
plateaus, ρ = 0.696 and 0.235, are the familiar published values.

## Decisions taken without asking

- **HLLD**, float, no fallback needed.
- **2.5-D field lines are in-plane only.** A guide-field-only bottle has none,
  correctly.
- **Wall is the default boundary**, and the coils' changes are carried
  through the vessel volumetrically. This goes beyond the spec, which asks
  for "re-imposed at the boundary"; Open does exactly that.
- **The background is ρ = 0.1, p = 0.005.** It is tenuous enough to emit
  ~0.2% of the ball's light, and dense enough that v_A (and the step) stays
  sane. The floors are ρ 10⁻⁴ and p 10⁻⁶.
- **Temperature is β against B_ref, not against Field**, so Field confines
  the same ball harder rather than rescaling everything.
- **Speed is Alfvén crossings (at B_ref) per second**, default 0.3.
- **The Green Orb is a preset file** (`docs/green-orb.preset`), not an option
  parameter. The defaults are the same bottle with the clip's own colours.
  Aurora is named for the oxygen line, and nothing refers to the inspiration.
- **Events land on frames**; Quench is a boolean (latching) and the coils
  come back when it is released.
- **Provisional About/ATTRIBUTIONS** hand copies with `guide = ""`, as in the
  rest of the unreleased tranche.
- **The display name is `SW Containment`** (14 characters), per the spec.

## What is genuinely verified, and what is assumed

Verified on this machine: every check above, all 13 negative controls
detected, the mutation detected by two checks, the sweep (31 parameters
live), `verify.sh` green on a fresh universal build (`lipo` shows both
slices, `oxbow probe` reads `SW Containment` / `CT01` / effect).

Not verified, or not done:

- **Never loaded into Resolume.** Never built or run on Windows, and never
  run on another GPU. The x86_64 slice has never executed.
- **Open is not indefinitely stable.** The default look under Open runs away
  after about 3 τ_A (about 10 s at the default Speed): field builds at the
  edges and the floors fire. Short events (a quench fireball, a few seconds
  of cusp jets) are fine. Wall is the default for this reason.
- **The dual-energy switch and EGLM are not conservative.** With a cusp, the
  total energy drifts 2×10⁻³ over 600 frames of Wall. The guide-field bottle,
  where neither engages, conserves to 3×10⁻⁷.
- **The GLM energy term, the carried coil change and the T-scaled curvature
  are modelling choices** stated above, not textbook ideal MHD.
- **Resolume's FFT bins** are assumed to be what rosette assumed.
- **The GPU-less CI runner** would need a software GL; `verify.sh` runs a few
  minutes here.
- **Stretch goals not attempted**: the tokamak view and the "Over"
  registration. OpenFX port and browser demo: not required for 0.1.0.
