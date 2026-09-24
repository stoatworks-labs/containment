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
  R_c ≥ 1.05× the grid's half-diagonal (Open's margin included, so no coil
  is ever inside the simulated plasma), normalised so |B| at radius 0.5 in
  the first gap is Field. Plus a uniform guide Bz. The closed form is written
  in `Physics.cpp` and in the shader, both marked `//= mirrored`.
- **Open**: a window onto a bigger bottle. The grid carries a margin all
  round the frame, 0.1 frame heights deep (`kMarginFraction`), simulated and
  never shown. In it the plasma is relaxed towards the ambient plasma at rest
  in the coils' own field, at a rate rising as depth² (`kSpongePower`), sized
  so the fastest wave (c_h) loses 6 e-foldings crossing it (`kSpongeEFolds`).
  The relaxation is exact for the step, so no rate is too stiff. Beyond the
  margin the ghost is zero-gradient but never below the ambient, its field
  B_vac(ghost) + (B − B_vac)(the cell it copies); ψ leaves. The coils' slow
  changes are carried through the volume, as under Wall. The ball's X / Y and
  the clip are the frame's; the margin reads the clip's edge.
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

**Default: Wall, guide field 1.4 × Field, Field 0.5, six poles, Fuel 0.3.** The guide
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
as they should. On the default look 56% of cell-steps take the entropy branch
(the background), and the floors then fire on none (7×10⁻⁶ before the 2026-09-24 build).

**The open boundary drains the background.** A zero-gradient ghost lets the
tenuous plasma stream out along field lines that cross the edge, leaving
β ≈ 0 holes that only the floors held up (millions of clamps, then a
runaway). The ghost then held no less than the ambient background, and Open
still ran away after about 3 τ_A as |B| climbed at the edges (28.6× the
coils' own field by 20 τ_A, `--negative`'s "open"). A pure-vacuum field ghost
was worse: NaN by t = 4.5. What holds is the absorbing margin (above): 20 τ_A
of the default look with |B| never past the coils' own field, and no floors.

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

**A released ball fingers, and the marker's half-contour is not its edge.**
This is what turned `--balance` red at Detail 512 (4.5 cells against a
bound of 3), and it was the check, not the physics. The top-hat ball is let
go out of balance and rings; while its edge decelerates it is
Rayleigh–Taylor unstable (a dense ball, a tenuous background, k across B, so
no tension holds it), seeded by the Cartesian grid on the axes. At 512 the
edge is resolved well enough to finger: the second moment of the marker's
> 0.5 region, over a disc's, grew steadily from 1.000 to 1.117 over 8 τ_A. In the
mixed layer a cell half ball by MASS is only ρ_out / (ρ_in + ρ_out) ≈ 13%
ball by volume, so the half-contour counts the layer and runs outwards. The
absorbing margin made it worse, correctly: the old ghost sent the first fast
wave back within 0.2 τ_A and cushioned the overshoot (the edge reached
0.1878 instead of 0.1944). The edge is now the ball's VOLUME, M / ρ_core: at
pressure balance the ball's plasma, in the core and in any finger, sits on
the core's adiabat with the core's Bz / ρ, so it has the core's density
everywhere. That lands 0.03 and 0.09 cells from flux conservation. The core
also heats by 2–3% in entropy during the ring-down (the compressions converge
on the axis), so an adiabatic prediction of the edge misses by 0.2–0.5
cells, past what the ringing allows it; the entropy is printed, the
adiabatic edge is not asserted.

**The absorbing layer relaxes towards the plugin's ambient plasma**, so any
check that loads a different plasma under Open gets a Riemann problem at the
frame's edge. `--resist` (ρ 100, p 1000) read a peak of −0.996 against 0.025;
it runs in a Wall now.

**|B| over the grid swings 2× as the coils turn past its corners.** The
margin puts the grid's corners within a few hundredths of the coil circle.
`--open` compared |B| with its value at the start and failed on the coils'
own field (5.4 against 2.6, periodically). It now compares with the coils'
own peak at the same moment: 1.00 throughout.

**A plane pulse is the wrong probe of an absorbing layer.** It runs ALONG the
top and bottom margins too, the layer damps it there as it should, and the
step that leaves diffracts into the frame (2.4% of it, sitting by the top and
bottom edges with a transverse velocity). Nothing the ball sends out runs
along an edge from outside the frame. The probe is a cylindrical wave from
the centre, differenced cell for cell against the same frame inside an
unbounded plasma.

**The layer echoes a low-frequency oblique wave by ~12%, and so did the
zero-gradient ghost.** A relaxation layer does not reflect a wave meeting it
square on (every characteristic is damped alike); at a slant its impedance
changes and it does. Depth barely helped (0.2 frame heights: 4–5%; 0.3:
2–4%) and costs cells. What differs is what happens next: the old ghost kept
9.3% of the wave in the frame two crossings later, the layer 0.74%.

**A bilinear resample conserves a sum only at integer upsampling.** `--glow`'s
raster half held to 10⁻⁴ at 512² on a 256² grid, and read 0.999874 at
320x180, where the grid is sampled DOWN and whole cells are skipped. Its
bound is now the resampler's own ripple where the light is, computed on the
CPU for the two sizes, with the grid no bigger than the raster.

**Apple's software renderer is ~50 s a frame at the default Detail.**
`CT_RENDERER=software` forces it here; `--state` took 2.5 minutes and the
Alfvén wave at Detail 128 read exactly what the GPU reads (0.70762, 0.054%).
The physics and the sweep cannot run in CI on it.

**`--dump-shaders` into a missing directory said "wrote 14 programs".**
`std::ofstream` fails silently; it now checks and exits 1.

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
| `--divb` | max 0.1, rms 0.01 of \|B\|/Δx | O(1) is a monopole a cell across; GLM holds it to truncation level. Measured max 0.028 / 0.024 at Detail 128 / 256; with GLM off 0.20 / 0.15 |
| `--balance` pressure, Bz | max ρv² / p_total + 0.5% | residual ringing is balanced by ρ dv/dt ~ ρv²/L; half a percent for the edge's own width |
| `--balance` edge | 1 cell | the ball's volume radius, √(M / π ρ_core), against flux conservation, r₀ √(Bz₀ / Bz_core). A fast wave across B moves Bz and ρ in proportion, so the ringing leaves Bz/ρ alone; the only plasma the core does not describe is what the ignition's tanh edge (half-width 0.75 cells) laid down part-mixed, less than a cell's width of the ball |
| `--balance` β = 1 | inside the edge layer ± 1 cell | β = 1 lies where 2p = B², inside the layer between marker 0.9 and 0.1 |
| `--rt` rate | (kL/2) gkA/γ² + g/(k c_s²) + 2%, L = 2 cells | a diffuse interface weakens only the buoyancy term; plus compressibility and the fit |
| `--rt` stable | peak < 2× the seed | the unstable cases reach > 20× in the same time |
| `--cusp` angles | half a bin + 2Δx/r | the histogram's bin and the grid at the sampling circle |
| `--cusp` contrast | 1.5× | leaks, not ripple (an eight-pole cusp's gaps are narrow: 2.6×) |
| `--frozen` | correlation > 0.99 | two advections of one mass flux, differing only in numerical diffusion at 32 cells a wavelength |
| `--quench` front | 2Δx / Δt | the front found to a cell at each end |
| `--quench` decay | 10⁻⁶ of the field | the plugin's product of per-frame decays is exp(−Σ dt / τ_q) exactly |
| `--resist` | 1% | the 5-point Laplacian on a Gaussian 20 cells wide is (Δx/w)²/12; the compression is 10⁻⁵ |
| `--floors` | < 10⁻⁴ of cell-steps | "rarely" made a number; the default look measures 0 |
| `--open` echo | what remains two frame crossings after the first echo left < (first echo)² | each further echo is at most the first's fraction of the one before; one factor of margin. The first echo is printed (11.6%), not bounded |
| `--open` long run | \|B\| < 2× the coils' own peak at that moment; floors < 10⁻⁴; ρ > 0.1 ρ_ambient; emission > ¼ of settled; 20 τ_A | a runaway passes 2× in a few τ_A (the old ghost: 28.6×); the rest are `--floors`' bound and "the ball is still there" |
| `--still` | bit for bit | `texelFetch` of the input at Mix 0 |
| `--glow` grid | 6 × 233 u | sums of up to 233 products, two passes a stage, three stages |
| `--glow` raster | the resampler's ripple + 20 u | for each texture, the pixels' summed bilinear weight on each cell against its mean, plus 2⁻⁹ a pixel where the filter's fraction is not a multiple of 1/256 (GL leaves sub-texel precision to the implementation), weighted by where the light is; computed on the CPU for the two sizes. 1.2×10⁻⁶ at 512² (exactly twice the grid), 0.35 at 320x180 -- still far inside the 0.6 an additive bloom adds |
| `--reference` (offline) | 4 cells of 8192 + 0.002 | a head found at 0.1% of the jump; a rarefaction head is a kink that minmod smears over √(c Δx t / 2) = 0.0033 |
| `--vacuum` (offline) | div B, J_z < Σ_k \|I_k\| (2h²/(d_k−h)⁴ + 2N u/(h d_k)); growth r^(N/2−1) to 4N(0.1/R_c)^N; \|B\|(0.5) = Field to 10⁻¹²; N cusps | the central difference's truncation (a line current's third derivative is 6\|I\|/d⁴) and rounding; the coils' next multipole term |
| `--names` (offline) | exact | unique parameter names; the display name ≤ 16 bytes with the fleet's prefix |

**The Brio–Wu reference.** The published solution exists as plots, not a
table, and there is no network here. The reference is computed by the
harness: 1-D ideal MHD in double precision, MUSCL–Hancock with minmod and the
**Rusanov** flux (a different Riemann solver from the plugin's), at 8192
cells. Its two fast-rarefaction heads are checked against the closed-form
fast speeds of the initial states (0.3195 vs 0.3208, 0.8686 vs 0.8684). Its
plateaus, ρ = 0.696 and 0.235, are the familiar published values.

## Would this hold on another rasteriser, at another raster?

`tools/verify.sh` runs every GL check twice: at the raster the check asks for,
and with every rig at 320x180, CI's (`--size 320x180`). Under `--size` the
grid keeps the aspect the check asked for, so a physics check runs on the
same cells at both and only the light's path to the raster changes; the two
passes print the same physics numbers. `CT_RENDERER=software` forces Apple's
software renderer, a different rasteriser, on this Mac.

| check | raster dependence | another rasteriser |
| --- | --- | --- |
| `--briowu` | none: the grid is Detail cells, the tube 1-D along a square box's axis; same numbers at 320x180 | texelFetch and float arithmetic only; no filtered sample in the solver |
| `--alfven` | none; same at 320x180 | **run on the software renderer**: Detail 128 read 0.70762 and lost 0.054%, digit for digit the GPU's |
| `--conserve` | none; same at 320x180 | sums of state texels read back in double; the bound is float rounding per substep, any IEEE float32 target |
| `--divb` | none | as `--conserve` |
| `--balance` | none; same at 320x180 | radial means of state texels; no raster enters |
| `--rt` | none | the fit is on state texels |
| `--cusp` | none | the histogram is of state texels on a circle, binned on the CPU |
| `--frozen` | none | correlation of state texels |
| `--quench` | none | the front from state texels; the coils' decay is CPU double |
| `--resist` | none | one state texel against a closed form |
| `--floors` | the default look at 320x180 either way | counts from the GPU's own flags; a rasteriser with other float slack could move the count, not past 10⁻⁴ |
| `--open` | none (the probe's grid is its own) | state texels differenced between two runs on the same GPU: any systematic of the rasteriser cancels |
| `--still` | 480x270 and 1280x720, and 320x180: bit for bit at all three | `texelFetch` of the input and a return: no arithmetic on the path, exact on any GPU |
| `--glow` grid | the emission grid's, not the raster's | normalised kernels in float: 6 × 233 u covers any IEEE float32 target |
| `--glow` raster | 512² (1.2×10⁻⁶) and 320x180 (0.35): the bound is computed per size | the filter's sub-texel precision is in the bound (2⁻⁹ per pixel off a 1/256 step) |
| `--state` | 320x180 both passes | **run on the software renderer**: passes; state, not pixels |
| `--mutation` | runs `--briowu` and `--alfven` | as those |
| `--reference`, `--vacuum`, `--names`, `--presets` | CPU only (`--offline`) | no GPU at all |

Deliberately not relied on: exact cancellation (`--glow`'s bound carries
rounding even where the resample is exact), `pow( 1.0, 1.0 ) == 1.0`, and
`tanh` past |x| ≈ 44 (clamped; see the traps).

## The recorded mutations

Two, both proving the harness drives the shaders the plugin compiles.

- **In the harness, every run** (`--mutation`): the monotonised-central
  slope's `0.5 * abs( l + r )` becomes `0.5 * abs( l - r )`, through the
  plugin's own shader assembly (`SetShaderMutationForTest`). `--briowu` fails
  10 assertions (waves 12–18 cells off, the contact 16 cells wide) and
  `--alfven` 1 (the amplitude converges at order 0.96). The check asserts the
  text is in 3 shipped programs, so a rename cannot make it vacuous.
- **By hand, once, 2026-09-24**: in `Shaders.cpp`'s `fastSpeed()`,
  `float a2  = Gamma * p / r;` became `Gamma * p * r` (one character), and
  the build ran the checks. `--quench` failed 2 (the front never left: 0.0000
  against the exact shock's 3.2154, and the ball's area went to 0). `--briowu`
  and `--alfven` PASSED: a wrong sound speed only widens HLLD's wave-speed
  estimates, which makes the scheme more diffusive but still consistent, and
  Brio–Wu's positions sit inside their 3 cells. Reverted by copying the file
  back; `git diff` empty; `--briowu` passes again.

## Decisions taken without asking

- **HLLD**, float, no fallback needed.
- **2.5-D field lines are in-plane only.** A guide-field-only bottle has none,
  correctly.
- **Wall stays the default boundary; Open is sound but dearer** (decided
  2026-09-24). The unfinished work had made Open the default. With the
  absorbing margin Open now holds the default look for 20 τ_A with the field
  never past the coils' own, but the default look under Open costs 13–16 ms
  a frame at 1080p against Wall's 8: the margin adds a third more cells, and
  its corners sit near the coils, where the field sets a shorter step (29
  substeps a frame against 22). The spec asks the default to hold 60 fps at
  1080p with headroom, and 83–95% of the frame is not headroom. So Clip Orb,
  Green Orb, Guide Field Bubble and Rayleigh-Taylor are Wall; Cusp Leak and
  Quench Fireball, where what leaves should leave, are Open. Under both the
  coils' changes are carried through the volume, which goes beyond the spec's
  "re-imposed at the boundary".
- **The margin is 0.1 frame heights, 6 e-foldings, rate ∝ depth².** Deeper
  margins echo less (0.2: ~5%, 0.3: ~3%) and cost cells on every frame; what
  the release needs is that the echo leaves, and it does at 0.1.
- **`--balance`'s edge is the ball's volume, not the marker's half-contour**,
  and its tolerance is one cell, derived. The half-contour is printed beside
  it so the mixed layer stays visible.
- **The Preset dropdown** (index 0, the fleet's override model; `Presets.h`):
  Custom, then Clip Orb (= the constructor's defaults, held by `--presets`),
  Green Orb, Guide Field Bubble, Cusp Leak, Rayleigh-Taylor, Quench Fireball.
  Choosing a row re-ignites.
- **Fuel** (0–2 per τ_A, default 0.3): a steady gas puff that tops the ball's
  footprint back up towards its ignition profile. It only adds, at rest, so it
  holds a leaking ball up without pinning it. It has no physics check of its
  own (it is a source, like Feed); the sweep holds it live, and `--open`'s
  long run is of the fuelled default.
- **CI runs what a GPU-less runner can**: `--offline`, glslc, and `--state` at
  320x180 under `--allow-no-gl`. The physics and the sweep run in `verify.sh`.
- **The background is ρ = 0.1, p = 0.005.** It is tenuous enough to emit
  ~0.2% of the ball's light, and dense enough that v_A (and the step) stays
  sane. The floors are ρ 10⁻⁴ and p 10⁻⁶.
- **Temperature is β against B_ref, not against Field**, so Field confines
  the same ball harder rather than rescaling everything.
- **Speed is Alfvén crossings (at B_ref) per second**, default 0.3.
- **The Green Orb is a Preset row and a preset file** (`docs/green-orb.preset`,
  for `cttest --preset`). The defaults are the same bottle with the clip's own
  colours.
  Aurora is named for the oxygen line, and nothing refers to the inspiration.
- **Events land on frames**; Quench is a boolean (latching) and the coils
  come back when it is released.
- **Provisional About/ATTRIBUTIONS** hand copies with `guide = ""`, as in the
  rest of the unreleased tranche.
- **The display name is `SW Containment`** (14 characters), per the spec,
  written once (`kDisplayName` in `Controls.h`) so `--names` can hold it to 16
  bytes; `oxbow probe` reads the bundle's copy.

## What is genuinely verified, and what is assumed

Verified on this machine (2026-09-24, M4 Max, macOS 26.4): every check above,
at each check's own raster and at 320x180; all 18 negative controls detected
(3 of them offline); the in-harness mutation caught by two checks and the
hand mutation by one; the sweep (33 parameters live); the pipe; `verify.sh` green (47 steps) on a fresh
universal build (`lipo` shows both slices, `oxbow probe` reads
`SW Containment` / `CT01` / effect). On Apple's software renderer: `--state`
and `--alfven` only.

Not verified, or not done:

- **Never loaded into Resolume.** Never built or run on Windows, and never
  run on another GPU. The x86_64 slice has never executed.
- **Open echoes ~12% of a low-frequency wave once**, as the zero-gradient
  ghost did; the echo then leaves. The long run is 20 τ_A (~70 s at the
  default Speed), not an evening.
- **The dual-energy switch and EGLM are not conservative.** With a cusp, the
  total energy drifts 2×10⁻³ over 600 frames of Wall. The guide-field bottle,
  where neither engages, conserves to 3×10⁻⁷.
- **The GLM energy term, the carried coil change, the T-scaled curvature, the
  absorbing margin and Fuel are modelling choices** stated above, not
  textbook ideal MHD.
- **Resolume's FFT bins** are assumed to be what rosette assumed.
- **CI has never run** (no remote). Its steps were run here by hand.
- **Stretch goals not attempted**: the tokamak view and the "Over"
  registration. OpenFX port, browser demo and user guide: not required for
  0.1.0.
