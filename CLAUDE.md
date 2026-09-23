# containment

A ball of hot plasma held in a magnetic bottle, for Resolume Arena/Avenue, as
an FFGL effect. The clip is what the plasma is made of, and the plasma obeys
2.5-D compressible ideal MHD on the GPU. C++/GLSL, CMake MODULE → universal
`.bundle` (macOS) + Windows `.dll`. MIT. ID `CT01`, display name
`SW Containment`, bundle id `com.stoatworks.ffgl.containment`.

Read `AGENTS.md` before touching the solver, the boundaries, the clock, the
floors or the dual-energy switch.

## Commands (CMake)
- Configure: `cmake -B build -DCMAKE_BUILD_TYPE=Release`
- Fast dev build: add `-DCMAKE_OSX_ARCHITECTURES=arm64`
- Build: `cmake --build build`
- Render offline: `./build/cttest --out /tmp/frame.png --frames 300`
- The green look: `./build/cttest --preset docs/green-orb.preset --out /tmp/orb.png --frames 400`
- List parameters: `./build/cttest --list`
- Set anything by name: `./build/cttest --set "Poles=1" --set "Curvature=0.5"`
- The state's ranges (development): `./build/cttest --stats --frames 600`
- Every program as compiled: `./build/cttest --dump-shaders DIR`
- Film: `./build/cttest --film 780 --size 1280x720 --preset docs/green-orb.preset --script docs/demo.cues | ffmpeg -f rawvideo -pix_fmt rgba -s 1280x720 -r 60 -i - out.mp4`
- Film a clip through it: `ffmpeg -i in.mov -f rawvideo -pix_fmt rgba - | ./build/cttest --pipe --size WxH [--script cues] | ffmpeg …`
- Do NOT `cmake --install` (it writes into Arena's Extra Effects).

## Verify
- Everything: `tools/verify.sh` (fresh universal build, dumped-shader compile,
  every check, the mutation, the negative controls, the sweep, lipo, plist,
  ad-hoc signature, `oxbow probe`, the bench)
- The solver: `--briowu`, `--alfven`, `--conserve`, `--divb`
- The bottle: `--balance`, `--rt`, `--cusp`, `--frozen`, `--quench`, `--resist`
- The model's guards: `--floors`
- The light and the host: `--still`, `--glow`, `--state`
- The harness drives the shipped shader: `--mutation`
- Every check against a wrong model: `--negative` (`CT_NEGATIVE=name` for one)
- No dead controls: `python3 tools/sweep.py`
- Cost: `./build/cttest --bench`

## Notes
- **The state is four RGBA32F textures** (A ρ,m; B E,B; C ψ, signal speed,
  ρK, flags; D ρ·colour, ρχ), ping-ponged. Never half floats.
- **A substep is: reduce → clock → predict → flux x → flux y → update.** dt
  lives on the GPU; the CPU plans the count from last frame's speed. Short
  plans make time run slow, never unstable (logged via Diag).
- **Each face's flux is computed once** in its own pass and read by both
  cells; do not fold it back into the update (latency and exact telescoping).
- **c_h is half the summed signal speed**; at the full sum GLM blows up.
- **GLM's B change does no work on the gas** (EGLM energy term); without it
  low-β pressure goes negative.
- **Dual-energy switch**: p from the entropy where p < 2% of kinetic+magnetic.
- **Wall = no-slip, line-tied, ghost field = coils + mirrored perturbation.**
  Coil changes are carried through the vessel volumetrically (Sources pass).
- **Open ghost never falls below the ambient background.** Open is stable for
  ~3 τ_A of the default look, not indefinitely (see AGENTS.md).
- **Curvature g_eff scales with T / T_ball**; the drive's stream function is
  windowed about the ball.
- **GLSL `tanh` returns NaN past |x| ≈ 44 here**; clamp its argument.
- **Every declared sampler needs a real texture bound**, even unused ones.
- Script tracks interpolate: a press is three keys.
- All host parameters are 0..1 and mapped in `Controls.cpp`; option
  parameters hold the element value; events act on the rising edge.
- `containment_core` is an OBJECT library; override `SetTextParameter`.
- Local repo only: no remote, no tag, not registered on the website.

## Not done yet
- Never loaded into Resolume (`oxbow probe` only). Never built on Windows. No
  OFX port, no browser demo, no user guide. `StoatworksAbout.h` and
  `ATTRIBUTIONS.md` are provisional hand copies.

## Diagnostics

`source/Diag.{h,cpp}` is a log file only, with no crash handler. It records
which shader failed to compile, the GL strings, the host clock's unit, and
when the substep cap bit.

    ~/Library/Logs/containment/containment.YYYY-MM-DD.log
