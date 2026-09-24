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
  (raw RGBA in and out; a partial frame at EOF ends the stream; a cue naming
  no parameter exits 2; a reader that hangs up gives exit 1, SIGPIPE ignored)
- Do NOT `cmake --install` (it writes into Arena's Extra Effects).

## Verify
- Everything: `tools/verify.sh` (the FFGL pin, a fresh universal build, the
  shaders, `--offline`, every GL check at its own raster AND at 320x180, the
  mutation, the negative controls, the pipe, the sweep, lipo, plist, ad-hoc
  signature, `oxbow probe`, the bench). It runs every step and counts.
- No GL (what CI runs): `./build/cttest --offline` (`--names`, `--presets`,
  `--reference`, `--vacuum` and their negative controls)
- Any check at another raster: `--size 320x180` (the grid keeps the check's
  aspect; only the light's path to the raster moves)
- On another rasteriser: `CT_RENDERER=software` (Apple's software renderer,
  ~50 s a frame -- `--state` and `--alfven` are what is practical)
- Shaders alone: `tools/check-shaders.sh build/cttest`
- The solver: `--briowu`, `--alfven`, `--conserve`, `--divb`
- The bottle: `--balance`, `--rt`, `--cusp`, `--frozen`, `--quench`, `--resist`
- The model's guards: `--floors`, `--open`
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
- **Open is a window onto a bigger bottle**: a margin 0.1 frame heights deep,
  simulated and never shown, holds an absorbing layer that relaxes the plasma
  towards the ambient at rest in the coils' field. Wall is the default (Open
  costs 13–16 ms at 1080p against Wall's 8). The
  layer's target is the plugin's ambient, so a check that loads any other
  plasma must use Wall (see AGENTS.md). Its ghost never falls below ambient.
- **Preset** is index 0 (the fleet's override model, `Presets.h`); row 1 is
  the constructor's defaults and `--presets` holds them together. **Fuel**
  tops the ball's footprint back up.
- **Curvature g_eff scales with T / T_ball**; the drive's stream function is
  windowed about the ball.
- **GLSL `tanh` returns NaN past |x| ≈ 44 here**; clamp its argument.
- **Every declared sampler needs a real texture bound**, even unused ones.
- Script tracks interpolate: a press is three keys.
- All host parameters are 0..1 and mapped in `Controls.cpp`; option
  parameters hold the element value; events act on the rising edge.
- `containment_core` is an OBJECT library; override `SetTextParameter`.
- The display name and id are `kDisplayName` / `kPluginCode` in `Controls.h`.
- Released: public `stoatworks-labs/containment`, v0.1.0 (2026-09-24), registered on
  the website with a user guide (`docs/USER-GUIDE.md` is the source; the PDF is generated
  by the website's `build_guides.py`).

## Not done yet
- Never loaded into Resolume on macOS (`oxbow probe` only). On Windows the
  fleet's Arena gate passes 9/9 on llvmpipe (`plugin-bench/arena/expect/containment.json`).
  No OFX port. `StoatworksAbout.h` and `ATTRIBUTIONS.md` are generated by
  stoatworks-backend; do not hand-edit them.

## Browser demo
- `demo/` is containment-demo.stoatworks-labs.com: a static-assets Worker
  (`wrangler.toml`, a Worker ROUTE on a proxied AAAA 100:: record, because the
  zone is at its 100 custom-domain limit), no build step. `deploy.yml` redeploys
  it on every push to main; by hand: `cf-run npx wrangler deploy`.
- `demo/plugin.js` carries every shader piece unedited plus `PROGRAM_PIECES`
  (SourceFor's joins). **Change a shader, copy it across**:
  `python3 demo/tools/check_shaders.py` (verify step "demo") fails otherwise.
- The CPU half (Controls, MakeCoils, ChooseGrid, Pcg, Drive, Presets/P(),
  ProcessOpenGL's sequence) is a hand port: change it there too.
- `demo/vendor/` is the shared kit: never edit it; re-vendor with
  `stoatworks-backend/resolume-demo/sync.sh containment`.
- The page defaults to Detail 128 (`DEMO_DETAIL`), not 256; said on the page.

## Diagnostics

`source/Diag.{h,cpp}` is a log file only, with no crash handler. It records
which shader failed to compile, the GL strings, the host clock's unit, and
when the substep cap bit.

    ~/Library/Logs/containment/containment.YYYY-MM-DD.log
