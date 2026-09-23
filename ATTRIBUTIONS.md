# Attributions

containment is built on other people's work. This file lists what that work is,
who did it, and what it is doing here.

> **Provisional.** Across the fleet this file is generated from master lists in
> `stoatworks-backend` by `scripts/sync-attributions.py`. containment is not
> registered there yet, so this copy is hand-written. Register it before release
> — and note that the script's `--only` flag truncates the file rather than
> filtering it.

## Third-party code this project uses

### Resolume FFGL SDK

<https://github.com/resolume/ffgl>
Licence: BSD-3-Clause
Copyright: FreeFrame

Vendored as a git submodule at `external/ffgl`, pinned to `b1afaf9`.

The plugin ABI itself. An FFGL effect is defined by this SDK's headers — there
is no other way to be loadable by Resolume Arena and Avenue.

### GLEW — the OpenGL Extension Wrangler Library

<https://github.com/nigels-com/glew>
Licence: BSD-3-Clause (with Mesa 3-D and Khronos components)
Copyright: Milan Ikits, Marcelo E. Magallon and Lev Povalahev

Windows only, from vcpkg, statically linked. The SDK's headers pull it in for
the OpenGL function pointers; macOS uses the system OpenGL framework instead.

### zlib

<https://zlib.net>
Licence: zlib
Copyright: Jean-loup Gailly and Mark Adler

Ships with macOS. The offline harness links it to deflate its PNG output.
Nothing in the shipped plugin uses it.

## Work from elsewhere in the fleet

### rosette, millpond — the audio analyser and the clock

<https://github.com/stoatworks-labs/rosette>
Licence: MIT
Copyright: Stoatworks Labs

`source/Audio.{h,cpp}` is rosette's analyser (from macroblock's), with its
primed first frame, by way of millpond. The host-clock unit vote in
`UpdateClock` is rosette's, unchanged.

### vectrix, tinsel, millpond

<https://github.com/stoatworks-labs/vectrix>
Licence: MIT
Copyright: Stoatworks Labs

`GLState.h` (put the host's state back however the frame ends) is vectrix's,
which took it from resolume-scopes. `PassBuffer` (`FFGLFBO` with the SDK's
colour-texture leak fixed), `Diag`, the CMake shape, the harness shape, the
`--pipe`/`--film`/`--script` format, `tools/sweep.py` and `tools/verify.sh`
come from tinsel and millpond.

## Method

Textbook numerics and physics, described in books and papers rather than
copied from anyone's source:

- Ideal MHD in conservative form; MUSCL–Hancock with a monotonised-central
  limiter — Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*.
- The HLLD Riemann solver — T. Miyoshi and K. Kusano, *J. Comput. Phys.* 208,
  315 (2005).
- Hyperbolic–parabolic divergence cleaning (GLM) and its EGLM energy source —
  A. Dedner et al., *J. Comput. Phys.* 175, 645 (2002); the damping as in
  A. Mignone and P. Tzeferacos, *J. Comput. Phys.* 229, 2117 (2010).
- The dual-energy (entropy) switch — D. Ryu et al., *ApJ* 414, 1 (1993);
  D. Balsara and D. Spicer, *J. Comput. Phys.* 148, 133 (1999).
- The Brio–Wu shock tube — M. Brio and C. C. Wu, *J. Comput. Phys.* 75, 400
  (1988); the circularly polarised Alfvén wave — G. Tóth, *J. Comput. Phys.*
  161, 605 (2000); the 1-2-3 problem — B. Einfeldt et al., *J. Comput. Phys.*
  92, 273 (1991); the exact Euler Riemann solver in the harness — Toro, ch. 4.
- The magnetic Rayleigh–Taylor growth rate — Chandrasekhar, *Hydrodynamic and
  Hydromagnetic Stability*, ch. X.
- Bremsstrahlung emissivity ∝ n²√T of an optically thin plasma.

## Getting this wrong

If your work is here and the description is inaccurate, the licence is wrong, or
you would rather not be listed — open an issue and it will be fixed.
