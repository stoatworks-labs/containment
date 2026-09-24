# Attributions

Containment is built on other people's work. This file lists what that work is, who did
it, and what it is doing here.

It is generated — the master lists live in the `stoatworks-backend` repo and are
pushed out by `scripts/sync-attributions.py`. Edit it there, not here.

## Code we derived from other people's work

Someone else solved this first, and this project would not exist in its current form without their work.

### Audio analyser and host-clock vote — Stoatworks rosette

<https://github.com/stoatworks-labs/rosette>  
Licence: MIT  
Copyright: Stoatworks Labs

source/Audio.{h,cpp} is rosette's analyser (itself from macroblock's), with its primed first frame, by way of millpond. The host-clock unit vote in UpdateClock is rosette's, unchanged.

### GL state guard, pass buffers, diagnostics and harness — Stoatworks vectrix

<https://github.com/stoatworks-labs/vectrix>  
Licence: MIT  
Copyright: Stoatworks Labs

GLState.h (put the host's state back however the frame ends) is vectrix's, which took it from resolume-scopes. PassBuffer (FFGLFBO with the SDK's colour-texture leak fixed), Diag, the CMake shape, the harness shape, the --pipe/--film/--script format, tools/sweep.py and tools/verify.sh come from tinsel and millpond; the override preset model is graticule's.

## Third-party code this project uses

Libraries, SDKs and frameworks the project is built on or bundles.

### Resolume FFGL SDK

<https://github.com/resolume/ffgl>  
Licence: BSD-3-Clause  
Copyright: FreeFrame

Vendored as a git submodule at external/ffgl (third_party/ffgl in oxbow).

The plugin ABI itself. An FFGL effect or source is defined by this SDK's headers — there is no other way to be loadable by Resolume Arena and Avenue.

### GLEW — the OpenGL Extension Wrangler Library

<https://github.com/nigels-com/glew>  
Licence: BSD-3-Clause (with Mesa 3-D and Khronos components)  
Copyright: Milan Ikits, Marcelo E. Magallon and Lev Povalahev

Arrives inside the FFGL submodule at external/ffgl/deps/glew-2.1.0. Not fetched separately.

Resolves OpenGL entry points on Windows, where the system headers stop at OpenGL 1.1.

### libpng

<http://www.libpng.org/pub/png/libpng.html>  
Licence: PNG Reference Library License (libpng)  
Copyright: the PNG Reference Library authors

Arrives inside the FFGL submodule, under the SDK's CustomThumbnail sample.

Part of the upstream SDK tree rather than something these plugins call directly — listed because it is present in the checkout.

## Inspirations

What this set out to be. No code, assets or binaries from any of these were used or examined — the debt is to the idea.

### Ideal MHD numerics

MUSCL-Hancock with a monotonised-central limiter (Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics); the HLLD Riemann solver, T. Miyoshi and K. Kusano, J. Comput. Phys. 208, 315 (2005); GLM divergence cleaning and its EGLM energy source, A. Dedner et al., J. Comput. Phys. 175, 645 (2002), damped as in A. Mignone and P. Tzeferacos, J. Comput. Phys. 229, 2117 (2010); the dual-energy switch, D. Ryu et al., ApJ 414, 1 (1993) and D. Balsara and D. Spicer, J. Comput. Phys. 148, 133 (1999). Implemented from the papers; nothing is copied from anyone's source.

## Standards and published specifications

What the implementation is measured against.

- **M. Brio and C. C. Wu, J. Comput. Phys. 75, 400 (1988); G. Toth, J. Comput. Phys. 161, 605 (2000); B. Einfeldt et al., J. Comput. Phys. 92, 273 (1991)** — The Brio-Wu shock tube, the circularly polarised Alfven wave and the 1-2-3 problem the harness checks against; the exact Euler Riemann solver in the harness is Toro's, chapter 4.
- **S. Chandrasekhar, Hydrodynamic and Hydromagnetic Stability, ch. X** — The magnetic Rayleigh-Taylor growth rate the harness measures the fingers against.
- **Bremsstrahlung emissivity of an optically thin plasma** — The light: rho^2 sqrt(T) times the plasma's colour.

## Getting this wrong

If your work is here and the description is inaccurate, the licence is wrong, or you would rather not be listed — open an issue and it will be fixed.
