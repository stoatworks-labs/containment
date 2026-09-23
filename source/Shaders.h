#pragma once

/**
    The passes, as GLSL source.

    ----------------------------------------------------------------- the state

    Four RGBA32F textures on the grid, ping-ponged:

        A  ( rho, rho u, rho v, rho w )
        B  ( E, Bx, By, Bz )
        C  ( psi, sum of signal speeds, max signal speed, floor flag )
        D  ( rho r, rho g, rho b, rho chi )       -- the picture, and the ball marker

    C's last three channels are not state: the update writes each cell's
    signal speeds (for the CFL step and GLM's c_h) and whether a floor fired,
    so the reduction that follows reads one texture rather than recomputing.

    ----------------------------------------------------------------- one step

        reduce    C -> 1x1: max( |u| + c_fx + |v| + c_fy ), max signal speed, floors fired
        clock     1x1: this substep's dt, the time done so far this frame, c_h
        predict   MUSCL-Hancock's half step, all four textures: U* = U - dt/2 div F( U +- dU/2 )
        flux x/y  every face's flux ONCE: U* reconstructed with MC-limited
                  slopes, GLM, HLLD, resistivity; the walls
        update    the flux differences, the sources, the dual-energy switch,
                  cooling, the floors

    The faces are their own passes for two reasons. Each is computed once and
    read by both its cells, so what one cell loses its neighbour gains bit for
    bit. And a single update doing four HLLD solves per cell was one long,
    register-heavy fragment the GPU could not keep enough of in flight: at the
    default grid a substep cost the latency of that shader, not its work.

    dt lives on the GPU. The CPU plans a number of substeps from last frame's
    speeds; each substep takes min( the CFL step, what is left of the frame /
    the substeps left ). If the plan was short the frame covers less time --
    simulated time runs slow, never unstable.

    ------------------------------------------------------------------ the rest

        ignite, sources   lay down a ball; feed the clip, heat, the pellet
        emission          rho^2 sqrt(T) times the tracer (or a debug view), on the grid
        glow              box reduce to <= 256 cells, then three cascaded
                          mirror-bounded Gaussians: it moves light, never makes any
        poisson           A_z, del^2 A = -J_z, by multigrid V-cycles, for Field Lines
        composite         onto the host raster: exposure, the shoulder, lines, mix

    Every program is assembled from the pieces in Shaders.cpp. `SourceFor`
    is what the plugin compiles and what `cttest --dump-shaders` writes out,
    so glslc checks the text that actually runs.
*/

#include <string>
#include <vector>

namespace containment
{

enum class Program
{
	Reduce = 0,
	Clock,
	Predict,
	Flux,
	Update,
	Ignite,
	Sources,
	Emission,
	Downsample,
	Blur,
	Residual,
	Smooth,
	Prolong,
	Composite,
	Count
};

struct ProgramSource
{
	const char* name;
	std::string vertex;
	std::string fragment;
};

/// The full source of one program, as compiled.
ProgramSource SourceFor( Program program );

/// Every program, in enum order.
std::vector< ProgramSource > AllSources();

/// The glow's Gaussians, in frame heights, and their share of the glare.
constexpr int kGlowStages              = 3;
constexpr float kGlowSigma[ 3 ]        = { 0.012f, 0.045f, 0.15f };
constexpr float kGlowShare[ 3 ]        = { 0.5f, 0.33f, 0.17f };
/// The glow runs on a copy of the emission reduced to at most this many cells
/// on the short side.
constexpr int kGlowCells = 256;
/// The widest blur, in taps to one side. Past this the Gaussian is truncated
/// -- never reached at kGlowCells with the sigmas above (3 sigma < 116).
constexpr int kMaxBlurTaps = 128;

} // namespace containment
