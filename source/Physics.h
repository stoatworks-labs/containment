#pragma once

#include <cstdint>
#include <vector>

/**
    The CPU side of the model: the coils, the stirring force's modes, and the
    grid. Everything the GPU is handed as uniforms is computed here, in double.

    The vacuum field of the coils is written twice: here (to normalise the
    currents and for the harness) and in the shaders (where every ghost cell
    and every Ignite evaluates it). Both are the closed-form field of infinite
    line currents; they are marked `//= mirrored`.
*/
namespace containment
{
constexpr int kMaxCoils = 12;

/// The coil set, in frame-height units, relative to the grid's origin.
struct Coils
{
	int count = 0;
	double x[ kMaxCoils ]       = {};
	double y[ kMaxCoils ]       = {};
	/// Signed current over 2 pi (mu0 = 1), so B = sum I' ( -(y - yk), x - xk ) / d^2.
	double current[ kMaxCoils ] = {};
	/// The uniform guide field along the axis.
	double guide = 0.0;
};

/**
    The bottle.

    `poles` line currents of alternating sign on a circle of `radius` about
    (cx, cy), the first at `angle`. Normalised so that |B| at radius 0.5 in the
    first gap between coils (the first cusp) is `field` -- before `strength`,
    which scales the whole set (1 in operation, falling to 0 in a quench).
    The guide field is `guide` (already in absolute units) times `strength`.
*/
Coils MakeCoils( int poles, double radius, double angle, double field, double guide, double strength, double cx,
                 double cy );

//= mirrored in Shaders.cpp, vacuumField()
/// The coils' field at (x, y): ( Bx, By, Bz ).
void VacuumField( const Coils& coils, double x, double y, double& bx, double& by, double& bz );

//= mirrored in Shaders.cpp, vacuumPotential()
/// The flux function A_z of the coils' in-plane field: B_perp = grad A_z x z.
double VacuumPotential( const Coils& coils, double x, double y );

/// The grid: square cells, `cells` on the frame's short side, the frame's long
/// side rounded to a multiple of 8 so the glow's three 2x reductions are
/// exact -- and, for Open, a margin of `margin` cells all round the frame that
/// is simulated but never shown (it holds the absorbing layer).
struct Grid
{
	int nx    = 0;  ///< simulated cells, margin included
	int ny    = 0;
	int fx    = 0;  ///< the frame's own cells
	int fy    = 0;
	int ox    = 0;  ///< where the frame starts, in cells (the margin)
	int oy    = 0;
	double dx = 0.0;///< cell size, frame heights
	double lx = 0.0;///< nx * dx
	double ly = 0.0;///< ny * dx: 1 on a landscape frame with no margin
};

Grid ChooseGrid( int width, int height, int cells, int margin = 0 );

/// A 32-bit PCG generator (O'Neill), for the drive's random walk and nothing
/// in a shader.
class Pcg
{
public:
	explicit Pcg( uint64_t seed = 0x853c49e6748fea9bULL );
	uint32_t Next();
	double Uniform();///< [0, 1)
	double Normal(); ///< Box-Muller, mean 0, variance 1

private:
	uint64_t state;
	uint64_t increment = 0xda3e39cb94b95bdbULL;
};

/**
    The stirring force: a sum of Fourier modes of a stream function, so the
    force is divergence-free by construction (f = curl( phi z )). Each mode's
    two amplitudes are Ornstein-Uhlenbeck processes with correlation time
    `kDriveTime`, so the stirring wanders rather than flickers.

    The wavevectors sit in a band around 2 pi / scale at random angles; the
    rms acceleration is `strength`.
*/
class Drive
{
public:
	struct Mode
	{
		float kx, ky, a, b;///< phi += a cos( k.x ) + b sin( k.x )
	};

	/// Advance by `dt` of simulated time, with the given scale and strength.
	void Advance( double dt, double scale, double strength );
	const std::vector< Mode >& Modes() const
	{
		return modes;
	}
	void Reset();

private:
	Pcg random;
	std::vector< Mode > modes;
	std::vector< double > state;///< 2 per mode, unit-variance OU
	std::vector< double > angle;///< the wavevectors' directions
	std::vector< double > band; ///< and their lengths, as a multiple of 2 pi / scale
};

} // namespace containment
