#include "Physics.h"

#include "Controls.h"

#include <algorithm>
#include <cmath>

namespace containment
{
namespace
{
constexpr double kPi = 3.14159265358979323846;
}

//---------------------------------------------------------------------------
void VacuumField( const Coils& coils, double x, double y, double& bx, double& by, double& bz )
{
	//= mirrored in Shaders.cpp, vacuumField()
	bx = 0.0;
	by = 0.0;
	for( int k = 0; k < coils.count; ++k )
	{
		const double dx = x - coils.x[ k ];
		const double dy = y - coils.y[ k ];
		const double d2 = dx * dx + dy * dy;
		bx += -coils.current[ k ] * dy / d2;
		by += coils.current[ k ] * dx / d2;
	}
	bz = coils.guide;
}

double VacuumPotential( const Coils& coils, double x, double y )
{
	//= mirrored in Shaders.cpp, vacuumPotential()
	//B = ( dA/dy, -dA/dx ) with A = -I' ln d, which is the field above.
	double a = 0.0;
	for( int k = 0; k < coils.count; ++k )
	{
		const double dx = x - coils.x[ k ];
		const double dy = y - coils.y[ k ];
		a -= coils.current[ k ] * 0.5 * std::log( dx * dx + dy * dy );
	}
	return a;
}

Coils MakeCoils( int poles, double radius, double angle, double field, double guide, double strength, double cx,
                 double cy )
{
	Coils coils;
	coils.count = std::clamp( poles, 0, kMaxCoils );
	coils.guide = guide * strength;
	if( coils.count == 0 )
		return coils;

	for( int k = 0; k < coils.count; ++k )
	{
		const double phi    = angle + 2.0 * kPi * k / coils.count;
		coils.x[ k ]        = cx + radius * std::cos( phi );
		coils.y[ k ]        = cy + radius * std::sin( phi );
		coils.current[ k ] = ( k % 2 == 0 ) ? 1.0 : -1.0;
	}

	//Normalise: |B| at radius 0.5, in the first gap, is `field`. A pure
	//multipole has the same |B| all round a circle; a finite coil radius
	//ripples it a little, and the gap is where the plasma leaks, so that is
	//where the confinement is quoted.
	const double gap = angle + kPi / coils.count;
	double bx = 0.0, by = 0.0, bz = 0.0;
	VacuumField( coils, cx + 0.5 * std::cos( gap ), cy + 0.5 * std::sin( gap ), bx, by, bz );
	const double here  = std::sqrt( bx * bx + by * by );
	const double scale = here > 0.0 ? field * strength / here : 0.0;
	for( int k = 0; k < coils.count; ++k )
		coils.current[ k ] *= scale;
	return coils;
}

//---------------------------------------------------------------------------
Grid ChooseGrid( int width, int height, int cells )
{
	Grid grid;
	const bool landscape = width >= height;
	const double aspect  = landscape ? static_cast< double >( width ) / std::max( height, 1 )
	                                 : static_cast< double >( height ) / std::max( width, 1 );
	const int shortSide  = cells;
	const int longSide   = std::max( 8, 8 * static_cast< int >( std::lround( shortSide * aspect / 8.0 ) ) );

	grid.nx = landscape ? longSide : shortSide;
	grid.ny = landscape ? shortSide : longSide;
	//The frame height is 1 whichever way round the frame is.
	grid.dx = 1.0 / grid.ny;
	grid.lx = grid.nx * grid.dx;
	grid.ly = 1.0;
	return grid;
}

//---------------------------------------------------------------------------
Pcg::Pcg( uint64_t seed ) : state( 0 )
{
	Next();
	state += seed;
	Next();
}

uint32_t Pcg::Next()
{
	const uint64_t old = state;
	state              = old * 6364136223846793005ULL + increment;
	const uint32_t xorshifted = static_cast< uint32_t >( ( ( old >> 18u ) ^ old ) >> 27u );
	const uint32_t rot        = static_cast< uint32_t >( old >> 59u );
	return ( xorshifted >> rot ) | ( xorshifted << ( ( -rot ) & 31 ) );
}

double Pcg::Uniform()
{
	return Next() / 4294967296.0;
}

double Pcg::Normal()
{
	const double u1 = std::max( Uniform(), 1e-12 );
	const double u2 = Uniform();
	return std::sqrt( -2.0 * std::log( u1 ) ) * std::cos( 2.0 * kPi * u2 );
}

//---------------------------------------------------------------------------
void Drive::Reset()
{
	random = Pcg( 0x5eedULL );
	modes.clear();
	state.clear();
	angle.clear();
	band.clear();
}

void Drive::Advance( double dt, double scale, double strength )
{
	if( state.empty() )
	{
		state.resize( 2 * kDriveModes );
		for( double& s : state )
			s = random.Normal();
		for( int m = 0; m < kDriveModes; ++m )
		{
			angle.push_back( 2.0 * kPi * random.Uniform() );
			band.push_back( 0.75 + 0.5 * random.Uniform() );
		}
	}

	//Ornstein-Uhlenbeck, exact for any step: x <- x e^{-dt/T} + sqrt(1 - e^{-2dt/T}) N.
	const double decay = std::exp( -std::max( dt, 0.0 ) / kDriveTime );
	const double kick  = std::sqrt( std::max( 0.0, 1.0 - decay * decay ) );
	for( double& s : state )
		s = s * decay + kick * random.Normal();

	//The force of mode m is |k| times its stream amplitude, so the rms
	//acceleration over the modes is strength when each amplitude is
	//strength / ( |k| sqrt( modes ) ): each of the 2M unit-variance
	//amplitudes contributes half its square on average (cos^2, sin^2).
	modes.resize( kDriveModes );
	for( int m = 0; m < kDriveModes; ++m )
	{
		const double k    = 2.0 * kPi / std::max( scale, 1e-3 ) * band[ m ];
		const double amp  = strength / ( k * std::sqrt( static_cast< double >( kDriveModes ) ) );
		modes[ m ].kx     = static_cast< float >( k * std::cos( angle[ m ] ) );
		modes[ m ].ky     = static_cast< float >( k * std::sin( angle[ m ] ) );
		modes[ m ].a      = static_cast< float >( amp * state[ 2 * m ] );
		modes[ m ].b      = static_cast< float >( amp * state[ 2 * m + 1 ] );
	}
}

} // namespace containment
